
library(tidyverse)
library(here)
library(readxl)
library(sf)
library(tigris)
library(spdep)        # mat2listw(), moran.test(), n.comp.nb()
library(rstan)        # brms car() only takes binary adjacency, so fit directly
library(posterior)    # summarise_draws()
library(loo)

options(tigris_use_cache = TRUE, mc.cores = parallel::detectCores())

# =============================================================================
# SPATIAL PANEL MODELS — commuting-weighted BYM2 (alternative to 06_)
# =============================================================================
# Same panel, formulas, and priors as 06_Models_Spatial.R, but the ICAR
# neighbourhood is functional rather than physical: counties are linked by
# ACS county-to-county commuting flows, weighted by flow intensity, so a
# county's effect is pulled toward the counties it exchanges workers with.
#
#   w_ij = (flow_ij / workers_i + flow_ji / workers_j) / 2
#
# i.e. the mean of the two directional out-commuting shares, which keeps W
# symmetric (required by the ICAR) and scale-free (LA-Orange is not 1000x a
# rural pair just because of population). Links with w_ij < 0.5% are dropped.
#
# Each wave uses the latest flows vintage that ends before it, so the graph
# changes by wave. The spatial (ICAR) component is therefore wave-specific,
# while the iid component stays per county as in 06_.
#
# brms car() converts any non-zero weight to 1, so the weighted BYM2 is written
# in Stan (code/stan/bym2_weighted.stan) with brms's default priors.
#
#   Model 0 — state FE + year FE + BYM2 (spatial baseline)
#   Model 1 — adds controls, excluding federal land share
#   Model 2 — adds USFS + DOI land share
# =============================================================================

data2 <- readRDS(here("output", "Cty_Agg.rds")) |>
  mutate(
    log_HH_TaxBurden  = log(HH_TaxBurden),
    log_GovtsPer10k   = log(GovtsPer10k),
    log_population    = log(population),
    log_med_hh_income = log(med_hh_income)
  )

geoids <- sort(unique(data2$GEOID))
waves  <- sort(unique(data2$year))

# =============================================================================
# COMMUTING WEIGHT MATRICES — one per wave
# =============================================================================

# Census releases county flows from non-overlapping 5-year ACS (2006-10, 2011-15,
# 2016-20; plus 2009-13, which ends after 2012 and is superseded by 2011-15 for
# 2017). Table 1 is residence county -> workplace county; column positions
# (res state, res county, work state, work county, workers) differ by vintage.
flow_vintages <- tribble(
  ~wave, ~vintage,    ~url,                                                                                                ~cols,
  2012,  "2006_2010", "https://www2.census.gov/programs-surveys/demo/tables/metro-micro/2010/commuting-employment-2010/table1.xlsx", c(1, 2, 3, 4, 5),
  2017,  "2011_2015", "https://www2.census.gov/programs-surveys/demo/tables/metro-micro/2015/commuting-flows-2015/table1.xlsx",      c(1, 2, 5, 6, 9),
  2022,  "2016_2020", "https://www2.census.gov/programs-surveys/demo/tables/metro-micro/2020/commuting-flows-2020/table1.xlsx",      c(1, 2, 5, 6, 9)
)
stopifnot(setequal(flow_vintages$wave, waves))

read_flows <- function(vintage, url, cols) {
  path <- here("input", paste0("commuting_flows_", vintage, ".xlsx"))
  if (!file.exists(path)) download.file(url, destfile = path, mode = "wb")

  # Workplace state FIPS is 3 digits in some vintages ("006"); header and
  # footnote rows drop out on the 5-digit residence code and non-numeric flow
  x <- read_excel(path, col_names = FALSE, col_types = "text")
  tibble(
    res     = paste0(x[[cols[1]]], x[[cols[2]]]),
    work    = paste0(str_sub(x[[cols[3]]], -2), x[[cols[4]]]),
    workers = suppressWarnings(as.numeric(x[[cols[5]]]))
  ) |>
    filter(str_detect(res, "^\\d{5}$"), !is.na(workers))
}

# Drop weak, noisy links. At 0.5% every wave's graph stays connected (968-999
# edges, min degree 1); at 1% the 2016-20 graph splits into 10 components, so
# the ICAR would need per-component sum-to-zero constraints and scaling
min_weight <- 0.005

commute_W <- function(flows) {
  stopifnot(all(geoids %in% flows$res))

  # Denominator is all resident workers, including those working outside the sample
  res_workers <- flows |> count(res, wt = workers, name = "res_total")

  share <- flows |>
    filter(res != work, res %in% geoids, work %in% geoids) |>
    left_join(res_workers, by = "res") |>
    mutate(share = workers / res_total)

  S <- matrix(0, length(geoids), length(geoids), dimnames = list(geoids, geoids))
  S[cbind(share$res, share$work)] <- share$share
  W <- (S + t(S)) / 2
  W[W < min_weight] <- 0

  nb <- mat2listw(W, style = "B", zero.policy = TRUE)$neighbours
  stopifnot(min(card(nb)) > 0, n.comp.nb(nb)$nc == 1)  # no islands, one graph
  W
}

W_list <- flow_vintages |>
  pmap(\(wave, vintage, url, cols) commute_W(read_flows(vintage, url, cols))) |>
  set_names(flow_vintages$wave)
W_list <- W_list[as.character(waves)]

# BYM2 scaling: geometric mean of the marginal variances of the (weighted)
# ICAR, i.e. diag of the generalised inverse of Q = D - W (Riebler et al. 2016)
scaling_factor <- map_dbl(W_list, \(W) {
  Q <- diag(rowSums(W)) - W
  exp(mean(log(diag(MASS::ginv(Q)))))
})

# Edges of all wave graphs, indexed into the stacked county-wave vector
# ((t - 1) * J + j), as the Stan model expects
J <- length(geoids)
edges <- imap(unname(W_list), \(W, t) {
  e <- which(upper.tri(W) & W > 0, arr.ind = TRUE)
  tibble(node1 = (t - 1) * J + e[, 1], node2 = (t - 1) * J + e[, 2], weight = W[e])
}) |>
  bind_rows()

# --- Motivation: spatial autocorrelation in non-spatial OLS residuals --------
# Moran's I on each wave's residuals under that wave's commuting weights
ols2 <- lm(
  log_HH_TaxBurden ~ usfs_doi_share_pct + metro + unincorp_share_pct +
                     log_population + log_GovtsPer10k + log_med_hh_income +
                     pct_bachelors_plus + factor(state_name) + factor(year),
  data = data2
)
ols2_resid <- tibble(GEOID = data2$GEOID, year = data2$year, resid = residuals(ols2))

iwalk(W_list, \(W, wave) {
  r <- ols2_resid |> filter(year == as.numeric(wave))
  stopifnot(identical(sort(r$GEOID), geoids))  # balanced panel: every county each wave
  print(moran.test(r$resid[match(geoids, r$GEOID)], mat2listw(W, style = "W")))
})

# =============================================================================
# MODELS
# =============================================================================

bym2_weighted <- stan_model(here("code", "stan", "bym2_weighted.stan"))

# Formula is the fixed-effect part only; the county effect is always added
fit_bym2_w <- function(formula) {
  X <- model.matrix(formula, data2)[, -1, drop = FALSE]
  X <- scale(X, center = TRUE, scale = FALSE)  # as brms does, so Intercept is at the mean
  y <- model.response(model.frame(formula, data2))

  fit <- sampling(
    bym2_weighted,
    data = list(
      N = nrow(X), K = ncol(X), X = X, y = y,
      J = J, T = length(waves),
      cw = (match(data2$year, waves) - 1) * J + match(data2$GEOID, geoids),
      N_edges = nrow(edges), node1 = edges$node1, node2 = edges$node2,
      weight = edges$weight, scaling_factor = as.array(scaling_factor),
      y_median = median(y)
    ),
    chains = 4, iter = 4000, seed = 2012,
    control = list(adapt_delta = 0.95)
  )
  attr(fit, "coef_names") <- colnames(X)
  fit
}

model0_cw <- fit_bym2_w(
  log_HH_TaxBurden ~ factor(state_name) + factor(year)
)

model1_cw <- fit_bym2_w(
  log_HH_TaxBurden ~ metro + unincorp_share_pct + log_population + log_GovtsPer10k +
                     log_med_hh_income + pct_bachelors_plus +
                     factor(state_name) + factor(year)
)

model2_cw <- fit_bym2_w(
  log_HH_TaxBurden ~ usfs_doi_share_pct + metro + unincorp_share_pct +
                     log_population + log_GovtsPer10k +
                     log_med_hh_income + pct_bachelors_plus +
                     factor(state_name) + factor(year)
)

models_cw <- list(
  "M0: FE + BYM2"   = model0_cw,
  "M1: + Controls"  = model1_cw,
  "M2: + Fed Share" = model2_cw
)

# --- Convergence -------------------------------------------------------------
walk(models_cw, \(m) {
  print(m, pars = c("Intercept", "b", "sdcar", "rhocar", "sigma"))  # Rhat < 1.01, ESS
  check_hmc_diagnostics(m)                                          # divergences, treedepth
})

# --- Model comparison (PSIS-LOO) ---------------------------------------------
loo_cw <- map(models_cw, loo)  # rstan's loo() method reads log_lik
loo_compare(loo_cw)

# =============================================================================
# RESULTS
# =============================================================================

coef_labels <- c(
  "usfs_doi_share_pct" = "USFS+DOI Share (%)",
  "metroNonMMetro"     = "Non-Metro",
  "unincorp_share_pct" = "Unincorporated Share (%)",
  "log_population"     = "log(Population)",
  "log_GovtsPer10k"    = "log(Govts per 10k)",
  "log_med_hh_income"  = "log(Median HH Income)",
  "pct_bachelors_plus" = "Bachelor's Degree or Higher (%)",
  "sdcar"              = "SD (BYM2 county effect)",
  "rhocar"             = "Rho (spatial share)",
  "sigma"              = "Sigma (residual)"
)

# Posterior median [95% CrI] per term, one column per model
tidy_bym2_w <- function(fit) {
  draws <- as_draws_df(as.array(fit, pars = c("b", "sdcar", "rhocar", "sigma")))
  summarise_draws(draws, median, ~quantile(.x, c(0.025, 0.975))) |>
    mutate(term = c(attr(fit, "coef_names"), "sdcar", "rhocar", "sigma")) |>
    filter(term %in% names(coef_labels)) |>
    transmute(term, est = sprintf("%.3f [%.3f, %.3f]", median, `2.5%`, `97.5%`))
}

imap(models_cw, \(m, nm) tidy_bym2_w(m) |> rename(!!nm := est)) |>
  reduce(full_join, by = "term") |>
  mutate(term = factor(term, levels = names(coef_labels), labels = coef_labels)) |>
  arrange(term) |>
  bind_rows(tibble(term = "ELPD (LOO)",
                   !!!map(loo_cw, \(l) sprintf("%.1f", l$estimates["elpd_loo", "Estimate"]))))
# Posterior medians with 95% credible intervals. State and year fixed effects
# included but not shown. BYM2 county effect on ACS commuting weights (2006-10, 2011-15, 2016-20 for the
# 2012, 2017, 2022 waves).

# --- Map of posterior mean county effect by wave (Model 2) -------------------
car_effect <- tibble(
  GEOID      = rep(geoids, times = length(waves)),
  year       = rep(waves, each = J),
  car_effect = colMeans(as.matrix(model2_cw, pars = "rcar"))
)

counties(state = unique(substr(geoids, 1, 2)), cb = TRUE, year = 2013) |>
  select(GEOID) |>
  filter(GEOID %in% geoids) |>
  inner_join(car_effect, by = "GEOID") |>
  ggplot() +
  geom_sf(aes(fill = car_effect), color = "white", linewidth = 0.1) +
  scale_fill_gradient2() +
  facet_wrap(~year) +
  labs(title = "Model 2 — Posterior Mean BYM2 County Effect (commuting weights)",
       fill  = "log scale")
