
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
# Flows are the 2006-2010 ACS vintage, the latest ending before the first wave
# (2012), held fixed across waves like PAD-US and RUCC so each county has a
# single spatial effect.
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

# =============================================================================
# COMMUTING WEIGHT MATRIX
# =============================================================================

# Table 1: residence county -> workplace county, 5-year ACS 2006-2010
flows_path <- here("input", "commuting_flows_2006_2010.xlsx")
if (!file.exists(flows_path)) {
  download.file(
    "https://www2.census.gov/programs-surveys/demo/tables/metro-micro/2010/commuting-employment-2010/table1.xlsx",
    destfile = flows_path, mode = "wb"
  )
}

# Columns: res state, res county, work state (3 digits, "006"), work county,
# workers. Header and footnote rows drop out on the 5-digit residence code
# and non-numeric flow
flows <- read_excel(flows_path, col_names = FALSE, col_types = "text") |>
  transmute(
    res     = paste0(...1, ...2),
    work    = paste0(str_sub(...3, -2), ...4),
    workers = suppressWarnings(as.numeric(...5))
  ) |>
  filter(str_detect(res, "^\\d{5}$"), !is.na(workers))

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

# Drop weak, noisy links. At 0.5% the graph stays connected (983 edges, min
# degree 1); at 1% (on 2016-20 flows) it splits into 10 components, so the
# ICAR would need per-component sum-to-zero constraints and scaling
min_weight <- 0.005
W[W < min_weight] <- 0

nb <- mat2listw(W, style = "B", zero.policy = TRUE)$neighbours
stopifnot(min(card(nb)) > 0, n.comp.nb(nb)$nc == 1)  # no islands, one graph

# BYM2 scaling: geometric mean of the marginal variances of the (weighted)
# ICAR, i.e. diag of the generalised inverse of Q = D - W (Riebler et al. 2016)
Q <- diag(rowSums(W)) - W
scaling_factor <- exp(mean(log(diag(MASS::ginv(Q)))))

edges <- which(upper.tri(W) & W > 0, arr.ind = TRUE)

# --- Motivation: spatial autocorrelation in non-spatial OLS residuals --------
# Moran's I on county-mean residuals under the commuting weights
ols2 <- lm(
  log_HH_TaxBurden ~ usfs_doi_share_pct + metro + unincorp_share_pct +
                     log_population + log_GovtsPer10k + log_med_hh_income +
                     pct_bachelors_plus + factor(state_name) + factor(year),
  data = data2
)
ols2_resid <- tibble(GEOID = data2$GEOID, resid = residuals(ols2)) |>
  group_by(GEOID) |>
  summarise(resid = mean(resid)) |>
  arrange(GEOID)

moran.test(ols2_resid$resid, mat2listw(W, style = "W"))

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
      J = length(geoids), county = match(data2$GEOID, geoids),
      N_edges = nrow(edges), node1 = edges[, 1], node2 = edges[, 2],
      weight = W[edges], scaling_factor = scaling_factor,
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
  "metroNon-Metro"     = "Non-Metro",
  "unincorp_share_pct" = "Unincorporated Share (%)",
  "log_population"     = "log(Population)",
  "log_GovtsPer10k"    = "log(Govts per 10k)",
  "log_med_hh_income"  = "log(Median HH Income)",
  "pct_bachelors_plus" = "Bachelor's Degree or Higher (%)",
  "sdcar"              = "SD (BYM2 county effect)",
  "rhocar"             = "Rho (spatial share)",
  "sigma"              = "Sigma (residual)"
)

# Posterior median and 95% CrI per term
tidy_bym2_w <- function(fit) {
  draws <- as_draws_df(as.array(fit, pars = c("b", "sdcar", "rhocar", "sigma")))
  summarise_draws(draws, median, ~quantile(.x, c(0.025, 0.975))) |>
    transmute(term = c(attr(fit, "coef_names"), "sdcar", "rhocar", "sigma"),
              median, lower = `2.5%`, upper = `97.5%`) |>
    filter(term %in% names(coef_labels))
}

coefs_cw <- imap(models_cw, \(m, nm) tidy_bym2_w(m) |> mutate(model = nm, .before = 1)) |>
  list_rbind() |>
  mutate(label = factor(coef_labels[term], levels = coef_labels),
         model = fct_inorder(model))

dir.create(here("output", "tables"),  showWarnings = FALSE)
dir.create(here("output", "figures"), showWarnings = FALSE)
write_csv(coefs_cw, here("output", "tables", "06b_coefs.csv"))

# --- Results table: median [95% CrI], one column per model -------------------
results_cw <- coefs_cw |>
  arrange(label) |>
  transmute(model, term = label, est = sprintf("%.3f [%.3f, %.3f]", median, lower, upper)) |>
  pivot_wider(names_from = model, values_from = est, names_sort = TRUE) |>  # model level order
  bind_rows(tibble(term = "ELPD (LOO)",
                   !!!map(loo_cw, \(l) sprintf("%.1f", l$estimates["elpd_loo", "Estimate"]))))
# Posterior medians with 95% credible intervals. State and year fixed effects
# included but not shown. BYM2 county effect on ACS 2006-2010 commuting weights.
results_cw
write_csv(results_cw, here("output", "tables", "06b_results.csv"))

# --- Coefficient plot: posterior median and 95% CrI --------------------------
coefplot_cw <- coefs_cw |>
  filter(!term %in% c("sdcar", "rhocar", "sigma")) |>
  ggplot(aes(x = median, xmin = lower, xmax = upper, y = fct_rev(label), color = model)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(position = position_dodge(width = 0.6, orientation = "y", reverse = TRUE),
                  size = 0.3) +
  labs(title = "BYM2 Models (commuting weights) — Coefficients with 95% CrIs",
       x = "Coefficient (log HH Tax Burden)", y = NULL, color = NULL)
coefplot_cw
ggsave(here("output", "figures", "06b_coefplot.png"), coefplot_cw, width = 8, height = 5)

# --- Map of posterior mean spatial effect (Model 2) --------------------------
counties(state = unique(substr(geoids, 1, 2)), cb = TRUE, year = 2013) |>
  select(GEOID) |>
  filter(GEOID %in% geoids) |>
  arrange(GEOID) |>
  mutate(car_effect = colMeans(as.matrix(model2_cw, pars = "rcar"))) |>
  ggplot() +
  geom_sf(aes(fill = car_effect), color = "white", linewidth = 0.1) +
  scale_fill_gradient2() +
  labs(title = "Model 2 — Posterior Mean BYM2 County Effect (commuting weights)",
       fill  = "log scale")
