
library(tidyverse)
library(here)
library(readxl)
library(sf)
library(tigris)
library(spdep)        # mat2listw(), moran.test(), n.comp.nb()
library(Matrix)

source(here("code", "inla_helpers.R"))  # INLA, shared priors, tidy_inla(), compare_cpo()

options(tigris_use_cache = TRUE)

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
# INLA's bym2 model uses only the non-zero pattern of the graph, so the
# weighted BYM2 is written as an rgeneric model below, with the same
# parameterisation (Riebler et al. 2016) and priors as bym2 in 06_.
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
data2$county_id <- match(data2$GEOID, geoids)

# INLA predicts NA responses and treats NA covariates as 0, rather than
# dropping the row, so check up front
model_vars <- c("log_HH_TaxBurden", "usfs_doi_share_pct", "metro", "unincorp_share_pct",
                "log_population", "log_GovtsPer10k", "log_med_hh_income",
                "pct_bachelors_plus", "state_name", "year")
stopifnot(all(map_lgl(data2[model_vars], \(x) all(!is.na(x) & !is.infinite(x)))))

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
# WEIGHTED BYM2 (rgeneric)
# =============================================================================
# Latent vector x = (v, u), as in INLA's bym2: v (indexed by the data) is the
# total county effect and u the scaled ICAR component, with
#
#   v | u ~ N( sqrt(phi / tau) u, (1 - phi) / tau I ),   u ~ ICAR(R)
#
# where R = D - W scaled so the geometric mean of the ICAR marginal variances
# is 1 (Riebler et al. 2016). Theta = (log tau, logit phi).

J  <- length(geoids)
Qw <- Matrix(diag(rowSums(W)) - W, sparse = TRUE)                # weighted ICAR structure
R  <- inla.scale.model(Qw, constr = list(A = matrix(1, 1, J), e = 0))
log_det_R <- sum(log(head(eigen(as.matrix(R), symmetric = TRUE, only.values = TRUE)$values, J - 1)))

# PC prior on logit phi for this graph, as bym2 computes it for a binary one.
# Table of (logit phi, log density) from INLA's internal helper.
phi_tab <- INLA:::inla.pc.bym.phi(Q = Qw, rankdef = 1, u = pc_phi_param[1],
                                  alpha = pc_phi_param[2], return.as.table = TRUE) |>
  str_remove("^table:") |> str_squish() |> str_split_1(" ") |> as.numeric() |>
  matrix(ncol = 2)

bym2_weighted <- function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                                  "log.prior", "quit"),
                          theta = NULL) {
  tau <- \() exp(theta[1L])
  phi <- \() 1 / (1 + exp(-theta[2L]))

  graph <- function() Q()
  Q <- function() {
    k <- tau(); p <- phi()
    I <- Matrix::Diagonal(J)
    a <- k / (1 - p)
    b <- -sqrt(p * k) / (1 - p)
    rbind(cbind(a * I, b * I),
          cbind(b * I, R + p / (1 - p) * I))
  }
  mu <- function() numeric(0)
  # v | u is proper (J dims); u is intrinsic with rank J - 1
  log.norm.const <- function() {
    -0.5 * (2 * J - 1) * log(2 * pi) + 0.5 * (J * log(tau() / (1 - phi())) + log_det_R)
  }
  # pc.prec(u, alpha) on tau, as a density for log tau; tabulated PC prior on logit phi
  log.prior <- function() {
    lambda <- -log(prec_param[2]) / prec_param[1]
    log(lambda / 2) - lambda * exp(-theta[1L] / 2) - theta[1L] / 2 +
      approx(phi_tab[, 1], phi_tab[, 2], theta[2L], rule = 2)$y
  }
  initial <- function() c(4, -3)  # bym2 defaults
  quit <- function() invisible()

  if (is.null(theta)) theta <- initial()
  do.call(match.arg(cmd), args = list())
}

bym2_w_model <- inla.rgeneric.define(
  bym2_weighted, J = J, R = R, log_det_R = log_det_R,
  phi_tab = phi_tab, prec_param = pc_prec$param
)

# =============================================================================
# MODELS
# =============================================================================

formulas <- list(
  "M0: FE + BYM2" = log_HH_TaxBurden ~ factor(state_name) + factor(year),
  "M1: + Controls" =
    log_HH_TaxBurden ~ metro + unincorp_share_pct + log_population + log_GovtsPer10k +
                       log_med_hh_income + pct_bachelors_plus +
                       factor(state_name) + factor(year),
  "M2: + Fed Share" =
    log_HH_TaxBurden ~ usfs_doi_share_pct + metro + unincorp_share_pct +
                       log_population + log_GovtsPer10k +
                       log_med_hh_income + pct_bachelors_plus +
                       factor(state_name) + factor(year)
)

# Formula is the fixed-effect part only; the county effect is always added.
# Sum-to-zero on u, the second half of the latent vector, as bym2 does.
fit_bym2_w <- function(formula) {
  fit_inla(
    update(formula, . ~ . + f(county_id, model = bym2_w_model, n = 2 * J,
                              extraconstr = list(A = cbind(matrix(0, 1, J), matrix(1, 1, J)),
                                                 e = 0))),
    data2
  )
}

models_cw <- map(formulas, fit_bym2_w)

# --- Fit checks --------------------------------------------------------------
# mode.status 0 = optimiser converged cleanly; also see CPO failures below
walk(models_cw, \(m) print(summary(m)))
map_int(models_cw, \(m) m$mode$mode.status)

# --- Model comparison (LOO via CPO) ------------------------------------------
compare_cpo(models_cw)

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
  "rhocar"             = "Phi (spatial share)",
  "sigma"              = "Sigma (residual)"
)

# rgeneric hyperparameters are named Theta1/Theta2 (log tau, logit phi)
hyper_bym2_w <- c(hyper_sigma, list(
  sdcar  = list(name = "Theta1 for county_id", fun = sd_from_log_prec),
  rhocar = list(name = "Theta2 for county_id", fun = plogis)
))

# Posterior median and 95% CrI per term
coefs_cw <- imap(models_cw, \(m, nm) tidy_inla(m, hyper_bym2_w) |> mutate(model = nm, .before = 1)) |>
  list_rbind() |>
  filter(term %in% names(coef_labels)) |>
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
  bind_rows(tibble(term = "ELPD (LOO, CPO)",
                   !!!map(models_cw, \(m) sprintf("%.1f", elpd_cpo(m)))))
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
# First half of the latent vector is the total county effect v
counties(state = unique(substr(geoids, 1, 2)), cb = TRUE, year = 2013) |>
  select(GEOID) |>
  filter(GEOID %in% geoids) |>
  arrange(GEOID) |>
  mutate(car_effect = models_cw[["M2: + Fed Share"]]$summary.random$county_id$mean[1:J]) |>
  ggplot() +
  geom_sf(aes(fill = car_effect), color = "white", linewidth = 0.1) +
  scale_fill_gradient2() +
  labs(title = "Model 2 — Posterior Mean BYM2 County Effect (commuting weights)",
       fill  = "log scale")
