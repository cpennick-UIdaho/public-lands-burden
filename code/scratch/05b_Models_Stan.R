
library(tidyverse)
library(here)
library(sf)
library(tigris)
library(rstan)        # same engine as 06b_, so only the county effect differs
library(posterior)    # summarise_draws()
library(loo)

options(tigris_use_cache = TRUE, mc.cores = parallel::detectCores())

# =============================================================================
# NON-SPATIAL PANEL MODELS — Stan counterpart to 06_ / 06b_
# =============================================================================
# Same panel (all waves), formulas, priors, and likelihood as the BYM2 models,
# fit in Stan (code/stan/lm_county.stan), so estimates and LOO are directly
# comparable with the spatial fits. Each formula is fit two ways:
#
#   Pooled — no county effect; Bayesian analogue of the OLS in 05_Models.R
#            (and of the ols2 used for Moran's I in 06_/06b_)
#   IID    — iid county effect, i.e. BYM2 with rho = 0. Against the BYM2 fits
#            this isolates the gain from spatial structure; against Pooled,
#            the gain from modelling repeated observations of a county.
#
#   Model 0 — state FE + year FE
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
# MODELS
# =============================================================================

lm_county <- stan_model(here("code", "scratch", "stan", "lm_county.stan"))

fit_lm <- function(formula, has_county) {
  X <- model.matrix(formula, data2)[, -1, drop = FALSE]
  X <- scale(X, center = TRUE, scale = FALSE)  # as brms does, so Intercept is at the mean
  y <- model.response(model.frame(formula, data2))

  fit <- sampling(
    lm_county,
    data = list(
      N = nrow(X), K = ncol(X), X = X, y = y,
      J = length(geoids), county = match(data2$GEOID, geoids),
      has_county = as.integer(has_county),
      y_median = median(y)
    ),
    chains = 4, iter = 4000, seed = 2012,
    control = list(adapt_delta = 0.95)
  )
  attr(fit, "coef_names") <- colnames(X)
  attr(fit, "has_county") <- has_county
  fit
}

formulas <- list(
  "M0: FE" = log_HH_TaxBurden ~ factor(state_name) + factor(year),
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

models_pool <- map(formulas, fit_lm, has_county = FALSE)
models_iid  <- map(formulas, fit_lm, has_county = TRUE)

models_ns <- c(
  set_names(models_pool, paste(names(formulas), "(Pooled)")),
  set_names(models_iid,  paste(names(formulas), "(IID)"))
)

# --- Convergence -------------------------------------------------------------
walk(models_ns, \(m) {
  print(m, pars = c("Intercept", "b", "sd_county", "sigma"))  # Rhat < 1.01, ESS
  check_hmc_diagnostics(m)                                     # divergences, treedepth
})

# --- Model comparison (PSIS-LOO) ---------------------------------------------
# log_lik rows follow data2, as in 06_ and 06b_, so these can be passed to
# loo_compare() together with loo_cw (06b_) or map(models_sp, loo) (06_)
loo_ns <- map(models_ns, loo)  # rstan's loo() method reads log_lik
loo_compare(loo_ns)

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
  "sd_county"          = "SD (iid county effect)",
  "sigma"              = "Sigma (residual)"
)

# Posterior median and 95% CrI per term
tidy_lm <- function(fit) {
  sd    <- if (attr(fit, "has_county")) "sd_county"  # zero-length when pooled
  draws <- as_draws_df(as.array(fit, pars = c("b", sd, "sigma")))
  summarise_draws(draws, median, ~quantile(.x, c(0.025, 0.975))) |>
    transmute(term = c(attr(fit, "coef_names"), sd, "sigma"),
              median, lower = `2.5%`, upper = `97.5%`) |>
    filter(term %in% names(coef_labels))
}

coefs_ns <- imap(models_ns, \(m, nm) tidy_lm(m) |> mutate(model = nm, .before = 1)) |>
  list_rbind() |>
  mutate(label = factor(coef_labels[term], levels = coef_labels),
         model = fct_inorder(model))

dir.create(here("output", "tables"),  showWarnings = FALSE)
dir.create(here("output", "figures"), showWarnings = FALSE)
write_csv(coefs_ns, here("output", "tables", "05b_coefs.csv"))

# --- Results table: median [95% CrI], one column per model -------------------
results_ns <- coefs_ns |>
  arrange(label) |>
  transmute(model, term = label, est = sprintf("%.3f [%.3f, %.3f]", median, lower, upper)) |>
  pivot_wider(names_from = model, values_from = est, names_sort = TRUE) |>  # model level order
  bind_rows(tibble(term = "ELPD (LOO)",
                   !!!map(loo_ns, \(l) sprintf("%.1f", l$estimates["elpd_loo", "Estimate"]))))
# Posterior medians with 95% credible intervals. State and year fixed effects
# included but not shown. IID = non-spatial county random effect.
results_ns
write_csv(results_ns, here("output", "tables", "05b_results.csv"))

# --- Coefficient plot: posterior median and 95% CrI --------------------------
coefplot_ns <- coefs_ns |>
  filter(!term %in% c("sd_county", "sigma")) |>
  ggplot(aes(x = median, xmin = lower, xmax = upper, y = fct_rev(label), color = model)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(position = position_dodge(width = 0.6, orientation = "y", reverse = TRUE),
                  size = 0.3) +
  labs(title = "Non-Spatial Models — Coefficients with 95% CrIs",
       x = "Coefficient (log HH Tax Burden)", y = NULL, color = NULL)
coefplot_ns
ggsave(here("output", "figures", "05b_coefplot.png"), coefplot_ns, width = 8, height = 5)

# --- Map of posterior mean iid county effect (Model 2) -----------------------
# Compare with the BYM2 maps in 06_ / 06b_: spatial clustering left here is
# what the ICAR component picks up
counties(state = unique(substr(geoids, 1, 2)), cb = TRUE, year = 2013) |>
  select(GEOID) |>
  filter(GEOID %in% geoids) |>
  arrange(GEOID) |>
  mutate(county_effect = colMeans(as.matrix(models_iid[["M2: + Fed Share"]], pars = "r"))) |>
  ggplot() +
  geom_sf(aes(fill = county_effect), color = "white", linewidth = 0.1) +
  scale_fill_gradient2() +
  labs(title = "Model 2 — Posterior Mean IID County Effect (non-spatial)",
       fill  = "log scale")
