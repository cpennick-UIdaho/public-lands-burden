
library(tidyverse)
library(here)
library(sf)
library(tigris)

source(here("code", "inla_helpers.R"))  # INLA, shared priors, tidy_inla(), compare_cpo()

options(tigris_use_cache = TRUE)

# =============================================================================
# NON-SPATIAL PANEL MODELS — INLA counterpart to 06_ / 06b_
# =============================================================================
# Same panel (all waves), formulas, priors, and likelihood as the BYM2 models,
# fit in INLA (code/inla_helpers.R), so estimates and LOO are directly
# comparable with the spatial fits. Each formula is fit two ways:
#
#   Pooled — no county effect; Bayesian analogue of the ols2 used for
#            Moran's I in 06_/06b_ (not of 05_Models.R, which is a 2012-only
#            cross-section with a different control set)
#   IID    — iid county effect, i.e. BYM2 with phi = 0. Against the BYM2 fits
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
    log_n_govts       = log(n_govts),
    log_population    = log(population),
    log_med_hh_income = log(med_hh_income)
  )

geoids <- sort(unique(data2$GEOID))
data2$county_id <- match(data2$GEOID, geoids)

# INLA predicts NA responses and treats NA covariates as 0, rather than
# dropping the row, so check up front
model_vars <- c("log_HH_TaxBurden", "usfs_doi_share_pct", "metro", "unincorp_share_pct",
                "log_population", "log_n_govts", "log_med_hh_income",
                "pct_bachelors_plus", "state_name", "year")
stopifnot(all(map_lgl(data2[model_vars], \(x) all(!is.na(x) & !is.infinite(x)))))

# =============================================================================
# MODELS
# =============================================================================

formulas <- list(
  "M0: FE" = log_HH_TaxBurden ~ factor(state_name) + factor(year),
  "M1: + Controls" =
    log_HH_TaxBurden ~ metro + unincorp_share_pct + log_population + log_n_govts +
                       log_med_hh_income + pct_bachelors_plus +
                       factor(state_name) + factor(year),
  "M2: + Fed Share" =
    log_HH_TaxBurden ~ usfs_doi_share_pct + metro + unincorp_share_pct +
                       log_population + log_n_govts +
                       log_med_hh_income + pct_bachelors_plus +
                       factor(state_name) + factor(year)
)

fit_iid <- function(formula) {
  fit_inla(update(formula, . ~ . + f(county_id, model = "iid",
                                     hyper = list(prec = pc_prec))),
           data2)
}

models_pool <- map(formulas, fit_inla, data = data2)
models_iid  <- map(formulas, fit_iid)

models_ns <- c(
  set_names(models_pool, paste(names(formulas), "(Pooled)")),
  set_names(models_iid,  paste(names(formulas), "(IID)"))
)

# --- Fit checks --------------------------------------------------------------
# mode.status 0 = optimiser converged cleanly; also see CPO failures below
walk(models_ns, \(m) print(summary(m)))
map_int(models_ns, \(m) m$mode$mode.status)

# --- Model comparison (LOO via CPO) ------------------------------------------
# CPO rows follow data2, as in 06_ and 06b_, so these fits can be passed to
# compare_cpo() together with models_sp (06_) or models_cw (06b_)
compare_cpo(models_ns)

# =============================================================================
# RESULTS
# =============================================================================

coef_labels <- c(
  "usfs_doi_share_pct" = "USFS+DOI Share (%)",
  "metroNon-Metro"     = "Non-Metro",
  "unincorp_share_pct" = "Unincorporated Share (%)",
  "log_population"     = "log(Population)",
  "log_n_govts"        = "log(# Govts)",
  "log_med_hh_income"  = "log(Median HH Income)",
  "pct_bachelors_plus" = "Bachelor's Degree or Higher (%)",
  "sd_county"          = "SD (iid county effect)",
  "sigma"              = "Sigma (residual)"
)

hyper_iid <- c(hyper_sigma,
               list(sd_county = list(name = "Log precision for county_id",
                                     fun  = sd_from_log_prec)))

# Posterior median and 95% CrI per term
tidy_ns <- function(fit) {
  hyper <- if (is.null(fit$summary.random$county_id)) hyper_sigma else hyper_iid
  tidy_inla(fit, hyper) |> filter(term %in% names(coef_labels))
}

coefs_ns <- imap(models_ns, \(m, nm) tidy_ns(m) |> mutate(model = nm, .before = 1)) |>
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
  bind_rows(tibble(term = "ELPD (LOO, CPO)",
                   !!!map(models_ns, \(m) sprintf("%.3f", elpd_cpo(m)))))
results_ns
write_html_table(results_ns, here("output", "tables", "05b_results.html"),
                 caption = "Posterior medians with 95% credible intervals. State and year fixed effects included but not shown. IID = non-spatial county random effect.")

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
  mutate(county_effect = models_iid[["M2: + Fed Share"]]$summary.random$county_id$mean) |>
  ggplot() +
  geom_sf(aes(fill = county_effect), color = "white", linewidth = 0.1) +
  scale_fill_gradient2() +
  labs(title = "Model 2 — Posterior Mean IID County Effect (non-spatial)",
       fill  = "log scale")
