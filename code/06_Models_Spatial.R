
library(tidyverse)
library(here)
library(sf)
library(tigris)
library(spdep)        # poly2nb(), moran.test()

source(here("code", "inla_helpers.R"))  # INLA, shared priors, tidy_inla(), compare_cpo()

options(tigris_use_cache = TRUE)

# =============================================================================
# SPATIAL PANEL MODELS — log(HH Tax Burden) with county BYM2 errors
# =============================================================================
# County x COG-F wave panel (2012, 2017, 2022). Each model adds a BYM2 county
# random effect: a spatially structured ICAR component (a county's effect is
# conditional on its neighbours') plus an unstructured iid component, mixed by
# phi. Repeated waves per county separate the iid county effect from the
# residual sigma, which a single cross-section cannot do.
#
# Fit in INLA with the built-in bym2 latent model (scaled ICAR, sum-to-zero).
#
# Federal land share uses a single PAD-US vintage, so it is constant within
# county: its coefficient is identified from between-county variation only.
#
#   Model 0 — state FE + year FE + BYM2 (spatial baseline)
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

# INLA predicts NA responses and treats NA covariates as 0, rather than
# dropping the row, so check up front
model_vars <- c("log_HH_TaxBurden", "usfs_doi_share_pct", "metro", "unincorp_share_pct",
                "log_population", "log_n_govts", "log_med_hh_income",
                "pct_bachelors_plus", "state_name", "year")
stopifnot(all(map_lgl(data2[model_vars], \(x) all(!is.na(x) & !is.infinite(x)))))

# =============================================================================
# ADJACENCY MATRIX — rook contiguity
# =============================================================================

counties_sf <- counties(state = unique(substr(data2$GEOID, 1, 2)), cb = TRUE, year = 2013) |>
  select(GEOID) |>
  filter(GEOID %in% data2$GEOID) |>
  arrange(GEOID)

stopifnot(nrow(counties_sf) == n_distinct(data2$GEOID))

# bym2 index: row of counties_sf (and of W)
data2$county_id <- match(data2$GEOID, counties_sf$GEOID)

nb <- poly2nb(counties_sf, queen = FALSE)

# ICAR needs every county to have a neighbour. Link islands (e.g. San Juan, WA)
# to their nearest county by centroid distance.
islands <- which(card(nb) == 0)
if (length(islands) > 0) {
  nearest <- knn2nb(knearneigh(st_centroid(st_geometry(counties_sf)), k = 1))
  for (i in islands) {
    j <- nearest[[i]]
    nb[[i]] <- as.integer(sort(unique(c(nb[[i]][nb[[i]] != 0], j))))
    nb[[j]] <- as.integer(sort(unique(c(nb[[j]][nb[[j]] != 0], i))))
  }
  message("Linked islands to nearest county: ",
          paste(counties_sf$GEOID[islands], collapse = ", "))
}
stopifnot(n.comp.nb(nb)$nc == 1)  # single connected graph

W <- nb2mat(nb, style = "B")
dimnames(W) <- list(counties_sf$GEOID, counties_sf$GEOID)

# --- Motivation: spatial autocorrelation in non-spatial OLS residuals --------
# Moran's I on county-mean residuals (one value per county, ordered by GEOID)
ols2 <- lm(
  log_HH_TaxBurden ~ usfs_doi_share_pct + metro + unincorp_share_pct +
                     log_population + log_n_govts + log_med_hh_income +
                     pct_bachelors_plus + factor(state_name) + factor(year),
  data = data2
)
ols2_resid <- tibble(GEOID = data2$GEOID, resid = residuals(ols2)) |>
  group_by(GEOID) |>
  summarise(resid = mean(resid)) |>
  arrange(GEOID)

moran.test(ols2_resid$resid, nb2listw(nb, style = "W"))

# =============================================================================
# MODELS
# =============================================================================

formulas <- list(
  "M0: FE + BYM2" = log_HH_TaxBurden ~ factor(state_name) + factor(year),
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

# Formula is the fixed-effect part only; the county effect is always added.
# INLA reads only the non-zero pattern of W, so the graph is binary.
fit_bym2 <- function(formula) {
  fit_inla(
    update(formula, . ~ . + f(county_id, model = "bym2", graph = W,
                              scale.model = TRUE, constr = TRUE,
                              hyper = list(prec = pc_prec,
                                           phi  = list(prior = "pc", param = pc_phi_param)))),
    data2
  )
}

models_sp <- map(formulas, fit_bym2)

# --- Fit checks --------------------------------------------------------------
# mode.status 0 = optimiser converged cleanly; also see CPO failures below
walk(models_sp, \(m) print(summary(m)))
map_int(models_sp, \(m) m$mode$mode.status)

# --- Model comparison (LOO via CPO) ------------------------------------------
compare_cpo(models_sp)

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
  "sdcar"              = "SD (BYM2 county effect)",
  "rhocar"             = "Phi (spatial share)",
  "sigma"              = "Sigma (residual)"
)

# Term names (sdcar, rhocar) kept from the brms version so 00_run_all.R
# lines up across scripts
hyper_bym2 <- c(hyper_sigma, list(
  sdcar  = list(name = "Log precision for county_id", fun = sd_from_log_prec),
  rhocar = list(name = "Logit phi for county_id",     fun = plogis)
))

# Posterior median and 95% CrI per term (M0 has no covariates)
coefs_sp <- imap(models_sp, \(m, nm) tidy_inla(m, hyper_bym2) |> mutate(model = nm, .before = 1)) |>
  list_rbind() |>
  filter(term %in% names(coef_labels)) |>
  mutate(label = factor(coef_labels[term], levels = coef_labels),
         model = fct_inorder(model))

dir.create(here("output", "tables"),  showWarnings = FALSE)
dir.create(here("output", "figures"), showWarnings = FALSE)
write_csv(coefs_sp, here("output", "tables", "06_coefs.csv"))

# --- Results table: median [95% CrI], one column per model -------------------
results_sp <- coefs_sp |>
  arrange(label) |>
  transmute(model, term = label, est = sprintf("%.3f [%.3f, %.3f]", median, lower, upper)) |>
  pivot_wider(names_from = model, values_from = est, names_sort = TRUE) |>  # model level order
  bind_rows(tibble(term = "ELPD (LOO, CPO)",
                   !!!map(models_sp, \(m) sprintf("%.3f", elpd_cpo(m)))))
results_sp
write_html_table(results_sp, here("output", "tables", "06_results.html"),
                 caption = "Posterior medians with 95% credible intervals. State and year fixed effects included but not shown. BYM2 county random effect on rook contiguity.")

# --- Coefficient plot: posterior median and 95% CrI --------------------------
coefplot_sp <- coefs_sp |>
  filter(!term %in% c("sdcar", "rhocar", "sigma")) |>
  ggplot(aes(x = median, xmin = lower, xmax = upper, y = fct_rev(label), color = model)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(position = position_dodge(width = 0.6, orientation = "y", reverse = TRUE),
                  size = 0.3) +
  labs(title = "BYM2 Models (rook contiguity) — Coefficients with 95% CrIs",
       x = "Coefficient (log HH Tax Burden)", y = NULL, color = NULL)
coefplot_sp
ggsave(here("output", "figures", "06_coefplot.png"), coefplot_sp, width = 8, height = 5)

# --- Map of posterior mean spatial effect (Model 2) --------------------------
# bym2 latent vector is (total county effect, scaled ICAR); map the first half
J <- nrow(counties_sf)
counties_sf |>
  mutate(car_effect = models_sp[["M2: + Fed Share"]]$summary.random$county_id$mean[1:J]) |>
  ggplot() +
  geom_sf(aes(fill = car_effect), color = "white", linewidth = 0.1) +
  scale_fill_gradient2() +
  labs(title = "Model 2 — Posterior Mean BYM2 County Effect",
       fill  = "log scale")
