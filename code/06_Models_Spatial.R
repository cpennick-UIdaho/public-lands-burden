
library(tidyverse)
library(here)
library(sf)
library(tigris)
library(spdep)        # poly2nb(), moran.test()
library(brms)         # Bayesian regression with car() spatial terms
library(modelsummary)

options(tigris_use_cache = TRUE, mc.cores = parallel::detectCores())

# =============================================================================
# SPATIAL PANEL MODELS — log(HH Tax Burden) with county BYM2 errors
# =============================================================================
# County x COG-F wave panel (2012, 2017, 2022). Each model adds a BYM2 county
# random effect: a spatially structured ICAR component (a county's effect is
# conditional on its neighbours') plus an unstructured iid component, mixed by
# rho. Repeated waves per county separate the iid county effect from the
# residual sigma, which a single cross-section cannot do.
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
    log_GovtsPer10k   = log(GovtsPer10k),
    log_population    = log(population),
    log_med_hh_income = log(med_hh_income)
  )

# =============================================================================
# ADJACENCY MATRIX — rook contiguity
# =============================================================================

counties_sf <- counties(state = unique(substr(data2$GEOID, 1, 2)), cb = TRUE, year = 2013) |>
  select(GEOID) |>
  filter(GEOID %in% data2$GEOID) |>
  arrange(GEOID)

stopifnot(nrow(counties_sf) == n_distinct(data2$GEOID))

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
                     log_population + log_GovtsPer10k + log_med_hh_income +
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

fit_bym2 <- function(formula) {
  brm(
    formula,
    data    = data2,
    data2   = list(W = W),
    family  = gaussian(),
    prior   = prior(normal(0, 1), class = "b"),
    chains  = 4, iter = 4000, seed = 2012,
    control = list(adapt_delta = 0.95)
  )
}

model0_sp <- fit_bym2(
  log_HH_TaxBurden ~ factor(state_name) + factor(year) +
                     car(W, gr = GEOID, type = "bym2")
)

model1_sp <- fit_bym2(
  log_HH_TaxBurden ~ metro + unincorp_share_pct + log_population + log_GovtsPer10k +
                     log_med_hh_income + pct_bachelors_plus +
                     factor(state_name) + factor(year) +
                     car(W, gr = GEOID, type = "bym2")
)

model2_sp <- fit_bym2(
  log_HH_TaxBurden ~ usfs_doi_share_pct + metro + unincorp_share_pct +
                     log_population + log_GovtsPer10k +
                     log_med_hh_income + pct_bachelors_plus +
                     factor(state_name) + factor(year) +
                     car(W, gr = GEOID, type = "bym2")
)

models_sp <- list(
  "M0: FE + BYM2"   = model0_sp,
  "M1: + Controls"  = model1_sp,
  "M2: + Fed Share" = model2_sp
)

# --- Convergence -------------------------------------------------------------
walk(models_sp, \(m) print(summary(m)))  # check Rhat < 1.01, ESS, divergences

# --- Model comparison (PSIS-LOO) ---------------------------------------------
loo_compare(map(models_sp, loo))

# =============================================================================
# RESULTS
# =============================================================================

coef_labels <- c(
  "b_usfs_doi_share_pct" = "USFS+DOI Share (%)",
  "b_metroNonMMetro"     = "Non-Metro",
  "b_unincorp_share_pct" = "Unincorporated Share (%)",
  "b_log_population"     = "log(Population)",
  "b_log_GovtsPer10k"    = "log(Govts per 10k)",
  "b_log_med_hh_income"  = "log(Median HH Income)",
  "b_pct_bachelors_plus" = "Bachelor's Degree or Higher (%)",
  "sdcar"                = "SD (BYM2 county effect)",
  "rhocar"               = "Rho (spatial share)",
  "sigma"                = "Sigma (residual)"
)

modelsummary(
  models_sp,
  statistic   = "conf.int",
  coef_map    = coef_labels,
  gof_map     = c("nobs", "elpd", "r.squared"),
  notes       = "Posterior medians with 95% credible intervals. State and year fixed effects included but not shown. BYM2 county random effect on rook contiguity.",
  title       = "Federal Land Share and log(HH Tax Burden) with County BYM2 Errors — Western US Counties, 2012–2022"
)

# --- Map of posterior mean spatial effect (Model 2) --------------------------
counties_sf |>
  mutate(car_effect = colMeans(as_draws_matrix(model2_sp, variable = "rcar"))) |>
  ggplot() +
  geom_sf(aes(fill = car_effect), color = "white", linewidth = 0.1) +
  scale_fill_gradient2() +
  labs(title = "Model 2 — Posterior Mean BYM2 County Effect",
       fill  = "log scale")
