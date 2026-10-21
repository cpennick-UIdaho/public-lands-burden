
library(tidyverse)
library(broom)
library(here)
library(moments)      # skewness(), kurtosis()
library(lmtest)       # bptest()
library(sandwich)     # vcovHC() for robust SEs
library(car)          # vif()
library(modelsummary) # regression tables and coefficient plots

# Note: re-run fed_share_acs.R then Cty_Agg.R before running this script to
# ensure pct_bachelors_plus and all updated variables are present in the .rds.
data2 <- readRDS(here("output", "Cty_Agg.rds"))

# =============================================================================
# PHASE 2 — VARIABLE DIAGNOSTICS
# =============================================================================

diag_vars <- c("HH_TaxBurden", "HH_Expenditures", "usfs_doi_share_pct",
               "unincorp_share_pct", "GovtsPer10k",
               "med_hh_income", "population")

# --- Skewness / kurtosis table -----------------------------------------------
data2 |>
  select(all_of(diag_vars)) |>
  summarise(across(everything(), list(
    mean     = \(x) mean(x, na.rm = TRUE),
    sd       = \(x) sd(x, na.rm = TRUE),
    skewness = \(x) skewness(x, na.rm = TRUE),
    kurtosis = \(x) kurtosis(x, na.rm = TRUE),
    n_miss   = \(x) sum(is.na(x))
  ))) |>
  pivot_longer(everything(), names_to = c("variable", "stat"), names_sep = "_(?=[^_]+$)") |>
  pivot_wider(names_from = stat, values_from = value)

# --- Histograms --------------------------------------------------------------
data2 |>
  select(all_of(diag_vars)) |>
  pivot_longer(everything(), names_to = "variable", values_to = "value") |>
  ggplot(aes(x = value)) +
  geom_histogram(bins = 30) +
  facet_wrap(~variable, scales = "free") +
  labs(title = "Variable Distributions", x = NULL, y = "Count")

# --- Log-transformed DVs (apply if skewness > 1 in table above) -------------
# Coefficients become semi-elasticities: % change in DV per unit change in IV.
data2 <- data2 |>
  mutate(
    log_HH_TaxBurden    = log(HH_TaxBurden),
    log_HH_Expenditures = log(HH_Expenditures),
    log_GovtsPer10k     = log(GovtsPer10k),
    log_median_hh_inc   = log(med_hh_income),
    log_population      = log(population)
  )

# =============================================================================
# PHASE 3 — DESCRIPTIVE STATISTICS
# =============================================================================

# --- Summary table: overall and by Metro/Non-Metro ---------------------------
model_vars <- c("HH_TaxBurden", "HH_Expenditures", "GovtsPer10k",
                "usfs_doi_share_pct", "unincorp_share_pct",
                "med_hh_income", "n_govts", "population")

datasummary(
  All(data2 |> select(all_of(model_vars))) ~ Mean + SD + Min + Max + N,
  data = data2,
  title = "Descriptive Statistics — All Counties"
)

datasummary(
  All(data2 |> select(all_of(model_vars))) ~ metro * (Mean + SD + N),
  data = data2,
  title = "Descriptive Statistics by Metro/Non-Metro"
)

# --- Correlation matrix (continuous variables) --------------------------------
cor_vars <- c("log_HH_TaxBurden", "log_HH_Expenditures", "log_GovtsPer10k",
              "usfs_doi_share_pct", "unincorp_share_pct",
              "log_median_hh_inc", "log_population")

cor_matrix <- data2 |>
  select(all_of(cor_vars)) |>
  cor(use = "pairwise.complete.obs") |>
  round(2)

cor_matrix

# Formatted correlation table (continuous variables only; metro excluded)
datasummary_correlation(
  data2 |> select(all_of(cor_vars[cor_vars != "metro"])),
  fmt    = 2,
  title  = "Correlation Matrix — Key Continuous Variables",
  notes  = "Pairwise complete observations. Values rounded to 2 decimal places."
)

# --- Scatter plots: federal land share vs. DVs, by metro status --------------
data2 |>
  select(usfs_doi_share_pct, log_HH_TaxBurden, log_HH_Expenditures, log_GovtsPer10k, metro) |>
  pivot_longer(c(log_HH_TaxBurden, log_HH_Expenditures, log_GovtsPer10k),
               names_to = "outcome", values_to = "value") |>
  ggplot(aes(x = usfs_doi_share_pct, y = value, color = metro)) +
  geom_point(alpha = 0.4, size = 0.9) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~outcome, scales = "free_y") +
  labs(title = "USFS + DOI Land Share vs. Outcome Variables",
       x = "USFS + DOI Land Share (%)", y = NULL, color = NULL)

# =============================================================================
# PHASE 4 — REGRESSION MODELS
# =============================================================================
# All models: OLS with state fixed effects.
# Year FE omitted — single 2012 cross-section (year is constant).
# Log-transformed DVs used if skewness diagnostics (Phase 2) support it.
#
# Swap log_HH_TaxBurden / log_HH_Expenditures / log_GovtsPer10k for the
# untransformed versions if transformation is not warranted.
# =============================================================================

# Model 1 — Household Tax Burden (total effect of federal land) ---------------
model1 <- lm(
  log_HH_TaxBurden ~ usfs_doi_share_pct + metro + unincorp_share_pct +
                     log_population + log_GovtsPer10k + factor(state_name),
  data = data2
)



# =============================================================================
# PHASE 5 — ASSUMPTION CHECKING
# =============================================================================

check_assumptions <- function(model, label) {
  message("\n--- ", label, " ---")

  # Breusch-Pagan heteroskedasticity test
  bp <- bptest(model)
  message("Breusch-Pagan: chi2 = ", round(bp$statistic, 2),
          ", p = ", round(bp$p.value, 4))

  # Cook's distance — flag high-influence observations
  cd    <- cooks.distance(model)
  n     <- length(cd)
  threshold <- 4 / n
  n_inf <- sum(cd > threshold, na.rm = TRUE)
  message("Cook's D > 4/n (", round(threshold, 4), "): ", n_inf, " counties")
  if (n_inf > 0) {
    influential <- data2[which(cd > threshold), ] |>
      select(state_name, county_name) |>
      mutate(cooks_d = cd[cd > threshold])
    print(influential)
  }

  # VIF — multicollinearity
  message("VIF:")
  vif_vals <- vif(model)
  print(round(vif_vals, 2))
}

check_assumptions(model1,     "Model 1 — HH Tax Burden")


# Q-Q plots for residual normality
par(mfrow = c(2, 3))
walk2(
  list(model1),
  c("Tax Burden"),
  \(m, lbl) qqnorm(residuals(m), main = lbl)
)
par(mfrow = c(1, 1))

# Robust SE correction (HC3) — apply if Breusch-Pagan is significant
rob_se <- \(m) sqrt(diag(vcovHC(m, type = "HC3")))

# =============================================================================
# PHASE 6 — RESULTS PRESENTATION
# =============================================================================

fe_pattern <- "factor\\(state_name\\)"

coef_labels <- c(
  "usfs_doi_share_pct" = "USFS+DOI Share (%)",
  "metroNon-Metro"     = "Non-Metro",
  "unincorp_share_pct" = "Unincorporated Share (%)",
  "log_population"     = "log(Population)",
  "log_GovtsPer10k"    = "log(Govts per 10k)"
)

# --- Regression table --------------------------------------------------------
modelsummary(
  list("log(HH Tax Burden)" = model1),
  vcov        = \(m) vcovHC(m, type = "HC3"),
  coef_omit   = fe_pattern,
  coef_rename = coef_labels,
  stars       = c("*" = 0.1, "**" = 0.05, "***" = 0.01),
  gof_map     = c("nobs", "r.squared", "adj.r.squared"),
  notes       = "HC3 robust standard errors. State fixed effects included but not shown. DV is log-transformed.",
  title       = "Federal Land Share and Household Tax Burden — Western US Counties, 2012"
)

# --- Coefficient plot --------------------------------------------------------
modelplot(
  list("log(HH Tax Burden)" = model1),
  vcov        = \(m) vcovHC(m, type = "HC3"),
  coef_omit   = paste0(fe_pattern, "|Intercept"),
  coef_rename = coef_labels
) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    title = "Model 1 — Coefficient Estimates with 95% CIs",
    x     = "Coefficient (log HH Tax Burden)",
    color = NULL
  )

# --- Predicted log tax burden over federal land share, by metro status -------
pred_grid <- expand_grid(
  usfs_doi_share_pct = seq(0, 100, by = 1),
  metro              = c("Metro", "Non-Metro"),
  unincorp_share_pct = mean(data2$unincorp_share_pct, na.rm = TRUE),
  log_population     = mean(data2$log_population,     na.rm = TRUE),
  log_GovtsPer10k    = mean(data2$log_GovtsPer10k,    na.rm = TRUE),
  state_name         = "Colorado"  # mid-range reference state; illustrative only
)

pred_grid |>
  bind_cols(
    predict(model1, newdata = pred_grid, interval = "confidence") |>
      as_tibble()
  ) |>
  ggplot(aes(x = usfs_doi_share_pct, y = fit, color = metro, fill = metro)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.8) +
  labs(
    title    = "Predicted Household Tax Burden by Federal Land Share and Metro Status",
    subtitle = "Other variables held at means; Colorado state FE applied (illustrative)",
    x        = "USFS + DOI Land Share (%)",
    y        = "Predicted log(HH Tax Burden)",
    color    = NULL, fill = NULL
  )

# =============================================================================
# PHASE 7 — NOMINAL TAX BURDEN MODEL
# =============================================================================
# DV: HH_TaxBurden (untransformed, nominal dollars per household)
# Key predictors: USFS + DOI land share (%) and median household income ($).
# Model 6a: parsimonious — federal land share + income + state FEs.
# Model 6b: full controls — adds metro, unincorporated share, education.
# State FEs included throughout for cross-state comparability.
# HC3 robust SEs applied (heteroskedasticity expected with a nominal DV).
# =============================================================================

# Model 6a — Parsimonious specification ----------------------------------------
model6a_tax <- lm(
  HH_TaxBurden ~ usfs_doi_share_pct + med_hh_income + factor(state_name),
  data = data2
)

# Model 6b — Full controls (parallel to Phase 4 models) -----------------------
model6b_tax <- lm(
  HH_TaxBurden ~ usfs_doi_share_pct + med_hh_income + metro +
                 unincorp_share_pct + pct_bachelors_plus + factor(state_name),
  data = data2
)

# --- Assumption checks -------------------------------------------------------
check_assumptions(model6a_tax, "Model 6a — Nominal Tax Burden (parsimonious)")
check_assumptions(model6b_tax, "Model 6b — Nominal Tax Burden (full controls)")

# --- Results table -----------------------------------------------------------
modelsummary(
  list(
    "Nominal Tax Burden (Parsimonious)" = model6a_tax,
    "Nominal Tax Burden (Full Controls)" = model6b_tax
  ),
  vcov      = \(m) vcovHC(m, type = "HC3"),
  coef_omit = fe_pattern,
  coef_rename = c(
    "usfs_doi_share_pct" = "USFS+DOI Share (%)",
    "med_hh_income"      = "Median HH Income ($)",
    "metroNon-Metro"     = "Non-Metro",
    "unincorp_share_pct" = "Unincorporated Share (%)",
    "pct_bachelors_plus" = "Bachelor's Degree or Higher (%)"
  ),
  stars   = c("*" = 0.1, "**" = 0.05, "***" = 0.01),
  gof_map = c("nobs", "r.squared", "adj.r.squared"),
  notes   = "HC3 robust standard errors. State fixed effects included but not shown. DV is nominal dollars per household.",
  title   = "Federal Land Share, Income, and Nominal Household Tax Burden — Western US Counties, 2012"
)

# --- Coefficient plot: Model 6a vs. 6b ---------------------------------------
modelplot(
  list(
    "Parsimonious" = model6a_tax,
    "Full Controls" = model6b_tax
  ),
  vcov      = \(m) vcovHC(m, type = "HC3"),
  coef_omit = fe_pattern,
  coef_rename = c(
    "usfs_doi_share_pct" = "USFS+DOI Share (%)",
    "med_hh_income"      = "Median HH Income ($)",
    "metroNon-Metro"     = "Non-Metro",
    "unincorp_share_pct" = "Unincorporated Share (%)",
    "pct_bachelors_plus" = "Bachelor's Degree or Higher (%)"
  )
) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    title = "Nominal Tax Burden: USFS+DOI Share and Median Income Coefficients",
    x     = "Coefficient (nominal $ per household)",
    color = NULL
  )

# --- Predicted nominal tax burden over federal land share, by income quartile -
income_quartiles <- quantile(data2$med_hh_income, probs = c(0.25, 0.5, 0.75),
                             na.rm = TRUE)

pred_grid6 <- expand_grid(
  usfs_doi_share_pct = seq(0, 100, by = 1),
  med_hh_income      = income_quartiles,
  state_name         = "Colorado"  # reference state; illustrative only
)

pred_grid6 |>
  bind_cols(
    predict(model6a_tax, newdata = pred_grid6, interval = "confidence") |>
      as_tibble() |> rename(fit = fit, lwr = lwr, upr = upr)
  ) |>
  mutate(income_label = paste0("$", scales::comma(med_hh_income),
                               " (", c("Q1", "Q2", "Q3")[match(med_hh_income, income_quartiles)],
                               ")")) |>
  ggplot(aes(x = usfs_doi_share_pct, y = fit,
             color = income_label, fill = income_label)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.8) +
  scale_y_continuous(labels = scales::dollar) +
  labs(
    title    = "Predicted Nominal Tax Burden by Federal Land Share and Income Level",
    subtitle = "Other variables held at means; Colorado state FE applied (illustrative)",
    x        = "USFS + DOI Land Share (%)",
    y        = "Predicted HH Tax Burden (nominal $)",
    color    = "Median HH Income", fill = "Median HH Income"
  )

# =============================================================================
# PHASE 8 — QUADRATIC SPECIFICATION AND SENSITIVITY CHECK
# =============================================================================
# Motivation: quartile analysis (Phase 7) revealed a U-shaped relationship
# between federal land share and nominal tax burden in non-metro counties.
# Model 7 formalises this with a quadratic term interacted with metro status.
#
# usfs_doi_share_pct is mean-centred before squaring to reduce collinearity
# between the linear and quadratic terms.
#
# Sensitivity check: refit Model 6b and Model 7 after dropping the 20 counties
# with Cook's D > 4/n from Model 7 to assess influence on the quadratic shape.
# =============================================================================

# --- Mean-centre federal land share ------------------------------------------
fed_mean <- mean(data2$usfs_doi_share_pct, na.rm = TRUE)

data2 <- data2 |>
  mutate(
    usfs_doi_c  = usfs_doi_share_pct - fed_mean,
    usfs_doi_c2 = usfs_doi_c^2
  )

# Model 7 — Quadratic × Metro interaction -------------------------------------
model7_tax <- lm(
  HH_TaxBurden ~ (usfs_doi_c + usfs_doi_c2) * metro +
                  med_hh_income + unincorp_share_pct +
                  pct_bachelors_plus + factor(state_name),
  data = data2
)

# Model fit comparison: linear (6b) vs. quadratic (7)
anova(model6b_tax, model7_tax)
cat("AIC — Model 6b:", AIC(model6b_tax), " | Model 7:", AIC(model7_tax), "\n")
cat("BIC — Model 6b:", BIC(model6b_tax), " | Model 7:", BIC(model7_tax), "\n")

# --- Predicted curves by metro status ----------------------------------------
pred_grid7 <- expand_grid(
  usfs_doi_share_pct = seq(0, 100, by = 0.5),
  metro              = c("Metro", "Non-Metro"),
  med_hh_income      = mean(data2$med_hh_income,      na.rm = TRUE),
  unincorp_share_pct = mean(data2$unincorp_share_pct, na.rm = TRUE),
  pct_bachelors_plus = mean(data2$pct_bachelors_plus, na.rm = TRUE),
  state_name         = "Colorado"  # reference state; illustrative only
) |>
  mutate(usfs_doi_c  = usfs_doi_share_pct - fed_mean,
         usfs_doi_c2 = usfs_doi_c^2)

pred_grid7 |>
  bind_cols(
    predict(model7_tax, newdata = pred_grid7, interval = "confidence") |>
      as_tibble()
  ) |>
  ggplot(aes(x = usfs_doi_share_pct, y = fit, color = metro, fill = metro)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.9) +
  geom_rug(data = data2, aes(x = usfs_doi_share_pct, color = metro),
           inherit.aes = FALSE, alpha = 0.3, length = unit(0.03, "npc")) +
  scale_y_continuous(labels = scales::dollar) +
  labs(
    title    = "Predicted Nominal Tax Burden — Quadratic × Metro Interaction (Model 7)",
    subtitle = "Other controls at means; Colorado state FE applied (illustrative)",
    x        = "USFS + DOI Land Share (%)",
    y        = "Predicted HH Tax Burden (nominal $)",
    color    = NULL, fill = NULL
  )

# =============================================================================
# PHASE 8b — SENSITIVITY CHECK: EXCLUDE HIGH-INFLUENCE COUNTIES
# =============================================================================
# Threshold: Cook's D > 4/n from Model 7. Counties flagged are removed from
# both the linear (6b) and quadratic (7) refits for direct comparison.
# =============================================================================

# --- Identify and remove high-influence observations -------------------------
cd7       <- cooks.distance(model7_tax)
threshold7 <- 4 / length(cd7)

influential_ids <- data2 |>
  mutate(cooks_d = cd7) |>
  filter(cooks_d > threshold7) |>
  select(state_name, county_name, metro, usfs_doi_share_pct, HH_TaxBurden, cooks_d) |>
  arrange(desc(cooks_d))

influential_ids  # inspect flagged counties

data2_trim <- data2 |>
  anti_join(influential_ids |> select(state_name, county_name),
            by = c("state_name", "county_name"))

cat("Counties removed:", nrow(data2) - nrow(data2_trim), "\n")

# --- Refit models on trimmed data --------------------------------------------
model6b_trim <- lm(
  HH_TaxBurden ~ usfs_doi_share_pct + med_hh_income + metro +
                 unincorp_share_pct + pct_bachelors_plus + factor(state_name),
  data = data2_trim
)

model7_trim <- lm(
  HH_TaxBurden ~ (usfs_doi_c + usfs_doi_c2) * metro +
                  med_hh_income + unincorp_share_pct +
                  pct_bachelors_plus + factor(state_name),
  data = data2_trim
)

# --- Coefficient comparison table --------------------------------------------
modelsummary(
  list(
    "Linear — Full"       = model6b_tax,
    "Linear — Trimmed"    = model6b_trim,
    "Quadratic — Full"    = model7_tax,
    "Quadratic — Trimmed" = model7_trim
  ),
  vcov      = \(m) vcovHC(m, type = "HC3"),
  coef_omit = fe_pattern,
  coef_rename = c(
    "usfs_doi_share_pct"         = "USFS+DOI Share (%)",
    "usfs_doi_c"                 = "USFS+DOI Share (centered)",
    "usfs_doi_c2"                = "USFS+DOI Share² (centered)",
    "med_hh_income"              = "Median HH Income ($)",
    "metroNon-Metro"             = "Non-Metro",
    "unincorp_share_pct"         = "Unincorporated Share (%)",
    "pct_bachelors_plus"         = "Bachelor's Degree or Higher (%)",
    "usfs_doi_c:metroNon-Metro"  = "USFS+DOI Share × Non-Metro",
    "usfs_doi_c2:metroNon-Metro" = "USFS+DOI Share² × Non-Metro"
  ),
  stars   = c("*" = 0.1, "**" = 0.05, "***" = 0.01),
  gof_map = c("nobs", "r.squared", "adj.r.squared"),
  notes   = "HC3 robust SEs. State FEs included but not shown. Trimmed = 20 high-Cook's-D counties excluded.",
  title   = "Sensitivity Check: Full vs. Influential-Observation-Trimmed Samples"
)

# --- Overlay predicted curves: full vs. trimmed ------------------------------
pred_trim <- expand_grid(
  usfs_doi_share_pct = seq(0, 100, by = 0.5),
  metro              = c("Metro", "Non-Metro"),
  med_hh_income      = mean(data2_trim$med_hh_income,      na.rm = TRUE),
  unincorp_share_pct = mean(data2_trim$unincorp_share_pct, na.rm = TRUE),
  pct_bachelors_plus = mean(data2_trim$pct_bachelors_plus, na.rm = TRUE),
  state_name         = "Colorado"
) |>
  mutate(usfs_doi_c  = usfs_doi_share_pct - fed_mean,
         usfs_doi_c2 = usfs_doi_c^2)

bind_rows(
  pred_grid7 |>
    bind_cols(predict(model7_tax,  newdata = pred_grid7, interval = "confidence") |>
                as_tibble()) |>
    mutate(sample = "Full"),
  pred_trim |>
    bind_cols(predict(model7_trim, newdata = pred_trim,  interval = "confidence") |>
                as_tibble()) |>
    mutate(sample = "Trimmed (−20 counties)")
) |>
  ggplot(aes(x = usfs_doi_share_pct, y = fit,
             color = sample, linetype = sample)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = sample), alpha = 0.10, color = NA) +
  geom_line(linewidth = 0.85) +
  scale_y_continuous(labels = scales::dollar) +
  facet_wrap(~metro) +
  labs(
    title    = "Sensitivity Check: Full vs. Trimmed Sample — Quadratic Model 7",
    subtitle = "Shaded = 95% CI; 20 high-Cook's-D counties excluded in trimmed version",
    x        = "USFS + DOI Land Share (%)",
    y        = "Predicted HH Tax Burden (nominal $)",
    color    = NULL, fill = NULL, linetype = NULL
  )
