
library(tidyverse)
library(here)

source(here("code", "inla_helpers.R"))  # compare_cpo()

# =============================================================================
# RUN ALL — sources the pipeline in order
# =============================================================================
# Each script runs in its own environment, so objects that share a name across
# scripts (data2, coef_labels, W, ...) don't leak between steps. Top-level
# values print as they would interactively; with Rscript, plots go to
# Rplots.pdf. The model scripts (05b_, 06_, 06b_) save their results tables to
# output/tables/ and coefficient plots to output/figures/; this script adds
# loo_compare_all.csv, coefs_all.csv, and coefplot_all.png across them.
#
#   Rscript code/00_run_all.R
#
# 02_ needs a Census API key: CENSUS_API_KEY in the environment, or else the
# census_api_key file in the folder above the repo (read at run time).
#
# Toggle steps below. The models (05b_, 06_, 06b_) are fit in INLA and take
# a few minutes.
# =============================================================================

run_data     <- TRUE   # 01_, 02_, 04_: build output/Cty_Agg.rds
run_unincorp <- FALSE  # 03_: slow spatial step; runs anyway if its CSV is missing
run_models   <- c("05", "05b", "06", "06b")
run_plots    <- TRUE   # 07_

# =============================================================================

run_step <- function(script) {
  message("\n==== ", script, " ", strrep("=", max(0, 70 - nchar(script))))
  t0  <- Sys.time()
  env <- new.env(parent = globalenv())
  source(here("code", script), local = env, print.eval = TRUE)
  message("==== ", script, " done in ",
          format(round(difftime(Sys.time(), t0, units = "mins"), 1)))
  invisible(env)
}

if (run_data && Sys.getenv("CENSUS_API_KEY") == "" &&
    !file.exists(here("..", "census_api_key")))
  stop("No Census API key: set CENSUS_API_KEY or add ../census_api_key (02_ needs it)")

# --- Data --------------------------------------------------------------------
if (run_data) {
  run_step("01_All_Gov_West.R")
  run_step("02_fed_share_acs.R")
}
if (run_unincorp || !file.exists(here("output", "unincorp_share.csv")))
  run_step("03_unincorp_share.R")
if (run_data)
  run_step("04_Cty_Agg.R")

# --- Models ------------------------------------------------------------------
model_scripts <- c(
  "05"  = "05_Models.R",
  "05b" = "05b_Models_INLA.R",
  "06"  = "06_Models_Spatial.R",
  "06b" = "06b_Models_Spatial_Commute.R"
)
envs <- map(model_scripts[run_models], run_step)

# --- LOO across non-spatial and spatial fits ---------------------------------
# CPO rows follow data2 in 05b_, 06_, and 06b_, so all fits compare directly
fits_all <- c(
  if ("05b" %in% run_models) envs[["05b"]]$models_ns,
  if ("06"  %in% run_models) envs[["06"]]$models_sp |>
                               set_names(\(nm) paste(nm, "(BYM2 rook)")),
  if ("06b" %in% run_models) envs[["06b"]]$models_cw |>
                               set_names(\(nm) paste(nm, "(BYM2 commute)"))
)
if (length(fits_all) > 1) {
  loo_all <- compare_cpo(fits_all)
  print(loo_all)
  write_csv(loo_all, here("output", "tables", "loo_compare_all.csv"))
}

# --- Coefficients across specifications --------------------------------------
# M1 and M2 from each script, one color per county-effect specification
coefs_all <- bind_rows(
  if ("05b" %in% run_models) envs[["05b"]]$coefs_ns |>
    mutate(spec = str_extract(model, "(?<=\\()[^)]+(?=\\)$)")),
  if ("06"  %in% run_models) envs[["06"]]$coefs_sp  |> mutate(spec = "BYM2 rook"),
  if ("06b" %in% run_models) envs[["06b"]]$coefs_cw |> mutate(spec = "BYM2 commute")
)
if (nrow(coefs_all) > 0) {
  coefs_all <- coefs_all |>
    mutate(model = str_remove(model, " \\(.*\\)$") |> str_remove(" \\+ BYM2$"),
           label = as.character(label),
           spec  = fct_inorder(spec)) |>
    filter(model != "M0: FE")
  write_csv(coefs_all, here("output", "tables", "coefs_all.csv"))

  # Covariates only; county-effect and residual SDs differ in meaning by spec
  coefplot_all <- coefs_all |>
    filter(!term %in% c("sd_county", "sdcar", "rhocar", "sigma")) |>
    mutate(label = fct_rev(fct_inorder(label))) |>
    ggplot(aes(x = median, xmin = lower, xmax = upper, y = label, color = spec)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_pointrange(position = position_dodge(width = 0.7, orientation = "y", reverse = TRUE),
                    size = 0.3) +
    facet_wrap(~model) +
    labs(title = "Coefficients with 95% CrIs by County-Effect Specification",
         x = "Coefficient (log HH Tax Burden)", y = NULL, color = NULL)
  print(coefplot_all)
  ggsave(here("output", "figures", "coefplot_all.png"), coefplot_all, width = 11, height = 5)
}

# --- Plots -------------------------------------------------------------------
if (run_plots)
  run_step("07_plots.R")
