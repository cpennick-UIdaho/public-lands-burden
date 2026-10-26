
library(INLA)

# =============================================================================
# INLA HELPERS — sourced by 05b_, 06_, 06b_, and 00_run_all.R
# =============================================================================
# Shared priors, fitting call, posterior summaries, and LOO comparison, so the
# non-spatial, rook BYM2, and commuting BYM2 fits differ only in the county
# effect. Not a pipeline step.
# =============================================================================

# PC prior on every precision (residual, iid, and BYM2): P(SD > 1) = 0.01
pc_prec <- list(prior = "pc.prec", param = c(1, 0.01))

# PC prior on the BYM2 mixing parameter: P(phi < 0.5) = 0.5 (INLA's default)
pc_phi_param <- c(0.5, 0.5)

# Slopes ~ N(0, 1), as in the earlier Stan/brms fits; flat intercept (INLA default)
fit_inla <- function(formula, data) {
  fit <- inla(
    formula,
    family            = "gaussian",
    data              = data,
    control.family    = list(hyper = list(prec = pc_prec)),
    control.fixed     = list(mean = 0, prec = 1),
    control.compute   = list(cpo = TRUE, waic = TRUE, dic = TRUE)
  )
  # Recompute any CPO values INLA flags as unreliable
  if (any(fit$cpo$failure > 0, na.rm = TRUE)) fit <- inla.cpo(fit)
  fit
}

# --- Posterior summaries -----------------------------------------------------
# Hyperparameters are summarised on the scale we report, by transforming their
# internal (log precision / logit phi) marginals
sd_from_log_prec <- \(x) exp(-x / 2)

hyper_sigma <- list(sigma = list(name = "Log precision for the Gaussian observations",
                                 fun  = sd_from_log_prec))

# Posterior median and 95% CrI per term. `hyper`: named list, one entry per
# reported term, giving the internal hyperparameter name and back-transform
tidy_inla <- function(fit, hyper = hyper_sigma) {
  fixed <- fit$summary.fixed |>
    as_tibble(rownames = "term") |>
    transmute(term, median = `0.5quant`, lower = `0.025quant`, upper = `0.975quant`)

  hyp <- imap(hyper, \(h, term) {
    marg <- inla.tmarginal(h$fun, fit$internal.marginals.hyperpar[[h$name]])
    tibble(term   = term,
           median = inla.qmarginal(0.5,   marg),
           lower  = inla.qmarginal(0.025, marg),
           upper  = inla.qmarginal(0.975, marg))
  }) |>
    list_rbind()

  bind_rows(fixed, hyp)
}

# --- Model comparison (LOO via CPO) ------------------------------------------
# sum(log CPO) is the leave-one-out ELPD. Differences and SEs are pointwise
# against the best model, as in loo::loo_compare(). Fits must share row order.
elpd_cpo <- \(fit) sum(log(fit$cpo$cpo))

compare_cpo <- function(fits) {
  ll   <- map(fits, \(m) log(m$cpo$cpo))
  stopifnot(n_distinct(lengths(ll)) == 1)
  best <- ll[[which.max(map_dbl(ll, sum))]]
  tibble(
    model        = names(fits),
    elpd_diff    = map_dbl(ll, \(x) sum(x - best)),
    se_diff      = map_dbl(ll, \(x) sqrt(length(x)) * sd(x - best)),
    elpd_loo     = map_dbl(ll, sum),
    se_elpd_loo  = map_dbl(ll, \(x) sqrt(length(x)) * sd(x)),
    waic         = map_dbl(fits, \(m) m$waic$waic),
    dic          = map_dbl(fits, \(m) m$dic$dic),
    cpo_failures = map_int(fits, \(m) sum(m$cpo$failure > 0, na.rm = TRUE))
  ) |>
    arrange(desc(elpd_diff))
}

# --- HTML tables -------------------------------------------------------------
# Numeric columns rounded to 3 decimals; preformatted strings pass through as
# is. Missing cells (terms not in a model) are left blank.
write_html_table <- function(df, path, caption = NULL) {
  old <- options(knitr.kable.NA = "")
  on.exit(options(old))
  style <- paste("<style>table { border-collapse: collapse; font-family: sans-serif; }",
                 "th, td { padding: 4px 10px; border-bottom: 1px solid #ccc; }",
                 "caption { caption-side: bottom; text-align: left; color: #555; }</style>")
  tbl <- knitr::kable(df, format = "html", digits = 3, caption = caption)
  writeLines(c(style, as.character(tbl)), path)
}
