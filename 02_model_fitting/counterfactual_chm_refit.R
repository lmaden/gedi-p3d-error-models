#!/usr/bin/env Rscript
# =====================================================================
# counterfactual_chm_refit.R
#
# Refits the Chapter 1 CHM Stage 2 Bayesian hierarchical
# model on the 3DEP-DTM-substituted error series (err_alt), to test
# whether the CHM model's 90%-site-random variance structure persists
# when the CNN-DTM-compensation layer is removed.
#
# Motivation. The DTM-substitution step established that CHM error and DTM error
# are inversely correlated at most sites (median r = -0.60), so the
# CNN DTM is acting as a partial compensator for DSM overshoot inside
# CHM = DSM - DTM. If that compensation is site-varying, some of what
# the Ch1 CHM model absorbed into site random intercepts may actually
# be the site-varying residual of CNN compensation. This refit tests
# that by fitting the same hierarchical structure to err_alt (CHM error
# with 3DEP DTM swapped in for CNN DTM), where the compensation layer
# has been removed.
#
# Design. Everything is held EXACTLY the same as the Ch1 Stage 2 CHM
# fit (formula, priors, Student-t family, sampler settings, sample
# fraction) except for the response variable and the site set.
#
#   Response:  err_alt = err_p3d + err_dtm
#              = chm_error_mean + dtm_error_mean
#              (algebraically identical to
#               (p3d_chm_mean + p3d_dtm_mean) - dep_dtm_mean - als_chm_mean)
#   Sample:    rows of mod_chm_s2 whose shot_number also appears in
#              mod_dtm_s2 AND whose site != 10 (CLBJ has no 3DEP)
#   Sites:     18 (Site 10 excluded; all others retained)
#
# Two fits are produced, in parallel with 16-site:
#   fit_alt_18 : all 18 sites with 3DEP coverage
#   fit_alt_15 : 15-site sensitivity (drops flagged 1, 2, 3)
#
# Comparison targets from 16-site (section_10b_chm_sensitivity_refit):
#   - fit_chm_full (19-site Ch1 CHM)     : bias -1.054, sd_site 4.069
#   - fit_chm_16   (16-site sensitivity) : bias -2.425, sd_site 1.516
#
# Expected output axis:
#   If fit_alt_18's sd_site is close to fit_chm_full's (~4), then the
#   CNN-DTM-compensation was NOT what the site intercepts were absorbing.
#   If fit_alt_18's sd_site drops substantially, then much of the Ch1
#   site-RE structure WAS the CNN compensation residual -- a major
# mechanistic finding.
#
# Inputs:
#   /gpfs/data1/vclgp/lmaden/chpt1/checkpoints/08_model_prep.rds
#   /gpfs/data1/vclgp/lmaden/chpt1/checkpoints/10_models_stage2.rds
#
# Outputs:
#   checkpoints/groundwork_task4_refit.rds
#   manuscript_tables/coupling_refit_comparison.csv
#   manuscript_tables/coupling_refit_fixef_comparison.csv
#
# Runtime expectation: 16-site's single 16-site refit took ~58 hours on
# the cluster. This script runs TWO refits; budget roughly 5 days
# wall-clock. Both refits are cached via brm(file=...) so interruptions
# can be resumed.
# =====================================================================

suppressPackageStartupMessages({
  library(brms)
  library(data.table)
  library(posterior)
  library(loo)
})

# Source project infrastructure.
if (!exists("log_progress"))  source("analysis_utils.R")
if (!exists("PROJECT_ROOT"))  source("analysis_config.R")

# Sanity-check required config variables up front.
.required_vars <- c("PROJECT_ROOT", "CPU_BUDGET", "RAM_FRAC", "MODELS_DIR")
.missing_vars  <- .required_vars[!vapply(.required_vars, exists, logical(1))]
if (length(.missing_vars) > 0) {
  stop("Missing required config variables: ",
       paste(.missing_vars, collapse = ", "),
       ". Check analysis_config.R.")
}

manuscript_tables_dir <- file.path(PROJECT_ROOT, "manuscript_tables")

# Site exclusion conventions (tracker-level character codes, consistent
# with mod_chm_s2$site which holds character site labels -- same
# convention as section_10b).
SITE_EXCLUDE_NO_3DEP <- c("10")                # Site 10 has no 3DEP DTM
SITE_FLAGGED         <- c("1", "2", "3")       # flagged ALS reference sites

log_progress("================================================================")
log_progress("CHM refit on 3DEP-DTM-substituted error (err_alt)")
log_progress(sprintf("  No-3DEP exclusion: %s",
                     paste(SITE_EXCLUDE_NO_3DEP, collapse = ", ")))
log_progress(sprintf("  Flagged (for sensitivity fit only): %s",
                     paste(SITE_FLAGGED, collapse = ", ")))
log_progress("================================================================")

# ---------------------------------------------------------------------
# Load upstream checkpoints
# ---------------------------------------------------------------------

log_progress("Loading stage-2 checkpoint (for formula, priors, full fit)...")
s2 <- load_checkpoint("10_models_stage2")
stopifnot(
  "fit_chm_s2"     %in% names(s2),
  "chm_formula_s2" %in% names(s2)
)
fit_chm_full   <- s2$fit_chm_s2
chm_formula_s2 <- s2$chm_formula_s2

log_progress("Loading model-prep checkpoint (for mod_chm_s2 and mod_dtm_s2)...")
mp <- load_checkpoint("08_model_prep")
stopifnot("mod_chm_s2" %in% names(mp), "mod_dtm_s2" %in% names(mp))
mod_chm_s2 <- mp$mod_chm_s2
mod_dtm_s2 <- mp$mod_dtm_s2

# ---------------------------------------------------------------------
# Build the err_alt response by joining CHM and DTM Stage 2 samples
# ---------------------------------------------------------------------
# mod_chm_s2 has chm_error_mean (err_p3d) but not the DTM-side info we
# need to reconstruct err_alt. We use the identity:
#
#   err_alt  =  err_p3d + err_dtm
#            =  chm_error_mean + dtm_error_mean
#
# Both are already in their respective mod_*_s2 tables. Inner-join on
# shot_number to get footprints present in BOTH Stage 2 samples. This
# will be smaller than mod_chm_s2 because (a) the 33% CHM and DTM
# samples were drawn independently, and (b) Site 10 is missing from
# mod_dtm_s2 by pipeline exclusion.

log_progress("Joining mod_chm_s2 and mod_dtm_s2 on shot_number...")
chm_dt <- as.data.table(mod_chm_s2)
dtm_dt <- as.data.table(mod_dtm_s2)

# Rename the DTM response before join to avoid a collision.
setnames(dtm_dt, "dtm_error_mean", ".dtm_error_mean")

# Inner join. All predictors come from mod_chm_s2 (the authoritative CHM
# sample); the only column pulled from the DTM side is .dtm_error_mean.
joined <- merge(
  chm_dt,
  dtm_dt[, .(shot_number, .dtm_error_mean)],
  by = "shot_number", all = FALSE
)

# Compute err_alt and overwrite the response column so the Ch1 formula
# object attaches cleanly. CRITICAL: chm_formula_s2 has LHS =
# chm_error_mean; we swap the underlying values while keeping the column
# name so the formula, contrasts, standata structure, and priors all
# map identically.
joined[, chm_error_mean_ORIG := chm_error_mean]
joined[, chm_error_mean      := chm_error_mean + .dtm_error_mean]

# Guard: |err_alt| <= 100 (same as section_01 does for err_p3d). Should
# be a near no-op; the DTM-substitution step already applied this guard at ingest.
n_before_guard <- nrow(joined)
joined <- joined[abs(chm_error_mean) <= 100]
n_after_guard  <- nrow(joined)
if (n_before_guard != n_after_guard) {
  log_progress(sprintf("  |err_alt| <= 100 guard removed %d rows",
                       n_before_guard - n_after_guard))
}

log_progress(sprintf("  mod_chm_s2        : %s rows",
                     format(nrow(chm_dt), big.mark = ",")))
log_progress(sprintf("  mod_dtm_s2        : %s rows",
                     format(nrow(dtm_dt), big.mark = ",")))
log_progress(sprintf("  Joined sample     : %s rows (err_alt series)",
                     format(nrow(joined), big.mark = ",")))
log_progress(sprintf("  err_alt: mean = %+.3f, sd = %.3f",
                     mean(joined$chm_error_mean),
                     sd(joined$chm_error_mean)))

# ---------------------------------------------------------------------
# Build the two fit datasets: 18-site and 15-site
# ---------------------------------------------------------------------

# 18-site: drop Site 10. Site 10 should already be absent from the join
# (because mod_dtm_s2 excluded it), but we assert this.
stopifnot(!"10" %in% as.character(joined$site))

mod_alt_18 <- copy(joined)
mod_alt_18[, site := droplevels(factor(site))]

sites_18 <- sort(unique(as.character(mod_alt_18$site)))
log_progress(sprintf("  18-site sample: %s rows, %d sites: %s",
                     format(nrow(mod_alt_18), big.mark = ","),
                     length(sites_18), paste(sites_18, collapse = ", ")))
stopifnot(length(sites_18) == 18L)

# 15-site: drop flagged sites. This is the err_alt analogue of 16-site.
keep_15 <- !as.character(mod_alt_18$site) %in% SITE_FLAGGED
mod_alt_15 <- mod_alt_18[keep_15, ]
mod_alt_15[, site := droplevels(factor(site))]

sites_15 <- sort(unique(as.character(mod_alt_15$site)))
log_progress(sprintf("  15-site sample: %s rows, %d sites: %s",
                     format(nrow(mod_alt_15), big.mark = ","),
                     length(sites_15), paste(sites_15, collapse = ", ")))
stopifnot(length(sites_15) == 15L,
          !any(SITE_FLAGGED %in% sites_15))

# ---------------------------------------------------------------------
# Priors and sampler settings (matched EXACTLY to section_10 / 10b)
# ---------------------------------------------------------------------

log_progress("Defining priors and control settings (matched to section_10)...")

priors_chm <- fit_chm_full$prior

ctrl_chm <- list(
  adapt_delta   = 0.98,
  max_treedepth = 15,
  step_size     = 0.0005
)

n_chains <- 4
n_iter   <- 3000
n_warmup <- 1500

parallel_chains   <- max(1L, min(n_chains, CPU_BUDGET))
threads_per_chain <- if (CPU_BUDGET >= 8) 2L else 1L

use_cmdstanr   <- requireNamespace("cmdstanr", quietly = TRUE)
backend_choice <- if (use_cmdstanr) "cmdstanr" else "rstan"

full_backend <- tryCatch(fit_chm_full$backend, error = function(e) "unknown")
log_progress(sprintf("  Backend: %s (full fit used: %s)",
                     backend_choice, full_backend))
if (!is.na(full_backend) && full_backend != "unknown" &&
    full_backend != backend_choice) {
  warning(sprintf(
    "Backend mismatch: full fit used '%s' but refit will use '%s'. Results are still valid but introduce a minor nuisance source of variation.",
    full_backend, backend_choice))
}
log_progress(sprintf("  Chains: %d parallel, %d thread(s) per chain",
                     parallel_chains, threads_per_chain))
log_progress(sprintf("  Iter: %d (warmup %d)", n_iter, n_warmup))
log_progress(sprintf("  Control: adapt_delta=%.2f, max_treedepth=%d, step_size=%.4f",
                     ctrl_chm$adapt_delta, ctrl_chm$max_treedepth,
                     ctrl_chm$step_size))

dir.create(MODELS_DIR, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------
# RAM estimation (mirrors section_10b)
# ---------------------------------------------------------------------

estimate_ram_needed <- function(n_rows, n_chains = 4) {
  base_mb <- 500
  per_row_mb <- 0.15
  per_chain_mb <- base_mb + (n_rows * per_row_mb)
  (per_chain_mb * n_chains) / 1024
}

ram_needed_gb <- estimate_ram_needed(nrow(mod_alt_18), n_chains = 4)
ram_avail_gb  <- tryCatch({
  as.numeric(system("free -g | awk '/^Mem:/ {print $2}'", intern = TRUE)) * RAM_FRAC
}, error = function(e) 16)

log_progress(sprintf("  RAM estimate: need ~%.1f GB, budget ~%.1f GB",
                     ram_needed_gb, ram_avail_gb))
if (ram_needed_gb > ram_avail_gb * 0.8) {
  warning(sprintf("Estimated RAM %.1f GB exceeds 80%% of %.1f GB budget",
                  ram_needed_gb, ram_avail_gb))
}

# ---------------------------------------------------------------------
# Refit 1: 18-site err_alt
# ---------------------------------------------------------------------

log_progress("")
log_progress("============  Refit 1: err_alt, 18 sites  ============")
log_progress("Refitting CHM Stage 2 formula on err_alt (18 sites)...")

fit_alt_18 <- brm(
  formula = chm_formula_s2,
  data    = mod_alt_18,
  prior   = priors_chm,
  family  = student(),
  warmup  = n_warmup,
  iter    = n_iter,
  chains  = n_chains,
  cores   = parallel_chains,
  threads = threading(threads_per_chain),
  control = ctrl_chm,
  backend = backend_choice,
  file    = file.path(MODELS_DIR, "fit_alt_s2_18site"),
  seed    = 2025
)

log_progress("  Refit 1 complete")

max_rhat_18 <- max(brms::rhat(fit_alt_18),       na.rm = TRUE)
min_ess_18  <- min(brms::neff_ratio(fit_alt_18), na.rm = TRUE)
log_progress(sprintf("  Max Rhat: %.3f, Min ESS ratio: %.3f",
                     max_rhat_18, min_ess_18))
if (max_rhat_18 > 1.01) {
  warning("Some 18-site err_alt parameters have Rhat > 1.01")
}

log_progress("Computing LOO for the 18-site refit...")
loo_alt_18 <- tryCatch({
  log_progress("  Attempting LOO with moment matching...")
  loo(fit_alt_18, moment_match = TRUE, cores = parallel_chains)
}, error = function(e) {
  log_progress(sprintf("  Moment matching failed: %s", e$message))
  tryCatch({
    log_progress("  Falling back to basic LOO...")
    loo(fit_alt_18, cores = parallel_chains)
  }, error = function(e2) {
    log_progress(sprintf("  LOO failed entirely: %s", e2$message))
    NULL
  })
})
if (!is.null(loo_alt_18)) {
  log_progress(sprintf("  18-site LOO ELPD: %.1f (SE: %.1f)",
                       loo_alt_18$estimates["elpd_loo", "Estimate"],
                       loo_alt_18$estimates["elpd_loo", "SE"]))
}

# ---------------------------------------------------------------------
# Refit 2: 15-site err_alt (sensitivity -- flagged 3 removed)
# ---------------------------------------------------------------------

log_progress("")
log_progress("============  Refit 2: err_alt, 15 sites (sensitivity)  ============")
log_progress("Refitting CHM Stage 2 formula on err_alt (15 sites, flagged dropped)...")

fit_alt_15 <- brm(
  formula = chm_formula_s2,
  data    = mod_alt_15,
  prior   = priors_chm,
  family  = student(),
  warmup  = n_warmup,
  iter    = n_iter,
  chains  = n_chains,
  cores   = parallel_chains,
  threads = threading(threads_per_chain),
  control = ctrl_chm,
  backend = backend_choice,
  file    = file.path(MODELS_DIR, "fit_alt_s2_15site"),
  seed    = 2025
)

log_progress("  Refit 2 complete")

max_rhat_15 <- max(brms::rhat(fit_alt_15),       na.rm = TRUE)
min_ess_15  <- min(brms::neff_ratio(fit_alt_15), na.rm = TRUE)
log_progress(sprintf("  Max Rhat: %.3f, Min ESS ratio: %.3f",
                     max_rhat_15, min_ess_15))
if (max_rhat_15 > 1.01) {
  warning("Some 15-site err_alt parameters have Rhat > 1.01")
}

log_progress("Computing LOO for the 15-site refit...")
loo_alt_15 <- tryCatch({
  log_progress("  Attempting LOO with moment matching...")
  loo(fit_alt_15, moment_match = TRUE, cores = parallel_chains)
}, error = function(e) {
  log_progress(sprintf("  Moment matching failed: %s", e$message))
  tryCatch({
    log_progress("  Falling back to basic LOO...")
    loo(fit_alt_15, cores = parallel_chains)
  }, error = function(e2) {
    log_progress(sprintf("  LOO failed entirely: %s", e2$message))
    NULL
  })
})
if (!is.null(loo_alt_15)) {
  log_progress(sprintf("  15-site LOO ELPD: %.1f (SE: %.1f)",
                       loo_alt_15$estimates["elpd_loo", "Estimate"],
                       loo_alt_15$estimates["elpd_loo", "SE"]))
}

# ---------------------------------------------------------------------
# Side-by-side global-parameter comparison (uses 16-site fit if available)
# ---------------------------------------------------------------------

extract_summary <- function(fit, label) {
  fx    <- fixef(fit, probs = c(0.025, 0.5, 0.975))
  b0    <- fx["Intercept", ]

  vc      <- VarCorr(fit, probs = c(0.025, 0.5, 0.975))
  sd_site <- vc$site$sd["Intercept", ]

  re       <- ranef(fit)$site[, "Estimate", "Intercept"]
  re_range <- range(re)
  re_sd    <- sd(re)

  pars      <- as_draws_df(fit, variable = c("nu", "b_sigma_Intercept"))
  nu_vec    <- if ("nu" %in% names(pars)) pars$nu else NA_real_
  sigma_vec <- if ("b_sigma_Intercept" %in% names(pars))
                 exp(pars$b_sigma_Intercept) else NA_real_

  data.table(
    fit                      = label,
    n_obs                    = nrow(fit$data),
    n_sites                  = length(unique(as.character(fit$data$site))),
    intercept_mean           = unname(b0["Estimate"]),
    intercept_lo             = unname(b0["Q2.5"]),
    intercept_hi             = unname(b0["Q97.5"]),
    sd_site_mean             = unname(sd_site["Estimate"]),
    sd_site_lo               = unname(sd_site["Q2.5"]),
    sd_site_hi               = unname(sd_site["Q97.5"]),
    re_site_min              = re_range[1],
    re_site_max              = re_range[2],
    re_site_sd_posteriormean = re_sd,
    nu_median                = if (is.numeric(nu_vec))    median(nu_vec)    else NA_real_,
    sigma_median             = if (is.numeric(sigma_vec)) median(sigma_vec) else NA_real_
  )
}

log_progress("")
log_progress("Building global-parameter comparison table...")

cmp_rows <- list(
  extract_summary(fit_chm_full, "Ch1 CHM (19-site, err_p3d)"),
  extract_summary(fit_alt_18,   "Counterfactual (18-site, err_alt)"),
  extract_summary(fit_alt_15,   "Counterfactual (15-site, err_alt, flagged dropped)")
)

# Pull in 16-site's fit_chm_16 if its checkpoint exists.
if (checkpoint_exists("10b_chm_sensitivity")) {
  trackB <- load_checkpoint("10b_chm_sensitivity")
  if (!is.null(trackB$fit_chm_16)) {
    cmp_rows <- c(cmp_rows,
                  list(extract_summary(trackB$fit_chm_16,
                                       "16-site (16-site, err_p3d, flagged dropped)")))
  }
}

cmp <- rbindlist(cmp_rows, use.names = TRUE)
print(cmp)

# ---------------------------------------------------------------------
# Fixed-effect coefficient comparison (18-site err_alt vs Ch1 full)
# ---------------------------------------------------------------------
# Mirrors section_10b exactly. Tells us which predictors' effects change
# direction or lose/gain credible non-zero status under err_alt.

log_progress("Comparing fixed-effect coefficients (err_alt 18-site vs Ch1 full)...")

fx_full <- as.data.table(fixef(fit_chm_full, probs = c(0.025, 0.975)),
                         keep.rownames = "predictor")
fx_alt  <- as.data.table(fixef(fit_alt_18,   probs = c(0.025, 0.975)),
                         keep.rownames = "predictor")

fx_join <- merge(
  fx_full[, .(predictor, est_full = Estimate, lo_full = Q2.5, hi_full = Q97.5)],
  fx_alt[,  .(predictor, est_alt  = Estimate, lo_alt  = Q2.5, hi_alt  = Q97.5)],
  by = "predictor", all = TRUE
)

fx_join[, shift := est_alt - est_full]

fx_join[, includes_zero_full   := (lo_full <= 0) & (hi_full >= 0)]
fx_join[, includes_zero_alt    := (lo_alt  <= 0) & (hi_alt  >= 0)]
fx_join[, status_preserved     := includes_zero_full == includes_zero_alt]
fx_join[, cis_overlap_bw_fits  := !(hi_full < lo_alt | hi_alt < lo_full)]

n_missing_rows <- sum(is.na(fx_join$status_preserved))
if (n_missing_rows > 0) {
  missing_pred <- fx_join[is.na(status_preserved), predictor]
  warning(sprintf(
    "%d predictor(s) present in only one fit: %s. Treating as status-not-preserved.",
    n_missing_rows, paste(missing_pred, collapse = ", ")))
  fx_join[is.na(status_preserved), status_preserved := FALSE]
}

fx_join <- fx_join[order(-abs(shift))]
print(fx_join[seq_len(min(.N, 12))])

# ---------------------------------------------------------------------
# Write CSV outputs
# ---------------------------------------------------------------------

dir.create(manuscript_tables_dir, showWarnings = FALSE, recursive = TRUE)
fwrite(cmp,     file.path(manuscript_tables_dir,
                          "coupling_refit_comparison.csv"))
fwrite(fx_join, file.path(manuscript_tables_dir,
                          "coupling_refit_fixef_comparison.csv"))
log_progress(sprintf("  Wrote: %s/coupling_refit_comparison.csv",
                     manuscript_tables_dir))
log_progress(sprintf("  Wrote: %s/coupling_refit_fixef_comparison.csv",
                     manuscript_tables_dir))

# ---------------------------------------------------------------------
# Save checkpoint
# ---------------------------------------------------------------------

sd_site_full <- cmp[fit == "Ch1 CHM (19-site, err_p3d)",       sd_site_mean]
sd_site_18   <- cmp[fit == "Counterfactual (18-site, err_alt)",         sd_site_mean]
sd_site_15   <- cmp[fit == "Counterfactual (15-site, err_alt, flagged dropped)",
                    sd_site_mean]
sd_site_16   <- if ("16-site (16-site, err_p3d, flagged dropped)" %in% cmp$fit) {
                  cmp[fit == "16-site (16-site, err_p3d, flagged dropped)",
                      sd_site_mean]
                } else NA_real_

shrink_18_vs_full <- 1 - (sd_site_18 / sd_site_full)
shrink_15_vs_16   <- if (!is.na(sd_site_16)) 1 - (sd_site_15 / sd_site_16) else NA_real_

checkpoint_data <- list(
  fit_alt_18         = fit_alt_18,
  fit_alt_15         = fit_alt_15,
  chm_formula_s2     = chm_formula_s2,
  sites_excluded_18  = SITE_EXCLUDE_NO_3DEP,
  sites_excluded_15  = union(SITE_EXCLUDE_NO_3DEP, SITE_FLAGGED),
  comparison_table   = cmp,
  fixef_comparison   = fx_join,
  loo_alt_18         = loo_alt_18,
  loo_alt_15         = loo_alt_15,
  max_rhat_18        = max_rhat_18,
  min_ess_18         = min_ess_18,
  max_rhat_15        = max_rhat_15,
  min_ess_15         = min_ess_15,
  sd_site_full       = sd_site_full,
  sd_site_alt_18     = sd_site_18,
  sd_site_alt_15     = sd_site_15,
  sd_site_trackB_16  = sd_site_16,
  shrink_18_vs_full  = shrink_18_vs_full,
  shrink_15_vs_16    = shrink_15_vs_16
)
save_checkpoint("groundwork_task4_refit", checkpoint_data)

# ---------------------------------------------------------------------
# Interpretive verdict
# ---------------------------------------------------------------------

log_progress("")
log_progress("================================================================")
log_progress("Counterfactual refit interpretive summary")
log_progress("================================================================")

log_progress(sprintf("  Ch1 full         sd_site = %.3f  (err_p3d, 19 sites)",
                     sd_site_full))
log_progress(sprintf("  Counterfactual -site   sd_site = %.3f  (err_alt, 18 sites)",
                     sd_site_18))
if (!is.na(sd_site_16)) {
  log_progress(sprintf("  16-site  sd_site = %.3f  (err_p3d, flagged dropped)",
                       sd_site_16))
}
log_progress(sprintf("  Counterfactual -site   sd_site = %.3f  (err_alt, flagged dropped)",
                     sd_site_15))

log_progress("")
log_progress(sprintf("  sd_site shrinkage (18-site err_alt vs 19-site err_p3d): %+.1f%%",
                     100 * shrink_18_vs_full))
if (!is.na(shrink_15_vs_16)) {
  log_progress(sprintf("  sd_site shrinkage (15-site err_alt vs 16-site err_p3d): %+.1f%%",
                       100 * shrink_15_vs_16))
}

log_progress("")
log_progress("  Interpretation cheat-sheet:")
log_progress("  - Large POSITIVE shrinkage (err_alt sd_site much smaller) ==>")
log_progress("    CNN-DTM compensation residuals were a major component of the")
log_progress("    Ch1 CHM site-random variance. Site effects absorbed compensation")
log_progress("    pattern. Key finding.")
log_progress("  - ~Zero shrinkage ==>")
log_progress("    CNN-DTM compensation residuals were NOT what site effects absorbed.")
log_progress("    Ch1 site intercepts reflect something else (DSM variation,")
log_progress("    phenology, acquisition geometry). Look to Tasks 1/2/3/6.")
log_progress("  - NEGATIVE shrinkage (sd_site grows under err_alt) ==>")
log_progress("    The CNN was homogenizing across sites; removing it reveals more")
log_progress("    site-level structure. Unlikely but possible.")

log_progress("")
log_progress("  The 15-site-vs-16-site comparison (shrink_15_vs_16) is the cleanest")
log_progress("  because both fits drop the same three flagged reference sites, so")
log_progress("  the only difference is the response variable (err_alt vs err_p3d).")
log_progress("  This isolates the CNN-DTM-compensation effect from the reference-")
log_progress("  quality effect.")

log_progress("================================================================")
