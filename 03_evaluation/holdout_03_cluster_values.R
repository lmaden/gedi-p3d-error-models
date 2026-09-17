#!/usr/bin/env Rscript
# =====================================================================
# holdout_03_cluster_values.R
#
# Step 03: the two remaining cluster operations.
#
#   S18 PPC: posterior predictive draws for CHM 16-site (overlay
#            against observed for density check) + 95% PI coverage
#            verification .
#
#   S12 DTM: rhat / neff_ratio extraction for DTM 18-site, to mirror
#            the CHM diagnostics already produced in holdout scoring.
#
# Not produced here:
#   S10 — per-site residual histograms (uses the holdout residuals CSV)
#   S13 — temporal-gap polish (data join at the plotting step)
#   S19 — coefficient comparison (uses an existing CSV)
#   S22 — variance proportions (uses holdout scoring variance partition CSV)
#
# Outputs:
#   section_L_s18_chm_ppc_summary.csv      - n, mean, sd of observed
#                                            and posterior-predictive
#                                            means; coverage rates
#   section_L_s18_chm_ppc_draws.csv        - subsampled (5k) observed
#                                            and 100 yrep draws for
#                                            density overlay plotting
#                                            at the plotting step
#   section_L_s18_dtm_ppc_summary.csv      - same for DTM (verifies v04
#                                            "94.0% coverage" claim)
#   section_L_s12_dtm_mcmc_diagnostics.csv - parameter / rhat /
#                                            neff_ratio for DTM
#
# Wall-time estimate: 10-20 min (PPC posterior draws are the slow
# step; DTM rhat/ess is seconds).
# =====================================================================

suppressPackageStartupMessages({
  library(brms)
  library(posterior)
  library(data.table)
})

Sys.setenv(DISPLAY = "")
options(device = pdf)

if (!exists("log_progress"))  source("analysis_utils.R")
if (!exists("PROJECT_ROOT"))  source("analysis_config.R")

manuscript_tables_dir <- file.path(PROJECT_ROOT, "manuscript_tables")

log_progress("================================================================")
log_progress("S18 PPC and S12 DTM diagnostics")
log_progress("================================================================")

# ---------------------------------------------------------------------
# Load both brmsfits
# ---------------------------------------------------------------------

log_progress("Loading CHM 16-site brmsfit")
cp_chm <- load_checkpoint("10b_chm_sensitivity")
fit_chm <- cp_chm$fit_chm_16

log_progress("Loading DTM 18-site brmsfit")
cp_dtm <- load_checkpoint("10_models_stage2")
fit_dtm <- cp_dtm$fit_dtm_s2

# ---------------------------------------------------------------------
# (S12 DTM) Convergence diagnostics — fast, do first
# ---------------------------------------------------------------------

log_subsection("S12 DTM: convergence diagnostics")

dtm_rhats <- brms::rhat(fit_dtm)
dtm_ess   <- brms::neff_ratio(fit_dtm)

log_progress(sprintf("  DTM parameters with rhat: %d", length(dtm_rhats)))
log_progress(sprintf("  rhat: max = %.4f, p95 = %.4f, p99 = %.4f",
                     max(dtm_rhats, na.rm = TRUE),
                     quantile(dtm_rhats, 0.95, na.rm = TRUE),
                     quantile(dtm_rhats, 0.99, na.rm = TRUE)))
log_progress(sprintf("  neff_ratio: min = %.4f, p05 = %.4f",
                     min(dtm_ess, na.rm = TRUE),
                     quantile(dtm_ess, 0.05, na.rm = TRUE)))
log_progress(sprintf("  N params with rhat > 1.01: %d",
                     sum(dtm_rhats > 1.01, na.rm = TRUE)))

dtm_mcmc <- data.table(
  parameter  = names(dtm_rhats),
  rhat       = unname(dtm_rhats),
  neff_ratio = unname(dtm_ess[names(dtm_rhats)])
)
fwrite(dtm_mcmc, file.path(manuscript_tables_dir,
                          "section_L_s12_dtm_mcmc_diagnostics.csv"))

# Quick CHM cross-check for caption rewrite
chm_param_n <- length(brms::rhat(fit_chm))
dtm_param_n <- length(dtm_rhats)
log_progress(sprintf("  S12 caption parameter counts: CHM = %d, DTM = %d",
                     chm_param_n, dtm_param_n))

# ---------------------------------------------------------------------
# (S18) Posterior predictive checks for both products
# ---------------------------------------------------------------------

log_subsection("S18 CHM 16-site: posterior predictive draws and coverage")

# Generate posterior predictive draws on a stratified subsample of the
# fit's training data. The PPC plot in v04 is a density overlay of
# observed vs many yrep draws, plus an implied 95% PI coverage stat.
set.seed(20260512)
chm_data <- as.data.table(fit_chm$data)
n_chm <- nrow(chm_data)
ppc_n <- min(5000L, n_chm)
ppc_idx <- sort(sample(seq_len(n_chm), ppc_n))
chm_ppc_data <- chm_data[ppc_idx]
log_progress(sprintf("  CHM PPC subsample: %d obs from %d training rows",
                     ppc_n, n_chm))

log_progress("  Computing posterior_predict on subsample (slow)")
n_draws <- 100  # 100 draws is plenty for density overlay
chm_yrep <- posterior_predict(fit_chm, newdata = chm_ppc_data,
                              ndraws = n_draws,
                              allow_new_levels = FALSE)
log_progress(sprintf("  yrep matrix: %d draws x %d obs",
                     nrow(chm_yrep), ncol(chm_yrep)))

# Save observed + a sample of yrep draws for density overlay plotting at the plotting step
chm_obs <- chm_ppc_data$chm_error_mean
# Store yrep in long format: row = obs, cols = obs_idx, draw, yrep_value
chm_ppc_long <- data.table(
  obs_idx   = rep(seq_len(ppc_n), each = n_draws),
  draw      = rep(seq_len(n_draws), times = ppc_n),
  yrep      = as.vector(t(chm_yrep)),
  observed  = rep(chm_obs, each = n_draws)
)
fwrite(chm_ppc_long, file.path(manuscript_tables_dir,
                              "section_L_s18_chm_ppc_draws.csv"))

# Compute PI coverage on all training data (NOT just the subsample,
# since coverage stat needs the full holdout fan)
log_progress("  Computing 95% PI coverage on full training data")
# Use predictive_interval to get per-obs 95% bands
ppc_intervals <- predictive_interval(fit_chm, newdata = chm_data,
                                     prob = 0.95,
                                     allow_new_levels = FALSE)
covered <- chm_data$chm_error_mean >= ppc_intervals[, 1] &
           chm_data$chm_error_mean <= ppc_intervals[, 2]
chm_coverage <- mean(covered, na.rm = TRUE)
log_progress(sprintf("  CHM 95%% PI coverage: %.4f (v04 caption: 94.0%%)",
                     chm_coverage))

chm_ppc_summary <- data.table(
  product           = "CHM_16_site",
  n_observed        = sum(!is.na(chm_data$chm_error_mean)),
  obs_mean          = mean(chm_data$chm_error_mean, na.rm = TRUE),
  obs_sd            = sd(chm_data$chm_error_mean, na.rm = TRUE),
  obs_q025          = quantile(chm_data$chm_error_mean, 0.025, na.rm = TRUE),
  obs_q975          = quantile(chm_data$chm_error_mean, 0.975, na.rm = TRUE),
  yrep_mean         = mean(chm_yrep, na.rm = TRUE),
  yrep_sd           = sd(chm_yrep, na.rm = TRUE),
  coverage_95pi     = chm_coverage,
  n_ppc_subsample   = ppc_n,
  n_draws           = n_draws
)
fwrite(chm_ppc_summary, file.path(manuscript_tables_dir,
                                 "section_L_s18_chm_ppc_summary.csv"))

# DTM PPC: coverage check only (faster, no draws stored — DTM half
# of S18 isn't being regenerated, but we verify the 94.0% claim)
log_subsection("S18 DTM 18-site: coverage verification only")
dtm_data <- as.data.table(fit_dtm$data)
log_progress("  Computing 95% PI coverage for DTM 18-site")
dtm_intervals <- predictive_interval(fit_dtm, newdata = dtm_data,
                                     prob = 0.95,
                                     allow_new_levels = FALSE)
dtm_covered <- dtm_data$dtm_error_mean >= dtm_intervals[, 1] &
               dtm_data$dtm_error_mean <= dtm_intervals[, 2]
dtm_coverage <- mean(dtm_covered, na.rm = TRUE)
log_progress(sprintf("  DTM 95%% PI coverage: %.4f (v04 caption: 94.0%%)",
                     dtm_coverage))

dtm_ppc_summary <- data.table(
  product       = "DTM_18_site",
  n_observed    = sum(!is.na(dtm_data$dtm_error_mean)),
  obs_mean      = mean(dtm_data$dtm_error_mean, na.rm = TRUE),
  obs_sd        = sd(dtm_data$dtm_error_mean, na.rm = TRUE),
  coverage_95pi = dtm_coverage,
  n_data        = nrow(dtm_data)
)
fwrite(dtm_ppc_summary, file.path(manuscript_tables_dir,
                                 "section_L_s18_dtm_ppc_summary.csv"))

# ---------------------------------------------------------------------
# Console summary
# ---------------------------------------------------------------------

log_progress("================================================================")
log_progress("SUMMARY")
log_progress("================================================================")
log_progress("")
log_progress("[S12 caption parameter counts]")
log_progress(sprintf("  CHM: %d  (v04 caption '193 (CHM)' -> update to %d)",
                     chm_param_n, chm_param_n))
log_progress(sprintf("  DTM: %d  (v04 caption probably reports a number for DTM)",
                     dtm_param_n))
log_progress("")
log_progress("[S12 DTM convergence]")
log_progress(sprintf("  Max rhat:        %.4f", max(dtm_rhats, na.rm = TRUE)))
log_progress(sprintf("  Min neff_ratio:  %.4f", min(dtm_ess, na.rm = TRUE)))
log_progress(sprintf("  N rhat > 1.01:   %d", sum(dtm_rhats > 1.01, na.rm = TRUE)))
log_progress("")
log_progress("[S18 PPC coverage — gates caption '94.0% coverage' claim]")
log_progress(sprintf("  CHM 16-site:  %.3f%%  (v04: 94.0%%; §2A target: 94.0%%)",
                     chm_coverage * 100))
log_progress(sprintf("  DTM 18-site:  %.3f%%  (v04: 94.0%%; §2B target: 94.0%%)",
                     dtm_coverage * 100))
log_progress("")
log_progress("Outputs written:")
log_progress("  section_L_s18_chm_ppc_summary.csv")
log_progress("  section_L_s18_chm_ppc_draws.csv         (5k obs x 100 yrep draws)")
log_progress("  section_L_s18_dtm_ppc_summary.csv       (coverage verification only)")
log_progress("  section_L_s12_dtm_mcmc_diagnostics.csv")
log_progress("================================================================")
log_progress("Complete. All cluster work for the holdout sequence is done.")
log_progress("================================================================")
