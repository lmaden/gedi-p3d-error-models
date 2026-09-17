# =============================================================================
# Chapter 1 (GEDI vs P3D) — v03 finalization
# Segment 4.0 pre-flight: full pipeline for placeholder extraction
# =============================================================================
#
# PURPOSE
#   Extract all numeric values needed to fill ~38 placeholders in
#   chapter1_v03.docx and supplemental_figures_v03.docx (segment 4.1 fill-in pass).
#
# SCOPE
#   - Item 9:  Min bulk-ESS per Stage 2 fit (Table 4)
#   - Item 11: PSIS-LOO elpd, SE, Pareto-k diagnostics (Table 4)
#              + within-CHM 16-vs-19 comparison on shared 16-site subset (§3.5)
#   - Item 12: In-sample CRPS per stratum (Table 6)
#   - Item 13: Posterior medians for Bayesian R² (full + decomposed) and
#              Student-t nu (Table 4)
#   - Item 14: Prior justification (option b) — uses prior predictive summary
#              + screening sensitivity in lieu of un-run prior-perturbation refits
#
# ENVIRONMENT
#   R 4.5.0, rh9 cluster
#   PROJ 25 / GDAL 3.11 / LAStools 250710
#   brms 2.23.0, posterior 1.6.1
#   scoringRules 1.1.3 (installed to /gpfs/data1/vclgp/lmaden/Rlib)
#
# CHECKPOINTS USED
#   10b_chm_sensitivity.rds    -> fit_chm_16 + cached loo_chm_16
#   10_models_stage2.rds       -> fit_chm_19 (sensitivity) + fit_dtm_18
#                                 + cached loo_chm_19 + loo_dtm_18
#   09_models_stage1.rds       -> Stage 1 baseline fits (no prior-perturbation refits)
#   manuscript_tables/         -> sensitivity_*.csv (all 16-vs-19 screening
#                                 comparisons, not prior comparisons)
#   tables/prior_predictive_summary.csv -> prior predictive check artifact
#
# RUNTIME
#   ~5 min total on rh9 cluster, dominated by CRPS posterior_predict
#   (~1 min CHM 16-site, ~2.5 min DTM 18-site)
#
# AUTHOR
#   L. Maden, May 2026, chat N+4 pre-flight
# =============================================================================


# -----------------------------------------------------------------------------
# Setup: library path + packages
# -----------------------------------------------------------------------------

# Project-local R library (cluster-specific)
.libPaths(c("/gpfs/data1/vclgp/lmaden/Rlib", .libPaths()))

library(brms)
library(posterior)
library(scoringRules)

# Checkpoint directory
CK <- "/gpfs/data1/vclgp/lmaden/chpt1/checkpoints"


# -----------------------------------------------------------------------------
# Load Stage 2 fit objects + cached LOO objects
# -----------------------------------------------------------------------------

ck_chm <- readRDS(file.path(CK, "10b_chm_sensitivity.rds"))
ck_dtm <- readRDS(file.path(CK, "10_models_stage2.rds"))

# Assignments per Stage 2 schema.
# Note: the 19-site CHM sensitivity fit lives in the DTM checkpoint file,
# not the CHM file. The CHM file holds only the 16-site primary plus diagnostics.
fit_chm_16 <- ck_chm$data$fit_chm_16     # 16-site CHM primary  (124966 obs)
fit_chm_19 <- ck_dtm$data$fit_chm_s2     # 19-site CHM sensitivity (204672 obs)
fit_dtm_18 <- ck_dtm$data$fit_dtm_s2     # 18-site DTM primary  (202278 obs)

loo_chm_16 <- ck_chm$data$loo_chm_16
loo_chm_19 <- ck_dtm$data$loo_chm_s2
loo_dtm_18 <- ck_dtm$data$loo_dtm_s2


# =============================================================================
# Item 11 — PSIS-LOO Table 4 rows (from cached loo objects, no recompute)
# =============================================================================

cat("\n##############################################\n")
cat("## Item 11: PSIS-LOO Table 4 rows\n")
cat("##############################################\n")

for (lbl in c("CHM 16-site primary", "DTM 18-site primary", "CHM 19-site sensitivity")) {
  loo_obj <- switch(lbl,
                    "CHM 16-site primary"       = loo_chm_16,
                    "DTM 18-site primary"       = loo_dtm_18,
                    "CHM 19-site sensitivity"   = loo_chm_19)
  est <- loo_obj$estimates["elpd_loo", "Estimate"]
  se  <- loo_obj$estimates["elpd_loo", "SE"]
  pk_bad_pct <- mean(loo_obj$diagnostics$pareto_k > 0.7) * 100
  cat(sprintf("  %-28s elpd_loo = %.1f (SE %.1f);  Pareto-k > 0.7: %.2f%%\n",
              lbl, est, se, pk_bad_pct))
}


# =============================================================================
# Item 9  — Min bulk-ESS (absolute)
# Item 13 — Bayesian R² (full + decomposed) and Student-t nu, posterior medians
# =============================================================================

cat("\n##############################################\n")
cat("## Items 9 + 13: bulk-ESS, R², nu (posterior medians)\n")
cat("##############################################\n")

share_pct <- function(part, whole) round(100 * part / whole, 1)

for (lbl in c("CHM 16-site primary", "DTM 18-site primary", "CHM 19-site sensitivity")) {
  ft <- switch(lbl,
               "CHM 16-site primary"       = fit_chm_16,
               "DTM 18-site primary"       = fit_dtm_18,
               "CHM 19-site sensitivity"   = fit_chm_19)
  
  # Min bulk-ESS
  bulk_min <- min(summarise_draws(as_draws_df(ft), "ess_bulk")$ess_bulk, na.rm = TRUE)
  
  # Bayesian R²: full and fixed-only (random share = full - fixed)
  r2_full  <- median(bayes_R2(ft, summary = FALSE)[, "R2"])
  r2_fixed <- median(bayes_R2(ft, re_formula = NA, summary = FALSE)[, "R2"])
  r2_rand  <- r2_full - r2_fixed
  
  # Student-t nu
  nu_med <- median(as_draws_df(ft)$nu)
  
  cat(sprintf("\n  %s:\n", lbl))
  cat(sprintf("    Min bulk-ESS:       %.0f\n", bulk_min))
  cat(sprintf("    R² full:            %.3f\n", r2_full))
  cat(sprintf("    R² fixed-effects:   %.3f  (%.1f%% of full)\n",
              r2_fixed, share_pct(r2_fixed, r2_full)))
  cat(sprintf("    R² random-effects:  %.3f  (%.1f%% of full)\n",
              r2_rand, share_pct(r2_rand, r2_full)))
  cat(sprintf("    Student-t nu:       %.2f\n", nu_med))
}


# =============================================================================
# Item 11 (§3.5 sentence) — within-CHM 16-vs-19 elpd comparison
#   loo_compare() errors on different-n fits (124966 vs 204672), so we
#   subset the 19-site pointwise elpd to the 16-site shared subset and
#   compute elpd_diff manually.
# =============================================================================

cat("\n##############################################\n")
cat("## Item 11 (§3.5): 16-vs-19 elpd comparison on shared subset\n")
cat("##############################################\n")

sites_16 <- sort(unique(fit_chm_16$data$site))     # {4-15, 17-20}
idx_19_to_16 <- which(fit_chm_19$data$site %in% sites_16)

stopifnot(length(idx_19_to_16) == nrow(fit_chm_16$data))  # both should be 124966

pw_16        <- loo_chm_16$pointwise[, "elpd_loo"]
pw_19_subset <- loo_chm_19$pointwise[idx_19_to_16, "elpd_loo"]

diff_pw   <- pw_19_subset - pw_16   # positive = 19-site better at that obs
elpd_diff <- sum(diff_pw)
se_diff   <- sqrt(length(diff_pw) * var(diff_pw))

cat(sprintf("  sum elpd, 16-site fit:                  %.1f\n", sum(pw_16)))
cat(sprintf("  sum elpd, 19-site fit (16-site subset): %.1f\n", sum(pw_19_subset)))
cat(sprintf("  elpd_diff (19-site − 16-site):          %.1f (SE %.1f)\n",
            elpd_diff, se_diff))
cat(sprintf("  z = elpd_diff / SE:                     %.2f\n", elpd_diff / se_diff))
cat(sprintf("  per-obs mean elpd, 16-site:             %.4f\n", mean(pw_16)))
cat(sprintf("  per-obs mean elpd, 19-site on subset:   %.4f\n", mean(pw_19_subset)))
cat("  Interpretation: negative elpd_diff => 16-site fit preferred\n")


# =============================================================================
# Item 12 — In-sample CRPS per stratum, sample-based (crps_sample)
#   Methodologically symmetric across CHM and DTM. In-sample to match
#   Table 6's convention (RMSE/MAE/R²/bias/correlation are all in-sample
#   per pred_metrics in 11_diagnostics.rds).
# =============================================================================

crps_by_stratum <- function(fit, label, n_draws = 500) {
  cat(sprintf("\n=== CRPS for %s ===\n", label))
  
  cat("  Drawing posterior predictive (n =", n_draws, ")...\n")
  t0 <- Sys.time()
  pp <- posterior_predict(fit, ndraws = n_draws)
  cat(sprintf("  Done in %.1f sec. Dim: %d × %d\n",
              as.numeric(difftime(Sys.time(), t0, units = "secs")),
              nrow(pp), ncol(pp)))
  
  # brms strips underscores in formula resp name (chmerrormean vs chm_error_mean);
  # match back to the actual fit$data column.
  resp_var <- as.character(formula(fit)$resp)
  if (!resp_var %in% names(fit$data)) {
    candidates <- names(fit$data)
    match_idx <- which(gsub("_", "", candidates) == resp_var)
    if (length(match_idx) == 1) {
      resp_var <- candidates[match_idx]
    } else {
      resp_var <- names(fit$data)[1]  # fallback: first col is the response in brms
    }
  }
  cat("  Response column used:", resp_var, "\n")
  y <- fit$data[[resp_var]]
  stopifnot(length(y) == ncol(pp))
  
  cat("  Computing crps_sample...\n")
  t0 <- Sys.time()
  crps_per_obs <- crps_sample(y = y, dat = t(pp))   # obs × draws
  cat(sprintf("  Done in %.1f sec\n",
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  
  d <- fit$data
  cat("\n  Per-stratum mean CRPS (lower is better):\n")
  
  cat(sprintf("    %-25s n=%6d   CRPS = %.3f\n",
              "All Forest Types", length(y), mean(crps_per_obs)))
  
  for (ft in c("BDF", "DNF", "EBF", "ENF")) {
    idx <- d$lc_l1_code == ft
    if (sum(idx) > 0) {
      cat(sprintf("    %-25s n=%6d   CRPS = %.3f\n",
                  ft, sum(idx), mean(crps_per_obs[idx])))
    } else {
      cat(sprintf("    %-25s n=%6d   CRPS = (no obs)\n", ft, 0))
    }
  }
  
  invisible(list(per_obs = crps_per_obs, mean = mean(crps_per_obs)))
}

cat("\n##############################################\n")
cat("## Item 12: In-sample CRPS per stratum\n")
cat("##############################################\n")

crps_chm16 <- crps_by_stratum(fit_chm_16, "CHM 16-site primary")
crps_dtm18 <- crps_by_stratum(fit_dtm_18, "DTM 18-site primary")


# =============================================================================
# Item 14 — Prior justification artifacts (option b, no refit needed)
#   Prior predictive check + screening sensitivity replace the un-run
#   prior-perturbation refits.
# =============================================================================

cat("\n##############################################\n")
cat("## Item 14: Prior predictive summary (for S2.5 ¶2)\n")
cat("##############################################\n")

prior_pp_path <- "/gpfs/data1/vclgp/lmaden/chpt1/tables/prior_predictive_summary.csv"
if (file.exists(prior_pp_path)) {
  pp_summary <- read.csv(prior_pp_path)
  print(pp_summary)
}

# Screening sensitivity assets already in hand from §3.5 elpd comparison above
# (16-vs-19 site contrast) plus 16-vs-19 CSVs under manuscript_tables/.

cat("\n##############################################\n")
cat("## End of Segment 4.0 pre-flight\n")
cat("##############################################\n")
