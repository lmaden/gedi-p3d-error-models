# =====================================================================
# section_13_summary.R
# Overall summary and final outputs
# =====================================================================
# Ensure dependencies loaded
if (!exists("log_progress")) source("analysis_utils.R")
if (!exists("PROJECT_ROOT")) source("analysis_config.R")

if (!exists("chm_df") || !exists("dtm_df")) {
  log_progress("⚠ Data not loaded. Loading from checkpoint...")
  if (checkpoint_exists("01_data_ingest")) {
    data <- load_checkpoint("01_data_ingest")
    chm_df <- data$chm_df
    dtm_df <- data$dtm_df
  } else {
    stop("Data not available. Run: source('section_01_ingest.R') first")
  }
}

log_progress("Creating overall summary...")

# Compute overall statistics
overall <- tibble(
  metric = c("CHM_bias","CHM_RMSE","CHM_MAE","CHM_NMAD","CHM_Q95abs",
             "DTM_bias","DTM_RMSE","DTM_MAE","DTM_NMAD","DTM_Q95abs",
             "PPC_coverage_CHM","PPC_coverage_DTM",
             "BayesR2_CHM_mean","BayesR2_DTM_mean"),
  value = c(
    mean(chm_df$chm_error_mean, na.rm=TRUE),
    sqrt(mean(chm_df$chm_error_mean^2, na.rm=TRUE)),
    mean(abs(chm_df$chm_error_mean), na.rm=TRUE),
    nm_ad(chm_df$chm_error_mean),
    q95_abs(chm_df$chm_error_mean),
    mean(dtm_df$dtm_error_mean, na.rm=TRUE),
    sqrt(mean(dtm_df$dtm_error_mean^2, na.rm=TRUE)),
    mean(abs(dtm_df$dtm_error_mean), na.rm=TRUE),
    nm_ad(dtm_df$dtm_error_mean),
    q95_abs(dtm_df$dtm_error_mean),
    if (exists("cov_chm")) cov_chm else NA_real_,
    if (exists("cov_dtm")) cov_dtm else NA_real_,
    if (exists("br2_chm")) mean(br2_chm) else NA_real_,
    if (exists("br2_dtm")) mean(br2_dtm) else NA_real_
  )
)

write_csv(overall, file.path(out_tables, "overall_summary.csv"))

# Print summary to console
log_subsection("Overall Performance Summary")
cat("\n")
cat("CHM Performance:\n")
cat(sprintf("  Bias:      %+.3f m\n", overall$value[overall$metric == "CHM_bias"]))
cat(sprintf("  RMSE:      %.3f m\n", overall$value[overall$metric == "CHM_RMSE"]))
cat(sprintf("  MAE:       %.3f m\n", overall$value[overall$metric == "CHM_MAE"]))
cat(sprintf("  NMAD:      %.3f m\n", overall$value[overall$metric == "CHM_NMAD"]))
cat(sprintf("  Q95(|e|):  %.3f m\n", overall$value[overall$metric == "CHM_Q95abs"]))
cat("\n")
cat("DTM Performance:\n")
cat(sprintf("  Bias:      %+.3f m\n", overall$value[overall$metric == "DTM_bias"]))
cat(sprintf("  RMSE:      %.3f m\n", overall$value[overall$metric == "DTM_RMSE"]))
cat(sprintf("  MAE:       %.3f m\n", overall$value[overall$metric == "DTM_MAE"]))
cat(sprintf("  NMAD:      %.3f m\n", overall$value[overall$metric == "DTM_NMAD"]))
cat(sprintf("  Q95(|e|):  %.3f m\n", overall$value[overall$metric == "DTM_Q95abs"]))
cat("\n")

if (exists("cov_chm") && exists("cov_dtm")) {
  cat("Model Calibration:\n")
  cat(sprintf("  CHM coverage: %.1f%%\n", 100*cov_chm))
  cat(sprintf("  DTM coverage: %.1f%%\n", 100*cov_dtm))
  cat("\n")
}

if (exists("br2_chm") && exists("br2_dtm")) {
  cat("Model Fit (Bayes R²):\n")
  cat(sprintf("  CHM: %.3f\n", mean(br2_chm)))
  cat(sprintf("  DTM: %.3f\n", mean(br2_dtm)))
  cat("\n")
}

log_progress("✓ Summary complete")

# Save checkpoint for interactive mode
if (!exists("BATCH_MODE") || !BATCH_MODE) {
  log_progress("Saving checkpoint for interactive mode...")
  save_checkpoint("13_summary", list(
    overall = overall,
    cov_chm = if(exists("cov_chm")) cov_chm else NULL,
    cov_dtm = if(exists("cov_dtm")) cov_dtm else NULL
  ))
}