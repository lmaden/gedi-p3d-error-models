# =====================================================================
# compute_definitive_r2.R
# 
# Definitive R² computation for CHM and DTM models
# Loads models fresh from checkpoint, computes all R² metrics cleanly
# No reliance on existing environment variables
# =====================================================================

library(brms)

# Clear any existing R² variables to avoid confusion
suppressWarnings(rm(list = c("br2_chm", "br2_dtm", "cov_chm", "cov_dtm"), envir = globalenv()))

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("DEFINITIVE R² COMPUTATION\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# =====================================================================
# Load configuration
# =====================================================================

if (!exists("log_progress")) source("analysis_utils.R")
if (!exists("out_tables")) source("analysis_config.R")

# =====================================================================
# Load models FRESH from checkpoint (not from environment)
# =====================================================================

cat("Loading models fresh from checkpoint...\n")

if (checkpoint_exists("10_models_stage2")) {
  cat("  Loading Stage 2 models from checkpoint...\n")
  ckpt <- load_checkpoint("10_models_stage2")
  model_chm <- ckpt$fit_chm_s2
  model_dtm <- ckpt$fit_dtm_s2
  model_stage <- "Stage 2"
  rm(ckpt)  # Clean up
  gc()
} else if (checkpoint_exists("09_models_stage1")) {
  cat("  Loading Stage 1 models from checkpoint...\n")
  ckpt <- load_checkpoint("09_models_stage1")
  model_chm <- ckpt$fit_chm_s1
  model_dtm <- ckpt$fit_dtm_s1
  model_stage <- "Stage 1"
  rm(ckpt)
  gc()
} else {
  stop("No model checkpoints found!")
}

cat(sprintf("  ✓ Loaded %s models\n", model_stage))
cat(sprintf("  CHM model: %d observations\n", nobs(model_chm)))
cat(sprintf("  DTM model: %d observations\n", nobs(model_dtm)))

# =====================================================================
# Compute Bayes R² - Full Model (Fixed + Random Effects)
# =====================================================================

cat("\n")
cat(paste(rep("-", 70), collapse = ""), "\n")
cat("BAYES R² - FULL MODEL (Fixed + Random Effects)\n")
cat(paste(rep("-", 70), collapse = ""), "\n")

cat("  Computing CHM Bayes R² (full model)...\n")
r2_chm_full <- bayes_R2(model_chm, re_formula = NULL)  # NULL = include all random effects
cat(sprintf("    Posterior draws: %d\n", nrow(r2_chm_full)))
cat(sprintf("    CHM R² = %.4f [%.4f, %.4f]\n", 
            median(r2_chm_full[,1]),
            quantile(r2_chm_full[,1], 0.025),
            quantile(r2_chm_full[,1], 0.975)))

cat("  Computing DTM Bayes R² (full model)...\n")
r2_dtm_full <- bayes_R2(model_dtm, re_formula = NULL)
cat(sprintf("    Posterior draws: %d\n", nrow(r2_dtm_full)))
cat(sprintf("    DTM R² = %.4f [%.4f, %.4f]\n", 
            median(r2_dtm_full[,1]),
            quantile(r2_dtm_full[,1], 0.025),
            quantile(r2_dtm_full[,1], 0.975)))

# =====================================================================
# Compute Bayes R² - Fixed Effects Only
# =====================================================================

cat("\n")
cat(paste(rep("-", 70), collapse = ""), "\n")
cat("BAYES R² - FIXED EFFECTS ONLY (Generalizable)\n")
cat(paste(rep("-", 70), collapse = ""), "\n")

cat("  Computing CHM Bayes R² (fixed only)...\n")
r2_chm_fixed <- bayes_R2(model_chm, re_formula = NA)  # NA = exclude random effects
cat(sprintf("    CHM R² (fixed) = %.4f [%.4f, %.4f]\n", 
            median(r2_chm_fixed[,1]),
            quantile(r2_chm_fixed[,1], 0.025),
            quantile(r2_chm_fixed[,1], 0.975)))

cat("  Computing DTM Bayes R² (fixed only)...\n")
r2_dtm_fixed <- bayes_R2(model_dtm, re_formula = NA)
cat(sprintf("    DTM R² (fixed) = %.4f [%.4f, %.4f]\n", 
            median(r2_dtm_fixed[,1]),
            quantile(r2_dtm_fixed[,1], 0.025),
            quantile(r2_dtm_fixed[,1], 0.975)))

# =====================================================================
# Compute Squared Correlation (for comparison with plots)
# =====================================================================

cat("\n")
cat(paste(rep("-", 70), collapse = ""), "\n")
cat("SQUARED CORRELATION (Point Estimate)\n")
cat(paste(rep("-", 70), collapse = ""), "\n")

cat("  Computing posterior predictions...\n")

# CHM
y_obs_chm <- model_chm$data$chm_error_mean
y_pred_chm <- colMeans(posterior_predict(model_chm, ndraws = 500))
cor2_chm <- cor(y_obs_chm, y_pred_chm)^2
cat(sprintf("    CHM cor² = %.4f\n", cor2_chm))

# DTM  
y_obs_dtm <- model_dtm$data$dtm_error_mean
y_pred_dtm <- colMeans(posterior_predict(model_dtm, ndraws = 500))
cor2_dtm <- cor(y_obs_dtm, y_pred_dtm)^2
cat(sprintf("    DTM cor² = %.4f\n", cor2_dtm))

# =====================================================================
# Compute PPC Coverage
# =====================================================================

cat("\n")
cat(paste(rep("-", 70), collapse = ""), "\n")
cat("PPC COVERAGE (95% Prediction Intervals)\n")
cat(paste(rep("-", 70), collapse = ""), "\n")

compute_coverage <- function(model, y_obs, name) {
  y_rep <- posterior_predict(model, ndraws = 500)
  pi_lower <- apply(y_rep, 2, quantile, probs = 0.025)
  pi_upper <- apply(y_rep, 2, quantile, probs = 0.975)
  coverage <- mean((y_obs >= pi_lower) & (y_obs <= pi_upper))
  cat(sprintf("    %s coverage = %.1f%%\n", name, coverage * 100))
  return(coverage)
}

cov_chm <- compute_coverage(model_chm, y_obs_chm, "CHM")
cov_dtm <- compute_coverage(model_dtm, y_obs_dtm, "DTM")

# =====================================================================
# Summary Table
# =====================================================================

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("DEFINITIVE SUMMARY\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

summary_df <- data.frame(
  Model = c("CHM", "DTM"),
  Stage = model_stage,
  N_obs = c(nobs(model_chm), nobs(model_dtm)),
  
  # Full model R²
  R2_full_median = c(median(r2_chm_full[,1]), median(r2_dtm_full[,1])),
  R2_full_lower = c(quantile(r2_chm_full[,1], 0.025), quantile(r2_dtm_full[,1], 0.025)),
  R2_full_upper = c(quantile(r2_chm_full[,1], 0.975), quantile(r2_dtm_full[,1], 0.975)),
  
  # Fixed effects only R²
  R2_fixed_median = c(median(r2_chm_fixed[,1]), median(r2_dtm_fixed[,1])),
  R2_fixed_lower = c(quantile(r2_chm_fixed[,1], 0.025), quantile(r2_dtm_fixed[,1], 0.025)),
  R2_fixed_upper = c(quantile(r2_chm_fixed[,1], 0.975), quantile(r2_dtm_fixed[,1], 0.975)),
  
  # Proportion from random effects
  Prop_from_random = c(
    1 - median(r2_chm_fixed[,1]) / median(r2_chm_full[,1]),
    1 - median(r2_dtm_fixed[,1]) / median(r2_dtm_full[,1])
  ),
  
  # Point estimates
  Cor_squared = c(cor2_chm, cor2_dtm),
  PPC_coverage = c(cov_chm, cov_dtm)
)

# Print summary
cat("Model Performance Summary:\n\n")

cat(sprintf("CHM (%s, n=%d):\n", model_stage, summary_df$N_obs[1]))
cat(sprintf("  Bayes R² (full model):    %.3f [%.3f, %.3f]\n",
            summary_df$R2_full_median[1], summary_df$R2_full_lower[1], summary_df$R2_full_upper[1]))
cat(sprintf("  Bayes R² (fixed effects): %.3f [%.3f, %.3f]\n",
            summary_df$R2_fixed_median[1], summary_df$R2_fixed_lower[1], summary_df$R2_fixed_upper[1]))
cat(sprintf("  Proportion from random:   %.1f%%\n", summary_df$Prop_from_random[1] * 100))
cat(sprintf("  Squared correlation:      %.3f\n", summary_df$Cor_squared[1]))
cat(sprintf("  PPC coverage (95%%):       %.1f%%\n", summary_df$PPC_coverage[1] * 100))

cat(sprintf("\nDTM (%s, n=%d):\n", model_stage, summary_df$N_obs[2]))
cat(sprintf("  Bayes R² (full model):    %.3f [%.3f, %.3f]\n",
            summary_df$R2_full_median[2], summary_df$R2_full_lower[2], summary_df$R2_full_upper[2]))
cat(sprintf("  Bayes R² (fixed effects): %.3f [%.3f, %.3f]\n",
            summary_df$R2_fixed_median[2], summary_df$R2_fixed_lower[2], summary_df$R2_fixed_upper[2]))
cat(sprintf("  Proportion from random:   %.1f%%\n", summary_df$Prop_from_random[2] * 100))
cat(sprintf("  Squared correlation:      %.3f\n", summary_df$Cor_squared[2]))
cat(sprintf("  PPC coverage (95%%):       %.1f%%\n", summary_df$PPC_coverage[2] * 100))

# =====================================================================
# Export
# =====================================================================

cat("\n")
output_file <- file.path(out_tables, "definitive_r2_metrics.csv")
write.csv(summary_df, output_file, row.names = FALSE)
cat(sprintf("✓ Saved to: %s\n", output_file))

# Also export the full posterior draws for reference
saveRDS(list(
  r2_chm_full = r2_chm_full,
  r2_dtm_full = r2_dtm_full,
  r2_chm_fixed = r2_chm_fixed,
  r2_dtm_fixed = r2_dtm_fixed,
  model_stage = model_stage
), file.path(out_tables, "r2_posterior_draws.rds"))
cat(sprintf("✓ Posterior draws saved to: %s\n", file.path(out_tables, "r2_posterior_draws.rds")))

# =====================================================================
# Interpretation Guide
# =====================================================================

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("INTERPRETATION GUIDE\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("Which R² to report:\n")
cat("  • For 'model explains X% of variance': Use Bayes R² (full model)\n")
cat("  • For 'generalizable predictive power': Use Bayes R² (fixed effects)\n")
cat("  • For obs vs pred plots: Squared correlation is fine for visual\n\n")

cat("Why Bayes R² differs from cor²:\n")
cat("  • Bayes R² accounts for posterior uncertainty in predictions\n")
cat("  • cor² is a point estimate that can be inflated by outliers\n")
cat("  • Bayes R² is generally more conservative and appropriate for reporting\n\n")

cat("✓ Definitive R² computation complete\n")
