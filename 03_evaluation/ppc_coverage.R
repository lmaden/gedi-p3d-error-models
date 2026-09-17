# =====================================================================
# ppc_coverage.R
# Computes PPC coverage (calibration check) for CHM and DTM models
# Run after section_11_diagnostics.R or section_10_models_stage2.R
# =====================================================================

library(brms)

# Ensure config is loaded
if (!exists("log_progress")) source("analysis_utils.R")
if (!exists("out_tables")) source("analysis_config.R")

cat("Computing PPC coverage and R² metrics...\n")

# =====================================================================
# Load models if not in environment
# =====================================================================

if (!exists("model_chm_final") || !exists("model_dtm_final")) {
  cat("  Loading models from checkpoint...\n")
  if (checkpoint_exists("10_models_stage2")) {
    data <- load_checkpoint("10_models_stage2")
    model_chm_final <- data$fit_chm_s2
    model_dtm_final <- data$fit_dtm_s2
    cat("  ✓ Stage 2 models loaded\n")
  } else if (checkpoint_exists("09_models_stage1")) {
    data <- load_checkpoint("09_models_stage1")
    model_chm_final <- data$fit_chm_s1
    model_dtm_final <- data$fit_dtm_s1
    cat("  ✓ Stage 1 models loaded\n")
  } else {
    stop("No models available. Run modeling sections first.")
  }
}

# =====================================================================
# Compute PPC Coverage
# =====================================================================

compute_ppc_coverage <- function(model, response_var, name, alpha = 0.95) {
  # Computes the proportion of observations within the (1-alpha)% prediction interval
  # For a well-calibrated model, this should be close to alpha (e.g., 95%)
  
  cat(sprintf("  Computing %s PPC coverage...\n", name))
  
  tryCatch({
    # Get posterior predictions (full distribution, not just mean)
    y_rep <- posterior_predict(model, ndraws = 500)
    
    # Get observed values
    y_obs <- model$data[[response_var]]
    
    # Compute prediction intervals for each observation
    lower <- (1 - alpha) / 2
    upper <- 1 - lower
    
    pi_lower <- apply(y_rep, 2, quantile, probs = lower)
    pi_upper <- apply(y_rep, 2, quantile, probs = upper)
    
    # Check how many observations fall within their prediction interval
    within_interval <- (y_obs >= pi_lower) & (y_obs <= pi_upper)
    coverage <- mean(within_interval)
    
    cat(sprintf("    %s %.0f%% PPC coverage: %.1f%%\n", 
                name, alpha * 100, coverage * 100))
    
    if (abs(coverage - alpha) < 0.05) {
      cat(sprintf("    → Well calibrated (within 5%% of target)\n"))
    } else if (coverage < alpha) {
      cat(sprintf("    → Undercoverage: prediction intervals too narrow\n"))
    } else {
      cat(sprintf("    → Overcoverage: prediction intervals too wide\n"))
    }
    
    return(coverage)
    
  }, error = function(e) {
    cat(sprintf("    ⚠ %s coverage failed: %s\n", name, e$message))
    return(NA_real_)
  })
}

# Compute coverage for both models
cov_chm <- compute_ppc_coverage(model_chm_final, "chm_error_mean", "CHM")
cov_dtm <- compute_ppc_coverage(model_dtm_final, "dtm_error_mean", "DTM")

# =====================================================================
# Compute and Compare R² Metrics
# =====================================================================

cat("\nComputing R² metrics...\n")

# Bayes R² (proper Bayesian measure)
cat("  Bayes R² (accounts for posterior uncertainty):\n")

br2_chm <- tryCatch({
  r2 <- bayes_R2(model_chm_final)
  cat(sprintf("    CHM: %.3f [%.3f, %.3f]\n", 
              median(r2[,1]), quantile(r2[,1], 0.025), quantile(r2[,1], 0.975)))
  r2
}, error = function(e) {
  cat(sprintf("    ⚠ CHM Bayes R² failed: %s\n", e$message))
  NULL
})

br2_dtm <- tryCatch({
  r2 <- bayes_R2(model_dtm_final)
  cat(sprintf("    DTM: %.3f [%.3f, %.3f]\n", 
              median(r2[,1]), quantile(r2[,1], 0.025), quantile(r2[,1], 0.975)))
  r2
}, error = function(e) {
  cat(sprintf("    ⚠ DTM Bayes R² failed: %s\n", e$message))
  NULL
})

# Squared correlation (for comparison)
cat("\n  Squared correlation (point estimate only):\n")

compute_cor2 <- function(model, response_var, name) {
  tryCatch({
    y_obs <- model$data[[response_var]]
    y_pred <- colMeans(posterior_predict(model, ndraws = 100))
    cor2 <- cor(y_obs, y_pred)^2
    cat(sprintf("    %s: %.3f\n", name, cor2))
    return(cor2)
  }, error = function(e) {
    cat(sprintf("    ⚠ %s cor² failed: %s\n", name, e$message))
    return(NA_real_)
  })
}

cor2_chm <- compute_cor2(model_chm_final, "chm_error_mean", "CHM")
cor2_dtm <- compute_cor2(model_dtm_final, "dtm_error_mean", "DTM")

# =====================================================================
# Variance Partitioning (Fixed vs Random Effects)
# =====================================================================

cat("\nVariance partitioning (Fixed vs Random effects):\n")

compute_variance_partition <- function(model, name) {
  tryCatch({
    # Full model R² (includes random effects)
    r2_full <- bayes_R2(model, re_formula = NULL)
    
    # Fixed effects only R²
    r2_fixed <- bayes_R2(model, re_formula = NA)
    
    full_median <- median(r2_full[,1])
    fixed_median <- median(r2_fixed[,1])
    random_contrib <- full_median - fixed_median
    
    cat(sprintf("  %s:\n", name))
    cat(sprintf("    Fixed effects R²:  %.3f (%.1f%% of total)\n", 
                fixed_median, 100 * fixed_median / full_median))
    cat(sprintf("    Random effects:    +%.3f (%.1f%% of explained variance)\n",
                random_contrib, 100 * random_contrib / full_median))
    cat(sprintf("    Full model R²:     %.3f\n", full_median))
    
    return(list(
      r2_fixed = fixed_median,
      r2_full = full_median,
      r2_random_contrib = random_contrib,
      prop_fixed = fixed_median / full_median,
      prop_random = random_contrib / full_median
    ))
    
  }, error = function(e) {
    cat(sprintf("  ⚠ %s variance partition failed: %s\n", name, e$message))
    return(NULL)
  })
}

var_part_chm <- compute_variance_partition(model_chm_final, "CHM")
var_part_dtm <- compute_variance_partition(model_dtm_final, "DTM")

# =====================================================================
# Export Results
# =====================================================================

cat("\nExporting results...\n")

# Create comprehensive R² comparison table
r2_comparison <- data.frame(
  model = c("CHM", "DTM"),
  bayes_r2_median = c(
    ifelse(!is.null(br2_chm), median(br2_chm[,1]), NA),
    ifelse(!is.null(br2_dtm), median(br2_dtm[,1]), NA)
  ),
  bayes_r2_lower = c(
    ifelse(!is.null(br2_chm), quantile(br2_chm[,1], 0.025), NA),
    ifelse(!is.null(br2_dtm), quantile(br2_dtm[,1], 0.025), NA)
  ),
  bayes_r2_upper = c(
    ifelse(!is.null(br2_chm), quantile(br2_chm[,1], 0.975), NA),
    ifelse(!is.null(br2_dtm), quantile(br2_dtm[,1], 0.975), NA)
  ),
  cor_squared = c(cor2_chm, cor2_dtm),
  ppc_coverage_95 = c(cov_chm, cov_dtm),
  r2_fixed_only = c(
    ifelse(!is.null(var_part_chm), var_part_chm$r2_fixed, NA),
    ifelse(!is.null(var_part_dtm), var_part_dtm$r2_fixed, NA)
  ),
  prop_variance_from_random = c(
    ifelse(!is.null(var_part_chm), var_part_chm$prop_random, NA),
    ifelse(!is.null(var_part_dtm), var_part_dtm$prop_random, NA)
  )
)

write.csv(r2_comparison, file.path(out_tables, "model_r2_comparison.csv"), row.names = FALSE)
cat(sprintf("  ✓ Saved to %s\n", file.path(out_tables, "model_r2_comparison.csv")))

# Print summary
cat("\n")
cat(paste(rep("=", 70), collapse = ""))
cat("\n")
cat("SUMMARY\n")
cat(paste(rep("=", 70), collapse = ""))
cat("\n\n")

cat("R² Metrics (use Bayes R² for reporting):\n")
cat(sprintf("  CHM Bayes R²: %.3f   (cor² = %.3f)\n", 
            ifelse(!is.null(br2_chm), median(br2_chm[,1]), NA), cor2_chm))
cat(sprintf("  DTM Bayes R²: %.3f   (cor² = %.3f)\n", 
            ifelse(!is.null(br2_dtm), median(br2_dtm[,1]), NA), cor2_dtm))

cat("\n95% PPC Coverage (target = 95%):\n")
cat(sprintf("  CHM: %.1f%%\n", cov_chm * 100))
cat(sprintf("  DTM: %.1f%%\n", cov_dtm * 100))

cat("\nVariance Partitioning:\n")
if (!is.null(var_part_chm)) {
  cat(sprintf("  CHM: %.1f%% from fixed effects, %.1f%% from random effects\n",
              var_part_chm$prop_fixed * 100, var_part_chm$prop_random * 100))
}
if (!is.null(var_part_dtm)) {
  cat(sprintf("  DTM: %.1f%% from fixed effects, %.1f%% from random effects\n",
              var_part_dtm$prop_fixed * 100, var_part_dtm$prop_random * 100))
}

cat("\n✓ PPC coverage and R² metrics complete\n")

# Make variables available in global environment for section_13
cat("\nVariables exported to global environment:\n")
cat("  cov_chm, cov_dtm, br2_chm, br2_dtm\n")
