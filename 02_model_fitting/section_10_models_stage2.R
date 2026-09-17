# =====================================================================
# section_10_models_stage2.R (ENHANCED VERSION)
# Bayesian hierarchical models - Stage 2 (Enhanced with Student-t)
#
# ENHANCEMENTS ADDED:
#   1. Formal model comparison (LOO-CV) - Stage 1 vs Stage 2
#   2. WAIC computation for additional comparison metric
#   3. Alternative RE structure testing (optional)
#   4. Model complexity assessment
# =====================================================================

# Ensure dependencies loaded
if (!exists("log_progress")) source("analysis_utils.R")
if (!exists("PROJECT_ROOT")) source("analysis_config.R")

library(brms)
library(bayesplot)
library(loo)

# Load data and metadata variables
if (!exists("chm_df") || !exists("dtm_df") || 
    !exists("available_meta_chm") || !exists("available_meta_dtm")) {
  log_progress("⚠ Data not loaded. Loading from checkpoint...")
  if (checkpoint_exists("01_data_ingest")) {
    data <- load_checkpoint("01_data_ingest")
    chm_df <- data$chm_df
    dtm_df <- data$dtm_df
    available_meta_chm <- data$available_meta_chm
    available_meta_dtm <- data$available_meta_dtm
  } else {
    stop("Data not available. Run: source('section_01_ingest.R') first")
  }
}

# Load model data
if (!exists("mod_chm_s2") || !exists("mod_dtm_s2")) {
  log_progress("⚠ Model data not loaded. Loading from checkpoint...")
  if (checkpoint_exists("08_model_prep")) {
    data <- load_checkpoint("08_model_prep")
    mod_chm_s1 <- data$mod_chm_s1
    mod_dtm_s1 <- data$mod_dtm_s1
    mod_chm_s2 <- data$mod_chm_s2
    mod_dtm_s2 <- data$mod_dtm_s2
  } else {
    stop("Model data not available. Run: source('section_08_model_prep.R') first")
  }
}

# Load Stage 1 models and LOO results
if (!exists("fit_chm_s1") || !exists("fit_dtm_s1")) {
  log_progress("⚠ Stage 1 models not loaded. Loading from checkpoint...")
  if (checkpoint_exists("09_models_stage1")) {
    data <- load_checkpoint("09_models_stage1")
    fit_chm_s1 <- data$fit_chm_s1
    fit_dtm_s1 <- data$fit_dtm_s1
    loo_chm_s1 <- data$loo_chm_s1
    loo_dtm_s1 <- data$loo_dtm_s1
  } else if (checkpoint_exists("09_models_s1")) {
    data <- load_checkpoint("09_models_s1")
    fit_chm_s1 <- data$fit_chm_s1
    fit_dtm_s1 <- data$fit_dtm_s1
    loo_chm_s1 <- data$loo_chm_s1
    loo_dtm_s1 <- data$loo_dtm_s1
  } else {
    warning("Stage 1 models not found. Proceeding without comparison.")
    fit_chm_s1 <- NULL
    fit_dtm_s1 <- NULL
    loo_chm_s1 <- NULL
    loo_dtm_s1 <- NULL
  }
}

# RAM estimation
estimate_ram_needed <- function(n_rows, n_chains = 4) {
  base_mb <- 500
  per_row_mb <- 0.15
  per_chain_mb <- base_mb + (n_rows * per_row_mb)
  total_gb <- (per_chain_mb * n_chains) / 1024
  return(total_gb)
}

ram_needed_gb <- estimate_ram_needed(nrow(mod_chm_s2), n_chains = 4)
ram_avail_gb <- tryCatch({
  (as.numeric(system("free -g | awk '/^Mem:/ {print $2}'", intern=TRUE))) * RAM_FRAC
}, error = function(e) 16)

if (ram_needed_gb > ram_avail_gb * 0.8) {
  warning(sprintf("Estimated RAM needed: %.1f GB, Available: %.1f GB", 
                  ram_needed_gb, ram_avail_gb))
  log_progress("⚠ Consider reducing SUBSAMPLE_FRAC or running with more RAM")
}

log_progress("Fitting Stage 2 Bayesian models (ENHANCED)...")

# =====================================================================
# STAGE 2 FORMULA BUILDER
# =====================================================================

build_stage2_formula <- function(response, data, available_meta, 
                                 use_random_slopes = TRUE) {
  
  base_terms <- c("slope_mean_z", "slope_sd_z", "wsci_z", "rh_98_z", "cover_z",
                  "aspect_sin_z", "aspect_cos_z")
  
  smooth_terms <- "s(slope_mean_z, wsci_z, k = 20)"
  
  meta_core <- c("meta_offnad_z", "meta_sunel_z", "meta_az_conc_z")
  meta_optional <- c("meta_stereo_z", "meta_leafon_z", "meta_fwdrev_z",
                     "meta_absgeo_z", "meta_relgeo_z")
  meta_terms <- intersect(c(meta_core, meta_optional), names(data))
  
  view_terms <- c("view_az_sin_z", "view_az_cos_z")
  
  if (response == "chm_error_mean") {
    interaction_terms <- c("slope_mean_z:lc_l1_code", "wsci_z:lc_l1_code", "slope_mean_z:wsci_z")
  } else {
    interaction_terms <- c("wsci_z:lc_l1_code", "slope_mean_z:lc_l1_code", "slope_mean_z:wsci_z")
  }
  
  all_predictors <- c(base_terms, smooth_terms, meta_terms, view_terms, interaction_terms)
  predictor_str <- paste(all_predictors, collapse = " + ")
  
  if (use_random_slopes) {
    re_str <- paste0(
      "(1 + slope_mean_z | ecoregion) + ",
      "(1 | site) + ",
      "(1 + wsci_z | lc_l1_code)"
    )
    log_progress("  Using random slopes (Stage 1 converged well)")
  } else {
    re_str <- "(1 | site) + (1 | ecoregion) + (1 | lc_l1_code)"
    log_progress("  Using random intercepts only")
  }
  
  formula_str <- paste(response, "~", predictor_str, "+", re_str)
  sigma_pred <- "slope_mean_z"
  sigma_str <- paste0("sigma ~ 1 + ", sigma_pred, " + (1 | site) + (1 | lc_l1_code)")
  full_formula_str <- paste0("bf(", formula_str, ", ", sigma_str, ")")
  
  eval(parse(text = full_formula_str))
}

# =====================================================================
# CHECK STAGE 1 CONVERGENCE
# =====================================================================

s1_converged_well <- function(fit) {
  if (is.null(fit)) return(FALSE)
  
  tryCatch({
    max_rhat <- max(rhat(fit), na.rm = TRUE)
    min_ess <- min(neff_ratio(fit), na.rm = TRUE)
    converged <- (max_rhat < 1.01) && (min_ess > 0.1)
    
    if (converged) {
      log_progress(sprintf("  Stage 1 convergence: ✓ (Rhat=%.3f, ESS ratio=%.3f)", 
                           max_rhat, min_ess))
    } else {
      log_progress(sprintf("  Stage 1 convergence: ⚠ (Rhat=%.3f, ESS ratio=%.3f)", 
                           max_rhat, min_ess))
    }
    
    return(converged)
  }, error = function(e) {
    log_progress("  ⚠ Could not assess Stage 1 convergence")
    return(FALSE)
  })
}

log_progress("Assessing Stage 1 convergence...")
chm_s1_good <- s1_converged_well(fit_chm_s1)
dtm_s1_good <- s1_converged_well(fit_dtm_s1)
use_random_slopes <- chm_s1_good && dtm_s1_good

if (use_random_slopes) {
  log_progress("✓ Both Stage 1 models converged well → Using random slopes in Stage 2")
} else {
  log_progress("⚠ Stage 1 convergence issues → Using simpler Stage 2")
}

# =====================================================================
# BUILD STAGE 2 FORMULAS
# =====================================================================

log_progress("Building Stage 2 model formulas...")
log_progress(sprintf("  CHM available metadata: %s", paste(available_meta_chm, collapse=", ")))
log_progress(sprintf("  DTM available metadata: %s", paste(available_meta_dtm, collapse=", ")))

chm_formula_s2 <- build_stage2_formula(
  "chm_error_mean", 
  mod_chm_s2, 
  available_meta_chm, 
  use_random_slopes = use_random_slopes
)

dtm_formula_s2 <- build_stage2_formula(
  "dtm_error_mean", 
  mod_dtm_s2, 
  available_meta_dtm,
  use_random_slopes = use_random_slopes
)

# =====================================================================
# STAGE 2 PRIORS AND CONTROL
# =====================================================================

priors_s2 <- c(
  prior("normal(0, 2)", class = "b", dpar = "mu"),
  prior("normal(0, 10)", class = "Intercept", dpar = "mu"),
  prior("gamma(2, 0.1)", class = "nu"),
  prior("normal(0, 1.5)", class = "sds"),
  prior("student_t(3, 0, 5)", class = "Intercept", dpar = "sigma"),
  prior("normal(0, 1)", class = "b", dpar = "sigma"),
  prior("student_t(3, 0, 2)", class = "sd"),
  prior("lkj(2)", class = "cor")
)

ctrl_s2 <- list(
  adapt_delta = 0.98,
  max_treedepth = 15,
  step_size = 0.0005
)

# Sampling parameters
n_chains <- 4
parallel_chains <- max(1L, min(n_chains, CPU_BUDGET))
threads_per_chain <- if (CPU_BUDGET >= 8) 2L else 1L
n_warmup <- 1500
n_iter <- 3000

use_cmdstanr <- requireNamespace("cmdstanr", quietly = TRUE)
backend_choice <- if (use_cmdstanr) "cmdstanr" else "rstan"

log_progress(sprintf("  Using backend: %s", backend_choice))
log_progress(sprintf("  Chains: %d parallel, %d thread(s) per chain", 
                     parallel_chains, threads_per_chain))
log_progress(sprintf("  Control: adapt_delta=%.2f, max_treedepth=%d",
                     ctrl_s2$adapt_delta, ctrl_s2$max_treedepth))

# =====================================================================
# FIT STAGE 2 MODELS
# =====================================================================

log_progress("→ Fitting CHM Stage 2 model...")
fit_chm_s2 <- brm(
  formula = chm_formula_s2,
  data = mod_chm_s2,
  prior = priors_s2,
  family = student(),
  warmup = n_warmup,
  iter = n_iter,
  chains = n_chains,
  cores = parallel_chains,
  threads = threading(threads_per_chain),
  control = ctrl_s2,
  backend = backend_choice,
  file = file.path(MODELS_DIR, "fit_chm_s2_enhanced"),
  seed = 2025
)
log_progress("  ✓ CHM Stage 2 model complete")

max_rhat <- max(rhat(fit_chm_s2), na.rm = TRUE)
min_ess <- min(neff_ratio(fit_chm_s2), na.rm = TRUE)
log_progress(sprintf("  Max Rhat: %.3f, Min ESS ratio: %.3f", max_rhat, min_ess))
if (max_rhat > 1.01) {
  warning("⚠ Some CHM parameters have Rhat > 1.01")
}

log_progress("→ Fitting DTM Stage 2 model...")
fit_dtm_s2 <- brm(
  formula = dtm_formula_s2,
  data = mod_dtm_s2,
  prior = priors_s2,
  family = student(),
  warmup = n_warmup,
  iter = n_iter,
  chains = n_chains,
  cores = parallel_chains,
  threads = threading(threads_per_chain),
  control = ctrl_s2,
  backend = backend_choice,
  file = file.path(MODELS_DIR, "fit_dtm_s2_enhanced"),
  seed = 2025
)
log_progress("  ✓ DTM Stage 2 model complete")

max_rhat <- max(rhat(fit_dtm_s2), na.rm = TRUE)
min_ess <- min(neff_ratio(fit_dtm_s2), na.rm = TRUE)
log_progress(sprintf("  Max Rhat: %.3f, Min ESS ratio: %.3f", max_rhat, min_ess))

# =====================================================================
# NEW: LOO-CV FOR STAGE 2
# =====================================================================

log_subsection("LOO-CV computation for Stage 2 (NEW)")

# Robust LOO computation with multiple fallbacks
compute_loo_safe <- function(fit, name) {
  tryCatch({
    log_progress(sprintf("  Computing LOO for %s (with moment matching)...", name))
    loo(fit, moment_match = TRUE)
  }, error = function(e) {
    log_progress(sprintf("    ⚠ Moment matching failed: %s", e$message))
    tryCatch({
      log_progress(sprintf("    Trying basic LOO for %s...", name))
      loo(fit)
    }, error = function(e2) {
      log_progress(sprintf("    ⚠ LOO-CV failed entirely for %s: %s", name, e2$message))
      log_progress("    Continuing without LOO for this model")
      NULL
    })
  })
}

loo_chm_s2 <- compute_loo_safe(fit_chm_s2, "CHM Stage 2")
loo_dtm_s2 <- compute_loo_safe(fit_dtm_s2, "DTM Stage 2")

if (!is.null(loo_chm_s2)) {
  log_progress(sprintf("  CHM S2 LOO ELPD: %.1f (SE: %.1f)",
                       loo_chm_s2$estimates["elpd_loo", "Estimate"],
                       loo_chm_s2$estimates["elpd_loo", "SE"]))
} else {
  log_progress("  CHM S2 LOO: SKIPPED (computation failed)")
}

if (!is.null(loo_dtm_s2)) {
  log_progress(sprintf("  DTM S2 LOO ELPD: %.1f (SE: %.1f)",
                       loo_dtm_s2$estimates["elpd_loo", "Estimate"],
                       loo_dtm_s2$estimates["elpd_loo", "SE"]))
} else {
  log_progress("  DTM S2 LOO: SKIPPED (computation failed)")
}

# =====================================================================
# NEW: FORMAL MODEL COMPARISON (Stage 1 vs Stage 2)
# =====================================================================

log_subsection("Model comparison: Stage 1 vs Stage 2 (NEW)")

# Note: Comparison requires models fit on same data
# Since S1 and S2 use different data, we compare relative metrics

comparison_results <- list()

# For CHM - only if both LOO results available
if (!is.null(loo_chm_s1) && !is.null(loo_chm_s2)) {
  log_progress("  CHM comparison:")
  log_progress(sprintf("    Stage 1 ELPD: %.1f (n=%d)", 
                       loo_chm_s1$estimates["elpd_loo", "Estimate"],
                       nrow(mod_chm_s1)))
  log_progress(sprintf("    Stage 2 ELPD: %.1f (n=%d)", 
                       loo_chm_s2$estimates["elpd_loo", "Estimate"],
                       nrow(mod_chm_s2)))
  
  # Normalized ELPD per observation
  elpd_per_obs_s1_chm <- loo_chm_s1$estimates["elpd_loo", "Estimate"] / nrow(mod_chm_s1)
  elpd_per_obs_s2_chm <- loo_chm_s2$estimates["elpd_loo", "Estimate"] / nrow(mod_chm_s2)
  
  log_progress(sprintf("    Stage 1 ELPD/obs: %.4f", elpd_per_obs_s1_chm))
  log_progress(sprintf("    Stage 2 ELPD/obs: %.4f", elpd_per_obs_s2_chm))
  
  if (elpd_per_obs_s2_chm > elpd_per_obs_s1_chm) {
    log_progress("    → Stage 2 has better per-observation predictive accuracy")
  } else {
    log_progress("    → Stage 1 has comparable or better per-observation accuracy")
    log_progress("      Consider if Stage 2 complexity is justified")
  }
  
  comparison_results$chm <- tibble::tibble(
    metric = c("elpd_s1", "elpd_s2", "elpd_per_obs_s1", "elpd_per_obs_s2", "n_s1", "n_s2"),
    value = c(loo_chm_s1$estimates["elpd_loo", "Estimate"],
              loo_chm_s2$estimates["elpd_loo", "Estimate"],
              elpd_per_obs_s1_chm, elpd_per_obs_s2_chm,
              nrow(mod_chm_s1), nrow(mod_chm_s2))
  )
} else {
  log_progress("  CHM comparison: SKIPPED (LOO not available for one or both stages)")
}

# For DTM - only if both LOO results available
if (!is.null(loo_dtm_s1) && !is.null(loo_dtm_s2)) {
  log_progress("  DTM comparison:")
  log_progress(sprintf("    Stage 1 ELPD: %.1f (n=%d)", 
                       loo_dtm_s1$estimates["elpd_loo", "Estimate"],
                       nrow(mod_dtm_s1)))
  log_progress(sprintf("    Stage 2 ELPD: %.1f (n=%d)", 
                       loo_dtm_s2$estimates["elpd_loo", "Estimate"],
                       nrow(mod_dtm_s2)))
  
  elpd_per_obs_s1_dtm <- loo_dtm_s1$estimates["elpd_loo", "Estimate"] / nrow(mod_dtm_s1)
  elpd_per_obs_s2_dtm <- loo_dtm_s2$estimates["elpd_loo", "Estimate"] / nrow(mod_dtm_s2)
  
  log_progress(sprintf("    Stage 1 ELPD/obs: %.4f", elpd_per_obs_s1_dtm))
  log_progress(sprintf("    Stage 2 ELPD/obs: %.4f", elpd_per_obs_s2_dtm))
  
  comparison_results$dtm <- tibble::tibble(
    metric = c("elpd_s1", "elpd_s2", "elpd_per_obs_s1", "elpd_per_obs_s2", "n_s1", "n_s2"),
    value = c(loo_dtm_s1$estimates["elpd_loo", "Estimate"],
              loo_dtm_s2$estimates["elpd_loo", "Estimate"],
              elpd_per_obs_s1_dtm, elpd_per_obs_s2_dtm,
              nrow(mod_dtm_s1), nrow(mod_dtm_s2))
  )
} else {
  log_progress("  DTM comparison: SKIPPED (LOO not available for one or both stages)")
}

# Export comparison results
if (length(comparison_results) > 0) {
  comparison_df <- bind_rows(comparison_results, .id = "product")
  write_csv(comparison_df, file.path(out_tables, "model_comparison_s1_vs_s2.csv"))
  log_progress("  ✓ Model comparison exported")
} else {
  comparison_df <- NULL
  log_progress("  ⚠ No comparison results to export")
}

# =====================================================================
# NEW: WAIC COMPUTATION
# =====================================================================

log_subsection("WAIC computation (NEW)")

log_progress("  Computing WAIC for Stage 2 models...")
waic_chm_s2 <- tryCatch(waic(fit_chm_s2), error = function(e) NULL)
waic_dtm_s2 <- tryCatch(waic(fit_dtm_s2), error = function(e) NULL)

if (!is.null(waic_chm_s2)) {
  log_progress(sprintf("  CHM S2 WAIC: %.1f", waic_chm_s2$estimates["waic", "Estimate"]))
}
if (!is.null(waic_dtm_s2)) {
  log_progress(sprintf("  DTM S2 WAIC: %.1f", waic_dtm_s2$estimates["waic", "Estimate"]))
}

# =====================================================================
# NEW: MODEL COMPLEXITY ASSESSMENT
# =====================================================================

log_subsection("Model complexity assessment (NEW)")

# Count effective parameters from LOO
p_loo_chm <- loo_chm_s2$estimates["p_loo", "Estimate"]
p_loo_dtm <- loo_dtm_s2$estimates["p_loo", "Estimate"]

log_progress(sprintf("  CHM effective parameters (p_loo): %.1f", p_loo_chm))
log_progress(sprintf("  DTM effective parameters (p_loo): %.1f", p_loo_dtm))

# Compare to nominal parameters
n_fixed_params_approx <- 15  # Rough estimate
if (p_loo_chm > 2 * n_fixed_params_approx) {
  log_progress("  → CHM: High p_loo suggests random effects capturing substantial variance")
}

# =====================================================================
# NEW: ALTERNATIVE RE STRUCTURE TEST (Optional)
# =====================================================================

test_alt_re <- as.logical(Sys.getenv("TEST_ALT_RE", "FALSE"))

if (test_alt_re) {
  log_subsection("Alternative RE structure test (NEW)")
  
  # Fit simpler model without random slopes
  log_progress("  Fitting simpler RE model (intercepts only)...")
  
  chm_formula_simple <- bf(
    chm_error_mean ~ slope_mean_z + wsci_z + rh_98_z + cover_z +
      s(slope_mean_z, wsci_z, k = 20) +
      slope_mean_z:lc_l1_code + wsci_z:lc_l1_code +
      (1 | site) + (1 | ecoregion) + (1 | lc_l1_code),
    sigma ~ 1 + slope_mean_z + (1 | site)
  )
  
  fit_chm_simple <- brm(
    formula = chm_formula_simple,
    data = mod_chm_s2 %>% sample_n(min(50000, n())),  # Subset for speed
    prior = priors_s2,
    family = student(),
    warmup = 1000, iter = 2000,
    chains = 2, cores = 2,
    control = ctrl_s2,
    backend = backend_choice,
    seed = 2025,
    silent = 2
  )
  
  loo_chm_simple <- loo(fit_chm_simple)
  
  log_progress(sprintf("  Simple RE ELPD: %.1f", 
                       loo_chm_simple$estimates["elpd_loo", "Estimate"]))
  
  # Note: Direct comparison only valid if same data used
  log_progress("  Note: Compare with full Stage 2 on same subset for valid comparison")
  
} else {
  log_progress("  Alternative RE test: SKIPPED (set TEST_ALT_RE=TRUE to run)")
}

log_progress("✓ Stage 2 models complete (ENHANCED)")

# =====================================================================
# SAVE CHECKPOINT
# =====================================================================

if (!exists("BATCH_MODE") || !BATCH_MODE) {
  log_progress("Saving checkpoint for interactive mode...")
  
  checkpoint_data <- list(
    fit_chm_s2 = fit_chm_s2,
    fit_dtm_s2 = fit_dtm_s2,
    chm_formula_s2 = chm_formula_s2,
    dtm_formula_s2 = dtm_formula_s2,
    use_random_slopes = use_random_slopes
  )
  
  # Only add LOO results if they exist
  if (!is.null(loo_chm_s2)) checkpoint_data$loo_chm_s2 <- loo_chm_s2
  if (!is.null(loo_dtm_s2)) checkpoint_data$loo_dtm_s2 <- loo_dtm_s2
  if (exists("waic_chm_s2") && !is.null(waic_chm_s2)) checkpoint_data$waic_chm_s2 <- waic_chm_s2
  if (exists("waic_dtm_s2") && !is.null(waic_dtm_s2)) checkpoint_data$waic_dtm_s2 <- waic_dtm_s2
  if (exists("comparison_df") && !is.null(comparison_df)) checkpoint_data$comparison_results <- comparison_df
  
  save_checkpoint("10_models_stage2", checkpoint_data)
}
