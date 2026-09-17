# =====================================================================
# section_09_models_stage1.R (ENHANCED VERSION)
# Bayesian hierarchical models - Stage 1 (Student-t)
#
# ENHANCEMENTS ADDED:
#   1. Prior sensitivity analysis - tests robustness to prior choices
#   2. K-fold cross-validation setup (optional, computationally expensive)
#   3. Enhanced MCMC diagnostics - trace plots, divergences, pairs
#   4. LOO-CV computation for model comparison
# =====================================================================

# Ensure dependencies loaded
if (!exists("log_progress")) source("analysis_utils.R")
if (!exists("PROJECT_ROOT")) source("analysis_config.R")

library(brms)
library(bayesplot)

if (!exists("chm_df") || !exists("dtm_df")) {
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
if (!exists("mod_chm_s1") || !exists("mod_dtm_s1")) {
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

# RAM estimation
estimate_ram_needed <- function(n_rows, n_chains = 4) {
  base_mb <- 500
  per_row_mb <- 0.15
  per_chain_mb <- base_mb + (n_rows * per_row_mb)
  total_gb <- (per_chain_mb * n_chains) / 1024
  return(total_gb)
}

ram_needed_gb <- estimate_ram_needed(nrow(mod_chm_s1), n_chains = 4)
ram_avail_gb <- tryCatch({
  (as.numeric(system("free -g | awk '/^Mem:/ {print $2}'", intern=TRUE))) * RAM_FRAC
}, error = function(e) 16)  # Default to 16GB if can't detect

if (ram_needed_gb > ram_avail_gb * 0.8) {
  warning(sprintf("Estimated RAM needed: %.1f GB, Available: %.1f GB", 
                  ram_needed_gb, ram_avail_gb))
  log_progress("⚠ Consider reducing SUBSAMPLE_FRAC or running with more RAM")
}

log_progress("Fitting Stage 1 Bayesian models (ENHANCED)...")

# =====================================================================
# FORMULA BUILDERS (from original)
# =====================================================================

build_enhanced_formula <- function(response, data, available_meta, 
                                   interaction_term = NULL,
                                   use_smooth = TRUE) {
  base_terms <- c("slope_mean_z", "slope_sd_z", "wsci_z", "rh_98_z", "cover_z",
                  "aspect_sin_z", "aspect_cos_z")
  
  if (use_smooth) {
    smooth_terms <- "s(slope_mean_z, k = 5)"
  } else {
    smooth_terms <- "I(slope_mean_z^2)"
  }
  
  meta_core <- c("meta_offnad_z", "meta_sunel_z", "meta_az_conc_z")
  meta_optional <- c("meta_stereo_z", "meta_leafon_z", "meta_fwdrev_z")
  meta_terms <- intersect(c(meta_core, meta_optional), names(data))
  
  view_terms <- c("view_az_sin_z", "view_az_cos_z")
  
  linear_predictors <- c(base_terms, meta_terms, view_terms)
  predictor_str <- paste(linear_predictors, collapse = " + ")
  predictor_str <- paste(predictor_str, "+", smooth_terms)
  
  if (!is.null(interaction_term)) {
    predictor_str <- paste(predictor_str, "+", interaction_term)
  }
  
  re_str <- "(1 | site) + (1 | gr(ecoregion, by = eco_on)) + (1 | lc_l1_code)"
  formula_str <- paste(response, "~", predictor_str, "+", re_str)
  
  sigma_predictors <- if (response == "chm_error_mean") "slope_mean_z" else "wsci_z"
  sigma_str <- paste0("sigma ~ 1 + ", sigma_predictors, " + (1 | site)")
  
  full_formula_str <- paste0("bf(", formula_str, ", ", sigma_str, ")")
  eval(parse(text = full_formula_str))
}

log_progress(sprintf("  CHM available metadata: %s", paste(available_meta_chm, collapse=", ")))
log_progress(sprintf("  DTM available metadata: %s", paste(available_meta_dtm, collapse=", ")))

# Build formulas
chm_formula_s1 <- build_enhanced_formula(
  response = "chm_error_mean",
  data = mod_chm_s1,
  available_meta = available_meta_chm,
  interaction_term = "slope_mean_z:lc_l1_code",
  use_smooth = TRUE
)

dtm_formula_s1 <- build_enhanced_formula(
  response = "dtm_error_mean",
  data = mod_dtm_s1,
  available_meta = available_meta_dtm,
  interaction_term = "wsci_z:lc_l1_code",
  use_smooth = TRUE
)

# =====================================================================
# DEFAULT PRIORS
# =====================================================================

priors_s1_default <- c(
  prior("normal(0, 2)", class = "b", dpar = "mu"),
  prior("normal(0, 10)", class = "Intercept", dpar = "mu"),
  prior("gamma(2, 0.1)", class = "nu"),
  prior("normal(0, 2)", class = "sds"),
  prior("student_t(3, 0, 5)", class = "Intercept", dpar = "sigma"),
  prior("normal(0, 1)", class = "b", dpar = "sigma"),
  prior("student_t(3, 0, 2)", class = "sd")
)

ctrl_s1 <- list(
  adapt_delta = 0.95,
  max_treedepth = 14,
  step_size = 0.001
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
log_progress(sprintf("  Iterations: %d warmup + %d sampling = %d total",
                     n_warmup, n_iter - n_warmup, n_iter))

# =====================================================================
# FIT MAIN MODELS
# =====================================================================

log_progress("→ Fitting CHM Stage 1 model...")
fit_chm_s1 <- brm(
  formula = chm_formula_s1,
  data = mod_chm_s1,
  prior = priors_s1_default,
  family = student(),
  warmup = n_warmup,   
  iter = n_iter,       
  chains = n_chains,
  cores = parallel_chains,
  threads = threading(threads_per_chain),
  control = ctrl_s1,  
  backend = backend_choice,
  file = file.path(MODELS_DIR, "fit_chm_s1_enhanced"),
  seed = 2025
)
log_progress("  ✓ CHM model complete")

# Quick convergence check
max_rhat_chm <- max(rhat(fit_chm_s1), na.rm = TRUE)
min_ess_chm <- min(neff_ratio(fit_chm_s1), na.rm = TRUE)
log_progress(sprintf("  Max Rhat: %.3f, Min ESS ratio: %.3f", max_rhat_chm, min_ess_chm))
if (max_rhat_chm > 1.01) {
  warning("⚠ Some parameters have Rhat > 1.01. Check convergence diagnostics.")
}

log_progress("→ Fitting DTM Stage 1 model...")
fit_dtm_s1 <- brm(
  formula = dtm_formula_s1,
  data = mod_dtm_s1,
  prior = priors_s1_default,
  family = student(),
  warmup = n_warmup,
  iter = n_iter,
  chains = n_chains,
  cores = parallel_chains,
  threads = threading(threads_per_chain),
  control = ctrl_s1,
  backend = backend_choice,
  file = file.path(MODELS_DIR, "fit_dtm_s1_enhanced"),
  seed = 2025
)
log_progress("  ✓ DTM model complete")

max_rhat_dtm <- max(rhat(fit_dtm_s1), na.rm = TRUE)
min_ess_dtm <- min(neff_ratio(fit_dtm_s1), na.rm = TRUE)
log_progress(sprintf("  Max Rhat: %.3f, Min ESS ratio: %.3f", max_rhat_dtm, min_ess_dtm))

# =====================================================================
# NEW: ENHANCED MCMC DIAGNOSTICS
# =====================================================================

log_subsection("Enhanced MCMC diagnostics (NEW)")

# Check for divergences
np_chm <- nuts_params(fit_chm_s1)
np_dtm <- nuts_params(fit_dtm_s1)

n_div_chm <- sum(np_chm$Value[np_chm$Parameter == "divergent__"])
n_div_dtm <- sum(np_dtm$Value[np_dtm$Parameter == "divergent__"])

log_progress(sprintf("  CHM divergent transitions: %d", n_div_chm))
log_progress(sprintf("  DTM divergent transitions: %d", n_div_dtm))

if (n_div_chm > 0 || n_div_dtm > 0) {
  log_progress("  ⚠ Divergent transitions detected - consider increasing adapt_delta")
}

# Check max treedepth
n_maxtd_chm <- sum(np_chm$Value[np_chm$Parameter == "treedepth__"] >= ctrl_s1$max_treedepth)
n_maxtd_dtm <- sum(np_dtm$Value[np_dtm$Parameter == "treedepth__"] >= ctrl_s1$max_treedepth)

if (n_maxtd_chm > 0 || n_maxtd_dtm > 0) {
  log_progress(sprintf("  ⚠ Iterations at max treedepth: CHM=%d, DTM=%d", n_maxtd_chm, n_maxtd_dtm))
}

# Trace plots for key parameters
log_progress("  Creating trace plots...")
key_params <- c("b_slope_mean_z", "b_wsci_z", "nu")

tryCatch({
  p_trace_chm <- mcmc_trace(fit_chm_s1, pars = key_params, np = np_chm) +
    ggtitle("CHM Stage 1: MCMC Trace Plots")
  ggsave(file.path(out_plots, "06_trace_chm_s1.pdf"), p_trace_chm,
         width = 10, height = 8, bg = "white")
  log_progress("    ✓ CHM trace plots saved")
}, error = function(e) {
  log_progress(sprintf("    ⚠ Trace plot failed: %s", e$message))
})

# Pairs plot to check for funnel shapes
log_progress("  Creating pairs plots...")
tryCatch({
  p_pairs <- mcmc_pairs(fit_chm_s1, pars = c("b_slope_mean_z", "b_wsci_z", "sigma"),
                        off_diag_fun = "hex", np = np_chm)
  ggsave(file.path(out_plots, "06_pairs_chm_s1.pdf"), p_pairs,
         width = 10, height = 10, bg = "white")
  log_progress("    ✓ Pairs plots saved")
}, error = function(e) {
  log_progress(sprintf("    ⚠ Pairs plot failed: %s", e$message))
})

# =====================================================================
# NEW: LOO-CV FOR MODEL COMPARISON
# =====================================================================

log_subsection("LOO-CV computation (NEW)")

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

loo_chm_s1 <- compute_loo_safe(fit_chm_s1, "CHM Stage 1")
loo_dtm_s1 <- compute_loo_safe(fit_dtm_s1, "DTM Stage 1")

# Report LOO results (only if successful)
if (!is.null(loo_chm_s1)) {
  log_progress(sprintf("  CHM LOO ELPD: %.1f (SE: %.1f)",
                       loo_chm_s1$estimates["elpd_loo", "Estimate"],
                       loo_chm_s1$estimates["elpd_loo", "SE"]))
  pk_chm <- loo_chm_s1$diagnostics$pareto_k
  log_progress(sprintf("  CHM Pareto k > 0.7: %d (%.2f%%)", 
                       sum(pk_chm > 0.7), 100 * mean(pk_chm > 0.7)))
} else {
  log_progress("  CHM LOO: SKIPPED (computation failed)")
  pk_chm <- NULL
}

if (!is.null(loo_dtm_s1)) {
  log_progress(sprintf("  DTM LOO ELPD: %.1f (SE: %.1f)",
                       loo_dtm_s1$estimates["elpd_loo", "Estimate"],
                       loo_dtm_s1$estimates["elpd_loo", "SE"]))
  pk_dtm <- loo_dtm_s1$diagnostics$pareto_k
  log_progress(sprintf("  DTM Pareto k > 0.7: %d (%.2f%%)", 
                       sum(pk_dtm > 0.7), 100 * mean(pk_dtm > 0.7)))
} else {
  log_progress("  DTM LOO: SKIPPED (computation failed)")
  pk_dtm <- NULL
}

# Export LOO summary (only for successful computations)
loo_rows <- list()
if (!is.null(loo_chm_s1)) {
  loo_rows$chm <- tibble::tibble(
    model = "CHM_S1",
    elpd_loo = loo_chm_s1$estimates["elpd_loo", "Estimate"],
    se_elpd = loo_chm_s1$estimates["elpd_loo", "SE"],
    p_loo = loo_chm_s1$estimates["p_loo", "Estimate"],
    n_high_pareto_k = sum(pk_chm > 0.7)
  )
}
if (!is.null(loo_dtm_s1)) {
  loo_rows$dtm <- tibble::tibble(
    model = "DTM_S1",
    elpd_loo = loo_dtm_s1$estimates["elpd_loo", "Estimate"],
    se_elpd = loo_dtm_s1$estimates["elpd_loo", "SE"],
    p_loo = loo_dtm_s1$estimates["p_loo", "Estimate"],
    n_high_pareto_k = sum(pk_dtm > 0.7)
  )
}

if (length(loo_rows) > 0) {
  loo_summary <- bind_rows(loo_rows)
  write_csv(loo_summary, file.path(out_tables, "loo_summary_stage1.csv"))
  log_progress("  ✓ LOO summary exported")
} else {
  log_progress("  ⚠ No LOO results to export")
  loo_summary <- NULL
}

# =====================================================================
# NEW: PRIOR SENSITIVITY ANALYSIS
# =====================================================================

# Only run if explicitly requested (computationally expensive)
run_prior_sensitivity <- as.logical(Sys.getenv("RUN_PRIOR_SENSITIVITY", "FALSE"))

if (run_prior_sensitivity) {
  log_subsection("Prior sensitivity analysis (NEW)")
  
  # Define alternative prior specifications
  prior_specs <- list(
    "weakly_informative" = c(
      prior("normal(0, 5)", class = "b", dpar = "mu"),
      prior("normal(0, 20)", class = "Intercept", dpar = "mu"),
      prior("gamma(2, 0.1)", class = "nu"),
      prior("normal(0, 3)", class = "sds"),
      prior("student_t(3, 0, 5)", class = "Intercept", dpar = "sigma"),
      prior("normal(0, 2)", class = "b", dpar = "sigma"),
      prior("student_t(3, 0, 3)", class = "sd")
    ),
    "more_regularizing" = c(
      prior("normal(0, 1)", class = "b", dpar = "mu"),
      prior("normal(0, 5)", class = "Intercept", dpar = "mu"),
      prior("gamma(2, 0.1)", class = "nu"),
      prior("normal(0, 1)", class = "sds"),
      prior("student_t(3, 0, 5)", class = "Intercept", dpar = "sigma"),
      prior("normal(0, 0.5)", class = "b", dpar = "sigma"),
      prior("student_t(3, 0, 1)", class = "sd")
    )
  )
  
  # Use smaller subset for sensitivity analysis
  sens_data <- mod_chm_s1 %>% sample_n(min(3000, n()))
  
  # Fit default model on subset for fair comparison
  log_progress("  Fitting default priors on subset...")
  fit_default_sens <- brm(
    chm_error_mean ~ slope_mean_z + wsci_z + rh_98_z + (1 | site),
    data = sens_data,
    family = student(),
    prior = c(
      prior("normal(0, 2)", class = "b"),
      prior("normal(0, 10)", class = "Intercept"),
      prior("gamma(2, 0.1)", class = "nu"),
      prior("student_t(3, 0, 2)", class = "sd")
    ),
    chains = 2, iter = 1500, warmup = 750,
    silent = 2, refresh = 0, seed = 2025
  )
  
  sens_results <- list()
  sens_results[["default"]] <- fixef(fit_default_sens)[, c("Estimate", "Est.Error")]
  
  for (prior_name in names(prior_specs)) {
    log_progress(sprintf("  Testing %s priors...", prior_name))
    
    sens_fit <- tryCatch({
      brm(
        chm_error_mean ~ slope_mean_z + wsci_z + rh_98_z + (1 | site),
        data = sens_data,
        family = student(),
        prior = c(
          prior_specs[[prior_name]][1:2],  # b and Intercept for mu
          prior_specs[[prior_name]][3],     # nu
          prior_specs[[prior_name]][7]      # sd
        ),
        chains = 2, iter = 1500, warmup = 750,
        silent = 2, refresh = 0, seed = 2025
      )
    }, error = function(e) {
      log_progress(sprintf("    ⚠ Failed: %s", e$message))
      NULL
    })
    
    if (!is.null(sens_fit)) {
      sens_results[[prior_name]] <- fixef(sens_fit)[, c("Estimate", "Est.Error")]
    }
  }
  
  # Compare coefficient estimates across priors
  if (length(sens_results) > 1) {
    log_progress("  Prior sensitivity summary:")
    
    # Create comparison table
    params <- rownames(sens_results[["default"]])
    sens_comparison <- tibble::tibble(parameter = params)
    
    for (prior_name in names(sens_results)) {
      sens_comparison[[paste0(prior_name, "_est")]] <- sens_results[[prior_name]][, "Estimate"]
      sens_comparison[[paste0(prior_name, "_se")]] <- sens_results[[prior_name]][, "Est.Error"]
    }
    
    write_csv(sens_comparison, file.path(out_tables, "prior_sensitivity_comparison.csv"))
    
    # Check for sensitivity (large changes across priors)
    for (param in params[-1]) {  # Skip intercept
      estimates <- sapply(sens_results, function(x) x[param, "Estimate"])
      range_est <- max(estimates) - min(estimates)
      mean_se <- mean(sapply(sens_results, function(x) x[param, "Est.Error"]))
      
      if (range_est > 2 * mean_se) {
        log_progress(sprintf("    ⚠ %s sensitive to priors (range=%.2f, avg SE=%.2f)", 
                             param, range_est, mean_se))
      }
    }
    
    log_progress("  ✓ Prior sensitivity analysis complete")
  }
} else {
  log_progress("  Prior sensitivity analysis: SKIPPED (set RUN_PRIOR_SENSITIVITY=TRUE to run)")
}

# =====================================================================
# NEW: K-FOLD CROSS-VALIDATION (Optional)
# =====================================================================

run_kfold <- as.logical(Sys.getenv("RUN_KFOLD", "FALSE"))

if (run_kfold) {
  log_subsection("K-fold cross-validation (NEW)")
  
  log_progress("  Running 5-fold CV for CHM (this will take a while)...")
  kfold_chm <- kfold(fit_chm_s1, K = 5, save_fits = TRUE)
  
  log_progress("  Running 5-fold CV for DTM...")
  kfold_dtm <- kfold(fit_dtm_s1, K = 5, save_fits = TRUE)
  
  log_progress(sprintf("  CHM k-fold ELPD: %.1f (SE: %.1f)",
                       kfold_chm$estimates["elpd_kfold", "Estimate"],
                       kfold_chm$estimates["elpd_kfold", "SE"]))
  log_progress(sprintf("  DTM k-fold ELPD: %.1f (SE: %.1f)",
                       kfold_dtm$estimates["elpd_kfold", "Estimate"],
                       kfold_dtm$estimates["elpd_kfold", "SE"]))
  
  # Save k-fold results
  kfold_summary <- tibble::tibble(
    model = c("CHM_S1", "DTM_S1"),
    elpd_kfold = c(kfold_chm$estimates["elpd_kfold", "Estimate"],
                   kfold_dtm$estimates["elpd_kfold", "Estimate"]),
    se_elpd_kfold = c(kfold_chm$estimates["elpd_kfold", "SE"],
                      kfold_dtm$estimates["elpd_kfold", "SE"])
  )
  write_csv(kfold_summary, file.path(out_tables, "kfold_summary_stage1.csv"))
  
} else {
  log_progress("  K-fold CV: SKIPPED (set RUN_KFOLD=TRUE to run)")
  kfold_chm <- NULL
  kfold_dtm <- NULL
}

log_progress("✓ Stage 1 models complete (ENHANCED)")

# =====================================================================
# SAVE CHECKPOINT
# =====================================================================

if (!exists("BATCH_MODE") || !BATCH_MODE) {
  log_progress("Saving checkpoint for interactive mode...")
  
  checkpoint_data <- list(
    fit_chm_s1 = fit_chm_s1,
    fit_dtm_s1 = fit_dtm_s1,
    chm_formula_s1 = chm_formula_s1,
    dtm_formula_s1 = dtm_formula_s1,
    loo_chm_s1 = loo_chm_s1,
    loo_dtm_s1 = loo_dtm_s1,
    n_divergences_chm = n_div_chm,
    n_divergences_dtm = n_div_dtm
  )
  
  if (exists("kfold_chm") && !is.null(kfold_chm)) {
    checkpoint_data$kfold_chm <- kfold_chm
    checkpoint_data$kfold_dtm <- kfold_dtm
  }
  
  if (exists("sens_results")) {
    checkpoint_data$sens_results <- sens_results
  }
  
  save_checkpoint("09_models_stage1", checkpoint_data)
}
