# =====================================================================
# section_11_diagnostics.R (ENHANCED VERSION)
# Comprehensive model diagnostics
#
# ENHANCEMENTS ADDED:
#   1. Predictive accuracy metrics (RMSE, MAE, correlation)
#   2. Stratified residual analysis (by slope, LC, site)
#   3. Effect size interpretation table
#   4. Probability of direction for effects
# =====================================================================

# Ensure dependencies loaded
if (!exists("log_progress")) source("analysis_utils.R")
if (!exists("PROJECT_ROOT")) source("analysis_config.R")

library(brms)
library(bayesplot)
library(loo)

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

# Load models
if (!exists("fit_chm_s2") && !exists("fit_chm_s1")) {
  log_progress("⚠ Models not loaded. Loading from checkpoint...")
  if (checkpoint_exists("10_models_stage2")) {
    data <- load_checkpoint("10_models_stage2")
    fit_chm_s2 <- data$fit_chm_s2
    fit_dtm_s2 <- data$fit_dtm_s2
    loo_chm_s2 <- data$loo_chm_s2
    loo_dtm_s2 <- data$loo_dtm_s2
    model_chm_final <- fit_chm_s2
    model_dtm_final <- fit_dtm_s2
    loo_chm_final <- loo_chm_s2
    loo_dtm_final <- loo_dtm_s2
  } else if (checkpoint_exists("09_models_stage1")) {
    data <- load_checkpoint("09_models_stage1")
    fit_chm_s1 <- data$fit_chm_s1
    fit_dtm_s1 <- data$fit_dtm_s1
    loo_chm_s1 <- data$loo_chm_s1
    loo_dtm_s1 <- data$loo_dtm_s1
    model_chm_final <- fit_chm_s1
    model_dtm_final <- fit_dtm_s1
    loo_chm_final <- loo_chm_s1
    loo_dtm_final <- loo_dtm_s1
  } else {
    stop("No models available. Run modeling sections first.")
  }
} else {
  model_chm_final <- if (exists("fit_chm_s2")) fit_chm_s2 else fit_chm_s1
  model_dtm_final <- if (exists("fit_dtm_s2")) fit_dtm_s2 else fit_dtm_s1
  loo_chm_final <- if (exists("loo_chm_s2")) loo_chm_s2 else if (exists("loo_chm_s1")) loo_chm_s1 else NULL
  loo_dtm_final <- if (exists("loo_dtm_s2")) loo_dtm_s2 else if (exists("loo_dtm_s1")) loo_dtm_s1 else NULL
}

log_progress("Running enhanced model diagnostics...")

model_stage <- if (exists("fit_chm_s2")) "Stage 2" else "Stage 1"
log_progress(sprintf("  Diagnosing %s models", model_stage))

# Define forest classes for focused analysis
FOREST_CLASSES <- c("EBF", "BDF", "ENF", "DNF")
FOREST_LABELS <- c("Broadleaf Evergreen", "Broadleaf Deciduous", 
                   "Needleleaf Evergreen", "Needleleaf Deciduous")
names(FOREST_LABELS) <- FOREST_CLASSES

# Create forest subsets
chm_forest <- chm_df %>% filter(lc_l1_code %in% FOREST_CLASSES)
dtm_forest <- dtm_df %>% filter(lc_l1_code %in% FOREST_CLASSES)
forest_present_chm <- intersect(FOREST_CLASSES, unique(chm_df$lc_l1_code))

log_progress(sprintf("  Forest classes for analysis: %s", paste(forest_present_chm, collapse = ", ")))

# Check model family
is_student_chm <- family(model_chm_final)$family == "student"
is_student_dtm <- family(model_dtm_final)$family == "student"

# =====================================================================
# 1. CONVERGENCE DIAGNOSTICS (Original + Enhanced)
# =====================================================================

log_subsection("Convergence diagnostics")

rhat_chm <- brms::rhat(model_chm_final)
rhat_dtm <- brms::rhat(model_dtm_final)

p_rhat_chm <- bayesplot::mcmc_rhat(rhat_chm) + 
  ggtitle("R-hat (CHM)") +
  geom_vline(xintercept = 1.01, linetype = "dashed", color = "red")

p_rhat_dtm <- bayesplot::mcmc_rhat(rhat_dtm) + 
  ggtitle("R-hat (DTM)") +
  geom_vline(xintercept = 1.01, linetype = "dashed", color = "red")

ggsave(file.path(out_plots, "06a_rhat_chm.pdf"), p_rhat_chm, 
       width = 7, height = 5, bg = "white")
ggsave(file.path(out_plots, "06a_rhat_dtm.pdf"), p_rhat_dtm, 
       width = 7, height = 5, bg = "white")

log_progress(sprintf("  CHM: max Rhat = %.4f, min ESS ratio = %.3f",
                     max(rhat_chm, na.rm = TRUE),
                     min(neff_ratio(model_chm_final), na.rm = TRUE)))
log_progress(sprintf("  DTM: max Rhat = %.4f, min ESS ratio = %.3f",
                     max(rhat_dtm, na.rm = TRUE),
                     min(neff_ratio(model_dtm_final), na.rm = TRUE)))

# =====================================================================
# 2. STUDENT-T ν DIAGNOSTICS (Original)
# =====================================================================

log_subsection("Student-t degrees of freedom (ν) diagnostics")

nu_summary_chm <- NULL
nu_summary_dtm <- NULL

if (is_student_chm) {
  nu_chm <- as_draws_df(model_chm_final, variable = "nu")$nu
  nu_summary_chm <- c(
    median = median(nu_chm),
    mean = mean(nu_chm),
    q025 = quantile(nu_chm, 0.025),
    q975 = quantile(nu_chm, 0.975)
  )
  log_progress(sprintf("  CHM ν: %.2f [%.2f, %.2f]",
                       nu_summary_chm["median"], nu_summary_chm["q025"], nu_summary_chm["q975"]))
  
  if (nu_summary_chm["median"] < 5) {
    log_progress("    → Heavy tails (ν < 5) - Student-t justified!")
  }
  
  p_nu_chm <- ggplot(data.frame(nu = nu_chm), aes(x = nu)) +
    geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
    geom_vline(xintercept = nu_summary_chm["median"], color = "red", linetype = "dashed") +
    labs(title = "CHM: Posterior distribution of ν",
         subtitle = sprintf("Median = %.2f", nu_summary_chm["median"])) +
    theme_cowplot()
  ggsave(file.path(out_plots, "06d_chm_nu_posterior.pdf"), p_nu_chm,
         width = 7, height = 5, bg = "white")
}

if (is_student_dtm) {
  nu_dtm <- as_draws_df(model_dtm_final, variable = "nu")$nu
  nu_summary_dtm <- c(
    median = median(nu_dtm),
    mean = mean(nu_dtm),
    q025 = quantile(nu_dtm, 0.025),
    q975 = quantile(nu_dtm, 0.975)
  )
  log_progress(sprintf("  DTM ν: %.2f [%.2f, %.2f]",
                       nu_summary_dtm["median"], nu_summary_dtm["q025"], nu_summary_dtm["q975"]))
  
  p_nu_dtm <- ggplot(data.frame(nu = nu_dtm), aes(x = nu)) +
    geom_histogram(bins = 50, fill = "darkorange", alpha = 0.7) +
    geom_vline(xintercept = nu_summary_dtm["median"], color = "red", linetype = "dashed") +
    labs(title = "DTM: Posterior distribution of ν",
         subtitle = sprintf("Median = %.2f", nu_summary_dtm["median"])) +
    theme_cowplot()
  ggsave(file.path(out_plots, "06d_dtm_nu_posterior.pdf"), p_nu_dtm,
         width = 7, height = 5, bg = "white")
}

# =====================================================================
# 3. POSTERIOR PREDICTIVE CHECKS (Original)
# =====================================================================

log_subsection("Posterior predictive checks")

ppc_chm <- pp_check(model_chm_final, ndraws = 200) + 
  ggtitle("CHM: Posterior predictive density overlay")
ppc_dtm <- pp_check(model_dtm_final, ndraws = 200) + 
  ggtitle("DTM: Posterior predictive density overlay")

ggsave(file.path(out_plots, "06b_ppc_chm.pdf"), ppc_chm, width = 6, height = 5, bg = "white")
ggsave(file.path(out_plots, "06b_ppc_dtm.pdf"), ppc_dtm, width = 6, height = 5, bg = "white")

# PPC error histograms (Original)
ppc_err_chm <- tryCatch({
  pp_check(model_chm_final, type = "error_hist", ndraws = 10) +
    ggtitle("CHM: Posterior predictive error histograms")
}, error = function(e) NULL)

ppc_err_dtm <- tryCatch({
  pp_check(model_dtm_final, type = "error_hist", ndraws = 10) +
    ggtitle("DTM: Posterior predictive error histograms")
}, error = function(e) NULL)

if (!is.null(ppc_err_chm)) {
  ggsave(file.path(out_plots, "06b_ppc_error_chm.pdf"), ppc_err_chm,
         width = 10, height = 6, bg = "white")
}
if (!is.null(ppc_err_dtm)) {
  ggsave(file.path(out_plots, "06b_ppc_error_dtm.pdf"), ppc_err_dtm,
         width = 10, height = 6, bg = "white")
}

# =====================================================================
# 3b. RANDOM EFFECTS EXAMINATION (Original)
# =====================================================================

log_subsection("Random effects summary")

# Extract random effects with error handling
ranef_chm <- tryCatch(ranef(model_chm_final), error = function(e) NULL)
ranef_dtm <- tryCatch(ranef(model_dtm_final), error = function(e) NULL)

# Ecoregion effects
if (!is.null(ranef_chm) && "ecoregion" %in% names(ranef_chm)) {
  eco_effects_chm <- ranef_chm$ecoregion[, "Estimate", "Intercept"]
  log_progress(sprintf("  CHM ecoregion effects: range = [%.2f, %.2f]",
                       min(eco_effects_chm), max(eco_effects_chm)))
  
  eco_table_chm <- as.data.frame(ranef_chm$ecoregion)
  eco_table_chm$ecoregion <- rownames(eco_table_chm)
  write_csv(eco_table_chm, file.path(out_tables, "random_effects_ecoregion_chm.csv"))
}

if (!is.null(ranef_dtm) && "ecoregion" %in% names(ranef_dtm)) {
  eco_effects_dtm <- ranef_dtm$ecoregion[, "Estimate", "Intercept"]
  log_progress(sprintf("  DTM ecoregion effects: range = [%.2f, %.2f]",
                       min(eco_effects_dtm), max(eco_effects_dtm)))
  
  eco_table_dtm <- as.data.frame(ranef_dtm$ecoregion)
  eco_table_dtm$ecoregion <- rownames(eco_table_dtm)
  write_csv(eco_table_dtm, file.path(out_tables, "random_effects_ecoregion_dtm.csv"))
}

# Land cover effects
if (!is.null(ranef_chm) && "lc_l1_code" %in% names(ranef_chm)) {
  lc_effects_chm <- ranef_chm$lc_l1_code[, "Estimate", "Intercept"]
  log_progress(sprintf("  CHM LC effects: range = [%.2f, %.2f]",
                       min(lc_effects_chm), max(lc_effects_chm)))
  
  lc_table_chm <- as.data.frame(ranef_chm$lc_l1_code)
  lc_table_chm$lc_l1_code <- rownames(lc_table_chm)
  write_csv(lc_table_chm, file.path(out_tables, "random_effects_lc_chm.csv"))
}

if (!is.null(ranef_dtm) && "lc_l1_code" %in% names(ranef_dtm)) {
  lc_effects_dtm <- ranef_dtm$lc_l1_code[, "Estimate", "Intercept"]
  log_progress(sprintf("  DTM LC effects: range = [%.2f, %.2f]",
                       min(lc_effects_dtm), max(lc_effects_dtm)))
  
  lc_table_dtm <- as.data.frame(ranef_dtm$lc_l1_code)
  lc_table_dtm$lc_l1_code <- rownames(lc_table_dtm)
  write_csv(lc_table_dtm, file.path(out_tables, "random_effects_lc_dtm.csv"))
}

# Forest-only LC effects
log_subsection("Random effects for forest types only")

if (!is.null(ranef_chm) && "lc_l1_code" %in% names(ranef_chm)) {
  lc_table_chm <- as.data.frame(ranef_chm$lc_l1_code)
  lc_table_chm$lc_l1_code <- rownames(lc_table_chm)
  
  lc_table_forest_chm <- lc_table_chm %>% filter(lc_l1_code %in% FOREST_CLASSES)
  
  if (nrow(lc_table_forest_chm) > 0) {
    log_progress("  CHM forest type random effects:")
    print(lc_table_forest_chm %>% select(lc_l1_code, Estimate.Intercept, Est.Error.Intercept))
    write_csv(lc_table_forest_chm, file.path(out_tables, "random_effects_lc_chm_FOREST.csv"))
  }
}

if (!is.null(ranef_dtm) && "lc_l1_code" %in% names(ranef_dtm)) {
  lc_table_dtm <- as.data.frame(ranef_dtm$lc_l1_code)
  lc_table_dtm$lc_l1_code <- rownames(lc_table_dtm)
  
  lc_table_forest_dtm <- lc_table_dtm %>% filter(lc_l1_code %in% FOREST_CLASSES)
  
  if (nrow(lc_table_forest_dtm) > 0) {
    log_progress("  DTM forest type random effects:")
    print(lc_table_forest_dtm %>% select(lc_l1_code, Estimate.Intercept, Est.Error.Intercept))
    write_csv(lc_table_forest_dtm, file.path(out_tables, "random_effects_lc_dtm_FOREST.csv"))
  }
}

# Site effects
if (!is.null(ranef_chm) && "site" %in% names(ranef_chm)) {
  site_effects_chm <- ranef_chm$site[, "Estimate", "Intercept"]
  log_progress(sprintf("  CHM site effects: range = [%.2f, %.2f]",
                       min(site_effects_chm), max(site_effects_chm)))
  
  site_table_chm <- as.data.frame(ranef_chm$site)
  site_table_chm$site <- rownames(site_table_chm)
  write_csv(site_table_chm, file.path(out_tables, "random_effects_site_chm.csv"))
}

if (!is.null(ranef_dtm) && "site" %in% names(ranef_dtm)) {
  site_effects_dtm <- ranef_dtm$site[, "Estimate", "Intercept"]
  log_progress(sprintf("  DTM site effects: range = [%.2f, %.2f]",
                       min(site_effects_dtm), max(site_effects_dtm)))
  
  site_table_dtm <- as.data.frame(ranef_dtm$site)
  site_table_dtm$site <- rownames(site_table_dtm)
  write_csv(site_table_dtm, file.path(out_tables, "random_effects_site_dtm.csv"))
}

# =====================================================================
# 4. NEW: PREDICTIVE ACCURACY METRICS
# =====================================================================

log_subsection("Predictive accuracy assessment (NEW)")

log_progress("  Computing posterior predictions...")

# Safe posterior prediction function
compute_predictions_safe <- function(model, name) {
  tryCatch({
    y_rep <- posterior_predict(model, ndraws = 100)
    y_obs <- model$data[[ifelse(grepl("chm", name, ignore.case = TRUE), 
                                 "chm_error_mean", "dtm_error_mean")]]
    y_pred <- colMeans(y_rep)
    list(y_rep = y_rep, y_obs = y_obs, y_pred = y_pred, success = TRUE)
  }, error = function(e) {
    log_progress(sprintf("    ⚠ Posterior prediction failed for %s: %s", name, e$message))
    list(success = FALSE)
  })
}

pred_chm <- compute_predictions_safe(model_chm_final, "CHM")
pred_dtm <- compute_predictions_safe(model_dtm_final, "DTM")

# Calculate metrics only if predictions succeeded
pred_metrics_list <- list()

if (pred_chm$success) {
  pred_metrics_list$chm <- tibble::tibble(
    product = "CHM",
    n_obs = length(pred_chm$y_obs),
    RMSE = sqrt(mean((pred_chm$y_obs - pred_chm$y_pred)^2)),
    MAE = mean(abs(pred_chm$y_obs - pred_chm$y_pred)),
    bias = mean(pred_chm$y_pred - pred_chm$y_obs),
    correlation = cor(pred_chm$y_obs, pred_chm$y_pred),
    R2 = cor(pred_chm$y_obs, pred_chm$y_pred)^2
  )
  
  # Observed vs Predicted plot
  p_obs_pred_chm <- ggplot(data.frame(observed = pred_chm$y_obs, predicted = pred_chm$y_pred),
                           aes(x = observed, y = predicted)) +
    geom_hex(bins = 50) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    scale_fill_viridis_c(trans = "log10") +
    coord_fixed() +
    labs(title = "CHM: Observed vs Predicted",
         subtitle = sprintf("R² = %.3f, RMSE = %.2f m", 
                            pred_metrics_list$chm$R2, pred_metrics_list$chm$RMSE),
         x = "Observed error (m)", y = "Predicted error (m)") +
    theme_cowplot()
  
  ggsave(file.path(out_plots, "06e_obs_vs_pred_chm.pdf"), p_obs_pred_chm,
         width = 7, height = 6, bg = "white")
} else {
  log_progress("  CHM predictions: SKIPPED (computation failed)")
}

if (pred_dtm$success) {
  pred_metrics_list$dtm <- tibble::tibble(
    product = "DTM",
    n_obs = length(pred_dtm$y_obs),
    RMSE = sqrt(mean((pred_dtm$y_obs - pred_dtm$y_pred)^2)),
    MAE = mean(abs(pred_dtm$y_obs - pred_dtm$y_pred)),
    bias = mean(pred_dtm$y_pred - pred_dtm$y_obs),
    correlation = cor(pred_dtm$y_obs, pred_dtm$y_pred),
    R2 = cor(pred_dtm$y_obs, pred_dtm$y_pred)^2
  )
} else {
  log_progress("  DTM predictions: SKIPPED (computation failed)")
}

if (length(pred_metrics_list) > 0) {
  pred_metrics <- bind_rows(pred_metrics_list)
  log_progress("  Predictive accuracy:")
  print(pred_metrics)
  write_csv(pred_metrics, file.path(out_tables, "predictive_accuracy.csv"))
} else {
  pred_metrics <- NULL
  log_progress("  ⚠ No predictive metrics computed")
}

# =====================================================================
# 4b. NEW: FOREST-FOCUSED PREDICTIVE ACCURACY
# =====================================================================

log_subsection("Predictive accuracy assessment (FOREST ONLY)")

if (pred_chm$success) {
  # Get model data and add predictions
  model_data_pred <- model_chm_final$data
  model_data_pred$y_obs <- pred_chm$y_obs
  model_data_pred$y_pred <- pred_chm$y_pred
  
  # Filter to forest classes
  forest_pred_data <- model_data_pred %>% filter(lc_l1_code %in% FOREST_CLASSES)
  
  if (nrow(forest_pred_data) > 100) {
    # Overall forest metrics
    pred_metrics_forest <- tibble::tibble(
      product = "CHM_FOREST",
      n_obs = nrow(forest_pred_data),
      RMSE = sqrt(mean((forest_pred_data$y_obs - forest_pred_data$y_pred)^2)),
      MAE = mean(abs(forest_pred_data$y_obs - forest_pred_data$y_pred)),
      bias = mean(forest_pred_data$y_pred - forest_pred_data$y_obs),
      correlation = cor(forest_pred_data$y_obs, forest_pred_data$y_pred),
      R2 = cor(forest_pred_data$y_obs, forest_pred_data$y_pred)^2
    )
    
    # Per-forest-type metrics
    pred_metrics_by_forest <- forest_pred_data %>%
      group_by(lc_l1_code) %>%
      summarise(
        n_obs = n(),
        RMSE = sqrt(mean((y_obs - y_pred)^2)),
        MAE = mean(abs(y_obs - y_pred)),
        bias = mean(y_pred - y_obs),
        correlation = cor(y_obs, y_pred),
        R2 = cor(y_obs, y_pred)^2,
        .groups = "drop"
      )
    
    log_progress("  Forest-only predictive accuracy:")
    print(pred_metrics_forest)
    log_progress("  By forest type:")
    print(pred_metrics_by_forest)
    
    write_csv(pred_metrics_forest, file.path(out_tables, "predictive_accuracy_FOREST.csv"))
    write_csv(pred_metrics_by_forest, file.path(out_tables, "predictive_accuracy_by_lc_FOREST.csv"))
    
    # Forest-only observed vs predicted plot
    tryCatch({
      p_obs_pred_forest <- ggplot(forest_pred_data, aes(x = y_obs, y = y_pred)) +
        geom_hex(bins = 40) +
        geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
        scale_fill_viridis_c(trans = "log10") +
        coord_fixed() +
        labs(title = "CHM: Observed vs Predicted — Forest Types Only",
             subtitle = sprintf("R² = %.3f, RMSE = %.2f m", 
                                pred_metrics_forest$R2, pred_metrics_forest$RMSE),
             x = "Observed error (m)", y = "Predicted error (m)") +
        theme_cowplot()
      
      ggsave(file.path(out_plots, "06e_obs_vs_pred_chm_FOREST.pdf"), p_obs_pred_forest,
             width = 7, height = 6, bg = "white")
    }, error = function(e) {
      log_progress(sprintf("    ⚠ Forest obs vs pred plot failed: %s", e$message))
    })
    
    # Faceted by forest type
    tryCatch({
      p_obs_pred_forest_facet <- ggplot(forest_pred_data, aes(x = y_obs, y = y_pred)) +
        geom_hex(bins = 30) +
        geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
        facet_wrap(~lc_l1_code, ncol = 2) +
        scale_fill_viridis_c(trans = "log10") +
        coord_fixed() +
        labs(title = "CHM: Observed vs Predicted by Forest Type",
             x = "Observed error (m)", y = "Predicted error (m)") +
        theme_cowplot()
      
      ggsave(file.path(out_plots, "06e_obs_vs_pred_chm_by_forest.pdf"), p_obs_pred_forest_facet,
             width = 10, height = 8, bg = "white")
    }, error = function(e) {
      log_progress(sprintf("    ⚠ Forest faceted obs vs pred plot failed: %s", e$message))
    })
    
  } else {
    log_progress("  ⚠ Insufficient forest observations for forest-specific metrics")
  }
} else {
  log_progress("  ⚠ Forest predictions: SKIPPED (CHM predictions failed)")
}

# =====================================================================
# 5. NEW: STRATIFIED RESIDUAL ANALYSIS
# =====================================================================

log_subsection("Stratified residual analysis (NEW)")

# Safe residual extraction
extract_residuals_safe <- function(model, name) {
  tryCatch({
    resid <- residuals(model, type = "response")[, "Estimate"]
    fitted_vals <- fitted(model)[, "Estimate"]
    list(residual = resid, fitted = fitted_vals, success = TRUE)
  }, error = function(e) {
    log_progress(sprintf("    ⚠ Residual extraction failed for %s: %s", name, e$message))
    list(success = FALSE)
  })
}

resid_result <- extract_residuals_safe(model_chm_final, "CHM")

if (resid_result$success) {
  model_data_chm <- model_chm_final$data
  model_data_chm$residual <- resid_result$residual
  model_data_chm$fitted <- resid_result$fitted
  
  # Residuals by slope bin
  log_progress("  Analyzing residuals by slope...")
  model_data_chm <- model_data_chm %>%
    mutate(slope_bin = cut(slope_mean_z, 
                           breaks = quantile(slope_mean_z, seq(0, 1, 0.2), na.rm = TRUE),
                           include.lowest = TRUE,
                           labels = c("Q1", "Q2", "Q3", "Q4", "Q5")))
  
  resid_by_slope <- model_data_chm %>%
    group_by(slope_bin) %>%
    summarise(
      n = n(),
      mean_resid = mean(residual, na.rm = TRUE),
      sd_resid = sd(residual, na.rm = TRUE),
      median_resid = median(residual, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(!is.na(slope_bin))
  
  log_progress("    Residuals by slope quintile:")
  print(resid_by_slope)
  
  # Check for systematic patterns
  max_bias_slope <- max(abs(resid_by_slope$mean_resid))
  if (max_bias_slope > 0.5) {
    log_progress(sprintf("    ⚠ Systematic residual bias by slope (max = %.2f m)", max_bias_slope))
  }
  
  # Residuals by LC
  log_progress("  Analyzing residuals by land cover...")
  resid_by_lc <- model_data_chm %>%
    group_by(lc_l1_code) %>%
    summarise(
      n = n(),
      mean_resid = mean(residual, na.rm = TRUE),
      sd_resid = sd(residual, na.rm = TRUE),
      median_resid = median(residual, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(n > 100)
  
  log_progress("    Residuals by LC (classes with n>100):")
  print(resid_by_lc)
  
  # Residuals by site
  log_progress("  Analyzing residuals by site...")
  resid_by_site <- model_data_chm %>%
    group_by(site) %>%
    summarise(
      n = n(),
      mean_resid = mean(residual, na.rm = TRUE),
      sd_resid = sd(residual, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Flag sites with high residual bias
  biased_sites <- resid_by_site %>% filter(abs(mean_resid) > 1)
  if (nrow(biased_sites) > 0) {
    log_progress(sprintf("    ⚠ %d sites with mean |residual| > 1m", nrow(biased_sites)))
  }
  
  # Export residual summaries
  write_csv(resid_by_slope, file.path(out_tables, "residuals_by_slope.csv"))
  write_csv(resid_by_lc, file.path(out_tables, "residuals_by_lc.csv"))
  write_csv(resid_by_site, file.path(out_tables, "residuals_by_site.csv"))
  
  # Plot residuals vs fitted by stratum
  tryCatch({
    p_resid_slope <- ggplot(model_data_chm %>% filter(!is.na(slope_bin)), 
                            aes(x = fitted, y = residual)) +
      geom_hex(bins = 30) +
      geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
      facet_wrap(~slope_bin, ncol = 3) +
      scale_fill_viridis_c(trans = "log10") +
      labs(title = "CHM Residuals vs Fitted by Slope Quintile",
           x = "Fitted value (m)", y = "Residual (m)") +
      theme_cowplot() +
      theme(strip.background = element_rect(fill = "lightblue"))
    
    ggsave(file.path(out_plots, "08_residuals_by_slope.pdf"), p_resid_slope,
           width = 12, height = 8, bg = "white")
    log_progress("    ✓ Residual by slope plot saved")
  }, error = function(e) {
    log_progress(sprintf("    ⚠ Residual by slope plot failed: %s", e$message))
  })
  
  # Top 6 LC classes
  tryCatch({
    top_lc <- resid_by_lc %>% arrange(desc(n)) %>% slice(1:6) %>% pull(lc_l1_code)
    
    p_resid_lc <- ggplot(model_data_chm %>% filter(lc_l1_code %in% top_lc),
                         aes(x = fitted, y = residual)) +
      geom_hex(bins = 30) +
      geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
      facet_wrap(~lc_l1_code, ncol = 3) +
      scale_fill_viridis_c(trans = "log10") +
      labs(title = "CHM Residuals vs Fitted by Land Cover",
           x = "Fitted value (m)", y = "Residual (m)") +
      theme_cowplot()
    
    ggsave(file.path(out_plots, "08_residuals_by_lc.pdf"), p_resid_lc,
           width = 12, height = 10, bg = "white")
    log_progress("    ✓ Residual by LC plot saved")
  }, error = function(e) {
    log_progress(sprintf("    ⚠ Residual by LC plot failed: %s", e$message))
  })
  
  # =====================================================================
  # 5b. NEW: FOREST-FOCUSED RESIDUAL ANALYSIS
  # =====================================================================
  
  log_subsection("Stratified residual analysis (FOREST ONLY)")
  
  # Filter to forest classes
  model_data_forest <- model_data_chm %>% 
    filter(lc_l1_code %in% FOREST_CLASSES)
  
  if (nrow(model_data_forest) > 100) {
    # Residuals by forest type
    resid_by_forest <- model_data_forest %>%
      group_by(lc_l1_code) %>%
      summarise(
        n = n(),
        mean_resid = mean(residual, na.rm = TRUE),
        sd_resid = sd(residual, na.rm = TRUE),
        median_resid = median(residual, na.rm = TRUE),
        rmse_resid = sqrt(mean(residual^2, na.rm = TRUE)),
        .groups = "drop"
      )
    
    log_progress("    Residuals by forest type:")
    print(resid_by_forest)
    write_csv(resid_by_forest, file.path(out_tables, "residuals_by_lc_FOREST.csv"))
    
    # Residuals by slope for forests
    resid_forest_slope <- model_data_forest %>%
      filter(!is.na(slope_bin)) %>%
      group_by(slope_bin) %>%
      summarise(
        n = n(),
        mean_resid = mean(residual, na.rm = TRUE),
        sd_resid = sd(residual, na.rm = TRUE),
        .groups = "drop"
      )
    write_csv(resid_forest_slope, file.path(out_tables, "residuals_by_slope_FOREST.csv"))
    
    # Plot residuals by forest type
    tryCatch({
      p_resid_forest <- ggplot(model_data_forest,
                               aes(x = fitted, y = residual)) +
        geom_hex(bins = 30) +
        geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
        facet_wrap(~lc_l1_code, ncol = 2) +
        scale_fill_viridis_c(trans = "log10") +
        labs(title = "CHM Residuals vs Fitted — Forest Types Only",
             x = "Fitted value (m)", y = "Residual (m)") +
        theme_cowplot()
      
      ggsave(file.path(out_plots, "08_residuals_by_lc_FOREST.pdf"), p_resid_forest,
             width = 10, height = 8, bg = "white")
      log_progress("    ✓ Forest residual plot saved")
    }, error = function(e) {
      log_progress(sprintf("    ⚠ Forest residual plot failed: %s", e$message))
    })
    
    # Residual distribution by forest type (density plot)
    tryCatch({
      p_resid_forest_dens <- ggplot(model_data_forest,
                                     aes(x = residual, fill = lc_l1_code)) +
        geom_density(alpha = 0.5) +
        geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
        scale_fill_brewer(palette = "Dark2") +
        coord_cartesian(xlim = c(-15, 15)) +
        labs(title = "CHM Residual Distribution by Forest Type",
             x = "Residual (m)", y = "Density",
             fill = "Forest Type") +
        theme_cowplot()
      
      ggsave(file.path(out_plots, "08_residual_density_FOREST.pdf"), p_resid_forest_dens,
             width = 10, height = 6, bg = "white")
    }, error = function(e) {
      log_progress(sprintf("    ⚠ Forest residual density plot failed: %s", e$message))
    })
    
    # Slope × forest type residual interaction
    tryCatch({
      resid_forest_slope_lc <- model_data_forest %>%
        filter(!is.na(slope_bin)) %>%
        group_by(lc_l1_code, slope_bin) %>%
        summarise(
          n = n(),
          mean_resid = mean(residual, na.rm = TRUE),
          se_resid = sd(residual, na.rm = TRUE) / sqrt(n()),
          .groups = "drop"
        ) %>%
        filter(n >= 30)
      
      p_resid_forest_interact <- ggplot(resid_forest_slope_lc,
                                         aes(x = slope_bin, y = mean_resid, 
                                             color = lc_l1_code, group = lc_l1_code)) +
        geom_line(linewidth = 1) +
        geom_point(size = 2) +
        geom_errorbar(aes(ymin = mean_resid - 1.96*se_resid,
                          ymax = mean_resid + 1.96*se_resid), width = 0.2) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
        scale_color_brewer(palette = "Dark2") +
        labs(title = "Mean Residual by Slope × Forest Type",
             subtitle = "Non-zero patterns suggest missing interactions",
             x = "Slope quintile", y = "Mean residual (m)",
             color = "Forest Type") +
        theme_cowplot() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      ggsave(file.path(out_plots, "08_residuals_slope_forest_interaction.pdf"), 
             p_resid_forest_interact, width = 10, height = 6, bg = "white")
      
      write_csv(resid_forest_slope_lc, file.path(out_tables, "residuals_by_slope_by_lc_FOREST.csv"))
    }, error = function(e) {
      log_progress(sprintf("    ⚠ Forest slope interaction plot failed: %s", e$message))
    })
    
  } else {
    log_progress("    ⚠ Insufficient forest observations for residual analysis")
  }
  
} else {
  log_progress("  Stratified residual analysis: SKIPPED (residual extraction failed)")
  resid_by_slope <- NULL
  resid_by_lc <- NULL
  resid_by_site <- NULL
  model_data_chm <- NULL
}

# =====================================================================
# 6. NEW: EFFECT SIZE INTERPRETATION
# =====================================================================

log_subsection("Effect size interpretation (NEW)")

# Safe effect size extraction
effect_table <- tryCatch({
  # Get fixed effects
  fe_chm <- fixef(model_chm_final)
  
  # Get posterior draws for probability of direction
  draws <- as_draws_df(model_chm_final)
  
  # Build effect interpretation table
  effect_tbl <- tibble::tibble(
    predictor = rownames(fe_chm),
    estimate = fe_chm[, "Estimate"],
    std_error = fe_chm[, "Est.Error"],
    lower_95 = fe_chm[, "Q2.5"],
    upper_95 = fe_chm[, "Q97.5"],
    effect_1sd_m = fe_chm[, "Estimate"],
    prob_positive = NA_real_,
    prob_negative = NA_real_
  )
  
  # Calculate probability of direction
  for (i in 1:nrow(effect_tbl)) {
    param_name <- paste0("b_", effect_tbl$predictor[i])
    param_name <- gsub(":", "", param_name)
    
    matching_cols <- grep(effect_tbl$predictor[i], names(draws), value = TRUE, fixed = TRUE)
    if (length(matching_cols) > 0) {
      param_draws <- draws[[matching_cols[1]]]
      effect_tbl$prob_positive[i] <- mean(param_draws > 0)
      effect_tbl$prob_negative[i] <- mean(param_draws < 0)
    }
  }
  
  # Add interpretation columns
  effect_tbl <- effect_tbl %>%
    mutate(
      magnitude = case_when(
        abs(estimate) < 0.25 ~ "Negligible",
        abs(estimate) < 0.5 ~ "Small",
        abs(estimate) < 1.0 ~ "Moderate",
        abs(estimate) < 2.0 ~ "Large",
        TRUE ~ "Very large"
      ),
      excludes_zero = (lower_95 > 0) | (upper_95 < 0),
      direction = case_when(
        prob_positive > 0.975 ~ "Positive (very certain)",
        prob_positive > 0.95 ~ "Positive (certain)",
        prob_positive > 0.90 ~ "Positive (likely)",
        prob_negative > 0.975 ~ "Negative (very certain)",
        prob_negative > 0.95 ~ "Negative (certain)",
        prob_negative > 0.90 ~ "Negative (likely)",
        TRUE ~ "Uncertain"
      )
    )
  
  effect_tbl
}, error = function(e) {
  log_progress(sprintf("  ⚠ Effect size extraction failed: %s", e$message))
  NULL
})

if (!is.null(effect_table)) {
  # Filter to main effects (not interactions or smooth terms)
  main_effects <- effect_table %>%
    filter(!grepl(":", predictor),
           !grepl("Intercept", predictor),
           !grepl("^s\\(", predictor))
  
  log_progress("  Main effect interpretation:")
  print(main_effects %>% select(predictor, estimate, std_error, magnitude, direction))
  
  write_csv(effect_table, file.path(out_tables, "effect_size_interpretation.csv"))
  
  # Create effect size plot
  tryCatch({
    p_effects <- ggplot(main_effects, aes(x = reorder(predictor, estimate), y = estimate)) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = lower_95, ymax = upper_95), width = 0.2) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      coord_flip() +
      labs(title = "CHM Model: Fixed Effect Estimates",
           subtitle = "Error bars show 95% credible intervals",
           x = "Predictor", y = "Effect on error (m per 1 SD)") +
      theme_cowplot()
    
    ggsave(file.path(out_plots, "09_effect_sizes.pdf"), p_effects,
           width = 10, height = 8, bg = "white")
    log_progress("  ✓ Effect size plot saved")
  }, error = function(e) {
    log_progress(sprintf("  ⚠ Effect size plot failed: %s", e$message))
  })
} else {
  log_progress("  Effect size interpretation: SKIPPED")
  main_effects <- NULL
}

# =====================================================================
# 7. LOO AND PARETO-k (Original + Enhanced)
# =====================================================================

log_subsection("LOO-CV and Pareto-k diagnostics")

pk_chm <- NULL
pk_dtm <- NULL

if (!is.null(loo_chm_final)) {
  tryCatch({
    pk_chm <- as.numeric(loo_chm_final$diagnostics$pareto_k)
    n_high_k_chm <- sum(pk_chm > 0.7, na.rm = TRUE)
    pct_high_k_chm <- 100 * n_high_k_chm / length(pk_chm)
    
    log_progress(sprintf("  CHM: %d (%.2f%%) observations with Pareto k > 0.7",
                         n_high_k_chm, pct_high_k_chm))
    
    # Identify high-k observations
    if (n_high_k_chm > 0 && !is.null(model_data_chm)) {
      high_k_indices <- which(pk_chm > 0.7)
      if (length(high_k_indices) <= nrow(model_data_chm)) {
        high_k_data <- model_data_chm[high_k_indices, ]
        log_progress("    High Pareto-k observations tend to have:")
        log_progress(sprintf("      Higher |residual|: %.2f vs %.2f overall",
                             mean(abs(high_k_data$residual), na.rm = TRUE),
                             mean(abs(model_data_chm$residual), na.rm = TRUE)))
      }
    }
  }, error = function(e) {
    log_progress(sprintf("  ⚠ CHM Pareto-k analysis failed: %s", e$message))
  })
} else {
  log_progress("  CHM LOO: Not available")
}

if (!is.null(loo_dtm_final)) {
  tryCatch({
    pk_dtm <- as.numeric(loo_dtm_final$diagnostics$pareto_k)
    n_high_k_dtm <- sum(pk_dtm > 0.7, na.rm = TRUE)
    pct_high_k_dtm <- 100 * n_high_k_dtm / length(pk_dtm)
    
    log_progress(sprintf("  DTM: %d (%.2f%%) observations with Pareto k > 0.7",
                         n_high_k_dtm, pct_high_k_dtm))
  }, error = function(e) {
    log_progress(sprintf("  ⚠ DTM Pareto-k analysis failed: %s", e$message))
  })
} else {
  log_progress("  DTM LOO: Not available")
}

# Pareto k histogram (only if we have data)
if (!is.null(pk_chm) || !is.null(pk_dtm)) {
  tryCatch({
    pdf(file.path(out_plots, "06c_pareto_k_hist.pdf"), width = 10, height = 4.2)
    par(mfrow = c(1, 2))
    if (!is.null(pk_chm)) {
      hist(pk_chm, breaks = 30, main = "Pareto-k (CHM)", xlab = "k", col = "steelblue")
      abline(v = c(0.5, 0.7), col = c("orange", "red"), lty = 2, lwd = 2)
    } else {
      plot.new()
      text(0.5, 0.5, "CHM LOO not available", cex = 1.5)
    }
    if (!is.null(pk_dtm)) {
      hist(pk_dtm, breaks = 30, main = "Pareto-k (DTM)", xlab = "k", col = "darkorange")
      abline(v = c(0.5, 0.7), col = c("orange", "red"), lty = 2, lwd = 2)
    } else {
      plot.new()
      text(0.5, 0.5, "DTM LOO not available", cex = 1.5)
    }
    dev.off()
    log_progress("  ✓ Pareto-k histogram saved")
  }, error = function(e) {
    log_progress(sprintf("  ⚠ Pareto-k histogram failed: %s", e$message))
    try(dev.off(), silent = TRUE)
  })
} else {
  log_progress("  Pareto-k histograms: SKIPPED (no LOO results available)")
}

# =====================================================================
# 8. BAYES R²
# =====================================================================

log_subsection("Bayes R² calculation")

br2_chm <- tryCatch({
  bayes_R2(model_chm_final)
}, error = function(e) {
  log_progress(sprintf("  ⚠ CHM Bayes R² failed: %s", e$message))
  NULL
})

br2_dtm <- tryCatch({
  bayes_R2(model_dtm_final)
}, error = function(e) {
  log_progress(sprintf("  ⚠ DTM Bayes R² failed: %s", e$message))
  NULL
})

if (!is.null(br2_chm)) {
  log_progress(sprintf("  CHM Bayes R²: %.3f [%.3f, %.3f]", 
                       median(br2_chm[,1]), 
                       quantile(br2_chm[,1], 0.025),
                       quantile(br2_chm[,1], 0.975)))
}

if (!is.null(br2_dtm)) {
  log_progress(sprintf("  DTM Bayes R²: %.3f [%.3f, %.3f]", 
                       median(br2_dtm[,1]), 
                       quantile(br2_dtm[,1], 0.025),
                       quantile(br2_dtm[,1], 0.975)))
}

# =====================================================================
# 9. CONDITIONAL EFFECTS (Original)
# =====================================================================

log_subsection("Conditional effects plots")

save_ce_plots <- function(fit, effects, dpar, outfile) {
  model_vars <- names(fit$data)
  effects_present <- character(0)
  
  for (eff in effects) {
    if (grepl(":", eff, fixed = TRUE)) {
      parts <- strsplit(eff, ":")[[1]]
      if (all(parts %in% model_vars)) {
        effects_present <- c(effects_present, eff)
      }
    } else {
      if (eff %in% model_vars) {
        effects_present <- c(effects_present, eff)
      }
    }
  }
  
  if (length(effects_present) == 0) {
    log_progress(sprintf("  ⚠ Skipping %s: no effects available", basename(outfile)))
    return(invisible(NULL))
  }
  
  tryCatch({
    ce <- conditional_effects(fit, effects = effects_present, dpar = dpar)
    plots <- plot(ce, ask = FALSE)
    g <- cowplot::plot_grid(plotlist = plots, ncol = 3)
    ggsave(outfile, g, width = 12, height = ceiling(length(plots)/3)*4, bg = "white")
    log_progress(sprintf("  ✓ Saved %s", basename(outfile)))
  }, error = function(e) {
    log_progress(sprintf("  ⚠ Error creating %s: %s", basename(outfile), e$message))
  })
}

base_effects <- c("slope_mean_z", "wsci_z", "rh_98_z", "cover_z", "lc_l1_code")
meta_effects_chm <- intersect(c("meta_offnad_z", "meta_sunel_z", "meta_stereo_z", "meta_leafon_z"), 
                              available_meta_chm)

chm_mu_effects <- c(base_effects, meta_effects_chm)
save_ce_plots(model_chm_final, chm_mu_effects, "mu", file.path(out_plots, "07a_ce_chm_mu.pdf"))
save_ce_plots(model_chm_final, c("slope_mean_z:lc_l1_code"), "mu", file.path(out_plots, "07b_ce_chm_mu_inter.pdf"))

if (is_student_chm) {
  save_ce_plots(model_chm_final, c("slope_mean_z", "wsci_z"), "sigma",
                file.path(out_plots, "07c_ce_chm_sigma.pdf"))
}

# =====================================================================
# 10. SUMMARY STATISTICS TABLE
# =====================================================================

log_subsection("Creating summary statistics table")

# Build summary stats with available data
summary_stats <- tryCatch({
  stats_df <- tibble::tibble(
    Model = c("CHM", "DTM"),
    Family = c(family(model_chm_final)$family, family(model_dtm_final)$family),
    N_obs = c(nobs(model_chm_final), nobs(model_dtm_final)),
    Max_Rhat = c(max(rhat_chm, na.rm = TRUE), max(rhat_dtm, na.rm = TRUE)),
    Min_ESS_ratio = c(min(neff_ratio(model_chm_final), na.rm = TRUE),
                      min(neff_ratio(model_dtm_final), na.rm = TRUE))
  )
  
  # Add predictive metrics if available
  if (!is.null(pred_metrics) && nrow(pred_metrics) >= 2) {
    stats_df$RMSE <- pred_metrics$RMSE
    stats_df$MAE <- pred_metrics$MAE
  } else {
    stats_df$RMSE <- NA_real_
    stats_df$MAE <- NA_real_
  }
  
  # Add Bayes R² if available
  if (!is.null(br2_chm)) {
    stats_df$Bayes_R2_median[1] <- median(br2_chm[,1])
    stats_df$Bayes_R2_lower[1] <- quantile(br2_chm[,1], 0.025)
    stats_df$Bayes_R2_upper[1] <- quantile(br2_chm[,1], 0.975)
  } else {
    stats_df$Bayes_R2_median[1] <- NA_real_
    stats_df$Bayes_R2_lower[1] <- NA_real_
    stats_df$Bayes_R2_upper[1] <- NA_real_
  }
  
  if (!is.null(br2_dtm)) {
    stats_df$Bayes_R2_median[2] <- median(br2_dtm[,1])
    stats_df$Bayes_R2_lower[2] <- quantile(br2_dtm[,1], 0.025)
    stats_df$Bayes_R2_upper[2] <- quantile(br2_dtm[,1], 0.975)
  } else {
    stats_df$Bayes_R2_median[2] <- NA_real_
    stats_df$Bayes_R2_lower[2] <- NA_real_
    stats_df$Bayes_R2_upper[2] <- NA_real_
  }
  
  # Add ν if Student-t
  if (is_student_chm && !is.null(nu_summary_chm)) {
    stats_df$nu_median[1] <- nu_summary_chm["median"]
    stats_df$nu_lower[1] <- nu_summary_chm["q025"]
    stats_df$nu_upper[1] <- nu_summary_chm["q975"]
  }
  
  if (is_student_dtm && !is.null(nu_summary_dtm)) {
    stats_df$nu_median[2] <- nu_summary_dtm["median"]
    stats_df$nu_lower[2] <- nu_summary_dtm["q025"]
    stats_df$nu_upper[2] <- nu_summary_dtm["q975"]
  }
  
  stats_df
}, error = function(e) {
  log_progress(sprintf("  ⚠ Summary stats creation failed: %s", e$message))
  NULL
})

if (!is.null(summary_stats)) {
  write_csv(summary_stats, file.path(out_tables, "model_summary_statistics.csv"))
  log_progress("  ✓ Summary statistics table saved")
} else {
  log_progress("  ⚠ Summary statistics table: SKIPPED")
}

log_progress("✓ Enhanced model diagnostics complete")

# =====================================================================
# SAVE CHECKPOINT
# =====================================================================

if (!exists("BATCH_MODE") || !BATCH_MODE) {
  log_progress("Saving checkpoint for interactive mode...")
  
  checkpoint_data <- list(
    diagnostics_complete = TRUE
  )
  
  # Only add items that exist and are not NULL
  if (!is.null(br2_chm)) checkpoint_data$br2_chm <- br2_chm
  if (!is.null(br2_dtm)) checkpoint_data$br2_dtm <- br2_dtm
  if (!is.null(loo_chm_final)) checkpoint_data$loo_chm_final <- loo_chm_final
  if (!is.null(loo_dtm_final)) checkpoint_data$loo_dtm_final <- loo_dtm_final
  if (!is.null(summary_stats)) checkpoint_data$summary_stats <- summary_stats
  if (!is.null(pred_metrics)) checkpoint_data$pred_metrics <- pred_metrics
  if (!is.null(effect_table)) checkpoint_data$effect_table <- effect_table
  if (exists("resid_by_slope") && !is.null(resid_by_slope)) checkpoint_data$resid_by_slope <- resid_by_slope
  if (exists("resid_by_lc") && !is.null(resid_by_lc)) checkpoint_data$resid_by_lc <- resid_by_lc
  if (exists("resid_by_site") && !is.null(resid_by_site)) checkpoint_data$resid_by_site <- resid_by_site
  if (!is.null(nu_summary_chm)) checkpoint_data$nu_summary_chm <- nu_summary_chm
  if (!is.null(nu_summary_dtm)) checkpoint_data$nu_summary_dtm <- nu_summary_dtm
  
  save_checkpoint("11_diagnostics", checkpoint_data)
}
