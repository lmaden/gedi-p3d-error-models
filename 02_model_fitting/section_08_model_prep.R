# =====================================================================
# section_08_model_prep.R (ENHANCED VERSION)
# Prepare modeling datasets (Stage 1 & 2)
#
# ENHANCEMENTS ADDED:
#   1. Prior predictive checks - validates priors before fitting
#   2. Sample balance verification - ensures stratification worked
#   3. Data summary statistics for modeling
#   4. Separate configurable sample fractions for Stage 1 vs Stage 2
#
# SAMPLE SIZE CONFIGURATION:
#   - Stage 1 (Section 09): 20% - for initial model comparison
#   - Stage 2 (Section 10): 33% - for final model fitting
#   - Can override via environment variables:
#       Sys.setenv(STAGE1_FRAC = "0.20", STAGE2_FRAC = "0.33")
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
    available_meta_chm <- data$available_meta_chm
    available_meta_dtm <- data$available_meta_dtm
  } else {
    stop("Data not available. Run: source('section_01_ingest.R') first")
  }
}

log_progress("Preparing modeling datasets (ENHANCED)...")

set.seed(2025)
include_ecoregion_re <- TRUE

# =====================================================================
# SAMPLE SIZE CONFIGURATION
# =====================================================================
# Stage 1 (Section 09): Initial model comparison - smaller sample for speed
# Stage 2 (Section 10): Final model fitting - larger sample for accuracy

stage1_frac <- as.numeric(Sys.getenv("STAGE1_FRAC", "0.20"))  # 20% default
stage2_frac <- as.numeric(Sys.getenv("STAGE2_FRAC", "0.33"))  # 33% default

log_progress(sprintf("  Sample configuration: Stage 1 = %.0f%%, Stage 2 = %.0f%%",
                     100 * stage1_frac, 100 * stage2_frac))

# =====================================================================
# 1. STAGE 1: SMALLER SAMPLE FOR MODEL COMPARISON (Section 09)
# =====================================================================

log_subsection("Stage 1 datasets (for section_09 model comparison)")

# Filter and sample
mod_chm_s1 <- chm_df %>%
  filter(is.finite(chm_error_mean), !is.na(lc_l1_code)) %>%
  group_by(site) %>%
  sample_frac(stage1_frac) %>%
  ungroup() %>%
  mutate(lc_l1_code = factor(lc_l1_code), 
         ecoregion = factor(ecoregion), 
         site = factor(site),
         eco_on = as.numeric(include_ecoregion_re))

mod_dtm_s1 <- dtm_df %>%
  filter(is.finite(dtm_error_mean), !is.na(lc_l1_code)) %>%
  group_by(site) %>%
  sample_frac(stage1_frac) %>%
  ungroup() %>%
  mutate(lc_l1_code = factor(lc_l1_code), 
         ecoregion = factor(ecoregion), 
         site = factor(site),
         eco_on = as.numeric(include_ecoregion_re))

log_progress(sprintf("  Stage 1 CHM: %s rows (%.0f%%)", 
                     format(nrow(mod_chm_s1), big.mark=","), 100 * stage1_frac))
log_progress(sprintf("  Stage 1 DTM: %s rows (%.0f%%)", 
                     format(nrow(mod_dtm_s1), big.mark=","), 100 * stage1_frac))

# =====================================================================
# 2. STAGE 2: LARGER SAMPLE FOR FINAL MODEL (Section 10)
# =====================================================================

log_subsection("Stage 2 datasets (for section_10 final model)")

mod_chm_s2 <- chm_df %>%
  filter(is.finite(chm_error_mean), !is.na(lc_l1_code)) %>%
  group_by(site) %>%
  sample_frac(stage2_frac) %>%
  ungroup() %>%
  mutate(lc_l1_code = factor(lc_l1_code), 
         ecoregion = factor(ecoregion),
         site = factor(site), 
         eco_on = as.numeric(include_ecoregion_re))

mod_dtm_s2 <- dtm_df %>%
  filter(is.finite(dtm_error_mean), !is.na(lc_l1_code)) %>%
  group_by(site) %>%
  sample_frac(stage2_frac) %>%
  ungroup() %>%
  mutate(lc_l1_code = factor(lc_l1_code), 
         ecoregion = factor(ecoregion),
         site = factor(site), 
         eco_on = as.numeric(include_ecoregion_re))

log_progress(sprintf("  Stage 2 CHM: %s rows (%.0f%%)", 
                     format(nrow(mod_chm_s2), big.mark=","), 100 * stage2_frac))
log_progress(sprintf("  Stage 2 DTM: %s rows (%.0f%%)", 
                     format(nrow(mod_dtm_s2), big.mark=","), 100 * stage2_frac))

# =====================================================================
# 3. SAMPLE BALANCE VERIFICATION
# =====================================================================

log_subsection("Sample balance verification (NEW)")

# Check balance by site × LC combination for Stage 1
balance_chm_s1 <- mod_chm_s1 %>%
  group_by(site, lc_l1_code) %>%
  summarise(n = n(), .groups = "drop")

balance_by_site <- balance_chm_s1 %>%
  group_by(site) %>%
  summarise(
    n_total = sum(n),
    n_lc_classes = n(),
    min_per_lc = min(n),
    max_per_lc = max(n),
    cv_per_lc = 100 * sd(n) / mean(n),
    .groups = "drop"
  )

log_progress("  Stage 1 CHM balance by site:")
log_progress(sprintf("    Total samples: %s", format(sum(balance_by_site$n_total), big.mark=",")))
log_progress(sprintf("    Sites: %d", nrow(balance_by_site)))
log_progress(sprintf("    LC classes per site: min=%d, max=%d", 
                     min(balance_by_site$n_lc_classes), max(balance_by_site$n_lc_classes)))
log_progress(sprintf("    Samples per site-LC: min=%d, max=%d",
                     min(balance_chm_s1$n), max(balance_chm_s1$n)))

# Flag problematic combinations
sparse_combos <- balance_chm_s1 %>% filter(n < 30)
if (nrow(sparse_combos) > 0) {
  log_progress(sprintf("  ⚠ %d site-LC combinations have <30 samples", nrow(sparse_combos)))
  log_progress("    Random effects may be unstable for these groups")
}

# Also check Stage 2 balance
balance_chm_s2 <- mod_chm_s2 %>%
  group_by(site, lc_l1_code) %>%
  summarise(n = n(), .groups = "drop")

balance_by_site_s2 <- balance_chm_s2 %>%
  group_by(site) %>%
  summarise(
    n_total = sum(n),
    n_lc_classes = n(),
    min_per_lc = min(n),
    max_per_lc = max(n),
    .groups = "drop"
  )

log_progress("  Stage 2 CHM balance by site:")
log_progress(sprintf("    Total samples: %s", format(sum(balance_by_site_s2$n_total), big.mark=",")))
log_progress(sprintf("    Sites: %d", nrow(balance_by_site_s2)))
log_progress(sprintf("    Samples per site-LC: min=%d, max=%d",
                     min(balance_chm_s2$n), max(balance_chm_s2$n)))

write_csv(balance_by_site, file.path(out_tables, "model_sample_balance_by_site.csv"))
write_csv(balance_chm_s1, file.path(out_tables, "model_sample_balance_sitexlc.csv"))
write_csv(balance_by_site_s2, file.path(out_tables, "model_sample_balance_by_site_s2.csv"))
write_csv(balance_chm_s2, file.path(out_tables, "model_sample_balance_sitexlc_s2.csv"))

# =====================================================================
# 4. MODELING DATA SUMMARY STATISTICS
# =====================================================================

log_subsection("Modeling data summary statistics (NEW)")

# Summary of response variable in model data
model_data_summary <- tibble::tibble(
  dataset = c("Stage 1 CHM", "Stage 1 DTM", "Stage 2 CHM", "Stage 2 DTM"),
  sample_pct = c(stage1_frac, stage1_frac, stage2_frac, stage2_frac) * 100,
  n = c(nrow(mod_chm_s1), nrow(mod_dtm_s1), nrow(mod_chm_s2), nrow(mod_dtm_s2)),
  n_sites = c(n_distinct(mod_chm_s1$site), n_distinct(mod_dtm_s1$site),
              n_distinct(mod_chm_s2$site), n_distinct(mod_dtm_s2$site)),
  n_ecoregions = c(n_distinct(mod_chm_s1$ecoregion), n_distinct(mod_dtm_s1$ecoregion),
                   n_distinct(mod_chm_s2$ecoregion), n_distinct(mod_dtm_s2$ecoregion)),
  n_lc_classes = c(n_distinct(mod_chm_s1$lc_l1_code), n_distinct(mod_dtm_s1$lc_l1_code),
                   n_distinct(mod_chm_s2$lc_l1_code), n_distinct(mod_dtm_s2$lc_l1_code)),
  error_mean = c(mean(mod_chm_s1$chm_error_mean, na.rm=T), mean(mod_dtm_s1$dtm_error_mean, na.rm=T),
                 mean(mod_chm_s2$chm_error_mean, na.rm=T), mean(mod_dtm_s2$dtm_error_mean, na.rm=T)),
  error_sd = c(sd(mod_chm_s1$chm_error_mean, na.rm=T), sd(mod_dtm_s1$dtm_error_mean, na.rm=T),
               sd(mod_chm_s2$chm_error_mean, na.rm=T), sd(mod_dtm_s2$dtm_error_mean, na.rm=T))
)

write_csv(model_data_summary, file.path(out_tables, "model_data_summary.csv"))
log_progress("  Model data summary:")
print(model_data_summary)

# =====================================================================
# 5. PRIOR PREDICTIVE CHECK
# =====================================================================
# CRITICAL: Validates priors before expensive model fitting

log_subsection("Prior predictive check (NEW)")

# Check if brms is available
if (requireNamespace("brms", quietly = TRUE)) {
  library(brms)
  
  # Define priors matching section_09
  priors_for_ppc <- c(
    prior("normal(0, 2)", class = "b"),
    prior("normal(0, 10)", class = "Intercept"),
    prior("gamma(2, 0.1)", class = "nu"),
    prior("student_t(3, 0, 2)", class = "sd")
  )
  
  # Use small subset for prior predictive check (fast)
  ppc_data <- mod_chm_s1 %>% 
    sample_n(min(2000, n())) %>%
    filter(is.finite(chm_error_mean))
  
  log_progress("  Fitting prior-only model for prior predictive check...")
  log_progress(sprintf("    Using %d samples", nrow(ppc_data)))
  
  # Fit prior-only model
  prior_only_fit <- tryCatch({
    brm(
      chm_error_mean ~ slope_mean_z + wsci_z + (1 | site),
      data = ppc_data,
      family = student(),
      prior = priors_for_ppc,
      sample_prior = "only",  # CRITICAL: only samples from prior
      chains = 2,
      iter = 1000,
      warmup = 500,
      seed = 2025,
      silent = 2,
      refresh = 0
    )
  }, error = function(e) {
    log_progress(sprintf("  ⚠ Prior predictive check failed: %s", e$message))
    NULL
  })
  
  if (!is.null(prior_only_fit)) {
    # Generate prior predictive samples
    pp_prior <- posterior_predict(prior_only_fit, ndraws = 100)
    
    # Summarize prior predictive range
    pp_range <- range(pp_prior)
    pp_mean <- mean(pp_prior)
    pp_sd <- sd(as.vector(pp_prior))
    
    log_progress(sprintf("  Prior predictive range: [%.1f, %.1f] m", pp_range[1], pp_range[2]))
    log_progress(sprintf("  Prior predictive mean: %.1f m, SD: %.1f m", pp_mean, pp_sd))
    
    # Check if prior predictions are in reasonable range
    # Errors > 100m are likely unrealistic
    reasonable_range <- c(-100, 100)
    pct_in_range <- mean(pp_prior > reasonable_range[1] & pp_prior < reasonable_range[2]) * 100
    
    log_progress(sprintf("  Predictions in reasonable range [%.0f, %.0f]m: %.1f%%",
                         reasonable_range[1], reasonable_range[2], pct_in_range))
    
    if (pct_in_range < 50) {
      log_progress("  ⚠ WARNING: Priors allow many unreasonable predictions")
      log_progress("    Consider tighter priors on coefficients")
    } else if (pct_in_range > 99) {
      log_progress("  ⚠ NOTE: Priors may be too informative")
      log_progress("    Consider slightly wider priors")
    } else {
      log_progress("  ✓ Priors produce reasonable prior predictive range")
    }
    
    # Compare to observed data range
    obs_range <- range(ppc_data$chm_error_mean)
    log_progress(sprintf("  Observed data range: [%.1f, %.1f] m", obs_range[1], obs_range[2]))
    
    # Save prior predictive check plot
    p_pp <- pp_check(prior_only_fit, ndraws = 50) + 
      ggtitle("Prior Predictive Check", 
              subtitle = "Samples from prior only - should bracket plausible data range") +
      theme_cowplot()
    
    ggsave(file.path(out_plots, "05_prior_predictive_check.pdf"), p_pp,
           width = 8, height = 5, bg = "white")
    
    log_progress("  ✓ Prior predictive check plot saved")
    
    # Export prior predictive summary
    ppc_summary <- tibble::tibble(
      metric = c("pp_min", "pp_max", "pp_mean", "pp_sd", 
                 "pct_in_reasonable_range", "obs_min", "obs_max"),
      value = c(pp_range[1], pp_range[2], pp_mean, pp_sd, 
                pct_in_range, obs_range[1], obs_range[2])
    )
    write_csv(ppc_summary, file.path(out_tables, "prior_predictive_summary.csv"))
    
  }
} else {
  log_progress("  ⚠ brms not available - skipping prior predictive check")
  prior_only_fit <- NULL
}

# =====================================================================
# 6. FACTOR LEVEL VERIFICATION
# =====================================================================

log_subsection("Factor level verification")

# Ensure factor levels are consistent between Stage 1 and Stage 2
sites_s1 <- levels(mod_chm_s1$site)
sites_s2 <- levels(factor(mod_chm_s2$site))
sites_missing_in_s2 <- setdiff(sites_s1, sites_s2)
sites_extra_in_s2 <- setdiff(sites_s2, sites_s1)

if (length(sites_missing_in_s2) > 0) {
  log_progress(sprintf("  ⚠ Sites in S1 but not S2: %s", paste(sites_missing_in_s2, collapse=", ")))
}
if (length(sites_extra_in_s2) > 0) {
  log_progress(sprintf("  ⚠ Sites in S2 but not S1: %s", paste(sites_extra_in_s2, collapse=", ")))
}

lc_s1 <- levels(mod_chm_s1$lc_l1_code)
lc_s2 <- levels(factor(mod_chm_s2$lc_l1_code))
log_progress(sprintf("  Stage 1 LC classes: %s", paste(lc_s1, collapse=", ")))
log_progress(sprintf("  Stage 2 LC classes: %s", paste(lc_s2, collapse=", ")))

log_progress("✓ Model datasets prepared (ENHANCED)")

# =====================================================================
# SAVE CHECKPOINT
# =====================================================================

if (!exists("BATCH_MODE") || !BATCH_MODE) {
  log_progress("Saving checkpoint for interactive mode...")
  
  checkpoint_data <- list(
    mod_chm_s1 = mod_chm_s1,
    mod_dtm_s1 = mod_dtm_s1,
    mod_chm_s2 = mod_chm_s2,
    mod_dtm_s2 = mod_dtm_s2,
    stage1_frac = stage1_frac,
    stage2_frac = stage2_frac,
    model_data_summary = model_data_summary,
    balance_by_site = balance_by_site,
    balance_by_site_s2 = balance_by_site_s2,
    prior_predictive_complete = !is.null(prior_only_fit)
  )
  
  if (exists("ppc_summary")) {
    checkpoint_data$ppc_summary <- ppc_summary
  }
  
  save_checkpoint("08_model_prep", checkpoint_data)
}
