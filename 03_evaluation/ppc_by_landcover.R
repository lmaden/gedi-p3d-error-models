# =====================================================================
# POSTERIOR PREDICTIVE CHECKS BY LAND COVER CLASS (UPDATED)
# =====================================================================
# Generate PPC density plots comparing observed vs predicted errors
# for each of the 4 main forested land cover classes
# 
# OUTPUT: Combined 2×4 panel figure (CHM top row, DTM bottom row)
#         suitable for publication
#
# UPDATES FROM LEGACY SCRIPT:
#   1. Uses current checkpoint paths (no subdirectory)
#   2. Uses 33% subsample fraction to match Stage 2 models
#   3. Integrates with analysis_config.R and analysis_utils.R
#   4. Creates single combined 2×4 panel figure for publication
# =====================================================================

# =====================================================================
# LOAD WORKFLOW CONFIGURATION
# =====================================================================

# Source the main config and utils (adjust path if running standalone)
config_path <- "analysis_config.R"
utils_path <- "analysis_utils.R"

if (file.exists(config_path)) {
  source(config_path)
} else {
  # Fallback configuration if running standalone
  suppressPackageStartupMessages({
    library(brms)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(bayesplot)
    library(cowplot)
  })
  
  PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/gpfs/data1/vclgp/lmaden/chpt1")
  out_plots <- file.path(PROJECT_ROOT, "plots")
  out_tables <- file.path(PROJECT_ROOT, "tables")
  dir.create(out_plots, showWarnings = FALSE, recursive = TRUE)
  dir.create(out_tables, showWarnings = FALSE, recursive = TRUE)
  theme_set(theme_cowplot())
}

if (file.exists(utils_path)) {
  source(utils_path)
} else {
  # Minimal fallback
  log_progress <- function(msg) cat(sprintf("[%s] %s\n", Sys.time(), msg))
  CHECKPOINT_DIR <- file.path(PROJECT_ROOT, "checkpoints")
}

# =====================================================================
# CONFIGURATION
# =====================================================================

# Define forested classes (excluding Mixed Forest due to low n)
forested_classes <- c("BDF", "DNF", "EBF", "ENF")

# Land cover full names for plot labels
lc_names <- c(

  "BDF" = "Broadleaf Deciduous",
  "DNF" = "Deciduous Needleleaf",
  "EBF" = "Evergreen Broadleaf",
  "ENF" = "Evergreen Needleleaf"
)

# Sample size per class for predictions (will use all available if fewer)
n_samples <- 5000

# Number of posterior draws for predictions
n_draws <- 500

# Number of posterior draws to plot (for visibility)
n_draws_plot <- 50

# Subsample fraction - MUST MATCH section_08_model_prep.R Stage 2
stage2_frac <- as.numeric(Sys.getenv("STAGE2_FRAC", "0.33"))

# Random seed for reproducibility
set.seed(2025)

log_progress("=" %>% rep(70) %>% paste(collapse = ""))
log_progress("POSTERIOR PREDICTIVE CHECKS BY LAND COVER (UPDATED)")
log_progress("=" %>% rep(70) %>% paste(collapse = ""))

# =====================================================================
# LOAD MODELS FROM CHECKPOINT
# =====================================================================

log_progress("Loading Stage 2 models from checkpoint...")

model_checkpoint <- file.path(CHECKPOINT_DIR, "10_models_stage2.rds")

if (file.exists(model_checkpoint)) {
  data_models <- readRDS(model_checkpoint)
  fit_chm_s2 <- data_models$data$fit_chm_s2
  fit_dtm_s2 <- data_models$data$fit_dtm_s2
  log_progress(sprintf("✓ Models loaded from: %s", model_checkpoint))
} else {
  stop("Model checkpoint not found at: ", model_checkpoint)
}

# =====================================================================
# LOAD AND PREPARE HOLDOUT DATA
# =====================================================================

log_progress("Loading full dataset for holdout creation...")

data_checkpoint <- file.path(CHECKPOINT_DIR, "01_data_ingest.rds")
if (!file.exists(data_checkpoint)) {
  stop("Data checkpoint not found at: ", data_checkpoint)
}

full_data <- readRDS(data_checkpoint)
chm_df_full_raw <- full_data$data$chm_df
dtm_df_full_raw <- full_data$data$dtm_df
log_progress(sprintf("✓ Loaded full dataset: %s CHM, %s DTM observations",
                     format(nrow(chm_df_full_raw), big.mark = ","),
                     format(nrow(dtm_df_full_raw), big.mark = ",")))

# Apply quality filters (same as section_08)
log_progress("Applying quality filters...")

chm_df_full <- chm_df_full_raw %>%
  filter(is.finite(chm_error_mean), !is.na(lc_l1_code), !is.na(site), !is.na(ecoregion)) %>%
  mutate(lc_l1_code = factor(lc_l1_code),
         ecoregion = factor(ecoregion),
         site = factor(site))

dtm_df_full <- dtm_df_full_raw %>%
  filter(is.finite(dtm_error_mean), !is.na(lc_l1_code), !is.na(site), !is.na(ecoregion)) %>%
  mutate(lc_l1_code = factor(lc_l1_code),
         ecoregion = factor(ecoregion),
         site = factor(site))

log_progress(sprintf("✓ Filtered dataset: %s CHM, %s DTM observations",
                     format(nrow(chm_df_full), big.mark = ","),
                     format(nrow(dtm_df_full), big.mark = ",")))

# Recreate training sample using SAME seed and fraction as section_08
log_progress(sprintf("Recreating training sample (seed=2025, frac=%.0f%%)...", 100 * stage2_frac))

set.seed(2025)

chm_df_full$original_row <- 1:nrow(chm_df_full)
dtm_df_full$original_row <- 1:nrow(dtm_df_full)

mod_chm_s2_train <- chm_df_full %>%
  group_by(site) %>%
  sample_frac(stage2_frac) %>%
  ungroup()

mod_dtm_s2_train <- dtm_df_full %>%
  group_by(site) %>%
  sample_frac(stage2_frac) %>%
  ungroup()

log_progress(sprintf("✓ Training data: %s CHM, %s DTM observations",
                     format(nrow(mod_chm_s2_train), big.mark = ","),
                     format(nrow(mod_dtm_s2_train), big.mark = ",")))

# Get factor levels from training
train_lc_levels <- levels(droplevels(mod_chm_s2_train$lc_l1_code))
train_site_levels <- levels(droplevels(mod_chm_s2_train$site))
train_eco_levels <- levels(droplevels(mod_chm_s2_train$ecoregion))

# Create holdout sets
log_progress("Creating holdout datasets...")

train_chm_rows <- mod_chm_s2_train$original_row
train_dtm_rows <- mod_dtm_s2_train$original_row

chm_holdout_raw <- chm_df_full %>%
  filter(!original_row %in% train_chm_rows) %>%
  select(-original_row)

dtm_holdout_raw <- dtm_df_full %>%
  filter(!original_row %in% train_dtm_rows) %>%
  select(-original_row)

# Filter holdout to match training factor levels
chm_holdout <- chm_holdout_raw %>%
  filter(lc_l1_code %in% train_lc_levels,
         site %in% train_site_levels,
         ecoregion %in% train_eco_levels) %>%
  mutate(lc_l1_code = factor(lc_l1_code, levels = train_lc_levels),
         site = factor(site, levels = train_site_levels),
         ecoregion = factor(ecoregion, levels = train_eco_levels))

dtm_holdout <- dtm_holdout_raw %>%
  filter(lc_l1_code %in% train_lc_levels,
         site %in% train_site_levels,
         ecoregion %in% train_eco_levels) %>%
  mutate(lc_l1_code = factor(lc_l1_code, levels = train_lc_levels),
         site = factor(site, levels = train_site_levels),
         ecoregion = factor(ecoregion, levels = train_eco_levels))

log_progress(sprintf("✓ CHM holdout: %s observations (%.1f%%)",
                     format(nrow(chm_holdout), big.mark = ","),
                     100 * nrow(chm_holdout) / nrow(chm_df_full)))
log_progress(sprintf("✓ DTM holdout: %s observations (%.1f%%)",
                     format(nrow(dtm_holdout), big.mark = ","),
                     100 * nrow(dtm_holdout) / nrow(dtm_df_full)))

# =====================================================================
# APPLY SITE EXCLUSIONS (must match model training!)
# =====================================================================
# Site 10 was excluded from DTM training due to +143m datum offset
# We must exclude it from holdout predictions too

SITES_EXCLUDE_CHM <- if (exists("SITES_EXCLUDE_CHM")) SITES_EXCLUDE_CHM else ""
SITES_EXCLUDE_DTM <- if (exists("SITES_EXCLUDE_DTM")) SITES_EXCLUDE_DTM else "10"

log_progress("Applying site exclusions to match model training...")

if (nzchar(SITES_EXCLUDE_CHM)) {
 sites_excl_chm <- trimws(unlist(strsplit(SITES_EXCLUDE_CHM, ",")))
  chm_holdout <- chm_holdout %>% filter(!site %in% sites_excl_chm)
  log_progress(sprintf("  CHM: Excluded sites %s", SITES_EXCLUDE_CHM))
}

if (nzchar(SITES_EXCLUDE_DTM)) {
  sites_excl_dtm <- trimws(unlist(strsplit(SITES_EXCLUDE_DTM, ",")))
  dtm_holdout <- dtm_holdout %>% filter(!site %in% sites_excl_dtm)
  log_progress(sprintf("  DTM: Excluded sites %s", SITES_EXCLUDE_DTM))
}

log_progress(sprintf("✓ After exclusions: CHM holdout = %s, DTM holdout = %s",
                     format(nrow(chm_holdout), big.mark = ","),
                     format(nrow(dtm_holdout), big.mark = ",")))

# Drop unused factor levels after exclusions
chm_holdout <- chm_holdout %>%
  mutate(site = droplevels(site),
         lc_l1_code = droplevels(lc_l1_code),
         ecoregion = droplevels(ecoregion))

dtm_holdout <- dtm_holdout %>%
  mutate(site = droplevels(site),
         lc_l1_code = droplevels(lc_l1_code),
         ecoregion = droplevels(ecoregion))

# Prepare holdout data with eco_on flag
mod_chm_s2 <- chm_holdout %>% mutate(eco_on = 1)
mod_dtm_s2 <- dtm_holdout %>% mutate(eco_on = 1)

# Clean up
chm_df_full <- chm_df_full %>% select(-original_row)
dtm_df_full <- dtm_df_full %>% select(-original_row)
mod_chm_s2_train <- mod_chm_s2_train %>% select(-original_row)
mod_dtm_s2_train <- mod_dtm_s2_train %>% select(-original_row)

log_progress("\n*** Using HOLDOUT data for out-of-sample validation ***")

# =====================================================================
# HELPER FUNCTION: STRATIFIED SAMPLING
# =====================================================================

stratified_sample <- function(data, group_var, n_per_group, seed = 2025) {
  set.seed(seed)
  
  sample_sizes <- data %>%
    group_by(!!sym(group_var)) %>%
    summarize(n_available = n(), .groups = "drop") %>%
    mutate(n_sample = pmin(n_per_group, n_available))
  
  sampled_data <- data %>%
    left_join(sample_sizes, by = group_var) %>%
    group_by(!!sym(group_var)) %>%
    group_split() %>%
    lapply(function(df) {
      n_to_sample <- df$n_sample[1]
      df %>%
        select(-n_available, -n_sample) %>%
        slice_sample(n = n_to_sample, replace = FALSE)
    }) %>%
    bind_rows()
  
  return(list(data = sampled_data, sizes = sample_sizes))
}

# =====================================================================
# CHM POSTERIOR PREDICTIVE CHECKS
# =====================================================================

log_progress("\n" %>% paste0(rep("=", 70) %>% paste(collapse = "")))
log_progress("CHM POSTERIOR PREDICTIVE CHECKS")
log_progress(rep("=", 70) %>% paste(collapse = ""))

# Filter to forested classes
chm_forest <- mod_chm_s2 %>%
  filter(lc_l1_code %in% forested_classes) %>%
  filter(is.finite(chm_error_mean))

log_progress(sprintf("CHM: %d total observations in forested classes", nrow(chm_forest)))

# Stratified sampling
chm_sampled <- stratified_sample(chm_forest, "lc_l1_code", n_samples)
chm_sample <- chm_sampled$data
log_progress(sprintf("CHM: Sampled %d observations", nrow(chm_sample)))

# Generate posterior predictions
log_progress("Generating CHM posterior predictions...")
chm_yrep <- posterior_predict(fit_chm_s2, newdata = chm_sample, ndraws = n_draws,
                               allow_new_levels = TRUE)
log_progress(sprintf("✓ Generated %d posterior draws for %d observations",
                     nrow(chm_yrep), ncol(chm_yrep)))

# Create CHM plots
chm_y <- chm_sample$chm_error_mean
chm_plots <- list()

for (lc in forested_classes) {
  lc_idx <- which(chm_sample$lc_l1_code == lc)
  lc_name <- lc_names[lc]
  
  y_obs <- chm_y[lc_idx]
  y_rep_matrix <- chm_yrep[, lc_idx, drop = FALSE]
  
  # Remove non-finite values
 finite_obs <- is.finite(y_obs)
  finite_rep <- apply(y_rep_matrix, 2, function(col) all(is.finite(col)))
  keep_idx <- finite_obs & finite_rep
  
  y_obs <- y_obs[keep_idx]
  y_rep_matrix <- y_rep_matrix[, keep_idx, drop = FALSE]
  
  if (length(y_obs) > 0 && ncol(y_rep_matrix) > 0) {
    n_plot <- min(n_draws_plot, nrow(y_rep_matrix))
    
    p <- bayesplot::ppc_dens_overlay(y_obs, y_rep_matrix[1:n_plot, ]) +
      ggtitle(paste0(lc_name, " (n=", format(length(y_obs), big.mark = ","), ")")) +
      coord_cartesian(xlim = c(-25, 25)) +
      labs(x = "CHM Error (m)", y = "Density") +
      theme_minimal(base_size = 10) +
      theme(
        legend.position = "none",
        plot.title = element_text(face = "bold", size = 9, hjust = 0.5),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        plot.margin = margin(5, 5, 5, 5)
      )
    
    chm_plots[[lc]] <- p
    log_progress(sprintf("  ✓ CHM %s: %d observations", lc, length(y_obs)))
  } else {
    log_progress(sprintf("  ⚠ Skipping CHM %s - insufficient data", lc))
  }
}

# =====================================================================
# DTM POSTERIOR PREDICTIVE CHECKS
# =====================================================================

log_progress("\n" %>% paste0(rep("=", 70) %>% paste(collapse = "")))
log_progress("DTM POSTERIOR PREDICTIVE CHECKS")
log_progress(rep("=", 70) %>% paste(collapse = ""))

# Filter to forested classes
dtm_forest <- mod_dtm_s2 %>%
  filter(lc_l1_code %in% forested_classes) %>%
  filter(is.finite(dtm_error_mean))

log_progress(sprintf("DTM: %d total observations in forested classes", nrow(dtm_forest)))

# Stratified sampling
dtm_sampled <- stratified_sample(dtm_forest, "lc_l1_code", n_samples)
dtm_sample <- dtm_sampled$data
log_progress(sprintf("DTM: Sampled %d observations", nrow(dtm_sample)))

# Generate posterior predictions
log_progress("Generating DTM posterior predictions...")
dtm_yrep <- posterior_predict(fit_dtm_s2, newdata = dtm_sample, ndraws = n_draws,
                               allow_new_levels = TRUE)
log_progress(sprintf("✓ Generated %d posterior draws for %d observations",
                     nrow(dtm_yrep), ncol(dtm_yrep)))

# Create DTM plots
dtm_y <- dtm_sample$dtm_error_mean
dtm_plots <- list()

for (lc in forested_classes) {
  lc_idx <- which(dtm_sample$lc_l1_code == lc)
  lc_name <- lc_names[lc]
  
  y_obs <- dtm_y[lc_idx]
  y_rep_matrix <- dtm_yrep[, lc_idx, drop = FALSE]
  
  # Remove non-finite values
  finite_obs <- is.finite(y_obs)
  finite_rep <- apply(y_rep_matrix, 2, function(col) all(is.finite(col)))
  keep_idx <- finite_obs & finite_rep
  
  y_obs <- y_obs[keep_idx]
  y_rep_matrix <- y_rep_matrix[, keep_idx, drop = FALSE]
  
  if (length(y_obs) > 0 && ncol(y_rep_matrix) > 0) {
    n_plot <- min(n_draws_plot, nrow(y_rep_matrix))
    
    p <- bayesplot::ppc_dens_overlay(y_obs, y_rep_matrix[1:n_plot, ]) +
      ggtitle(paste0(lc_name, " (n=", format(length(y_obs), big.mark = ","), ")")) +
      coord_cartesian(xlim = c(-25, 25)) +
      labs(x = "DTM Error (m)", y = "Density") +
      theme_minimal(base_size = 10) +
      theme(
        legend.position = "none",
        plot.title = element_text(face = "bold", size = 9, hjust = 0.5),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        plot.margin = margin(5, 5, 5, 5)
      )
    
    dtm_plots[[lc]] <- p
    log_progress(sprintf("  ✓ DTM %s: %d observations", lc, length(y_obs)))
  } else {
    log_progress(sprintf("  ⚠ Skipping DTM %s - insufficient data", lc))
  }
}

# =====================================================================
# CREATE COMBINED 2×4 PANEL FIGURE
# =====================================================================

log_progress("\nCreating combined 2×4 panel figure...")

# Arrange CHM plots in order
chm_row <- plot_grid(
  chm_plots[["BDF"]], chm_plots[["DNF"]], 
  chm_plots[["EBF"]], chm_plots[["ENF"]],
  nrow = 1,
  labels = c("A", "B", "C", "D"),
  label_size = 10,
  label_fontface = "bold"
)

# Arrange DTM plots in order
dtm_row <- plot_grid(
  dtm_plots[["BDF"]], dtm_plots[["DNF"]], 
  dtm_plots[["EBF"]], dtm_plots[["ENF"]],
  nrow = 1,
  labels = c("E", "F", "G", "H"),
  label_size = 10,
  label_fontface = "bold"
)

# Add row labels
chm_row_labeled <- plot_grid(
  ggdraw() + draw_label("CHM", fontface = "bold", size = 11, angle = 90),
  chm_row,
  nrow = 1,
  rel_widths = c(0.03, 1)
)

dtm_row_labeled <- plot_grid(
  ggdraw() + draw_label("DTM", fontface = "bold", size = 11, angle = 90),
  dtm_row,
  nrow = 1,
  rel_widths = c(0.03, 1)
)

# Combine rows
combined_plot <- plot_grid(
  chm_row_labeled,
  dtm_row_labeled,
  nrow = 2,
  rel_heights = c(1, 1)
)

# Add overall title
final_plot <- plot_grid(
  ggdraw() + draw_label(
    "Posterior Predictive Checks by Forest Type (Holdout Data)",
    fontface = "bold", size = 12
  ),
  combined_plot,
  ncol = 1,
  rel_heights = c(0.05, 1)
)

# Save combined figure
output_file <- file.path(out_plots, "ppc_by_forest_type_combined.pdf")
ggsave(output_file, final_plot, width = 12, height = 6, bg = "white")
log_progress(sprintf("✓ Saved combined figure: %s", output_file))

# Also save as PNG for quick viewing
output_png <- file.path(out_plots, "ppc_by_forest_type_combined.png")
ggsave(output_png, final_plot, width = 12, height = 6, dpi = 300, bg = "white")
log_progress(sprintf("✓ Saved PNG: %s", output_png))

# =====================================================================
# SAVE INDIVIDUAL PLOTS (for flexibility)
# =====================================================================

# CHM individual
if (length(chm_plots) > 0) {
  p_chm <- plot_grid(plotlist = chm_plots, ncol = 2)
  ggsave(file.path(out_plots, "ppc_chm_by_forest_type.pdf"),
         p_chm, width = 10, height = 8, bg = "white")
  log_progress("✓ Saved individual CHM plots")
}

# DTM individual
if (length(dtm_plots) > 0) {
  p_dtm <- plot_grid(plotlist = dtm_plots, ncol = 2)
  ggsave(file.path(out_plots, "ppc_dtm_by_forest_type.pdf"),
         p_dtm, width = 10, height = 8, bg = "white")
  log_progress("✓ Saved individual DTM plots")
}

# =====================================================================
# SUMMARY
# =====================================================================

log_progress("\n" %>% paste0(rep("=", 70) %>% paste(collapse = "")))
log_progress("POSTERIOR PREDICTIVE CHECKS COMPLETED")
log_progress(rep("=", 70) %>% paste(collapse = ""))

log_progress("\nKey outputs:")
log_progress(sprintf("  Combined figure: %s", output_file))
log_progress(sprintf("  Individual CHM:  %s", file.path(out_plots, "ppc_chm_by_forest_type.pdf")))
log_progress(sprintf("  Individual DTM:  %s", file.path(out_plots, "ppc_dtm_by_forest_type.pdf")))

log_progress("\nConfiguration used:")
log_progress(sprintf("  Stage 2 subsample: %.0f%%", 100 * stage2_frac))
log_progress(sprintf("  Holdout fraction:  %.0f%%", 100 * (1 - stage2_frac)))
log_progress(sprintf("  Posterior draws:   %d (plotted: %d)", n_draws, n_draws_plot))
log_progress(sprintf("  Samples per class: up to %d", n_samples))

log_progress("\nSuggested figure caption:")
cat("\n")
cat("Figure X. Posterior predictive checks stratified by forest functional type\n")
cat("for CHM (top row, A-D) and DTM (bottom row, E-F) models. Observed holdout\n")
cat("error distributions (dark line) are compared against replicated datasets\n")
cat("drawn from the posterior predictive distribution (light lines). Forest\n")
cat("types: BDF = Broadleaf Deciduous, DNF = Deciduous Needleleaf, EBF =\n")
cat("Evergreen Broadleaf, ENF = Evergreen Needleleaf. Sample sizes shown in\n")
cat("parentheses represent the holdout validation set.\n")
cat("\n")

log_progress(rep("=", 70) %>% paste(collapse = ""))
