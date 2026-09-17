# =====================================================================
# generate_ppc_overall.R
#
# PURPOSE:
#   Generate standalone PPC density overlay figure showing observed
#   vs. posterior predictive distributions for CHM and DTM overall
#   (not stratified by land cover). This is the "Figure 8" replacement.
#
# OUTPUT:
#   - figure_ppc_overall.pdf/png (1×2 panel: CHM left, DTM right)
#
# RUN: source("generate_ppc_overall.R")
# =====================================================================

Sys.setenv(DISPLAY = "")
options(device = pdf)

# =====================================================================
# LOAD DEPENDENCIES
# =====================================================================

config_path <- "analysis_config.R"
utils_path <- "analysis_utils.R"

if (file.exists(config_path)) {
  source(config_path)
} else {
  suppressPackageStartupMessages({
    library(brms)
    library(dplyr)
    library(ggplot2)
    library(bayesplot)
    library(cowplot)
  })
  
  PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/gpfs/data1/vclgp/lmaden/chpt1")
  out_plots <- file.path(PROJECT_ROOT, "plots")
  dir.create(out_plots, showWarnings = FALSE, recursive = TRUE)
  theme_set(theme_cowplot())
}

if (file.exists(utils_path)) {
  source(utils_path)
} else {
  log_progress <- function(msg) cat(sprintf("[%s] %s\n", Sys.time(), msg))
  CHECKPOINT_DIR <- file.path(PROJECT_ROOT, "checkpoints")
}

# =====================================================================
# CONFIGURATION
# =====================================================================

# Number of holdout observations to use for PPC
n_samples <- 10000

# Number of posterior draws for predictions
n_draws <- 500

# Number of posterior draws to plot (for visibility)
n_draws_plot <- 50

# Stage 2 subsample fraction (must match model training)
stage2_frac <- as.numeric(Sys.getenv("STAGE2_FRAC", "0.33"))

# Seed for reproducibility
set.seed(2025)

log_progress(strrep("=", 70))
log_progress("OVERALL POSTERIOR PREDICTIVE CHECK FIGURE")
log_progress(strrep("=", 70))

# =====================================================================
# LOAD MODELS
# =====================================================================

log_progress("Loading Stage 2 models...")

model_checkpoint <- file.path(CHECKPOINT_DIR, "10_models_stage2.rds")

if (file.exists(model_checkpoint)) {
  data_models <- readRDS(model_checkpoint)
  fit_chm_s2 <- data_models$data$fit_chm_s2
  fit_dtm_s2 <- data_models$data$fit_dtm_s2
  log_progress(sprintf("  Models loaded from: %s", model_checkpoint))
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
log_progress(sprintf("  Full dataset: %s CHM, %s DTM",
                     format(nrow(chm_df_full_raw), big.mark = ","),
                     format(nrow(dtm_df_full_raw), big.mark = ",")))

# Quality filters (same as section_08)
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

# Get factor levels from training
train_lc_levels_chm <- levels(droplevels(mod_chm_s2_train$lc_l1_code))
train_site_levels_chm <- levels(droplevels(mod_chm_s2_train$site))
train_eco_levels_chm <- levels(droplevels(mod_chm_s2_train$ecoregion))

train_lc_levels_dtm <- levels(droplevels(mod_dtm_s2_train$lc_l1_code))
train_site_levels_dtm <- levels(droplevels(mod_dtm_s2_train$site))
train_eco_levels_dtm <- levels(droplevels(mod_dtm_s2_train$ecoregion))

# Create holdout sets
train_chm_rows <- mod_chm_s2_train$original_row
train_dtm_rows <- mod_dtm_s2_train$original_row

chm_holdout <- chm_df_full %>%
  filter(!original_row %in% train_chm_rows) %>%
  select(-original_row) %>%
  filter(lc_l1_code %in% train_lc_levels_chm,
         site %in% train_site_levels_chm,
         ecoregion %in% train_eco_levels_chm) %>%
  mutate(lc_l1_code = factor(lc_l1_code, levels = train_lc_levels_chm),
         site = factor(site, levels = train_site_levels_chm),
         ecoregion = factor(ecoregion, levels = train_eco_levels_chm))

dtm_holdout <- dtm_df_full %>%
  filter(!original_row %in% train_dtm_rows) %>%
  select(-original_row) %>%
  filter(lc_l1_code %in% train_lc_levels_dtm,
         site %in% train_site_levels_dtm,
         ecoregion %in% train_eco_levels_dtm) %>%
  mutate(lc_l1_code = factor(lc_l1_code, levels = train_lc_levels_dtm),
         site = factor(site, levels = train_site_levels_dtm),
         ecoregion = factor(ecoregion, levels = train_eco_levels_dtm))

log_progress(sprintf("  CHM holdout: %s observations", format(nrow(chm_holdout), big.mark = ",")))
log_progress(sprintf("  DTM holdout: %s observations", format(nrow(dtm_holdout), big.mark = ",")))

# =====================================================================
# APPLY SITE EXCLUSIONS
# =====================================================================

# Site 10 excluded from DTM model
SITES_EXCLUDE_DTM <- "10"

sites_excl_dtm <- trimws(unlist(strsplit(SITES_EXCLUDE_DTM, ",")))
dtm_holdout <- dtm_holdout %>% filter(!site %in% sites_excl_dtm)
log_progress(sprintf("  DTM holdout after excluding Site 10: %s observations",
                     format(nrow(dtm_holdout), big.mark = ",")))

# =====================================================================
# SUBSAMPLE HOLDOUT FOR PPC
# =====================================================================

log_progress(sprintf("Subsampling %s observations per product...",
                     format(n_samples, big.mark = ",")))

# Stratified by site to ensure representation
chm_sample <- chm_holdout %>%
  group_by(site) %>%
  slice_sample(n = ceiling(n_samples / n_distinct(chm_holdout$site)),
               replace = FALSE) %>%
  ungroup() %>%
  slice_sample(n = min(n_samples, nrow(.)))

dtm_sample <- dtm_holdout %>%
  group_by(site) %>%
  slice_sample(n = ceiling(n_samples / n_distinct(dtm_holdout$site)),
               replace = FALSE) %>%
  ungroup() %>%
  slice_sample(n = min(n_samples, nrow(.)))

log_progress(sprintf("  CHM sample: %d observations from %d sites",
                     nrow(chm_sample), n_distinct(chm_sample$site)))
log_progress(sprintf("  DTM sample: %d observations from %d sites",
                     nrow(dtm_sample), n_distinct(dtm_sample$site)))

# =====================================================================
# GENERATE POSTERIOR PREDICTIONS
# =====================================================================

log_progress("Generating CHM posterior predictions...")
chm_yrep <- posterior_predict(fit_chm_s2, newdata = chm_sample,
                               ndraws = n_draws, allow_new_levels = TRUE)
log_progress(sprintf("  CHM: %d draws x %d observations", nrow(chm_yrep), ncol(chm_yrep)))

log_progress("Generating DTM posterior predictions...")
dtm_yrep <- posterior_predict(fit_dtm_s2, newdata = dtm_sample,
                               ndraws = n_draws, allow_new_levels = TRUE)
log_progress(sprintf("  DTM: %d draws x %d observations", nrow(dtm_yrep), ncol(dtm_yrep)))

# =====================================================================
# COMPUTE COVERAGE STATISTICS
# =====================================================================

log_progress("Computing 95% PPC coverage...")

compute_coverage <- function(y_obs, y_rep, alpha = 0.95) {
  lower <- (1 - alpha) / 2
  upper <- 1 - lower
  pi_lower <- apply(y_rep, 2, quantile, probs = lower, na.rm = TRUE)
  pi_upper <- apply(y_rep, 2, quantile, probs = upper, na.rm = TRUE)
  within <- (y_obs >= pi_lower) & (y_obs <= pi_upper)
  return(mean(within, na.rm = TRUE))
}

chm_y_obs <- chm_sample$chm_error_mean
dtm_y_obs <- dtm_sample$dtm_error_mean

cov_chm <- compute_coverage(chm_y_obs, chm_yrep)
cov_dtm <- compute_coverage(dtm_y_obs, dtm_yrep)

log_progress(sprintf("  CHM 95%% coverage: %.1f%%", cov_chm * 100))
log_progress(sprintf("  DTM 95%% coverage: %.1f%%", cov_dtm * 100))

# =====================================================================
# CREATE PPC DENSITY OVERLAY PLOTS
# =====================================================================

log_progress("Creating PPC density overlay plots...")

# Clean non-finite values
chm_finite <- is.finite(chm_y_obs) & apply(chm_yrep, 2, function(col) all(is.finite(col)))
dtm_finite <- is.finite(dtm_y_obs) & apply(dtm_yrep, 2, function(col) all(is.finite(col)))

chm_y_clean <- chm_y_obs[chm_finite]
chm_yrep_clean <- chm_yrep[, chm_finite, drop = FALSE]

dtm_y_clean <- dtm_y_obs[dtm_finite]
dtm_yrep_clean <- dtm_yrep[, dtm_finite, drop = FALSE]

n_plot_chm <- min(n_draws_plot, nrow(chm_yrep_clean))
n_plot_dtm <- min(n_draws_plot, nrow(dtm_yrep_clean))

# CHM PPC plot
p_chm <- ppc_dens_overlay(chm_y_clean, chm_yrep_clean[1:n_plot_chm, ]) +
  coord_cartesian(xlim = c(-30, 30)) +
  labs(
    title = sprintf("CHM (n = %s)", format(length(chm_y_clean), big.mark = ",")),
    subtitle = sprintf("95%% coverage: %.1f%%", cov_chm * 100),
    x = "CHM Error (m)",
    y = "Density"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40"),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    panel.grid.minor = element_blank()
  )

# DTM PPC plot
p_dtm <- ppc_dens_overlay(dtm_y_clean, dtm_yrep_clean[1:n_plot_dtm, ]) +
  coord_cartesian(xlim = c(-30, 30)) +
  labs(
    title = sprintf("DTM (n = %s)", format(length(dtm_y_clean), big.mark = ",")),
    subtitle = sprintf("95%% coverage: %.1f%%", cov_dtm * 100),
    x = "DTM Error (m)",
    y = "Density"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40"),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    panel.grid.minor = element_blank()
  )

# =====================================================================
# COMBINE INTO 1×2 PANEL FIGURE
# =====================================================================

log_progress("Assembling combined figure...")

combined <- plot_grid(
  p_chm, p_dtm,
  nrow = 1,
  labels = c("(a)", "(b)"),
  label_size = 12,
  label_fontface = "bold"
)

# Add overall title
final_plot <- plot_grid(
  ggdraw() + draw_label(
    "Posterior Predictive Checks",
    fontface = "bold", size = 14
  ) +
  draw_label(
    "Observed error distribution (dark line) vs. posterior predictive draws (light lines)",
    size = 10, y = 0.2, color = "gray40"
  ),
  combined,
  ncol = 1,
  rel_heights = c(0.08, 1)
)

# =====================================================================
# SAVE
# =====================================================================

log_progress("Saving figures...")

ggsave(file.path(out_plots, "figure_ppc_overall.pdf"),
       final_plot, width = 10, height = 5, bg = "white", device = pdf)
log_progress(sprintf("  Saved: %s/figure_ppc_overall.pdf", out_plots))

ggsave(file.path(out_plots, "figure_ppc_overall.png"),
       final_plot, width = 10, height = 5, dpi = 300, bg = "white")
log_progress(sprintf("  Saved: %s/figure_ppc_overall.png", out_plots))

# =====================================================================
# SUMMARY
# =====================================================================

log_progress(strrep("=", 70))
log_progress("PPC OVERALL FIGURE COMPLETE")
log_progress(strrep("=", 70))

log_progress(sprintf("\n  CHM: %d observations, %.1f%% coverage",
                     length(chm_y_clean), cov_chm * 100))
log_progress(sprintf("  DTM: %d observations, %.1f%% coverage",
                     length(dtm_y_clean), cov_dtm * 100))

cat("\n")
cat("SUGGESTED CAPTION:\n\n")
cat("Figure 8. Posterior predictive checks for the CHM (a) and DTM (b) models.\n")
cat("Observed holdout error distributions (dark line) are compared against 50\n")
cat("replicated datasets drawn from the posterior predictive distribution (light\n")
cat("lines). The close agreement between observed and predicted distributions\n")
cat("confirms adequate model fit. Coverage of 95% posterior predictive intervals\n")
cat("is annotated on each panel.\n")

log_progress(strrep("=", 70))
