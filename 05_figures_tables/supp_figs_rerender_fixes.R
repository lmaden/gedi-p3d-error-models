#!/usr/bin/env Rscript
# =====================================================================
# supp_figs_rerender_fixes.R
#
# Patches three issues from supp_figs_rerender.R:
#
#   (1) figS06: data source corrected
#       Original used fit_chm$data (training-only); v04 used the full
#       chm_df sampled to 50,000 rows. Mirror that so the R^2 values
#       are comparable apples-to-apples with v04. If they still come
#       in much lower than v04 (~0.5-0.6), that's a real finding
#       about the 16-site refit; if they come up, my prior data-source
#       diagnosis was right.
#
#   (2) figS18: validate_newdata "New factor levels not allowed" fix
#       The DTM holdout contains 'LMS' (and possibly other) land cover
#       levels not in fit_dtm$data. Likely cause: 01_data_ingest was
#       updated post-model-fit. Filter PPC samples to only factor
#       levels present in the model BEFORE posterior_predict.
#
#   (3) figS20: re-run (was blocked by figS18 error in v1)
#
# figS12 ran cleanly in v1 - SKIPPED here, keep existing PDF.
#
# RUN: from /gpfs/data1/vclgp/lmaden/chpt1
#   Rscript supp_figs_rerender_fixes.R
#
# =====================================================================

Sys.setenv(DISPLAY = "")
options(device = pdf)

suppressPackageStartupMessages({
  library(brms)
  library(bayesplot)
  library(ggplot2)
  library(cowplot)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(scales)
  library(viridis)
  library(data.table)
  library(sf)
  library(gstat)
})

bayesplot::bayesplot_theme_set(bayesplot::theme_default())

if (!exists("log_progress")) source("analysis_utils.R")
if (!exists("PROJECT_ROOT")) source("analysis_config.R")

stopifnot(exists("PROJECT_ROOT"), exists("CHECKPOINT_DIR"))

out_dir <- file.path(PROJECT_ROOT, "plots", "groundwork")
manuscript_tables_dir <- file.path(PROJECT_ROOT, "manuscript_tables")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(manuscript_tables_dir, showWarnings = FALSE, recursive = TRUE)

log_progress(strrep("=", 70))
log_progress("Supplementary figure fixes: figS06 + figS18 + figS20")
log_progress(strrep("=", 70))

# ---------------------------------------------------------------------
# Load checkpoints (same as v1)
# ---------------------------------------------------------------------

log_progress("Loading 16-site CHM refit...")
cp_chm <- load_checkpoint("10b_chm_sensitivity")
fit_chm <- cp_chm$fit_chm_16
sites_excluded_chm <- cp_chm$sites_excluded
log_progress(sprintf("  CHM model: %d obs across %d sites; excludes %s",
                     nrow(fit_chm$data),
                     length(unique(as.character(fit_chm$data$site))),
                     paste(sites_excluded_chm, collapse = ", ")))

log_progress("Loading 18-site DTM...")
cp_s2 <- load_checkpoint("10_models_stage2")
fit_dtm <- cp_s2$fit_dtm_s2
log_progress(sprintf("  DTM model: %d obs across %d sites",
                     nrow(fit_dtm$data),
                     length(unique(as.character(fit_dtm$data$site)))))

log_progress("Loading full dataset (01_data_ingest)...")
data_chk <- load_checkpoint("01_data_ingest")
chm_df_full_raw <- data_chk$chm_df
dtm_df_full_raw <- data_chk$dtm_df

FOREST_CLASSES <- c("BDF", "DNF", "EBF", "ENF")
FOREST_LABELS <- c(
  "BDF" = "Broadleaf Deciduous",
  "DNF" = "Deciduous Needleleaf",
  "EBF" = "Evergreen Broadleaf",
  "ENF" = "Evergreen Needleleaf"
)

theme_pub <- function(base_size = 10) {
  theme_cowplot(font_size = base_size) +
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      axis.line  = element_line(color = "black", linewidth = 0.5),
      axis.text  = element_text(color = "black", size = base_size - 1),
      axis.title = element_text(color = "black", size = base_size),
      plot.title    = element_text(size = base_size + 1, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size - 1, hjust = 0, color = "gray40")
    )
}

# ---------------------------------------------------------------------
# Helper: filter newdata to factor levels present in the model
# Mirrors the factor-compatibility approach used in
# section_12_spatial.R's extract_brms_residuals (lines 136-202).
# ---------------------------------------------------------------------

filter_to_model_levels <- function(newdata, model, name) {
  model_data <- model$data
  factor_vars <- names(model_data)[sapply(model_data, is.factor)]
  valid_idx <- rep(TRUE, nrow(newdata))
  for (fvar in factor_vars) {
    if (fvar %in% names(newdata)) {
      model_levels <- levels(model_data[[fvar]])
      bad <- setdiff(unique(as.character(newdata[[fvar]])), model_levels)
      if (length(bad) > 0) {
        n_dropped <- sum(as.character(newdata[[fvar]]) %in% bad)
        log_progress(sprintf(
          "    %s: dropping %s rows with unknown %s level(s): %s",
          name, format(n_dropped, big.mark = ","),
          fvar, paste(bad, collapse = ", ")))
      }
      valid_idx <- valid_idx & (as.character(newdata[[fvar]]) %in% model_levels)
    }
  }
  out <- newdata[valid_idx, ]
  for (fvar in factor_vars) {
    if (fvar %in% names(out)) {
      out[[fvar]] <- factor(out[[fvar]], levels = levels(model_data[[fvar]]))
    }
  }
  out
}

# =====================================================================
# figS06 v2 — Predicted vs Observed, FULL-data sampling approach
# =====================================================================

log_progress("")
log_progress(strrep("-", 70))
log_progress("figS06 v2 — full-data sampling (mirrors section_15 approach)")
log_progress(strrep("-", 70))

# Prep CHM full dataset; filter to 16 retained sites
chm_df_full <- chm_df_full_raw %>%
  filter(is.finite(chm_error_mean),
         !is.na(lc_l1_code), !is.na(site), !is.na(ecoregion)) %>%
  filter(!as.character(site) %in% sites_excluded_chm) %>%
  mutate(lc_l1_code = factor(lc_l1_code),
         ecoregion  = factor(ecoregion),
         site       = factor(site))

log_progress(sprintf("  Full 16-site CHM data: %s rows",
                     format(nrow(chm_df_full), big.mark = ",")))

# Sample 50,000 rows (matches section_15 n_sample = 50000, seed = 42)
set.seed(42)
chm_sampled <- chm_df_full %>%
  sample_n(min(50000, nrow(.)))

log_progress(sprintf("  Sampled to %s rows",
                     format(nrow(chm_sampled), big.mark = ",")))

# Filter to model-compatible factor levels (defensive)
chm_sampled <- filter_to_model_levels(chm_sampled, fit_chm, "CHM (figS06)")

# Restrict to forest classes
chm_forest <- as.data.table(chm_sampled)[lc_l1_code %in% FOREST_CLASSES]
log_progress(sprintf("  Forest-class subset: %s rows",
                     format(nrow(chm_forest), big.mark = ",")))

log_progress("  Computing fitted() (this is the slow step)...")
fitted_full <- fitted(fit_chm, newdata = chm_forest,
                      summary = TRUE, ndraws = 100)
chm_forest[, predicted := fitted_full[, "Estimate"]]
chm_forest[, observed  := chm_error_mean]

# Per-forest-type anchors, BOTH R^2 flavors for transparency
r2_anchors <- chm_forest[, .(
  n      = .N,
  r2_cor = cor(observed, predicted, use = "complete.obs")^2,
  r2_ss  = 1 - sum((observed - predicted)^2, na.rm = TRUE) /
               sum((observed - mean(observed, na.rm = TRUE))^2, na.rm = TRUE),
  rmse   = sqrt(mean((observed - predicted)^2, na.rm = TRUE)),
  mae    = mean(abs(observed - predicted), na.rm = TRUE),
  bias   = mean(predicted - observed, na.rm = TRUE)
), by = lc_l1_code]
setorder(r2_anchors, lc_l1_code)

log_progress("  Per-forest-type anchors (v04 reference in parentheses):")
v04_ref <- list(BDF = 0.58, DNF = 0.57, EBF = 0.54, ENF = 0.29)
for (rr in seq_len(nrow(r2_anchors))) {
  cls <- as.character(r2_anchors$lc_l1_code[rr])
  log_progress(sprintf(
    "    %s: cor^2 = %.3f (v04 ~%.2f) | SS R^2 = %.3f | RMSE = %.2f m | n = %s",
    cls,
    r2_anchors$r2_cor[rr],
    if (!is.null(v04_ref[[cls]])) v04_ref[[cls]] else NA_real_,
    r2_anchors$r2_ss[rr],
    r2_anchors$rmse[rr],
    format(r2_anchors$n[rr], big.mark = ",")))
}

fwrite(r2_anchors,
       file.path(manuscript_tables_dir,
                 "section_L5_figS06_r2_anchors_v2_fulldata.csv"))

# Per-panel label: v04 format
label_df <- r2_anchors[, .(
  lc_l1_code,
  lab = sprintf("R\u00b2 = %.2f\nRMSE = %.1f m\nn = %s",
                r2_cor, rmse, format(n, big.mark = ","))
)]

lim_val <- quantile(abs(c(chm_forest$observed, chm_forest$predicted)),
                    0.995, na.rm = TRUE)
lim_val <- ceiling(lim_val)
log_progress(sprintf("  Axis limits: [-%d, %d] m", lim_val, lim_val))

label_df[, x := -lim_val * 0.9]
label_df[, y :=  lim_val * 0.85]

p_s06 <- ggplot(chm_forest, aes(x = observed, y = predicted)) +
  geom_hex(bins = 50, aes(fill = after_stat(count))) +
  scale_fill_viridis_c(trans = "log10", name = "Count",
                       labels = scales::label_comma()) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",
              color = "red", linewidth = 0.6) +
  geom_text(data = label_df,
            aes(x = x, y = y, label = lab),
            hjust = 0, vjust = 1, size = 3, inherit.aes = FALSE) +
  facet_wrap(~ lc_l1_code, ncol = 2,
             labeller = labeller(lc_l1_code = FOREST_LABELS)) +
  coord_fixed(xlim = c(-lim_val, lim_val), ylim = c(-lim_val, lim_val)) +
  labs(
    title = "CHM Predicted vs Observed Error by Forest Type",
    x     = "Observed error (m)",
    y     = "Predicted error (m)"
  ) +
  theme_pub() +
  theme(legend.position = "right")

ggsave(file.path(out_dir, "figS06_pred_vs_obs_by_forest_16site_v2.pdf"),
       p_s06, width = 10, height = 10, bg = "white", device = pdf)
log_progress(sprintf("  Saved: %s/figS06_pred_vs_obs_by_forest_16site_v2.pdf",
                     out_dir))

# =====================================================================
# figS18 v2 — PPC with model-level factor filter
# =====================================================================

log_progress("")
log_progress(strrep("-", 70))
log_progress("figS18 v2 — filter PPC samples to model factor levels")
log_progress(strrep("-", 70))

n_samples       <- 10000
n_draws_ppc     <- 500
n_draws_plot    <- 50
stage2_frac     <- as.numeric(Sys.getenv("STAGE2_FRAC", "0.33"))

chm_full <- chm_df_full_raw %>%
  filter(is.finite(chm_error_mean), !is.na(lc_l1_code),
         !is.na(site), !is.na(ecoregion)) %>%
  mutate(lc_l1_code = factor(lc_l1_code),
         ecoregion  = factor(ecoregion),
         site       = factor(site))
dtm_full <- dtm_df_full_raw %>%
  filter(is.finite(dtm_error_mean), !is.na(lc_l1_code),
         !is.na(site), !is.na(ecoregion)) %>%
  mutate(lc_l1_code = factor(lc_l1_code),
         ecoregion  = factor(ecoregion),
         site       = factor(site))

# CHM filter to 16 sites BEFORE training split
chm_full <- chm_full %>%
  filter(!as.character(site) %in% sites_excluded_chm) %>%
  mutate(site = droplevels(factor(site)))

chm_full$original_row <- seq_len(nrow(chm_full))
dtm_full$original_row <- seq_len(nrow(dtm_full))

set.seed(2025)
chm_train_smp <- chm_full %>% group_by(site) %>%
  sample_frac(stage2_frac) %>% ungroup()
dtm_train_smp <- dtm_full %>% group_by(site) %>%
  sample_frac(stage2_frac) %>% ungroup()

chm_holdout <- chm_full %>%
  filter(!original_row %in% chm_train_smp$original_row) %>%
  select(-original_row)
dtm_holdout <- dtm_full %>%
  filter(!original_row %in% dtm_train_smp$original_row) %>%
  select(-original_row) %>%
  filter(!as.character(site) %in% "10")  # Site 10 excluded from DTM

log_progress(sprintf("  CHM holdout (raw): %s | DTM holdout (raw): %s",
                     format(nrow(chm_holdout), big.mark = ","),
                     format(nrow(dtm_holdout), big.mark = ",")))

# Apply MODEL-LEVEL filter (this is the fix vs v1)
log_progress("  Filtering CHM holdout to model factor levels:")
chm_holdout <- filter_to_model_levels(chm_holdout, fit_chm, "CHM (figS18)")
log_progress("  Filtering DTM holdout to model factor levels:")
dtm_holdout <- filter_to_model_levels(dtm_holdout, fit_dtm, "DTM (figS18)")
log_progress(sprintf("  CHM holdout (model-compatible): %s | DTM: %s",
                     format(nrow(chm_holdout), big.mark = ","),
                     format(nrow(dtm_holdout), big.mark = ",")))

# Per-site stratified subsample to n_samples
chm_ppc_sample <- chm_holdout %>%
  group_by(site) %>%
  slice_sample(n = ceiling(n_samples / n_distinct(chm_holdout$site)),
               replace = FALSE) %>% ungroup() %>%
  slice_sample(n = min(n_samples, nrow(.)))
dtm_ppc_sample <- dtm_holdout %>%
  group_by(site) %>%
  slice_sample(n = ceiling(n_samples / n_distinct(dtm_holdout$site)),
               replace = FALSE) %>% ungroup() %>%
  slice_sample(n = min(n_samples, nrow(.)))

log_progress(sprintf("  CHM PPC sample: %d obs | DTM PPC sample: %d obs",
                     nrow(chm_ppc_sample), nrow(dtm_ppc_sample)))

log_progress("  Generating CHM posterior predictions...")
chm_yrep <- posterior_predict(fit_chm, newdata = chm_ppc_sample,
                              ndraws = n_draws_ppc, allow_new_levels = TRUE)
log_progress("  Generating DTM posterior predictions...")
dtm_yrep <- posterior_predict(fit_dtm, newdata = dtm_ppc_sample,
                              ndraws = n_draws_ppc, allow_new_levels = TRUE)

compute_coverage <- function(y_obs, y_rep, alpha = 0.95) {
  lower <- (1 - alpha) / 2; upper <- 1 - lower
  pi_lo <- apply(y_rep, 2, quantile, probs = lower, na.rm = TRUE)
  pi_hi <- apply(y_rep, 2, quantile, probs = upper, na.rm = TRUE)
  mean((y_obs >= pi_lo) & (y_obs <= pi_hi), na.rm = TRUE)
}

chm_y_obs <- chm_ppc_sample$chm_error_mean
dtm_y_obs <- dtm_ppc_sample$dtm_error_mean
cov_chm <- compute_coverage(chm_y_obs, chm_yrep)
cov_dtm <- compute_coverage(dtm_y_obs, dtm_yrep)

log_progress(sprintf("  CHM 95%% PI coverage: %.2f%%", cov_chm * 100))
log_progress(sprintf("  DTM 95%% PI coverage: %.2f%%", cov_dtm * 100))

fwrite(data.table(
  product       = c("CHM", "DTM"),
  n_sample      = c(length(chm_y_obs), length(dtm_y_obs)),
  coverage_95pi = c(cov_chm, cov_dtm)),
  file.path(manuscript_tables_dir, "section_L5_figS18_coverage_v2.csv"))

chm_finite <- is.finite(chm_y_obs) &
  apply(chm_yrep, 2, function(col) all(is.finite(col)))
dtm_finite <- is.finite(dtm_y_obs) &
  apply(dtm_yrep, 2, function(col) all(is.finite(col)))
chm_y_clean    <- chm_y_obs[chm_finite]
chm_yrep_clean <- chm_yrep[, chm_finite, drop = FALSE]
dtm_y_clean    <- dtm_y_obs[dtm_finite]
dtm_yrep_clean <- dtm_yrep[, dtm_finite, drop = FALSE]

n_plot_chm <- min(n_draws_plot, nrow(chm_yrep_clean))
n_plot_dtm <- min(n_draws_plot, nrow(dtm_yrep_clean))

p_ppc_chm <- ppc_dens_overlay(chm_y_clean,
                              chm_yrep_clean[seq_len(n_plot_chm), ]) +
  coord_cartesian(xlim = c(-30, 30)) +
  labs(
    title    = sprintf("CHM (n = %s)", format(length(chm_y_clean), big.mark = ",")),
    subtitle = sprintf("95%% PI coverage: %.1f%% (full holdout)", cov_chm * 100),
    x = "CHM Error (m)", y = "Density"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position    = "none",
    plot.title         = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle      = element_text(size = 10, hjust = 0.5, color = "gray40"),
    axis.title         = element_text(size = 11),
    axis.text          = element_text(size = 10),
    panel.grid.minor   = element_blank()
  )

p_ppc_dtm <- ppc_dens_overlay(dtm_y_clean,
                              dtm_yrep_clean[seq_len(n_plot_dtm), ]) +
  coord_cartesian(xlim = c(-30, 30)) +
  labs(
    title    = sprintf("DTM (n = %s)", format(length(dtm_y_clean), big.mark = ",")),
    subtitle = sprintf("95%% PI coverage: %.1f%% (full holdout)", cov_dtm * 100),
    x = "DTM Error (m)", y = "Density"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position    = "none",
    plot.title         = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle      = element_text(size = 10, hjust = 0.5, color = "gray40"),
    axis.title         = element_text(size = 11),
    axis.text          = element_text(size = 10),
    panel.grid.minor   = element_blank()
  )

ppc_combined <- plot_grid(
  p_ppc_chm, p_ppc_dtm, nrow = 1,
  labels = c("(a)", "(b)"),
  label_size = 12, label_fontface = "bold"
)

p_s18 <- plot_grid(
  ggdraw() +
    draw_label("Posterior Predictive Checks",
               fontface = "bold", size = 14) +
    draw_label("Observed error distribution (dark line) vs. posterior predictive draws (light lines)",
               size = 10, y = 0.2, color = "gray40"),
  ppc_combined,
  ncol = 1,
  rel_heights = c(0.08, 1)
)

ggsave(file.path(out_dir, "figS18_ppc_overall_16site_v2.pdf"),
       p_s18, width = 10, height = 5, bg = "white", device = pdf)
log_progress(sprintf("  Saved: %s/figS18_ppc_overall_16site_v2.pdf", out_dir))

# =====================================================================
# figS20 v2 — variograms (CHM 16-site + DTM 18-site)
# Same methodology as v1 but with the filter_to_model_levels guard
# baked into the residual-extraction step.
# =====================================================================

log_progress("")
log_progress(strrep("-", 70))
log_progress("figS20 v2 — semi-variograms")
log_progress(strrep("-", 70))

TARGET_CRS       <- 5070
SOURCE_CRS       <- 4326
GEDI_SWATH_WIDTH <- 6000
MIN_LAG_DISTANCE <- 100
MAX_LAG_DISTANCE <- 50000
N_LAGS           <- 20

chm_df_all <- chm_df_full_raw
dtm_df_all <- dtm_df_full_raw
stopifnot(all(c("x", "y") %in% names(chm_df_all)))
stopifnot(all(c("x", "y") %in% names(dtm_df_all)))

# Filter CHM to 16 sites
chm_df_all <- chm_df_all %>%
  filter(!as.character(site) %in% sites_excluded_chm)

set.seed(42)
chm_s <- chm_df_all %>%
  tidyr::drop_na(x, y, chm_error_mean) %>%
  group_by(site) %>%
  group_modify(~ slice_sample(.x, n = min(5000, nrow(.x)), replace = FALSE)) %>%
  ungroup()
dtm_s <- dtm_df_all %>%
  tidyr::drop_na(x, y, dtm_error_mean) %>%
  group_by(site) %>%
  group_modify(~ slice_sample(.x, n = min(5000, nrow(.x)), replace = FALSE)) %>%
  ungroup()
log_progress(sprintf("  CHM spatial sample: %s | DTM spatial sample: %s",
                     format(nrow(chm_s), big.mark = ","),
                     format(nrow(dtm_s), big.mark = ",")))

# Project to Albers
chm_sf <- st_transform(
  st_as_sf(chm_s, coords = c("x", "y"), crs = SOURCE_CRS, remove = FALSE),
  TARGET_CRS)
dtm_sf <- st_transform(
  st_as_sf(dtm_s, coords = c("x", "y"), crs = SOURCE_CRS, remove = FALSE),
  TARGET_CRS)

chm_s <- chm_s %>%
  mutate(x_proj = st_coordinates(chm_sf)[, 1],
         y_proj = st_coordinates(chm_sf)[, 2])
dtm_s <- dtm_s %>%
  mutate(x_proj = st_coordinates(dtm_sf)[, 1],
         y_proj = st_coordinates(dtm_sf)[, 2])

# Mutate to factors as the model expects (so filter_to_model_levels works)
chm_s <- chm_s %>%
  mutate(lc_l1_code = factor(lc_l1_code),
         site       = factor(site),
         ecoregion  = factor(ecoregion))
dtm_s <- dtm_s %>%
  mutate(lc_l1_code = factor(lc_l1_code),
         site       = factor(site),
         ecoregion  = factor(ecoregion))

# Filter to model-compatible levels BEFORE fitted()
log_progress("  Filtering CHM spatial sample to model levels:")
chm_s_v <- filter_to_model_levels(chm_s, fit_chm, "CHM (figS20)")
log_progress("  Filtering DTM spatial sample to model levels:")
dtm_s_v <- filter_to_model_levels(dtm_s, fit_dtm, "DTM (figS20)")

log_progress(sprintf("  CHM after filter: %s | DTM after filter: %s",
                     format(nrow(chm_s_v), big.mark = ","),
                     format(nrow(dtm_s_v), big.mark = ",")))

log_progress("  Computing CHM fitted() for residuals...")
pred_chm <- fitted(fit_chm, newdata = chm_s_v, summary = TRUE, ndraws = 100)
chm_s_v$residual <- chm_s_v$chm_error_mean - pred_chm[, "Estimate"]

log_progress("  Computing DTM fitted() for residuals...")
pred_dtm <- fitted(fit_dtm, newdata = dtm_s_v, summary = TRUE, ndraws = 100)
dtm_s_v$residual <- dtm_s_v$dtm_error_mean - pred_dtm[, "Estimate"]

# Variograms (50 km cutoff, 20 bins)
vg_calc_projected <- function(df_with_coords, res_col,
                              max_dist = MAX_LAG_DISTANCE,
                              n_lags = N_LAGS,
                              min_lag = MIN_LAG_DISTANCE) {
  df <- df_with_coords %>%
    mutate(res = .data[[res_col]]) %>%
    filter(is.finite(res))
  sf_obj <- st_as_sf(df, coords = c("x_proj", "y_proj"), crs = TARGET_CRS)
  g <- gstat::gstat(id = "res", formula = res ~ 1, data = sf_obj)
  lag_width <- (max_dist - min_lag) / n_lags
  gstat::variogram(g, cutoff = max_dist, width = lag_width)
}

log_progress("  Computing empirical variograms...")
vg_chm <- vg_calc_projected(chm_s_v, "residual")
vg_dtm <- vg_calc_projected(dtm_s_v, "residual")

fit_variogram_models <- function(vg_empirical, name) {
  if (is.null(vg_empirical)) return(NULL)
  results <- list()
  max_gamma <- max(vg_empirical$gamma)
  max_dist  <- max(vg_empirical$dist)
  init_nugget <- vg_empirical$gamma[1] * 0.8
  init_sill   <- max_gamma - init_nugget
  gamma_target <- init_nugget + 0.63 * init_sill
  init_range  <- vg_empirical$dist[which.min(abs(vg_empirical$gamma - gamma_target))]
  if (init_range < MIN_LAG_DISTANCE) init_range <- max_dist / 4
  for (m in c("Exp", "Sph", "Gau", "Mat")) {
    tryCatch({
      vg_model <- if (m == "Mat") {
        gstat::vgm(psill = init_sill, model = m, range = init_range,
                   nugget = init_nugget, kappa = 1.5)
      } else {
        gstat::vgm(psill = init_sill, model = m, range = init_range,
                   nugget = init_nugget)
      }
      fit <- gstat::fit.variogram(vg_empirical, vg_model)
      predicted <- gstat::variogramLine(fit, maxdist = max_dist,
                                        n = nrow(vg_empirical))
      pred_at_emp <- approx(predicted$dist, predicted$gamma,
                            vg_empirical$dist)$y
      weights <- vg_empirical$np / sum(vg_empirical$np)
      sse <- sum(weights * (vg_empirical$gamma - pred_at_emp)^2, na.rm = TRUE)
      nugget       <- fit$psill[1]
      total_sill   <- sum(fit$psill)
      range_param  <- fit$range[2]
      results[[m]] <- list(
        model = fit, sse = sse, nugget = nugget,
        sill = total_sill, range = range_param,
        nugget_sill_ratio = nugget / total_sill
      )
      log_progress(sprintf("    %s %s: nugget=%.2f, sill=%.2f, range=%.0f m, SSE=%.4f",
                           name, m, nugget, total_sill, range_param, sse))
    }, error = function(e) {
      log_progress(sprintf("    %s %s: fit failed (%s)", name, m, e$message))
    })
  }
  if (length(results) == 0) return(NULL)
  sse_vals <- sapply(results, function(x) x$sse)
  best <- names(which.min(sse_vals))
  results$best     <- best
  results$best_fit <- results[[best]]
  log_progress(sprintf("    %s best model: %s", name, best))
  results
}

vg_fits_chm <- fit_variogram_models(vg_chm, "CHM")
vg_fits_dtm <- fit_variogram_models(vg_dtm, "DTM")

params_df <- data.table(
  product           = c("CHM", "DTM"),
  best_model        = c(vg_fits_chm$best, vg_fits_dtm$best),
  nugget            = c(vg_fits_chm$best_fit$nugget,
                        vg_fits_dtm$best_fit$nugget),
  sill              = c(vg_fits_chm$best_fit$sill,
                        vg_fits_dtm$best_fit$sill),
  range_m           = c(vg_fits_chm$best_fit$range,
                        vg_fits_dtm$best_fit$range),
  nugget_sill_ratio = c(vg_fits_chm$best_fit$nugget_sill_ratio,
                        vg_fits_dtm$best_fit$nugget_sill_ratio)
)
fwrite(params_df,
       file.path(manuscript_tables_dir,
                 "section_L5_figS20_variogram_params_v2.csv"))

build_variogram_panel <- function(vg_data, vg_fit, title_prefix) {
  vg_plot_data <- vg_data
  fitted_line <- gstat::variogramLine(vg_fit$best_fit$model,
                                       maxdist = max(vg_data$dist), n = 100)
  ggplot() +
    geom_point(data = vg_plot_data,
               aes(x = dist / 1000, y = gamma, size = np), alpha = 0.7) +
    geom_line(data = fitted_line,
              aes(x = dist / 1000, y = gamma),
              color = "red", linewidth = 1) +
    geom_hline(yintercept = vg_fit$best_fit$sill,
               linetype = "dashed", color = "blue", alpha = 0.7) +
    geom_vline(xintercept = vg_fit$best_fit$range / 1000,
               linetype = "dashed", color = "forestgreen", alpha = 0.7) +
    geom_vline(xintercept = GEDI_SWATH_WIDTH / 1000,
               linetype = "dotted", color = "orange", alpha = 0.5) +
    annotate("text", x = GEDI_SWATH_WIDTH / 1000 + 1,
             y = max(vg_plot_data$gamma) * 0.95,
             label = "GEDI swath", color = "orange",
             hjust = 0, size = 3) +
    labs(
      title    = sprintf("%s: Semi-variogram with Fitted Model", title_prefix),
      subtitle = sprintf("%s model: nugget=%.2f m\u00b2, sill=%.2f m\u00b2, range=%.1f km\nNugget/Sill ratio=%.3f",
                         vg_fit$best, vg_fit$best_fit$nugget,
                         vg_fit$best_fit$sill, vg_fit$best_fit$range / 1000,
                         vg_fit$best_fit$nugget_sill_ratio),
      x = "Distance (km)", y = "Semivariance (m\u00b2)", size = "N pairs"
    ) +
    theme_cowplot() +
    theme(plot.subtitle = element_text(size = 9))
}

p_vg_chm <- build_variogram_panel(vg_chm, vg_fits_chm, "CHM")
p_vg_dtm <- build_variogram_panel(vg_dtm, vg_fits_dtm, "DTM")

p_s20 <- plot_grid(p_vg_chm, p_vg_dtm, ncol = 2,
                   labels = c("(a)", "(b)"),
                   label_size = 12, label_fontface = "bold")

ggsave(file.path(out_dir, "figS20_variogram_2panel_16site_v2.pdf"),
       p_s20, width = 16, height = 7, bg = "white", device = pdf)
log_progress(sprintf("  Saved: %s/figS20_variogram_2panel_16site_v2.pdf",
                     out_dir))

# ---------------------------------------------------------------------
# Done
# ---------------------------------------------------------------------

log_progress("")
log_progress(strrep("=", 70))
log_progress("Supplementary figure fixes complete. New outputs:")
log_progress(sprintf("  %s/figS06_pred_vs_obs_by_forest_16site_v2.pdf", out_dir))
log_progress(sprintf("  %s/figS18_ppc_overall_16site_v2.pdf",           out_dir))
log_progress(sprintf("  %s/figS20_variogram_2panel_16site_v2.pdf",      out_dir))
log_progress("")
log_progress("(figS12 from v1 run unchanged - keep that PDF.)")
log_progress(strrep("=", 70))
