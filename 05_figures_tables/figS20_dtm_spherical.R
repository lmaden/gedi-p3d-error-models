#!/usr/bin/env Rscript
# =====================================================================
# figS20_dtm_spherical.R
#
# Re-runs figS20 only, with DTM model selection restricted to Spherical.
# CHM keeps the full Exp/Sph/Gau/Mat selection (which picked Sph anyway,
# matching v04 — so this preserves the figS20 v2 CHM panel unchanged
# up to seed-level Monte Carlo noise in fitted() draws).
#
# Rationale: DTM model and data are unchanged from v04. The empirical
# variogram differs slightly because the spatial sample inherits the
# updated 01_data_ingest. Letting model selection switch from Sph to
# Gau on a 0.5% SSE margin produces a misleadingly different reported
# range (1.0 km vs ~4 km) for an unchanged underlying model. Forcing
# Sph for DTM preserves visual and numerical comparability with v04.
#
# RUN: from /gpfs/data1/vclgp/lmaden/chpt1
#   Rscript figS20_dtm_spherical.R
#
# Date: 2026-05-14
# =====================================================================

Sys.setenv(DISPLAY = "")
options(device = pdf)

suppressPackageStartupMessages({
  library(brms)
  library(ggplot2)
  library(cowplot)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(sf)
  library(gstat)
})

if (!exists("log_progress")) source("analysis_utils.R")
if (!exists("PROJECT_ROOT")) source("analysis_config.R")

out_dir <- file.path(PROJECT_ROOT, "plots", "groundwork")
manuscript_tables_dir <- file.path(PROJECT_ROOT, "manuscript_tables")

log_progress(strrep("=", 70))
log_progress("figS20: force the spherical model for DTM")
log_progress(strrep("=", 70))

# ---------------------------------------------------------------------
# Load checkpoints
# ---------------------------------------------------------------------

cp_chm <- load_checkpoint("10b_chm_sensitivity")
fit_chm <- cp_chm$fit_chm_16
sites_excluded_chm <- cp_chm$sites_excluded
log_progress(sprintf("  CHM 16-site refit loaded; excludes: %s",
                     paste(sites_excluded_chm, collapse = ", ")))

cp_s2 <- load_checkpoint("10_models_stage2")
fit_dtm <- cp_s2$fit_dtm_s2
log_progress("  DTM 18-site model loaded")

data_chk <- load_checkpoint("01_data_ingest")
chm_df_full_raw <- data_chk$chm_df
dtm_df_full_raw <- data_chk$dtm_df

# ---------------------------------------------------------------------
# Constants (match patch_v1)
# ---------------------------------------------------------------------

TARGET_CRS       <- 5070
SOURCE_CRS       <- 4326
GEDI_SWATH_WIDTH <- 6000
MIN_LAG_DISTANCE <- 100
MAX_LAG_DISTANCE <- 50000
N_LAGS           <- 20

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
        log_progress(sprintf("    %s: dropping %s rows with unknown %s: %s",
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

# ---------------------------------------------------------------------
# Spatial sample + residual extraction (same as patch_v1)
# ---------------------------------------------------------------------

chm_df_all <- chm_df_full_raw %>%
  filter(!as.character(site) %in% sites_excluded_chm)
dtm_df_all <- dtm_df_full_raw

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
log_progress(sprintf("  CHM spatial sample: %s | DTM: %s",
                     format(nrow(chm_s), big.mark = ","),
                     format(nrow(dtm_s), big.mark = ",")))

chm_sf <- st_transform(
  st_as_sf(chm_s, coords = c("x", "y"), crs = SOURCE_CRS, remove = FALSE),
  TARGET_CRS)
dtm_sf <- st_transform(
  st_as_sf(dtm_s, coords = c("x", "y"), crs = SOURCE_CRS, remove = FALSE),
  TARGET_CRS)

chm_s <- chm_s %>%
  mutate(x_proj = st_coordinates(chm_sf)[, 1],
         y_proj = st_coordinates(chm_sf)[, 2],
         lc_l1_code = factor(lc_l1_code),
         site       = factor(site),
         ecoregion  = factor(ecoregion))
dtm_s <- dtm_s %>%
  mutate(x_proj = st_coordinates(dtm_sf)[, 1],
         y_proj = st_coordinates(dtm_sf)[, 2],
         lc_l1_code = factor(lc_l1_code),
         site       = factor(site),
         ecoregion  = factor(ecoregion))

log_progress("  Filtering to model factor levels...")
chm_s_v <- filter_to_model_levels(chm_s, fit_chm, "CHM")
dtm_s_v <- filter_to_model_levels(dtm_s, fit_dtm, "DTM")

log_progress("  Computing CHM fitted() residuals...")
pred_chm <- fitted(fit_chm, newdata = chm_s_v, summary = TRUE, ndraws = 100)
chm_s_v$residual <- chm_s_v$chm_error_mean - pred_chm[, "Estimate"]

log_progress("  Computing DTM fitted() residuals...")
pred_dtm <- fitted(fit_dtm, newdata = dtm_s_v, summary = TRUE, ndraws = 100)
dtm_s_v$residual <- dtm_s_v$dtm_error_mean - pred_dtm[, "Estimate"]

# ---------------------------------------------------------------------
# Empirical variograms
# ---------------------------------------------------------------------

vg_calc_projected <- function(df, res_col) {
  df2 <- df %>% mutate(res = .data[[res_col]]) %>% filter(is.finite(res))
  sf_obj <- st_as_sf(df2, coords = c("x_proj", "y_proj"), crs = TARGET_CRS)
  g <- gstat::gstat(id = "res", formula = res ~ 1, data = sf_obj)
  lag_width <- (MAX_LAG_DISTANCE - MIN_LAG_DISTANCE) / N_LAGS
  gstat::variogram(g, cutoff = MAX_LAG_DISTANCE, width = lag_width)
}

vg_chm <- vg_calc_projected(chm_s_v, "residual")
vg_dtm <- vg_calc_projected(dtm_s_v, "residual")

# ---------------------------------------------------------------------
# Model fitting:
#   CHM: full Exp/Sph/Gau/Mat selection (preserves patch_v1 behavior)
#   DTM: Sph only (THIS IS THE PATCH)
# ---------------------------------------------------------------------

fit_one_model <- function(vg_empirical, m, name) {
  max_gamma <- max(vg_empirical$gamma)
  max_dist  <- max(vg_empirical$dist)
  init_nugget <- vg_empirical$gamma[1] * 0.8
  init_sill   <- max_gamma - init_nugget
  gamma_target <- init_nugget + 0.63 * init_sill
  init_range  <- vg_empirical$dist[which.min(abs(vg_empirical$gamma - gamma_target))]
  if (init_range < MIN_LAG_DISTANCE) init_range <- max_dist / 4

  vg_model <- if (m == "Mat") {
    gstat::vgm(psill = init_sill, model = m, range = init_range,
               nugget = init_nugget, kappa = 1.5)
  } else {
    gstat::vgm(psill = init_sill, model = m, range = init_range,
               nugget = init_nugget)
  }
  fit <- gstat::fit.variogram(vg_empirical, vg_model)
  predicted <- gstat::variogramLine(fit, maxdist = max_dist, n = nrow(vg_empirical))
  pred_at_emp <- approx(predicted$dist, predicted$gamma, vg_empirical$dist)$y
  weights <- vg_empirical$np / sum(vg_empirical$np)
  sse <- sum(weights * (vg_empirical$gamma - pred_at_emp)^2, na.rm = TRUE)
  nugget <- fit$psill[1]; sill <- sum(fit$psill); rng <- fit$range[2]
  log_progress(sprintf("    %s %s: nugget=%.2f, sill=%.2f, range=%.0f m, SSE=%.4f",
                       name, m, nugget, sill, rng, sse))
  list(model = fit, sse = sse, nugget = nugget,
       sill = sill, range = rng, nugget_sill_ratio = nugget / sill)
}

# CHM: full selection
log_progress("  Fitting CHM variogram models (full selection):")
chm_results <- list()
for (m in c("Exp", "Sph", "Gau", "Mat")) {
  tryCatch({
    chm_results[[m]] <- fit_one_model(vg_chm, m, "CHM")
  }, error = function(e) log_progress(sprintf("    CHM %s: failed (%s)", m, e$message)))
}
chm_best <- names(which.min(sapply(chm_results, function(x) x$sse)))
vg_fits_chm <- list(best = chm_best, best_fit = chm_results[[chm_best]])
log_progress(sprintf("    CHM best: %s", chm_best))

# DTM: SPH ONLY (the patch)
log_progress("  Fitting DTM variogram model (Sph forced for v04 continuity):")
dtm_sph <- fit_one_model(vg_dtm, "Sph", "DTM")
vg_fits_dtm <- list(best = "Sph", best_fit = dtm_sph)
log_progress("    DTM model: Sph (forced)")

# ---------------------------------------------------------------------
# Save params CSV
# ---------------------------------------------------------------------

params_df <- data.table(
  product           = c("CHM", "DTM"),
  best_model        = c(vg_fits_chm$best, vg_fits_dtm$best),
  selection_mode    = c("full",           "Sph-forced"),
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
                 "section_L5_figS20_variogram_params_v3_sph_dtm.csv"))

# ---------------------------------------------------------------------
# Build panels
# ---------------------------------------------------------------------

build_variogram_panel <- function(vg_data, vg_fit, title_prefix) {
  fitted_line <- gstat::variogramLine(vg_fit$best_fit$model,
                                       maxdist = max(vg_data$dist), n = 100)
  ggplot() +
    geom_point(data = vg_data,
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
             y = max(vg_data$gamma) * 0.95,
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

ggsave(file.path(out_dir, "figS20_variogram_2panel_16site_v3_sph.pdf"),
       p_s20, width = 16, height = 7, bg = "white", device = pdf)

log_progress("")
log_progress(strrep("=", 70))
log_progress("figS20 patch complete. Output:")
log_progress(sprintf("  %s/figS20_variogram_2panel_16site_v3_sph.pdf", out_dir))
log_progress(strrep("=", 70))
