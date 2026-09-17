# =====================================================================
# fig6_spatial_residuals.R
# Figure 6: spatial distribution of model residuals across four
# representative sites (CHM 16-site refit; DTM 18-site unchanged).
#
# v04 -> v10 changes:
#   - CHM residuals computed from fit_chm_16 (was 19-site fit_chm_s2);
#     numbers shift; spatial structure shouldn't change qualitatively
#     because all four sites (manuscript IDs 4, 12, 18, 19) are inside
#     the 16-site frame.
#   - DTM residuals unchanged.
#   - Site selection: manuscript IDs 4 (tracker 4), 12 (tracker 12),
#     18 (tracker 19, ltbmu_20180710), 19 (tracker 20, neon_sjer)
#     .
#   - In-figure titles use manuscript IDs.
#
# Sources:
#   - load_checkpoint("10b_chm_sensitivity")$fit_chm_16
#   - load_checkpoint("10_models_stage2")$fit_dtm_s2
#   - load_checkpoint("08_model_prep")$mod_chm_s2, mod_dtm_s2
# =====================================================================

source("fig_common.R")
fig_banner("Figure 6",
                 "Spatial residuals 4-site (manuscript 4, 12, 18, 19)")

suppressPackageStartupMessages({
  library(brms); library(ggplot2); library(scales)
  library(data.table); library(dplyr); library(patchwork); library(sf)
})

# ---------------------------------------------------------------------
# Site selection (manuscript IDs; map to tracker via SITE_LOOKUP)
# ---------------------------------------------------------------------
SITE_SELECTION_MANUSCRIPT <- c(4L, 12L, 18L, 19L)
SITE_SELECTION_TRACKER    <- manuscript_to_tracker(SITE_SELECTION_MANUSCRIPT)
log_progress(sprintf("  Manuscript IDs: %s -> tracker IDs: %s",
                     paste(SITE_SELECTION_MANUSCRIPT, collapse = ","),
                     paste(SITE_SELECTION_TRACKER,    collapse = ",")))

# ---------------------------------------------------------------------
# 1. Load fits and data
# ---------------------------------------------------------------------
log_subsection("Loading fits and data")

chm_ck <- load_checkpoint("10b_chm_sensitivity")
fit_chm_16 <- chm_ck$fit_chm_16

s2 <- load_checkpoint("10_models_stage2")
fit_dtm <- s2$fit_dtm_s2

mp <- load_checkpoint("08_model_prep")
mod_chm <- mp$mod_chm_s2
mod_dtm <- mp$mod_dtm_s2

# Filter to selection (tracker ID basis).
mod_chm_sub <- mod_chm[as.character(mod_chm$site) %in% SITE_SELECTION_TRACKER, ]
mod_dtm_sub <- mod_dtm[as.character(mod_dtm$site) %in% SITE_SELECTION_TRACKER, ]
log_progress(sprintf("  CHM rows: %d  DTM rows: %d",
                     nrow(mod_chm_sub), nrow(mod_dtm_sub)))

# ---------------------------------------------------------------------
# 2. Identify lon/lat columns. Footprints carry lon/lat (or x/y).
# ---------------------------------------------------------------------
ALBERS_CRS <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"

# Coordinate detection. Tier 1: canonical projected names (per
# section_12_spatial.R / publication_figures.R, which use
# x_proj / y_proj in meters). Tier 2: explicit lon/lat names. Tier 3:
# ambiguous x / y - validate by VALUE RANGE since the model-prep data
# uses 'x' and 'y' for lon/lat in degrees (NOT projected meters); we
# need to detect that and transform rather than naively divide by 1000.
detect_xy <- function(df) {
  # Tier 1: explicit projected names
  for (pair in list(c("x_proj", "y_proj"),
                    c("x_albers", "y_albers"),
                    c("easting", "northing"))) {
    if (all(pair %in% names(df))) {
      return(list(cols = pair, kind = "projected_m"))
    }
  }
  # Tier 2: explicit lon/lat names
  for (lon_c in c("lon", "longitude", "shot_lon", "x_lon")) {
    if (lon_c %in% names(df)) {
      for (lat_c in c("lat", "latitude", "shot_lat", "y_lat")) {
        if (lat_c %in% names(df)) {
          return(list(cols = c(lon_c, lat_c), kind = "lonlat"))
        }
      }
    }
  }
  # Tier 3: ambiguous x / y - detect by value range
  if ("x" %in% names(df) && "y" %in% names(df)) {
    x_range <- range(df$x, na.rm = TRUE)
    y_range <- range(df$y, na.rm = TRUE)
    if (max(abs(x_range)) <= 180 && max(abs(y_range)) <= 90) {
      return(list(cols = c("x", "y"), kind = "lonlat"))
    } else {
      return(list(cols = c("x", "y"), kind = "projected_m"))
    }
  }
  stop("Could not detect coordinate columns. Have: ",
       paste(names(df), collapse = ", "))
}

xy_chm <- detect_xy(mod_chm_sub)
xy_dtm <- detect_xy(mod_dtm_sub)
log_progress(sprintf("  CHM coords: %s, %s (kind: %s)",
                     xy_chm$cols[1], xy_chm$cols[2], xy_chm$kind))
log_progress(sprintf("  DTM coords: %s, %s (kind: %s)",
                     xy_dtm$cols[1], xy_dtm$cols[2], xy_dtm$kind))

# ---------------------------------------------------------------------
# 3. Compute residuals on holdout split.
# Residual = observed - posterior predictive mean.
# ---------------------------------------------------------------------
log_subsection("Computing residuals")

# Sample within each site to keep plotting tractable.
N_PER_SITE <- 5000L
set.seed(42L)

subsample_per_site <- function(df) {
  out <- df[, {
    n <- .N
    if (n <= N_PER_SITE) .SD else .SD[sample.int(n, N_PER_SITE)]
  }, by = "site"]
  out
}

chm_sub <- subsample_per_site(as.data.table(mod_chm_sub))
dtm_sub <- subsample_per_site(as.data.table(mod_dtm_sub))

predict_mean <- function(fit, df, ndraws = 100L) {
  pm <- posterior_predict(fit, newdata = df, ndraws = ndraws,
                          allow_new_levels = TRUE, re_formula = NULL)
  colMeans(pm)
}

chm_sub[, y_pred := predict_mean(fit_chm_16, chm_sub)]
chm_sub[, resid  := y_pred - chm_error_mean]
dtm_sub[, y_pred := predict_mean(fit_dtm,    dtm_sub)]
dtm_sub[, resid  := y_pred - dtm_error_mean]

# Map tracker -> manuscript for axis titles.
chm_sub[, manuscript_id := tracker_to_manuscript(site)]
dtm_sub[, manuscript_id := tracker_to_manuscript(site)]

# Use detected coordinate columns. Convert lon/lat to Albers if so.
to_albers <- function(df, xy_info) {
  if (xy_info$kind == "lonlat") {
    sf_pts <- st_as_sf(df, coords = xy_info$cols, crs = 4326, remove = FALSE)
    sf_pts <- st_transform(sf_pts, ALBERS_CRS)
    coords <- st_coordinates(sf_pts)
    df[, x_km := coords[, 1] / 1000]
    df[, y_km := coords[, 2] / 1000]
  } else if (xy_info$kind == "projected_m") {
    df[, x_km := df[[xy_info$cols[1]]] / 1000]
    df[, y_km := df[[xy_info$cols[2]]] / 1000]
  } else {
    stop("Unknown coord kind: ", xy_info$kind)
  }
  df
}

chm_sub <- to_albers(chm_sub, xy_chm)
dtm_sub <- to_albers(dtm_sub, xy_dtm)

# ---------------------------------------------------------------------
# 4. Build panels (CHM/DTM x 4 sites = 8 panels) with shared color scale
# ---------------------------------------------------------------------
log_subsection("Composing panels")

# Cap color scale at +-10 m (matches v04).
LIM <- 10
chm_sub[, resid_c := pmin(pmax(resid, -LIM), LIM)]
dtm_sub[, resid_c := pmin(pmax(resid, -LIM), LIM)]

build_resid_panel <- function(df, manu_id, product, prod_label) {
  d <- df[manuscript_id == manu_id]
  if (!nrow(d)) return(plot_spacer())
  ggplot(d, aes(x = x_km, y = y_km, color = resid_c)) +
    geom_point(size = 0.4, alpha = 0.9) +
    scale_color_distiller(palette = "RdBu", direction = -1,
                          limits = c(-LIM, LIM),
                          oob = scales::squish,
                          name = "Residual (m)") +
    coord_equal() +
    labs(
      title = sprintf("%s %s Site %d", prod_label, MINUS, manu_id),
      x = "Easting (km)", y = "Northing (km)"
    ) +
    theme_section_G(base_size = 9) +
    theme(plot.title = element_text(size = 9, face = "bold"),
          legend.position = "right",
          axis.title  = element_text(size = 8),
          axis.text   = element_text(size = 7))
}

build_resid_panel_clean <- function(...) {
  p <- build_resid_panel(...)
  p + theme(legend.position = "none")
}

chm_panels <- lapply(SITE_SELECTION_MANUSCRIPT,
                     function(m) build_resid_panel_clean(chm_sub, m, "CHM", "CHM"))
dtm_panels <- lapply(SITE_SELECTION_MANUSCRIPT,
                     function(m) build_resid_panel_clean(dtm_sub, m, "DTM", "DTM"))

# Use one panel with legend (rightmost CHM) to extract a shared legend.
legend_p <- build_resid_panel(chm_sub, SITE_SELECTION_MANUSCRIPT[1], "CHM", "CHM")
legend_obj <- cowplot::get_legend(legend_p +
                                    theme(legend.position = "right"))

grid_top <- patchwork::wrap_plots(chm_panels, nrow = 1)
grid_bot <- patchwork::wrap_plots(dtm_panels, nrow = 1)
grid_all <- (grid_top / grid_bot)

final_plot <- cowplot::plot_grid(
  cowplot::plot_grid(grid_all, NULL, ncol = 1, rel_heights = c(1, 0.02)),
  legend_obj,
  nrow = 1, rel_widths = c(1, 0.10)
)

# ---------------------------------------------------------------------
# 5. Save and verify
# ---------------------------------------------------------------------
save_figure(final_plot, "fig06_spatial_residuals",
            width_in = 11.0, height_in = 5.4)

log_subsection("VERIFICATION (console anchors)")
cat("\nSite selection :\n")
for (i in seq_along(SITE_SELECTION_MANUSCRIPT)) {
  m <- SITE_SELECTION_MANUSCRIPT[i]
  tr <- SITE_SELECTION_TRACKER[i]
  ck <- chm_sub[manuscript_id == m]
  dk <- dtm_sub[manuscript_id == m]
  cat(sprintf("  Manuscript %d (tracker %d):  CHM n=%s  DTM n=%s\n",
              m, tr,
              formatC(nrow(ck), big.mark = ",", format = "d"),
              formatC(nrow(dk), big.mark = ",", format = "d")))
}
cat("\nResidual = predicted err - observed err (m); convention matches\n")
cat("canonical verify_predictive_accuracy.R (Bias = mean(pred - obs)).\n")
cat("Color scale capped at +/- 10 m to preserve mid-range visibility.\n")
cat("Red regions: model over-predicts observed error (predicted more positive).\n")
cat("Blue regions: model under-predicts observed error (predicted more negative).\n")

log_progress("fig6_spatial_residuals.R complete.")
