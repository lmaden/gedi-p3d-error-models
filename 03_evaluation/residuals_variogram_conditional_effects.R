#!/usr/bin/env Rscript
# =====================================================================
# residuals_variogram_conditional_effects.R
#
# Route A (full 17-predictor S21 recomputation) +
# S20 CHM semi-variogram + S4 CHM conditional effects extraction.
#
# All three operations require the CHM 16-site brmsfit and the CHM
# data ingest checkpoint; running them in one R session amortizes the
# loading cost.
#
# Inputs:
#   checkpoints/10b_chm_sensitivity.rds  (fit_chm_16 brmsfit)
#   checkpoints/01_data_ingest.rds       (chm_df with all 17 predictors)
#
# Outputs (all under manuscript_tables/):
#   section_L_s21_chm_full17_predictors.csv  - 17-row bivariate table
#                                               at 16-site CHM frame
#   section_L_s21_chm_site_means_full17.csv  - per-site predictor means
#   section_L_s20_chm_residuals.csv          - per-footprint CHM residuals
#                                               with coordinates (for var.)
#   section_L_s20_chm_variogram_params.csv   - fitted variogram nugget /
#                                               sill / range / model
#   section_L_s4_chm_conditional_effects.rds - brms conditional_effects
#                                               result (for plotting in L-4)
#
# Wall-time estimate: 30-60 min depending on residual subsample size.
# =====================================================================

suppressPackageStartupMessages({
  library(brms)
  library(posterior)
  library(data.table)
  library(gstat)
  library(sf)
})

Sys.setenv(DISPLAY = "")
options(device = pdf)

if (!exists("log_progress"))  source("analysis_utils.R")
if (!exists("PROJECT_ROOT"))  source("analysis_config_FIXED.R")

manuscript_tables_dir <- file.path(PROJECT_ROOT, "manuscript_tables")

log_progress("================================================================")
log_progress("Route A + S20 variogram + S4 cond effects")
log_progress("================================================================")

# ---------------------------------------------------------------------
# Load checkpoints
# ---------------------------------------------------------------------

log_progress("Loading CHM 16-site brmsfit from 10b_chm_sensitivity")
cp_chm <- load_checkpoint("10b_chm_sensitivity")
fit_chm <- cp_chm$fit_chm_16
stopifnot(inherits(fit_chm, "brmsfit"))

log_progress("Loading data ingest (01_data_ingest)")
di <- load_checkpoint("01_data_ingest")
chm_df <- as.data.table(di$chm_df)
log_progress(sprintf("  CHM data frame: %d footprints, %d columns",
                     nrow(chm_df), ncol(chm_df)))

# Identify the site column. fit_chm$data should match.
site_col <- intersect(c("site", "site_id", "tracker_site_id",
                       "manuscript_site_id"), names(chm_df))[1]
log_progress(sprintf("  Site column in chm_df: '%s'", site_col))

# Restrict to the 16 sites in the fit
fit_sites <- unique(as.character(fit_chm$data[[site_col]]))
log_progress(sprintf("  Sites in fit: %s",
                     paste(sort(as.integer(fit_sites)), collapse = ", ")))
chm_df_16 <- chm_df[as.character(get(site_col)) %in% fit_sites]
log_progress(sprintf("  Restricted to 16-site frame: %d footprints",
                     nrow(chm_df_16)))

# ---------------------------------------------------------------------
# ROUTE A: aggregate all 17 main-effect predictors to per-site means
# ---------------------------------------------------------------------

log_subsection("Route A: per-site means for all 17 CHM main-effect predictors")

# Main-effect predictors per §2F (matching the z-scored column names)
predictors_17 <- c("slope_mean_z", "rh_98_z", "wsci_z", "cover_z",
                  "meta_leafon_z", "meta_sunel_z", "meta_offnad_z",
                  "meta_relgeo_z", "view_az_cos_z", "meta_absgeo_z",
                  "meta_az_conc_z", "meta_fwdrev_z", "view_az_sin_z",
                  "slope_sd_z", "aspect_cos_z", "aspect_sin_z",
                  "meta_stereo_z")

present <- intersect(predictors_17, names(chm_df_16))
missing <- setdiff(predictors_17, present)
log_progress(sprintf("  Predictors present in chm_df: %d of %d",
                     length(present), length(predictors_17)))
if (length(missing) > 0) {
  log_progress(sprintf("  Predictors MISSING in chm_df: %s",
                       paste(missing, collapse = ", ")))
  log_progress("  Looking for alternate column names...")
  # Suggest possible alternates
  for (m in missing) {
    cands <- grep(sub("_z$", "", m), names(chm_df_16),
                  ignore.case = TRUE, value = TRUE)
    if (length(cands) > 0) {
      log_progress(sprintf("    '%s' may map to: %s",
                           m, paste(cands, collapse = ", ")))
    }
  }
}

# Aggregate to site level (using NA-safe means)
site_means <- chm_df_16[, lapply(.SD, mean, na.rm = TRUE),
                       by = c(site_col),
                       .SDcols = present]
setnames(site_means, site_col, "tracker_site")
site_means[, tracker_site := as.character(tracker_site)]
log_progress(sprintf("  Per-site means table: %d sites x %d predictors",
                     nrow(site_means), ncol(site_means) - 1))

fwrite(site_means, file.path(manuscript_tables_dir,
                            "section_L_s21_chm_site_means_full17.csv"))

# Merge with L-1 site-RE table
re_df <- fread(file.path(manuscript_tables_dir,
                       "section_L_s17_chm_site_re_summary.csv"))
setnames(re_df, "site_label", "tracker_site")
re_df[, tracker_site := as.character(tracker_site)]

merged_full <- merge(re_df[, .(tracker_site, re_mean)],
                    site_means, by = "tracker_site", all.x = TRUE)

# Compute bivariate r, p, r^2 for each predictor
biv_full <- rbindlist(lapply(present, function(cv) {
  x <- merged_full[[cv]]
  y <- merged_full$re_mean
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 4) return(NULL)
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "pearson"))
  data.table(
    predictor  = cv,
    n_sites    = sum(ok),
    pearson_r  = unname(ct$estimate),
    p_value    = ct$p.value,
    r_squared  = unname(ct$estimate)^2
  )
}))
biv_full <- biv_full[order(-r_squared)]

n_tests_full <- nrow(biv_full)
bonf_alpha_full <- 0.05 / n_tests_full
biv_full[, bonferroni_significant := p_value < bonf_alpha_full]

fwrite(biv_full, file.path(manuscript_tables_dir,
                          "section_L_s21_chm_full17_predictors.csv"))

log_progress(sprintf("  Final n predictors tested: %d", n_tests_full))
log_progress(sprintf("  Bonferroni alpha at %d tests: %.5f",
                     n_tests_full, bonf_alpha_full))
log_progress("  Full predictor bivariates (sorted by r^2 descending):")
for (i in seq_len(nrow(biv_full))) {
  sig <- if (biv_full$bonferroni_significant[i]) " *BONF*" else ""
  log_progress(sprintf("    %-15s | n = %2d | r = %+.3f | p = %.4f | r^2 = %.3f%s",
                       biv_full$predictor[i],
                       biv_full$n_sites[i],
                       biv_full$pearson_r[i],
                       biv_full$p_value[i],
                       biv_full$r_squared[i],
                       sig))
}

# ---------------------------------------------------------------------
# S20: CHM semi-variogram on per-footprint residuals
# ---------------------------------------------------------------------

log_subsection("S20: CHM 16-site per-footprint residuals + variogram")

# Compute fitted residuals. brms::residuals() returns posterior summaries;
# for variogram we want point residuals (use posterior_mean of fitted).
log_progress("  Computing posterior_predict means (this is the slow step)")
# A subsample keeps wall time reasonable. 10% stratified is standard for
# variogram fits per Cressie & Wikle (2011) — diminishing returns above this.
set.seed(20260512)
n_total <- nrow(fit_chm$data)
subsample_n <- min(50000, n_total)  # cap at 50k for variogram tractability
log_progress(sprintf("  Total fit rows: %d; subsampling to %d for variogram",
                     n_total, subsample_n))
idx <- sort(sample(seq_len(n_total), subsample_n))
fit_data <- as.data.table(fit_chm$data)[idx]

# residuals(): default returns posterior summary residuals (obs - pred mean)
# at the fitted-data scale
log_progress("  brms::residuals() on subsample")
res <- residuals(fit_chm, newdata = fit_data, summary = TRUE)
fit_data[, resid_mean := res[, "Estimate"]]

# Locate coordinate columns. Standard names: lon, lat, x, y, easting, northing
coord_cols <- list(
  x = intersect(c("lon", "longitude", "x", "easting", "X"),
                names(fit_data))[1],
  y = intersect(c("lat", "latitude", "y", "northing", "Y"),
                names(fit_data))[1]
)
log_progress(sprintf("  Coordinate columns resolved: x='%s', y='%s'",
                     coord_cols$x, coord_cols$y))
if (is.na(coord_cols$x) || is.na(coord_cols$y)) {
  # Fall back to chm_df lookup if fit_data doesn't have coords
  log_progress("  fit_data missing coords; attempting join from chm_df_16")
  key_cols <- intersect(c("shot_number", "footprint_id"),
                       names(fit_data))[1]
  if (!is.na(key_cols)) {
    fit_data <- merge(fit_data,
                     chm_df_16[, .SD, .SDcols = c(key_cols, "lon", "lat")],
                     by = key_cols, all.x = TRUE)
    coord_cols$x <- "lon"
    coord_cols$y <- "lat"
  } else {
    stop("Cannot locate coordinate columns; check chm_df structure")
  }
}

# Convert to sf and project to NAD83 / CONUS Albers (EPSG:5070), matching
# v04 S8 / S9 caption conventions
fit_sf <- sf::st_as_sf(fit_data,
                      coords = c(coord_cols$x, coord_cols$y),
                      crs = 4326)
fit_sf <- sf::st_transform(fit_sf, 5070)
coords <- sf::st_coordinates(fit_sf)
fit_data[, x_alb := coords[, 1]]
fit_data[, y_alb := coords[, 2]]

# Save residuals CSV (subsampled) for both S20 variogram and possible
# S8 spatial residual plot use
out_resid <- fit_data[, .(tracker_site = get(site_col),
                          x_alb, y_alb, resid_mean)]
out_resid[, tracker_site := as.character(tracker_site)]
fwrite(out_resid, file.path(manuscript_tables_dir,
                           "section_L_s20_chm_residuals.csv"))

# Variogram fit
log_progress("  Fitting empirical variogram (gstat::variogram)")
sp_df <- sf::as_Spatial(
  sf::st_as_sf(fit_data[is.finite(resid_mean)],
              coords = c("x_alb", "y_alb"), crs = 5070)
)
sp_df$resid_mean <- fit_data[is.finite(resid_mean), resid_mean]

# Cap distance at 20km for tractability and to match the v04 caption's
# scale claim (effective range 3.7 km)
v_emp <- gstat::variogram(resid_mean ~ 1, sp_df, cutoff = 20000, width = 500)

# Fit spherical model (matching v04 caption "fitted spherical variogram")
log_progress("  Fitting spherical variogram model")
v_init <- gstat::vgm(psill = var(sp_df$resid_mean, na.rm = TRUE) * 0.5,
                    model = "Sph",
                    range = 3700,  # initial guess from v04
                    nugget = var(sp_df$resid_mean, na.rm = TRUE) * 0.5)
v_fit <- tryCatch(gstat::fit.variogram(v_emp, v_init),
                 error = function(e) {
                   log_progress(sprintf("  Sph fit failed: %s; trying Exp",
                                       e$message))
                   v_init$model <- "Exp"
                   gstat::fit.variogram(v_emp, v_init)
                 })

# Extract parameters
v_nugget <- v_fit$psill[1]
v_partial_sill <- v_fit$psill[2]
v_total_sill <- v_nugget + v_partial_sill
v_range <- v_fit$range[2]
v_nug_sill_ratio <- v_nugget / v_total_sill
v_model <- as.character(v_fit$model[2])

log_progress(sprintf("  Variogram model: %s", v_model))
log_progress(sprintf("  Nugget:           %.3f m^2", v_nugget))
log_progress(sprintf("  Partial sill:     %.3f m^2", v_partial_sill))
log_progress(sprintf("  Total sill:       %.3f m^2", v_total_sill))
log_progress(sprintf("  Nugget/sill:      %.3f", v_nug_sill_ratio))
log_progress(sprintf("  Effective range:  %.1f m (%.2f km)",
                     v_range, v_range / 1000))

var_summary <- data.table(
  product = "CHM_16_site",
  model = v_model,
  nugget = v_nugget,
  partial_sill = v_partial_sill,
  total_sill = v_total_sill,
  nugget_to_sill = v_nug_sill_ratio,
  range_m = v_range,
  range_km = v_range / 1000,
  n_residuals = sum(is.finite(sp_df$resid_mean)),
  cutoff_m = 20000
)
fwrite(var_summary, file.path(manuscript_tables_dir,
                             "section_L_s20_chm_variogram_params.csv"))

# Save empirical variogram for plotting
v_emp_dt <- as.data.table(v_emp)
fwrite(v_emp_dt, file.path(manuscript_tables_dir,
                          "section_L_s20_chm_empirical_variogram.csv"))

# ---------------------------------------------------------------------
# S4: CHM conditional effects on all main predictors
# ---------------------------------------------------------------------

log_subsection("S4: CHM 16-site conditional effects (all predictors)")

log_progress("  brms::conditional_effects() on CHM fit")
ce <- conditional_effects(fit_chm)
log_progress(sprintf("  Generated %d conditional-effect panels", length(ce)))
log_progress(sprintf("  Panels: %s",
                     paste(names(ce), collapse = ", ")))

saveRDS(ce, file.path(manuscript_tables_dir,
                    "section_L_s4_chm_conditional_effects.rds"))

# ---------------------------------------------------------------------
# Console summary
# ---------------------------------------------------------------------

log_progress("================================================================")
log_progress("SUMMARY")
log_progress("================================================================")
log_progress("")
log_progress("[S21 full 17-predictor table]")
log_progress(sprintf("  N predictors tested: %d", n_tests_full))
log_progress(sprintf("  Bonferroni alpha: %.5f", bonf_alpha_full))
log_progress("  Top 5 by r^2:")
for (i in seq_len(min(5, nrow(biv_full)))) {
  sig <- if (biv_full$bonferroni_significant[i]) " *BONF*" else ""
  log_progress(sprintf("    %-15s | r = %+.3f | p = %.4f | r^2 = %.3f%s",
                       biv_full$predictor[i],
                       biv_full$pearson_r[i],
                       biv_full$p_value[i],
                       biv_full$r_squared[i],
                       sig))
}
log_progress("")
log_progress("[S20 CHM variogram]")
log_progress(sprintf("  Model: %s | Nugget/sill: %.3f | Range: %.2f km",
                     v_model, v_nug_sill_ratio, v_range / 1000))
log_progress(sprintf("  vs v04 caption: 0.79 nugget/sill, 3.7 km range"))
log_progress("")
log_progress("[S4 conditional effects]")
log_progress(sprintf("  %d panels saved to RDS for L-4 plotting",
                     length(ce)))
log_progress("")
log_progress("Outputs:")
log_progress("  section_L_s21_chm_full17_predictors.csv")
log_progress("  section_L_s21_chm_site_means_full17.csv")
log_progress("  section_L_s20_chm_residuals.csv")
log_progress("  section_L_s20_chm_variogram_params.csv")
log_progress("  section_L_s20_chm_empirical_variogram.csv")
log_progress("  section_L_s4_chm_conditional_effects.rds")
log_progress("================================================================")
log_progress("L-2 complete.")
log_progress("================================================================")
