#!/usr/bin/env Rscript
# =====================================================================
# holdout_04_dtm_panel_data.R
#
# Step 04: DTM-side data for the dual-panel figures S18, S20, S21.
#
# The DTM model is unchanged under 16-site primary (DTM analysis uses
# all 18 DTM-OK sites independent of CHM screening). This script:
#   - Reproduces v04 DTM-side numerics from the current fit_dtm
#     checkpoint (validates v04 caption claims against current data)
#   - Saves the data needed for the plotting step to render the DTM panels of
#     S18 (PPC), S20 (variogram), S21 (bivariates)
#
# v04 DTM-side caption claims to verify:
#   S18 panel b:  94.0% coverage              (step 03 already got 94.44%)
#   S20 right:    0.76 nugget/sill, 3.3 km range
#   S21 right:    canopy cover r^2 = 0.24
#
# Applies the same factor-level filter-then-recast pattern that worked
# for CHM in holdout_02_scoring.R.
#
# Verified column names (from holdout_diagnostic.R):
#   Response: dtm_error_mean
#   Coords:   x, y in WGS84 (EPSG:4326) -> reproject to EPSG:5070
#   Site:     'site' (18 sites for DTM analysis)
#
# Wall time estimate: 15-25 min.
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
if (!exists("PROJECT_ROOT"))  source("analysis_config.R")

manuscript_tables_dir <- file.path(PROJECT_ROOT, "manuscript_tables")

# Verified DTM column names
RESPONSE_COL  <- "dtm_error_mean"
COORD_X       <- "x"
COORD_Y       <- "y"
SITE_COL      <- "site"
COORD_CRS_IN  <- 4326
COORD_CRS_OUT <- 5070

log_progress("================================================================")
log_progress("DTM-side dual-panel data (S18, S20, S21)")
log_progress("================================================================")

# ---------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------

log_progress("Loading DTM 18-site brmsfit")
cp_dtm <- load_checkpoint("10_models_stage2")
fit_dtm <- cp_dtm$fit_dtm_s2

log_progress("Loading data ingest")
di <- load_checkpoint("01_data_ingest")
dtm_df <- as.data.table(di$dtm_df)
log_progress(sprintf("  DTM data frame: %d footprints, %d columns",
                     nrow(dtm_df), ncol(dtm_df)))

# Sanity-check required columns are present in dtm_df
required <- c(RESPONSE_COL, COORD_X, COORD_Y, SITE_COL,
              "lc_l1_code", "ecoregion")
missing <- setdiff(required, names(dtm_df))
if (length(missing) > 0) {
  log_progress(sprintf("  MISSING required columns: %s",
                       paste(missing, collapse = ", ")))
  log_progress(sprintf("  Available columns: %s",
                       paste(sort(names(dtm_df)), collapse = ", ")))
  stop("Required columns missing.")
}
log_progress(sprintf("  Required cols present: %s",
                     paste(required, collapse = ", ")))

# Quick sanity stats on response and coords
log_progress(sprintf("  %s: n_finite = %d, range [%+.2f, %+.2f]",
                     RESPONSE_COL,
                     sum(is.finite(dtm_df[[RESPONSE_COL]])),
                     min(dtm_df[[RESPONSE_COL]], na.rm = TRUE),
                     max(dtm_df[[RESPONSE_COL]], na.rm = TRUE)))

# Restrict to fit's 18 sites
fit_sites <- unique(as.character(fit_dtm$data[[SITE_COL]]))
log_progress(sprintf("  Fit DTM sites (%d): %s", length(fit_sites),
                     paste(sort(fit_sites), collapse = ", ")))
dtm_df_18 <- dtm_df[as.character(get(SITE_COL)) %in% fit_sites]
log_progress(sprintf("  18-site DTM frame: %d footprints", nrow(dtm_df_18)))

# Helper: get level set of a column
get_levels <- function(x) {
  if (is.factor(x)) levels(x) else sort(unique(as.character(x)))
}

lc_levels_fit   <- get_levels(fit_dtm$data$lc_l1_code)
ec_levels_fit   <- get_levels(fit_dtm$data$ecoregion)
site_levels_fit <- get_levels(fit_dtm$data$site)
log_progress(sprintf("  Fit lc_l1_code levels (%d): %s",
                     length(lc_levels_fit),
                     paste(sort(lc_levels_fit), collapse = ", ")))
log_progress(sprintf("  Fit ecoregion levels (%d)", length(ec_levels_fit)))

# Report extras (defensive)
lc_extra <- setdiff(get_levels(dtm_df_18$lc_l1_code), lc_levels_fit)
ec_extra <- setdiff(get_levels(dtm_df_18$ecoregion),  ec_levels_fit)
log_progress(sprintf("  lc levels in dtm_df_18 not in fit: %s",
                     if (length(lc_extra) == 0) "(none)"
                     else paste(lc_extra, collapse = ", ")))
log_progress(sprintf("  ecoregion levels in dtm_df_18 not in fit: %s",
                     if (length(ec_extra) == 0) "(none)"
                     else paste(ec_extra, collapse = ", ")))

# Filter to fit's factor levels (drops any rows with LMS-style new levels)
n_before <- nrow(dtm_df_18)
dtm_filtered <- dtm_df_18[
  as.character(lc_l1_code) %in% lc_levels_fit &
  as.character(ecoregion)  %in% ec_levels_fit
]
n_after <- nrow(dtm_filtered)
log_progress(sprintf("  After filtering to fit's factor levels: %d rows (%d dropped, %.2f%%)",
                     n_after, n_before - n_after,
                     100 * (n_before - n_after) / n_before))

# ---------------------------------------------------------------------
# (S21 DTM) Route A: site REs + per-site predictor bivariates
# FAST — do first
# ---------------------------------------------------------------------

log_subsection("S21 DTM: site REs + per-site predictor bivariates")

# Per-site REs
dtm_site_re <- ranef(fit_dtm, summary = TRUE,
                     probs = c(0.025, 0.5, 0.975))[[SITE_COL]]
dtm_re_df <- data.table(
  tracker_site = dimnames(dtm_site_re)[[1]],
  re_mean      = dtm_site_re[, "Estimate",  "Intercept"],
  re_sd        = dtm_site_re[, "Est.Error", "Intercept"],
  re_q025      = dtm_site_re[, "Q2.5",      "Intercept"],
  re_q50       = dtm_site_re[, "Q50",       "Intercept"],
  re_q975      = dtm_site_re[, "Q97.5",     "Intercept"]
)
dtm_re_df[, sort_key := suppressWarnings(as.integer(tracker_site))]
setorder(dtm_re_df, sort_key, na.last = TRUE)
dtm_re_df[, sort_key := NULL]
fwrite(dtm_re_df, file.path(manuscript_tables_dir,
                            "section_L_s21_dtm_site_re_summary.csv"))
log_progress(sprintf("  DTM site RE range (posterior mean): [%+.3f, %+.3f]",
                     min(dtm_re_df$re_mean), max(dtm_re_df$re_mean)))

# 17 main-effect predictors (same set as CHM per section_10 formula construction)
predictors_17 <- c("slope_mean_z", "rh_98_z", "wsci_z", "cover_z",
                   "meta_leafon_z", "meta_sunel_z", "meta_offnad_z",
                   "meta_relgeo_z", "view_az_cos_z", "meta_absgeo_z",
                   "meta_az_conc_z", "meta_fwdrev_z", "view_az_sin_z",
                   "slope_sd_z", "aspect_cos_z", "aspect_sin_z",
                   "meta_stereo_z")

present <- intersect(predictors_17, names(dtm_df_18))
missing_preds <- setdiff(predictors_17, names(dtm_df_18))
log_progress(sprintf("  Predictors present in dtm_df: %d of %d",
                     length(present), length(predictors_17)))
if (length(missing_preds) > 0) {
  log_progress(sprintf("  Predictors MISSING from dtm_df: %s",
                       paste(missing_preds, collapse = ", ")))
}

# Aggregate per-site means over 18-site DTM data
site_means <- dtm_df_18[, lapply(.SD, mean, na.rm = TRUE),
                        by = c(SITE_COL),
                        .SDcols = present]
setnames(site_means, SITE_COL, "tracker_site")
site_means[, tracker_site := as.character(tracker_site)]
fwrite(site_means, file.path(manuscript_tables_dir,
                             "section_L_s21_dtm_site_means_full17.csv"))

# Bivariate r^2 against site REs
dtm_re_df[, tracker_site := as.character(tracker_site)]
merged <- merge(dtm_re_df[, .(tracker_site, re_mean)],
                site_means, by = "tracker_site", all.x = TRUE)

biv <- rbindlist(lapply(present, function(cv) {
  x <- merged[[cv]]
  y <- merged$re_mean
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 4) return(NULL)
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "pearson"))
  data.table(
    predictor = cv,
    n_sites   = sum(ok),
    pearson_r = unname(ct$estimate),
    p_value   = ct$p.value,
    r_squared = unname(ct$estimate)^2
  )
}))
biv <- biv[order(-r_squared)]
n_tests    <- nrow(biv)
bonf_alpha <- 0.05 / n_tests
biv[, bonferroni_significant := p_value < bonf_alpha]
fwrite(biv, file.path(manuscript_tables_dir,
                      "section_L_s21_dtm_full17_predictors.csv"))

log_progress(sprintf("  DTM N predictors tested: %d", n_tests))
log_progress(sprintf("  Bonferroni alpha at %d tests: %.5f",
                     n_tests, bonf_alpha))
log_progress("  Top 5 DTM predictors by r^2:")
for (i in seq_len(min(5, nrow(biv)))) {
  sig <- if (biv$bonferroni_significant[i]) " *BONF*" else ""
  log_progress(sprintf("    %-15s | n = %2d | r = %+.3f | p = %.4f | r^2 = %.3f%s",
                       biv$predictor[i], biv$n_sites[i],
                       biv$pearson_r[i], biv$p_value[i],
                       biv$r_squared[i], sig))
}

# Specifically pull cover_z to compare against v04's claim
cover_row <- biv[predictor == "cover_z"]
if (nrow(cover_row) > 0) {
  log_progress(sprintf("  cover_z (v04 caption: r^2 = 0.24):"))
  log_progress(sprintf("    r = %+.3f, r^2 = %.3f, p = %.4f",
                       cover_row$pearson_r, cover_row$r_squared,
                       cover_row$p_value))
}

# ---------------------------------------------------------------------
# (S18 DTM) PPC draws — moderate cost
# ---------------------------------------------------------------------

log_subsection("S18 DTM: posterior predictive draws")

set.seed(20260512)
dtm_train <- as.data.table(fit_dtm$data)
ppc_n <- min(5000L, nrow(dtm_train))
ppc_idx <- sort(sample(seq_len(nrow(dtm_train)), ppc_n))
dtm_ppc_data <- dtm_train[ppc_idx]
log_progress(sprintf("  DTM PPC subsample: %d obs from %d training rows",
                     ppc_n, nrow(dtm_train)))

n_draws <- 100
log_progress("  Computing posterior_predict (slow ~5 min)")
dtm_yrep <- posterior_predict(fit_dtm, newdata = dtm_ppc_data,
                              ndraws = n_draws,
                              allow_new_levels = FALSE)
log_progress(sprintf("  yrep matrix: %d draws x %d obs",
                     nrow(dtm_yrep), ncol(dtm_yrep)))

dtm_obs <- dtm_ppc_data$dtm_error_mean
dtm_ppc_long <- data.table(
  obs_idx  = rep(seq_len(ppc_n), each = n_draws),
  draw     = rep(seq_len(n_draws), times = ppc_n),
  yrep     = as.vector(t(dtm_yrep)),
  observed = rep(dtm_obs, each = n_draws)
)
fwrite(dtm_ppc_long, file.path(manuscript_tables_dir,
                              "section_L_s18_dtm_ppc_draws.csv"))

# ---------------------------------------------------------------------
# (S20 DTM) Residuals + variogram — slowest, last
# ---------------------------------------------------------------------

log_subsection("S20 DTM: residuals + variogram")

subsample_n <- min(50000L, nrow(dtm_filtered))
idx <- sort(sample(seq_len(nrow(dtm_filtered)), subsample_n))
sub <- dtm_filtered[idx]
log_progress(sprintf("  Subsample: %d rows", nrow(sub)))

# Recast factor columns to fit's level sets
sub[, lc_l1_code := factor(as.character(lc_l1_code), levels = lc_levels_fit)]
sub[, ecoregion  := factor(as.character(ecoregion),  levels = ec_levels_fit)]
sub[, site       := factor(as.character(site),       levels = site_levels_fit)]

# Defensive NA-drop after recast
n_na <- sum(is.na(sub$lc_l1_code) | is.na(sub$ecoregion) | is.na(sub$site))
if (n_na > 0) {
  log_progress(sprintf("  WARNING: %d NAs after factor recast; dropping", n_na))
  sub <- sub[!is.na(lc_l1_code) & !is.na(ecoregion) & !is.na(site)]
  log_progress(sprintf("  %d rows remain", nrow(sub)))
}

log_progress("  brms::residuals() on subsample (slow step ~3-5 min)")
res <- residuals(fit_dtm, newdata = sub, summary = TRUE,
                 allow_new_levels = FALSE)
sub[, resid_mean := res[, "Estimate"]]

log_progress(sprintf("  Residuals: mean = %+.3f m, sd = %.3f m, range [%+.2f, %+.2f]",
                     mean(sub$resid_mean, na.rm = TRUE),
                     sd(sub$resid_mean, na.rm = TRUE),
                     min(sub$resid_mean, na.rm = TRUE),
                     max(sub$resid_mean, na.rm = TRUE)))

# Reproject WGS84 -> Albers
log_progress(sprintf("  Reprojecting from EPSG:%d to EPSG:%d",
                     COORD_CRS_IN, COORD_CRS_OUT))
sf_sub <- sf::st_as_sf(sub, coords = c(COORD_X, COORD_Y),
                       crs = COORD_CRS_IN)
sf_sub <- sf::st_transform(sf_sub, COORD_CRS_OUT)
alb <- sf::st_coordinates(sf_sub)
sub[, x_alb := alb[, 1]]
sub[, y_alb := alb[, 2]]

out_resid <- sub[is.finite(resid_mean) & is.finite(x_alb) & is.finite(y_alb),
                 .(site       = as.character(site),
                   ecoregion  = as.character(ecoregion),
                   lc_l1_code = as.character(lc_l1_code),
                   x_alb, y_alb, resid_mean)]
fwrite(out_resid, file.path(manuscript_tables_dir,
                            "section_L_s20_dtm_residuals.csv"))
log_progress(sprintf("  Saved %d residual rows", nrow(out_resid)))

# Variogram
log_progress("  Fitting empirical variogram (cutoff 20 km, width 500 m)")
sp_df <- sf::as_Spatial(
  sf::st_as_sf(out_resid, coords = c("x_alb", "y_alb"),
              crs = COORD_CRS_OUT)
)
sp_df$resid_mean <- out_resid$resid_mean
v_emp <- gstat::variogram(resid_mean ~ 1, sp_df,
                         cutoff = 20000, width = 500)
log_progress(sprintf("  Empirical variogram: %d distance bins", nrow(v_emp)))

log_progress("  Fitting spherical variogram model")
v_init <- gstat::vgm(
  psill  = var(sp_df$resid_mean, na.rm = TRUE) * 0.5,
  model  = "Sph",
  range  = 3300,  # v04 reported 3.3 km for DTM
  nugget = var(sp_df$resid_mean, na.rm = TRUE) * 0.5
)
v_fit <- tryCatch(
  gstat::fit.variogram(v_emp, v_init),
  error = function(e) {
    log_progress(sprintf("  Sph fit failed: %s; falling back to Exp",
                         e$message))
    v_init$model <- "Exp"
    gstat::fit.variogram(v_emp, v_init)
  }
)

v_nugget       <- v_fit$psill[1]
v_partial_sill <- v_fit$psill[2]
v_total_sill   <- v_nugget + v_partial_sill
v_range        <- v_fit$range[2]
v_nug_sill     <- v_nugget / v_total_sill
v_model        <- as.character(v_fit$model[2])

log_progress("  --- DTM variogram fit ---")
log_progress(sprintf("    Model:           %s", v_model))
log_progress(sprintf("    Nugget:          %.4f m^2", v_nugget))
log_progress(sprintf("    Partial sill:    %.4f m^2", v_partial_sill))
log_progress(sprintf("    Total sill:      %.4f m^2", v_total_sill))
log_progress(sprintf("    Nugget/sill:     %.4f", v_nug_sill))
log_progress(sprintf("    Range:           %.1f m (%.3f km)",
                     v_range, v_range / 1000))

var_summary <- data.table(
  product        = "DTM_18_site",
  model          = v_model,
  nugget         = v_nugget,
  partial_sill   = v_partial_sill,
  total_sill     = v_total_sill,
  nugget_to_sill = v_nug_sill,
  range_m        = v_range,
  range_km       = v_range / 1000,
  n_residuals    = length(sp_df$resid_mean),
  cutoff_m       = 20000
)
fwrite(var_summary, file.path(manuscript_tables_dir,
                              "section_L_s20_dtm_variogram_params.csv"))
fwrite(as.data.table(v_emp),
       file.path(manuscript_tables_dir,
                "section_L_s20_dtm_empirical_variogram.csv"))

# ---------------------------------------------------------------------
# Console summary
# ---------------------------------------------------------------------

log_progress("================================================================")
log_progress("SUMMARY — DTM-side dual-panel data verified")
log_progress("================================================================")
log_progress("")
log_progress("[S21 DTM right column — verifies v04 'canopy cover r^2 = 0.24']")
if (nrow(cover_row) > 0) {
  log_progress(sprintf("  cover_z bivariate: r = %+.3f, r^2 = %.3f, p = %.4f",
                       cover_row$pearson_r, cover_row$r_squared,
                       cover_row$p_value))
}
log_progress(sprintf("  Top DTM predictor by r^2: %s (r^2 = %.3f)",
                     biv$predictor[1], biv$r_squared[1]))
log_progress("")
log_progress("[S20 DTM right panel — verifies v04 nugget/sill 0.76, range 3.3 km]")
log_progress(sprintf("  Model:        %s   (v04 Sph)", v_model))
log_progress(sprintf("  Nugget/sill:  %.3f   (v04: 0.76)", v_nug_sill))
log_progress(sprintf("  Range:        %.2f km   (v04: 3.3 km)", v_range / 1000))
log_progress("")
log_progress("[S18 DTM panel b — data prepared for density overlay]")
log_progress(sprintf("  PPC subsample: %d obs x %d draws", ppc_n, n_draws))
log_progress(sprintf("  Coverage already verified: 94.44%%"))
log_progress("")
log_progress("Outputs written:")
log_progress("  section_L_s18_dtm_ppc_draws.csv          (density overlay data)")
log_progress("  section_L_s20_dtm_residuals.csv          (for plotting)")
log_progress("  section_L_s20_dtm_variogram_params.csv")
log_progress("  section_L_s20_dtm_empirical_variogram.csv")
log_progress("  section_L_s21_dtm_site_re_summary.csv")
log_progress("  section_L_s21_dtm_site_means_full17.csv")
log_progress("  section_L_s21_dtm_full17_predictors.csv")
log_progress("================================================================")
log_progress("Complete. All holdout cluster work done.")
log_progress("================================================================")
