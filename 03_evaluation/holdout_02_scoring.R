#!/usr/bin/env Rscript
# =====================================================================
# holdout_02_scoring.R
#
# Final holdout scoring. Fixes the validate_newdata() factor-level error from
# section_L2_final.R by filtering chm_df_16 to factor levels present
# in fit_chm$data before passing to brms::residuals().
#
# Specific fixes vs v1:
#   1. Filter rows where lc_l1_code or ecoregion has a level not in
#      fit_chm$data (LMS is the case in point; also defensive for
#      ecoregion).
#   2. Recast all three grouping-factor columns of the subsample
#      explicitly to factor() with fit's level set, so the factor
#      `levels()` attribute matches what brms expects (not just the
#      unique values).
#   3. Variance-partition pattern excludes `sigma_Intercept$` so the
#      distributional random intercepts on sigma aren't double-counted
#      in the 3-way response decomposition.
#   4. Run cheap deterministic steps (S12 rhat/ess, S22 var partition,
#      S4 conditional effects) BEFORE the slow risky residuals step,
#      so partial failure still saves usable outputs.
#
# Verified column names from holdout_diagnostic.R:
#   Response: chm_error_mean
#   Coords:   x, y in WGS84 (EPSG:4326), reproject to EPSG:5070
#   Site:     'site' (16 tracker IDs: 4-15, 17-20)
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

RESPONSE_COL  <- "chm_error_mean"
COORD_X       <- "x"
COORD_Y       <- "y"
SITE_COL      <- "site"
COORD_CRS_IN  <- 4326   # WGS84 geographic
COORD_CRS_OUT <- 5070   # NAD83 / CONUS Albers

log_progress("================================================================")
log_progress("Filtered residuals and variogram")
log_progress("  (Order: S12 -> S22 -> S4 -> S20; slow step last)")
log_progress("================================================================")

# ---------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------

log_progress("Loading CHM 16-site brmsfit")
cp_chm <- load_checkpoint("10b_chm_sensitivity")
fit_chm <- cp_chm$fit_chm_16

log_progress("Loading CHM data ingest")
di <- load_checkpoint("01_data_ingest")
chm_df <- as.data.table(di$chm_df)

# ---------------------------------------------------------------------
# (S12) MCMC convergence diagnostics — no newdata needed, do first
# ---------------------------------------------------------------------

log_subsection("S12 prep: MCMC convergence diagnostics")

rhats      <- brms::rhat(fit_chm)
ess_ratios <- brms::neff_ratio(fit_chm)

log_progress(sprintf("  Total parameters with rhat: %d", length(rhats)))
log_progress(sprintf("  rhat: max = %.4f, p95 = %.4f, p99 = %.4f",
                     max(rhats, na.rm = TRUE),
                     quantile(rhats, 0.95, na.rm = TRUE),
                     quantile(rhats, 0.99, na.rm = TRUE)))
log_progress(sprintf("  neff_ratio: min = %.4f, p05 = %.4f",
                     min(ess_ratios, na.rm = TRUE),
                     quantile(ess_ratios, 0.05, na.rm = TRUE)))
log_progress(sprintf("  N params with rhat > 1.01: %d",
                     sum(rhats > 1.01, na.rm = TRUE)))

mcmc_dt <- data.table(
  parameter  = names(rhats),
  rhat       = unname(rhats),
  neff_ratio = unname(ess_ratios[names(rhats)])
)
fwrite(mcmc_dt, file.path(manuscript_tables_dir,
                         "section_L_s12_chm_mcmc_diagnostics.csv"))

# ---------------------------------------------------------------------
# (S22) Variance partition posterior — no newdata needed, do early
# ---------------------------------------------------------------------

log_subsection("S22 prep: variance partition posterior")

draws <- as_draws_df(fit_chm)
sd_all_cols <- grep("^sd_.*__Intercept$", names(draws), value = TRUE)
log_progress(sprintf("  All sd__Intercept-suffixed parameters: %s",
                     paste(sd_all_cols, collapse = ", ")))

# Exclude sigma distributional random intercepts (sd_*__sigma_Intercept)
sd_int_cols <- sd_all_cols[!grepl("sigma_Intercept$", sd_all_cols)]
log_progress(sprintf("  Response intercept sds (used for §2D-style partition): %s",
                     paste(sd_int_cols, collapse = ", ")))

if (length(sd_int_cols) == 0) {
  log_progress("  WARNING: no response intercept sd parameters found; skipping S22")
} else {
  v_draws_list <- list()
  for (cn in sd_int_cols) {
    grp_name <- gsub("^sd_(.+)__Intercept$", "\\1", cn)
    v_draws_list[[grp_name]] <- draws[[cn]]^2
  }
  v_draws <- as.data.table(v_draws_list)
  total_re_var <- rowSums(as.matrix(v_draws))

  prop_mat <- as.matrix(v_draws) / total_re_var
  prop_summary <- data.table(
    grouping   = colnames(prop_mat),
    share_mean = apply(prop_mat, 2, mean),
    share_q025 = apply(prop_mat, 2, function(x) quantile(x, 0.025)),
    share_q50  = apply(prop_mat, 2, function(x) quantile(x, 0.50)),
    share_q975 = apply(prop_mat, 2, function(x) quantile(x, 0.975))
  )
  fwrite(prop_summary, file.path(manuscript_tables_dir,
                                "section_L_s22_chm_variance_partition.csv"))

  log_progress("  Variance share posterior (cross-check vs §2D 56/42/2):")
  for (i in seq_len(nrow(prop_summary))) {
    log_progress(sprintf("    %-14s | mean = %.3f | 95%% CI [%.3f, %.3f]",
                         prop_summary$grouping[i],
                         prop_summary$share_mean[i],
                         prop_summary$share_q025[i],
                         prop_summary$share_q975[i]))
  }
}

# ---------------------------------------------------------------------
# (S4) Conditional effects — no newdata needed
# ---------------------------------------------------------------------

log_subsection("S4: CHM 16-site conditional effects")
log_progress("  brms::conditional_effects() on CHM fit")
ce <- conditional_effects(fit_chm)
log_progress(sprintf("  Generated %d panels", length(ce)))
saveRDS(ce, file.path(manuscript_tables_dir,
                    "section_L_s4_chm_conditional_effects.rds"))

log_progress("")
log_progress("  >>> Fast deterministic outputs saved.")
log_progress("  >>> Proceeding to slow residuals + variogram step.")
log_progress("")

# ---------------------------------------------------------------------
# (S20) Filtered residuals + variogram — the slow risky step
# ---------------------------------------------------------------------

log_subsection("S20: filtered residuals + variogram")

# Subset to 16 sites
fit_sites <- unique(as.character(fit_chm$data[[SITE_COL]]))
chm_df_16 <- chm_df[as.character(get(SITE_COL)) %in% fit_sites]
log_progress(sprintf("  16-site frame: %d rows", nrow(chm_df_16)))

# Helper: get level set of a column (works for factor or character)
get_levels <- function(x) {
  if (is.factor(x)) levels(x)
  else sort(unique(as.character(x)))
}

# Extract level sets from fit
lc_levels_fit   <- get_levels(fit_chm$data$lc_l1_code)
ec_levels_fit   <- get_levels(fit_chm$data$ecoregion)
site_levels_fit <- get_levels(fit_chm$data$site)

log_progress(sprintf("  Fit lc_l1_code levels (%d): %s",
                     length(lc_levels_fit),
                     paste(sort(lc_levels_fit), collapse = ", ")))
log_progress(sprintf("  Fit ecoregion levels (%d total)", length(ec_levels_fit)))
log_progress(sprintf("  Fit site levels (%d): %s",
                     length(site_levels_fit),
                     paste(sort(site_levels_fit), collapse = ", ")))

# Identify and report any extra levels in chm_df_16
lc_in_chm <- get_levels(chm_df_16$lc_l1_code)
ec_in_chm <- get_levels(chm_df_16$ecoregion)
lc_extra <- setdiff(lc_in_chm, lc_levels_fit)
ec_extra <- setdiff(ec_in_chm, ec_levels_fit)
log_progress(sprintf("  lc levels in chm_df_16 not in fit: %s",
                     if (length(lc_extra) == 0) "(none)"
                     else paste(lc_extra, collapse = ", ")))
log_progress(sprintf("  ecoregion levels in chm_df_16 not in fit: %s",
                     if (length(ec_extra) == 0) "(none)"
                     else paste(ec_extra, collapse = ", ")))

# Filter rows whose values aren't in fit's level sets
n_before <- nrow(chm_df_16)
chm_filtered <- chm_df_16[
  as.character(lc_l1_code) %in% lc_levels_fit &
  as.character(ecoregion)  %in% ec_levels_fit
]
n_after <- nrow(chm_filtered)
log_progress(sprintf("  After filtering to fit's factor levels: %d rows (%d dropped, %.2f%%)",
                     n_after, n_before - n_after,
                     100 * (n_before - n_after) / n_before))

# Subsample 50k
set.seed(20260512)
subsample_n <- min(50000L, nrow(chm_filtered))
idx <- sort(sample(seq_len(nrow(chm_filtered)), subsample_n))
sub <- chm_filtered[idx]
log_progress(sprintf("  Subsample for residuals + variogram: %d rows", nrow(sub)))

# Recast factor columns with fit's level set (defensive — the `levels()`
# attribute can carry orphan levels even when no rows hold them)
sub[, lc_l1_code := factor(as.character(lc_l1_code), levels = lc_levels_fit)]
sub[, ecoregion  := factor(as.character(ecoregion),  levels = ec_levels_fit)]
sub[, site       := factor(as.character(site),       levels = site_levels_fit)]

# Final sanity: confirm no NAs introduced
n_na_lc <- sum(is.na(sub$lc_l1_code))
n_na_ec <- sum(is.na(sub$ecoregion))
n_na_site <- sum(is.na(sub$site))
if (n_na_lc + n_na_ec + n_na_site > 0) {
  log_progress(sprintf("  WARNING: NAs after recasting: lc=%d, ec=%d, site=%d",
                       n_na_lc, n_na_ec, n_na_site))
  sub <- sub[!is.na(lc_l1_code) & !is.na(ecoregion) & !is.na(site)]
  log_progress(sprintf("  Dropped NAs; %d rows remain", nrow(sub)))
}

# Compute residuals on the cleaned subsample
log_progress("  brms::residuals() on subsample (slow step ~3-10 min)")
res <- residuals(fit_chm, newdata = sub, summary = TRUE,
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

# Save residuals (reusable for step 03 S10 per-site histograms)
out_resid <- sub[is.finite(resid_mean) & is.finite(x_alb) & is.finite(y_alb),
                .(site       = as.character(site),
                  ecoregion  = as.character(ecoregion),
                  lc_l1_code = as.character(lc_l1_code),
                  x_alb, y_alb, resid_mean)]
fwrite(out_resid, file.path(manuscript_tables_dir,
                           "section_L_s20_chm_residuals.csv"))
log_progress(sprintf("  Saved %d residual rows for S20 + S10 use",
                     nrow(out_resid)))

# Empirical variogram
log_progress("  Fitting empirical variogram (cutoff 20 km, width 500 m)")
sp_df <- sf::as_Spatial(
  sf::st_as_sf(out_resid, coords = c("x_alb", "y_alb"),
              crs = COORD_CRS_OUT)
)
sp_df$resid_mean <- out_resid$resid_mean
v_emp <- gstat::variogram(resid_mean ~ 1, sp_df,
                         cutoff = 20000, width = 500)
log_progress(sprintf("  Empirical variogram: %d distance bins", nrow(v_emp)))

# Spherical fit (fall back to exponential)
log_progress("  Fitting spherical variogram model")
v_init <- gstat::vgm(
  psill  = var(sp_df$resid_mean, na.rm = TRUE) * 0.5,
  model  = "Sph",
  range  = 3700,
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

log_progress("  --- variogram fit ---")
log_progress(sprintf("    Model:           %s", v_model))
log_progress(sprintf("    Nugget:          %.4f m^2", v_nugget))
log_progress(sprintf("    Partial sill:    %.4f m^2", v_partial_sill))
log_progress(sprintf("    Total sill:      %.4f m^2", v_total_sill))
log_progress(sprintf("    Nugget/sill:     %.4f", v_nug_sill))
log_progress(sprintf("    Range:           %.1f m (%.3f km)",
                     v_range, v_range / 1000))

var_summary <- data.table(
  product        = "CHM_16_site",
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
                             "section_L_s20_chm_variogram_params.csv"))
fwrite(as.data.table(v_emp),
       file.path(manuscript_tables_dir,
                "section_L_s20_chm_empirical_variogram.csv"))

# ---------------------------------------------------------------------
# Console summary
# ---------------------------------------------------------------------

log_progress("================================================================")
log_progress("SUMMARY")
log_progress("================================================================")
log_progress("")
log_progress("[S20 CHM variogram — GATES CAPTION REWRITE]")
log_progress(sprintf("  Model:           %s   (v04 reported Spherical)", v_model))
log_progress(sprintf("  Nugget/sill:     %.3f   (v04 reported 0.79)", v_nug_sill))
log_progress(sprintf("  Range:           %.2f km   (v04 reported 3.7 km)",
                     v_range / 1000))
log_progress(sprintf("  Residuals used:  %d (16-site, fit-compatible LC levels)",
                     length(sp_df$resid_mean)))
log_progress("")
log_progress("[S4 conditional effects]")
log_progress(sprintf("  %d panels saved for plotting", length(ce)))
log_progress("")
log_progress("[S12 MCMC diagnostics]")
log_progress(sprintf("  Max rhat:        %.4f", max(rhats, na.rm = TRUE)))
log_progress(sprintf("  Min neff_ratio:  %.4f", min(ess_ratios, na.rm = TRUE)))
log_progress(sprintf("  N rhat > 1.01:   %d", sum(rhats > 1.01, na.rm = TRUE)))
log_progress("")
if (exists("prop_summary")) {
  log_progress("[S22 variance partition — cross-check vs §2D]")
  for (i in seq_len(nrow(prop_summary))) {
    log_progress(sprintf("  %-14s share: %.1f%%  [%.1f%%, %.1f%%]",
                         prop_summary$grouping[i],
                         prop_summary$share_mean[i] * 100,
                         prop_summary$share_q025[i] * 100,
                         prop_summary$share_q975[i] * 100))
  }
}
log_progress("")
log_progress("Outputs written:")
log_progress("  section_L_s20_chm_residuals.csv          (reusable for step 03 S10)")
log_progress("  section_L_s20_chm_variogram_params.csv")
log_progress("  section_L_s20_chm_empirical_variogram.csv")
log_progress("  section_L_s4_chm_conditional_effects.rds")
log_progress("  section_L_s12_chm_mcmc_diagnostics.csv")
log_progress("  section_L_s22_chm_variance_partition.csv")
log_progress("================================================================")
log_progress("All five critical-figure caption numerics in hand.")
log_progress("Scoring complete.")
log_progress("================================================================")
