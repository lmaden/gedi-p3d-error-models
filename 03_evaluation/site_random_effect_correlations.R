#!/usr/bin/env Rscript
# =====================================================================
# site_random_effect_correlations.R
#
# Fixes two bugs in the first version of this script:
#
# Bug 1 (silent correctness). The CHM ranef site labels are TRACKER
# IDs (4-15, 17-20) per how the CHM stage 2 fit was specified. The
# site-cov CSV's manuscript_site column has MANUSCRIPT IDs (4-9,
# 11-19). Merging directly produced 14 matches with 3 INCORRECT
# pairings (CHM tracker 17/18/19 silently matched to site-cov
# manuscript 17/18/19, which are tracker 18/19/20). Fix: merge on
# tracker_site, which both tables have via the lookup.
#
# Bug 2. The patch tested all 35 numeric columns including descriptive
# error statistics (p3d_median_bias, alt_mean_bias, var_err_p3d, etc.)
# These are dependent-variable-derived and trivially correlate with
# site RE intercepts. v04 S21 caption framed the figure as testing the
# 17 main-effect predictors. Fix: restrict to predictor site-level
# means (cover_avg, rh_98_avg, wsci_avg, slope_mean_avg, and the
# meta_* covariates that have site means).
#
# Output: writes section_L_s21_chm_site_re_bivariate_r2.csv with the
# corrected and restricted bivariate table (overwrites v1 file).
# =====================================================================

suppressPackageStartupMessages({
  library(data.table)
})

if (!exists("log_progress"))  source("analysis_utils.R")
if (!exists("PROJECT_ROOT"))  source("analysis_config_FIXED.R")

manuscript_tables_dir <- file.path(PROJECT_ROOT, "manuscript_tables")

log_progress("================================================================")
log_progress("S21 patch v2: correct keying + predictor-only restriction")
log_progress("================================================================")

# ---------------------------------------------------------------------
# Load the L-1 CHM site-RE table (its site_label = tracker_id)
# ---------------------------------------------------------------------

re_path <- file.path(manuscript_tables_dir,
                   "section_L_s17_chm_site_re_summary.csv")
chm_re_df <- fread(re_path)
setnames(chm_re_df, "site_label", "tracker_site")
chm_re_df[, tracker_site := as.character(tracker_site)]
log_progress(sprintf("Loaded CHM site-RE table: %d rows", nrow(chm_re_df)))
log_progress(sprintf("  Tracker IDs in CHM ranef: %s",
                     paste(sort(as.integer(chm_re_df$tracker_site)),
                           collapse = ", ")))

# ---------------------------------------------------------------------
# Load the per-site covariates CSV
# ---------------------------------------------------------------------

cov_path <- file.path(manuscript_tables_dir,
                   "groundwork_task4_phase25_site_covariates.csv")
site_cov <- fread(cov_path)
log_progress(sprintf("Loaded site-cov table: %d rows", nrow(site_cov)))
log_progress(sprintf("  Tracker IDs in site-cov: %s",
                     paste(sort(as.integer(site_cov$tracker_site)),
                           collapse = ", ")))
site_cov[, tracker_site := as.character(tracker_site)]

# ---------------------------------------------------------------------
# Merge on tracker_site (the correct key)
# ---------------------------------------------------------------------

merged <- merge(chm_re_df[, .(tracker_site, re_mean)],
                site_cov, by = "tracker_site", all.x = TRUE)
log_progress(sprintf("Merged table: %d rows (expected 16 CHM sites; missing rows are sites not in Phase 2.5 frame)",
                     nrow(merged)))

n_matched <- sum(!is.na(merged$manuscript_site))
log_progress(sprintf("  Matched rows (have site-cov data): %d", n_matched))
unmatched <- merged[is.na(manuscript_site), tracker_site]
if (length(unmatched) > 0) {
  log_progress(sprintf("  Unmatched tracker IDs (CHM ranef without site-cov data): %s",
                       paste(unmatched, collapse = ", ")))
}

# ---------------------------------------------------------------------
# Restrict to PREDICTOR site-level means
# Per §2F, the model's 17 main-effect continuous predictors are:
#   slope_mean_z, rh_98_z, wsci_z, cover_z, meta_leafon_z,
#   meta_sunel_z, meta_offnad_z, meta_relgeo_z, view_az_cos_z,
#   meta_absgeo_z, meta_az_conc_z, meta_fwdrev_z, view_az_sin_z,
#   slope_sd_z, aspect_cos_z, aspect_sin_z, meta_stereo_z.
#
# The per-site CSV has site-level means for a subset. Map by name:
# ---------------------------------------------------------------------

predictor_to_csv <- c(
  slope_mean_z   = "slope_mean_avg",
  rh_98_z        = "rh_98_avg",
  wsci_z         = "wsci_avg",
  cover_z        = "cover_avg",
  meta_leafon_z  = "meta_leaf_on_ratio_avg",
  meta_sunel_z   = "meta_sun_elev_avg_avg",
  meta_offnad_z  = "meta_off_nadir_avg_avg",
  meta_relgeo_z  = "meta_rel_geoacc_avg_avg",
  meta_absgeo_z  = "meta_abs_geoacc_avg_avg",
  meta_az_conc_z = "meta_az_concentration_avg",
  meta_fwdrev_z  = "meta_fwd_rev_ratio_avg",
  meta_stereo_z  = "meta_stereo_ratio_avg"
  # NOT in CSV at site-level: view_az_cos_z, view_az_sin_z,
  # slope_sd_z, aspect_cos_z, aspect_sin_z
)

predictors_present <- predictor_to_csv[predictor_to_csv %in% names(merged)]
predictors_missing <- setdiff(predictor_to_csv,
                             names(merged)[names(merged) %in% predictor_to_csv])

log_progress(sprintf("Predictor site-means present in CSV: %d of %d (%s)",
                     length(predictors_present),
                     length(predictor_to_csv),
                     paste(predictors_present, collapse = ", ")))
if (length(predictors_missing) > 0) {
  log_progress(sprintf("Predictors NOT in CSV (would need aggregation from data ingest): %s",
                       paste(predictors_missing, collapse = ", ")))
}

# ---------------------------------------------------------------------
# Compute bivariate r, p, r^2 for predictor covariates only
# ---------------------------------------------------------------------

biv <- rbindlist(lapply(seq_along(predictors_present), function(i) {
  csv_name      <- predictors_present[i]
  predictor_lbl <- names(predictors_present)[i]
  x <- merged[[csv_name]]
  y <- merged$re_mean
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 4) return(NULL)
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "pearson"))
  data.table(
    predictor    = predictor_lbl,
    csv_column   = csv_name,
    n_sites      = sum(ok),
    pearson_r    = unname(ct$estimate),
    p_value      = ct$p.value,
    r_squared    = unname(ct$estimate)^2
  )
}))
biv <- biv[order(-r_squared)]

fwrite(biv, file.path(manuscript_tables_dir,
                     "section_L_s21_chm_site_re_bivariate_r2.csv"))

# Bonferroni-adjusted significance: alpha / n_tests
n_tests <- nrow(biv)
bonf_alpha <- 0.05 / n_tests
biv[, bonferroni_significant := p_value < bonf_alpha]

# ---------------------------------------------------------------------
# Save the merged metadata table (corrected) too
# ---------------------------------------------------------------------

# Build a clean output table: tracker_site, manuscript_site, RE, plus
# only the predictor columns and corr_errP3D_errDTM (relevant to §2H).
keep_cols <- c("tracker_site", "manuscript_site", "re_mean",
              "corr_errP3D_errDTM", unname(predictors_present))
keep_cols <- intersect(keep_cols, names(merged))
out_merged <- merged[, ..keep_cols]
fwrite(out_merged, file.path(manuscript_tables_dir,
                            "section_L_s21_chm_site_re_metadata.csv"))

# ---------------------------------------------------------------------
# Console summary
# ---------------------------------------------------------------------

log_progress("================================================================")
log_progress("SUMMARY")
log_progress("================================================================")
log_progress(sprintf("N predictors tested: %d (out of v04's 17)", n_tests))
log_progress(sprintf("Bonferroni alpha at %d tests: %.4f", n_tests, bonf_alpha))
log_progress("")
log_progress("All predictor bivariates (sorted by r^2 descending):")
for (i in seq_len(nrow(biv))) {
  sig <- if (biv$bonferroni_significant[i]) " *BONF*" else ""
  log_progress(sprintf("  %-15s | csv: %-25s | n = %2d | r = %+.3f | p = %.4f | r^2 = %.3f%s",
                       biv$predictor[i],
                       biv$csv_column[i],
                       biv$n_sites[i],
                       biv$pearson_r[i],
                       biv$p_value[i],
                       biv$r_squared[i],
                       sig))
}
log_progress("")
log_progress("Specific v04 caption reference (CHM: fwd/rev ratio, r^2 = 0.23):")
row <- biv[predictor == "meta_fwdrev_z"]
if (nrow(row) > 0) {
  log_progress(sprintf("  meta_fwdrev_z (meta_fwd_rev_ratio_avg): r = %+.3f, r^2 = %.3f, p = %.4f",
                       row$pearson_r, row$r_squared, row$p_value))
} else {
  log_progress("  meta_fwdrev_z not in result table")
}
log_progress("")
log_progress("Outputs written:")
log_progress(sprintf("  %s/section_L_s21_chm_site_re_metadata.csv",
                     manuscript_tables_dir))
log_progress(sprintf("  %s/section_L_s21_chm_site_re_bivariate_r2.csv (OVERWRITTEN with corrected values)",
                     manuscript_tables_dir))
log_progress("================================================================")
log_progress("S21 patch v2 complete.")
log_progress("================================================================")
