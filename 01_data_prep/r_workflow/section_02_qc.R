# =====================================================================
# section_02_qc.R (ENHANCED VERSION)
# Quality control audits and validation
#
# ENHANCEMENTS ADDED:
#   1. Distribution shape tests (kurtosis, skewness, Shapiro-Wilk)
#      - Documents justification for Student-t likelihood
#   2. Predictor range validation
#      - Catches unrealistic values early
#   3. Early outlier detection (flags, not removes)
# =====================================================================

# Ensure dependencies loaded
if (!exists("log_progress")) source("analysis_utils.R")
if (!exists("PROJECT_ROOT")) source("analysis_config.R")

# Load moments package for skewness/kurtosis
if (!requireNamespace("moments", quietly = TRUE)) {
  log_progress("Installing moments package...")
  install.packages("moments", repos = "https://cloud.r-project.org")
}
library(moments)

# Load or define metadata column lists
if (!exists("available_meta_chm") || !exists("available_meta_dtm")) {
  log_progress("⚠ Metadata lists not in environment. Attempting to load from checkpoint...")
  if (checkpoint_exists("01_data_ingest")) {
    data <- load_checkpoint("01_data_ingest")
    available_meta_chm <- data$available_meta_chm
    available_meta_dtm <- data$available_meta_dtm
    if (!exists("chm_df")) chm_df <- data$chm_df
    if (!exists("dtm_df")) dtm_df <- data$dtm_df
  } else {
    log_progress("⚠ No checkpoint found. Using default metadata columns...")
    available_meta_chm <- c("meta_offnad_z", "meta_sunel_z", "meta_az_conc_z", 
                            "meta_absgeo_z", "meta_relgeo_z", "meta_tot_z", 
                            "meta_stereo_z", "meta_fwdrev_z", "meta_leafon_z")
    available_meta_dtm <- available_meta_chm
  }
}

# Ensure data is loaded
if (!exists("chm_df") || !exists("dtm_df")) {
  if (checkpoint_exists("01_data_ingest")) {
    data <- load_checkpoint("01_data_ingest")
    chm_df <- data$chm_df
    dtm_df <- data$dtm_df
  } else {
    stop("Data not available. Run: source('section_01_ingest.R') first")
  }
}

log_progress("Running QC audits (ENHANCED)...")

# =====================================================================
# 1. SCHEMA AUDIT (Original)
# =====================================================================

audit_schema <- function(df, name, cols_types) {
  miss <- setdiff(names(cols_types), names(df))
  if (length(miss)) warning(sprintf("[%s] missing columns: %s", name, paste(miss, collapse=", ")))
  out <- tibble::tibble(
    column = names(cols_types),
    expected = unname(cols_types),
    present  = column %in% names(df),
    class    = vapply(column, function(c) if (c %in% names(df)) paste(class(df[[c]]), collapse="|") else NA_character_, character(1))
  )
  readr::write_csv(out, file.path(out_tables, sprintf("qc_schema_%s.csv", name)))
  invisible(out)
}

core_cols_chm <- c(
  chm_error_mean="numeric", w_chm="numeric",
  slope_mean="numeric", slope_mean_z="numeric", slope_sd_z="numeric",
  wsci_z="numeric", rh_98_z="numeric", cover_z="numeric",
  aspect_sin_z="numeric", aspect_cos_z="numeric",
  view_az_sin="numeric", view_az_cos="numeric",
  lc_l1_code="factor", ecoregion="factor", site="factor",
  als_chm_valid_frac="numeric", p3d_chm_valid_frac="numeric"
)

core_cols_dtm <- c(
  dtm_error_mean="numeric", w_dtm="numeric",
  slope_mean="numeric", slope_mean_z="numeric", slope_sd_z="numeric",
  wsci_z="numeric", rh_98_z="numeric", cover_z="numeric",
  aspect_sin_z="numeric", aspect_cos_z="numeric",
  view_az_sin="numeric", view_az_cos="numeric",
  lc_l1_code="factor", ecoregion="factor", site="factor",
  dep_dtm_valid_frac="numeric", p3d_dtm_valid_frac="numeric"
)

meta_cols_check_chm <- rep("numeric", length(available_meta_chm))
names(meta_cols_check_chm) <- available_meta_chm

meta_cols_check_dtm <- rep("numeric", length(available_meta_dtm))
names(meta_cols_check_dtm) <- available_meta_dtm

cols_types_chm <- c(core_cols_chm, meta_cols_check_chm)
cols_types_dtm <- c(core_cols_dtm, meta_cols_check_dtm)

audit_schema(chm_df, "chm", cols_types_chm)
audit_schema(dtm_df, "dtm", cols_types_dtm)

# =====================================================================
# 2. MISSINGNESS AUDIT (Original)
# =====================================================================

audit_missingness <- function(df, name) {
  nn <- sapply(df, function(x) sum(is.na(x)))
  frac <- nn / nrow(df)
  out <- tibble::tibble(var = names(df), n_na = as.integer(nn), frac = as.numeric(frac)) %>%
    dplyr::arrange(desc(frac), var)
  readr::write_csv(out, file.path(out_tables, sprintf("qc_missingness_%s.csv", name)))
  invisible(out)
}
audit_missingness(chm_df, "chm")
audit_missingness(dtm_df, "dtm")

# =====================================================================
# 3. DUPLICATE AUDIT (Original)
# =====================================================================

dup_chm <- chm_df %>% dplyr::count(site, shot_number) %>% dplyr::filter(n > 1)
dup_dtm <- dtm_df %>% dplyr::count(site, shot_number) %>% dplyr::filter(n > 1)
write_csv(dup_chm, file.path(out_tables, "qc_duplicates_chm.csv"))
write_csv(dup_dtm, file.path(out_tables, "qc_duplicates_dtm.csv"))

# =====================================================================
# 4. VALID FRACTION HISTOGRAMS (Original)
# =====================================================================

vf_chm <- chm_df %>% 
  dplyr::transmute(als=als_chm_valid_frac, p3d=p3d_chm_valid_frac, w=w_chm) %>% 
  tidyr::pivot_longer(everything())

vf_dtm <- dtm_df %>% 
  dplyr::transmute(dep=dep_dtm_valid_frac, p3d=p3d_dtm_valid_frac, w=w_dtm) %>% 
  tidyr::pivot_longer(everything())

p_vf_chm <- ggplot(vf_chm, aes(value)) + 
  geom_histogram(bins=60) + 
  facet_wrap(~name, nrow=1) + 
  labs(title="CHM: valid fractions & weights")

p_vf_dtm <- ggplot(vf_dtm, aes(value)) + 
  geom_histogram(bins=60) + 
  facet_wrap(~name, nrow=1) + 
  labs(title="DTM: valid fractions & weights")

ggsave(file.path(out_plots,"00_valid_frac_chm.pdf"), p_vf_chm, 
       width=10, height=3.6, bg="white")
ggsave(file.path(out_plots,"00_valid_frac_dtm.pdf"), p_vf_dtm, 
       width=10, height=3.6, bg="white")

# =====================================================================
# 5. NEW: DISTRIBUTION SHAPE DIAGNOSTICS
# =====================================================================
# CRITICAL FOR BAYESIAN MODEL SELECTION
# Documents justification for Student-t likelihood

log_subsection("Distribution shape diagnostics (NEW)")

# Sample for normality tests (Shapiro-Wilk limited to 5000)
set.seed(42)
chm_errors <- chm_df$chm_error_mean[is.finite(chm_df$chm_error_mean)]
dtm_errors <- dtm_df$dtm_error_mean[is.finite(dtm_df$dtm_error_mean)]

sample_chm <- sample(chm_errors, min(5000, length(chm_errors)))
sample_dtm <- sample(dtm_errors, min(5000, length(dtm_errors)))

# Shapiro-Wilk normality tests
shapiro_chm <- shapiro.test(sample_chm)
shapiro_dtm <- shapiro.test(sample_dtm)

log_progress(sprintf("  CHM Shapiro-Wilk: W=%.4f, p=%.2e", 
                     shapiro_chm$statistic, shapiro_chm$p.value))
log_progress(sprintf("  DTM Shapiro-Wilk: W=%.4f, p=%.2e", 
                     shapiro_dtm$statistic, shapiro_dtm$p.value))

# Kurtosis check (CRITICAL for Student-t justification)
# Gaussian kurtosis = 3; values > 3 indicate heavy tails
kurtosis_chm <- moments::kurtosis(chm_errors)
kurtosis_dtm <- moments::kurtosis(dtm_errors)

skewness_chm <- moments::skewness(chm_errors)
skewness_dtm <- moments::skewness(dtm_errors)

log_progress(sprintf("  CHM kurtosis: %.2f (Gaussian=3)", kurtosis_chm))
log_progress(sprintf("  DTM kurtosis: %.2f (Gaussian=3)", kurtosis_dtm))
log_progress(sprintf("  CHM skewness: %.2f", skewness_chm))
log_progress(sprintf("  DTM skewness: %.2f", skewness_dtm))

# Interpretation and recommendations
if (kurtosis_chm > 5 || kurtosis_dtm > 5) {
  log_progress("  → HEAVY TAILS DETECTED: Student-t family RECOMMENDED")
  log_progress("    Kurtosis > 5 indicates significant departure from Gaussian")
}

if (kurtosis_chm > 7 || kurtosis_dtm > 7) {
  log_progress("  → VERY HEAVY TAILS: Kurtosis > 7 confirms Student-t essential")
}

# Create distribution diagnostics table
dist_diagnostics <- tibble::tibble(
  product = c("CHM", "DTM"),
  n_valid = c(length(chm_errors), length(dtm_errors)),
  mean = c(mean(chm_errors), mean(dtm_errors)),
  median = c(median(chm_errors), median(dtm_errors)),
  sd = c(sd(chm_errors), sd(dtm_errors)),
  nmad = c(nm_ad(chm_errors), nm_ad(dtm_errors)),
  shapiro_W = c(shapiro_chm$statistic, shapiro_dtm$statistic),
  shapiro_p = c(shapiro_chm$p.value, shapiro_dtm$p.value),
  skewness = c(skewness_chm, skewness_dtm),
  kurtosis = c(kurtosis_chm, kurtosis_dtm),
  q01 = c(quantile(chm_errors, 0.01), quantile(dtm_errors, 0.01)),
  q99 = c(quantile(chm_errors, 0.99), quantile(dtm_errors, 0.99)),
  heavy_tails = c(kurtosis_chm > 5, kurtosis_dtm > 5),
  student_t_recommended = c(kurtosis_chm > 5 || shapiro_chm$p.value < 0.05,
                            kurtosis_dtm > 5 || shapiro_dtm$p.value < 0.05)
)

write_csv(dist_diagnostics, file.path(out_tables, "qc_distribution_diagnostics.csv"))
log_progress("  ✓ Distribution diagnostics saved")

# Print summary table
log_progress("\n  Distribution Summary:")
print(dist_diagnostics %>% select(product, kurtosis, skewness, shapiro_p, student_t_recommended))

# =====================================================================
# 6. NEW: PREDICTOR RANGE VALIDATION
# =====================================================================
# Catches unrealistic values that could indicate data issues

log_subsection("Predictor range validation (NEW)")

# Define expected ranges for key predictors
predictor_ranges <- list(
  slope_mean = c(0, 90),           # Degrees - cannot exceed 90
  cover = c(0, 1),                 # Fraction
  rh_98 = c(0, 100),               # Meters (reasonable canopy height)
  wsci = c(0, 200)                 # Typical WSCI range
)

# Also check z-scored predictors for extreme values (>5 SD is suspicious)
z_score_threshold <- 5

range_violations <- list()

for (pred in names(predictor_ranges)) {
  if (pred %in% names(chm_df)) {
    vals <- chm_df[[pred]]
    expected_range <- predictor_ranges[[pred]]
    out_of_range <- sum(vals < expected_range[1] | vals > expected_range[2], na.rm = TRUE)
    
    if (out_of_range > 0) {
      log_progress(sprintf("  ⚠ %s: %d values outside expected range [%.1f, %.1f]",
                           pred, out_of_range, expected_range[1], expected_range[2]))
      range_violations[[pred]] <- out_of_range
    }
  }
}

# Check z-scored predictors for extreme values
z_cols <- grep("_z$", names(chm_df), value = TRUE)
extreme_z <- sapply(z_cols, function(col) {
  if (col %in% names(chm_df)) {
    vals <- chm_df[[col]]
    sum(abs(vals) > z_score_threshold, na.rm = TRUE)
  } else {
    0
  }
})

extreme_z <- extreme_z[extreme_z > 0]
if (length(extreme_z) > 0) {
  log_progress("  Z-scored predictors with |z| > 5:")
  for (col in names(extreme_z)) {
    pct <- 100 * extreme_z[col] / nrow(chm_df)
    log_progress(sprintf("    %s: %d (%.3f%%)", col, extreme_z[col], pct))
  }
}

if (length(range_violations) == 0 && length(extreme_z) == 0) {
  log_progress("  ✓ All predictor ranges within expected bounds")
}

# =====================================================================
# 7. NEW: EARLY OUTLIER DETECTION (FLAG, DON'T REMOVE)
# =====================================================================
# Identifies statistical outliers for later investigation

log_subsection("Outlier detection (NEW)")

# IQR-based outlier detection
chm_df <- chm_df %>%
  mutate(
    outlier_iqr = abs(chm_error_mean - median(chm_error_mean, na.rm = TRUE)) > 
                 3 * IQR(chm_error_mean, na.rm = TRUE),
    outlier_extreme = abs(chm_error_mean) > 50  # >50m is likely problematic
  )

dtm_df <- dtm_df %>%
  mutate(
    outlier_iqr = abs(dtm_error_mean - median(dtm_error_mean, na.rm = TRUE)) > 
                 3 * IQR(dtm_error_mean, na.rm = TRUE),
    outlier_extreme = abs(dtm_error_mean) > 50
  )

n_outliers_chm <- sum(chm_df$outlier_iqr, na.rm = TRUE)
n_outliers_dtm <- sum(dtm_df$outlier_iqr, na.rm = TRUE)
n_extreme_chm <- sum(chm_df$outlier_extreme, na.rm = TRUE)
n_extreme_dtm <- sum(dtm_df$outlier_extreme, na.rm = TRUE)

log_progress(sprintf("  CHM: %d IQR outliers (%.2f%%), %d extreme (>50m)",
                     n_outliers_chm, 100 * n_outliers_chm / nrow(chm_df), n_extreme_chm))
log_progress(sprintf("  DTM: %d IQR outliers (%.2f%%), %d extreme (>50m)",
                     n_outliers_dtm, 100 * n_outliers_dtm / nrow(dtm_df), n_extreme_dtm))

# Export outlier summary by site
outlier_by_site <- chm_df %>%
  group_by(site) %>%
  summarise(
    n = n(),
    n_outliers = sum(outlier_iqr, na.rm = TRUE),
    pct_outliers = 100 * n_outliers / n,
    n_extreme = sum(outlier_extreme, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(pct_outliers))

write_csv(outlier_by_site, file.path(out_tables, "qc_outliers_by_site.csv"))

# Flag sites with high outlier rates
high_outlier_sites <- outlier_by_site %>% filter(pct_outliers > 5)
if (nrow(high_outlier_sites) > 0) {
  log_progress("  ⚠ Sites with >5% outlier rate:")
  for (i in 1:nrow(high_outlier_sites)) {
    log_progress(sprintf("    Site %s: %.1f%% outliers", 
                         high_outlier_sites$site[i], 
                         high_outlier_sites$pct_outliers[i]))
  }
}

# =====================================================================
# 8. QC SUMMARY REPORT
# =====================================================================

log_subsection("QC Summary Report")

qc_summary <- tibble::tibble(
  Check = c(
    "Schema completeness (CHM)",
    "Schema completeness (DTM)",
    "Duplicate records (CHM)",
    "Duplicate records (DTM)",
    "Kurtosis - Student-t needed (CHM)",
    "Kurtosis - Student-t needed (DTM)",
    "Outliers flagged (CHM)",
    "Outliers flagged (DTM)"
  ),
  Status = c(
    if (all(cols_types_chm %in% names(chm_df))) "✓ PASS" else "⚠ WARN",
    if (all(cols_types_dtm %in% names(dtm_df))) "✓ PASS" else "⚠ WARN",
    if (nrow(dup_chm) == 0) "✓ PASS" else sprintf("⚠ %d duplicates", nrow(dup_chm)),
    if (nrow(dup_dtm) == 0) "✓ PASS" else sprintf("⚠ %d duplicates", nrow(dup_dtm)),
    if (kurtosis_chm > 5) "✓ YES (heavy tails)" else "○ NO",
    if (kurtosis_dtm > 5) "✓ YES (heavy tails)" else "○ NO",
    sprintf("%d (%.2f%%)", n_outliers_chm, 100 * n_outliers_chm / nrow(chm_df)),
    sprintf("%d (%.2f%%)", n_outliers_dtm, 100 * n_outliers_dtm / nrow(dtm_df))
  )
)

write_csv(qc_summary, file.path(out_tables, "qc_summary_report.csv"))
log_progress("\n  QC Summary:")
print(qc_summary)

log_progress("✓ Enhanced QC audits complete")

# Save updated data with outlier flags
if (!exists("BATCH_MODE") || !BATCH_MODE) {
  log_progress("Saving QC checkpoint...")
  save_checkpoint("02_qc", list(
    qc_complete = TRUE,
    dist_diagnostics = dist_diagnostics,
    outlier_by_site = outlier_by_site,
    kurtosis_chm = kurtosis_chm,
    kurtosis_dtm = kurtosis_dtm,
    student_t_recommended = (kurtosis_chm > 5 || kurtosis_dtm > 5)
  ))
}
