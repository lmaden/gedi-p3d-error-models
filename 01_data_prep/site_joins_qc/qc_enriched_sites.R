#!/usr/bin/env Rscript
# qc_enriched_sites_comprehensive.R
# 
# Comprehensive QC for enriched GEDI site data with zonal sampling
# Designed to identify issues before Bayesian analysis
#
# Major QC Categories:
#   1. Data Integrity (duplicates, formats, constants)
#   2. Value Range Validation (physical bounds)
#   3. Cross-Variable Consistency (logical relationships)
#   4. Zonal Statistics Internal Consistency
#   5. Error Distribution Diagnostics
#   6. Site-Level Anomaly Detection
#   7. Covariate Quality for Bayesian Modeling
#   8. Missing Data Pattern Analysis

suppressPackageStartupMessages({
  library(data.table)
  library(future)
  library(future.apply)
  library(parallelly)
})

# ============================================================================
# CONFIGURATION
# ============================================================================

ENRICHED_DIR <- "/gpfs/data1/vclgp/lmaden/chpt1/data/enriched_by_site"
LOG_DIR      <- "/gpfs/data1/vclgp/lmaden/chpt1/logs"
dir.create(LOG_DIR, showWarnings = FALSE, recursive = TRUE)

SITES <- 1:20
MAX_WORKERS <- 10L
DT_THREADS  <- 2L
data.table::setDTthreads(DT_THREADS)

TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")

# Output files
REPORT_FILE    <- file.path(LOG_DIR, sprintf("qc_comprehensive_report_%s.txt", TIMESTAMP))
SUMMARY_CSV    <- file.path(LOG_DIR, sprintf("qc_site_summary_%s.csv", TIMESTAMP))
VARIABLE_CSV   <- file.path(LOG_DIR, sprintf("qc_variable_summary_%s.csv", TIMESTAMP))
ISSUES_CSV     <- file.path(LOG_DIR, sprintf("qc_issues_%s.csv", TIMESTAMP))
FLAGGED_CSV    <- file.path(LOG_DIR, sprintf("qc_flagged_observations_%s.csv", TIMESTAMP))
EXCLUSIONS_CSV <- file.path(LOG_DIR, sprintf("qc_recommended_exclusions_%s.csv", TIMESTAMP))

# ============================================================================
# COLUMN DEFINITIONS
# ============================================================================

# Core GEDI columns (should always be present and valid)
GEDI_CORE <- c("shot_number", "site", "delta_time")

# GEDI metrics
GEDI_METRICS <- c("rh_98", "elev_low", "cover", "pft", "wsci")

# Zonal product columns
DTM_REF_COLS <- c("dep_dtm_mean", "dep_dtm_median", "dep_dtm_p25", "dep_dtm_p75", 
                  "dep_dtm_p90", "dep_dtm_valid_frac", "dep_dtm")
DTM_TEST_COLS <- c("p3d_dtm_mean", "p3d_dtm_median", "p3d_dtm_p25", "p3d_dtm_p75",
                   "p3d_dtm_p90", "p3d_dtm_valid_frac", "mhrsi_dtm")
CHM_REF_COLS <- c("als_chm_mean", "als_chm_median", "als_chm_sd", "als_chm_p25", 
                  "als_chm_p75", "als_chm_p90", "als_chm_valid_frac", "als_chm")
CHM_TEST_COLS <- c("p3d_chm_mean", "p3d_chm_median", "p3d_chm_sd", "p3d_chm_p25",
                   "p3d_chm_p75", "p3d_chm_p90", "p3d_chm_valid_frac", "p3d_chm")

# Error columns
ERROR_COLS <- c("error_mean", "error_median", "dtm_error_mean", "dtm_error_median",
                "chm_rmse", "dtm_rmse")

# Terrain columns
TERRAIN_COLS <- c("slope_mean", "slope_sd", "slope_valid_frac",
                  "aspect_sin_mean", "aspect_cos_mean", "aspect_valid_frac")

# Ancillary
LC_COLS <- c("lc2022_id_center", "lc2022_id_mode", "lc2022_mode_frac", 
             "lc2022_l1_code", "lc2022_l1_name", "ecoregion")

# Physical bounds for validation
BOUNDS <- list(
  # Elevation bounds (meters) - reasonable global range
  elev_min = -500,
  elev_max = 9000,
  # Canopy height bounds (meters)
  chm_min = -5,      # Allow small negative due to noise
  chm_max = 100,     # Tallest trees ~115m
  # Error bounds (meters) - flag extreme discrepancies
  error_warn = 10,
  error_critical = 25,
  # Terrain
  slope_min = 0,
  slope_max = 90,
  trig_min = -1,
  trig_max = 1,
  # GEDI
  rh98_min = -5,
  rh98_max = 120,
  cover_min = 0,
  cover_max = 1
)

# Thresholds for QC flags
THRESHOLDS <- list(
  # Missing data
  missing_warn_pct = 10,
  missing_critical_pct = 50,
  # Outliers
  outlier_warn_pct = 5,
  outlier_critical_pct = 15,
  # Bias
  bias_warn = 3,      # meters
  bias_critical = 5,
  # Site anomaly (z-score)
  site_zscore_warn = 2,
  site_zscore_critical = 3,
  # Correlation for consistency
  percentile_order_tolerance = 0.01,  # Allow tiny numeric imprecision
  # Skewness/kurtosis for normality
  skew_warn = 2,
  kurt_warn = 7
)

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

determine_workers <- function(max_target = MAX_WORKERS) {
  n <- suppressWarnings(parallelly::availableCores())
  if (is.na(n)) n <- 1L

  max(1L, min(as.integer(n), max_target))
}

find_enriched_csv <- function(site_id) {
  patterns <- c(
    sprintf("site_%02d_enriched.csv.gz", site_id),
    sprintf("site_%02d_enriched.csv", site_id),
    sprintf("site_%d_enriched.csv.gz", site_id),
    sprintf("site_%d_enriched.csv", site_id)
  )
  for (p in patterns) {
    fp <- file.path(ENRICHED_DIR, p)
    if (file.exists(fp)) return(fp)
  }
  NA_character_
}

safe_quantile <- function(x, probs, na.rm = TRUE) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(rep(NA_real_, length(probs)))
  quantile(x, probs = probs, na.rm = na.rm)
}

# IQR-based outlier detection
count_outliers_iqr <- function(x, k = 1.5) {
  x <- x[is.finite(x)]
  if (length(x) < 4) return(list(n_low = 0L, n_high = 0L, pct = 0, 
                                  lower_fence = NA, upper_fence = NA))
  q <- quantile(x, c(0.25, 0.75))
  iqr <- q[2] - q[1]
  lower <- q[1] - k * iqr
  upper <- q[2] + k * iqr
  n_low <- sum(x < lower)
  n_high <- sum(x > upper)
  list(n_low = n_low, n_high = n_high, 
       pct = 100 * (n_low + n_high) / length(x),
       lower_fence = lower, upper_fence = upper)
}

# MAD-based outlier detection (more robust)
count_outliers_mad <- function(x, k = 3) {
  x <- x[is.finite(x)]
  if (length(x) < 4) return(list(n_outliers = 0L, pct = 0))
  med <- median(x)
  mad_val <- mad(x, constant = 1.4826)
  if (mad_val == 0) return(list(n_outliers = 0L, pct = 0))
  n_out <- sum(abs(x - med) > k * mad_val)
  list(n_outliers = n_out, pct = 100 * n_out / length(x))
}

# Comprehensive column statistics
compute_col_stats <- function(x, col_name) {
  x_finite <- x[is.finite(x)]
  n_total <- length(x)
  n_valid <- length(x_finite)
  n_na <- sum(is.na(x))
  n_nan <- sum(is.nan(x))
  n_inf <- sum(is.infinite(x))
  
  if (n_valid == 0) {
    return(data.table(
      column = col_name, n_total = n_total, n_valid = 0L,
      n_na = n_na, n_nan = n_nan, n_inf = n_inf,
      pct_valid = 0, pct_missing = 100,
      mean = NA_real_, median = NA_real_, sd = NA_real_,
      min = NA_real_, max = NA_real_, range = NA_real_,
      p01 = NA_real_, p05 = NA_real_, p10 = NA_real_,
      p25 = NA_real_, p75 = NA_real_, 
      p90 = NA_real_, p95 = NA_real_, p99 = NA_real_,
      iqr = NA_real_, cv = NA_real_,
      skewness = NA_real_, kurtosis = NA_real_,
      n_negative = NA_integer_, pct_negative = NA_real_,
      n_zero = NA_integer_, pct_zero = NA_real_,
      n_unique = NA_integer_,
      outliers_iqr_n = NA_integer_, outliers_iqr_pct = NA_real_,
      outliers_mad_n = NA_integer_, outliers_mad_pct = NA_real_
    ))
  }
  
  pcts <- safe_quantile(x_finite, c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99))
  mu <- mean(x_finite)
  sigma <- sd(x_finite)
  
  # Skewness and kurtosis
  skew <- if (sigma > 0 && n_valid > 2) mean(((x_finite - mu) / sigma)^3) else NA_real_
  kurt <- if (sigma > 0 && n_valid > 3) mean(((x_finite - mu) / sigma)^4) - 3 else NA_real_
  
  outliers_iqr <- count_outliers_iqr(x_finite)
  outliers_mad <- count_outliers_mad(x_finite)
  
  data.table(
    column = col_name,
    n_total = n_total,
    n_valid = n_valid,
    n_na = n_na,
    n_nan = n_nan,
    n_inf = n_inf,
    pct_valid = 100 * n_valid / n_total,
    pct_missing = 100 * (n_total - n_valid) / n_total,
    mean = mu,
    median = pcts[5],
    sd = sigma,
    min = min(x_finite),
    max = max(x_finite),
    range = max(x_finite) - min(x_finite),
    p01 = pcts[1], p05 = pcts[2], p10 = pcts[3],
    p25 = pcts[4], p75 = pcts[6],
    p90 = pcts[7], p95 = pcts[8], p99 = pcts[9],
    iqr = pcts[6] - pcts[4],
    cv = if (abs(mu) > 1e-10) abs(sigma / mu) else NA_real_,
    skewness = skew,
    kurtosis = kurt,
    n_negative = sum(x_finite < 0),
    pct_negative = 100 * sum(x_finite < 0) / n_valid,
    n_zero = sum(x_finite == 0),
    pct_zero = 100 * sum(x_finite == 0) / n_valid,
    n_unique = uniqueN(x_finite),
    outliers_iqr_n = outliers_iqr$n_low + outliers_iqr$n_high,
    outliers_iqr_pct = outliers_iqr$pct,
    outliers_mad_n = outliers_mad$n_outliers,
    outliers_mad_pct = outliers_mad$pct
  )
}

# ============================================================================
# QC CHECK FUNCTIONS
# ============================================================================

# 1. Data Integrity Checks
check_data_integrity <- function(dt, site_id) {
  issues <- list()
  n_rows <- nrow(dt)
  
  # Check for duplicate shot_numbers
  if ("shot_number" %in% names(dt)) {
    n_dup_shots <- n_rows - uniqueN(dt$shot_number)
    if (n_dup_shots > 0) {
      issues$duplicate_shots <- list(
        severity = "CRITICAL",
        n = n_dup_shots,
        pct = 100 * n_dup_shots / n_rows,
        msg = sprintf("%d duplicate shot_numbers (%.2f%%)", n_dup_shots, 100 * n_dup_shots / n_rows)
      )
    }
    
    # Check shot_number format (should be numeric string)
    invalid_shots <- sum(!grepl("^\\d+$", dt$shot_number))
    if (invalid_shots > 0) {
      issues$invalid_shot_format <- list(
        severity = "WARNING",
        n = invalid_shots,
        msg = sprintf("%d shot_numbers with invalid format", invalid_shots)
      )
    }
  }
  
  # Check for completely duplicate rows
  n_dup_rows <- n_rows - uniqueN(dt)
  if (n_dup_rows > 0) {
    issues$duplicate_rows <- list(
      severity = "CRITICAL",
      n = n_dup_rows,
      pct = 100 * n_dup_rows / n_rows,
      msg = sprintf("%d completely duplicate rows", n_dup_rows)
    )
  }
  
  # Check for constant columns (zero variance)
  numeric_cols <- names(dt)[sapply(dt, is.numeric)]
  constant_cols <- c()
  for (col in numeric_cols) {
    vals <- dt[[col]][is.finite(dt[[col]])]
    if (length(vals) > 0 && length(unique(vals)) == 1) {
      constant_cols <- c(constant_cols, col)
    }
  }
  if (length(constant_cols) > 0) {
    issues$constant_columns <- list(
      severity = "WARNING",
      columns = constant_cols,
      n = length(constant_cols),
      msg = sprintf("%d constant numeric columns: %s", 
                    length(constant_cols), paste(constant_cols, collapse = ", "))
    )
  }
  
  # Check site column consistency
  if ("site" %in% names(dt)) {
    unique_sites <- unique(dt$site)
    if (length(unique_sites) > 1) {
      issues$mixed_sites <- list(
        severity = "CRITICAL",
        sites = unique_sites,
        msg = sprintf("Multiple site IDs in file: %s", paste(unique_sites, collapse = ", "))
      )
    }
    if (!site_id %in% unique_sites) {
      issues$site_mismatch <- list(
        severity = "CRITICAL",
        expected = site_id,
        found = unique_sites,
        msg = sprintf("Expected site %d but found %s", site_id, paste(unique_sites, collapse = ", "))
      )
    }
  }
  
  issues
}

# 2. Value Range Validation
check_value_ranges <- function(dt, site_id) {
  issues <- list()
  flagged <- data.table()
  
  add_flags <- function(idx, flag_name, values) {
    if (length(idx) > 0) {
      flagged <<- rbind(flagged, data.table(
        site = site_id,
        row_idx = idx[1:min(500, length(idx))],
        shot_number = if ("shot_number" %in% names(dt)) dt$shot_number[idx[1:min(500, length(idx))]] else NA,
        flag = flag_name,
        value = values[1:min(500, length(idx))]
      ), fill = TRUE)
    }
  }
  
  # Elevation checks
  for (col in c("dep_dtm_mean", "p3d_dtm_mean", "dep_dtm", "mhrsi_dtm", "elev_low")) {
    if (col %in% names(dt)) {
      x <- dt[[col]]
      idx_low <- which(x < BOUNDS$elev_min & is.finite(x))
      idx_high <- which(x > BOUNDS$elev_max & is.finite(x))
      if (length(idx_low) > 0) {
        issues[[paste0(col, "_below_min")]] <- list(
          severity = "CRITICAL", n = length(idx_low),
          msg = sprintf("%s: %d values below %dm", col, length(idx_low), BOUNDS$elev_min)
        )
        add_flags(idx_low, paste0(col, "_below_min"), x[idx_low])
      }
      if (length(idx_high) > 0) {
        issues[[paste0(col, "_above_max")]] <- list(
          severity = "CRITICAL", n = length(idx_high),
          msg = sprintf("%s: %d values above %dm", col, length(idx_high), BOUNDS$elev_max)
        )
        add_flags(idx_high, paste0(col, "_above_max"), x[idx_high])
      }
    }
  }
  
  # CHM checks
  for (col in c("als_chm_mean", "p3d_chm_mean", "als_chm", "p3d_chm")) {
    if (col %in% names(dt)) {
      x <- dt[[col]]
      idx_low <- which(x < BOUNDS$chm_min & is.finite(x))
      idx_high <- which(x > BOUNDS$chm_max & is.finite(x))
      if (length(idx_low) > 0) {
        issues[[paste0(col, "_excessive_negative")]] <- list(
          severity = "WARNING", n = length(idx_low),
          msg = sprintf("%s: %d values below %dm", col, length(idx_low), BOUNDS$chm_min)
        )
        add_flags(idx_low, paste0(col, "_excessive_negative"), x[idx_low])
      }
      if (length(idx_high) > 0) {
        issues[[paste0(col, "_excessive_positive")]] <- list(
          severity = "WARNING", n = length(idx_high),
          msg = sprintf("%s: %d values above %dm", col, length(idx_high), BOUNDS$chm_max)
        )
        add_flags(idx_high, paste0(col, "_excessive_positive"), x[idx_high])
      }
    }
  }
  
  # RH98 checks
  if ("rh_98" %in% names(dt)) {
    x <- dt$rh_98
    idx_neg <- which(x < BOUNDS$rh98_min & is.finite(x))
    idx_high <- which(x > BOUNDS$rh98_max & is.finite(x))
    if (length(idx_neg) > 0) {
      issues$rh98_negative <- list(
        severity = "WARNING", n = length(idx_neg),
        msg = sprintf("rh_98: %d negative values", length(idx_neg))
      )
      add_flags(idx_neg, "rh98_negative", x[idx_neg])
    }
    if (length(idx_high) > 0) {
      issues$rh98_extreme <- list(
        severity = "WARNING", n = length(idx_high),
        msg = sprintf("rh_98: %d values > %dm", length(idx_high), BOUNDS$rh98_max)
      )
      add_flags(idx_high, "rh98_extreme", x[idx_high])
    }
  }
  
  # Slope checks
  if ("slope_mean" %in% names(dt)) {
    x <- dt$slope_mean
    idx_invalid <- which((x < BOUNDS$slope_min | x > BOUNDS$slope_max) & is.finite(x))
    if (length(idx_invalid) > 0) {
      issues$slope_out_of_range <- list(
        severity = "WARNING", n = length(idx_invalid),
        msg = sprintf("slope_mean: %d values outside [0, 90]", length(idx_invalid))
      )
      add_flags(idx_invalid, "slope_out_of_range", x[idx_invalid])
    }
  }
  
  # Trig function checks (aspect_sin, aspect_cos)
  for (col in c("aspect_sin_mean", "aspect_cos_mean")) {
    if (col %in% names(dt)) {
      x <- dt[[col]]
      idx_invalid <- which((x < BOUNDS$trig_min | x > BOUNDS$trig_max) & is.finite(x))
      if (length(idx_invalid) > 0) {
        issues[[paste0(col, "_out_of_range")]] <- list(
          severity = "WARNING", n = length(idx_invalid),
          msg = sprintf("%s: %d values outside [-1, 1]", col, length(idx_invalid))
        )
        add_flags(idx_invalid, paste0(col, "_out_of_range"), x[idx_invalid])
      }
    }
  }
  
  # Valid fraction checks (should be in [0, 1])
  frac_cols <- names(dt)[grepl("_frac$|_valid_frac$|_ratio$", names(dt))]
  for (col in frac_cols) {
    if (is.numeric(dt[[col]])) {
      x <- dt[[col]]
      idx_invalid <- which((x < 0 | x > 1) & is.finite(x))
      if (length(idx_invalid) > 0) {
        issues[[paste0(col, "_invalid")]] <- list(
          severity = "WARNING", n = length(idx_invalid),
          msg = sprintf("%s: %d values outside [0, 1]", col, length(idx_invalid))
        )
        add_flags(idx_invalid, paste0(col, "_invalid"), x[idx_invalid])
      }
    }
  }
  
  # Cover check
  if ("cover" %in% names(dt)) {
    x <- dt$cover
    idx_invalid <- which((x < BOUNDS$cover_min | x > BOUNDS$cover_max) & is.finite(x))
    if (length(idx_invalid) > 0) {
      issues$cover_out_of_range <- list(
        severity = "WARNING", n = length(idx_invalid),
        msg = sprintf("cover: %d values outside [0, 1]", length(idx_invalid))
      )
      add_flags(idx_invalid, "cover_out_of_range", x[idx_invalid])
    }
  }
  
  list(issues = issues, flagged = flagged)
}

# 3. Cross-Variable Consistency
check_cross_variable_consistency <- function(dt, site_id) {
  issues <- list()
  flagged <- data.table()
  n_rows <- nrow(dt)
  
  add_flags <- function(idx, flag_name, values) {
    if (length(idx) > 0) {
      flagged <<- rbind(flagged, data.table(
        site = site_id,
        row_idx = idx[1:min(200, length(idx))],
        shot_number = if ("shot_number" %in% names(dt)) dt$shot_number[idx[1:min(200, length(idx))]] else NA,
        flag = flag_name,
        value = values[1:min(200, length(idx))]
      ), fill = TRUE)
    }
  }
  
  # CHM vs rh_98 consistency (P3D CHM should correlate with GEDI rh_98)
  if ("p3d_chm_mean" %in% names(dt) && "rh_98" %in% names(dt)) {
    valid_idx <- which(is.finite(dt$p3d_chm_mean) & is.finite(dt$rh_98))
    if (length(valid_idx) > 100) {
      cor_val <- cor(dt$p3d_chm_mean[valid_idx], dt$rh_98[valid_idx], method = "spearman")
      if (!is.na(cor_val) && cor_val < 0.3) {
        issues$low_chm_rh98_correlation <- list(
          severity = "WARNING",
          correlation = cor_val,
          msg = sprintf("Low correlation between p3d_chm_mean and rh_98: %.3f", cor_val)
        )
      }
    }
  }
  
  # DTM reference vs test: extreme differences
  if ("dep_dtm_mean" %in% names(dt) && "p3d_dtm_mean" %in% names(dt)) {
    dtm_diff <- dt$p3d_dtm_mean - dt$dep_dtm_mean
    valid_diff <- dtm_diff[is.finite(dtm_diff)]
    
    if (length(valid_diff) > 0) {
      # Extreme individual differences
      idx_extreme <- which(abs(dtm_diff) > BOUNDS$error_critical & is.finite(dtm_diff))
      if (length(idx_extreme) > 0) {
        issues$dtm_extreme_differences <- list(
          severity = "WARNING",
          n = length(idx_extreme),
          pct = 100 * length(idx_extreme) / length(valid_diff),
          msg = sprintf("DTM: %d footprints with |diff| > %dm (%.2f%%)", 
                        length(idx_extreme), BOUNDS$error_critical,
                        100 * length(idx_extreme) / length(valid_diff))
        )
        add_flags(idx_extreme, "dtm_extreme_diff", dtm_diff[idx_extreme])
      }
    }
  }
  
  # CHM reference vs test: extreme differences
  if ("als_chm_mean" %in% names(dt) && "p3d_chm_mean" %in% names(dt)) {
    chm_diff <- dt$p3d_chm_mean - dt$als_chm_mean
    valid_diff <- chm_diff[is.finite(chm_diff)]
    
    if (length(valid_diff) > 0) {
      idx_extreme <- which(abs(chm_diff) > BOUNDS$error_critical & is.finite(chm_diff))
      if (length(idx_extreme) > 0) {
        issues$chm_extreme_differences <- list(
          severity = "WARNING",
          n = length(idx_extreme),
          pct = 100 * length(idx_extreme) / length(valid_diff),
          msg = sprintf("CHM: %d footprints with |diff| > %dm (%.2f%%)", 
                        length(idx_extreme), BOUNDS$error_critical,
                        100 * length(idx_extreme) / length(valid_diff))
        )
        add_flags(idx_extreme, "chm_extreme_diff", chm_diff[idx_extreme])
      }
    }
  }
  
  # GEDI elev_low vs DTM consistency
  if ("elev_low" %in% names(dt) && "dep_dtm_mean" %in% names(dt)) {
    elev_diff <- dt$elev_low - dt$dep_dtm_mean
    valid_diff <- elev_diff[is.finite(elev_diff)]
    if (length(valid_diff) > 100) {
      # Large systematic offset could indicate geolocation issues
      mean_diff <- mean(valid_diff)
      if (abs(mean_diff) > 10) {
        issues$elev_dtm_offset <- list(
          severity = "WARNING",
          mean_offset = mean_diff,
          msg = sprintf("Mean elev_low - dep_dtm_mean offset: %.2fm (potential geolocation issue)", mean_diff)
        )
      }
    }
  }
  
  # Error correlation (CHM and DTM errors should be somewhat independent)
  if ("error_mean" %in% names(dt) && "dtm_error_mean" %in% names(dt)) {
    valid_idx <- which(is.finite(dt$error_mean) & is.finite(dt$dtm_error_mean))
    if (length(valid_idx) > 100) {
      cor_val <- cor(dt$error_mean[valid_idx], dt$dtm_error_mean[valid_idx])
      if (!is.na(cor_val)) {
        issues$chm_dtm_error_correlation <- list(
          severity = "INFO",
          correlation = cor_val,
          msg = sprintf("CHM-DTM error correlation: %.3f", cor_val)
        )
      }
    }
  }
  
  list(issues = issues, flagged = flagged)
}

# 4. Zonal Statistics Internal Consistency
check_zonal_consistency <- function(dt, site_id) {
  issues <- list()
  flagged <- data.table()
  
  # Check percentile ordering (p25 <= median <= p75 <= p90)
  check_percentile_order <- function(prefix) {
    cols <- paste0(prefix, c("_p25", "_median", "_p75", "_p90"))
    cols_present <- cols[cols %in% names(dt)]
    if (length(cols_present) < 2) return(NULL)
    
    # Get column values
    vals <- lapply(cols_present, function(c) dt[[c]])
    names(vals) <- cols_present
    
    violations <- 0
    for (i in 1:(length(cols_present) - 1)) {
      lower <- vals[[i]]
      upper <- vals[[i + 1]]
      valid_idx <- which(is.finite(lower) & is.finite(upper))
      if (length(valid_idx) > 0) {
        violations <- violations + sum(lower[valid_idx] > upper[valid_idx] + THRESHOLDS$percentile_order_tolerance)
      }
    }
    
    if (violations > 0) {
      list(
        severity = "WARNING",
        n = violations,
        msg = sprintf("%s percentile ordering violated in %d rows", prefix, violations)
      )
    } else NULL
  }
  
  for (prefix in c("dep_dtm", "p3d_dtm", "als_chm", "p3d_chm")) {
    result <- check_percentile_order(prefix)
    if (!is.null(result)) {
      issues[[paste0(prefix, "_percentile_order")]] <- result
    }
  }
  
  # Check mean vs median relationship (extreme skew indicator)
  check_mean_median_skew <- function(mean_col, median_col, name) {
    if (!all(c(mean_col, median_col) %in% names(dt))) return(NULL)
    
    x_mean <- dt[[mean_col]]
    x_median <- dt[[median_col]]
    valid_idx <- which(is.finite(x_mean) & is.finite(x_median))
    
    if (length(valid_idx) < 100) return(NULL)
    
    # Compute relative difference
    denom <- pmax(abs(x_mean[valid_idx]), abs(x_median[valid_idx]), 0.1)
    rel_diff <- abs(x_mean[valid_idx] - x_median[valid_idx]) / denom
    
    extreme_skew <- sum(rel_diff > 0.5)  # 50% relative difference
    if (extreme_skew > 0.01 * length(valid_idx)) {
      list(
        severity = "INFO",
        n = extreme_skew,
        pct = 100 * extreme_skew / length(valid_idx),
        msg = sprintf("%s: %d footprints (%.1f%%) with large mean-median discrepancy (high skew within footprint)",
                      name, extreme_skew, 100 * extreme_skew / length(valid_idx))
      )
    } else NULL
  }
  
  for (prefix in c("dep_dtm", "p3d_dtm", "als_chm", "p3d_chm")) {
    result <- check_mean_median_skew(paste0(prefix, "_mean"), paste0(prefix, "_median"), prefix)
    if (!is.null(result)) {
      issues[[paste0(prefix, "_mean_median_skew")]] <- result
    }
  }
  
  # Check valid_frac consistency (if valid_frac is 0, mean should be NA)
  for (prefix in c("dep_dtm", "p3d_dtm", "als_chm", "p3d_chm", "slope")) {
    frac_col <- paste0(prefix, "_valid_frac")
    mean_col <- paste0(prefix, "_mean")
    if (frac_col %in% names(dt) && mean_col %in% names(dt)) {
      zero_frac_idx <- which(dt[[frac_col]] == 0)
      if (length(zero_frac_idx) > 0) {
        has_value <- sum(!is.na(dt[[mean_col]][zero_frac_idx]))
        if (has_value > 0) {
          issues[[paste0(prefix, "_frac_mean_inconsistency")]] <- list(
            severity = "WARNING",
            n = has_value,
            msg = sprintf("%s: %d rows with valid_frac=0 but non-NA mean", prefix, has_value)
          )
        }
      }
    }
  }
  
  list(issues = issues, flagged = flagged)
}

# 5. Error Distribution Diagnostics
check_error_distributions <- function(dt, site_id) {
  issues <- list()
  distributions <- list()
  
  analyze_error <- function(col_name, display_name) {
    if (!col_name %in% names(dt)) return(NULL)
    x <- dt[[col_name]]
    x_valid <- x[is.finite(x)]
    n <- length(x_valid)
    if (n < 50) return(NULL)
    
    # Basic stats
    mu <- mean(x_valid)
    med <- median(x_valid)
    sigma <- sd(x_valid)
    skew <- if (sigma > 0) mean(((x_valid - mu) / sigma)^3) else NA
    kurt <- if (sigma > 0) mean(((x_valid - mu) / sigma)^4) - 3 else NA
    
    # Normality test (sample if too large)
    shapiro_p <- NA
    if (n >= 3 && n <= 5000) {
      shapiro_p <- tryCatch(shapiro.test(x_valid)$p.value, error = function(e) NA)
    } else if (n > 5000) {
      sample_idx <- sample(n, 5000)
      shapiro_p <- tryCatch(shapiro.test(x_valid[sample_idx])$p.value, error = function(e) NA)
    }
    
    # Bimodality check using dip test (if available)
    bimodal_indicator <- NA
    if (requireNamespace("diptest", quietly = TRUE) && n >= 10) {
      dip_p <- tryCatch(diptest::dip.test(x_valid)$p.value, error = function(e) NA)
      bimodal_indicator <- if (!is.na(dip_p) && dip_p < 0.05) TRUE else FALSE
    }
    
    # RMSE
    rmse <- sqrt(mean(x_valid^2))
    
    # Store distribution info
    dist_info <- list(
      column = col_name,
      n = n,
      mean = mu,
      median = med,
      sd = sigma,
      rmse = rmse,
      skewness = skew,
      kurtosis = kurt,
      shapiro_p = shapiro_p,
      bimodal = bimodal_indicator
    )
    
    # Flag issues
    issue_list <- list()
    
    if (!is.na(skew) && abs(skew) > THRESHOLDS$skew_warn) {
      issue_list$high_skewness <- list(
        severity = "INFO",
        value = skew,
        msg = sprintf("%s: High skewness (%.2f), consider robust methods", display_name, skew)
      )
    }
    
    if (!is.na(kurt) && kurt > THRESHOLDS$kurt_warn) {
      issue_list$heavy_tails <- list(
        severity = "INFO",
        value = kurt,
        msg = sprintf("%s: Heavy tails (kurtosis=%.2f), outliers may impact analysis", display_name, kurt)
      )
    }
    
    if (!is.na(bimodal_indicator) && bimodal_indicator) {
      issue_list$potential_bimodality <- list(
        severity = "WARNING",
        msg = sprintf("%s: Potential bimodality detected", display_name)
      )
    }
    
    list(dist = dist_info, issues = issue_list)
  }
  
  # Analyze error columns
  error_analyses <- list(
    chm_error = analyze_error("error_mean", "CHM error"),
    dtm_error = analyze_error("dtm_error_mean", "DTM error")
  )
  
  for (name in names(error_analyses)) {
    if (!is.null(error_analyses[[name]])) {
      distributions[[name]] <- error_analyses[[name]]$dist
      for (issue_name in names(error_analyses[[name]]$issues)) {
        issues[[paste0(name, "_", issue_name)]] <- error_analyses[[name]]$issues[[issue_name]]
      }
    }
  }
  
  list(issues = issues, distributions = distributions)
}

# 6. Missing Data Analysis
analyze_missing_patterns <- function(dt, site_id) {
  n_rows <- nrow(dt)
  all_cols <- names(dt)
  
  # Column-wise missing
  col_missing <- data.table(
    column = all_cols,
    n_missing = sapply(dt, function(x) sum(is.na(x) | is.nan(x) | (is.numeric(x) & is.infinite(x)))),
    n_total = n_rows
  )
  col_missing[, pct_missing := 100 * n_missing / n_total]
  col_missing[, severity := fcase(
    pct_missing == 0, "OK",
    pct_missing < THRESHOLDS$missing_warn_pct, "LOW",
    pct_missing < THRESHOLDS$missing_critical_pct, "MODERATE",
    pct_missing < 100, "HIGH",
    default = "COMPLETE"
  )]
  setorder(col_missing, -pct_missing)
  
  # Row-wise missing
  missing_per_row <- rowSums(sapply(dt, function(x) is.na(x) | is.nan(x) | (is.numeric(x) & is.infinite(x))))
  
  # Co-occurrence analysis: which columns are missing together?
  # Sample for efficiency
  sample_size <- min(10000, n_rows)
  sample_idx <- if (n_rows > sample_size) sample(n_rows, sample_size) else 1:n_rows
  dt_sample <- dt[sample_idx]
  
  # Create missing indicator matrix for high-missing columns
  high_missing_cols <- col_missing[pct_missing > 10, column]
  if (length(high_missing_cols) >= 2 && length(high_missing_cols) <= 20) {
    missing_matrix <- sapply(high_missing_cols, function(col) {
      is.na(dt_sample[[col]]) | is.nan(dt_sample[[col]]) | 
        (is.numeric(dt_sample[[col]]) & is.infinite(dt_sample[[col]]))
    })
    
    # Compute co-occurrence
    co_occur <- crossprod(missing_matrix)
    diag(co_occur) <- 0
    
    # Find highly co-occurring pairs
    high_cooccur <- which(co_occur > 0.9 * sample_size, arr.ind = TRUE)
    high_cooccur <- high_cooccur[high_cooccur[,1] < high_cooccur[,2], , drop = FALSE]
  } else {
    high_cooccur <- matrix(nrow = 0, ncol = 2)
  }
  
  # Critical column check
  critical_status <- data.table(
    column = c(GEDI_CORE, GEDI_METRICS),
    expected = TRUE
  )
  critical_status[, present := column %in% all_cols]
  critical_status[, n_missing := sapply(column, function(c) {
    if (c %in% all_cols) sum(is.na(dt[[c]])) else NA_integer_
  })]
  critical_status[, pct_missing := 100 * n_missing / n_rows]
  
  list(
    col_missing = col_missing,
    n_complete_rows = sum(missing_per_row == 0),
    pct_complete_rows = 100 * sum(missing_per_row == 0) / n_rows,
    mean_missing_per_row = mean(missing_per_row),
    median_missing_per_row = median(missing_per_row),
    max_missing_per_row = max(missing_per_row),
    critical_status = critical_status,
    n_high_cooccur_pairs = nrow(high_cooccur)
  )
}

# 7. Covariate Quality for Bayesian Modeling
check_covariate_quality <- function(dt, site_id) {
  issues <- list()
  quality <- list()
  
  # Key covariates for modeling
  covariates <- c("slope_mean", "slope_sd", "aspect_sin_mean", "aspect_cos_mean",
                  "rh_98", "cover", "elev_low")
  covariates <- covariates[covariates %in% names(dt)]
  
  if (length(covariates) < 2) return(list(issues = issues, quality = quality))
  
  # Check for sufficient variability
  for (cov in covariates) {
    x <- dt[[cov]][is.finite(dt[[cov]])]
    if (length(x) > 100) {
      cv <- sd(x) / abs(mean(x) + 0.01)
      iqr_val <- IQR(x)
      range_val <- diff(range(x))
      
      quality[[cov]] <- list(
        n_valid = length(x),
        mean = mean(x),
        sd = sd(x),
        cv = cv,
        iqr = iqr_val,
        range = range_val
      )
      
      # Flag low variability
      if (cv < 0.05 && range_val < 1) {
        issues[[paste0(cov, "_low_variability")]] <- list(
          severity = "WARNING",
          cv = cv,
          range = range_val,
          msg = sprintf("%s: Low variability (CV=%.3f, range=%.2f) - may not be informative", 
                        cov, cv, range_val)
        )
      }
    }
  }
  
  # Check for multicollinearity among numeric covariates
  numeric_covs <- covariates[sapply(covariates, function(c) is.numeric(dt[[c]]))]
  if (length(numeric_covs) >= 2) {
    cov_data <- na.omit(dt[, ..numeric_covs])
    if (nrow(cov_data) > 100) {
      cor_matrix <- cor(cov_data)
      
      # Find high correlations (exclude diagonal)
      high_cor_pairs <- which(abs(cor_matrix) > 0.85 & upper.tri(cor_matrix), arr.ind = TRUE)
      if (nrow(high_cor_pairs) > 0) {
        pair_names <- apply(high_cor_pairs, 1, function(idx) {
          paste(numeric_covs[idx[1]], "-", numeric_covs[idx[2]], 
                sprintf("(r=%.2f)", cor_matrix[idx[1], idx[2]]))
        })
        issues$high_covariate_correlation <- list(
          severity = "INFO",
          pairs = pair_names,
          msg = sprintf("High correlations among covariates: %s", paste(pair_names, collapse = "; "))
        )
      }
    }
  }
  
  list(issues = issues, quality = quality)
}

# ============================================================================
# MAIN QC FUNCTION PER SITE
# ============================================================================

qc_site_comprehensive <- function(site_id) {
  csv_path <- find_enriched_csv(site_id)
  
  if (is.na(csv_path)) {
    return(list(
      site = site_id,
      status = "NOT_FOUND",
      all_issues = list(file_not_found = list(severity = "CRITICAL", 
                                               msg = "Enriched CSV file not found")),
      summary = data.table(site = site_id, status = "NOT_FOUND", n_rows = 0L)
    ))
  }
  
  # Read data
  dt <- tryCatch(
    fread(csv_path, colClasses = c(shot_number = "character")),
    error = function(e) NULL
  )
  
  if (is.null(dt) || nrow(dt) == 0) {
    return(list(
      site = site_id,
      status = if (is.null(dt)) "READ_ERROR" else "EMPTY",
      all_issues = list(data_read_error = list(severity = "CRITICAL", 
                                                msg = "Failed to read or empty file")),
      summary = data.table(site = site_id, status = "READ_ERROR", n_rows = 0L)
    ))
  }
  
  n_rows <- nrow(dt)
  n_cols <- ncol(dt)
  all_issues <- list()
  all_flagged <- data.table()
  
  # Run all QC checks
  # 1. Data Integrity
  integrity <- check_data_integrity(dt, site_id)
  all_issues <- c(all_issues, integrity)
  
  # 2. Value Ranges
  ranges <- check_value_ranges(dt, site_id)
  all_issues <- c(all_issues, ranges$issues)
  all_flagged <- rbind(all_flagged, ranges$flagged, fill = TRUE)
  
  # 3. Cross-Variable Consistency
  consistency <- check_cross_variable_consistency(dt, site_id)
  all_issues <- c(all_issues, consistency$issues)
  all_flagged <- rbind(all_flagged, consistency$flagged, fill = TRUE)
  
  # 4. Zonal Consistency
  zonal <- check_zonal_consistency(dt, site_id)
  all_issues <- c(all_issues, zonal$issues)
  all_flagged <- rbind(all_flagged, zonal$flagged, fill = TRUE)
  
  # 5. Error Distributions
  errors <- check_error_distributions(dt, site_id)
  all_issues <- c(all_issues, errors$issues)
  
  # 6. Missing Data
  missing <- analyze_missing_patterns(dt, site_id)
  
  # 7. Covariate Quality
  covariates <- check_covariate_quality(dt, site_id)
  all_issues <- c(all_issues, covariates$issues)
  
  # Compute column statistics for all numeric columns
  numeric_cols <- names(dt)[sapply(dt, is.numeric)]
  col_stats <- rbindlist(lapply(numeric_cols, function(col) {
    compute_col_stats(dt[[col]], col)
  }), fill = TRUE)
  col_stats[, site := site_id]
  
  # Compute bias metrics
  chm_bias <- dtm_bias <- list(mean = NA_real_, median = NA_real_, sd = NA_real_, n = 0L)
  
  if ("als_chm_mean" %in% names(dt) && "p3d_chm_mean" %in% names(dt)) {
    chm_diff <- dt$p3d_chm_mean - dt$als_chm_mean
    valid <- chm_diff[is.finite(chm_diff)]
    if (length(valid) > 0) {
      chm_bias <- list(mean = mean(valid), median = median(valid), 
                       sd = sd(valid), n = length(valid))
    }
  }
  
  if ("dep_dtm_mean" %in% names(dt) && "p3d_dtm_mean" %in% names(dt)) {
    dtm_diff <- dt$p3d_dtm_mean - dt$dep_dtm_mean
    valid <- dtm_diff[is.finite(dtm_diff)]
    if (length(valid) > 0) {
      dtm_bias <- list(mean = mean(valid), median = median(valid), 
                       sd = sd(valid), n = length(valid))
    }
  }
  
  # Count issues by severity
  n_critical <- sum(sapply(all_issues, function(x) x$severity == "CRITICAL"))
  n_warning <- sum(sapply(all_issues, function(x) x$severity == "WARNING"))
  n_info <- sum(sapply(all_issues, function(x) x$severity == "INFO"))
  
  # Compute quality score (0-100)
  # Deduct points for issues
  quality_score <- 100
  quality_score <- quality_score - (n_critical * 15)
  quality_score <- quality_score - (n_warning * 3)
  quality_score <- quality_score - (n_info * 0.5)
  # Deduct for missing data
  quality_score <- quality_score - max(0, (100 - missing$pct_complete_rows) * 0.2)
  quality_score <- max(0, min(100, quality_score))
  
  # Create summary
  site_summary <- data.table(
    site = site_id,
    status = "OK",
    n_rows = n_rows,
    n_cols = n_cols,
    quality_score = round(quality_score, 1),
    n_issues_critical = n_critical,
    n_issues_warning = n_warning,
    n_issues_info = n_info,
    n_flagged_rows = nrow(all_flagged),
    n_complete_rows = missing$n_complete_rows,
    pct_complete_rows = missing$pct_complete_rows,
    mean_missing_per_row = missing$mean_missing_per_row,
    chm_bias_mean = chm_bias$mean,
    chm_bias_median = chm_bias$median,
    chm_bias_sd = chm_bias$sd,
    chm_n_valid = chm_bias$n,
    dtm_bias_mean = dtm_bias$mean,
    dtm_bias_median = dtm_bias$median,
    dtm_bias_sd = dtm_bias$sd,
    dtm_n_valid = dtm_bias$n
  )
  
  # Add error distribution stats if available
  if (!is.null(errors$distributions$chm_error)) {
    site_summary[, `:=`(
      chm_error_mean = errors$distributions$chm_error$mean,
      chm_error_sd = errors$distributions$chm_error$sd,
      chm_error_rmse = errors$distributions$chm_error$rmse,
      chm_error_skew = errors$distributions$chm_error$skewness,
      chm_error_kurt = errors$distributions$chm_error$kurtosis
    )]
  }
  if (!is.null(errors$distributions$dtm_error)) {
    site_summary[, `:=`(
      dtm_error_mean = errors$distributions$dtm_error$mean,
      dtm_error_sd = errors$distributions$dtm_error$sd,
      dtm_error_rmse = errors$distributions$dtm_error$rmse,
      dtm_error_skew = errors$distributions$dtm_error$skewness,
      dtm_error_kurt = errors$distributions$dtm_error$kurtosis
    )]
  }
  
  list(
    site = site_id,
    status = "OK",
    n_rows = n_rows,
    all_issues = all_issues,
    flagged_rows = all_flagged,
    col_stats = col_stats,
    missing = missing,
    error_distributions = errors$distributions,
    covariate_quality = covariates$quality,
    summary = site_summary
  )
}

# ============================================================================
# CROSS-SITE ANOMALY DETECTION
# ============================================================================

detect_site_anomalies <- function(summaries) {
  anomalies <- list()
  
  numeric_cols <- names(summaries)[sapply(summaries, is.numeric)]
  # Focus on key metrics
  key_metrics <- c("chm_bias_mean", "dtm_bias_mean", "chm_error_rmse", "dtm_error_rmse",
                   "pct_complete_rows", "quality_score")
  key_metrics <- intersect(key_metrics, numeric_cols)
  
  for (metric in key_metrics) {
    x <- summaries[[metric]]
    valid_idx <- which(is.finite(x))
    if (length(valid_idx) < 3) next
    
    x_valid <- x[valid_idx]
    mu <- mean(x_valid)
    sigma <- sd(x_valid)
    
    if (sigma > 0) {
      z_scores <- (x_valid - mu) / sigma
      
      # Find outlier sites
      for (i in seq_along(valid_idx)) {
        site_idx <- valid_idx[i]
        z <- z_scores[i]
        
        if (abs(z) >= THRESHOLDS$site_zscore_critical) {
          anomalies[[paste0("site_", summaries$site[site_idx], "_", metric)]] <- list(
            site = summaries$site[site_idx],
            metric = metric,
            value = x[site_idx],
            z_score = z,
            severity = "CRITICAL",
            msg = sprintf("Site %d: %s = %.3f (z = %.2f, %s extreme)", 
                          summaries$site[site_idx], metric, x[site_idx], z,
                          if (z > 0) "high" else "low")
          )
        } else if (abs(z) >= THRESHOLDS$site_zscore_warn) {
          anomalies[[paste0("site_", summaries$site[site_idx], "_", metric)]] <- list(
            site = summaries$site[site_idx],
            metric = metric,
            value = x[site_idx],
            z_score = z,
            severity = "WARNING",
            msg = sprintf("Site %d: %s = %.3f (z = %.2f, %s outlier)", 
                          summaries$site[site_idx], metric, x[site_idx], z,
                          if (z > 0) "high" else "low")
          )
        }
      }
    }
  }
  
  anomalies
}

# ============================================================================
# REPORT GENERATION
# ============================================================================

write_report <- function(results, summaries, site_anomalies, filepath) {
  con <- file(filepath, "w")
  on.exit(close(con))
  
  write_section <- function(title) {
    cat("\n", file = con)
    cat(paste(rep("=", 80), collapse = ""), "\n", file = con)
    cat(title, "\n", file = con)
    cat(paste(rep("=", 80), collapse = ""), "\n", file = con)
  }
  
  write_subsection <- function(title) {
    cat("\n", file = con)
    cat(paste(rep("-", 60), collapse = ""), "\n", file = con)
    cat(title, "\n", file = con)
    cat(paste(rep("-", 60), collapse = ""), "\n", file = con)
  }
  
  # Header
  cat("╔══════════════════════════════════════════════════════════════════════════════╗\n", file = con)
  cat("║         COMPREHENSIVE QC REPORT - ENRICHED GEDI SITE DATA                    ║\n", file = con)
  cat("║                    Zonal Statistics Validation                               ║\n", file = con)
  cat("╚══════════════════════════════════════════════════════════════════════════════╝\n", file = con)
  cat(sprintf("\nGenerated: %s\n", Sys.time()), file = con)
  cat(sprintf("Sites analyzed: %d\n", nrow(summaries)), file = con)
  cat(sprintf("Total observations: %s\n", format(sum(summaries$n_rows), big.mark = ",")), file = con)
  
  # ==================== EXECUTIVE SUMMARY ====================
  write_section("EXECUTIVE SUMMARY")
  
  # Overall statistics
  cat("\n┌─ OVERALL STATISTICS ─────────────────────────────────────────────────────────\n", file = con)
  cat(sprintf("│ Total sites:              %d\n", nrow(summaries)), file = con)
  cat(sprintf("│ Total observations:       %s\n", format(sum(summaries$n_rows), big.mark = ",")), file = con)
  cat(sprintf("│ Mean quality score:       %.1f / 100\n", mean(summaries$quality_score, na.rm = TRUE)), file = con)
  cat(sprintf("│ Mean complete rows:       %.1f%%\n", mean(summaries$pct_complete_rows, na.rm = TRUE)), file = con)
  cat("└──────────────────────────────────────────────────────────────────────────────\n", file = con)
  
  # Issue summary
  total_critical <- sum(summaries$n_issues_critical, na.rm = TRUE)
  total_warning <- sum(summaries$n_issues_warning, na.rm = TRUE)
  total_info <- sum(summaries$n_issues_info, na.rm = TRUE)
  
  cat("\n┌─ ISSUE SUMMARY ──────────────────────────────────────────────────────────────\n", file = con)
  cat(sprintf("│ CRITICAL issues:          %d\n", total_critical), file = con)
  cat(sprintf("│ WARNING issues:           %d\n", total_warning), file = con)
  cat(sprintf("│ INFO issues:              %d\n", total_info), file = con)
  cat(sprintf("│ Site-level anomalies:     %d\n", length(site_anomalies)), file = con)
  cat("└──────────────────────────────────────────────────────────────────────────────\n", file = con)
  
  # Bias summary
  cat("\n┌─ BIAS SUMMARY (P3D - Reference) ─────────────────────────────────────────────\n", file = con)
  cat(sprintf("│ CHM Bias (vs ALS):        Mean = %.3f m, SD = %.3f m\n", 
              mean(summaries$chm_bias_mean, na.rm = TRUE),
              mean(summaries$chm_bias_sd, na.rm = TRUE)), file = con)
  cat(sprintf("│ DTM Bias (vs 3DEP):       Mean = %.3f m, SD = %.3f m\n", 
              mean(summaries$dtm_bias_mean, na.rm = TRUE),
              mean(summaries$dtm_bias_sd, na.rm = TRUE)), file = con)
  cat("└──────────────────────────────────────────────────────────────────────────────\n", file = con)
  
  # ==================== SITE QUALITY RANKING ====================
  write_section("SITE QUALITY RANKING")
  
  summaries_sorted <- summaries[order(-quality_score)]
  cat(sprintf("\n%-6s %8s %12s %10s %10s %12s %12s\n",
              "Site", "Score", "N Rows", "Complete%", "Critical", "CHM Bias", "DTM Bias"), file = con)
  cat(paste(rep("-", 80), collapse = ""), "\n", file = con)
  
  for (i in seq_len(nrow(summaries_sorted))) {
    s <- summaries_sorted[i]
    status <- if (s$quality_score >= 80) "✓" else if (s$quality_score >= 50) "~" else "✗"
    cat(sprintf("%s %-4d %8.1f %12s %9.1f%% %10d %12.3f %12.3f\n",
                status, s$site, s$quality_score, format(s$n_rows, big.mark = ","),
                s$pct_complete_rows, s$n_issues_critical,
                ifelse(is.na(s$chm_bias_mean), NA, s$chm_bias_mean),
                ifelse(is.na(s$dtm_bias_mean), NA, s$dtm_bias_mean)), file = con)
  }
  
  # ==================== CRITICAL ISSUES BY SITE ====================
  write_section("DETAILED ISSUES BY SITE")
  
  for (site_name in names(results)) {
    r <- results[[site_name]]
    if (is.null(r$all_issues) || length(r$all_issues) == 0) next
    
    # Filter to show warnings and critical only
    important_issues <- Filter(function(x) x$severity %in% c("CRITICAL", "WARNING"), r$all_issues)
    if (length(important_issues) == 0) next
    
    write_subsection(sprintf("SITE %d", r$site))
    
    for (issue_name in names(important_issues)) {
      issue <- important_issues[[issue_name]]
      icon <- if (issue$severity == "CRITICAL") "✗" else "⚠"
      cat(sprintf("  %s [%s] %s\n", icon, issue$severity, issue$msg), file = con)
    }
  }
  
  # ==================== SITE ANOMALIES ====================
  if (length(site_anomalies) > 0) {
    write_section("CROSS-SITE ANOMALIES")
    
    cat("\nSites with metrics significantly different from the population:\n\n", file = con)
    
    for (anom_name in names(site_anomalies)) {
      anom <- site_anomalies[[anom_name]]
      icon <- if (anom$severity == "CRITICAL") "✗" else "⚠"
      cat(sprintf("  %s %s\n", icon, anom$msg), file = con)
    }
  }
  
  # ==================== DATA COMPLETENESS ====================
  write_section("DATA COMPLETENESS BY SITE")
  
  cat(sprintf("\n%-6s %12s %12s %15s %15s %15s\n",
              "Site", "N Rows", "Complete", "Complete %", "Mean Miss/Row", "DTM Valid %"), file = con)
  cat(paste(rep("-", 85), collapse = ""), "\n", file = con)
  
  for (i in seq_len(nrow(summaries))) {
    s <- summaries[i]
    dtm_valid_pct <- if (!is.na(s$dtm_n_valid)) 100 * s$dtm_n_valid / s$n_rows else NA
    cat(sprintf("%-6d %12s %12s %14.1f%% %15.1f %14.1f%%\n",
                s$site, format(s$n_rows, big.mark = ","),
                format(s$n_complete_rows, big.mark = ","),
                s$pct_complete_rows, s$mean_missing_per_row,
                ifelse(is.na(dtm_valid_pct), NA, dtm_valid_pct)), file = con)
  }
  
  # ==================== ERROR STATISTICS ====================
  write_section("ERROR DISTRIBUTION STATISTICS")
  
  cat("\n┌─ CHM Error (P3D - ALS) ──────────────────────────────────────────────────────\n", file = con)
  cat(sprintf("%-6s %10s %10s %10s %10s %10s\n", 
              "Site", "Mean", "SD", "RMSE", "Skewness", "Kurtosis"), file = con)
  cat(paste(rep("-", 60), collapse = ""), "\n", file = con)
  
  for (i in seq_len(nrow(summaries))) {
    s <- summaries[i]
    if ("chm_error_mean" %in% names(s) && !is.na(s$chm_error_mean)) {
      cat(sprintf("%-6d %10.3f %10.3f %10.3f %10.2f %10.2f\n",
                  s$site, s$chm_error_mean, s$chm_error_sd, s$chm_error_rmse,
                  ifelse(is.na(s$chm_error_skew), NA, s$chm_error_skew),
                  ifelse(is.na(s$chm_error_kurt), NA, s$chm_error_kurt)), file = con)
    } else {
      cat(sprintf("%-6d %10s %10s %10s %10s %10s\n", s$site, "NA", "NA", "NA", "NA", "NA"), file = con)
    }
  }
  
  cat("\n┌─ DTM Error (P3D - 3DEP) ─────────────────────────────────────────────────────\n", file = con)
  cat(sprintf("%-6s %10s %10s %10s %10s %10s\n", 
              "Site", "Mean", "SD", "RMSE", "Skewness", "Kurtosis"), file = con)
  cat(paste(rep("-", 60), collapse = ""), "\n", file = con)
  
  for (i in seq_len(nrow(summaries))) {
    s <- summaries[i]
    if ("dtm_error_mean" %in% names(s) && !is.na(s$dtm_error_mean)) {
      cat(sprintf("%-6d %10.3f %10.3f %10.3f %10.2f %10.2f\n",
                  s$site, s$dtm_error_mean, s$dtm_error_sd, s$dtm_error_rmse,
                  ifelse(is.na(s$dtm_error_skew), NA, s$dtm_error_skew),
                  ifelse(is.na(s$dtm_error_kurt), NA, s$dtm_error_kurt)), file = con)
    } else {
      cat(sprintf("%-6d %10s %10s %10s %10s %10s\n", s$site, "NA", "NA", "NA", "NA", "NA"), file = con)
    }
  }
  
  # ==================== RECOMMENDATIONS ====================
  write_section("RECOMMENDATIONS")
  
  cat("\nBased on this comprehensive QC analysis:\n\n", file = con)
  
  # Sites to consider excluding
  low_quality_sites <- summaries[quality_score < 50, site]
  if (length(low_quality_sites) > 0) {
    cat(sprintf("⚠ SITES TO REVIEW FOR EXCLUSION (quality score < 50):\n"), file = con)
    cat(sprintf("   Sites: %s\n\n", paste(low_quality_sites, collapse = ", ")), file = con)
  }
  
  # Sites with critical issues
  critical_sites <- summaries[n_issues_critical > 0, site]
  if (length(critical_sites) > 0) {
    cat(sprintf("✗ SITES WITH CRITICAL ISSUES:\n"), file = con)
    cat(sprintf("   Sites: %s\n", paste(critical_sites, collapse = ", ")), file = con)
    cat("   Review the detailed issues above before including in analysis.\n\n", file = con)
  }
  
  # High bias sites
  high_chm_bias <- summaries[!is.na(chm_bias_mean) & abs(chm_bias_mean) > THRESHOLDS$bias_critical, site]
  high_dtm_bias <- summaries[!is.na(dtm_bias_mean) & abs(dtm_bias_mean) > THRESHOLDS$bias_critical, site]
  
  if (length(high_chm_bias) > 0) {
    cat(sprintf("⚠ SITES WITH HIGH CHM BIAS (|bias| > %dm):\n", THRESHOLDS$bias_critical), file = con)
    cat(sprintf("   Sites: %s\n", paste(high_chm_bias, collapse = ", ")), file = con)
    cat("   Consider site-specific bias correction or exclusion.\n\n", file = con)
  }
  
  if (length(high_dtm_bias) > 0) {
    cat(sprintf("⚠ SITES WITH HIGH DTM BIAS (|bias| > %dm):\n", THRESHOLDS$bias_critical), file = con)
    cat(sprintf("   Sites: %s\n", paste(high_dtm_bias, collapse = ", ")), file = con)
    cat("   Consider site-specific bias correction or exclusion.\n\n", file = con)
  }
  
  # Low data completeness
  low_complete <- summaries[pct_complete_rows < 50, site]
  if (length(low_complete) > 0) {
    cat(sprintf("⚠ SITES WITH LOW DATA COMPLETENESS (<50%% complete rows):\n"), file = con)
    cat(sprintf("   Sites: %s\n", paste(low_complete, collapse = ", ")), file = con)
    cat("   Missing data may bias results; consider imputation or exclusion.\n\n", file = con)
  }
  
  # Good sites
  good_sites <- summaries[quality_score >= 80 & n_issues_critical == 0, site]
  if (length(good_sites) > 0) {
    cat(sprintf("✓ HIGH QUALITY SITES (score >= 80, no critical issues):\n"), file = con)
    cat(sprintf("   Sites: %s\n\n", paste(good_sites, collapse = ", ")), file = con)
  }
  
  # Bayesian modeling notes
  cat("─────────────────────────────────────────────────────────────────────────────────\n", file = con)
  cat("NOTES FOR BAYESIAN ANALYSIS:\n\n", file = con)
  cat("1. Error distributions: Check skewness and kurtosis values above.\n", file = con)
  cat("   High skewness (>2) or kurtosis (>7) suggest using robust likelihoods\n", file = con)
  cat("   (e.g., Student-t instead of Gaussian).\n\n", file = con)
  cat("2. Site-level random effects: Consider hierarchical structure given\n", file = con)
  cat("   the between-site variability in bias and error metrics.\n\n", file = con)
  cat("3. Missing data: The 3DEP DTM has limited coverage (~43% in some sites).\n", file = con)
  cat("   Consider modeling this missingness if non-random.\n\n", file = con)
  cat("4. Covariates: Check for multicollinearity warnings above.\n", file = con)
  cat("   Consider centering/scaling predictors for stable MCMC.\n\n", file = con)
  
  cat("\n", file = con)
  cat(paste(rep("=", 80), collapse = ""), "\n", file = con)
  cat("END OF COMPREHENSIVE QC REPORT\n", file = con)
  cat(paste(rep("=", 80), collapse = ""), "\n", file = con)
}

# ============================================================================
# MAIN EXECUTION
# ============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════════════╗\n")
cat("║     COMPREHENSIVE QC ANALYSIS - ENRICHED GEDI SITE DATA                      ║\n")
cat("╚══════════════════════════════════════════════════════════════════════════════╝\n")
cat(sprintf("Processing %d sites\n", length(SITES)))
cat(sprintf("Timestamp: %s\n", TIMESTAMP))
cat(sprintf("Workers: %d\n", determine_workers()))
cat("────────────────────────────────────────────────────────────────────────────────\n")

# Setup parallel processing
n_workers <- determine_workers()
plan(multisession, workers = n_workers)

# Process all sites
cat("\nProcessing sites:\n")
results <- future_lapply(SITES, function(s) {
  res <- qc_site_comprehensive(s)
  cat(sprintf("  Site %02d: %s (n=%s, score=%.0f)\n", 
              s, res$status,
              if (!is.null(res$n_rows)) format(res$n_rows, big.mark = ",") else "NA",
              if (!is.null(res$summary$quality_score)) res$summary$quality_score else 0))
  res
}, future.seed = TRUE)

names(results) <- sprintf("site_%02d", SITES)
plan(sequential)

# Extract summaries
summaries <- rbindlist(lapply(results, function(r) r$summary), fill = TRUE)

# Detect cross-site anomalies
site_anomalies <- detect_site_anomalies(summaries)

# Collect all column statistics
all_col_stats <- rbindlist(lapply(results, function(r) {
  if (!is.null(r$col_stats)) r$col_stats else NULL
}), fill = TRUE)

# Collect all flagged rows
all_flagged <- rbindlist(lapply(results, function(r) {
  if (!is.null(r$flagged_rows) && nrow(r$flagged_rows) > 0) r$flagged_rows else NULL
}), fill = TRUE)

# Collect all issues for CSV export
all_issues_dt <- rbindlist(lapply(names(results), function(site_name) {
  r <- results[[site_name]]
  if (is.null(r$all_issues) || length(r$all_issues) == 0) return(NULL)
  rbindlist(lapply(names(r$all_issues), function(issue_name) {
    issue <- r$all_issues[[issue_name]]
    data.table(
      site = r$site,
      issue_name = issue_name,
      severity = issue$severity,
      message = issue$msg
    )
  }), fill = TRUE)
}), fill = TRUE)

# Add site anomalies to issues
if (length(site_anomalies) > 0) {
  anom_dt <- rbindlist(lapply(names(site_anomalies), function(anom_name) {
    anom <- site_anomalies[[anom_name]]
    data.table(
      site = anom$site,
      issue_name = anom_name,
      severity = anom$severity,
      message = anom$msg
    )
  }), fill = TRUE)
  all_issues_dt <- rbind(all_issues_dt, anom_dt, fill = TRUE)
}

# Generate recommended exclusions
exclusions <- data.table()

# Low quality sites
low_quality <- summaries[quality_score < 50]
if (nrow(low_quality) > 0) {
  exclusions <- rbind(exclusions, data.table(
    site = low_quality$site,
    reason = "Low quality score",
    value = low_quality$quality_score,
    recommendation = "REVIEW"
  ))
}

# Critical issues
critical <- summaries[n_issues_critical > 0]
if (nrow(critical) > 0) {
  exclusions <- rbind(exclusions, data.table(
    site = critical$site,
    reason = "Critical issues present",
    value = critical$n_issues_critical,
    recommendation = "REVIEW"
  ))
}

# High bias
high_bias <- summaries[abs(chm_bias_mean) > THRESHOLDS$bias_critical | 
                        abs(dtm_bias_mean) > THRESHOLDS$bias_critical]
if (nrow(high_bias) > 0) {
  exclusions <- rbind(exclusions, data.table(
    site = high_bias$site,
    reason = "High bias",
    value = pmax(abs(high_bias$chm_bias_mean), abs(high_bias$dtm_bias_mean), na.rm = TRUE),
    recommendation = "REVIEW"
  ))
}

# Write outputs
cat("\n────────────────────────────────────────────────────────────────────────────────\n")
cat("Writing outputs...\n")

# Site summary CSV
fwrite(summaries, SUMMARY_CSV)
cat(sprintf("  Site summary:      %s\n", SUMMARY_CSV))

# Variable statistics CSV
if (nrow(all_col_stats) > 0) {
  fwrite(all_col_stats, VARIABLE_CSV)
  cat(sprintf("  Variable stats:    %s\n", VARIABLE_CSV))
}

# Issues CSV
if (nrow(all_issues_dt) > 0) {
  fwrite(all_issues_dt, ISSUES_CSV)
  cat(sprintf("  All issues:        %s\n", ISSUES_CSV))
}

# Flagged observations CSV
if (nrow(all_flagged) > 0) {
  fwrite(all_flagged, FLAGGED_CSV)
  cat(sprintf("  Flagged obs:       %s\n", FLAGGED_CSV))
}

# Exclusions CSV
if (nrow(exclusions) > 0) {
  exclusions <- unique(exclusions)
  fwrite(exclusions, EXCLUSIONS_CSV)
  cat(sprintf("  Exclusions:        %s\n", EXCLUSIONS_CSV))
}

# Write detailed report
write_report(results, summaries, site_anomalies, REPORT_FILE)
cat(sprintf("  Full report:       %s\n", REPORT_FILE))

# ============================================================================
# CONSOLE SUMMARY
# ============================================================================

cat("\n════════════════════════════════════════════════════════════════════════════════\n")
cat("QC ANALYSIS COMPLETE\n")
cat("════════════════════════════════════════════════════════════════════════════════\n")

cat("\nKEY FINDINGS:\n")
cat(sprintf("  Total sites:          %d\n", nrow(summaries)))
cat(sprintf("  Total observations:   %s\n", format(sum(summaries$n_rows), big.mark = ",")))
cat(sprintf("  Mean quality score:   %.1f / 100\n", mean(summaries$quality_score, na.rm = TRUE)))
cat(sprintf("  Sites with critical:  %d\n", sum(summaries$n_issues_critical > 0)))

cat("\nBIAS (P3D - Reference):\n")
cat(sprintf("  CHM (vs ALS):   %.3f ± %.3f m\n", 
            mean(summaries$chm_bias_mean, na.rm = TRUE),
            sd(summaries$chm_bias_mean, na.rm = TRUE)))
cat(sprintf("  DTM (vs 3DEP):  %.3f ± %.3f m\n", 
            mean(summaries$dtm_bias_mean, na.rm = TRUE),
            sd(summaries$dtm_bias_mean, na.rm = TRUE)))

if (nrow(exclusions) > 0) {
  cat(sprintf("\n⚠ %d site(s) flagged for review - see %s\n", 
              uniqueN(exclusions$site), basename(EXCLUSIONS_CSV)))
}

cat("\n════════════════════════════════════════════════════════════════════════════════\n")
