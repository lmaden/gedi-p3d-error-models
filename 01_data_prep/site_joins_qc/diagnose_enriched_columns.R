#!/usr/bin/env Rscript
# diagnose_enriched_columns.R
# 
# Comprehensive diagnostic for enriched GEDI site files
# Adapted for zonal sampling workflow (footprint-level statistics)
#
# This script:
#   1. Scans all available site files
#   2. Groups columns by source/product type
#   3. Identifies zonal statistics computed (mean, std, count, etc.)
#   4. Checks cross-site consistency
#   5. Reports data completeness

suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================================
# CONFIGURATION
# ============================================================================

ENRICHED_DIR <- "/gpfs/data1/vclgp/lmaden/chpt1/data/enriched_by_site"
SAMPLE_ROWS <- 1000  # Rows to sample for data completeness checks

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

find_all_site_files <- function(dir) {
  # Find all enriched CSV files in directory
  if (!dir.exists(dir)) {
    return(data.table(site_id = integer(), path = character(), exists = logical()))
  }
  
  files <- list.files(dir, pattern = "site_.*_enriched\\.csv(\\.gz)?$", full.names = TRUE)
  
  if (length(files) == 0) {
    return(data.table(site_id = integer(), path = character(), exists = logical()))
  }
  
  # Extract site IDs from filenames
  site_ids <- as.integer(gsub(".*site_0*(\\d+)_enriched.*", "\\1", basename(files)))
  
  data.table(
    site_id = site_ids,
    path = files,
    exists = TRUE
  )[order(site_id)]
}

categorize_columns <- function(col_names) {
  # Categorize columns by source/product type
  categories <- list()
  
  # GEDI core attributes
  gedi_patterns <- c("^shot_", "^beam", "^lat", "^lon", "^elev", "^rh", "^quality",
                     "^sensitivity", "^degrade", "^solar", "^land_cover", "^pft",
                     "^region", "^urban", "^water", "^snow", "^night", "^granule",
                     "^orbit", "^track", "^delta_time", "^selected_algorithm")
  categories$gedi_core <- col_names[Reduce(`|`, lapply(gedi_patterns, function(p) grepl(p, col_names, ignore.case = TRUE)))]
  
  # 3DEP DTM columns
  categories$dep_dtm <- col_names[grepl("dep.*dtm|dtm.*dep|3dep|dep_", col_names, ignore.case = TRUE)]
  
  # P3D / MHRSI DTM columns
  categories$p3d_dtm <- col_names[grepl("p3d.*dtm|mhrsi.*dtm|dtm.*(p3d|mhrsi)", col_names, ignore.case = TRUE)]
  
  # ALS CHM columns
  categories$als_chm <- col_names[grepl("als.*chm|chm.*als|als_", col_names, ignore.case = TRUE)]
  
  # P3D / MHRSI CHM columns
  categories$p3d_chm <- col_names[grepl("p3d.*chm|mhrsi.*chm|chm.*(p3d|mhrsi)", col_names, ignore.case = TRUE)]
  
  # Error/difference columns
  categories$error <- col_names[grepl("error|diff|bias|residual", col_names, ignore.case = TRUE)]
  
  # Slope/terrain columns
  categories$terrain <- col_names[grepl("slope|aspect|tpi|tri|rough|curv|terrain", col_names, ignore.case = TRUE)]
  
  # Canopy metrics (non-CHM)
  categories$canopy <- col_names[grepl("canopy|cover|fhd|pai|agbd|biomass|height|rh\\d", col_names, ignore.case = TRUE)]
  
  # Site/location metadata
  categories$metadata <- col_names[grepl("^site|^cell|^tile|^file|^path|^src|^date|^time|^year|^month|^day", col_names, ignore.case = TRUE)]
  
  # Identify uncategorized columns
  all_categorized <- unique(unlist(categories))
  categories$other <- setdiff(col_names, all_categorized)
  
  categories
}

extract_zonal_stats <- function(col_names) {
  # Identify what zonal statistics were computed
  stat_suffixes <- c("mean", "median", "std", "sd", "min", "max", "sum", 
                     "count", "range", "p05", "p10", "p25", "p50", "p75", "p90", "p95",
                     "q1", "q3", "iqr", "var", "cv", "n", "valid")
  
  detected <- list()
  for (stat in stat_suffixes) {
    pattern <- sprintf("_%s$|_%s_|\\b%s\\b", stat, stat, stat)
    matching <- col_names[grepl(pattern, col_names, ignore.case = TRUE)]
    if (length(matching) > 0) {
      detected[[stat]] <- matching
    }
  }
  detected
}

get_column_summary <- function(dt, cols, max_cols = 50) {
  # Get summary statistics for specified columns
  cols <- intersect(cols, names(dt))
  if (length(cols) == 0) return(NULL)
  
  # Limit columns if too many
  if (length(cols) > max_cols) {
    cols <- cols[1:max_cols]
  }
  
  summaries <- lapply(cols, function(col) {
    x <- dt[[col]]
    data.table(
      column = col,
      class = class(x)[1],
      n_total = length(x),
      n_valid = sum(!is.na(x)),
      pct_valid = round(100 * sum(!is.na(x)) / length(x), 1),
      n_unique = uniqueN(x, na.rm = TRUE),
      min = if (is.numeric(x)) round(min(x, na.rm = TRUE), 4) else NA_real_,
      max = if (is.numeric(x)) round(max(x, na.rm = TRUE), 4) else NA_real_,
      mean = if (is.numeric(x)) round(mean(x, na.rm = TRUE), 4) else NA_real_
    )
  })
  rbindlist(summaries)
}

# ============================================================================
# MAIN DIAGNOSTIC
# ============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════════╗\n")
cat("║           ENRICHED SITE COLUMNS DIAGNOSTIC (Zonal Sampling)              ║\n")
cat("╚══════════════════════════════════════════════════════════════════════════╝\n")

cat(sprintf("\nData directory: %s\n", ENRICHED_DIR))
cat(sprintf("Sample size for completeness checks: %d rows\n", SAMPLE_ROWS))

# Find all site files
cat("\n─────────────────────────────────────────────────────────────────────────────\n")
cat("SCANNING FOR SITE FILES\n")
cat("─────────────────────────────────────────────────────────────────────────────\n")

site_files <- find_all_site_files(ENRICHED_DIR)

if (nrow(site_files) == 0) {
  cat("\n⚠ No enriched site files found in directory!\n")
  cat("  Check that ENRICHED_DIR path is correct.\n\n")
  quit(status = 1)
}

cat(sprintf("\nFound %d site files:\n", nrow(site_files)))
for (i in seq_len(nrow(site_files))) {
  cat(sprintf("  [%02d] %s\n", site_files$site_id[i], basename(site_files$path[i])))
}

# ============================================================================
# CROSS-SITE COLUMN CONSISTENCY CHECK
# ============================================================================

cat("\n─────────────────────────────────────────────────────────────────────────────\n")
cat("CROSS-SITE COLUMN CONSISTENCY\n")
cat("─────────────────────────────────────────────────────────────────────────────\n")

# Read column names from all sites
site_columns <- list()
site_nrows <- c()

for (i in seq_len(nrow(site_files))) {
  tryCatch({
    dt <- fread(site_files$path[i], nrows = 5)
    site_columns[[as.character(site_files$site_id[i])]] <- names(dt)
    
    # Get row count (read full file for accurate count)
    dt_count <- fread(site_files$path[i], select = 1L)
    site_nrows[as.character(site_files$site_id[i])] <- nrow(dt_count)
  }, error = function(e) {
    cat(sprintf("  ⚠ Error reading site %d: %s\n", site_files$site_id[i], e$message))
  })
}

# Find common and unique columns
all_columns <- unique(unlist(site_columns))
column_presence <- sapply(site_columns, function(cols) all_columns %in% cols)

if (is.vector(column_presence)) {
  # Only one site
  column_presence <- matrix(column_presence, ncol = 1)
  colnames(column_presence) <- names(site_columns)[1]
}

rownames(column_presence) <- all_columns
n_sites_per_col <- rowSums(column_presence)

# Columns in all sites vs some sites
cols_in_all <- names(n_sites_per_col[n_sites_per_col == length(site_columns)])
cols_in_some <- names(n_sites_per_col[n_sites_per_col < length(site_columns) & n_sites_per_col > 0])

cat(sprintf("\nTotal unique columns across all sites: %d\n", length(all_columns)))
cat(sprintf("Columns present in ALL sites: %d\n", length(cols_in_all)))
cat(sprintf("Columns in SOME sites only: %d\n", length(cols_in_some)))

if (length(cols_in_some) > 0 && length(cols_in_some) <= 30) {
  cat("\nColumns with inconsistent presence:\n")
  for (col in cols_in_some) {
    sites_with <- names(site_columns)[sapply(site_columns, function(x) col %in% x)]
    cat(sprintf("  • %s (in sites: %s)\n", col, paste(sites_with, collapse = ", ")))
  }
} else if (length(cols_in_some) > 30) {
  cat(sprintf("\n  (Too many inconsistent columns to list - %d total)\n", length(cols_in_some)))
}

# Row counts per site
cat("\nRow counts per site:\n")
for (site in names(site_nrows)) {
  cat(sprintf("  Site %s: %s rows\n", site, format(site_nrows[site], big.mark = ",")))
}
cat(sprintf("\nTotal rows across all sites: %s\n", format(sum(site_nrows), big.mark = ",")))

# ============================================================================
# COLUMN CATEGORIZATION (using first available site as reference)
# ============================================================================

cat("\n─────────────────────────────────────────────────────────────────────────────\n")
cat("COLUMN CATEGORIZATION\n")
cat("─────────────────────────────────────────────────────────────────────────────\n")

# Use columns common to all sites for categorization
reference_cols <- cols_in_all
categories <- categorize_columns(reference_cols)

cat("\nColumns grouped by source/product type:\n\n")

for (cat_name in names(categories)) {
  cols <- categories[[cat_name]]
  if (length(cols) > 0) {
    cat_display <- switch(cat_name,
      "gedi_core" = "GEDI Core Attributes",
      "dep_dtm" = "3DEP DTM Products",
      "p3d_dtm" = "P3D/MHRSI DTM Products",
      "als_chm" = "ALS CHM Products",
      "p3d_chm" = "P3D/MHRSI CHM Products",
      "error" = "Error/Difference Metrics",
      "terrain" = "Terrain Variables",
      "canopy" = "Canopy Metrics",
      "metadata" = "Site/Location Metadata",
      "other" = "Other/Uncategorized"
    )
    
    cat(sprintf("┌─ %s (%d columns)\n", cat_display, length(cols)))
    
    # Show columns (limit display for large categories)
    max_show <- 20
    cols_to_show <- if (length(cols) > max_show) cols[1:max_show] else cols
    for (col in cols_to_show) {
      cat(sprintf("│    %s\n", col))
    }
    if (length(cols) > max_show) {
      cat(sprintf("│    ... and %d more\n", length(cols) - max_show))
    }
    cat("│\n")
  }
}

# ============================================================================
# ZONAL STATISTICS DETECTION
# ============================================================================

cat("\n─────────────────────────────────────────────────────────────────────────────\n")
cat("ZONAL STATISTICS DETECTED\n")
cat("─────────────────────────────────────────────────────────────────────────────\n")

zonal_stats <- extract_zonal_stats(reference_cols)

if (length(zonal_stats) > 0) {
  cat("\nZonal statistics suffixes found in column names:\n")
  for (stat in names(zonal_stats)) {
    n_cols <- length(zonal_stats[[stat]])
    cat(sprintf("  ✓ _%s : %d columns\n", stat, n_cols))
  }
} else {
  cat("\n⚠ No standard zonal statistics suffixes detected (mean, std, etc.)\n")
  cat("  This might indicate pixel-level data or non-standard naming.\n")
}

# ============================================================================
# DTM/CHM VARIABLE MAPPING CHECK
# ============================================================================

cat("\n─────────────────────────────────────────────────────────────────────────────\n")
cat("DTM/CHM VARIABLE MAPPING\n")
cat("─────────────────────────────────────────────────────────────────────────────\n")

# Check for expected variable pairs for bias calculation
cat("\nChecking for variable pairs needed for error calculation:\n")

# DTM pairs
dtm_ref_patterns <- c("dep_dtm", "3dep_dtm", "dep.*dtm.*mean", "dtm.*dep.*mean")
dtm_test_patterns <- c("p3d_dtm", "mhrsi_dtm", "p3d.*dtm.*mean", "dtm.*p3d.*mean")

dtm_ref_found <- reference_cols[Reduce(`|`, lapply(dtm_ref_patterns, function(p) grepl(p, reference_cols, ignore.case = TRUE)))]
dtm_test_found <- reference_cols[Reduce(`|`, lapply(dtm_test_patterns, function(p) grepl(p, reference_cols, ignore.case = TRUE)))]

cat("\n  DTM Reference (3DEP):\n")
if (length(dtm_ref_found) > 0) {
  for (col in dtm_ref_found) cat(sprintf("    ✓ %s\n", col))
} else {
  cat("    ✗ No 3DEP DTM columns found\n")
}

cat("\n  DTM Test (P3D/MHRSI):\n")
if (length(dtm_test_found) > 0) {
  for (col in dtm_test_found) cat(sprintf("    ✓ %s\n", col))
} else {
  cat("    ✗ No P3D/MHRSI DTM columns found\n")
}

# CHM pairs
chm_ref_patterns <- c("als_chm", "als.*chm.*mean", "chm.*als.*mean")
chm_test_patterns <- c("p3d_chm", "mhrsi_chm", "p3d.*chm.*mean", "chm.*p3d.*mean")

chm_ref_found <- reference_cols[Reduce(`|`, lapply(chm_ref_patterns, function(p) grepl(p, reference_cols, ignore.case = TRUE)))]
chm_test_found <- reference_cols[Reduce(`|`, lapply(chm_test_patterns, function(p) grepl(p, reference_cols, ignore.case = TRUE)))]

cat("\n  CHM Reference (ALS):\n")
if (length(chm_ref_found) > 0) {
  for (col in chm_ref_found) cat(sprintf("    ✓ %s\n", col))
} else {
  cat("    ✗ No ALS CHM columns found\n")
}

cat("\n  CHM Test (P3D/MHRSI):\n")
if (length(chm_test_found) > 0) {
  for (col in chm_test_found) cat(sprintf("    ✓ %s\n", col))
} else {
  cat("    ✗ No P3D/MHRSI CHM columns found\n")
}

# ============================================================================
# DATA COMPLETENESS CHECK (sample from first site)
# ============================================================================

cat("\n─────────────────────────────────────────────────────────────────────────────\n")
cat("DATA COMPLETENESS CHECK\n")
cat("─────────────────────────────────────────────────────────────────────────────\n")

# Read sample from first site
first_site <- site_files$path[1]
cat(sprintf("\nSampling %d rows from: %s\n", SAMPLE_ROWS, basename(first_site)))

dt_sample <- fread(first_site, nrows = SAMPLE_ROWS)

# Focus on key product columns
key_cols <- unique(c(dtm_ref_found, dtm_test_found, chm_ref_found, chm_test_found))
key_cols <- intersect(key_cols, names(dt_sample))

if (length(key_cols) > 0) {
  cat("\nCompleteness for key DTM/CHM columns:\n")
  summary_dt <- get_column_summary(dt_sample, key_cols)
  if (!is.null(summary_dt)) {
    for (i in seq_len(nrow(summary_dt))) {
      row <- summary_dt[i]
      status <- if (row$pct_valid >= 95) "✓" else if (row$pct_valid >= 50) "~" else "✗"
      cat(sprintf("  %s %s: %.1f%% valid (%d/%d), range [%.2f, %.2f]\n",
                  status, row$column, row$pct_valid, row$n_valid, row$n_total,
                  row$min, row$max))
    }
  }
}

# Check for pre-computed error columns
error_cols <- reference_cols[grepl("error|diff|bias", reference_cols, ignore.case = TRUE)]
if (length(error_cols) > 0) {
  cat("\nPre-computed error/difference columns:\n")
  error_cols_sample <- intersect(error_cols, names(dt_sample))
  summary_dt <- get_column_summary(dt_sample, error_cols_sample)
  if (!is.null(summary_dt)) {
    for (i in seq_len(nrow(summary_dt))) {
      row <- summary_dt[i]
      cat(sprintf("  %s: %.1f%% valid, mean=%.4f, range [%.2f, %.2f]\n",
                  row$column, row$pct_valid, row$mean, row$min, row$max))
    }
  }
} else {
  cat("\n⚠ No pre-computed error columns found.\n")
  cat("  Errors will need to be calculated from DTM/CHM pairs.\n")
}

# ============================================================================
# FULL COLUMN LIST EXPORT
# ============================================================================

cat("\n─────────────────────────────────────────────────────────────────────────────\n")
cat("FULL COLUMN LIST\n")
cat("─────────────────────────────────────────────────────────────────────────────\n")

cat(sprintf("\nAll %d columns present in all sites:\n\n", length(cols_in_all)))

# Print in 3 columns for readability
n_cols <- length(cols_in_all)
cols_sorted <- sort(cols_in_all)

for (i in seq_along(cols_sorted)) {
  cat(sprintf("  %3d. %s\n", i, cols_sorted[i]))
}

# ============================================================================
# SUMMARY & RECOMMENDATIONS
# ============================================================================

cat("\n═══════════════════════════════════════════════════════════════════════════\n")
cat("SUMMARY & RECOMMENDATIONS\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")

cat("\n")

# DTM analysis feasibility
if (length(dtm_ref_found) > 0 && length(dtm_test_found) > 0) {
  cat("✓ DTM error analysis: FEASIBLE\n")
  cat(sprintf("  Reference: %s\n", paste(dtm_ref_found, collapse = ", ")))
  cat(sprintf("  Test: %s\n", paste(dtm_test_found, collapse = ", ")))
} else {
  cat("✗ DTM error analysis: NOT FEASIBLE\n")
  if (length(dtm_ref_found) == 0) cat("  Missing: 3DEP DTM reference\n")
  if (length(dtm_test_found) == 0) cat("  Missing: P3D/MHRSI DTM test\n")
}

cat("\n")

# CHM analysis feasibility
if (length(chm_ref_found) > 0 && length(chm_test_found) > 0) {
  cat("✓ CHM error analysis: FEASIBLE\n")
  cat(sprintf("  Reference: %s\n", paste(chm_ref_found, collapse = ", ")))
  cat(sprintf("  Test: %s\n", paste(chm_test_found, collapse = ", ")))
} else {
  cat("✗ CHM error analysis: NOT FEASIBLE\n")
  if (length(chm_ref_found) == 0) cat("  Missing: ALS CHM reference\n")
  if (length(chm_test_found) == 0) cat("  Missing: P3D/MHRSI CHM test\n")
}

cat("\n")

# Zonal stats check
if (length(zonal_stats) > 0) {
  cat(sprintf("✓ Zonal statistics present: %s\n", paste(names(zonal_stats), collapse = ", ")))
} else {
  cat("⚠ No zonal statistics suffixes detected - verify sampling approach\n")
}

cat("\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("Diagnostic complete.\n\n")
