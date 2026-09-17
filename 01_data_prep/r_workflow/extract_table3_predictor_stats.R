#!/usr/bin/env Rscript
# =====================================================================
# extract_table3_predictor_stats_v2.R
# Extract summary statistics for all predictor variables used in modeling
# FIXED: Corrected column names for aspect_sin_mean and aspect_cos_mean
#
# Run on cluster: Rscript extract_table3_predictor_stats_v2.R
# =====================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# =====================================================================
# Configuration
# =====================================================================

PROJECT_ROOT <- "/gpfs/data1/vclgp/lmaden/chpt1"
ENRICHED_DIR <- file.path(PROJECT_ROOT, "data", "enriched_by_site")
OUTPUT_DIR <- file.path(PROJECT_ROOT, "manuscript_tables")

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Filtering parameters (matching your analysis)
VALID_FRAC_THRESH <- 0.5
FOREST_HEIGHT_THRESH <- 2.0
SITES_EXCLUDE_CHM <- c(16)

cat("======================================================================\n")
cat("Table 3: Predictor Variable Statistics Extraction (v2 - Fixed)\n")
cat("======================================================================\n\n")

# =====================================================================
# Define predictor variables - CORRECTED COLUMN NAMES
# =====================================================================

TERRAIN_PREDICTORS <- c(
  "slope_mean",        # Mean terrain slope (degrees)
  "slope_sd",          # Slope standard deviation
  "aspect_sin_mean",   # FIXED: was aspect_sin

"aspect_cos_mean"    # FIXED: was aspect_cos
)

VEGETATION_PREDICTORS <- c(
  "als_chm_mean",      # Reference canopy height - mean
  "als_chm_p90",       # Reference canopy height - 90th percentile  
  "cover",             # Canopy cover fraction (GEDI)
  "wsci",              # Waveform structural complexity index
  "rh_98"              # GEDI relative height 98 (canopy top)
)

ACQUISITION_PREDICTORS <- c(
  "meta_stereo_ratio",       # Stereo pair ratio
  "meta_off_nadir_avg",      # Mean off-nadir view angle
  "meta_sun_elev_avg",       # Mean sun elevation
  "meta_leaf_on_ratio",      # Leaf-on image ratio
  "meta_az_concentration"    # Azimuth concentration
)

# Combine all predictors
ALL_PREDICTORS <- c(TERRAIN_PREDICTORS, VEGETATION_PREDICTORS, ACQUISITION_PREDICTORS)

cat("Looking for these predictor columns:\n")
cat(paste(" -", ALL_PREDICTORS, collapse = "\n"), "\n\n")

# =====================================================================
# Discover and load site files
# =====================================================================

cat("Discovering enriched site files...\n")

site_files <- list.files(ENRICHED_DIR, 
                         pattern = "^site_[0-9]+_enriched\\.csv(\\.gz)?$",
                         full.names = TRUE)

get_site_num <- function(path) {
  bn <- basename(path)
  m <- regmatches(bn, regexec("site_([0-9]+)_enriched", bn))[[1]]
  if (length(m) == 2) as.integer(m[2]) else NA_integer_
}

site_nums <- sapply(site_files, get_site_num)
site_files <- site_files[order(site_nums)]
site_nums <- sort(site_nums)

cat(sprintf("Found %d site files\n\n", length(site_files)))

# =====================================================================
# Load all data
# =====================================================================

cat("Loading data from all sites...\n")

all_data <- list()

for (i in seq_along(site_files)) {
  site_id <- site_nums[i]
  
  # Skip excluded sites
  if (site_id %in% SITES_EXCLUDE_CHM) {
    cat(sprintf("  Site %d: Skipped (excluded)\n", site_id))
    next
  }
  
  cat(sprintf("  Site %d: Loading...", site_id))
  
  tryCatch({
    DT <- fread(site_files[i], showProgress = FALSE)
    
    # Check what columns exist
    available_cols <- intersect(names(DT), ALL_PREDICTORS)
    
    if (i == 1) {
      cat("\n  Available predictor columns in data:\n")
      cat(paste("   -", available_cols, collapse = "\n"), "\n")
      missing_cols <- setdiff(ALL_PREDICTORS, names(DT))
      if (length(missing_cols) > 0) {
        cat("  Missing columns:\n")
        cat(paste("   -", missing_cols, collapse = "\n"), "\n")
      }
    }
    
    # Apply CHM filtering
    # Check which filter columns exist
    filter_cols <- c("p3d_chm_mean", "als_chm_mean", "als_chm_valid_frac", 
                     "p3d_chm_valid_frac", "slope_valid_frac", "als_chm_p90")
    has_filter_cols <- all(filter_cols %in% names(DT))
    
    if (has_filter_cols) {
      DT <- DT[
        is.finite(p3d_chm_mean) & is.finite(als_chm_mean) &
        is.finite(als_chm_valid_frac) & is.finite(p3d_chm_valid_frac) &
        als_chm_valid_frac >= VALID_FRAC_THRESH & 
        p3d_chm_valid_frac >= VALID_FRAC_THRESH &
        is.finite(slope_valid_frac) & slope_valid_frac >= VALID_FRAC_THRESH &
        is.finite(als_chm_p90) & als_chm_p90 >= FOREST_HEIGHT_THRESH
      ]
    }
    
    # Select predictor columns
    if (length(available_cols) > 0) {
      DT_subset <- DT[, ..available_cols]
      all_data[[as.character(site_id)]] <- DT_subset
    }
    
    cat(sprintf(" %s rows, %d predictors\n", 
                format(nrow(DT), big.mark = ","), length(available_cols)))
    
  }, error = function(e) {
    cat(sprintf(" ERROR: %s\n", e$message))
  })
}

# Combine all data
cat("\nCombining data from all sites...\n")
combined_data <- rbindlist(all_data, fill = TRUE)
cat(sprintf("Total observations: %s\n\n", format(nrow(combined_data), big.mark = ",")))

# =====================================================================
# Calculate summary statistics
# =====================================================================

cat("Calculating summary statistics...\n\n")

calc_stats <- function(x, var_name) {
  x_valid <- x[is.finite(x)]
  
  if (length(x_valid) == 0) {
    return(data.frame(
      Variable = var_name,
      n = 0,
      n_missing = length(x),
      Mean = NA,
      SD = NA,
      Min = NA,
      Max = NA
    ))
  }
  
  data.frame(
    Variable = var_name,
    n = length(x_valid),
    n_missing = sum(!is.finite(x)),
    Mean = round(mean(x_valid), 3),
    SD = round(sd(x_valid), 3),
    Min = round(min(x_valid), 3),
    Max = round(max(x_valid), 3)
  )
}

# Calculate stats for all predictors
stats_list <- list()

for (var in ALL_PREDICTORS) {
  if (var %in% names(combined_data)) {
    stats_list[[var]] <- calc_stats(combined_data[[var]], var)
  } else {
    cat(sprintf("  Warning: %s not found in data\n", var))
    stats_list[[var]] <- data.frame(
      Variable = var,
      n = NA, n_missing = NA, Mean = NA, SD = NA, Min = NA, Max = NA
    )
  }
}

stats_df <- rbindlist(stats_list)

# Add category labels
stats_df$Category <- case_when(
  stats_df$Variable %in% TERRAIN_PREDICTORS ~ "Terrain",
  stats_df$Variable %in% VEGETATION_PREDICTORS ~ "Vegetation",
  stats_df$Variable %in% ACQUISITION_PREDICTORS ~ "Acquisition",
  TRUE ~ "Other"
)

# Add descriptions - UPDATED for correct column names
var_descriptions <- c(
  slope_mean = "Mean terrain slope",
  slope_sd = "Slope variability (SD)",
  aspect_sin_mean = "Aspect sine (E-W)",
  aspect_cos_mean = "Aspect cosine (N-S)",
  als_chm_mean = "Reference canopy height (mean)",
  als_chm_p90 = "Reference canopy height (P90)",
  cover = "Canopy cover fraction",
  wsci = "Waveform structural complexity",
  rh_98 = "Relative height 98th percentile",
  meta_stereo_ratio = "Stereo pair ratio",
  meta_off_nadir_avg = "Off-nadir angle (mean)",
  meta_sun_elev_avg = "Sun elevation (mean)",
  meta_leaf_on_ratio = "Leaf-on image ratio",
  meta_az_concentration = "Azimuth concentration"
)

var_units <- c(
  slope_mean = "degrees",
  slope_sd = "degrees",
  aspect_sin_mean = "unitless",
  aspect_cos_mean = "unitless",
  als_chm_mean = "m",
  als_chm_p90 = "m",
  cover = "proportion",
  wsci = "unitless",
  rh_98 = "m",
  meta_stereo_ratio = "proportion",
  meta_off_nadir_avg = "degrees",
  meta_sun_elev_avg = "degrees",
  meta_leaf_on_ratio = "proportion",
  meta_az_concentration = "unitless"
)

var_sources <- c(
  slope_mean = "3DEP DTM",
  slope_sd = "3DEP DTM",
  aspect_sin_mean = "3DEP DTM",
  aspect_cos_mean = "3DEP DTM",
  als_chm_mean = "ALS CHM",
  als_chm_p90 = "ALS CHM",
  cover = "GEDI L2A",
  wsci = "GEDI L2A",
  rh_98 = "GEDI L2A",
  meta_stereo_ratio = "P3D metadata",
  meta_off_nadir_avg = "P3D metadata",
  meta_sun_elev_avg = "P3D metadata",
  meta_leaf_on_ratio = "P3D metadata",
  meta_az_concentration = "P3D metadata"
)

stats_df$Description <- var_descriptions[stats_df$Variable]
stats_df$Units <- var_units[stats_df$Variable]
stats_df$Source <- var_sources[stats_df$Variable]

# Reorder columns
stats_df <- stats_df[, c("Category", "Variable", "Description", "Units", "Source",
                          "n", "Mean", "SD", "Min", "Max")]

# =====================================================================
# Output results
# =====================================================================

output_file <- file.path(OUTPUT_DIR, "Table3_predictor_stats_v2.csv")
write.csv(stats_df, output_file, row.names = FALSE)

cat("======================================================================\n")
cat("Table 3 Predictor Summary Statistics\n")
cat("======================================================================\n\n")

print(stats_df, row.names = FALSE)

cat("\n")
cat(sprintf("Results saved to: %s\n", output_file))

# =====================================================================
# Generate formatted output
# =====================================================================

txt_output <- file.path(OUTPUT_DIR, "Table3_formatted_v2.txt")

sink(txt_output)
cat("Table 3. Predictor variables used in Bayesian hierarchical error models.\n")
cat("Variables are grouped by category with units, data source, and summary statistics.\n")
cat(sprintf("Statistics calculated from %s GEDI footprints across 19 sites.\n\n", 
            format(nrow(combined_data), big.mark = ",")))

for (cat_name in c("Terrain", "Vegetation", "Acquisition")) {
  cat_data <- stats_df[stats_df$Category == cat_name, ]
  cat(sprintf("\n%s Predictors:\n", cat_name))
  cat(strrep("-", 60), "\n")
  
  for (i in 1:nrow(cat_data)) {
    row <- cat_data[i, ]
    if (!is.na(row$Mean)) {
      range_str <- sprintf("%.2f–%.2f", row$Min, row$Max)
      cat(sprintf("  %s (%s): Mean=%.2f, SD=%.2f, Range=%s\n",
                  row$Description, row$Units, row$Mean, row$SD, range_str))
    } else {
      cat(sprintf("  %s (%s): N/A\n", row$Description, row$Units))
    }
  }
}
sink()

cat(sprintf("Formatted table saved to: %s\n", txt_output))

cat("\nDone!\n")
