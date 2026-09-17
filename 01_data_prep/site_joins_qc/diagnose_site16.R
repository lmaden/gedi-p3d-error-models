# =====================================================================
# diagnose_site16.R
# Investigate why Site 16 has 0 footprints after QC filtering
# 
# Run: Rscript diagnose_site16.R
# =====================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(sf)
})

PROJECT_ROOT <- "/gpfs/data1/vclgp/lmaden/chpt1"
ENRICHED_DIR <- file.path(PROJECT_ROOT, "data", "enriched_by_site")
GEDI_DIR     <- file.path(PROJECT_ROOT, "gedi")

cat("=" |> rep(70) |> paste(collapse=""), "\n")
cat("Site 16 Diagnostic Report\n")
cat("=" |> rep(70) |> paste(collapse=""), "\n\n")

# =====================================================================
# 1. Load Site 16 enriched CSV
# =====================================================================

site16_path <- file.path(ENRICHED_DIR, "site_16_enriched.csv.gz")
if (!file.exists(site16_path)) {
  site16_path <- file.path(ENRICHED_DIR, "site_16_enriched.csv")
}

cat("Loading Site 16 data from:", site16_path, "\n\n")

DT <- fread(site16_path, showProgress = FALSE)
cat(sprintf("Total rows in raw file: %s\n\n", format(nrow(DT), big.mark=",")))

# =====================================================================
# 2. Check column availability
# =====================================================================

cat("=" |> rep(50) |> paste(collapse=""), "\n")
cat("COLUMN AVAILABILITY CHECK\n")
cat("=" |> rep(50) |> paste(collapse=""), "\n")

key_cols <- c(
  "p3d_chm_mean", "als_chm_mean", "als_chm_p90",
  "als_chm_valid_frac", "p3d_chm_valid_frac",
  "p3d_dtm_mean", "dep_dtm_mean",
  "dep_dtm_valid_frac", "p3d_dtm_valid_frac",
  "slope_valid_frac", "slope_mean",
  "ecoregion", "lc2022_l1_code", "shot_number"
)

for (col in key_cols) {
  present <- col %in% names(DT)
  if (present) {
    n_valid <- sum(is.finite(DT[[col]]) | !is.na(DT[[col]]))
    n_na <- sum(is.na(DT[[col]]))
    cat(sprintf("  %-25s: PRESENT (valid: %s, NA: %s)\n", 
                col, format(n_valid, big.mark=","), format(n_na, big.mark=",")))
  } else {
    cat(sprintf("  %-25s: MISSING!\n", col))
  }
}

# =====================================================================
# 3. Check value distributions for key filtering variables
# =====================================================================

cat("\n")
cat("=" |> rep(50) |> paste(collapse=""), "\n")
cat("VALUE DISTRIBUTIONS (filtering variables)\n")
cat("=" |> rep(50) |> paste(collapse=""), "\n")

check_distribution <- function(DT, col) {
  if (!(col %in% names(DT))) {
    cat(sprintf("\n%s: Column not present\n", col))
    return()
  }
  
  vals <- DT[[col]]
  n_total <- length(vals)
  n_na <- sum(is.na(vals))
  n_finite <- sum(is.finite(vals))
  
  cat(sprintf("\n%s:\n", col))
  cat(sprintf("  Total: %s | NA: %s | Finite: %s\n", 
              format(n_total, big.mark=","),
              format(n_na, big.mark=","),
              format(n_finite, big.mark=",")))
  
  if (n_finite > 0) {
    finite_vals <- vals[is.finite(vals)]
    cat(sprintf("  Min: %.4f | Max: %.4f | Mean: %.4f | Median: %.4f\n",
                min(finite_vals), max(finite_vals), 
                mean(finite_vals), median(finite_vals)))
    
    # Show quantiles
    qs <- quantile(finite_vals, probs = c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99))
    cat(sprintf("  Quantiles: 1%%=%.2f, 5%%=%.2f, 25%%=%.2f, 50%%=%.2f, 75%%=%.2f, 95%%=%.2f, 99%%=%.2f\n",
                qs[1], qs[2], qs[3], qs[4], qs[5], qs[6], qs[7]))
  }
}

check_distribution(DT, "als_chm_p90")
check_distribution(DT, "als_chm_mean")
check_distribution(DT, "p3d_chm_mean")
check_distribution(DT, "als_chm_valid_frac")
check_distribution(DT, "p3d_chm_valid_frac")
check_distribution(DT, "slope_valid_frac")
check_distribution(DT, "dep_dtm_valid_frac")
check_distribution(DT, "p3d_dtm_valid_frac")

# =====================================================================
# 4. Step-by-step filtering to identify where data is lost
# =====================================================================

cat("\n")
cat("=" |> rep(50) |> paste(collapse=""), "\n")
cat("STEP-BY-STEP CHM FILTERING\n")
cat("=" |> rep(50) |> paste(collapse=""), "\n")

n0 <- nrow(DT)
cat(sprintf("\nStarting rows: %s\n", format(n0, big.mark=",")))

# Step 1: finite p3d_chm_mean
if ("p3d_chm_mean" %in% names(DT)) {
  n1 <- sum(is.finite(DT$p3d_chm_mean))
  cat(sprintf("After is.finite(p3d_chm_mean): %s (lost %s)\n", 
              format(n1, big.mark=","), format(n0 - n1, big.mark=",")))
} else {
  cat("p3d_chm_mean: COLUMN MISSING - this is the problem!\n")
  n1 <- 0
}

# Step 2: finite als_chm_mean
if ("als_chm_mean" %in% names(DT)) {
  n2 <- sum(is.finite(DT$p3d_chm_mean) & is.finite(DT$als_chm_mean), na.rm = TRUE)
  cat(sprintf("After is.finite(als_chm_mean): %s (lost %s)\n", 
              format(n2, big.mark=","), format(n1 - n2, big.mark=",")))
} else {
  cat("als_chm_mean: COLUMN MISSING - this is the problem!\n")
  n2 <- 0
}

# Step 3: valid fractions for CHM
if (all(c("als_chm_valid_frac", "p3d_chm_valid_frac") %in% names(DT))) {
  n3 <- sum(
    is.finite(DT$p3d_chm_mean) & is.finite(DT$als_chm_mean) &
    is.finite(DT$als_chm_valid_frac) & is.finite(DT$p3d_chm_valid_frac) &
    DT$als_chm_valid_frac >= 0.5 & DT$p3d_chm_valid_frac >= 0.5,
    na.rm = TRUE
  )
  cat(sprintf("After valid_frac >= 0.5: %s (lost %s)\n", 
              format(n3, big.mark=","), format(n2 - n3, big.mark=",")))
} else {
  cat("CHM valid_frac columns: MISSING\n")
  n3 <- 0
}

# Step 4: slope valid fraction
if ("slope_valid_frac" %in% names(DT)) {
  n4 <- sum(
    is.finite(DT$p3d_chm_mean) & is.finite(DT$als_chm_mean) &
    is.finite(DT$als_chm_valid_frac) & is.finite(DT$p3d_chm_valid_frac) &
    DT$als_chm_valid_frac >= 0.5 & DT$p3d_chm_valid_frac >= 0.5 &
    is.finite(DT$slope_valid_frac) & DT$slope_valid_frac >= 0.5,
    na.rm = TRUE
  )
  cat(sprintf("After slope_valid_frac >= 0.5: %s (lost %s)\n", 
              format(n4, big.mark=","), format(n3 - n4, big.mark=",")))
} else {
  cat("slope_valid_frac: COLUMN MISSING\n")
  n4 <- 0
}

# Step 5: forest filter (als_chm_p90 >= 2)
if ("als_chm_p90" %in% names(DT)) {
  n5 <- sum(
    is.finite(DT$p3d_chm_mean) & is.finite(DT$als_chm_mean) &
    is.finite(DT$als_chm_valid_frac) & is.finite(DT$p3d_chm_valid_frac) &
    DT$als_chm_valid_frac >= 0.5 & DT$p3d_chm_valid_frac >= 0.5 &
    is.finite(DT$slope_valid_frac) & DT$slope_valid_frac >= 0.5 &
    is.finite(DT$als_chm_p90) & DT$als_chm_p90 >= 2,
    na.rm = TRUE
  )
  cat(sprintf("After als_chm_p90 >= 2m (forest): %s (lost %s)\n", 
              format(n5, big.mark=","), format(n4 - n5, big.mark=",")))
} else {
  cat("als_chm_p90: COLUMN MISSING\n")
  n5 <- 0
}

cat(sprintf("\nFINAL CHM FOOTPRINTS: %s\n", format(n5, big.mark=",")))

# =====================================================================
# 5. Check land cover distribution
# =====================================================================

cat("\n")
cat("=" |> rep(50) |> paste(collapse=""), "\n")
cat("LAND COVER DISTRIBUTION\n")
cat("=" |> rep(50) |> paste(collapse=""), "\n")

lc_col <- if ("lc2022_l1_code" %in% names(DT)) "lc2022_l1_code" else 
          if ("lc2022_mode_l1_code" %in% names(DT)) "lc2022_mode_l1_code" else NULL

if (!is.null(lc_col)) {
  lc_table <- sort(table(DT[[lc_col]]), decreasing = TRUE)
  cat("\nLand cover classes (raw data):\n")
  for (i in seq_along(lc_table)) {
    pct <- 100 * lc_table[i] / sum(lc_table)
    cat(sprintf("  %s: %s (%.1f%%)\n", 
                names(lc_table)[i], 
                format(lc_table[i], big.mark=","),
                pct))
  }
} else {
  cat("No land cover column found\n")
}

# =====================================================================
# 6. Check ecoregion
# =====================================================================

cat("\n")
cat("=" |> rep(50) |> paste(collapse=""), "\n")
cat("ECOREGION\n")
cat("=" |> rep(50) |> paste(collapse=""), "\n")

if ("ecoregion" %in% names(DT)) {
  eco_table <- table(DT$ecoregion)
  cat(sprintf("Ecoregion: %s\n", names(eco_table)[1]))
} else {
  cat("Ecoregion column not found\n")
}

# =====================================================================
# 7. Check geopackage
# =====================================================================

cat("\n")
cat("=" |> rep(50) |> paste(collapse=""), "\n")
cat("GEOPACKAGE CHECK\n")
cat("=" |> rep(50) |> paste(collapse=""), "\n")

gpkg_path <- file.path(GEDI_DIR, "16", "GEDI_site16_hq_ALL.gpkg")
cat(sprintf("Expected path: %s\n", gpkg_path))
cat(sprintf("File exists: %s\n", file.exists(gpkg_path)))

if (file.exists(gpkg_path)) {
  gdf <- tryCatch({
    st_read(gpkg_path, quiet = TRUE)
  }, error = function(e) {
    cat(sprintf("Error reading gpkg: %s\n", e$message))
    NULL
  })
  
  if (!is.null(gdf)) {
    cat(sprintf("Rows in gpkg: %s\n", format(nrow(gdf), big.mark=",")))
    cat(sprintf("CRS: %s\n", st_crs(gdf)$input))
    
    # Get bounds
    bbox <- st_bbox(gdf)
    cat(sprintf("Bounding box: [%.4f, %.4f] to [%.4f, %.4f]\n",
                bbox["xmin"], bbox["ymin"], bbox["xmax"], bbox["ymax"]))
  }
}

# =====================================================================
# 8. Sample of raw data
# =====================================================================

cat("\n")
cat("=" |> rep(50) |> paste(collapse=""), "\n")
cat("SAMPLE OF RAW DATA (first 5 rows, key columns)\n")
cat("=" |> rep(50) |> paste(collapse=""), "\n\n")

sample_cols <- intersect(
  c("shot_number", "als_chm_p90", "als_chm_mean", "p3d_chm_mean", 
    "als_chm_valid_frac", "p3d_chm_valid_frac", "slope_valid_frac",
    "lc2022_l1_code", "ecoregion"),
  names(DT)
)

print(DT[1:min(5, nrow(DT)), ..sample_cols])

cat("\n")
cat("=" |> rep(70) |> paste(collapse=""), "\n")
cat("END OF DIAGNOSTIC REPORT\n")
cat("=" |> rep(70) |> paste(collapse=""), "\n")
