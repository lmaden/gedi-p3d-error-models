# =====================================================================
# extract_section2_1_study_sites_v2.R
# Extract data for Section 2.1 Study Sites (Table 1 + prose)
# 
# v2: Added elev_q25 and elev_q75 for DTM IQR column
# Run on cluster with: Rscript extract_section2_1_study_sites_v2.R
# =====================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(sf)
  library(dplyr)
  library(tidyr)
})

# =====================================================================
# Configuration
# =====================================================================

PROJECT_ROOT <- "/gpfs/data1/vclgp/lmaden/chpt1"
ENRICHED_DIR <- file.path(PROJECT_ROOT, "data", "enriched_by_site")
GEDI_DIR     <- file.path(PROJECT_ROOT, "gedi")
OUTPUT_DIR   <- file.path(PROJECT_ROOT, "manuscript_tables")

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Filtering parameters (matching section_01_ingest.R)
VALID_FRAC_THRESH   <- 0.5
FOREST_HEIGHT_THRESH <- 2.0  # als_chm_p90 >= 2m
ERROR_ABS_MAX       <- 100   # |error| <= 100m

# Site exclusions (matching analysis_config.R)
# Site 10: +143m DTM offset (vertical datum issue) - exclude from DTM only
# Site 16: Non-forest site (grassland/cropland, no slope data) - exclude from both
SITES_EXCLUDE_CHM <- c(16)    # Site 16 excluded (non-forest)
SITES_EXCLUDE_DTM <- c(10, 16) # Site 10 (datum issue) + Site 16 (non-forest)

cat("=" |> rep(70) |> paste(collapse=""), "\n")
cat("Extracting Section 2.1 Study Sites Data (v2 - with DTM IQR)\n")
cat("=" |> rep(70) |> paste(collapse=""), "\n\n")

# =====================================================================
# Pure R State/Province lookup - NO EXTERNAL PACKAGES NEEDED
# Uses bounding boxes for approximate location (sufficient for site centroids)
# =====================================================================

cat("Setting up state/province lookup (pure R, no external packages)...\n")

# US States with approximate bounding boxes
# Format: name, abbrev, lon_min, lon_max, lat_min, lat_max
# Ordered roughly by priority/likelihood for GEDI forest sites
us_states_bb <- list(
  # Pacific Northwest / West Coast
  list("Washington", "WA", -124.8, -116.9, 45.5, 49.0),
  list("Oregon", "OR", -124.6, -116.5, 41.9, 46.3),
  list("California", "CA", -124.5, -114.1, 32.5, 42.0),
  # Northern Rockies
  list("Idaho", "ID", -117.3, -111.0, 42.0, 49.0),
  list("Montana", "MT", -116.1, -104.0, 44.4, 49.0),
  list("Wyoming", "WY", -111.1, -104.1, 41.0, 45.0),
  # Southwest
  list("Colorado", "CO", -109.1, -102.0, 37.0, 41.0),
  list("Utah", "UT", -114.1, -109.0, 37.0, 42.0),
  list("Nevada", "NV", -120.0, -114.0, 35.0, 42.0),
  list("Arizona", "AZ", -114.8, -109.0, 31.3, 37.0),
  list("New Mexico", "NM", -109.1, -103.0, 31.3, 37.0),
  # Northern Great Plains
  list("North Dakota", "ND", -104.1, -96.6, 45.9, 49.0),
  list("South Dakota", "SD", -104.1, -96.4, 42.5, 46.0),
  list("Nebraska", "NE", -104.1, -95.3, 40.0, 43.0),
  list("Minnesota", "MN", -97.3, -89.5, 43.5, 49.4),
  list("Wisconsin", "WI", -92.9, -86.2, 42.5, 47.1),
  list("Iowa", "IA", -96.7, -90.1, 40.4, 43.5),
  # Great Lakes / Upper Midwest
  list("Michigan", "MI", -90.5, -82.1, 41.7, 48.3),
  list("Illinois", "IL", -91.5, -87.0, 37.0, 42.5),
  list("Indiana", "IN", -88.1, -84.8, 37.8, 41.8),
  list("Ohio", "OH", -84.9, -80.5, 38.4, 42.0),
  # Northeast
  list("Maine", "ME", -71.1, -66.9, 43.0, 47.5),
  list("New Hampshire", "NH", -72.6, -70.6, 42.7, 45.3),
  list("Vermont", "VT", -73.5, -71.5, 42.7, 45.0),
  list("Massachusetts", "MA", -73.5, -69.9, 41.2, 42.9),
  list("New York", "NY", -79.8, -71.9, 40.5, 45.1),
  list("Pennsylvania", "PA", -80.6, -74.7, 39.7, 42.3),
  list("New Jersey", "NJ", -75.6, -73.9, 38.9, 41.4),
  list("Connecticut", "CT", -73.8, -71.8, 40.9, 42.1),
  list("Rhode Island", "RI", -71.9, -71.1, 41.1, 42.0),
  # Mid-Atlantic / Southeast
  list("Maryland", "MD", -79.5, -75.0, 37.9, 39.7),
  list("Delaware", "DE", -75.8, -75.0, 38.4, 39.8),
  list("Virginia", "VA", -83.7, -75.2, 36.5, 39.5),
  list("West Virginia", "WV", -82.7, -77.7, 37.2, 40.6),
  list("North Carolina", "NC", -84.3, -75.5, 33.8, 36.6),
  list("South Carolina", "SC", -83.4, -78.5, 32.0, 35.2),
  list("Georgia", "GA", -85.6, -80.8, 30.4, 35.0),
  list("Florida", "FL", -87.6, -80.0, 24.5, 31.0),
  list("Tennessee", "TN", -90.3, -81.6, 35.0, 36.7),
  list("Kentucky", "KY", -89.6, -81.9, 36.5, 39.2),
  list("Alabama", "AL", -88.5, -84.9, 30.2, 35.0),
  list("Mississippi", "MS", -91.7, -88.1, 30.2, 35.0),
  list("Louisiana", "LA", -94.1, -88.8, 29.0, 33.0),
  list("Arkansas", "AR", -94.6, -89.6, 33.0, 36.5),
  list("Missouri", "MO", -95.8, -89.1, 36.0, 40.6),
  list("Kansas", "KS", -102.1, -94.6, 37.0, 40.0),
  list("Oklahoma", "OK", -103.0, -94.4, 33.6, 37.0),
  list("Texas", "TX", -106.7, -93.5, 25.8, 36.5),
  # Alaska
  list("Alaska", "AK", -180.0, -129.0, 51.0, 71.5)
)

# Canadian provinces with bounding boxes
ca_provinces_bb <- list(
  list("British Columbia", "BC", -139.1, -114.0, 48.3, 60.0),
  list("Alberta", "AB", -120.0, -110.0, 49.0, 60.0),
  list("Saskatchewan", "SK", -110.0, -101.4, 49.0, 60.0),
  list("Manitoba", "MB", -102.0, -88.9, 49.0, 60.0),
  list("Ontario", "ON", -95.2, -74.3, 41.7, 56.9),
  list("Quebec", "QC", -79.8, -57.0, 45.0, 62.6),
  list("New Brunswick", "NB", -69.1, -63.7, 44.5, 48.1),
  list("Nova Scotia", "NS", -66.4, -59.7, 43.4, 47.1),
  list("Prince Edward Island", "PE", -64.5, -62.0, 45.9, 47.1),
  list("Newfoundland and Labrador", "NL", -67.8, -52.6, 46.6, 60.4)
)

# Function to get state/province from coordinates (pure R, no external packages)
get_state_province <- function(lon, lat) {
  if (is.na(lon) || is.na(lat)) {
    return(list(name = NA_character_, code = NA_character_))
  }
  
  # Check US states first
  for (st in us_states_bb) {
    if (lon >= st[[3]] && lon <= st[[4]] && 
        lat >= st[[5]] && lat <= st[[6]]) {
      return(list(name = st[[1]], code = paste0("US-", st[[2]])))
    }
  }
  
  # Check Canadian provinces
  for (prov in ca_provinces_bb) {
    if (lon >= prov[[3]] && lon <= prov[[4]] && 
        lat >= prov[[5]] && lat <= prov[[6]]) {
      return(list(name = prov[[1]], code = paste0("CA-", prov[[2]])))
    }
  }
  
  # Not found - return coordinates for manual lookup
  return(list(name = NA_character_, code = NA_character_))
}

cat("  Ready (pure R bounding box lookup - no compilation needed)\n\n")

# =====================================================================
# Discover and load all site files
# =====================================================================

cat("Discovering site files...\n")

site_files <- list.files(ENRICHED_DIR, 
                         pattern = "^site_[0-9]+_enriched\\.csv(\\.gz)?$",
                         full.names = TRUE)

# Extract site numbers from filenames
get_site_num <- function(path) {
  bn <- basename(path)
  m <- regmatches(bn, regexec("site_([0-9]+)_enriched", bn))[[1]]
  if (length(m) == 2) as.integer(m[2]) else NA_integer_
}

site_nums <- sapply(site_files, get_site_num)
site_files <- site_files[order(site_nums)]
site_nums <- sort(site_nums)

cat(sprintf("Found %d site files: %s\n\n", 
            length(site_files), 
            paste(site_nums, collapse=", ")))

# =====================================================================
# Process each site
# =====================================================================

results <- list()

for (i in seq_along(site_files)) {
  site_id <- site_nums[i]
  site_path <- site_files[i]
  
  cat(sprintf("[Site %02d] Loading...", site_id))
  
  # Read the enriched CSV
  DT <- fread(site_path, showProgress = FALSE)
  n_raw <- nrow(DT)
  
  # =========================================
  # Apply CHM filtering (same as section_01)
  # =========================================
  if (site_id %in% SITES_EXCLUDE_CHM) {
    CHM <- data.table()
  } else {
    CHM <- DT[
      is.finite(p3d_chm_mean) & is.finite(als_chm_mean) &
      is.finite(als_chm_valid_frac) & is.finite(p3d_chm_valid_frac) &
      als_chm_valid_frac >= VALID_FRAC_THRESH & 
      p3d_chm_valid_frac >= VALID_FRAC_THRESH &
      is.finite(slope_valid_frac) & slope_valid_frac >= VALID_FRAC_THRESH &
      is.finite(als_chm_p90) & als_chm_p90 >= FOREST_HEIGHT_THRESH
    ]
    # Calculate error and filter extremes
    CHM[, chm_error := p3d_chm_mean - als_chm_mean]
    CHM <- CHM[abs(chm_error) <= ERROR_ABS_MAX | !is.finite(chm_error)]
  }
  
  # =========================================
  # Apply DTM filtering (same as section_01)
  # =========================================
  if (site_id %in% SITES_EXCLUDE_DTM) {
    DTM <- data.table()
  } else {
    DTM <- DT[
      is.finite(p3d_dtm_mean) & is.finite(dep_dtm_mean) &
      is.finite(dep_dtm_valid_frac) & is.finite(p3d_dtm_valid_frac) &
      dep_dtm_valid_frac >= VALID_FRAC_THRESH & 
      p3d_dtm_valid_frac >= VALID_FRAC_THRESH &
      is.finite(slope_valid_frac) & slope_valid_frac >= VALID_FRAC_THRESH &
      is.finite(als_chm_p90) & als_chm_p90 >= FOREST_HEIGHT_THRESH
    ]
    # Calculate error and filter extremes
    DTM[, dtm_error := p3d_dtm_mean - dep_dtm_mean]
    DTM <- DTM[abs(dtm_error) <= ERROR_ABS_MAX | !is.finite(dtm_error)]
  }
  
  # =========================================
  # Load geopackage for coordinates
  # =========================================
  gpkg_path <- file.path(GEDI_DIR, site_id, 
                         sprintf("GEDI_site%d_hq_ALL.gpkg", site_id))
  
  coords_df <- NULL
  if (file.exists(gpkg_path)) {
    gdf <- tryCatch({
      st_read(gpkg_path, quiet = TRUE) |>
        st_transform(4326)
    }, error = function(e) NULL)
    
    if (!is.null(gdf) && nrow(gdf) > 0) {
      # Extract coordinates properly - st_coordinates on the whole sf object
      coord_matrix <- st_coordinates(gdf)
      coords_df <- data.frame(
        shot_number = as.character(gdf$shot_number),
        lon = coord_matrix[, 1],
        lat = coord_matrix[, 2],
        stringsAsFactors = FALSE
      )
    }
  }
  
  # =========================================
  # Calculate site-level summaries
  # =========================================
  
  # Get ecoregion (most common value from raw data - this is a site characteristic)
  ecoregion <- if (nrow(DT) > 0 && "ecoregion" %in% names(DT)) {
    ec <- DT$ecoregion[!is.na(DT$ecoregion)]
    if (length(ec) > 0) names(sort(table(ec), decreasing=TRUE))[1] else NA_character_
  } else NA_character_
  
  # Get dominant land cover from CHM-FILTERED data (not raw!)
  # This represents the actual forest composition in the analysis
  lc_col <- if ("lc2022_l1_code" %in% names(CHM)) "lc2022_l1_code" else 
            if ("lc2022_mode_l1_code" %in% names(CHM)) "lc2022_mode_l1_code" else 
            if ("lc2022_l1_code" %in% names(DT)) "lc2022_l1_code" else
            if ("lc2022_mode_l1_code" %in% names(DT)) "lc2022_mode_l1_code" else NULL
  
  # Use CHM data if available, otherwise fall back to raw
  lc_source <- if (nrow(CHM) > 0 && !is.null(lc_col) && lc_col %in% names(CHM)) CHM else DT
  
  dominant_lc <- if (!is.null(lc_col) && nrow(lc_source) > 0 && lc_col %in% names(lc_source)) {
    lc <- lc_source[[lc_col]][!is.na(lc_source[[lc_col]]) & lc_source[[lc_col]] != ""]
    if (length(lc) > 0) names(sort(table(lc), decreasing=TRUE))[1] else NA_character_
  } else NA_character_
  
  # Merge coordinates with CHM data for site extent calculation
  CHM_with_coords <- NULL
  if (!is.null(coords_df) && nrow(CHM) > 0) {
    CHM[, shot_number := as.character(shot_number)]
    CHM_with_coords <- merge(CHM, coords_df, by = "shot_number", all.x = TRUE)
  }
  
  # Calculate geographic extent
  if (!is.null(CHM_with_coords) && nrow(CHM_with_coords) > 0 &&
      sum(!is.na(CHM_with_coords$lon)) > 0) {
    lon_range <- range(CHM_with_coords$lon, na.rm = TRUE)
    lat_range <- range(CHM_with_coords$lat, na.rm = TRUE)
    centroid_lon <- mean(lon_range)
    centroid_lat <- mean(lat_range)
    
    # Approximate area using bounding box (in km²)
    # Using simple spherical approximation
    lat_km <- 111.32  # km per degree latitude
    lon_km <- 111.32 * cos(centroid_lat * pi / 180)  # km per degree longitude
    width_km <- abs(diff(lon_range)) * lon_km
    height_km <- abs(diff(lat_range)) * lat_km
    area_km2 <- width_km * height_km
  } else {
    lon_range <- c(NA, NA)
    lat_range <- c(NA, NA)
    centroid_lon <- NA_real_
    centroid_lat <- NA_real_
    area_km2 <- NA_real_
  }
  
  # =========================================
  # Reverse geocode to get state/province
  # =========================================
  state_info <- get_state_province(centroid_lon, centroid_lat)
  
  # Footprint counts
  n_chm <- nrow(CHM)
  n_dtm <- nrow(DTM)
  
  # Slope statistics (from CHM dataset)
  slope_stats <- if (nrow(CHM) > 0 && "slope_mean" %in% names(CHM)) {
    sl <- CHM$slope_mean[is.finite(CHM$slope_mean)]
    if (length(sl) > 10) {
      qs <- quantile(sl, probs = c(0.25, 0.75), na.rm = TRUE)
      c(
        slope_q25 = as.numeric(qs[1]),
        slope_q75 = as.numeric(qs[2]),
        slope_min = min(sl),
        slope_max = max(sl),
        slope_median = median(sl)
      )
    } else {
      c(slope_q25 = NA_real_, slope_q75 = NA_real_, slope_min = NA_real_, 
        slope_max = NA_real_, slope_median = NA_real_)
    }
  } else {
    c(slope_q25 = NA_real_, slope_q75 = NA_real_, slope_min = NA_real_, 
      slope_max = NA_real_, slope_median = NA_real_)
  }
  
  # Canopy height statistics (als_chm_p90 from CHM dataset)
  chm_stats <- if (nrow(CHM) > 0 && "als_chm_p90" %in% names(CHM)) {
    ch <- CHM$als_chm_p90[is.finite(CHM$als_chm_p90)]
    if (length(ch) > 10) {
      qs <- quantile(ch, probs = c(0.25, 0.75), na.rm = TRUE)
      c(
        chm_q25 = as.numeric(qs[1]),
        chm_q75 = as.numeric(qs[2]),
        chm_min = min(ch),
        chm_max = max(ch),
        chm_median = median(ch)
      )
    } else {
      c(chm_q25 = NA_real_, chm_q75 = NA_real_, chm_min = NA_real_, 
        chm_max = NA_real_, chm_median = NA_real_)
    }
  } else {
    c(chm_q25 = NA_real_, chm_q75 = NA_real_, chm_min = NA_real_, 
      chm_max = NA_real_, chm_median = NA_real_)
  }
  
  # =========================================
  # Elevation statistics - NOW WITH Q25/Q75!
  # (dep_dtm_mean from DTM dataset)
  # =========================================
  elev_stats <- if (nrow(DTM) > 0 && "dep_dtm_mean" %in% names(DTM)) {
    el <- DTM$dep_dtm_mean[is.finite(DTM$dep_dtm_mean)]
    if (length(el) > 10) {
      qs <- quantile(el, probs = c(0.25, 0.75), na.rm = TRUE)
      c(
        elev_q25 = as.numeric(qs[1]),
        elev_q75 = as.numeric(qs[2]),
        elev_min = min(el),
        elev_max = max(el),
        elev_median = median(el)
      )
    } else {
      c(elev_q25 = NA_real_, elev_q75 = NA_real_, elev_min = NA_real_, 
        elev_max = NA_real_, elev_median = NA_real_)
    }
  } else {
    c(elev_q25 = NA_real_, elev_q75 = NA_real_, elev_min = NA_real_, 
      elev_max = NA_real_, elev_median = NA_real_)
  }
  
  # Store results
  results[[i]] <- data.frame(
    site_id = site_id,
    state_name = state_info$name,
    state_code = state_info$code,
    n_raw = n_raw,
    n_chm = n_chm,
    n_dtm = n_dtm,
    ecoregion = ecoregion,
    dominant_lc = dominant_lc,
    centroid_lat = centroid_lat,
    centroid_lon = centroid_lon,
    lat_min = lat_range[1],
    lat_max = lat_range[2],
    lon_min = lon_range[1],
    lon_max = lon_range[2],
    area_km2 = area_km2,
    slope_q25 = slope_stats[["slope_q25"]],
    slope_q75 = slope_stats[["slope_q75"]],
    slope_min = slope_stats[["slope_min"]],
    slope_max = slope_stats[["slope_max"]],
    slope_median = slope_stats[["slope_median"]],
    chm_q25 = chm_stats[["chm_q25"]],
    chm_q75 = chm_stats[["chm_q75"]],
    chm_min = chm_stats[["chm_min"]],
    chm_max = chm_stats[["chm_max"]],
    chm_median = chm_stats[["chm_median"]],
    elev_q25 = elev_stats[["elev_q25"]],
    elev_q75 = elev_stats[["elev_q75"]],
    elev_min = elev_stats[["elev_min"]],
    elev_max = elev_stats[["elev_max"]],
    elev_median = elev_stats[["elev_median"]],
    stringsAsFactors = FALSE
  )
  
  cat(sprintf(" Done (CHM: %s, DTM: %s)\n", 
              format(n_chm, big.mark=","),
              format(n_dtm, big.mark=",")))
}

# Combine results
site_summary <- bind_rows(results)

# =====================================================================
# Calculate overall statistics (forest sites only)
# =====================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse=""), "\n")
cat("Calculating overall statistics (forest sites only)...\n")
cat("=" |> rep(70) |> paste(collapse=""), "\n")

# Filter to forest sites (exclude Site 16)
forest_sites <- site_summary |> filter(!(site_id %in% c(16)))

overall_stats <- list(
  # Geographic extent
  lat_range = c(min(forest_sites$lat_min, na.rm=TRUE), 
                max(forest_sites$lat_max, na.rm=TRUE)),
  lon_range = c(min(forest_sites$lon_min, na.rm=TRUE),
                max(forest_sites$lon_max, na.rm=TRUE)),
  
  # Elevation range
  elev_range = c(min(forest_sites$elev_min, na.rm=TRUE),
                 max(forest_sites$elev_max, na.rm=TRUE)),
  
  # Canopy height range
  chm_range = c(min(forest_sites$chm_min, na.rm=TRUE),
                max(forest_sites$chm_max, na.rm=TRUE)),
  
  # Slope range
  slope_range = c(min(forest_sites$slope_min, na.rm=TRUE),
                  max(forest_sites$slope_max, na.rm=TRUE)),
  
  # Sample sizes
  total_chm = sum(forest_sites$n_chm, na.rm=TRUE),
  total_dtm = sum(forest_sites$n_dtm, na.rm=TRUE),
  total_raw = sum(forest_sites$n_raw, na.rm=TRUE),
  
  # Total area
  total_area_km2 = sum(forest_sites$area_km2, na.rm=TRUE),
  
  # Count of sites (forest sites only)
  n_sites_total = nrow(forest_sites),
  n_sites_chm = sum(forest_sites$n_chm > 0),
  n_sites_dtm = sum(forest_sites$n_dtm > 0),
  
  # Ecoregions (forest sites only)
  n_ecoregions = length(unique(forest_sites$ecoregion[!is.na(forest_sites$ecoregion)])),
  ecoregions = unique(forest_sites$ecoregion[!is.na(forest_sites$ecoregion)]),
  
  # Forest types (forest sites only)
  forest_types = unique(forest_sites$dominant_lc[!is.na(forest_sites$dominant_lc)]),
  
  # States/provinces (forest sites only)
  states = unique(forest_sites$state_name[!is.na(forest_sites$state_name)]),
  
  # Note about excluded sites
  excluded_sites = "Site 16 (non-forest: grassland/cropland with no slope data)"
)

# =====================================================================
# Output results
# =====================================================================

# Save full summary table
write.csv(site_summary, 
          file.path(OUTPUT_DIR, "section2_1_site_summary_full_v2.csv"),
          row.names = FALSE)

# Create formatted Table 1 for manuscript - NOW WITH DTM IQR
# Exclude Site 16 (non-forest site) from the table entirely
table1 <- site_summary |>
  filter(!(site_id %in% c(16))) |>  # Exclude non-forest sites
  mutate(
    # Create state abbreviation (extract from ISO code, e.g., "US-ME" -> "ME")
    state_abbrev = ifelse(!is.na(state_code), 
                          sub("^[A-Z]{2}-", "", state_code), 
                          NA_character_)
  ) |>
  select(
    site_id,
    state_abbrev,
    centroid_lat,
    centroid_lon,
    ecoregion,
    dominant_lc,
    area_km2,
    n_chm,
    slope_q25, slope_q75,
    chm_q25, chm_q75,
    elev_q25, elev_q75
  ) |>
  mutate(
    centroid_lat = round(centroid_lat, 2),
    centroid_lon = round(centroid_lon, 2),
    area_km2 = round(area_km2, 0),
    slope_iqr = sprintf("%.1f–%.1f", slope_q25, slope_q75),
    chm_iqr = sprintf("%.1f–%.1f", chm_q25, chm_q75),
    elev_iqr = ifelse(!is.na(elev_q25) & !is.na(elev_q75),
                      sprintf("%.0f–%.0f", elev_q25, elev_q75),
                      NA_character_),
    n_chm_fmt = format(n_chm, big.mark = ",")
  ) |>
  select(
    `Site ID` = site_id,
    `State` = state_abbrev,
    `Lat (°N)` = centroid_lat,
    `Lon (°W)` = centroid_lon,
    `Ecoregion` = ecoregion,
    `Forest Type` = dominant_lc,
    `Area (km²)` = area_km2,
    `n Footprints` = n_chm_fmt,
    `Slope IQR (°)` = slope_iqr,
    `CHM IQR (m)` = chm_iqr,
    `Elev IQR (m)` = elev_iqr
  )

write.csv(table1, 
          file.path(OUTPUT_DIR, "Table1_site_characteristics_v2.csv"),
          row.names = FALSE)

# Save overall stats as text summary
sink(file.path(OUTPUT_DIR, "section2_1_prose_values_v2.txt"))
cat("=" |> rep(70) |> paste(collapse=""), "\n")
cat("Section 2.1 Prose Values (v2 - with DTM IQR)\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("=" |> rep(70) |> paste(collapse=""), "\n\n")

cat("GEOGRAPHIC EXTENT:\n")
cat(sprintf("  Latitude range: %.2f°N to %.2f°N (span: %.2f°)\n",
            overall_stats$lat_range[1], overall_stats$lat_range[2],
            diff(overall_stats$lat_range)))
cat(sprintf("  Longitude range: %.2f°W to %.2f°W (span: %.2f°)\n",
            abs(overall_stats$lon_range[2]), abs(overall_stats$lon_range[1]),
            abs(diff(overall_stats$lon_range))))
cat(sprintf("  Total mapped extent: ~%.0f km²\n\n", overall_stats$total_area_km2))

cat("STATES/PROVINCES REPRESENTED:\n")
cat(sprintf("  Count: %d\n", length(overall_stats$states)))
cat(sprintf("  Names: %s\n\n", paste(sort(overall_stats$states), collapse = ", ")))

cat("ELEVATION:\n")
cat(sprintf("  Range: %.0f m to %.0f m\n\n", 
            overall_stats$elev_range[1], overall_stats$elev_range[2]))

cat("CANOPY HEIGHT (ALS reference, P90):\n")
cat(sprintf("  Range: %.1f m to %.1f m\n\n",
            overall_stats$chm_range[1], overall_stats$chm_range[2]))

cat("SLOPE:\n")
cat(sprintf("  Range: %.1f° to %.1f°\n\n",
            overall_stats$slope_range[1], overall_stats$slope_range[2]))

cat("SAMPLE SIZES:\n")
cat(sprintf("  Total raw footprints: %s\n", format(overall_stats$total_raw, big.mark=",")))
cat(sprintf("  CHM analysis footprints (after QC): %s\n", format(overall_stats$total_chm, big.mark=",")))
cat(sprintf("  DTM analysis footprints (after QC): %s\n\n", format(overall_stats$total_dtm, big.mark=",")))

cat("SITES:\n")
cat(sprintf("  Total forest sites: %d (Site 16 excluded: non-forest)\n", overall_stats$n_sites_total))
cat(sprintf("  Sites in CHM analysis: %d\n", overall_stats$n_sites_chm))
cat(sprintf("  Sites in DTM analysis: %d (Site 10 excluded: datum offset)\n\n", overall_stats$n_sites_dtm))

cat("ECOREGIONS:\n")
cat(sprintf("  Number of EPA Level III ecoregions: %d\n", overall_stats$n_ecoregions))
cat("  Ecoregion names:\n")
for (eco in sort(overall_stats$ecoregions)) {
  cat(sprintf("    - %s\n", eco))
}

cat("\nFOREST TYPES:\n")
cat(sprintf("  Types represented: %s\n", paste(sort(overall_stats$forest_types), collapse=", ")))

cat("\n")
cat("=" |> rep(70) |> paste(collapse=""), "\n")
cat("FORMATTED VALUES FOR PROSE (copy-paste ready)\n")
cat("=" |> rep(70) |> paste(collapse=""), "\n\n")

cat("Paragraph 2 template:\n")
cat(sprintf('  "The %d sites span approximately %.0f degrees of latitude from %.2f°N to %.2f°N\n',
            overall_stats$n_sites_total,
            diff(overall_stats$lat_range),
            overall_stats$lat_range[1], overall_stats$lat_range[2]))
cat(sprintf('   and %.0f degrees of longitude from %.2f°W to %.2f°W, encompassing a total\n',
            abs(diff(overall_stats$lon_range)),
            abs(overall_stats$lon_range[2]), abs(overall_stats$lon_range[1])))
cat(sprintf('   mapped extent of approximately %.0f km². Sites represent %d EPA Level III\n',
            overall_stats$total_area_km2, overall_stats$n_ecoregions))
cat(sprintf('   ecoregions... Elevation ranges from %.0f m to %.0f m across all sites."\n\n',
            overall_stats$elev_range[1], overall_stats$elev_range[2]))

cat("Paragraph 3 template:\n")
cat(sprintf('  "Reference canopy heights measured from airborne LiDAR range from\n'))
cat(sprintf('   approximately %.0f m in recently disturbed stands to over %.0f m in mature\n',
            overall_stats$chm_range[1], overall_stats$chm_range[2]))
cat(sprintf('   Pacific Northwest forests..."\n\n'))

sink()

# Print summary to console
cat("\n")
cat("=" |> rep(70) |> paste(collapse=""), "\n")
cat("SUMMARY (19 forest sites; Site 16 excluded as non-forest)\n")
cat("=" |> rep(70) |> paste(collapse=""), "\n")
cat(sprintf("Total CHM footprints: %s\n", format(overall_stats$total_chm, big.mark=",")))
cat(sprintf("Total DTM footprints: %s\n", format(overall_stats$total_dtm, big.mark=",")))
cat(sprintf("Latitude range: %.2f°N to %.2f°N\n", 
            overall_stats$lat_range[1], overall_stats$lat_range[2]))
cat(sprintf("Longitude range: %.2f°W to %.2f°W\n",
            abs(overall_stats$lon_range[2]), abs(overall_stats$lon_range[1])))
cat(sprintf("States/provinces: %s\n", paste(sort(overall_stats$states), collapse=", ")))
cat(sprintf("Ecoregions: %d\n", overall_stats$n_ecoregions))
cat(sprintf("Elevation range: %.0f m to %.0f m\n", 
            overall_stats$elev_range[1], overall_stats$elev_range[2]))
cat(sprintf("Max canopy height: %.1f m\n", overall_stats$chm_range[2]))
cat(sprintf("Max slope: %.1f°\n\n", overall_stats$slope_range[2]))

cat("Output files saved to:", OUTPUT_DIR, "\n")
cat("  - section2_1_site_summary_full_v2.csv (complete data with elev Q25/Q75)\n")
cat("  - Table1_site_characteristics_v2.csv (formatted for manuscript with DTM IQR)\n")
cat("  - section2_1_prose_values_v2.txt (values for prose)\n\n")

cat("Done!\n")
