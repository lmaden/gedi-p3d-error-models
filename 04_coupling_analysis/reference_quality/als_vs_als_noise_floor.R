#!/usr/bin/env Rscript
# =============================================================================
# als_vs_als_noise_floor.R
#
# The ALS-versus-ALS comparison.
# Paired-year ALS-vs-ALS comparison: compute the irreducible reference-noise
# floor and the "reducible fraction" of P3D-ALS CHM error.
#
# Inputs:
#   - Per-site enriched CSVs (data/enriched_by_site/site_NN_enriched.csv.gz)
#   - Per-site GEDI gpkg footprint geometry (gedi/NN/GEDI_siteNN_hq_ALL.gpkg)
#   - manuscript_tables/site_id_lookup.csv (tracker # <-> manuscript # mapping)
#   - manuscript_tables/als_repeat_acquisition_inventory.csv (multi-year CHM dirs)
#   - ALS CHM rasters on cluster (/gpfs/data1/vclgp/data/gedi/imported/usa/)
#
# Outputs (manuscript_tables/):
#   - als_noise_extraction_summary.csv — per-footprint × per-year p90
#   - als_noise_pairwise_summary.csv   — per-site × per-pair stats
#   - als_noise_reducible_fraction.csv  — the headline table
#
# Outputs (plots/groundwork/):
#   - task5_als_als_distributions.pdf
#   - task5_reducible_fraction.pdf
#   - task5_howland_detail.pdf
#   - task5_rmse_vs_gap.pdf
#
# Design:
#   - Full 25 m-radius p90 extraction matching pipeline kernel method
#   - Chunked by site with per-site checkpointing for resumability
#   - Forest-filtered footprints only (same set as P3D-ALS manuscript analysis)
#   - Growth correction with/without as sensitivity
#   - All pairwise year combinations computed per site
#
# Run: source("als_vs_als_noise_floor.R") from Pane 1 (R console)
# =============================================================================

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 0: Setup
# ─────────────────────────────────────────────────────────────────────────────

source("analysis_config.R")
source("analysis_utils.R")

log_section("The ALS-versus-ALS comparison: Paired-year ALS-vs-ALS comparison")

set.seed(2025)

# Paths
ALS_BASE      <- Sys.getenv("ALS_BASE", "/gpfs/data1/vclgp/data/gedi/imported/usa")
MANUSCRIPT_TB <- file.path(PROJECT_ROOT, "manuscript_tables")
GEDI_BASE     <- file.path(PROJECT_ROOT, "gedi")
ENRICHED_DIR  <- file.path(PROJECT_ROOT, "data", "enriched_by_site")
PLOT_DIR      <- file.path(PROJECT_ROOT, "plots", "groundwork")
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(MANUSCRIPT_TB, recursive = TRUE, showWarnings = FALSE)

# Extraction parameters
BUFFER_M         <- 25      # Match pipeline 25 m GEDI footprint radius
MIN_VALID_PIXELS <- 500     # ~25% of theoretical max (pi*25^2 = 1963 at 1m)
CHUNK_SIZE       <- 1000    # Footprints per tile-scan chunk
FOREST_THRESH_M  <- 2       # als_chm_p90 >= 2 m (matches chm_forest_thresh_m)
VALID_FRAC_THRESH <- 0.5    # Matches pipeline QC threshold
ERROR_OUTLIER_M  <- 100     # abs(chm_error) <= 100 m (matches pipeline)

# Checkpoint name
CKPT_NAME <- "als_noise_extraction"

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: Load reference tables
# ─────────────────────────────────────────────────────────────────────────────

log_subsection("Loading reference tables")

# Site ID lookup
lookup_path <- file.path(MANUSCRIPT_TB, "site_id_lookup.csv")
stopifnot(file.exists(lookup_path))
site_lookup <- fread(lookup_path)
# Keep only the 19 manuscript sites (drop omitted rows with NA manuscript_site)
site_lookup <- site_lookup[!is.na(manuscript_site) & manuscript_site != ""]
site_lookup[, manuscript_site := as.integer(manuscript_site)]
site_lookup[, tracker_site := as.integer(tracker_site)]
log_progress(sprintf("  Site lookup: %d manuscript sites loaded", nrow(site_lookup)))

# ALS inventory
inv_path <- file.path(MANUSCRIPT_TB, "als_repeat_acquisition_inventory.csv")
stopifnot(file.exists(inv_path))
als_inv <- fread(inv_path)
log_progress(sprintf("  ALS inventory: %d rows (site × sibling dir combinations)", nrow(als_inv)))

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Identify multi-year sites and assign years
# ─────────────────────────────────────────────────────────────────────────────

log_subsection("Identifying multi-year sites and assigning acquisition years")

# Robust year extraction from sibling directory name and date_hint.
# Priority: directory-name year (most reliable), then date_hint.
# Handles the known anomaly: konz2017 has date_hint "2000" (regex matched
# coordinate digits, not the actual acquisition year).
extract_year <- function(sibling_dir, date_hint, pipeline_dir) {
  # Strategy 1: Extract trailing 4-digit year from directory name
  # e.g. "neon_jerc2021" → 2021, "neon_harv" → NA
  m <- regmatches(sibling_dir, regexpr("(20[0-2][0-9])$", sibling_dir))
  if (length(m) == 1 && nzchar(m)) {
    yr <- as.integer(m)
    if (!is.na(yr) && yr >= 2008 && yr <= 2025) return(yr)
  }

  # Strategy 2: Extract year embedded mid-name with underscore prefix
  # e.g. "ltbmu_20180710" → 2018 (from full date), "nasa_howland_2009" → 2009
  m2 <- regmatches(sibling_dir, regexpr("_(20[0-2][0-9])", sibling_dir))
  if (length(m2) == 1 && nzchar(m2)) {
    yr <- as.integer(gsub("^_", "", m2))
    if (!is.na(yr) && yr >= 2008 && yr <= 2025) return(yr)
  }

  # Strategy 3: Use date_hint (first 4 chars)
  if (!is.na(date_hint) && nzchar(date_hint)) {
    yr <- as.integer(substr(date_hint, 1, 4))
    # Sanity check: reject implausible years (the konz "2000" case)
    if (!is.na(yr) && yr >= 2008 && yr <= 2025) return(yr)
  }

  # Fallback: NA (will be logged as a warning)
  return(NA_integer_)
}

# Filter inventory to CHM-present rows only
als_inv <- als_inv[chm_subdir_present == "yes"]

# For site 12 (neon_stei2022): drop CHEQ rows (CRS 32615, different site)
# These belong to site 14 (neon_steicheq2022)
als_inv <- als_inv[!(manuscript_site == 12 & grepl("steicheq", sibling_dir))]

# Assign years
als_inv[, acq_year := mapply(extract_year, sibling_dir, date_hint, site_dir_used_in_pipeline)]

# Log any NA years
na_years <- als_inv[is.na(acq_year)]
if (nrow(na_years) > 0) {
  log_progress("  WARNING: Could not assign year to these directories:")
  for (i in seq_len(nrow(na_years))) {
    log_progress(sprintf("    Site %d, dir=%s, date_hint=%s",
                         na_years$manuscript_site[i],
                         na_years$sibling_dir[i],
                         na_years$date_hint[i]))
  }
}

# Count years per site
site_year_count <- als_inv[!is.na(acq_year), .(
  n_years = uniqueN(acq_year),
  years = paste(sort(unique(acq_year)), collapse = ",")
), by = manuscript_site]

multi_year_sites <- site_year_count[n_years >= 2]$manuscript_site
log_progress(sprintf("  Multi-year sites (>= 2 distinct years): %d of 19",
                     length(multi_year_sites)))

# Display site × year matrix
for (ms in sort(multi_year_sites)) {
  yrs <- sort(unique(als_inv[manuscript_site == ms & !is.na(acq_year)]$acq_year))
  lk  <- site_lookup[manuscript_site == ms]
  log_progress(sprintf("  Site %2d (%s): %d years [%s]",
                       ms, lk$site_dir_used_in_pipeline, length(yrs),
                       paste(yrs, collapse = ", ")))
}

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: Literature-based growth rates
# ─────────────────────────────────────────────────────────────────────────────

log_subsection("Growth rate lookup table")

# Conservative midpoint annual height increment rates by dominant forest type.
# Sources: FIA national-scale estimates (Domke et al. 2020 USDA GTR);
# NEON AOP flux tower biomass increments; Pan et al. 2011 Science (US temperate).
# Values are deliberately central estimates; the sensitivity analysis reports
# results at +/- 50% of these rates.
growth_rate_table <- data.table(
  lc_code  = c("ENF", "BDF", "DNF", "EBF", "MFT", "IWL", "GRS", "SHR", "UNK"),
  rate_m_yr = c(0.35,  0.50,  0.35,  0.55,  0.40,  0.30,  0.00,  0.10,  0.00)
)
log_progress("  Growth rates (m/yr):")
for (i in seq_len(nrow(growth_rate_table))) {
  log_progress(sprintf("    %s: %.2f", growth_rate_table$lc_code[i],
                       growth_rate_table$rate_m_yr[i]))
}

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4: Helper functions
# ─────────────────────────────────────────────────────────────────────────────

log_subsection("Defining helper functions")

# Pre-compute circular offset grid (done ONCE, reused for all extractions).
OFFSET_GRID <- {
  og <- expand.grid(dx = seq(-BUFFER_M, BUFFER_M, by = 1),
                    dy = seq(-BUFFER_M, BUFFER_M, by = 1))
  og <- og[sqrt(og$dx^2 + og$dy^2) <= BUFFER_M, ]
  log_progress(sprintf("  Offset grid: %d sample points per footprint (r=%d m, 1 m spacing)",
                       nrow(og), BUFFER_M))
  og
}

#' Read all pixel values from a GeoTIFF tile using command-line GDAL.
#'
#' Converts the tile to raw ENVI binary via gdal_translate, then reads
#' with readBin(). This completely bypasses terra's C++ pixel-reading
#' backend, which segfaults in this environment.
#'
#' @param tif_path Path to the GeoTIFF file.
#' @param n_cells  Expected number of cells (ncol × nrow).
#' @return Numeric vector of length n_cells, or NULL on failure.
read_tile_gdal <- function(tif_path, n_cells) {
  tmp <- tempfile(fileext = ".bin", tmpdir = file.path(PROJECT_ROOT, "tmp"))
  on.exit(unlink(c(tmp, paste0(tmp, ".hdr"), paste0(tmp, ".aux.xml")),
                 force = TRUE), add = TRUE)

  ret <- system2("gdal_translate",
                  c("-q", "-of", "ENVI", "-ot", "Float32",
                    shQuote(tif_path), shQuote(tmp)),
                  stdout = FALSE, stderr = FALSE)

  if (ret != 0 || !file.exists(tmp)) {
    warning("gdal_translate failed for: ", basename(tif_path))
    return(NULL)
  }

  v <- readBin(tmp, what = "double", n = n_cells, size = 4, endian = "little")

  # Clean nodata / sentinel values (CHM valid range: roughly -10 to 100 m)
  v[!is.finite(v) | v < -100 | v > 200] <- NA
  v
}

#' Extract ALS CHM p90 within a 25 m radius at GEDI footprint locations.
#'
#' TILE-BY-TILE CELL-INDEX approach — completely bypasses terra::extract()
#' and terra::vrt(), both of which segfault in this GDAL/terra environment.
#'
#' For each chunk of footprints:
#'   1. Generate all sample-point coordinates on the 1 m offset grid.
#'   2. For each CHM tile: find sample points inside the tile, convert
#'      coordinates to cell indices via cellFromXY(), read values via r[cells].
#'   3. After all tiles: aggregate per footprint to compute p90.
#'
#' @param chm_dir    Path to the directory containing CHM .tif tiles.
#' @param pts_xy     Matrix (n x 2) of footprint centroid coordinates in the
#'                   ALS CRS (same projection as the tiles).
#' @param offset_grid Data.frame with columns dx, dy.
#' @param min_valid  Minimum non-NA pixel count to compute p90.
#' @param chunk_sz   Number of footprints per processing chunk.
#' @return A data.table with columns: idx, p90, n_valid.
extract_chm_p90 <- function(chm_dir, pts_xy, offset_grid = OFFSET_GRID,
                            min_valid = 500, chunk_sz = 1000) {

  tifs <- list.files(chm_dir, pattern = "\\.(tif|tiff)$",
                     full.names = TRUE, ignore.case = TRUE)
  if (length(tifs) == 0) stop("No TIF files in: ", chm_dir)

  n_pts     <- nrow(pts_xy)
  n_offsets <- nrow(offset_grid)
  n_chunks  <- ceiling(n_pts / chunk_sz)
  n_tiles   <- length(tifs)
  results   <- vector("list", n_chunks)

  # Pre-load tile metadata: extent, resolution, dimensions (header reads only)
  log_progress(sprintf("    Reading %d tile headers...", n_tiles))
  tile_meta <- lapply(tifs, function(f) {
    r <- terra::rast(f)
    e <- as.vector(terra::ext(r))
    list(path = f, xmin = e[1], xmax = e[2], ymin = e[3], ymax = e[4],
         res_x = terra::res(r)[1], res_y = terra::res(r)[2],
         ncol = terra::ncol(r), nrow = terra::nrow(r))
  })
  tile_xmin <- vapply(tile_meta, `[[`, 0, "xmin")
  tile_xmax <- vapply(tile_meta, `[[`, 0, "xmax")
  tile_ymin <- vapply(tile_meta, `[[`, 0, "ymin")
  tile_ymax <- vapply(tile_meta, `[[`, 0, "ymax")

  for (ci in seq_len(n_chunks)) {
    start_i    <- (ci - 1L) * chunk_sz + 1L
    end_i      <- min(ci * chunk_sz, n_pts)
    n_in_chunk <- end_i - start_i + 1L
    chunk_xy   <- pts_xy[start_i:end_i, , drop = FALSE]

    # Generate sample-point coordinates for this chunk
    sx    <- rep(chunk_xy[, 1], each = n_offsets) + rep(offset_grid$dx, n_in_chunk)
    sy    <- rep(chunk_xy[, 2], each = n_offsets) + rep(offset_grid$dy, n_in_chunk)
    fp_id <- rep(seq_len(n_in_chunk), each = n_offsets)
    chm_v <- rep(NA_real_, length(fp_id))

    # Scan tiles: read pixel data via gdal_translate + readBin (command-line
    # GDAL, which works reliably), then index in pure R. This completely
    # bypasses terra's C++ pixel-reading backend which segfaults here.
    for (ti in seq_len(n_tiles)) {
      in_tile <- which(
        is.na(chm_v) &
        sx >= tile_xmin[ti] & sx <= tile_xmax[ti] &
        sy >= tile_ymin[ti] & sy <= tile_ymax[ti]
      )
      if (length(in_tile) == 0L) next

      # Read tile via gdal_translate to raw ENVI binary + readBin
      tm <- tile_meta[[ti]]
      all_v <- read_tile_gdal(tm$path, tm$ncol * tm$nrow)
      if (is.null(all_v)) next

      # Pure R coordinate → cell-index conversion
      col_i <- as.integer(floor((sx[in_tile] - tm$xmin) / tm$res_x)) + 1L
      row_i <- as.integer(floor((tm$ymax - sy[in_tile]) / tm$res_y)) + 1L

      ok <- col_i >= 1L & col_i <= tm$ncol & row_i >= 1L & row_i <= tm$nrow
      if (any(ok)) {
        cell_i <- (row_i[ok] - 1L) * tm$ncol + col_i[ok]
        chm_v[in_tile[ok]] <- all_v[cell_i]
      }
      rm(all_v)
    }

    # Aggregate per footprint
    dt <- data.table(fp = fp_id, chm = chm_v)
    agg <- dt[, {
      x  <- chm[!is.na(chm)]
      nv <- length(x)
      list(p90 = if (nv >= min_valid) quantile(x, 0.9, names = FALSE) else NA_real_,
           n_valid = nv)
    }, by = fp]
    agg[, idx := fp + start_i - 1L]
    results[[ci]] <- agg[, .(idx, p90, n_valid)]

    rm(sx, sy, fp_id, chm_v, dt, agg)
    if (ci %% 5 == 0) gc(verbose = FALSE)
    if (ci %% 2 == 0 || ci == n_chunks) {
      log_progress(sprintf("    Chunk %d/%d (%d%%, footprints %d-%d)",
                           ci, n_chunks, round(100 * end_i / n_pts),
                           start_i, end_i))
    }
  }

  rbindlist(results)
}


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5: Main extraction loop — per site, with checkpointing
# ─────────────────────────────────────────────────────────────────────────────

log_section("Main extraction: per-site ALS CHM p90 from all available years")

# Verify gdal_translate is available (required for pixel reading)
if (system2("gdal_translate", "--version", stdout = FALSE, stderr = FALSE) != 0) {
  stop("gdal_translate not found on PATH. Source the tmux-quad-chpt1-v2.sh module stack.")
}
log_progress("  gdal_translate available ✓")

# Expected runtime: ~5-15 min per site × year (I/O bound by VRT tile reads),
# ~60-90 min per site with 5-6 years, total ~6-12 hours for all 15 sites.
# Each site checkpoints on completion, so the loop is fully resumable.

# Load or initialize checkpoint
if (checkpoint_exists(CKPT_NAME)) {
  log_progress("  Resuming from existing checkpoint...")
  ckpt_data <- load_checkpoint(CKPT_NAME)
  extraction_results <- ckpt_data$extraction_results
  completed_sites    <- ckpt_data$completed_sites
  footprint_meta     <- ckpt_data$footprint_meta
  log_progress(sprintf("  Already completed: sites %s",
                       paste(sort(completed_sites), collapse = ", ")))
} else {
  extraction_results <- list()
  completed_sites    <- integer(0)
  footprint_meta     <- list()
}

# Process each multi-year site
for (ms in sort(multi_year_sites)) {

  if (ms %in% completed_sites) {
    log_progress(sprintf("[Site %d] Already completed — skipping", ms))
    next
  }

  lk <- site_lookup[manuscript_site == ms]
  tracker_num  <- lk$tracker_site
  pipeline_dir <- lk$site_dir_used_in_pipeline
  dominant_lc  <- lk$dominant_lc

  log_subsection(sprintf("Site %d (%s) — tracker %d, LC = %s",
                         ms, pipeline_dir, tracker_num, dominant_lc))

  # ── 5a: Load enriched CSV and apply forest filter ──────────────────────

  # Enriched CSVs use zero-padded site numbers: site_02_enriched.csv.gz
  enriched_path <- file.path(ENRICHED_DIR,
                             sprintf("site_%02d_enriched.csv.gz", tracker_num))
  if (!file.exists(enriched_path)) {
    enriched_path <- file.path(ENRICHED_DIR,
                               sprintf("site_%02d_enriched.csv", tracker_num))
  }
  # Fallback: non-padded (site_2 instead of site_02)
  if (!file.exists(enriched_path)) {
    enriched_path <- file.path(ENRICHED_DIR,
                               sprintf("site_%d_enriched.csv.gz", tracker_num))
  }
  if (!file.exists(enriched_path)) {
    enriched_path <- file.path(ENRICHED_DIR,
                               sprintf("site_%d_enriched.csv", tracker_num))
  }
  # Last resort: glob for any matching pattern
  if (!file.exists(enriched_path)) {
    candidates <- list.files(ENRICHED_DIR,
                             pattern = sprintf("site_0?%d_enriched\\.csv", tracker_num),
                             full.names = TRUE)
    if (length(candidates) > 0) enriched_path <- candidates[1]
  }
  if (!file.exists(enriched_path)) {
    log_progress(sprintf("  WARNING: Enriched CSV not found for tracker %d — skipping",
                         tracker_num))
    log_progress(sprintf("           Searched in: %s", ENRICHED_DIR))
    # On first miss, list directory contents to help diagnose
    if (!exists(".enriched_dir_listed")) {
      avail <- list.files(ENRICHED_DIR, pattern = "site_.*enriched")
      log_progress(sprintf("           Available files (%d): %s",
                           length(avail),
                           paste(head(avail, 10), collapse = ", ")))
      if (length(avail) > 10) log_progress("           ...")
      .enriched_dir_listed <<- TRUE
    }
    next
  }

  # Read only the columns we need
  need_cols <- c("shot_number", "site",
                 "p3d_chm_mean", "als_chm_mean", "als_chm_p90",
                 "als_chm_valid_frac", "p3d_chm_valid_frac",
                 "slope_valid_frac", "error_mean",
                 "lc2022_l1_code", "lc2022_mode_l1_code",
                 "rh_98", "cover", "slope_mean")
  hdr <- names(fread(enriched_path, nrows = 0))
  use_cols <- intersect(need_cols, hdr)

  dt_raw <- fread(enriched_path, select = use_cols,
                  colClasses = list(character = "shot_number"),
                  showProgress = FALSE)
  log_progress(sprintf("  Loaded enriched CSV: %s rows", format(nrow(dt_raw), big.mark = ",")))

  # Harmonize LC column name
  if (!("lc2022_l1_code" %in% names(dt_raw)) && "lc2022_mode_l1_code" %in% names(dt_raw)) {
    setnames(dt_raw, "lc2022_mode_l1_code", "lc2022_l1_code")
  }

  # Apply the SAME forest filter as section_01_ingest.R
  dt_forest <- dt_raw[
    is.finite(p3d_chm_mean) & is.finite(als_chm_mean) &
    is.finite(als_chm_valid_frac) & is.finite(p3d_chm_valid_frac) &
    als_chm_valid_frac >= VALID_FRAC_THRESH &
    p3d_chm_valid_frac >= VALID_FRAC_THRESH &
    is.finite(slope_valid_frac) & slope_valid_frac >= VALID_FRAC_THRESH &
    is.finite(als_chm_p90) & als_chm_p90 >= FOREST_THRESH_M
  ]

  # Compute CHM error (P3D - ALS) matching pipeline logic
  dt_forest[, chm_error := fifelse(is.finite(error_mean),
                                   error_mean,
                                   p3d_chm_mean - als_chm_mean)]
  dt_forest <- dt_forest[abs(chm_error) <= ERROR_OUTLIER_M | !is.finite(chm_error)]

  n_forest <- nrow(dt_forest)
  log_progress(sprintf("  After forest filter: %s footprints",
                       format(n_forest, big.mark = ",")))

  if (n_forest < 30) {
    log_progress("  Too few forest footprints (< 30) — skipping site")
    completed_sites <- c(completed_sites, ms)
    next
  }

  # ── 5b: Load GEDI footprint geometry ───────────────────────────────────

  gpkg_path <- file.path(GEDI_BASE, as.character(tracker_num),
                         sprintf("GEDI_site%s_hq_ALL.gpkg", tracker_num))
  if (!file.exists(gpkg_path)) {
    log_progress(sprintf("  WARNING: GEDI gpkg not found: %s — skipping", gpkg_path))
    next
  }

  gedi_sf <- tryCatch(
    suppressMessages(st_read(gpkg_path, quiet = TRUE, int64_as_string = TRUE)),
    error = function(e) suppressMessages(st_read(gpkg_path, quiet = TRUE))
  )

  # Match footprints by shot_number
  gedi_sf$shot_key <- if (is.character(gedi_sf$shot_number)) {
    gedi_sf$shot_number
  } else {
    format(gedi_sf$shot_number, scientific = FALSE, trim = TRUE)
  }
  gedi_sf <- gedi_sf[!is.na(gedi_sf$shot_key), ]
  gedi_sf <- gedi_sf[!duplicated(gedi_sf$shot_key), ]

  dt_forest[, shot_key := as.character(shot_number)]
  matched <- dt_forest$shot_key %in% gedi_sf$shot_key
  n_matched <- sum(matched)
  log_progress(sprintf("  Matched to GEDI geometry: %d / %d footprints",
                       n_matched, n_forest))

  if (n_matched < 30) {
    log_progress("  Too few matched footprints (< 30) — skipping")
    completed_sites <- c(completed_sites, ms)
    next
  }

  dt_forest <- dt_forest[matched]
  gedi_matched <- gedi_sf[match(dt_forest$shot_key, gedi_sf$shot_key), ]

  # Store footprint metadata for later (P3D-ALS error, LC, etc.)
  footprint_meta[[as.character(ms)]] <- dt_forest[, .(
    shot_number, shot_key, als_chm_p90_pipeline = als_chm_p90,
    chm_error, p3d_chm_mean, als_chm_mean,
    lc_code = lc2022_l1_code, rh_98, cover, slope_mean
  )]

  # ── 5c: Determine ALS CRS and reproject footprints ────────────────────

  # Get the CRS from the inventory (pipeline's acquisition)
  site_inv <- als_inv[manuscript_site == ms & !is.na(acq_year)]
  pipeline_row <- site_inv[is_pipeline_dir == "yes"]
  if (nrow(pipeline_row) == 0) pipeline_row <- site_inv[1]

  als_epsg <- as.integer(pipeline_row$crs_epsg[1])
  if (is.na(als_epsg) || als_epsg == 0) {
    log_progress(sprintf("  WARNING: No valid EPSG for site %d — skipping", ms))
    completed_sites <- c(completed_sites, ms)
    next
  }

  # Reproject GEDI footprints to ALS CRS (UTM) and extract coordinate matrix
  pts_sf <- st_transform(gedi_matched, crs = als_epsg)
  pts_xy <- sf::st_coordinates(pts_sf)  # n × 2 matrix (X, Y) in ALS CRS
  log_progress(sprintf("  Footprints reprojected to EPSG:%d", als_epsg))

  # ── 5d: Extract p90 from each available CHM year ──────────────────────

  years_available <- sort(unique(site_inv$acq_year))
  site_extractions <- data.table()

  for (yr in years_available) {
    yr_row <- site_inv[acq_year == yr]
    # If multiple dirs map to same year, pick the one with most tiles
    if (nrow(yr_row) > 1) {
      yr_row <- yr_row[order(-n_tif)][1]
    }

    chm_dir <- file.path(ALS_BASE, yr_row$sibling_dir, "chm")
    log_progress(sprintf("  Year %d: dir=%s (%d tiles)",
                         yr, yr_row$sibling_dir, yr_row$n_tif))

    if (!dir.exists(chm_dir)) {
      log_progress(sprintf("    WARNING: CHM dir not found: %s — skipping year", chm_dir))
      next
    }

    # Check CRS compatibility
    yr_epsg <- as.integer(yr_row$crs_epsg)
    if (!is.na(yr_epsg) && yr_epsg != als_epsg) {
      log_progress(sprintf("    WARNING: CRS mismatch (year=%d, pipeline=%d) — skipping",
                           yr_epsg, als_epsg))
      next
    }

    # Extract p90 via tile-by-tile cell-index lookup (no VRT, no terra::extract)
    t0 <- Sys.time()
    extr <- extract_chm_p90(chm_dir, pts_xy,
                            offset_grid = OFFSET_GRID,
                            min_valid   = MIN_VALID_PIXELS,
                            chunk_sz    = CHUNK_SIZE)
    t_ext <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    n_valid_extr <- sum(!is.na(extr$p90))
    log_progress(sprintf("    Extraction done in %.1f sec: %d / %d footprints valid",
                         t_ext, n_valid_extr, nrow(extr)))

    # Attach metadata
    extr[, `:=`(
      manuscript_site = ms,
      acq_year = yr,
      sibling_dir = yr_row$sibling_dir,
      is_pipeline = yr_row$is_pipeline_dir == "yes",
      shot_key = dt_forest$shot_key[extr$idx]
    )]

    site_extractions <- rbind(site_extractions, extr[, .(
      manuscript_site, shot_key, acq_year, sibling_dir, is_pipeline, p90, n_valid
    )])

    gc(verbose = FALSE)
  }

  # Store results
  extraction_results[[as.character(ms)]] <- site_extractions
  completed_sites <- c(completed_sites, ms)

  # ── 5e: Checkpoint after each site ────────────────────────────────────

  save_checkpoint(CKPT_NAME, list(
    extraction_results = extraction_results,
    completed_sites    = completed_sites,
    footprint_meta     = footprint_meta
  ))

  log_progress(sprintf("[Site %d] Complete. %d year × footprint rows saved.",
                       ms, nrow(site_extractions)))

  # Memory cleanup
  rm(dt_raw, dt_forest, gedi_sf, gedi_matched, pts_sf, pts_xy,
     site_extractions, site_inv)
  gc(verbose = FALSE)
}

log_progress(sprintf("Extraction complete. Sites processed: %s",
                     paste(sort(completed_sites), collapse = ", ")))


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 6: Combine extractions and compute pairwise ALS-ALS differences
# ─────────────────────────────────────────────────────────────────────────────

log_section("Computing pairwise ALS-ALS differences")

# Reload checkpoint if needed (e.g., running sections independently)
if (!exists("extraction_results") || length(extraction_results) == 0) {
  ckpt_data <- load_checkpoint(CKPT_NAME)
  extraction_results <- ckpt_data$extraction_results
  completed_sites    <- ckpt_data$completed_sites
  footprint_meta     <- ckpt_data$footprint_meta
}

# Combine all site extractions
all_extr <- rbindlist(extraction_results)
log_progress(sprintf("  Total extraction rows: %s",
                     format(nrow(all_extr), big.mark = ",")))

# Write extraction summary CSV
fwrite(all_extr, file.path(MANUSCRIPT_TB, "als_noise_extraction_summary.csv"))
log_progress("  Wrote extraction summary CSV")

# For each site, compute all pairwise year differences
pairwise_results <- list()

for (ms in sort(unique(all_extr$manuscript_site))) {

  site_extr <- all_extr[manuscript_site == ms & !is.na(p90)]
  years <- sort(unique(site_extr$acq_year))

  if (length(years) < 2) {
    log_progress(sprintf("  Site %d: < 2 years with valid data — skipping pairwise", ms))
    next
  }

  # Pivot to wide: one row per footprint, one column per year
  wide <- dcast(site_extr, shot_key ~ acq_year, value.var = "p90")
  year_cols <- as.character(years)

  # Get footprint metadata (LC, P3D-ALS error)
  meta <- footprint_meta[[as.character(ms)]]
  wide <- merge(wide, meta, by = "shot_key", all.x = TRUE)

  # Get growth rate for dominant LC
  lk <- site_lookup[manuscript_site == ms]
  site_lc <- lk$dominant_lc

  # Per-footprint LC for growth correction (use site-level if footprint LC missing)
  wide[, growth_lc := fifelse(
    lc_code %in% growth_rate_table$lc_code,
    lc_code,
    site_lc
  )]
  wide <- merge(wide, growth_rate_table, by.x = "growth_lc", by.y = "lc_code",
                all.x = TRUE, sort = FALSE)
  wide[is.na(rate_m_yr), rate_m_yr := 0]

  # Generate all ordered year pairs
  for (i in seq_along(years)[-length(years)]) {
    for (j in (i+1):length(years)) {
      yr_a <- years[i]
      yr_b <- years[j]
      col_a <- as.character(yr_a)
      col_b <- as.character(yr_b)

      if (!(col_a %in% names(wide)) || !(col_b %in% names(wide))) next

      # Footprints with valid p90 from BOTH years
      pair_dt <- wide[!is.na(get(col_a)) & !is.na(get(col_b))]
      n_pair <- nrow(pair_dt)

      if (n_pair < 10) next

      gap_years <- yr_b - yr_a
      raw_diff  <- pair_dt[[col_b]] - pair_dt[[col_a]]

      # Growth correction: expected_growth = rate × gap
      expected_growth <- pair_dt$rate_m_yr * gap_years
      corrected_diff  <- raw_diff - expected_growth

      # Also compute corrected at ±50% growth rate (sensitivity)
      corrected_lo <- raw_diff - expected_growth * 0.5
      corrected_hi <- raw_diff - expected_growth * 1.5

      # Summary statistics
      rmse_raw  <- sqrt(mean(raw_diff^2))
      rmse_corr <- sqrt(mean(corrected_diff^2))
      rmse_lo   <- sqrt(mean(corrected_lo^2))
      rmse_hi   <- sqrt(mean(corrected_hi^2))
      mae_raw   <- mean(abs(raw_diff))
      mae_corr  <- mean(abs(corrected_diff))
      bias_raw  <- mean(raw_diff)
      bias_corr <- mean(corrected_diff)
      sd_raw    <- sd(raw_diff)
      sd_corr   <- sd(corrected_diff)

      # P3D-ALS RMSE on the SAME footprint set (apples-to-apples)
      p3d_als_err <- pair_dt$chm_error
      p3d_als_err <- p3d_als_err[is.finite(p3d_als_err)]
      p3d_als_rmse_same <- if (length(p3d_als_err) >= 10) sqrt(mean(p3d_als_err^2)) else NA_real_
      p3d_als_mae_same  <- if (length(p3d_als_err) >= 10) mean(abs(p3d_als_err)) else NA_real_
      p3d_als_bias_same <- if (length(p3d_als_err) >= 10) mean(p3d_als_err) else NA_real_

      # Reducible fraction (growth-corrected)
      reducible <- if (!is.na(p3d_als_rmse_same) && p3d_als_rmse_same > 0 && rmse_corr > 0) {
        1 - (rmse_corr / p3d_als_rmse_same)^2
      } else NA_real_

      pairwise_results[[length(pairwise_results) + 1]] <- data.table(
        manuscript_site = ms,
        site_name = lk$site_dir_used_in_pipeline,
        dominant_lc = site_lc,
        flag_status = lk$flag_status,
        year_a = yr_a,
        year_b = yr_b,
        gap_years = gap_years,
        n_footprints = n_pair,
        als_als_rmse_raw = round(rmse_raw, 3),
        als_als_mae_raw = round(mae_raw, 3),
        als_als_bias_raw = round(bias_raw, 3),
        als_als_sd_raw = round(sd_raw, 3),
        als_als_rmse_corr = round(rmse_corr, 3),
        als_als_mae_corr = round(mae_corr, 3),
        als_als_bias_corr = round(bias_corr, 3),
        als_als_sd_corr = round(sd_corr, 3),
        als_als_rmse_corr_lo = round(rmse_lo, 3),
        als_als_rmse_corr_hi = round(rmse_hi, 3),
        growth_rate_m_yr = round(unique(pair_dt$rate_m_yr)[1], 2),
        p3d_als_rmse_same_footprints = round(p3d_als_rmse_same, 3),
        p3d_als_mae_same_footprints = round(p3d_als_mae_same, 3),
        p3d_als_bias_same_footprints = round(p3d_als_bias_same, 3),
        reducible_fraction = round(reducible, 3)
      )
    }
  }
}

pairwise_dt <- rbindlist(pairwise_results)
log_progress(sprintf("  Total year-pairs computed: %d across %d sites",
                     nrow(pairwise_dt), uniqueN(pairwise_dt$manuscript_site)))

# Write pairwise summary
fwrite(pairwise_dt, file.path(MANUSCRIPT_TB, "als_noise_pairwise_summary.csv"))
log_progress("  Wrote pairwise summary CSV")


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 7: Build the reducible-fraction summary table
# ─────────────────────────────────────────────────────────────────────────────

log_section("Building reducible-fraction summary table")

# For each site, select the SHORTEST-GAP pair as the primary noise-floor
# estimate (minimizes growth-correction uncertainty), and also report
# the full range across all pairs.

summary_rows <- list()

for (ms in sort(unique(pairwise_dt$manuscript_site))) {
  site_pairs <- pairwise_dt[manuscript_site == ms]
  shortest   <- site_pairs[gap_years == min(gap_years)]
  # If multiple pairs have the same shortest gap, take the one with most footprints
  shortest   <- shortest[order(-n_footprints)][1]

  summary_rows[[length(summary_rows) + 1]] <- data.table(
    manuscript_site = ms,
    site_name = shortest$site_name,
    dominant_lc = shortest$dominant_lc,
    flag_status = shortest$flag_status,
    n_pairs_total = nrow(site_pairs),
    primary_year_a = shortest$year_a,
    primary_year_b = shortest$year_b,
    primary_gap_years = shortest$gap_years,
    primary_n_footprints = shortest$n_footprints,
    primary_als_als_rmse_raw = shortest$als_als_rmse_raw,
    primary_als_als_rmse_corr = shortest$als_als_rmse_corr,
    primary_als_als_rmse_corr_range = sprintf("[%.2f, %.2f]",
      shortest$als_als_rmse_corr_lo, shortest$als_als_rmse_corr_hi),
    p3d_als_rmse_same = shortest$p3d_als_rmse_same_footprints,
    reducible_fraction = shortest$reducible_fraction,
    all_pairs_rmse_corr_range = sprintf("[%.2f, %.2f]",
      min(site_pairs$als_als_rmse_corr), max(site_pairs$als_als_rmse_corr)),
    all_pairs_reducible_range = sprintf("[%.2f, %.2f]",
      min(site_pairs$reducible_fraction, na.rm = TRUE),
      max(site_pairs$reducible_fraction, na.rm = TRUE))
  )
}

summary_dt <- rbindlist(summary_rows)

# Print summary table to console
log_progress("")
log_progress("╔══════════════════════════════════════════════════════════════════════════╗")
log_progress("║  REDUCIBLE FRACTION SUMMARY (shortest-gap pair per site)               ║")
log_progress("╚══════════════════════════════════════════════════════════════════════════╝")
log_progress(sprintf("  %-4s %-20s %-4s %-8s %-7s %-8s %-8s %-8s %-8s",
                     "Site", "Name", "LC", "Flag", "Gap(yr)", "ALS-ALS", "P3D-ALS",
                     "Red.Fr.", "N"))
log_progress(paste(rep("─", 85), collapse = ""))

for (i in seq_len(nrow(summary_dt))) {
  r <- summary_dt[i]
  log_progress(sprintf("  %-4d %-20s %-4s %-8s %-7d %-8.2f %-8.2f %-8.2f %-8d",
    r$manuscript_site, substr(r$site_name, 1, 20), r$dominant_lc,
    substr(r$flag_status, 1, 8), r$primary_gap_years,
    r$primary_als_als_rmse_corr, r$p3d_als_rmse_same,
    r$reducible_fraction, r$primary_n_footprints))
}
log_progress(paste(rep("─", 85), collapse = ""))

# Write summary CSV
fwrite(summary_dt, file.path(MANUSCRIPT_TB, "als_noise_reducible_fraction.csv"))
log_progress("  Wrote reducible-fraction summary CSV")


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 8: Diagnostic figures
# ─────────────────────────────────────────────────────────────────────────────

log_section("Generating diagnostic figures")

png_dev <- get_png_device()

# ── Figure 1: ALS-ALS RMSE (growth-corrected) vs P3D-ALS RMSE bar chart ─────

p_bars <- ggplot(summary_dt, aes(x = reorder(site_name, manuscript_site))) +
  geom_col(aes(y = p3d_als_rmse_same, fill = "P3D-ALS RMSE"), alpha = 0.7, width = 0.6) +
  geom_col(aes(y = primary_als_als_rmse_corr, fill = "ALS-ALS RMSE\n(growth-corrected)"),
           alpha = 0.9, width = 0.4) +
  geom_text(aes(y = p3d_als_rmse_same + 0.3,
                label = sprintf("%.0f%%", 100 * reducible_fraction)),
            size = 2.8, hjust = 0.5) +
  scale_fill_manual(values = c("P3D-ALS RMSE" = "#2166ac",
                               "ALS-ALS RMSE\n(growth-corrected)" = "#b2182b")) +
  coord_flip() +
  labs(x = NULL, y = "RMSE (m)", fill = NULL,
       title = "ALS-ALS noise floor vs P3D-ALS error",
       subtitle = "Shortest-gap pair per site. Percentage = reducible fraction.") +
  theme_cowplot(12) +
  theme(legend.position = "top")

ggsave(file.path(PLOT_DIR, "task5_reducible_fraction.pdf"),
       p_bars, width = 10, height = 8, bg = "white")
log_progress("  Wrote task5_reducible_fraction.pdf")

# ── Figure 2: RMSE vs time gap (all pairs) ──────────────────────────────────

p_gap <- ggplot(pairwise_dt, aes(x = gap_years, y = als_als_rmse_corr)) +
  geom_point(aes(color = dominant_lc, size = n_footprints), alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "grey30", linetype = "dashed") +
  facet_wrap(~dominant_lc, scales = "free_y") +
  scale_size_continuous(range = c(1, 5), guide = "none") +
  labs(x = "Time gap (years)", y = "ALS-ALS RMSE (growth-corrected, m)",
       color = "Forest type",
       title = "ALS-ALS RMSE vs acquisition time gap") +
  theme_cowplot(11) +
  theme(legend.position = "top")

ggsave(file.path(PLOT_DIR, "task5_rmse_vs_gap.pdf"),
       p_gap, width = 10, height = 7, bg = "white")
log_progress("  Wrote task5_rmse_vs_gap.pdf")

# ── Figure 3: ALS-ALS difference distributions (shortest-gap, all sites) ────

# Rebuild per-footprint differences for shortest-gap pairs
diff_for_plot <- list()

for (ms in sort(unique(all_extr$manuscript_site))) {
  site_extr <- all_extr[manuscript_site == ms & !is.na(p90)]
  years <- sort(unique(site_extr$acq_year))
  if (length(years) < 2) next

  meta <- footprint_meta[[as.character(ms)]]
  lk <- site_lookup[manuscript_site == ms]

  # Find shortest gap
  best_gap <- Inf; best_a <- NA; best_b <- NA
  for (i in seq_along(years)[-length(years)]) {
    for (j in (i+1):length(years)) {
      g <- years[j] - years[i]
      if (g < best_gap) { best_gap <- g; best_a <- years[i]; best_b <- years[j] }
    }
  }

  wide <- dcast(site_extr, shot_key ~ acq_year, value.var = "p90")
  col_a <- as.character(best_a)
  col_b <- as.character(best_b)
  if (!(col_a %in% names(wide)) || !(col_b %in% names(wide))) next

  pair <- wide[!is.na(get(col_a)) & !is.na(get(col_b))]
  if (nrow(pair) < 10) next

  rate <- growth_rate_table[lc_code == lk$dominant_lc]$rate_m_yr
  if (length(rate) == 0) rate <- 0
  raw_d <- pair[[col_b]] - pair[[col_a]]
  corr_d <- raw_d - rate * best_gap

  diff_for_plot[[length(diff_for_plot) + 1]] <- data.table(
    manuscript_site = ms,
    site_label = sprintf("Site %d (%s)", ms, substr(lk$site_dir_used_in_pipeline, 1, 14)),
    diff_corr = corr_d,
    diff_raw = raw_d,
    gap = best_gap
  )
}

diff_all <- rbindlist(diff_for_plot)

if (nrow(diff_all) > 0) {
  p_dist <- ggplot(diff_all, aes(y = site_label, x = diff_corr)) +
    geom_violin(fill = "#4393c3", alpha = 0.5, scale = "width") +
    stat_summary(fun = median, geom = "point", size = 2, color = "black") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    coord_cartesian(xlim = c(-10, 10)) +
    labs(x = "ALS year-B minus ALS year-A (growth-corrected, m)",
         y = NULL,
         title = "ALS-ALS difference distributions (shortest-gap pair per site)",
         subtitle = "Black dot = median. Dashed line = zero.") +
    theme_cowplot(11)

  ggsave(file.path(PLOT_DIR, "task5_als_als_distributions.pdf"),
         p_dist, width = 10, height = max(6, 0.5 * length(diff_for_plot)), bg = "white")
  log_progress("  Wrote task5_als_als_distributions.pdf")
} else {
  log_progress("  WARNING: No data available for distribution plot")
}

# ── Figure 4: NASA Howland detail (if available) ─────────────────────────────

if (2 %in% unique(pairwise_dt$manuscript_site)) {
  howl <- pairwise_dt[manuscript_site == 2]
  log_progress("  NASA Howland (FLAGGED, site 2) detail:")
  for (i in seq_len(nrow(howl))) {
    log_progress(sprintf("    %d vs %d: ALS-ALS RMSE = %.2f m (raw), %.2f m (corr), P3D-ALS RMSE = %.2f m, reducible = %.1f%%, n = %d",
      howl$year_a[i], howl$year_b[i],
      howl$als_als_rmse_raw[i], howl$als_als_rmse_corr[i],
      howl$p3d_als_rmse_same_footprints[i],
      100 * howl$reducible_fraction[i],
      howl$n_footprints[i]))
  }

  # Scatter: pipeline-year p90 vs alternate-year p90 at Howland
  howl_extr <- all_extr[manuscript_site == 2 & !is.na(p90)]
  if (nrow(howl_extr) > 0 && length(unique(howl_extr$acq_year)) >= 2) {
    howl_wide <- dcast(howl_extr, shot_key ~ acq_year, value.var = "p90")
    yr_cols <- setdiff(names(howl_wide), "shot_key")
    if (length(yr_cols) == 2) {
      setnames(howl_wide, yr_cols, c("year_a", "year_b"))
      howl_wide <- howl_wide[!is.na(year_a) & !is.na(year_b)]

      p_howl <- ggplot(howl_wide, aes(x = year_a, y = year_b)) +
        geom_point(alpha = 0.4, size = 1) +
        geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
        coord_equal() +
        labs(x = sprintf("ALS CHM p90 — %s (m)", yr_cols[1]),
             y = sprintf("ALS CHM p90 — %s (m)", yr_cols[2]),
             title = "NASA Howland: ALS p90 year-to-year comparison",
             subtitle = sprintf("(FLAGGED site, offset = -7.92 m, n = %d)", nrow(howl_wide))) +
        theme_cowplot(11)

      ggsave(file.path(PLOT_DIR, "task5_howland_detail.pdf"),
             p_howl, width = 7, height = 7, bg = "white")
      log_progress("  Wrote task5_howland_detail.pdf")
    }
  }
} else {
  log_progress("  NASA Howland (site 2) not in results — no detail plot")
}


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 9: Final summary and next steps
# ─────────────────────────────────────────────────────────────────────────────

log_section("The ALS-versus-ALS comparison — Complete")

log_progress("Output files:")
log_progress(sprintf("  %s/als_noise_extraction_summary.csv", MANUSCRIPT_TB))
log_progress(sprintf("  %s/als_noise_pairwise_summary.csv", MANUSCRIPT_TB))
log_progress(sprintf("  %s/als_noise_reducible_fraction.csv", MANUSCRIPT_TB))
log_progress(sprintf("  %s/task5_reducible_fraction.pdf", PLOT_DIR))
log_progress(sprintf("  %s/task5_rmse_vs_gap.pdf", PLOT_DIR))
log_progress(sprintf("  %s/task5_als_als_distributions.pdf", PLOT_DIR))
log_progress(sprintf("  %s/task5_howland_detail.pdf", PLOT_DIR))
log_progress("")
log_progress("Key outputs:")
log_progress("  1. als_noise_reducible_fraction.csv")
log_progress("  2. als_noise_pairwise_summary.csv")
log_progress("  3. The four figures above")
log_progress("")
log_progress("Interpretation guidance (read AFTER results are in):")
log_progress("  - reducible_fraction > 0.80: ALS noise floor is small relative to")
log_progress("    P3D error → most P3D-ALS CHM error is genuinely P3D.")
log_progress("  - reducible_fraction < 0.50: ALS noise floor is a major component →")
log_progress("    site-level CHM variance partly reflects reference instability.")
log_progress("  - If Howland's ALS-ALS RMSE is large (>3 m): reference instability")
log_progress("    at flagged sites is persistent across acquisitions.")
log_progress("  - If ALS-ALS RMSE is small at non-flagged sites but large at Howland:")
log_progress("    confirms the flagging diagnostic and strengthens the 16-site narrative.")
