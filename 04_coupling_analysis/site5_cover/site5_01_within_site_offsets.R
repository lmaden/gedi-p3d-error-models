#!/usr/bin/env Rscript
# =============================================================================
# site5_01_within_site_offsets.R
#
# Step 01.
# Within-site spatial structure of the ALS reference quality offset
# (als_chm_p90 - rh_98) at the three flagged sites, plus an unflagged
# control (HARV) for visual contrast.
#
# Goal:
#   step 01 asks whether the systematic offsets at the flagged 3
#   (Site 1 usda_me, Site 2 nasa_howland, Site 3 neon_sawb) are uniform
#   across each site (= global ALS pipeline issue) or patchy (= dense-
#   canopy interpolation, edge effects, slope/LC dependence).
#
#   Within-site analog of Track A's Supp Fig S23 cross-site diagnostic.
#
#   step 04 (16-site atlas) is gated on step 01 finding spatial structure
#   worth generalizing to non-flagged sites for the residual ~12% V_site.
#
# Approach: footprint-only. No raster reads. Per-footprint als_chm_p90
# (pipeline 25 m kernel mean) and rh_98 are already in the enriched CSVs;
# we tile, aggregate, test for spatial structure.
#
# Inputs:
#   - data/enriched_by_site/site_NN_enriched.csv.gz (tracker NN)
#   - gedi/N/GEDI_siteN_hq_ALL.gpkg (footprint geometry)
#   - manuscript_tables/site_id_lookup.csv
#   - manuscript_tables/als_repeat_acquisition_inventory.csv (for ALS CRS)
#
# Outputs (manuscript_tables/):
#   - groundwork_task6_phase1_tile_summary.csv      per-site x per-tile
#   - groundwork_task6_phase1_site_summary.csv      per-site spatial metrics
#   - groundwork_task6_phase1_regression.csv        per-site OLS coefs
#
# Outputs (plots/groundwork/):
#   - task6_phase1_siteNN_<short>_detail.pdf        per-site 4-panel
#   - task6_phase1_offset_maps_combined.pdf         side-by-side maps
#   - task6_phase1_offset_distributions.pdf         tile-offset density
#
# Run: source("site5_01_within_site_offsets.R") from Pane 1.
#
# Tunables via env:
#   TASK6_TILE_M=250            (default; try 500 for sensitivity)
#   TASK6_INCLUDE_CONTROL=TRUE  (default; FALSE skips HARV)
#   TASK6_MIN_FP=5              (min footprints per tile)
# =============================================================================

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 0: Setup
# ─────────────────────────────────────────────────────────────────────────────

source("analysis_config.R")
source("analysis_utils.R")

log_section("Step 01: Within-site spatial structure of the ALS offset")

set.seed(2026)

suppressPackageStartupMessages({
  library(data.table)
  library(sf)
  library(ggplot2)
  library(cowplot)
  library(spdep)
})

# Configuration
TILE_SIZE_M             <- as.numeric(Sys.getenv("TASK6_TILE_M", "250"))
MIN_FOOTPRINTS_PER_TILE <- as.integer(Sys.getenv("TASK6_MIN_FP", "5"))
INCLUDE_CONTROL         <- toupper(Sys.getenv("TASK6_INCLUDE_CONTROL", "TRUE")) %in%
                            c("TRUE", "T", "1", "YES")
N_PERM                  <- 999
FOREST_THRESH_M         <- 2
VALID_FRAC_THRESH       <- 0.5
ERROR_OUTLIER_M         <- 100

# Site list (tracker numbers; manuscript == tracker for sites 1-4)
FLAGGED_SITES <- c(1L, 2L, 3L)
CONTROL_SITES <- if (INCLUDE_CONTROL) 4L else integer(0)
SITES_PHASE1  <- c(FLAGGED_SITES, CONTROL_SITES)

log_progress(sprintf("  Tile size:           %g m", TILE_SIZE_M))
log_progress(sprintf("  Min footprints/tile: %d",   MIN_FOOTPRINTS_PER_TILE))
log_progress(sprintf("  Sites (tracker #):   %s",
                     paste(SITES_PHASE1, collapse = ", ")))
log_progress(sprintf("  Permutations for Moran's I: %d", N_PERM))

# Paths
ALS_BASE      <- Sys.getenv("ALS_BASE", "/gpfs/data1/vclgp/data/gedi/imported/usa")
MANUSCRIPT_TB <- file.path(PROJECT_ROOT, "manuscript_tables")
GEDI_BASE     <- file.path(PROJECT_ROOT, "gedi")
ENRICHED_DIR  <- file.path(PROJECT_ROOT, "data", "enriched_by_site")
PLOT_DIR      <- file.path(PROJECT_ROOT, "plots", "groundwork")
dir.create(PLOT_DIR,      recursive = TRUE, showWarnings = FALSE)
dir.create(MANUSCRIPT_TB, recursive = TRUE, showWarnings = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: Load reference tables
# ─────────────────────────────────────────────────────────────────────────────

log_subsection("Loading reference tables")

site_lookup <- fread(file.path(MANUSCRIPT_TB, "site_id_lookup.csv"))
site_lookup <- site_lookup[!is.na(manuscript_site) & manuscript_site != ""]
site_lookup[, manuscript_site := as.integer(manuscript_site)]
site_lookup[, tracker_site    := as.integer(tracker_site)]

als_inv <- fread(file.path(MANUSCRIPT_TB, "als_repeat_acquisition_inventory.csv"))
als_inv <- als_inv[chm_subdir_present == "yes"]

log_progress(sprintf("  Site lookup: %d sites", nrow(site_lookup)))

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Helper functions
# ─────────────────────────────────────────────────────────────────────────────

log_subsection("Defining helper functions")

# Locate enriched CSV with zero-padded preference (same fallback ladder as the ALS inventory)
find_enriched <- function(tracker_num) {
  candidates <- c(
    sprintf("site_%02d_enriched.csv.gz", tracker_num),
    sprintf("site_%02d_enriched.csv",    tracker_num),
    sprintf("site_%d_enriched.csv.gz",   tracker_num),
    sprintf("site_%d_enriched.csv",      tracker_num)
  )
  for (cand in candidates) {
    p <- file.path(ENRICHED_DIR, cand)
    if (file.exists(p)) return(p)
  }
  hits <- list.files(ENRICHED_DIR,
                     pattern = sprintf("site_0?%d_enriched\\.csv", tracker_num),
                     full.names = TRUE)
  if (length(hits) > 0) return(hits[1])
  stop("Enriched CSV not found for tracker site ", tracker_num)
}

# Read enriched CSV with forest filter; return offset-bearing data.table
read_enriched_forest <- function(tracker_num) {
  enriched_path <- find_enriched(tracker_num)
  log_progress(sprintf("  Reading: %s", basename(enriched_path)))

  need_cols <- c("shot_number",
                 "p3d_chm_mean", "als_chm_mean", "als_chm_p90",
                 "als_chm_valid_frac", "p3d_chm_valid_frac",
                 "slope_valid_frac", "error_mean",
                 "lc2022_l1_code", "lc2022_mode_l1_code",
                 "rh_98", "cover", "slope_mean")
  hdr <- names(fread(enriched_path, nrows = 0))
  use_cols <- intersect(need_cols, hdr)
  missing_critical <- setdiff(c("shot_number", "als_chm_p90", "rh_98"), use_cols)
  if (length(missing_critical) > 0) {
    stop("Missing critical columns in ", basename(enriched_path), ": ",
         paste(missing_critical, collapse = ", "))
  }

  dt <- fread(enriched_path, select = use_cols,
              colClasses = list(character = "shot_number"),
              showProgress = FALSE)
  log_progress(sprintf("    %s rows loaded", format(nrow(dt), big.mark = ",")))

  # Harmonize LC column. NB: lc2022_l1_code uses lowercase-L plus digit-1
  # for "Level 1"; visually indistinguishable from lc2022_ll_code (double-L).
  if (!("lc2022_l1_code" %in% names(dt)) && "lc2022_mode_l1_code" %in% names(dt)) {
    setnames(dt, "lc2022_mode_l1_code", "lc2022_l1_code")
  }
  if (!("lc2022_l1_code" %in% names(dt))) dt[, lc2022_l1_code := NA_character_]

  # Forest filter (matches the pipeline conventions)
  dt_f <- dt[
    is.finite(p3d_chm_mean) & is.finite(als_chm_mean) &
    is.finite(als_chm_valid_frac) & is.finite(p3d_chm_valid_frac) &
    is.finite(slope_valid_frac) &
    als_chm_valid_frac >= VALID_FRAC_THRESH &
    p3d_chm_valid_frac >= VALID_FRAC_THRESH &
    slope_valid_frac >= VALID_FRAC_THRESH &
    is.finite(als_chm_p90) & als_chm_p90 >= FOREST_THRESH_M &
    is.finite(rh_98)
  ]

  # CHM error: prefer error_mean (the unprefixed pipeline column for CHM, not
  # chm_error_mean), fall back to direct subtraction.
  dt_f[, chm_err := fifelse(is.finite(error_mean), error_mean,
                             p3d_chm_mean - als_chm_mean)]
  dt_f <- dt_f[!is.finite(chm_err) | abs(chm_err) <= ERROR_OUTLIER_M]

  # Offset = ALS p90 - GEDI rh_98 (Track A convention; negative = ALS underestimates)
  dt_f[, offset := als_chm_p90 - rh_98]

  dt_f[, shot_key := as.character(shot_number)]
  log_progress(sprintf("    %s after forest filter", format(nrow(dt_f), big.mark = ",")))
  dt_f
}

# Load GEDI footprint geometry, project to ALS CRS, return matched xy + shot_key
load_gedi_xy <- function(tracker_num, als_epsg, target_keys) {
  gpkg_path <- file.path(GEDI_BASE, as.character(tracker_num),
                         sprintf("GEDI_site%s_hq_ALL.gpkg", tracker_num))
  if (!file.exists(gpkg_path)) stop("GEDI gpkg not found: ", gpkg_path)

  gedi_sf <- tryCatch(
    suppressMessages(st_read(gpkg_path, quiet = TRUE, int64_as_string = TRUE)),
    error = function(e) suppressMessages(st_read(gpkg_path, quiet = TRUE))
  )

  gedi_sf$shot_key <- if (is.character(gedi_sf$shot_number)) {
    gedi_sf$shot_number
  } else {
    format(gedi_sf$shot_number, scientific = FALSE, trim = TRUE)
  }
  gedi_sf <- gedi_sf[!is.na(gedi_sf$shot_key) & !duplicated(gedi_sf$shot_key), ]

  # Restrict to keys in target (forest-filtered) set
  gedi_sf <- gedi_sf[gedi_sf$shot_key %in% target_keys, ]
  if (nrow(gedi_sf) == 0L) stop("No matching GEDI footprints for tracker ", tracker_num)

  pts_proj <- st_transform(gedi_sf, crs = als_epsg)
  xy <- sf::st_coordinates(pts_proj)
  data.table(shot_key = pts_proj$shot_key, x = xy[, 1], y = xy[, 2])
}

# Assign tile (col, row, id) given an origin and tile size.
# Origin is the SW corner snapped down to the nearest tile-size multiple.
assign_tile <- function(x, y, tile_m) {
  origin_x <- floor(min(x) / tile_m) * tile_m
  origin_y <- floor(min(y) / tile_m) * tile_m
  col <- as.integer(floor((x - origin_x) / tile_m))
  row <- as.integer(floor((y - origin_y) / tile_m))
  list(
    tile_col = col,
    tile_row = row,
    tile_id  = sprintf("c%04d_r%04d", col, row),
    origin_x = origin_x,
    origin_y = origin_y
  )
}

# Per-tile distance from centroid to the convex hull boundary of the
# full footprint cloud at this site. Negative if the centroid falls
# outside the hull (shouldn't happen if tiles only exist where footprints do).
tile_edge_distance <- function(tile_dt, all_xy_sf) {
  hull <- st_convex_hull(st_union(all_xy_sf))
  hull_line <- st_cast(hull, "LINESTRING")
  centroids <- st_as_sf(tile_dt, coords = c("tile_centroid_x", "tile_centroid_y"),
                         crs = st_crs(all_xy_sf))
  d <- as.numeric(st_distance(centroids, hull_line))
  inside <- as.logical(st_intersects(centroids, hull, sparse = FALSE))
  ifelse(inside, d, -d)
}

# Moran's I via spdep using distance-band rook contiguity on tile centroids.
moran_on_tiles <- function(tile_dt, value_col, tile_m, n_perm = 999) {
  v <- tile_dt[[value_col]]
  ok <- is.finite(v)
  if (sum(ok) < 4) return(list(I = NA_real_, p_perm = NA_real_, n_used = sum(ok)))

  sub <- tile_dt[ok]
  v_sub <- v[ok]

  coords <- as.matrix(sub[, .(tile_centroid_x, tile_centroid_y)])
  # Rook contiguity: neighbors share an edge => centroids ~ tile_m apart.
  # Use 1.05*tile_m to allow floating-point slop; exclude diagonals
  # (which are tile_m * sqrt(2) ≈ 1.414 * tile_m apart).
  nb <- tryCatch(
    dnearneigh(coords, d1 = 0, d2 = 1.05 * tile_m),
    error = function(e) NULL
  )
  if (is.null(nb)) return(list(I = NA_real_, p_perm = NA_real_, n_used = nrow(sub)))

  # Drop isolates
  card <- card(nb)
  keep <- card > 0
  if (sum(keep) < 4) return(list(I = NA_real_, p_perm = NA_real_, n_used = sum(keep)))

  nb2     <- subset(nb, keep)
  v_keep  <- v_sub[keep]
  lw      <- tryCatch(nb2listw(nb2, style = "W", zero.policy = TRUE),
                      error = function(e) NULL)
  if (is.null(lw)) return(list(I = NA_real_, p_perm = NA_real_, n_used = sum(keep)))

  mc <- tryCatch(
    moran.mc(v_keep, lw, nsim = n_perm, zero.policy = TRUE,
             alternative = "two.sided"),
    error = function(e) NULL
  )
  if (is.null(mc)) return(list(I = NA_real_, p_perm = NA_real_, n_used = sum(keep)))

  list(I = unname(mc$statistic), p_perm = mc$p.value, n_used = sum(keep))
}

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: Per-site processing
# ─────────────────────────────────────────────────────────────────────────────

log_section("Per-site tile aggregation and spatial-structure tests")

all_tile_rows  <- list()
all_site_rows  <- list()
all_reg_rows   <- list()
all_plot_data  <- list()

for (tnum in SITES_PHASE1) {

  lk <- site_lookup[tracker_site == tnum]
  if (nrow(lk) == 0L) {
    log_progress(sprintf("  Tracker %d not in lookup; skipping", tnum))
    next
  }
  ms       <- lk$manuscript_site
  pdir     <- lk$site_dir_used_in_pipeline
  flag     <- lk$flag_status
  lc       <- lk$dominant_lc

  log_subsection(sprintf("Manuscript site %d (tracker %d) — %s — flag=%s — LC=%s",
                         ms, tnum, pdir, flag, lc))

  # ── 3a: Read enriched CSV with forest filter ─────────────────────────────
  dt_f <- read_enriched_forest(tnum)

  # ── 3b: Determine ALS CRS ────────────────────────────────────────────────
  inv_row <- als_inv[manuscript_site == ms & is_pipeline_dir == "yes"]
  if (nrow(inv_row) == 0L) inv_row <- als_inv[manuscript_site == ms][1]
  als_epsg <- as.integer(inv_row$crs_epsg[1])
  if (is.na(als_epsg) || als_epsg == 0L) {
    log_progress(sprintf("  WARNING: no valid EPSG for tracker %d; skipping", tnum))
    next
  }
  log_progress(sprintf("  ALS CRS: EPSG:%d", als_epsg))

  # ── 3c: Load GEDI footprint geometry, project to ALS UTM ────────────────
  gxy <- load_gedi_xy(tnum, als_epsg, dt_f$shot_key)
  log_progress(sprintf("  Matched GEDI geometry: %s footprints",
                       format(nrow(gxy), big.mark = ",")))

  fp <- merge(dt_f, gxy, by = "shot_key", all = FALSE)
  log_progress(sprintf("  Joined footprints: %s rows", format(nrow(fp), big.mark = ",")))

  if (nrow(fp) < 50) {
    log_progress("  Too few joined footprints (<50); skipping site")
    next
  }

  # ── 3d: Assign tiles ─────────────────────────────────────────────────────
  tinfo <- assign_tile(fp$x, fp$y, TILE_SIZE_M)
  fp[, `:=`(tile_col = tinfo$tile_col,
            tile_row = tinfo$tile_row,
            tile_id  = tinfo$tile_id)]
  origin_x <- tinfo$origin_x
  origin_y <- tinfo$origin_y

  log_progress(sprintf("  Tile origin (UTM): (%.0f, %.0f); tile_size = %g m",
                       origin_x, origin_y, TILE_SIZE_M))

  # ── 3e: Aggregate per tile ───────────────────────────────────────────────
  # Mode of dominant LC: most-common code among the tile's footprints
  mode_chr <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0L) return(NA_character_)
    tbl <- sort(table(x), decreasing = TRUE)
    names(tbl)[1]
  }

  tile_dt <- fp[, .(
    n_footprints     = .N,
    offset_mean      = mean(offset, na.rm = TRUE),
    offset_sd        = sd(offset, na.rm = TRUE),
    chm_err_mean     = mean(chm_err, na.rm = TRUE),
    als_p90_mean     = mean(als_chm_p90, na.rm = TRUE),
    rh_98_mean       = mean(rh_98, na.rm = TRUE),
    cover_mean       = mean(cover, na.rm = TRUE),
    slope_mean       = mean(slope_mean, na.rm = TRUE),
    dominant_lc      = mode_chr(lc2022_l1_code),
    tile_col         = first(tile_col),
    tile_row         = first(tile_row)
  ), by = tile_id]

  # Tile centroids (UTM)
  tile_dt[, tile_centroid_x := origin_x + (tile_col + 0.5) * TILE_SIZE_M]
  tile_dt[, tile_centroid_y := origin_y + (tile_row + 0.5) * TILE_SIZE_M]

  # Apply minimum-footprint filter
  n_tiles_total <- nrow(tile_dt)
  tile_dt[, valid_tile := n_footprints >= MIN_FOOTPRINTS_PER_TILE]
  n_tiles_valid <- sum(tile_dt$valid_tile)
  log_progress(sprintf("  Tiles total: %d; valid (>=%d footprints): %d",
                       n_tiles_total, MIN_FOOTPRINTS_PER_TILE, n_tiles_valid))

  if (n_tiles_valid < 4) {
    log_progress("  Too few valid tiles (<4); skipping spatial tests")
    next
  }

  # ── 3f: Edge distance per tile (relative to convex hull of all footprints)
  fp_sf <- st_as_sf(fp[, .(x, y)], coords = c("x", "y"), crs = als_epsg)
  tile_dt[, edge_distance_m := tile_edge_distance(tile_dt, fp_sf)]

  # Annotate
  tile_dt[, `:=`(
    manuscript_site = ms,
    tracker_site    = tnum,
    site_short      = pdir,
    flag_status     = flag,
    dominant_lc_site = lc,
    tile_size_m     = TILE_SIZE_M,
    als_epsg        = als_epsg
  )]

  # ── 3g: Spatial-structure tests on valid tiles only ──────────────────────
  vt <- tile_dt[valid_tile == TRUE]

  # (i) Within-site SD of tile-mean offset
  ws_sd     <- sd(vt$offset_mean, na.rm = TRUE)
  ws_mean   <- mean(vt$offset_mean, na.rm = TRUE)
  ws_iqr    <- diff(quantile(vt$offset_mean, c(.25, .75), na.rm = TRUE))
  ws_q05    <- quantile(vt$offset_mean, .05, na.rm = TRUE)
  ws_q95    <- quantile(vt$offset_mean, .95, na.rm = TRUE)

  # (ii) Moran's I
  m_off <- moran_on_tiles(vt, "offset_mean", TILE_SIZE_M, n_perm = N_PERM)
  m_err <- moran_on_tiles(vt, "chm_err_mean", TILE_SIZE_M, n_perm = N_PERM)

  # (iii) Tile-level OLS: offset ~ cover + slope + edge_distance + lc
  # Build design data; drop tiles with NA in any predictor
  reg_dt <- vt[
    is.finite(offset_mean) & is.finite(cover_mean) &
    is.finite(slope_mean) & is.finite(edge_distance_m)
  ]
  reg_summary <- NULL
  reg_n       <- 0L
  reg_r2adj   <- NA_real_

  if (nrow(reg_dt) >= 8) {
    # Convert dominant_lc to factor with levels present at this site
    if (uniqueN(reg_dt$dominant_lc) > 1) {
      reg_dt[, dominant_lc := factor(dominant_lc)]
      fit <- tryCatch(
        lm(offset_mean ~ cover_mean + slope_mean + edge_distance_m + dominant_lc,
           data = reg_dt),
        error = function(e) NULL
      )
    } else {
      fit <- tryCatch(
        lm(offset_mean ~ cover_mean + slope_mean + edge_distance_m, data = reg_dt),
        error = function(e) NULL
      )
    }
    if (!is.null(fit)) {
      smry <- summary(fit)
      cf   <- as.data.table(smry$coefficients, keep.rownames = "term")
      setnames(cf, c("term", "estimate", "std_error", "t_value", "p_value"))
      cf[, `:=`(manuscript_site = ms, tracker_site = tnum,
                site_short = pdir, flag_status = flag,
                n_obs = nrow(reg_dt),
                r2 = smry$r.squared, r2_adj = smry$adj.r.squared)]
      reg_summary <- cf
      reg_n       <- nrow(reg_dt)
      reg_r2adj   <- smry$adj.r.squared
    }
  }

  # ── 3h: Per-site summary row ────────────────────────────────────────────
  site_row <- data.table(
    manuscript_site  = ms,
    tracker_site     = tnum,
    site_short       = pdir,
    flag_status      = flag,
    dominant_lc      = lc,
    tile_size_m      = TILE_SIZE_M,
    n_footprints     = nrow(fp),
    n_tiles_total    = n_tiles_total,
    n_tiles_valid    = n_tiles_valid,
    offset_mean_pooled  = mean(fp$offset, na.rm = TRUE),
    offset_median_pooled = median(fp$offset, na.rm = TRUE),
    offset_sd_pooled    = sd(fp$offset, na.rm = TRUE),
    tile_offset_mean    = ws_mean,
    tile_offset_sd      = ws_sd,
    tile_offset_iqr     = ws_iqr,
    tile_offset_q05     = ws_q05,
    tile_offset_q95     = ws_q95,
    moran_I_offset      = m_off$I,
    moran_p_offset      = m_off$p_perm,
    moran_n_offset      = m_off$n_used,
    moran_I_chm_err     = m_err$I,
    moran_p_chm_err     = m_err$p_perm,
    reg_n               = reg_n,
    reg_r2_adj          = reg_r2adj
  )

  log_progress(sprintf(
    "  -> tile_offset_mean = %+.2f m, sd = %.2f m, n_valid = %d",
    ws_mean, ws_sd, n_tiles_valid))
  log_progress(sprintf(
    "     Moran I (offset) = %+.3f, p_perm = %.3f (n=%d)",
    m_off$I, m_off$p_perm, m_off$n_used))
  if (!is.na(reg_r2adj)) {
    log_progress(sprintf("     OLS adj-R^2 = %.3f (n=%d)", reg_r2adj, reg_n))
  }

  # ── 3i: Stash for combined outputs ──────────────────────────────────────
  all_tile_rows[[as.character(tnum)]] <- tile_dt
  all_site_rows[[as.character(tnum)]] <- site_row
  if (!is.null(reg_summary)) {
    all_reg_rows[[as.character(tnum)]] <- reg_summary
  }
  all_plot_data[[as.character(tnum)]] <- list(
    tile_dt = vt, fp = fp, ms = ms, pdir = pdir,
    flag = flag, lc = lc, als_epsg = als_epsg,
    tile_size_m = TILE_SIZE_M
  )

  # cleanup
  rm(dt_f, gxy, fp, fp_sf, tile_dt, vt, reg_dt)
  gc(verbose = FALSE)
}

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4: Write CSV outputs
# ─────────────────────────────────────────────────────────────────────────────

log_section("Writing CSV outputs")

if (length(all_tile_rows) == 0L) stop("No site processed successfully.")

tile_summary_dt <- rbindlist(all_tile_rows, use.names = TRUE, fill = TRUE)
site_summary_dt <- rbindlist(all_site_rows, use.names = TRUE, fill = TRUE)
reg_summary_dt  <- if (length(all_reg_rows) > 0L) {
  rbindlist(all_reg_rows, use.names = TRUE, fill = TRUE)
} else {
  data.table()
}

fwrite(tile_summary_dt, file.path(MANUSCRIPT_TB, "groundwork_task6_phase1_tile_summary.csv"))
fwrite(site_summary_dt, file.path(MANUSCRIPT_TB, "groundwork_task6_phase1_site_summary.csv"))
fwrite(reg_summary_dt,  file.path(MANUSCRIPT_TB, "groundwork_task6_phase1_regression.csv"))

log_progress(sprintf("  tile rows: %d", nrow(tile_summary_dt)))
log_progress(sprintf("  site rows: %d", nrow(site_summary_dt)))
log_progress(sprintf("  regression rows: %d", nrow(reg_summary_dt)))

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5: Per-site detail figures
# ─────────────────────────────────────────────────────────────────────────────

log_section("Generating per-site detail figures")

# Common offset color scale: clip to robust 5/95 across all sites
all_tile_offsets <- tile_summary_dt[valid_tile == TRUE]$offset_mean
off_lim <- quantile(all_tile_offsets, c(0.02, 0.98), na.rm = TRUE)
off_lim <- c(-1, 1) * max(abs(off_lim))   # symmetric around zero
log_progress(sprintf("  Common offset color scale: %+.1f m to %+.1f m", off_lim[1], off_lim[2]))

err_lim <- quantile(tile_summary_dt[valid_tile == TRUE]$chm_err_mean,
                    c(0.02, 0.98), na.rm = TRUE)
err_lim <- c(-1, 1) * max(abs(err_lim))

make_site_detail <- function(pd) {
  vt <- pd$tile_dt
  ms <- pd$ms; pdir <- pd$pdir; flag <- pd$flag

  # Convert tile centroids back to a relative grid for nicer axes (km from SW)
  vt[, x_km := (tile_centroid_x - min(tile_centroid_x)) / 1000]
  vt[, y_km := (tile_centroid_y - min(tile_centroid_y)) / 1000]

  title_main <- sprintf("Site %d (%s)%s",
                        ms, pdir,
                        if (flag == "FLAGGED") " — FLAGGED" else "")

  # Panel A: tile offset map
  pA <- ggplot(vt, aes(x_km, y_km, fill = offset_mean)) +
    geom_tile(width = pd$tile_size_m / 1000, height = pd$tile_size_m / 1000) +
    scale_fill_gradient2(low = "#b2182b", mid = "#f7f7f7", high = "#2166ac",
                         midpoint = 0, limits = off_lim, oob = scales::squish,
                         name = "ALS p90 −\nGEDI rh_98 (m)") +
    coord_equal() +
    labs(x = "x (km from SW corner)", y = "y (km)", title = "Tile offset (negative = ALS underestimates)") +
    theme_cowplot(11) +
    theme(legend.position = "right")

  # Panel B: tile CHM error map
  pB <- ggplot(vt, aes(x_km, y_km, fill = chm_err_mean)) +
    geom_tile(width = pd$tile_size_m / 1000, height = pd$tile_size_m / 1000) +
    scale_fill_gradient2(low = "#762a83", mid = "#f7f7f7", high = "#1b7837",
                         midpoint = 0, limits = err_lim, oob = scales::squish,
                         name = "P3D − ALS\nCHM error (m)") +
    coord_equal() +
    labs(x = "x (km)", y = "y (km)", title = "Tile CHM error (P3D − ALS)") +
    theme_cowplot(11) +
    theme(legend.position = "right")

  # Panel C: offset vs cover scatter
  pC <- ggplot(vt, aes(cover_mean, offset_mean)) +
    geom_point(aes(size = n_footprints), alpha = 0.5, color = "#2166ac") +
    geom_smooth(method = "lm", se = TRUE, color = "grey30", linetype = "dashed") +
    geom_hline(yintercept = 0, color = "grey40", linetype = "dotted") +
    scale_size_continuous(range = c(1, 4), guide = "none") +
    labs(x = "Tile mean canopy cover", y = "Tile mean offset (m)",
         title = "Offset vs canopy cover") +
    theme_cowplot(11)

  # Panel D: offset vs edge distance scatter
  pD <- ggplot(vt, aes(edge_distance_m, offset_mean)) +
    geom_point(aes(size = n_footprints), alpha = 0.5, color = "#762a83") +
    geom_smooth(method = "lm", se = TRUE, color = "grey30", linetype = "dashed") +
    geom_hline(yintercept = 0, color = "grey40", linetype = "dotted") +
    scale_size_continuous(range = c(1, 4), guide = "none") +
    labs(x = "Distance to convex hull boundary (m)", y = "Tile mean offset (m)",
         title = "Offset vs edge distance") +
    theme_cowplot(11)

  combined <- plot_grid(pA, pB, pC, pD, ncol = 2, labels = c("A", "B", "C", "D"))
  title_grob <- ggdraw() +
    draw_label(title_main, fontface = "bold", x = 0, hjust = 0, size = 13) +
    theme(plot.margin = margin(t = 4, l = 8))
  plot_grid(title_grob, combined, ncol = 1, rel_heights = c(0.05, 1))
}

for (k in names(all_plot_data)) {
  pd <- all_plot_data[[k]]
  fn <- sprintf("task6_phase1_site%02d_%s_detail.pdf", pd$ms,
                substr(pd$pdir, 1, 16))
  fp_out <- file.path(PLOT_DIR, fn)
  fig <- make_site_detail(pd)
  ggsave(fp_out, fig, width = 12, height = 10, bg = "white")
  log_progress(sprintf("  Wrote %s", fn))
}

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 6: Combined comparison figures
# ─────────────────────────────────────────────────────────────────────────────

log_section("Generating combined comparison figures")

# (a) Side-by-side offset maps with common color scale
make_offset_map <- function(pd) {
  vt <- pd$tile_dt
  vt[, x_km := (tile_centroid_x - min(tile_centroid_x)) / 1000]
  vt[, y_km := (tile_centroid_y - min(tile_centroid_y)) / 1000]
  ttl <- sprintf("Site %d — %s%s", pd$ms, pd$pdir,
                 if (pd$flag == "FLAGGED") " (FLAGGED)" else " (control)")
  ggplot(vt, aes(x_km, y_km, fill = offset_mean)) +
    geom_tile(width = pd$tile_size_m / 1000, height = pd$tile_size_m / 1000) +
    scale_fill_gradient2(low = "#b2182b", mid = "#f7f7f7", high = "#2166ac",
                         midpoint = 0, limits = off_lim, oob = scales::squish,
                         name = "Offset (m)") +
    coord_equal() +
    labs(x = "x (km)", y = "y (km)", title = ttl) +
    theme_cowplot(10) +
    theme(legend.position = "right",
          plot.title = element_text(size = 11))
}

ordered_keys <- names(all_plot_data)[order(sapply(all_plot_data, function(p) p$ms))]
maps <- lapply(ordered_keys, function(k) make_offset_map(all_plot_data[[k]]))

combined_maps <- plot_grid(plotlist = maps, ncol = 2, align = "hv")
ggsave(file.path(PLOT_DIR, "task6_phase1_offset_maps_combined.pdf"),
       combined_maps,
       width = 14, height = max(8, 4 * ceiling(length(maps) / 2)),
       bg = "white")
log_progress("  Wrote task6_phase1_offset_maps_combined.pdf")

# (b) Tile-offset distribution per site (violin/density)
plot_dist_dt <- tile_summary_dt[valid_tile == TRUE,
  .(manuscript_site, site_short, flag_status, offset_mean)]
plot_dist_dt[, site_label := sprintf("Site %d (%s)%s",
                                      manuscript_site, site_short,
                                      ifelse(flag_status == "FLAGGED", " *", ""))]
plot_dist_dt[, site_label := factor(site_label,
  levels = unique(site_label[order(manuscript_site)]))]

p_dist <- ggplot(plot_dist_dt, aes(x = site_label, y = offset_mean,
                                   fill = flag_status)) +
  geom_violin(scale = "width", alpha = 0.6) +
  geom_jitter(width = 0.1, alpha = 0.3, size = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("FLAGGED" = "#d6604d", "ok" = "#4393c3"),
                    name = "Flag status") +
  labs(x = NULL, y = "Tile mean offset (ALS p90 − GEDI rh_98, m)",
       title = "Within-site distribution of tile-mean offsets",
       subtitle = sprintf("Flagged sites marked with *; tile size = %g m; min %d footprints/tile",
                          TILE_SIZE_M, MIN_FOOTPRINTS_PER_TILE)) +
  theme_cowplot(11) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 25, hjust = 1))

ggsave(file.path(PLOT_DIR, "task6_phase1_offset_distributions.pdf"),
       p_dist, width = 10, height = 6, bg = "white")
log_progress("  Wrote task6_phase1_offset_distributions.pdf")

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 7: Final summary
# ─────────────────────────────────────────────────────────────────────────────

log_section("Step 01 — Complete")

cat("\nSite-level summary:\n")
print(site_summary_dt[, .(
  ms = manuscript_site, site_short, flag_status,
  n_tiles_valid,
  off_mean_pool   = round(offset_mean_pooled,    2),
  tile_off_mean   = round(tile_offset_mean,      2),
  tile_off_sd     = round(tile_offset_sd,        2),
  moran_I         = round(moran_I_offset,        3),
  moran_p         = round(moran_p_offset,        3),
  reg_R2_adj      = round(reg_r2_adj,            3)
)])

if (nrow(reg_summary_dt) > 0L) {
  cat("\nKey OLS slopes (excluding intercept and LC factor levels):\n")
  reg_print <- reg_summary_dt[term %in% c("cover_mean", "slope_mean", "edge_distance_m"),
    .(ms = manuscript_site, site_short,
      term, est = round(estimate, 4), se = round(std_error, 4),
      p = round(p_value, 4))]
  print(reg_print)
}

log_progress("")
log_progress("Output files:")
log_progress(sprintf("  %s/groundwork_task6_phase1_tile_summary.csv", MANUSCRIPT_TB))
log_progress(sprintf("  %s/groundwork_task6_phase1_site_summary.csv", MANUSCRIPT_TB))
log_progress(sprintf("  %s/groundwork_task6_phase1_regression.csv",   MANUSCRIPT_TB))
for (k in ordered_keys) {
  pd <- all_plot_data[[k]]
  log_progress(sprintf("  %s/task6_phase1_site%02d_%s_detail.pdf",
                       PLOT_DIR, pd$ms, substr(pd$pdir, 1, 16)))
}
log_progress(sprintf("  %s/task6_phase1_offset_maps_combined.pdf",   PLOT_DIR))
log_progress(sprintf("  %s/task6_phase1_offset_distributions.pdf",   PLOT_DIR))

log_progress("")
log_progress("step 04 trigger heuristic (16-site atlas worth running if any are TRUE):")
log_progress("  (a) within-site tile_offset_sd >= 0.5 * between-site SD of medians")
log_progress("  (b) Moran I p_perm < 0.05 with |I| >= 0.2 at any flagged site")
log_progress("  (c) any tile-level OLS slope (cover_mean, slope_mean,")
log_progress("      edge_distance_m) clearing |t| >= 2 at any flagged site")
log_progress("")
log_progress("Outputs:")
log_progress("  - The three CSVs above")
log_progress("  - The PDF detail figures")
log_progress("  - This script's stdout/stderr log")
