#!/usr/bin/env Rscript
# =============================================================================
# site5_04_cross_site_atlas.R
#
# Step 04.
# 19-site atlas of within-site ALS offset spatial structure, plus
# cross-site analysis of per-site cover-slope (b_cover) against the
# three site-RE frames (Ch1, 16-site, Decomposition 15-site err_alt).
#
# Hypothesis (from step 01):
#   Per-site variation in the cover-slope b_cover (the OLS coefficient
#   on tile-mean canopy cover in the within-site offset regression) is
#   the structural mechanism for the residual ~12% V_site (sigma_site
#   ~ 1.34 m) remainder after 16-site's flagging removal and CNN-DTM
#   compensation are accounted for.
#
#   Prediction: more steeply negative b_cover ~ more positive site RE.
#   Sign of correlation should be NEGATIVE.
#
# Design changes vs step 01:
#   1. SITES expanded to all 19 manuscript sites.
#   2. edge_distance_m DROPPED from the tile-level OLS (closed in step 01).
#   3. Cross-site b_cover summary block added; correlations against
#      Ch1, 16-site, decomposition site REs across multiple frames.
#   4. CRS fallback: if the ALS inventory has no valid EPSG (Site 5),
#      compute UTM zone from GEDI footprint mean longitude.
#
# Inputs:
#   - data/enriched_by_site/site_NN_enriched.csv.gz
#   - gedi/N/GEDI_siteN_hq_ALL.gpkg
#   - manuscript_tables/site_id_lookup.csv
#   - manuscript_tables/als_repeat_acquisition_inventory.csv
#   - manuscript_tables/site_random_effects_chm_dtm.csv  (Ch1 REs)
#   - manuscript_tables/track_b_re_intercepts.csv         (16-site REs)
#   - manuscript_tables/coupling_random_effect_shifts.csv (decomposition step)
#
# Outputs (manuscript_tables/):
#   - groundwork_task6_phase2_tile_summary.csv     per (site x tile)
#   - groundwork_task6_phase2_site_summary.csv     per site (incl. b_cover)
#   - groundwork_task6_phase2_regression.csv       per-site OLS coefs
#   - groundwork_task6_phase2_cross_site.csv       per site b_cover + RE
#   - groundwork_task6_phase2_correlations.csv     cross-site r tests
#
# Outputs (plots/groundwork/):
#   - task6_phase2_offset_maps_atlas.pdf           19-site offset map atlas
#   - task6_phase2_bcover_bar.pdf                  per-site b_cover bar chart
#   - task6_phase2_bcover_vs_re_panels.pdf         3-panel scatter vs REs
#
# Run: source("site5_04_cross_site_atlas.R") from Pane 1.
#
# Tunables via env:
#   TASK6_TILE_M=250            (default; matches step 01)
#   TASK6_MIN_FP=5              (min footprints per tile)
# =============================================================================

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 0: Setup
# ─────────────────────────────────────────────────────────────────────────────

source("analysis_config.R")
source("analysis_utils.R")

log_section("Step 04: 19-site atlas + cross-site b_cover analysis")

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
N_PERM                  <- 999
FOREST_THRESH_M         <- 2
VALID_FRAC_THRESH       <- 0.5
ERROR_OUTLIER_M         <- 100
MIN_TILES_FOR_OLS       <- 8
MIN_TILES_FOR_MORAN     <- 4

ALL_SITES <- 1:19  # Manuscript site IDs

log_progress(sprintf("  Tile size:           %g m", TILE_SIZE_M))
log_progress(sprintf("  Min footprints/tile: %d",   MIN_FOOTPRINTS_PER_TILE))
log_progress(sprintf("  Sites:               manuscript 1-19 (all)"))
log_progress(sprintf("  Permutations for Moran's I: %d", N_PERM))

# Paths
ALS_BASE      <- Sys.getenv("ALS_BASE", "/gpfs/data1/vclgp/data/gedi/imported/usa")
MANUSCRIPT_TB <- file.path(PROJECT_ROOT, "manuscript_tables")
GEDI_BASE     <- file.path(PROJECT_ROOT, "gedi")
ENRICHED_DIR  <- file.path(PROJECT_ROOT, "data", "enriched_by_site")
PLOT_DIR      <- file.path(PROJECT_ROOT, "plots", "groundwork")
dir.create(PLOT_DIR,      recursive = TRUE, showWarnings = FALSE)
dir.create(MANUSCRIPT_TB, recursive = TRUE, showWarnings = FALSE)

# RE file paths
RE_CH1_PATH      <- file.path(MANUSCRIPT_TB, "site_random_effects_chm_dtm.csv")
RE_TRACKB_PATH   <- file.path(MANUSCRIPT_TB, "track_b_re_intercepts.csv")
RE_PHASE3_PATH   <- file.path(MANUSCRIPT_TB, "coupling_random_effect_shifts.csv")

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
# SECTION 2: RE CSV discovery — defensive column-name resolution
# ─────────────────────────────────────────────────────────────────────────────

log_subsection("Discovering RE CSV column names")

#' Read an RE CSV with defensive column-name discovery.
#' Returns a normalized data.table with columns (tracker_site, re_estimate)
#' or NULL if the file/columns can't be resolved.
read_re_csv <- function(path, label,
                        key_candidates,
                        value_candidates,
                        prefer_match = NULL) {
  if (!file.exists(path)) {
    log_progress(sprintf("  %s: NOT FOUND at %s", label, path))
    return(NULL)
  }
  dt <- tryCatch(fread(path), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0L) {
    log_progress(sprintf("  %s: read failed or empty", label))
    return(NULL)
  }
  cols <- names(dt)
  log_progress(sprintf("  %s: columns = [%s]", label, paste(cols, collapse = ", ")))

  # Resolve key
  key_col <- NA_character_
  for (cand in key_candidates) {
    if (cand %in% cols) { key_col <- cand; break }
  }
  if (is.na(key_col)) {
    log_progress(sprintf("    -> WARNING: no recognized key column. Tried: %s",
                         paste(key_candidates, collapse = " | ")))
    return(NULL)
  }

  # Resolve value: prefer columns matching `prefer_match` regex first
  val_col <- NA_character_
  if (!is.null(prefer_match)) {
    matches <- grep(prefer_match, cols, value = TRUE, ignore.case = TRUE)
    matches <- intersect(matches, value_candidates)
    if (length(matches) > 0) val_col <- matches[1]
  }
  if (is.na(val_col)) {
    for (cand in value_candidates) {
      if (cand %in% cols) { val_col <- cand; break }
    }
  }
  if (is.na(val_col)) {
    log_progress(sprintf("    -> WARNING: no recognized value column. Tried: %s",
                         paste(value_candidates, collapse = " | ")))
    return(NULL)
  }

  out <- data.table(
    tracker_site = suppressWarnings(as.integer(dt[[key_col]])),
    re_estimate  = suppressWarnings(as.numeric(dt[[val_col]]))
  )
  out <- out[!is.na(tracker_site) & !is.na(re_estimate)]
  log_progress(sprintf("    -> resolved: key=%s, value=%s; %d sites loaded",
                       key_col, val_col, nrow(out)))
  out
}

ch1_re <- read_re_csv(
  RE_CH1_PATH, "Ch1 site REs",
  key_candidates = c("tracker_site_id", "tracker_site", "site",
                     "tracker_site_num", "ms_site_id", "manuscript_site"),
  value_candidates = c("chm_intercept_estimate", "chm_intercept", "chm_re",
                       "chm_site_re", "re_chm", "Estimate_CHM",
                       "Estimate", "estimate", "mean", "intercept"),
  prefer_match = "chm"
)

trackb_re <- read_re_csv(
  RE_TRACKB_PATH, "16-site REs",
  key_candidates = c("tracker_site", "site", "tracker_site_id",
                     "ms_site_id", "manuscript_site"),
  value_candidates = c("intercept_estimate", "intercept", "Estimate",
                       "estimate", "mean", "re_intercept", "chm_re",
                       "chm_intercept")
)

phase3_re <- read_re_csv(
  RE_PHASE3_PATH, "Decomposition 15-site err_alt REs",
  key_candidates = c("site", "tracker_site", "tracker_site_id",
                     "ms_site_id", "manuscript_site"),
  value_candidates = c("re_alt_15", "re_alt", "re_15site", "re_phase3_15",
                       "re_alt_15site", "re_phase3", "phase3_re",
                       "alt_re", "Estimate", "estimate", "mean",
                       "re_shift", "delta_re"),
  prefer_match = "alt|phase3|15"
)

if (is.null(ch1_re) && is.null(trackb_re) && is.null(phase3_re)) {
  log_progress("  WARNING: NO RE CSVs resolved. Cross-site correlations will all be NA.")
  log_progress("           step 04 still produces tile/site/cross-site CSVs;")
  log_progress("           inspect column names manually if you need correlations.")
}

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: Helper functions (duplicated from step 01 + new fallbacks)
# ─────────────────────────────────────────────────────────────────────────────

log_subsection("Defining helper functions")

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

  if (!("lc2022_l1_code" %in% names(dt)) && "lc2022_mode_l1_code" %in% names(dt)) {
    setnames(dt, "lc2022_mode_l1_code", "lc2022_l1_code")
  }
  if (!("lc2022_l1_code" %in% names(dt))) dt[, lc2022_l1_code := NA_character_]

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

  dt_f[, chm_err := fifelse(is.finite(error_mean), error_mean,
                             p3d_chm_mean - als_chm_mean)]
  dt_f <- dt_f[!is.finite(chm_err) | abs(chm_err) <= ERROR_OUTLIER_M]

  dt_f[, offset := als_chm_p90 - rh_98]
  dt_f[, shot_key := as.character(shot_number)]
  log_progress(sprintf("    %s after forest filter", format(nrow(dt_f), big.mark = ",")))
  dt_f
}

# Reproject GEDI footprints to a target CRS (with fallback if als_epsg is NA)
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
  gedi_sf <- gedi_sf[gedi_sf$shot_key %in% target_keys, ]
  if (nrow(gedi_sf) == 0L) stop("No matching GEDI footprints for tracker ", tracker_num)

  # CRS fallback: compute UTM zone from mean longitude in WGS84
  if (is.na(als_epsg) || als_epsg == 0L) {
    gedi_wgs <- st_transform(gedi_sf, 4326)
    coords_wgs <- sf::st_coordinates(gedi_wgs)
    mean_lon <- mean(coords_wgs[, 1], na.rm = TRUE)
    mean_lat <- mean(coords_wgs[, 2], na.rm = TRUE)
    utm_zone <- floor((mean_lon + 180) / 6) + 1L
    als_epsg <- 32600L + utm_zone  # northern hemisphere
    log_progress(sprintf(
      "    CRS fallback: site centroid lon=%.2f lat=%.2f -> EPSG:%d (UTM %dN)",
      mean_lon, mean_lat, als_epsg, utm_zone))
  }

  pts_proj <- st_transform(gedi_sf, crs = als_epsg)
  xy <- sf::st_coordinates(pts_proj)
  list(
    xy = data.table(shot_key = pts_proj$shot_key, x = xy[, 1], y = xy[, 2]),
    epsg = als_epsg
  )
}

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

moran_on_tiles <- function(tile_dt, value_col, tile_m, n_perm = 999) {
  v <- tile_dt[[value_col]]
  ok <- is.finite(v)
  if (sum(ok) < MIN_TILES_FOR_MORAN) return(list(I = NA_real_, p_perm = NA_real_, n_used = sum(ok)))

  sub <- tile_dt[ok]
  v_sub <- v[ok]

  coords <- as.matrix(sub[, .(tile_centroid_x, tile_centroid_y)])
  nb <- tryCatch(
    dnearneigh(coords, d1 = 0, d2 = 1.05 * tile_m),
    error = function(e) NULL
  )
  if (is.null(nb)) return(list(I = NA_real_, p_perm = NA_real_, n_used = nrow(sub)))

  card <- card(nb)
  keep <- card > 0
  if (sum(keep) < MIN_TILES_FOR_MORAN) return(list(I = NA_real_, p_perm = NA_real_, n_used = sum(keep)))

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
# SECTION 4: Per-site processing loop (all 19 sites)
# ─────────────────────────────────────────────────────────────────────────────

log_section("Per-site tile aggregation across all 19 sites")

all_tile_rows <- list()
all_site_rows <- list()
all_reg_rows  <- list()
all_plot_data <- list()

for (ms in ALL_SITES) {

  lk <- site_lookup[manuscript_site == ms]
  if (nrow(lk) == 0L) {
    log_progress(sprintf("  Manuscript site %d not in lookup; skipping", ms))
    next
  }
  tnum     <- lk$tracker_site
  pdir     <- lk$site_dir_used_in_pipeline
  flag     <- lk$flag_status
  lc       <- lk$dominant_lc

  log_subsection(sprintf("Manuscript site %d (tracker %d) — %s — flag=%s — LC=%s",
                         ms, tnum, pdir, flag, lc))

  # ── 4a: Read enriched CSV ─────────────────────────────────────────────────
  dt_f <- tryCatch(read_enriched_forest(tnum), error = function(e) {
    log_progress(sprintf("  ERROR reading enriched CSV: %s", conditionMessage(e)))
    NULL
  })
  if (is.null(dt_f) || nrow(dt_f) < 50L) {
    log_progress(sprintf("  Too few forest footprints (%d); skipping site",
                         if (is.null(dt_f)) 0 else nrow(dt_f)))
    # Still write a stub site row so the cross-site analysis sees it
    all_site_rows[[as.character(ms)]] <- data.table(
      manuscript_site = ms, tracker_site = tnum,
      site_short = pdir, flag_status = flag, dominant_lc = lc,
      tile_size_m = TILE_SIZE_M,
      n_footprints = if (is.null(dt_f)) 0L else nrow(dt_f),
      n_tiles_total = NA_integer_, n_tiles_valid = NA_integer_,
      pooled_offset_mean = NA_real_, pooled_offset_median = NA_real_,
      pooled_offset_sd = NA_real_,
      tile_offset_mean = NA_real_, tile_offset_sd = NA_real_,
      moran_I_offset = NA_real_, moran_p_offset = NA_real_,
      b_cover = NA_real_, se_b_cover = NA_real_,
      t_b_cover = NA_real_, p_b_cover = NA_real_,
      b_slope = NA_real_, se_b_slope = NA_real_,
      t_b_slope = NA_real_, p_b_slope = NA_real_,
      reg_n = 0L, reg_r2_adj = NA_real_, skip_reason = "too_few_footprints"
    )
    next
  }

  # ── 4b: Determine ALS CRS (with fallback) ────────────────────────────────
  inv_row <- als_inv[manuscript_site == ms & is_pipeline_dir == "yes"]
  if (nrow(inv_row) == 0L) inv_row <- als_inv[manuscript_site == ms][1]
  als_epsg <- if (nrow(inv_row) > 0L) {
    suppressWarnings(as.integer(inv_row$crs_epsg[1]))
  } else NA_integer_

  if (is.na(als_epsg) || als_epsg == 0L) {
    log_progress("  No valid CRS in inventory; will use UTM-zone fallback")
  } else {
    log_progress(sprintf("  ALS CRS: EPSG:%d", als_epsg))
  }

  # ── 4c: Load + project GEDI footprints (with fallback CRS if needed) ────
  gedi_loaded <- tryCatch(
    load_gedi_xy(tnum, als_epsg, dt_f$shot_key),
    error = function(e) {
      log_progress(sprintf("  ERROR loading GEDI: %s", conditionMessage(e)))
      NULL
    }
  )
  if (is.null(gedi_loaded)) {
    all_site_rows[[as.character(ms)]] <- data.table(
      manuscript_site = ms, tracker_site = tnum,
      site_short = pdir, flag_status = flag, dominant_lc = lc,
      tile_size_m = TILE_SIZE_M,
      n_footprints = nrow(dt_f),
      n_tiles_total = NA_integer_, n_tiles_valid = NA_integer_,
      pooled_offset_mean = NA_real_, pooled_offset_median = NA_real_,
      pooled_offset_sd = NA_real_,
      tile_offset_mean = NA_real_, tile_offset_sd = NA_real_,
      moran_I_offset = NA_real_, moran_p_offset = NA_real_,
      b_cover = NA_real_, se_b_cover = NA_real_,
      t_b_cover = NA_real_, p_b_cover = NA_real_,
      b_slope = NA_real_, se_b_slope = NA_real_,
      t_b_slope = NA_real_, p_b_slope = NA_real_,
      reg_n = 0L, reg_r2_adj = NA_real_, skip_reason = "gedi_load_failed"
    )
    next
  }
  used_epsg <- gedi_loaded$epsg
  fp <- merge(dt_f, gedi_loaded$xy, by = "shot_key", all = FALSE)
  log_progress(sprintf("  Joined footprints: %s rows", format(nrow(fp), big.mark = ",")))

  # ── 4d: Tile assignment ──────────────────────────────────────────────────
  tinfo <- assign_tile(fp$x, fp$y, TILE_SIZE_M)
  fp[, `:=`(tile_col = tinfo$tile_col,
            tile_row = tinfo$tile_row,
            tile_id  = tinfo$tile_id)]
  origin_x <- tinfo$origin_x
  origin_y <- tinfo$origin_y

  # ── 4e: Aggregate per tile ──────────────────────────────────────────────
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

  tile_dt[, tile_centroid_x := origin_x + (tile_col + 0.5) * TILE_SIZE_M]
  tile_dt[, tile_centroid_y := origin_y + (tile_row + 0.5) * TILE_SIZE_M]

  n_tiles_total <- nrow(tile_dt)
  tile_dt[, valid_tile := n_footprints >= MIN_FOOTPRINTS_PER_TILE]
  n_tiles_valid <- sum(tile_dt$valid_tile)
  log_progress(sprintf("  Tiles total: %d; valid (>=%d footprints): %d",
                       n_tiles_total, MIN_FOOTPRINTS_PER_TILE, n_tiles_valid))

  # Annotate
  tile_dt[, `:=`(
    manuscript_site = ms,
    tracker_site    = tnum,
    site_short      = pdir,
    flag_status     = flag,
    dominant_lc_site = lc,
    tile_size_m     = TILE_SIZE_M,
    als_epsg        = used_epsg
  )]

  # ── 4f: Spatial-structure tests (Moran's I) ─────────────────────────────
  vt <- tile_dt[valid_tile == TRUE]
  if (n_tiles_valid >= MIN_TILES_FOR_MORAN) {
    m_off <- moran_on_tiles(vt, "offset_mean", TILE_SIZE_M, n_perm = N_PERM)
  } else {
    m_off <- list(I = NA_real_, p_perm = NA_real_, n_used = n_tiles_valid)
  }

  # ── 4g: Tile-level OLS — DROPS edge_distance_m vs step 01 ───────────────
  reg_dt <- vt[is.finite(offset_mean) & is.finite(cover_mean) & is.finite(slope_mean)]

  reg_summary <- NULL
  reg_n       <- 0L
  reg_r2adj   <- NA_real_
  b_cover     <- NA_real_; se_b_cover <- NA_real_
  t_b_cover   <- NA_real_; p_b_cover  <- NA_real_
  b_slope     <- NA_real_; se_b_slope <- NA_real_
  t_b_slope   <- NA_real_; p_b_slope  <- NA_real_

  if (nrow(reg_dt) >= MIN_TILES_FOR_OLS) {
    if (uniqueN(reg_dt$dominant_lc) > 1) {
      reg_dt[, dominant_lc := factor(dominant_lc)]
      fit <- tryCatch(
        lm(offset_mean ~ cover_mean + slope_mean + dominant_lc, data = reg_dt),
        error = function(e) NULL
      )
    } else {
      fit <- tryCatch(
        lm(offset_mean ~ cover_mean + slope_mean, data = reg_dt),
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

      cov_row <- cf[term == "cover_mean"]
      if (nrow(cov_row) > 0L) {
        b_cover    <- cov_row$estimate
        se_b_cover <- cov_row$std_error
        t_b_cover  <- cov_row$t_value
        p_b_cover  <- cov_row$p_value
      }
      slp_row <- cf[term == "slope_mean"]
      if (nrow(slp_row) > 0L) {
        b_slope    <- slp_row$estimate
        se_b_slope <- slp_row$std_error
        t_b_slope  <- slp_row$t_value
        p_b_slope  <- slp_row$p_value
      }
    }
  }

  # ── 4h: Per-site summary ────────────────────────────────────────────────
  site_row <- data.table(
    manuscript_site      = ms,
    tracker_site         = tnum,
    site_short           = pdir,
    flag_status          = flag,
    dominant_lc          = lc,
    tile_size_m          = TILE_SIZE_M,
    n_footprints         = nrow(fp),
    n_tiles_total        = n_tiles_total,
    n_tiles_valid        = n_tiles_valid,
    pooled_offset_mean   = mean(fp$offset, na.rm = TRUE),
    pooled_offset_median = median(fp$offset, na.rm = TRUE),
    pooled_offset_sd     = sd(fp$offset, na.rm = TRUE),
    tile_offset_mean     = if (nrow(vt) > 0) mean(vt$offset_mean, na.rm = TRUE) else NA_real_,
    tile_offset_sd       = if (nrow(vt) > 0) sd(vt$offset_mean, na.rm = TRUE) else NA_real_,
    moran_I_offset       = m_off$I,
    moran_p_offset       = m_off$p_perm,
    b_cover              = b_cover,
    se_b_cover           = se_b_cover,
    t_b_cover            = t_b_cover,
    p_b_cover            = p_b_cover,
    b_slope              = b_slope,
    se_b_slope           = se_b_slope,
    t_b_slope            = t_b_slope,
    p_b_slope            = p_b_slope,
    reg_n                = reg_n,
    reg_r2_adj           = reg_r2adj,
    skip_reason          = NA_character_
  )

  log_progress(sprintf(
    "  -> tile_off_mean = %+.2f m, sd = %.2f m, n_valid = %d, b_cover = %s, R2adj = %s",
    site_row$tile_offset_mean, site_row$tile_offset_sd, n_tiles_valid,
    if (is.finite(b_cover))  sprintf("%+.2f", b_cover)  else "NA",
    if (is.finite(reg_r2adj)) sprintf("%.3f", reg_r2adj) else "NA"))

  # Stash
  all_tile_rows[[as.character(ms)]] <- tile_dt
  all_site_rows[[as.character(ms)]] <- site_row
  if (!is.null(reg_summary)) {
    all_reg_rows[[as.character(ms)]] <- reg_summary
  }
  all_plot_data[[as.character(ms)]] <- list(
    tile_dt = vt, ms = ms, pdir = pdir, flag = flag, lc = lc,
    tile_size_m = TILE_SIZE_M
  )

  rm(dt_f, gedi_loaded, fp, tile_dt, vt, reg_dt)
  gc(verbose = FALSE)
}

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5: Write per-site CSV outputs
# ─────────────────────────────────────────────────────────────────────────────

log_section("Writing per-site CSV outputs")

if (length(all_tile_rows) == 0L) stop("No site processed successfully.")

tile_summary_dt <- rbindlist(all_tile_rows, use.names = TRUE, fill = TRUE)
site_summary_dt <- rbindlist(all_site_rows, use.names = TRUE, fill = TRUE)
reg_summary_dt  <- if (length(all_reg_rows) > 0L) {
  rbindlist(all_reg_rows, use.names = TRUE, fill = TRUE)
} else data.table()

fwrite(tile_summary_dt, file.path(MANUSCRIPT_TB, "groundwork_task6_phase2_tile_summary.csv"))
fwrite(site_summary_dt, file.path(MANUSCRIPT_TB, "groundwork_task6_phase2_site_summary.csv"))
fwrite(reg_summary_dt,  file.path(MANUSCRIPT_TB, "groundwork_task6_phase2_regression.csv"))

log_progress(sprintf("  tile rows: %d", nrow(tile_summary_dt)))
log_progress(sprintf("  site rows: %d", nrow(site_summary_dt)))
log_progress(sprintf("  regression rows: %d", nrow(reg_summary_dt)))

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 6: Cross-site b_cover analysis
# ─────────────────────────────────────────────────────────────────────────────

log_section("Cross-site b_cover analysis")

# Build the cross-site table: per-site b_cover joined with each RE source.
cross_site <- site_summary_dt[, .(
  manuscript_site, tracker_site, site_short, flag_status, dominant_lc,
  n_tiles_valid, tile_offset_mean, tile_offset_sd,
  b_cover, se_b_cover, t_b_cover, p_b_cover, reg_r2_adj
)]

# Merge in REs (each may be NULL if file/columns failed to resolve)
if (!is.null(ch1_re)) {
  cross_site <- merge(cross_site,
                       ch1_re[, .(tracker_site, ch1_chm_re = re_estimate)],
                       by = "tracker_site", all.x = TRUE)
} else {
  cross_site[, ch1_chm_re := NA_real_]
}
if (!is.null(trackb_re)) {
  cross_site <- merge(cross_site,
                       trackb_re[, .(tracker_site, track_b_re = re_estimate)],
                       by = "tracker_site", all.x = TRUE)
} else {
  cross_site[, track_b_re := NA_real_]
}
if (!is.null(phase3_re)) {
  cross_site <- merge(cross_site,
                       phase3_re[, .(tracker_site, phase3_re = re_estimate)],
                       by = "tracker_site", all.x = TRUE)
} else {
  cross_site[, phase3_re := NA_real_]
}

setorder(cross_site, manuscript_site)

# Frames
make_frame <- function(name, mask) {
  list(name = name, mask = mask)
}
frames <- list(
  make_frame("19_sites",
             rep(TRUE, nrow(cross_site))),
  make_frame("16_non_flagged",
             cross_site$flag_status != "FLAGGED"),
  make_frame("15_err_alt",
             cross_site$flag_status != "FLAGGED" & cross_site$flag_status != "DTM_excluded"),
  make_frame("forest_only_no_flagged",
             cross_site$flag_status != "FLAGGED" &
               cross_site$dominant_lc %in% c("BDF", "ENF", "DNF", "EBF", "MFT", "IWL"))
)

# Predictor / response combinations
predictors  <- c("b_cover", "tile_offset_sd", "tile_offset_mean")
responses   <- c("ch1_chm_re", "track_b_re", "phase3_re")

corr_rows <- list()
for (fr in frames) {
  sub <- cross_site[fr$mask]
  for (p in predictors) {
    for (r in responses) {
      x <- sub[[p]]; y <- sub[[r]]
      ok <- is.finite(x) & is.finite(y)
      if (sum(ok) < 4L) {
        corr_rows[[length(corr_rows) + 1]] <- data.table(
          frame = fr$name, predictor = p, response = r,
          n = sum(ok),
          r_pearson = NA_real_, p_pearson = NA_real_,
          r_spearman = NA_real_, p_spearman = NA_real_
        )
        next
      }
      ct_p <- tryCatch(cor.test(x[ok], y[ok], method = "pearson"),
                        error = function(e) NULL)
      ct_s <- tryCatch(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE),
                        error = function(e) NULL)
      corr_rows[[length(corr_rows) + 1]] <- data.table(
        frame = fr$name, predictor = p, response = r,
        n = sum(ok),
        r_pearson = if (!is.null(ct_p)) unname(ct_p$estimate) else NA_real_,
        p_pearson = if (!is.null(ct_p)) ct_p$p.value         else NA_real_,
        r_spearman = if (!is.null(ct_s)) unname(ct_s$estimate) else NA_real_,
        p_spearman = if (!is.null(ct_s)) ct_s$p.value         else NA_real_
      )
    }
  }
}
corr_dt <- rbindlist(corr_rows)

fwrite(cross_site, file.path(MANUSCRIPT_TB, "groundwork_task6_phase2_cross_site.csv"))
fwrite(corr_dt,    file.path(MANUSCRIPT_TB, "groundwork_task6_phase2_correlations.csv"))
log_progress(sprintf("  Wrote cross_site (%d rows) and correlations (%d rows)",
                     nrow(cross_site), nrow(corr_dt)))

# Print headline table
cat("\nCross-site b_cover and tile_offset_sd correlations (Pearson):\n")
print(corr_dt[predictor %in% c("b_cover", "tile_offset_sd"),
              .(frame, predictor, response, n,
                r = round(r_pearson, 3), p = round(p_pearson, 4))])

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 7: Atlas figure — 19-site offset maps
# ─────────────────────────────────────────────────────────────────────────────

log_section("Generating atlas + cross-site figures")

all_tile_offsets <- tile_summary_dt[valid_tile == TRUE]$offset_mean
off_lim <- quantile(all_tile_offsets, c(0.02, 0.98), na.rm = TRUE)
off_lim <- c(-1, 1) * max(abs(off_lim))
log_progress(sprintf("  Common offset color scale: %+.1f to %+.1f m",
                     off_lim[1], off_lim[2]))

make_atlas_map <- function(pd) {
  vt <- pd$tile_dt
  if (nrow(vt) == 0L) {
    return(ggplot() + theme_void() +
             labs(title = sprintf("Site %d — %s (no valid tiles)", pd$ms, pd$pdir)))
  }
  vt[, x_km := (tile_centroid_x - min(tile_centroid_x)) / 1000]
  vt[, y_km := (tile_centroid_y - min(tile_centroid_y)) / 1000]
  ttl_flag <- if (pd$flag == "FLAGGED") "*" else ""
  ttl <- sprintf("Site %d — %s%s", pd$ms, pd$pdir, ttl_flag)
  ggplot(vt, aes(x_km, y_km, fill = offset_mean)) +
    geom_tile(width = pd$tile_size_m / 1000, height = pd$tile_size_m / 1000) +
    scale_fill_gradient2(low = "#b2182b", mid = "#f7f7f7", high = "#2166ac",
                         midpoint = 0, limits = off_lim, oob = scales::squish,
                         name = "Offset (m)", guide = "none") +
    coord_equal() +
    labs(x = NULL, y = NULL, title = ttl) +
    theme_cowplot(8) +
    theme(plot.title = element_text(size = 9, face = "plain"),
          axis.text = element_text(size = 6),
          axis.ticks = element_line(linewidth = 0.2))
}

ordered_keys <- names(all_plot_data)[order(as.integer(names(all_plot_data)))]
atlas_panels <- lapply(ordered_keys, function(k) make_atlas_map(all_plot_data[[k]]))

# Add a shared legend strip
legend_panel <- ggplot(data.frame(x = 1, y = 1), aes(x, y, fill = x)) +
  geom_tile() +
  scale_fill_gradient2(low = "#b2182b", mid = "#f7f7f7", high = "#2166ac",
                       midpoint = 0, limits = off_lim, oob = scales::squish,
                       name = "ALS p90 -\nGEDI rh_98 (m)") +
  theme_void() + theme(legend.position = "right",
                       legend.title = element_text(size = 9))
shared_legend <- get_legend(legend_panel)

atlas_grid <- plot_grid(plotlist = atlas_panels, ncol = 4, align = "hv")
atlas_combined <- plot_grid(atlas_grid, shared_legend, ncol = 2, rel_widths = c(1, 0.10))

ggsave(file.path(PLOT_DIR, "task6_phase2_offset_maps_atlas.pdf"),
       atlas_combined,
       width = 16, height = 4 * ceiling(length(atlas_panels) / 4),
       bg = "white")
log_progress("  Wrote task6_phase2_offset_maps_atlas.pdf")

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 8: Per-site b_cover bar chart
# ─────────────────────────────────────────────────────────────────────────────

bcover_dt <- cross_site[is.finite(b_cover)]
bcover_dt[, site_label := sprintf("%d %s", manuscript_site, site_short)]
bcover_dt[, site_label := factor(site_label,
  levels = site_label[order(b_cover)])]

p_bar <- ggplot(bcover_dt, aes(x = site_label, y = b_cover, fill = flag_status)) +
  geom_col(width = 0.75) +
  geom_errorbar(aes(ymin = b_cover - 1.96 * se_b_cover,
                    ymax = b_cover + 1.96 * se_b_cover),
                width = 0.25, linewidth = 0.4) +
  geom_hline(yintercept = 0, color = "grey30") +
  scale_fill_manual(values = c("FLAGGED" = "#d6604d", "ok" = "#4393c3",
                                "DTM_excluded" = "#bf812d"),
                    name = "Flag status") +
  coord_flip() +
  labs(x = NULL, y = "b_cover (m of offset per unit canopy cover)",
       title = "Per-site cover slope (b_cover) sorted ascending",
       subtitle = sprintf(
         "Tile-level OLS: offset ~ cover_mean + slope_mean + LC. Tile size = %g m.",
         TILE_SIZE_M)) +
  theme_cowplot(11) +
  theme(legend.position = "top")

ggsave(file.path(PLOT_DIR, "task6_phase2_bcover_bar.pdf"),
       p_bar, width = 9, height = 7, bg = "white")
log_progress("  Wrote task6_phase2_bcover_bar.pdf")

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 9: 3-panel scatter: b_cover vs each RE
# ─────────────────────────────────────────────────────────────────────────────

mk_scatter <- function(dt, re_col, re_label, frame_label) {
  d <- dt[is.finite(b_cover) & is.finite(get(re_col))]
  if (nrow(d) < 3L) {
    return(ggplot() + theme_void() +
             labs(title = sprintf("%s\n(insufficient data)", re_label)))
  }
  ct <- cor.test(d$b_cover, d[[re_col]])
  r_pe <- unname(ct$estimate); p_pe <- ct$p.value
  ggplot(d, aes(x = b_cover, y = get(re_col), color = flag_status)) +
    geom_point(size = 2.5) +
    geom_text(aes(label = manuscript_site), nudge_x = 0.5, size = 3,
              color = "grey20") +
    geom_smooth(method = "lm", se = TRUE, color = "grey30",
                linetype = "dashed", inherit.aes = FALSE,
                aes(x = b_cover, y = get(re_col))) +
    geom_hline(yintercept = 0, color = "grey50", linetype = "dotted") +
    geom_vline(xintercept = 0, color = "grey50", linetype = "dotted") +
    scale_color_manual(values = c("FLAGGED" = "#d6604d", "ok" = "#4393c3",
                                   "DTM_excluded" = "#bf812d"),
                       name = "Flag status") +
    labs(x = "b_cover (m / unit cover)",
         y = re_label,
         title = sprintf("%s\n%s: r = %+.2f, p = %.3f, n = %d",
                          re_label, frame_label, r_pe, p_pe, nrow(d))) +
    theme_cowplot(10) +
    theme(legend.position = "none",
          plot.title = element_text(size = 10))
}

# Frames: 16-site, 19-site (published CHM fit) and 15-site (decomposition)
panel_ch1   <- mk_scatter(cross_site, "ch1_chm_re",   "Ch1 CHM site RE",         "19 sites")
panel_trkb  <- mk_scatter(cross_site[flag_status != "FLAGGED"],
                           "track_b_re", "16-site RE",      "16 non-flagged")
panel_ph3   <- mk_scatter(cross_site[flag_status != "FLAGGED" & flag_status != "DTM_excluded"],
                           "phase3_re",  "Decomposition 15-site err_alt RE", "15 sites")

# Add a shared legend
legend_dummy <- ggplot(cross_site, aes(x = b_cover, y = b_cover, color = flag_status)) +
  geom_point() +
  scale_color_manual(values = c("FLAGGED" = "#d6604d", "ok" = "#4393c3",
                                 "DTM_excluded" = "#bf812d"),
                     name = "Flag status") +
  theme_cowplot(10) + theme(legend.position = "top")
leg <- get_legend(legend_dummy)

panels_row <- plot_grid(panel_ch1, panel_trkb, panel_ph3, ncol = 3, align = "hv")
panels_combined <- plot_grid(leg, panels_row, ncol = 1, rel_heights = c(0.07, 1))

ggsave(file.path(PLOT_DIR, "task6_phase2_bcover_vs_re_panels.pdf"),
       panels_combined, width = 15, height = 5.5, bg = "white")
log_progress("  Wrote task6_phase2_bcover_vs_re_panels.pdf")

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 10: Final summary
# ─────────────────────────────────────────────────────────────────────────────

log_section("Step 04 — Complete")

cat("\nPer-site b_cover summary (sorted by manuscript_site):\n")
print(site_summary_dt[, .(
  ms = manuscript_site, site = site_short, flag = flag_status, lc = dominant_lc,
  n_valid_tiles = n_tiles_valid,
  tile_off_mean = round(tile_offset_mean, 2),
  tile_off_sd   = round(tile_offset_sd,   2),
  b_cover       = round(b_cover,          2),
  t_bcov        = round(t_b_cover,        1),
  R2_adj        = round(reg_r2_adj,       3)
)])

cat("\nCross-site b_cover correlation against site REs (Pearson, all frames):\n")
print(corr_dt[predictor == "b_cover",
              .(frame, response, n, r = round(r_pearson, 3),
                p = round(p_pearson, 4))])

log_progress("")
log_progress("Output files:")
log_progress(sprintf("  %s/groundwork_task6_phase2_tile_summary.csv",  MANUSCRIPT_TB))
log_progress(sprintf("  %s/groundwork_task6_phase2_site_summary.csv",  MANUSCRIPT_TB))
log_progress(sprintf("  %s/groundwork_task6_phase2_regression.csv",    MANUSCRIPT_TB))
log_progress(sprintf("  %s/groundwork_task6_phase2_cross_site.csv",    MANUSCRIPT_TB))
log_progress(sprintf("  %s/groundwork_task6_phase2_correlations.csv",  MANUSCRIPT_TB))
log_progress(sprintf("  %s/task6_phase2_offset_maps_atlas.pdf",        PLOT_DIR))
log_progress(sprintf("  %s/task6_phase2_bcover_bar.pdf",               PLOT_DIR))
log_progress(sprintf("  %s/task6_phase2_bcover_vs_re_panels.pdf",      PLOT_DIR))

log_progress("")
log_progress("Hypothesis test:")
log_progress("  Predicted sign: r(b_cover, site_RE) NEGATIVE.")
log_progress("  (More steeply negative b_cover -> stronger ALS underestimation")
log_progress("   -> P3D appears to overestimate -> positive site RE.)")
log_progress("")
log_progress("Decision rule for the variance puzzle remainder:")
log_progress("  * r negative & significant in 16-site frame:")
log_progress("    -> b_cover variation explains residual ~12% V_site;")
log_progress("       variance puzzle closed.")
log_progress("  * r non-significant in 16-site frame:")
log_progress("    -> remainder is something else (DSM-side phenology,")
log_progress("       image vintage, sun-canopy geometry, footprint geoloc).")
log_progress("")
log_progress("Outputs:")
log_progress("  - The five CSVs above")
log_progress("  - The three PDFs above")
log_progress("  - This script's stdout/stderr log")
