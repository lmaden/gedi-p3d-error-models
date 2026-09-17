#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(sf); library(dplyr); library(tidyr); library(readr)
  library(future); library(future.apply); library(parallelly)
})
sf::sf_use_s2(FALSE)

# ---------- user paths ----------
SITE_FOOT_DIR   <- "/gpfs/data1/vclgp/lmaden/chpt1/data/footprints/sites"  # site_XX_footprints.parquet/gpkg
GEDI_ROOT       <- "/gpfs/data1/vclgp/lmaden/chpt1/gedi"                   # <id>/GEDI_site<id>_hq_ALL.gpkg
OUT_LOOKUPS_DIR <- "/gpfs/data1/vclgp/lmaden/chpt1/data/footprints/lookups"
dir.create(OUT_LOOKUPS_DIR, recursive = TRUE, showWarnings = FALSE)

# Options
USE_BUFFER      <- FALSE         # FALSE: PIP; TRUE: buffer 12.5 m and use intersects
JOIN_PREDICATE  <- "within"      # "within" (default) or "intersects"
BUFFER_M        <- 12.5
CRS_METERS      <- 5070          # For buffering only
CHUNK_SIZE      <- 50000         # GEDI shots per chunk
MAX_WORKERS     <- 30L
DEFAULT_WORKERS <- 12L

# ---- helpers ----
determine_workers <- function(max_target = MAX_WORKERS) {
  o <- Sys.getenv("IMG_META_WORKERS", "")
  if (nzchar(o)) {
    w <- suppressWarnings(as.integer(o))
    if (!is.na(w) && w > 0) return(max(1L, min(w, max_target)))
  }
  n <- NA_integer_
  if (requireNamespace("parallelly", quietly = TRUE)) n <- suppressWarnings(parallelly::availableCores())
  if (is.na(n)) {
    n <- suppressWarnings(try(parallel::detectCores(), silent = TRUE))
    n <- if (inherits(n, "try-error") || is.na(n)) 1L else n
  }
  max(1L, min(as.integer(n), max_target))
}

SHOT_ID_CANDIDATES <- c("shot_number","shot_num","shot_id","Shot_Number","shotID")
choose_shot_id <- function(gedi) {
  hit <- SHOT_ID_CANDIDATES[SHOT_ID_CANDIDATES %in% names(gedi)][1]
  if (is.na(hit)) stop("Cannot find a GEDI shot id column. Present: ", paste(names(gedi), collapse=", "))
  hit
}

agg_pairs_to_features <- function(pairs, shot_col = "shot_key") {
  if (!nrow(pairs)) return(tibble(!!shot_col := character(0)))
  
  s <- pairs |>
    group_by(.data[[shot_col]]) |>
    summarise(
      meta_tot_ct         = dplyr::n(),
      meta_stereo_ct      = sum(replace_na(as.integer(stereo_ind), 0L), na.rm = TRUE),
      meta_fwd_ct         = sum(replace_na(as.integer(is_forward), 0L), na.rm = TRUE),
      meta_rev_ct         = sum(replace_na(as.integer(is_reverse), 0L), na.rm = TRUE),
      meta_leaf_on_ct     = sum(replace_na(as.logical(leaf_on),  FALSE), na.rm = TRUE),
      meta_leaf_off_ct    = sum(replace_na(as.logical(leaf_off), FALSE), na.rm = TRUE),
      
      offnad_sum          = sum(off_nadir_avg, na.rm = TRUE),
      offnad_n            = sum(!is.na(off_nadir_avg)),
      sunel_sum           = sum(sun_elevation_avg, na.rm = TRUE),
      sunel_n             = sum(!is.na(sun_elevation_avg)),
      absgeo_sum          = sum(absolute_geolocation_accuracy, na.rm = TRUE),
      absgeo_n            = sum(!is.na(absolute_geolocation_accuracy)),
      relgeo_sum          = sum(relative_geolocation_accuracy, na.rm = TRUE),
      relgeo_n            = sum(!is.na(relative_geolocation_accuracy)),
      az_sin_sum          = sum(az_sin, na.rm = TRUE),
      az_cos_sum          = sum(az_cos, na.rm = TRUE),
      az_n                = sum(!is.na(target_azimuth_avg)),
      .groups = "drop"
    )
  
  veh_ct <- pairs |>
    filter(!is.na(veh) & veh != "") |>
    count(.data[[shot_col]], veh, name = "ct") |>
    pivot_wider(names_from = veh, values_from = ct,
                names_prefix = "meta_veh_ct_", values_fill = 0)
  
  res <- left_join(s, veh_ct, by = shot_col)
  
  veh_cols <- grep("^meta_veh_ct_", names(res), value = TRUE)
  for (v in veh_cols) {
    res[[sub("^meta_veh_ct_", "meta_veh_ratio_", v)]] <-
      ifelse(res$meta_tot_ct > 0, res[[v]] / res$meta_tot_ct, NA_real_)
  }
  
  res |>
    mutate(
      meta_stereo_any       = as.integer(meta_stereo_ct > 0),
      meta_stereo_ratio     = if_else(meta_tot_ct > 0, meta_stereo_ct / meta_tot_ct, NA_real_),
      meta_fwd_ratio        = if_else((meta_fwd_ct + meta_rev_ct) > 0, meta_fwd_ct / (meta_fwd_ct + meta_rev_ct), NA_real_),
      meta_rev_ratio        = if_else((meta_fwd_ct + meta_rev_ct) > 0, meta_rev_ct / (meta_fwd_ct + meta_rev_ct), NA_real_),
      meta_leaf_on_ratio    = if_else((meta_leaf_on_ct + meta_leaf_off_ct) > 0, meta_leaf_on_ct / (meta_leaf_on_ct + meta_leaf_off_ct), NA_real_),
      
      meta_off_nadir_avg    = if_else(offnad_n > 0, offnad_sum / offnad_n, NA_real_),
      meta_sun_elev_avg     = if_else(sunel_n  > 0,  sunel_sum  / sunel_n,  NA_real_),
      meta_abs_geoacc_avg   = if_else(absgeo_n > 0,  absgeo_sum / absgeo_n, NA_real_),
      meta_rel_geoacc_avg   = if_else(relgeo_n > 0,  relgeo_sum / relgeo_n, NA_real_),
      
      meta_target_azimuth_avg = if_else(az_n > 0, (atan2(az_sin_sum, az_cos_sum) * 180/pi) %% 360, NA_real_),
      meta_az_concentration   = if_else(az_n > 0, sqrt(az_sin_sum^2 + az_cos_sum^2) / az_n, NA_real_)
    ) |>
    select(
      all_of(shot_col),
      starts_with("meta_veh_ct_"), starts_with("meta_veh_ratio_"),
      meta_tot_ct, meta_stereo_ct, meta_stereo_any, meta_stereo_ratio,
      meta_fwd_ct, meta_rev_ct, meta_fwd_ratio, meta_rev_ratio,
      meta_leaf_on_ct, meta_leaf_off_ct, meta_leaf_on_ratio,
      meta_off_nadir_avg, meta_sun_elev_avg,
      meta_abs_geoacc_avg, meta_rel_geoacc_avg, meta_target_azimuth_avg, meta_az_concentration
    )
}

process_one_site <- function(site_id) {
  # inputs
  dsn_pq   <- file.path(SITE_FOOT_DIR, sprintf("site_%02d_footprints.parquet", site_id))
  dsn_gpkg <- sub("\\.parquet$", ".gpkg", dsn_pq)
  fpd      <- if (file.exists(dsn_pq)) dsn_pq else if (file.exists(dsn_gpkg)) dsn_gpkg else NA
  if (is.na(fpd)) { message(sprintf("[site %d] no per-site footprints found", site_id)); return(invisible(NULL)) }
  
  gedi_gpkg <- file.path(GEDI_ROOT, as.character(site_id), sprintf("GEDI_site%d_hq_ALL.gpkg", site_id))
  if (!file.exists(gedi_gpkg)) { message(sprintf("[site %d] missing GEDI file", site_id)); return(invisible(NULL)) }
  
  message(sprintf("[site %d] reading GEDI + footprints …", site_id))
  fp   <- suppressWarnings(st_read(fpd, quiet = TRUE))
  gedi <- suppressWarnings(st_read(gedi_gpkg, quiet = TRUE,
                                   options = "PROMOTE_TO_MULTI=NO",
                                   int64_as_string = TRUE))
  if (!nrow(gedi)) { message(sprintf("[site %d] GEDI empty", site_id)); return(invisible(NULL)) }
  
  nm_shot <- choose_shot_id(gedi)
  # keep only the shot id column; geometry is preserved automatically
  gedi <- gedi[, nm_shot, drop = FALSE]
  names(gedi)[1] <- "shot_key"
  # ensure character (avoids 64-bit int issues)
  gedi$shot_key <- as.character(gedi$shot_key)
  
  # CRS alignment
  gedi_x <- st_transform(gedi, st_crs(fp))
  
  # choose join
  join_fun <- switch(tolower(JOIN_PREDICATE),
                     within = st_within,
                     intersects = st_intersects,
                     contains = st_contains,
                     st_within)
  
  # optional 12.5 m buffer
  if (USE_BUFFER) {
    gedi_m <- st_transform(gedi_x, CRS_METERS)
    gbuf   <- st_buffer(gedi_m, BUFFER_M)
    gedi_x <- st_transform(gbuf, st_crs(fp))
    # st_join uses intersects with polygons anyway when buffered
    pairs  <- suppressWarnings(st_join(st_set_geometry(gedi, st_geometry(gedi_x)),
                                       fp, join = st_intersects, left = FALSE))
  } else {
    pairs  <- suppressWarnings(st_join(gedi_x, fp, join = join_fun, left = FALSE))
  }
  
  if (!nrow(pairs)) {
    message(sprintf("[site %d] no overlaps → writing empty lookup", site_id))
    out <- file.path(OUT_LOOKUPS_DIR, sprintf("meta_lookup_site%02d.parquet", site_id))
    if (requireNamespace("arrow", quietly = TRUE) && ("write_parquet" %in% getNamespaceExports("arrow")))
      arrow::write_parquet(tibble(shot_key=character()), out)
    else
      readr::write_csv(tibble(shot_key=character()), sub("\\.parquet$", ".csv.gz", out))
    return(invisible(NULL))
  }
  
  # Chunk by GEDI to limit memory (pairs are computed per-chunk)
  # Build chunk indices
  n <- nrow(gedi_x)
  starts <- seq(1, n, by = CHUNK_SIZE)
  out_chunks <- vector("list", length(starts))
  for (i in seq_along(starts)) {
    a <- starts[i]; b <- min(n, a + CHUNK_SIZE - 1L)
    gi <- gedi_x[a:b, , drop = FALSE]
    pr <- suppressWarnings(st_join(gi, fp, join = join_fun, left = FALSE))
    if (nrow(pr)) {
      st_geometry(pr) <- NULL
      out_chunks[[i]] <- agg_pairs_to_features(pr, "shot_key")
    } else {
      out_chunks[[i]] <- tibble(shot_key = character(0))
    }
    if (i %% 5 == 0) message(sprintf("[site %d] chunks: %d/%d", site_id, i, length(starts)))
  }
  
  res <- bind_rows(out_chunks)
  
  out <- file.path(OUT_LOOKUPS_DIR, sprintf("meta_lookup_site%02d.parquet", site_id))
  if (requireNamespace("arrow", quietly = TRUE) && ("write_parquet" %in% getNamespaceExports("arrow"))) {
    arrow::write_parquet(res, out)
  } else {
    readr::write_csv(res, sub("\\.parquet$", ".csv.gz", out))
  }
  message(sprintf("[site %d] wrote lookup (%d rows): %s", site_id, nrow(res), out))
  invisible(NULL)
}

# ---------- main ----------
files <- list.files(SITE_FOOT_DIR, pattern = "^site_([0-9]{2})_footprints\\.(parquet|gpkg)$",
                    full.names = TRUE, ignore.case = TRUE)
if (!length(files)) stop("No per-site footprint files found in: ", SITE_FOOT_DIR)

site_ids <- as.integer(sub("^site_0*([0-9]+)_.*$", "\\1", basename(files)))
site_ids <- sort(unique(site_ids))

workers <- min(determine_workers(), length(site_ids))
if (workers >= 2L) {
  message(sprintf("Using %d worker(s) (cap %d).", workers, MAX_WORKERS)); plan(multisession, workers = workers)
} else {
  message("Running sequentially."); plan(sequential)
}
on.exit(plan(sequential), add = TRUE)

future.apply::future_lapply(site_ids, process_one_site,
                            future.chunk.size = 1L, future.seed = TRUE)

message("✔ Step 2 done: per‑site metadata lookups created in: ", OUT_LOOKUPS_DIR)
