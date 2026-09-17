#!/usr/bin/env Rscript
# Normalize imagery footprints PER SITE and write per-site cleaned outputs.
# Input pattern: footprints_dir/cats_<SITEID>_within_aoi.geojson
# Paired GEDI:   gedi_root/<SITEID>/GEDI_site<SITEID>_hq_ALL.gpkg
# Output:        out_dir/sites/site_<SITEID>_footprints.parquet (or .gpkg fallback)

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(lubridate)
  library(future)
  library(future.apply)
  library(parallelly)
})

sf::sf_use_s2(FALSE)

# -------- user paths --------
FOOTPRINTS_DIR <- "/gpfs/data1/vclgp/lmaden/chpt1/data/footprints"   # cats_<id>_within_aoi.geojson
GEDI_ROOT      <- "/gpfs/data1/vclgp/lmaden/chpt1/gedi"              # <id>/GEDI_site<id>_hq_ALL.gpkg
OUT_DIR        <- "/gpfs/data1/vclgp/lmaden/chpt1/data/footprints"   # outputs under OUT_DIR/sites/

# Options
CLIP_TO_GEDI_HULL      <- TRUE
MAX_TARGET_WORKERS     <- 30L
DEFAULT_WORKERS        <- 10L

# -------- helpers --------
determine_workers <- function(min_target = DEFAULT_WORKERS, max_target = MAX_TARGET_WORKERS) {
  o <- Sys.getenv("IMG_META_WORKERS", "")
  if (nzchar(o)) {
    w <- suppressWarnings(as.integer(o))
    if (!is.na(w) && w > 0) return(max(1L, min(w, max_target)))
  }
  n <- NA_integer_
  if (requireNamespace("parallelly", quietly = TRUE)) {
    n <- suppressWarnings(parallelly::availableCores(ensureLoadAverage = TRUE))
  }
  if (is.na(n)) {
    n <- suppressWarnings(try(parallel::detectCores(), silent = TRUE))
    n <- if (inherits(n, "try-error") || is.na(n)) 1L else n
  }
  max(1L, min(as.integer(n), max_target))
}

has_write_geoparquet <- function() {
  if (!requireNamespace("arrow", quietly = TRUE)) return(FALSE)
  "write_geoparquet" %in% getNamespaceExports("arrow")
}

write_geo <- function(g, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  if (tools::file_ext(path) == "parquet" && has_write_geoparquet()) {
    arrow::write_geoparquet(g, path)
    return(path)
  }
  alt <- sub("\\.parquet$", ".gpkg", path)
  message(sprintf("GeoParquet not available → writing GPKG: %s", alt))
  if (file.exists(alt)) unlink(alt)
  sf::st_write(g, alt, quiet = TRUE)
  alt
}

norm_vehicle <- function(v) {
  if (is.na(v) || !nzchar(v)) return("UNK")
  vv <- toupper(gsub("[ -]", "", v))
  vv <- gsub("WORLDVIEW", "WV", vv)
  vv <- gsub("GEOEYE", "GE", vv)
  vv <- gsub("QUICKBIRD", "QB", vv)
  for (k in c("WV01","WV02","WV03","WV04","GE01","QB02")) {
    if (grepl(k, vv, fixed = TRUE)) return(k)
  }
  vv
}

is_leaf_off <- function(ts) {
  if (is.na(ts)) return(FALSE)
  m <- month(ts); d <- mday(ts)
  (m > 10L) || (m == 10L && d >= 14L) || (m < 3L) || (m == 3L && d <= 8L)
}

detect_cols <- function(fp) {
  col_has <- function(nm) nm %in% names(fp)
  nm_vehicle <- if (col_has("vehicle_name")) "vehicle_name" else names(fp)[grepl("vehicle", names(fp), ignore.case=TRUE)][1]
  nm_time    <- if (col_has("collect_time_end")) "collect_time_end" else names(fp)[grepl("collect.*time", names(fp), ignore.case=TRUE)][1]
  nm_stereo  <- if (col_has("has_stereo_pair")) "has_stereo_pair" else names(fp)[grepl("stereo", names(fp), ignore.case=TRUE)][1]
  nm_scan    <- if (col_has("scan_direction")) "scan_direction" else names(fp)[grepl("scan.*direction", names(fp), ignore.case=TRUE)][1]
  list(
    vehicle = nm_vehicle %||% NA_character_,
    time    = nm_time    %||% NA_character_,
    stereo  = nm_stereo  %||% NA_character_,
    scan    = nm_scan    %||% NA_character_,
    offnad  = "off_nadir_avg",
    az      = "target_azimuth_avg",
    sunel   = "sun_elevation_avg",
    absacc  = "absolute_geolocation_accuracy",
    relacc  = "relative_geolocation_accuracy"
  )
}
`%||%` <- function(a,b) if (!is.null(a)) a else b
coerce_numeric <- function(x) suppressWarnings(as.numeric(x))

normalize_footprints <- function(fp) {
  cols <- detect_cols(fp)
  for (nm in unlist(cols)) if (!is.na(nm) && !(nm %in% names(fp))) fp[[nm]] <- NA
  if (!("geometry" %in% names(fp))) stop("Footprints layer has no geometry column")
  
  ts <- suppressWarnings(lubridate::ymd_hms(fp[[cols$time]], tz = "UTC"))
  if (all(is.na(ts))) {
    ts <- suppressWarnings(lubridate::parse_date_time(fp[[cols$time]],
                                                      orders = c("Ymd HMS", "Ymd HM", "Ymd", "YmdTz", "Ymd HMSz"), tz = "UTC"))
  }
  
  fp |>
    mutate(
      veh         = vapply(.data[[cols$vehicle]], norm_vehicle, character(1)),
      ts_utc      = ts,
      leaf_off    = vapply(ts, is_leaf_off, logical(1)),
      leaf_on     = !leaf_off,
      stereo_ind  = as.integer(tidyr::replace_na(as.integer(.data[[cols$stereo]]), 0L)),
      scan_dir    = tolower(trimws(as.character(.data[[cols$scan]]))),
      is_forward  = as.integer(scan_dir == "forward"),
      is_reverse  = as.integer(scan_dir == "reverse"),
      off_nadir_avg = coerce_numeric(.data[[cols$offnad]]),
      target_azimuth_avg = coerce_numeric(.data[[cols$az]]),
      sun_elevation_avg  = coerce_numeric(.data[[cols$sunel]]),
      absolute_geolocation_accuracy = coerce_numeric(.data[[cols$absacc]]),
      relative_geolocation_accuracy = coerce_numeric(.data[[cols$relacc]])
    ) |>
    mutate(
      az_sin = sin(pi * target_azimuth_avg / 180),
      az_cos = cos(pi * target_azimuth_avg / 180)
    ) |>
    dplyr::select(
      veh, ts_utc, leaf_off, leaf_on, stereo_ind, scan_dir,
      is_forward, is_reverse,
      off_nadir_avg, target_azimuth_avg, sun_elevation_avg,
      absolute_geolocation_accuracy, relative_geolocation_accuracy,
      az_sin, az_cos, geometry
    )
}

extract_site_id <- function(path) {
  b <- basename(path)
  m <- regexec("^cats_([0-9]+)_within_aoi\\.geojson$", b, ignore.case = TRUE)
  r <- regmatches(b, m)[[1]]
  if (length(r) < 2) return(NA_integer_)
  as.integer(r[2])
}

process_one_site_file <- function(geojson_path) {
  site_id <- extract_site_id(geojson_path)
  if (is.na(site_id)) {
    message("Skipping (cannot parse site id): ", geojson_path); return(invisible(NULL))
  }
  
  gedi_gpkg <- file.path(GEDI_ROOT, as.character(site_id), sprintf("GEDI_site%d_hq_ALL.gpkg", site_id))
  if (!file.exists(gedi_gpkg)) {
    message(sprintf("[site %d] Missing GEDI file → skip (%s)", site_id, gedi_gpkg))
    return(invisible(NULL))
  }
  
  message(sprintf("[site %d] Reading footprints: %s", site_id, geojson_path))
  fp_raw <- suppressWarnings(sf::st_read(geojson_path, quiet = TRUE))
  if (!nrow(fp_raw)) { message(sprintf("[site %d] Footprints empty → skip", site_id)); return(invisible(NULL)) }
  
  if (any(!sf::st_is_valid(fp_raw))) {
    message(sprintf("[site %d] Fixing invalid geometries …", site_id))
    if (requireNamespace("lwgeom", quietly = TRUE)) fp_raw <- lwgeom::st_make_valid(fp_raw) else fp_raw$geometry <- sf::st_buffer(fp_raw$geometry, 0)
  }
  
  fp_clean <- normalize_footprints(fp_raw)
  
  # Optional: clip to GEDI convex hull
  if (CLIP_TO_GEDI_HULL) {
    gedi <- suppressWarnings(sf::st_read(gedi_gpkg, quiet = TRUE, options = "PROMOTE_TO_MULTI=NO"))
    if (nrow(gedi)) {
      hull <- gedi |>
        st_transform(st_crs(fp_clean)) |>
        st_union() |>
        st_convex_hull()
      # ---- FIXED: build logical vector over ALL footprints ----
      idx <- lengths(st_intersects(fp_clean, hull)) > 0
      kept <- sum(idx); total <- nrow(fp_clean)
      fp_clean <- fp_clean[idx, , drop = FALSE]
      message(sprintf("[site %d] Clip kept %d / %d footprints", site_id, kept, total))
    }
  }
  
  out_sites_dir <- file.path(OUT_DIR, "sites")
  dir.create(out_sites_dir, showWarnings = FALSE, recursive = TRUE)
  out_pq <- file.path(out_sites_dir, sprintf("site_%02d_footprints.parquet", site_id))
  message(sprintf("[site %d] Writing %s (n=%d)", site_id, out_pq, nrow(fp_clean)))
  write_geo(fp_clean, out_pq)
  
  if (nrow(fp_clean)) {
    topv <- head(sort(table(fp_clean$veh), decreasing = TRUE), 5)
    message(sprintf("[site %d] Vehicles top: %s", site_id,
                    paste(paste(names(topv), as.integer(topv), sep=":"), collapse=", ")))
  }
  invisible(NULL)
}

# -------- main --------
files <- list.files(FOOTPRINTS_DIR,
                    pattern = "^cats_[0-9]+_within_aoi\\.geojson$",
                    full.names = TRUE, ignore.case = TRUE)
if (!length(files)) stop("No site GeoJSONs found in: ", FOOTPRINTS_DIR)

site_ids <- vapply(files, extract_site_id, integer(1))
ord <- order(site_ids, na.last = NA)
files <- files[ord]; site_ids <- site_ids[ord]

workers <- min(determine_workers(), length(files))
if (workers >= 2L) {
  message(sprintf("Using %d worker(s) via multisession (cap %d).", workers, MAX_TARGET_WORKERS))
  plan(multisession, workers = workers)
} else {
  message("Running sequentially.")
  plan(sequential)
}
on.exit(plan(sequential), add = TRUE)

future.apply::future_lapply(files, process_one_site_file,
                            future.chunk.size = 1L, future.seed = TRUE)

message("✔ Per‑site normalization complete. Outputs in: ", file.path(OUT_DIR, "sites"))

