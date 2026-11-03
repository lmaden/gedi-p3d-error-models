#!/usr/bin/env Rscript
# =====================================================================
# lidar_products_rebuild.R  --  lidar ground and canopy at every footprint
# Paper: the lidar references (Section 2.2.2, Text S1), behind every error in the tables and figures.
#
# Every reference is built from each site's own point cloud with one method
# at all 18 sites:
#   ground  TIN of the class-2 returns, noise removed (lidR ivf), 1 m grid
#   CHM     all returns normalized on that ground raster, heights outside
#           -2..120 m dropped, pit-free 1 m CHM
#           (Khosravipour et al. 2014; lidR pitfree), negatives set to 0
#   surface ground + CHM (built downstream)
# Both are averaged over the footprint disc: the 529 cells of 1 m whose centres
# lie within 13 cells of the footprint's cell, on a grid aligned to whole
# metres. Output per footprint: lg_mean, lg_valid, chm_mean, chm_p90, chm_valid.
#
# Checks written next to them, per site (never used as references):
#   ground minus 3DEP (dep_dtm_mean); rebuilt CHM against the delivered CHM
#   (als_chm_mean, als_chm_p90); rebuilt CHM p90 against GEDI rh_98.
#
# Footprints: the enriched footprint tables of the sampling stage
# (site_NN_enriched.csv.gz), rows with P3D CHM and slope valid
# on at least half the disc (the frame needs both). shot_number is text.
#
# Output: $CH1/manuscript_tables/lidar_products/lp_site_NN.csv.gz, lp_summary.csv,
#         lp_log.txt
# Run
#   module purge; module load rh9/R/4.5.0; module load rh9/gdal/3.11.0
#   export R_LIBS=<R_LIBRARY> CH1=<CLUSTER_ROOT>/chpt1
#   cd $CH1/scripts/reviewed
#   Rscript lidar_products_rebuild.R --limit 200 --sites 4 1            (smoke)
#   nohup Rscript lidar_products_rebuild.R --sites ... --cores 24 > lp_A.out 2>&1 &
#   Options: --sites (tracker ids), --limit N, --cores N, --overwrite
# =====================================================================

suppressPackageStartupMessages({
  library(data.table); library(sf); library(lidR); library(terra); library(parallel)
})
set_lidr_threads(1L); setDTthreads(1L)

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default = NULL) {
  i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[i + 1]
}
getints <- function(flag) {
  i <- match(flag, args); if (is.na(i)) return(NULL)
  j <- i + 1; out <- integer(0)
  while (j <= length(args) && grepl("^[0-9]+$", args[j])) { out <- c(out, as.integer(args[j])); j <- j + 1 }
  out
}

CH1       <- Sys.getenv("CH1", "<CLUSTER_ROOT>/chpt1")
IMPORTED  <- "<LIDAR_ARCHIVE>"   # one folder per flight, each with LAZ_ground/
LIMIT     <- as.integer(getarg("--limit", "0"))
CORES     <- as.integer(getarg("--cores", "24"))
OVERWRITE <- "--overwrite" %in% args
SEED      <- 20260923L
BUF       <- 30      # m around a cell's footprints so the TINs cover every disc
R_DISC    <- 13L     # the footprint disc at 1 m (helper_speed._metric_kernel of the sampling stage)

# tracker id -> point-cloud folder. Sites 6 and 20: the flight whose CHM matches
# the delivered one (flight comparison). Sites 1-3: the only flights.
SITE_CODE <- c(`1` = "usda_me", `2` = "nasa_howland", `3` = "neon_sawb",
               `4` = "neon_harv2019", `5` = "usda_sc", `6` = "neon_jerc2019",
               `7` = "neon_tall2021", `8` = "neon_dela2021", `9` = "neon_leno2021",
               `11` = "neon_konz2020", `12` = "neon_stei2022", `13` = "neon_unde2022",
               `14` = "neon_steicheq2022", `15` = "neon_wood2021", `17` = "neon_ster2022",
               `18` = "neon_cper2021", `19` = "ltbmu_20180710", `20` = "neon_sjer2019")
CRS_FALLBACK <- c(`5` = 26917L)

SITES <- getints("--sites"); if (is.null(SITES)) SITES <- as.integer(names(SITE_CODE))
OUTDIR <- file.path(CH1, "manuscript_tables", "lidar_products")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
LOG <- file(file.path(OUTDIR, "lp_log.txt"), open = "at")
msg <- function(...) {
  line <- paste0(format(Sys.time(), "%H:%M:%S "), sprintf(...))
  cat(line, "\n", sep = ""); writeLines(line, LOG); flush(LOG)
}
if (!requireNamespace("geometry", quietly = TRUE))
  msg("NOTE: R package 'geometry' is not installed; some lidR steps need it (without it, cells failed at site 1)")
msg("lidar_products_rebuild.R | sites %s | limit %d | cores %d%s",
    paste(SITES, collapse = " "), LIMIT, CORES, if (LIMIT > 0) " | SMOKE" else "")

disc <- CJ(dx = -R_DISC:R_DISC, dy = -R_DISC:R_DISC)[dx^2 + dy^2 <= R_DISC^2]
stopifnot(nrow(disc) == 529L)
PITFREE <- pitfree(thresholds = c(0, 2, 5, 10, 15), max_edge = c(0, 1.5))

footprints <- function(sid) {
  en <- fread(file.path(CH1, "data", "enriched_by_site", sprintf("site_%02d_enriched.csv.gz", sid)),
              select = c("shot_number", "dep_dtm_mean", "als_chm_mean", "als_chm_p90", "rh_98",
                         "p3d_chm_valid_frac", "slope_valid_frac"),
              colClasses = c(shot_number = "character"))
  en <- unique(en[p3d_chm_valid_frac >= 0.5 & slope_valid_frac >= 0.5], by = "shot_number")
  gp <- file.path(CH1, "gedi", sid, sprintf("GEDI_site%d_hq_ALL.gpkg", sid))
  fp <- suppressWarnings(st_read(gp, quiet = TRUE, int64_as_string = TRUE,
          query = sprintf("select shot_number, lon_lowestmode, lat_lowestmode from GEDI_site%d_hq_ALL", sid)))
  if (inherits(fp, "sf")) fp <- st_drop_geometry(fp)
  fp <- as.data.table(fp)
  if (!is.character(fp$shot_number)) stop("gpkg shot_number was not read as text")
  fp <- unique(fp[!is.na(lon_lowestmode) & !is.na(lat_lowestmode)], by = "shot_number")
  d <- merge(en, fp, by = "shot_number")
  if (LIMIT > 0L && nrow(d) > LIMIT) { set.seed(SEED); d <- d[sample.int(.N, LIMIT)] }
  d
}

one_cell <- function(fp, ctg, crs_str) {
  fx <- floor(fp$x); fy <- floor(fp$y)
  out <- data.table(shot_number = fp$shot_number, lg_mean = NA_real_, lg_valid = 0L, lg_sd = NA_real_,
                    chm_mean = NA_real_, chm_p90 = NA_real_, chm_valid = 0L,
                    n_ground = 0L, n_noise = 0L)
  a <- clip_rectangle(ctg, min(fx) - R_DISC - BUF, min(fy) - R_DISC - BUF,
                      max(fx) + R_DISC + 1 + BUF, max(fy) + R_DISC + 1 + BUF)
  if (is.null(a) || is.empty(a)) return(out)
  if (is.na(st_crs(a))) st_crs(a) <- st_crs(crs_str)
  a <- filter_poi(a, Classification != 7L & Classification != 18L)   # delivered noise classes
  # ground: class-2 returns, isolated points removed (ivf on the ground returns
  # only; ivf on all returns fails in filter_with_grid)
  g <- filter_poi(a, Classification == 2L)
  if (npoints(g) < 10L) return(out)
  n0 <- npoints(g)
  g <- classify_noise(g, ivf(res = 5, n = 6))
  g <- filter_poi(g, Classification != LASNOISE)
  out[, `:=`(n_ground = npoints(g), n_noise = n0 - npoints(g))]
  if (npoints(g) < 10L) return(out)
  # 1 m grid on whole metres, 10 m wider than the discs; the ground TIN uses
  # returns out to BUF so it is complete over the whole grid
  E <- R_DISC + 10L
  tmpl <- rast(xmin = min(fx) - E, xmax = max(fx) + E + 1, ymin = min(fy) - E, ymax = max(fy) + E + 1,
               resolution = 1, crs = crs_str)
  dtm <- rasterize_terrain(g, res = tmpl, algorithm = tin())
  # canopy: returns inside the grid (pitfree's filter_with_grid rejects any point
  # outside it), normalized on that ground; heights outside -2..120 m are noise
  ex <- ext(tmpl)
  a  <- filter_poi(a, X >= ex[1] & X < ex[2] & Y >= ex[3] & Y < ex[4])
  if (is.empty(a)) return(out)
  nl <- normalize_height(a, dtm, na.rm = TRUE)
  nl <- filter_poi(nl, Z > -2 & Z < 120)
  if (is.empty(nl)) return(out)
  chm <- rasterize_canopy(nl, res = tmpl, algorithm = PITFREE)
  chm <- terra::clamp(chm, lower = 0, values = TRUE)
  rm(a, g, nl)
  for (k in seq_len(nrow(fp))) {
    pts <- cbind(fx[k] + 0.5 + disc$dx, fy[k] + 0.5 + disc$dy)
    e <- terra::extract(dtm, pts); v <- e[[ncol(e)]]; ok <- is.finite(v)
    set(out, k, "lg_valid", sum(ok))
    if (any(ok)) { set(out, k, "lg_mean", mean(v[ok])); set(out, k, "lg_sd", sd(v[ok])) }
    e <- terra::extract(chm, pts); h <- e[[ncol(e)]]; okh <- is.finite(h)
    set(out, k, "chm_valid", sum(okh))
    if (any(okh)) {
      set(out, k, "chm_mean", mean(h[okh]))
      set(out, k, "chm_p90", as.numeric(quantile(h[okh], 0.9, type = 7)))
    }
  }
  out
}

summ <- list()
for (sid in SITES) {
  code <- SITE_CODE[[as.character(sid)]]
  of <- file.path(OUTDIR, sprintf("lp_site_%02d%s.csv.gz", sid, if (LIMIT > 0) "_smoke" else ""))
  if (file.exists(of) && !OVERWRITE) { msg("site %d (%s): output exists, skipped", sid, code); next }
  t0 <- Sys.time()
  files <- list.files(file.path(IMPORTED, code, "LAZ_ground"), pattern = "\\.la[sz]$",
                      full.names = TRUE, ignore.case = TRUE)
  if (!length(files)) { msg("site %d (%s): NO LAZ, skipped", sid, code); next }
  d <- footprints(sid)
  ctg <- suppressWarnings(readLAScatalog(files, progress = FALSE))
  opt_progress(ctg) <- FALSE
  if (is.na(st_crs(ctg))) {
    fb <- CRS_FALLBACK[as.character(sid)]
    if (is.na(fb)) { msg("site %d (%s): LAZ carry no CRS and no fallback, skipped", sid, code); next }
    st_crs(ctg) <- fb
  }
  crs_h <- st_crs(ctg); crs_str <- crs_h$wkt
  xy <- st_coordinates(st_transform(st_as_sf(d, coords = c("lon_lowestmode", "lat_lowestmode"), crs = 4326), crs_h))
  d[, `:=`(x = xy[, 1], y = xy[, 2])]
  hb <- as.data.table(st_drop_geometry(ctg@data))[, .(Min.X, Max.X, Min.Y, Max.Y)]
  inside <- hb[d, on = .(Min.X <= x, Max.X >= x, Min.Y <= y, Max.Y >= y), .(shot_number = i.shot_number),
               mult = "first", nomatch = NULL]
  d <- d[shot_number %in% inside$shot_number]          # only footprints the point cloud can cover
  d[, cell := as.integer(factor(paste(floor(x / 1000), floor(y / 1000))))]
  cells <- sort(unique(d$cell))
  msg("site %d (%s): %s footprints inside the point cloud, %s cells of 1 km, %s LAZ files, CRS %s",
      sid, code, format(nrow(d), big.mark = ","), format(length(cells), big.mark = ","),
      format(length(files), big.mark = ","), if (!is.na(crs_h$epsg)) crs_h$epsg else crs_h$input)

  res <- mclapply(cells, function(cid) {
    setDTthreads(1L)
    fp <- d[cell == cid]
    r <- tryCatch(one_cell(fp, ctg, crs_str), error = function(e) conditionMessage(e))
    if (is.character(r)) data.table(shot_number = fp$shot_number, err = TRUE, errmsg = r)
    else r[, `:=`(err = FALSE, errmsg = NA_character_)]
  }, mc.cores = CORES, mc.preschedule = FALSE)
  R <- rbindlist(res, fill = TRUE)
  n_err <- R[err == TRUE, uniqueN(shot_number)]
  if (n_err) for (m in head(unique(R[err == TRUE, errmsg]), 3)) msg("  cell error: %s", m)
  out <- merge(d[, .(shot_number, x, y, cell, dep_dtm_mean, als_chm_mean, als_chm_p90, rh_98)],
               R[, -c("err", "errmsg")], by = "shot_number", all.x = TRUE)
  out[, tracker := sid]
  fwrite(out, of)

  g  <- out[lg_valid >= 265L & is.finite(dep_dtm_mean), lg_mean - dep_dtm_mean]
  cc <- out[chm_valid >= 265L & is.finite(als_chm_p90) & (chm_p90 >= 2 | als_chm_p90 >= 2)]
  gg <- out[chm_valid >= 265L & chm_p90 >= 2 & is.finite(rh_98)]
  md <- function(x) if (length(x)) median(x) else NA_real_
  rr <- function(a, b) if (length(a) > 30) cor(a, b) else NA_real_
  s <- data.table(tracker = sid, code = code, n_fp = nrow(out),
                  n_ground_ok = sum(out$lg_valid >= 265L, na.rm = TRUE),
                  n_chm_ok = sum(out$chm_valid >= 265L, na.rm = TRUE), n_cell_err = n_err,
                  gap_3dep_med = md(g), gap_3dep_within05 = if (length(g)) mean(abs(g - md(g)) < 0.5) else NA,
                  chm_mean_minus_delivered = md(cc$chm_mean - cc$als_chm_mean),
                  chm_p90_minus_delivered = md(cc$chm_p90 - cc$als_chm_p90),
                  r_chm_delivered = rr(cc$chm_p90, cc$als_chm_p90),
                  p90_minus_gedi = md(gg$chm_p90 - gg$rh_98), r_p90_gedi = rr(gg$chm_p90, gg$rh_98),
                  delivered_p90_minus_gedi = md(gg$als_chm_p90 - gg$rh_98),
                  r_delivered_gedi = rr(gg$als_chm_p90[is.finite(gg$als_chm_p90)], gg$rh_98[is.finite(gg$als_chm_p90)]),
                  minutes = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1))
  summ[[length(summ) + 1L]] <- s
  msg("site %d: ground %s, CHM %s of %s | ground-3DEP %+.2f (%.2f within 0.5 m) | CHM mean/p90 minus delivered %+.2f/%+.2f r %.3f | p90 minus GEDI %+.2f r %.3f (delivered %+.2f r %.3f) | errors %d | %.1f min",
      sid, format(s$n_ground_ok, big.mark = ","), format(s$n_chm_ok, big.mark = ","), format(s$n_fp, big.mark = ","),
      s$gap_3dep_med, s$gap_3dep_within05, s$chm_mean_minus_delivered, s$chm_p90_minus_delivered, s$r_chm_delivered,
      s$p90_minus_gedi, s$r_p90_gedi, s$delivered_p90_minus_gedi, s$r_delivered_gedi, n_err, s$minutes)
}
if (length(summ)) {
  sf_ <- file.path(OUTDIR, if (LIMIT > 0) "lp_summary_smoke.csv" else "lp_summary.csv")
  fwrite(rbindlist(summ), sf_, append = file.exists(sf_))
}
msg("done")
close(LOG)
