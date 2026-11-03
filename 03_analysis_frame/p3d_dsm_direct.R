#!/usr/bin/env Rscript
# =====================================================================
# p3d_dsm_direct.R  --  the delivered P3D DSM at every frame footprint
# Paper: the directly sampled DSM error (Eq. 1) of Table 2, Table 3b, Tables S2 to S4, Figs. 2 to 4 and Fig. S2.
#
# The frame's DSM error is the CHM error plus the DTM error (Eq. 4,
# rework_model_prep.R STEP 2), because the frame never sampled the
# delivered DSM. This script samples it on the same disc as
# the P3D CHM and DTM, so the paper can compute the DSM error directly (Eq. 1):
#   disc   the 1,961 cells of 0.5 m whose centres lie within 25 cells (12.5 m)
#          of the centre of the cell holding the footprint, snapped to the
#          raster's own grid (fig5_profiles.R, which reproduced the frame's
#          P3D DTM exactly at 142 footprints)
#   place  x, y of lp_site_NN.csv.gz (point-cloud CRS, from
#          lidar_products_rebuild.R) projected to the P3D CRS, as in fig5_profiles.R
#
# Per frame footprint it writes
#   p3d_dsm_mean, p3d_dsm_valid        delivered DSM: mean, valid cells of 1,961
#   p3d_dtm_mean_re, p3d_dtm_valid_re  delivered DTM on the same disc (alignment gate,
#                                      at every 10th 1 km cell; NA elsewhere)
#   p3d_dtm_mean, p3d_chm_mean         the enriched values the frame used
#   delta    p3d_dsm_mean - p3d_dtm_mean - p3d_chm_mean
#            the direct DSM error is the frame's dsm_error_mean + delta. The
#            height model is the DSM minus the DTM, set to 0 where negative, so
#            delta is 0 or negative when the cells line up.
# Gates, per site, in the log and the summary
#   DTM    recomputed minus enriched: median |diff| near 0 (a wrong CRS, tile
#          set or disc shows at once)
#   delta  share above +0.05 m near 0 (a positive delta cannot come from
#          the clipping, so it means misaligned cells)
#
# Input   $CH1/manuscript_tables/rework/frame_footprints.csv.gz (as written by rework_model_prep.R,
#         whose dsm_error_mean is the Eq. 4 error; not the frame updated by merge_dsm_direct.R)
#         $CH1/manuscript_tables/lidar_products/lp_site_NN.csv.gz
#         $CH1/data/enriched_by_site/site_NN_enriched.csv.gz
#         <LIDAR_ARCHIVE>/<code>/LAZ_ground  (CRS only)
#         $CH2/data/<site>/vricon_raster_50cm/{dsm,dtm}/data/*.tif     (0.5 m)
# Output  $CH1/manuscript_tables/dsm_direct/dsm_site_NN.csv.gz,
#         dsm_direct_summary_v2.csv (per-site gates, with the position path used),
#         dsm_direct_log.txt; merge_dsm_direct.R stacks the per-site files
# Run (16 cores, under 16 GB; reads only, writes only the output folder)
#   module purge; module load rh9/R/4.5.0; module load rh9/gdal/3.11.0
#   export R_LIBS=<R_LIBRARY> CH1=<CLUSTER_ROOT>/chpt1 CH2=<CLUSTER_ROOT>/chpt2 IMPORTED=<LIDAR_ARCHIVE>
#   If terra pixel reads segfault ("memory not mapped") while header reads work, drop a
#   LAStools lib folder from the library path for the run:
#   export LD_LIBRARY_PATH=$(echo "$LD_LIBRARY_PATH" | tr ':' '\n' | grep -v lastools | paste -sd:)
#   export TMPDIR=$CH1/tmp
#   cd $CH1/scripts/reviewed
#   Rscript p3d_dsm_direct.R --sites 4 --limit 300 --cores 4      (smoke, a few minutes)
#   nohup Rscript p3d_dsm_direct.R --cores 16 > dsm_direct.out 2>&1 &
#   Options: --sites (tracker ids), --limit N, --cores N, --cells-per-task N,
#            --dtm-every N, --method window|points, --overwrite,
#            --positions auto|lidar|lidar_sf|gedi_terra|gedi_sf|gedi_nad (default auto)
#   --method points reads point by point as fig5_profiles.R does (slow fallback)
# The smoke run also checks that the windowed read equals fig5_profiles.R's
# point-by-point read at 20 footprints, and stops if it does not.
# =====================================================================

Sys.setenv(GDAL_CACHEMAX = "512")   # MB of raster cache per worker
suppressPackageStartupMessages({ library(data.table); library(sf); library(terra); library(parallel) })
setDTthreads(1L)

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
CH2       <- Sys.getenv("CH2", "<CLUSTER_ROOT>/chpt2")
IMPORTED  <- Sys.getenv("IMPORTED", "<LIDAR_ARCHIVE>")
LIMIT     <- as.integer(getarg("--limit", "0"))
CORES     <- as.integer(getarg("--cores", "16"))
CELLS_PER_TASK <- as.integer(getarg("--cells-per-task", "20"))
DTM_EVERY <- as.integer(getarg("--dtm-every", "10"))   # re-read the DTM at every Nth 1 km cell
METHOD    <- getarg("--method", "window")   # window (fast) or points (fig5_profiles.R's reads, slow)
if (!METHOD %in% c("window", "points")) stop("--method must be window or points")
# How footprints are placed on the P3D grid. "auto" tries every path below on 300
# footprints and keeps the one whose recomputed P3D DTM matches the frame's (the
# lidar path matched exactly at 15 sites but not at the three USFS point clouds,
# trackers 1, 5 and 19, which use older NAD83 coordinate systems).
POSITIONS <- getarg("--positions", "auto")
POS_PATHS <- c("lidar", "lidar_sf", "gedi_terra", "gedi_sf", "gedi_nad")
if (!POSITIONS %in% c("auto", POS_PATHS)) stop("--positions must be auto or one of ", paste(POS_PATHS, collapse = ", "))
OVERWRITE <- "--overwrite" %in% args
SEED      <- 20261014L

# tracker id -> point-cloud folder (as in lidar_products_rebuild.R)
SITE_CODE <- c(`1` = "usda_me", `2` = "nasa_howland", `3` = "neon_sawb",
               `4` = "neon_harv2019", `5` = "usda_sc", `6` = "neon_jerc2019",
               `7` = "neon_tall2021", `8` = "neon_dela2021", `9` = "neon_leno2021",
               `11` = "neon_konz2020", `12` = "neon_stei2022", `13` = "neon_unde2022",
               `14` = "neon_steicheq2022", `15` = "neon_wood2021", `17` = "neon_ster2022",
               `18` = "neon_cper2021", `19` = "ltbmu_20180710", `20` = "neon_sjer2019")
CRS_FALLBACK <- c(`5` = 26917L)

SITES <- getints("--sites"); if (is.null(SITES)) SITES <- as.integer(names(SITE_CODE))
OUTDIR <- file.path(CH1, "manuscript_tables", "dsm_direct")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
LOG <- file(file.path(OUTDIR, "dsm_direct_log.txt"), open = "at")
msg <- function(...) {
  line <- paste0(format(Sys.time(), "%H:%M:%S "), sprintf(...))
  cat(line, "\n", sep = ""); writeLines(line, LOG); flush(LOG)
}
msg("p3d_dsm_direct.R | sites %s | limit %d | cores %d | %d cells per task | DTM gate every %d cells | %s reads%s",
    paste(SITES, collapse = " "), LIMIT, CORES, CELLS_PER_TASK, DTM_EVERY, METHOD, if (LIMIT > 0) " | SMOKE" else "")

disc <- CJ(k = -25:25, l = -25:25)[k^2 + l^2 <= 625]
stopifnot(nrow(disc) == 1961L)
ND <- nrow(disc)

frame <- fread(file.path(CH1, "manuscript_tables", "rework", "frame_footprints.csv.gz"),
               select = c("site", "shot_number", "dsm_error_mean"),
               colClasses = c(shot_number = "character"))
msg("frame: %s footprints at %d sites", format(nrow(frame), big.mark = ","), uniqueN(frame$site))

p3d_files <- function(sid, kind) {
  d <- file.path(CH2, "data", sid, "vricon_raster_50cm", kind, "data")
  f <- list.files(d, pattern = "\\.tif$", full.names = TRUE)
  if (!length(f)) stop("no P3D ", kind, " tiles in ", d)
  f
}

lidar_crs <- function(sid, code) {
  files <- list.files(file.path(IMPORTED, code, "LAZ_ground"), pattern = "\\.la[sz]$",
                      full.names = TRUE, ignore.case = TRUE)
  if (!length(files)) stop("no LAZ files for site ", sid, " (", code, ")")
  ctg <- suppressWarnings(lidR::readLAScatalog(files, progress = FALSE))
  cr <- sf::st_crs(ctg)
  if (is.na(cr)) {
    fb <- CRS_FALLBACK[as.character(sid)]
    if (is.na(fb)) stop("LAZ carry no CRS and no fallback for site ", sid)
    cr <- sf::st_crs(as.integer(fb))
  }
  cr$wkt
}

clean <- function(v) { v[!is.finite(v) | v < -1000 | v > 10000] <- NA_real_; v }
summarise_disc <- function(v) {
  M <- matrix(v, nrow = ND)                       # one column per footprint
  cnt <- colSums(!is.na(M))
  mn <- colMeans(M, na.rm = TRUE); mn[cnt == 0L] <- NA_real_
  list(mean = mn, valid = as.integer(cnt))
}

# Disc means of one raster for footprints at px, py (raster CRS), all within one
# 1 km cell. The centre cell and the disc are those of fig5_profiles.R: the cell
# holding the footprint on the raster's own grid, and the cells k, l away from it.
# One window covering every disc is read once with readValues (no new raster is
# made; terra's crop() segfaulted on the cluster VRT) and indexed in memory.
MAX_WIN <- 4200L    # cells per side; a 1 km cell plus the discs needs about 2,052
disc_means <- function(R, px, py) {
  e <- ext(R); x0 <- e[1]; y0 <- e[3]; rs <- res(R)
  if (any(abs(rs - 0.5) > 1e-6)) stop("raster resolution is not 0.5 m: ", paste(rs, collapse = " "))
  n  <- length(px)
  NR <- nrow(R); NC <- ncol(R)
  ci <- floor((px - x0) / 0.5)                    # 0-based column from the west edge
  ri <- floor((py - y0) / 0.5)                    # 0-based row from the south edge
  c0 <- max(0L, min(ci) - 25L); c1 <- min(NC - 1L, max(ci) + 25L)
  r0 <- max(0L, min(ri) - 25L); r1 <- min(NR - 1L, max(ri) + 25L)
  if (c1 < c0 || r1 < r0) return(list(mean = rep(NA_real_, n), valid = integer(n)))
  if (c1 - c0 + 1L > MAX_WIN || r1 - r0 + 1L > MAX_WIN)
    stop(sprintf("window of %d x %d cells is too large", c1 - c0 + 1L, r1 - r0 + 1L))
  nrw <- r1 - r0 + 1L; ncw <- c1 - c0 + 1L
  readStart(R)
  v <- readValues(R, row = NR - r1, nrows = nrw, col = c0 + 1L, ncols = ncw)
  readStop(R)
  m  <- matrix(v, nrow = nrw, ncol = ncw, byrow = TRUE)   # row 1 is the northern edge
  cc <- rep(ci, each = ND) + rep(disc$k, times = n)       # 0-based, whole raster
  rr <- rep(ri, each = ND) + rep(disc$l, times = n)
  inb <- cc >= c0 & cc <= c1 & rr >= r0 & rr <= r1
  out <- rep(NA_real_, length(cc))
  out[inb] <- m[cbind(r1 - rr[inb] + 1L, cc[inb] - c0 + 1L)]
  summarise_disc(clean(out))
}

# the same means read point by point, exactly as fig5_profiles.R (smoke-test reference)
disc_means_points <- function(R, px, py) {
  e <- ext(R); x0 <- e[1]; y0 <- e[3]
  cx <- x0 + (floor((px - x0) / 0.5) + 0.5) * 0.5
  cy <- y0 + (floor((py - y0) / 0.5) + 0.5) * 0.5
  n <- length(px)
  X <- rep(cx, each = ND) + rep(disc$k * 0.5, times = n)
  Y <- rep(cy, each = ND) + rep(disc$l * 0.5, times = n)
  summarise_disc(clean(terra::extract(R, cbind(X, Y), method = "simple")[, 1]))
}

# Tile extents of one product at one site (header reads only).
tile_index <- function(files) rbindlist(lapply(files, function(f) {
  r <- rast(f); e <- as.vector(ext(r))
  data.table(f = f, xmin = e[1], xmax = e[2], ymin = e[3], ymax = e[4], crs = crs(r))
}))
# The tiles under one 1 km cell's discs, opened as fig5_profiles.R does: one tile
# directly, or the few tiles it needs as a small VRT. A VRT of all of a site's
# tiles segfaulted in terra on the cluster (crop and readValues).
cell_means <- function(TI, px, py, fun) {
  n <- length(px)
  sel <- TI[xmin <= max(px) + 15 & xmax >= min(px) - 15 & ymin <= max(py) + 15 & ymax >= min(py) - 15, f]
  if (!length(sel)) return(list(mean = rep(NA_real_, n), valid = integer(n)))
  if (length(sel) == 1L) return(fun(rast(sel), px, py))
  tf <- tempfile(fileext = ".vrt"); on.exit(unlink(tf))
  fun(vrt(sel, filename = tf, overwrite = TRUE), px, py)
}

# one task = a batch of 1 km cells; the DTM is re-read for every DTM_EVERY-th cell (gate)
one_task <- function(D, TI_dsm, TI_dtm) {
  dm <- if (METHOD == "window") disc_means else disc_means_points
  out <- list()
  for (cid in unique(D$cell)) {
    ch <- D[cell == cid]
    a <- cell_means(TI_dsm, ch$px, ch$py, dm)
    b <- if (cid %% DTM_EVERY == 0L) cell_means(TI_dtm, ch$px, ch$py, dm) else
      list(mean = rep(NA_real_, nrow(ch)), valid = rep(NA_integer_, nrow(ch)))
    out[[length(out) + 1L]] <- data.table(shot_number = ch$shot_number, p3d_dsm_mean = a$mean,
                                          p3d_dsm_valid = a$valid, p3d_dtm_mean_re = b$mean,
                                          p3d_dtm_valid_re = b$valid)
  }
  rbindlist(out)
}

# GEDI's own footprint positions (lowest mode), as lidar_products_rebuild.R reads them
gedi_lonlat <- function(sid) {
  gp <- file.path(CH1, "gedi", sid, sprintf("GEDI_site%d_hq_ALL.gpkg", sid))
  if (!file.exists(gp)) return(NULL)
  fp <- suppressWarnings(sf::st_read(gp, quiet = TRUE, int64_as_string = TRUE,
          query = sprintf("select shot_number, lon_lowestmode, lat_lowestmode from GEDI_site%d_hq_ALL", sid)))
  if (inherits(fp, "sf")) fp <- sf::st_drop_geometry(fp)
  fp <- as.data.table(fp)
  if (!is.character(fp$shot_number)) stop("gpkg shot_number was not read as text")
  unique(fp[!is.na(lon_lowestmode) & !is.na(lat_lowestmode),
            .(shot_number, lon = as.numeric(lon_lowestmode), lat = as.numeric(lat_lowestmode))], by = "shot_number")
}

# candidate positions on the P3D grid, one matrix (x, y) per path
candidates <- function(d, crsL, crsP) {
  tsf <- function(df, from) sf::st_coordinates(sf::st_transform(
    sf::st_as_sf(df, coords = names(df), crs = from), sf::st_crs(crsP)))[, 1:2]
  out <- list(lidar    = crds(project(vect(as.matrix(d[, .(x, y)]), crs = crsL), crsP)),
              lidar_sf = tsf(as.data.frame(d[, .(x, y)]), sf::st_crs(crsL)))
  if (all(c("lon", "lat") %in% names(d)) && any(is.finite(d$lon))) {
    ll <- as.matrix(d[, .(lon, lat)])
    out$gedi_terra <- crds(project(vect(ll, crs = "EPSG:4326"), crsP))
    out$gedi_sf    <- tsf(as.data.frame(d[, .(lon, lat)]), 4326)
    out$gedi_nad   <- crds(project(vect(ll, crs = "EPSG:6318"), crsP))   # lon, lat read as NAD83(2011)
  }
  lapply(out, function(m) unname(as.matrix(m)[, 1:2, drop = FALSE]))
}

# DTM gate for each candidate on a sample: share of footprints whose recomputed
# P3D DTM is within 1 cm of the frame's (the enriched value)
score_positions <- function(d, cand, TI_dtm, n = 300L) {
  set.seed(SEED)
  k <- if (nrow(d) > n) sample.int(nrow(d), n) else seq_len(nrow(d))
  rbindlist(lapply(names(cand), function(nm) {
    P <- cand[[nm]][k, , drop = FALSE]
    s <- data.table(px = P[, 1], py = P[, 2], enr = d$p3d_dtm_mean[k])[is.finite(px) & is.finite(py)]
    if (!nrow(s)) return(data.table(path = nm, n = 0L, within_1cm = NA_real_, median_absdiff = NA_real_))
    s[, cell := as.integer(factor(paste(floor(px / 1000), floor(py / 1000))))]
    m <- rbindlist(lapply(split(s, by = "cell"), function(u)
      data.table(enr = u$enr, re = cell_means(TI_dtm, u$px, u$py, disc_means)$mean)))
    dd <- abs(m$re - m$enr)
    data.table(path = nm, n = sum(is.finite(dd)), within_1cm = mean(dd < 0.01, na.rm = TRUE),
               median_absdiff = median(dd, na.rm = TRUE))
  }))
}

# pre-flight: every input of every site exists before any site starts
pf <- rbindlist(lapply(SITES, function(sid) {
  code <- SITE_CODE[[as.character(sid)]]
  nt <- function(kind) length(list.files(file.path(CH2, "data", sid, "vricon_raster_50cm", kind, "data"), pattern = "\\.tif$"))
  data.table(site = sid, code = if (is.null(code)) NA_character_ else code,
             frame = nrow(frame[site == sid]),
             lp = file.exists(file.path(CH1, "manuscript_tables", "lidar_products", sprintf("lp_site_%02d.csv.gz", sid))),
             enriched = file.exists(file.path(CH1, "data", "enriched_by_site", sprintf("site_%02d_enriched.csv.gz", sid))),
             laz = !is.null(code) && dir.exists(file.path(IMPORTED, code, "LAZ_ground")),
             dsm_tiles = nt("dsm"), dtm_tiles = nt("dtm"))
}))
print(pf)
bad <- pf[is.na(code) | frame == 0L | !lp | !enriched | !laz | dsm_tiles == 0L | dtm_tiles == 0L]
if (nrow(bad)) { msg("PRE-FLIGHT FAILED at site(s) %s", paste(bad$site, collapse = " ")); close(LOG); quit(status = 1) }
msg("pre-flight: all inputs present for %d site(s)", nrow(pf))

summ <- list()
for (sid in SITES) {
  code <- SITE_CODE[[as.character(sid)]]
  if (is.null(code)) { msg("site %d: not a tracker id of the frame, skipped", sid); next }
  of <- file.path(OUTDIR, sprintf("dsm_site_%02d%s.csv.gz", sid, if (LIMIT > 0) "_smoke" else ""))
  if (file.exists(of) && !OVERWRITE) { msg("site %d (%s): output exists, skipped", sid, code); next }
  t0 <- Sys.time()
  fr <- frame[site == sid]
  lp <- fread(file.path(CH1, "manuscript_tables", "lidar_products", sprintf("lp_site_%02d.csv.gz", sid)),
              select = c("shot_number", "x", "y"), colClasses = c(shot_number = "character"))
  en <- fread(file.path(CH1, "data", "enriched_by_site", sprintf("site_%02d_enriched.csv.gz", sid)),
              select = c("shot_number", "p3d_chm_mean", "p3d_dtm_mean"),
              colClasses = c(shot_number = "character"))
  d <- merge(merge(fr, unique(lp, by = "shot_number"), by = "shot_number"),
             unique(en, by = "shot_number"), by = "shot_number")
  if (nrow(d) != nrow(fr))
    msg("site %d: WARNING %d frame footprints lack x, y or enriched values", sid, nrow(fr) - nrow(d))
  if (LIMIT > 0L && nrow(d) > LIMIT) { set.seed(SEED); d <- d[sample.int(.N, LIMIT)] }
  ll <- gedi_lonlat(sid)
  if (!is.null(ll)) d <- merge(d, ll, by = "shot_number", all.x = TRUE)

  crsL <- lidar_crs(sid, code)
  f_dsm <- p3d_files(sid, "dsm"); f_dtm <- p3d_files(sid, "dtm")
  TI_dsm <- tile_index(f_dsm); TI_dtm <- tile_index(f_dtm)
  if (uniqueN(c(TI_dsm$crs, TI_dtm$crs)) > 1L) stop("site ", sid, ": P3D tiles carry more than one CRS")
  crsP <- TI_dsm$crs[1]

  cand <- candidates(d, crsL, crsP)
  if (POSITIONS == "auto") {
    sc <- score_positions(d, cand, TI_dtm)
    msg("site %d: position paths, share of 300 footprints whose recomputed DTM is within 1 cm of the frame's: %s",
        sid, paste(sprintf("%s %.3f (median |diff| %.3f m)", sc$path, sc$within_1cm, sc$median_absdiff), collapse = " | "))
    chosen <- sc$path[which.max(sc$within_1cm)]
    if (max(sc$within_1cm, na.rm = TRUE) < 0.95) msg("site %d: WARNING no path reproduces the frame's DTM (best %.3f)", sid, max(sc$within_1cm, na.rm = TRUE))
  } else chosen <- POSITIONS
  if (is.null(cand[[chosen]])) stop("site ", sid, ": position path ", chosen, " is not available")
  msg("site %d: footprints placed with the %s path", sid, chosen)
  d[, `:=`(px = cand[[chosen]][, 1], py = cand[[chosen]][, 2])]
  d <- d[is.finite(px) & is.finite(py)]
  d[, cell := as.integer(factor(paste(floor(px / 1000), floor(py / 1000))))]
  d <- d[order(cell, py, px)]
  cells <- sort(unique(d$cell))
  d[, task := match(cell, cells) %/% CELLS_PER_TASK]
  msg("site %d (%s): %s footprints in %s cells of 1 km, %d DSM and %d DTM tiles, P3D CRS %s",
      sid, code, format(nrow(d), big.mark = ","), format(length(cells), big.mark = ","),
      length(f_dsm), length(f_dtm), crs(rast(f_dsm[1]), describe = TRUE)$code)

  if (LIMIT > 0L && METHOD == "window") {   # smoke: the windowed read must equal fig5_profiles.R's point reads
    t <- d[cell %in% head(unique(d$cell), 5L)]
    for (kind in c("dsm", "dtm")) {
      TIk <- if (kind == "dsm") TI_dsm else TI_dtm
      p <- rbindlist(lapply(split(t, by = "cell"), function(u) as.data.table(cell_means(TIk, u$px, u$py, disc_means_points))))
      msg("  smoke check %s: point reads done for %d footprints in %d cells", kind, nrow(t), uniqueN(t$cell))
      w <- rbindlist(lapply(split(t, by = "cell"), function(u) as.data.table(cell_means(TIk, u$px, u$py, disc_means))))
      dm <- max(abs(w$mean - p$mean), na.rm = TRUE)
      msg("  smoke check %s: windowed minus point reads, max |diff| %.2e m, valid cells equal %s",
          kind, dm, identical(w$valid, p$valid))
      if (!is.finite(dm) || dm > 1e-6 || !identical(w$valid, p$valid))
        stop("smoke check failed for the ", kind, " at site ", sid)
    }
  }

  res <- mclapply(split(d, by = "task"), function(D) {
    setDTthreads(1L)
    r <- tryCatch(one_task(D, TI_dsm, TI_dtm), error = function(e) conditionMessage(e))
    if (is.character(r)) data.table(shot_number = D$shot_number, err = r) else r[, err := NA_character_]
  }, mc.cores = CORES, mc.preschedule = FALSE)
  crashed <- !vapply(res, is.data.table, logical(1))
  if (any(crashed)) msg("  %d task(s) crashed: %s", sum(crashed), paste(head(unique(as.character(res[crashed])), 2), collapse = " | "))
  R <- rbindlist(res[!crashed], fill = TRUE)
  if (!"err" %in% names(R)) R[, err := NA_character_]
  n_err <- R[!is.na(err), .N] + sum(vapply(split(d, by = "task")[crashed], nrow, integer(1)))
  if (R[!is.na(err), .N]) for (m in head(unique(R[!is.na(err), err]), 3)) msg("  task error: %s", m)
  out <- merge(d[, .(shot_number, dsm_error_mean, p3d_dtm_mean, p3d_chm_mean)],
               R[, -"err"], by = "shot_number", all.x = TRUE)
  out[, `:=`(site = sid, delta = p3d_dsm_mean - p3d_dtm_mean - p3d_chm_mean)]
  out[, dsm_error_direct := dsm_error_mean + delta]
  fwrite(out, of)

  ok   <- out[is.finite(delta)]
  ddtm <- out[is.finite(p3d_dtm_mean_re) & is.finite(p3d_dtm_mean), p3d_dtm_mean_re - p3d_dtm_mean]
  q    <- function(x, p) if (length(x)) as.numeric(quantile(x, p)) else NA_real_
  s <- data.table(
    site = sid, code = code, n_frame = nrow(fr), n_sampled = nrow(out), n_chunk_err = n_err,
    n_dsm_valid_half = sum(out$p3d_dsm_valid >= 981L, na.rm = TRUE),
    dtm_absdiff_median = q(abs(ddtm), 0.5), dtm_absdiff_p99 = q(abs(ddtm), 0.99),
    dtm_within_1cm = if (length(ddtm)) mean(abs(ddtm) < 0.01) else NA_real_,
    delta_min = q(ok$delta, 0), delta_p01 = q(ok$delta, 0.01), delta_median = q(ok$delta, 0.5),
    delta_p99 = q(ok$delta, 0.99), delta_max = q(ok$delta, 1), delta_mean = mean(ok$delta),
    delta_share_above_5cm = mean(ok$delta > 0.05),
    dsm_bias_eq4 = mean(ok$dsm_error_mean), dsm_bias_direct = mean(ok$dsm_error_direct),
    dsm_rmse_eq4 = sqrt(mean(ok$dsm_error_mean^2)), dsm_rmse_direct = sqrt(mean(ok$dsm_error_direct^2)),
    positions = chosen,
    minutes = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1))
  summ[[length(summ) + 1L]] <- s
  msg("site %d: %s sampled, %s with DSM on half the disc | DTM gate median |diff| %.3f m, %.0f%% within 1 cm | delta min %.2f median %.3f max %.2f, %.2f%% above +5 cm | DSM bias %+.2f -> %+.2f, RMSE %.2f -> %.2f | errors %d | %.1f min",
      sid, format(s$n_sampled, big.mark = ","), format(s$n_dsm_valid_half, big.mark = ","),
      s$dtm_absdiff_median, 100 * s$dtm_within_1cm, s$delta_min, s$delta_median, s$delta_max,
      100 * s$delta_share_above_5cm, s$dsm_bias_eq4, s$dsm_bias_direct, s$dsm_rmse_eq4,
      s$dsm_rmse_direct, n_err, s$minutes)
}
if (length(summ)) {
  sf_ <- file.path(OUTDIR, if (LIMIT > 0) "dsm_direct_summary_v2_smoke.csv" else "dsm_direct_summary_v2.csv")   # with the positions column
  fwrite(rbindlist(summ), sf_, append = file.exists(sf_))
}
msg("done")
close(LOG)
