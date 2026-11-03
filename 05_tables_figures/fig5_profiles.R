# fig5_profiles.R  --  Fig. 5, step 2 of 3: cut the transect profiles (runs on the cluster)
# Input  fig5_footprints.csv (same folder; from fig5_select_transects.R):
#        the footprints of the top 3 windows at sites 6 (Jones Center), 4 (Harvard
#        Forest) and 14 (Chequamegon), x/y in each point cloud's native CRS.
# Rasters
#   lidar  $CH2/data/als_rebuild/site_NN/tiles/tile_<x0>_<y0>.tif (1 km tiles at 1 m from a wall-to-wall
#          run of the lidar_products_rebuild.R recipe, not in this release; band 1 ground_m, band 2 chm_m,
#          native CRS; checked against lidar_products_rebuild.R to the millimeter)
#   P3D    $CH2/data/<sid>/vricon_raster_50cm/{dsm,dtm}/data/*.tif (0.5 m)
# Output $CH1/manuscript_tables/fig5/
#   fig5_profile_siteNN.csv    every 0.5 m along each window (+25 m each end):
#                              P3D DSM, P3D DTM, lidar ground, lidar CHM
#   fig5_fpcheck_siteNN.csv    footprint errors recomputed here vs the frame's
# Gate: recomputed footprint errors must match the frame (median |diff| printed per
# product; a mismatch means wrong tiles, CRS or sign). One core, < 2 GB, ~2 min.
suppressPackageStartupMessages({ library(data.table); library(terra) })

CH1  <- Sys.getenv("CH1_ROOT", "<CLUSTER_ROOT>/chpt1")
CH2  <- Sys.getenv("CH2_ROOT", "<CLUSTER_ROOT>/chpt2")
OUT  <- file.path(CH1, "manuscript_tables", "fig5")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
here <- tryCatch(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))),
                 error = function(e) getwd())
fp_all <- fread(file.path(here, "fig5_footprints.csv"), colClasses = list(character = "shot_number"))

STEP <- 0.5; EXT <- 25
lid_disc <- CJ(dx = -13:13, dy = -13:13)[dx^2 + dy^2 <= 169]          # 529 cells at 1 m
stopifnot(nrow(lid_disc) == 529L)
p3d_disc <- CJ(k = -25:25, l = -25:25)[k^2 + l^2 <= 625]               # 1,961 cells at 0.5 m
stopifnot(nrow(p3d_disc) == 1961L)

p3d_files <- function(sid, kind) {
  d <- file.path(CH2, "data", sid, "vricon_raster_50cm", kind, "data")
  f <- list.files(d, pattern = "\\.tif$", full.names = TRUE)
  if (!length(f)) stop("no P3D ", kind, " tiles in ", d)
  f
}
# tiles whose extent intersects a bbox (header reads only)
pick <- function(files, bb) {
  keep <- vapply(files, function(f) {
    e <- as.vector(ext(rast(f)))
    e[1] <= bb[2] && e[2] >= bb[1] && e[3] <= bb[4] && e[4] >= bb[3]
  }, logical(1))
  files[keep]
}
val <- function(r, xy) extract(r, xy, method = "simple")[, 1]

for (s in unique(fp_all$site)) {
  cat(sprintf("\n==== site %d ====\n", s))
  F   <- fp_all[site == s]
  bbL <- c(min(F$x) - 200, max(F$x) + 200, min(F$y) - 200, max(F$y) + 200)
  # lidar tiles (whole-km grid, native CRS)
  kx  <- seq(floor(bbL[1] / 1000) * 1000, floor(bbL[2] / 1000) * 1000, 1000)
  ky  <- seq(floor(bbL[3] / 1000) * 1000, floor(bbL[4] / 1000) * 1000, 1000)
  lt  <- file.path(CH2, "data", "als_rebuild", sprintf("site_%02d", s), "tiles",
                   sprintf("tile_%d_%d.tif", rep(kx, length(ky)), rep(ky, each = length(kx))))
  lt  <- lt[file.exists(lt)]
  cat(sprintf("lidar tiles: %d\n", length(lt)))
  L   <- if (length(lt) > 1) vrt(lt) else rast(lt)
  crsL <- crs(L)
  # P3D tiles, in the P3D CRS
  f_dsm <- p3d_files(s, "dsm"); f_dtm <- p3d_files(s, "dtm")
  cat(sprintf("P3D files in dirs: dsm %d (e.g. %s), dtm %d (e.g. %s)\n", length(f_dsm),
              basename(f_dsm[1]), length(f_dtm), basename(f_dtm[1])))
  crsP  <- crs(rast(f_dsm[1]))
  bbv   <- vect(cbind(c(bbL[1], bbL[2], bbL[2], bbL[1]), c(bbL[3], bbL[3], bbL[4], bbL[4])),
                crs = crsL)
  bbP   <- as.vector(ext(project(bbv, crsP)))
  f_dsm <- pick(f_dsm, bbP); f_dtm <- pick(f_dtm, bbP)
  cat(sprintf("P3D tiles: dsm %d, dtm %d; P3D CRS %s; lidar CRS %s\n", length(f_dsm), length(f_dtm),
              crs(rast(f_dsm[1]), describe = TRUE)$code, crs(L, describe = TRUE)$code))
  DSM <- if (length(f_dsm) > 1) vrt(f_dsm) else rast(f_dsm)
  DTM <- if (length(f_dtm) > 1) vrt(f_dtm) else rast(f_dtm)
  toP <- function(xy) crds(project(vect(xy, crs = crsL), crsP))

  prof <- list(); chk <- list()
  for (k in sort(unique(F$rank))) {
    W  <- F[rank == k][order(shot_number)]
    a  <- c(W$x[1], W$y[1]); b <- c(W$x[nrow(W)], W$y[nrow(W)])
    len <- sqrt(sum((b - a)^2)); u <- (b - a) / len
    d  <- seq(-EXT, len + EXT, by = STEP)
    xy <- cbind(a[1] + d * u[1], a[2] + d * u[2])
    pp <- toP(xy)
    prof[[k]] <- data.table(site = s, rank = k, dist = d, x = xy[, 1], y = xy[, 2],
                            p3d_dsm = val(DSM, pp), p3d_dtm = val(DTM, pp),
                            l_ground = val(L[[1]], xy), l_chm = val(L[[2]], xy))
    # footprint check: the analysis discs (529 cells of 1 m for lidar, 1,961 cells of 0.5 m for P3D)
    for (j in seq_len(nrow(W))) {
      fx <- floor(W$x[j]); fy <- floor(W$y[j])
      lxy <- cbind(fx + 0.5 + lid_disc$dx, fy + 0.5 + lid_disc$dy)
      g <- val(L[[1]], lxy); h <- val(L[[2]], lxy)
      pc <- toP(cbind(W$x[j], W$y[j]))
      e  <- ext(DSM); x0 <- e[1]; y0 <- e[3]
      cx <- x0 + (floor((pc[1] - x0) / 0.5) + 0.5) * 0.5
      cy <- y0 + (floor((pc[2] - y0) / 0.5) + 0.5) * 0.5
      pxy <- cbind(cx + p3d_disc$k * 0.5, cy + p3d_disc$l * 0.5)
      sv <- val(DSM, pxy); tv <- val(DTM, pxy)
      lg <- mean(g, na.rm = TRUE); lc <- mean(h, na.rm = TRUE)
      ps <- mean(sv, na.rm = TRUE); pt <- mean(tv, na.rm = TRUE)
      chk[[length(chk) + 1]] <- data.table(
        site = s, rank = k, shot_number = W$shot_number[j],
        dist = sum((c(W$x[j], W$y[j]) - a) * u),
        off_line = abs(-(W$x[j] - a[1]) * u[2] + (W$y[j] - a[2]) * u[1]),
        lg_valid = sum(is.finite(g)), p3d_valid = sum(is.finite(sv) & is.finite(tv)),
        dtm_err = pt - lg, chm_err = (ps - pt) - lc, dsm_err = ps - (lg + lc),
        frame_dtm = W$dtm_error_mean[j], frame_chm = W$chm_error_mean[j], frame_dsm = W$dsm_error_mean[j])
    }
  }
  prof <- rbindlist(prof); chk <- rbindlist(chk)
  fwrite(prof, file.path(OUT, sprintf("fig5_profile_site%02d.csv", s)))
  fwrite(chk,  file.path(OUT, sprintf("fig5_fpcheck_site%02d.csv", s)))
  cat(sprintf("profile NA share: p3d_dsm %.3f p3d_dtm %.3f l_ground %.3f l_chm %.3f\n",
              mean(is.na(prof$p3d_dsm)), mean(is.na(prof$p3d_dtm)),
              mean(is.na(prof$l_ground)), mean(is.na(prof$l_chm))))
  cat(sprintf("max footprint distance off the line: %.1f m\n", max(chk$off_line)))
  cat("GATE recomputed minus frame, median |diff| (m):\n")
  print(chk[, .(n = .N, dtm = round(median(abs(dtm_err - frame_dtm)), 3),
                chm = round(median(abs(chm_err - frame_chm)), 3),
                dsm = round(median(abs(dsm_err - frame_dsm)), 3),
                max_any = round(max(abs(c(dtm_err - frame_dtm, chm_err - frame_chm))), 2)), by = rank])
}
cat("\ndone\n")
quit(save = "no")
