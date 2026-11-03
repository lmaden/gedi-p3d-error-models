# fig5_select_transects.R  --  Fig. 5, step 1 of 3: pick the three transects
# Picks one ~1 km stretch of GEDI track per site by a fixed rule, so each panel
# shows what is typical at its site:
#   window = consecutive shots of one beam spanning at most 18 shot steps (~1 km),
#   at least 12 footprints in the analysed sample, no gap over 3 steps,
#   all footprints in forest classes;
#   score  = sum over DSM, DTM, CHM of ((window mean error - site mean error) / site SD)^2;
#   lowest score wins. The top 5 per site are written so a fallback exists if a
#   window has holes in the rasters.
# Sites (tracker ids): 6 Jones Center, 4 Harvard Forest, 14 Chequamegon (paper Sites 6, 4, 13).
# Inputs  (P3D_DATA_ROOT) rework/frame_footprints.csv.gz;
#                         lidar_ground/lg_site_{06_neon_jerc2019,04,14}.csv.gz (x, y in the point cloud's native CRS)
# The coordinates come from the lidar ground files of the earlier ground-only rebuild (lg_site), which cover
# slightly fewer frame footprints than lidar_products/lp_site_NN.csv.gz (x, y identical where both exist); the
# published windows reproduce only with this footprint set (with the lp files the Site 13 rank-1 window shifts by
# one shot and the Site 6 rank-2 window gains a footprint).
# Outputs (P3D_RESULTS)   fig5/fig5_windows.csv (top 5 windows per site),
#                         fig5/fig5_footprints_siteNN.csv (their footprints),
#                         fig5/fig5_footprints.csv (ranks 1 to 3 of the three sites; input of fig5_profiles.R,
#                         which is run on the cluster with this file in its folder)
suppressPackageStartupMessages(library(data.table))

DATA_ROOT <- Sys.getenv("P3D_DATA_ROOT", "<LOCAL_DATA_ROOT>")
OUT_ROOT  <- Sys.getenv("P3D_RESULTS", "results")
FRAME <- file.path(DATA_ROOT, "rework", "frame_footprints.csv.gz")
LG    <- c(`6`  = "lidar_ground/lg_site_06_neon_jerc2019.csv.gz",
           `4`  = "lidar_ground/lg_site_04.csv.gz",
           `14` = "lidar_ground/lg_site_14.csv.gz")
OUT   <- file.path(OUT_ROOT, "fig5")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

FOREST   <- c("ENF", "DNF", "EBF", "BDF", "MFT")
MAX_SPAN <- 18L
MIN_N    <- 12L
MAX_GAP  <- 3L

fr <- fread(FRAME, colClasses = list(character = "shot_number"))
# The published transects were selected with the DSM error of Eq. 4, before the DSM was sampled directly; with the directly sampled error the same windows rank first at Sites 6 and 13 and second at Site 4.
fr[, dsm_error_mean := dsm_error_eq4]
fr <- fr[site %in% as.integer(names(LG))]

res <- list(); fp3 <- list()
for (s in names(LG)) {
  f  <- fr[site == as.integer(s)]
  xy <- fread(file.path(DATA_ROOT, LG[[s]]), select = c("shot_number", "x", "y"),
              colClasses = list(character = "shot_number"))
  f  <- merge(f, xy, by = "shot_number")
  cat(sprintf("site %s: %d of %d frame footprints have coordinates\n",
              s, nrow(f), nrow(fr[site == as.integer(s)])))
  ref <- fr[site == as.integer(s), .(b_dsm = mean(dsm_error_mean), s_dsm = sd(dsm_error_mean),
                                     b_dtm = mean(dtm_error_mean), s_dtm = sd(dtm_error_mean),
                                     b_chm = mean(chm_error_mean), s_chm = sd(chm_error_mean))]
  # beam track = all digits but the last 8; along-track index = last 8 digits
  nc <- nchar(f$shot_number)
  f[, track := substr(shot_number, 1, nc - 8)]
  f[, idx := as.integer(substr(shot_number, nc - 7, nc))]
  setorder(f, track, idx)
  cand <- list()
  for (tk in unique(f$track)) {
    g <- f[track == tk]
    if (nrow(g) < MIN_N) next
    for (i in seq_len(nrow(g))) {
      w <- g[idx >= g$idx[i] & idx <= g$idx[i] + MAX_SPAN]
      if (nrow(w) < MIN_N || any(diff(w$idx) > MAX_GAP) || !all(w$lc_l1_code %in% FOREST)) next
      step_m <- sqrt(diff(range(w$x))^2 + diff(range(w$y))^2) / max(1, diff(range(w$idx)))
      if (step_m < 50 || step_m > 70) next   # guards the along-track reading of shot_number
      sc <- ((mean(w$dsm_error_mean) - ref$b_dsm) / ref$s_dsm)^2 +
            ((mean(w$dtm_error_mean) - ref$b_dtm) / ref$s_dtm)^2 +
            ((mean(w$chm_error_mean) - ref$b_chm) / ref$s_chm)^2
      cand[[length(cand) + 1]] <- data.table(
        site = as.integer(s), track = tk, idx0 = g$idx[i], idx1 = max(w$idx), n = nrow(w),
        step_m = round(step_m, 1),
        x0 = w$x[1], y0 = w$y[1], x1 = w$x[nrow(w)], y1 = w$y[nrow(w)],
        dsm = mean(w$dsm_error_mean), dtm = mean(w$dtm_error_mean), chm = mean(w$chm_error_mean),
        h_p90 = mean(w$chm_p90), site_dsm = ref$b_dsm, site_dtm = ref$b_dtm, site_chm = ref$b_chm,
        score = sc)
    }
  }
  cand <- rbindlist(cand)
  cat(sprintf("site %s: %d candidate windows\n", s, nrow(cand)))
  # keep non-overlapping windows only, best first
  setorder(cand, score)
  keep <- cand[0]
  for (j in seq_len(nrow(cand))) {
    r <- cand[j]
    clash <- keep[track == r$track & idx0 <= r$idx1 & idx1 >= r$idx0]
    if (nrow(clash) == 0) keep <- rbind(keep, r)
    if (nrow(keep) == 5) break
  }
  keep[, rank := seq_len(.N)]
  res[[s]] <- keep
  # footprints of the top 5 windows, for the cluster profile script
  fp <- rbindlist(lapply(seq_len(nrow(keep)), function(k) {
    f[track == keep$track[k] & idx >= keep$idx0[k] & idx <= keep$idx1[k],
      .(site, rank = k, shot_number, x, y, chm_p90, dsm_error_mean, dtm_error_mean, chm_error_mean)]
  }))
  fwrite(fp, file.path(OUT, sprintf("fig5_footprints_site%02d.csv", as.integer(s))))
  fp3[[s]] <- fp[rank <= 3]
}
res <- rbindlist(res)
fwrite(res, file.path(OUT, "fig5_windows.csv"))
fwrite(rbindlist(fp3), file.path(OUT, "fig5_footprints.csv"))
print(res[, .(site, rank, track, idx0, idx1, n, step_m, dsm = round(dsm, 2), dtm = round(dtm, 2),
              chm = round(chm, 2), site_dsm = round(site_dsm, 2), site_dtm = round(site_dtm, 2),
              site_chm = round(site_chm, 2), h = round(h_p90, 1), score = round(score, 3))])
