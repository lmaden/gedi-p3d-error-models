#!/usr/bin/env Rscript
# =====================================================================
# rework_model_prep.R  --  one frame, one training sample, three responses, folds
# Paper: the analysis frame, training sample and folds (Sections 2.3 to 2.5, Text S1, Table S1 constants).
#
# Every reference comes from the site's own point cloud
# (lidar_products_rebuild.R), one method at every site:
#   chm_error_mean  P3D CHM disc mean - rebuilt lidar CHM disc mean
#   dtm_error_mean  P3D DTM disc mean - rebuilt lidar ground disc mean
#   dsm_error_mean  chm + dtm = P3D DSM - (lidar ground + lidar CHM) (Eq. 4; the
#                   error models use it; merge_dsm_direct.R later replaces it in the
#                   local copy of the frame with the directly sampled DSM error)
# Vertical datums (standard conversions, not calibrated to 3DEP):
#   site 19 (LTBMU) ellipsoidal NAD83 heights -> NAVD88 with GEOID18 (H = h - N)
#   site 3  (SAWB)  heights on the ITRF2014 ellipsoid -> NAD83(2011): + the
#                   ellipsoid-height shift at each footprint (+1.21 m there;
#                   gap vs 3DEP before the fix -1.13 m)
#
# Frame (per footprint): P3D CHM and slope valid on >= half the disc; rebuilt
# CHM and ground on >= 265 of 529 cells; rebuilt CHM p90 >= 2 m; P3D DTM valid
# >= 0.5; all 17 predictors present; |each error| <= 100 m; lidar ground within
# 5 m of 3DEP around the site's median gap (datum screen; footprints
# without 3DEP cannot be screened and are dropped).
# Predictors: standardized ONCE across all sites of the frame (standardizing
# within each site would hide between-site differences from the fixed effects).
# Training sample: up to CAP footprints per site at random (all of a smaller
# site), so every site weighs about equally; the same rows are scored when a
# site is held out.
# Folds: sites sorted by training rows, dealt into folds of three (snake).
#
# OUTPUT
#   checkpoints/09_rework_prep.rds   list(dat, folds, zconst, frame_n, ...)
#   manuscript_tables/rework/prep_persite.csv, folds.csv, zconst.csv,
#   frame_footprints.csv.gz, prep_YYYYMMDD.log
# RUN (minutes, no Stan)
#   module purge; module load rh9/R/4.5.0; module load rh9/gdal/3.11.0
#   export PROJECT_ROOT=<CLUSTER_ROOT>/chpt1 R_LIBS=<R_LIBRARY>
#   cd $PROJECT_ROOT/scripts/reviewed && Rscript rework_model_prep.R
#   SITES="4,5,..." restricts the site list (e.g. after a site fails its checks)
# =====================================================================

suppressPackageStartupMessages({ library(data.table); library(sf); library(terra) })

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "<CLUSTER_ROOT>/chpt1")
CHKPT  <- file.path(PROJECT_ROOT, "checkpoints")
LPDIR  <- file.path(PROJECT_ROOT, "manuscript_tables", "lidar_products")
ENRDIR <- file.path(PROJECT_ROOT, "data", "enriched_by_site")
OUTTBL <- file.path(PROJECT_ROOT, "manuscript_tables", "rework")
GEOID  <- file.path(PROJECT_ROOT, "data", "geoid", "us_noaa_g2018u0.tif")   # NOAA GEOID18 grid (PROJ data)
CAP    <- as.integer(Sys.getenv("CAP", "5000"))
SEED   <- 2026L
GAP_MAX <- 5     # m; datum screen, lidar ground vs 3DEP around the site median
dir.create(OUTTBL, recursive = TRUE, showWarnings = FALSE)

SITES <- if (nzchar(Sys.getenv("SITES"))) as.integer(strsplit(Sys.getenv("SITES"), ",")[[1]]) else
  c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 11L, 12L, 13L, 14L, 15L, 17L, 18L, 19L, 20L)
# point-cloud CRS of the lp_site files (lidar_products_rebuild.R log) for the datum steps
VERT <- list(`19` = list(type = "geoid18", crs = 3310L),
             `3`  = list(type = "itrf2014", crs = 32619L))

LOG <- file.path(OUTTBL, sprintf("prep_%s.log", format(Sys.Date(), "%Y%m%d")))
con <- file(LOG, open = "wt"); sink(con, split = TRUE)
hr  <- function() cat(strrep("=", 92), "\n")
sec <- function(x) { cat("\n"); hr(); cat("## ", x, "\n"); hr() }
fatal <- function(...) { cat("\nFATAL: ", ..., "\n", sep = ""); sink(); close(con); quit(status = 1) }
cat("rework_model_prep.R\nstarted: ", format(Sys.time()), "\nsites: ", paste(SITES, collapse = " "),
    " | cap ", CAP, " per site\n", sep = "")

RAW <- c("slope_mean", "slope_sd", "wsci", "rh_98", "cover", "aspect_sin_mean", "aspect_cos_mean",
         "meta_off_nadir_avg", "meta_sun_elev_avg", "meta_az_concentration", "meta_stereo_ratio",
         "meta_leaf_on_ratio", "meta_fwd_ratio", "meta_rev_ratio", "meta_abs_geoacc_avg",
         "meta_rel_geoacc_avg", "meta_target_azimuth_avg")

# =====================================================================
sec("STEP 1  --  FOOTPRINTS, LIDAR PRODUCTS, DATUMS")
# =====================================================================
parts <- list()
for (s in SITES) {
  lf <- file.path(LPDIR, sprintf("lp_site_%02d.csv.gz", s))
  if (!file.exists(lf)) fatal("missing ", lf, ". Run lidar_products_rebuild.R first.")
  lp <- fread(lf, select = c("shot_number", "x", "y", "lg_mean", "lg_valid", "chm_mean", "chm_p90", "chm_valid"),
              colClasses = c(shot_number = "character"))
  v <- VERT[[as.character(s)]]
  if (!is.null(v)) {
    ok <- is.finite(lp$x)
    ll <- st_coordinates(st_transform(st_as_sf(lp[ok], coords = c("x", "y"), crs = v$crs), 6318))
    if (v$type == "geoid18") {
      if (!file.exists(GEOID)) fatal("missing ", GEOID)
      e <- terra::extract(rast(GEOID), cbind(ll[, 1] %% 360, ll[, 2])); N <- e[[ncol(e)]]
      if (abs(median(N, na.rm = TRUE) + 24) > 3) fatal("GEOID18 undulation at site ", s, " not near -24 m")
      lp[ok, lg_mean := lg_mean - N]
      cat(sprintf("site %d: GEOID18 applied, N %.2f to %.2f m\n", s, min(N, na.rm = TRUE), max(N, na.rm = TRUE)))
    } else if (v$type == "itrf2014") {
      sh <- sf_project("EPSG:7912", "EPSG:6319", cbind(ll[, 1], ll[, 2], 0))[, 3]
      if (!all(is.finite(sh)) || abs(median(sh) - 1.2) > 0.5) fatal("ITRF2014 to NAD83(2011) shift at site ", s, " not near +1.2 m")
      lp[ok, lg_mean := lg_mean + sh]
      cat(sprintf("site %d: ITRF2014 -> NAD83(2011) ellipsoid shift applied, %+.3f to %+.3f m\n", s, min(sh), max(sh)))
    }
  }
  en <- fread(file.path(ENRDIR, sprintf("site_%02d_enriched.csv.gz", s)),
              select = c("shot_number", "ecoregion", "lc2022_l1_code", RAW, "p3d_chm_mean", "p3d_chm_valid_frac",
                         "slope_valid_frac", "p3d_dtm_mean", "p3d_dtm_valid_frac", "dep_dtm_mean"),
              colClasses = c(shot_number = "character", lc2022_l1_code = "character"))
  en <- unique(en, by = "shot_number")
  x <- merge(en, lp[, -c("x", "y")], by = "shot_number")
  x[, site := s]
  cat(sprintf("site %2d: %s enriched rows, %s with lidar products\n", s,
              format(nrow(en), big.mark = ","), format(nrow(x), big.mark = ",")))
  parts[[length(parts) + 1L]] <- x
}
F <- rbindlist(parts); rm(parts)

# =====================================================================
sec("STEP 2  --  FRAME AND RESPONSES")
# =====================================================================
F[, `:=`(chm_error_mean = p3d_chm_mean - chm_mean, dtm_error_mean = p3d_dtm_mean - lg_mean)]
F[, dsm_error_mean := chm_error_mean + dtm_error_mean]
# datum screen: part of the site 1 delivery (usda_me, x < ~400 km) has its
# ground a flat +25.7 m above 3DEP while P3D agrees with 3DEP there; a blunder
# screen on the departure from the site's median gap, not a calibration to 3DEP
F[, ground_gap := lg_mean - dep_dtm_mean]
F[, gap_med := median(ground_gap, na.rm = TRUE), by = site]
cat("ground minus 3DEP beyond", GAP_MAX, "m of the site median, by site:",
    paste(F[, .(k = sum(abs(ground_gap - gap_med) > GAP_MAX, na.rm = TRUE)), by = site][k > 0, sprintf("%d:%d", site, k)], collapse = " "), "\n")
F[, `:=`(lc_l1_code = fifelse(is.na(lc2022_l1_code) | lc2022_l1_code == "", "UNK", lc2022_l1_code),
         fwdrev = meta_fwd_ratio - meta_rev_ratio,
         view_az_sin = sin(pi * meta_target_azimuth_avg / 180),
         view_az_cos = cos(pi * meta_target_azimuth_avg / 180))]
PRED <- c("slope_mean", "slope_sd", "wsci", "rh_98", "cover", "aspect_sin_mean", "aspect_cos_mean",
          "meta_off_nadir_avg", "meta_sun_elev_avg", "meta_az_concentration", "meta_stereo_ratio",
          "meta_leaf_on_ratio", "fwdrev", "meta_abs_geoacc_avg", "meta_rel_geoacc_avg", "view_az_sin", "view_az_cos")
clauses <- list(
  p3d_chm_valid = F[, p3d_chm_valid_frac >= 0.5 & slope_valid_frac >= 0.5],
  lidar_chm     = F[, chm_valid >= 265 & is.finite(chm_p90) & chm_p90 >= 2],
  lidar_ground  = F[, lg_valid >= 265 & is.finite(lg_mean)],
  ground_datum  = F[, is.finite(ground_gap) & abs(ground_gap - gap_med) <= GAP_MAX],
  p3d_dtm_valid = F[, p3d_dtm_valid_frac >= 0.5 & is.finite(p3d_dtm_mean)],
  predictors    = complete.cases(F[, ..PRED]) & !is.na(F$ecoregion) & F$ecoregion != "",
  errors_100    = F[, is.finite(chm_error_mean) & is.finite(dtm_error_mean) &
                      abs(chm_error_mean) <= 100 & abs(dtm_error_mean) <= 100 & abs(dsm_error_mean) <= 100])
keep <- Reduce(`&`, lapply(clauses, function(z) { z[is.na(z)] <- FALSE; z }))
for (nm in names(clauses)) cat(sprintf("  fails %-14s %9s\n", nm, format(sum(!clauses[[nm]] %in% TRUE), big.mark = ",")))
F <- F[keep]
cat(sprintf("frame: %s footprints at %d sites\n", format(nrow(F), big.mark = ","), uniqueN(F$site)))

# =====================================================================
sec("STEP 3  --  GATES")
# =====================================================================
cat(sprintf("identity max |dsm - chm - dtm|: %g\n", F[, max(abs(dsm_error_mean - chm_error_mean - dtm_error_mean))]))
persite <- F[, .(n_frame = .N,
                 chm_bias = mean(chm_error_mean), dtm_bias = mean(dtm_error_mean), dsm_bias = mean(dsm_error_mean),
                 ground_minus_3dep_med = median(lg_mean - dep_dtm_mean, na.rm = TRUE),
                 ground_minus_3dep_within05 = mean(abs((lg_mean - dep_dtm_mean) - median(lg_mean - dep_dtm_mean, na.rm = TRUE)) < 0.5, na.rm = TRUE),
                 rebuilt_p90_minus_gedi = median(chm_p90 - rh_98), r_p90_gedi = cor(chm_p90, rh_98)),
             by = site][order(site)]
print(persite[, lapply(.SD, function(z) if (is.numeric(z)) round(z, 3) else z)])
fwrite(persite, file.path(OUTTBL, "prep_persite.csv"))
bad <- persite[abs(ground_minus_3dep_med) > 1]
if (nrow(bad) && Sys.getenv("ALLOW_DATUM") != "1")
  fatal("datum gate: lidar ground more than 1 m from 3DEP (median) at site(s) ", paste(bad$site, collapse = " "))
if (uniqueN(F$site) != length(SITES)) fatal("sites lost to the frame: ", paste(setdiff(SITES, unique(F$site)), collapse = " "))

# =====================================================================
sec("STEP 4  --  PREDICTORS STANDARDIZED ACROSS ALL SITES")
# =====================================================================
ZN <- c(slope_mean = "slope_mean_z", slope_sd = "slope_sd_z", wsci = "wsci_z", rh_98 = "rh_98_z", cover = "cover_z",
        aspect_sin_mean = "aspect_sin_z", aspect_cos_mean = "aspect_cos_z", meta_off_nadir_avg = "meta_offnad_z",
        meta_sun_elev_avg = "meta_sunel_z", meta_az_concentration = "meta_az_conc_z", meta_stereo_ratio = "meta_stereo_z",
        meta_leaf_on_ratio = "meta_leafon_z", fwdrev = "meta_fwdrev_z", meta_abs_geoacc_avg = "meta_absgeo_z",
        meta_rel_geoacc_avg = "meta_relgeo_z", view_az_sin = "view_az_sin_z", view_az_cos = "view_az_cos_z")
zconst <- data.table(variable = names(ZN), z = unname(ZN),
                     mean = sapply(names(ZN), function(v) mean(F[[v]])), sd = sapply(names(ZN), function(v) sd(F[[v]])))
for (i in seq_len(nrow(zconst))) F[, (zconst$z[i]) := (get(zconst$variable[i]) - zconst$mean[i]) / zconst$sd[i]]
print(zconst[, .(variable, mean = round(mean, 3), sd = round(sd, 3))])
fwrite(zconst, file.path(OUTTBL, "zconst.csv"))
rz <- F[, .(m = mean(rh_98_z)), by = site][order(site)]
cat("per-site mean of rh_98_z (differs between sites now):", paste(sprintf("%d:%+.2f", rz$site, rz$m), collapse = " "), "\n")

# =====================================================================
sec("STEP 5  --  TRAINING SAMPLE, UP TO CAP FOOTPRINTS PER SITE")
# =====================================================================
set.seed(SEED)
dat <- F[F[, .I[if (.N > CAP) sample.int(.N, CAP) else seq_len(.N)], by = site]$V1]
dat[, `:=`(site = factor(site), ecoregion = factor(as.character(ecoregion)), lc_l1_code = factor(lc_l1_code))]
print(dat[, .(n_train = .N), by = site][order(as.integer(as.character(site)))])
cat(sprintf("training rows: %s\n", format(nrow(dat), big.mark = ",")))
pooled <- dat[, lapply(.SD, function(z) c(mean(z), sd(z))), .SDcols = c("chm_error_mean", "dtm_error_mean", "dsm_error_mean")]
cat("training bias / SD:\n"); print(round(pooled, 3))

# =====================================================================
sec("STEP 6  --  FOLDS OF THREE SITES, AND THE LAND-COVER PRE-FLIGHT")
# =====================================================================
ns <- dat[, .N, by = site][order(-N)]
K <- max(2L, as.integer(round(nrow(ns) / 3)))
snake <- rep(c(seq_len(K), rev(seq_len(K))), length.out = nrow(ns))
ns[, fold := snake]
folds <- ns[, .(site = as.integer(as.character(site)), n = N, fold)][order(fold, -n)]
print(folds[, .(sites = paste(site, collapse = " "), n = sum(n)), by = fold])
fwrite(folds, file.path(OUTTBL, "folds.csv"))
KEEP <- c("ENF", "IWL", "BDF", "RCP", "GRS", "UNK", "DNF", "WTR", "IMP", "EBF")
dat[, lc_chk := fifelse(as.character(lc_l1_code) %in% KEEP, as.character(lc_l1_code), "OTHER")]
cat("\ntraining rows per land-cover group (all sites):\n"); print(dat[, .N, by = lc_chk][order(-N)])
lc_bad <- FALSE
for (k in sort(unique(folds$fold))) {
  hs <- as.character(folds[fold == k]$site)
  miss <- setdiff(unique(dat[as.character(site) %in% hs]$lc_chk), unique(dat[!as.character(site) %in% hs]$lc_chk))
  thin <- dat[!as.character(site) %in% hs, .N, by = lc_chk][N < 30]
  cat(sprintf("  fold %d (sites %s): missing %s | training groups under 30 rows: %s\n", k, paste(hs, collapse = " "),
              if (length(miss)) paste(miss, collapse = " ") else "none",
              if (nrow(thin)) paste(sprintf("%s(%d)", thin$lc_chk, thin$N), collapse = " ") else "none"))
  if (length(miss)) lc_bad <- TRUE
}
dat[, lc_chk := NULL]
if (lc_bad) fatal("a fold would hold out a land-cover group its training lacks")

# =====================================================================
sec("STEP 7  --  SAVE")
# =====================================================================
out <- file.path(CHKPT, "09_rework_prep.rds")
saveRDS(list(dat = dat, folds = folds, zconst = zconst, frame_n = F[, .N, by = site], cap = CAP,
             sites = SITES, built = Sys.time(),
             source = "enriched + lidar_products_rebuild.R (ground, pit-free CHM); global z; capped sample"), out)
fwrite(F[, .(site, shot_number, lc_l1_code, chm_p90, rh_98, chm_error_mean, dtm_error_mean, dsm_error_mean)],
       file.path(OUTTBL, "frame_footprints.csv.gz"))
cat("saved ", out, " (", format(nrow(dat), big.mark = ","), " training rows) and frame_footprints.csv.gz\n", sep = "")
hr(); cat("finished: ", format(Sys.time()), "\nlog: ", LOG, "\n", sep = ""); hr()
sink(); close(con)
