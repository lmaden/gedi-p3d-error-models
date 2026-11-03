# merge_dsm_direct.R  --  the directly sampled DSM error into the analysis frame
# Stacks the 18 per-site files of p3d_dsm_direct.R into dsm_direct/dsm_direct_frame.csv.gz and updates the frame:
# dsm_error_mean <- dsm_error_direct (the delivered DSM on the footprint disc minus the lidar surface, Eq. 1); the
# Eq. 4 value (CHM error plus DTM error) is kept as dsm_error_eq4. The error models (04_error_models) were fit to the
# Eq. 4 error; the tables and figures (05_tables_figures) read the updated frame.
# Run once, after p3d_dsm_direct.R, on the local copy of the cluster outputs (P3D_DATA_ROOT, laid out as the cluster's
# $CH1/manuscript_tables). The frame as written by rework_model_prep.R is kept as rework/frame_footprints_eq4dsm.csv.gz.
# Inputs  (P3D_DATA_ROOT) rework/frame_footprints.csv.gz (rework_model_prep.R; DSM error of Eq. 4)
#                         dsm_direct/dsm_site_NN.csv.gz (p3d_dsm_direct.R, one file per site)
# Outputs (P3D_DATA_ROOT) dsm_direct/dsm_direct_frame.csv.gz (site, shot_number, p3d_dsm_mean, p3d_dsm_valid, delta,
#                         dsm_error_direct); rework/frame_footprints.csv.gz (updated in place)
# Checks  one per-site row for every frame footprint and no other; the Eq. 4 error in the per-site files equals the
#         frame's; every footprint has a direct value.
suppressPackageStartupMessages(library(data.table))
ROOT  <- Sys.getenv("P3D_DATA_ROOT", "<LOCAL_DATA_ROOT>")
FRAME <- file.path(ROOT, "rework", "frame_footprints.csv.gz")
SITES <- c(1:9, 11:15, 17:20)                     # tracker ids
fr <- fread(FRAME, colClasses = list(character = "shot_number"))
if ("dsm_error_eq4" %in% names(fr)) stop("the frame already holds dsm_error_eq4; start from the frame of rework_model_prep.R")
dd <- rbindlist(lapply(SITES, function(s) fread(file.path(ROOT, "dsm_direct", sprintf("dsm_site_%02d.csv.gz", s)),
                                                colClasses = list(character = "shot_number"))))
m <- merge(fr[, .(site, shot_number, e4 = dsm_error_mean)],
           dd[, .(site, shot_number, e4_site = dsm_error_mean, dsm_error_direct)], by = c("site", "shot_number"))
stopifnot(nrow(dd) == nrow(fr), nrow(m) == nrow(fr), !anyDuplicated(dd, by = c("site", "shot_number")),
          max(abs(m$e4 - m$e4_site)) == 0, !anyNA(m$dsm_error_direct))
fwrite(dd[, .(site, shot_number, p3d_dsm_mean, p3d_dsm_valid, delta, dsm_error_direct)],
       file.path(ROOT, "dsm_direct", "dsm_direct_frame.csv.gz"))
invisible(file.copy(FRAME, file.path(ROOT, "rework", "frame_footprints_eq4dsm.csv.gz"), overwrite = FALSE))
fr[, dsm_error_eq4 := dsm_error_mean]
fr[dd, on = .(site, shot_number), dsm_error_mean := i.dsm_error_direct]
fwrite(fr, FRAME)
dd[, w := 1 / .N, by = site]
cat(sprintf("frame updated: %s footprints at %d sites | DSM bias %+.3f -> %+.3f m | Eq. 4 sides within 0.1 m at %.2f%% of footprints (%.2f%% with sites weighted equally)\n",
            format(nrow(fr), big.mark = ","), uniqueN(fr$site), mean(fr$dsm_error_eq4), mean(fr$dsm_error_mean),
            100 * mean(abs(dd$delta) < 0.1), 100 * dd[, sum(w * (abs(delta) < 0.1)) / sum(w)]))
