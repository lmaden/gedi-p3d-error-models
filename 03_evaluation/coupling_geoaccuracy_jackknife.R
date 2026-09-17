#!/usr/bin/env Rscript
# =====================================================================
# coupling_geoaccuracy_jackknife.R
#
# Section 3.4: the cross-site association between within-site CHM-DTM error
# coupling and the absolute horizontal geolocation accuracy of the source
# imagery, with the robustness checks the paper reports.
#
# Runs entirely off the deposited per-site table, so it reproduces the
# published values without cluster access:
#   data_derived/coupling_site_level/coupling_site_covariates.csv
#     corr_errP3D_errDTM      within-site Pearson r between CHM and DTM error
#     meta_abs_geoacc_avg_avg site-mean absolute geoaccuracy of the imagery (m)
#     flag_status             "ok" for the 15 sites common to both analyses;
#                             "FLAGGED" for Sites 1 to 3, which the canopy-height
#                             analysis excludes on reference quality
#
# Published values (Section 3.4): r = -0.70, jackknife -0.62 to -0.77, and
# r = -0.70 with the two smallest sites removed.
# =====================================================================

csv <- Sys.getenv("COUPLING_COVARIATES_CSV",
                  "data_derived/coupling_site_level/coupling_site_covariates.csv")
stopifnot(file.exists(csv))

d <- read.csv(csv, stringsAsFactors = FALSE)
d <- d[d$flag_status == "ok" &
       !is.na(d$corr_errP3D_errDTM) & !is.na(d$meta_abs_geoacc_avg_avg), ]
stopifnot(nrow(d) == 15)

r_full <- cor(d$corr_errP3D_errDTM, d$meta_abs_geoacc_avg_avg)
cat(sprintf("cross-site r over %d sites: %.3f\n", nrow(d), r_full))

# site-level jackknife: drop one site, refit
loo <- vapply(seq_len(nrow(d)), function(i)
  cor(d$corr_errP3D_errDTM[-i], d$meta_abs_geoacc_avg_avg[-i]), numeric(1))
cat(sprintf("site-level jackknife: %.3f to %.3f (sign stable: %s)\n",
            min(loo), max(loo), all(sign(loo) == sign(r_full))))
cat(sprintf("  most influential site: %s (dropping it gives r = %.3f)\n",
            d$manuscript_site[which.min(abs(loo))], loo[which.min(abs(loo))]))

# repeat without the two smallest sites by footprint count
small <- order(d$n)[1:2]
cat(sprintf("two smallest sites: %s (n = %d) and %s (n = %d)\n",
            d$manuscript_site[small[1]], d$n[small[1]],
            d$manuscript_site[small[2]], d$n[small[2]]))
cat(sprintf("excluding them: r = %.3f over %d sites\n",
            cor(d$corr_errP3D_errDTM[-small], d$meta_abs_geoacc_avg_avg[-small]),
            nrow(d) - 2L))
