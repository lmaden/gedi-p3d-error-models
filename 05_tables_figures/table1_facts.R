# table1_facts.R  --  Table 1, two facts only the cluster holds (runs on the cluster)
# 1. Level II ecoregion of every footprint in the analysed sample (frame_footprints.csv.gz
#    joined to the enriched tables on shot_number, read as text): per site the dominant
#    ecoregion and its share, the second one if any, and the count of distinct ecoregions.
# 2. Lidar acquisition dates from the point clouds themselves: 10 random LAZ files per flight,
#    header (GPS time type, file creation date) and the GPS time of every 50th point.
#    Adjusted standard GPS time converts to a calendar date; GPS week time cannot, so for
#    those files only the header creation date is reported.
# Inputs  $CH1/manuscript_tables/rework/frame_footprints.csv.gz (rework_model_prep.R; site and shot_number only)
#         $CH1/data/enriched_by_site/site_NN_enriched.csv.gz (footprint sampling stage)
#         <LIDAR_ARCHIVE>/<flight>/LAZ_ground/ (point clouds)
# Output  $CH1/manuscript_tables/table1/{ecoregion_by_site,lidar_dates}.csv
# One core, < 4 GB, about 5-15 minutes (the LAZ reads dominate).
# RUN: Rscript table1_facts.R > table1_facts.out 2>&1
suppressPackageStartupMessages({ library(data.table); library(lidR) })

CH1 <- Sys.getenv("PROJECT_ROOT", "<CLUSTER_ROOT>/chpt1")
OUT <- file.path(CH1, "manuscript_tables", "table1"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
IMPORTED <- "<LIDAR_ARCHIVE>"
SITE_CODE <- c(`1` = "usda_me", `2` = "nasa_howland", `3` = "neon_sawb",
               `4` = "neon_harv2019", `5` = "usda_sc", `6` = "neon_jerc2019",
               `7` = "neon_tall2021", `8` = "neon_dela2021", `9` = "neon_leno2021",
               `11` = "neon_konz2020", `12` = "neon_stei2022", `13` = "neon_unde2022",
               `14` = "neon_steicheq2022", `15` = "neon_wood2021", `17` = "neon_ster2022",
               `18` = "neon_cper2021", `19` = "ltbmu_20180710", `20` = "neon_sjer2019")

# ---- 1. ecoregions ----
fr <- fread(file.path(CH1, "manuscript_tables", "rework", "frame_footprints.csv.gz"),
            select = c("site", "shot_number"), colClasses = list(character = "shot_number"))
cat(sprintf("frame: %d footprints at %d sites\n", nrow(fr), uniqueN(fr$site)))
eco <- rbindlist(lapply(sort(unique(fr$site)), function(s) {
  en <- fread(file.path(CH1, "data", "enriched_by_site", sprintf("site_%02d_enriched.csv.gz", s)),
              select = c("shot_number", "ecoregion"), colClasses = list(character = c("shot_number", "ecoregion")))
  m <- merge(fr[site == s], unique(en, by = "shot_number"), by = "shot_number")
  if (nrow(m) != nrow(fr[site == s])) cat(sprintf("  site %d: %d of %d joined\n", s, nrow(m), nrow(fr[site == s])))
  m[, .N, by = .(site, ecoregion)]
}))
eco[, share := N / sum(N), by = site]
setorder(eco, site, -N)
fwrite(eco, file.path(OUT, "ecoregion_by_site.csv"))
cat("\nLevel II ecoregion by site (share of footprints):\n")
print(eco[, .(ecoregions = paste(sprintf("%s %.0f%%", ecoregion, 100 * share), collapse = ", ")), by = site])
tot <- eco[, .(N = sum(N)), by = ecoregion][order(-N)]
cat(sprintf("\ndistinct ecoregions over the sample: %d\n", nrow(tot))); print(tot)
cat(sprintf("distinct ecoregions that are some site's dominant one: %d\n",
            uniqueN(eco[, .SD[1], by = site]$ecoregion)))

# ---- 2. lidar dates ----
set.seed(2026)
gps0 <- as.POSIXct("1980-01-06", tz = "UTC")
dates <- rbindlist(lapply(names(SITE_CODE), function(s) {
  d <- file.path(IMPORTED, SITE_CODE[[s]], "LAZ_ground")
  f <- list.files(d, pattern = "\\.la[sz]$", full.names = TRUE, recursive = TRUE)
  if (!length(f)) return(data.table(site = as.integer(s), flight = SITE_CODE[[s]], n_files = 0L))
  pick <- sample(f, min(10L, length(f)))
  rbindlist(lapply(pick, function(p) {
    h <- tryCatch(readLASheader(p), error = function(e) NULL)
    if (is.null(h)) return(data.table(site = as.integer(s), flight = SITE_CODE[[s]], n_files = length(f), file = basename(p), err = "header"))
    phb <- h@PHB
    adj <- isTRUE(as.logical(phb[["Global Encoding"]][["GPS Time Type"]]))
    t <- tryCatch(readLAS(p, select = "t", filter = "-keep_every_nth 50")$gpstime, error = function(e) NULL)
    t0 <- if (length(t)) min(t) else NA_real_; t1 <- if (length(t)) max(t) else NA_real_
    data.table(site = as.integer(s), flight = SITE_CODE[[s]], n_files = length(f), file = basename(p),
               gps_adjusted = adj, created_year = phb[["File Creation Year"]],
               created_doy = phb[["File Creation Day of Year"]],
               date_first = if (adj && is.finite(t0)) format(gps0 + t0 + 1e9, "%Y-%m-%d") else NA_character_,
               date_last  = if (adj && is.finite(t1)) format(gps0 + t1 + 1e9, "%Y-%m-%d") else NA_character_,
               gps_min = t0, gps_max = t1, err = "")
  }), fill = TRUE)
}), fill = TRUE)
fwrite(dates, file.path(OUT, "lidar_dates.csv"))
cat("\nlidar acquisition dates (10 random LAZ files per flight):\n")
print(dates[, .(n_files = n_files[1], adjusted = paste(unique(gps_adjusted), collapse = "/"),
                first = suppressWarnings(min(date_first, na.rm = TRUE)),
                last = suppressWarnings(max(date_last, na.rm = TRUE)),
                created = paste(unique(sprintf("%s-%03d", created_year, created_doy)), collapse = " ")),
            by = .(site, flight)])
cat("\ndone\n")
quit(save = "no")
