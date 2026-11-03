# table1_sites.R  --  Table 1, the 18 sites
# Columns: site (number + formal name), state, Level II ecoregion (dominant over the site's
# footprints), forest share of footprints, lidar (source, year), P3D imagery years,
# footprints, median lidar canopy height.
# Inputs  (P3D_DERIVED)   site_names.csv; ecoregion_by_site.csv (table1_facts.R, cluster);
#                         p3d_year_summary.csv (first and last image year per site and product).
#                         The lidar years are typed in below, as read from lidar_dates.csv (table1_facts.R).
#         (P3D_RESULTS)   persite/persite.csv (persite_class_tables.R: n, median canopy height)
#         (P3D_DATA_ROOT) rework/frame_footprints.csv.gz (forest share)
# Outputs (P3D_RESULTS)   table1/table1_sites.csv and table1_sites.md
# Checks  footprints sum to 630,390; ecoregion counts equal persite n; 9 dominant ecoregions.
suppressPackageStartupMessages(library(data.table))
DATA_ROOT <- Sys.getenv("P3D_DATA_ROOT", "<LOCAL_DATA_ROOT>")
DERIVED   <- Sys.getenv("P3D_DERIVED", "data_derived")
OUT_ROOT  <- Sys.getenv("P3D_RESULTS", "results")
OUT  <- file.path(OUT_ROOT, "table1"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
nm  <- fread(file.path(DERIVED, "site_names.csv"))
eco <- fread(file.path(DERIVED, "ecoregion_by_site.csv"))
ps  <- fread(file.path(OUT_ROOT, "persite", "persite.csv"))
p3d <- fread(file.path(DERIVED, "p3d_year_summary.csv"))[product == "Vricon DSM"]
fr  <- fread(file.path(DATA_ROOT, "rework", "frame_footprints.csv.gz"),
             select = c("site", "lc_l1_code"))

# lidar year per flight, read from lidar_dates.csv (table1_facts.R): GPS dates where
# adjusted GPS time is stored (usda_me Aug-Sep 2014, usda_sc Jun-Jul 2014, ltbmu Jul-Aug 2018;
# ltbmu's 2011-09-14 rows are files with GPS time 0), else the flight date in NEON file names
# (YYYYMMDDHH) or the Howland file names (Jun2012); neon_jerc2019 from the folder year.
lid <- data.table(tracker = c(1:9, 11:15, 17:20),
                  source  = c("USFS", "NASA", "NEON", "NEON", "USFS", rep("NEON", 11), "USFS", "NEON"),
                  year    = c(2014, 2012, 2014, 2019, 2014, 2019, 2021, 2021, 2021, 2020, 2022, 2022, 2022, 2021, 2022, 2021, 2018, 2019))

FOREST <- c("BDF", "ENF", "DNF", "EBF", "MFT")
fs <- fr[, .(forest_pct = 100 * mean(lc_l1_code %in% FOREST)), by = site]
dom <- eco[order(site, -N)][, .SD[1], by = site][, .(site, ecoregion, eco_share = share)]
tt  <- function(x) gsub("\\b(And|Usa)\\b", "\\L\\1", tools::toTitleCase(tolower(x)), perl = TRUE)
dom[, ecoregion := sub("usa", "USA", tt(ecoregion))]

T1 <- nm[, .(paper_site, formal_name_for_table1, state, tracker)]
T1 <- merge(T1, dom, by.x = "tracker", by.y = "site")
T1 <- merge(T1, fs, by.x = "tracker", by.y = "site")
T1 <- merge(T1, lid, by = "tracker")
T1 <- merge(T1, p3d[, .(tracker = site_id, img_first = first_year, img_last = last_year)], by = "tracker")
T1 <- merge(T1, ps[, .(tracker = site, n, median_height)], by = "tracker")
setorder(T1, paper_site)

# checks
g1 <- sum(T1$n) == 630390L
g2 <- all(eco[, .(N = sum(N)), by = site][ps, on = "site"][, N == n])
g3 <- uniqueN(dom$ecoregion) == 9L
cat(sprintf("GATE footprints sum 630,390: %s | ecoregion counts = persite n: %s | 9 dominant ecoregions: %s\n",
            g1, g2, g3))
cat(sprintf("forest share of all footprints: %.1f%%\n", 100 * mean(fr$lc_l1_code %in% FOREST)))
cat("sites whose dominant ecoregion holds < 90% of footprints:\n"); print(T1[eco_share < 0.9, .(paper_site, ecoregion, eco_share)])

fwrite(T1, file.path(OUT, "table1_sites.csv"))
md <- c("| Site | State | Level II ecoregion | Forest (% of footprints) | Lidar (source, year) | P3D imagery (years) | Footprints | Median lidar canopy height (m) |",
        "|---|---|---|---|---|---|---|---|")
md <- c(md, T1[, sprintf("| %d %s | %s | %s | %.0f | %s %d | %d to %d | %s | %.1f |", paper_site, formal_name_for_table1,
                         state, ecoregion, forest_pct, source, year, img_first, img_last,
                         format(n, big.mark = ","), median_height)])
writeLines(md, file.path(OUT, "table1_sites.md"))
cat(md, sep = "\n")
