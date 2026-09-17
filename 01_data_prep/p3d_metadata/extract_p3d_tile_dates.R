# =====================================================================
# extract_p3d_tile_dates.R
#
# PURPOSE:
#   Read each site_XX_footprints.gpkg, extract ts_utc (P3D tile
#   acquisition dates), and produce a site-level summary CSV.
#
# INPUT:
#   /gpfs/data1/vclgp/lmaden/chpt1/data/footprints/sites/site_XX_footprints.gpkg
#
# OUTPUT:
#   /gpfs/data1/vclgp/lmaden/chpt1/tables/p3d_tile_dates_by_site.csv
#   /gpfs/data1/vclgp/lmaden/chpt1/tables/p3d_tile_dates_all.csv  (every tile)
#
# RUN ON CLUSTER: Rscript extract_p3d_tile_dates.R
# =====================================================================

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(lubridate)
})

# =====================================================================
# CONFIGURATION
# =====================================================================

PROJECT_ROOT  <- Sys.getenv("PROJECT_ROOT", "/gpfs/data1/vclgp/lmaden/chpt1")
FOOTPRINT_DIR <- file.path(PROJECT_ROOT, "data", "footprints", "sites")
out_tables    <- file.path(PROJECT_ROOT, "tables")
dir.create(out_tables, showWarnings = FALSE, recursive = TRUE)

cat("=== EXTRACT P3D TILE DATES ===\n\n")
cat(sprintf("Footprint dir: %s\n", FOOTPRINT_DIR))

# =====================================================================
# DISCOVER GEOPACKAGE FILES
# =====================================================================

gpkg_files <- list.files(FOOTPRINT_DIR,
                          pattern = "^site_[0-9]+_footprints\\.gpkg$",
                          full.names = TRUE)

if (length(gpkg_files) == 0) {
  stop("No geopackage files found in: ", FOOTPRINT_DIR)
}

# Extract site numbers from filenames
get_site_num <- function(path) {
  m <- regmatches(basename(path), regexec("site_([0-9]+)_footprints", basename(path)))
  as.integer(m[[1]][2])
}

gpkg_df <- data.frame(
  path = gpkg_files,
  site = sapply(gpkg_files, get_site_num),
  stringsAsFactors = FALSE
) %>% arrange(site)

cat(sprintf("Found %d geopackage files: sites %s\n\n",
            nrow(gpkg_df),
            paste(gpkg_df$site, collapse = ", ")))

# =====================================================================
# EXTRACT DATES FROM EACH GEOPACKAGE
# =====================================================================

all_tiles <- list()

for (i in seq_len(nrow(gpkg_df))) {
  site_id <- gpkg_df$site[i]
  gpkg_path <- gpkg_df$path[i]

  cat(sprintf("Site %2d: Reading %s... ", site_id, basename(gpkg_path)))

  tryCatch({
    # Read the gpkg — only need ts_utc column + geometry for count
    # st_read with quiet to suppress CRS messages
    layer_name <- st_layers(gpkg_path)$name[1]
    tiles <- st_read(gpkg_path, layer = layer_name, quiet = TRUE)

    if (!"ts_utc" %in% names(tiles)) {
      cat("⚠ No ts_utc column found. Columns: ",
          paste(names(tiles), collapse = ", "), "\n")
      next
    }

    # Parse dates
    tiles$ts_utc_parsed <- ymd_hms(tiles$ts_utc, quiet = TRUE)

    n_total <- nrow(tiles)
    n_valid <- sum(!is.na(tiles$ts_utc_parsed))
    n_na    <- n_total - n_valid

    # Store all tile records (without geometry to save space)
    tile_record <- st_drop_geometry(tiles) %>%
      mutate(
        site = site_id,
        ts_date = as.Date(ts_utc_parsed),
        ts_year = year(ts_utc_parsed)
      ) %>%
      select(site, ts_utc, ts_utc_parsed, ts_date, ts_year,
             any_of(c("veh", "leaf_off", "leaf_on",
                       "off_nadir_avg", "sun_elevation_avg")))

    all_tiles[[as.character(site_id)]] <- tile_record

    cat(sprintf("%d tiles (%d valid dates), range: %s to %s\n",
                n_total, n_valid,
                format(min(tiles$ts_utc_parsed, na.rm = TRUE), "%Y-%m-%d"),
                format(max(tiles$ts_utc_parsed, na.rm = TRUE), "%Y-%m-%d")))

  }, error = function(e) {
    cat(sprintf("ERROR: %s\n", conditionMessage(e)))
  })
}

# =====================================================================
# COMBINE AND SUMMARIZE
# =====================================================================

cat("\n--- Summarizing ---\n\n")

all_tiles_df <- bind_rows(all_tiles)

cat(sprintf("Total tiles across all sites: %d\n", nrow(all_tiles_df)))

# Site-level summary
site_summary <- all_tiles_df %>%
  filter(!is.na(ts_utc_parsed)) %>%
  group_by(site) %>%
  summarise(
    n_tiles        = n(),
    p3d_date_min   = min(ts_date, na.rm = TRUE),
    p3d_date_max   = max(ts_date, na.rm = TRUE),
    p3d_date_median = median(ts_date, na.rm = TRUE),
    p3d_date_mean  = as.Date(mean(ts_utc_parsed, na.rm = TRUE)),
    p3d_year_min   = min(ts_year, na.rm = TRUE),
    p3d_year_max   = max(ts_year, na.rm = TRUE),
    p3d_span_years = as.numeric(difftime(max(ts_utc_parsed, na.rm = TRUE),
                                          min(ts_utc_parsed, na.rm = TRUE),
                                          units = "days")) / 365.25,
    .groups = "drop"
  ) %>%
  arrange(site)

cat("\nSite-level P3D date summary:\n")
print(as.data.frame(site_summary), row.names = FALSE)

# =====================================================================
# SAVE (before optional diagnostics so outputs are guaranteed)
# =====================================================================

# Site-level summary
out_summary <- file.path(out_tables, "p3d_tile_dates_by_site.csv")
write.csv(site_summary, out_summary, row.names = FALSE)
cat(sprintf("\n✓ Site summary saved: %s\n", out_summary))

# All individual tile records (drop parsed datetime and any list columns)
out_all <- file.path(out_tables, "p3d_tile_dates_all.csv")
tiles_out <- all_tiles_df %>%
  select(site, ts_utc, ts_date, ts_year) %>%
  mutate(across(everything(), unlist))
write.csv(tiles_out, out_all, row.names = FALSE)
cat(sprintf("✓ All tiles saved:   %s\n", out_all))

# Year distribution per site (cosmetic — wrapped in tryCatch)
tryCatch({
  cat("\nP3D tile year distribution:\n")
  year_dist <- data.frame(
    site = unlist(all_tiles_df$site),
    ts_year = unlist(all_tiles_df$ts_year)
  ) %>%
    filter(!is.na(ts_year)) %>%
    count(site, ts_year) %>%
    tidyr::pivot_wider(names_from = ts_year, values_from = n, values_fill = 0)
  print(as.data.frame(year_dist), row.names = FALSE)
}, error = function(e) {
  cat(sprintf("  (Year distribution table skipped: %s)\n", conditionMessage(e)))
})
cat(sprintf("✓ All tiles saved:   %s\n", out_all))

cat("\n=== DONE ===\n")
