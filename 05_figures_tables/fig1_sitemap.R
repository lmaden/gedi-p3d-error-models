# =====================================================================
# fig1_sitemap.R
# Figure 1: study area map of the 19 GEDI cal/val sites across CONUS.
#
# The background is the EPA CEC Level II ecoregion shapefile. The nine study
# ecoregions are filled from the Spectral palette at alpha 0.75; ecoregions
# adjacent to the study set are uniform gray at alpha 0.7.
#
# Marker shape encodes which analyses a site enters:
#   filled square        15 sites in both the canopy-height and terrain analyses
#   square with X        Site 10, whose terrain reference is defective
#   open square outline  Sites 1, 2 and 3, outside the canopy-height analysis
#
# Manuscript site IDs 1 to 19 are labelled beside each marker, and the two
# legends stack on the right. Site centroids come from the site summary CSV;
# where one is missing the script falls back to a hardcoded lat/lon for that
# site (Site 16 at 40.49, -103.04, matching neon_ster2022).
#
# Source data:
#   - section2_1_site_summary_full_v2.csv (lat/lon per site)
#   - site_id_lookup.csv (tracker -> manuscript ID; loaded by fig_common.R)
#   - NA_CEC_Eco_Level2.shp (ecoregion polygons)
# =====================================================================

source("fig_common.R")
fig_banner("Figure 1",
                 "Site map with ecoregions, frame-class markers, ID labels")

suppressPackageStartupMessages({
  library(sf); library(ggplot2); library(dplyr); library(data.table)
  library(ggrepel)
})

# ---------------------------------------------------------------------
# Hardcoded fallback centroids (G-6: section2_1_site_summary_full_v2.csv
# has NA centroid for Site 16; values from v09 notebook hardcoded list).
# ---------------------------------------------------------------------
G6_FALLBACK_CENTROIDS <- data.table(
  manuscript_id = c(16L),
  fallback_lat  = c(40.49),
  fallback_lon  = c(-103.04)  # Colorado; matches neon_ster2022 (tracker 17)
)

# ---------------------------------------------------------------------
# Canonical 9 study ecoregion names (uppercase to match notebook's
# NA_L2NAME field values; toupper() applied to shapefile values before
# matching so case differences are absorbed).
# ---------------------------------------------------------------------
STUDY_ECOREGIONS <- c(
  "ATLANTIC HIGHLANDS",
  "MIXED WOOD PLAINS",
  "MISSISSIPPI ALLUVIAL AND SOUTHEAST USA COASTAL PLAINS",
  "SOUTHEASTERN USA PLAINS",
  "SOUTH CENTRAL SEMIARID PRAIRIES",
  "MIXED WOOD SHIELD",
  "WEST-CENTRAL SEMIARID PRAIRIES",
  "WESTERN CORDILLERA",
  "MEDITERRANEAN CALIFORNIA"
)

# ---------------------------------------------------------------------
# State-outline source resolution 
# ---------------------------------------------------------------------
get_states_sf <- function() {
  if (requireNamespace("rnaturalearth",      quietly = TRUE) &&
      requireNamespace("rnaturalearthdata",  quietly = TRUE) &&
      requireNamespace("rnaturalearthhires", quietly = TRUE)) {
    log_progress("  Using 'rnaturalearth' for CONUS state outlines (preferred)")
    s <- try(rnaturalearth::ne_states(country = "United States of America",
                                       returnclass = "sf"),
             silent = TRUE)
    if (!inherits(s, "try-error") && !is.null(s) && nrow(s) > 0) return(s)
    log_progress("  rnaturalearth call failed; trying maps")
  } else if (requireNamespace("rnaturalearth", quietly = TRUE)) {
    log_progress("  rnaturalearth installed but rnaturalearthhires is not; using 'maps' instead")
  }
  if (requireNamespace("maps", quietly = TRUE)) {
    log_progress("  Using 'maps' package for CONUS state outlines")
    states_df <- ggplot2::map_data("state")
    states_sf_list <- split(states_df, states_df$group) |>
      lapply(function(g) {
        m <- as.matrix(g[, c("long", "lat")])
        m <- rbind(m, m[1, ])
        sf::st_polygon(list(m))
      })
    states_sfc <- sf::st_sfc(states_sf_list, crs = 4326)
    out <- sf::st_sf(
      group = unique(states_df$group),
      region = vapply(split(states_df$region, states_df$group),
                      `[`, character(1), 1),
      geometry = states_sfc
    )
    return(out)
  }
  for (fn in c("us_states.shp", "cb_2018_us_state_20m.shp",
               "tl_2020_us_state.shp")) {
    p <- file.path(PROJECT_ROOT, "rasters", fn)
    if (file.exists(p)) {
      log_progress(sprintf("  Using project shapefile: %s", p))
      s <- sf::st_read(p, quiet = TRUE)
      if (!is.na(sf::st_crs(s)) && sf::st_crs(s) != sf::st_crs(4326)) {
        s <- sf::st_transform(s, 4326)
      }
      return(s)
    }
  }
  log_progress("  WARNING: no state-outline data available; drawing points only")
  NULL
}

# ---------------------------------------------------------------------
# 1. Load site coordinates (with G-6 fallback)
# ---------------------------------------------------------------------
log_subsection("Loading site coordinates")

site_summary_path <- file.path(manuscript_tables_dir,
                               "section2_1_site_summary_full_v2.csv")
if (!file.exists(site_summary_path)) {
  stop("Required site summary not found: ", site_summary_path)
}
sites <- as.data.table(fread(site_summary_path))

nms <- tolower(names(sites))
setnames(sites, names(sites), nms)

manuscript_col <- intersect(c("site_id", "manuscript_id", "manuscript_site"),
                            names(sites))[1]
lat_col        <- intersect(c("centroid_lat", "lat", "latitude", "site_lat", "y"),
                            names(sites))[1]
lon_col        <- intersect(c("centroid_lon", "lon", "longitude", "site_lon", "x"),
                            names(sites))[1]

if (any(is.na(c(manuscript_col, lat_col, lon_col)))) {
  stop("Could not locate manuscript/lat/lon columns in ", site_summary_path,
       "\n  Columns present: ", paste(names(sites), collapse = ", "))
}

# v122 FIX (7 Sep 2026): the site_id column of section2_1_site_summary_full_v2.csv
# is the TRACKER id (row 16 = omitted neon_nogp; rows 17-20 = manuscript 16-19).
# v117 assigned it directly as manuscript_id, which drew manuscript sites 17, 18
# and 19 at the tracker 17, 18, 19 locations (CO, CO, Tahoe) and never plotted SJER.
sites[, source_tracker_id := suppressWarnings(as.integer(get(manuscript_col)))]
sites[, manuscript_id     := tracker_to_manuscript(source_tracker_id)]
sites[, lat           := as.numeric(get(lat_col))]
sites[, lon           := as.numeric(get(lon_col))]

sites <- sites[!is.na(manuscript_id) & manuscript_id %in% 1:19]
sites[, tracker_id := manuscript_to_tracker(manuscript_id)]
stopifnot(nrow(sites) == 19L)
stopifnot(all(sites$tracker_id == sites$source_tracker_id))

# G-6 fallback: fill any NA centroids from the hardcoded notebook lookup.
na_rows_before <- sites[is.na(lat) | is.na(lon), .(manuscript_id, tracker_id)]
if (nrow(na_rows_before) > 0L) {
  log_progress(sprintf("  %d site(s) with NA centroid in CSV; attempting G-6 fallback...",
                       nrow(na_rows_before)))
  for (i in seq_len(nrow(G6_FALLBACK_CENTROIDS))) {
    mid <- G6_FALLBACK_CENTROIDS$manuscript_id[i]
    if (any(sites$manuscript_id == mid & (is.na(sites$lat) | is.na(sites$lon)))) {
      sites[manuscript_id == mid & is.na(lat),
            lat := G6_FALLBACK_CENTROIDS$fallback_lat[i]]
      sites[manuscript_id == mid & is.na(lon),
            lon := G6_FALLBACK_CENTROIDS$fallback_lon[i]]
      log_progress(sprintf("    Site %d: filled from G6_FALLBACK_CENTROIDS (lat=%.2f, lon=%.2f)",
                           mid,
                           G6_FALLBACK_CENTROIDS$fallback_lat[i],
                           G6_FALLBACK_CENTROIDS$fallback_lon[i]))
    }
  }
}

# Any sites still missing centroids after fallback get dropped.
still_missing <- sites[is.na(lat) | is.na(lon)]
if (nrow(still_missing) > 0L) {
  log_progress(sprintf("  %d site(s) still missing centroid after fallback; will be dropped:",
                       nrow(still_missing)))
  for (i in seq_len(nrow(still_missing))) {
    log_progress(sprintf("    Manuscript %d (tracker %d)",
                         still_missing$manuscript_id[i],
                         still_missing$tracker_id[i]))
  }
  sites <- sites[!is.na(lat) & !is.na(lon)]
}

# Assign frame class from manuscript ID.
sites[, frame := manuscript_frame_class(manuscript_id)]
sites[, class_label := factor(
  fcase(
    frame == "both",          "Both analyses (15 sites)",
    frame == "dtm_only",      "Terrain analysis only (3 sites)",
    frame == "chm_only",      "Canopy-height analysis only (Site 10)",
    frame == "chm_screened",  "Terrain analysis only (3 sites)",
    frame == "dtm_screened",  "Canopy-height analysis only (Site 10)",
    default = "Other"),
  levels = c("Both analyses (15 sites)",
             "Terrain analysis only (3 sites)",
             "Canopy-height analysis only (Site 10)"))]

log_progress(sprintf("  Loaded %d sites (%s)",
                     nrow(sites),
                     paste(sort(sites$manuscript_id), collapse = ",")))
print(sites[, .(manuscript_id, tracker_id, lat, lon, frame, class_label)])

# ---------------------------------------------------------------------
# 2. Load ecoregion shapefile
# ---------------------------------------------------------------------
log_subsection("Loading ecoregion shapefile")

ECO_PATHS_TO_TRY <- c(
  "/gpfs/data1/vclgp/lmaden/chpt1/data/NA_CEC_Eco_Level2.shp",
  file.path(PROJECT_ROOT, "data", "NA_CEC_Eco_Level2.shp"),
  file.path(PROJECT_ROOT, "rasters", "us_eco_l2.shp"),
  Sys.getenv("CEC_ECO_L2_SHP", "NA_CEC_Eco_Level2.shp")  # local copy of the CEC level-2 ecoregions
)
eco_path <- NULL
for (p in ECO_PATHS_TO_TRY) {
  if (file.exists(p)) {
    eco_path <- p
    break
  }
}
if (is.null(eco_path)) {
  stop("Could not find ecoregion shapefile. Tried:\n  ",
       paste(ECO_PATHS_TO_TRY, collapse = "\n  "))
}
log_progress(sprintf("  Using ecoregion shapefile: %s", eco_path))
eco_sf <- sf::st_read(eco_path, quiet = TRUE)
log_progress(sprintf("  Loaded %d ecoregion polygons", nrow(eco_sf)))

# Identify L2 name column (notebook uses NA_L2NAME; older variants use
# US_L2NAME or L2_KEY).
l2_col <- intersect(c("NA_L2NAME", "US_L2NAME", "L2_NAME", "L2_KEY"),
                    names(eco_sf))[1]
if (is.na(l2_col)) {
  stop("Cannot find an L2 ecoregion name column in ", eco_path,
       "\n  Available columns: ", paste(names(eco_sf), collapse = ", "))
}
log_progress(sprintf("  L2 name column: %s", l2_col))

# Ensure ecoregion sf is in WGS84 for buffer math (5° buffer in lat/lon).
if (!is.na(st_crs(eco_sf)) && st_crs(eco_sf) != st_crs(4326)) {
  eco_sf <- st_transform(eco_sf, 4326)
}

# ---------------------------------------------------------------------
# 3. Filter ecoregions to those near study sites (5° buffer per notebook)
# ---------------------------------------------------------------------
sites_sf_wgs84 <- st_as_sf(sites, coords = c("lon", "lat"), crs = 4326)
sites_buffered <- st_buffer(st_union(sites_sf_wgs84), 5.0)

eco_near_idx <- lengths(st_intersects(eco_sf, sites_buffered)) > 0L
eco_near <- eco_sf[eco_near_idx, ]
log_progress(sprintf("  Ecoregions within 5° of any study site: %d polygons",
                     nrow(eco_near)))

# Split into study (Spectral-colored) vs other (gray).
eco_near$.l2_upper <- toupper(eco_near[[l2_col]])
eco_main  <- eco_near[eco_near$.l2_upper %in%  STUDY_ECOREGIONS, ]
eco_other <- eco_near[!eco_near$.l2_upper %in% STUDY_ECOREGIONS, ]

log_progress(sprintf("  Study ecoregions (Spectral fill): %d polygons (%d unique names)",
                     nrow(eco_main), length(unique(eco_main$.l2_upper))))
log_progress(sprintf("  Near-study ecoregions (gray fill): %d polygons",
                     nrow(eco_other)))

# Convert to a factor with the canonical 9-level ordering for consistent
# Spectral palette assignment across runs.
if (nrow(eco_main) > 0L) {
  eco_main$study_ecoregion <- factor(eco_main$.l2_upper,
                                      levels = STUDY_ECOREGIONS)
}

# ---------------------------------------------------------------------
# 4. State outlines (basemap)
# ---------------------------------------------------------------------
log_subsection("Building base map")
us_states_sf <- get_states_sf()

NON_CONUS <- c("Alaska", "Hawaii", "Puerto Rico",
               "United States Virgin Islands", "Guam", "American Samoa",
               "Northern Mariana Islands",
               "Commonwealth of the Northern Mariana Islands",
               "U.S. Minor Outlying Islands",
               "United States Minor Outlying Islands")
conus_bbox <- c(xmin = -125, xmax = -66, ymin = 24, ymax = 50)
if (!is.null(us_states_sf)) {
  name_col <- intersect(c("name", "name_en", "NAME", "STATE_NAME",
                          "name_long", "admin"),
                        names(us_states_sf))[1]
  if (!is.na(name_col)) {
    us_states_sf <- us_states_sf[!us_states_sf[[name_col]] %in% NON_CONUS, ]
  }
  us_states_sf <- us_states_sf[!sf::st_is_empty(us_states_sf), ]
  us_states_sf <- sf::st_make_valid(us_states_sf)
  cropped <- tryCatch(
    st_crop(us_states_sf,
            st_bbox(c(xmin = unname(conus_bbox["xmin"]),
                      ymin = unname(conus_bbox["ymin"]),
                      xmax = unname(conus_bbox["xmax"]),
                      ymax = unname(conus_bbox["ymax"])),
                    crs = st_crs(us_states_sf))),
    error = function(e) {
      log_progress(sprintf("  st_crop skipped (%s); using name-subset only",
                           conditionMessage(e)))
      NULL
    }
  )
  if (!is.null(cropped)) us_states_sf <- cropped
}

# ---------------------------------------------------------------------
# 5. Project everything to Albers Equal Area Conic for plotting
# ---------------------------------------------------------------------
albers_crs <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"

if (!is.null(us_states_sf)) us_states_sf <- st_transform(us_states_sf, albers_crs)
if (nrow(eco_other) > 0L)    eco_other    <- st_transform(eco_other,    albers_crs)
if (nrow(eco_main) > 0L)     eco_main     <- st_transform(eco_main,     albers_crs)
sites_sf <- st_transform(sites_sf_wgs84, albers_crs)
sites_sf$x <- st_coordinates(sites_sf)[, 1]
sites_sf$y <- st_coordinates(sites_sf)[, 2]

# ---------------------------------------------------------------------
# 6. Build plot (layered: basemap > other-ecos > study-ecos > state
#    outlines on top > site markers > ID labels)
# ---------------------------------------------------------------------
log_subsection("Composing plot")

shape_map <- c(
  "Both analyses (15 sites)"     = 15L,   # solid square
  "Terrain analysis only (3 sites)"   = 0L,    # open square
  "Canopy-height analysis only (Site 10)"   = 7L     # square with X
)
color_map <- c(
  "Both analyses (15 sites)"     = "#222222",
  "Terrain analysis only (3 sites)"   = "#7f7f7f",
  "Canopy-height analysis only (Site 10)"   = "#c0c0c0"
)

p <- ggplot()

# Layer 1: light-gray basemap of CONUS states.
if (!is.null(us_states_sf)) {
  p <- p + geom_sf(data = us_states_sf, fill = "gray97",
                   color = "gray70", linewidth = 0.2)
}

# Layer 2: non-study ecoregions (uniform gray).
if (nrow(eco_other) > 0L) {
  p <- p + geom_sf(data = eco_other,
                   fill = "gray85", color = "gray70",
                   linewidth = 0.15, alpha = 0.7)
}

# Layer 3: study ecoregions (Spectral palette).
if (nrow(eco_main) > 0L) {
  p <- p +
    geom_sf(data = eco_main,
            aes(fill = study_ecoregion),
            color = "gray55", linewidth = 0.18, alpha = 0.75) +
    scale_fill_brewer(palette = "Spectral",
                      name   = "Level II ecoregion",
                      labels = function(x) sub("Usa", "USA", tools::toTitleCase(tolower(x))),
                      drop   = FALSE,
                      guide  = guide_legend(order = 1,
                                            ncol = 1,
                                            byrow = FALSE,
                                            keywidth = unit(10, "pt"),
                                            keyheight = unit(10, "pt")))
}

# Layer 4: state outlines re-drawn on top of ecoregion fills.
if (!is.null(us_states_sf)) {
  p <- p + geom_sf(data = us_states_sf, fill = NA,
                   color = "gray45", linewidth = 0.25)
}

# Layer 5: site markers and ID labels.
p <- p +
  geom_point(data = sites_sf,
             aes(x = x, y = y, shape = class_label, color = class_label),
             size = 3.5, stroke = 1.0) +
  ggrepel::geom_text_repel(
    data = sites_sf,
    aes(x = x, y = y, label = manuscript_id),
    size = 3.0, fontface = "bold",
    color = "black",
    bg.color = "white", bg.r = 0.12,
    min.segment.length = 0,
    segment.color = "gray30", segment.size = 0.30,
    box.padding = 0.30, point.padding = 0.20,
    max.overlaps = Inf,
    seed = 42L
  ) +
  scale_shape_manual(values = shape_map,
                     name = "Analyses",
                     guide = guide_legend(order = 2,
                                          ncol = 1,
                                          keywidth = unit(10, "pt"),
                                          keyheight = unit(10, "pt"))) +
  scale_color_manual(values = color_map,
                     name = "Analyses",
                     guide = guide_legend(order = 2,
                                          ncol = 1,
                                          keywidth = unit(10, "pt"),
                                          keyheight = unit(10, "pt"))) +
  coord_sf(crs = albers_crs, expand = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_section_G(base_size = 11) +
  theme(
    panel.grid       = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position  = "right",
    legend.direction = "vertical",
    legend.box       = "vertical",
    legend.box.just  = "left",
    legend.title     = element_text(size = 9, face = "bold"),
    legend.text      = element_text(size = 8),
    legend.spacing.y = unit(8, "pt"),
    legend.margin    = margin(0, 0, 0, 6),
    axis.text        = element_blank(),
    axis.ticks       = element_blank(),
    plot.margin      = margin(4, 4, 4, 4)
  )

# ---------------------------------------------------------------------
# 7. Save (larger dimensions to accommodate two right-side legends)
# ---------------------------------------------------------------------
save_figure(p, "fig01_sitemap_v122", width_in = 9.5, height_in = 5.5)

# ---------------------------------------------------------------------
# 8. Console verification
# ---------------------------------------------------------------------
log_subsection("VERIFICATION (console anchors)")

cat("\nSite frame class breakdown (must match story-lock v2.3 §0):\n")
print(sites[, .N, by = class_label][order(class_label)])

cat("\nExpected counts:\n")
cat("  Dual-OK (CHM 16 \u2229 DTM 18):     15 sites\n")
cat("  DTM-only (CHM-screened 1,2,3): 3 sites\n")
cat("  CHM-only (DTM-screened 10):    1 site\n")
cat(sprintf("  Total plotted:                 %d sites (target: 19; one omitted if G-6 fallback unavailable)\n",
            nrow(sites)))

cat("\nFull manuscript-ID -> tracker-ID -> frame -> coords mapping:\n")
print(sites[order(manuscript_id),
            .(manuscript_id, tracker_id, frame,
              lat = round(lat, 2), lon = round(lon, 2))])

cat(sprintf("\nEcoregion layer:\n"))
cat(sprintf("  Source: %s\n", eco_path))
cat(sprintf("  Polygons within 5\u00B0 of any study site: %d\n", nrow(eco_near)))
cat(sprintf("  Study ecoregions (Spectral-colored): %d polygons covering %d unique names\n",
            nrow(eco_main), length(unique(eco_main$.l2_upper))))
cat(sprintf("  Other near-study ecoregions (gray): %d polygons\n",
            nrow(eco_other)))

if (nrow(eco_main) > 0L) {
  found <- sort(unique(eco_main$.l2_upper))
  missing_names <- setdiff(STUDY_ECOREGIONS, found)
  extra_names   <- setdiff(found, STUDY_ECOREGIONS)
  if (length(missing_names) == 0L && length(extra_names) == 0L) {
    cat("  All 9 expected study ecoregion names matched, no extras.\n")
  } else {
    if (length(missing_names) > 0L) {
      cat(sprintf("  WARNING: %d expected name(s) NOT found in shapefile:\n",
                  length(missing_names)))
      for (nm in missing_names) cat(sprintf("    - %s\n", nm))
    }
    if (length(extra_names) > 0L) {
      cat(sprintf("  WARNING: %d unexpected name(s) matched STUDY_ECOREGIONS:\n",
                  length(extra_names)))
      for (nm in extra_names) cat(sprintf("    + %s\n", nm))
    }
  }
}

log_progress("section_G_fig1_sitemap.R complete.")
