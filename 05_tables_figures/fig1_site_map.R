# fig1_site_map.R  --  Fig. 1, the 18 sites on the Level II ecoregions
# CEC Level II ecoregions in an Albers projection, the 9 study ecoregions in the Paul Tol muted palette
# (colour-blind safe) and the others grey; one marker per site, sized by footprint count in three classes;
# labels are the 1-18 paper numbers (site_names.csv); map cropped to the lower 48; legend below the map,
# 180 mm wide. Site markers are the site centroids of the footprint sampling stage (section2_1 site summary).
# Inputs  (P3D_DERIVED)   site_names.csv
#         (P3D_DATA_ROOT) external/section2_1_site_summary_full_v2.csv (centroids by TRACKER id),
#                         external/na_cec_eco_l2/NA_CEC_Eco_Level2.shp (CEC Level II ecoregions)
#         (P3D_RESULTS)   table1/table1_sites.csv (table1_sites.R: footprints, ecoregion)
# Outputs (P3D_RESULTS)   fig1/fig1_site_map.{png,pdf}
# Checks  18 sites with centroids; prints any site whose plotted ecoregion differs from its Table 1 ecoregion.
suppressPackageStartupMessages({ library(data.table); library(sf); library(ggplot2); library(ggrepel) })
sf_use_s2(FALSE)
DATA_ROOT <- Sys.getenv("P3D_DATA_ROOT", "<LOCAL_DATA_ROOT>")
DERIVED   <- Sys.getenv("P3D_DERIVED", "data_derived")
OUT_ROOT  <- Sys.getenv("P3D_RESULTS", "results")
OUT  <- file.path(OUT_ROOT, "fig1"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

nm  <- fread(file.path(DERIVED, "site_names.csv"))
cen <- fread(file.path(DATA_ROOT, "external", "section2_1_site_summary_full_v2.csv"))
st  <- merge(nm[, .(paper_site, tracker)], cen[, .(tracker = site_id, lat = centroid_lat, lon = centroid_lon)], by = "tracker")
stopifnot(nrow(st) == 18L, !anyNA(st$lat))
t1  <- fread(file.path(OUT_ROOT, "table1", "table1_sites.csv"))

STUDY <- c("ATLANTIC HIGHLANDS", "MIXED WOOD PLAINS", "MISSISSIPPI ALLUVIAL AND SOUTHEAST USA COASTAL PLAINS",
           "SOUTHEASTERN USA PLAINS", "SOUTH CENTRAL SEMIARID PRAIRIES", "MIXED WOOD SHIELD",
           "WEST-CENTRAL SEMIARID PRAIRIES", "WESTERN CORDILLERA", "MEDITERRANEAN CALIFORNIA")
eco <- st_read(file.path(DATA_ROOT, "external", "na_cec_eco_l2", "NA_CEC_Eco_Level2.shp"), quiet = TRUE)
eco <- st_transform(st_make_valid(eco), 4326)
eco <- suppressWarnings(st_crop(eco, st_bbox(c(xmin = -126, ymin = 23.5, xmax = -65.5, ymax = 50.5), crs = 4326)))
eco$l2 <- toupper(eco$NA_L2NAME)

# lower-48 outline and state lines from 'maps'
sdf <- map_data("state")
states <- st_sf(geometry = st_sfc(lapply(split(sdf, sdf$group), function(g) {
  m <- as.matrix(g[, c("long", "lat")]); st_polygon(list(rbind(m, m[1, ]))) }), crs = 4326))
states <- st_make_valid(states)
us <- st_union(states)

eco <- suppressWarnings(st_intersection(eco, us))
eco <- eco[!st_is_empty(eco), ]
st <- merge(st, t1[, .(paper_site, n)], by = "paper_site")
st[, size_cl := factor(ifelse(n < 5000, "Under 5,000", ifelse(n < 50000, "5,000 to 50,000", "Over 50,000")), levels = c("Under 5,000", "5,000 to 50,000", "Over 50,000"))]  # no site has 50,000 to 200,000 footprints, so the top class reads "Over 50,000"
stopifnot(st[n >= 50000, all(n > 200000)])
sites <- st_as_sf(st, coords = c("lon", "lat"), crs = 4326)

# check: ecoregion under each centroid vs Table 1 (dominant over footprints)
hit <- st_join(sites, eco[, "l2"], join = st_intersects)
hit <- as.data.table(hit)[, .(paper_site, l2)]
chk <- merge(hit, t1[, .(paper_site, ecoregion = toupper(ecoregion))], by = "paper_site")
chk[, same := l2 == ecoregion]
cat("centroid ecoregion vs Table 1:\n"); print(chk[same == FALSE | is.na(same)])

alb <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
eco <- st_transform(eco, alb); states <- st_transform(states, alb); sites <- st_transform(sites, alb)
sites$x <- st_coordinates(sites)[, 1]; sites$y <- st_coordinates(sites)[, 2]
eco_main  <- eco[eco$l2 %in% STUDY, ]; eco_main$f <- factor(eco_main$l2, levels = STUDY)
eco_other <- eco[!eco$l2 %in% STUDY, ]
# Paul Tol "muted" qualitative palette (colour-blind safe), assigned so neighbouring ecoregions differ
TOL <- setNames(c("#CC6677", "#88CCEE", "#117733", "#DDCC77", "#AA4499", "#44AA99", "#999933", "#882255", "#332288"), STUDY)
lab <- function(x) sub("Usa", "USA", sub(" And ", " and ", tools::toTitleCase(tolower(x))))

p <- ggplot() +
  geom_sf(data = eco_other, fill = "gray96", color = NA) +
  geom_sf(data = eco_main, aes(fill = f), color = NA, alpha = 0.8) +
  geom_sf(data = states, fill = NA, color = "gray45", linewidth = 0.2) +
  geom_point(data = sites, aes(x, y, size = size_cl), shape = 22, fill = "black", color = "white", stroke = 0.5) +
  scale_size_manual(values = c(1.6, 2.6, 5.2), name = "Footprints", guide = guide_legend(override.aes = list(colour = "black"), order = 2)) +
  geom_text_repel(data = sites, aes(x, y, label = paper_site), size = 3, fontface = "bold",
                  bg.color = "white", bg.r = 0.12, min.segment.length = 0, segment.color = "gray30",
                  segment.size = 0.3, box.padding = 0.35, point.padding = 0.2, max.overlaps = Inf, seed = 42) +
  scale_fill_manual(values = TOL, name = "Level II ecoregion", labels = lab, drop = FALSE,
                    guide = guide_legend(keywidth = unit(9, "pt"), keyheight = unit(9, "pt"), ncol = 2, order = 1, title.position = "top")) +
  coord_sf(crs = alb, expand = FALSE) + labs(x = NULL, y = NULL) +
  theme_void(base_size = 10) +
  theme(legend.title = element_text(size = 8.5, face = "bold"), legend.text = element_text(size = 8),
        legend.position = "bottom", legend.box = "vertical", legend.box.just = "left", legend.spacing.y = unit(1, "mm"),
        plot.margin = margin(4, 4, 4, 4))

ggsave(file.path(OUT, "fig1_site_map.png"), p, width = 180, height = 145, units = "mm", dpi = 600, bg = "white")
ggsave(file.path(OUT, "fig1_site_map.pdf"), p, width = 180, height = 145, units = "mm", bg = "white")
cat("written", file.path(OUT, "fig1_site_map.png"), "\n")
