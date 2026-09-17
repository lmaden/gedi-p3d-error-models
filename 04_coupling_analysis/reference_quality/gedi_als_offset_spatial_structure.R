# =====================================================================
# gedi_als_offset_spatial_structure.R
#
# The spatial-structure step: within-site spatial structure of per-footprint
# GEDI-vs-ALS p90 offset distributions.
#
# the offset-distribution step (gedi_als_offset_distributions.R) characterized
# per-site offset distributions and showed that:
#   - per-site median offset correlates with Ch1 CHM site RE at r = 0.90
#     (19 sites) / 0.71 (16 non-flagged)
#   - per-site offset SD does NOT correlate with residual REs
# the offset-distribution step ruled out within-site GEDI-ALS dispersion as a candidate for
# the unexplained ~12% sigma_site remainder.
#
# the spatial-structure step is descriptive / Chapter-2-design-input rather than a probe of
# the residual variance. It characterizes whether the within-site offset
# field is spatially structured (clustered / autocorrelated) or random
# in space, which matters for Chapter 2's GWR anchor model design and
# matches the analytic depth of the Site 5 analysis.
#
# Three deliverables:
#   1. Per-site Moran's I on the offset_p90 field (k-NN weights, k=8)
#   2. Per-site empirical variogram + spherical model fit
#   3. Per-site point maps for the 4 Task-1 sites + 3 flagged sites
#
# Scope:
#   - Forest footprints (mirrors the offset-distribution step's forest subset)
#   - All 19 manuscript sites (subject to GEDI gpkg availability)
#   - Per-site coordinate attachment from gedi/GEDI_site<TR>_hq_ALL.gpkg
#     (joined on shot_number)
#
# Outputs (CSV) -> $PROJECT_ROOT/manuscript_tables/:
#   groundwork_task3_phase2_morans_i.csv
#   groundwork_task3_phase2_variogram_params.csv
#   groundwork_task3_phase2_summary.csv
# Outputs (PDF) -> $PROJECT_ROOT/plots/groundwork/:
#   groundwork_task3_phase2_morans_i_summary.pdf
#   groundwork_task3_phase2_variograms.pdf
#   groundwork_task3_phase2_point_map_site<NN>_<short>.pdf  (one per
#                                                            POINT_MAP site)
# Checkpoint: groundwork_task3_phase2 (via save_checkpoint())
#
# Runtime estimate: ~10-25 minutes (per-site read + spdep + gstat at
# ~8k subsample per site is fast; full ingest of 19 enriched CSVs
# dominates).
# =====================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(sf)
  library(spdep)
  library(gstat)
  library(sp)
})

if (!exists("log_progress")) source("analysis_utils.R")
if (!exists("PROJECT_ROOT")) source("analysis_config.R")

log_progress("=== The spatial-structure step: within-site spatial structure ===")

MS_TABLES_DIR   <- file.path(PROJECT_ROOT, "manuscript_tables")
GROUNDWORK_PLOT <- file.path(PROJECT_ROOT, "plots", "groundwork")
GEDI_DIR        <- gedi_base_dir   # from analysis_config.R

dir.create(GROUNDWORK_PLOT, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("  PROJECT_ROOT = %s\n", PROJECT_ROOT))
cat(sprintf("  ENRICHED_DIR = %s\n", ENRICHED_DIR))
cat(sprintf("  GEDI_DIR     = %s\n", GEDI_DIR))
cat(sprintf("  getwd()      = %s\n", getwd()))

.required_dirs <- c(
  "PROJECT_ROOT" = PROJECT_ROOT,
  "ENRICHED_DIR" = ENRICHED_DIR,
  "GEDI_DIR"     = GEDI_DIR,
  "manuscript_tables" = MS_TABLES_DIR
)
.missing <- .required_dirs[!vapply(.required_dirs, dir.exists, logical(1))]
if (length(.missing)) {
  stop("Required directories missing:\n  ",
       paste(names(.missing), "->", .missing, collapse = "\n  "))
}

# ---- Configuration ----

# Subsample size per site for Moran's I + variogram. Full N (up to ~600k
# at the largest sites) is intractable for distance-matrix construction.
# 8000 random points per site is large enough to characterize spatial
# structure cleanly and small enough for spdep::knearneigh + gstat::variogram
# to run in seconds.
N_SAMPLE_PER_SITE <- 8000L

# k-NN for Moran's I weights. k=8 is the standard choice for GEDI at this
# spacing (~60 m between adjacent shots within a track, ~600 m between
# tracks). k=8 captures both within-track and adjacent-track neighbors.
KNN_K <- 8L

# Variogram parameters. 5 km cutoff covers within-site structure for
# every site; 100 m bins resolve the GEDI footprint geometry.
VARIO_CUTOFF_M <- 5000
VARIO_WIDTH_M  <- 100

# Sites for point maps: 4 Task-1 anchors + 3 flagged sites (with overlap
# at site 3 = neon_sawb, which is both flagged AND a Task-1 anchor for
# Mechanism A). Six unique sites total.
POINT_MAP_SITES <- c(1L, 2L, 3L, 4L, 5L, 6L)
# 1 = usda_me (FLAGGED), 2 = nasa_howland (FLAGGED), 3 = neon_sawb
# (FLAGGED + Task-1), 4 = HARV (Task-1), 5 = usda_sc (Task-1),
# 6 = JERC (Task-1)

# Forest LC codes (mirrors the offset-distribution step)
FOREST_LCS <- c("BDF", "DNF", "EBF", "ENF")

# ---- Load site lookup ----

lookup <- data.table::fread(file.path(MS_TABLES_DIR, "site_id_lookup.csv"))
log_progress(sprintf("Loaded site_id_lookup.csv (%d rows)", nrow(lookup)))

lookup_ms <- lookup[!is.na(manuscript_site)]

# ---- Helper: read enriched CSV + coords for one site ----
#
# Mirrors the filter cascade from the offset-distribution step to keep the subset comparable.
# Adds coord attachment from the per-site GEDI gpkg.

.read_site_phase2 <- function(tracker_id) {
  ms_row <- lookup_ms[tracker_site == tracker_id]
  if (!nrow(ms_row)) {
    log_progress(sprintf("  WARNING: tracker %d not in lookup; skipping",
                         tracker_id))
    return(NULL)
  }
  ms_id <- ms_row$manuscript_site[1]
  short <- ms_row$site_short_name[1]
  flag  <- ms_row$flag_status[1]
  dlc   <- ms_row$dominant_lc[1]

  # Locate enriched CSV
  patterns <- c(
    file.path(ENRICHED_DIR, sprintf("site_%02d_enriched.csv.gz", tracker_id)),
    file.path(ENRICHED_DIR, sprintf("site_%d_enriched.csv.gz",   tracker_id)),
    file.path(ENRICHED_DIR, sprintf("site_%02d_enriched.csv",    tracker_id)),
    file.path(ENRICHED_DIR, sprintf("site_%d_enriched.csv",      tracker_id))
  )
  path <- patterns[file.exists(patterns)][1]
  if (is.na(path)) {
    log_progress(sprintf("  WARNING: no enriched CSV for tracker %d",
                         tracker_id))
    return(NULL)
  }

  hdr <- names(data.table::fread(path, nrows = 0))
  needed <- c(
    "shot_number", "site",
    "rh_98", "als_chm_p90", "cover", "error_mean",
    "lc2022_l1_code", "lc2022_l1_name",
    "lc2022_mode_l1_code", "lc2022_mode_l1_name",
    "als_chm_valid_frac"
  )
  have <- intersect(needed, hdr)
  DT <- data.table::fread(path, select = have,
                          colClasses = list(character = "shot_number"),
                          nThread = min(4L, CPU_BUDGET),
                          showProgress = FALSE)

  # LC column normalization (matches section_01_ingest.R)
  if (!("lc2022_l1_code" %in% names(DT)) &&
      "lc2022_mode_l1_code" %in% names(DT)) {
    setnames(DT, "lc2022_mode_l1_code", "lc2022_l1_code")
  }
  if (!("lc2022_l1_name" %in% names(DT)) &&
      "lc2022_mode_l1_name" %in% names(DT)) {
    setnames(DT, "lc2022_mode_l1_name", "lc2022_l1_name")
  }

  # the offset-distribution step offset definition
  DT[, offset_p90 := rh_98 - als_chm_p90]

  # the offset-distribution step filter cascade
  DT_clean <- DT[
    is.finite(rh_98) & is.finite(als_chm_p90) & is.finite(offset_p90) &
    is.finite(als_chm_valid_frac) & als_chm_valid_frac >= 0.5
  ]

  # Forest subset
  if ("lc2022_l1_code" %in% names(DT_clean)) {
    DT_forest <- DT_clean[lc2022_l1_code %in% FOREST_LCS]
  } else {
    DT_forest <- DT_clean[als_chm_p90 >= chm_forest_thresh_m]
  }

  if (nrow(DT_forest) == 0) {
    log_progress(sprintf("  [tracker %d / ms %s] %s: 0 forest footprints; skipping",
                         tracker_id, as.character(ms_id), short))
    return(NULL)
  }

  # ---- Coord attachment from GEDI gpkg ----

  # Canonical layout (per section_01_ingest.R):
  #   <GEDI_DIR>/<tracker_id>/GEDI_site<tracker_id>_hq_ALL.gpkg
  # Fallbacks below cover zero-padded subdir, flat layout (no subdir),
  # and zero-padded filename - in case gpkg files were ever staged
  # differently for some sites.
  gpkg_candidates <- c(
    file.path(GEDI_DIR, as.character(tracker_id),
              sprintf(gedi_file_name, as.character(tracker_id))),
    file.path(GEDI_DIR, sprintf("%02d", tracker_id),
              sprintf(gedi_file_name, as.character(tracker_id))),
    file.path(GEDI_DIR, sprintf(gedi_file_name, as.character(tracker_id))),
    file.path(GEDI_DIR, sprintf("GEDI_site%02d_hq_ALL.gpkg", tracker_id))
  )
  gpkg_path <- gpkg_candidates[file.exists(gpkg_candidates)][1]
  if (is.na(gpkg_path)) {
    log_progress(sprintf("  WARNING: no GEDI gpkg for tracker %d (tried: %s)",
                         tracker_id, paste(basename(gpkg_candidates),
                                           collapse = ", ")))
    return(NULL)
  }

  gp <- tryCatch(sf::st_read(gpkg_path, quiet = TRUE),
                 error = function(e) {
                   log_progress(sprintf("  ERROR reading %s: %s",
                                        gpkg_path, e$message))
                   NULL
                 })
  if (is.null(gp) || !nrow(gp)) {
    log_progress(sprintf("  WARNING: empty/unreadable gpkg for tracker %d",
                         tracker_id))
    return(NULL)
  }

  # Find shot_number column in gpkg (case may vary)
  shot_col <- intersect(c("shot_number", "shotNumber", "SHOT_NUMBER",
                          "shot_no", "ShotNumber"), names(gp))
  if (!length(shot_col)) {
    log_progress(sprintf("  ERROR: no shot_number column in gpkg %s; cols: %s",
                         basename(gpkg_path),
                         paste(names(gp), collapse = ", ")))
    return(NULL)
  }
  gp$shot_number <- as.character(gp[[shot_col[1]]])

  # Project to per-site UTM derived from data centroid (works for any
  # source CRS as long as one is set; if missing, assume WGS84).
  if (is.na(sf::st_crs(gp))) {
    sf::st_crs(gp) <- 4326
  }
  if (sf::st_crs(gp)$epsg != 4326L) {
    gp_4326 <- sf::st_transform(gp, 4326)
  } else {
    gp_4326 <- gp
  }
  cent <- suppressWarnings(sf::st_centroid(sf::st_union(gp_4326)))
  cent_xy <- sf::st_coordinates(cent)
  utm_zone <- floor((cent_xy[1, "X"] + 180) / 6) + 1
  utm_north <- cent_xy[1, "Y"] >= 0
  utm_epsg <- if (utm_north) 32600L + utm_zone else 32700L + utm_zone
  gp_utm <- sf::st_transform(gp_4326, crs = utm_epsg)

  utm_coords <- sf::st_coordinates(gp_utm)
  coords_dt <- data.table(
    shot_number = gp_utm$shot_number,
    x_utm = utm_coords[, "X"],
    y_utm = utm_coords[, "Y"]
  )

  # Join coords onto filtered footprints
  setkey(DT_forest, shot_number)
  setkey(coords_dt, shot_number)
  DT_join <- coords_dt[DT_forest, on = "shot_number", nomatch = 0]

  n_joined <- nrow(DT_join)
  if (n_joined < 100L) {
    log_progress(sprintf("  WARNING: tracker %d only %d footprints after coord join; skipping",
                         tracker_id, n_joined))
    return(NULL)
  }

  DT_join[, `:=`(
    tracker_site    = tracker_id,
    manuscript_site = ms_id,
    site_short_name = short,
    flag_status     = flag,
    dominant_lc     = dlc
  )]

  log_progress(sprintf("  tracker %d / ms %02d %s: n=%s forest footprints with coords (utm_epsg=%d)",
                       tracker_id, ms_id, short,
                       format(n_joined, big.mark = ","),
                       utm_epsg))

  list(data = DT_join, utm_epsg = utm_epsg)
}

# ---- Helper: Moran's I on a per-site subsample ----

.compute_morans_i <- function(d) {
  if (nrow(d) > N_SAMPLE_PER_SITE) {
    set.seed(2025)
    d <- d[sample(.N, N_SAMPLE_PER_SITE)]
  }
  k_use <- min(KNN_K, nrow(d) - 1L)
  if (k_use < 2L) {
    return(list(I = NA_real_, expected_I = NA_real_,
                var_I = NA_real_, z = NA_real_, p = NA_real_,
                n = nrow(d)))
  }

  pts <- as.matrix(d[, .(x_utm, y_utm)])
  knn <- spdep::knearneigh(pts, k = k_use)
  nb  <- spdep::knn2nb(knn)
  lw  <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)

  mt <- tryCatch(
    spdep::moran.test(d$offset_p90, lw, zero.policy = TRUE,
                      randomisation = TRUE),
    error = function(e) NULL
  )
  if (is.null(mt)) {
    return(list(I = NA_real_, expected_I = NA_real_,
                var_I = NA_real_, z = NA_real_, p = NA_real_,
                n = nrow(d)))
  }

  list(
    I          = unname(mt$estimate["Moran I statistic"]),
    expected_I = unname(mt$estimate["Expectation"]),
    var_I      = unname(mt$estimate["Variance"]),
    z          = unname(mt$statistic),
    p          = mt$p.value,
    n          = nrow(d)
  )
}

# ---- Helper: empirical variogram + spherical fit on a per-site subsample ----

.compute_variogram <- function(d) {
  if (nrow(d) > N_SAMPLE_PER_SITE) {
    set.seed(2025)
    d <- d[sample(.N, N_SAMPLE_PER_SITE)]
  }
  if (nrow(d) < 100L) {
    return(list(empirical = NULL, fit = NULL))
  }
  sp_pts <- sp::SpatialPointsDataFrame(
    coords = as.matrix(d[, .(x_utm, y_utm)]),
    data = data.frame(offset = d$offset_p90)
  )
  vgm_emp <- tryCatch(
    gstat::variogram(offset ~ 1, data = sp_pts,
                     cutoff = VARIO_CUTOFF_M, width = VARIO_WIDTH_M),
    error = function(e) NULL
  )
  if (is.null(vgm_emp) || !nrow(vgm_emp)) {
    return(list(empirical = NULL, fit = NULL))
  }

  # Initialize spherical model from empirical curve
  init_sill   <- max(vgm_emp$gamma, na.rm = TRUE)
  valid_g     <- vgm_emp$gamma[vgm_emp$np > 0]
  init_nugget <- if (length(valid_g)) min(valid_g) else 0
  init_range  <- VARIO_CUTOFF_M / 4

  vgm_init <- gstat::vgm(psill  = max(0, init_sill - init_nugget),
                         model  = "Sph",
                         range  = init_range,
                         nugget = init_nugget)
  vgm_fit <- tryCatch(
    gstat::fit.variogram(vgm_emp, vgm_init),
    error   = function(e) NULL,
    warning = function(w) NULL
  )

  list(empirical = vgm_emp, fit = vgm_fit)
}

# =====================================================================
# Per-site loop
# =====================================================================

log_progress("Reading per-site data + coords + computing spatial structure...")

all_trackers <- sort(unique(lookup_ms$tracker_site))
per_site_results <- list()

for (t in all_trackers) {
  log_progress(sprintf("--- tracker site %d ---", t))
  res <- .read_site_phase2(t)
  if (is.null(res)) next

  d <- res$data
  morans <- .compute_morans_i(d)
  vario  <- .compute_variogram(d)

  per_site_results[[as.character(t)]] <- list(
    data     = d,
    utm_epsg = res$utm_epsg,
    morans   = morans,
    vario    = vario
  )

  log_progress(sprintf("  Moran's I = %s (p = %s, n = %d)",
                       formatC(morans$I, digits = 3, format = "f"),
                       formatC(morans$p, digits = 3, format = "g"),
                       morans$n))
}

if (!length(per_site_results)) {
  stop("No sites successfully processed; check enriched CSV and gpkg paths.")
}

log_progress(sprintf("Successfully processed %d sites", length(per_site_results)))

# =====================================================================
# Build summary tables
# =====================================================================

# ---- Moran's I table ----

morans_rows <- lapply(per_site_results, function(r) {
  s <- r$morans
  data.table(
    manuscript_site = r$data$manuscript_site[1],
    tracker_site    = r$data$tracker_site[1],
    site_short_name = r$data$site_short_name[1],
    flag_status     = r$data$flag_status[1],
    dominant_lc     = r$data$dominant_lc[1],
    morans_I        = s$I,
    expected_I      = s$expected_I,
    var_I           = s$var_I,
    z_score         = s$z,
    p_value         = s$p,
    n_subsample     = s$n,
    n_total         = nrow(r$data),
    utm_epsg        = r$utm_epsg
  )
})
morans_dt <- rbindlist(morans_rows, fill = TRUE)
morans_dt <- morans_dt[order(manuscript_site)]
fwrite(morans_dt, file.path(MS_TABLES_DIR, "groundwork_task3_phase2_morans_i.csv"))
log_progress("Wrote groundwork_task3_phase2_morans_i.csv")

# ---- Variogram parameters table ----

vgm_rows <- lapply(per_site_results, function(r) {
  fit <- r$vario$fit
  ms_id <- r$data$manuscript_site[1]
  base <- data.table(
    manuscript_site = ms_id,
    tracker_site    = r$data$tracker_site[1],
    site_short_name = r$data$site_short_name[1],
    flag_status     = r$data$flag_status[1],
    nugget          = NA_real_,
    partial_sill    = NA_real_,
    total_sill      = NA_real_,
    range_m         = NA_real_,
    n_subsample     = r$morans$n
  )
  if (is.null(fit) || nrow(fit) < 1L) {
    return(base)
  }
  nug_idx <- which(fit$model == "Nug")
  sph_idx <- which(fit$model == "Sph")
  nug <- if (length(nug_idx)) fit$psill[nug_idx[1]] else 0
  sph <- if (length(sph_idx)) fit$psill[sph_idx[1]] else NA_real_
  rng <- if (length(sph_idx)) fit$range[sph_idx[1]] else NA_real_
  base[, `:=`(
    nugget       = nug,
    partial_sill = sph,
    total_sill   = nug + ifelse(is.na(sph), 0, sph),
    range_m      = rng
  )]
  base
})
vgm_dt <- rbindlist(vgm_rows, fill = TRUE)
vgm_dt <- vgm_dt[order(manuscript_site)]
fwrite(vgm_dt, file.path(MS_TABLES_DIR, "groundwork_task3_phase2_variogram_params.csv"))
log_progress("Wrote groundwork_task3_phase2_variogram_params.csv")

# ---- Combined summary ----

summary_dt <- merge(
  morans_dt[, .(manuscript_site, tracker_site, site_short_name,
                flag_status, dominant_lc, morans_I, p_value,
                n_subsample, n_total)],
  vgm_dt[, .(manuscript_site, nugget, partial_sill, total_sill, range_m)],
  by = "manuscript_site", all.x = TRUE
)
summary_dt[, nugget_to_sill := nugget / total_sill]
summary_dt <- summary_dt[order(manuscript_site)]
fwrite(summary_dt, file.path(MS_TABLES_DIR, "groundwork_task3_phase2_summary.csv"))
log_progress("Wrote groundwork_task3_phase2_summary.csv")

# ---- Save checkpoint ----

save_checkpoint("groundwork_task3_phase2", list(
  per_site_results = per_site_results,
  morans_dt        = morans_dt,
  vgm_dt           = vgm_dt,
  summary_dt       = summary_dt,
  config = list(
    n_sample_per_site = N_SAMPLE_PER_SITE,
    knn_k             = KNN_K,
    vario_cutoff_m    = VARIO_CUTOFF_M,
    vario_width_m     = VARIO_WIDTH_M,
    forest_lcs        = FOREST_LCS,
    point_map_sites   = POINT_MAP_SITES
  ),
  timestamp = Sys.time()
))
log_progress("Saved checkpoint groundwork_task3_phase2")

# =====================================================================
# FIGURES
# =====================================================================

# ---- Figure 1: cross-site Moran's I summary ----

log_progress("Figure 1: Moran's I cross-site summary...")

plot_morans <- copy(morans_dt)
plot_morans[, label := sprintf("Site %02d %s%s",
                                manuscript_site, site_short_name,
                                ifelse(flag_status == "FLAGGED", " ★", ""))]
plot_morans[, label := factor(label, levels = label[order(morans_I)])]
plot_morans[, fill_grp := ifelse(flag_status == "FLAGGED", "FLAGGED", "OK")]

p_morans <- ggplot(plot_morans, aes(x = morans_I, y = label, fill = fill_grp)) +
  geom_col() +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  geom_text(aes(label = sprintf("I = %.3f, p = %.2g (n = %s)",
                                 morans_I, p_value,
                                 format(n_subsample, big.mark = ","))),
            hjust = -0.05, size = 3) +
  scale_fill_manual(values = c("FLAGGED" = "#B06DB3", "OK" = "#2C6E9B"),
                    name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.50))) +
  labs(
    x = sprintf("Moran's I on offset_p90 (rh_98 − als_chm_p90, k = %d NN weights)", KNN_K),
    y = NULL,
    title = "The spatial-structure step Fig 1: Within-site spatial autocorrelation of GEDI − ALS offset",
    subtitle = sprintf("Per-site %s-footprint subsample. Positive I = clustered in space; near zero = spatially random.",
                       format(N_SAMPLE_PER_SITE, big.mark = ","))
  ) +
  theme_cowplot(font_size = 10) +
  theme(legend.position = "top")

ggsave(file.path(GROUNDWORK_PLOT, "groundwork_task3_phase2_morans_i_summary.pdf"),
       p_morans, width = 11, height = 9)
log_progress("  wrote groundwork_task3_phase2_morans_i_summary.pdf")

# ---- Figure 2: empirical variogram facet ----

log_progress("Figure 2: empirical variograms by site...")

vario_long <- rbindlist(lapply(per_site_results, function(r) {
  emp <- r$vario$empirical
  if (is.null(emp) || !nrow(emp)) return(NULL)
  data.table(
    manuscript_site = r$data$manuscript_site[1],
    site_short_name = r$data$site_short_name[1],
    flag_status     = r$data$flag_status[1],
    dist_m          = emp$dist,
    gamma           = emp$gamma,
    np              = emp$np
  )
}), fill = TRUE)

# Spherical model fit overlays
fit_long <- rbindlist(lapply(per_site_results, function(r) {
  fit <- r$vario$fit
  if (is.null(fit)) return(NULL)
  d_eval <- seq(0, VARIO_CUTOFF_M, length.out = 100)
  pred <- tryCatch(
    gstat::variogramLine(fit, dist_vector = d_eval),
    error = function(e) NULL
  )
  if (is.null(pred)) return(NULL)
  data.table(
    manuscript_site = r$data$manuscript_site[1],
    site_short_name = r$data$site_short_name[1],
    flag_status     = r$data$flag_status[1],
    dist_m          = pred$dist,
    gamma_fit       = pred$gamma
  )
}), fill = TRUE)

# Build panel labels and apply consistent ordering
vario_long[, panel_label := sprintf("Site %02d %s%s",
                                     manuscript_site, site_short_name,
                                     ifelse(flag_status == "FLAGGED", " ★", ""))]
panel_levels <- unique(vario_long[order(manuscript_site), panel_label])
vario_long[, panel_label := factor(panel_label, levels = panel_levels)]

if (nrow(fit_long)) {
  fit_long[, panel_label := sprintf("Site %02d %s%s",
                                     manuscript_site, site_short_name,
                                     ifelse(flag_status == "FLAGGED", " ★", ""))]
  fit_long[, panel_label := factor(panel_label, levels = panel_levels)]
}

p_vario <- ggplot(vario_long, aes(x = dist_m, y = gamma)) +
  geom_point(size = 1.2, alpha = 0.8, colour = "#2C6E9B") +
  geom_line(linewidth = 0.4, alpha = 0.5, colour = "#2C6E9B")

if (nrow(fit_long)) {
  p_vario <- p_vario +
    geom_line(data = fit_long, aes(x = dist_m, y = gamma_fit),
              colour = "#D45E5E", linewidth = 0.6, inherit.aes = FALSE)
}

p_vario <- p_vario +
  facet_wrap(~ panel_label, scales = "free_y", ncol = 5) +
  scale_x_continuous(labels = label_number(scale = 1e-3, suffix = " km")) +
  labs(
    x = "Lag distance (km)",
    y = "Semivariance of offset_p90 (m²)",
    title = "The spatial-structure step Fig 2: Empirical variograms of GEDI − ALS offset",
    subtitle = sprintf("Per-site subsample (n ≤ %s). Blue = empirical; red = spherical model fit. Spherical params in groundwork_task3_phase2_variogram_params.csv.",
                       format(N_SAMPLE_PER_SITE, big.mark = ","))
  ) +
  theme_cowplot(font_size = 10) +
  theme(
    strip.background = element_rect(fill = "grey95", colour = NA),
    strip.text = element_text(size = 8, lineheight = 0.9),
    panel.spacing = unit(0.3, "lines")
  )

ggsave(file.path(GROUNDWORK_PLOT, "groundwork_task3_phase2_variograms.pdf"),
       p_vario, width = 14, height = 10)
log_progress("  wrote groundwork_task3_phase2_variograms.pdf")

# ---- Figure 3: per-site point maps for selected sites ----

log_progress("Figure 3: per-site point maps for selected sites...")

selected_trackers <- lookup_ms[manuscript_site %in% POINT_MAP_SITES]$tracker_site
n_maps_written <- 0L

for (t in selected_trackers) {
  if (is.null(per_site_results[[as.character(t)]])) {
    log_progress(sprintf("  point map skipped for tracker %d (no data)", t))
    next
  }
  r <- per_site_results[[as.character(t)]]
  d <- r$data
  ms_id <- d$manuscript_site[1]
  short <- d$site_short_name[1]
  flag  <- d$flag_status[1]

  # Symmetric percentile clamp on color limits
  q_lim <- quantile(abs(d$offset_p90), 0.98, na.rm = TRUE)
  q_lim <- max(q_lim, 1)  # ensure non-zero

  p_map <- ggplot(d, aes(x = x_utm, y = y_utm, colour = offset_p90)) +
    geom_point(size = 0.4, alpha = 0.6) +
    coord_equal() +
    scale_colour_distiller(type = "div", palette = "RdBu",
                           limits = c(-q_lim, q_lim),
                           oob = scales::squish,
                           name = "rh_98 − als_p90 (m)") +
    labs(
      x = sprintf("Easting (m, EPSG %s)", r$utm_epsg),
      y = "Northing (m)",
      title = sprintf("The spatial-structure step Fig 3: Site %02d %s%s",
                      ms_id, short,
                      ifelse(flag == "FLAGGED", " ★ FLAGGED", "")),
      subtitle = sprintf("n = %s footprints; Moran's I = %.3f (p = %.2g, k=%d NN, n_sub = %s)",
                         format(nrow(d), big.mark = ","),
                         r$morans$I, r$morans$p,
                         KNN_K,
                         format(r$morans$n, big.mark = ","))
    ) +
    theme_cowplot(font_size = 10)

  out_pdf <- file.path(GROUNDWORK_PLOT,
                       sprintf("groundwork_task3_phase2_point_map_site%02d_%s.pdf",
                               ms_id, short))
  ggsave(out_pdf, p_map, width = 8, height = 7)
  log_progress(sprintf("  wrote %s", basename(out_pdf)))
  n_maps_written <- n_maps_written + 1L
}

# =====================================================================
# CONSOLE SUMMARY
# =====================================================================

cat("\n\n====================  SUMMARY  ====================\n\n")

cat("Per-site Moran's I + variogram parameters (manuscript-numbered):\n\n")
print_cols <- c("manuscript_site", "site_short_name", "flag_status",
                "morans_I", "p_value", "range_m", "nugget", "total_sill",
                "nugget_to_sill", "n_subsample", "n_total")
print(summary_dt[, ..print_cols], row.names = FALSE, digits = 3)

cat("\n\nReading guide:\n")
cat("  - morans_I near 0    : within-site offset is spatially random.\n")
cat("  - morans_I > 0       : offset clusters in space (positive autocorrelation).\n")
cat("  - range_m            : distance over which spatial dependence operates.\n")
cat("  - nugget_to_sill     : share of variance NOT explained by spatial structure\n")
cat("                         (1.0 = pure noise; 0 = fully spatially structured).\n")

cat("\nArtifacts written:\n")
cat(sprintf("  %s\n", file.path(MS_TABLES_DIR, "groundwork_task3_phase2_morans_i.csv")))
cat(sprintf("  %s\n", file.path(MS_TABLES_DIR, "groundwork_task3_phase2_variogram_params.csv")))
cat(sprintf("  %s\n", file.path(MS_TABLES_DIR, "groundwork_task3_phase2_summary.csv")))
cat(sprintf("  %s (checkpoint)\n",
            file.path(Sys.getenv("PROJECT_ROOT", PROJECT_ROOT),
                      "checkpoints", "groundwork_task3_phase2.rds")))
cat(sprintf("  %s\n", file.path(GROUNDWORK_PLOT, "groundwork_task3_phase2_morans_i_summary.pdf")))
cat(sprintf("  %s\n", file.path(GROUNDWORK_PLOT, "groundwork_task3_phase2_variograms.pdf")))
cat(sprintf("  groundwork_task3_phase2_point_map_site<NN>_<short>.pdf  (%d files)\n",
            n_maps_written))

cat("\n===========================================================\n\n")

log_progress("the spatial-structure step complete.")
