# =====================================================================
# section_01_ingest.R (FIXED VERSION)
# Data ingest and enrichment from per-site CSVs
#
# FIXES APPLIED:
#   1. Product-specific site exclusion (SITES_EXCLUDE_CHM, SITES_EXCLUDE_DTM)
#   2. Z-scaled azimuth variables (view_az_sin_z, view_az_cos_z)
# =====================================================================

# Ensure dependencies loaded
if (!exists("log_progress")) source("analysis_utils.R")
if (!exists("PROJECT_ROOT")) source("analysis_config.R")
log_progress("Starting data ingest...")

# Parse site specifications
.parse_sites <- function(s) {
  s <- trimws(s); if (!nzchar(s)) return(integer())
  parts <- strsplit(s, ",", fixed = TRUE)[[1]]
  out <- integer(0)
  for (p in parts) {
    p <- trimws(p)
    if (grepl("-", p)) {
      ab <- as.integer(strsplit(p, "-", fixed = TRUE)[[1]])
      if (length(ab) == 2 && !any(is.na(ab))) out <- c(out, seq.int(ab[1], ab[2]))
    } else {
      v <- suppressWarnings(as.integer(p))
      if (!is.na(v)) out <- c(out, v)
    }
  }
  unique(out)
}

# Discover per-site files
log_progress("Discovering per-site enriched CSVs...")
site_file_tbl <- (function() {
  stopifnot(dir.exists(ENRICHED_DIR))
  pats <- c("^site_([0-9]{1,2})_enriched\\.csv\\.gz$",
            "^site_([0-9]{1,2})_enriched\\.csv$")
  f <- list.files(ENRICHED_DIR, pattern = "site_.*_enriched\\.csv(\\.gz)?$",
                  full.names = TRUE)
  if (!length(f)) stop("No per-site enriched CSVs under: ", ENRICHED_DIR)
  get_site <- function(x) {
    b <- basename(x)
    for (pat in pats) {
      m <- regexec(pat, b); r <- regmatches(b, m)[[1]]
      if (length(r) == 2) return(as.integer(r[2]))
    }
    NA_integer_
  }
  sid <- vapply(f, get_site, integer(1)); keep <- !is.na(sid)
  data.frame(site = sid[keep], path = f[keep], stringsAsFactors = FALSE)
})()

# =====================================================================
# FIXED: Product-specific site exclusion
# =====================================================================
# Parse product-specific exclusions from config
excl_chm <- if (exists("SITES_EXCLUDE_CHM")) .parse_sites(SITES_EXCLUDE_CHM) else integer()
excl_dtm <- if (exists("SITES_EXCLUDE_DTM")) .parse_sites(SITES_EXCLUDE_DTM) else integer()

# Legacy global exclusion (backward compatibility)
incl <- .parse_sites(SITES_INCLUDE)
excl_global <- .parse_sites(SITES_EXCLUDE)

# Apply global include filter to file table (affects both products)
if (length(incl)) site_file_tbl <- subset(site_file_tbl, site %in% incl)
# Note: Global exclude is now deprecated - product-specific exclusions applied in read function

site_file_tbl <- site_file_tbl[order(site_file_tbl$site), , drop = FALSE]
stopifnot(nrow(site_file_tbl) > 0)

log_progress(sprintf("Found %d sites to process", nrow(site_file_tbl)))
log_progress(sprintf("CHM exclusions: %s", 
                     ifelse(length(excl_chm), paste(excl_chm, collapse=","), "(none)")))
log_progress(sprintf("DTM exclusions: %s", 
                     ifelse(length(excl_dtm), paste(excl_dtm, collapse=","), "(none)")))

# Per-site read function
.read_site_politely <- function(site_path, site_id) {
  hdr <- names(data.table::fread(site_path, nrows = 0))
  use <- intersect(sel_cols, hdr)
  cc  <- list(shot_number = "character")
  DT  <- data.table::fread(site_path, select = use,
                           colClasses = cc, nThread = min(4L, CPU_BUDGET),
                           showProgress = FALSE)
  if (!nrow(DT)) return(list(chm = data.table::data.table(), dtm = data.table::data.table()))
  
  need_na <- setdiff(sel_cols, names(DT))
  if (length(need_na)) for (nm in need_na) DT[, (nm) := NA_real_]
  
  # Robust LC column names
  if (!("lc2022_l1_code" %in% names(DT)) && "lc2022_mode_l1_code" %in% names(DT))
    setnames(DT, "lc2022_mode_l1_code", "lc2022_l1_code")
  if (!("lc2022_l1_name" %in% names(DT)) && "lc2022_mode_l1_name" %in% names(DT))
    setnames(DT, "lc2022_mode_l1_name", "lc2022_l1_name")
  
  # Vendor spelling fix
  if ("meta_veh_ct_GEO1" %in% names(DT)) setnames(DT, "meta_veh_ct_GEO1", "meta_veh_ct_GE01")
  if ("meta_veh_ratio_GEO1" %in% names(DT)) setnames(DT, "meta_veh_ratio_GEO1", "meta_veh_ratio_GE01")
  
  DT[, `:=`(
    lc_l1_code = factor(fifelse(is.na(lc2022_l1_code) | lc2022_l1_code == "", "UNK", lc2022_l1_code)),
    lc_l1_name = factor(fifelse(is.na(lc2022_l1_name) | lc2022_l1_name == "", "Unknown", lc2022_l1_name)),
    site       = factor(site),
    ecoregion  = factor(ecoregion)
  )]
  
  # =====================================================================
  # CHM subset (FIXED: product-specific exclusion)
  # =====================================================================
  # Check if this site should be excluded from CHM analysis
  if (length(excl_chm) && site_id %in% excl_chm) {
    log_progress(sprintf("  Site %d excluded from CHM analysis", site_id))
    CHM <- data.table::data.table()
  } else {
    CHM <- DT[is.finite(p3d_chm_mean) & is.finite(als_chm_mean) &
                is.finite(als_chm_valid_frac) & is.finite(p3d_chm_valid_frac) &
                als_chm_valid_frac >= 0.5 & p3d_chm_valid_frac >= 0.5 &
                is.finite(slope_valid_frac) & slope_valid_frac >= 0.5]
    if (apply_chm_forest_filter)
      CHM <- CHM[is.finite(als_chm_p90) & als_chm_p90 >= chm_forest_thresh_m]
    
    if (nrow(CHM)) {
      CHM[, chm_error_mean := fifelse(is.finite(error_mean),
                                      error_mean, p3d_chm_mean - als_chm_mean)]
      CHM <- CHM[abs(chm_error_mean) <= 100 | !is.finite(chm_error_mean)]
      # FIXED: Added z-scaled azimuth variables
      CHM[, `:=`(
        view_az_sin     = sin(pi * meta_target_azimuth_avg / 180),
        view_az_cos     = cos(pi * meta_target_azimuth_avg / 180),
        view_az_sin_z   = zscale(sin(pi * meta_target_azimuth_avg / 180)),
        view_az_cos_z   = zscale(cos(pi * meta_target_azimuth_avg / 180)),
        meta_offnad_z   = zscale(meta_off_nadir_avg),
        meta_sunel_z    = zscale(meta_sun_elev_avg),
        meta_az_conc_z  = zscale(meta_az_concentration),
        meta_absgeo_z   = zscale(meta_abs_geoacc_avg),
        meta_relgeo_z   = zscale(meta_rel_geoacc_avg),
        meta_tot_z      = zscale(meta_tot_ct),
        meta_stereo_z   = zscale(meta_stereo_ratio),
        meta_fwdrev_z   = zscale(meta_fwd_ratio - meta_rev_ratio),
        meta_leafon_z   = zscale(meta_leaf_on_ratio),
        slope_mean_z    = zscale(slope_mean),
        slope_sd_z      = zscale(slope_sd),
        wsci_z          = zscale(wsci),
        rh_98_z         = zscale(rh_98),
        cover_z         = zscale(cover),
        aspect_sin_z    = zscale(aspect_sin_mean),
        aspect_cos_z    = zscale(aspect_cos_mean),
        w_chm           = pmin(als_chm_valid_frac, p3d_chm_valid_frac, na.rm = TRUE),
        v_GE01 = meta_veh_ratio_GE01,
        v_WV01 = meta_veh_ratio_WV01,
        v_WV02 = meta_veh_ratio_WV02,
        v_WV03 = meta_veh_ratio_WV03
      )]
    }
  }
  
  # =====================================================================
  # DTM subset (FIXED: product-specific exclusion)
  # =====================================================================
  # Check if this site should be excluded from DTM analysis
  if (length(excl_dtm) && site_id %in% excl_dtm) {
    log_progress(sprintf("  Site %d excluded from DTM analysis", site_id))
    DTM <- data.table::data.table()
  } else {
    DTM <- DT[is.finite(p3d_dtm_mean) & is.finite(dep_dtm_mean) &
                is.finite(dep_dtm_valid_frac) & is.finite(p3d_dtm_valid_frac) &
                dep_dtm_valid_frac >= 0.5 & p3d_dtm_valid_frac >= 0.5 &
                is.finite(slope_valid_frac) & slope_valid_frac >= 0.5]
    if (apply_dtm_forest_filter && dtm_forest_proxy == "als_chm_p90")
      DTM <- DTM[is.finite(als_chm_p90) & als_chm_p90 >= chm_forest_thresh_m]
    if (apply_dtm_forest_filter && dtm_forest_proxy == "cover")
      DTM <- DTM[is.finite(cover) & cover >= dtm_cover_threshold]
    
    if (nrow(DTM)) {
      DTM[, dtm_error_mean := fifelse(is.finite(dtm_error_mean),
                                      dtm_error_mean, p3d_dtm_mean - dep_dtm_mean)]
      DTM <- DTM[abs(dtm_error_mean) <= 100 | !is.finite(dtm_error_mean)]
      # FIXED: Added z-scaled azimuth variables
      DTM[, `:=`(
        view_az_sin     = sin(pi * meta_target_azimuth_avg / 180),
        view_az_cos     = cos(pi * meta_target_azimuth_avg / 180),
        view_az_sin_z   = zscale(sin(pi * meta_target_azimuth_avg / 180)),
        view_az_cos_z   = zscale(cos(pi * meta_target_azimuth_avg / 180)),
        meta_offnad_z   = zscale(meta_off_nadir_avg),
        meta_sunel_z    = zscale(meta_sun_elev_avg),
        meta_az_conc_z  = zscale(meta_az_concentration),
        meta_absgeo_z   = zscale(meta_abs_geoacc_avg),
        meta_relgeo_z   = zscale(meta_rel_geoacc_avg),
        meta_tot_z      = zscale(meta_tot_ct),
        meta_stereo_z   = zscale(meta_stereo_ratio),
        meta_fwdrev_z   = zscale(meta_fwd_ratio - meta_rev_ratio),
        meta_leafon_z   = zscale(meta_leaf_on_ratio),
        slope_mean_z    = zscale(slope_mean),
        slope_sd_z      = zscale(slope_sd),
        wsci_z          = zscale(wsci),
        rh_98_z         = zscale(rh_98),
        cover_z         = zscale(cover),
        aspect_sin_z    = zscale(aspect_sin_mean),
        aspect_cos_z    = zscale(aspect_cos_mean),
        w_dtm           = pmin(dep_dtm_valid_frac, p3d_dtm_valid_frac, na.rm = TRUE),
        v_GE01 = meta_veh_ratio_GE01,
        v_WV01 = meta_veh_ratio_WV01,
        v_WV02 = meta_veh_ratio_WV02,
        v_WV03 = meta_veh_ratio_WV03
      )]
    }
  }
  
  # FIXED: Added z-scaled azimuth to keep lists
  chm_keep <- c(
    "site","ecoregion","shot_number","lc_l1_code","lc_l1_name",
    "slope_mean","slope_mean_z","slope_sd_z","wsci_z","rh_98_z","cover_z",
    "aspect_sin_z","aspect_cos_z",
    "meta_stereo_any","meta_off_nadir_avg","meta_sun_elev_avg",
    "meta_target_azimuth_avg",
    "view_az_sin","view_az_cos","view_az_sin_z","view_az_cos_z",  # Added z-scaled
    "meta_offnad_z","meta_sunel_z","meta_az_conc_z",
    "meta_absgeo_z","meta_relgeo_z","meta_tot_z","meta_stereo_z",
    "meta_fwdrev_z","meta_leafon_z",
    "v_GE01","v_WV01","v_WV02","v_WV03",
    "chm_error_mean","w_chm","als_chm_p90",
    "als_chm_valid_frac","p3d_chm_valid_frac"
  )
  
  dtm_keep <- c(
    "site","ecoregion","shot_number","lc_l1_code","lc_l1_name",
    "slope_mean","slope_mean_z","slope_sd_z","wsci_z","rh_98_z","cover_z",
    "aspect_sin_z","aspect_cos_z",
    "meta_stereo_any","meta_off_nadir_avg","meta_sun_elev_avg",
    "meta_target_azimuth_avg",
    "view_az_sin","view_az_cos","view_az_sin_z","view_az_cos_z",  # Added z-scaled
    "meta_offnad_z","meta_sunel_z","meta_az_conc_z",
    "meta_absgeo_z","meta_relgeo_z","meta_tot_z","meta_stereo_z",
    "meta_fwdrev_z","meta_leafon_z",
    "v_GE01","v_WV01","v_WV02","v_WV03",
    "dtm_error_mean","w_dtm",
    "dep_dtm_valid_frac","p3d_dtm_valid_frac"
  )
  
  CHM <- if (nrow(CHM)) CHM[, intersect(chm_keep, names(CHM)), with = FALSE] else CHM
  DTM <- if (nrow(DTM)) DTM[, intersect(dtm_keep, names(DTM)), with = FALSE] else DTM
  list(chm = CHM, dtm = DTM)
}

# Main per-site loop
log_progress(sprintf("Reading %d sites...", nrow(site_file_tbl)))
pb <- progress_bar$new(format = "  [:bar] :current/:total (:percent) eta: :eta",
                       total = nrow(site_file_tbl), clear = FALSE, width = 68)
chm_parts <- vector("list", nrow(site_file_tbl))
dtm_parts <- vector("list", nrow(site_file_tbl))

for (i in seq_len(nrow(site_file_tbl))) {
  s  <- site_file_tbl$site[i]
  fp <- site_file_tbl$path[i]
  res <- .read_site_politely(fp, s)
  chm_parts[[i]] <- res$chm
  dtm_parts[[i]] <- res$dtm
  pb$tick()
  Sys.sleep(0.05)
}

chm_df <- data.table::rbindlist(chm_parts, use.names = TRUE, fill = TRUE)
dtm_df <- data.table::rbindlist(dtm_parts, use.names = TRUE, fill = TRUE)

log_progress(sprintf("CHM: %s rows | DTM: %s rows", 
                     format(nrow(chm_df), big.mark=","),
                     format(nrow(dtm_df), big.mark=",")))
log_progress(sprintf("RAM (CHM) ≈ %s | RAM (DTM) ≈ %s", 
                     format(object.size(chm_df), units = "auto"),
                     format(object.size(dtm_df), units = "auto")))
gc()

# Attach coordinates if requested
if (attach_coords) {
  log_progress("Attaching coordinates via GEDI gpkg...")
  add_coords <- function(df) {
    stopifnot("shot_number" %in% names(df))
    sites <- sort(unique(df$site))
    out   <- df
    out$shot_key <- as.character(out$shot_number)
    if (!("x" %in% names(out))) out$x <- NA_real_
    if (!("y" %in% names(out))) out$y <- NA_real_
    for (s in sites) {
      gpkg <- file.path(gedi_base_dir, as.character(s),
                        sprintf(gedi_file_name, as.character(s)))
      if (!file.exists(gpkg)) next
      g <- tryCatch(
        suppressMessages(sf::st_read(gpkg, quiet = TRUE, int64_as_string = TRUE)),
        error = function(e) suppressMessages(sf::st_read(gpkg, quiet = TRUE))
      )
      g <- sf::st_transform(g, coord_crs_out)
      g_xy <- g %>%
        dplyr::select(shot_number) %>%
        dplyr::mutate(
          shot_key = if (is.character(shot_number)) shot_number else
            format(shot_number, scientific = FALSE, trim = TRUE),
          x = sf::st_coordinates(.)[,1],
          y = sf::st_coordinates(.)[,2]
        ) %>%
        sf::st_drop_geometry() %>%
        dplyr::filter(!is.na(shot_key)) %>%
        dplyr::distinct(shot_key, .keep_all = TRUE)
      idx  <- which(out$site == s)
      sk   <- out$shot_key[idx]
      pos  <- match(sk, g_xy$shot_key)
      out$x[idx] <- g_xy$x[pos]; out$y[idx] <- g_xy$y[pos]
      n_miss <- sum(is.na(out$x[idx]) | is.na(out$y[idx]))
      log_progress(sprintf("  Site %s: %d/%d rows (missed %d)",
                           s, length(idx) - n_miss, length(idx), n_miss))
    }
    out$shot_key <- NULL
    out
  }
  chm_df <- add_coords(chm_df)
  dtm_df <- add_coords(dtm_df)
}

# Convert to data.frame for modeling
chm_df <- as.data.frame(chm_df)
dtm_df <- as.data.frame(dtm_df)

# Remove duplicate columns
chm_df <- chm_df[, !duplicated(names(chm_df))]
dtm_df <- dtm_df[, !duplicated(names(dtm_df))]

# Verify metadata availability (FIXED: includes z-scaled azimuth)
available_meta_chm <- intersect(
  c("meta_offnad_z", "meta_sunel_z", "meta_az_conc_z", 
    "meta_absgeo_z", "meta_relgeo_z", "meta_tot_z", 
    "meta_stereo_z", "meta_fwdrev_z", "meta_leafon_z",
    "view_az_sin_z", "view_az_cos_z"),  # Added z-scaled azimuth
  names(chm_df)
)

available_meta_dtm <- intersect(
  c("meta_offnad_z", "meta_sunel_z", "meta_az_conc_z", 
    "meta_absgeo_z", "meta_relgeo_z", "meta_tot_z", 
    "meta_stereo_z", "meta_fwdrev_z", "meta_leafon_z",
    "view_az_sin_z", "view_az_cos_z"),  # Added z-scaled azimuth
  names(dtm_df)
)

log_progress(sprintf("CHM metadata available: %s", 
                     paste(available_meta_chm, collapse=", ")))
log_progress(sprintf("DTM metadata available: %s", 
                     paste(available_meta_dtm, collapse=", ")))

log_progress("✓ Data ingest complete")

# Save checkpoint for interactive mode (if not in batch mode)
if (!exists("BATCH_MODE") || !BATCH_MODE) {
  log_progress("Saving checkpoint for interactive mode...")
  save_checkpoint("01_data_ingest", list(
    chm_df = chm_df,
    dtm_df = dtm_df,
    available_meta_chm = available_meta_chm,
    available_meta_dtm = available_meta_dtm,
    chm_rows = nrow(chm_df),
    dtm_rows = nrow(dtm_df),
    sites_excluded_chm = excl_chm,
    sites_excluded_dtm = excl_dtm
  ))
}
