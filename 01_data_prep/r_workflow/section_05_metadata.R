# =====================================================================
# section_05_metadata.R (ENHANCED VERSION)
# MHRSI metadata relationship analysis
#
# ENHANCEMENTS ADDED:
#   1. Forest-focused metadata analysis (EBF, BDF, ENF, DNF)
#   2. Leaf-on × forest type interaction (key for phenology hypothesis)
#   3. Deciduous vs Evergreen metadata comparisons
#   4. Forest-specific metadata summary statistics
# =====================================================================

# Ensure dependencies loaded
if (!exists("log_progress")) source("analysis_utils.R")
if (!exists("PROJECT_ROOT")) source("analysis_config.R")

if (!exists("chm_df") || !exists("dtm_df")) {
  log_progress("⚠ Data not loaded. Loading from checkpoint...")
  if (checkpoint_exists("01_data_ingest")) {
    data <- load_checkpoint("01_data_ingest")
    chm_df <- data$chm_df
    dtm_df <- data$dtm_df
    available_meta_chm <- data$available_meta_chm
    available_meta_dtm <- data$available_meta_dtm
  } else {
    stop("Data not available. Run: source('section_01_ingest.R') first")
  }
}

if (!exists("n_bins_common")) {
  n_chm_slope <- sum(is.finite(chm_df$slope_mean) & is.finite(chm_df$chm_error_mean))
  n_dtm_slope <- sum(is.finite(dtm_df$slope_mean) & is.finite(dtm_df$dtm_error_mean))
  n_bins_common <- max(choose_hex_bins(n_chm_slope), choose_hex_bins(n_dtm_slope))
}

log_progress("Analyzing MHRSI metadata relationships (ENHANCED)...")

# =====================================================================
# FOREST CLASS DEFINITIONS
# =====================================================================

FOREST_CLASSES <- c("EBF", "BDF", "ENF", "DNF")
FOREST_LABELS <- c("Broadleaf Evergreen", "Broadleaf Deciduous", 
                   "Needleleaf Evergreen", "Needleleaf Deciduous")
names(FOREST_LABELS) <- FOREST_CLASSES

# Create forest subsets with phenology classification
chm_forest <- chm_df %>% 
  filter(lc_l1_code %in% FOREST_CLASSES) %>%
  mutate(
    phenology = case_when(
      lc_l1_code %in% c("EBF", "ENF") ~ "Evergreen",
      lc_l1_code %in% c("BDF", "DNF") ~ "Deciduous",
      TRUE ~ NA_character_
    ),
    leaf_structure = case_when(
      lc_l1_code %in% c("EBF", "BDF") ~ "Broadleaf",
      lc_l1_code %in% c("ENF", "DNF") ~ "Needleleaf",
      TRUE ~ NA_character_
    )
  )

dtm_forest <- dtm_df %>% 
  filter(lc_l1_code %in% FOREST_CLASSES) %>%
  mutate(
    phenology = case_when(
      lc_l1_code %in% c("EBF", "ENF") ~ "Evergreen",
      lc_l1_code %in% c("BDF", "DNF") ~ "Deciduous",
      TRUE ~ NA_character_
    ),
    leaf_structure = case_when(
      lc_l1_code %in% c("EBF", "BDF") ~ "Broadleaf",
      lc_l1_code %in% c("ENF", "DNF") ~ "Needleleaf",
      TRUE ~ NA_character_
    )
  )

forest_present_chm <- intersect(FOREST_CLASSES, unique(chm_df$lc_l1_code))
log_progress(sprintf("  Forest classes present: %s", paste(forest_present_chm, collapse = ", ")))
log_progress(sprintf("  Forest observations: CHM=%s (%.1f%%)", 
                     format(nrow(chm_forest), big.mark=","),
                     100 * nrow(chm_forest) / nrow(chm_df)))

# Helper: check if column has enough data to plot
can_plot_meta <- function(df, col, err_col, w_col, min_n = 1000) {
  if (!col %in% names(df)) return(FALSE)
  # Check for sufficient non-NA combinations of x, y, and weight
  valid <- is.finite(df[[col]]) & is.finite(df[[err_col]]) & is.finite(df[[w_col]])
  sum(valid) >= min_n
}

# Helper: safe plot wrapper
try_plot <- function(expr, name) {
  tryCatch(expr, error = function(e) {
    log_progress(sprintf("  ⚠ Skipping %s: insufficient data after filtering", name))
    NULL
  })
}

# =====================================================================
# 1. VIEWING GEOMETRY (Original)
# =====================================================================

log_subsection("Viewing geometry (off-nadir, sun elevation, azimuth)")
plot_list_chm <- list()
plot_list_dtm <- list()

if (can_plot_meta(chm_df, "meta_offnad_z", "chm_error_mean", "w_chm")) {
  p <- try_plot(
    make_hex_trend_plot(chm_df, "meta_offnad_z","chm_error_mean","w_chm",
                        "Off-nadir (z)","CHM error (m)",
                        "CHM error vs off-nadir — hex + trend",
                        n_bins_x = n_bins_common, n_bins_y = n_bins_common),
    "CHM off-nadir plot")
  if (!is.null(p)) plot_list_chm$offnad <- p
}

if (can_plot_meta(chm_df, "meta_sunel_z", "chm_error_mean", "w_chm")) {
  p <- try_plot(
    make_hex_trend_plot(chm_df, "meta_sunel_z","chm_error_mean","w_chm",
                        "Sun elevation (z)","CHM error (m)",
                        "CHM error vs sun elevation — hex + trend",
                        n_bins_x = n_bins_common, n_bins_y = n_bins_common),
    "CHM sun elevation plot")
  if (!is.null(p)) plot_list_chm$sunel <- p
}

if (can_plot_meta(chm_df, "meta_az_conc_z", "chm_error_mean", "w_chm")) {
  p <- try_plot(
    make_hex_trend_plot(chm_df, "meta_az_conc_z","chm_error_mean","w_chm",
                        "Azimuth concentration (z)","CHM error (m)",
                        "CHM error vs azimuth concentration — hex + trend",
                        n_bins_x = n_bins_common, n_bins_y = n_bins_common),
    "CHM azimuth concentration plot")
  if (!is.null(p)) plot_list_chm$azconc <- p
}

if (can_plot_meta(dtm_df, "meta_offnad_z", "dtm_error_mean", "w_dtm")) {
  p <- try_plot(
    make_hex_trend_plot(dtm_df, "meta_offnad_z","dtm_error_mean","w_dtm",
                        "Off-nadir (z)","DTM error (m)",
                        "DTM error vs off-nadir — hex + trend",
                        n_bins_x = n_bins_common, n_bins_y = n_bins_common),
    "DTM off-nadir plot")
  if (!is.null(p)) plot_list_dtm$offnad <- p
}

if (can_plot_meta(dtm_df, "meta_sunel_z", "dtm_error_mean", "w_dtm")) {
  p <- try_plot(
    make_hex_trend_plot(dtm_df, "meta_sunel_z","dtm_error_mean","w_dtm",
                        "Sun elevation (z)","DTM error (m)",
                        "DTM error vs sun elevation — hex + trend",
                        n_bins_x = n_bins_common, n_bins_y = n_bins_common),
    "DTM sun elevation plot")
  if (!is.null(p)) plot_list_dtm$sunel <- p
}

if (can_plot_meta(dtm_df, "meta_az_conc_z", "dtm_error_mean", "w_dtm")) {
  p <- try_plot(
    make_hex_trend_plot(dtm_df, "meta_az_conc_z","dtm_error_mean","w_dtm",
                        "Azimuth concentration (z)","DTM error (m)",
                        "DTM error vs azimuth concentration — hex + trend",
                        n_bins_x = n_bins_common, n_bins_y = n_bins_common),
    "DTM azimuth concentration plot")
  if (!is.null(p)) plot_list_dtm$azconc <- p
}

# Save only if we have plots
if (length(plot_list_chm) > 0) {
  combined_chm <- patchwork::wrap_plots(plot_list_chm, ncol = length(plot_list_chm))
  ggsave(file.path(out_plots, "02b_chm_meta_panels.pdf"),
         combined_chm, width=14, height=4.5, bg="white")
}

if (length(plot_list_dtm) > 0) {
  combined_dtm <- patchwork::wrap_plots(plot_list_dtm, ncol = length(plot_list_dtm))
  ggsave(file.path(out_plots, "02b_dtm_meta_panels.pdf"),
         combined_dtm, width=14, height=4.5, bg="white")
}

# =====================================================================
# 1b. NEW: VIEWING GEOMETRY (FOREST ONLY)
# =====================================================================

log_subsection("Viewing geometry (FOREST ONLY)")
plot_list_chm_forest <- list()

if (can_plot_meta(chm_forest, "meta_offnad_z", "chm_error_mean", "w_chm", min_n = 500)) {
  p <- try_plot(
    make_hex_trend_plot(chm_forest, "meta_offnad_z","chm_error_mean","w_chm",
                        "Off-nadir (z)","CHM error (m)",
                        "CHM error vs off-nadir — Forest Types Only",
                        n_bins_x = n_bins_common, n_bins_y = n_bins_common),
    "CHM forest off-nadir plot")
  if (!is.null(p)) plot_list_chm_forest$offnad <- p
}

if (can_plot_meta(chm_forest, "meta_sunel_z", "chm_error_mean", "w_chm", min_n = 500)) {
  p <- try_plot(
    make_hex_trend_plot(chm_forest, "meta_sunel_z","chm_error_mean","w_chm",
                        "Sun elevation (z)","CHM error (m)",
                        "CHM error vs sun elevation — Forest Types Only",
                        n_bins_x = n_bins_common, n_bins_y = n_bins_common),
    "CHM forest sun elevation plot")
  if (!is.null(p)) plot_list_chm_forest$sunel <- p
}

if (length(plot_list_chm_forest) > 0) {
  combined_forest <- patchwork::wrap_plots(plot_list_chm_forest, ncol = length(plot_list_chm_forest))
  ggsave(file.path(out_plots, "02b_chm_meta_panels_FOREST.pdf"),
         combined_forest, width=12, height=4.5, bg="white")
}

# =====================================================================
# 2. GEOLOCATION ACCURACY (Original)
# =====================================================================

log_subsection("Geolocation accuracy")
plot_list_geo <- list()

if (can_plot_meta(chm_df, "meta_absgeo_z", "chm_error_mean", "w_chm")) {
  p <- try_plot(
    make_hex_trend_plot(chm_df, "meta_absgeo_z","chm_error_mean","w_chm",
                        "Abs geoacc (z)","CHM error (m)",
                        "CHM error vs absolute geolocation accuracy (z)"),
    "CHM absolute geolocation accuracy plot")
  if (!is.null(p)) plot_list_geo$absgeo <- p
}

if (can_plot_meta(chm_df, "meta_relgeo_z", "chm_error_mean", "w_chm")) {
  p <- try_plot(
    make_hex_trend_plot(chm_df, "meta_relgeo_z","chm_error_mean","w_chm",
                        "Rel geoacc (z)","CHM error (m)",
                        "CHM error vs relative geolocation accuracy (z)"),
    "CHM relative geolocation accuracy plot")
  if (!is.null(p)) plot_list_geo$relgeo <- p
}

if (length(plot_list_geo) > 0) {
  combined_geo <- patchwork::wrap_plots(plot_list_geo, ncol = length(plot_list_geo))
  ggsave(file.path(out_plots, "02c_chm_hex_geoacc.pdf"),
         combined_geo, width=12, height=5, bg="white")
}

# =====================================================================
# 3. STEREO, FORWARD/REVERSE, LEAF-ON (Original)
# =====================================================================

log_subsection("Stereo, forward/reverse, seasonality")
plot_list_other <- list()

if (can_plot_meta(chm_df, "meta_stereo_z", "chm_error_mean", "w_chm")) {
  p <- try_plot(
    make_hex_trend_plot(chm_df, "meta_stereo_z","chm_error_mean","w_chm",
                        "Stereo ratio (z)","CHM error (m)",
                        "CHM error vs stereo ratio (z)"),
    "CHM stereo ratio plot")
  if (!is.null(p)) plot_list_other$stereo <- p
}

if (can_plot_meta(chm_df, "meta_fwdrev_z", "chm_error_mean", "w_chm")) {
  p <- try_plot(
    make_hex_trend_plot(chm_df, "meta_fwdrev_z","chm_error_mean","w_chm",
                        "Fwd−Rev ratio (z)","CHM error (m)",
                        "CHM error vs forward−reverse (z)"),
    "CHM forward-reverse ratio plot")
  if (!is.null(p)) plot_list_other$fwdrev <- p
}

if (can_plot_meta(chm_df, "meta_leafon_z", "chm_error_mean", "w_chm")) {
  p <- try_plot(
    make_hex_trend_plot(chm_df, "meta_leafon_z","chm_error_mean","w_chm",
                        "Leaf-on ratio (z)","CHM error (m)",
                        "CHM error vs leaf-on ratio (z)"),
    "CHM leaf-on ratio plot")
  if (!is.null(p)) plot_list_other$leafon <- p
}

if (length(plot_list_other) > 0) {
  combined_other <- patchwork::wrap_plots(plot_list_other, ncol = length(plot_list_other))
  ggsave(file.path(out_plots, "02d_chm_hex_stereo_fwdrev_leafon.pdf"),
         combined_other, width=16, height=5, bg="white")
}

# =====================================================================
# 3b. NEW: LEAF-ON × FOREST TYPE INTERACTION (Key for phenology hypothesis)
# =====================================================================

log_subsection("Leaf-on × Forest type interaction (NEW - Key analysis)")

if ("meta_leafon_z" %in% names(chm_forest) && 
    sum(is.finite(chm_forest$meta_leafon_z)) > 1000) {
  
  # Hex plots faceted by forest type
  tryCatch({
    p_leafon_forest_hex <- ggplot(chm_forest %>% filter(is.finite(meta_leafon_z)),
                                   aes(x = meta_leafon_z, y = chm_error_mean)) +
      geom_hex(bins = 30) +
      geom_smooth(method = "gam", color = "red", linewidth = 1) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "white") +
      facet_wrap(~lc_l1_code, ncol = 2) +
      scale_fill_viridis_c(trans = "log10", name = "Count") +
      labs(title = "CHM Error vs Leaf-on Ratio by Forest Type",
           subtitle = "Testing phenology hypothesis: deciduous forests should show stronger leaf-on effect",
           x = "Leaf-on ratio (z-scaled)", y = "CHM error (m)") +
      theme_cowplot()
    
    ggsave(file.path(out_plots, "02d_leafon_by_forest_type.pdf"), p_leafon_forest_hex,
           width = 10, height = 8, bg = "white")
    log_progress("  ✓ Leaf-on × forest type hex plot saved")
  }, error = function(e) {
    log_progress(sprintf("  ⚠ Leaf-on × forest hex plot failed: %s", e$message))
  })
  
  # Deciduous vs Evergreen comparison
  tryCatch({
    p_leafon_phenology <- ggplot(chm_forest %>% filter(is.finite(meta_leafon_z), !is.na(phenology)),
                                  aes(x = meta_leafon_z, y = chm_error_mean, color = phenology)) +
      geom_smooth(method = "gam", linewidth = 1.2) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      scale_color_manual(values = c("Deciduous" = "#E69F00", "Evergreen" = "#009E73")) +
      labs(title = "Leaf-on Effect: Deciduous vs Evergreen Forests",
           subtitle = "Deciduous forests expected to show stronger relationship with leaf-on ratio",
           x = "Leaf-on ratio (z-scaled)", y = "CHM error (m)",
           color = "Phenology") +
      theme_cowplot()
    
    ggsave(file.path(out_plots, "02d_leafon_deciduous_vs_evergreen.pdf"), p_leafon_phenology,
           width = 9, height = 6, bg = "white")
    log_progress("  ✓ Deciduous vs Evergreen leaf-on comparison saved")
  }, error = function(e) {
    log_progress(sprintf("  ⚠ Deciduous vs Evergreen plot failed: %s", e$message))
  })
  
  # Interaction statistics: binned mean error by leaf-on × phenology
  tryCatch({
    leafon_interaction <- chm_forest %>%
      filter(is.finite(meta_leafon_z), !is.na(phenology)) %>%
      mutate(leafon_bin = cut(meta_leafon_z, breaks = seq(-3, 3, 1), include.lowest = TRUE)) %>%
      filter(!is.na(leafon_bin)) %>%
      group_by(phenology, leafon_bin) %>%
      summarise(
        n = n(),
        mean_error = mean(chm_error_mean, na.rm = TRUE),
        sd_error = sd(chm_error_mean, na.rm = TRUE),
        se_error = sd_error / sqrt(n),
        .groups = "drop"
      ) %>%
      filter(n >= 30)
    
    p_leafon_interact <- ggplot(leafon_interaction, 
                                 aes(x = leafon_bin, y = mean_error, 
                                     color = phenology, group = phenology)) +
      geom_line(linewidth = 1) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = mean_error - 1.96*se_error,
                        ymax = mean_error + 1.96*se_error), width = 0.2) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      scale_color_manual(values = c("Deciduous" = "#E69F00", "Evergreen" = "#009E73")) +
      labs(title = "Leaf-on × Phenology Interaction",
           subtitle = "Mean CHM error in leaf-on bins with 95% CI",
           x = "Leaf-on ratio bin (z-scaled)", y = "Mean CHM error (m)",
           color = "Phenology") +
      theme_cowplot() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(file.path(out_plots, "02d_leafon_phenology_interaction.pdf"), p_leafon_interact,
           width = 9, height = 6, bg = "white")
    
    write_csv(leafon_interaction, file.path(out_tables, "leafon_phenology_interaction.csv"))
    log_progress("  ✓ Leaf-on × phenology interaction saved")
    
    # Check for differential effect
    decid_range <- leafon_interaction %>% 
      filter(phenology == "Deciduous") %>%
      summarise(range = max(mean_error) - min(mean_error)) %>%
      pull(range)
    
    everg_range <- leafon_interaction %>% 
      filter(phenology == "Evergreen") %>%
      summarise(range = max(mean_error) - min(mean_error)) %>%
      pull(range)
    
    log_progress(sprintf("  Leaf-on effect range: Deciduous=%.2fm, Evergreen=%.2fm", 
                         decid_range, everg_range))
    
    if (decid_range > everg_range * 1.5) {
      log_progress("  ✓ Phenology hypothesis SUPPORTED: Deciduous shows stronger leaf-on effect")
    } else if (everg_range > decid_range * 1.5) {
      log_progress("  ⚠ Unexpected: Evergreen shows stronger leaf-on effect")
    } else {
      log_progress("  ~ Similar leaf-on effect for both phenology types")
    }
    
  }, error = function(e) {
    log_progress(sprintf("  ⚠ Leaf-on interaction statistics failed: %s", e$message))
  })
  
  # Individual forest type leaf-on interaction
  tryCatch({
    leafon_by_lc <- chm_forest %>%
      filter(is.finite(meta_leafon_z)) %>%
      mutate(leafon_bin = cut(meta_leafon_z, breaks = seq(-3, 3, 1), include.lowest = TRUE)) %>%
      filter(!is.na(leafon_bin)) %>%
      group_by(lc_l1_code, leafon_bin) %>%
      summarise(
        n = n(),
        mean_error = mean(chm_error_mean, na.rm = TRUE),
        se_error = sd(chm_error_mean, na.rm = TRUE) / sqrt(n()),
        .groups = "drop"
      ) %>%
      filter(n >= 20)
    
    p_leafon_lc <- ggplot(leafon_by_lc, 
                           aes(x = leafon_bin, y = mean_error, 
                               color = lc_l1_code, group = lc_l1_code)) +
      geom_line(linewidth = 1) +
      geom_point(size = 2) +
      geom_errorbar(aes(ymin = mean_error - 1.96*se_error,
                        ymax = mean_error + 1.96*se_error), width = 0.2) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      scale_color_brewer(palette = "Dark2") +
      labs(title = "Leaf-on Effect by Forest Type",
           x = "Leaf-on ratio bin (z-scaled)", y = "Mean CHM error (m)",
           color = "Forest Type") +
      theme_cowplot() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(file.path(out_plots, "02d_leafon_by_lc_interaction.pdf"), p_leafon_lc,
           width = 10, height = 6, bg = "white")
    
    write_csv(leafon_by_lc, file.path(out_tables, "leafon_by_forest_type.csv"))
    
  }, error = function(e) {
    log_progress(sprintf("  ⚠ Leaf-on by forest type failed: %s", e$message))
  })
  
} else {
  log_progress("  ⚠ Leaf-on × forest analysis: SKIPPED (insufficient data)")
}

# =====================================================================
# 3c. NEW: OFF-NADIR × FOREST TYPE INTERACTION
# =====================================================================

log_subsection("Off-nadir × Forest type interaction (NEW)")

if ("meta_offnad_z" %in% names(chm_forest) && 
    sum(is.finite(chm_forest$meta_offnad_z)) > 1000) {
  
  tryCatch({
    offnad_by_phenology <- chm_forest %>%
      filter(is.finite(meta_offnad_z), !is.na(phenology)) %>%
      mutate(offnad_bin = cut(meta_offnad_z, breaks = seq(-2, 3, 1), include.lowest = TRUE)) %>%
      filter(!is.na(offnad_bin)) %>%
      group_by(phenology, offnad_bin) %>%
      summarise(
        n = n(),
        mean_error = mean(chm_error_mean, na.rm = TRUE),
        se_error = sd(chm_error_mean, na.rm = TRUE) / sqrt(n()),
        .groups = "drop"
      ) %>%
      filter(n >= 30)
    
    p_offnad_phenology <- ggplot(offnad_by_phenology, 
                                  aes(x = offnad_bin, y = mean_error, 
                                      color = phenology, group = phenology)) +
      geom_line(linewidth = 1) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = mean_error - 1.96*se_error,
                        ymax = mean_error + 1.96*se_error), width = 0.2) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      scale_color_manual(values = c("Deciduous" = "#E69F00", "Evergreen" = "#009E73")) +
      labs(title = "Off-nadir Effect: Deciduous vs Evergreen Forests",
           x = "Off-nadir bin (z-scaled)", y = "Mean CHM error (m)",
           color = "Phenology") +
      theme_cowplot() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(file.path(out_plots, "02b_offnadir_phenology_interaction.pdf"), p_offnad_phenology,
           width = 9, height = 6, bg = "white")
    
    write_csv(offnad_by_phenology, file.path(out_tables, "offnadir_phenology_interaction.csv"))
    log_progress("  ✓ Off-nadir × phenology interaction saved")
    
  }, error = function(e) {
    log_progress(sprintf("  ⚠ Off-nadir × phenology failed: %s", e$message))
  })
}

# =====================================================================
# 4. VENDOR COMPOSITION (Original)
# =====================================================================

log_subsection("Vendor composition")
vd_chm <- vendor_long(chm_df)
vd_dtm <- vendor_long(dtm_df)
ggsave(file.path(out_plots, "02e_vendor_composition_chm.pdf"),
       p_vendors(vd_chm, "Vendor composition (CHM)"), 
       width=6, height=4, bg="white")
ggsave(file.path(out_plots, "02e_vendor_composition_dtm.pdf"),
       p_vendors(vd_dtm, "Vendor composition (DTM)"), 
       width=6, height=4, bg="white")

# =====================================================================
# 4b. NEW: VENDOR COMPOSITION BY FOREST TYPE
# =====================================================================

log_subsection("Vendor composition (FOREST ONLY)")

tryCatch({
  vd_chm_forest <- vendor_long(chm_forest)
  ggsave(file.path(out_plots, "02e_vendor_composition_chm_FOREST.pdf"),
         p_vendors(vd_chm_forest, "Vendor composition — Forest Types Only"), 
         width=6, height=4, bg="white")
}, error = function(e) {
  log_progress(sprintf("  ⚠ Forest vendor composition failed: %s", e$message))
})

# =====================================================================
# 5. AZIMUTH ANALYSIS (Original)
# =====================================================================

log_subsection("Azimuth relationships")
if (isTRUE(any(is.finite(chm_df$meta_target_azimuth_avg)))) {
  p_rose_chm <- try_plot(
    rose_error_plot(chm_df, "chm_error_mean", "CHM |error| by viewing azimuth"),
    "CHM azimuth rose plot")
  if (!is.null(p_rose_chm)) {
    ggsave(file.path(out_plots, "02f_chm_azimuth_rose_abs_error.pdf"),
           p_rose_chm, width=6, height=6, bg="white")
  }
}
if (isTRUE(any(is.finite(dtm_df$meta_target_azimuth_avg)))) {
  p_rose_dtm <- try_plot(
    rose_error_plot(dtm_df, "dtm_error_mean", "DTM |error| by viewing azimuth"),
    "DTM azimuth rose plot")
  if (!is.null(p_rose_dtm)) {
    ggsave(file.path(out_plots, "02f_dtm_azimuth_rose_abs_error.pdf"),
           p_rose_dtm, width=6, height=6, bg="white")
  }
}

# Cyclic azimuth trends
p_cyc_chm <- try_plot(
  cyclic_trend_plot(chm_df, "meta_target_azimuth_avg", 
                    "chm_error_mean", "w_chm",
                    "CHM error vs azimuth (cyclic GAM)"),
  "CHM cyclic azimuth plot")

p_cyc_dtm <- try_plot(
  cyclic_trend_plot(dtm_df, "meta_target_azimuth_avg", 
                    "dtm_error_mean", "w_dtm",
                    "DTM error vs azimuth (cyclic GAM)"),
  "DTM cyclic azimuth plot")

if (!is.null(p_cyc_chm)) {
  ggsave(file.path(out_plots, "02g_chm_cyclic_azimuth.pdf"), p_cyc_chm, 
         width=8, height=4.5, bg="white")
}
if (!is.null(p_cyc_dtm)) {
  ggsave(file.path(out_plots, "02g_dtm_cyclic_azimuth.pdf"), p_cyc_dtm, 
         width=8, height=4.5, bg="white")
}

# =====================================================================
# 6. MISSINGNESS HEATMAPS (Original)
# =====================================================================

log_subsection("Metadata missingness")

# Build list of columns to check - only those that exist
meta_check_cols_chm <- intersect(
  c("meta_off_nadir_avg","meta_sun_elev_avg","meta_absgeo_z","meta_relgeo_z",
    "meta_tot_z","meta_stereo_z","meta_fwdrev_z","meta_leafon_z",
    "v_GE01","v_WV01","v_WV02","v_WV03"),
  names(chm_df)
)

meta_check_cols_dtm <- intersect(
  c("meta_off_nadir_avg","meta_sun_elev_avg","meta_absgeo_z","meta_relgeo_z",
    "meta_tot_z","meta_stereo_z","meta_fwdrev_z","meta_leafon_z",
    "v_GE01","v_WV01","v_WV02","v_WV03"),
  names(dtm_df)
)

if (length(meta_check_cols_chm) > 0) {
  missing_heat(chm_df, meta_check_cols_chm,
               "Missingness of key MHRSI metadata (CHM)",
               file.path(out_plots, "02h_missingness_mhrsi_chm.pdf"))
}

if (length(meta_check_cols_dtm) > 0) {
  missing_heat(dtm_df, meta_check_cols_dtm,
               "Missingness of key MHRSI metadata (DTM)",
               file.path(out_plots, "02h_missingness_mhrsi_dtm.pdf"))
}

# =====================================================================
# 7. NEW: FOREST-SPECIFIC METADATA SUMMARY STATISTICS
# =====================================================================

log_subsection("Metadata summary statistics by forest type (NEW)")

# Compute metadata effect statistics by forest type
meta_cols_to_check <- c("meta_offnad_z", "meta_sunel_z", "meta_leafon_z", 
                        "meta_stereo_z", "meta_fwdrev_z")
meta_cols_present <- intersect(meta_cols_to_check, names(chm_forest))

if (length(meta_cols_present) > 0) {
  
  # Overall forest metadata stats
  forest_meta_stats <- chm_forest %>%
    group_by(lc_l1_code) %>%
    summarise(
      n = n(),
      across(all_of(meta_cols_present), 
             list(
               mean = ~mean(.x, na.rm = TRUE),
               sd = ~sd(.x, na.rm = TRUE),
               n_valid = ~sum(is.finite(.x))
             ),
             .names = "{.col}_{.fn}"),
      chm_error_mean = mean(chm_error_mean, na.rm = TRUE),
      chm_error_sd = sd(chm_error_mean, na.rm = TRUE),
      .groups = "drop"
    )
  
  write_csv(forest_meta_stats, file.path(out_tables, "metadata_stats_by_forest_type.csv"))
  log_progress("  ✓ Metadata statistics by forest type saved")
  
  # Phenology-level metadata stats
  phenology_meta_stats <- chm_forest %>%
    filter(!is.na(phenology)) %>%
    group_by(phenology) %>%
    summarise(
      n = n(),
      across(all_of(meta_cols_present), 
             list(
               mean = ~mean(.x, na.rm = TRUE),
               sd = ~sd(.x, na.rm = TRUE)
             ),
             .names = "{.col}_{.fn}"),
      chm_error_mean = mean(chm_error_mean, na.rm = TRUE),
      chm_error_sd = sd(chm_error_mean, na.rm = TRUE),
      .groups = "drop"
    )
  
  write_csv(phenology_meta_stats, file.path(out_tables, "metadata_stats_by_phenology.csv"))
  
  log_progress("  Metadata coverage by phenology:")
  print(phenology_meta_stats %>% select(phenology, n, starts_with("meta_leafon")))
}

# =====================================================================
# 8. NEW: COMPREHENSIVE FOREST METADATA SUMMARY TABLE
# =====================================================================

log_subsection("Comprehensive metadata-error relationship table (FOREST)")

# For each metadata variable, compute correlation with error by forest type
if (length(meta_cols_present) > 0) {
  
  meta_error_corr <- chm_forest %>%
    group_by(lc_l1_code) %>%
    summarise(
      n = n(),
      across(all_of(meta_cols_present),
             ~cor(.x, chm_error_mean, use = "pairwise.complete.obs"),
             .names = "cor_{.col}"),
      .groups = "drop"
    )
  
  write_csv(meta_error_corr, file.path(out_tables, "metadata_error_correlations_FOREST.csv"))
  log_progress("  ✓ Metadata-error correlations by forest type saved")
  
  # Also by phenology
  meta_error_corr_phenology <- chm_forest %>%
    filter(!is.na(phenology)) %>%
    group_by(phenology) %>%
    summarise(
      n = n(),
      across(all_of(meta_cols_present),
             ~cor(.x, chm_error_mean, use = "pairwise.complete.obs"),
             .names = "cor_{.col}"),
      .groups = "drop"
    )
  
  write_csv(meta_error_corr_phenology, file.path(out_tables, "metadata_error_correlations_by_phenology.csv"))
  
  log_progress("  Metadata-error correlations by phenology:")
  print(meta_error_corr_phenology)
}

log_progress("✓ Enhanced metadata analysis complete")

# =====================================================================
# SAVE CHECKPOINT
# =====================================================================

if (!exists("BATCH_MODE") || !BATCH_MODE) {
  log_progress("Saving checkpoint for interactive mode...")
  save_checkpoint("05_metadata", list(
    metadata_complete = TRUE,
    forest_classes = FOREST_CLASSES,
    forest_present = forest_present_chm
  ))
}
