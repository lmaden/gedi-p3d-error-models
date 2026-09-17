# =====================================================================
# section_04_terrain_lc.R (ENHANCED VERSION)
# Terrain and land cover relationship analysis
#
# ENHANCEMENTS ADDED:
#   1. Boxplots by site (in addition to by LC)
#   2. Per-LC error statistics table
#   3. Per-site error statistics table
#   4. Slope distribution by LC
#   5. Interaction visualization (slope effect by LC)
#   6. Forest-focused versions of all plots (EBF, BDF, ENF, DNF)
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
  } else {
    stop("Data not available. Run: source('section_01_ingest.R') first")
  }
}

if (!exists("n_bins_common")) {
  n_chm_slope <- sum(is.finite(chm_df$slope_mean) & is.finite(chm_df$chm_error_mean))
  n_dtm_slope <- sum(is.finite(dtm_df$slope_mean) & is.finite(dtm_df$dtm_error_mean))
  n_bins_common <- max(choose_hex_bins(n_chm_slope), choose_hex_bins(n_dtm_slope))
}

log_progress("Analyzing terrain & land cover relationships (ENHANCED)...")

# Define forest classes for focused analysis
FOREST_CLASSES <- c("EBF", "BDF", "ENF", "DNF")
FOREST_LABELS <- c("Broadleaf Evergreen", "Broadleaf Deciduous", 
                   "Needleleaf Evergreen", "Needleleaf Deciduous")
names(FOREST_LABELS) <- FOREST_CLASSES

log_progress(sprintf("  Forest classes for focused analysis: %s", 
                     paste(FOREST_CLASSES, collapse = ", ")))

# Check which forest classes are present
forest_present_chm <- intersect(FOREST_CLASSES, unique(chm_df$lc_l1_code))
forest_present_dtm <- intersect(FOREST_CLASSES, unique(dtm_df$lc_l1_code))
log_progress(sprintf("  CHM forest classes present: %s", paste(forest_present_chm, collapse = ", ")))
log_progress(sprintf("  DTM forest classes present: %s", paste(forest_present_dtm, collapse = ", ")))

# Top land cover classes (all types)
top_lc_chm <- top_levels(chm_df$lc_l1_code, 6L)
top_lc_dtm <- top_levels(dtm_df$lc_l1_code, 6L)

log_progress(sprintf("  Top CHM classes (all): %s", paste(top_lc_chm, collapse=", ")))
log_progress(sprintf("  Top DTM classes (all): %s", paste(top_lc_dtm, collapse=", ")))

# =====================================================================
# 1. HEX PLOTS BY LC (Original - All Top 6)
# =====================================================================

log_subsection("Hex trends by land cover (all top 6)")
p_hex_slope_chm_by_lc <- make_hex_trend_plot(
  chm_df %>% filter(lc_l1_code %in% top_lc_chm),
  "slope_mean", "chm_error_mean", "w_chm",
  "Slope [degrees]", "CHM error (m)", "CHM error vs slope — by LC (top 6)",
  n_bins_x = n_bins_common, n_bins_y = n_bins_common
) + facet_wrap(~ lc_l1_code, ncol = 3)

p_hex_slope_dtm_by_lc <- make_hex_trend_plot(
  dtm_df %>% filter(lc_l1_code %in% top_lc_dtm),
  "slope_mean", "dtm_error_mean", "w_dtm",
  "Slope [degrees]", "DTM error (m)", "DTM error vs slope — by LC (top 6)",
  n_bins_x = n_bins_common, n_bins_y = n_bins_common
) + facet_wrap(~ lc_l1_code, ncol = 3)

ggsave(file.path(out_plots, "03a_chm_hex_slope_by_lc.pdf"), 
       p_hex_slope_chm_by_lc, width = 12, height = 10, bg = "white")
ggsave(file.path(out_plots, "03a_dtm_hex_slope_by_lc.pdf"), 
       p_hex_slope_dtm_by_lc, width = 12, height = 10, bg = "white")

# =====================================================================
# 1b. NEW: HEX PLOTS BY LC (Forest Classes Only)
# =====================================================================

log_subsection("Hex trends by land cover (FOREST ONLY)")

if (length(forest_present_chm) >= 2) {
  p_hex_slope_chm_forest <- make_hex_trend_plot(
    chm_df %>% filter(lc_l1_code %in% forest_present_chm),
    "slope_mean", "chm_error_mean", "w_chm",
    "Slope [degrees]", "CHM error (m)", 
    "CHM error vs slope — Forest Types (EBF, BDF, ENF, DNF)",
    n_bins_x = n_bins_common, n_bins_y = n_bins_common
  ) + facet_wrap(~ lc_l1_code, ncol = 2)
  
  ggsave(file.path(out_plots, "03a_chm_hex_slope_by_lc_FOREST.pdf"), 
         p_hex_slope_chm_forest, width = 10, height = 8, bg = "white")
}

if (length(forest_present_dtm) >= 2) {
  p_hex_slope_dtm_forest <- make_hex_trend_plot(
    dtm_df %>% filter(lc_l1_code %in% forest_present_dtm),
    "slope_mean", "dtm_error_mean", "w_dtm",
    "Slope [degrees]", "DTM error (m)", 
    "DTM error vs slope — Forest Types (EBF, BDF, ENF, DNF)",
    n_bins_x = n_bins_common, n_bins_y = n_bins_common
  ) + facet_wrap(~ lc_l1_code, ncol = 2)
  
  ggsave(file.path(out_plots, "03a_dtm_hex_slope_by_lc_FOREST.pdf"), 
         p_hex_slope_dtm_forest, width = 10, height = 8, bg = "white")
}

# =====================================================================
# 2. BOXPLOTS BY LC (Original - All Classes)
# =====================================================================

log_subsection("Boxplots by land cover (all classes)")
L_clip <- quantile(abs(c(chm_df$chm_error_mean, dtm_df$dtm_error_mean)), 
                   probs = 0.995, na.rm = TRUE)

plot_box_by_lc <- function(df, err_col, lc_col, title, n_per_class = 20000, 
                           lc_filter = NULL, order_by_median = FALSE) {
  set.seed(42)
  dd <- df %>%
    dplyr::select(all_of(c(err_col, lc_col))) %>%
    dplyr::filter(is.finite(.data[[err_col]]))
  
  if (!is.null(lc_filter)) {
    dd <- dd %>% filter(.data[[lc_col]] %in% lc_filter)
  }
  
  dd <- dd %>%
    dplyr::group_by(.data[[lc_col]]) %>%
    dplyr::group_modify(~ dplyr::slice_sample(.x,
                                              n = min(n_per_class, nrow(.x)),
                                              replace = FALSE)) %>%
    dplyr::ungroup()
  
  if (order_by_median) {
    lc_order <- dd %>%
      group_by(.data[[lc_col]]) %>%
      summarise(med = median(.data[[err_col]], na.rm = TRUE), .groups = "drop") %>%
      arrange(med) %>%
      pull(.data[[lc_col]])
    dd[[lc_col]] <- factor(dd[[lc_col]], levels = lc_order)
  }
  
  ggplot(dd, aes(x = .data[[lc_col]], y = .data[[err_col]])) +
    geom_boxplot(outlier.alpha = 0.1, width = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    coord_cartesian(ylim = c(-L_clip, L_clip)) +
    labs(x = "LC Level-1", y = "Error (m)", title = title) +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

p_box_chm <- plot_box_by_lc(chm_df, "chm_error_mean", "lc_l1_code", 
                            "CHM error by LC (boxplots)")
p_box_dtm <- plot_box_by_lc(dtm_df, "dtm_error_mean", "lc_l1_code", 
                            "DTM error by LC (boxplots)")
ggsave(file.path(out_plots, "03b_chm_box_by_lc.pdf"), p_box_chm, 
       width = 10, height = 6, bg = "white")
ggsave(file.path(out_plots, "03b_dtm_box_by_lc.pdf"), p_box_dtm, 
       width = 10, height = 6, bg = "white")

# =====================================================================
# 2b. NEW: BOXPLOTS BY LC (Forest Classes Only)
# =====================================================================

log_subsection("Boxplots by land cover (FOREST ONLY)")

p_box_chm_forest <- plot_box_by_lc(chm_df, "chm_error_mean", "lc_l1_code", 
                                    "CHM error by Forest Type",
                                    lc_filter = forest_present_chm,
                                    order_by_median = TRUE)
p_box_dtm_forest <- plot_box_by_lc(dtm_df, "dtm_error_mean", "lc_l1_code", 
                                    "DTM error by Forest Type",
                                    lc_filter = forest_present_dtm,
                                    order_by_median = TRUE)
ggsave(file.path(out_plots, "03b_chm_box_by_lc_FOREST.pdf"), p_box_chm_forest, 
       width = 8, height = 6, bg = "white")
ggsave(file.path(out_plots, "03b_dtm_box_by_lc_FOREST.pdf"), p_box_dtm_forest, 
       width = 8, height = 6, bg = "white")

# =====================================================================
# 3. BOXPLOTS BY SITE (Original)
# =====================================================================

log_subsection("Boxplots by site")

plot_box_by_site <- function(df, err_col, site_col, title, n_per_site = 10000) {
  set.seed(42)
  dd <- df %>%
    dplyr::select(all_of(c(err_col, site_col))) %>%
    dplyr::filter(is.finite(.data[[err_col]])) %>%
    dplyr::group_by(.data[[site_col]]) %>%
    dplyr::group_modify(~ dplyr::slice_sample(.x,
                                              n = min(n_per_site, nrow(.x)),
                                              replace = FALSE)) %>%
    dplyr::ungroup()
  
  # Order sites by median error
  site_order <- dd %>%
    group_by(.data[[site_col]]) %>%
    summarise(med = median(.data[[err_col]], na.rm = TRUE), .groups = "drop") %>%
    arrange(med) %>%
    pull(.data[[site_col]])
  
  dd[[site_col]] <- factor(dd[[site_col]], levels = site_order)
  
  ggplot(dd, aes(x = .data[[site_col]], y = .data[[err_col]])) +
    geom_boxplot(outlier.alpha = 0.1, width = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    coord_cartesian(ylim = c(-L_clip, L_clip)) +
    labs(x = "Site (ordered by median error)", y = "Error (m)", title = title) +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

p_box_chm_site <- plot_box_by_site(chm_df, "chm_error_mean", "site", 
                                    "CHM error by site (boxplots)")
p_box_dtm_site <- plot_box_by_site(dtm_df, "dtm_error_mean", "site", 
                                    "DTM error by site (boxplots)")
ggsave(file.path(out_plots, "03b_chm_box_by_site.pdf"), p_box_chm_site, 
       width = 12, height = 6, bg = "white")
ggsave(file.path(out_plots, "03b_dtm_box_by_site.pdf"), p_box_dtm_site, 
       width = 12, height = 6, bg = "white")

# =====================================================================
# 4. ECDF BY LC (Original - Top 6)
# =====================================================================

log_subsection("Empirical CDFs by land cover (top 6)")
p_ecdf_chm <- plot_ecdf_by_lc(chm_df, "chm_error_mean", "lc_l1_code", 
                              top_lc_chm, "ECDF of |CHM error| by LC (top 6)")
p_ecdf_dtm <- plot_ecdf_by_lc(dtm_df, "dtm_error_mean", "lc_l1_code", 
                              top_lc_dtm, "ECDF of |DTM error| by LC (top 6)")
ggsave(file.path(out_plots, "03c_chm_ecdf_by_lc.pdf"), p_ecdf_chm, 
       width = 8, height = 6, bg = "white")
ggsave(file.path(out_plots, "03c_dtm_ecdf_by_lc.pdf"), p_ecdf_dtm, 
       width = 8, height = 6, bg = "white")

# =====================================================================
# 4b. NEW: ECDF BY LC (Forest Classes Only)
# =====================================================================

log_subsection("Empirical CDFs by land cover (FOREST ONLY)")

p_ecdf_chm_forest <- plot_ecdf_by_lc(chm_df, "chm_error_mean", "lc_l1_code", 
                                      forest_present_chm, 
                                      "ECDF of |CHM error| by Forest Type")
p_ecdf_dtm_forest <- plot_ecdf_by_lc(dtm_df, "dtm_error_mean", "lc_l1_code", 
                                      forest_present_dtm, 
                                      "ECDF of |DTM error| by Forest Type")
ggsave(file.path(out_plots, "03c_chm_ecdf_by_lc_FOREST.pdf"), p_ecdf_chm_forest, 
       width = 8, height = 6, bg = "white")
ggsave(file.path(out_plots, "03c_dtm_ecdf_by_lc_FOREST.pdf"), p_ecdf_dtm_forest, 
       width = 8, height = 6, bg = "white")

# =====================================================================
# 5. PER-LC ERROR STATISTICS TABLE
# =====================================================================

log_subsection("Per-LC error statistics")

compute_error_stats_by_group <- function(df, err_col, group_col) {
  df %>%
    filter(is.finite(.data[[err_col]])) %>%
    group_by(.data[[group_col]]) %>%
    summarise(
      n = n(),
      mean = mean(.data[[err_col]], na.rm = TRUE),
      median = median(.data[[err_col]], na.rm = TRUE),
      sd = sd(.data[[err_col]], na.rm = TRUE),
      nmad = 1.4826 * median(abs(.data[[err_col]] - median(.data[[err_col]], na.rm = TRUE)), na.rm = TRUE),
      mae = mean(abs(.data[[err_col]]), na.rm = TRUE),
      rmse = sqrt(mean(.data[[err_col]]^2, na.rm = TRUE)),
      q05 = quantile(.data[[err_col]], 0.05, na.rm = TRUE),
      q95 = quantile(.data[[err_col]], 0.95, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(n))
}

lc_stats_chm <- compute_error_stats_by_group(chm_df, "chm_error_mean", "lc_l1_code")
lc_stats_dtm <- compute_error_stats_by_group(dtm_df, "dtm_error_mean", "lc_l1_code")

write_csv(lc_stats_chm, file.path(out_tables, "error_stats_by_lc_chm.csv"))
write_csv(lc_stats_dtm, file.path(out_tables, "error_stats_by_lc_dtm.csv"))

log_progress("  CHM error stats by LC:")
print(lc_stats_chm %>% select(lc_l1_code, n, mean, sd, rmse))

# Forest-only stats
lc_stats_chm_forest <- lc_stats_chm %>% filter(lc_l1_code %in% FOREST_CLASSES)
lc_stats_dtm_forest <- lc_stats_dtm %>% filter(lc_l1_code %in% FOREST_CLASSES)

write_csv(lc_stats_chm_forest, file.path(out_tables, "error_stats_by_lc_chm_FOREST.csv"))
write_csv(lc_stats_dtm_forest, file.path(out_tables, "error_stats_by_lc_dtm_FOREST.csv"))

log_progress("  CHM error stats by FOREST type:")
print(lc_stats_chm_forest %>% select(lc_l1_code, n, mean, sd, rmse))

# =====================================================================
# 6. PER-SITE ERROR STATISTICS TABLE
# =====================================================================

log_subsection("Per-site error statistics")

site_stats_chm <- compute_error_stats_by_group(chm_df, "chm_error_mean", "site")
site_stats_dtm <- compute_error_stats_by_group(dtm_df, "dtm_error_mean", "site")

write_csv(site_stats_chm, file.path(out_tables, "error_stats_by_site_chm.csv"))
write_csv(site_stats_dtm, file.path(out_tables, "error_stats_by_site_dtm.csv"))

log_progress("  CHM error stats by site:")
print(site_stats_chm %>% select(site, n, mean, sd, rmse))

# Flag sites with high absolute mean error
high_bias_sites <- site_stats_chm %>% filter(abs(mean) > 5)
if (nrow(high_bias_sites) > 0) {
  log_progress(sprintf("  ⚠ %d sites with |mean error| > 5m:", nrow(high_bias_sites)))
  for (i in 1:nrow(high_bias_sites)) {
    log_progress(sprintf("    Site %s: mean=%.2f, sd=%.2f, n=%d",
                         high_bias_sites$site[i], high_bias_sites$mean[i],
                         high_bias_sites$sd[i], high_bias_sites$n[i]))
  }
}

# =====================================================================
# 7. SLOPE DISTRIBUTION BY LC (Original - Top 6)
# =====================================================================

log_subsection("Slope distribution by LC")

p_slope_by_lc <- chm_df %>%
  filter(lc_l1_code %in% top_lc_chm) %>%
  ggplot(aes(x = slope_mean, fill = lc_l1_code)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(0, 50)) +
  labs(x = "Slope (degrees)", y = "Density", 
       title = "Slope distribution by land cover class (top 6)",
       fill = "LC Class") +
  theme_cowplot()

ggsave(file.path(out_plots, "03d_slope_density_by_lc.pdf"), p_slope_by_lc,
       width = 10, height = 6, bg = "white")

# =====================================================================
# 7b. NEW: SLOPE DISTRIBUTION BY LC (Forest Classes Only)
# =====================================================================

log_subsection("Slope distribution by LC (FOREST ONLY)")

p_slope_by_lc_forest <- chm_df %>%
  filter(lc_l1_code %in% forest_present_chm) %>%
  ggplot(aes(x = slope_mean, fill = lc_l1_code)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(0, 50)) +
  scale_fill_brewer(palette = "Dark2", labels = FOREST_LABELS[forest_present_chm]) +
  labs(x = "Slope (degrees)", y = "Density", 
       title = "Slope distribution by Forest Type",
       fill = "Forest Type") +
  theme_cowplot()

ggsave(file.path(out_plots, "03d_slope_density_by_lc_FOREST.pdf"), p_slope_by_lc_forest,
       width = 10, height = 6, bg = "white")

# =====================================================================
# 8. SLOPE × LC INTERACTION VISUALIZATION (Original - Top 6)
# =====================================================================

log_subsection("Slope × LC interaction visualization")

# Compute mean error in slope bins for each LC
interaction_data <- chm_df %>%
  filter(lc_l1_code %in% top_lc_chm, is.finite(slope_mean), is.finite(chm_error_mean)) %>%
  mutate(slope_bin = cut(slope_mean, breaks = seq(0, 40, 5), include.lowest = TRUE)) %>%
  filter(!is.na(slope_bin)) %>%
  group_by(lc_l1_code, slope_bin) %>%
  summarise(
    n = n(),
    mean_error = weighted.mean(chm_error_mean, w_chm, na.rm = TRUE),
    se_error = sd(chm_error_mean, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  filter(n >= 50)  # Only keep bins with enough data

p_interaction <- ggplot(interaction_data, 
                        aes(x = slope_bin, y = mean_error, color = lc_l1_code, group = lc_l1_code)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_error - 1.96 * se_error, 
                    ymax = mean_error + 1.96 * se_error), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Slope bin (degrees)", y = "Mean CHM error (m)",
       title = "CHM error vs slope by land cover class (top 6)",
       subtitle = "Lines show how slope effect differs by LC (interaction)",
       color = "LC Class") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(out_plots, "03e_slope_lc_interaction.pdf"), p_interaction,
       width = 10, height = 6, bg = "white")

# Check if interaction is substantial
log_progress("  Checking for slope × LC interaction...")
slope_effect_by_lc <- interaction_data %>%
  group_by(lc_l1_code) %>%
  summarise(
    slope_range = diff(range(mean_error)),
    .groups = "drop"
  )

max_slope_effect <- max(slope_effect_by_lc$slope_range)
min_slope_effect <- min(slope_effect_by_lc$slope_range)

if (max_slope_effect / min_slope_effect > 2) {
  log_progress(sprintf("  ⚠ Slope effect varies substantially by LC (range: %.2f to %.2f m)",
                       min_slope_effect, max_slope_effect))
  log_progress("    Consider including slope:lc_l1_code interaction in model")
} else {
  log_progress(sprintf("  Slope effect relatively consistent across LC (range: %.2f to %.2f m)",
                       min_slope_effect, max_slope_effect))
}

# =====================================================================
# 8b. NEW: SLOPE × LC INTERACTION (Forest Classes Only)
# =====================================================================

log_subsection("Slope × LC interaction visualization (FOREST ONLY)")

interaction_data_forest <- chm_df %>%
  filter(lc_l1_code %in% forest_present_chm, is.finite(slope_mean), is.finite(chm_error_mean)) %>%
  mutate(slope_bin = cut(slope_mean, breaks = seq(0, 40, 5), include.lowest = TRUE)) %>%
  filter(!is.na(slope_bin)) %>%
  group_by(lc_l1_code, slope_bin) %>%
  summarise(
    n = n(),
    mean_error = weighted.mean(chm_error_mean, w_chm, na.rm = TRUE),
    se_error = sd(chm_error_mean, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  filter(n >= 50)

p_interaction_forest <- ggplot(interaction_data_forest, 
                                aes(x = slope_bin, y = mean_error, color = lc_l1_code, group = lc_l1_code)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_error - 1.96 * se_error, 
                    ymax = mean_error + 1.96 * se_error), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_brewer(palette = "Dark2", labels = FOREST_LABELS[forest_present_chm]) +
  labs(x = "Slope bin (degrees)", y = "Mean CHM error (m)",
       title = "CHM error vs slope by Forest Type",
       subtitle = "Comparing Evergreen vs Deciduous, Broadleaf vs Needleleaf",
       color = "Forest Type") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(out_plots, "03e_slope_lc_interaction_FOREST.pdf"), p_interaction_forest,
       width = 10, height = 6, bg = "white")

# Forest-specific interaction check
log_progress("  Checking for slope × forest type interaction...")
slope_effect_forest <- interaction_data_forest %>%
  group_by(lc_l1_code) %>%
  summarise(
    slope_range = diff(range(mean_error)),
    min_error = min(mean_error),
    max_error = max(mean_error),
    .groups = "drop"
  )
log_progress("  Slope effect range by forest type:")
for (i in 1:nrow(slope_effect_forest)) {
  log_progress(sprintf("    %s: %.2f to %.2f m (range: %.2f m)",
                       slope_effect_forest$lc_l1_code[i],
                       slope_effect_forest$min_error[i],
                       slope_effect_forest$max_error[i],
                       slope_effect_forest$slope_range[i]))
}

# =====================================================================
# 9. ECOREGION ANALYSIS (Original)
# =====================================================================

log_subsection("Error by ecoregion")

if ("ecoregion" %in% names(chm_df)) {
  eco_stats_chm <- compute_error_stats_by_group(chm_df, "chm_error_mean", "ecoregion")
  write_csv(eco_stats_chm, file.path(out_tables, "error_stats_by_ecoregion_chm.csv"))
  
  # Boxplot by ecoregion
  p_box_eco <- chm_df %>%
    filter(is.finite(chm_error_mean)) %>%
    ggplot(aes(x = reorder(ecoregion, chm_error_mean, FUN = median), y = chm_error_mean)) +
    geom_boxplot(outlier.alpha = 0.05) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    coord_flip(ylim = c(-L_clip, L_clip)) +
    labs(x = "Ecoregion", y = "CHM error (m)",
         title = "CHM error by ecoregion") +
    theme_cowplot()
  
  ggsave(file.path(out_plots, "03f_chm_box_by_ecoregion.pdf"), p_box_eco,
         width = 10, height = 8, bg = "white")
  
  log_progress("  CHM error stats by ecoregion:")
  print(eco_stats_chm %>% select(ecoregion, n, mean, sd, rmse))
}

# =====================================================================
# 10. NEW: DECIDUOUS VS EVERGREEN COMPARISON
# =====================================================================

log_subsection("Deciduous vs Evergreen comparison (NEW)")

# Create leaf type variable
chm_forest <- chm_df %>%
  filter(lc_l1_code %in% FOREST_CLASSES) %>%
  mutate(
    leaf_type = case_when(
      lc_l1_code %in% c("EBF", "ENF") ~ "Evergreen",
      lc_l1_code %in% c("BDF", "DNF") ~ "Deciduous",
      TRUE ~ NA_character_
    ),
    needle_type = case_when(
      lc_l1_code %in% c("ENF", "DNF") ~ "Needleleaf",
      lc_l1_code %in% c("EBF", "BDF") ~ "Broadleaf",
      TRUE ~ NA_character_
    )
  )

# Boxplot: Deciduous vs Evergreen
p_decid_ever <- chm_forest %>%
  filter(!is.na(leaf_type), is.finite(chm_error_mean)) %>%
  ggplot(aes(x = leaf_type, y = chm_error_mean, fill = leaf_type)) +
  geom_boxplot(outlier.alpha = 0.05, width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  coord_cartesian(ylim = c(-L_clip, L_clip)) +
  scale_fill_manual(values = c("Deciduous" = "#E69F00", "Evergreen" = "#009E73")) +
  labs(x = "", y = "CHM error (m)",
       title = "CHM error: Deciduous vs Evergreen Forests") +
  theme_cowplot() +
  theme(legend.position = "none")

# Boxplot: Broadleaf vs Needleleaf
p_broad_needle <- chm_forest %>%
  filter(!is.na(needle_type), is.finite(chm_error_mean)) %>%
  ggplot(aes(x = needle_type, y = chm_error_mean, fill = needle_type)) +
  geom_boxplot(outlier.alpha = 0.05, width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  coord_cartesian(ylim = c(-L_clip, L_clip)) +
  scale_fill_manual(values = c("Broadleaf" = "#56B4E9", "Needleleaf" = "#D55E00")) +
  labs(x = "", y = "CHM error (m)",
       title = "CHM error: Broadleaf vs Needleleaf Forests") +
  theme_cowplot() +
  theme(legend.position = "none")

# Combined plot
p_forest_comparison <- p_decid_ever | p_broad_needle
ggsave(file.path(out_plots, "03g_forest_type_comparison.pdf"), p_forest_comparison,
       width = 12, height = 6, bg = "white")

# Summary statistics
log_progress("  Forest type comparison:")
decid_ever_stats <- chm_forest %>%
  filter(!is.na(leaf_type), is.finite(chm_error_mean)) %>%
  group_by(leaf_type) %>%
  summarise(
    n = n(),
    mean = mean(chm_error_mean),
    median = median(chm_error_mean),
    sd = sd(chm_error_mean),
    .groups = "drop"
  )
print(decid_ever_stats)

broad_needle_stats <- chm_forest %>%
  filter(!is.na(needle_type), is.finite(chm_error_mean)) %>%
  group_by(needle_type) %>%
  summarise(
    n = n(),
    mean = mean(chm_error_mean),
    median = median(chm_error_mean),
    sd = sd(chm_error_mean),
    .groups = "drop"
  )
print(broad_needle_stats)

# =====================================================================
# 11. NEW: SLOPE × FOREST TYPE FACETED HEX (2x2 grid)
# =====================================================================

log_subsection("Slope × Forest type faceted hex (2x2 grid)")

if (length(forest_present_chm) >= 2) {
  # Create ordered factor for nice facet labels
  chm_forest_hex <- chm_df %>%
    filter(lc_l1_code %in% forest_present_chm) %>%
    mutate(forest_label = factor(lc_l1_code, 
                                  levels = c("BDF", "EBF", "DNF", "ENF"),
                                  labels = c("Broadleaf\nDeciduous", "Broadleaf\nEvergreen",
                                            "Needleleaf\nDeciduous", "Needleleaf\nEvergreen")))
  
  p_hex_forest_2x2 <- ggplot(chm_forest_hex, aes(x = slope_mean, y = chm_error_mean)) +
    geom_hex(bins = 40) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "white", linewidth = 0.5) +
    geom_smooth(method = "gam", color = "black", linewidth = 1) +
    scale_fill_viridis_c(trans = "log10", name = "Count") +
    facet_wrap(~ forest_label, ncol = 2) +
    coord_cartesian(xlim = c(0, 40), ylim = c(-15, 15)) +
    labs(x = "Slope (degrees)", y = "CHM error (m)",
         title = "CHM Error vs Slope by Forest Type",
         subtitle = "Comparing phenology (rows) and leaf structure (columns)") +
    theme_cowplot() +
    theme(strip.text = element_text(size = 11))
  
  ggsave(file.path(out_plots, "03h_chm_hex_forest_2x2.pdf"), p_hex_forest_2x2,
         width = 10, height = 8, bg = "white")
}

log_progress("✓ Enhanced terrain & LC analysis complete")

# =====================================================================
# SAVE CHECKPOINT
# =====================================================================

if (!exists("BATCH_MODE") || !BATCH_MODE) {
  log_progress("Saving checkpoint for interactive mode...")
  save_checkpoint("04_terrain_lc", list(
    terrain_complete = TRUE,
    top_lc_chm = top_lc_chm,
    top_lc_dtm = top_lc_dtm,
    forest_classes = FOREST_CLASSES,
    forest_present_chm = forest_present_chm,
    forest_present_dtm = forest_present_dtm,
    lc_stats_chm = lc_stats_chm,
    lc_stats_dtm = lc_stats_dtm,
    lc_stats_chm_forest = lc_stats_chm_forest,
    lc_stats_dtm_forest = lc_stats_dtm_forest,
    site_stats_chm = site_stats_chm,
    site_stats_dtm = site_stats_dtm
  ))
}
