# =====================================================================
# section_03_eda.R (ENHANCED VERSION)
# Core exploratory data analysis
#
# ENHANCEMENTS ADDED:
#   1. Comprehensive error distribution statistics table
#   2. Bivariate distribution checks (predictor joint distributions)
#   3. Leverage point identification (hat values, Cook's distance)
#   4. Per-site balance diagnostics
#   5. Forest-focused versions of key plots (EBF, BDF, ENF, DNF)
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

# Load moments if not already loaded
if (!requireNamespace("moments", quietly = TRUE)) {
  install.packages("moments", repos = "https://cloud.r-project.org")
}
library(moments)

log_progress("Running core EDA (ENHANCED)...")

# Define forest classes for focused analysis
FOREST_CLASSES <- c("EBF", "BDF", "ENF", "DNF")
FOREST_LABELS <- c("Broadleaf Evergreen", "Broadleaf Deciduous", 
                   "Needleleaf Evergreen", "Needleleaf Deciduous")
names(FOREST_LABELS) <- FOREST_CLASSES

# Check which forest classes are present
forest_present_chm <- intersect(FOREST_CLASSES, unique(chm_df$lc_l1_code))
forest_present_dtm <- intersect(FOREST_CLASSES, unique(dtm_df$lc_l1_code))
log_progress(sprintf("  Forest classes present: %s", paste(forest_present_chm, collapse = ", ")))

# Create forest-only subsets
chm_forest <- chm_df %>% filter(lc_l1_code %in% FOREST_CLASSES)
dtm_forest <- dtm_df %>% filter(lc_l1_code %in% FOREST_CLASSES)
log_progress(sprintf("  Forest observations: CHM=%s, DTM=%s", 
                     format(nrow(chm_forest), big.mark=","),
                     format(nrow(dtm_forest), big.mark=",")))

# =====================================================================
# 1. ERROR DISTRIBUTIONS (Original + Enhanced)
# =====================================================================

log_subsection("Error distributions")
p_chm_histqq <- plot_hist_qq(chm_df$chm_error_mean, "CHM error (P3D − ALS)")
p_dtm_histqq <- plot_hist_qq(dtm_df$dtm_error_mean, "DTM error (P3D − 3DEP)", bins = 40)
ggsave(file.path(out_plots, "01_chm_hist_qq.pdf"), p_chm_histqq, 
       width = 10, height = 4, bg = "white")
ggsave(file.path(out_plots, "01_dtm_hist_qq.pdf"), p_dtm_histqq, 
       width = 10, height = 4, bg = "white")

# =====================================================================
# 1b. NEW: FOREST-ONLY ERROR DISTRIBUTIONS
# =====================================================================

log_subsection("Error distributions (FOREST ONLY)")
p_chm_histqq_forest <- plot_hist_qq(chm_forest$chm_error_mean, 
                                     "CHM error (P3D − ALS) — Forest Types Only")
p_dtm_histqq_forest <- plot_hist_qq(dtm_forest$dtm_error_mean, 
                                     "DTM error (P3D − 3DEP) — Forest Types Only", bins = 40)
ggsave(file.path(out_plots, "01_chm_hist_qq_FOREST.pdf"), p_chm_histqq_forest, 
       width = 10, height = 4, bg = "white")
ggsave(file.path(out_plots, "01_dtm_hist_qq_FOREST.pdf"), p_dtm_histqq_forest, 
       width = 10, height = 4, bg = "white")

# =====================================================================
# 2. COMPREHENSIVE ERROR DISTRIBUTION STATISTICS
# =====================================================================

log_subsection("Comprehensive error distribution analysis")

compute_error_stats <- function(errors, label) {
  errors <- errors[is.finite(errors)]
  tibble::tibble(
    product = label,
    n = length(errors),
    mean = mean(errors),
    median = median(errors),
    sd = sd(errors),
    nmad = nm_ad(errors),
    iqr = IQR(errors),
    q01 = quantile(errors, 0.01),
    q05 = quantile(errors, 0.05),
    q25 = quantile(errors, 0.25),
    q75 = quantile(errors, 0.75),
    q95 = quantile(errors, 0.95),
    q99 = quantile(errors, 0.99),
    skewness = moments::skewness(errors),
    kurtosis = moments::kurtosis(errors),
    mae = mean(abs(errors)),
    rmse = sqrt(mean(errors^2)),
    q95_abs = q95_abs(errors)
  )
}

# All data stats
error_stats <- bind_rows(
  compute_error_stats(chm_df$chm_error_mean, "CHM"),
  compute_error_stats(dtm_df$dtm_error_mean, "DTM")
)

write_csv(error_stats, file.path(out_tables, "error_distribution_statistics.csv"))
log_progress("  Error distribution statistics (all data):")
print(error_stats %>% select(product, mean, median, sd, nmad, kurtosis, rmse))

# Forest-only stats
error_stats_forest <- bind_rows(
  compute_error_stats(chm_forest$chm_error_mean, "CHM_FOREST"),
  compute_error_stats(dtm_forest$dtm_error_mean, "DTM_FOREST")
)

write_csv(error_stats_forest, file.path(out_tables, "error_distribution_statistics_FOREST.csv"))
log_progress("  Error distribution statistics (FOREST ONLY):")
print(error_stats_forest %>% select(product, mean, median, sd, nmad, kurtosis, rmse))

# =====================================================================
# 3. HEX BIN SIZING (Original)
# =====================================================================

log_subsection("Computing optimal hex bin sizes")
n_chm_slope <- sum(is.finite(chm_df$slope_mean) & is.finite(chm_df$chm_error_mean))
n_dtm_slope <- sum(is.finite(dtm_df$slope_mean) & is.finite(dtm_df$dtm_error_mean))
n_chm_wsci  <- sum(is.finite(chm_df$wsci_z) & is.finite(chm_df$chm_error_mean))
n_dtm_wsci  <- sum(is.finite(dtm_df$wsci_z) & is.finite(dtm_df$dtm_error_mean))
n_bins_common <- max(choose_hex_bins(n_chm_slope), choose_hex_bins(n_dtm_slope),
                     choose_hex_bins(n_chm_wsci),  choose_hex_bins(n_dtm_wsci))

log_progress(sprintf("  Using %d hex bins for consistency", n_bins_common))

# =====================================================================
# 4. TERRAIN RELATIONSHIPS (Original)
# =====================================================================

log_subsection("Terrain relationships (slope, WSCI)")
p_hex_slope_chm <- make_hex_trend_plot(
  chm_df, "slope_mean","chm_error_mean","w_chm",
  "Slope [degrees]","CHM error (m)", "CHM error vs slope — hex + trend",
  n_bins_x = n_bins_common, n_bins_y = n_bins_common)

p_hex_wsci_chm <- make_hex_trend_plot(
  chm_df, "wsci_z","chm_error_mean","w_chm",
  "WSCI (z)","CHM error (m)", "CHM error vs WSCI — hex + trend",
  n_bins_x = n_bins_common, n_bins_y = n_bins_common)

p_hex_slope_dtm <- make_hex_trend_plot(
  dtm_df, "slope_mean","dtm_error_mean","w_dtm",
  "Slope [degrees]","DTM error (m)", "DTM error vs slope — hex + trend",
  n_bins_x = n_bins_common, n_bins_y = n_bins_common)

p_hex_wsci_dtm <- make_hex_trend_plot(
  dtm_df, "wsci_z","dtm_error_mean","w_dtm",
  "WSCI (z)","DTM error (m)", "DTM error vs WSCI — hex + trend",
  n_bins_x = n_bins_common, n_bins_y = n_bins_common)

# Sync hex scales
hex_global_max <- max(attr(p_hex_slope_chm, "hex_max"), attr(p_hex_wsci_chm, "hex_max"),
                      attr(p_hex_slope_dtm, "hex_max"), attr(p_hex_wsci_dtm, "hex_max"),
                      na.rm = TRUE)
sync_hex_scale <- scale_fill_viridis_c(
  trans  = "log10",
  breaks = scales::trans_breaks("log10", function(x) 10^x),
  labels = scales::label_number(accuracy = 1),
  limits = c(1, hex_global_max),
  name   = "count (log10)"
)
p_hex_slope_chm <- p_hex_slope_chm + sync_hex_scale
p_hex_wsci_chm  <- p_hex_wsci_chm  + sync_hex_scale
p_hex_slope_dtm <- p_hex_slope_dtm + sync_hex_scale
p_hex_wsci_dtm  <- p_hex_wsci_dtm  + sync_hex_scale

ggsave(file.path(out_plots, "02_chm_hex_slope_wsci.pdf"),
       (p_chm_histqq / (p_hex_slope_chm | p_hex_wsci_chm)), 
       width = 12, height = 10, bg = "white")
ggsave(file.path(out_plots, "02_dtm_hex_slope_wsci.pdf"),
       (p_dtm_histqq / (p_hex_slope_dtm | p_hex_wsci_dtm)), 
       width = 12, height = 10, bg = "white")

# =====================================================================
# 4b. NEW: TERRAIN RELATIONSHIPS (FOREST ONLY)
# =====================================================================

log_subsection("Terrain relationships (FOREST ONLY)")

p_hex_slope_chm_forest <- make_hex_trend_plot(
  chm_forest, "slope_mean","chm_error_mean","w_chm",
  "Slope [degrees]","CHM error (m)", "CHM error vs slope — Forest Types Only",
  n_bins_x = n_bins_common, n_bins_y = n_bins_common)

p_hex_wsci_chm_forest <- make_hex_trend_plot(
  chm_forest, "wsci_z","chm_error_mean","w_chm",
  "WSCI (z)","CHM error (m)", "CHM error vs WSCI — Forest Types Only",
  n_bins_x = n_bins_common, n_bins_y = n_bins_common)

p_hex_slope_dtm_forest <- make_hex_trend_plot(
  dtm_forest, "slope_mean","dtm_error_mean","w_dtm",
  "Slope [degrees]","DTM error (m)", "DTM error vs slope — Forest Types Only",
  n_bins_x = n_bins_common, n_bins_y = n_bins_common)

p_hex_wsci_dtm_forest <- make_hex_trend_plot(
  dtm_forest, "wsci_z","dtm_error_mean","w_dtm",
  "WSCI (z)","DTM error (m)", "DTM error vs WSCI — Forest Types Only",
  n_bins_x = n_bins_common, n_bins_y = n_bins_common)

ggsave(file.path(out_plots, "02_chm_hex_slope_wsci_FOREST.pdf"),
       (p_chm_histqq_forest / (p_hex_slope_chm_forest | p_hex_wsci_chm_forest)), 
       width = 12, height = 10, bg = "white")
ggsave(file.path(out_plots, "02_dtm_hex_slope_wsci_FOREST.pdf"),
       (p_dtm_histqq_forest / (p_hex_slope_dtm_forest | p_hex_wsci_dtm_forest)), 
       width = 12, height = 10, bg = "white")

# =====================================================================
# 5. BIVARIATE DISTRIBUTION CHECKS
# =====================================================================

log_subsection("Bivariate predictor distributions")

# Sample for faster plotting
set.seed(42)
chm_sample <- chm_df %>% 
  filter(is.finite(slope_mean) & is.finite(wsci_z)) %>%
  sample_n(min(50000, n()))

# Slope × WSCI (key interaction in model)
p_bivar_slope_wsci <- ggplot(chm_sample, aes(x = slope_mean, y = wsci_z)) +
  geom_hex(bins = 50) +
  scale_fill_viridis_c(trans = "log10", name = "Count") +
  labs(title = "Joint distribution: Slope × WSCI",
       subtitle = "Check for sparse regions that may cause extrapolation",
       x = "Slope (degrees)", y = "WSCI (z-scaled)") +
  theme_cowplot()

ggsave(file.path(out_plots, "01e_bivar_slope_wsci.pdf"), p_bivar_slope_wsci,
       width = 8, height = 6, bg = "white")

# Slope × Cover (another important pair)
p_bivar_slope_cover <- ggplot(chm_sample, aes(x = slope_mean, y = cover_z)) +
  geom_hex(bins = 50) +
  scale_fill_viridis_c(trans = "log10", name = "Count") +
  labs(title = "Joint distribution: Slope × Canopy Cover",
       x = "Slope (degrees)", y = "Cover (z-scaled)") +
  theme_cowplot()

ggsave(file.path(out_plots, "01e_bivar_slope_cover.pdf"), p_bivar_slope_cover,
       width = 8, height = 6, bg = "white")

# Check for sparse corners (potential extrapolation zones)
slope_q90 <- quantile(chm_df$slope_mean, 0.90, na.rm = TRUE)
wsci_q90 <- quantile(chm_df$wsci_z, 0.90, na.rm = TRUE)
n_high_slope_high_wsci <- sum(chm_df$slope_mean > slope_q90 & chm_df$wsci_z > wsci_q90, na.rm = TRUE)
pct_corner <- 100 * n_high_slope_high_wsci / nrow(chm_df)

log_progress(sprintf("  High slope + high WSCI (>90th pctile both): %d (%.2f%%)",
                     n_high_slope_high_wsci, pct_corner))
if (pct_corner < 1) {
  log_progress("  ⚠ Sparse corner detected - model may extrapolate in steep+complex terrain")
}

# =====================================================================
# 5b. NEW: BIVARIATE DISTRIBUTION (FOREST ONLY)
# =====================================================================

log_subsection("Bivariate predictor distributions (FOREST ONLY)")

chm_forest_sample <- chm_forest %>% 
  filter(is.finite(slope_mean) & is.finite(wsci_z)) %>%
  sample_n(min(50000, n()))

p_bivar_slope_wsci_forest <- ggplot(chm_forest_sample, aes(x = slope_mean, y = wsci_z)) +
  geom_hex(bins = 50) +
  scale_fill_viridis_c(trans = "log10", name = "Count") +
  labs(title = "Joint distribution: Slope × WSCI (Forest Types Only)",
       x = "Slope (degrees)", y = "WSCI (z-scaled)") +
  theme_cowplot()

ggsave(file.path(out_plots, "01e_bivar_slope_wsci_FOREST.pdf"), p_bivar_slope_wsci_forest,
       width = 8, height = 6, bg = "white")

# =====================================================================
# 6. CORRELATIONS (Original)
# =====================================================================

log_subsection("Correlation matrices")
corr_cols_chm <- c("chm_error_mean","slope_mean_z","slope_sd_z","wsci_z",
                   "rh_98_z","cover_z","aspect_sin_z","aspect_cos_z",
                   available_meta_chm,
                   "view_az_sin","view_az_cos")
corr_cols_dtm <- c("dtm_error_mean","slope_mean_z","slope_sd_z","wsci_z",
                   "rh_98_z","cover_z","aspect_sin_z","aspect_cos_z",
                   available_meta_dtm,
                   "view_az_sin","view_az_cos")

corr_plot(chm_df, corr_cols_chm, "CHM correlation matrix (Pearson)", 
          file.path(out_plots,"04a_chm_corr_pearson.pdf"))
corr_plot(dtm_df, corr_cols_dtm, "DTM correlation matrix (Pearson)", 
          file.path(out_plots,"04a_dtm_corr_pearson.pdf"))

# =====================================================================
# 6b. NEW: CORRELATIONS (FOREST ONLY)
# =====================================================================

log_subsection("Correlation matrices (FOREST ONLY)")
corr_plot(chm_forest, corr_cols_chm, "CHM correlation matrix (Pearson) — Forest Types Only", 
          file.path(out_plots,"04a_chm_corr_pearson_FOREST.pdf"))
corr_plot(dtm_forest, corr_cols_dtm, "DTM correlation matrix (Pearson) — Forest Types Only", 
          file.path(out_plots,"04a_dtm_corr_pearson_FOREST.pdf"))

# =====================================================================
# 7. PAIRS PLOTS (Original)
# =====================================================================

log_subsection("Pairs plots (this may take a few minutes)")
make_pairs(chm_df, corr_cols_chm, file.path(out_plots, "04b_chm_pairs.pdf"))
make_pairs(dtm_df, corr_cols_dtm, file.path(out_plots, "04b_dtm_pairs.pdf"))

# =====================================================================
# 8. PCA (Original)
# =====================================================================

log_subsection("Principal component analysis")
driver_cols <- c("slope_mean_z","slope_sd_z","wsci_z","rh_98_z","cover_z",
                 "aspect_sin_z","aspect_cos_z",
                 available_meta_chm,
                 "view_az_sin","view_az_cos")
do_pca_plot(chm_df, driver_cols, "CHM driver PCA", 
            file.path(out_plots, "04c_chm_driver_pca.pdf"))
do_pca_plot(dtm_df, driver_cols, "DTM driver PCA", 
            file.path(out_plots, "04c_dtm_driver_pca.pdf"))

# =====================================================================
# 9. BIVARIATE MAPS (Original)
# =====================================================================

log_subsection("Bivariate mean error maps")
L_bivar <- quantile(abs(c(chm_df$chm_error_mean, dtm_df$dtm_error_mean)), 
                    probs = 0.98, na.rm = TRUE)

p_bivar_chm_sw <- bivar_mean_map(
  chm_df, "slope_mean", "wsci_z", "w_chm", "chm_error_mean",
  nx = 80, ny = 80, min_n = 50,
  xlab = "Slope [degrees]", ylab = "WSCI (z)",
  title = "CHM: weighted mean error in (slope × WSCI)",
  fixed_limits = c(-L_bivar, L_bivar)
)

p_bivar_dtm_sw <- bivar_mean_map(
  dtm_df, "slope_mean", "wsci_z", "w_dtm", "dtm_error_mean",
  nx = 80, ny = 80, min_n = 50,
  xlab = "Slope [degrees]", ylab = "WSCI (z)",
  title = "DTM: weighted mean error in (slope × WSCI)",
  fixed_limits = c(-L_bivar, L_bivar)
)

ggsave(file.path(out_plots, "05a_chm_bivar_slope_wsci.pdf"), 
       p_bivar_chm_sw, width = 6, height = 5, bg = "white")
ggsave(file.path(out_plots, "05a_dtm_bivar_slope_wsci.pdf"), 
       p_bivar_dtm_sw, width = 6, height = 5, bg = "white")

# =====================================================================
# 9b. NEW: BIVARIATE MAPS (FOREST ONLY)
# =====================================================================

log_subsection("Bivariate mean error maps (FOREST ONLY)")

p_bivar_chm_sw_forest <- bivar_mean_map(
  chm_forest, "slope_mean", "wsci_z", "w_chm", "chm_error_mean",
  nx = 80, ny = 80, min_n = 30,
  xlab = "Slope [degrees]", ylab = "WSCI (z)",
  title = "CHM: weighted mean error (slope × WSCI) — Forest Types Only",
  fixed_limits = c(-L_bivar, L_bivar)
)

p_bivar_dtm_sw_forest <- bivar_mean_map(
  dtm_forest, "slope_mean", "wsci_z", "w_dtm", "dtm_error_mean",
  nx = 80, ny = 80, min_n = 30,
  xlab = "Slope [degrees]", ylab = "WSCI (z)",
  title = "DTM: weighted mean error (slope × WSCI) — Forest Types Only",
  fixed_limits = c(-L_bivar, L_bivar)
)

ggsave(file.path(out_plots, "05a_chm_bivar_slope_wsci_FOREST.pdf"), 
       p_bivar_chm_sw_forest, width = 6, height = 5, bg = "white")
ggsave(file.path(out_plots, "05a_dtm_bivar_slope_wsci_FOREST.pdf"), 
       p_bivar_dtm_sw_forest, width = 6, height = 5, bg = "white")

# =====================================================================
# 10. LEVERAGE POINT IDENTIFICATION
# =====================================================================

log_subsection("Leverage point identification")

# Prepare data for leverage analysis
predictors <- c("slope_mean_z", "wsci_z", "rh_98_z", "cover_z")
predictors_present <- intersect(predictors, names(chm_df))

chm_complete <- chm_df %>% 
  select(chm_error_mean, all_of(predictors_present)) %>% 
  na.omit()

# Subsample for computational efficiency
set.seed(42)
if (nrow(chm_complete) > 100000) {
  chm_complete <- chm_complete %>% sample_n(100000)
}

log_progress(sprintf("  Using %s observations for leverage analysis", 
                     format(nrow(chm_complete), big.mark = ",")))

# Fit simple OLS to get hat values
ols_fit <- lm(chm_error_mean ~ ., data = chm_complete)
hat_vals <- hatvalues(ols_fit)
p <- length(coef(ols_fit))
n <- nrow(chm_complete)

# Identify high leverage points (hat > 2p/n)
high_leverage_threshold <- 2 * p / n
high_leverage <- hat_vals > high_leverage_threshold

log_progress(sprintf("  High leverage points (hat > 2p/n = %.4f): %d (%.2f%%)",
                     high_leverage_threshold,
                     sum(high_leverage), 
                     100 * mean(high_leverage)))

# Cook's distance for influential points
cooks_d <- cooks.distance(ols_fit)
influential_threshold <- 4 / n
influential <- cooks_d > influential_threshold

log_progress(sprintf("  Influential points (Cook's D > 4/n = %.4f): %d (%.2f%%)",
                     influential_threshold,
                     sum(influential), 
                     100 * mean(influential)))

# Very influential points
very_influential <- cooks_d > 1
if (sum(very_influential) > 0) {
  log_progress(sprintf("  ⚠ Very influential points (Cook's D > 1): %d", 
                       sum(very_influential)))
}

# Create leverage diagnostics dataframe
leverage_diag <- tibble::tibble(
  row_id = 1:n,
  hat_value = hat_vals,
  cooks_d = cooks_d,
  high_leverage = high_leverage,
  influential = influential
)

# Summary by characteristics of high-leverage points
high_lev_data <- chm_complete[high_leverage, ]
if (nrow(high_lev_data) > 0) {
  log_progress("  High-leverage point characteristics:")
  log_progress(sprintf("    Mean slope_z: %.2f (vs overall %.2f)", 
                       mean(high_lev_data$slope_mean_z), 
                       mean(chm_complete$slope_mean_z)))
  log_progress(sprintf("    Mean wsci_z: %.2f (vs overall %.2f)", 
                       mean(high_lev_data$wsci_z), 
                       mean(chm_complete$wsci_z)))
}

# Plot leverage diagnostics
p_leverage <- ggplot(leverage_diag, aes(x = hat_value, y = cooks_d)) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_hline(yintercept = influential_threshold, color = "orange", linetype = "dashed") +
  geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
  geom_vline(xintercept = high_leverage_threshold, color = "blue", linetype = "dashed") +
  scale_y_log10() +
  labs(title = "Leverage vs Influence Diagnostics",
       subtitle = sprintf("High leverage (blue): %d | Influential (orange): %d",
                          sum(high_leverage), sum(influential)),
       x = "Hat value (leverage)", y = "Cook's distance (influence, log scale)") +
  theme_cowplot()

ggsave(file.path(out_plots, "04d_leverage_influence.pdf"), p_leverage,
       width = 8, height = 6, bg = "white")

# Export leverage diagnostics summary
leverage_summary <- tibble::tibble(
  metric = c("n_observations", "n_predictors", "high_leverage_threshold", 
             "n_high_leverage", "pct_high_leverage",
             "influential_threshold", "n_influential", "pct_influential",
             "n_very_influential"),
  value = c(n, p, high_leverage_threshold, sum(high_leverage), 100*mean(high_leverage),
            influential_threshold, sum(influential), 100*mean(influential),
            sum(very_influential))
)
write_csv(leverage_summary, file.path(out_tables, "leverage_diagnostics.csv"))

# =====================================================================
# 11. PER-SITE BALANCE DIAGNOSTICS
# =====================================================================

log_subsection("Per-site balance diagnostics")

site_balance <- chm_df %>%
  group_by(site) %>%
  summarise(
    n = n(),
    n_lc_classes = n_distinct(lc_l1_code),
    slope_mean = mean(slope_mean, na.rm = TRUE),
    slope_sd = sd(slope_mean, na.rm = TRUE),
    slope_range = max(slope_mean, na.rm = TRUE) - min(slope_mean, na.rm = TRUE),
    pct_steep = 100 * mean(slope_mean > 30, na.rm = TRUE),
    wsci_mean = mean(wsci_z, na.rm = TRUE),
    cover_mean = mean(cover_z, na.rm = TRUE),
    error_mean = mean(chm_error_mean, na.rm = TRUE),
    error_sd = sd(chm_error_mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(n)

write_csv(site_balance, file.path(out_tables, "site_balance_diagnostics.csv"))

log_progress(sprintf("  Site sample sizes: min=%d, max=%d, CV=%.1f%%",
                     min(site_balance$n), max(site_balance$n),
                     100 * sd(site_balance$n) / mean(site_balance$n)))

log_progress(sprintf("  LC classes per site: min=%d, max=%d",
                     min(site_balance$n_lc_classes), max(site_balance$n_lc_classes)))

# Flag under-represented sites
small_sites <- site_balance %>% filter(n < 1000)
if (nrow(small_sites) > 0) {
  log_progress(sprintf("  ⚠ %d sites have <1000 observations", nrow(small_sites)))
}

# =====================================================================
# SUMMARY
# =====================================================================

log_progress("✓ Enhanced Core EDA complete")

# Save checkpoint for interactive mode
if (!exists("BATCH_MODE") || !BATCH_MODE) {
  log_progress("Saving checkpoint for interactive mode...")
  save_checkpoint("03_eda", list(
    eda_complete = TRUE,
    n_bins_common = if(exists("n_bins_common")) n_bins_common else NULL,
    error_stats = error_stats,
    error_stats_forest = error_stats_forest,
    site_balance = site_balance,
    leverage_summary = leverage_summary,
    forest_classes = FOREST_CLASSES,
    forest_present_chm = forest_present_chm
  ))
}
