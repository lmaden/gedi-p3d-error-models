# =====================================================================
# publication_figures.R
# Publication-ready figures for manuscript
#
# Creates three key figures identified as missing/important:
#   1. Predicted vs Observed scatterplot (with posterior uncertainty)
#   2. Conditional effects plots for key predictors
#   3. Spatial residual maps across study sites
#
# Follows conventions from:
#   - Gabry et al. (2019) "Visualization in Bayesian workflow"
#   - Remote Sensing of Environment publication standards
# =====================================================================

# =====================================================================
# CRITICAL: Force non-interactive PDF output for headless servers
# This must be set BEFORE loading any graphics packages
# =====================================================================
Sys.setenv(DISPLAY = "")
options(device = pdf)
options(bitmapType = "Xlib")  # Fallback, won't be used with PDF device

# Prevent any X11 device attempts
if (capabilities("X11")) {
  grDevices::X11.options(type = "Xlib")
}

# Ensure dependencies loaded
if (!exists("log_progress")) source("analysis_utils.R")
if (!exists("PROJECT_ROOT")) source("analysis_config.R")

library(brms)
library(bayesplot)
library(tidybayes)
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)
library(sf)
library(patchwork)
library(scales)
library(viridis)

# Disable interactive graphics in bayesplot
bayesplot::bayesplot_theme_set(bayesplot::theme_default())

log_progress("Creating publication-ready figures...")

# =====================================================================
# DATA AND MODEL LOADING
# =====================================================================

log_subsection("Loading data and models")

# Load data
if (!exists("chm_df") || !exists("dtm_df")) {
 log_progress("  Loading data from checkpoint...")
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

# Load models
model_chm <- NULL
model_dtm <- NULL
model_source <- "none"

if (checkpoint_exists("10_models_stage2")) {
 log_progress("  Loading Stage 2 models...")
 model_data <- load_checkpoint("10_models_stage2")
 model_chm <- model_data$fit_chm_s2
 model_dtm <- model_data$fit_dtm_s2
 model_source <- "Stage 2"
} else if (checkpoint_exists("09_models_stage1")) {
 log_progress("  Loading Stage 1 models...")
 model_data <- load_checkpoint("09_models_stage1")
 model_chm <- model_data$fit_chm_s1
 model_dtm <- model_data$fit_dtm_s1
 model_source <- "Stage 1"
} else {
 stop("No models available. Run modeling sections first.")
}

log_progress(sprintf("  Models loaded: %s", model_source))

# Define forest classes and publication color palette
FOREST_CLASSES <- c("BDF", "DNF", "EBF", "ENF")
FOREST_LABELS <- c(
 "BDF" = "Broadleaf Deciduous",
 "DNF" = "Deciduous Needleleaf", 
 "EBF" = "Evergreen Broadleaf",
 "ENF" = "Evergreen Needleleaf"
)

# Publication-quality color palette (colorblind-friendly)
FOREST_COLORS <- c(
 "BDF" = "#E69F00",  # Orange
 "DNF" = "#56B4E9",  # Sky blue
 "EBF" = "#009E73",  # Bluish green
 "ENF" = "#CC79A7"   # Reddish purple
)

# =====================================================================
# PUBLICATION THEME
# =====================================================================

theme_publication <- function(base_size = 11) {
 theme_cowplot(font_size = base_size) +
   theme(
     # Clean white background
     plot.background = element_rect(fill = "white", color = NA),
     panel.background = element_rect(fill = "white", color = NA),
     
     # Refined axis styling
     axis.line = element_line(color = "black", linewidth = 0.5),
     axis.ticks = element_line(color = "black", linewidth = 0.4),
     axis.text = element_text(color = "black", size = base_size - 1),
     axis.title = element_text(color = "black", size = base_size, face = "plain"),
     
     # Legend styling
     legend.background = element_rect(fill = "white", color = NA),
     legend.key = element_rect(fill = "white", color = NA),
     legend.text = element_text(size = base_size - 1),
     legend.title = element_text(size = base_size, face = "plain"),
     
     # Title styling
     plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0),
     plot.subtitle = element_text(size = base_size, hjust = 0, color = "gray30"),
     
     # Panel spacing for facets
     panel.spacing = unit(0.8, "lines"),
     strip.background = element_rect(fill = "gray95", color = NA),
     strip.text = element_text(size = base_size, face = "plain")
   )
}

# =====================================================================
# FIGURE 1: PREDICTED VS OBSERVED SCATTERPLOT
# =====================================================================

log_subsection("Creating Predicted vs Observed plots")

create_pred_vs_obs_plot <- function(model, df, error_col, product_name, 
                                    n_sample = 50000, n_draws = 100) {
 # Subsample for computational efficiency
 set.seed(42)
 
 # Get model data to ensure compatibility
 model_data <- model$data
 
 # Identify which observations can be predicted
 # (must have all predictor values matching model training data factor levels)
 df_subset <- df %>%
   filter(is.finite(.data[[error_col]])) %>%
   sample_n(min(n_sample, n()))
 
 # Get factor variables from model
 factor_vars <- names(model_data)[sapply(model_data, is.factor)]
 
 # Filter to observations with compatible factor levels
 for (fvar in factor_vars) {
   if (fvar %in% names(df_subset)) {
     model_levels <- levels(model_data[[fvar]])
     df_subset <- df_subset %>%
       filter(as.character(.data[[fvar]]) %in% model_levels)
   }
 }
 
 # Ensure factor levels match exactly
 for (fvar in factor_vars) {
   if (fvar %in% names(df_subset)) {
     df_subset[[fvar]] <- factor(df_subset[[fvar]], 
                                  levels = levels(model_data[[fvar]]))
   }
 }
 
 log_progress(sprintf("    Predicting for %d observations...", nrow(df_subset)))
 
 # Get posterior predictions (expected values)
 pred_summary <- tryCatch({
   # Use fitted() for expected values (posterior mean of mu)
   fitted_vals <- fitted(model, newdata = df_subset, 
                         summary = TRUE, ndraws = n_draws)
   
   tibble(
     observed = df_subset[[error_col]],
     predicted = fitted_vals[, "Estimate"],
     pred_lower = fitted_vals[, "Q2.5"],
     pred_upper = fitted_vals[, "Q97.5"],
     site = df_subset$site,
     lc_l1_code = df_subset$lc_l1_code
   )
 }, error = function(e) {
   log_progress(sprintf("    ⚠ Prediction failed: %s", e$message))
   return(NULL)
 })
 
 if (is.null(pred_summary)) return(NULL)
 
 # Calculate accuracy metrics
 rmse <- sqrt(mean((pred_summary$predicted - pred_summary$observed)^2))
 mae <- mean(abs(pred_summary$predicted - pred_summary$observed))
 r2 <- cor(pred_summary$predicted, pred_summary$observed)^2
 bias <- mean(pred_summary$predicted - pred_summary$observed)
 
 # Determine axis limits (symmetric around zero, capture 99% of data)
 lim_val <- quantile(abs(c(pred_summary$observed, pred_summary$predicted)), 
                     0.995, na.rm = TRUE)
 lim_val <- ceiling(lim_val)
 
 # Create main plot
 p <- ggplot(pred_summary, aes(x = observed, y = predicted)) +
   # Hex bins for density
   geom_hex(bins = 80, aes(fill = after_stat(count))) +
   scale_fill_viridis_c(
     trans = "log10",
     name = "Count",
     breaks = scales::trans_breaks("log10", function(x) 10^x),
     labels = scales::label_comma()
   ) +
   # 1:1 reference line
   geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
               color = "red", linewidth = 0.8) +
   # Zero reference lines
   geom_hline(yintercept = 0, linetype = "dotted", color = "gray50", linewidth = 0.4) +
   geom_vline(xintercept = 0, linetype = "dotted", color = "gray50", linewidth = 0.4) +
   # Axis settings
   coord_fixed(xlim = c(-lim_val, lim_val), ylim = c(-lim_val, lim_val)) +
   # Labels
   labs(
     title = sprintf("%s: Predicted vs Observed Error", product_name),
     subtitle = sprintf("R² = %.3f | RMSE = %.2f m | MAE = %.2f m | Bias = %.2f m",
                        r2, rmse, mae, bias),
     x = "Observed error (m)",
     y = "Predicted error (m)"
   ) +
   theme_publication() +
   theme(
     legend.position = c(0.02, 0.98),
     legend.justification = c(0, 1),
     legend.background = element_rect(fill = alpha("white", 0.9), color = NA)
   )
 
 # Store metrics as attribute
 attr(p, "metrics") <- list(r2 = r2, rmse = rmse, mae = mae, bias = bias)
 attr(p, "data") <- pred_summary
 
 return(p)
}

# Create CHM predicted vs observed
p_pred_obs_chm <- create_pred_vs_obs_plot(
 model_chm, chm_df, "chm_error_mean", "CHM",
 n_sample = 50000, n_draws = 100
)

if (!is.null(p_pred_obs_chm)) {
 ggsave(file.path(out_plots, "pub_01_pred_vs_obs_chm.pdf"), 
        p_pred_obs_chm, width = 7, height = 7, bg = "white", device = pdf)
 log_progress("  ✓ CHM predicted vs observed plot saved")
}

# Create DTM predicted vs observed
p_pred_obs_dtm <- create_pred_vs_obs_plot(
 model_dtm, dtm_df, "dtm_error_mean", "DTM",
 n_sample = 50000, n_draws = 100
)

if (!is.null(p_pred_obs_dtm)) {
 ggsave(file.path(out_plots, "pub_01_pred_vs_obs_dtm.pdf"), 
        p_pred_obs_dtm, width = 7, height = 7, bg = "white", device = pdf)
 log_progress("  ✓ DTM predicted vs observed plot saved")
}

# Combined panel figure
if (!is.null(p_pred_obs_chm) && !is.null(p_pred_obs_dtm)) {
 p_pred_obs_combined <- (p_pred_obs_chm + labs(title = "A) CHM")) | 
                        (p_pred_obs_dtm + labs(title = "B) DTM"))
 
 ggsave(file.path(out_plots, "pub_01_pred_vs_obs_combined.pdf"), 
        p_pred_obs_combined, width = 14, height = 7, bg = "white", device = pdf)
 log_progress("  ✓ Combined predicted vs observed plot saved")
}

# Create stratified version by forest type
if (!is.null(p_pred_obs_chm)) {
 pred_data_chm <- attr(p_pred_obs_chm, "data")
 
 if ("lc_l1_code" %in% names(pred_data_chm)) {
   pred_forest <- pred_data_chm %>%
     filter(lc_l1_code %in% FOREST_CLASSES)
   
   if (nrow(pred_forest) > 1000) {
     # Calculate per-forest-type metrics
     forest_metrics <- pred_forest %>%
       group_by(lc_l1_code) %>%
       summarise(
         r2 = cor(predicted, observed)^2,
         rmse = sqrt(mean((predicted - observed)^2)),
         n = n(),
         .groups = "drop"
       ) %>%
       mutate(label = sprintf("R² = %.2f\nRMSE = %.1f m\nn = %s", 
                              r2, rmse, format(n, big.mark = ",")))
     
     lim_val <- quantile(abs(c(pred_forest$observed, pred_forest$predicted)), 
                         0.995, na.rm = TRUE)
     lim_val <- ceiling(lim_val)
     
     p_pred_obs_forest <- ggplot(pred_forest, aes(x = observed, y = predicted)) +
       geom_hex(bins = 50, aes(fill = after_stat(count))) +
       scale_fill_viridis_c(trans = "log10", name = "Count") +
       geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
                   color = "red", linewidth = 0.6) +
       geom_text(data = forest_metrics, 
                 aes(x = -lim_val * 0.9, y = lim_val * 0.85, label = label),
                 hjust = 0, vjust = 1, size = 3, inherit.aes = FALSE) +
       facet_wrap(~ lc_l1_code, ncol = 2, 
                  labeller = labeller(lc_l1_code = FOREST_LABELS)) +
       coord_fixed(xlim = c(-lim_val, lim_val), ylim = c(-lim_val, lim_val)) +
       labs(
         title = "CHM Predicted vs Observed Error by Forest Type",
         x = "Observed error (m)",
         y = "Predicted error (m)"
       ) +
       theme_publication() +
       theme(legend.position = "right")
     
     ggsave(file.path(out_plots, "pub_01_pred_vs_obs_by_forest.pdf"), 
            p_pred_obs_forest, width = 10, height = 10, bg = "white", device = pdf)
     log_progress("  ✓ Forest-stratified predicted vs observed plot saved")
   }
 }
}

# =====================================================================
# FIGURE 2: CONDITIONAL EFFECTS PLOTS
# =====================================================================

log_subsection("Creating Conditional Effects plots")

create_conditional_effects_panel <- function(model, effects, product_name,
                                             resolution = 100) {
 # Check which effects are available in the model
 model_vars <- names(model$data)
 effects_present <- character(0)
 
 for (eff in effects) {
   if (grepl(":", eff, fixed = TRUE)) {
     parts <- strsplit(eff, ":")[[1]]
     if (all(parts %in% model_vars)) {
       effects_present <- c(effects_present, eff)
     }
   } else {
     if (eff %in% model_vars) {
       effects_present <- c(effects_present, eff)
     }
   }
 }
 
 if (length(effects_present) == 0) {
   log_progress(sprintf("    ⚠ No effects available for %s", product_name))
   return(NULL)
 }
 
 log_progress(sprintf("    Computing conditional effects for: %s", 
                      paste(effects_present, collapse = ", ")))
 
 # Compute conditional effects
 ce <- tryCatch({
   conditional_effects(model, effects = effects_present, 
                       resolution = resolution, method = "fitted")
 }, error = function(e) {
   log_progress(sprintf("    ⚠ conditional_effects failed: %s", e$message))
   return(NULL)
 })
 
 if (is.null(ce)) return(NULL)
 
 # Create individual plots with consistent styling
 plot_list <- list()
 
 # Nice labels for predictors
 predictor_labels <- c(
   "slope_mean_z" = "Terrain Slope (z-score)",
   "wsci_z" = "Waveform Structural Complexity (z-score)",
   "rh_98_z" = "Canopy Height RH98 (z-score)",
   "cover_z" = "Canopy Cover (z-score)",
   "meta_offnad_z" = "Off-Nadir Angle (z-score)",
   "meta_sunel_z" = "Sun Elevation (z-score)",
   "meta_leafon_z" = "Leaf-On Fraction (z-score)",
   "meta_stereo_z" = "Stereo Ratio (z-score)",
   "meta_az_conc_z" = "Azimuth Concentration (z-score)"
 )
 
 for (eff_name in names(ce)) {
   eff_data <- ce[[eff_name]]
   
   # Get the predictor variable name (first column that's not estimate__)
   pred_var <- setdiff(names(eff_data), 
                       c("estimate__", "se__", "lower__", "upper__", "cond__", "effect1__", "effect2__"))[1]
   
   x_label <- ifelse(pred_var %in% names(predictor_labels),
                     predictor_labels[pred_var], pred_var)
   
   p <- ggplot(eff_data, aes(x = .data[[pred_var]])) +
     # Credible interval ribbon
     geom_ribbon(aes(ymin = lower__, ymax = upper__), 
                 fill = "steelblue", alpha = 0.3) +
     # Mean line
     geom_line(aes(y = estimate__), color = "steelblue", linewidth = 1) +
     # Zero reference
     geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.4) +
     labs(
       x = x_label,
       y = sprintf("%s Error (m)", product_name)
     ) +
     theme_publication(base_size = 10)
   
   plot_list[[eff_name]] <- p
 }
 
 return(plot_list)
}

# Define key effects to visualize
# Terrain and vegetation structure effects (most interpretable)
key_effects <- c("slope_mean_z", "wsci_z", "rh_98_z", "cover_z")

# Check for available metadata effects
meta_effects <- intersect(
 c("meta_offnad_z", "meta_sunel_z", "meta_leafon_z", "meta_stereo_z"),
 names(model_chm$data)
)

# Create CHM conditional effects
ce_plots_chm <- create_conditional_effects_panel(
 model_chm, c(key_effects, meta_effects), "CHM"
)

if (!is.null(ce_plots_chm) && length(ce_plots_chm) > 0) {
 # Combine into panel figure
 n_plots <- length(ce_plots_chm)
 ncol <- min(3, n_plots)
 nrow <- ceiling(n_plots / ncol)
 
 p_ce_chm <- wrap_plots(ce_plots_chm, ncol = ncol) +
   plot_annotation(
     title = "CHM Error: Conditional Effects of Key Predictors",
     subtitle = "Shaded regions show 95% credible intervals",
     theme = theme_publication()
   )
 
 ggsave(file.path(out_plots, "pub_02_conditional_effects_chm.pdf"), 
        p_ce_chm, width = 4 * ncol, height = 3.5 * nrow, bg = "white", device = pdf)
 log_progress("  ✓ CHM conditional effects plot saved")
}

# Create DTM conditional effects
ce_plots_dtm <- create_conditional_effects_panel(
 model_dtm, c(key_effects, meta_effects), "DTM"
)

if (!is.null(ce_plots_dtm) && length(ce_plots_dtm) > 0) {
 n_plots <- length(ce_plots_dtm)
 ncol <- min(3, n_plots)
 nrow <- ceiling(n_plots / ncol)
 
 p_ce_dtm <- wrap_plots(ce_plots_dtm, ncol = ncol) +
   plot_annotation(
     title = "DTM Error: Conditional Effects of Key Predictors",
     subtitle = "Shaded regions show 95% credible intervals",
     theme = theme_publication()
   )
 
 ggsave(file.path(out_plots, "pub_02_conditional_effects_dtm.pdf"), 
        p_ce_dtm, width = 4 * ncol, height = 3.5 * nrow, bg = "white", device = pdf)
 log_progress("  ✓ DTM conditional effects plot saved")
}

# Create comparative panel (CHM vs DTM for key predictors)
if (!is.null(ce_plots_chm) && !is.null(ce_plots_dtm)) {
 # Find common effects
 common_effects <- intersect(names(ce_plots_chm), names(ce_plots_dtm))
 
 if (length(common_effects) >= 2) {
   # Select top 4 most important predictors for comparison
   compare_effects <- head(common_effects, 4)
   
   comparison_plots <- list()
   for (eff in compare_effects) {
     # Add product labels to titles
     p_chm <- ce_plots_chm[[eff]] + 
       labs(title = "CHM") +
       theme(plot.title = element_text(size = 10, hjust = 0.5))
     p_dtm <- ce_plots_dtm[[eff]] + 
       labs(title = "DTM") +
       theme(plot.title = element_text(size = 10, hjust = 0.5))
     
     comparison_plots[[eff]] <- p_chm | p_dtm
   }
   
   p_ce_compare <- wrap_plots(comparison_plots, ncol = 1) +
     plot_annotation(
       title = "Conditional Effects: CHM vs DTM Comparison",
       subtitle = "Comparing how predictors influence error in each product",
       theme = theme_publication()
     )
   
   ggsave(file.path(out_plots, "pub_02_conditional_effects_comparison.pdf"), 
          p_ce_compare, width = 10, height = 3 * length(compare_effects), bg = "white", device = pdf)
   log_progress("  ✓ CHM vs DTM conditional effects comparison saved")
 }
}

# =====================================================================
# FIGURE 3: SPATIAL RESIDUAL MAPS
# =====================================================================

log_subsection("Creating Spatial Residual Maps")

# Check if coordinates are available
has_coords <- all(c("x", "y") %in% names(chm_df))

if (!has_coords) {
 log_progress("  ⚠ Coordinates not available, skipping spatial residual maps")
} else {
 
 # Target CRS for mapping (Albers Equal Area for CONUS)
 TARGET_CRS <- 5070
 SOURCE_CRS <- 4326
 
 create_spatial_residual_map <- function(model, df, error_col, product_name,
                                          n_per_site = 3000, n_draws = 50) {
   
   # Sample data stratified by site
   set.seed(42)
   df_sample <- df %>%
     filter(is.finite(.data[[error_col]]), is.finite(x), is.finite(y)) %>%
     group_by(site) %>%
     mutate(.rand = runif(dplyr::n())) %>%
     arrange(.rand, .by_group = TRUE) %>%
     slice_head(n = n_per_site) %>%
     select(-.rand) %>%
     ungroup()
   
   log_progress(sprintf("    Sampled %d observations across %d sites", 
                        nrow(df_sample), n_distinct(df_sample$site)))
   
   # Get factor variables and filter compatible observations
   model_data <- model$data
   factor_vars <- names(model_data)[sapply(model_data, is.factor)]
   
   for (fvar in factor_vars) {
     if (fvar %in% names(df_sample)) {
       model_levels <- levels(model_data[[fvar]])
       df_sample <- df_sample %>%
         filter(as.character(.data[[fvar]]) %in% model_levels)
       df_sample[[fvar]] <- factor(df_sample[[fvar]], 
                                    levels = levels(model_data[[fvar]]))
     }
   }
   
   # Compute predictions
   log_progress("    Computing posterior predictions...")
   pred_vals <- tryCatch({
     fitted(model, newdata = df_sample, summary = TRUE, ndraws = n_draws)
   }, error = function(e) {
     log_progress(sprintf("    ⚠ Prediction failed: %s", e$message))
     return(NULL)
   })
   
   if (is.null(pred_vals)) return(NULL)
   
   # Compute residuals
   df_sample$predicted <- pred_vals[, "Estimate"]
   df_sample$residual <- df_sample[[error_col]] - df_sample$predicted
   
   # Project coordinates
   log_progress("    Projecting coordinates...")
   df_sf <- st_as_sf(df_sample, coords = c("x", "y"), crs = SOURCE_CRS)
   df_sf <- st_transform(df_sf, crs = TARGET_CRS)
   
   coords_proj <- st_coordinates(df_sf)
   df_sample$x_proj <- coords_proj[, 1]
   df_sample$y_proj <- coords_proj[, 2]
   
   # Ensure df_sample is a clean tibble (sf/matrix operations can affect structure)
   df_sample <- as_tibble(df_sample)
   
   # Calculate site-level summaries for labeling
   site_summary <- df_sample %>%
     group_by(site) %>%
     summarise(
       x_center = mean(x_proj),
       y_center = mean(y_proj),
       mean_resid = mean(residual),
       sd_resid = sd(residual),
       n = n(),
       .groups = "drop"
     )
   
   # Determine residual color limits (symmetric)
   resid_lim <- quantile(abs(df_sample$residual), 0.98, na.rm = TRUE)
   resid_lim <- ceiling(resid_lim)
   
   # Select sites to map (top 6 by sample size for visibility)
   sites_to_map <- df_sample %>%
     count(site) %>%
     arrange(desc(n)) %>%
     slice_head(n = 6) %>%
     pull(site)
   
   # Main multi-site panel figure
   map_data <- df_sample %>% filter(site %in% sites_to_map)
   
   p_map <- ggplot(map_data, aes(x = x_proj/1000, y = y_proj/1000, color = residual)) +
     geom_point(size = 0.3, alpha = 0.6) +
     scale_color_gradient2(
       low = "#2166AC", mid = "white", high = "#B2182B",
       midpoint = 0, 
       limits = c(-resid_lim, resid_lim),
       oob = scales::squish,
       name = "Residual (m)"
     ) +
     facet_wrap(~ site, scales = "free", ncol = 3) +
     coord_equal() +
     labs(
       title = sprintf("%s Model Residuals: Spatial Distribution", product_name),
       subtitle = "Residual = Observed error − Predicted error | Projected coordinates (Albers Equal Area)",
       x = "Easting (km)",
       y = "Northing (km)"
     ) +
     theme_publication(base_size = 10) +
     theme(
       axis.text = element_text(size = 7),
       strip.text = element_text(size = 9),
       legend.position = "right"
     )
   
   # Store data as attribute for further analysis
   attr(p_map, "data") <- df_sample
   attr(p_map, "site_summary") <- site_summary
   
   return(p_map)
 }
 
 # Create CHM spatial residual map
 p_spatial_chm <- create_spatial_residual_map(
   model_chm, chm_df, "chm_error_mean", "CHM",
   n_per_site = 3000, n_draws = 50
 )
 
 if (!is.null(p_spatial_chm)) {
   ggsave(file.path(out_plots, "pub_03_spatial_residuals_chm.pdf"), 
          p_spatial_chm, width = 12, height = 10, bg = "white", device = pdf)
   log_progress("  ✓ CHM spatial residual map saved")
   
   # Also save site summary
   site_summary_chm <- attr(p_spatial_chm, "site_summary")
   if (!is.null(site_summary_chm)) {
     write_csv(site_summary_chm, 
               file.path(out_tables, "spatial_residual_summary_chm.csv"))
   }
 }
 
 # Create DTM spatial residual map
 p_spatial_dtm <- create_spatial_residual_map(
   model_dtm, dtm_df, "dtm_error_mean", "DTM",
   n_per_site = 3000, n_draws = 50
 )
 
 if (!is.null(p_spatial_dtm)) {
   ggsave(file.path(out_plots, "pub_03_spatial_residuals_dtm.pdf"), 
          p_spatial_dtm, width = 12, height = 10, bg = "white", device = pdf)
   log_progress("  ✓ DTM spatial residual map saved")
   
   site_summary_dtm <- attr(p_spatial_dtm, "site_summary")
   if (!is.null(site_summary_dtm)) {
     write_csv(site_summary_dtm, 
               file.path(out_tables, "spatial_residual_summary_dtm.csv"))
   }
 }
 
 # Create residual histogram by site (supplementary diagnostic)
 if (!is.null(p_spatial_chm)) {
   resid_data_chm <- attr(p_spatial_chm, "data")
   
   p_resid_hist <- ggplot(resid_data_chm, aes(x = residual)) +
     geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7, color = "white") +
     geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
     facet_wrap(~ site, scales = "free_y", ncol = 4) +
     labs(
       title = "CHM Model Residual Distributions by Site",
       x = "Residual (m)",
       y = "Count"
     ) +
     theme_publication(base_size = 9) +
     theme(strip.text = element_text(size = 8))
   
   ggsave(file.path(out_plots, "pub_03_residual_histograms_by_site.pdf"), 
          p_resid_hist, width = 14, height = 12, bg = "white", device = pdf)
   log_progress("  ✓ Residual histograms by site saved")
 }
}

# =====================================================================
# SUMMARY OUTPUT
# =====================================================================

log_subsection("Summary")

# Collect all generated files
generated_files <- list.files(out_plots, pattern = "^pub_", full.names = FALSE)

log_progress(sprintf("Generated %d publication figure files:", length(generated_files)))
for (f in generated_files) {
 log_progress(sprintf("  - %s", f))
}

# Create summary metrics table
if (exists("p_pred_obs_chm") && !is.null(p_pred_obs_chm) &&
   exists("p_pred_obs_dtm") && !is.null(p_pred_obs_dtm)) {
 
 metrics_chm <- attr(p_pred_obs_chm, "metrics")
 metrics_dtm <- attr(p_pred_obs_dtm, "metrics")
 
 pred_accuracy_summary <- tibble::tibble(
   Product = c("CHM", "DTM"),
   R2 = c(metrics_chm$r2, metrics_dtm$r2),
   RMSE_m = c(metrics_chm$rmse, metrics_dtm$rmse),
   MAE_m = c(metrics_chm$mae, metrics_dtm$mae),
   Bias_m = c(metrics_chm$bias, metrics_dtm$bias)
 )
 
 write_csv(pred_accuracy_summary, 
           file.path(out_tables, "publication_predictive_accuracy.csv"))
 
 log_progress("\nPredictive Accuracy Summary:")
 print(pred_accuracy_summary)
}

log_progress("\n✓ Publication figures complete!")
log_progress(sprintf("  Output directory: %s", out_plots))

# =====================================================================
# SAVE CHECKPOINT
# =====================================================================

if (!exists("BATCH_MODE") || !BATCH_MODE) {
 checkpoint_data <- list(
   pub_figures_complete = TRUE,
   generated_files = generated_files
 )
 
 if (exists("pred_accuracy_summary")) {
   checkpoint_data$pred_accuracy_summary <- pred_accuracy_summary
 }
 
 save_checkpoint("15_publication_figures", checkpoint_data)
 log_progress("Checkpoint saved: 15_publication_figures")
}
