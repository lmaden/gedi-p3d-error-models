# =====================================================================
# section_12_spatial_FIXED.R
# Spatial autocorrelation diagnostics with CRITICAL FIXES
#
# FIXES IN THIS VERSION:
#   1. PROPER RESIDUAL EXTRACTION from brms models using posterior_predict
#   2. CRS TRANSFORMATION to projected coordinates before variogram analysis
#   3. GEDI TRANSECT AWARENESS - accounts for along-track sampling bias
#   4. Improved variogram lag spacing based on GEDI shot spacing (~60m)
#   5. Distance-band Moran's I to examine scale-dependent autocorrelation
# =====================================================================

# Ensure dependencies loaded
if (!exists("log_progress")) source("analysis_utils.R")
if (!exists("PROJECT_ROOT")) source("analysis_config.R")

library(sf)
library(spdep)
library(gstat)
library(brms)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)

# =====================================================================
# CONFIGURATION
# =====================================================================

# Target CRS for spatial analysis (Albers Equal Area for CONUS)
# This ensures distance calculations are in METERS
TARGET_CRS <- 5070  # EPSG:5070 = NAD83 / Conus Albers

# Source CRS (WGS84 lat/lon from GEDI)
SOURCE_CRS <- 4326

# GEDI-specific parameters
GEDI_ALONG_TRACK_SPACING <- 60     # meters between shots along track
GEDI_SWATH_WIDTH <- 6000           # meters (~6 km swath)
MIN_LAG_DISTANCE <- 100            # meters - minimum variogram lag
MAX_LAG_DISTANCE <- 50000          # meters - maximum lag for variogram

# =====================================================================
# DATA LOADING
# =====================================================================

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

# =====================================================================
# LOAD MODELS - CRITICAL FOR PROPER RESIDUAL EXTRACTION
# =====================================================================

log_progress("Loading Bayesian models for residual extraction...")

model_chm_final <- NULL
model_dtm_final <- NULL
model_source <- "none"

# Try Stage 2 models first
if (checkpoint_exists("10_models_stage2")) {
  log_progress("  Loading Stage 2 models...")
  data <- load_checkpoint("10_models_stage2")
  model_chm_final <- data$fit_chm_s2
  model_dtm_final <- data$fit_dtm_s2
  model_source <- "stage2"
  log_progress("  ✓ Stage 2 models loaded")
} else if (checkpoint_exists("09_models_stage1")) {
  log_progress("  Loading Stage 1 models...")
  data <- load_checkpoint("09_models_stage1")
  model_chm_final <- data$fit_chm_s1
  model_dtm_final <- data$fit_dtm_s1
  model_source <- "stage1"
  log_progress("  ✓ Stage 1 models loaded")
} else {
  log_progress("  ⚠ No models found - will use raw errors instead of residuals")
  log_progress("    NOTE: Spatial analysis on raw errors shows data structure,")
  log_progress("          not unexplained model variance. Results should be interpreted accordingly.")
}

# =====================================================================
# CHECK COORDINATES
# =====================================================================

do_spatial_diag <- exists("attach_coords") && attach_coords && 
  all(c("x","y") %in% names(chm_df)) && 
  all(c("x","y") %in% names(dtm_df))

if (!do_spatial_diag) {
  log_progress("⊘ Coordinates not available, skipping spatial diagnostics")
  if (!exists("BATCH_MODE") || !BATCH_MODE) {
    save_checkpoint("12_spatial", list(spatial_complete = FALSE))
  }
  return(invisible(NULL))
}

log_progress("Running spatial autocorrelation diagnostics (FIXED VERSION)...")
log_progress(sprintf("  Model source: %s", model_source))

# =====================================================================
# 1. SAMPLE FOR SPATIAL ANALYSIS
# =====================================================================

log_subsection("Preparing spatial samples")

set.seed(42)
chm_s <- chm_df %>% 
  tidyr::drop_na(x, y, chm_error_mean) %>% 
  group_by(site) %>%
  group_modify(~ slice_sample(.x, n = min(5000, nrow(.x)), replace = FALSE)) %>% 
  ungroup()

dtm_s <- dtm_df %>% 
  tidyr::drop_na(x, y, dtm_error_mean) %>% 
  group_by(site) %>%
  group_modify(~ slice_sample(.x, n = min(5000, nrow(.x)), replace = FALSE)) %>% 
  ungroup()

log_progress(sprintf("  Spatial sample CHM: %s", format(nrow(chm_s), big.mark=",")))
log_progress(sprintf("  Spatial sample DTM: %s", format(nrow(dtm_s), big.mark=",")))

# =====================================================================
# 2. PROPER RESIDUAL EXTRACTION FROM BRMS MODELS
# =====================================================================

log_subsection("Extracting model residuals (FIXED)")

extract_brms_residuals <- function(model, newdata, response_col, name) {
  # This function properly extracts residuals from a brms model
  # by computing: observed - E[predicted]
  #
  # CRITICAL: Handles factor level mismatches by filtering to compatible observations
  
  if (is.null(model)) {
    log_progress(sprintf("    %s: No model available, using raw errors", name))
    return(list(residuals = newdata[[response_col]], valid_idx = seq_len(nrow(newdata))))
  }
  
  tryCatch({
    log_progress(sprintf("    %s: Computing posterior predictions...", name))
    
    # Get the model's training data to check factor levels
    model_data <- model$data
    
    # Identify factor variables in the model
    factor_vars <- names(model_data)[sapply(model_data, is.factor)]
    log_progress(sprintf("      Factor variables in model: %s", 
                         paste(factor_vars, collapse = ", ")))
    
    # For each factor variable, check which observations in newdata have valid levels
    valid_idx <- rep(TRUE, nrow(newdata))
    
    for (fvar in factor_vars) {
      if (fvar %in% names(newdata)) {
        model_levels <- levels(model_data[[fvar]])
        newdata_values <- as.character(newdata[[fvar]])
        
        # Find observations with levels not in model
        invalid_levels <- setdiff(unique(newdata_values), model_levels)
        
        if (length(invalid_levels) > 0) {
          log_progress(sprintf("      ⚠ %s has new levels not in model: %s", 
                               fvar, paste(invalid_levels, collapse = ", ")))
          
          # Mark these observations as invalid for prediction
          valid_idx <- valid_idx & (newdata_values %in% model_levels)
        }
      }
    }
    
    n_valid <- sum(valid_idx)
    n_total <- nrow(newdata)
    n_excluded <- n_total - n_valid
    
    if (n_excluded > 0) {
      log_progress(sprintf("      Excluding %d/%d (%.1f%%) observations with incompatible factor levels",
                           n_excluded, n_total, 100 * n_excluded / n_total))
    }
    
    if (n_valid < 100) {
      log_progress("      ⚠ Too few valid observations for prediction, using raw errors")
      return(list(residuals = newdata[[response_col]], valid_idx = seq_len(nrow(newdata))))
    }
    
    # Subset to valid observations
    newdata_valid <- newdata[valid_idx, ]
    
    # Ensure factor levels match exactly
    for (fvar in factor_vars) {
      if (fvar %in% names(newdata_valid)) {
        newdata_valid[[fvar]] <- factor(as.character(newdata_valid[[fvar]]), 
                                         levels = levels(model_data[[fvar]]))
      }
    }
    
    # Compute posterior mean predictions on valid subset
    log_progress(sprintf("      Computing predictions for %d observations...", n_valid))
    pred <- fitted(model, newdata = newdata_valid, summary = TRUE, 
                   re_formula = NULL,  # Include all random effects
                   allow_new_levels = TRUE,
                   sample_new_levels = "uncertainty")
    
    # Residuals = observed - predicted (posterior mean)
    observed_valid <- newdata_valid[[response_col]]
    predicted <- pred[, "Estimate"]
    residuals_valid <- observed_valid - predicted
    
    # Create full residual vector (NA for excluded observations)
    residuals_full <- rep(NA_real_, n_total)
    residuals_full[valid_idx] <- residuals_valid
    
    # Report summary
    log_progress(sprintf("      ✓ Residuals computed: mean=%.3f, sd=%.3f", 
                         mean(residuals_valid, na.rm = TRUE), 
                         sd(residuals_valid, na.rm = TRUE)))
    log_progress(sprintf("      Raw error sd=%.3f vs Residual sd=%.3f (%.1f%% variance reduction)",
                         sd(observed_valid, na.rm = TRUE),
                         sd(residuals_valid, na.rm = TRUE),
                         100 * (1 - var(residuals_valid, na.rm = TRUE) / var(observed_valid, na.rm = TRUE))))
    
    return(list(residuals = residuals_full, valid_idx = which(valid_idx)))
    
  }, error = function(e) {
    log_progress(sprintf("      ⚠ Residual extraction failed: %s", e$message))
    log_progress("      Falling back to raw errors")
    return(list(residuals = newdata[[response_col]], valid_idx = seq_len(nrow(newdata))))
  })
}

# Extract residuals for CHM
chm_result <- extract_brms_residuals(
  model_chm_final, 
  chm_s, 
  "chm_error_mean",
  "CHM"
)
res_chm_s <- chm_result$residuals

# Extract residuals for DTM
dtm_result <- extract_brms_residuals(
  model_dtm_final, 
  dtm_s, 
  "dtm_error_mean",
  "DTM"
)
res_dtm_s <- dtm_result$residuals

# Add residuals to data
chm_s$residual <- res_chm_s
dtm_s$residual <- res_dtm_s

# Filter to only observations with valid residuals for spatial analysis
chm_s <- chm_s %>% filter(!is.na(residual))
dtm_s <- dtm_s %>% filter(!is.na(residual))

log_progress(sprintf("  After factor filtering - CHM: %s, DTM: %s observations",
                     format(nrow(chm_s), big.mark = ","),
                     format(nrow(dtm_s), big.mark = ",")))

# =====================================================================
# 3. CRS TRANSFORMATION - CRITICAL FOR PROPER DISTANCE CALCULATIONS
# =====================================================================

log_subsection("Transforming to projected CRS (FIXED)")

log_progress(sprintf("  Source CRS: EPSG:%d (WGS84 lat/lon)", SOURCE_CRS))
log_progress(sprintf("  Target CRS: EPSG:%d (Albers Equal Area)", TARGET_CRS))

# Convert to sf, transform, extract projected coordinates
chm_sf <- st_as_sf(chm_s, coords = c("x", "y"), crs = SOURCE_CRS) %>%
  st_transform(crs = TARGET_CRS)

dtm_sf <- st_as_sf(dtm_s, coords = c("x", "y"), crs = SOURCE_CRS) %>%
  st_transform(crs = TARGET_CRS)

# Extract projected coordinates
chm_coords <- st_coordinates(chm_sf)
chm_s$x_proj <- chm_coords[, 1]
chm_s$y_proj <- chm_coords[, 2]

dtm_coords <- st_coordinates(dtm_sf)
dtm_s$x_proj <- dtm_coords[, 1]
dtm_s$y_proj <- dtm_coords[, 2]

log_progress("  ✓ Coordinates transformed to meters")

# Calculate approximate spatial extent
x_range <- diff(range(chm_s$x_proj, na.rm = TRUE))
y_range <- diff(range(chm_s$y_proj, na.rm = TRUE))
log_progress(sprintf("  CHM spatial extent: %.1f km × %.1f km", 
                     x_range / 1000, y_range / 1000))

# =====================================================================
# 4. MORAN'S I WITH DISTANCE BANDS (accounts for GEDI sampling)
# =====================================================================

log_subsection("Moran's I test (with distance bands)")

# Function for distance-banded Moran's I
moran_distance_bands <- function(xy_df, res_vec, 
                                  bands = c(100, 500, 1000, 5000, 10000, 25000)) {
  coords <- as.matrix(xy_df[, c("x_proj", "y_proj")])
  results <- list()
  
  for (i in 1:(length(bands) - 1)) {
    d_lower <- bands[i]
    d_upper <- bands[i + 1]
    
    tryCatch({
      # Create distance-based neighbors within band
      nb <- spdep::dnearneigh(coords, d1 = d_lower, d2 = d_upper)
      
      # Check if we have any neighbors
      n_links <- sum(card(nb))
      if (n_links < 10) {
        results[[paste0(d_lower, "-", d_upper)]] <- list(
          band = paste0(d_lower/1000, "-", d_upper/1000, " km"),
          moran_i = NA,
          p_value = NA,
          n_links = n_links,
          note = "Too few neighbor links"
        )
        next
      }
      
      lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
      moran_result <- spdep::moran.test(res_vec, lw, zero.policy = TRUE)
      
      results[[paste0(d_lower, "-", d_upper)]] <- list(
        band = paste0(d_lower/1000, "-", d_upper/1000, " km"),
        moran_i = moran_result$estimate[1],
        p_value = moran_result$p.value,
        n_links = n_links
      )
    }, error = function(e) {
      results[[paste0(d_lower, "-", d_upper)]] <- list(
        band = paste0(d_lower/1000, "-", d_upper/1000, " km"),
        moran_i = NA,
        p_value = NA,
        note = e$message
      )
    })
  }
  
  return(results)
}

# Also compute traditional k-NN Moran's I for comparison
moran_calc <- function(xy_df, res_vec, k = 8) {
  coords <- as.matrix(xy_df[, c("x_proj", "y_proj")])
  nb <- spdep::knn2nb(spdep::knearneigh(coords, k = k))
  lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
  spdep::moran.test(res_vec, lw, zero.policy = TRUE)
}

# k-NN Moran's I
log_progress("  Computing k-NN Moran's I (k=8)...")

moran_chm <- tryCatch({
  moran_calc(chm_s, chm_s$residual)
}, error = function(e) {
  log_progress(sprintf("  ⚠ CHM Moran's I failed: %s", e$message))
  NULL
})

moran_dtm <- tryCatch({
  moran_calc(dtm_s, dtm_s$residual)
}, error = function(e) {
  log_progress(sprintf("  ⚠ DTM Moran's I failed: %s", e$message))
  NULL
})

if (!is.null(moran_chm)) {
  capture.output(moran_chm, file = file.path(out_tables, "moran_chm.txt"))
  log_progress(sprintf("  CHM Moran's I (k=8 NN): %.4f (p=%.4g)", 
                       moran_chm$estimate[1], moran_chm$p.value))
}

if (!is.null(moran_dtm)) {
  capture.output(moran_dtm, file = file.path(out_tables, "moran_dtm.txt"))
  log_progress(sprintf("  DTM Moran's I (k=8 NN): %.4f (p=%.4g)", 
                       moran_dtm$estimate[1], moran_dtm$p.value))
}

# Distance-banded Moran's I
log_progress("  Computing distance-banded Moran's I...")

# Define bands relevant to GEDI sampling:
# 0-100m: Within-track (adjacent shots)
# 100-500m: Along-track (same swath)
# 500-2000m: Cross-track within swath
# 2000-10000m: Between adjacent tracks
# 10000-50000m: Regional scale
distance_bands <- c(0, 100, 500, 2000, 10000, 50000)

moran_bands_chm <- moran_distance_bands(chm_s, chm_s$residual, bands = distance_bands)
moran_bands_dtm <- moran_distance_bands(dtm_s, dtm_s$residual, bands = distance_bands)

# Report distance-banded results
log_progress("  CHM Moran's I by distance band:")
for (band_result in moran_bands_chm) {
  if (!is.na(band_result$moran_i)) {
    log_progress(sprintf("    %s: I=%.4f (p=%.4g, n_links=%d)", 
                         band_result$band, band_result$moran_i, 
                         band_result$p_value, band_result$n_links))
  } else {
    log_progress(sprintf("    %s: %s", band_result$band, 
                         ifelse(is.null(band_result$note), "NA", band_result$note)))
  }
}

# Export distance-band results
moran_bands_df <- tryCatch({
  chm_rows <- lapply(names(moran_bands_chm), function(nm) {
    x <- moran_bands_chm[[nm]]
    data.frame(
      product = "CHM", 
      band = as.character(x$band), 
      moran_i = as.numeric(x$moran_i), 
      p_value = as.numeric(x$p_value), 
      n_links = as.numeric(ifelse(is.null(x$n_links), NA, x$n_links)),
      stringsAsFactors = FALSE
    )
  })
  
  dtm_rows <- lapply(names(moran_bands_dtm), function(nm) {
    x <- moran_bands_dtm[[nm]]
    data.frame(
      product = "DTM", 
      band = as.character(x$band), 
      moran_i = as.numeric(x$moran_i), 
      p_value = as.numeric(x$p_value), 
      n_links = as.numeric(ifelse(is.null(x$n_links), NA, x$n_links)),
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, c(chm_rows, dtm_rows))
}, error = function(e) {
  log_progress(sprintf("  ⚠ Could not create distance-band dataframe: %s", e$message))
  NULL
})

if (!is.null(moran_bands_df)) {
  write.csv(moran_bands_df, file.path(out_tables, "moran_by_distance_band.csv"), row.names = FALSE)
  log_progress("  ✓ Distance-banded Moran's I exported")
}

# =====================================================================
# 5. VARIOGRAM ANALYSIS WITH PROPER CRS AND LAG SPACING
# =====================================================================

log_subsection("Empirical variogram analysis (FIXED CRS)")

vg_calc_projected <- function(df_with_coords, res_col, max_dist = MAX_LAG_DISTANCE, 
                               n_lags = 20, min_lag = MIN_LAG_DISTANCE) {
  # Create sf object with projected coordinates
  df <- df_with_coords %>%
    mutate(res = .data[[res_col]]) %>%
    filter(is.finite(res))
  
  sf_obj <- st_as_sf(df, coords = c("x_proj", "y_proj"), crs = TARGET_CRS)
  
  # Create gstat object
  g <- gstat::gstat(id = "res", formula = res ~ 1, data = sf_obj)
  
  # Calculate variogram with specified lag spacing
  # Width = lag bin size, cutoff = maximum distance
  lag_width <- (max_dist - min_lag) / n_lags
  
  gstat::variogram(g, cutoff = max_dist, width = lag_width)
}

# Calculate empirical variograms with projected coordinates
log_progress("  Computing CHM variogram (projected CRS)...")
vg_chm <- tryCatch({
  vg_calc_projected(chm_s, "residual")
}, error = function(e) {
  log_progress(sprintf("  ⚠ CHM variogram failed: %s", e$message))
  NULL
})

log_progress("  Computing DTM variogram (projected CRS)...")
vg_dtm <- tryCatch({
  vg_calc_projected(dtm_s, "residual")
}, error = function(e) {
  log_progress(sprintf("  ⚠ DTM variogram failed: %s", e$message))
  NULL
})

if (!is.null(vg_chm)) {
  saveRDS(vg_chm, file.path(out_tables, "variogram_chm_projected.rds"))
  log_progress("  ✓ CHM empirical variogram computed (distances in meters)")
  log_progress(sprintf("    Distance range: %.0f m to %.0f m", 
                       min(vg_chm$dist), max(vg_chm$dist)))
}

if (!is.null(vg_dtm)) {
  saveRDS(vg_dtm, file.path(out_tables, "variogram_dtm_projected.rds"))
  log_progress("  ✓ DTM empirical variogram computed (distances in meters)")
}

# =====================================================================
# 6. SEMI-VARIOGRAM MODEL FITTING
# =====================================================================

log_subsection("Semi-variogram model fitting")

fit_variogram_models <- function(vg_empirical, name) {
  if (is.null(vg_empirical)) return(NULL)
  
  results <- list()
  
  # Get initial parameter estimates
  max_gamma <- max(vg_empirical$gamma)
  max_dist <- max(vg_empirical$dist)
  
  # Better initial estimates
  # Nugget: semivariance at shortest lag
  init_nugget <- vg_empirical$gamma[1] * 0.8
  # Sill: maximum semivariance
  init_sill <- max_gamma - init_nugget
  # Range: distance where gamma reaches ~63% of sill (for exponential)
  # or ~95% of sill (for spherical/Gaussian)
  gamma_target <- init_nugget + 0.63 * init_sill
  init_range <- vg_empirical$dist[which.min(abs(vg_empirical$gamma - gamma_target))]
  if (init_range < MIN_LAG_DISTANCE) init_range <- max_dist / 4
  
  log_progress(sprintf("    Initial estimates: nugget=%.2f, psill=%.2f, range=%.0f m",
                       init_nugget, init_sill, init_range))
  
  # Try different models
  model_types <- c("Exp", "Sph", "Gau", "Mat")  # Added Matern
  
  for (model_type in model_types) {
    tryCatch({
      # Initial variogram model
      if (model_type == "Mat") {
        vg_model <- gstat::vgm(
          psill = init_sill,
          model = model_type,
          range = init_range,
          nugget = init_nugget,
          kappa = 1.5  # Matern smoothness
        )
      } else {
        vg_model <- gstat::vgm(
          psill = init_sill,
          model = model_type,
          range = init_range,
          nugget = init_nugget
        )
      }
      
      # Fit model to empirical variogram
      fit <- gstat::fit.variogram(vg_empirical, vg_model)
      
      # Calculate sum of squared errors (weighted by number of pairs)
      predicted <- gstat::variogramLine(fit, maxdist = max_dist, n = nrow(vg_empirical))
      pred_at_emp <- approx(predicted$dist, predicted$gamma, vg_empirical$dist)$y
      
      # Weight by number of pairs (np)
      weights <- vg_empirical$np / sum(vg_empirical$np)
      sse <- sum(weights * (vg_empirical$gamma - pred_at_emp)^2, na.rm = TRUE)
      
      # Extract parameters
      nugget <- fit$psill[1]
      partial_sill <- fit$psill[2]
      total_sill <- sum(fit$psill)
      range_param <- fit$range[2]
      
      results[[model_type]] <- list(
        model = fit,
        sse = sse,
        nugget = nugget,
        partial_sill = partial_sill,
        sill = total_sill,
        range = range_param,
        nugget_sill_ratio = nugget / total_sill
      )
      
      log_progress(sprintf("    %s %s: nugget=%.2f, sill=%.2f, range=%.0f m, SSE=%.4f",
                           name, model_type, 
                           nugget, total_sill, range_param, sse))
      
    }, error = function(e) {
      log_progress(sprintf("    ⚠ %s %s fit failed: %s", name, model_type, e$message))
    })
  }
  
  # Select best model by SSE
  if (length(results) > 0) {
    sse_values <- sapply(results, function(x) x$sse)
    best_model <- names(which.min(sse_values))
    log_progress(sprintf("    Best model for %s: %s", name, best_model))
    
    results$best <- best_model
    results$best_fit <- results[[best_model]]
  }
  
  return(results)
}

log_progress("  Fitting variogram models to CHM...")
vg_fits_chm <- fit_variogram_models(vg_chm, "CHM")

log_progress("  Fitting variogram models to DTM...")
vg_fits_dtm <- fit_variogram_models(vg_dtm, "DTM")

# Export variogram fit parameters
if (!is.null(vg_fits_chm) && !is.null(vg_fits_chm$best_fit)) {
  vg_params_chm <- tibble::tibble(
    product = "CHM",
    best_model = vg_fits_chm$best,
    nugget = vg_fits_chm$best_fit$nugget,
    partial_sill = vg_fits_chm$best_fit$partial_sill,
    sill = vg_fits_chm$best_fit$sill,
    range_m = vg_fits_chm$best_fit$range,
    nugget_sill_ratio = vg_fits_chm$best_fit$nugget_sill_ratio,
    sse = vg_fits_chm$best_fit$sse
  )
  
  log_progress(sprintf("  CHM nugget/sill ratio: %.3f", vg_params_chm$nugget_sill_ratio))
  log_progress(sprintf("  CHM effective range: %.0f meters (%.1f km)", 
                       vg_params_chm$range_m, vg_params_chm$range_m / 1000))
  
  if (vg_params_chm$nugget_sill_ratio > 0.5) {
    log_progress("    → High nugget ratio: majority of variance is micro-scale/measurement error")
  } else if (vg_params_chm$nugget_sill_ratio < 0.25) {
    log_progress("    → Low nugget ratio: strong spatial structure in residuals")
  }
} else {
  vg_params_chm <- NULL
}

if (!is.null(vg_fits_dtm) && !is.null(vg_fits_dtm$best_fit)) {
  vg_params_dtm <- tibble::tibble(
    product = "DTM",
    best_model = vg_fits_dtm$best,
    nugget = vg_fits_dtm$best_fit$nugget,
    partial_sill = vg_fits_dtm$best_fit$partial_sill,
    sill = vg_fits_dtm$best_fit$sill,
    range_m = vg_fits_dtm$best_fit$range,
    nugget_sill_ratio = vg_fits_dtm$best_fit$nugget_sill_ratio,
    sse = vg_fits_dtm$best_fit$sse
  )
} else {
  vg_params_dtm <- NULL
}

vg_params <- bind_rows(vg_params_chm, vg_params_dtm)
if (nrow(vg_params) > 0) {
  write_csv(vg_params, file.path(out_tables, "variogram_parameters_projected.csv"))
  log_progress("  ✓ Variogram parameters exported (distances in meters)")
}

# =====================================================================
# 7. VARIOGRAM PLOT WITH FITTED MODEL
# =====================================================================

if (!is.null(vg_chm) && !is.null(vg_fits_chm$best_fit)) {
  vg_plot_data <- vg_chm
  vg_plot_data$type <- "Empirical"
  
  # Add fitted line
  fitted_line <- gstat::variogramLine(vg_fits_chm$best_fit$model, 
                                       maxdist = max(vg_chm$dist), n = 100)
  fitted_line$type <- paste0("Fitted (", vg_fits_chm$best, ")")
  fitted_line$np <- NA
  
  # Add GEDI-relevant distance markers
  p_vg_chm <- ggplot() +
    geom_point(data = vg_plot_data, aes(x = dist/1000, y = gamma, size = np), alpha = 0.7) +
    geom_line(data = fitted_line, aes(x = dist/1000, y = gamma), color = "red", linewidth = 1) +
    geom_hline(yintercept = vg_fits_chm$best_fit$sill, linetype = "dashed", color = "blue", alpha = 0.7) +
    geom_vline(xintercept = vg_fits_chm$best_fit$range/1000, linetype = "dashed", color = "forestgreen", alpha = 0.7) +
    # Add GEDI swath width reference
    geom_vline(xintercept = GEDI_SWATH_WIDTH/1000, linetype = "dotted", color = "orange", alpha = 0.5) +
    annotate("text", x = GEDI_SWATH_WIDTH/1000 + 1, y = max(vg_plot_data$gamma) * 0.95, 
             label = "GEDI swath", color = "orange", hjust = 0, size = 3) +
    labs(title = "CHM: Semi-variogram with Fitted Model",
         subtitle = sprintf("%s model: nugget=%.2f m², sill=%.2f m², range=%.1f km\nNugget/Sill ratio=%.3f",
                            vg_fits_chm$best, vg_fits_chm$best_fit$nugget,
                            vg_fits_chm$best_fit$sill, vg_fits_chm$best_fit$range/1000,
                            vg_fits_chm$best_fit$nugget_sill_ratio),
         x = "Distance (km)", y = "Semivariance (m²)", size = "N pairs") +
    theme_cowplot() +
    theme(plot.subtitle = element_text(size = 9))
  
  ggsave(file.path(out_plots, "10_variogram_chm_projected.pdf"), p_vg_chm,
         width = 10, height = 7, bg = "white")
  log_progress("  ✓ CHM variogram plot saved (with proper distance units)")
}

# =====================================================================
# 8. SPATIAL RESIDUAL MAPS (with projected coordinates)
# =====================================================================

log_subsection("Spatial residual maps")

sites_to_map <- chm_s %>%
  group_by(site) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(desc(n)) %>%
  slice_head(n = 6) %>%
  pull(site)

if (length(sites_to_map) > 0) {
  map_data <- chm_s %>% filter(site %in% sites_to_map)
  
  # Use projected coordinates for proper aspect ratio
  p_spatial_resid <- ggplot(map_data, aes(x = x_proj/1000, y = y_proj/1000, color = residual)) +
    geom_point(size = 0.5, alpha = 0.5) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", 
                          midpoint = 0, limits = c(-10, 10), oob = scales::squish,
                          name = "Residual (m)") +
    facet_wrap(~site, scales = "free", ncol = 3) +
    labs(title = sprintf("Spatial Distribution of CHM %s by Site",
                         ifelse(model_source == "none", "Raw Errors", "Model Residuals")),
         subtitle = "Projected coordinates (Albers Equal Area) - units in km",
         x = "Easting (km)", y = "Northing (km)") +
    coord_equal() +  # Ensures proper aspect ratio
    theme_cowplot() +
    theme(axis.text = element_text(size = 6))
  
  ggsave(file.path(out_plots, "10_spatial_residual_maps_projected.pdf"), p_spatial_resid,
         width = 12, height = 10, bg = "white")
  log_progress("  ✓ Spatial residual maps saved (projected coordinates)")
}

# =====================================================================
# 9. PER-SITE SPATIAL AUTOCORRELATION
# =====================================================================

log_subsection("Per-site spatial autocorrelation")

moran_by_site <- chm_s %>%
  group_by(site) %>%
  group_modify(function(site_data, keys) {
    if (nrow(site_data) < 100) {
      return(tibble::tibble(moran_i = NA_real_, p_value = NA_real_, 
                            n = nrow(site_data), mean_nn_dist_m = NA_real_))
    }
    
    tryCatch({
      coords <- as.matrix(site_data[, c("x_proj", "y_proj")])
      
      # Calculate mean nearest neighbor distance for context
      nn <- spdep::knearneigh(coords, k = 1)
      nn_dists <- sqrt(rowSums((coords - coords[nn$nn[,1], ])^2))
      mean_nn_dist <- mean(nn_dists)
      
      moran_result <- moran_calc(site_data, site_data$residual, k = 8)
      tibble::tibble(
        moran_i = moran_result$estimate[1],
        p_value = moran_result$p.value,
        n = nrow(site_data),
        mean_nn_dist_m = mean_nn_dist
      )
    }, error = function(e) {
      tibble::tibble(moran_i = NA_real_, p_value = NA_real_, 
                     n = nrow(site_data), mean_nn_dist_m = NA_real_)
    })
  }) %>%
  ungroup()

write_csv(moran_by_site, file.path(out_tables, "moran_by_site_projected.csv"))

# Summarize
sig_autocorr <- moran_by_site %>% 
  filter(!is.na(moran_i), p_value < 0.05, moran_i > 0.1)

log_progress(sprintf("  Sites with significant positive autocorrelation: %d / %d",
                     nrow(sig_autocorr), 
                     sum(!is.na(moran_by_site$moran_i))))

log_progress(sprintf("  Mean nearest-neighbor distance across sites: %.0f m",
                     mean(moran_by_site$mean_nn_dist_m, na.rm = TRUE)))

if (nrow(sig_autocorr) > 0) {
  log_progress("  Top 5 sites by Moran's I:")
  top_sites <- sig_autocorr %>% arrange(desc(moran_i)) %>% slice(1:5)
  for (j in 1:nrow(top_sites)) {
    log_progress(sprintf("    Site %s: I=%.3f (p=%.4g, mean NN dist=%.0f m)", 
                         top_sites$site[j], top_sites$moran_i[j], 
                         top_sites$p_value[j], top_sites$mean_nn_dist_m[j]))
  }
}

# =====================================================================
# 10. INTERPRETATION SUMMARY
# =====================================================================

log_subsection("Interpretation Summary")

log_progress("=" %>% strrep(70))
log_progress("SPATIAL DIAGNOSTICS SUMMARY")
log_progress("=" %>% strrep(70))

log_progress(sprintf("\nAnalysis type: %s", 
                     ifelse(model_source == "none", 
                            "RAW ERRORS (model not available)", 
                            paste0("MODEL RESIDUALS (", model_source, ")"))))

if (!is.null(moran_chm)) {
  log_progress(sprintf("\nCHM Moran's I (k=8 NN): %.4f", moran_chm$estimate[1]))
  if (moran_chm$estimate[1] > 0.5) {
    log_progress("  INTERPRETATION: Very strong positive spatial autocorrelation")
    if (model_source != "none") {
      log_progress("  → Model is NOT fully capturing spatial structure")
      log_progress("  → Consider adding spatial random effects (GP) or additional covariates")
    } else {
      log_progress("  → Expected for raw errors - indicates spatial structure to model")
    }
  } else if (moran_chm$estimate[1] > 0.2) {
    log_progress("  INTERPRETATION: Moderate positive spatial autocorrelation")
  } else {
    log_progress("  INTERPRETATION: Weak spatial autocorrelation - model captures structure well")
  }
}

if (!is.null(vg_params_chm)) {
  log_progress(sprintf("\nCHM Variogram:"))
  log_progress(sprintf("  Effective range: %.1f km", vg_params_chm$range_m / 1000))
  log_progress(sprintf("  Nugget/Sill ratio: %.3f", vg_params_chm$nugget_sill_ratio))
  
  if (vg_params_chm$nugget_sill_ratio < 0.25) {
    log_progress("  → Strong spatial dependence: %.0f%% of variance is spatially structured",
                 (1 - vg_params_chm$nugget_sill_ratio) * 100)
  } else if (vg_params_chm$nugget_sill_ratio > 0.75) {
    log_progress("  → Weak spatial dependence: most variance is noise/micro-scale")
  }
}

log_progress("\nNOTE ON GEDI SAMPLING:")
log_progress("  Moran's I may be inflated due to GEDI's along-track sampling design.")
log_progress("  Points within ~60m are from same acquisition, sharing atmospheric/geometric conditions.")
log_progress("  Distance-banded Moran's I (in moran_by_distance_band.csv) provides scale-specific insights.")

log_progress("=" %>% strrep(70))

log_progress("\n✓ Enhanced spatial diagnostics complete (FIXED VERSION)")

# =====================================================================
# SAVE CHECKPOINT
# =====================================================================

if (!exists("BATCH_MODE") || !BATCH_MODE) {
  log_progress("Saving checkpoint...")
  
  checkpoint_data <- list(
    spatial_complete = TRUE,
    model_source = model_source,
    moran_chm = if(exists("moran_chm")) moran_chm else NULL,
    moran_dtm = if(exists("moran_dtm")) moran_dtm else NULL,
    moran_by_site = moran_by_site,
    moran_bands_chm = moran_bands_chm,
    moran_bands_dtm = moran_bands_dtm,
    vg_params = vg_params,
    target_crs = TARGET_CRS
  )
  
  if (!is.null(vg_fits_chm)) checkpoint_data$vg_fits_chm <- vg_fits_chm
  if (!is.null(vg_fits_dtm)) checkpoint_data$vg_fits_dtm <- vg_fits_dtm
  
  save_checkpoint("12_spatial_fixed", checkpoint_data)
}
