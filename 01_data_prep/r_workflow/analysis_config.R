# =====================================================================
# analysis_config.R (FIXED VERSION)
# Global configuration and environment setup
# 
# FIXES APPLIED:
#   1. Product-specific site exclusion (Site 10 DTM only)
#   2. Corrected sprintf format string (lines 134-149)
# =====================================================================

# This file is sourced by all section scripts
# Do not put heavy computation here - only configuration!

suppressPackageStartupMessages({
  library(dplyr); library(data.table); library(readr)
  library(ggplot2); library(cowplot); library(patchwork)
  library(GGally); library(corrplot)
  library(brms);   library(loo);      library(bayesplot)
  library(tidybayes); library(posterior); library(matrixStats)
  library(performance); library(broom); library(broom.mixed)
  library(bayestestR)
  library(terra); library(sf)
  library(spdep); library(gstat); library(mgcv)
  library(progress)
})

# =====================================================================
# Site Exclusion Configuration (FIXED!)
# =====================================================================
# Site 10: +143m DTM offset (vertical datum issue) - exclude from DTM only
# Sites 5, 6: Normal bias values (QC scores 63.5, 68.5) - NO exclusion needed
# Sites 12, 14: High bias but NOT datum issues - handle separately if needed

SITES_EXCLUDE_CHM <- ""        # No CHM exclusions needed
SITES_EXCLUDE_DTM <- "10"      # Only Site 10 has datum issue affecting DTM

# Legacy environment variables for backward compatibility
# These are now DEPRECATED - use product-specific variables above
Sys.setenv(SITES_INCLUDE = "")
Sys.setenv(SITES_EXCLUDE = "")  # Cleared - product-specific vars take precedence

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/gpfs/data1/vclgp/lmaden/chpt1")
OUT_ROOT     <- PROJECT_ROOT
dir.create(file.path(OUT_ROOT, "plots"),   recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_ROOT, "tables"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_ROOT, "rasters"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_ROOT, "tmp"),     recursive = TRUE, showWarnings = FALSE)

out_plots   <- file.path(OUT_ROOT, "plots")
out_tables  <- file.path(OUT_ROOT, "tables")
out_rasters <- file.path(OUT_ROOT, "rasters")

# Graphics setup
options(bitmapType = "cairo")
theme_set(theme_cowplot())

# Resource configuration
get_total_cores <- function() {
  s <- Sys.getenv("SLURM_CPUS_PER_TASK", "")
  if (nzchar(s)) as.integer(s) else parallel::detectCores()
}

CPU_FRAC      <- as.numeric(Sys.getenv("CPU_FRAC", "0.10"))
RAM_FRAC      <- as.numeric(Sys.getenv("RAM_FRAC", "0.10"))
CPU_TOTAL     <- max(1L, get_total_cores())
CPU_BUDGET    <- max(1L, floor(CPU_TOTAL * CPU_FRAC))
SKIP_HEAVY    <- as.logical(Sys.getenv("SKIP_HEAVY", "FALSE"))

# PNG device setup
get_png_device <- function() {
  if (requireNamespace("ragg", quietly = TRUE)) return(ragg::agg_png)
  if (isTRUE(capabilities("cairo"))) {
    return(function(filename, width, height, units, res, ...) {
      grDevices::png(filename, width = width, height = height,
                     units = units, res = res, type = "cairo", ...)
    })
  }
  if (requireNamespace("Cairo", quietly = TRUE)) return(Cairo::CairoPNG)
  return(function(filename, width, height, units, res, ...) {
    grDevices::png(filename, width = width, height = height,
                   units = units, res = res, ...)
  })
}
PNG_DEV <- get_png_device()

# EDA toggles
EDA_ENABLE           <- as.logical(Sys.getenv("EDA_ENABLE", "TRUE"))
EDA_SAMPLE_FRAC      <- as.numeric(Sys.getenv("EDA_SAMPLE_FRAC", "0.25"))
EDA_MAX_N            <- as.integer(Sys.getenv("EDA_MAX_N", "750000"))
EDA_MAX_PER_GROUP    <- as.integer(Sys.getenv("EDA_MAX_PER_GROUP", "60000"))
EDA_GROUP_HINTS      <- strsplit(Sys.getenv("EDA_GROUP_HINTS", "lc_l1_code,site"), ",")[[1]] |> trimws()

# Threading limits
Sys.setenv(
  OMP_NUM_THREADS       = as.character(CPU_BUDGET),
  OPENBLAS_NUM_THREADS  = as.character(max(1L, min(2L, CPU_BUDGET))),
  MKL_NUM_THREADS       = as.character(max(1L, min(2L, CPU_BUDGET)))
)

if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(max(1L, min(2L, CPU_BUDGET)))
  RhpcBLASctl::omp_set_num_threads(CPU_BUDGET)
}

options(mc.cores = CPU_BUDGET)

# terra / GDAL tuning
terraOptions(
  tempdir  = file.path(OUT_ROOT, "tmp"),
  memfrac  = min(0.40, RAM_FRAC),
  progress = 1
)
sf::sf_use_s2(FALSE)

# CmdStan path
if (requireNamespace("cmdstanr", quietly = TRUE)) {
  cmdstan_path_env <- Sys.getenv("CMDSTAN", "")
  if (nzchar(cmdstan_path_env)) cmdstanr::set_cmdstan_path(cmdstan_path_env)
  options(brms.backend = "cmdstanr")
}

# data.table threading
data.table::setDTthreads(min(8L, CPU_BUDGET))
op <- options(datatable.integer64 = "character")

MGCV_THREADS <- max(1L, min(CPU_BUDGET, 8L))

# =====================================================================
# Model Fitting Configuration
# =====================================================================

# Create models directory
MODELS_DIR <- file.path(OUT_ROOT, "models")
dir.create(MODELS_DIR, recursive = TRUE, showWarnings = FALSE)

# MCMC sampling parameters - increased for Student-t models
WARMUP <- as.integer(Sys.getenv("MCMC_WARMUP", "1500"))  # Increased from 1000
ITER   <- as.integer(Sys.getenv("MCMC_ITER", "3000"))    # Increased from 2000

# Control parameters - increased for complex models
ADAPT_DELTA <- as.numeric(Sys.getenv("ADAPT_DELTA", "0.95"))        # Increased from 0.90
MAX_TREEDEPTH <- as.integer(Sys.getenv("MAX_TREEDEPTH", "14"))     # Increased from 12

# Model family - changed to Student-t
MODEL_FAMILY <- Sys.getenv("MODEL_FAMILY", "student")  # Changed from "gaussian"

# Number of chains
N_CHAINS <- as.integer(Sys.getenv("N_CHAINS", "4"))

# Log configuration (FIXED sprintf!)
cat(sprintf(
  paste0(
    "\n[BAYESIAN MODEL CONFIG]\n",
    "  Family: %s\n",
    "  MCMC: %d warmup + %d sampling = %d total iterations\n",
    "  Chains: %d\n",
    "  Control: adapt_delta=%.2f, max_treedepth=%d\n",
    "  Models dir: %s\n\n"
  ),
  MODEL_FAMILY, 
  WARMUP, 
  ITER - WARMUP, 
  ITER,
  N_CHAINS,
  ADAPT_DELTA, 
  MAX_TREEDEPTH,
  MODELS_DIR
))

# =====================================================================
# Data Ingest Configuration
# =====================================================================

ENRICHED_DIR <- Sys.getenv("ENRICHED_DIR", file.path(PROJECT_ROOT, "data", "enriched_by_site"))

# Legacy site inclusion/exclusion (DEPRECATED - kept for compatibility)
SITES_INCLUDE <- Sys.getenv("SITES_INCLUDE", "")
SITES_EXCLUDE <- Sys.getenv("SITES_EXCLUDE", "")

# Filtering parameters
apply_chm_forest_filter <- TRUE
chm_forest_thresh_m     <- 2
apply_dtm_forest_filter <- TRUE
dtm_forest_proxy        <- "als_chm_p90"
dtm_cover_threshold     <- 0.20

# Coordinate attachment
attach_coords   <- TRUE
gedi_base_dir   <- file.path(PROJECT_ROOT, "gedi")
gedi_file_name  <- "GEDI_site%s_hq_ALL.gpkg"
coord_crs_out   <- 4326

# Column definitions
common_cols <- c(
  "site","ecoregion","shot_number",
  "lc2022_l1_code","lc2022_l1_name",
  "lc2022_mode_l1_code","lc2022_mode_l1_name",
  "slope_valid_frac","slope_mean","slope_sd",
  "wsci","rh_98","cover",
  "aspect_sin_mean","aspect_cos_mean"
)

meta_cols <- c(
  "meta_stereo_any",
  "meta_off_nadir_avg","meta_sun_elev_avg",
  "meta_target_azimuth_avg","meta_az_concentration",
  "meta_abs_geoacc_avg","meta_rel_geoacc_avg",
  "meta_tot_ct","meta_stereo_ct","meta_stereo_ratio",
  "meta_fwd_ct","meta_rev_ct","meta_fwd_ratio","meta_rev_ratio",
  "meta_leaf_on_ct","meta_leaf_off_ct","meta_leaf_on_ratio",
  "meta_veh_ct_GE01","meta_veh_ratio_GE01","meta_veh_ct_GEO1","meta_veh_ratio_GEO1",
  "meta_veh_ct_WV01","meta_veh_ratio_WV01",
  "meta_veh_ct_WV02","meta_veh_ratio_WV02",
  "meta_veh_ct_WV03","meta_veh_ratio_WV03"
)

chm_cols <- c(
  "p3d_chm_mean","als_chm_mean",
  "als_chm_valid_frac","p3d_chm_valid_frac",
  "als_chm_p90","error_mean"
)

dtm_cols <- c(
  "p3d_dtm_mean","dep_dtm_mean",
  "dep_dtm_valid_frac","p3d_dtm_valid_frac",
  "dtm_error_mean"
)

sel_cols <- unique(c(common_cols, meta_cols, chm_cols, dtm_cols))

# =====================================================================
# Banner
# =====================================================================

try({
  gdal_v <- tryCatch(sf::gdal_version(), error = function(e) NA_character_)
  proj_v <- tryCatch(sf::proj_info()$version, error = function(e) NA_character_)
  cat(sprintf(
    "[chpt1 MODULAR - ENHANCED] root=%s | cores=%d (budget=%d) | GDAL=%s | PROJ=%s\n",
    PROJECT_ROOT, CPU_TOTAL, CPU_BUDGET, gdal_v, proj_v
  ))
})

# Log site exclusion configuration
cat(sprintf(
  paste0(
    "\n[SITE EXCLUSION CONFIG]\n",
    "  CHM exclusions: %s\n",
    "  DTM exclusions: %s\n\n"
  ),
  ifelse(nzchar(SITES_EXCLUDE_CHM), SITES_EXCLUDE_CHM, "(none)"),
  ifelse(nzchar(SITES_EXCLUDE_DTM), SITES_EXCLUDE_DTM, "(none)")
))

set.seed(2025)
