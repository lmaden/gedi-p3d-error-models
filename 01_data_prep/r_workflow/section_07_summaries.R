# =====================================================================
# section_07_summaries.R (ENHANCED VERSION)
# Stratified performance summaries
#
# ENHANCEMENTS ADDED:
#   1. Forest-focused stratified summaries (EBF, BDF, ENF, DNF only)
#   2. Deciduous vs Evergreen comparison tables
#   3. Broadleaf vs Needleleaf comparison tables
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

log_progress("Computing stratified performance summaries (ENHANCED)...")

# Define forest classes
FOREST_CLASSES <- c("EBF", "BDF", "ENF", "DNF")
log_progress(sprintf("  Forest classes for focused analysis: %s", 
                     paste(FOREST_CLASSES, collapse = ", ")))

# Create forest subsets
chm_forest <- chm_df %>% filter(lc_l1_code %in% FOREST_CLASSES)
dtm_forest <- dtm_df %>% filter(lc_l1_code %in% FOREST_CLASSES)

log_progress(sprintf("  Forest observations: CHM=%s (%.1f%%), DTM=%s (%.1f%%)",
                     format(nrow(chm_forest), big.mark=","),
                     100 * nrow(chm_forest) / nrow(chm_df),
                     format(nrow(dtm_forest), big.mark=","),
                     100 * nrow(dtm_forest) / nrow(dtm_df)))

# Add phenology/structure classification
chm_forest <- chm_forest %>%
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

dtm_forest <- dtm_forest %>%
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

# Helper functions for binning
mk_slope_bins <- function(slope, breaks = c(0,5,15,30,45,Inf)) {
  cut(slope, breaks = breaks, right = FALSE,
      labels = c("[0,5)","[5,15)","[15,30)","[30,45)","[45,+)"))
}

mk_qbin <- function(x, probs = c(0, .25, .5, .75, 1)) {
  q <- quantile(x, probs = probs, na.rm=TRUE)
  cut(x, unique(q), include.lowest = TRUE, ordered_result = TRUE)
}

summarize_error <- function(df, err_col, group_vars) {
  df %>%
    group_by(across(all_of(group_vars))) %>%
    summarize(
      n      = dplyr::n(),
      bias   = mean(.data[[err_col]], na.rm = TRUE),
      med    = median(.data[[err_col]], na.rm = TRUE),
      NMAD   = nm_ad(.data[[err_col]]),
      MAE    = mean(abs(.data[[err_col]]), na.rm = TRUE),
      RMSE   = sqrt(mean(.data[[err_col]]^2, na.rm = TRUE)),
      Q95abs = q95_abs(.data[[err_col]]),
      .groups = "drop"
    ) %>% arrange(desc(n))
}

# =====================================================================
# 1. ADD BINS TO DATAFRAMES (Original)
# =====================================================================

log_subsection("Creating stratification bins")

# Always add slope bin
chm_df <- chm_df %>% mutate(slope_bin = mk_slope_bins(slope_mean))
dtm_df <- dtm_df %>% mutate(slope_bin = mk_slope_bins(slope_mean))
chm_forest <- chm_forest %>% mutate(slope_bin = mk_slope_bins(slope_mean))
dtm_forest <- dtm_forest %>% mutate(slope_bin = mk_slope_bins(slope_mean))

# Conditionally add metadata bins
meta_bin_cols <- character(0)

if ("meta_offnad_z" %in% names(chm_df)) {
  chm_df <- chm_df %>% mutate(offnad_bin = mk_qbin(meta_offnad_z))
  dtm_df <- dtm_df %>% mutate(offnad_bin = mk_qbin(meta_offnad_z))
  chm_forest <- chm_forest %>% mutate(offnad_bin = mk_qbin(meta_offnad_z))
  dtm_forest <- dtm_forest %>% mutate(offnad_bin = mk_qbin(meta_offnad_z))
  meta_bin_cols <- c(meta_bin_cols, "offnad_bin")
}

if ("meta_sunel_z" %in% names(chm_df)) {
  chm_df <- chm_df %>% mutate(sunel_bin = mk_qbin(meta_sunel_z))
  dtm_df <- dtm_df %>% mutate(sunel_bin = mk_qbin(meta_sunel_z))
  chm_forest <- chm_forest %>% mutate(sunel_bin = mk_qbin(meta_sunel_z))
  dtm_forest <- dtm_forest %>% mutate(sunel_bin = mk_qbin(meta_sunel_z))
  meta_bin_cols <- c(meta_bin_cols, "sunel_bin")
}

if ("meta_leafon_z" %in% names(chm_df)) {
  chm_df <- chm_df %>% mutate(leafon_bin = mk_qbin(meta_leafon_z))
  dtm_df <- dtm_df %>% mutate(leafon_bin = mk_qbin(meta_leafon_z))
  chm_forest <- chm_forest %>% mutate(leafon_bin = mk_qbin(meta_leafon_z))
  dtm_forest <- dtm_forest %>% mutate(leafon_bin = mk_qbin(meta_leafon_z))
  meta_bin_cols <- c(meta_bin_cols, "leafon_bin")
}

if ("meta_stereo_z" %in% names(chm_df)) {
  chm_df <- chm_df %>% mutate(stereo_bin = mk_qbin(meta_stereo_z))
  dtm_df <- dtm_df %>% mutate(stereo_bin = mk_qbin(meta_stereo_z))
  chm_forest <- chm_forest %>% mutate(stereo_bin = mk_qbin(meta_stereo_z))
  dtm_forest <- dtm_forest %>% mutate(stereo_bin = mk_qbin(meta_stereo_z))
  meta_bin_cols <- c(meta_bin_cols, "stereo_bin")
}

log_progress(sprintf("  Created %d metadata bin columns: %s", 
                     length(meta_bin_cols), paste(meta_bin_cols, collapse=", ")))

# =====================================================================
# 2. TERRAIN/LC SUMMARIES - ALL DATA (Original)
# =====================================================================

log_subsection("Terrain and land cover stratifications (ALL DATA)")
tab_chm_slope     <- summarize_error(chm_df, "chm_error_mean", c("slope_bin"))
tab_chm_lc        <- summarize_error(chm_df, "chm_error_mean", c("lc_l1_code"))
tab_chm_ecoregion <- summarize_error(chm_df, "chm_error_mean", c("ecoregion"))
tab_chm_lc_slope  <- summarize_error(chm_df, "chm_error_mean", c("lc_l1_code","slope_bin"))

tab_dtm_slope     <- summarize_error(dtm_df, "dtm_error_mean", c("slope_bin"))
tab_dtm_lc        <- summarize_error(dtm_df, "dtm_error_mean", c("lc_l1_code"))
tab_dtm_ecoregion <- summarize_error(dtm_df, "dtm_error_mean", c("ecoregion"))
tab_dtm_lc_slope  <- summarize_error(dtm_df, "dtm_error_mean", c("lc_l1_code","slope_bin"))

# =====================================================================
# 3. NEW: TERRAIN/LC SUMMARIES - FOREST ONLY
# =====================================================================

log_subsection("Terrain and land cover stratifications (FOREST ONLY)")

# Forest types by slope
tab_chm_forest_slope <- summarize_error(chm_forest, "chm_error_mean", c("slope_bin"))
tab_dtm_forest_slope <- summarize_error(dtm_forest, "dtm_error_mean", c("slope_bin"))

# Individual forest types
tab_chm_forest_lc <- summarize_error(chm_forest, "chm_error_mean", c("lc_l1_code"))
tab_dtm_forest_lc <- summarize_error(dtm_forest, "dtm_error_mean", c("lc_l1_code"))

# Forest types crossed with slope
tab_chm_forest_lc_slope <- summarize_error(chm_forest, "chm_error_mean", c("lc_l1_code", "slope_bin"))
tab_dtm_forest_lc_slope <- summarize_error(dtm_forest, "dtm_error_mean", c("lc_l1_code", "slope_bin"))

# Forest types by ecoregion
tab_chm_forest_eco <- summarize_error(chm_forest, "chm_error_mean", c("ecoregion"))
tab_dtm_forest_eco <- summarize_error(dtm_forest, "dtm_error_mean", c("ecoregion"))

# Forest type × ecoregion
tab_chm_forest_lc_eco <- summarize_error(chm_forest, "chm_error_mean", c("lc_l1_code", "ecoregion"))
tab_dtm_forest_lc_eco <- summarize_error(dtm_forest, "dtm_error_mean", c("lc_l1_code", "ecoregion"))

log_progress("  Forest-only summary statistics:")
print(tab_chm_forest_lc %>% select(lc_l1_code, n, bias, RMSE))

# =====================================================================
# 4. NEW: PHENOLOGY COMPARISON (Deciduous vs Evergreen)
# =====================================================================

log_subsection("Phenology comparison: Deciduous vs Evergreen")

tab_chm_phenology <- summarize_error(chm_forest, "chm_error_mean", c("phenology"))
tab_dtm_phenology <- summarize_error(dtm_forest, "dtm_error_mean", c("phenology"))

tab_chm_phenology_slope <- summarize_error(chm_forest, "chm_error_mean", c("phenology", "slope_bin"))
tab_dtm_phenology_slope <- summarize_error(dtm_forest, "dtm_error_mean", c("phenology", "slope_bin"))

log_progress("  Deciduous vs Evergreen (CHM):")
print(tab_chm_phenology)

# =====================================================================
# 5. NEW: LEAF STRUCTURE COMPARISON (Broadleaf vs Needleleaf)
# =====================================================================

log_subsection("Leaf structure comparison: Broadleaf vs Needleleaf")

tab_chm_leafstruct <- summarize_error(chm_forest, "chm_error_mean", c("leaf_structure"))
tab_dtm_leafstruct <- summarize_error(dtm_forest, "dtm_error_mean", c("leaf_structure"))

tab_chm_leafstruct_slope <- summarize_error(chm_forest, "chm_error_mean", c("leaf_structure", "slope_bin"))
tab_dtm_leafstruct_slope <- summarize_error(dtm_forest, "dtm_error_mean", c("leaf_structure", "slope_bin"))

log_progress("  Broadleaf vs Needleleaf (CHM):")
print(tab_chm_leafstruct)

# =====================================================================
# 6. NEW: FULL FOREST TYPE MATRIX (Phenology × Leaf Structure)
# =====================================================================

log_subsection("Full forest type matrix")

tab_chm_forest_matrix <- summarize_error(chm_forest, "chm_error_mean", c("phenology", "leaf_structure"))
tab_dtm_forest_matrix <- summarize_error(dtm_forest, "dtm_error_mean", c("phenology", "leaf_structure"))

log_progress("  Phenology × Leaf Structure (CHM):")
print(tab_chm_forest_matrix)

# =====================================================================
# 7. METADATA SUMMARIES (Original)
# =====================================================================

log_subsection("Metadata stratifications")
meta_tabs_chm <- list()
meta_tabs_dtm <- list()
meta_tabs_chm_forest <- list()
meta_tabs_dtm_forest <- list()

if ("offnad_bin" %in% names(chm_df)) {
  meta_tabs_chm$offnad <- summarize_error(chm_df, "chm_error_mean", "offnad_bin")
  meta_tabs_dtm$offnad <- summarize_error(dtm_df, "dtm_error_mean", "offnad_bin")
  meta_tabs_chm_forest$offnad <- summarize_error(chm_forest, "chm_error_mean", "offnad_bin")
  meta_tabs_dtm_forest$offnad <- summarize_error(dtm_forest, "dtm_error_mean", "offnad_bin")
}

if ("sunel_bin" %in% names(chm_df)) {
  meta_tabs_chm$sunel <- summarize_error(chm_df, "chm_error_mean", "sunel_bin")
  meta_tabs_dtm$sunel <- summarize_error(dtm_df, "dtm_error_mean", "sunel_bin")
  meta_tabs_chm_forest$sunel <- summarize_error(chm_forest, "chm_error_mean", "sunel_bin")
  meta_tabs_dtm_forest$sunel <- summarize_error(dtm_forest, "dtm_error_mean", "sunel_bin")
}

if ("leafon_bin" %in% names(chm_df)) {
  meta_tabs_chm$leafon <- summarize_error(chm_df, "chm_error_mean", "leafon_bin")
  meta_tabs_dtm$leafon <- summarize_error(dtm_df, "dtm_error_mean", "leafon_bin")
  meta_tabs_chm_forest$leafon <- summarize_error(chm_forest, "chm_error_mean", "leafon_bin")
  meta_tabs_dtm_forest$leafon <- summarize_error(dtm_forest, "dtm_error_mean", "leafon_bin")
}

if ("stereo_bin" %in% names(chm_df)) {
  meta_tabs_chm$stereo <- summarize_error(chm_df, "chm_error_mean", "stereo_bin")
  meta_tabs_dtm$stereo <- summarize_error(dtm_df, "dtm_error_mean", "stereo_bin")
  meta_tabs_chm_forest$stereo <- summarize_error(chm_forest, "chm_error_mean", "stereo_bin")
  meta_tabs_dtm_forest$stereo <- summarize_error(dtm_forest, "dtm_error_mean", "stereo_bin")
}

# =====================================================================
# 8. WRITE ALL SUMMARY TABLES
# =====================================================================

log_subsection("Writing summary tables")

# Original tables (all data)
write_csv(tab_chm_slope,      file.path(out_tables, "chm_by_slope.csv"))
write_csv(tab_chm_lc,         file.path(out_tables, "chm_by_lc.csv"))
write_csv(tab_chm_ecoregion,  file.path(out_tables, "chm_by_ecoregion.csv"))
write_csv(tab_chm_lc_slope,   file.path(out_tables, "chm_by_lc_by_slope.csv"))
write_csv(tab_dtm_slope,      file.path(out_tables, "dtm_by_slope.csv"))
write_csv(tab_dtm_lc,         file.path(out_tables, "dtm_by_lc.csv"))
write_csv(tab_dtm_ecoregion,  file.path(out_tables, "dtm_by_ecoregion.csv"))
write_csv(tab_dtm_lc_slope,   file.path(out_tables, "dtm_by_lc_by_slope.csv"))

# Forest-only tables
write_csv(tab_chm_forest_slope,    file.path(out_tables, "chm_by_slope_FOREST.csv"))
write_csv(tab_chm_forest_lc,       file.path(out_tables, "chm_by_lc_FOREST.csv"))
write_csv(tab_chm_forest_lc_slope, file.path(out_tables, "chm_by_lc_by_slope_FOREST.csv"))
write_csv(tab_chm_forest_eco,      file.path(out_tables, "chm_by_ecoregion_FOREST.csv"))
write_csv(tab_chm_forest_lc_eco,   file.path(out_tables, "chm_by_lc_by_ecoregion_FOREST.csv"))

write_csv(tab_dtm_forest_slope,    file.path(out_tables, "dtm_by_slope_FOREST.csv"))
write_csv(tab_dtm_forest_lc,       file.path(out_tables, "dtm_by_lc_FOREST.csv"))
write_csv(tab_dtm_forest_lc_slope, file.path(out_tables, "dtm_by_lc_by_slope_FOREST.csv"))
write_csv(tab_dtm_forest_eco,      file.path(out_tables, "dtm_by_ecoregion_FOREST.csv"))
write_csv(tab_dtm_forest_lc_eco,   file.path(out_tables, "dtm_by_lc_by_ecoregion_FOREST.csv"))

# Phenology tables
write_csv(tab_chm_phenology,       file.path(out_tables, "chm_by_phenology_FOREST.csv"))
write_csv(tab_chm_phenology_slope, file.path(out_tables, "chm_by_phenology_by_slope_FOREST.csv"))
write_csv(tab_dtm_phenology,       file.path(out_tables, "dtm_by_phenology_FOREST.csv"))
write_csv(tab_dtm_phenology_slope, file.path(out_tables, "dtm_by_phenology_by_slope_FOREST.csv"))

# Leaf structure tables
write_csv(tab_chm_leafstruct,       file.path(out_tables, "chm_by_leafstructure_FOREST.csv"))
write_csv(tab_chm_leafstruct_slope, file.path(out_tables, "chm_by_leafstructure_by_slope_FOREST.csv"))
write_csv(tab_dtm_leafstruct,       file.path(out_tables, "dtm_by_leafstructure_FOREST.csv"))
write_csv(tab_dtm_leafstruct_slope, file.path(out_tables, "dtm_by_leafstructure_by_slope_FOREST.csv"))

# Forest type matrix tables
write_csv(tab_chm_forest_matrix, file.path(out_tables, "chm_forest_type_matrix.csv"))
write_csv(tab_dtm_forest_matrix, file.path(out_tables, "dtm_forest_type_matrix.csv"))

# Metadata tables (all data)
if (!is.null(meta_tabs_chm$offnad)) write_csv(meta_tabs_chm$offnad, file.path(out_tables, "chm_by_offnadir_quartiles.csv"))
if (!is.null(meta_tabs_chm$sunel))  write_csv(meta_tabs_chm$sunel,  file.path(out_tables, "chm_by_sunelev_quartiles.csv"))
if (!is.null(meta_tabs_chm$leafon)) write_csv(meta_tabs_chm$leafon, file.path(out_tables, "chm_by_leafon_quartiles.csv"))
if (!is.null(meta_tabs_chm$stereo)) write_csv(meta_tabs_chm$stereo, file.path(out_tables, "chm_by_stereo_quartiles.csv"))
if (!is.null(meta_tabs_dtm$offnad)) write_csv(meta_tabs_dtm$offnad, file.path(out_tables, "dtm_by_offnadir_quartiles.csv"))
if (!is.null(meta_tabs_dtm$sunel))  write_csv(meta_tabs_dtm$sunel,  file.path(out_tables, "dtm_by_sunelev_quartiles.csv"))
if (!is.null(meta_tabs_dtm$leafon)) write_csv(meta_tabs_dtm$leafon, file.path(out_tables, "dtm_by_leafon_quartiles.csv"))
if (!is.null(meta_tabs_dtm$stereo)) write_csv(meta_tabs_dtm$stereo, file.path(out_tables, "dtm_by_stereo_quartiles.csv"))

# Metadata tables (forest only)
if (!is.null(meta_tabs_chm_forest$offnad)) write_csv(meta_tabs_chm_forest$offnad, file.path(out_tables, "chm_by_offnadir_quartiles_FOREST.csv"))
if (!is.null(meta_tabs_chm_forest$sunel))  write_csv(meta_tabs_chm_forest$sunel,  file.path(out_tables, "chm_by_sunelev_quartiles_FOREST.csv"))
if (!is.null(meta_tabs_chm_forest$leafon)) write_csv(meta_tabs_chm_forest$leafon, file.path(out_tables, "chm_by_leafon_quartiles_FOREST.csv"))
if (!is.null(meta_tabs_chm_forest$stereo)) write_csv(meta_tabs_chm_forest$stereo, file.path(out_tables, "chm_by_stereo_quartiles_FOREST.csv"))
if (!is.null(meta_tabs_dtm_forest$offnad)) write_csv(meta_tabs_dtm_forest$offnad, file.path(out_tables, "dtm_by_offnadir_quartiles_FOREST.csv"))
if (!is.null(meta_tabs_dtm_forest$sunel))  write_csv(meta_tabs_dtm_forest$sunel,  file.path(out_tables, "dtm_by_sunelev_quartiles_FOREST.csv"))
if (!is.null(meta_tabs_dtm_forest$leafon)) write_csv(meta_tabs_dtm_forest$leafon, file.path(out_tables, "dtm_by_leafon_quartiles_FOREST.csv"))
if (!is.null(meta_tabs_dtm_forest$stereo)) write_csv(meta_tabs_dtm_forest$stereo, file.path(out_tables, "dtm_by_stereo_quartiles_FOREST.csv"))

log_progress("✓ Stratified summaries complete (ENHANCED)")

# Save checkpoint for interactive mode
if (!exists("BATCH_MODE") || !BATCH_MODE) {
  log_progress("Saving checkpoint for interactive mode...")
  checkpoint_data <- list(
    # Original tables
    tab_chm_slope = tab_chm_slope,
    tab_chm_lc = tab_chm_lc,
    tab_chm_ecoregion = tab_chm_ecoregion,
    tab_chm_lc_slope = tab_chm_lc_slope,
    tab_dtm_slope = tab_dtm_slope,
    tab_dtm_lc = tab_dtm_lc,
    tab_dtm_ecoregion = tab_dtm_ecoregion,
    tab_dtm_lc_slope = tab_dtm_lc_slope,
    meta_tabs_chm = meta_tabs_chm,
    meta_tabs_dtm = meta_tabs_dtm,
    # Forest tables
    tab_chm_forest_slope = tab_chm_forest_slope,
    tab_chm_forest_lc = tab_chm_forest_lc,
    tab_chm_forest_lc_slope = tab_chm_forest_lc_slope,
    tab_dtm_forest_slope = tab_dtm_forest_slope,
    tab_dtm_forest_lc = tab_dtm_forest_lc,
    tab_dtm_forest_lc_slope = tab_dtm_forest_lc_slope,
    tab_chm_phenology = tab_chm_phenology,
    tab_dtm_phenology = tab_dtm_phenology,
    tab_chm_leafstruct = tab_chm_leafstruct,
    tab_dtm_leafstruct = tab_dtm_leafstruct,
    forest_classes = FOREST_CLASSES
  )
  save_checkpoint("07_summaries", checkpoint_data)
}
