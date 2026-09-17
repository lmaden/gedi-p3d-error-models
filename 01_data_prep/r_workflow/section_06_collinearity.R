# =====================================================================
# section_06_collinearity.R
# Collinearity diagnostics (VIF checks)
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

log_progress("Running collinearity diagnostics...")

safe_check_collinearity <- function(df, cols, out_csv) {
  # Only use columns that actually exist in the dataframe
  cols_present <- intersect(cols, names(df))
  if (length(cols_present) < 2) {
    log_progress(sprintf("  ⚠ Skipping VIF: insufficient columns present (%d < 2)", length(cols_present)))
    return(invisible(NULL))
  }
  
  dd <- df %>% 
    dplyr::select(all_of(cols_present)) %>% 
    dplyr::mutate(across(everything(), as.numeric)) %>% 
    na.omit()
  
  if (nrow(dd) < 100) {
    log_progress(sprintf("  ⚠ Skipping VIF: insufficient data after na.omit (%d rows)", nrow(dd)))
    return(invisible(NULL))
  }
  
  if (nrow(dd) > 100000) dd <- dplyr::slice_sample(dd, n = 100000)
  fake_y <- rnorm(nrow(dd))
  fit <- lm(fake_y ~ ., data = dd)
  colinfo <- performance::check_collinearity(fit)
  readr::write_csv(as.data.frame(colinfo), out_csv)
  log_progress(sprintf("  VIF saved to: %s", basename(out_csv)))
}

# Build driver columns using only available metadata
base_drivers <- c("slope_mean_z","slope_sd_z","wsci_z","rh_98_z","cover_z")
meta_drivers_chm <- intersect(c("meta_offnad_z","meta_sunel_z","meta_az_conc_z"), 
                              available_meta_chm)
meta_drivers_dtm <- intersect(c("meta_offnad_z","meta_sunel_z","meta_az_conc_z"), 
                              available_meta_dtm)
view_drivers <- c("view_az_sin","view_az_cos")

driver_cols_chm <- c(base_drivers, meta_drivers_chm, view_drivers)
driver_cols_dtm <- c(base_drivers, meta_drivers_dtm, view_drivers)

log_progress(sprintf("  CHM drivers (%d): %s", length(driver_cols_chm), paste(driver_cols_chm, collapse=", ")))
log_progress(sprintf("  DTM drivers (%d): %s", length(driver_cols_dtm), paste(driver_cols_dtm, collapse=", ")))

safe_check_collinearity(chm_df, driver_cols_chm, 
                        file.path(out_tables, "collinearity_drivers_chm.csv"))
safe_check_collinearity(dtm_df, driver_cols_dtm, 
                        file.path(out_tables, "collinearity_drivers_dtm.csv"))

log_progress("✓ Collinearity diagnostics complete")

# Save checkpoint for interactive mode
if (!exists("BATCH_MODE") || !BATCH_MODE) {
  log_progress("Saving checkpoint for interactive mode...")
  save_checkpoint("06_collinearity", list(
    collin_complete = TRUE,
    driver_cols_chm = driver_cols_chm,
    driver_cols_dtm = driver_cols_dtm
  ))
}