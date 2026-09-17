# =====================================================================
# run_interactive.R
# Helper for running sections interactively
# =====================================================================

# Load core dependencies
source("analysis_config.R")
source("analysis_utils.R")

# Set interactive mode flag (enables checkpoint saves in section scripts)
BATCH_MODE <- FALSE


log_progress("Interactive mode loaded. Available functions:")
cat("\n")
cat("  list_checkpoints()              - Show saved checkpoints\n")
cat("  clear_checkpoint('name')        - Clear a checkpoint\n")
cat("  clear_all_checkpoints()         - Clear all checkpoints\n")
cat("\n")
cat("To run sections:\n")
cat("  source('section_01_ingest.R')   - Data ingest\n")
cat("  source('section_02_qc.R')       - Quality control\n")
cat("  source('section_03_eda.R')      - Core EDA\n")
cat("  source('section_04_terrain_lc.R') - Terrain & LC\n")
cat("  source('section_05_metadata.R') - Metadata analysis\n")
cat("  source('section_06_collinearity.R') - Collinearity\n")
cat("  source('section_07_summaries.R') - Summaries\n")
cat("  source('section_08_model_prep.R') - Model prep\n")
cat("  source('section_09_models_stage1.R') - Models S1\n")
cat("  source('section_10_models_stage2.R') - Models S2\n")
cat("  source('section_11_diagnostics.R') - Diagnostics\n")
cat("  source('section_12_spatial.R')  - Spatial analysis\n")
cat("  source('section_13_summary.R')  - Final summary\n")
cat("\n")
cat("To load data from checkpoints:\n")
cat("  data <- load_checkpoint('01_data_ingest')\n")
cat("  chm_df <- data$chm_df\n")
cat("  dtm_df <- data$dtm_df\n")
cat("\n")

# Helper function to load data if available
load_data <- function() {
  if (checkpoint_exists("01_data_ingest")) {
    log_progress("Loading data from checkpoint...")
    data <- load_checkpoint("01_data_ingest")
    assign("chm_df", data$chm_df, envir = .GlobalEnv)
    assign("dtm_df", data$dtm_df, envir = .GlobalEnv)
    log_progress(sprintf("✓ CHM: %s rows | DTM: %s rows",
                         format(nrow(chm_df), big.mark=","),
                         format(nrow(dtm_df), big.mark=",")))
  } else {
    log_progress("No data checkpoint found. Run: source('section_01_ingest.R')")
  }
}

load_models <- function() {
  if (checkpoint_exists("08_model_prep")) {
    log_progress("Loading model datasets from checkpoint...")
    mod_data <- load_checkpoint("08_model_prep")
    assign("mod_chm_s1", mod_data$mod_chm_s1, envir = .GlobalEnv)
    assign("mod_dtm_s1", mod_data$mod_dtm_s1, envir = .GlobalEnv)
    assign("mod_chm_s2", mod_data$mod_chm_s2, envir = .GlobalEnv)
    assign("mod_dtm_s2", mod_data$mod_dtm_s2, envir = .GlobalEnv)
    log_progress("✓ Model datasets loaded")
  } else {
    log_progress("No model prep checkpoint. Run: source('section_08_model_prep.R')")
  }
}

cat("Helper functions:\n")
cat("  load_data()    - Load chm_df, dtm_df from checkpoint\n")
cat("  load_models()  - Load modeling datasets from checkpoint\n")
cat("\n")