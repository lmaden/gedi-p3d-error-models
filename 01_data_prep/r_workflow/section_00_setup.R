# =====================================================================
# section_00_setup.R
# Environment setup and validation
# =====================================================================
# Ensure dependencies loaded
if (!exists("log_progress")) source("analysis_utils.R")
if (!exists("PROJECT_ROOT")) source("analysis_config.R")
log_progress("Validating environment configuration...")

# Verify directories exist
stopifnot(dir.exists(PROJECT_ROOT))
stopifnot(dir.exists(ENRICHED_DIR))

log_progress(sprintf("  Project root: %s", PROJECT_ROOT))
log_progress(sprintf("  CPU budget: %d/%d cores (%.0f%%)", 
                     CPU_BUDGET, CPU_TOTAL, 100*CPU_FRAC))
log_progress(sprintf("  RAM budget: %.0f%%", 100*RAM_FRAC))
log_progress(sprintf("  EDA enabled: %s (sample=%.0f%%, max=%s)", 
                     EDA_ENABLE, 100*EDA_SAMPLE_FRAC, 
                     format(EDA_MAX_N, big.mark=",")))

# Verify output directories
log_progress("Creating output directories...")
for (d in c(out_plots, out_tables, out_rasters)) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE)
    log_progress(sprintf("  Created: %s", d))
  }
}

log_progress("✓ Environment setup complete")

# Save checkpoint for interactive mode
if (!exists("BATCH_MODE") || !BATCH_MODE) {
  log_progress("Saving checkpoint for interactive mode...")
  save_checkpoint("00_environment", list(
    project_root = PROJECT_ROOT,
    cpu_budget = CPU_BUDGET,
    eda_enable = EDA_ENABLE
  ))
}