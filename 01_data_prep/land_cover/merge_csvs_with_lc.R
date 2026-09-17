#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

# Configuration
CSV_DIR <- "/gpfs/data1/vclgp/lmaden/chpt1/data/csvs_with_lc"
OUTPUT  <- "/gpfs/data1/vclgp/lmaden/chpt1/data/sampled_data_FINAL_with_lc.csv"

# Sites to merge (exclude 1 and 5)
SITES <- c(2:4, 6:20)

cat("Merging site CSVs with landcover...\n")

# Read and combine
dfs <- list()
for (site_id in SITES) {
  csv_path <- file.path(CSV_DIR, sprintf("site_%02d_with_lc.csv", site_id))
  
  if (!file.exists(csv_path)) {
    warning(sprintf("Site %d not found: %s", site_id, csv_path))
    next
  }
  
  cat(sprintf("Reading site %d...\n", site_id))
  df <- read_csv(csv_path, show_col_types = FALSE)
  dfs[[length(dfs) + 1]] <- df
}

# Combine
cat("\nCombining...\n")
merged <- bind_rows(dfs)

# Save
cat(sprintf("Writing merged CSV (%s rows)...\n", format(nrow(merged), big.mark=",")))
write_csv(merged, OUTPUT)

cat(sprintf("\n✔ Done! Wrote %s rows → %s\n", format(nrow(merged), big.mark=","), OUTPUT))
cat(sprintf("  Sites included: %s\n", paste(SITES, collapse=", ")))