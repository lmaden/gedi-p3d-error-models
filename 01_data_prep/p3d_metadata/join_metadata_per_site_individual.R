#!/usr/bin/env Rscript
# join_metadata_per_site_individual.R
# 
# Joins footprint metadata lookups to individual site CSVs (no merge required).
# Input:  per-site CSVs with LC data
# Output: per-site enriched CSVs with metadata columns attached

suppressPackageStartupMessages({
  library(data.table)
  library(future)
  library(future.apply)
  library(parallelly)
})

# ---------------- resource policy ----------------
DT_THREADS <- 2L
data.table::setDTthreads(DT_THREADS)

# ---------------- paths ----------------
CSV_WITH_LC_DIR <- "/gpfs/data1/vclgp/lmaden/chpt1/data/csvs_with_lc"
LOOKUPS_DIR     <- "/gpfs/data1/vclgp/lmaden/chpt1/data/footprints/lookups"
OUT_DIR         <- "/gpfs/data1/vclgp/lmaden/chpt1/data/enriched_by_site"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---------------- options ----------------
SHOT_COL    <- "shot_number"
GZIP_OUTPUT <- TRUE
MAX_WORKERS <- 10L

# Sites to process (will skip if input doesn't exist)
SITES <- 1:20

# ---------------- helpers ----------------
determine_workers <- function(max_target = MAX_WORKERS) {
  n <- NA_integer_
  if (requireNamespace("parallelly", quietly = TRUE)) {
    n <- suppressWarnings(parallelly::availableCores())
  }
  if (is.na(n)) {
    n <- suppressWarnings(try(parallel::detectCores(), silent = TRUE))
    n <- if (inherits(n, "try-error") || is.na(n)) 1L else n
  }
  max(1L, min(as.integer(n), max_target))
}

# Find the lookup file for a site (handles .csv.gz, .csv, .parquet)
find_lookup <- function(site_id) {
  patterns <- c(
    sprintf("meta_lookup_site%02d.csv.gz",  site_id),
    sprintf("meta_lookup_site%02d.parquet", site_id),
    sprintf("meta_lookup_site%02d.csv",     site_id),
    sprintf("meta_lookup_site%d.csv.gz",    site_id),
    sprintf("meta_lookup_site%d.parquet",   site_id),
    sprintf("meta_lookup_site%d.csv",       site_id)
  )
  for (p in patterns) {
    fp <- file.path(LOOKUPS_DIR, p)
    if (file.exists(fp)) return(fp)
  }
  NA_character_
}

# Find the input CSV for a site
find_input_csv <- function(site_id) {
  patterns <- c(
    sprintf("site_%02d_with_lc.csv", site_id),
    sprintf("site_%d_with_lc.csv",   site_id)
  )
  for (p in patterns) {
    fp <- file.path(CSV_WITH_LC_DIR, p)
    if (file.exists(fp)) return(fp)
  }
  NA_character_
}

# Load lookup table
load_lookup <- function(path) {
  if (grepl("\\.parquet$", path) && requireNamespace("arrow", quietly = TRUE)) {
    dt <- as.data.table(arrow::read_parquet(path))
  } else {
    dt <- data.table::fread(path, colClasses = c(shot_key = "character"))
  }
  setkey(dt, "shot_key")
  dt
}

# ---------------- main processing function ----------------
process_site <- function(site_id) {
  # Find input CSV
  csv_path <- find_input_csv(site_id)
  if (is.na(csv_path)) {
    message(sprintf("[site %02d] Input CSV not found, skipping.", site_id))
    return(invisible(NULL))
  }
  

  # Find lookup
  lk_path <- find_lookup(site_id)
  if (is.na(lk_path)) {
    message(sprintf("[site %02d] Lookup not found, skipping.", site_id))
    return(invisible(NULL))
  }
  
  message(sprintf("[site %02d] Processing: %s", site_id, basename(csv_path)))
  
  # Read input CSV
  dt <- data.table::fread(csv_path, colClasses = c(shot_number = "character"))
  
  # Create join key

dt[, shot_key := as.character(get(SHOT_COL))]
  
  # Load and join lookup
  lk <- load_lookup(lk_path)
  meta_cols <- setdiff(names(lk), "shot_key")
  
  setkey(dt, "shot_key")
  dt <- lk[dt]  # left join: all rows from dt, add lookup columns
  
  # Clean up join key
  dt[, shot_key := NULL]
  
  # Write output
  out_path <- file.path(OUT_DIR, sprintf("site_%02d_enriched.csv", site_id))
  data.table::fwrite(dt, out_path, sep = ",", quote = FALSE, na = "")
  
  # Optionally gzip
  if (GZIP_OUTPUT) {
    system2("gzip", c("-f", out_path), stdout = FALSE, stderr = FALSE)
    out_path <- paste0(out_path, ".gz")
  }
  
  message(sprintf("[site %02d] Wrote %s rows -> %s", site_id, format(nrow(dt), big.mark = ","), basename(out_path)))
  
  # Clean up
 rm(dt, lk)
  gc()
  
  invisible(out_path)
}

# ---------------- main ----------------
cat("============================================================\n")
cat("Joining footprint metadata to per-site CSVs\n")
cat("============================================================\n")
cat(sprintf("Input dir:   %s\n", CSV_WITH_LC_DIR))
cat(sprintf("Lookups dir: %s\n", LOOKUPS_DIR))
cat(sprintf("Output dir:  %s\n", OUT_DIR))
cat("============================================================\n\n")

# Check which sites have both inputs
valid_sites <- c()
for (s in SITES) {
  csv_ok <- !is.na(find_input_csv(s))
  lk_ok  <- !is.na(find_lookup(s))
  if (csv_ok && lk_ok) {
    valid_sites <- c(valid_sites, s)
  } else {
    if (!csv_ok) message(sprintf("[site %02d] Missing input CSV", s))
    if (!lk_ok)  message(sprintf("[site %02d] Missing lookup", s))
  }
}

cat(sprintf("\nFound %d sites with both inputs: %s\n\n", 
            length(valid_sites), paste(valid_sites, collapse = ", ")))

if (length(valid_sites) == 0) {
  stop("No valid sites found. Check paths.")
}

# Set up parallel processing
workers <- min(determine_workers(), length(valid_sites))
if (workers >= 2L) {
  message(sprintf("Using %d parallel workers.\n", workers))
  plan(multisession, workers = workers)
} else {
  message("Running sequentially.\n")
  plan(sequential)
}
on.exit(plan(sequential), add = TRUE)

# Process all sites
results <- future.apply::future_lapply(valid_sites, process_site,
                                        future.seed = TRUE)

cat("\n============================================================\n")
cat("COMPLETE\n")
cat(sprintf("Output files in: %s\n", OUT_DIR))
cat("============================================================\n")
