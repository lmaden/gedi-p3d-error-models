#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(data.table)
})

# ---------------- resource policy (polite defaults) ----------------
DT_THREADS <- 1L   # you can set to 2L or 3L later
options(mc.cores = 1L)
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
           VECLIB_MAXIMUM_THREADS="1", NUMEXPR_NUM_THREADS="1")
data.table::setDTthreads(DT_THREADS)

# ---------------- paths ----------------
MASTER_IN   <- "/gpfs/data1/vclgp/lmaden/chpt1/data/sampled_data_FINAL_with_lc.csv"
LOOKUPS_DIR <- "/gpfs/data1/vclgp/lmaden/chpt1/data/footprints/lookups"
OUT_DIR     <- "/gpfs/data1/vclgp/lmaden/chpt1/data/enriched_by_site"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Clean start: remove any previous outputs so we don't append to stale files
CLEAN_START <- TRUE
if (CLEAN_START) {
  old <- list.files(OUT_DIR, pattern = "^site_[0-9]{2}_enriched\\.csv(\\.gz)?$", full.names = TRUE)
  if (length(old)) file.remove(old)
}

# ---------------- schema in your master ----------------
SITE_COL <- "site"           # confirmed from your summary
SHOT_COL <- "shot_number"    # confirmed from your summary

# ---------------- helper: output path per site ----------------
out_path_for_site <- function(site_id) {
  file.path(OUT_DIR, sprintf("site_%02d_enriched.csv", as.integer(site_id)))
}

# ---------------- union of lookup columns (stable header) ----------------
lk_paths <- list.files(LOOKUPS_DIR,
                       pattern="^meta_lookup_site([0-9]{2})\\.(csv\\.gz|parquet)$",
                       full.names=TRUE)
if (!length(lk_paths)) stop("No lookup files found in: ", LOOKUPS_DIR)

lookup_cols_union <- unique(unlist(lapply(lk_paths, function(p) {
  if (grepl("\\.parquet$", p) && requireNamespace("arrow", quietly = TRUE)) {
    names(arrow::read_parquet(p, as_data_frame = TRUE, skip = 0, n_rows = 0))
  } else {
    names(data.table::fread(p, nrows = 0))
  }
})))
NEW_COLS <- setdiff(lookup_cols_union, "shot_key")  # meta_ columns to add

# Final column order: master columns first, then union of meta columns
master_header <- names(readr::read_csv(MASTER_IN, n_max = 0, show_col_types = FALSE))
FINAL_COL_ORDER <- c(master_header, setdiff(NEW_COLS, master_header))

# ---------------- lookup loading with tiny LRU cache ----------------
path_for_site <- function(site_id) {
  csv <- file.path(LOOKUPS_DIR, sprintf("meta_lookup_site%02d.csv.gz", as.integer(site_id)))
  pq  <- file.path(LOOKUPS_DIR, sprintf("meta_lookup_site%02d.parquet",  as.integer(site_id)))
  if (file.exists(pq)  && requireNamespace("arrow", quietly = TRUE)) return(pq)
  if (file.exists(csv)) return(csv)
  NA_character_
}

.lk_cache <- new.env(parent = emptyenv())
.lk_order <- character(0)
cache_limit <- 2L  # keep at most 2 sites resident

get_lookup_dt <- function(site_id) {
  sid <- as.integer(site_id); key <- as.character(sid)
  if (exists(key, envir = .lk_cache, inherits = FALSE)) {
    .lk_order <<- c(setdiff(.lk_order, key), key)  # MRU
    return(get(key, envir = .lk_cache, inherits = FALSE))
  }
  p <- path_for_site(sid)
  if (is.na(p)) return(NULL)
  if (grepl("\\.parquet$", p) && requireNamespace("arrow", quietly = TRUE)) {
    dt <- as.data.table(arrow::read_parquet(p))
  } else {
    dt <- data.table::fread(p, colClasses = c(shot_key = "character"))
  }
  setkey(dt, "shot_key")
  assign(key, dt, envir = .lk_cache)
  .lk_order <<- c(.lk_order, key)
  if (length(.lk_order) > cache_limit) {
    rm(list = .lk_order[1], envir = .lk_cache); .lk_order <<- .lk_order[-1]
  }
  dt
}

# ---------------- per-site header tracking ----------------
.header_written <- new.env(parent = emptyenv())

write_site_chunk <- function(dt_out, site_id) {
  path <- out_path_for_site(site_id)
  append_mode <- isTRUE(get0(as.character(site_id), envir = .header_written, inherits = FALSE))
  # Ensure stable schema for this chunk
  missing <- setdiff(FINAL_COL_ORDER, names(dt_out))
  if (length(missing)) for (nm in missing) dt_out[, (nm) := NA_real_]
  setcolorder(dt_out, FINAL_COL_ORDER)
  # Write (CSV uncompressed so we can append safely). We can gzip afterward.
  data.table::fwrite(dt_out, path,
                     append = append_mode,
                     col.names = !append_mode,
                     sep = ",", quote = FALSE, na = "", bom = FALSE)
  if (!append_mode) assign(as.character(site_id), TRUE, envir = .header_written)
}

# ---------------- streaming join ----------------
chunk_size <- 200000L   # tune down to 100k for even smaller RAM

cb <- function(chunk, pos) {
  DT <- as.data.table(chunk)
  
  # Keep shot_number as-is; create a TEMP shot_key for joining
  if (!(SHOT_COL %in% names(DT))) stop("Expected column ", SHOT_COL, " in master.")
  DT[, shot_key := as.character(get(SHOT_COL))]
  
  # Split by site, join, and append to that site's file
  sites <- sort(unique(DT[[SITE_COL]]))
  for (s in sites) {
    xi <- DT[get(SITE_COL) == s]
    lk <- get_lookup_dt(s)
    if (!is.null(lk)) {
      setkey(xi, "shot_key")
      xo <- lk[xi]   # left join: preserve xi rows; add lookup cols
    } else {
      xo <- xi
      # Ensure missing meta columns exist so schema stays stable
      missing <- setdiff(NEW_COLS, names(xo))
      if (length(missing)) for (nm in missing) xo[, (nm) := NA_real_]
    }
    # Drop the TEMP join key; keep your original shot_number column
    if ("shot_key" %in% names(xo)) xo[, shot_key := NULL]
    
    write_site_chunk(xo, s)
  }
  
  rm(DT); gc()
  invisible(NULL)
}

readr::read_csv_chunked(
  MASTER_IN,
  callback   = readr::SideEffectChunkCallback$new(cb),
  chunk_size = chunk_size,
  # ensure shot_number is read as character to avoid 64-bit issues
  col_types  = cols(.default = col_guess(), !!SHOT_COL := col_character()),
  progress   = TRUE
)

cat("✔ Wrote per-site enriched CSVs to: ", OUT_DIR, "\n")
cat("   Example file: ", out_path_for_site(1), "\n")

# ---------------- optional: gzip the outputs at the end ----------------
GZIP_AFTER <- TRUE
if (GZIP_AFTER) {
  outs <- list.files(OUT_DIR, pattern = "^site_[0-9]{2}_enriched\\.csv$", full.names = TRUE)
  for (f in outs) system2("gzip", c("-f", f))
  cat("✔ Compressed per-site CSVs to .csv.gz\n")
}