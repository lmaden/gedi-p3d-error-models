#!/usr/bin/env Rscript
# -------------------------------------------------------------------
# 02_sync_user_library.R
#
# Re-install (sync) all packages you previously had under
#   /gpfs/data1/vclgp/lmaden/Rlib
# into the *current* R user library (R 4.5.0), without ever
# exceeding 30 total CPU threads on a shared node.
# -------------------------------------------------------------------

suppressPackageStartupMessages({
  library(utils)
})

# ----------- USER PATHS (gsapp22) -----------
ROOT_RLIB   <- "/gpfs/data1/vclgp/lmaden/Rlib"
TMPDIR_PATH <- "/gpfs/data1/vclgp/lmaden/Rtmp"
LOG_DIR     <- "/gpfs/data1/vclgp/lmaden/chpt1/logs"

dir.create(TMPDIR_PATH, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR,     recursive = TRUE, showWarnings = FALSE)

# ----------- ENV + REPOS (keep compiles & caches on GPFS) -----------
Sys.setenv(TMPDIR = TMPDIR_PATH)
options(repos = c(
  PPM = "https://packagemanager.posit.co/cran/latest",
  CRAN = "https://cloud.r-project.org"
))

# ----------- CAP PARALLELISM ≤ 30 TOTAL THREADS -----------
if (!requireNamespace("parallelly", quietly = TRUE)) {
  install.packages("parallelly", quiet = TRUE)
}
CAP   <- 30L
AVAIL <- parallelly::availableCores(ensureLoadAverage = TRUE)  # respects cgroups/load
WORKERS <- max(1L, min(CAP, AVAIL))                            # never exceed CAP or AVAIL
# Split: PKG_PAR packages in parallel, each compiled with MAKE_PAR threads
PKG_PAR  <- min(6L, WORKERS)         # up to 6 packages concurrently
MAKE_PAR <- max(1L, WORKERS %/% PKG_PAR)
stopifnot(PKG_PAR * MAKE_PAR <= WORKERS)

options(Ncpus = PKG_PAR)
Sys.setenv(
  MAKEFLAGS = paste0("-j", MAKE_PAR),
  # also cap threaded math libs to avoid oversubscription
  OMP_NUM_THREADS        = as.character(MAKE_PAR),
  OPENBLAS_NUM_THREADS   = as.character(MAKE_PAR),
  MKL_NUM_THREADS        = as.character(MAKE_PAR),
  BLIS_NUM_THREADS       = as.character(MAKE_PAR),
  VECLIB_MAXIMUM_THREADS = as.character(MAKE_PAR),
  NUMEXPR_NUM_THREADS    = as.character(MAKE_PAR)
)

message(sprintf("CPU cap: total ≤ %d | installing %d package(s) at a time with make -j%d",
                WORKERS, PKG_PAR, MAKE_PAR))

# ----------- DETERMINE TARGET (CURRENT) USER LIB -----------
target_lib <- .libPaths()[1]
dir.create(target_lib, recursive = TRUE, showWarnings = FALSE)
# Ensure it's first on the path for this session
.libPaths(unique(c(target_lib, .libPaths())))

# ----------- DISCOVER SOURCE LIBRARIES TO SYNC FROM -----------
# Candidates: any "*-library/*" under ROOT_RLIB + literal macro path (if it exists)
macro_lit <- file.path(ROOT_RLIB, "%p-library", "%v")
cands <- unique(c(
  list.dirs(ROOT_RLIB, full.names = TRUE, recursive = TRUE),
  macro_lit
))
# Keep only directories that look like versioned libs and are not the target
is_valid_lib <- function(p) {
  if (!dir.exists(p)) return(FALSE)
  # heuristic: path contains "-library/" and has package folders inside
  grepl("-library", p, fixed = TRUE) && length(list.files(p, pattern = "^[A-Za-z]", no.. = TRUE)) > 0
}
source_libs <- unique(Filter(is_valid_lib, cands))
source_libs <- setdiff(normalizePath(source_libs, winslash = "/", mustWork = FALSE),
                       normalizePath(target_lib, winslash = "/", mustWork = FALSE))

if (length(source_libs) == 0L) {
  message("No source libraries found under: ", ROOT_RLIB,
          "\nNothing to sync. Exiting cleanly.")
  quit(status = 0L)
}

message("Target lib:   ", target_lib)
message("Source libs:\n  - ", paste(source_libs, collapse = "\n  - "))

# ----------- BUILD THE SET OF PACKAGES TO INSTALL -----------
# Helper: installed packages in a given library (names only)
pkgs_in <- function(lib) {
  ip <- tryCatch(utils::installed.packages(lib.loc = lib), error = function(e) NULL)
  if (is.null(ip) || nrow(ip) == 0L) character(0) else rownames(ip)
}

# Union of packages in all source libs
src_pkgs <- unique(unlist(lapply(source_libs, pkgs_in), use.names = FALSE))
# Drop base/recommended; keep user-installed only
base_recommended <- rownames(installed.packages(priority = c("base","recommended"),
                                                lib.loc = .Library.site))
src_pkgs <- setdiff(src_pkgs, base_recommended)

# What is already present in the target?
tgt_pkgs <- pkgs_in(target_lib)

need <- setdiff(src_pkgs, tgt_pkgs)
message(sprintf("Discovered %d unique package(s) in source libs; %d already in target; %d need install.",
                length(src_pkgs), length(src_pkgs) - length(need), length(need)))

# Write a manifest of what we plan to install (for records)
ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
manifest <- file.path(LOG_DIR, sprintf("pkg_sync_manifest_%s.txt", ts))
writeLines(sort(need), con = manifest)
message("→ Manifest of packages to install: ", manifest)

if (length(need) == 0L) {
  message("Nothing to install. Environment already synced.")
  quit(status = 0L)
}

# ----------- INSTALL, WITH ROBUST FALLBACK FOR FAILURES -----------
log_file <- file.path(LOG_DIR, sprintf("pkg_sync_install_%s.log", ts))
zz <- file(log_file, open = "wt"); sink(zz, type = "output"); sink(zz, type = "message")
on.exit({ try(sink(type = "message")); try(sink()); try(close(zz), silent = TRUE) }, add = TRUE)

message("Starting installation into: ", target_lib)
message(sprintf("Parallelism: Ncpus=%d, MAKEFLAGS=%s, TMPDIR=%s",
                getOption("Ncpus"), Sys.getenv("MAKEFLAGS"), Sys.getenv("TMPDIR")))

# Chunk to keep log readable; installation parallelism is controlled by Ncpus/MAKEFLAGS
chunk_size <- 40L
chunks <- split(need, ceiling(seq_along(need) / chunk_size))
success <- character(0)
failed  <- character(0)

for (i in seq_along(chunks)) {
  pk <- chunks[[i]]
  message(sprintf("\n[%d/%d] Installing %d package(s): %s",
                  i, length(chunks), length(pk), paste(pk, collapse = ", ")))
  try({
    install.packages(pk, lib = target_lib,
                     dependencies = c("Depends","Imports","LinkingTo"),
                     Ncpus = getOption("Ncpus"))
    success <- c(success, pk)
  }, silent = TRUE)
  
  # Verify each package actually landed; if not, try individually to capture errors
  landed <- intersect(pk, pkgs_in(target_lib))
  if (length(landed) < length(pk)) {
    retry <- setdiff(pk, landed)
    message(sprintf("  → %d package(s) need individual retry: %s",
                    length(retry), paste(retry, collapse = ", ")))
    for (p in retry) {
      ok <- FALSE
      try({
        install.packages(p, lib = target_lib,
                         dependencies = c("Depends","Imports","LinkingTo"),
                         Ncpus = 1L)    # retry single-threaded to reduce flakiness
        ok <- p %in% pkgs_in(target_lib)
      }, silent = TRUE)
      if (ok) success <- c(success, p) else failed <- c(failed, p)
    }
  }
  flush.console()
  Sys.sleep(0.1)  # polite pacing
}

# ----------- SUMMARY & SANITY CHECKS -----------
sink(type = "message"); sink()
close(zz)

# Re-open log tail to console
cat(readLines(log_file, n = 50L), sep = "\n")  # last 50 lines (context)

success <- sort(unique(success))
failed  <- sort(unique(setdiff(failed, success)))
post_tgt <- pkgs_in(target_lib)

cat("\n================ SUMMARY ================\n")
cat("Target library: ", target_lib, "\n", sep = "")
cat("Installed now:  ", length(post_tgt), " packages total\n", sep = "")
cat("Newly added:    ", length(success), "\n", sep = "")
if (length(failed)) {
  cat("FAILED:         ", length(failed), " → see log:\n", log_file, "\n", sep = "")
  cat("Failed names:   ", paste(failed, collapse = ", "), "\n", sep = "")
} else {
  cat("FAILED:         0\n")
}
cat("Manifest:       ", manifest, "\n", sep = "")
cat("Install log:    ", log_file, "\n", sep = "")
cat("R version:      ", as.character(getRversion()), "\n", sep = "")
cat("Platform:       ", R.version$platform, "\n", sep = "")
cat("Ncpus:          ", getOption("Ncpus"), " | MAKEFLAGS=", Sys.getenv("MAKEFLAGS"), "\n", sep = "")
cat("TMPDIR:         ", Sys.getenv("TMPDIR"), "\n", sep = "")
cat("=========================================\n\n")

# Final check: load a few core packages if present
for (p in c("dplyr","lubridate","sf","terra","cmdstanr")) {
  if (p %in% post_tgt) {
    cat("Testing library(", p, ") … ", sep = "")
    ok <- try(suppressPackageStartupMessages(require(p, character.only = TRUE)), silent = TRUE)
    cat(if (isTRUE(ok)) "OK\n" else "FAILED (load)\n")
  }
}

# If geospatial packages failed to install, print a helpful hint
if (length(failed) && any(failed %in% c("sf","terra","rgdal","stars"))) {
  cat("\nHint: geospatial packages need GDAL/PROJ/GEOS. Load modules, then re-run for failed ones:\n",
      "  module load gdal proj geos\n",
      "  install.packages(c('sf','terra'), lib='", target_lib, "')\n", sep = "")
}
