####### DIAGNOSTICS: read-only checks #######

# Helper to print tidy sections
hdr <- function(x) {
  cat("\n", paste(rep("=", nchar(x)), collapse = ""), "\n", x, "\n", 
      paste(rep("=", nchar(x)), collapse = ""), "\n", sep = "")
}
kv <- function(key, val) cat(sprintf("  %-28s : %s\n", key, if (length(val) == 0 || is.na(val)) "<NA>" else as.character(val)))

hdr("A) BASIC SYSTEM INFO")
kv("Sys.info()[\"nodename\"]", Sys.info()[["nodename"]])
kv("R.version.string", R.version.string)
kv("Platform (R.version$platform)", R.version$platform)
kv("getRversion()", as.character(getRversion()))
kv("R.home()", R.home())
kv("R.home('etc')", R.home("etc"))
kv("Which R (Sys.which('R'))", Sys.which("R"))
kv("Which Rscript", Sys.which("Rscript"))
kv("Interactive?", interactive())
kv("X11 DISPLAY", Sys.getenv("DISPLAY", unset = "<unset>"))
kv("Running in RStudio?", Sys.getenv("RSTUDIO", unset = "no"))

hdr("B) KEY ENVIRONMENT VARIABLES")
envs <- c("HOME","PWD","TMPDIR","R_LIBS_USER","R_LIBS","R_LIBS_SITE",
          "R_ENVIRON","R_ENVIRON_USER","R_PROFILE","R_PROFILE_USER",
          "CMDSTAN","RENV_PATHS_CACHE","XDG_CACHE_HOME",
          "OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS",
          "BLIS_NUM_THREADS","VECLIB_MAXIMUM_THREADS","NUMEXPR_NUM_THREADS",
          "IMG_META_WORKERS")
vals <- Sys.getenv(envs, unset = "<unset>")
for (i in seq_along(envs)) kv(envs[i], vals[i])

hdr("C) LIBRARY SEARCH PATHS & INSTALL TARGET")
print(.libPaths())
kv("Primary install lib (lib.loc used by install.packages)", .libPaths()[1])

hdr("D) WHERE WILL TEMP FILES GO?")
kv("tempdir()", normalizePath(tempdir(), mustWork = FALSE))
kv("tempfile(example)", normalizePath(dirname(tempfile()), mustWork = FALSE))

hdr("E) WORKING DIRECTORY")
kv("getwd()", getwd())

hdr("F) CRAN REPOS & COMPILE CORES")
print(getOption("repos"))
kv("options('Ncpus')", getOption("Ncpus"))
kv("parallel::detectCores()", tryCatch(parallel::detectCores(), error = function(e) NA_integer_))
if (requireNamespace("parallelly", quietly = TRUE)) {
  kv("parallelly::availableCores()", parallelly::availableCores())
} else {
  cat("  parallelly not installed; skipping availableCores() check.\n")
}

hdr("G) IMPORTANT DOTFILES IN PLAY")
# What R will try to read at startup (not editing anything; just showing presence + a peek)
rp_site <- file.path(R.home("etc"), "Rprofile.site")
rv_site <- file.path(R.home("etc"), "Renviron.site")
rp_user <- if (nzchar(Sys.getenv("R_PROFILE_USER"))) Sys.getenv("R_PROFILE_USER") else file.path(Sys.getenv("HOME"), ".Rprofile")
rv_user <- if (nzchar(Sys.getenv("R_ENVIRON_USER"))) Sys.getenv("R_ENVIRON_USER") else file.path(Sys.getenv("HOME"), ".Renviron")
files_to_show <- c("Rprofile.site" = rp_site, "Renviron.site" = rv_site, "User .Rprofile" = rp_user, "User .Renviron" = rv_user)
for (nm in names(files_to_show)) {
  f <- files_to_show[[nm]]
  kv(paste0(nm, " exists?"), file.exists(f))
  if (file.exists(f)) {
    cat("  --- first 30 lines of", nm, "(", f, ") ---\n")
    txt <- tryCatch(readLines(f, n = 30), error = function(e) paste("  [cannot read:", e$message, "]"))
    cat(paste0("  ", txt), sep = "\n"); cat("  --- end ---\n")
  }
}

hdr("H) VERIFY YOUR GPFS ROOTS EXIST")
gpfs_base <- "/gpfs/data1/vclgp/lmaden"
kv("GPFS base exists?", dir.exists(gpfs_base))
kv("GPFS write-test directory", file.path(gpfs_base, "_sanity"))

hdr("I) WRITE/READ TESTS (SAFE, TINY FILES)")
# We test writing to GPFS (allowed) and confirm we *don't* write to HOME unintentionally.
safe_dir <- file.path(gpfs_base, "_sanity")
ok_make <- tryCatch(dir.create(safe_dir, showWarnings = FALSE, recursive = TRUE), error = function(e) e)
kv("dir.create(GPFS/_sanity) ok?", if (isTRUE(ok_make)) "TRUE" else paste("FALSE:", conditionMessage(ok_make)))
test_file <- file.path(safe_dir, sprintf("write_test_%s.txt", Sys.getpid()))
ok_write <- tryCatch({ writeLines("hello-from-R", test_file); TRUE }, error = function(e) e)
kv("writeLines GPFS ok?", if (isTRUE(ok_write)) "TRUE" else paste("FALSE:", conditionMessage(ok_write)))
ok_read  <- tryCatch(readLines(test_file), error = function(e) e)
kv("readLines GPFS ok?", if (inherits(ok_read, "error")) paste("FALSE:", conditionMessage(ok_read)) else ok_read)
# Clean up tiny file (leave _sanity dir so you can inspect later)
try(unlink(test_file), silent = TRUE)

hdr("J) EXPAND NEW-STYLE USER LIB PATH (what %p-library/%v would become)")
# This shows *what* your R_LIBS_USER macro path would be, without requiring it to be set.
majmin <- paste0(R.version$major, ".", sub("\\..*$", "", R.version$minor))  # e.g., "4.4"
macro_expanded_userlib <- file.path(gpfs_base, "Rlib", sprintf("%s-library", R.version$platform), majmin)
kv("Macro-expanded userlib", macro_expanded_userlib)
kv("Legacy userlib (flat)",  file.path(gpfs_base, "Rlib"))

hdr("K) PACKAGES YOU CARE ABOUT: WHERE ARE THEY (IF INSTALLED)?")
pkgs_needed <- c(
  # original list
  "dplyr","data.table","readr","ggplot2","cowplot","patchwork",
  "GGally","corrplot","brms","loo","bayesplot","tidybayes","posterior",
  "matrixStats","performance","broom","broom.mixed","bayestestR",
  "terra","sf","spdep","gstat","mgcv",
  # second script adds:
  "stringr","lubridate","future","future.apply","parallelly","arrow","lwgeom"
)
ip <- installed.packages()[, c("Package","Version","LibPath")]
have <- ip[ip[,"Package"] %in% pkgs_needed, , drop = FALSE]
missing <- setdiff(pkgs_needed, have[, "Package"])
if (nrow(have)) {
  have <- have[order(have[, "Package"]), , drop = FALSE]
  print(as.data.frame(have), row.names = FALSE)
} else {
  cat("  None of the listed packages appear installed in the current library paths.\n")
}
kv("Missing packages count", length(missing))
if (length(missing)) print(sort(missing))

hdr("L) HEAVY GIS BACKENDS (IF sf/terra INSTALLED)")
if (requireNamespace("sf", quietly = TRUE)) {
  cat("sf extSoftVersion():\n"); print(sf::sf_extSoftVersion())
} else {
  cat("sf not installed or not in .libPaths(); skipping.\n")
}
if (requireNamespace("terra", quietly = TRUE)) {
  cat("terra::version():\n"); print(try(terra::version(), silent = TRUE))
} else {
  cat("terra not installed or not in .libPaths(); skipping.\n")
}

hdr("M) CMDSTAN (IF cmdstanr INSTALLED)")
if (requireNamespace("cmdstanr", quietly = TRUE)) {
  kv("cmdstanr::cmdstan_path()", tryCatch(cmdstanr::cmdstan_path(), error = function(e) paste("ERROR:", e$message)))
  kv("cmdstanr::cmdstan_version()", tryCatch(as.character(cmdstanr::cmdstan_version()), error = function(e) paste("ERROR:", e$message)))
} else {
  cat("cmdstanr not installed; skipping.\n")
}

hdr("N) THREAD / CORE CAP DETECTION (you said: NEVER > 30 cores)")
# What R thinks you can use (and what we will cap later)
detected <- tryCatch({
  if (requireNamespace("parallelly", quietly = TRUE)) parallelly::availableCores(ensureLoadAverage = TRUE)
  else parallel::detectCores()
}, error = function(e) NA_integer_)
kv("Detected cores (polite)", detected)
kv("Cap to", 30L)
cat("\nDiagnostics complete.\n")
