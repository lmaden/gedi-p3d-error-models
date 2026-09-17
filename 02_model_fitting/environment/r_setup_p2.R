## 2A) Put BOTH libraries on the search path (new first, old second)
new_style <- "/gpfs/data1/vclgp/lmaden/Rlib/x86_64-pc-linux-gnu-library/4.5"
old_style <- "/gpfs/data1/vclgp/lmaden/Rlib"

dir.create(new_style, recursive = TRUE, showWarnings = FALSE)  # harmless if exists
.libPaths(unique(c(new_style, old_style, .libPaths())))
.libPaths()  # verify new order (new_style should be first)

## 2B) See what you now have (re-check just your needed packages)
pkgs_needed <- c(
  "dplyr","data.table","readr","ggplot2","cowplot","patchwork",
  "GGally","corrplot","brms","loo","bayesplot","tidybayes","posterior",
  "matrixStats","performance","broom","broom.mixed","bayestestR",
  "terra","sf","spdep","gstat","mgcv",
  "stringr","lubridate","future","future.apply","parallelly","arrow","lwgeom"
)
ip <- installed.packages()[, c("Package","Version","LibPath","Built")]
have <- ip[ip[,"Package"] %in% pkgs_needed, , drop = FALSE]
missing <- setdiff(pkgs_needed, have[, "Package"])
if (nrow(have)) have[order(have[, "Package"]), , drop = FALSE] else "None found"

## 2C) Try loading a few common packages now visible via old_style
suppressPackageStartupMessages({
  cat("\nTrying to load core tidyverse bits…\n")
  for (p in c("dplyr","ggplot2","readr","stringr")) {
    ok <- requireNamespace(p, quietly = TRUE)
    cat(sprintf("  %-12s : %s\n", p, if (ok) "OK (found)" else "MISSING"))
  }
})

## 2D) Install any missing *lightweight* packages into the new-style lib
##     (skip sf/terra/lwgeom for now—they need system GIS libs)
light <- intersect(missing, c(
  "dplyr","data.table","readr","ggplot2","cowplot","patchwork","GGally",
  "corrplot","brms","loo","bayesplot","tidybayes","posterior","matrixStats",
  "performance","broom","broom.mixed","bayestestR","gstat","spdep",
  "stringr","lubridate","future","future.apply","parallelly","posterior"
))
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (length(light)) {
  cat("\nInstalling lightweight CRAN packages to:", .libPaths()[1], "\n")
  install.packages(light, lib = .libPaths()[1], Ncpus = 1L)
}

## 2E) Leave the GIS-heavy ones for after we verify GDAL/GEOS/PROJ
heavy <- intersect(missing, c("sf","terra","lwgeom","arrow"))
cat("\nGIS/system-dependent packages still to handle later:", paste(heavy, collapse=", "), "\n")

## 2F) Use a robust core cap (≤ 30) even if availableCores() is weird
determine_workers <- function(min_target = 10L, max_target = 30L) {
  env_override <- Sys.getenv("IMG_META_WORKERS", unset = "")
  if (nzchar(env_override)) {
    w <- suppressWarnings(as.integer(env_override))
    if (!is.na(w) && w > 0L) return(max(1L, min(w, max_target)))
  }
  ac <- NA_integer_
  if (requireNamespace("parallelly", quietly = TRUE)) {
    ac <- suppressWarnings(parallelly::availableCores())
  }
  if (!is.na(ac)) return(max(1L, min(ac, max_target)))
  # Fallbacks
  dc <- suppressWarnings(try(parallel::detectCores(), silent = TRUE))
  dc <- if (inherits(dc, "try-error") || is.na(dc)) 1L else dc
  max(1L, min(as.integer(dc), max_target))
}
workers <- determine_workers(10L, 30L)
workers  # sanity: should be between 1 and 30

