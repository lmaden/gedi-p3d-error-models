#!/usr/bin/env Rscript
# Local driver for the Fig. 1 site map.
#
# fig1_sitemap.R was written for the cluster. This driver supplies the
# few globals it expects and runs it on a workstation with R 4.4.0. Point FIG1_WORK_DIR
# at a directory holding that script and the two inputs listed in docs/MANIFEST.md
# (the site summary CSV and the CEC level-2 ecoregion shapefile); output lands there too.
WORK_DIR <- Sys.getenv("FIG1_WORK_DIR", ".")
stopifnot(dir.exists(WORK_DIR))

PROJECT_ROOT <- WORK_DIR
OUT_ROOT <- PROJECT_ROOT
log_progress    <- function(msg) cat(format(Sys.time(), "%H:%M:%S"), msg, "\n")
log_section     <- function(t) log_progress(t)
log_subsection  <- function(t) log_progress(paste("--", t))
suppressMessages(sf::sf_use_s2(FALSE))  # planar ops, matching the cluster build

setwd(WORK_DIR)
source("fig1_sitemap.R")
