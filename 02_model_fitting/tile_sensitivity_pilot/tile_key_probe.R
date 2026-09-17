#!/usr/bin/env Rscript
# =====================================================================
# tile_key_probe.R  --  what is the natural tile key?
#
# Follows tile_premise_check.R (16 Aug): the ten modeled acquisition
# predictors vary within site, but the joint value-combination key over
# all meta columns has 19,450 levels (CHM) with median 3 footprints per
# level. Either a true tile/strip ID exists in the data, or the key
# must be chosen from the value partitions.
#
# Three read-only questions:
#   1. Does mod_*_s2 carry an ID-like column (tile, strip, granule,
#      scene, catalog id, ...)? Every column is printed with its class
#      and distinct count so nothing hides behind a name.
#   2. Do the six same-granularity "geometry" columns (median 166
#      distinct per site) share ONE partition per site? If the joint
#      distinct count equals the max single-column count at every
#      site, they are one grid.
#   3. Same for the coarser "stack" columns, and for the ten modeled
#      predictors together (the candidate refit key).
#
# USAGE (login node fine; single core, about a minute):
#   export PROJECT_ROOT=/gpfs/data1/vclgp/lmaden/chpt1
#   export R_LIBS=/gpfs/data1/vclgp/lmaden/Rlib
#   Rscript $PROJECT_ROOT/scripts/reviewed/tile_key_probe.R
# =====================================================================
suppressPackageStartupMessages(library(data.table))

unwrap <- function(p) if (is.list(p) && !is.null(p$data)) p$data else p
mp <- unwrap(readRDS(file.path(
  Sys.getenv("PROJECT_ROOT", "/gpfs/data1/vclgp/lmaden/chpt1"),
  "checkpoints", "08_model_prep.rds")))
cat("elements in 08_model_prep:", paste(names(mp), collapse = " "), "\n")

GEO   <- c("meta_offnad_z", "meta_sunel_z", "meta_az_conc_z", "meta_absgeo_z",
           "view_az_sin_z", "view_az_cos_z")
STACK <- c("meta_relgeo_z", "meta_stereo_z", "meta_fwdrev_z", "meta_leafon_z")
MODELED10 <- c(GEO, STACK)   # the ten acquisition predictors in the fitted models

partition_report <- function(dat, cols, label) {
  cols <- intersect(cols, names(dat))
  j  <- dat[, .(joint = uniqueN(.SD)), by = site, .SDcols = cols]
  s  <- dat[, lapply(.SD, uniqueN), by = site, .SDcols = cols]
  mx <- melt(s, id.vars = "site")[, .(max_single = max(value)), by = site]
  x  <- merge(j, mx, by = "site")
  setorder(x, -joint)
  cat("\n--", label, "(", length(cols),
      "columns ): joint vs max single-column distinct, per site\n")
  print(x)
  cat(label, ": TOTAL joint =", sum(x$joint),
      "| joint == max_single at", sum(x$joint == x$max_single),
      "of", nrow(x), "sites\n")
  nfp <- dat[, .N, by = c("site", cols)]
  cat(label, ": footprints per joint level: min", min(nfp$N),
      "| median", median(nfp$N), "| max", max(nfp$N), "\n")
}

for (nm in c("mod_chm_s2", "mod_dtm_s2")) {
  dat <- as.data.table(mp[[nm]])
  cat("\n====================", nm, ":", format(nrow(dat), big.mark = ","),
      "rows,", ncol(dat), "columns\n")

  if (nm == "mod_chm_s2") {
    cat("\n-- every column: class | distinct values over all rows\n")
    for (cn in names(dat))
      cat(sprintf("  %-30s %-10s %s\n", cn, class(dat[[cn]])[1],
                  format(uniqueN(dat[[cn]]), big.mark = ",")))
  }

  idp <- grep("id$|^id|tile|strip|granule|quad|scene|image|cat|shot|orbit|beam",
              names(dat), ignore.case = TRUE, value = TRUE)
  cat("\nname-matched ID candidates:",
      if (length(idp)) paste(idp, collapse = " ") else "NONE", "\n")
  for (cn in idp)
    cat(sprintf("  %-30s distinct overall %s | median distinct per site %s\n", cn,
                format(uniqueN(dat[[cn]]), big.mark = ","),
                format(median(dat[, uniqueN(get(cn)), by = site]$V1),
                       big.mark = ",")))

  partition_report(dat, GEO,       "GEOMETRY six")
  partition_report(dat, STACK,     "STACK four")
  partition_report(dat, MODELED10, "MODELED ten")
}
