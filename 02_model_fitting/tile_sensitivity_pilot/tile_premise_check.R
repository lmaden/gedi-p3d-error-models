#!/usr/bin/env Rscript
# =====================================================================
# tile_premise_check.R  --  do the meta_* acquisition columns vary within
# site, and what tile key do they define?
#
# Read-only. Gates M2 (tile-level grouping). v113 sec 2.2.4 predicts
# 31-166 distinct values per site, each shared by a median of 17-68
# footprints. Two failure modes to catch:
#   - site-constant (1 distinct value per site): tile key collapses
#     into the site random effect; M2 as planned cannot be built
#   - footprint-unique (distinct values ~ rows): no replication
#     structure; exact-value grouping is wrong, a true tile/strip ID
#     or rounded key is needed
#
# USAGE (login node is fine; single core, well under a minute):
#   export PROJECT_ROOT=/gpfs/data1/vclgp/lmaden/chpt1
#   export R_LIBS=/gpfs/data1/vclgp/lmaden/Rlib
#   Rscript $PROJECT_ROOT/scripts/reviewed/tile_premise_check.R
# =====================================================================
suppressPackageStartupMessages(library(data.table))

unwrap <- function(p) if (is.list(p) && !is.null(p$data)) p$data else p
mp <- unwrap(readRDS(file.path(
  Sys.getenv("PROJECT_ROOT", "/gpfs/data1/vclgp/lmaden/chpt1"),
  "checkpoints", "08_model_prep.rds")))

for (nm in c("mod_chm_s2", "mod_dtm_s2")) {

  dat <- as.data.table(mp[[nm]])          # same access pattern as loo_fold_chm.R
  meta_cols <- grep("^meta_|^view_az", names(dat), value = TRUE)

  cat("\n====", nm, ":", format(nrow(dat), big.mark = ","), "rows,",
      uniqueN(as.character(dat$site)), "sites\n")
  cat("detected acquisition columns:", paste(meta_cols, collapse = " "), "\n")
  if (!length(meta_cols)) { cat("FATAL: no meta_*/view_az columns detected\n"); next }

  dv   <- dat[, lapply(.SD, uniqueN), by = site, .SDcols = meta_cols]
  long <- melt(dv, id.vars = "site", variable.name = "column",
               value.name = "n_distinct")
  cat("\n-- distinct values per site, summarised over sites, per column:\n")
  print(long[, .(min = min(n_distinct), median = as.double(median(n_distinct)),
                 max = max(n_distinct)), by = column])

  flat <- long[, all(n_distinct == 1), by = column][V1 == TRUE, as.character(column)]
  cat("\ncolumns constant within EVERY site:",
      if (length(flat)) paste(flat, collapse = " ") else "NONE", "\n")

  key_n <- dat[, .(tiles = uniqueN(.SD)), by = site, .SDcols = meta_cols]
  setorder(key_n, -tiles)
  cat("\n-- tile-key levels (unique meta_* combinations) per site:\n")
  print(key_n)
  cat("TOTAL tile levels:", sum(key_n$tiles), "\n")

  tk <- dat[, .N, by = c("site", meta_cols)]
  cat("footprints per tile: min", min(tk$N), "| median", median(tk$N),
      "| max", max(tk$N), "\n")
}
