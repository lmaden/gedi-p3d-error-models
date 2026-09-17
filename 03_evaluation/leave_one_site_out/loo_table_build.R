#!/usr/bin/env Rscript
# =====================================================================
# loo_table_build.R  --  read-only.  Supersedes loo_table_build.R.
#
# Fix to v1, v1's fault: the pooling helper took a parameter named
# `rung`, which is ALSO a column in the metrics tables. Inside
# DT[.rung == rung] data.table resolved `rung` to the column (character
# labels like "2_new_site_eco_lc_known"), not the argument, so every
# subset was empty and every cell came back NaN. Row selectors are now
# built as plain logical vectors outside data.table's scoping, and no
# helper argument shares a name with a column.
#
# Already verified by hand from v1's per-fold output, so these are
# expected to reproduce and the run is a formality for the other cells:
#   DTM all-18 rung 2  sqrt(2,286,599 / 202,278) = 3.3622  -> 3.362
#   CHM all-16 rung 2  sqrt(2,916,178 / 124,966) = 4.8307  -> 4.831
#   fold n sums 202,278 / 112,860 / 124,966 / 31,992, all exact.
#
# Pooling, established by loo_fold_gates.R:
#   RMSE  -> sqrt( sum(n_i * RMSE_i^2) / sum(n_i) )   SQUARED scale
#   others-> sum(n_i * x_i) / sum(n_i)                they are means
#
# WRITES NOTHING. Safe alongside the M2 DTM refit.
#
# USAGE:
#   export PROJECT_ROOT=/gpfs/data1/vclgp/lmaden/chpt1
#   export R_LIBS=/gpfs/data1/vclgp/lmaden/Rlib
#   nice -n 15 Rscript $PROJECT_ROOT/scripts/reviewed/loo_table_build.R \
#     2>&1 | tee $PROJECT_ROOT/manuscript_tables/loo_table_build.txt
# =====================================================================

suppressPackageStartupMessages(library(data.table))

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/gpfs/data1/vclgp/lmaden/chpt1")
LOO_DIR <- file.path(PROJECT_ROOT, "manuscript_tables", "loo")

hr  <- function() cat(strrep("=", 92), "\n")
sec <- function(x) { cat("\n"); hr(); cat("## ", x, "\n"); hr() }

cat("loo_table_build.R -- read-only\n")
cat("started: ", format(Sys.time()), "\n", sep = "")

recorded <- rbindlist(list(
  data.table(product="dtm", basis="all",   k=1, n=202278, rmse=3.003, bias=-0.339, crps=1.531, width90=8.910,  cov50=0.504, cov90=0.893, cov95=0.943),
  data.table(product="dtm", basis="all",   k=2, n=202278, rmse=3.362, bias=-0.182, crps=1.774, width90=10.533, cov50=0.529, cov90=0.892, cov95=0.942),
  data.table(product="dtm", basis="all",   k=3, n=202278, rmse=3.383, bias=-0.961, crps=1.763, width90=10.827, cov50=0.579, cov90=0.895, cov95=0.940),
  data.table(product="dtm", basis="valid", k=1, n=112860, rmse=3.106, bias=-0.229, crps=1.601, width90=9.510,  cov50=0.501, cov90=0.895, cov95=0.951),
  data.table(product="dtm", basis="valid", k=2, n=112860, rmse=3.534, bias= 0.392, crps=1.902, width90=10.331, cov50=0.448, cov90=0.878, cov95=0.940),
  data.table(product="dtm", basis="valid", k=3, n=112860, rmse=3.567, bias=-0.976, crps=1.878, width90=10.807, cov50=0.534, cov90=0.883, cov95=0.936),
  data.table(product="chm", basis="all",   k=1, n=124966, rmse=4.445, bias= 0.454, crps=2.220, width90=12.807, cov50=0.507, cov90=0.891, cov95=0.942),
  data.table(product="chm", basis="all",   k=2, n=124966, rmse=4.831, bias=-0.384, crps=2.608, width90=13.432, cov50=0.431, cov90=0.859, cov95=0.931),
  data.table(product="chm", basis="all",   k=3, n=124966, rmse=4.829, bias=-0.504, crps=2.600, width90=13.783, cov50=0.447, cov90=0.869, cov95=0.935),
  data.table(product="chm", basis="valid", k=1, n= 31992, rmse=4.079, bias=-0.163, crps=2.065, width90=12.858, cov50=0.486, cov90=0.910, cov95=0.963),
  data.table(product="chm", basis="valid", k=2, n= 31992, rmse=4.557, bias= 0.220, crps=2.476, width90=13.635, cov50=0.453, cov90=0.863, cov95=0.948),
  data.table(product="chm", basis="valid", k=3, n= 31992, rmse=4.527, bias=-0.172, crps=2.444, width90=14.135, cov50=0.488, cov90=0.879, cov95=0.951)
))

load_product <- function(product, suffix) {
  pat_any  <- sprintf("^loo_metrics_%s_site[0-9]+(_%s)?\\.csv$", product, suffix)
  pat_orig <- sprintf("^loo_metrics_%s_site[0-9]+\\.csv$",       product)
  pat_rr   <- sprintf("^loo_metrics_%s_site[0-9]+_%s\\.csv$",    product, suffix)
  files <- list.files(LOO_DIR, pattern = pat_any)
  orig  <- grep(pat_orig, files, value = TRUE)
  rr    <- grep(pat_rr,   files, value = TRUE)
  sid <- function(f) as.integer(sub(sprintf("^loo_metrics_%s_site([0-9]+).*$", product), "\\1", f))
  keep <- c(rr, orig[!sid(orig) %in% sid(rr)])
  cat(sprintf("\n%s: %d folds | %d replaced by _%s: %s\n", toupper(product),
              length(keep), length(rr), suffix, paste(sort(sid(rr)), collapse = ", ")))
  rbindlist(lapply(keep, function(f) {
    d <- fread(file.path(LOO_DIR, f)); d[, `:=`(.site = sid(f), .file = f)]; d
  }), fill = TRUE, use.names = TRUE)
}

pick <- function(d, pats, what) {
  for (p in pats) { h <- grep(p, names(d), ignore.case = TRUE, value = TRUE); if (length(h)) return(h[1]) }
  cat("  !! no column for ", what, "\n", sep = ""); NA_character_
}

pool_product <- function(DT, product) {
  C <- list(rung = pick(DT, c("^rung$"), "rung"), n = pick(DT, c("^n$"), "n"),
            rmse = pick(DT, c("^rmse$"), "rmse"), bias = pick(DT, c("^bias$"), "bias"),
            crps = pick(DT, c("^crps$"), "crps"), width90 = pick(DT, c("^width90$"), "width90"),
            cov50 = pick(DT, c("^cov50$"), "cov50"), cov90 = pick(DT, c("^cov90$"), "cov90"),
            cov95 = pick(DT, c("^cov95$"), "cov95"), eco = pick(DT, c("eco_known"), "eco_known"))

  # Plain vectors -- no data.table scoping anywhere below this line.
  rung_i <- as.integer(sub("^\\s*([0-9]).*$", "\\1", as.character(DT[[C$rung]])))
  site_i <- DT$.site
  nv     <- DT[[C$n]]
  ecov   <- as.logical(DT[[C$eco]])

  valid_sites <- sort(unique(site_i[rung_i == 2L & ecov %in% TRUE]))
  cat("  valid (eco_known) sites: ", paste(valid_sites, collapse = ", "),
      "  [", length(valid_sites), " of ", length(unique(site_i)), "]\n", sep = "")

  cat("\n  per-fold rung-2 rows:\n")
  s2 <- rung_i == 2L
  print(data.table(site = site_i[s2], n = nv[s2],
                   rmse = round(DT[[C$rmse]][s2], 4),
                   eco_known = ecov[s2], file = DT$.file[s2])[order(site)])

  one <- function(bas, kk) {
    sel <- rung_i == kk
    if (bas == "valid") sel <- sel & site_i %in% valid_sites
    n <- nv[sel]; N <- sum(n)
    wm <- function(cn) if (is.na(cn)) NA_real_ else sum(n * DT[[cn]][sel]) / N
    data.table(product = product, basis = bas, k = kk, folds = sum(sel), n = N,
               rmse    = sqrt(sum(n * DT[[C$rmse]][sel]^2) / N),
               bias = wm(C$bias), crps = wm(C$crps), width90 = wm(C$width90),
               cov50 = wm(C$cov50), cov90 = wm(C$cov90), cov95 = wm(C$cov95))
  }
  rbindlist(lapply(c("all", "valid"), function(b) rbindlist(lapply(1:3, function(k) one(b, k)))))
}

sec("RECOMPUTED POOLED METRICS")
got <- rbind(pool_product(load_product("dtm", "ad099"), "dtm"),
             pool_product(load_product("chm", "init0"), "chm"))

cat("\n")
print(got[, .(product, basis, rung = k, folds, n,
              rmse = round(rmse, 4), bias = round(bias, 4), crps = round(crps, 4),
              width90 = round(width90, 4), cov50 = round(cov50, 4),
              cov90 = round(cov90, 4), cov95 = round(cov95, 4))])

sec("DIFF AGAINST CHECKPOINT sec 6.2 / sec 6.3")
cmp <- merge(got, recorded, by = c("product", "basis", "k"), suffixes = c(".got", ".rec"))
setorder(cmp, product, basis, k)
metrics <- c("n", "rmse", "bias", "crps", "width90", "cov50", "cov90", "cov95")
tol <- c(n = 0, rmse = 5e-4, bias = 5e-4, crps = 5e-4, width90 = 5e-4,
         cov50 = 5e-4, cov90 = 5e-4, cov95 = 5e-4)

nbad <- 0L
for (i in seq_len(nrow(cmp))) {
  r <- cmp[i]
  cat(sprintf("\n%s  %-5s  rung %d   (folds %2d)\n", toupper(r$product), r$basis, r$k, r$folds))
  for (m in metrics) {
    g <- r[[paste0(m, ".got")]]; e <- r[[paste0(m, ".rec")]]
    ok <- isTRUE(abs(g - e) <= tol[[m]] + 1e-9)
    if (!ok) nbad <- nbad + 1L
    cat(sprintf("    %-8s recomputed %11.4f   recorded %11.4f   %s\n",
                m, g, e, if (ok) "OK" else "<-- DIFF"))
  }
}

sec(sprintf("SUMMARY: %d of 96 cell(s) differ beyond tolerance", nbad))
if (nbad == 0L) cat("All cells reproduce. Build the main-text LOO table from this run.\n")
cat("\nfinished: ", format(Sys.time()), "\n", sep = "")
