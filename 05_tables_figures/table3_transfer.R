# table3_transfer.R  --  Table 3: what the models explain and whether a correction transfers
# Per product (DSM, DTM, CHM):
#   R2 fixed (share of error variance the variables explain), R2 total minus R2 fixed (share the
#   group offsets add), from the full fits (r2.csv);
#   at held-out sites, pooled over all held-out footprints, as the median site and with each site
#   weighted equally ("Sites weighted equally": root of the mean site MSE, mean site coverage and
#   width; the frame Table 3b prints): RMSE with no correction (observed error itself), after a
#   single average correction (training-site mean), after the model's correction (new site), after
#   subtracting the site's own mean error; and the share of true errors inside the model's 90%
#   prediction interval. pct_cut_own and within_share_mse give the drop in RMSE from removing the
#   site's own mean and the within-site share of the mean squared error.
# Inputs  (P3D_DATA_ROOT) rework/folds/rowpreds_{dsm,dtm,chm}_fold{1..6}.csv.gz (rework_fold.R);
#                         dsm_direct/dsm_direct_frame.csv.gz (through dsm_direct_rowpreds.R)
#         (P3D_DERIVED)   folds/metrics_{dsm,dtm,chm}_fold{1..6}.csv (rework_fold.R);
#                         full_fits/{dsm,dtm,chm}/r2.csv (rework_fullfit.R)
# Outputs (P3D_RESULTS)   table3/table3_transfer.{csv,md}; says PARTIAL until all 6 folds x 3 products.
# Checks  fold-scope metrics recomputed from rowpreds must equal metrics_*.csv (1e-9); every site
#         held out once per product; no site in two folds.
suppressPackageStartupMessages(library(data.table))
DATA_ROOT <- Sys.getenv("P3D_DATA_ROOT", "<LOCAL_DATA_ROOT>")
DERIVED   <- Sys.getenv("P3D_DERIVED", "data_derived")
OUT_ROOT  <- Sys.getenv("P3D_RESULTS", "results")
OUT  <- file.path(OUT_ROOT, "table3"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
ARGF <- grep("^--file=", commandArgs(FALSE), value = TRUE)   # folder of this script (Rscript), for source()
HERE <- if (length(ARGF)) dirname(normalizePath(sub("^--file=", "", ARGF[1]))) else "05_tables_figures"
PRODUCTS <- c("DSM", "DTM", "CHM")

# ---- fold outputs: plain tags only (a rerun with non-default settings gets a suffix; flag it)
rp <- list.files(file.path(DATA_ROOT, "rework", "folds"), "^rowpreds_(dsm|dtm|chm)_fold[0-9]+.*\\.csv\\.gz$", full.names = TRUE)
tag <- sub("^rowpreds_(.*)\\.csv\\.gz$", "\\1", basename(rp))
plain <- grepl("^(dsm|dtm|chm)_fold[0-9]+$", tag)
if (any(!plain)) cat("NOTE non-default fold runs present, check which to use:", tag[!plain], "\n")
rp <- rp[plain]; tag <- tag[plain]
stopifnot(length(rp) > 0)
R <- rbindlist(lapply(seq_along(rp), function(i) {
  x <- fread(rp[i], colClasses = list(character = c("site", "shot_number")))
  x[, `:=`(product = toupper(sub("_fold.*", "", tag[i])), fold = as.integer(sub(".*_fold", "", tag[i])))]
}))
have <- unique(R[, .(product, fold)])
cat("fold outputs found:\n"); print(dcast(have, fold ~ product, fun.aggregate = length, value.var = "fold"))
COMPLETE <- nrow(have) == 18L

# ---- checks
rmse <- function(p, o) sqrt(mean((p - o)^2))
M <- rbindlist(lapply(file.path(DERIVED, "folds", sprintf("metrics_%s.csv", tag)), fread), fill = TRUE)
re <- R[, .(new_site = rmse(pred_new_site, observed), constant_global = rmse(pred_const_global, observed),
            constant_site = rmse(pred_const_site, observed),
            cov90 = mean(observed >= lo90 & observed <= hi90)), by = .(product, fold)]
mm <- dcast(M[scope == "fold"], product + fold ~ rung, value.var = "rmse")
cmp <- merge(re, mm, by = c("product", "fold"))
d <- cmp[, max(abs(c(new_site.x - new_site.y, constant_global.x - constant_global.y, constant_site.x - constant_site.y)))]
cat(sprintf("GATE rowpreds reproduce metrics files, max |diff| = %.1e  %s\n", d, ifelse(d < 1e-9, "PASS", "FAIL")))
dup <- R[, uniqueN(fold), by = .(product, site)][V1 > 1]
cat(sprintf("GATE no site held out in two folds: %s\n", ifelse(nrow(dup) == 0, "PASS", "FAIL")))
sites_per <- R[, uniqueN(site), by = product]
cat("held-out sites so far per product:\n"); print(sites_per)
diag <- unique(M[scope == "fold", .(product, fold, divergences, max_rhat, min_ess_ratio, pct_at_ceiling, fold_diagnostics_ok)])
cat("fold diagnostics (a few divergences are reported, not rerun):\n"); print(diag[order(product, fold)])

# ---- DSM rows switched to the directly sampled DSM error (the checks above use the stored Eq. 4 values)
source(file.path(HERE, "dsm_direct_rowpreds.R"))
R <- use_direct_dsm(R, "fold", DATA_ROOT)

# ---- transfer rows
st <- function(x) list(n = nrow(x),
  rmse_none   = sqrt(mean(x$observed^2)),
  rmse_global = rmse(x$pred_const_global, x$observed),
  rmse_model  = rmse(x$pred_new_site, x$observed),
  rmse_own    = rmse(x$pred_const_site, x$observed),
  cov90       = mean(x$observed >= x$lo90 & x$observed <= x$hi90),
  width90     = mean(x$hi90 - x$lo90))
P <- R[, st(.SD), by = product][, frame := "All held-out footprints"]
S <- R[, st(.SD), by = .(product, site)]
Md <- S[, lapply(.SD, median), by = product, .SDcols = c("rmse_none", "rmse_global", "rmse_model", "rmse_own", "cov90", "width90")]
Md[, `:=`(n = S[, .N, by = product]$N, frame = "Median site")]
# each held-out site weighted equally, as in Table 2:
# RMSE = root of the mean of the 18 site MSEs; coverage and width = means of the site values.
Sw <- S[, .(n = sum(n), rmse_none = sqrt(mean(rmse_none^2)), rmse_global = sqrt(mean(rmse_global^2)),
            rmse_model = sqrt(mean(rmse_model^2)), rmse_own = sqrt(mean(rmse_own^2)),
            cov90 = mean(cov90), width90 = mean(width90)), by = product][, frame := "Sites weighted equally"]
TR <- rbind(P, Md, Sw, use.names = TRUE)
# drop in RMSE from removing the site's own mean, and the within-site share of the mean squared error
TR[, `:=`(pct_cut_own = 100 * (1 - rmse_own / rmse_none), within_share_mse = 100 * rmse_own^2 / rmse_none^2)]
TR[frame == "Median site", `:=`(pct_cut_own = NA_real_, within_share_mse = NA_real_)]

# ---- variance shares from the full fits
r2 <- rbindlist(lapply(PRODUCTS, function(p) {
  f <- file.path(DERIVED, "full_fits", tolower(p), "r2.csv")
  if (!file.exists(f)) return(data.table(product = p, r2_fixed = NA_real_, r2_total = NA_real_))
  x <- fread(f); data.table(product = p, r2_fixed = x[quantity == "R2_fixed_only"]$median,
                            r2_total = x[quantity == "R2_total"]$median)
}))
r2[, r2_groups := r2_total - r2_fixed]
T3 <- merge(TR, r2, by = "product", all.x = TRUE)
T3[, product := factor(product, levels = PRODUCTS)]; setorder(T3, frame, product)
T3[, status := ifelse(COMPLETE, "FINAL", "PARTIAL")]
fwrite(T3, file.path(OUT, "table3_transfer.csv"))

f2 <- function(x) ifelse(is.na(x), "", sub("-", "\u2212", sprintf("%.2f", x), fixed = TRUE))
md <- c(if (!COMPLETE) "**PARTIAL (not all folds done), do not cite**\n",
        "| | Product | Explained by variables (%) | Added by group offsets (%) | RMSE, no correction (m) | RMSE, single average correction (m) | RMSE, model correction (m) | RMSE, site's own mean removed (m) | Inside 90% interval (%) | Width of 90% interval (m) |",
        "|---|---|---|---|---|---|---|---|---|---|")
for (fm in c("All held-out footprints", "Median site", "Sites weighted equally")) for (p in PRODUCTS) {
  r <- T3[frame == fm & product == p]; if (!nrow(r)) next
  md <- c(md, sprintf("| %s | %s | %s | %s | %s | %s | %s | %s | %.0f | %.1f |", ifelse(p == "DSM", fm, ""), p,
                      ifelse(fm != "All held-out footprints", "", sprintf("%.1f", 100 * r$r2_fixed)),
                      ifelse(fm != "All held-out footprints", "", sprintf("%.1f", 100 * r$r2_groups)),
                      f2(r$rmse_none), f2(r$rmse_global), f2(r$rmse_model), f2(r$rmse_own), 100 * r$cov90, r$width90))
}
writeLines(md, file.path(OUT, "table3_transfer.md"))
cat(md, sep = "\n")
