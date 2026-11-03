# tableS2_transfer_groups.R  --  Table S2, correction at each held-out site, with the groups held out together
# One table, three panels (DSM, DTM, CHM): per held-out site, its group, footprints, RMSE with no correction,
# after the single average correction, after the model's correction, after removing the site's own mean
# error, and the share of true errors inside the model's 90% prediction interval. DSM cells use the directly
# sampled DSM error (dsm_direct_rowpreds.R).
# Inputs  (P3D_DATA_ROOT) rework/folds/rowpreds_{dsm,dtm,chm}_fold{1..6}.csv.gz (group = fold number in the
#                         file name); dsm_direct/dsm_direct_frame.csv.gz (through dsm_direct_rowpreds.R)
#         (P3D_DERIVED)   site_names.csv
#         (P3D_RESULTS)   table3/table3_transfer.csv (table3_transfer.R, for the check)
# Outputs (P3D_RESULTS)   tableS2/tableS2_transfer_groups.{csv,md}
# Checks  no-correction pooled RMSE equals Table 3 rmse_none; medians of no-correction equal the Table 3
#         median-site rmse_none.
suppressPackageStartupMessages(library(data.table))
DATA_ROOT <- Sys.getenv("P3D_DATA_ROOT", "<LOCAL_DATA_ROOT>")
DERIVED   <- Sys.getenv("P3D_DERIVED", "data_derived")
OUT_ROOT  <- Sys.getenv("P3D_RESULTS", "results")
OUT  <- file.path(OUT_ROOT, "tableS2"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
ARGF <- grep("^--file=", commandArgs(FALSE), value = TRUE)   # folder of this script (Rscript), for source()
HERE <- if (length(ARGF)) dirname(normalizePath(sub("^--file=", "", ARGF[1]))) else "05_tables_figures"
nm   <- fread(file.path(DERIVED, "site_names.csv"))
rp <- list.files(file.path(DATA_ROOT, "rework", "folds"), "^rowpreds_(dsm|dtm|chm)_fold[0-9]+\\.csv\\.gz$", full.names = TRUE)
stopifnot(length(rp) == 18L)
R <- rbindlist(lapply(rp, function(f) {
  x <- fread(f, colClasses = list(character = c("site", "shot_number")))
  x[, product := toupper(sub("^rowpreds_(dsm|dtm|chm)_fold.*", "\\1", basename(f)))]
  x[, group := as.integer(sub(".*_fold([0-9]+)\\.csv\\.gz$", "\\1", basename(f)))]
}), fill = TRUE)
# DSM rows switched to the directly sampled DSM error
source(file.path(HERE, "dsm_direct_rowpreds.R"))
R <- use_direct_dsm(R, "group", DATA_ROOT)
rmse <- function(p, o) sqrt(mean((p - o)^2))
S <- R[, .(n = .N, rmse_none = sqrt(mean(observed^2)), rmse_global = rmse(pred_const_global, observed),
           rmse_model = rmse(pred_new_site, observed), rmse_own = rmse(pred_const_site, observed),
           cov90 = mean(observed >= lo90 & observed <= hi90)), by = .(product, group, site)]
S[, paper_site := nm$paper_site[match(as.integer(site), nm$tracker)]]
S[, short_name := nm$short_name[match(as.integer(site), nm$tracker)]]
stopifnot(!anyNA(S$paper_site), S[, .N, by = product]$N == 18L)

t3 <- fread(file.path(OUT_ROOT, "table3", "table3_transfer.csv"))
pool <- R[, .(none = sqrt(mean(observed^2))), by = product]
g3 <- merge(pool, t3[frame == "All held-out footprints"], by = "product")[, max(abs(none - rmse_none))]
med <- merge(S[, .(m = median(rmse_none)), by = product], t3[frame == "Median site"], by = "product")[, max(abs(m - rmse_none))]
cat(sprintf("GATE no-correction pooled = Table 3 (%.1e), median site = Table 3 (%.1e)  %s\n", g3, med, ifelse(g3 < 1e-9 & med < 1e-9, "PASS", "FAIL")))
cat("CHM: average raised RMSE at", S[product == "CHM" & rmse_global > rmse_none, paste(paper_site, short_name)], "\n")
cat("CHM: model raised RMSE at", S[product == "CHM" & rmse_model > rmse_none, paste(paper_site, short_name)], "\n")
seven <- S[product == "CHM"][order(rmse_none)][1:7]
cat("seven smallest CHM no-correction:", seven$paper_site, "| average raised", seven[rmse_global > rmse_none, .N], "| model raised", seven[rmse_model > rmse_none, .N], "\n")
setorder(S, product, paper_site); fwrite(S, file.path(OUT, "tableS2_transfer_groups.csv"))
f2 <- function(x) sprintf("%.2f", x); pc <- function(x) sprintf("%.0f", 100 * x)
bf <- function(txt, flag) ifelse(flag, paste0("**", txt, "**"), txt)
md <- character()
for (p in c("DSM", "DTM", "CHM")) {
  first <- p == "DSM"
  hdr <- if (first) "| Site | Held-out group | Held-out footprints | RMSE, no correction (m) | RMSE, single average correction (m) | RMSE, model correction (m) | RMSE, site's own mean removed (m) | True error inside the model's nominal 90% prediction interval (%) |" else
                    "| Site | RMSE, no correction (m) | RMSE, single average correction (m) | RMSE, model correction (m) | RMSE, site's own mean removed (m) | True error inside the model's nominal 90% prediction interval (%) |"
  sep <- if (first) "|---|---|---|---|---|---|---|---|" else "|---|---|---|---|---|---|"
  rows <- S[product == p, {
    a <- bf(f2(rmse_global), round(rmse_global, 2) > round(rmse_none, 2)); m <- bf(f2(rmse_model), round(rmse_model, 2) > round(rmse_none, 2))
    cv <- bf(pc(cov90), cov90 < 0.75)
    if (first) sprintf("| %d %s | %d | %s | %s | %s | %s | %s | %s |", paper_site, short_name, group, format(n, big.mark = ",", trim = TRUE), f2(rmse_none), a, m, f2(rmse_own), cv)
    else sprintf("| %d %s | %s | %s | %s | %s | %s |", paper_site, short_name, f2(rmse_none), a, m, f2(rmse_own), cv) }]
  md <- c(md, sprintf("**(%s) %s**", c(DSM = "a", DTM = "b", CHM = "c")[p], p), "", hdr, sep, rows, "")
}
writeLines(md, file.path(OUT, "tableS2_transfer_groups.md")); cat(md, sep = "\n")
