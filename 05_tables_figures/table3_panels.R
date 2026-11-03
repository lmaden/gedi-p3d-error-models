# table3_panels.R  --  Table 3 as printed, two panels
# (a) shares of error variance in the models fit to all training footprints, whole percentages, with the unexplained share
# (b) RMSE at held-out sites, interval coverage and width, each held-out site weighted equally
#     (frame "Sites weighted equally" of table3_transfer.R)
# Input  (P3D_RESULTS) table3/table3_transfer.csv (table3_transfer.R)
# Output (P3D_RESULTS) table3/table3_panels.md
suppressPackageStartupMessages(library(data.table))
OUT_ROOT <- Sys.getenv("P3D_RESULTS", "results")
TT <- fread(file.path(OUT_ROOT, "table3", "table3_transfer.csv"))
t3 <- TT[frame == "All held-out footprints"]   # panel (a) R2 shares (full fits)
t3[, product := factor(product, levels = c("DSM", "DTM", "CHM"))]; setorder(t3, product)
stopifnot(nrow(t3) == 3L, all(abs(t3$r2_total - t3$r2_fixed - t3$r2_groups) < 1e-9))
# panel (b) weights each held-out site equally, as in Table 2
tb <- TT[frame == "Sites weighted equally"]
tb[, product := factor(product, levels = c("DSM", "DTM", "CHM"))]; setorder(tb, product)
stopifnot(nrow(tb) == 3L)
pc <- function(x) sprintf("%.0f", 100 * x); f2 <- function(x) sprintf("%.2f", x)
md <- c("*(a) Share of the error variance in the models fit at all 18 sites (%)*", "",
        "| Product | Tested variables | Site, ecoregion and land-cover offsets | Unexplained |",
        "|---|---|---|---|",
        t3[, sprintf("| %s | %s | %s | %s |", product, pc(r2_fixed), pc(r2_groups), pc(1 - r2_total))], "",
        "*(b) Correction at held-out sites*", "",
        "| Product | RMSE, no correction (m) | RMSE, single average correction (m) | RMSE, model correction (m) | RMSE, site's own mean removed (m) | True error inside the model's nominal 90% prediction interval (%) | Mean width of the model's nominal 90% prediction interval (m) |",
        "|---|---|---|---|---|---|---|",
        tb[, sprintf("| %s | %s | %s | %s | %s | %s | %.1f |", product, f2(rmse_none), f2(rmse_global), f2(rmse_model), f2(rmse_own), pc(cov90), width90)])
writeLines(md, file.path(OUT_ROOT, "table3", "table3_panels.md")); cat(md, sep = "\n")
