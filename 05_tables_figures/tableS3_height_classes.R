# tableS3_height_classes.R  --  Table S3 class means: mean error by canopy height class, lidar and GEDI heights
# Same rule as Fig. 3b: site means per class, averaged over sites with >= 30 footprints in the class,
# so each site counts equally. Classes by lidar canopy height (90th percentile, as Fig. 3b) and by
# GEDI RH98 (Table S3b).
# Inputs  (P3D_DATA_ROOT) rework/frame_footprints.csv.gz
#         (P3D_RESULTS)   fig3/fig3b_class_means.csv (fig3_error_distributions.R, for the check)
# Outputs (P3D_RESULTS)   tableS3/tableS3_height_classes.{csv,md}
# Checks  lidar-class means reproduce the Fig. 3b means; footprints below 2 m GEDI height counted.
suppressPackageStartupMessages(library(data.table))
DATA_ROOT <- Sys.getenv("P3D_DATA_ROOT", "<LOCAL_DATA_ROOT>")
OUT_ROOT  <- Sys.getenv("P3D_RESULTS", "results")
OUT  <- file.path(OUT_ROOT, "tableS3"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
fr <- fread(file.path(DATA_ROOT, "rework", "frame_footprints.csv.gz"),
            select = c("site", "chm_p90", "rh_98", "dsm_error_mean", "dtm_error_mean", "chm_error_mean"))
stopifnot(nrow(fr) == 630390L)
PROD <- c("DSM", "DTM", "CHM")
br <- c(2, 5, 10, 15, 20, 30, Inf); lab <- c("2\u20135", "5\u201310", "10\u201315", "15\u201320", "20\u201330", ">30")

L <- melt(fr, id.vars = c("site", "chm_p90", "rh_98"),
          measure.vars = c("dsm_error_mean", "dtm_error_mean", "chm_error_mean"),
          variable.name = "product", value.name = "err")
L[, product := factor(toupper(sub("_error_mean", "", product)), levels = PROD)]
cm <- function(h) {
  L[, hc := cut(get(h), br, lab, right = FALSE)]
  S <- L[!is.na(hc), .(m = mean(err), n = .N), by = .(product, hc, site)][n >= 30]
  S[, .(mean_error = mean(m), n_sites = .N), by = .(product, hc)][order(product, hc)]
}
Lid <- cm("chm_p90")[, height := "Lidar"]; Ged <- cm("rh_98")[, height := "GEDI"]

f3b <- fread(file.path(OUT_ROOT, "fig3", "fig3b_class_means.csv"), encoding = "UTF-8")[, .(product, hc, mean_error, n_sites)]
g <- merge(Lid, f3b, by = c("product", "hc"))
d <- g[, max(abs(mean_error.x - mean_error.y))]; k <- g[, all(n_sites.x == n_sites.y)]
cat(sprintf("GATE lidar classes reproduce Fig. 3b: max |diff| = %.1e, site counts equal %s  %s\n",
            d, k, ifelse(nrow(g) == 18L && d < 1e-9 && k, "PASS", "FAIL")))
cat(sprintf("footprints with GEDI RH98 below 2 m (not in any GEDI class): %s (%.1f%%)\n",
            format(fr[rh_98 < 2, .N], big.mark = ","), 100 * fr[, mean(rh_98 < 2)]))
cat("GEDI-class ranges: "); print(Ged[, .(lo = min(mean_error), hi = max(mean_error)), by = product])

T3 <- rbind(Lid, Ged)
fwrite(T3, file.path(OUT, "tableS3_height_classes.csv"))
f2 <- function(x) sub("-", "\u2212", sprintf("%+.1f", x), fixed = TRUE)
stopifnot(T3[, uniqueN(n_sites), by = .(height, hc)][, all(V1 == 1L)])
W <- dcast(T3, hc ~ height + product, value.var = "mean_error")
N <- dcast(T3[product == "DSM"], hc ~ height, value.var = "n_sites")
W <- merge(W, N, by = "hc")
md <- c("| Canopy height class (m) | DSM, lidar classes | DTM, lidar classes | CHM, lidar classes | Sites | DSM, GEDI classes | DTM, GEDI classes | CHM, GEDI classes | Sites |",
        "|---|---|---|---|---|---|---|---|---|",
        W[, sprintf("| %s | %s | %s | %s | %d | %s | %s | %s | %d |", hc, f2(Lidar_DSM), f2(Lidar_DTM), f2(Lidar_CHM), Lidar,
                    f2(GEDI_DSM), f2(GEDI_DTM), f2(GEDI_CHM), GEDI)])
writeLines(md, file.path(OUT, "tableS3_height_classes.md")); cat(md, sep = "\n")
