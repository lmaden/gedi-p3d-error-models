# tableS4_persite.R  --  Table S4 site rows, accuracy of the three products at each site
# Columns per product: bias, RMSE, % of footprints with |error| > 5 m (the Fig. 4 statistics,
# extended to all three products; footprints and canopy height are in Table 1).
# Inputs  (P3D_DATA_ROOT) rework/frame_footprints.csv.gz
#         (P3D_DERIVED)   site_names.csv
#         (P3D_RESULTS)   persite/persite.csv (persite_class_tables.R, for the check)
# Outputs (P3D_RESULTS)   tableS4/tableS4_persite.{csv,md}
# Checks  bias and RMSE equal persite/persite.csv (source of Fig. 4) to 1e-9;
#         CHM bias = DSM bias - DTM bias at every site within 0.1 m (Eq. 4); the counts the text cites.
suppressPackageStartupMessages(library(data.table))
DATA_ROOT <- Sys.getenv("P3D_DATA_ROOT", "<LOCAL_DATA_ROOT>")
DERIVED   <- Sys.getenv("P3D_DERIVED", "data_derived")
OUT_ROOT  <- Sys.getenv("P3D_RESULTS", "results")
OUT  <- file.path(OUT_ROOT, "tableS4"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
fr <- fread(file.path(DATA_ROOT, "rework", "frame_footprints.csv.gz"),
            select = c("site", "dsm_error_mean", "dtm_error_mean", "chm_error_mean"))
nm <- fread(file.path(DERIVED, "site_names.csv"))
ps <- fread(file.path(OUT_ROOT, "persite", "persite.csv"))
stopifnot(nrow(fr) == 630390L)

st <- function(e) list(bias = mean(e), rmse = sqrt(mean(e^2)), over5 = 100 * mean(abs(e) > 5))
S <- fr[, c(setNames(st(dsm_error_mean), paste0("dsm_", c("bias", "rmse", "over5"))),
            setNames(st(dtm_error_mean), paste0("dtm_", c("bias", "rmse", "over5"))),
            setNames(st(chm_error_mean), paste0("chm_", c("bias", "rmse", "over5")))), by = site]
S <- merge(nm[, .(site = tracker, paper_site, short_name)], S, by = "site")
setorder(S, paper_site)

# checks
g <- merge(S, ps, by = "site")
d <- g[, max(abs(c(dsm_bias.x - dsm_bias.y, dsm_rmse.x - dsm_rmse.y, dtm_bias.x - dtm_bias.y,
                   dtm_rmse.x - dtm_rmse.y, chm_bias.x - chm_bias.y, chm_rmse.x - chm_rmse.y,
                   chm_over5 - chm_pct_over5)))]
cat(sprintf("GATE vs persite.csv (Fig. 4 source), max |diff| = %.1e  %s\n", d, ifelse(d < 1e-9, "PASS", "FAIL")))
e3 <- S[, max(abs(chm_bias - (dsm_bias - dtm_bias)))]
cat(sprintf("GATE CHM bias = DSM bias - DTM bias (Eq. 4) within 0.1 m (the DSM is sampled directly and the delivered CHM is 0 where the DSM is below the DTM), max |diff| = %.1e  %s\n", e3, ifelse(e3 < 0.1, "PASS", "FAIL")))
cat("text checks: DSM RMSE < CHM RMSE at", S[, sum(dsm_rmse < chm_rmse)], "of 18 (text 14);",
    "DSM RMSE lowest of the three at", S[, sum(dsm_rmse < pmin(dtm_rmse, chm_rmse))], "sites\n")
cat("sites where DSM RMSE >= CHM RMSE:", S[dsm_rmse >= chm_rmse, paste(paper_site, short_name, collapse = "; ")], "\n")

fwrite(S, file.path(OUT, "tableS4_persite.csv"))
f2 <- function(x) sub("-", "\u2212", sprintf("%.2f", x), fixed = TRUE)
md <- c("| Site | DSM bias (m) | DSM RMSE (m) | DSM > 5 m (%) | DTM bias (m) | DTM RMSE (m) | DTM > 5 m (%) | CHM bias (m) | CHM RMSE (m) | CHM > 5 m (%) |",
        "|---|---|---|---|---|---|---|---|---|---|",
        S[, sprintf("| %d %s | %s | %s | %.1f | %s | %s | %.1f | %s | %s | %.1f |", paper_site, short_name,
                    f2(dsm_bias), f2(dsm_rmse), dsm_over5, f2(dtm_bias), f2(dtm_rmse), dtm_over5,
                    f2(chm_bias), f2(chm_rmse), chm_over5)])
writeLines(md, file.path(OUT, "tableS4_persite.md"))
cat(md, sep = "\n")
