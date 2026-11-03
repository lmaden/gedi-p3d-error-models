# tableS3_rmse_over5.R  --  Table S3a RMSE and share of footprints with |error| > 5 m by lidar canopy height class
# Same rule as Fig. 3b and the Table S3 class means: statistic per site and class, sites with >= 30 footprints in
# the class, then the mean over sites so each site counts equally. rmse_rms is the root of the mean of the site MSEs
# (the Table 2 rule, printed in Table S3a); rmse is the mean of the site RMSEs (for reference). Also writes the pooled
# values for comparison with the pooled class table of persite_class_tables.R.
# Inputs  (P3D_DATA_ROOT) rework/frame_footprints.csv.gz
#         (P3D_RESULTS)   persite/height_class.csv (persite_class_tables.R) and
#                         tableS3/tableS3_height_classes.csv (tableS3_height_classes.R), for the checks
# Output  (P3D_RESULTS)   tableS3/tableS3_rmse_over5.csv
# Checks  pooled values reproduce persite/height_class.csv; site-weighted means and site counts reproduce
#         the Table S3 class means.
suppressPackageStartupMessages(library(data.table))
DATA_ROOT <- Sys.getenv("P3D_DATA_ROOT", "<LOCAL_DATA_ROOT>")
OUT_ROOT  <- Sys.getenv("P3D_RESULTS", "results")
OUT  <- file.path(OUT_ROOT, "tableS3"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
fr <- fread(file.path(DATA_ROOT, "rework", "frame_footprints.csv.gz"),
            select = c("site", "chm_p90", "dsm_error_mean", "dtm_error_mean", "chm_error_mean"))
stopifnot(nrow(fr) == 630390L)
br <- c(2, 5, 10, 15, 20, 30, Inf); lab <- c("2-5", "5-10", "10-15", "15-20", "20-30", ">30")
L <- melt(fr, id.vars = c("site", "chm_p90"), variable.name = "product", value.name = "err")
L[, product := toupper(sub("_error_mean", "", product))]
L[, hc := cut(chm_p90, br, lab, right = FALSE)]
S <- L[!is.na(hc), .(m = mean(err), rmse = sqrt(mean(err^2)), over5 = 100 * mean(abs(err) > 5), n = .N),
       by = .(product, hc, site)][n >= 30]
E <- S[, .(mean_error = mean(m), rmse = mean(rmse), rmse_rms = sqrt(mean(rmse^2)), over5 = mean(over5), n_sites = .N), by = .(product, hc)]
P <- L[!is.na(hc), .(pooled_rmse = sqrt(mean(err^2)), pooled_over5 = 100 * mean(abs(err) > 5), n = .N), by = .(product, hc)]
R <- merge(E, P, by = c("product", "hc"))[order(product, hc)]
# check against the pooled class table of persite_class_tables.R
pcl <- fread(file.path(OUT_ROOT, "persite", "height_class.csv"))
pcl[, hc := sub("-", "-", height_class)]
g <- merge(R, pcl[, .(product, hc = height_class, rmse_pcl = rmse, over5_pcl = pct_over5)], by = c("product", "hc"))
dmax <- g[, max(abs(pooled_rmse - rmse_pcl), abs(pooled_over5 - over5_pcl))]
cat(sprintf("GATE pooled values reproduce persite/height_class.csv: max |diff| = %.2e  %s\n", dmax, ifelse(nrow(g) == 18 && dmax < 1e-6, "PASS", "FAIL")))
# check against the Table S3 lidar-class means
s3 <- fread(file.path(OUT_ROOT, "tableS3", "tableS3_height_classes.csv"))[height == "Lidar"]
s3[, hc := gsub("\u2013", "-", hc)]
g2 <- merge(R, s3[, .(product, hc, me_s3 = mean_error, ns_s3 = n_sites)], by = c("product", "hc"))
cat(sprintf("GATE site-weighted means reproduce Table S3: max |diff| = %.2e, sites equal %s\n",
            g2[, max(abs(mean_error - me_s3))], g2[, all(n_sites == ns_s3)]))
fwrite(R, file.path(OUT, "tableS3_rmse_over5.csv"))
print(R[, .(product, hc, n_sites, mean_error = round(mean_error, 1), rmse = round(rmse, 1), rmse_rms = round(rmse_rms, 1), over5 = round(over5),
            pooled_rmse = round(pooled_rmse, 1), pooled_over5 = round(pooled_over5))])
