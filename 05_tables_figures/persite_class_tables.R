# =====================================================================
# persite_class_tables.R  --  accuracy per site, per canopy height class
#                             and per forest class (18 sites)
#
# Produces the per-site table behind Table 1 (footprints, median canopy
# height), Table S4 (checks) and Fig. 4, and the pooled class tables that
# Tables S3 and S4 are checked against. Errors per footprint, from the
# analysis frame (630,390 footprints at 18 sites):
#   CHM error = P3D CHM disc mean - rebuilt pit-free lidar CHM disc mean
#   DTM error = P3D DTM disc mean - rebuilt lidar ground disc mean
#   DSM error = P3D DSM disc mean - (lidar ground + lidar CHM), sampled
#               directly (p3d_dsm_direct.R, merge_dsm_direct.R)
# Sites are tracker numbers (1-9, 11-15, 17-20); site_names.csv gives the
# paper numbers.
#
# Inputs  (P3D_DATA_ROOT) rework/frame_footprints.csv.gz,
#                         lidar_products/lp_summary.csv
# Outputs (P3D_RESULTS)   persite/
#   persite.csv            per site, three products
#   height_class.csv       by lidar canopy height class (2-5 ... >30 m), pooled
#   forest_class.csv       by forest class, pooled
#   accuracy_summary.csv   three products: bias ... Q95, % |err| > 5 and 10 m;
#                          pooled, site-median, and pooled without site 1 rows
#   lp_checks.csv          per site lidar reference checks from lp_summary.csv
#   error_densities.png / .pdf   overview of the error densities (not in the paper)
#   persite_class_tables_log.txt
# =====================================================================

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

DATA_ROOT <- Sys.getenv("P3D_DATA_ROOT", "<LOCAL_DATA_ROOT>")
OUT_ROOT  <- Sys.getenv("P3D_RESULTS", "results")
OUT_DIR   <- file.path(OUT_ROOT, "persite")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
LOG <- file.path(OUT_DIR, "persite_class_tables_log.txt"); if (file.exists(LOG)) file.remove(LOG)
say <- function(...) { s <- sprintf(...); cat(s, "\n"); cat(s, "\n", file = LOG, append = TRUE) }

FOREST <- c("ENF", "BDF", "EBF", "DNF")
HCLASS <- c(2, 5, 10, 15, 20, 30, Inf)

F <- fread(file.path(DATA_ROOT, "rework", "frame_footprints.csv.gz"),
           colClasses = list(character = c("shot_number", "lc_l1_code")))
setnames(F, c("chm_error_mean", "dtm_error_mean", "dsm_error_mean"), c("chm_err", "dtm_err", "dsm_err"))
say("per-site and class accuracy tables, run %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
say("frame: %s footprints at %d sites; site 1 holds %s (%.1f%%)", format(nrow(F), big.mark = ","),
    uniqueN(F$site), format(F[site == 1, .N], big.mark = ","), 100 * F[, mean(site == 1)])
# not zero: the DSM error is sampled directly, while the delivered CHM is the DSM minus the DTM set to 0 where negative
say("identity max |dsm - chm - dtm| %.1e", F[, max(abs(dsm_err - chm_err - dtm_err))])

err_stats <- function(x) {
  x <- x[is.finite(x)]; n <- length(x); m <- mean(x); d <- x - m
  m2 <- mean(d^2); m3 <- mean(d^3); m4 <- mean(d^4)
  data.table(n = n, mean = m, median = median(x), sd = sd(x),
             nmad = 1.4826 * median(abs(x - median(x))), iqr = IQR(x),
             skew = m3 / m2^1.5, kurt = m4 / m2^2 - 3, mae = mean(abs(x)), rmse = sqrt(mean(x^2)),
             q95 = as.numeric(quantile(abs(x), 0.95)),
             pct_over5 = 100 * mean(abs(x) > 5), pct_over10 = 100 * mean(abs(x) > 10))
}
three <- function(d, ...) rbindlist(list(cbind(..., product = "CHM", err_stats(d$chm_err)),
                                         cbind(..., product = "DSM", err_stats(d$dsm_err)),
                                         cbind(..., product = "DTM", err_stats(d$dtm_err))))

# ---- accuracy summary (pooled, site median, and pooled without site 1)
ps <- F[, .(chm_mean = mean(chm_err), chm_sd = sd(chm_err), chm_rmse = sqrt(mean(chm_err^2)),
            chm_o5 = 100 * mean(abs(chm_err) > 5), chm_o10 = 100 * mean(abs(chm_err) > 10),
            dsm_mean = mean(dsm_err), dsm_sd = sd(dsm_err), dsm_rmse = sqrt(mean(dsm_err^2)),
            dsm_o5 = 100 * mean(abs(dsm_err) > 5), dsm_o10 = 100 * mean(abs(dsm_err) > 10),
            dtm_mean = mean(dtm_err), dtm_sd = sd(dtm_err), dtm_rmse = sqrt(mean(dtm_err^2)),
            dtm_o5 = 100 * mean(abs(dtm_err) > 5), dtm_o10 = 100 * mean(abs(dtm_err) > 10)), by = site]
sm <- rbindlist(lapply(c("chm", "dsm", "dtm"), function(p)
  data.table(frame = "site median", product = toupper(p), n = nrow(ps),
             mean = median(ps[[paste0(p, "_mean")]]), sd = median(ps[[paste0(p, "_sd")]]),
             rmse = median(ps[[paste0(p, "_rmse")]]), pct_over5 = median(ps[[paste0(p, "_o5")]]),
             pct_over10 = median(ps[[paste0(p, "_o10")]]))))
t3 <- rbind(three(F, frame = "pooled"), sm, three(F[site != 1], frame = "pooled without site 1"), fill = TRUE)
fwrite(t3, file.path(OUT_DIR, "accuracy_summary.csv"))
fmt <- function(x, f) ifelse(is.na(x), "-", sprintf(f, x))
say("\nACCURACY SUMMARY (18 sites, rebuilt lidar references)")
for (i in seq_len(nrow(t3))) with(t3[i], say("%-22s %s  n %9s  bias %+.2f  SD %.2f  RMSE %.2f  NMAD %s  Q95 %s  >5 m %s  >10 m %s",
  frame, product, format(n, big.mark = ","), mean, sd, rmse, fmt(nmad, "%.2f"), fmt(q95, "%.2f"),
  fmt(pct_over5, "%.1f%%"), fmt(pct_over10, "%.1f%%")))

# ---- height classes
F[, hclass := cut(chm_p90, HCLASS, right = FALSE, labels = c("2-5", "5-10", "10-15", "15-20", "20-30", ">30"))]
hc <- rbindlist(lapply(levels(F$hclass), function(h) three(F[hclass == h], height_class = h)))
hc_sites <- F[, .(n_sites = uniqueN(site), site1_share = mean(site == 1)), by = hclass]
hc[hc_sites, `:=`(n_sites = i.n_sites, site1_share = i.site1_share), on = c(height_class = "hclass")]
fwrite(hc, file.path(OUT_DIR, "height_class.csv"))
say("\nBY CANOPY HEIGHT CLASS (rebuilt lidar p90): bias / RMSE / %%>5 m")
for (h in levels(F$hclass)) {
  x <- hc[height_class == h]
  say("%-6s n %7d (%2d sites, site 1 %2.0f%%) | CHM %+.2f / %.2f / %4.1f%% | DSM %+.2f / %.2f / %4.1f%% | DTM %+.2f / %.2f / %4.1f%%",
      h, x$n[1], x$n_sites[1], 100 * x$site1_share[1],
      x[product == "CHM"]$mean, x[product == "CHM"]$rmse, x[product == "CHM"]$pct_over5,
      x[product == "DSM"]$mean, x[product == "DSM"]$rmse, x[product == "DSM"]$pct_over5,
      x[product == "DTM"]$mean, x[product == "DTM"]$rmse, x[product == "DTM"]$pct_over5)
}

# ---- per site
persite <- F[, {
  vc <- var(chm_err); vs <- var(dsm_err)
  .(n = .N, chm_bias = mean(chm_err), chm_sd = sd(chm_err), chm_rmse = sqrt(mean(chm_err^2)),
    dsm_bias = mean(dsm_err), dsm_sd = sd(dsm_err), dsm_rmse = sqrt(mean(dsm_err^2)),
    dtm_bias = mean(dtm_err), dtm_sd = sd(dtm_err), dtm_rmse = sqrt(mean(dtm_err^2)),
    chm_pct_over5 = 100 * mean(abs(chm_err) > 5), median_height = median(chm_p90),
    r_dsm_dtm = cor(dsm_err, dtm_err), share_removed_by_perfect_dtm = 1 - vs / vc,
    forest_class_share = mean(lc_l1_code %in% FOREST))
}, by = site][order(site)]
fwrite(persite, file.path(OUT_DIR, "persite.csv"))
say("\nPER SITE (tracker): n | CHM bias/SD | DSM bias/SD | DTM bias/SD | CHM %%>5 m | median height | perfect-DTM share")
for (i in seq_len(nrow(persite))) with(persite[i], say("%4d %7d | %+6.2f/%5.2f | %+6.2f/%5.2f | %+6.2f/%5.2f | %5.1f%% | %5.1f m | %+6.1f%%",
  site, n, chm_bias, chm_sd, dsm_bias, dsm_sd, dtm_bias, dtm_sd, chm_pct_over5, median_height,
  100 * share_removed_by_perfect_dtm))

# ---- forest classes
tf <- rbindlist(lapply(c("forest classes (4)", FOREST, "non-forest with canopy"), function(nm) {
  d <- if (nm == "forest classes (4)") F[lc_l1_code %in% FOREST] else if (nm == "non-forest with canopy") F[!lc_l1_code %in% FOREST] else F[lc_l1_code == nm]
  if (!nrow(d)) return(NULL)
  three(d, stratum = nm)
}))
fwrite(tf, file.path(OUT_DIR, "forest_class.csv"))

# ---- reference checks (lp_summary.csv)
lp <- fread(file.path(DATA_ROOT, "lidar_products", "lp_summary.csv"))
lpc <- lp[, .(site = tracker, flight = code, gap_3dep_med, gap_3dep_within05, chm_mean_minus_delivered,
              r_chm_delivered, p90_minus_gedi, r_p90_gedi, r_delivered_gedi)][order(site)]
fwrite(lpc, file.path(OUT_DIR, "lp_checks.csv"))
say("\nREFERENCE CHECKS (lp_summary; ground gap is raw, before the prep's datum fixes at sites 3 and 19)")
say("site flight              | ground-3DEP med, within 0.5 m | rebuilt-delivered CHM mean, r | rebuilt p90-GEDI, r (delivered r)")
for (i in seq_len(nrow(lpc))) with(lpc[i], say("%4d %-19s | %+7.2f %4.2f | %+6.2f %5.3f | %+6.2f %5.3f (%5.3f)",
  site, flight, gap_3dep_med, gap_3dep_within05, chm_mean_minus_delivered, r_chm_delivered, p90_minus_gedi, r_p90_gedi, r_delivered_gedi))

# ---- overview figure (not in the paper): error densities, pooled and each site weighted equally
SERIES <- c(CHM = "#2a78d6", DSM = "#eb6834", DTM = "#1baf7a")
XL <- c(-20, 20)
long <- melt(F[, .(site, CHM = chm_err, DSM = dsm_err, DTM = dtm_err)], id.vars = "site",
             variable.name = "product", value.name = "err")
long[, w_site := 1 / .N, by = .(site, product)]
dens <- function(d, w, panel) {
  rbindlist(lapply(names(SERIES), function(p) {
    x <- d[product == p]; ww <- if (w) x$w_site / sum(x$w_site) else rep(1 / nrow(x), nrow(x))
    k <- density(x$err, weights = ww, bw = 0.3, from = XL[1], to = XL[2], n = 801)
    data.table(panel = panel, product = p, x = k$x, y = k$y)
  }))
}
dd <- rbind(dens(long, FALSE, "All footprints"), dens(long, TRUE, "Each site weighted equally"))
dd[, product := factor(product, levels = names(SERIES))]
lab <- dd[, .SD[which.max(y)], by = .(panel, product)]
g <- ggplot(dd, aes(x, y, colour = product)) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.4) +
  geom_line(linewidth = 0.8) +
  geom_text(data = lab, aes(label = product), vjust = -0.6, size = 3.4, colour = "grey15", show.legend = FALSE) +
  geom_point(data = lab, size = 1.6, show.legend = FALSE) +
  facet_wrap(~panel, nrow = 1) +
  scale_colour_manual(values = SERIES, name = NULL) +
  scale_x_continuous(breaks = seq(-20, 20, 5), expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Error, P3D minus lidar (m)", y = "Density") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
        strip.text = element_text(face = "bold", hjust = 0), legend.position = "top",
        legend.justification = "left", panel.spacing = unit(1.2, "lines"))
ggsave(file.path(OUT_DIR, "error_densities.png"), g, width = 7.2, height = 3.2, dpi = 300, bg = "white")
ggsave(file.path(OUT_DIR, "error_densities.pdf"), g, width = 7.2, height = 3.2)
say("\nshare of errors outside the plotted -20 to +20 m: CHM %.2f%%, DSM %.2f%%, DTM %.2f%%",
    100 * F[, mean(abs(chm_err) > 20)], 100 * F[, mean(abs(dsm_err) > 20)], 100 * F[, mean(abs(dtm_err) > 20)])
say("\noutputs in %s", OUT_DIR)
