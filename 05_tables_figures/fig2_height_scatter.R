# fig2_height_scatter.R  --  Fig. 2, P3D canopy height against lidar canopy height
# Lidar canopy height = footprint mean of the lidar CHM (H_L in Eq. 3), all 630,390 footprints, hexagonal bins on a
# log count scale (raw footprint counts), dashed 1:1 line. Panel statistics weight each site equally (weight 1/n_site).
# (a) P3D CHM = lidar CHM mean + CHM error. (b) P3D DSM minus the lidar ground = lidar CHM mean + DSM error, with the
# directly sampled DSM error.
# Same join as the Fig. 2 panel statistics of site_weighted_accuracy.py.
# Inputs  (P3D_DATA_ROOT) rework/frame_footprints.csv.gz, lidar_products/lp_site_NN.csv.gz (chm_mean)
# Outputs (P3D_RESULTS)   fig2/fig2_height_scatter.{png,pdf}, fig2/fig2_scatter_stats.csv (panel statistics and the
#                         number of footprints outside the plotted range)
# Checks  panel statistics equal the published values (stopifnot below).
# SCATTER_PAL=viridis draws a viridis version (files with the suffix _viridis); the default is spectral.
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
DATA_ROOT <- Sys.getenv("P3D_DATA_ROOT", "<LOCAL_DATA_ROOT>")
OUT_ROOT  <- Sys.getenv("P3D_RESULTS", "results")
OUT  <- file.path(OUT_ROOT, "fig2"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
fr <- fread(file.path(DATA_ROOT, "rework", "frame_footprints.csv.gz"),
            colClasses = list(character = "shot_number"))
lp <- rbindlist(lapply(Sys.glob(file.path(DATA_ROOT, "lidar_products", "lp_site_*.csv.gz")),
                       function(f) fread(f, select = c("shot_number", "chm_mean"), colClasses = list(character = "shot_number"))))
d <- merge(fr[, .(shot_number, site, chm_error_mean, dsm_error_mean)], lp[!is.na(chm_mean)], by = "shot_number")
stopifnot(nrow(d) == 630390L, uniqueN(d$shot_number) == 630390L)
d[, `:=`(p3d_chm = chm_mean + chm_error_mean, p3d_dsm_lg = chm_mean + dsm_error_mean)]
d[, w := 1 / .N, by = site]
st <- function(y, x, w) { e <- y - x; mx <- weighted.mean(x, w); my <- weighted.mean(y, w)
  data.table(n = length(x), bias = weighted.mean(e, w), rmse = sqrt(weighted.mean(e^2, w)),
             r = weighted.mean((x - mx) * (y - my), w) / sqrt(weighted.mean((x - mx)^2, w) * weighted.mean((y - my)^2, w))) }
S <- rbind(cbind(panel = "a", st(d$p3d_chm, d$chm_mean, d$w)), cbind(panel = "b", st(d$p3d_dsm_lg, d$chm_mean, d$w)))
LIM <- c(-5, 40)
S[, outside := c(d[, sum(chm_mean < LIM[1] | chm_mean > LIM[2] | p3d_chm < LIM[1] | p3d_chm > LIM[2])],
                 d[, sum(chm_mean < LIM[1] | chm_mean > LIM[2] | p3d_dsm_lg < LIM[1] | p3d_dsm_lg > LIM[2])])]
# Published values, sites weighted equally (also written by site_weighted_accuracy.py to table2/fig2_panel_stats.csv):
# CHM -2.918/5.365/0.796; DSM minus lidar ground -1.261/3.632/0.898 (directly sampled DSM error)
stopifnot(abs(S$bias[1] + 2.918) < 0.001, abs(S$rmse[1] - 5.365) < 0.001, abs(S$r[1] - 0.796) < 0.001,
          abs(S$bias[2] + 1.261) < 0.001, abs(S$rmse[2] - 3.632) < 0.001, abs(S$r[2] - 0.898) < 0.001)
fwrite(S, file.path(OUT, "fig2_scatter_stats.csv")); print(S)

minus <- function(s) gsub("-", "\u2212", s)
lab <- function(i) minus(sprintf("Bias %.2f m, RMSE %.2f m, r = %.2f", S$bias[i], S$rmse[i], S$r[i]))
INK2 <- "#52514e"
PAL <- Sys.getenv("SCATTER_PAL", "spectral")  # the published figure uses spectral; viridis is an option
FILL <- switch(PAL,
  viridis  = scale_fill_viridis_c(option = "viridis", trans = "log10", name = "Footprints
(log scale)",
                                  breaks = c(1, 10, 100, 1000), labels = c("1", "10", "100", "1,000"), limits = c(1, NA)),
  spectral = scale_fill_distiller(palette = "Spectral", direction = -1, trans = "log10", name = "Footprints
(log scale)",
                                  breaks = c(1, 10, 100, 1000), labels = c("1", "10", "100", "1,000"), limits = c(1, NA)))
SUF <- if (PAL == "spectral") "" else paste0("_", PAL)
th <- theme_minimal(base_size = 9, base_family = "sans") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "#e6e5e1", linewidth = 0.25),
        plot.title = element_text(face = "bold", size = 9.5), plot.subtitle = element_text(size = 8.5, colour = "#1f1f1e"), axis.text = element_text(colour = INK2),
        legend.position = "right", legend.title = element_text(size = 8), legend.text = element_text(size = 7.5))
one <- function(y, ttl, i) {
  ggplot(d[chm_mean >= LIM[1] & chm_mean <= LIM[2] & get(y) >= LIM[1] & get(y) <= LIM[2]], aes(chm_mean, .data[[y]])) +
    geom_hex(binwidth = c(0.75, 0.75)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "black", linewidth = 0.4) +
    FILL +
    coord_equal(xlim = LIM, ylim = LIM, expand = FALSE) +
    scale_x_continuous(breaks = seq(0, 40, 10)) + scale_y_continuous(breaks = seq(0, 40, 10), labels = function(v) minus(v)) +
    labs(title = ttl, subtitle = lab(i), x = "Lidar canopy height, footprint mean (m)", y = "Canopy height from P3D (m)") + th
}
g <- one("p3d_chm", "(a) P3D CHM", 1) + one("p3d_dsm_lg", "(b) P3D DSM minus the lidar ground surface", 2) +
  plot_layout(guides = "collect")
ggsave(file.path(OUT, paste0("fig2_height_scatter", SUF, ".png")), g, width = 180, height = 90, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(OUT, paste0("fig2_height_scatter", SUF, ".pdf")), g, width = 180, height = 90, units = "mm", device = cairo_pdf)
