# figS2_variogram.R  --  Fig. S2, spatial correlation of the P3D DSM, DTM and CHM errors within sites
# Correlation at distance h = 1 - gamma(h) / within-site variance, from the classical semivariogram of each site's
# errors minus the site mean (within_site_variogram.py). Line and points: median over sites; band: middle
# half of sites (25th to 75th percentile). Distance classes with at least 100 footprint pairs at a site count for it.
# Input  (P3D_RESULTS) figS2/variogram_summary.csv (within_site_variogram.py)
# Output (P3D_RESULTS) figS2/figS2_variogram.{png,pdf}
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
OUT_ROOT <- Sys.getenv("P3D_RESULTS", "results")
OUT <- file.path(OUT_ROOT, "figS2"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
S <- fread(file.path(OUT, "variogram_summary.csv"))
S[, product := factor(toupper(product), levels = c("DSM", "DTM", "CHM"))]
S[, `:=`(r = 1 - g_med, r_lo = 1 - g_q3, r_hi = 1 - g_q1)]
stopifnot(S[, .N, by = product]$N == 12L, S[bin_lo == 30, round(r[product == "CHM"], 2)] == 0.57)
lab <- S[bin_lo == 30, .(product, h_med, r, l = sprintf("%.2f at %.0f m", r, h_med))]
COL <- c(DSM = "#2a78d6", DTM = "#eb6834", CHM = "#1baf7a"); INK <- "#1f1f1e"; INK2 <- "#52514e"
num <- function(x) gsub("-", "\u2212", formatC(x, format = "fg", big.mark = ","))
g <- ggplot(S, aes(h_med)) +
  geom_hline(yintercept = 0, colour = INK2, linewidth = 0.3) +
  geom_ribbon(aes(ymin = r_lo, ymax = r_hi, fill = product), alpha = 0.2) +
  geom_line(aes(y = r, colour = product), linewidth = 0.6) +
  geom_point(aes(y = r, colour = product), size = 1.4, shape = 21, fill = "white", stroke = 0.6) +
  geom_text(data = lab, aes(x = h_med * 1.18, y = r + 0.06, label = l), hjust = 0, vjust = 0, size = 2.6, colour = INK) +
  facet_wrap(~product, nrow = 1) +
  scale_colour_manual(values = COL, guide = "none") + scale_fill_manual(values = COL, guide = "none") +
  scale_x_log10(limits = c(50, 3000), breaks = c(60, 100, 300, 1000, 3000), labels = num) +
  scale_y_continuous(limits = c(-0.05, 1), breaks = seq(0, 1, 0.25), labels = function(x) num(x)) +
  labs(x = "Distance between footprints (m)", y = "Correlation of errors") +
  theme_minimal(base_size = 9, base_family = "sans") +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "#e6e5e1", linewidth = 0.25), axis.text = element_text(colour = INK2),
        strip.text = element_text(face = "bold", size = 9, hjust = 0), panel.spacing.x = unit(5, "mm"))
ggsave(file.path(OUT, "figS2_variogram.png"), g, width = 180, height = 70, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(OUT, "figS2_variogram.pdf"), g, width = 180, height = 70, units = "mm", device = cairo_pdf)
print(S[, .(product, bin = sprintf("%g-%g", bin_lo, bin_hi), n_sites, h_med = round(h_med), r = round(r, 2),
            r_lo = round(r_lo, 2), r_hi = round(r_hi, 2), same_track = round(same_track_med, 2))])
