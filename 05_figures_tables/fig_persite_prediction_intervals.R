# =====================================================================
# section_J_fig_persite_pi_ggplot.R   (self-contained: embedded verified data)
# J4a: CHM per-site prediction-interval width + calibration (16-site fit)
#   Panel A: per-site median 95% PI width (lollipop), sites ascending; CHM overall
#            median (17.6 m) and DTM overall median (11.0 m) reference lines.
#   Panel B: per-site empirical 95% PI coverage, points sized by footprint count;
#            nominal 0.95 and CHM overall 0.934 reference lines; under-covering
#            thin site(s) (cov < 0.90) highlighted.
#
# Data: embedded below, emitted verbatim from the VERIFIED
#       prediction_interval_widths_chm16.csv (section_H_pi_width_chm16.R output).
#       Site labels are MANUSCRIPT numbers (tracker->manuscript via s9b).
#       (Set PIW_CSV to a file path to override with a CSV instead of the embed.)
# No model load, no fabrication.
# =====================================================================

Sys.setenv(DISPLAY = ""); options(device = pdf)
suppressPackageStartupMessages({ library(ggplot2) })
have_patchwork <- requireNamespace("patchwork", quietly = TRUE)

CHPT1_ROOT <- Sys.getenv("CHPT1_ROOT", "/gpfs/data1/vclgp/lmaden/chpt1")
OUT_PNG <- file.path(CHPT1_ROOT, "fig_J4a_persite_PI_chm16_v117.png")
OUT_PDF <- file.path(CHPT1_ROOT, "fig_J4a_persite_PI_chm16_v117.pdf")
CHM <- "#2D5F3F"; DTM <- "#C7522A"; FLAG <- "#B8442E"; GREY <- "#6B6B6B"
OVR_MED <- 17.60; OVR_COV <- 0.9338; DTM_MED <- 11.0; NOMINAL <- 0.95

# ---- embedded verified per-site data (ascending median_width) ----
ms           <- c(17, 16, 19, 10, 15, 11, 18, 8, 12, 6, 5, 13, 9, 14, 4, 7)
median_width <- c(3.86525253729, 4.11608374882, 4.21190318562, 4.94592556874,
                  5.96122078148, 8.45132997657, 9.06703244819, 15.3823064255,
                  15.8514036102, 16.2560205798, 18.372393904, 18.4787030034,
                  21.1490533261, 23.5948730326, 26.0742113037, 27.9617377502)
coverage_95  <- c(0.808510638298, 1, 0.952504638219, 0.948670662364, 0.900406504065,
                  0.9405360134, 0.948618618619, 0.959243085881, 0.956998683633,
                  0.949124467582, 0.918761289961, 0.957251610385, 0.940290620872,
                  0.905683192261, 0.97234920281, 0.977216443784)
n            <- c(47, 54, 2695, 6507, 984, 8358, 33300, 4809, 6837, 12678,
                  142826, 5123, 7570, 827, 8969, 12114)

# optional CSV override
if (nzchar(Sys.getenv("PIW_CSV")) && file.exists(Sys.getenv("PIW_CSV"))) {
  suppressPackageStartupMessages(library(data.table))
  d <- fread(Sys.getenv("PIW_CSV")); ps <- d[product=="CHM" & level=="site"]
  trk <- as.integer(sub("site=","",ps$stratum))
  map <- c(`4`=4,`5`=5,`6`=6,`7`=7,`8`=8,`9`=9,`10`=10,`11`=11,`12`=12,`13`=13,
           `14`=14,`15`=15,`17`=16,`18`=17,`19`=18,`20`=19)
  o <- order(ps$median_width)
  ms <- unname(map[as.character(trk[o])]); median_width <- ps$median_width[o]
  coverage_95 <- ps$coverage_95[o]; n <- ps$n[o]
  cat("using PIW_CSV override\n")
}

df <- data.frame(ms = ms, median_width = median_width, coverage_95 = coverage_95, n = n)
df$lab <- factor(df$ms, levels = df$ms)              # x order = width order
df$under <- df$coverage_95 < 0.90
cat(sprintf("sites=%d  width %.2f-%.2f (%.1fx)\n", nrow(df),
            min(df$median_width), max(df$median_width),
            max(df$median_width)/min(df$median_width)))

pA <- ggplot(df, aes(lab, median_width)) +
  geom_hline(yintercept = OVR_MED, linetype = "dashed", color = GREY, linewidth = .4) +
  geom_hline(yintercept = DTM_MED, linetype = "dotted", color = DTM,  linewidth = .5) +
  geom_segment(aes(xend = lab, y = 0, yend = median_width), color = CHM, linewidth = .9) +
  geom_point(color = CHM, size = 2.4) +
  annotate("text", x = 1, y = OVR_MED + 0.7,
           label = sprintf("CHM overall median %.1f m", OVR_MED), hjust = 0, vjust = 0, size = 3, color = GREY) +
  annotate("text", x = 1, y = DTM_MED - 0.7,
           label = sprintf("DTM overall median %.1f m", DTM_MED), hjust = 0, vjust = 1, size = 3, color = DTM) +
  annotate("text", x = 1, y = max(df$median_width)*1.02,
           label = sprintf("%.1f\u00d7 range across sites", max(df$median_width)/min(df$median_width)),
           hjust = 0, vjust = 1, fontface = "italic", size = 3.1, color = "grey20") +
  scale_y_continuous(limits = c(0, max(df$median_width)*1.08), expand = c(0,0)) +
  labs(y = "Per-site median 95%\nPI width (m)") +
  theme_classic(base_size = 11) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), plot.margin = margin(6,8,2,6))

pB <- ggplot(df, aes(lab, coverage_95)) +
  geom_hline(yintercept = NOMINAL, linetype = "dashed", color = GREY, linewidth = .4) +
  geom_hline(yintercept = OVR_COV, color = CHM, alpha = .55, linewidth = .4) +
  geom_point(aes(size = n, color = under)) +
  scale_color_manual(values = c(`FALSE` = CHM, `TRUE` = FLAG), guide = "none") +
  scale_size_continuous(trans = "log10", range = c(1.6, 6), breaks = c(1e3,1e4,1e5),
                        labels = expression(10^3, 10^4, 10^5), name = "footprints") +
  geom_text(data = df[df$under, ], aes(label = sprintf("Site %d (n=%d)", ms, n)),
            color = FLAG, size = 2.7, vjust = 1.9, hjust = 0.2) +
  annotate("text", x = nrow(df)-0.2, y = NOMINAL + 0.004, label = "nominal 0.95",
           hjust = 1, vjust = 0, size = 3, color = GREY) +
  annotate("text", x = 0.7, y = OVR_COV - 0.006, label = "CHM overall 0.934",
           hjust = 0, vjust = 1, size = 3, color = CHM) +
  scale_y_continuous(limits = c(0.78, 1.02)) +
  labs(x = "Site (ordered by CHM median PI width)", y = "Per-site empirical\n95% PI coverage") +
  theme_classic(base_size = 11) +
  theme(legend.position = c(0.92, 0.16), legend.background = element_blank(),
        legend.key.size = unit(.4,"lines"), legend.title = element_text(size = 8),
        legend.text = element_text(size = 7.5), plot.margin = margin(2,8,6,6))

ttl <- "CHM prediction-interval width and calibration across sites (16-site fit)"
if (have_patchwork) {
  g <- patchwork::wrap_plots(pA, pB, ncol = 1, heights = c(1.12, 1)) +
       patchwork::plot_annotation(title = ttl,
         theme = theme(plot.title = element_text(size = 11, face = "bold")))
  ggsave(OUT_PDF, g, width = 7.4, height = 6.7); ggsave(OUT_PNG, g, width = 7.4, height = 6.7, dpi = 600)
} else {
  suppressPackageStartupMessages(library(gridExtra))
  g <- arrangeGrob(pA, pB, ncol = 1, heights = c(1.12,1), top = ttl)
  ggsave(OUT_PDF, g, width = 7.4, height = 6.7); ggsave(OUT_PNG, g, width = 7.4, height = 6.7, dpi = 600)
}
cat(sprintf("WROTE %s\n      %s\nDone.\n", OUT_PNG, OUT_PDF))
