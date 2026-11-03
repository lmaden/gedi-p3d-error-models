# fig4_site_accuracy.R  --  Fig. 4, accuracy at each site
# Rows are sites sorted by median lidar canopy height, and the site labels carry the footprint counts.
# (a) median lidar canopy height (the sort key, shown first)
# (b) bias of the two surfaces, DSM and DTM (the CHM bias is close to their difference, Eq. 4)
# (c) CHM RMSE  (d) share of footprints with a CHM error over 5 m (CHM green in (c) and (d))
# Per-site DSM and DTM RMSE go to Table S4.
# Inputs  (P3D_RESULTS) persite/persite.csv (persite_class_tables.R)
#         (P3D_DERIVED) site_names.csv
# Outputs (P3D_RESULTS) fig4/fig4_site_accuracy.{png,pdf}
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })

DERIVED  <- Sys.getenv("P3D_DERIVED", "data_derived")
OUT_ROOT <- Sys.getenv("P3D_RESULTS", "results")
OUT  <- file.path(OUT_ROOT, "fig4"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
ps <- fread(file.path(OUT_ROOT, "persite", "persite.csv"))
nm <- fread(file.path(DERIVED, "site_names.csv"))
ps <- merge(ps, nm[, .(site = tracker, paper_site, short_name)], by = "site")
stopifnot(nrow(ps) == 18L)
setorder(ps, median_height)
ps[, yi := .I]                                     # tallest canopy on the top row
ps[, lab := sprintf("%d %s (%s)", paper_site, short_name, format(n, big.mark = ",", trim = TRUE))]

COL  <- c(DSM = "#2a78d6", DTM = "#eb6834", CHM = "#1baf7a")
SHP  <- c(DSM = 21, DTM = 24)
INK2 <- "#52514e"; GRID <- "#e6e5e1"; BAR <- "#9c9b96"
minus <- function(x) gsub("-", "\u2212", format(x, trim = TRUE))
ys <- scale_y_continuous(breaks = ps$yi, labels = ps$lab, limits = c(0.5, 18.5), expand = c(0, 0))
stripes <- annotate("rect", xmin = -Inf, xmax = Inf, ymin = seq(2, 18, 2) - 0.5, ymax = seq(2, 18, 2) + 0.5,
                    fill = "#f3f2ef")
th <- theme_minimal(base_size = 9, base_family = "sans") +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = GRID, linewidth = 0.25),
        axis.text = element_text(colour = INK2), axis.text.y = element_text(colour = "#1f1f1e", size = 8),
        axis.title.x = element_text(size = 8.5, colour = "#1f1f1e", margin = margin(t = 3)),
        plot.title = element_text(face = "bold", size = 9, margin = margin(b = 14)),
        axis.title.y = element_blank(), legend.position = "none")
noy <- theme(axis.text.y = element_blank())

B <- rbind(ps[, .(yi, value = dsm_bias, product = "DSM")], ps[, .(yi, value = dtm_bias, product = "DTM")])
B[, y := yi + ifelse(product == "DSM", 0.15, -0.15)]

pa <- ggplot(ps) + stripes +
  geom_rect(aes(xmin = 0, xmax = median_height, ymin = yi - 0.3, ymax = yi + 0.3), fill = BAR) + ys +
  scale_x_continuous(limits = c(0, 28), breaks = seq(0, 25, 5), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "(a) Canopy height", x = "Median lidar canopy\nheight (m)") + th
pb <- ggplot(B, aes(value, y, colour = product, shape = product)) + stripes +
  geom_vline(xintercept = 0, colour = "#8f8e89", linewidth = 0.35) +
  geom_point(size = 2, fill = "white", stroke = 0.9) + ys +
  scale_colour_manual(values = COL, name = NULL, breaks = c("DSM", "DTM")) +
  scale_shape_manual(values = SHP, name = NULL, breaks = c("DSM", "DTM")) +
  scale_x_continuous(limits = c(-6, 7), breaks = seq(-6, 6, 2), labels = minus) +
  labs(title = "(b) DSM and DTM bias", x = "Mean error, P3D minus\nlidar (m)") + th + noy +
  theme(legend.position = "inside", legend.position.inside = c(0, 1.005), legend.justification = c(0, 0),
        legend.direction = "horizontal", legend.text = element_text(size = 8.5),
        legend.margin = margin(0, 0, 0, 0), legend.background = element_blank(),
        legend.key.size = unit(3.5, "mm")) +
  coord_cartesian(clip = "off")
pc <- ggplot(ps) + stripes +
  geom_rect(aes(xmin = 0, xmax = chm_rmse, ymin = yi - 0.3, ymax = yi + 0.3), fill = COL[["CHM"]]) + ys +
  scale_x_continuous(limits = c(0, 10.5), breaks = seq(0, 10, 2), expand = expansion(mult = c(0, 0.03))) +
  labs(title = "(c) CHM RMSE", x = "CHM RMSE (m)\n") + th + noy
pd <- ggplot(ps) + stripes +
  geom_rect(aes(xmin = 0, xmax = chm_pct_over5, ymin = yi - 0.3, ymax = yi + 0.3), fill = COL[["CHM"]]) + ys +
  scale_x_continuous(limits = c(0, 70), breaks = seq(0, 60, 20), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "(d) Large errors", x = "Footprints with\n|CHM error| > 5 m (%)") + th + noy

g <- pa + pb + pc + pd + plot_layout(nrow = 1, widths = c(0.95, 1.3, 0.9, 0.8))
ggsave(file.path(OUT, "fig4_site_accuracy.png"), g, width = 180, height = 125, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(OUT, "fig4_site_accuracy.pdf"), g, width = 180, height = 125, units = "mm", device = cairo_pdf)
cat("wrote fig4_site_accuracy.{png,pdf}\n")
