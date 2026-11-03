# fig5_transects_plot.R  --  Fig. 5, step 3 of 3: draw the three transects
# The same metres per millimetre in every panel (panel heights proportional to each elevation range), so canopy
# heights compare across panels. Rank 1 = the window chosen by fig5_select_transects.R.
# Inputs  (P3D_DATA_ROOT) fig5/fig5_profile_siteNN.csv, fig5/fig5_fpcheck_siteNN.csv (cluster, fig5_profiles.R)
# Outputs (P3D_RESULTS)   fig5/fig5_transects.{png,pdf}
# Panels  (a) Jones Center (tracker 6, paper Site 6), (b) Harvard Forest (4, Site 4),
#         (c) Chequamegon (14, Site 13). Colours: P3D DSM blue, P3D DTM orange,
#         lidar in grays.
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })

DATA_ROOT <- Sys.getenv("P3D_DATA_ROOT", "<LOCAL_DATA_ROOT>")
OUT_ROOT  <- Sys.getenv("P3D_RESULTS", "results")
IN   <- file.path(DATA_ROOT, "fig5")
OUT  <- file.path(OUT_ROOT, "fig5")
RANK <- 1L

panels <- list(
  list(tr = 6,  lab = "(a) Jones Center (Site 6)"),
  list(tr = 4,  lab = "(b) Harvard Forest (Site 4)"),
  list(tr = 14, lab = "(c) Chequamegon (Site 13)"))

COL <- c("P3D DSM" = "#2a78d6", "P3D DTM" = "#eb6834", "Lidar ground surface" = "#1f1f1e")
FILL <- c("Lidar canopy" = "#d4d3cf")
INK2 <- "#52514e"
num <- function(v) sub("-", "\u2212", sprintf("%+.1f m", v), fixed = TRUE)

XMAX <- 1100
SPAN <- {  # common vertical span: the largest panel range plus padding
  rr <- sapply(panels, function(P) { p <- fread(file.path(IN, sprintf("fig5_profile_site%02d.csv", P$tr)))[rank == RANK]
    diff(range(c(p$l_ground, p$p3d_dtm, p$l_ground + p$l_chm, p$p3d_dsm), na.rm = TRUE)) })
  cat("panel ranges (m):", round(rr, 1), "
"); RR <<- rr; max(rr) * 1.12 }
mk <- function(P) {
  p <- fread(file.path(IN, sprintf("fig5_profile_site%02d.csv", P$tr)))[rank == RANK]
  f <- fread(file.path(IN, sprintf("fig5_fpcheck_site%02d.csv", P$tr)),
             colClasses = list(character = "shot_number"))[rank == RANK]
  # distance from the first footprint; the profile starts 25 m before it
  p[, d := dist + 25]; f[, d := dist + 25]
  p[, surf := l_ground + l_chm]
  lines <- rbind(p[, .(d, z = p3d_dsm, what = "P3D DSM")],
                 p[, .(d, z = p3d_dtm, what = "P3D DTM")],
                 p[, .(d, z = l_ground, what = "Lidar ground surface")])
  lines[, what := factor(what, levels = names(COL))]
  ylo <- min(p$l_ground, p$p3d_dtm, na.rm = TRUE)
  yhi <- max(p$surf, p$p3d_dsm, na.rm = TRUE)
  pad <- 0.06 * (yhi - ylo)
  tick_lo <- ylo - 1.6 * pad; tick_hi <- ylo - 0.7 * pad
  sub <- sprintf("Mean error at the %d footprints   DSM %s   DTM %s   CHM %s", nrow(f),
                 num(mean(f$dsm_err)), num(mean(f$dtm_err)), num(mean(f$chm_err)))
  ggplot() +
    geom_ribbon(data = p, aes(x = d, ymin = l_ground, ymax = surf, fill = "Lidar canopy"),
                colour = NA) +
    geom_line(data = lines, aes(x = d, y = z, colour = what), linewidth = 0.35) +
    geom_segment(data = f, aes(x = d, xend = d, y = tick_lo, yend = tick_hi,
                               linetype = "GEDI footprint center"), colour = INK2, linewidth = 0.5,
                 key_glyph = "vpath") +
    scale_colour_manual(values = COL, name = NULL, drop = FALSE) +
    scale_fill_manual(values = FILL, name = NULL, labels = "Lidar surface") +
    scale_linetype_manual(values = c("GEDI footprint center" = "solid"), name = NULL) +
    scale_x_continuous(limits = c(0, XMAX), breaks = seq(0, 1000, 200), expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(0, 2000, 10), expand = c(0, 0)) +
    coord_cartesian(ylim = c(ylo - 2 * pad, yhi + pad)) +
    labs(title = P$lab, subtitle = sub, x = NULL, y = "Elevation (m)") +
    guides(colour = guide_legend(order = 1, override.aes = list(linewidth = 0.9)),
           fill = guide_legend(order = 2), linetype = guide_legend(order = 3)) +
    theme_minimal(base_size = 9, base_family = "sans") +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(colour = "#e6e5e1", linewidth = 0.25),
          axis.line.x = element_line(colour = INK2, linewidth = 0.3),
          axis.ticks.x = element_line(colour = INK2, linewidth = 0.3),
          axis.text = element_text(colour = INK2),
          plot.title = element_text(face = "bold", size = 9.5, margin = margin(b = 1)),
          plot.subtitle = element_text(colour = INK2, size = 8.5, margin = margin(b = 3)),
          legend.position = "none", legend.key.width = unit(9, "mm"),
          legend.key.height = unit(4, "mm"), legend.spacing.y = unit(0, "mm"),
          legend.margin = margin(0, 0, 0, 0),
          plot.margin = margin(4, 6, 4, 4))
}

# the legend sits in the empty right part of panel (c), whose transect is shorter
pl <- lapply(panels, mk)
pl[[3]] <- pl[[3]] + labs(x = "Distance along the GEDI track (m)") +
  theme(legend.position = "inside", legend.position.inside = c(0.99, 0.5),
        legend.justification = c(1, 0.5), legend.box = "vertical",
        legend.box.just = "left", legend.text = element_text(size = 8.5))
g <- wrap_plots(pl, ncol = 1, heights = RR * 1.09)

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(OUT, "fig5_transects.png"), g, width = 180, height = 190, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(OUT, "fig5_transects.pdf"), g, width = 180, height = 190, units = "mm", device = cairo_pdf)
cat("wrote", file.path(OUT, "fig5_transects.{png,pdf}"), "\n")
