# fig3_error_distributions.R  --  Fig. 3, error distributions and error by canopy height class
# (a) error distributions with each site weighted equally (weight 1/n_site): share of footprints per 0.5 m bin
#     on a log scale, with a normal curve of the same weighted mean and SD (dashed) so the heavy tails show.
#     One small panel per product, order DSM, DTM, CHM.
# (b) mean error by lidar canopy height class (90th percentile), each site weighted equally (site means
#     averaged over sites with >= 30 footprints in the class), with a 10th to 90th percentile bar per class
#     (each site weighted equally).
# Inputs  (P3D_DATA_ROOT) rework/frame_footprints.csv.gz (630,390 footprints)
#         (P3D_RESULTS)   table2/error_tails.csv (site_weighted_accuracy.py, for the check)
# Outputs (P3D_RESULTS)   fig3/fig3_error_distributions.{png,pdf}
#                         fig3/fig3b_class_means.csv (the plotted numbers of (b); Results 3.1)
# Checks  the beyond-3-SD shares of (a) equal error_tails.csv.
# Colours: DSM #2a78d6, DTM #eb6834, CHM #1baf7a. The CHM colour has less than 3:1 contrast,
# so every series is labelled directly and (b) adds point shapes.
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })

DATA_ROOT <- Sys.getenv("P3D_DATA_ROOT", "<LOCAL_DATA_ROOT>")
OUT_ROOT  <- Sys.getenv("P3D_RESULTS", "results")
fr <- fread(file.path(DATA_ROOT, "rework", "frame_footprints.csv.gz"),
            colClasses = list(character = "shot_number"))
stopifnot(nrow(fr) == 630390L)
OUTF <- file.path(OUT_ROOT, "fig3"); OUTT <- OUTF
dir.create(OUTT, showWarnings = FALSE, recursive = TRUE)

PROD <- c("DSM", "DTM", "CHM")
COL  <- c(DSM = "#2a78d6", DTM = "#eb6834", CHM = "#1baf7a")
INK2 <- "#52514e"; GRID <- "#e6e5e1"
minus <- function(x) gsub("-", "−", format(x, trim = TRUE))
L <- melt(fr, id.vars = c("site", "chm_p90"),
          measure.vars = c("dsm_error_mean", "dtm_error_mean", "chm_error_mean"),
          variable.name = "product", value.name = "err")
L[, product := factor(toupper(sub("_error_mean", "", product)), levels = PROD)]

th <- theme_minimal(base_size = 9, base_family = "sans") +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = GRID, linewidth = 0.25),
        axis.text = element_text(colour = INK2), axis.title = element_text(colour = "#1f1f1e"),
        plot.title = element_text(face = "bold", size = 9), legend.position = "none",
        strip.text = element_text(face = "bold", size = 9, hjust = 0))

# ---- (a) distributions ----
BW <- 0.5; XL <- c(-30, 30)
edges <- seq(XL[1], XL[2], by = BW)
L[, w := 1 / .N, by = .(product, site)]
H <- L[, {
  k <- findInterval(err, edges, rightmost.closed = TRUE)
  ok <- k >= 1 & k < length(edges)
  sw <- tapply(w[ok], factor(k[ok], levels = seq_len(length(edges) - 1)), sum); sw[is.na(sw)] <- 0
  .(x = head(edges, -1) + BW / 2, share = 100 * as.numeric(sw) / sum(w))
}, by = product]
N <- L[, {m <- weighted.mean(err, w); .(m = m, s = sqrt(weighted.mean((err - m)^2, w)))}, by = product]
Nc <- N[, .(x = seq(XL[1], XL[2], by = 0.1)), by = .(product, m, s)]
Nc[, share := 100 * BW * dnorm(x, m, s)]
FLOOR <- 1e-3                                   # 0.001% of footprints
H <- H[share > 0]; Nc <- Nc[share >= FLOOR]
tails <- L[, {m <- weighted.mean(err, w); s <- sqrt(weighted.mean((err - m)^2, w)); p <- weighted.mean(abs(err - m) > 3 * s, w)
              .(p = 100 * p, lab = sprintf("%.1f%% beyond 3 SD
(%.0f times the normal share)", 100 * p, p / (2 * pnorm(-3))))},
           by = product]
chk <- fread(file.path(OUT_ROOT, "table2", "error_tails.csv"))[frame == "Sites weighted equally"]
stopifnot(all(abs(tails[match(chk$product, as.character(product)), p] - chk$beyond3sd_pct) < 1e-6))
off <- L[abs(err) > 30, .N, by = product]; print(off)
pa <- ggplot() +
  annotate("segment", x = 0, xend = 0, y = FLOOR, yend = 20, colour = "#b9b8b3", linewidth = 0.3) +
  geom_line(data = Nc, aes(x, share), colour = INK2, linetype = "dashed", linewidth = 0.35) +
  geom_step(data = H, aes(x - BW / 2, share, colour = product), linewidth = 0.45, direction = "hv") +
  geom_text(data = tails, aes(x = XL[1] + 1, y = 250, label = lab), hjust = 0, vjust = 1,
            size = 2.5, colour = INK2, lineheight = 0.95) +
  facet_wrap(~product, nrow = 1) +
  scale_colour_manual(values = COL) +
  scale_y_log10(limits = c(FLOOR, 300), breaks = c(0.001, 0.01, 0.1, 1, 10),
                labels = c("0.001", "0.01", "0.1", "1", "10")) +
  scale_x_continuous(limits = XL, breaks = seq(-20, 20, 10), labels = minus, expand = c(0, 0)) +
  labs(title = "(a) Error distributions, each site weighted equally",
       x = "Error, P3D minus lidar (m)", y = "Share of footprints (%)") +
  th + theme(panel.spacing.x = unit(5, "mm"))

# ---- (b) mean error by canopy height class, sites weighted equally ----
br <- c(2, 5, 10, 15, 20, 30, Inf); lab <- c("2\u20135", "5\u201310", "10\u201315", "15\u201320", "20\u201330", ">30")
L[, hc := cut(chm_p90, br, lab, right = FALSE)]
S <- L[, .(m = mean(err), n = .N), by = .(product, hc, site)][n >= 30]
E <- S[, .(mean_error = mean(m), n_sites = .N), by = .(product, hc)][order(product, hc)]
wq <- function(x, w, p) { o <- order(x); x <- x[o]; cw <- cumsum(w[o]) / sum(w); x[which(cw >= p)[1]] }
keep <- S[, .(product, hc, site, n)]
LQ <- merge(L[!is.na(hc)], keep, by = c("product", "hc", "site"))
Q <- LQ[, .(q10 = wq(err, 1 / n, 0.10), q90 = wq(err, 1 / n, 0.90)), by = .(product, hc)]
E <- merge(E, Q, by = c("product", "hc"))[order(product, hc)]

E[, xi := as.integer(hc) + c(DSM = -0.18, DTM = 0, CHM = 0.18)[as.character(product)]]
lab_b <- E[hc == ">30"]
pb <- ggplot(E, aes(xi, mean_error, colour = product, shape = product)) +
  geom_hline(yintercept = 0, colour = "#b9b8b3", linewidth = 0.3) +
  geom_linerange(aes(ymin = q10, ymax = q90), linewidth = 0.5, alpha = 0.55, show.legend = FALSE) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2.2, fill = "white", stroke = 0.9) +
  geom_text(data = lab_b, aes(label = product), nudge_x = 0.5, hjust = 0, size = 3,
            fontface = "bold", show.legend = FALSE) +
  scale_colour_manual(values = COL) +
  scale_shape_manual(values = c(DSM = 21, DTM = 24, CHM = 22)) +
  scale_x_continuous(breaks = seq_along(lab), labels = lab, limits = c(0.6, 7.1)) +
  scale_y_continuous(breaks = seq(-15, 10, 5), labels = minus) +
  labs(title = "(b) Error by canopy height, each site weighted equally",
       x = "Lidar canopy height, 90th percentile (m)", y = "Error, P3D minus lidar (m)") +
  th

g <- pa / pb + plot_layout(heights = c(1, 1.15))
ggsave(file.path(OUTF, "fig3_error_distributions.png"), g, width = 180, height = 150, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(OUTF, "fig3_error_distributions.pdf"), g, width = 180, height = 150, units = "mm", device = cairo_pdf)
fwrite(E, file.path(OUTT, "fig3b_class_means.csv")); print(tails); print(E)
cat("wrote fig3_error_distributions.{png,pdf}\n")
