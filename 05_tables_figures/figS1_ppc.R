# figS1_ppc.R  --  Fig. S1, posterior predictive check of the three error models
# Share of the 68,786 training footprints per 0.5 m error bin (log scale, Fig. 3a style), observed (product colour)
# against 200 replicate error sets drawn from each fitted model (grey band = 5th to 95th percentile over replicates,
# dashed = replicate median). In-panel text: share of footprints with an error beyond 10 m, observed and replicated
# (median and 5th to 95th percentile), two decimals; figS1_numbers.csv also holds the 5 m shares and the 1st and 99th
# percentiles.
# Inputs  (P3D_DERIVED) ppc/{bins,stats}_{dsm,dtm,chm}.csv (rework_ppc.R on the cluster)
# Outputs (P3D_RESULTS) figS1/figS1_ppc.{png,pdf}, figS1/figS1_numbers.csv
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
DERIVED  <- Sys.getenv("P3D_DERIVED", "data_derived")
OUT_ROOT <- Sys.getenv("P3D_RESULTS", "results")
OUT <- file.path(OUT_ROOT, "figS1"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
IN <- file.path(DERIVED, "ppc")
B <- rbindlist(lapply(c("dsm", "dtm", "chm"), function(r) fread(file.path(IN, sprintf("bins_%s.csv", r)))))
S <- rbindlist(lapply(c("dsm", "dtm", "chm"), function(r) fread(file.path(IN, sprintf("stats_%s.csv", r)))), fill = TRUE)
B[, product := factor(toupper(response), levels = c("DSM", "DTM", "CHM"))]
S[, product := factor(toupper(response), levels = c("DSM", "DTM", "CHM"))]
stopifnot(S[set == "replicate", .N, by = product]$N == 200L, B[, .N, by = product]$N == 120L)
stopifnot(B[, abs(sum(obs) - 100) < 0.5, by = product]$V1)          # bins hold (nearly) all footprints
NUM <- S[, .(obs5 = pct_over5[set == "observed"], rep5 = median(pct_over5[set == "replicate"]),
             rep5_lo = quantile(pct_over5[set == "replicate"], 0.05), rep5_hi = quantile(pct_over5[set == "replicate"], 0.95),
             obs10 = pct_over10[set == "observed"], rep10 = median(pct_over10[set == "replicate"]),
             rep10_lo = quantile(pct_over10[set == "replicate"], 0.05), rep10_hi = quantile(pct_over10[set == "replicate"], 0.95),
             obs_q01 = q01[set == "observed"], rep_q01 = median(q01[set == "replicate"]),
             obs_q99 = q99[set == "observed"], rep_q99 = median(q99[set == "replicate"])), by = product]
fwrite(NUM, file.path(OUT, "figS1_numbers.csv")); print(NUM)
COL <- c(DSM = "#2a78d6", DTM = "#eb6834", CHM = "#1baf7a"); INK2 <- "#52514e"
FLOOR <- 1e-3; BW <- 0.5
B2 <- copy(B)[, `:=`(lo = pmax(rep_q05, FLOOR), hi = pmax(rep_q95, FLOOR), med = pmax(rep_median, FLOOR))]
lab <- NUM[, .(product, l = sprintf("Error beyond ±10 m\nobserved %.2f%%\nreplicates %.2f%% (%.2f to %.2f)", obs10, rep10, rep10_lo, rep10_hi))]
minus <- function(x) gsub("-", "\u2212", format(x, trim = TRUE))
g <- ggplot(B2) +
  geom_rect(aes(xmin = x - BW / 2, xmax = x + BW / 2, ymin = lo, ymax = hi), fill = "#cfcec9") +
  geom_step(aes(x - BW / 2, med), linetype = "dashed", colour = INK2, linewidth = 0.35, direction = "hv") +
  geom_step(data = B2[obs > 0], aes(x - BW / 2, obs, colour = product), linewidth = 0.45, direction = "hv") +
  geom_text(data = lab, aes(x = -29, y = 120, label = l), hjust = 0, vjust = 1, size = 2.6, colour = INK2, lineheight = 0.95) +
  facet_wrap(~product, nrow = 1) + scale_colour_manual(values = COL, guide = "none") +
  scale_y_log10(limits = c(FLOOR, 150), breaks = c(0.001, 0.01, 0.1, 1, 10), labels = c("0.001", "0.01", "0.1", "1", "10")) +
  scale_x_continuous(limits = c(-30, 30), breaks = seq(-20, 20, 10), labels = minus, expand = c(0, 0)) +
  labs(x = "Error, P3D minus lidar (m)", y = "Share of footprints (%)") +
  theme_minimal(base_size = 9, base_family = "sans") +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "#e6e5e1", linewidth = 0.25), axis.text = element_text(colour = INK2),
        strip.text = element_text(face = "bold", size = 9, hjust = 0), panel.spacing.x = unit(5, "mm"))
ggsave(file.path(OUT, "figS1_ppc.png"), g, width = 180, height = 75, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(OUT, "figS1_ppc.pdf"), g, width = 180, height = 75, units = "mm", device = cairo_pdf)
