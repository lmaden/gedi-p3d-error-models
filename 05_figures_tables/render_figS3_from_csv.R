# render_figS3_from_csv.R (14 Sep 2026): rebuild Supplementary Fig. S3 (predictor correlation matrix + VIF)
# from the two CSVs exported on the cluster by s3_vif_corr_probdir_16site.R. Run locally with R 4.4.0.
# Usage: Rscript render_figS3_from_csv.R <s3_corr_matrix_chm16.csv> <s3_vif_chm16.csv> <out_png>
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
args <- commandArgs(trailingOnly = TRUE)
corr_csv <- args[1]; vif_csv <- args[2]; out_png <- args[3]
lab <- c(slope_mean_z = "Slope (mean)", slope_sd_z = "Slope variability (SD)", aspect_sin_z = "Aspect (sin)",
         aspect_cos_z = "Aspect (cos)", wsci_z = "WSCI", rh_98_z = "RH98", cover_z = "Canopy cover",
         meta_offnad_z = "Off-nadir angle", meta_sunel_z = "Sun elevation", meta_az_conc_z = "Azimuth concentration",
         view_az_sin_z = "View azimuth (sin)", view_az_cos_z = "View azimuth (cos)", meta_stereo_z = "Stereo ratio",
         meta_fwdrev_z = "Forward-reverse ratio", meta_absgeo_z = "Absolute geoaccuracy",
         meta_relgeo_z = "Relative geoaccuracy", meta_leafon_z = "Leaf-on fraction")
R <- as.matrix(read.csv(corr_csv, row.names = 1, check.names = FALSE))
stopifnot(all(rownames(R) %in% names(lab)))
ord <- names(lab)[names(lab) %in% rownames(R)]
R <- R[ord, ord]
long <- data.frame(x = factor(rep(lab[ord], each = length(ord)), levels = lab[ord]),
                   y = factor(rep(lab[ord], times = length(ord)), levels = rev(lab[ord])),
                   r = as.vector(t(R)),
                   i = rep(seq_along(ord), each = length(ord)),
                   j = rep(seq_along(ord), times = length(ord)))
long <- long[long$j >= long$i, ]          # lower triangle only; the matrix is symmetric
# labels drop the leading zero so a two-decimal r fits inside one tile
fmt_r <- function(r) sub("^(-?)0\\.", "\\1.", sprintf("%.2f", r))
p1 <- ggplot(long, aes(x, y, fill = r)) + geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(abs(r) >= 0.3 & r < 1, fmt_r(r), "")), size = 1.8) +
  scale_fill_gradient2(low = "#b2182b", mid = "white", high = "#2166ac", limits = c(-1, 1), name = "Pearson r") +
  coord_equal() + labs(x = NULL, y = NULL, title = "(a) Predictor correlation matrix") +
  theme_minimal(base_size = 9) + theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7),
                                       axis.text.y = element_text(size = 7), panel.grid = element_blank(),
                                       legend.position = c(0.82, 0.78), legend.direction = "horizontal",
                                       legend.title = element_text(size = 7), legend.text = element_text(size = 6),
                                       legend.key.height = grid::unit(0.30, "cm"),
                                       legend.key.width = grid::unit(0.55, "cm"),
                                       plot.title = element_text(face = "bold", size = 9))
v <- read.csv(vif_csv); v$label <- factor(lab[v$predictor], levels = rev(lab[ord]))
p2 <- ggplot(v, aes(x = vif, y = label)) + geom_col(fill = "#4d4d4d", width = 0.7) +
  geom_vline(xintercept = 10, linetype = "dashed", colour = "#b2182b", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.1f", vif)), hjust = -0.15, size = 2.4) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = NULL, y = NULL, title = "(b) VIF") +
  theme_minimal(base_size = 9) + theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 7),
                                       plot.title = element_text(face = "bold", size = 9),
                                       plot.margin = grid::unit(c(4, 12, 4, 4), "pt"))
p <- p1 + p2 + plot_layout(widths = c(2.1, 1))
ggsave(out_png, p, width = 6.5, height = 4.2, units = "in", dpi = 400, bg = "white")
cat(sprintf("wrote %s | max VIF %.2f (%s)\n", out_png, max(v$vif), as.character(v$label[which.max(v$vif)])))
