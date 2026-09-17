# Supplementary Fig. S11: CHM site-level random intercepts and the random-effect
# variance partition, for the 16-site CHM fit.
#
# Shares the plotting path of supp_figs_plotting.R. Values are the published ones,
# Table S6b RE_CHM medians with 95% credible intervals and the Table S9 CHM variance
# shares (56.9 / 40.5 / 2.1), so the in-image legend matches the tables.
# Output directory is taken from FIG_OUT_DIR, defaulting to the working directory.
suppressPackageStartupMessages({ library(ggplot2); library(data.table); library(patchwork); library(scales) })
COL_CHM <- "#2D5F3F"; COL_POOL <- "#555555"; COL_REF <- "#A02020"; COL_EREG <- "#D8B044"; COL_LC <- "#88B0D8"
theme_supp <- function(base_size = 9) {
  theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
          strip.background = element_rect(fill = "grey95", color = NA),
          strip.text = element_text(size = base_size - 1, face = "bold"),
          plot.title = element_text(size = base_size + 1, face = "bold"),
          legend.title = element_text(size = base_size - 1),
          legend.text = element_text(size = base_size - 1))
}
re_df <- fread("
ms_site,re_mean,re_q025,re_q975
4,-1.46,-3.56,0.77
5,1.56,0.46,2.69
6,2.09,1.00,3.24
7,-2.19,-3.28,-1.03
8,-1.12,-2.21,0.03
9,-1.70,-2.78,-0.55
10,0.99,-0.32,2.32
11,-0.24,-1.54,1.09
12,0.27,-1.15,1.67
13,-0.58,-1.98,0.81
14,-1.36,-2.83,0.05
15,0.85,-0.82,2.39
16,0.94,-0.38,2.30
17,0.62,-0.72,1.95
18,-0.01,-1.77,1.78
19,1.23,-0.60,2.94")
setorder(re_df, ms_site)
re_df[, site_lab := factor(as.character(ms_site), levels = as.character(ms_site))]
p_top <- ggplot(re_df, aes(x = site_lab, y = re_mean)) +
  geom_errorbar(aes(ymin = re_q025, ymax = re_q975), width = 0, color = COL_POOL, linewidth = 0.4) +
  geom_point(color = COL_CHM, size = 2) +
  geom_hline(yintercept = 0, color = COL_REF, linetype = "dashed", linewidth = 0.4) +
  labs(x = "Site (manuscript numbering)", y = "Site random intercept (m)") +
  theme_supp() + theme(axis.text.x = element_text(size = 8))
vp <- data.table(grouping = c("site", "ecoregion", "lc_l1_code"), mean = c(0.569, 0.405, 0.021))
vp <- vp[order(-mean)]
nice <- c("site" = "Site", "ecoregion" = "Ecoregion", "lc_l1_code" = "Land cover")
vp[, nice_label := nice[grouping]]
vp[, legend_label := sprintf("%s (%.1f%%)", nice_label, mean * 100)]
vp[, legend_label := factor(legend_label, levels = legend_label)]
vp[, x := 1]
color_map <- setNames(c(COL_CHM, COL_EREG, COL_LC), vp$legend_label)
p_bot <- ggplot(vp, aes(x = x, y = mean, fill = legend_label)) +
  geom_col(width = 0.4, color = "white", linewidth = 0.5) +
  scale_fill_manual(values = color_map, name = NULL) +
  coord_flip() +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  labs(x = NULL, y = "Share of total variance") +
  theme_supp() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major.y = element_blank(),
        legend.position = "bottom", legend.direction = "horizontal", legend.text = element_text(size = 9))
p <- p_top / p_bot + plot_layout(heights = c(3, 1.2))
ggsave(file.path(Sys.getenv("FIG_OUT_DIR", "."), "figS11.png"), p, width = 6.5, height = 6.5, units = "in", dpi = 300, bg = "white")
cat("legend labels:", paste(levels(vp$legend_label), collapse = " | "), "\n")
