# =====================================================================
# suppfig_25_site5_cover_decile.R
# Supp Figure S25 (NEW): Site 5 cover-decile mean-offset bar chart.
#
# Site 5 supp deep-dive.
#
# Layout: single-panel vertical bar chart, D01-D10 on x-axis, mean
# ALS-vs-P3D canopy-height offset (m) on y-axis. Each bar shows a
# 95% CI from offset_se, the offset value as a text label, and the
# per-decile n_tiles below the axis. Headline annotation reports the
# ALL_POOLED b_cover slope from the LC-strata companion file.
#
# Notes:
# The earlier draft of this script was broken in three ways:
#   1. It coerced the character `decile` column ("D01"-"D10") with
#      as.integer(), giving NA for every row and a "DNA" x-axis.
#   2. It did not filter to Site 5 (the CSV contains 30 rows = 3 sites
#      x 10 deciles each; sites are 3, 4, 5 per diagnostic).
#   3. It did not detect the parallel integer column `decile_num`.
# The fixes here:
#   - Filter `manuscript_site == 5L` (confirmed Site 5 has the
#     ALL_POOLED-across-LC pattern that matches §2J anchors).
#   - Use `decile_num` integer column for ordering; label as D01-D10.
#   - Add per-decile 95% CI bars from offset_se (delete the
#     geom_errorbar block to revert to means-only if desired).
#   - Add a stopifnot-style verification block at the end.
# Note: the 30 rows are 3 sites x 10 deciles, not 3 land-cover strata x 10 deciles;
# LC stratification of Site 5 lives in the separate
# site5_lc_strata.csv consumed by S26.
#
# Source:
#   - manuscript_tables/groundwork_task6_phase1_5b_cover_decile_comparison.csv
#
# Story-lock §2J anchors (Site 5 ALL_POOLED):
#   D01 mean offset = +9.40 m (n ~ 1,826)
#   D04 mean offset = -0.78 m
#   D10 mean offset = -3.21 m
#   b_cover         = -17.18 m per unit cover (r_adj^2 = 0.366)
#   Sum n_tiles     = 18,255
# =====================================================================

source("fig_common.R")
fig_banner("Supp Figure S25",
                 "Site 5 cover-decile bar chart")

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(dplyr)
})

# ---------------------------------------------------------------------
# 1. Load CSV and filter to Site 5
# ---------------------------------------------------------------------
path <- file.path(manuscript_tables_dir,
                  "groundwork_task6_phase1_5b_cover_decile_comparison.csv")
if (!file.exists(path)) stop("Required: ", path)
dt <- as.data.table(fread(path))
log_progress(sprintf("  Loaded %d rows, %d cols", nrow(dt), ncol(dt)))
cat("  Columns: ", paste(names(dt), collapse = ", "), "\n")

required_cols <- c("decile_num", "manuscript_site", "n_tiles",
                   "offset_mean", "offset_se", "cover_med")
missing_cols <- setdiff(required_cols, names(dt))
if (length(missing_cols) > 0L) {
  stop("Required columns missing from ", basename(path), ": ",
       paste(missing_cols, collapse = ", "),
       "\nFound: ", paste(names(dt), collapse = ", "))
}

dt5 <- dt[manuscript_site == 5L]
log_progress(sprintf("  Site 5 rows: %d (expected 10)", nrow(dt5)))
if (nrow(dt5) != 10L) {
  stop("Expected 10 Site 5 rows, got ", nrow(dt5),
       ". Check manuscript_site column type/coding.")
}
if (!identical(sort(unique(dt5$decile_num)), 1:10)) {
  stop("decile_num at Site 5 does not span 1:10. Got: ",
       paste(sort(unique(dt5$decile_num)), collapse = ", "))
}

plot_df <- data.table(
  decile        = as.integer(dt5$decile_num),
  mean_offset   = as.numeric(dt5$offset_mean),
  offset_se     = as.numeric(dt5$offset_se),
  n_tiles       = as.integer(dt5$n_tiles),
  cover_median  = as.numeric(dt5$cover_med)
)
setorder(plot_df, decile)
plot_df[, dec_label := sprintf("D%02d", decile)]
plot_df[, dec_label := factor(dec_label, levels = dec_label)]
plot_df[, ci_lo := mean_offset - 1.96 * offset_se]
plot_df[, ci_hi := mean_offset + 1.96 * offset_se]
# Text label sits just past the CI end (above for positive bars,
# below for negative bars). vjust offsets by ~0.5 text heights.
plot_df[, label_y    := ifelse(mean_offset >= 0, ci_hi, ci_lo)]
plot_df[, label_vjust := ifelse(mean_offset >= 0, -0.6, 1.5)]

log_progress("Site 5 cover-decile summary:")
print(plot_df[, .(decile, dec_label,
                  mean_offset = round(mean_offset, 3),
                  offset_se   = round(offset_se,   3),
                  ci_lo       = round(ci_lo,       2),
                  ci_hi       = round(ci_hi,       2),
                  n_tiles,
                  cover_median = round(cover_median, 3))])

# ---------------------------------------------------------------------
# 2. Build plot
# ---------------------------------------------------------------------
bar_fill <- ifelse(plot_df$mean_offset >= 0, "#7baad8", "#e08a8a")

p <- ggplot(plot_df, aes(x = dec_label, y = mean_offset)) +
  geom_hline(yintercept = 0, color = COLOR_ZERO_LINE,
             linewidth = 0.4) +
  geom_col(fill = bar_fill, color = "#444444", linewidth = 0.3) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                width = 0.25, linewidth = 0.4, color = "#222222") +
  geom_text(aes(y = label_y, vjust = label_vjust,
                label = sprintf("%s%.2f",
                                ifelse(mean_offset < 0, MINUS, "+"),
                                abs(mean_offset))),
            size = 3.0) +
  geom_text(aes(y = -Inf,
                label = sprintf("n = %s",
                                formatC(n_tiles, big.mark = ",", format = "d"))),
            vjust = -0.6, size = 2.6, color = "gray30") +
  labs(
    title = "Site 5 ALS-vs-P3D canopy-height offset by canopy-cover decile",
    x = "Canopy-cover decile (D01 = lowest 10% cover; D10 = highest)",
    y = "Mean within-tile offset (m)"
  ) +
  annotate("text",
           x = 1, y = Inf,
           label = sprintf("Site 5, ALL_POOLED tile-level regression:\nb_cover = %s17.18 m per unit cover (r_adj^2 = 0.366)",
                           MINUS),
           hjust = 0, vjust = 1.6, size = 3.0, lineheight = 1.05) +
  scale_y_continuous(expand = expansion(mult = c(0.18, 0.22))) +
  theme_section_G(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    plot.title         = element_text(size = 11, face = "bold")
  )

save_figure(p, "suppfig_S25_site5_cover_decile",
            width_in = 7.5, height_in = 4.8)

# ---------------------------------------------------------------------
# 3. Verification against story-lock §2J anchors
# ---------------------------------------------------------------------
log_subsection("VERIFICATION (console anchors)")
cat("\nStory-lock §2J targets (Site 5, ALL_POOLED):\n")
cat("  D01 mean offset: +9.40 m (n ~ 1,826)\n")
cat("  D04 mean offset: -0.78 m\n")
cat("  D10 mean offset: -3.21 m\n")
cat("  b_cover (ALL_POOLED): -17.18 m per unit cover\n")
cat("  Sum n_tiles across deciles: 18,255\n\n")

chk <- function(label, observed, expected, tol = 0.05) {
  diff <- abs(observed - expected)
  flag <- if (diff <= tol) "OK " else "OFF"
  cat(sprintf("  [%s] %s: observed %.2f, expected %.2f, diff %.3f (tol %.2f)\n",
              flag, label, observed, expected, diff, tol))
}

chk("D01 mean offset", plot_df[decile == 1L,  mean_offset],  9.40)
chk("D04 mean offset", plot_df[decile == 4L,  mean_offset], -0.78)
chk("D10 mean offset", plot_df[decile == 10L, mean_offset], -3.21)
chk("n_tiles sum",     sum(plot_df$n_tiles),                18255, tol = 5)

log_progress("suppfig_25_site5_cover_decile.R complete.")
