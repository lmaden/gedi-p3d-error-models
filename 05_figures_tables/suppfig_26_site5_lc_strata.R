# =====================================================================
# suppfig_26_site5_lc_strata.R
# Supp Figure S26 (NEW): Site 5 b_cover forest plot stratified by
# dominant land cover. Layered companion to Supp Figure S25.
#
# Site 5 LC-stratified deep dive.
#
# Strata per story-lock §2J (5 of 10 strata in source CSV; 5 smaller
# strata exist - WTR, IMP, CWL, GRS, and an 843-tile unclassified set -
# noted in the docx caption but below §2J reporting scope):
#   ALL_POOLED:  n = 18,255   b_cover = -17.18   r_adj^2 = 0.366
#   ENF:         n =  9,869   b_cover = -16.83   r_adj^2 = 0.364
#   IWL:         n =  3,419   b_cover = -18.85   r_adj^2 = 0.396
#   RCP:         n =  1,856   b_cover = -15.46   r_adj^2 = 0.338
#   BDF:         n =    373   b_cover = -14.93   r_adj^2 = 0.330
#
# Notes:
# - The tile-level CSV `groundwork_task6_phase1_tile_summary.csv` has
#   no Site 5 rows (manuscript_site values 1..4 only). The two
#   Site-5-specific candidate filenames (`*_site5_tile_level.csv` and
#   `*_site5_tiles.csv`) also do not exist on cluster. So this figure
#   renders the summary-only forest-plot path; the hexbin-scatter
#   branch in the prior draft is removed for clarity.
# - The summary CSV has 10 strata; we filter explicitly to the 5 §2J
#   strata rather than relying on factor-level NA filtering.
# - 95% CIs on b_cover added from se_b_cover (column present in CSV;
#   delete the geom_errorbarh block to revert to point-only display).
# - Color cue for LC category added (forest = ENF/BDF in dark green;
#   non-forest = IWL/RCP in burnt orange; ALL_POOLED in dark gray).
#   Rationale: the manuscript's forest-functional-type framework
#   (BDF/ENF/DNF/EBF; Figure 4, Table 6) maps onto only ENF and BDF
#   at Site 5 (DNF and EBF are absent in temperate CONUS South
#   Carolina); the color cue makes the forest-vs-non-forest contrast
#   visible without dropping the 5,275 non-forest tiles from view.
#
# Source:
#   - groundwork_task6_phase1_5b_site5_lc_strata.csv
# =====================================================================

source("fig_common.R")
fig_banner("Supp Figure S26",
                 "Site 5 LC-stratified b_cover forest plot")

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(dplyr)
})

# ---------------------------------------------------------------------
# 1. Load summary CSV
# ---------------------------------------------------------------------
summary_path <- file.path(manuscript_tables_dir,
                          "groundwork_task6_phase1_5b_site5_lc_strata.csv")
if (!file.exists(summary_path)) stop("Required: ", summary_path)
summ <- as.data.table(fread(summary_path))
log_progress(sprintf("  Loaded %d rows, %d cols", nrow(summ), ncol(summ)))
cat("  Columns: ", paste(names(summ), collapse = ", "), "\n")

required_cols <- c("dominant_lc", "n_tiles", "b_cover", "se_b_cover", "r2_adj")
missing_cols <- setdiff(required_cols, names(summ))
if (length(missing_cols) > 0L) {
  stop("Required columns missing from ", basename(summary_path), ": ",
       paste(missing_cols, collapse = ", "),
       "\nFound: ", paste(names(summ), collapse = ", "))
}

summ_std <- data.table(
  stratum = as.character(summ$dominant_lc),
  b       = as.numeric(summ$b_cover),
  se_b    = as.numeric(summ$se_b_cover),
  n       = as.integer(summ$n_tiles),
  r2_adj  = as.numeric(summ$r2_adj)
)

log_progress("All strata in source CSV (sorted by n, descending):")
print(summ_std[order(-n)][, .(stratum,
                              b = round(b, 2),
                              se_b = round(se_b, 3),
                              r2_adj = round(r2_adj, 3),
                              n)])

# ---------------------------------------------------------------------
# 2. Filter to 5 §2J strata, assign LC category, compute CI
# ---------------------------------------------------------------------
target_strata <- c("ALL_POOLED", "ENF", "IWL", "RCP", "BDF")
summ_plot <- copy(summ_std)[stratum %in% target_strata]
if (nrow(summ_plot) != 5L) {
  stop("Expected all 5 §2J strata (",
       paste(target_strata, collapse = ", "),
       ") in CSV; found ", nrow(summ_plot),
       " (", paste(summ_plot$stratum, collapse = ", "), ").")
}

# Order ALL_POOLED at top of y-axis: ggplot draws the first factor
# level at the bottom, so reverse the §2J order.
summ_plot[, stratum := factor(stratum, levels = rev(target_strata))]

# LC category for the forest-vs-non-forest color cue.
#   Forest:     ENF (Evergreen Needleleaf), BDF (Broadleaf Deciduous).
#               At Site 5 (USDA South Carolina), these are the only
#               two of the four GEDI forest functional types present
#               (DNF and EBF have zero tiles here).
#   Non-forest: IWL (Inland Wetland), RCP (Row Crop / Pasture).
#   Pooled:     ALL_POOLED (reference; aggregates all LC classes).
summ_plot[, lc_category := fcase(
  as.character(stratum) %in% c("ENF", "BDF"),  "Forest",
  as.character(stratum) %in% c("IWL", "RCP"),  "Non-forest",
  as.character(stratum) == "ALL_POOLED",        "Pooled (all LC)"
)]
summ_plot[, lc_category := factor(lc_category,
                                  levels = c("Pooled (all LC)",
                                             "Forest", "Non-forest"))]

summ_plot[, ci_lo := b - 1.96 * se_b]
summ_plot[, ci_hi := b + 1.96 * se_b]

# Vertical-align the row-end text labels at a constant x slightly past
# the widest CI bar (avoids per-row label drift).
label_x <- max(summ_plot$ci_hi, na.rm = TRUE) + 0.6

log_progress("Strata plotted (5 of 10, §2J-consistent):")
print(summ_plot[order(stratum)][, .(stratum, lc_category,
                                    b = round(b, 2),
                                    ci_lo = round(ci_lo, 2),
                                    ci_hi = round(ci_hi, 2),
                                    r2_adj = round(r2_adj, 3),
                                    n)])

# ---------------------------------------------------------------------
# 3. Forest plot
# ---------------------------------------------------------------------
LC_CATEGORY_PALETTE <- c(
  "Pooled (all LC)" = "#444444",   # dark gray; existing convention
  "Forest"          = "#2e7d32",   # dark green; ENF + BDF
  "Non-forest"      = "#c46210"    # burnt orange; IWL + RCP
)

p <- ggplot(summ_plot, aes(x = b, y = stratum)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = COLOR_ZERO_LINE, linewidth = 0.4) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi, color = lc_category),
                 height = 0.18, linewidth = 0.5) +
  geom_point(aes(color = lc_category), size = 3.5) +
  geom_text(aes(label = sprintf("b = %s%.2f   r_adj^2 = %.3f   n = %s",
                                ifelse(b < 0, MINUS, ""), abs(b),
                                r2_adj,
                                formatC(n, big.mark = ",", format = "d"))),
            x = label_x, hjust = 0,
            size = 3.2, color = "black") +
  scale_color_manual(values = LC_CATEGORY_PALETTE,
                     name = NULL,
                     drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.80))) +
  labs(
    title = "Site 5 cover-vs-offset slope by land-cover stratum",
    x = "b_cover (m per unit canopy cover, with 95% CI)",
    y = NULL
  ) +
  theme_section_G(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    plot.title         = element_text(size = 11, face = "bold"),
    legend.position    = "top",
    legend.box.margin  = margin(0, 0, -4, 0),
    legend.key.width   = unit(0.7, "lines"),
    legend.text        = element_text(size = 9)
  )

save_figure(p, "suppfig_S26_site5_lc_strata",
            width_in = 8.0, height_in = 4.8)

# ---------------------------------------------------------------------
# 4. Verification against story-lock §2J anchors
# ---------------------------------------------------------------------
log_subsection("VERIFICATION (console anchors)")
cat("\nStory-lock §2J targets (b_cover m/unit cover, r_adj^2, n_tiles):\n")
cat("  ALL_POOLED: -17.18, 0.366, 18,255\n")
cat("  ENF:        -16.83, 0.364,  9,869\n")
cat("  IWL:        -18.85, 0.396,  3,419\n")
cat("  RCP:        -15.46, 0.338,  1,856\n")
cat("  BDF:        -14.93, 0.330,    373\n\n")

chk <- function(label, observed, expected, tol = 0.05) {
  diff <- abs(observed - expected)
  flag <- if (diff <= tol) "OK " else "OFF"
  cat(sprintf("  [%s] %s: observed %.2f, expected %.2f, diff %.3f (tol %.2f)\n",
              flag, label, observed, expected, diff, tol))
}

chk("ALL_POOLED b_cover", summ_plot[stratum == "ALL_POOLED", b], -17.18)
chk("ENF        b_cover", summ_plot[stratum == "ENF",        b], -16.83)
chk("IWL        b_cover", summ_plot[stratum == "IWL",        b], -18.85)
chk("RCP        b_cover", summ_plot[stratum == "RCP",        b], -15.46)
chk("BDF        b_cover", summ_plot[stratum == "BDF",        b], -14.93)

# LC-category tile-count breakdown (informational, for caption).
cat("\nLC-category breakdown of tiles plotted:\n")
cat_summary <- summ_plot[, .(n_tiles_total = sum(n), n_strata = .N), by = lc_category]
cat_summary <- cat_summary[order(lc_category)]
for (i in seq_len(nrow(cat_summary))) {
  cat(sprintf("  %-18s  n_strata = %d   n_tiles = %s\n",
              cat_summary$lc_category[i],
              cat_summary$n_strata[i],
              formatC(cat_summary$n_tiles_total[i], big.mark = ",", format = "d")))
}

# Note the 5 smaller strata excluded from the figure (for caption).
other_strata <- setdiff(summ_std$stratum, target_strata)
named_other  <- setdiff(other_strata, "")
empty_n      <- sum(summ_std$stratum == "")
cat(sprintf("\n  Note: source CSV contains %d strata; 5 below §2J reporting scope.\n",
            nrow(summ_std)))
if (length(named_other) > 0L) {
  ns <- summ_std[stratum %in% named_other][order(-n)]
  cat("    Named smaller strata (excluded from figure):\n")
  for (i in seq_len(nrow(ns))) {
    cat(sprintf("      %-4s  n=%6d  b_cover=%s%.2f  r_adj^2=%.3f\n",
                ns$stratum[i], ns$n[i],
                ifelse(ns$b[i] < 0, MINUS, ""), abs(ns$b[i]),
                ns$r2_adj[i]))
  }
}
if (empty_n > 0L) {
  emp <- summ_std[stratum == ""]
  cat(sprintf("    Unclassified-LC stratum (dominant_lc empty): n=%d  b_cover=%s%.2f  r_adj^2=%.3f\n",
              emp$n[1], ifelse(emp$b[1] < 0, MINUS, ""), abs(emp$b[1]),
              emp$r2_adj[1]))
}

log_progress("suppfig_26_site5_lc_strata.R complete.")
