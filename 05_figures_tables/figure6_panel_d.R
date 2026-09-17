#!/usr/bin/env Rscript
# ============================================================================
# figure6_panel_d.R
#
# Figure 6 panel (d): RH98 fixed-effect coefficient under the operational
# 16-site primary CHM fit ("P3D CHM") vs the 15-site counterfactual reference-
# substitution refit ("Counterfactual"). Visualizes the sign-flip that is the
# chapter's clearest single piece of mechanistic evidence for the DSM-DTM
# compensation channel (§3.4, §4.2, Supplementary Text S2.3).
#
# Numerics: brms::fixef(fit, robust = TRUE)["rh_98_z", ] for each fit.
#   - fit_chm_16:  16-site primary CHM, loaded from 10_models_stage2.rds.
#   - fit_alt_15:  15-site err_alt CHM, loaded from
#                  groundwork_task4_refit.rds$data$fit_alt_15.
#
# Cross-checked 2026-05-21 against chapter1_v06.docx Table 5 row 1 (RH98).
#
# Cluster convention (mirrors supp_figs_rerender_fixes.R):
#   - Set DISPLAY = "" and options(device = pdf) to bypass the X11-dependent
#     png() device.
#   - Output is a PDF written to /gpfs/data1/vclgp/lmaden/chpt1/plots.
#   - Convert to PNG downstream for Word embedding via
#     `pdftoppm -png -r 300 figure6_panel_d.pdf figure6_panel_d`
#     (writes figure6_panel_d-1.png).
#
# Notes on prior-version fixes:
#   - Y-axis now uses explicit breaks + a 1-decimal formatter (was rendering
#     in scientific notation with 1.11e-16 at the zero tick).
#   - Subtitle is ASCII-only (R's default PDF Type 1 fonts lack glyphs for
#     U+0394 Delta and U+2248 approx).
#   - X-axis tick label "err_alt" renamed to "Counterfactual" to mirror
#     chapter prose and pair cleanly with "P3D CHM".
#
# RUN: from /gpfs/data1/vclgp/lmaden/chpt1
#   Rscript figure6_panel_d.R
# ============================================================================

Sys.setenv(DISPLAY = "")
options(device = pdf)

suppressPackageStartupMessages({
  library(ggplot2)
  library(tibble)
})

# ----------------------------------------------------------------------------
# 1. Output directory
# ----------------------------------------------------------------------------

out_dir  <- "/gpfs/data1/vclgp/lmaden/chpt1/plots"
out_path <- file.path(out_dir, "figure6_panel_d.pdf")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------------------------------------------------------
# 2. Posterior summaries
# ----------------------------------------------------------------------------
# Option A (preferred, fully reproducible): extract directly from loaded fits.
# Uncomment if fit_chm_16 and fit_alt_15 are in the current R session.
#
# library(brms)
# operational <- fixef(fit_chm_16, robust = TRUE)["rh_98_z",
#                                                 c("Estimate", "Q2.5", "Q97.5")]
# counterfact <- fixef(fit_alt_15, robust = TRUE)["rh_98_z",
#                                                 c("Estimate", "Q2.5", "Q97.5")]
#
# Option B (default): hard-coded constants from the audit above. These match
# the cluster-side fixef() output to four decimals and the chapter Table 5
# row 1 to two decimals.

operational <- c(Estimate = -0.4406, Q2.5 = -0.4766, Q97.5 = -0.4060)
counterfact <- c(Estimate = +0.3934, Q2.5 = +0.3472, Q97.5 = +0.4390)

OP_LABEL <- "P3D CHM\n(uses CNN DTM)"
CF_LABEL <- "Counterfactual\n(uses 3DEP DTM)"

panel_data <- tibble(
  fit_label = c(OP_LABEL, CF_LABEL),
  fit_order = c(1L, 2L),
  median    = c(operational[["Estimate"]], counterfact[["Estimate"]]),
  q025      = c(operational[["Q2.5"]],     counterfact[["Q2.5"]]),
  q975      = c(operational[["Q97.5"]],    counterfact[["Q97.5"]])
)

shift <- counterfact[["Estimate"]] - operational[["Estimate"]]

# ----------------------------------------------------------------------------
# 3. Plot
# ----------------------------------------------------------------------------
# Color palette aligned with progress-report Figure 2C (P3D blue,
# counterfactual green).

bar_fills <- c(
  "P3D CHM\n(uses CNN DTM)"        = "#4A6FA5",
  "Counterfactual\n(uses 3DEP DTM)" = "#3F8657"
)

# Outside-bar label positions: above positive bars (vjust = 0), below
# negative bars (vjust = 1), offset by 0.045 m/SD beyond each error bar.

panel_data$label_y <- with(panel_data,
                           ifelse(median > 0, q975 + 0.045, q025 - 0.045))
panel_data$label_v <- with(panel_data,
                           ifelse(median > 0, 0,           1))

# Subtitle: ASCII-only so it renders correctly in the default PDF device.
subtitle_txt <- sprintf("Sign-flipped (shift: %+0.2f m per SD)", shift)

# Y-axis: explicit clean breaks AND an explicit 1-decimal formatter so the
# zero tick reads "0.0" and no tick renders in scientific notation.
y_breaks <- c(-0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6)
y_fmt    <- function(x) sprintf("%.1f", x)

panel_d <- ggplot(panel_data,
                  aes(x = reorder(fit_label, fit_order),
                      y = median,
                      fill = fit_label)) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
  geom_col(width = 0.55, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = q025, ymax = q975),
                width = 0.18, color = "black", linewidth = 0.4) +
  geom_text(aes(y = label_y,
                label = sprintf("%+0.2f", median),
                vjust = label_v),
            size = 4.2, fontface = "bold") +
  scale_fill_manual(values = bar_fills) +
  scale_y_continuous(limits = c(-0.65, 0.65),
                     breaks = y_breaks,
                     labels = y_fmt,
                     expand = c(0, 0)) +
  labs(
    title    = NULL,
    subtitle = subtitle_txt,
    tag      = "(d)",
    x        = NULL,
    y        = "RH98 coefficient (m per SD)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.tag           = element_text(face = "bold", size = 14),
    plot.tag.position  = c(0.05, 0.97),
    plot.subtitle      = element_text(face = "italic",
                                      color = "#B22222",
                                      size = 10,
                                      hjust = 0.5,
                                      margin = margin(t = 0, b = 8)),
    plot.margin        = margin(t = 10, r = 14, b = 8, l = 14),
    axis.title.y       = element_text(margin = margin(r = 6)),
    axis.text.x        = element_text(size = 10, lineheight = 0.9),
    axis.text.y        = element_text(size = 10),
    axis.line          = element_line(linewidth = 0.4),
    axis.ticks         = element_line(linewidth = 0.4),
    legend.position    = "none"
  )

# ----------------------------------------------------------------------------
# 4. Save
# ----------------------------------------------------------------------------

ggsave(
  filename = out_path,
  plot     = panel_d,
  width    = 4.0,
  height   = 4.5,
  bg       = "white",
  device   = pdf            # no X11 dependency
)

message("Wrote ", normalizePath(out_path))
message("Operational rh_98_z: ",
        sprintf("%+0.4f [%+0.4f, %+0.4f]",
                operational[["Estimate"]],
                operational[["Q2.5"]],
                operational[["Q97.5"]]))
message("Counterfactual rh_98_z: ",
        sprintf("%+0.4f [%+0.4f, %+0.4f]",
                counterfact[["Estimate"]],
                counterfact[["Q2.5"]],
                counterfact[["Q97.5"]]))
message("Shift:               ",
        sprintf("%+0.4f m per SD (sign-flipped)", shift))
message("")
message("Next step (downstream, for Word embedding):")
message("  cd /gpfs/data1/vclgp/lmaden/chpt1/plots")
message("  pdftoppm -png -r 300 figure6_panel_d.pdf figure6_panel_d")
