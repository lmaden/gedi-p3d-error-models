# =====================================================================
# fig3_nu_posteriors.R
# Figure 3: Student-t degrees-of-freedom posterior distributions.
#
# v04 -> v10 changes:
#   - CHM panel uses 16-site fit (ν median 2.80) instead of
#     19-site (4.49).
#   - DTM panel unchanged from 18-site (ν median 3.16).
#   - APA33 (id=68: "Break the distribution into many more intervals
#     and shorten the domain of the x-axis"): bin count raised to
#     ~50 bins; x-axis range tightened to [2.0, 4.5] which covers both
#     posteriors comfortably without padding.
#   - In-panel annotation: median + 95% CI (lab convention).
#
# Sources:
#   - load_checkpoint("10b_chm_sensitivity")$fit_chm_16
#   - load_checkpoint("10_models_stage2")$fit_dtm_s2
# =====================================================================

source("fig_common.R")
fig_banner("Figure 3",
                 "Student-t nu posteriors (CHM 16-site, DTM 18-site)")

suppressPackageStartupMessages({
  library(brms); library(posterior); library(ggplot2)
  library(data.table); library(dplyr); library(patchwork)
})

# ---------------------------------------------------------------------
# 1. Load fits and extract nu posterior draws
# ---------------------------------------------------------------------
log_subsection("Loading fits and extracting nu posteriors")

chm_ck <- load_checkpoint("10b_chm_sensitivity")
if (!"fit_chm_16" %in% names(chm_ck)) {
  stop("fit_chm_16 not found in 10b_chm_sensitivity checkpoint.")
}
fit_chm_16 <- chm_ck$fit_chm_16

s2 <- load_checkpoint("10_models_stage2")
if (!"fit_dtm_s2" %in% names(s2)) {
  stop("fit_dtm_s2 not found in 10_models_stage2 checkpoint.")
}
fit_dtm <- s2$fit_dtm_s2

extract_nu <- function(fit, label) {
  d <- as_draws_df(fit, variable = "nu")
  if (!"nu" %in% names(d)) {
    # Older brms versions name it differently.
    cand <- grep("^nu", names(d), value = TRUE)
    if (!length(cand)) stop("Could not find 'nu' draws in fit (", label, ")")
    setnames(d, cand[1], "nu")
  }
  out <- data.table(product = label, nu = as.numeric(d$nu))
  q   <- as.numeric(quantile(out$nu, c(0.025, 0.5, 0.975)))
  attr(out, "summary") <- data.table(product = label,
                                     median  = q[2],
                                     lo      = q[1],
                                     hi      = q[3])
  out
}

nu_chm <- extract_nu(fit_chm_16, "CHM (16-site)")
nu_dtm <- extract_nu(fit_dtm,    "DTM (18-site)")
nu     <- rbind(nu_chm, nu_dtm)

summary_dt <- rbindlist(list(attr(nu_chm, "summary"),
                             attr(nu_dtm, "summary")))
log_progress("nu posterior summaries:")
print(summary_dt)

# ---------------------------------------------------------------------
# 2. Build plot
# ---------------------------------------------------------------------
log_subsection("Composing plot")

# Tight x-axis to APA33 spec. Bins of width 0.025 over [2.0, 4.5]
# gives ~100 bins; we use 60 for readability.
x_lo <- 2.0
x_hi <- 4.5
n_bins <- 60

# Annotation text per panel.
annot <- summary_dt[, .(
  product,
  label = sprintf("Median = %.2f\n95%% CI: [%.2f, %.2f]",
                  median, lo, hi)
)]

# Map product -> annotation x position (just inside the rightmost
# fraction of the panel).
annot[, x_pos := x_hi - 0.05]

pal <- c("CHM (16-site)" = COLOR_CHM,
         "DTM (18-site)" = COLOR_DTM)

p <- ggplot(nu, aes(x = nu, fill = product, color = product)) +
  geom_histogram(bins = n_bins, alpha = 0.85,
                 boundary = x_lo, closed = "left",
                 linewidth = 0.15) +
  geom_vline(data = summary_dt,
             aes(xintercept = median),
             color = "black", linetype = "dashed", linewidth = 0.45) +
  geom_text(data = annot,
            aes(x = x_pos, y = Inf, label = label),
            hjust = 1, vjust = 1.4, size = 3.3,
            color = "black",
            inherit.aes = FALSE,
            lineheight = 1.05) +
  facet_wrap(~ product, nrow = 1, scales = "fixed") +
  scale_fill_manual(values = pal, guide = "none") +
  scale_color_manual(values = pal, guide = "none") +
  scale_x_continuous(limits = c(x_lo, x_hi),
                     breaks = seq(x_lo, x_hi, by = 0.5),
                     expand = c(0, 0)) +
  labs(
    x = expression(paste("Student-", italic("t"), " degrees of freedom (", nu, ")")),
    y = "Posterior draws (count)"
  ) +
  theme_section_G(base_size = 11) +
  theme(
    legend.position = "none",
    strip.text      = element_text(size = 11, face = "bold"),
    panel.spacing.x = unit(1.5, "lines")
  )

# ---------------------------------------------------------------------
# 3. Save and verify
# ---------------------------------------------------------------------
save_figure(p, "fig03_nu_posteriors", width_in = 7.0, height_in = 3.6)

log_subsection("VERIFICATION (console anchors)")
cat("\nMust match story-lock §2A (CHM 16-site) and §2B (DTM 18-site):\n")
cat("  CHM 16-site nu: median 2.796, 95% CI [2.740, 2.853]\n")
cat("  DTM 18-site nu: median 3.16,  95% CI [3.11,  3.22]\n\n")
cat("Observed:\n")
summary_dt[, .(product,
               median = round(median, 3),
               lo     = round(lo, 3),
               hi     = round(hi, 3))] |> print()

log_progress("fig3_nu_posteriors.R complete.")
