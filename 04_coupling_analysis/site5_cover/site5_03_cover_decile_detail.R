#!/usr/bin/env Rscript
# =============================================================================
# site5_03_cover_decile_detail.R
#
# step 03 — investigate whether Site 5 (usda_sc, ms=5)'s extreme
# b_cover (-17.23) is a real cover-slope signal or an artifact.
#
# Motivation:
#   step 02 LOO showed that removing Site 5 from the 16-non-flagged
#   frame drops r(b_cover, RE_alt18) from -0.570 to -0.182. Site 5
#   single-handedly carries the apparent post-flagging mechanism
#   signal. Before settling on a Discussion-paragraph framing we want
#   to know whether Site 5's b_cover is:
#     (a) A real cover-driven ALS underestimation that escapes flagging
#         because pooled offset cancels (steep slope + balanced cover
#         distribution -> small mean), or
#     (b) An artifact of LC heterogeneity, extreme cover leverage, or
#         acquisition-specific weirdness.
#
# Diagnostics (footprint-level work all already done by step 04; this
# script reads the tile_summary CSV and analyses):
#
#   1. Site 5 detail figure (step 01 style: tile offset map, tile CHM
#      error map, offset vs cover, offset vs slope) + cover-distribution
#      histogram.
#   2. LC-stratified b_cover at Site 5: fit offset ~ cover separately
#      within each dominant_lc class; report b_cover per stratum.
#      If the slope holds within ENF alone (Site 5's site-LC), it's
#      a real within-LC cover effect. If b_cover is heterogeneous
#      across strata, it's confounded with LC.
#   3. Cover-decile aggregation comparison: Site 5 vs Site 3 (flagged,
#      |b_cover|=19.4, our "real cover-slope template") vs Site 4
#      (HARV, |b_cover|=4.88, our "weak mechanism template"). Tile-mean
#      offset by cover decile with SE bars. Real signal -> monotonic
#      decile pattern matching Site 3's shape. Artifact -> noisy or
#      leverage-driven decile pattern.
#
# Inputs:
#   - manuscript_tables/groundwork_task6_phase2_tile_summary.csv
#
# Outputs (manuscript_tables/):
#   - groundwork_task6_phase1_5b_site5_lc_strata.csv
#   - groundwork_task6_phase1_5b_cover_decile_comparison.csv
#
# Outputs (plots/groundwork/):
#   - task6_phase1_5b_site5_detail.pdf          5-panel Site 5 detail
#   - task6_phase1_5b_site5_lc_strata.pdf       LC-stratified scatter
#   - task6_phase1_5b_cover_decile_comparison.pdf  3-site decile compare
#
# Run: source("site5_03_cover_decile_detail.R")
# Wall-clock: ~15 seconds.
# =============================================================================

source("analysis_config.R")
source("analysis_utils.R")

log_section("Step 03: Site 5 (usda_sc) deep-dive")

set.seed(2026)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
})

MANUSCRIPT_TB <- file.path(PROJECT_ROOT, "manuscript_tables")
PLOT_DIR      <- file.path(PROJECT_ROOT, "plots", "groundwork")

TARGET_MS         <- 5L
COMPARISON_MS_HI  <- 3L  # neon_sawb — flagged, b_cover ≈ -19.4 (real cover-slope template)
COMPARISON_MS_LO  <- 4L  # neon_harv2019 — ok, b_cover ≈ -4.88 (weak-mechanism template)

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: Read tile summary, filter to valid tiles
# ─────────────────────────────────────────────────────────────────────────────

log_subsection("Reading tile summary")

ts <- fread(file.path(MANUSCRIPT_TB, "groundwork_task6_phase2_tile_summary.csv"))
log_progress(sprintf("  Tile summary: %s rows total",
                     format(nrow(ts), big.mark = ",")))

ts_v <- ts[valid_tile == TRUE]
log_progress(sprintf("  Valid tiles (>= min footprints): %s",
                     format(nrow(ts_v), big.mark = ",")))

s5 <- ts_v[manuscript_site == TARGET_MS]
s3 <- ts_v[manuscript_site == COMPARISON_MS_HI]
s4 <- ts_v[manuscript_site == COMPARISON_MS_LO]

log_progress(sprintf("  Site 5 valid tiles: %s", format(nrow(s5), big.mark = ",")))
log_progress(sprintf("  Site 3 valid tiles: %s (flagged, steep b_cover)",
                     format(nrow(s3), big.mark = ",")))
log_progress(sprintf("  Site 4 valid tiles: %s (HARV control)",
                     format(nrow(s4), big.mark = ",")))

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Site 5 detail figure (5 panels)
# ─────────────────────────────────────────────────────────────────────────────

log_subsection("Generating Site 5 detail figure")

s5[, x_km := (tile_centroid_x - min(tile_centroid_x)) / 1000]
s5[, y_km := (tile_centroid_y - min(tile_centroid_y)) / 1000]

# Symmetric color scales clipped to robust 2/98 quantiles
off_lim <- c(-1, 1) * max(abs(quantile(s5$offset_mean,  c(0.02, 0.98), na.rm = TRUE)))
err_lim <- c(-1, 1) * max(abs(quantile(s5$chm_err_mean, c(0.02, 0.98), na.rm = TRUE)))

# Tile width in km (assume tile_size_m is consistent; otherwise default 0.25)
tile_w_km <- if ("tile_size_m" %in% names(s5)) s5$tile_size_m[1] / 1000 else 0.25

pA <- ggplot(s5, aes(x_km, y_km, fill = offset_mean)) +
  geom_tile(width = tile_w_km, height = tile_w_km) +
  scale_fill_gradient2(low = "#b2182b", mid = "#f7f7f7", high = "#2166ac",
                       midpoint = 0, limits = off_lim, oob = scales::squish,
                       name = "ALS p90 -\nGEDI rh_98 (m)") +
  coord_equal() +
  labs(x = "x (km from SW corner)", y = "y (km)",
       title = "A. Tile offset (negative = ALS underestimates)") +
  theme_cowplot(10) +
  theme(plot.title = element_text(size = 11),
        legend.position = "right")

pB <- ggplot(s5, aes(x_km, y_km, fill = chm_err_mean)) +
  geom_tile(width = tile_w_km, height = tile_w_km) +
  scale_fill_gradient2(low = "#762a83", mid = "#f7f7f7", high = "#1b7837",
                       midpoint = 0, limits = err_lim, oob = scales::squish,
                       name = "P3D - ALS\nCHM err (m)") +
  coord_equal() +
  labs(x = "x (km)", y = "y (km)",
       title = "B. Tile CHM error (P3D - ALS)") +
  theme_cowplot(10) +
  theme(plot.title = element_text(size = 11),
        legend.position = "right")

pC <- ggplot(s5, aes(x = cover_mean, y = offset_mean)) +
  geom_hex(bins = 60) +
  scale_fill_continuous(type = "viridis", trans = "log10", name = "tiles") +
  geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 0.6,
              linetype = "dashed") +
  geom_hline(yintercept = 0, color = "grey40", linetype = "dotted") +
  labs(x = "Tile mean canopy cover", y = "Tile mean offset (m)",
       title = "C. Offset vs canopy cover (hex density)") +
  theme_cowplot(10) +
  theme(plot.title = element_text(size = 11))

pD <- ggplot(s5, aes(x = slope_mean, y = offset_mean)) +
  geom_hex(bins = 60) +
  scale_fill_continuous(type = "viridis", trans = "log10", name = "tiles") +
  geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 0.6,
              linetype = "dashed") +
  geom_hline(yintercept = 0, color = "grey40", linetype = "dotted") +
  labs(x = "Tile mean slope (deg)", y = "Tile mean offset (m)",
       title = "D. Offset vs terrain slope (hex density)") +
  theme_cowplot(10) +
  theme(plot.title = element_text(size = 11))

pE <- ggplot(s5, aes(x = cover_mean)) +
  geom_histogram(bins = 50, fill = "#4393c3", color = "grey20", linewidth = 0.2) +
  labs(x = "Tile mean canopy cover", y = "Number of tiles",
       title = "E. Cover distribution at Site 5") +
  theme_cowplot(10) +
  theme(plot.title = element_text(size = 11))

# Build the layout: row1 = maps, row2 = scatter+hist
maps_row    <- plot_grid(pA, pB, ncol = 2, align = "hv")
scatter_row <- plot_grid(pC, pD, pE, ncol = 3, align = "hv", rel_widths = c(1, 1, 0.8))

site5_title <- ggdraw() +
  draw_label(sprintf("Site %d (usda_sc) — flag=ok — LC=ENF — n_tiles=%s — pooled offset=%+.2f m, sd=%.2f m, b_cover=%+.2f",
                     TARGET_MS, format(nrow(s5), big.mark = ","),
                     mean(s5$offset_mean), sd(s5$offset_mean), -17.23),
             fontface = "bold", x = 0, hjust = 0, size = 12) +
  theme(plot.margin = margin(t = 4, l = 8))

site5_fig <- plot_grid(site5_title, maps_row, scatter_row,
                        ncol = 1, rel_heights = c(0.05, 1, 0.7))

ggsave(file.path(PLOT_DIR, "task6_phase1_5b_site5_detail.pdf"),
       site5_fig, width = 16, height = 11, bg = "white")
log_progress("  Wrote task6_phase1_5b_site5_detail.pdf")

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: LC-stratified b_cover at Site 5
# ─────────────────────────────────────────────────────────────────────────────

log_subsection("LC-stratified b_cover within Site 5")

# Tabulate dominant_lc presence at Site 5
lc_counts <- s5[, .(n_tiles = .N), by = dominant_lc][order(-n_tiles)]
log_progress("  LC composition of Site 5's valid tiles:")
for (i in seq_len(nrow(lc_counts))) {
  log_progress(sprintf("    %s: %s tiles (%.1f%%)",
                       lc_counts$dominant_lc[i],
                       format(lc_counts$n_tiles[i], big.mark = ","),
                       100 * lc_counts$n_tiles[i] / nrow(s5)))
}

# Fit offset ~ cover within each LC class with >= 50 tiles
strata_rows <- list()
for (i in seq_len(nrow(lc_counts))) {
  lc_i <- lc_counts$dominant_lc[i]
  if (is.na(lc_i)) next
  if (lc_counts$n_tiles[i] < 50L) next
  sub <- s5[dominant_lc == lc_i & is.finite(cover_mean) & is.finite(offset_mean)]
  if (nrow(sub) < 10L) next
  fit <- tryCatch(
    lm(offset_mean ~ cover_mean, data = sub),
    error = function(e) NULL
  )
  if (is.null(fit)) next
  smry <- summary(fit)
  cf   <- smry$coefficients
  cover_idx <- which(rownames(cf) == "cover_mean")
  if (length(cover_idx) == 0L) next
  strata_rows[[length(strata_rows) + 1]] <- data.table(
    manuscript_site = TARGET_MS,
    site_short      = "usda_sc",
    dominant_lc     = lc_i,
    n_tiles         = nrow(sub),
    cover_mean_avg  = mean(sub$cover_mean,  na.rm = TRUE),
    cover_mean_sd   = sd(sub$cover_mean,   na.rm = TRUE),
    cover_mean_min  = min(sub$cover_mean,  na.rm = TRUE),
    cover_mean_max  = max(sub$cover_mean,  na.rm = TRUE),
    offset_mean_avg = mean(sub$offset_mean, na.rm = TRUE),
    intercept       = cf[1, "Estimate"],
    b_cover         = cf[cover_idx, "Estimate"],
    se_b_cover      = cf[cover_idx, "Std. Error"],
    t_b_cover       = cf[cover_idx, "t value"],
    p_b_cover       = cf[cover_idx, "Pr(>|t|)"],
    r2_adj          = smry$adj.r.squared
  )
}

# Pooled fit at Site 5 for reference (without LC stratification)
pooled_fit <- lm(offset_mean ~ cover_mean,
                 data = s5[is.finite(cover_mean) & is.finite(offset_mean)])
pooled_smry <- summary(pooled_fit)
pooled_cf   <- pooled_smry$coefficients

strata_dt <- rbindlist(strata_rows, use.names = TRUE, fill = TRUE)
strata_dt <- rbind(
  strata_dt,
  data.table(
    manuscript_site = TARGET_MS,
    site_short      = "usda_sc",
    dominant_lc     = "ALL_POOLED",
    n_tiles         = nrow(s5),
    cover_mean_avg  = mean(s5$cover_mean,  na.rm = TRUE),
    cover_mean_sd   = sd(s5$cover_mean,   na.rm = TRUE),
    cover_mean_min  = min(s5$cover_mean,  na.rm = TRUE),
    cover_mean_max  = max(s5$cover_mean,  na.rm = TRUE),
    offset_mean_avg = mean(s5$offset_mean, na.rm = TRUE),
    intercept       = pooled_cf[1, "Estimate"],
    b_cover         = pooled_cf["cover_mean", "Estimate"],
    se_b_cover      = pooled_cf["cover_mean", "Std. Error"],
    t_b_cover       = pooled_cf["cover_mean", "t value"],
    p_b_cover       = pooled_cf["cover_mean", "Pr(>|t|)"],
    r2_adj          = pooled_smry$adj.r.squared
  ),
  use.names = TRUE, fill = TRUE
)

setorder(strata_dt, -n_tiles)
fwrite(strata_dt,
       file.path(MANUSCRIPT_TB, "groundwork_task6_phase1_5b_site5_lc_strata.csv"))
log_progress("  Wrote site5_lc_strata.csv")

cat("\nSite 5 LC-stratified b_cover:\n")
print(strata_dt[, .(LC = dominant_lc, n_tiles,
                    cover_range = sprintf("%.2f-%.2f", cover_mean_min, cover_mean_max),
                    cover_mean = round(cover_mean_avg, 3),
                    offset_mean = round(offset_mean_avg, 2),
                    b_cover = round(b_cover, 2),
                    se_bcov = round(se_b_cover, 3),
                    t_bcov = round(t_b_cover, 1),
                    R2adj = round(r2_adj, 3))])

# Generate LC-stratified scatter figure
strata_fits <- strata_dt[dominant_lc != "ALL_POOLED" & !is.na(b_cover)]
plot_lcs <- strata_fits$dominant_lc

s5_for_plot <- s5[dominant_lc %in% plot_lcs]
s5_for_plot[, lc_label := sprintf("%s (n=%d, b_cover=%+.2f)",
                                  dominant_lc,
                                  strata_fits$n_tiles[match(dominant_lc, strata_fits$dominant_lc)],
                                  strata_fits$b_cover[match(dominant_lc, strata_fits$dominant_lc)])]

p_lc <- ggplot(s5_for_plot, aes(x = cover_mean, y = offset_mean)) +
  geom_hex(bins = 50) +
  scale_fill_continuous(type = "viridis", trans = "log10", name = "tiles") +
  geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 0.6,
              linetype = "dashed") +
  geom_hline(yintercept = 0, color = "grey40", linetype = "dotted") +
  facet_wrap(~ lc_label, scales = "free_x") +
  labs(x = "Tile mean canopy cover", y = "Tile mean offset (m)",
       title = sprintf("Site 5 (usda_sc) cover-offset relationship by dominant LC class"),
       subtitle = sprintf("Pooled b_cover at Site 5 = %+.2f. If real cover effect, b_cover should be similar across strata.",
                          pooled_cf["cover_mean", "Estimate"])) +
  theme_cowplot(11) +
  theme(strip.background = element_rect(fill = "grey90"),
        strip.text = element_text(size = 10))

ggsave(file.path(PLOT_DIR, "task6_phase1_5b_site5_lc_strata.pdf"),
       p_lc, width = 14, height = 4 * ceiling(length(plot_lcs) / 3), bg = "white")
log_progress("  Wrote task6_phase1_5b_site5_lc_strata.pdf")

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4: Cover-decile comparison: Site 5 vs Site 3 (flagged) vs Site 4 (HARV)
# ─────────────────────────────────────────────────────────────────────────────

log_subsection("Cover-decile aggregation comparison")

# Compute decile aggregation per site, including SE on offset mean per decile
compute_deciles <- function(dt, ms, label) {
  d <- dt[is.finite(cover_mean) & is.finite(offset_mean)]
  d <- d[, decile := cut(cover_mean,
                          breaks = quantile(cover_mean, probs = seq(0, 1, 0.1), na.rm = TRUE),
                          include.lowest = TRUE,
                          labels = sprintf("D%02d", 1:10))]
  d_agg <- d[!is.na(decile),
             .(n_tiles    = .N,
               cover_med  = median(cover_mean, na.rm = TRUE),
               cover_lo   = quantile(cover_mean, 0.25, na.rm = TRUE),
               cover_hi   = quantile(cover_mean, 0.75, na.rm = TRUE),
               offset_mean = mean(offset_mean, na.rm = TRUE),
               offset_se   = sd(offset_mean,   na.rm = TRUE) / sqrt(.N)),
             by = decile]
  d_agg[, decile_num := as.integer(sub("D", "", as.character(decile)))]
  d_agg[, manuscript_site := ms]
  d_agg[, site_label := label]
  setorder(d_agg, decile_num)
  d_agg
}

dec_s5 <- compute_deciles(s5, TARGET_MS, "Site 5: usda_sc (ok, b_cover=-17.2)")
dec_s3 <- compute_deciles(s3, COMPARISON_MS_HI, "Site 3: neon_sawb (FLAGGED, b_cover=-19.4)")
dec_s4 <- compute_deciles(s4, COMPARISON_MS_LO, "Site 4: neon_harv2019 (ok, b_cover=-4.9)")

dec_combined <- rbind(dec_s5, dec_s3, dec_s4, use.names = TRUE)

fwrite(dec_combined,
       file.path(MANUSCRIPT_TB, "groundwork_task6_phase1_5b_cover_decile_comparison.csv"))
log_progress("  Wrote cover_decile_comparison.csv")

# Make site_label a factor with deliberate order
dec_combined[, site_label := factor(site_label,
  levels = c("Site 3: neon_sawb (FLAGGED, b_cover=-19.4)",
             "Site 5: usda_sc (ok, b_cover=-17.2)",
             "Site 4: neon_harv2019 (ok, b_cover=-4.9)"))]

p_dec <- ggplot(dec_combined,
                 aes(x = cover_med, y = offset_mean,
                     ymin = offset_mean - 1.96 * offset_se,
                     ymax = offset_mean + 1.96 * offset_se,
                     color = site_label, group = site_label)) +
  geom_hline(yintercept = 0, color = "grey50", linetype = "dotted") +
  geom_line(linewidth = 0.6) +
  geom_errorbar(width = 0.015, linewidth = 0.4) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c(
    "Site 3: neon_sawb (FLAGGED, b_cover=-19.4)" = "#b2182b",
    "Site 5: usda_sc (ok, b_cover=-17.2)"        = "#2166ac",
    "Site 4: neon_harv2019 (ok, b_cover=-4.9)"   = "#4d9221"
  ), name = "Site (flag, b_cover)") +
  labs(x = "Tile mean canopy cover (decile median)",
       y = "Tile mean offset (m), 95% CI",
       title = "Cover-decile aggregation: is Site 5's slope shape like a flagged site or a control?",
       subtitle = "Real cover-driven mechanism -> Site 5 (blue) tracks Site 3 (red) shape.\nLeverage artifact -> Site 5 looks like Site 4 (green) with extreme deciles pulling the regression.") +
  theme_cowplot(11) +
  theme(legend.position = "top",
        legend.text = element_text(size = 9))

ggsave(file.path(PLOT_DIR, "task6_phase1_5b_cover_decile_comparison.pdf"),
       p_dec, width = 11, height = 7, bg = "white")
log_progress("  Wrote cover_decile_comparison.pdf")

# Print decile summary
cat("\nCover-decile summary (offset_mean by decile):\n")
print(dcast(dec_combined,
            decile_num ~ site_label,
            value.var = "offset_mean",
            fun.aggregate = function(x) round(mean(x), 2)))

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5: Verdict heuristic
# ─────────────────────────────────────────────────────────────────────────────

log_subsection("Verdict heuristic")

# Test 1: Is Site 5's b_cover within ENF (its dominant LC) similar in magnitude
# to the pooled b_cover? If yes (within ~30%), the slope holds within LC.
enf_row <- strata_dt[dominant_lc == "ENF"]
pooled_b <- strata_dt[dominant_lc == "ALL_POOLED"]$b_cover
if (nrow(enf_row) > 0L) {
  enf_b <- enf_row$b_cover[1]
  ratio <- enf_b / pooled_b
  cat(sprintf("\nTest 1 (LC robustness): b_cover within ENF = %+.2f, pooled = %+.2f, ratio = %.2f\n",
              enf_b, pooled_b, ratio))
  if (abs(ratio - 1) < 0.3) {
    cat("  -> ENF-only slope is similar to pooled. LC heterogeneity is NOT the explanation.\n")
  } else if (abs(ratio) < 0.5) {
    cat("  -> ENF-only slope is much weaker than pooled. LC heterogeneity DRIVES the slope.\n")
  } else {
    cat("  -> ENF-only slope differs from pooled but in the same direction. Mixed signal.\n")
  }
}

# Test 2: Decile pattern shape — is Site 5's last-decile drop (D9->D10) larger
# than the median drop, suggesting leverage? Compare to Site 3's pattern.
last_drop <- function(d) {
  v <- d$offset_mean[order(d$decile_num)]
  v[length(v)] - v[length(v) - 1]
}
all_drops <- function(d) {
  v <- d$offset_mean[order(d$decile_num)]
  diff(v)
}
s5_last <- last_drop(dec_s5)
s5_drops <- all_drops(dec_s5)
s5_drop_ratio <- s5_last / median(s5_drops, na.rm = TRUE)

s3_last <- last_drop(dec_s3)
s3_drops <- all_drops(dec_s3)
s3_drop_ratio <- s3_last / median(s3_drops, na.rm = TRUE)

cat(sprintf("\nTest 2 (decile shape): Site 5 last-decile drop / median drop = %.2f\n", s5_drop_ratio))
cat(sprintf("                       Site 3 last-decile drop / median drop = %.2f (template)\n", s3_drop_ratio))
if (abs(s5_drop_ratio) < 2 * abs(s3_drop_ratio)) {
  cat("  -> Site 5's decile pattern is comparable to Site 3 (flagged template). Slope shape is consistent with a real cover effect.\n")
} else {
  cat("  -> Site 5's last-decile drop is much larger than Site 3's. Leverage from extreme cover.\n")
}

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 6: Final summary
# ─────────────────────────────────────────────────────────────────────────────

log_section("step 03 — Complete")

log_progress("Output files:")
log_progress(sprintf("  %s/groundwork_task6_phase1_5b_site5_lc_strata.csv",       MANUSCRIPT_TB))
log_progress(sprintf("  %s/groundwork_task6_phase1_5b_cover_decile_comparison.csv", MANUSCRIPT_TB))
log_progress(sprintf("  %s/task6_phase1_5b_site5_detail.pdf",                     PLOT_DIR))
log_progress(sprintf("  %s/task6_phase1_5b_site5_lc_strata.pdf",                  PLOT_DIR))
log_progress(sprintf("  %s/task6_phase1_5b_cover_decile_comparison.pdf",          PLOT_DIR))
log_progress("")
log_progress("Interpretation guide:")
log_progress("  * If LC-stratified b_cover within ENF (Site 5's dominant LC) is similar")
log_progress("    in magnitude to the pooled b_cover -> real cover-slope signal.")
log_progress("  * If decile-aggregation shape at Site 5 looks like Site 3 (flagged)")
log_progress("    rather than Site 4 (HARV) -> real cover-slope mechanism that")
log_progress("    escapes Track A flagging because pooled offset cancels.")
log_progress("    -> Discussion framing: 'partial mechanism plus a flag-rule")
log_progress("       refinement opportunity for slope-based criteria'.")
log_progress("  * If both tests fail -> Site 5's high-leverage role is artifactual,")
log_progress("    drop Site 5 from the cover-slope summary, narrative becomes")
log_progress("    'cover-slope mechanism applies to flagged sites, not to the")
log_progress("    residual ~12% V_site at non-flagged sites'.")
log_progress("")
log_progress("Outputs:")
log_progress("  - The two CSVs above")
log_progress("  - The three PDFs above")
log_progress("  - This script's stdout/stderr log")
