# =====================================================================
# temporal_mismatch_analysis.R (v3)
#
# PURPOSE:
#   Compute temporal gaps and correlate with error characteristics
#   for BOTH study products:
#     CHM: P3D imagery vs ALS reference (from tracker.xlsx)
#     DTM: P3D imagery vs 3DEP reference (from query_3dep_dates.py)
#
# PREREQUISITES:
#   1. Run extract_p3d_tile_dates.R  -> p3d_tile_dates_by_site.csv
#   2. Run query_3dep_dates.py       -> 3dep_date_summary.csv
#   3. Checkpoints: 10_models_stage2.rds, 01_data_ingest.rds
#
# OUTPUTS:
#   - tables/temporal_gap_chm.csv    (CHM site-level analysis)
#   - tables/temporal_gap_dtm.csv    (DTM site-level analysis)
#   - plots/pub_figSXX_temporal_gap_chm.pdf
#   - plots/pub_figSXX_temporal_gap_dtm.pdf
#   - Console: manuscript-ready values for Sections 2.2, 4.1, 4.8
#
# RUN ON CLUSTER: Rscript temporal_mismatch_analysis.R
# =====================================================================

Sys.setenv(DISPLAY = "")
options(device = pdf)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(brms)
  library(patchwork)
})

# =====================================================================
# CONFIGURATION
# =====================================================================

PROJECT_ROOT   <- Sys.getenv("PROJECT_ROOT", "/gpfs/data1/vclgp/lmaden/chpt1")
CHECKPOINT_DIR <- file.path(PROJECT_ROOT, "checkpoints")
out_plots      <- file.path(PROJECT_ROOT, "plots")
out_tables     <- file.path(PROJECT_ROOT, "tables")
dir.create(out_plots, showWarnings = FALSE, recursive = TRUE)
dir.create(out_tables, showWarnings = FALSE, recursive = TRUE)

safe_pdf <- if (capabilities("cairo")) cairo_pdf else pdf

# Site lists
SITES_CHM <- c(1:9, 11:15, 17:20)          # All sites except 16 (non-forest)
SITES_DTM <- c(1:9, 11:15, 17:20)          # Same, but site 10 also excluded
SITES_EXCLUDE_DTM <- 10                     # Datum offset

# Paper site renumbering: site 16 dropped -> 17=16, 18=17, 19=18, 20=19
paper_site_label <- function(site_num) {
  s <- as.integer(site_num)
  ifelse(s <= 15, s, s - 1L)
}

cat(strrep("=", 65), "\n")
cat("  TEMPORAL MISMATCH ANALYSIS (CHM + DTM)\n")
cat(strrep("=", 65), "\n\n")


# =====================================================================
# REFERENCE DATES
# =====================================================================

# --- ALS dates (CHM reference) from tracker.xlsx ---
# Site 3 corrected from 2104-07-07 to 2014-07-07

als_dates <- data.frame(
  site = 1:20,
  als_date = as.Date(c(
    "2015-02-18",  # 1  usda_me
    "2013-10-21",  # 2  nasa_howland
    "2014-07-07",  # 3  neon_sawb  *** CORRECTED ***
    "2019-10-06",  # 4  neon_harv
    "2015-04-06",  # 5  usda_sc
    "2021-10-26",  # 6  neon_jerc
    "2021-06-06",  # 7  neon_tall
    "2021-06-07",  # 8  neon_dela
    "2021-06-09",  # 9  neon_leno
    "2021-07-30",  # 10 neon_clbj
    "2020-08-16",  # 11 neon_konz
    "2022-07-09",  # 12 neon_stei
    "2022-07-23",  # 13 neon_unde
    "2022-07-20",  # 14 neon_steicheq
    "2021-07-21",  # 15 neon_wood
    "2021-07-14",  # 16 neon_nogp (excluded)
    "2022-08-19",  # 17 neon_ster
    "2021-07-05",  # 18 neon_cper
    "2019-08-02",  # 19 ltbmu
    "2013-08-03"   # 20 neon_sjer
  )),
  site_name = c(
    "usda_me", "nasa_howland", "neon_sawb", "neon_harv", "usda_sc",
    "neon_jerc", "neon_tall", "neon_dela", "neon_leno", "neon_clbj",
    "neon_konz", "neon_stei", "neon_unde", "neon_steicheq", "neon_wood",
    "neon_nogp", "neon_ster", "neon_cper", "ltbmu", "neon_sjer"
  ),
  stringsAsFactors = FALSE
)


# --- 3DEP dates (DTM reference) from query_3dep_dates.py ---

dep_csv <- file.path(out_tables, "3dep_date_summary.csv")

# Also check working directory (if user copied it there)
if (!file.exists(dep_csv)) {
  dep_csv_alt <- "3dep_date_summary.csv"
  if (file.exists(dep_csv_alt)) dep_csv <- dep_csv_alt
}

has_3dep_dates <- file.exists(dep_csv)

if (has_3dep_dates) {
  dep_raw <- read.csv(dep_csv, stringsAsFactors = FALSE)

  # The script outputs: site_id, n_tiles, pub_earliest, pub_latest,
  #                     source_earliest, source_latest
  # "source" dates = lidar collection dates (preferred for temporal gap)
  # "pub" dates    = USGS publication dates (fallback)

  dep_dates <- dep_raw %>%
    rename(site = site_id) %>%
    mutate(
      # Prefer source dates (actual lidar collection); fall back to pub dates
      dep_date_early = as.Date(ifelse(source_earliest != "N/A",
                                       source_earliest, pub_earliest)),
      dep_date_late  = as.Date(ifelse(source_latest != "N/A",
                                       source_latest, pub_latest)),
      # Midpoint as representative date
      dep_date_mid   = dep_date_early + (dep_date_late - dep_date_early) / 2,
      dep_date_source = ifelse(source_earliest != "N/A", "lidar_collection", "publication")
    )

  cat(sprintf("3DEP dates loaded for %d sites (source: %s)\n",
              nrow(dep_dates),
              paste(unique(dep_dates$dep_date_source), collapse = "/")))
} else {
  cat("WARNING: 3dep_date_summary.csv not found.\n")
  cat("  Run query_3dep_dates.py on your laptop first.\n")
  cat("  Looked in: ", file.path(out_tables, "3dep_date_summary.csv"), "\n")
  cat("  DTM temporal analysis will be SKIPPED.\n\n")
}


# =====================================================================
# P3D TILE DATES (shared by both CHM and DTM)
# =====================================================================

p3d_csv <- file.path(out_tables, "p3d_tile_dates_by_site.csv")
if (!file.exists(p3d_csv)) {
  stop("P3D tile date summary not found at:\n  ", p3d_csv,
       "\n\nRun extract_p3d_tile_dates.R first.")
}

p3d_dates <- read.csv(p3d_csv, stringsAsFactors = FALSE) %>%
  mutate(across(starts_with("p3d_date_"), as.Date))

cat(sprintf("P3D dates loaded for %d sites\n\n", nrow(p3d_dates)))


# =====================================================================
# LOAD MODELS AND DATA
# =====================================================================

cat("--- Loading models and data ---\n")

model_cp <- readRDS(file.path(CHECKPOINT_DIR, "10_models_stage2.rds"))
fit_chm  <- model_cp$data$fit_chm_s2
fit_dtm  <- model_cp$data$fit_dtm_s2

data_cp  <- readRDS(file.path(CHECKPOINT_DIR, "01_data_ingest.rds"))
chm_df   <- data_cp$data$chm_df
dtm_df   <- data_cp$data$dtm_df

FOREST_CLASSES <- c("BDF", "DNF", "EBF", "ENF")

cat("  Models and data loaded.\n\n")


# =====================================================================
# HELPER FUNCTIONS
# =====================================================================

extract_site_re <- function(fit) {
  re <- ranef(fit)$site
  data.frame(
    site         = as.integer(rownames(re[, , "Intercept"])),
    re_intercept = re[, "Estimate", "Intercept"],
    re_lower     = re[, "Q2.5", "Intercept"],
    re_upper     = re[, "Q97.5", "Intercept"],
    stringsAsFactors = FALSE
  )
}

compute_site_errors <- function(df, error_col) {
  df %>%
    filter(lc_l1_code %in% FOREST_CLASSES, is.finite(.data[[error_col]])) %>%
    group_by(site) %>%
    summarise(
      n_footprints   = n(),
      mean_error     = mean(.data[[error_col]], na.rm = TRUE),
      mean_abs_error = mean(abs(.data[[error_col]]), na.rm = TRUE),
      rmse           = sqrt(mean(.data[[error_col]]^2, na.rm = TRUE)),
      sd_error       = sd(.data[[error_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(site = as.integer(as.character(site)))
}

run_correlations <- function(df) {
  list(
    re    = cor.test(df$abs_gap_years, df$re_intercept),
    mae   = cor.test(df$abs_gap_years, df$mean_abs_error),
    bias  = cor.test(df$abs_gap_years, df$mean_error),
    signed = cor.test(df$gap_years, df$mean_error),
    rmse  = cor.test(df$abs_gap_years, df$rmse)
  )
}

print_correlations <- function(cors, product, n) {
  cat(sprintf("\n  %s RESULTS (n = %d sites):\n\n", product, n))
  cat(sprintf("  |gap| vs %s random intercept:    r = %+.3f  (p = %.3f)\n",
              product, cors$re$estimate, cors$re$p.value))
  cat(sprintf("  |gap| vs site mean |%s error|:   r = %+.3f  (p = %.3f)\n",
              product, cors$mae$estimate, cors$mae$p.value))
  cat(sprintf("  |gap| vs site mean %s bias:      r = %+.3f  (p = %.3f)\n",
              product, cors$bias$estimate, cors$bias$p.value))
  cat(sprintf("  Signed gap vs %s bias:           r = %+.3f  (p = %.3f)  *** directional ***\n",
              product, cors$signed$estimate, cors$signed$p.value))
  cat(sprintf("  |gap| vs site %s RMSE:           r = %+.3f  (p = %.3f)\n",
              product, cors$rmse$estimate, cors$rmse$p.value))
}

make_figure <- function(df, cors, product, ref_label) {

  theme_rse <- function(base_size = 9) {
    theme_minimal(base_size = base_size) %+replace%
      theme(
        panel.background  = element_rect(fill = "white", color = NA),
        plot.background   = element_rect(fill = "white", color = NA),
        panel.grid.major  = element_line(color = "gray92", linewidth = 0.3),
        panel.grid.minor  = element_blank(),
        axis.line         = element_line(color = "black", linewidth = 0.4),
        axis.ticks        = element_line(color = "black", linewidth = 0.3),
        axis.text         = element_text(color = "black", size = base_size - 1),
        axis.title        = element_text(color = "black", size = base_size),
        plot.tag          = element_text(size = 10, face = "bold"),
        plot.margin       = margin(6, 10, 6, 6)
      )
  }

  nudge <- 0.2

  cor_label <- function(ct) {
    sprintf("italic(r) == %+.2f ~~ (italic(p) == %.3f)",
            ct$estimate, ct$p.value)
  }

  # (a) |gap| vs random intercept
  p_a <- ggplot(df, aes(x = abs_gap_years, y = re_intercept)) +
    geom_errorbar(aes(ymin = re_lower, ymax = re_upper),
                  width = 0, linewidth = 0.3, color = "gray60") +
    geom_point(size = 2.5, color = "#0072B2") +
    geom_text(aes(label = paper_site),
              nudge_x = nudge, size = 2.5, color = "gray30") +
    geom_smooth(method = "lm", se = TRUE, color = "#D55E00",
                linewidth = 0.6, alpha = 0.1, linetype = "dashed") +
    annotate("text", x = Inf, y = Inf,
             label = cor_label(cors$re), parse = TRUE,
             hjust = 1.1, vjust = 1.3, size = 2.8) +
    labs(x = sprintf("|Temporal gap| (years)\n[P3D vs %s]", ref_label),
         y = "Site random intercept (m)",
         tag = "(a)") +
    theme_rse()

  # (b) |gap| vs MAE
  p_b <- ggplot(df, aes(x = abs_gap_years, y = mean_abs_error)) +
    geom_point(size = 2.5, color = "#0072B2") +
    geom_text(aes(label = paper_site),
              nudge_x = nudge, size = 2.5, color = "gray30") +
    geom_smooth(method = "lm", se = TRUE, color = "#D55E00",
                linewidth = 0.6, alpha = 0.1, linetype = "dashed") +
    annotate("text", x = Inf, y = Inf,
             label = cor_label(cors$mae), parse = TRUE,
             hjust = 1.1, vjust = 1.3, size = 2.8) +
    labs(x = sprintf("|Temporal gap| (years)\n[P3D vs %s]", ref_label),
         y = sprintf("Site mean |%s error| (m)", product),
         tag = "(b)") +
    theme_rse()

  # (c) Signed gap vs signed bias
  p_c <- ggplot(df, aes(x = gap_years, y = mean_error)) +
    geom_point(size = 2.5, color = "#0072B2") +
    geom_text(aes(label = paper_site),
              nudge_x = nudge, size = 2.5, color = "gray30") +
    geom_smooth(method = "lm", se = TRUE, color = "#D55E00",
                linewidth = 0.6, alpha = 0.1, linetype = "dashed") +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.3, color = "gray50") +
    annotate("text", x = Inf, y = Inf,
             label = cor_label(cors$signed), parse = TRUE,
             hjust = 1.1, vjust = 1.3, size = 2.8) +
    labs(x = sprintf("Signed gap (years)\n(+ve = P3D after %s)", ref_label),
         y = sprintf("Site mean %s bias (m)", product),
         tag = "(c)") +
    theme_rse()

  p_a + p_b + p_c + plot_layout(nrow = 1) +
    plot_annotation(theme = theme(
      plot.background = element_rect(fill = "white", color = NA)
    ))
}

print_gap_summary <- function(temporal, product, ref_label) {
  cat(sprintf("\n  %s temporal gaps (P3D vs %s):\n", product, ref_label))
  cat(sprintf("  Median |gap|:  %.1f years\n", median(temporal$abs_gap_years)))
  cat(sprintf("  IQR:           %.1f - %.1f years\n",
              quantile(temporal$abs_gap_years, 0.25),
              quantile(temporal$abs_gap_years, 0.75)))
  cat(sprintf("  Range:         %.1f - %.1f years\n",
              min(temporal$abs_gap_years), max(temporal$abs_gap_years)))
  cat(sprintf("  Sign: %d P3D before ref, %d after, %d ~same year\n",
              sum(temporal$gap_years < -0.5),
              sum(temporal$gap_years > 0.5),
              sum(abs(temporal$gap_years) <= 0.5)))
  cat(sprintf("  Median within-site P3D span: %.1f years\n",
              median(temporal$p3d_span_years)))
}


# =====================================================================
# =====================================================================
#  CHM ANALYSIS: P3D vs ALS
# =====================================================================
# =====================================================================

cat(strrep("=", 65), "\n")
cat("  CHM: P3D vs ALS temporal mismatch\n")
cat(strrep("=", 65), "\n")

# Compute gaps
chm_temporal <- als_dates %>%
  filter(site %in% SITES_CHM) %>%
  inner_join(p3d_dates %>% select(site, n_tiles, p3d_date_median,
                                    p3d_date_mean, p3d_date_min,
                                    p3d_date_max, p3d_span_years),
             by = "site") %>%
  mutate(
    ref_date      = als_date,
    gap_days      = as.numeric(difftime(p3d_date_median, als_date, units = "days")),
    gap_years     = gap_days / 365.25,
    abs_gap_years = abs(gap_years),
    paper_site    = paper_site_label(site)
  )

cat(sprintf("\nMatched %d CHM sites\n", nrow(chm_temporal)))
print_gap_summary(chm_temporal, "CHM", "ALS")

# Site random intercepts and error stats
chm_re     <- extract_site_re(fit_chm)
chm_errors <- compute_site_errors(chm_df, "chm_error_mean")

chm_analysis <- chm_temporal %>%
  inner_join(chm_re, by = "site") %>%
  inner_join(chm_errors, by = "site")

cat(sprintf("\n  Complete CHM data for %d sites\n", nrow(chm_analysis)))

# Correlations
chm_cors <- run_correlations(chm_analysis)
print_correlations(chm_cors, "CHM", nrow(chm_analysis))

# Figure
p_chm <- make_figure(chm_analysis, chm_cors, "CHM", "ALS")

chm_fig_path <- file.path(out_plots, "pub_figSXX_temporal_gap_chm.pdf")
ggsave(chm_fig_path, p_chm,
       width = 190, height = 70, units = "mm",
       bg = "white", device = safe_pdf)
ggsave(sub("\\.pdf$", ".png", chm_fig_path), p_chm,
       width = 190, height = 70, units = "mm", dpi = 300, bg = "white")
cat(sprintf("\n  CHM figure saved: %s\n", chm_fig_path))

# Save
write.csv(chm_analysis, file.path(out_tables, "temporal_gap_chm.csv"),
          row.names = FALSE)


# =====================================================================
# =====================================================================
#  DTM ANALYSIS: P3D vs 3DEP
# =====================================================================
# =====================================================================

if (has_3dep_dates) {

  cat("\n", strrep("=", 65), "\n", sep = "")
  cat("  DTM: P3D vs 3DEP temporal mismatch\n")
  cat(strrep("=", 65), "\n")

  # Sites for DTM (exclude site 10 and site 16)
  dtm_sites <- setdiff(SITES_CHM, SITES_EXCLUDE_DTM)

  dtm_temporal <- dep_dates %>%
    filter(site %in% dtm_sites) %>%
    select(site, dep_date_mid, dep_date_early, dep_date_late, dep_date_source) %>%
    inner_join(als_dates %>% select(site, site_name), by = "site") %>%
    inner_join(p3d_dates %>% select(site, n_tiles, p3d_date_median,
                                      p3d_date_mean, p3d_date_min,
                                      p3d_date_max, p3d_span_years),
               by = "site") %>%
    mutate(
      ref_date      = dep_date_mid,
      gap_days      = as.numeric(difftime(p3d_date_median, dep_date_mid, units = "days")),
      gap_years     = gap_days / 365.25,
      abs_gap_years = abs(gap_years),
      paper_site    = paper_site_label(site)
    )

  cat(sprintf("\nMatched %d DTM sites\n", nrow(dtm_temporal)))
  print_gap_summary(dtm_temporal, "DTM", "3DEP")

  # Site random intercepts and error stats
  dtm_re     <- extract_site_re(fit_dtm)
  dtm_errors <- compute_site_errors(dtm_df, "dtm_error_mean")

  dtm_analysis <- dtm_temporal %>%
    inner_join(dtm_re, by = "site") %>%
    inner_join(dtm_errors, by = "site")

  cat(sprintf("\n  Complete DTM data for %d sites\n", nrow(dtm_analysis)))

  # Correlations
  dtm_cors <- run_correlations(dtm_analysis)
  print_correlations(dtm_cors, "DTM", nrow(dtm_analysis))

  # Figure
  p_dtm <- make_figure(dtm_analysis, dtm_cors, "DTM", "3DEP")

  dtm_fig_path <- file.path(out_plots, "pub_figSXX_temporal_gap_dtm.pdf")
  ggsave(dtm_fig_path, p_dtm,
         width = 190, height = 70, units = "mm",
         bg = "white", device = safe_pdf)
  ggsave(sub("\\.pdf$", ".png", dtm_fig_path), p_dtm,
         width = 190, height = 70, units = "mm", dpi = 300, bg = "white")
  cat(sprintf("\n  DTM figure saved: %s\n", dtm_fig_path))

  # Save
  write.csv(dtm_analysis, file.path(out_tables, "temporal_gap_dtm.csv"),
            row.names = FALSE)

} else {
  cat("\n  DTM temporal analysis SKIPPED (no 3DEP dates).\n")
  cat("  Run query_3dep_dates.py, place output in: ",
      file.path(out_tables, "3dep_date_summary.csv"), "\n")
}


# =====================================================================
# =====================================================================
#  MANUSCRIPT-READY VALUES
# =====================================================================
# =====================================================================

cat("\n\n")
cat(strrep("=", 65), "\n")
cat("  MANUSCRIPT-READY VALUES\n")
cat(strrep("=", 65), "\n")

# --- CHM values ---
chm_med   <- median(chm_temporal$abs_gap_years)
chm_iqr   <- quantile(chm_temporal$abs_gap_years, c(0.25, 0.75))
chm_range <- range(chm_temporal$abs_gap_years)
chm_max   <- max(chm_temporal$abs_gap_years)

cat("\n--- Section 2.2 (P22) ---\n\n")
cat(sprintf(paste0(
  "The temporal gap between P3D source imagery and ALS reference\n",
  "acquisitions was computed for each site as the absolute difference\n",
  "between the median P3D tile acquisition date and the site-level\n",
  "ALS collection date. Across the %d CHM study sites, the median\n",
  "temporal gap was %.1f years (IQR: %.1f-%.1f years; range:\n",
  "%.1f-%.1f years)."),
  nrow(chm_temporal), chm_med, chm_iqr[1], chm_iqr[2],
  chm_range[1], chm_range[2]))

if (has_3dep_dates) {
  dtm_med   <- median(dtm_temporal$abs_gap_years)
  dtm_iqr   <- quantile(dtm_temporal$abs_gap_years, c(0.25, 0.75))
  dtm_range <- range(dtm_temporal$abs_gap_years)
  dtm_max   <- max(dtm_temporal$abs_gap_years)

  cat(sprintf(paste0(
    " For DTM, the gap between\n",
    "P3D imagery and 3DEP reference had a median of %.1f years\n",
    "(IQR: %.1f-%.1f; range: %.1f-%.1f years).\n"),
    dtm_med, dtm_iqr[1], dtm_iqr[2], dtm_range[1], dtm_range[2]))
} else {
  cat("\n[DTM gap values: run query_3dep_dates.py first]\n")
}


cat("\n\n--- Section 4.1 (P194): CHM temporal gap correlations ---\n\n")
cat(sprintf(paste0(
  "The temporal gap between P3D source imagery and ALS reference\n",
  "data varied from %.1f to %.1f years across sites (median %.1f\n",
  "years). [...] the Pearson correlation between site-level temporal\n",
  "gap and CHM random intercepts was r = %+.2f (p = %.3f), and the\n",
  "correlation with site-level mean absolute CHM error was r = %+.2f\n",
  "(p = %.3f) (Supplementary Figure SXX).\n"),
  chm_range[1], chm_range[2], chm_med,
  chm_cors$re$estimate, chm_cors$re$p.value,
  chm_cors$mae$estimate, chm_cors$mae$p.value))

cat(sprintf(
  "Directional test (signed gap vs CHM bias): r = %+.2f (p = %.3f)\n",
  chm_cors$signed$estimate, chm_cors$signed$p.value))

# Which continuation?
chm_weak <- abs(chm_cors$re$estimate) < 0.4 && abs(chm_cors$mae$estimate) < 0.4
cat(sprintf("-> %s: Use '%s correlation' continuation.\n",
            ifelse(chm_weak, "WEAK", "MODERATE/STRONG"),
            ifelse(chm_weak, "weak", "moderate")))


if (has_3dep_dates) {
  cat("\n\n--- Section 4.1 (P194): DTM temporal gap correlations ---\n\n")
  cat(sprintf(paste0(
    "For DTM, the temporal gap between P3D imagery and 3DEP reference\n",
    "ranged from %.1f to %.1f years (median %.1f years). The correlation\n",
    "between site-level temporal gap and DTM random intercepts was\n",
    "r = %+.2f (p = %.3f), and with site-level mean absolute DTM error\n",
    "was r = %+.2f (p = %.3f).\n"),
    dtm_range[1], dtm_range[2], dtm_med,
    dtm_cors$re$estimate, dtm_cors$re$p.value,
    dtm_cors$mae$estimate, dtm_cors$mae$p.value))

  cat(sprintf(
    "Directional test (signed gap vs DTM bias): r = %+.2f (p = %.3f)\n",
    dtm_cors$signed$estimate, dtm_cors$signed$p.value))

  dtm_weak <- abs(dtm_cors$re$estimate) < 0.4 && abs(dtm_cors$mae$estimate) < 0.4
  cat(sprintf("-> %s\n", ifelse(dtm_weak,
    "WEAK: Temporal mismatch not a major factor for DTM (terrain is stable).",
    "MODERATE/STRONG: Unexpected for DTM — investigate further.")))
}


cat("\n\n--- Section 4.8 (P231): Bounding argument ---\n\n")
cat("CHM bounding (growth rate 0.2-0.5 m/yr, temperate mature forest):\n")
cat(sprintf("  Max gap:    %.1f yr -> max contamination: %.1f m (at 0.5 m/yr)\n",
            chm_max, 0.5 * chm_max))
cat(sprintf("  Median gap: %.1f yr -> median contamination: %.1f m (at 0.5 m/yr)\n",
            chm_med, 0.5 * chm_med))
cat(sprintf("  As %% of CHM RMSE (5.61 m): max %.0f%%, median %.0f%%\n",
            100 * (0.5 * chm_max) / 5.61,
            100 * (0.5 * chm_med) / 5.61))
cat(sprintf("  As %% of 95%% PI width (~15.1 m): max %.0f%%, median %.0f%%\n",
            100 * (0.5 * chm_max) / 15.1,
            100 * (0.5 * chm_med) / 15.1))

if (has_3dep_dates) {
  cat(sprintf(paste0(
    "\nDTM bounding:\n",
    "  Terrain surfaces are geomorphically stable at these sites.\n",
    "  Max |gap|: %.1f years. No height growth signal applies to DTM.\n",
    "  Any 3DEP updates within this window reflect improved processing,\n",
    "  not terrain change.\n"), dtm_max))
}


# --- Supplementary figure captions ---
cat("\n\n--- Supplementary Figure Captions ---\n\n")

cat(sprintf(paste0(
  "Supplementary Figure SXX (CHM). Relationship between temporal gap\n",
  "(|P3D median acquisition date - ALS collection date|) and site-level\n",
  "CHM error characteristics. (a) Absolute temporal gap versus estimated\n",
  "site random intercept (95%% credible intervals). (b) Absolute temporal\n",
  "gap versus site-level mean absolute CHM error. (c) Signed temporal gap\n",
  "(positive = P3D after ALS) versus site-level mean CHM bias. n = %d sites.\n"),
  nrow(chm_analysis)))

if (has_3dep_dates) {
  cat(sprintf(paste0(
    "\nSupplementary Figure SXX (DTM). As above but for DTM, with temporal\n",
    "gap computed as |P3D median acquisition date - 3DEP reference midpoint\n",
    "date|. n = %d sites (Site 10 excluded due to datum offset).\n"),
    nrow(dtm_analysis)))
}


# =====================================================================
# SUMMARY TABLE
# =====================================================================

cat("\n\n--- Compact summary for quick reference ---\n\n")
cat(sprintf("%-6s %-10s %8s %8s %10s %10s\n",
            "Prod.", "Reference", "Med|gap|", "Max|gap|", "r(RE)", "r(MAE)"))
cat(strrep("-", 55), "\n")
cat(sprintf("%-6s %-10s %7.1f yr %7.1f yr %9.3f %9.3f\n",
            "CHM", "ALS", chm_med, chm_max,
            chm_cors$re$estimate, chm_cors$mae$estimate))
if (has_3dep_dates) {
  cat(sprintf("%-6s %-10s %7.1f yr %7.1f yr %9.3f %9.3f\n",
              "DTM", "3DEP", dtm_med, dtm_max,
              dtm_cors$re$estimate, dtm_cors$mae$estimate))
}


cat("\n", strrep("=", 65), "\n", sep = "")
cat("  TEMPORAL MISMATCH ANALYSIS COMPLETE\n")
cat(strrep("=", 65), "\n")
