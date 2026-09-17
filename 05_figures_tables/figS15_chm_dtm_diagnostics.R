# =====================================================================
# figS15_chm_dtm_diagnostics.R
#
# Generator for Supplementary Figure S15: error distribution diagnostics
# for CHM (P3D - ALS) and DTM (P3D - 3DEP) products.
#
# 8-panel composite (2 rows x 4 columns):
#   Row 1 - CHM, 16-site forest subset (drops flagged sites 1, 2, 3)
#     (a) Histogram + KDE + zero-line + mean-line + summary stats inset
#     (b) Q-Q plot vs theoretical normal with 95% pointwise envelope
#     (c) Hexbin error vs terrain slope, with LOESS conditional mean
#     (d) Hexbin error vs z-scaled WSCI, with LOESS conditional mean
#   Row 2 - DTM, 18-site forest subset (Site 10 already excluded upstream)
#     (e) Histogram + ...   (f) Q-Q ...   (g) Hexbin slope ...   (h) Hexbin WSCI ...
#
# Inputs:
#   checkpoints/01_data_ingest.rds  (chm_df, dtm_df with raw + z-scored predictors)
#
# Output:
#   plots/section_L_supp_v05/figS15.pdf  (PDF only; cluster lacks cairo/X11)
#
# Console:
#   Per-panel descriptive stats for caption sign-off.
#
# K-14 Item F (v2.7): replaces the ad-hoc figS15 embedded in v05c
# (which used the 19-site CHM forest subset and had no committed generator).
#
# Wall-time estimate: 1-3 minutes (mostly hexbin LOESS smoothing).
# =====================================================================

# Headless graphics setup
Sys.setenv(DISPLAY = "")
options(device = pdf)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(viridis)
  library(scales)
  library(moments)
})

# Optional package check
HAS_QQPLOTR <- requireNamespace("qqplotr", quietly = TRUE)
if (!HAS_QQPLOTR) {
  message("Note: qqplotr not available; Q-Q panels will use stat_qq/stat_qq_line ",
          "without confidence envelope.")
}

# ---------------------------------------------------------------------
# CONFIG
# ---------------------------------------------------------------------
if (!nzchar(Sys.getenv("PROJECT_ROOT"))) {
  Sys.setenv(PROJECT_ROOT = "/gpfs/data1/vclgp/lmaden/chpt1")
}
PROJECT_ROOT   <- Sys.getenv("PROJECT_ROOT")
checkpoint_dir <- file.path(PROJECT_ROOT, "checkpoints")
out_dir        <- file.path(PROJECT_ROOT, "plots", "section_L_supp_v05")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

FOREST_CLASSES    <- c("EBF", "BDF", "ENF", "DNF")
FLAGGED_SITES_CHM <- c("1", "2", "3")   # K-14 Item A: reference-quality screening
X_CLAMP_M         <- 30                  # histogram / hexbin y-axis clamp for visualization

# Colors (match supp_figs_plotting.R for visual consistency across supp)
COL_CHM    <- "#2D5F3F"
COL_DTM    <- "#C7522A"
COL_FIT    <- "#FF8800"   # KDE / mean / LOESS overlays
COL_ZERO   <- "black"     # zero-line reference (caption specifies black)
COL_QQ_FILL <- "#4477AA"  # Q-Q confidence envelope blue (matches caption)

# ---------------------------------------------------------------------
# 1. LOAD DATA AND APPLY FILTERS
# ---------------------------------------------------------------------
cat("\n[1/5] Loading data and applying filters...\n")
data <- readRDS(file.path(checkpoint_dir, "01_data_ingest.rds"))

chm_16 <- data$data$chm_df %>%
  filter(is.finite(chm_error_mean), !is.na(lc_l1_code),
         lc_l1_code %in% FOREST_CLASSES,
         !(site %in% FLAGGED_SITES_CHM))

dtm_18 <- data$data$dtm_df %>%
  filter(is.finite(dtm_error_mean), !is.na(lc_l1_code),
         lc_l1_code %in% FOREST_CLASSES)

cat(sprintf("  CHM 16-site forest subset: %s rows\n", format(nrow(chm_16), big.mark=",")))
cat(sprintf("  DTM 18-site forest subset: %s rows\n", format(nrow(dtm_18), big.mark=",")))

# Detect slope and WSCI column names (raw preferred for slope; z-scored for WSCI per caption)
detect_col <- function(df, candidates) {
  for (nm in candidates) if (nm %in% names(df)) return(nm)
  stop(sprintf("None of %s found in data.frame columns: %s",
               paste(candidates, collapse=","),
               paste(head(names(df), 20), collapse=",")))
}
slope_col <- detect_col(chm_16, c("slope_mean", "slope", "slope_mean_z"))
wsci_col  <- detect_col(chm_16, c("wsci_z", "wsci"))
slope_label <- if (slope_col == "slope_mean_z") "Terrain slope (z-scaled)" else "Terrain slope (degrees)"
wsci_label  <- if (wsci_col  == "wsci")         "WSCI"                     else "WSCI (z-scaled)"
cat(sprintf("  Slope column: %s  | label: %s\n", slope_col, slope_label))
cat(sprintf("  WSCI column : %s  | label: %s\n", wsci_col,  wsci_label))

# ---------------------------------------------------------------------
# 2. COMPUTE SUMMARY STATISTICS
# ---------------------------------------------------------------------
cat("\n[2/5] Computing summary statistics...\n")
compute_stats <- function(x) {
  list(n = length(x), mean = mean(x), sd = sd(x), median = median(x),
       skew = moments::skewness(x), kurt = moments::kurtosis(x))
}
s_chm <- compute_stats(chm_16$chm_error_mean)
s_dtm <- compute_stats(dtm_18$dtm_error_mean)

print_stats <- function(s, label) {
  cat(sprintf("  %s:\n", label))
  cat(sprintf("    n=%s, bias=%+.4f m, sd=%.4f m, median=%+.4f m, skew=%+.4f, kurtosis=%.4f\n",
              format(s$n, big.mark=","), s$mean, s$sd, s$median, s$skew, s$kurt))
}
print_stats(s_chm, "CHM 16-site")
print_stats(s_dtm, "DTM 18-site")

# ---------------------------------------------------------------------
# 3. PANEL BUILDERS
# ---------------------------------------------------------------------
theme_supp <- function() {
  theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
          plot.title       = element_text(face = "bold", size = 10),
          legend.position  = "right",
          legend.key.height = unit(0.4, "cm"),
          legend.key.width  = unit(0.3, "cm"),
          legend.title      = element_text(size = 7),
          legend.text       = element_text(size = 7))
}

# Histogram + KDE + stats inset
make_hist_panel <- function(df, err_col, color_main, s, product_label, panel_label) {
  df_p <- df[abs(df[[err_col]]) <= X_CLAMP_M, , drop = FALSE]
  inset <- sprintf("n = %s\nbias = %+.2f m\nskew = %+.2f\nkurtosis = %.2f",
                   format(s$n, big.mark=","), s$mean, s$skew, s$kurt)
  ggplot(df_p, aes(x = .data[[err_col]])) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = 80, fill = color_main, alpha = 0.55, color = NA) +
    geom_density(color = COL_FIT, linewidth = 0.7) +
    geom_vline(xintercept = 0,      linetype = "dashed", color = COL_ZERO, linewidth = 0.4) +
    geom_vline(xintercept = s$mean, linetype = "dotted", color = COL_FIT,  linewidth = 0.55) +
    annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.2, size = 2.6,
             family = "mono", label = inset) +
    coord_cartesian(xlim = c(-X_CLAMP_M, X_CLAMP_M)) +
    labs(x = sprintf("%s error (m)", product_label), y = "Density",
         title = sprintf("(%s) %s error distribution", panel_label, product_label)) +
    theme_supp()
}

# Q-Q plot
make_qq_panel <- function(df, err_col, color_main, product_label, panel_label,
                          sample_max = 50000) {
  x <- df[[err_col]]
  if (length(x) > sample_max) x <- sample(x, sample_max)
  df_s <- data.frame(err = x)
  p <- ggplot(df_s, aes(sample = err))
  if (HAS_QQPLOTR) {
    p <- p +
      qqplotr::stat_qq_band(distribution = "norm",
                            fill = COL_QQ_FILL, alpha = 0.30) +
      qqplotr::stat_qq_line(distribution = "norm",
                            color = color_main, linewidth = 0.5) +
      qqplotr::stat_qq_point(distribution = "norm",
                             color = color_main, alpha = 0.45, size = 0.4)
  } else {
    p <- p +
      stat_qq(color = color_main, alpha = 0.45, size = 0.4) +
      stat_qq_line(color = color_main, linewidth = 0.5)
  }
  p +
    labs(x = "Theoretical normal quantile",
         y = sprintf("%s sample quantile (m)", product_label),
         title = sprintf("(%s) %s Q-Q plot", panel_label, product_label)) +
    theme_supp()
}

# Hexbin error vs covariate
# NOTE: LOESS via geom_smooth has O(N^2) memory cost; full forest subsets
# (~250k-460k rows) exceed predLoess workspace. We pass a random subsample
# to geom_smooth while keeping full data for geom_hex density.
SMOOTH_N <- 20000

make_hex_panel <- function(df, err_col, cov_col, color_main,
                           product_label, panel_label, x_axis_label) {
  cols_keep <- c(err_col, cov_col)
  df_p <- df[is.finite(df[[cov_col]]), cols_keep, drop = FALSE]
  names(df_p) <- c("err", "x")
  df_p <- df_p[abs(df_p$err) <= X_CLAMP_M, , drop = FALSE]

  set.seed(2025)  # reproducible LOESS subsample
  df_smooth <- if (nrow(df_p) > SMOOTH_N) {
    df_p[sample(nrow(df_p), SMOOTH_N), , drop = FALSE]
  } else {
    df_p
  }

  ggplot(df_p, aes(x = x, y = err)) +
    geom_hex(bins = 50) +
    scale_fill_viridis(option = "plasma", trans = "log10",
                       name = "count",
                       labels = scales::trans_format("identity", scales::comma)) +
    geom_smooth(data = df_smooth, method = "loess", span = 0.3, se = TRUE,
                color = COL_FIT, fill = COL_FIT, alpha = 0.25, linewidth = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = COL_ZERO, linewidth = 0.4) +
    labs(x = x_axis_label,
         y = sprintf("%s error (m)", product_label),
         title = sprintf("(%s) %s vs %s",
                         panel_label, product_label,
                         tolower(sub("\\s*\\([^)]*\\)", "", x_axis_label)))) +
    theme_supp()
}

# ---------------------------------------------------------------------
# 4. BUILD ALL 8 PANELS
# ---------------------------------------------------------------------
cat("\n[3/5] Building panels...\n")

cat("  CHM (16-site)...\n")
p_a <- make_hist_panel(chm_16, "chm_error_mean", COL_CHM, s_chm, "CHM", "a")
p_b <- make_qq_panel  (chm_16, "chm_error_mean", COL_CHM,         "CHM", "b")
p_c <- make_hex_panel (chm_16, "chm_error_mean", slope_col, COL_CHM,
                       "CHM", "c", slope_label)
p_d <- make_hex_panel (chm_16, "chm_error_mean", wsci_col,  COL_CHM,
                       "CHM", "d", wsci_label)

cat("  DTM (18-site)...\n")
p_e <- make_hist_panel(dtm_18, "dtm_error_mean", COL_DTM, s_dtm, "DTM", "e")
p_f <- make_qq_panel  (dtm_18, "dtm_error_mean", COL_DTM,         "DTM", "f")
p_g <- make_hex_panel (dtm_18, "dtm_error_mean", slope_col, COL_DTM,
                       "DTM", "g", slope_label)
p_h <- make_hex_panel (dtm_18, "dtm_error_mean", wsci_col,  COL_DTM,
                       "DTM", "h", wsci_label)

# ---------------------------------------------------------------------
# 5. COMPOSE AND SAVE
# ---------------------------------------------------------------------
cat("\n[4/5] Composing 2 x 4 layout...\n")
fig <- (p_a | p_b | p_c | p_d) /
       (p_e | p_f | p_g | p_h)

cat("\n[5/5] Saving output...\n")
pdf_path <- file.path(out_dir, "figS15.pdf")

## Default PDF device (cairo not available on this cluster; not required
## for this figure since no Unicode characters appear in labels/annotations).
## PNG is intentionally not generated here because this cluster has neither
## cairo nor an X11 display. If a PNG preview is needed, convert from PDF
## locally (e.g., pdftools::pdf_render_page, ImageMagick, or open in a viewer).
ggsave(pdf_path, fig, width = 14, height = 7, units = "in")

cat(sprintf("  PDF: %s  (%.1f KB)\n", pdf_path, file.info(pdf_path)$size / 1024))

# ---------------------------------------------------------------------
# CAPTION SIGN-OFF SUMMARY
# ---------------------------------------------------------------------
cat("\n",
    strrep("=", 70), "\n",
    "figS15 v2.7 regeneration complete -- caption sign-off summary\n",
    strrep("=", 70), "\n", sep = "")
cat(sprintf("CHM (panels a-d, 16-site forest subset):\n"))
cat(sprintf("  n=%s, bias=%+.2f m, sd=%.2f m, skew=%+.2f, kurtosis=%.2f\n",
            format(s_chm$n, big.mark=","), s_chm$mean, s_chm$sd,
            s_chm$skew, s_chm$kurt))
cat(sprintf("DTM (panels e-h, 18-site forest subset):\n"))
cat(sprintf("  n=%s, bias=%+.2f m, sd=%.2f m, skew=%+.2f, kurtosis=%.2f\n",
            format(s_dtm$n, big.mark=","), s_dtm$mean, s_dtm$sd,
            s_dtm$skew, s_dtm$kurt))
cat(strrep("=", 70), "\n", sep = "")
