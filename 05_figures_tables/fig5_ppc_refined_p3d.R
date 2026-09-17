# =====================================================================
# fig5_ppc_refined_p3d.R   [HOUSE-STYLE GENERATOR — Figure 5]
#
# House-style (palette.R + theme_p3d) re-render of manuscript FIGURE 5
# (PPC by forest functional type). Save as a NEW file in scripts/reviewed/;
# does NOT overwrite fig4_ppc_refined.R.
#
# PRESENTATION-ONLY change vs fig4_ppc_refined.R:
#   - colors : COLOR_CHM/COLOR_DTM (#1f78b4/#ff7f00, from fig_common)
#              -> p3d_col$chm (#0072B2) / p3d_col$dtm (#D55E00)  [Okabe-Ito, locked]
#   - theme  : theme_section_G(9) -> theme_p3d(11, grid = "xy"), y-title re-enabled
#   - saver  : save_figure -> save_p3d(dpi = 600)  [ragg, cairo-free]
#   - titles : shortened to the LC abbreviation + n (full names live in the
#              caption glossary); inner-column y-axes suppressed (Density + ticks
#              on the leftmost panel of each row only)
#   - tags   : bold (a)-(h) drawn as a data-space annotation in each panel's
#              empty top-left corner (density ~ 0 at x = -25) so placement is
#              identical across ggplot 3.4.x and 4.0 (no plot.tag.location dep)
#   - x ticks: ASCII hyphen via label_number(style_negative = "hyphen")
#
# DATA LAYER (Sections 1-3) and the Section-5 FINGERPRINT are COPIED VERBATIM
#   from fig4_ppc_refined.R — no sampling / seed / exclusion / anchor
#   logic changed, so per-panel n is identical and the fingerprint still guards
#   the swap.
#
# RUN from scripts/reviewed/:
#   cd scripts/reviewed && Rscript fig5_ppc_refined_p3d.R
# Output: plots/section_G/fig05_ppc_refined_p3d.{png,pdf}
#   rendered 11.0 x 5.6 in (aspect 1.964 == docx wp:extent for image5.png) at
#   600 dpi; Word scales it into the 6.5-in column.
# =====================================================================

source("fig_common.R")   # load_checkpoint, log_*, SECTION_G_PLOTS, fig_banner
source("palette.R")               # house colors: p3d_col, p3d_gray  (source AFTER common)
source("theme_p3d.R")             # house theme + saver: theme_p3d(), save_p3d(), p3d_type
fig_banner("Figure 5 (PPC refined, house style)",
                 "PPC ribbon + observed, shared axes; house palette + theme_p3d (CHM 16-site, DTM 18-site)")

suppressPackageStartupMessages({
  library(brms); library(posterior); library(ggplot2)
  library(data.table); library(dplyr); library(patchwork)
  library(cowplot); library(scales)
})

# =====================================================================
# SECTIONS 1-3 BELOW ARE COPIED VERBATIM FROM
# fig4_ppc_by_forest_type.R (verified data layer), with one
# benign addition: n_avail (pre-cap holdout pool size) recorded per
# stratum for the fingerprint. No sampling/seed/exclusion logic changed.
# =====================================================================

# ---- 1. Load fits and reconstruct canonical holdout ----
log_subsection("Loading fits")
chm_ck <- load_checkpoint("10b_chm_sensitivity")
fit_chm_16 <- chm_ck$fit_chm_16
s2 <- load_checkpoint("10_models_stage2")
fit_dtm <- s2$fit_dtm_s2

log_subsection("Reconstructing canonical holdout")
stage2_frac <- as.numeric(Sys.getenv("STAGE2_FRAC", "0.33"))
SITES_EXCLUDE_CHM_16 <- c("1", "2", "3")
SITES_EXCLUDE_DTM_18 <- c("10")

full_ck <- load_checkpoint("01_data_ingest")
if (!all(c("chm_df", "dtm_df") %in% names(full_ck)))
  stop("01_data_ingest checkpoint missing chm_df/dtm_df. Names: ",
       paste(names(full_ck), collapse = ", "))

chm_df_full <- as.data.frame(full_ck$chm_df) %>%
  filter(is.finite(chm_error_mean), !is.na(lc_l1_code), !is.na(site), !is.na(ecoregion)) %>%
  mutate(lc_l1_code = factor(lc_l1_code), ecoregion = factor(ecoregion), site = factor(site))
dtm_df_full <- as.data.frame(full_ck$dtm_df) %>%
  filter(is.finite(dtm_error_mean), !is.na(lc_l1_code), !is.na(site), !is.na(ecoregion)) %>%
  mutate(lc_l1_code = factor(lc_l1_code), ecoregion = factor(ecoregion), site = factor(site))
chm_df_full$original_row <- seq_len(nrow(chm_df_full))
dtm_df_full$original_row <- seq_len(nrow(dtm_df_full))

set.seed(2025)
chm_train <- chm_df_full %>% group_by(site) %>% sample_frac(stage2_frac) %>% ungroup()
dtm_train <- dtm_df_full %>% group_by(site) %>% sample_frac(stage2_frac) %>% ungroup()
train_lc_chm <- levels(droplevels(chm_train$lc_l1_code)); train_site_chm <- levels(droplevels(chm_train$site)); train_eco_chm <- levels(droplevels(chm_train$ecoregion))
train_lc_dtm <- levels(droplevels(dtm_train$lc_l1_code)); train_site_dtm <- levels(droplevels(dtm_train$site)); train_eco_dtm <- levels(droplevels(dtm_train$ecoregion))

chm_holdout <- chm_df_full %>%
  filter(!original_row %in% chm_train$original_row) %>% select(-original_row) %>%
  filter(lc_l1_code %in% train_lc_chm, site %in% train_site_chm, ecoregion %in% train_eco_chm) %>%
  mutate(lc_l1_code = factor(lc_l1_code, levels = train_lc_chm),
         site = factor(site, levels = train_site_chm),
         ecoregion = factor(ecoregion, levels = train_eco_chm))
dtm_holdout <- dtm_df_full %>%
  filter(!original_row %in% dtm_train$original_row) %>% select(-original_row) %>%
  filter(lc_l1_code %in% train_lc_dtm, site %in% train_site_dtm, ecoregion %in% train_eco_dtm) %>%
  mutate(lc_l1_code = factor(lc_l1_code, levels = train_lc_dtm),
         site = factor(site, levels = train_site_dtm),
         ecoregion = factor(ecoregion, levels = train_eco_dtm))

mod_chm_16 <- chm_holdout %>% filter(!site %in% SITES_EXCLUDE_CHM_16) %>%
  mutate(site = droplevels(site), lc_l1_code = droplevels(lc_l1_code), ecoregion = droplevels(ecoregion))
mod_dtm <- dtm_holdout %>% filter(!site %in% SITES_EXCLUDE_DTM_18) %>%
  mutate(site = droplevels(site), lc_l1_code = droplevels(lc_l1_code), ecoregion = droplevels(ecoregion))
if ("eco_on" %in% all.vars(formula(fit_chm_16))) mod_chm_16$eco_on <- 1
if ("eco_on" %in% all.vars(formula(fit_dtm)))    mod_dtm$eco_on    <- 1
log_progress(sprintf("  CHM 16-site holdout rows: %s", formatC(nrow(mod_chm_16), big.mark = ",", format = "d")))
log_progress(sprintf("  DTM 18-site holdout rows: %s", formatC(nrow(mod_dtm), big.mark = ",", format = "d")))

# ---- 2. Forest functional type strata ----
LC_LABELS <- c(BDF = "Broadleaf Deciduous", DNF = "Deciduous Needleleaf",
               EBF = "Evergreen Broadleaf", ENF = "Evergreen Needleleaf")
LC_ORDER <- c("BDF", "DNF", "EBF", "ENF")
LC_COL_CANDIDATES <- c("lc_l1_code", "lc2022_l1_code", "lc_l1", "lc")
pick_lc_col <- function(df, label = "df") {
  for (c in LC_COL_CANDIDATES) if (c %in% names(df)) {
    vals <- unique(as.character(df[[c]])); n_match <- sum(LC_ORDER %in% vals)
    if (n_match >= 2L) { log_progress(sprintf("  [%s] Picked LC column '%s' (%d/%d expected forest codes present)", label, c, n_match, length(LC_ORDER))); return(c) }
  }
  stop("Could not find a usable LC column in ", label, ". Available: ", paste(names(df), collapse = ", "))
}
lc_chm <- pick_lc_col(mod_chm_16, "CHM 16-site")
lc_dtm <- pick_lc_col(mod_dtm,    "DTM 18-site")

# ---- 3. Stratified sampling + PPC draws (n_avail added for fingerprint) ----
log_subsection("Stratified sampling and PPC draws")
set.seed(42L)
N_PER_STRATUM <- 5000L
N_PPC_DRAWS   <- 50L
panel_data_for <- function(df, lc_col, fit, product_label) {
  df_f <- df[df[[lc_col]] %in% LC_ORDER, ]; out_panels <- list()
  for (lc in LC_ORDER) {
    sub <- df_f[df_f[[lc_col]] == lc, ]; n_avail <- nrow(sub)
    if (n_avail == 0) { log_progress(sprintf("  %s %s: no rows; skipping", product_label, lc)); next }
    n_take <- min(N_PER_STRATUM, n_avail); idx <- sample.int(n_avail, n_take); sub2 <- sub[idx, , drop = FALSE]
    yrep <- posterior_predict(fit, newdata = sub2, ndraws = N_PPC_DRAWS, allow_new_levels = TRUE, re_formula = NULL)
    response_col <- if (product_label == "CHM") "chm_error_mean" else "dtm_error_mean"
    y <- as.numeric(sub2[[response_col]])
    ok <- (!apply(is.na(yrep), 2, any)) & (!is.na(y)); n_drop <- sum(!ok)
    if (n_drop > 0L) { log_progress(sprintf("  %s %s: dropping %d/%d NA rows", product_label, lc, n_drop, length(ok))); yrep <- yrep[, ok, drop = FALSE]; y <- y[ok] }
    if (length(y) < 50L) { log_progress(sprintf("  %s %s: <50 usable rows; skipping", product_label, lc)); next }
    out_panels[[lc]] <- list(y = y, yrep = yrep, n = length(y), n_avail = n_avail, lc = lc, product = product_label)
  }
  out_panels
}
chm_panels <- panel_data_for(mod_chm_16, lc_chm, fit_chm_16, "CHM")
dtm_panels <- panel_data_for(mod_dtm,    lc_dtm, fit_dtm,    "DTM")

# =====================================================================
# SECTION 4 (NEW): refined panels -- ribbon + observed line, shared axes.
# =====================================================================
log_subsection("Composing refined panels (ribbon + observed, shared axes)")
XR <- c(-25, 25); GX_N <- 512L

compute_panel <- function(pd) {
  y <- pd$y; yrep <- pd$yrep
  bw <- density(y)$bw                                   # common bandwidth for comparability
  do <- density(y, bw = bw, from = XR[1], to = XR[2], n = GX_N)
  dmat <- t(apply(yrep, 1, function(v) density(v, bw = bw, from = XR[1], to = XR[2], n = GX_N)$y))
  list(gx = do$x, obs = do$y,
       q05 = apply(dmat, 2, quantile, 0.05),
       q95 = apply(dmat, 2, quantile, 0.95),
       n = pd$n, lc = pd$lc, product = pd$product)
}

chm_cp <- lapply(LC_ORDER, function(lc) if (!is.null(chm_panels[[lc]])) compute_panel(chm_panels[[lc]]))
dtm_cp <- lapply(LC_ORDER, function(lc) if (!is.null(dtm_panels[[lc]])) compute_panel(dtm_panels[[lc]]))
row_ymax <- function(cps) max(unlist(lapply(cps, function(c) if (is.null(c)) 0 else max(c$q95, c$obs))))
chm_ylim <- c(0, row_ymax(chm_cp) * 1.08)
dtm_ylim <- c(0, row_ymax(dtm_cp) * 1.08)

# RESTYLED (presentation only): house palette + theme_p3d; abbreviation + n
# title; corner (a)-(h) tag as a data-space annotation in the empty top-left;
# inner-column y-axes suppressed via show_y (leftmost panel of each row = TRUE).
render_panel <- function(cp, label, ylim, color, show_y = TRUE) {
  if (is.null(cp)) return(patchwork::plot_spacer())
  df <- data.frame(x = cp$gx, lo = cp$q05, hi = cp$q95, obs = cp$obs)
  ggplot(df, aes(x = x)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = scales::alpha(color, 0.22)) +
    geom_line(aes(y = obs), color = color, linewidth = 0.7) +
    annotate("text", x = XR[1] + 1.2, y = ylim[2] * 0.94, label = label,
             hjust = 0, vjust = 1, fontface = "bold",
             size = p3d_type$panel_label / .pt, color = p3d_gray$text_hi) +
    coord_cartesian(xlim = XR, ylim = ylim, expand = FALSE) +
    scale_x_continuous(breaks = seq(-20, 20, 10),
                       labels = scales::label_number(style_negative = "hyphen")) +
    scale_y_continuous(labels = scales::label_number(accuracy = 0.01)) +
    labs(title = sprintf("%s  (n = %s)", cp$lc,
                         formatC(cp$n, big.mark = ",", format = "d")),
         x = sprintf("%s error (m)", cp$product), y = "Density") +
    theme_p3d(base_size = 11, grid = "xy") +
    theme(
      plot.title   = element_text(size = p3d_type$strip, face = "bold", hjust = 0.5,
                                  color = p3d_gray$strip, margin = margin(b = 5)),
      axis.title.x = element_text(size = p3d_type$axis_title, color = p3d_gray$text_hi),
      axis.title.y = if (show_y) element_text(size = p3d_type$axis_title,
                                              color = p3d_gray$text_hi, angle = 90) else element_blank(),
      axis.text.y  = if (show_y) element_text(size = p3d_type$tick, color = p3d_gray$text_lo) else element_blank(),
      axis.ticks.y = if (show_y) element_line(color = p3d_gray$grid, linewidth = 0.3) else element_blank(),
      legend.position = "none",
      plot.margin  = margin(9, 9, 4, 7)
    )
}

labels_chm <- c("(a)", "(b)", "(c)", "(d)"); labels_dtm <- c("(e)", "(f)", "(g)", "(h)")
chm_plots <- lapply(seq_along(LC_ORDER), function(i) render_panel(chm_cp[[i]], labels_chm[i], chm_ylim, p3d_col$chm, show_y = (i == 1L)))
dtm_plots <- lapply(seq_along(LC_ORDER), function(i) render_panel(dtm_cp[[i]], labels_dtm[i], dtm_ylim, p3d_col$dtm, show_y = (i == 1L)))

panel_layout <- (chm_plots[[1]] | chm_plots[[2]] | chm_plots[[3]] | chm_plots[[4]]) /
                (dtm_plots[[1]] | dtm_plots[[2]] | dtm_plots[[3]] | dtm_plots[[4]])
final_plot <- cowplot::ggdraw(panel_layout) +
  cowplot::draw_label("CHM", x = 0.007, y = 0.75, angle = 90, size = 11, fontface = "bold", colour = p3d_gray$strip) +
  cowplot::draw_label("DTM", x = 0.007, y = 0.27, angle = 90, size = 11, fontface = "bold", colour = p3d_gray$strip)

# =====================================================================
# SECTION 5: save + FINGERPRINT (per-panel n vs current Fig 4; §2C pool)
# =====================================================================
save_p3d(final_plot, "fig05_ppc_refined_p3d", width_in = 11.0, height_in = 5.6,
         dpi = 600, dir = SECTION_G_PLOTS)

log_subsection("FINGERPRINT vs current Figure 4 (image5.png) per-panel n")
# Recorded n from the current manuscript Fig 4 raster. Tolerance allows
# the small NA-drop jitter; a large miss means wrong fit/data -> STOP.
N_ANCHOR <- list(CHM = c(BDF = 5000, DNF = 5000, EBF = 1520, ENF = 4998),
                 DTM = c(BDF = 5000, DNF = 5000, EBF = 5000, ENF = 4999))
check_n <- function(panels, product, tol = 5) {
  for (lc in LC_ORDER) {
    pd <- panels[[lc]]; got <- if (is.null(pd)) NA_integer_ else pd$n; exp <- N_ANCHOR[[product]][[lc]]
    if (is.na(got)) stop(sprintf("ANCHOR MISMATCH: %s %s panel missing (expected n=%d)", product, lc, exp))
    if (abs(got - exp) > tol)
      stop(sprintf("ANCHOR MISMATCH: %s %s n=%d, expected ~%d (tol %d) -> WRONG fit/data; do NOT use.",
                   product, lc, got, exp, tol))
    cat(sprintf("  OK  %s %s  n=%d (anchor %d)  holdout pool=%s\n",
                product, lc, got, exp, formatC(pd$n_avail, big.mark = ",", format = "d")))
  }
}
check_n(chm_panels, "CHM"); check_n(dtm_panels, "DTM")
cat("\nStory-lock §2C 16-site forest-subset FULL counts: ENF 193,453; BDF 44,394; DNF 9,035; EBF 2,266.\n")
cat("Approx 67% holdout pool: ENF ~129,514; BDF ~29,744; DNF ~6,054; EBF ~1,518.\n")
cat("Compare the printed 'holdout pool' above against these for the §2C sanity check.\n")
cat("Fingerprint PASSED.\n")
log_progress("fig4_ppc_refined.R complete.")
