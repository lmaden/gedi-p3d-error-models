# =====================================================================
# suppfig_S16_15panel_hexbin.R  (24 Aug 2026: 15-site basis)
# Supp Figure S16: 15-panel array of per-footprint within-site
# CHM-DTM error scatter, one panel per DTM-OK site.
#
# 
#
# Layout: 6 cols x 3 rows = 18 panels. Each panel is a hexbin of
# err_chm vs err_dtm at one site (joined on shot_number), with the
# within-site r annotated. Panels ordered by manuscript ID.
#
# Sources:
#   - load_checkpoint("08_model_prep")$mod_chm_s2, mod_dtm_s2
#   - coupling_rh98_per_site_slopes.csv (per-site r for annotation)
#
# Notes:
# (P1) `tracker_to_manuscript()` from fig_common.R errored
#      with "Error in .(tracker_id) : could not find function \".\""
#      when called from inside lapply() in this script, even though
#      it returns 5 for scalar 5L at the R prompt. Workaround: use an
#      inline mapping per story-lock §0 (tracker 1-15 -> manuscript
#      1-15; tracker 17 -> 16; 18 -> 17; 19 -> 18; 20 -> 19). Flagged
#      for the helper audit.
# (P2) In-panel annotation switched from annotate("text", ...) to
#      annotate("label", ...) with opaque white fill and no border.
#      High-density panels (Sites 1, 5, 18) have hexbins reaching the
#      panel corners, which was occluding the plain-text annotation.
#      The label-with-white-fill version sits on top of the data
#      regardless of corner density and remains readable.
# =====================================================================

source("fig_common.R")
fig_banner("Supp Figure S16",
                 "15-panel within-site CHM-DTM error hexbin array")

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(dplyr); library(patchwork)
})

# ---------------------------------------------------------------------
# 1. Load data
# ---------------------------------------------------------------------
mp <- load_checkpoint("08_model_prep")
mod_chm <- as.data.table(mp$mod_chm_s2)
mod_dtm <- as.data.table(mp$mod_dtm_s2)

ps <- as.data.table(fread(file.path(manuscript_tables_dir,
                                    "coupling_rh98_per_site_slopes.csv")))
site_col_ps <- intersect(c("tracker_site", "tracker_site_id",
                           "tracker_id", "site", "site_id"),
                         names(ps))[1]
setnames(ps, site_col_ps, "tracker_id_str")
ps[, tracker_id := suppressWarnings(as.integer(as.character(tracker_id_str)))]
r_col <- intersect(c("corr_errP3D_errDTM", "corr_err", "r_within"),
                   names(ps))[1]

# The 15 sites in both analyses (from common).
target_trackers <- intersect(TRACKER_PHASE25_15, unique(c(
  as.integer(as.character(mod_chm$site)),
  as.integer(as.character(mod_dtm$site)))))
log_progress(sprintf("  Target sites (intersect with available data): %s",
                     paste(sort(target_trackers), collapse = ",")))

# Shot-number key.
join_key_candidates <- c("shot_number", "shot_id", "footprint_id",
                          "gedi_shot_number")
key_chm <- intersect(join_key_candidates, names(mod_chm))[1]
key_dtm <- intersect(join_key_candidates, names(mod_dtm))[1]

# ---------------------------------------------------------------------
# 1a. Inline tracker -> manuscript mapping (bypasses scalar-broken helper)
# ---------------------------------------------------------------------
tracker_to_manuscript_local <- function(tid) {
  tid <- suppressWarnings(as.integer(tid))
  if (length(tid) == 0L || is.na(tid)) return(NA_integer_)
  if (tid >= 1L && tid <= 15L)          return(tid)
  if (tid == 17L)                       return(16L)
  if (tid == 18L)                       return(17L)
  if (tid == 19L)                       return(18L)
  if (tid == 20L)                       return(19L)
  return(NA_integer_)
}

# ---------------------------------------------------------------------
# 2. Per-site joined error vectors
# ---------------------------------------------------------------------
build_site_panel <- function(tid) {
  m <- tracker_to_manuscript_local(tid)
  chm_s <- mod_chm[as.character(site) == as.character(tid),
                   .(shot = get(key_chm), err_chm = chm_error_mean)]
  dtm_s <- mod_dtm[as.character(site) == as.character(tid),
                   .(shot = get(key_dtm), err_dtm = dtm_error_mean)]
  j <- merge(chm_s, dtm_s, by = "shot", all = FALSE)
  if (!nrow(j)) return(NULL)
  r_csv <- ps[tracker_id == tid, get(r_col)][1]
  r_emp <- if (nrow(j) >= 3L) cor(j$err_chm, j$err_dtm) else NA_real_
  list(manuscript_id = m, tracker_id = tid,
       data = j, r_csv = r_csv, r_emp = r_emp, n = nrow(j))
}

panel_data_list <- lapply(target_trackers, build_site_panel)
panel_data_list <- Filter(Negate(is.null), panel_data_list)
# Order by manuscript ID.
panel_data_list <- panel_data_list[order(sapply(panel_data_list,
                                                function(x) x$manuscript_id))]

log_progress(sprintf("  Built %d site panels", length(panel_data_list)))

# ---------------------------------------------------------------------
# 3. Build a hexbin panel function
# ---------------------------------------------------------------------
mk_panel <- function(pd) {
  r_use <- if (is.finite(pd$r_csv)) pd$r_csv else pd$r_emp
  ggplot(pd$data, aes(x = err_dtm, y = err_chm)) +
    geom_hex(bins = 40, show.legend = FALSE) +
    geom_smooth(method = "lm", se = FALSE,
                color = COLOR_REF_LINE, linewidth = 0.4) +
    geom_hline(yintercept = 0, color = COLOR_ZERO_LINE,
               linewidth = 0.15, linetype = "dotted") +
    geom_vline(xintercept = 0, color = COLOR_ZERO_LINE,
               linewidth = 0.15, linetype = "dotted") +
    scale_fill_viridis_c(trans = "log10", option = "viridis") +
    annotate("label",
             x = -Inf, y = Inf,
             label = sprintf("Site %d\nr = %s%.2f\nn = %s",
                             pd$manuscript_id,
                             ifelse(is.finite(r_use) && r_use < 0, MINUS, ""),
                             abs(r_use),
                             formatC(pd$n, big.mark = ",", format = "d")),
             hjust = 0, vjust = 1,
             size = 2.5, lineheight = 1.0,
             color = "black",
             fill = "white",
             label.size = NA,
             label.padding = unit(0.18, "lines"),
             label.r = unit(0, "lines")) +
    coord_cartesian(xlim = c(-15, 15), ylim = c(-15, 15)) +
    labs(x = NULL, y = NULL) +
    theme_section_G(base_size = 8) +
    theme(
      axis.text.x      = element_text(size = 7),
      axis.text.y      = element_text(size = 7),
      panel.grid       = element_blank(),
      plot.margin      = margin(2, 2, 2, 2),
      panel.spacing    = unit(0.3, "lines")
    )
}

panels <- lapply(panel_data_list, mk_panel)

# ---------------------------------------------------------------------
# 4. Compose 6-col x 3-row grid (or use ceiling if !=18)
# ---------------------------------------------------------------------
n_pan <- length(panels)
ncol <- 5L
nrow <- ceiling(n_pan / ncol)

combined <- patchwork::wrap_plots(panels, ncol = ncol, nrow = nrow) +
  patchwork::plot_annotation(
    title = "Per-footprint within-site CHM vs DTM error scatter",
    subtitle = "The 15 sites in both analyses; manuscript site ID, within-site r, and n shown per panel",
    theme = theme(plot.title    = element_text(size = 11, face = "bold"),
                  plot.subtitle = element_text(size = 9,  color = "gray30"))
  ) &
  labs(x = "DTM error (m)", y = "CHM error (m)") &
  theme(plot.tag = element_blank())

save_figure(combined, "suppfig_S16_15panel_hexbin",
            width_in = 9.6, height_in = 6.5)

# ---------------------------------------------------------------------
# 5. Verification
# ---------------------------------------------------------------------
log_subsection("VERIFICATION (console)")
cat("\nPer-site r values:\n")
for (pd in panel_data_list) {
  cat(sprintf("  Site %2d (tracker %2d):  r_csv = %s%5.3f   r_emp = %s%5.3f   n = %s\n",
              pd$manuscript_id, pd$tracker_id,
              ifelse(is.finite(pd$r_csv) && pd$r_csv < 0, MINUS, " "),
              abs(pd$r_csv),
              ifelse(is.finite(pd$r_emp) && pd$r_emp < 0, MINUS, " "),
              abs(pd$r_emp),
              formatC(pd$n, big.mark = ",", format = "d")))
}
r_vec <- sapply(panel_data_list, function(x) x$r_csv)
cat(sprintf("\n  Median r across panels: %.3f  (expected -0.663)\n",
            median(r_vec, na.rm = TRUE)))
cat(sprintf("  Sites with negative r:  %d of %d (expected 14 of 15)\n",
            sum(r_vec < 0, na.rm = TRUE),
            sum(is.finite(r_vec))))

log_progress("suppfig_S16_15panel_hexbin.R complete.")
