# =============================================================================
# gedi_als_offset_distributions.R
# -----------------------------------------------------------------------------
# The offset-distribution step -- GEDI-vs-ALS offset distributions and cross-site connection
# to the Ch1 sigma_site decomposition.
#
# the offset-distribution step covers everything not requiring GEDI footprint geometry.
# the spatial-structure step (within-site spatial structure: Moran's I, variograms, point maps)
# is deferred pending the offset-distribution step results.
#
# =============================================================================

# ---- 1. Setup --------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(moments)
  library(patchwork)
})

# Pane 1 cwd should be scripts/reviewed/
script_dir <- getwd()
if (!grepl("scripts/reviewed", script_dir)) {
  warning("CWD is not scripts/reviewed; paths assume PROJECT_ROOT anchoring.")
}

source("analysis_config.R")   # PROJECT_ROOT and friends
source("analysis_utils.R")          # save_checkpoint, log_msg

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", unset = "/gpfs/data1/vclgp/lmaden/chpt1")
DATA_DIR     <- file.path(PROJECT_ROOT, "data", "enriched_by_site")
LOOKUP_PATH  <- file.path(PROJECT_ROOT, "manuscript_tables", "site_id_lookup.csv")
OUT_TBL      <- file.path(PROJECT_ROOT, "manuscript_tables")
OUT_PLOT     <- file.path(PROJECT_ROOT, "plots", "groundwork")
dir.create(OUT_TBL,  recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_PLOT, recursive = TRUE, showWarnings = FALSE)

# RE intercept files -- pinned to the actual filenames in manuscript_tables/
# (verified column names; see comments in section 8).
CH1_RE_PATH    <- file.path(OUT_TBL, "site_random_effects_chm_dtm.csv")
PHASE3_RE_PATH <- file.path(OUT_TBL, "coupling_random_effect_shifts.csv")
TRACK_B_RE_PATH<- file.path(OUT_TBL, "track_b_re_intercepts.csv")

FOREST_LCS <- c("BDF", "DNF", "EBF", "ENF")

# ---- 2. Banner -------------------------------------------------------------
log_msg("=== gedi_als_offset_distributions.R ===")
log_msg("PROJECT_ROOT = ", PROJECT_ROOT)
log_msg("CPU_FRAC = ", Sys.getenv("CPU_FRAC", "unset"),
        " | RAM_FRAC = ", Sys.getenv("RAM_FRAC", "unset"))

# ---- 3. Site lookup --------------------------------------------------------
stopifnot(file.exists(LOOKUP_PATH))
site_lookup <- fread(LOOKUP_PATH)
sites <- site_lookup[!is.na(manuscript_site)][order(manuscript_site)]
log_msg("Loaded ", nrow(sites), " manuscript sites from lookup.")

# ---- 4. Read enriched CSVs -------------------------------------------------
# Enriched-CSV column conventions (verified via str() + charToRaw() on
# site_04_enriched.csv.gz):
#   - CHM error is `error_mean` (unprefixed; CHM was the original product).
#   - DTM error is `dtm_error_mean` (DTM was added later).
#   - LC Level-1 code with BDF/DNF/EBF/ENF/IWL/GRS values lives in
#     `lc2022_l1_code`.  NOTE: that's a lowercase-L followed by a digit-1
#     (Level-1), NOT a double-L; in most terminal fonts the two are visually
#     indistinguishable.  Same applies to `lc2022_l1_name`.
# We rename to canonical names (`chm_error_mean`, `lc_class`) right after the
# subset so the rest of the script body stays symmetric and readable.
required_cols <- c("rh_98", "als_chm_p90", "als_chm_mean",
                   "error_mean", "dtm_error_mean",
                   "cover", "wsci", "slope_mean", "lc2022_l1_code",
                   "meta_abs_geoacc_avg")

read_one_site <- function(tracker_id, manuscript_id, site_short) {
  fp_gz <- file.path(DATA_DIR, sprintf("site_%02d_enriched.csv.gz", tracker_id))
  fp    <- file.path(DATA_DIR, sprintf("site_%02d_enriched.csv",    tracker_id))
  path  <- if (file.exists(fp_gz)) fp_gz else if (file.exists(fp)) fp else NA_character_
  if (is.na(path)) {
    warning("Enriched CSV missing for tracker site ", tracker_id,
            " (manuscript ", manuscript_id, ", ", site_short, ")")
    return(NULL)
  }
  dt <- fread(path)

  # Aggressive column-name normalization: keep only [A-Za-z0-9_].  Catches
  # BOMs, NBSPs, zero-width spaces, and any other invisible-character cruft
  # without needing to enumerate variants.
  raw_names   <- names(dt)
  clean_names <- gsub("[^A-Za-z0-9_]+", "", raw_names)
  req_set     <- c(required_cols, "shot_number")
  match_idx   <- which(clean_names %in% req_set)
  miss        <- setdiff(req_set, clean_names)
  if (length(miss) > 0) {
    warning("Site ", manuscript_id, " (", site_short, ") missing cols: ",
            paste(miss, collapse = ", "))
  }

  keep_orig <- raw_names[match_idx]
  dt <- dt[, ..keep_orig]
  setnames(dt, keep_orig, clean_names[match_idx])

  dt[, `:=`(tracker_site = tracker_id,
            manuscript_site = manuscript_id,
            site_short = site_short)]
  dt
}

# Pre-loop probe: read site 1's header only and report any column names that
# need normalization.  Logs once so we can see the smoking gun if the issue
# was a hidden character.
{
  probe_path <- file.path(DATA_DIR,
                          sprintf("site_%02d_enriched.csv.gz",
                                  sites$tracker_site[1]))
  if (file.exists(probe_path)) {
    probe_raw   <- names(fread(probe_path, nrows = 0))
    probe_clean <- gsub("[^A-Za-z0-9_]+", "", probe_raw)
    diff_idx <- which(probe_raw != probe_clean)
    if (length(diff_idx) > 0) {
      log_msg("Column-name normalization probe (site_",
              sprintf("%02d", sites$tracker_site[1]), "):")
      for (i in diff_idx) {
        log_msg("  '", probe_raw[i], "' (bytes: ",
                paste(charToRaw(probe_raw[i]), collapse = " "),
                ") -> '", probe_clean[i], "'")
      }
    } else {
      log_msg("Column-name normalization probe: all column names are clean ASCII.")
    }
  }
}

log_msg("Reading per-site enriched CSVs...")
fp <- rbindlist(
  mapply(read_one_site,
         sites$tracker_site, sites$manuscript_site, sites$site_short_name,
         SIMPLIFY = FALSE),
  fill = TRUE, use.names = TRUE)
log_msg("  Loaded ", nrow(fp), " footprints across ",
        uniqueN(fp$manuscript_site), " sites.")

# Rename to canonical names used by the rest of the script.
# Done post-rbind on the final fp object so it always runs once on the full
# dataset (more robust than per-site rename for diagnosing schema surprises).
setnames(fp,
         old = c("error_mean",     "lc2022_l1_code"),
         new = c("chm_error_mean", "lc_class"),
         skip_absent = TRUE)

# Belt-and-suspenders fallbacks if the exact-name rename above missed the
# columns due to a column-name oddity that survived the trimws() in
# read_one_site (e.g. embedded zero-width chars, encoding mismatch).
if (!"chm_error_mean" %in% names(fp)) {
  cand <- grep("^error_mean$|^chm_error_mean$|^err.*chm.*mean$",
               names(fp), value = TRUE, ignore.case = TRUE)
  if (length(cand) >= 1) {
    log_msg("Fallback CHM-error rename: ", cand[1], " -> chm_error_mean")
    setnames(fp, cand[1], "chm_error_mean")
  }
}
if (!"lc_class" %in% names(fp)) {
  cand <- grep("^lc.*l[1l]_code$|^lc.*class$|^lc2022_l[1l]_code$",
               names(fp), value = TRUE, ignore.case = TRUE)
  if (length(cand) >= 1) {
    log_msg("Fallback LC-class rename: ", cand[1], " -> lc_class")
    setnames(fp, cand[1], "lc_class")
  }
}

# Assert canonical cols are present so we fail fast with a clear message
# rather than crashing in section 5 with .checkTypos.
canon <- c("chm_error_mean", "lc_class")
miss_canon <- setdiff(canon, names(fp))
if (length(miss_canon) > 0) {
  stop("After read+rename, fp is missing canonical column(s): ",
       paste(miss_canon, collapse = ", "),
       ". fp columns are: ", paste(sort(names(fp)), collapse = ", "))
}

# Sanity: print observed lc_class levels
lc_levels <- sort(unique(as.character(fp$lc_class)))
log_msg("  lc_class levels observed: ", paste(lc_levels, collapse = ", "))
if (!any(FOREST_LCS %in% lc_levels)) {
  warning("No FOREST_LCS codes found in lc_class. Adjust FOREST_LCS to match data.")
}

# Attach flag_status / dominant_lc / ecoregion
fp <- merge(fp,
            sites[, .(manuscript_site, dominant_lc, flag_status, ecoregion_id)],
            by = "manuscript_site", all.x = TRUE)

# ---- 5. Compute offsets ----------------------------------------------------
fp[, offset_p90  := rh_98 - als_chm_p90]
fp[, offset_mean := rh_98 - als_chm_mean]

n_pre <- nrow(fp)
fp <- fp[!is.na(offset_p90) & !is.na(offset_mean)]
log_msg("After NA filter on offsets: ", nrow(fp), "/", n_pre, " footprints retained.")

fp_for <- fp[lc_class %in% FOREST_LCS]
log_msg("Forest subset (lc_class in {", paste(FOREST_LCS, collapse=","),"}): ",
        nrow(fp_for), " footprints.")

# ---- 6. Per-site summary ---------------------------------------------------
make_site_summary <- function(dt, label) {
  dt[, .(
    pipeline             = label,
    n                    = .N,
    median_p90           = median(offset_p90),
    mean_p90             = mean(offset_p90),
    sd_p90               = sd(offset_p90),
    iqr_p90              = IQR(offset_p90),
    q05_p90              = quantile(offset_p90, 0.05),
    q95_p90              = quantile(offset_p90, 0.95),
    skew_p90             = skewness(offset_p90),
    kurt_p90             = kurtosis(offset_p90),
    median_mean          = median(offset_mean),
    sd_mean              = sd(offset_mean),
    r_offset_p90_chm_err = if (sum(!is.na(chm_error_mean)) > 30)
                             cor(offset_p90, chm_error_mean,
                                 use = "pairwise.complete.obs") else NA_real_,
    r_offset_p90_dtm_err = if (sum(!is.na(dtm_error_mean)) > 30)
                             cor(offset_p90, dtm_error_mean,
                                 use = "pairwise.complete.obs") else NA_real_,
    r_offset_p90_geoacc  = if (sum(!is.na(meta_abs_geoacc_avg)) > 30)
                             cor(offset_p90, meta_abs_geoacc_avg,
                                 use = "pairwise.complete.obs") else NA_real_,
    cover_avg            = mean(cover, na.rm = TRUE),
    rh_98_avg            = mean(rh_98, na.rm = TRUE),
    slope_avg            = mean(slope_mean, na.rm = TRUE),
    wsci_avg             = mean(wsci, na.rm = TRUE),
    meta_abs_geoacc_avg  = mean(meta_abs_geoacc_avg, na.rm = TRUE)
  ), by = .(manuscript_site, site_short, dominant_lc, flag_status)]
}

per_site_all <- make_site_summary(fp,     "all_lc")
per_site_for <- make_site_summary(fp_for, "forest_only")
per_site <- rbindlist(list(per_site_all, per_site_for))
setorder(per_site, pipeline, manuscript_site)
fwrite(per_site, file.path(OUT_TBL, "groundwork_task3_phase1_per_site_summary.csv"))
log_msg("Wrote per-site summary: ", nrow(per_site), " rows.")

# ---- 7. Stratified summaries ----------------------------------------------
# Global quartile breaks so bins are comparable across sites
q_break <- function(x) quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE)
slope_q <- q_break(fp$slope_mean)
cover_q <- q_break(fp$cover)

assign_bin <- function(x, q, lab) {
  cut(x, breaks = c(-Inf, q, Inf),
      labels = paste0(lab, c("_q1", "_q2", "_q3", "_q4")),
      include.lowest = TRUE)
}
fp[,     slope_bin := assign_bin(slope_mean, slope_q, "slope")]
fp[,     cover_bin := assign_bin(cover,      cover_q, "cover")]
fp_for[, slope_bin := assign_bin(slope_mean, slope_q, "slope")]
fp_for[, cover_bin := assign_bin(cover,      cover_q, "cover")]

stratify <- function(dt, by_var, stratum_label, pipeline_label) {
  if (!by_var %in% names(dt)) return(NULL)
  out <- dt[, .(
    n          = .N,
    median_p90 = median(offset_p90),
    sd_p90     = sd(offset_p90),
    iqr_p90    = IQR(offset_p90)
  ), by = c("manuscript_site", "site_short", by_var)]
  setnames(out, by_var, "stratum")
  out[, stratum := as.character(stratum)]
  out[, `:=`(pipeline = pipeline_label, stratum_var = stratum_label)]
  setcolorder(out, c("pipeline", "stratum_var", "manuscript_site",
                     "site_short", "stratum", "n",
                     "median_p90", "sd_p90", "iqr_p90"))
  out
}

strat_long <- rbindlist(list(
  stratify(fp,     "slope_bin", "slope", "all_lc"),
  stratify(fp,     "cover_bin", "cover", "all_lc"),
  stratify(fp,     "lc_class",  "lc",    "all_lc"),
  stratify(fp_for, "slope_bin", "slope", "forest_only"),
  stratify(fp_for, "cover_bin", "cover", "forest_only"),
  stratify(fp_for, "lc_class",  "lc",    "forest_only")
), fill = TRUE)
fwrite(strat_long, file.path(OUT_TBL, "groundwork_task3_phase1_stratified.csv"))
log_msg("Wrote stratified summary: ", nrow(strat_long), " rows.")

# ---- 8. Read RE intercept files (pinned columns) --------------------------
# Verified file structure:
#
#   site_random_effects_chm_dtm.csv
#     ms_site_id (= manuscript_site, 1-19 contiguous), tracker_site_id,
#     site_name, chm_re_mean (...se, lo, med, hi),
#     dtm_re_mean (...se, lo, med, hi)
#
#   coupling_random_effect_shifts.csv
#     site (= TRACKER numbering, skips 16), RE_full, RE_full_lo, RE_full_hi,
#     RE_alt18, RE_alt18_lo, RE_alt18_hi, delta, abs_delta
#     -- RE_full equals chm_re_mean from the file above (cross-checked);
#        we therefore take only RE_alt18 from this file.
#
#   track_b_re_intercepts.csv (just produced via ranef() on fit_chm_16)
#     tracker_site, track_b_chm_re, track_b_chm_re_q025, track_b_chm_re_q975
#     -- 16 rows; flagged sites 1-3 absent.

read_re_csv <- function(path, label) {
  if (!file.exists(path)) {
    log_msg("[", label, "] not found at ", path)
    return(NULL)
  }
  dt <- fread(path)
  log_msg("[", label, "] ", nrow(dt), " rows; cols: ",
          paste(names(dt), collapse = ", "))
  dt
}

ch1_re     <- read_re_csv(CH1_RE_PATH,     "Ch1 site REs")
phase3_re  <- read_re_csv(PHASE3_RE_PATH,  "Decomposition RE shifts")
track_b_re <- read_re_csv(TRACK_B_RE_PATH, "16-site RE")

# Build a per-site RE table keyed on manuscript_site.
re_long <- data.table(manuscript_site = sites$manuscript_site)

if (!is.null(ch1_re)) {
  ch1_slim <- ch1_re[, .(manuscript_site = ms_site_id,
                         ch1_chm_re      = chm_re_mean,
                         ch1_dtm_re      = dtm_re_mean)]
  re_long <- merge(re_long, ch1_slim, by = "manuscript_site", all.x = TRUE)
}

if (!is.null(phase3_re)) {
  # `site` column in this file is tracker_site -> map to manuscript_site
  ph3_slim <- merge(phase3_re[, .(tracker_site = site,
                                  phase3_18site_re_alt = RE_alt18)],
                    sites[, .(tracker_site, manuscript_site)],
                    by = "tracker_site", all.x = TRUE)
  ph3_slim <- ph3_slim[!is.na(manuscript_site),
                       .(manuscript_site, phase3_18site_re_alt)]
  re_long <- merge(re_long, ph3_slim, by = "manuscript_site", all.x = TRUE)
}

if (!is.null(track_b_re)) {
  tb_slim <- merge(track_b_re[, .(tracker_site, track_b_chm_re)],
                   sites[, .(tracker_site, manuscript_site)],
                   by = "tracker_site", all.x = TRUE)
  tb_slim <- tb_slim[!is.na(manuscript_site),
                     .(manuscript_site, track_b_chm_re)]
  re_long <- merge(re_long, tb_slim, by = "manuscript_site", all.x = TRUE)
}

log_msg("Per-site RE table assembled with columns: ",
        paste(names(re_long), collapse = ", "))

# ---- 9. Cross-site correlations -------------------------------------------
flagged_sites <- sites[flag_status == "FLAGGED"]$manuscript_site
dtm_excluded  <- sites[flag_status == "DTM_excluded"]$manuscript_site

xs <- merge(per_site_all, re_long, by = "manuscript_site", all.x = TRUE)

frames <- list(
  `19_all`               = xs,
  `16_nonflagged`        = xs[!manuscript_site %in% flagged_sites],
  `18_dtm_ok`            = xs[!manuscript_site %in% dtm_excluded],
  `15_nonflagged_dtm_ok` = xs[!manuscript_site %in% c(flagged_sites, dtm_excluded)]
)

targets <- c("median_p90", "sd_p90", "iqr_p90",
             "r_offset_p90_chm_err", "r_offset_p90_dtm_err")
predictors <- intersect(
  c("meta_abs_geoacc_avg", "cover_avg", "rh_98_avg", "slope_avg", "wsci_avg",
    "ch1_chm_re", "ch1_dtm_re", "track_b_chm_re", "phase3_18site_re_alt"),
  names(xs))

cor_one <- function(dt, x, y, frame_label) {
  if (!all(c(x, y) %in% names(dt))) return(NULL)
  z <- dt[, .(x = get(x), y = get(y))][!is.na(x) & !is.na(y)]
  if (nrow(z) < 4) return(NULL)
  data.table(
    frame        = frame_label,
    n            = nrow(z),
    target       = x,
    predictor    = y,
    pearson_r    = suppressWarnings(cor(z$x, z$y, method = "pearson")),
    spearman_rho = suppressWarnings(cor(z$x, z$y, method = "spearman")),
    pearson_p    = tryCatch(cor.test(z$x, z$y, method = "pearson")$p.value,
                            error = function(e) NA_real_),
    spearman_p   = tryCatch(cor.test(z$x, z$y, method = "spearman")$p.value,
                            error = function(e) NA_real_)
  )
}

cor_rows <- list()
for (fname in names(frames)) {
  f <- frames[[fname]]
  for (t in targets) for (p in predictors) {
    cor_rows[[length(cor_rows) + 1]] <- cor_one(f, t, p, fname)
  }
}
cross_cor <- rbindlist(cor_rows, fill = TRUE)
fwrite(cross_cor, file.path(OUT_TBL, "groundwork_task3_phase1_cross_site_correlations.csv"))
log_msg("Wrote cross-site correlations: ", nrow(cross_cor), " rows.")

# ---- 10. Plots ------------------------------------------------------------
flag_pal <- c(ok = "#2C7FB8", FLAGGED = "#D7301F", DTM_excluded = "#FDB863")

make_panel <- function(dt, title_suffix, file) {
  dt2 <- copy(dt)
  dt2[, panel_label := sprintf("%02d %s", manuscript_site, site_short)]
  panel_order <- dt2[, .(o = manuscript_site[1]), by = panel_label][
    order(o)]$panel_label
  dt2[, panel_label := factor(panel_label, levels = panel_order)]

  p <- ggplot(dt2, aes(offset_p90, fill = flag_status)) +
    geom_histogram(bins = 80, alpha = 0.85) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    facet_wrap(~ panel_label, ncol = 4, scales = "free_y") +
    scale_fill_manual(values = flag_pal, na.value = "grey60") +
    coord_cartesian(xlim = c(-25, 15)) +
    labs(title = paste("rh_98 - als_chm_p90 by site:", title_suffix),
         x = "Offset (m)", y = "Footprint count", fill = "Site flag") +
    theme_minimal(base_size = 9) +
    theme(strip.text = element_text(size = 8))
  ggsave(file, p, width = 11, height = 13)
}

make_panel(fp,     "all LC",
           file.path(OUT_PLOT, "groundwork_task3_phase1_offset_distributions_alllc.pdf"))
make_panel(fp_for, "forest only",
           file.path(OUT_PLOT, "groundwork_task3_phase1_offset_distributions_forest.pdf"))

make_xs_scatter <- function(dt, x, y, frame_label) {
  if (!all(c(x, y) %in% names(dt))) return(NULL)
  z <- dt[, .(x = get(x), y = get(y),
              site = manuscript_site, flag = flag_status)][!is.na(x) & !is.na(y)]
  if (nrow(z) < 4) return(NULL)
  ggplot(z, aes(x, y, label = site, colour = flag)) +
    geom_point(size = 2.4) +
    geom_text(nudge_y = 0.04 * diff(range(z$y, na.rm = TRUE)), size = 3) +
    scale_colour_manual(values = flag_pal, na.value = "grey60") +
    labs(x = x, y = y,
         title = sprintf("%s  (n=%d)", frame_label, nrow(z))) +
    theme_minimal(base_size = 10)
}

xs_predictors <- intersect(
  c("meta_abs_geoacc_avg", "cover_avg",
    "ch1_chm_re", "ch1_dtm_re", "track_b_chm_re", "phase3_18site_re_alt"),
  names(xs))
xs_plots <- list()
for (yvar in c("median_p90", "sd_p90", "r_offset_p90_chm_err")) {
  for (xvar in xs_predictors) {
    xs_plots[[paste(yvar, xvar, sep = "__")]] <-
      make_xs_scatter(xs, xvar, yvar, sprintf("19-site: %s vs %s", yvar, xvar))
  }
}
xs_plots <- xs_plots[!sapply(xs_plots, is.null)]
if (length(xs_plots) > 0) {
  composite <- wrap_plots(xs_plots, ncol = 3) +
    plot_annotation(title = "The offset-distribution step - cross-site offset summaries")
  ggsave(file.path(OUT_PLOT, "groundwork_task3_phase1_cross_site_summary.pdf"),
         composite,
         width = 14,
         height = 4 * ceiling(length(xs_plots) / 3),
         limitsize = FALSE)
}

# ---- 11. Checkpoint -------------------------------------------------------
save_checkpoint("groundwork_task3_phase1", list(
  per_site         = per_site,
  strat_long       = strat_long,
  cross_cor        = cross_cor,
  re_long          = re_long,
  fp_offsets_slim  = fp[, .(manuscript_site, shot_number,
                             offset_p90, offset_mean,
                             chm_error_mean, dtm_error_mean,
                             cover, slope_mean, lc_class,
                             meta_abs_geoacc_avg)]
))

log_msg("=== The offset-distribution step complete ===")
log_msg("CSVs:  ", OUT_TBL,
        "/groundwork_task3_phase1_{per_site_summary,stratified,cross_site_correlations}.csv")
log_msg("PDFs:  ", OUT_PLOT,
        "/groundwork_task3_phase1_{offset_distributions_alllc,offset_distributions_forest,cross_site_summary}.pdf")
log_msg("Checkpoint: groundwork_task3_phase1")
