#!/usr/bin/env Rscript
# =====================================================================
# holdout_01_posterior_values.R
#
# Step 01 (cheap pre-flight): extract posterior-summary values
# needed to sign off on captions for S17, S9b (absorbed), and S21
# BEFORE the slow step 03 PPC and variogram passes run.
#
# v2 change: probe checkpoint contents for the brmsfit object rather
# than hard-coding key names. v1 assumed key 'fit_chm_s2_sensitivity'
# but section_10b actually saves under 'fit_chm_16'; v2 finds either.
#
# Inputs (in priority order):
#   1. checkpoints/10b_chm_sensitivity.rds  (16-site CHM)
#   2. checkpoints/10_models_stage2.rds     (18-site DTM)
#
# Outputs (all under manuscript_tables/):
#   section_L_s17_chm_site_re_summary.csv
#   section_L_s9b_chm_grouping_levels.csv
#   section_L_s9b_dtm_grouping_levels.csv
#   section_L_s9b_groupwise_sigma.csv
#   section_L_s21_chm_site_re_metadata.csv
#   section_L_s21_chm_site_re_bivariate_r2.csv
# =====================================================================

suppressPackageStartupMessages({
  library(brms)
  library(posterior)
  library(data.table)
})

Sys.setenv(DISPLAY = "")
options(device = pdf)

if (!exists("log_progress"))  source("analysis_utils.R")
if (!exists("PROJECT_ROOT"))  source("analysis_config.R")

manuscript_tables_dir <- file.path(PROJECT_ROOT, "manuscript_tables")
if (!dir.exists(manuscript_tables_dir)) dir.create(manuscript_tables_dir, recursive = TRUE)

log_progress("================================================================")
log_progress("Step 01 (v2): pre-flight posterior summaries")
log_progress("  Targets: S17 r_site, S9b ranef levels, S21 RE-metadata bivariates")
log_progress("================================================================")

# ---------------------------------------------------------------------
# Defensive brmsfit discovery: given a loaded checkpoint list, find
# the first element that is a brmsfit object regardless of key name
# ---------------------------------------------------------------------

find_brmsfit_in_list <- function(x, label = "<checkpoint>") {
  if (inherits(x, "brmsfit")) return(x)
  if (!is.list(x)) {
    stop(sprintf("%s: not a brmsfit and not a list (class: %s)",
                 label, paste(class(x), collapse = "/")))
  }
  hits <- vapply(x, inherits, logical(1), what = "brmsfit")
  if (sum(hits) == 0) {
    stop(sprintf("%s: no brmsfit object found among keys: %s",
                 label, paste(names(x), collapse = ", ")))
  }
  if (sum(hits) > 1) {
    log_progress(sprintf("  %s: multiple brmsfits found (%s); using first: %s",
                         label,
                         paste(names(x)[hits], collapse = ", "),
                         names(x)[hits][1]))
  } else {
    log_progress(sprintf("  %s: brmsfit key = '%s'",
                         label, names(x)[hits][1]))
  }
  x[[ which(hits)[1] ]]
}

# ---------------------------------------------------------------------
# Load CHM 16-site (from 10b_chm_sensitivity checkpoint)
# ---------------------------------------------------------------------

log_progress("Loading 16-site CHM fit from 10b_chm_sensitivity checkpoint")
cp_chm <- load_checkpoint("10b_chm_sensitivity")
log_progress(sprintf("  Checkpoint keys: %s",
                     paste(names(cp_chm), collapse = ", ")))
fit_chm <- find_brmsfit_in_list(cp_chm, "10b_chm_sensitivity")
stopifnot(inherits(fit_chm, "brmsfit"))

# ---------------------------------------------------------------------
# Load DTM 18-site (from 10_models_stage2 checkpoint)
# ---------------------------------------------------------------------

log_progress("Loading 18-site DTM fit from 10_models_stage2 checkpoint")
cp_dtm <- load_checkpoint("10_models_stage2")
log_progress(sprintf("  Checkpoint keys: %s",
                     paste(names(cp_dtm), collapse = ", ")))
# The Stage 2 checkpoint contains both fit_chm_s2 and fit_dtm_s2; pick DTM
fit_dtm <- NULL
if ("fit_dtm_s2" %in% names(cp_dtm)) {
  fit_dtm <- cp_dtm$fit_dtm_s2
} else {
  dtm_named <- grep("dtm", names(cp_dtm), ignore.case = TRUE, value = TRUE)
  for (nm in dtm_named) {
    if (inherits(cp_dtm[[nm]], "brmsfit")) {
      fit_dtm <- cp_dtm[[nm]]
      log_progress(sprintf("  Resolved DTM fit via key '%s'", nm))
      break
    }
  }
}
stopifnot(inherits(fit_dtm, "brmsfit"))

# Confirm grouping levels
chm_grp <- names(ranef(fit_chm, summary = FALSE))
dtm_grp <- names(ranef(fit_dtm, summary = FALSE))
log_progress(sprintf("CHM grouping levels: %s", paste(chm_grp, collapse = ", ")))
log_progress(sprintf("DTM grouping levels: %s", paste(dtm_grp, collapse = ", ")))

guess_site_group <- function(grp_names) {
  cands <- c("site_id", "site", "siteID", "manuscript_site_id", "tracker_site_id")
  hit <- intersect(cands, grp_names)
  if (length(hit) > 0) return(hit[1])
  sg <- grep("site", grp_names, ignore.case = TRUE, value = TRUE)
  if (length(sg) > 0) return(sg[1])
  stop(sprintf("Cannot locate site grouping factor among: %s",
               paste(grp_names, collapse = ", ")))
}
chm_site_grp <- guess_site_group(chm_grp)
dtm_site_grp <- guess_site_group(dtm_grp)
log_progress(sprintf("CHM site grouping resolved to: '%s'", chm_site_grp))
log_progress(sprintf("DTM site grouping resolved to: '%s'", dtm_site_grp))

# ---------------------------------------------------------------------
# S17: CHM 16-site per-site RE summary
# ---------------------------------------------------------------------

log_subsection("S17: CHM 16-site per-site random intercepts")

chm_site_re <- ranef(fit_chm, summary = TRUE,
                     probs = c(0.025, 0.5, 0.975))[[chm_site_grp]]
chm_re_df <- data.table(
  site_label = dimnames(chm_site_re)[[1]],
  re_mean   = chm_site_re[, "Estimate",  "Intercept"],
  re_sd     = chm_site_re[, "Est.Error", "Intercept"],
  re_q025   = chm_site_re[, "Q2.5",      "Intercept"],
  re_q50    = chm_site_re[, "Q50",       "Intercept"],
  re_q975   = chm_site_re[, "Q97.5",     "Intercept"]
)
chm_re_df[, sort_key := suppressWarnings(as.integer(site_label))]
setorder(chm_re_df, sort_key, na.last = TRUE)
chm_re_df[, sort_key := NULL]

fwrite(chm_re_df,
       file.path(manuscript_tables_dir,
                 "section_L_s17_chm_site_re_summary.csv"))

s17_range <- range(chm_re_df$re_mean)
s17_q_range <- c(min(chm_re_df$re_q025), max(chm_re_df$re_q975))

log_progress(sprintf("  N sites in CHM ranef table: %d", nrow(chm_re_df)))
log_progress(sprintf("  CHM 16-site RE posterior-mean range: [%+.3f, %+.3f] m",
                     s17_range[1], s17_range[2]))
log_progress(sprintf("  CHM 16-site RE 95%% CI envelope:      [%+.3f, %+.3f] m",
                     s17_q_range[1], s17_q_range[2]))
log_progress(sprintf("  Site label contributing min RE: %s (mean %+.3f m)",
                     chm_re_df$site_label[which.min(chm_re_df$re_mean)],
                     min(chm_re_df$re_mean)))
log_progress(sprintf("  Site label contributing max RE: %s (mean %+.3f m)",
                     chm_re_df$site_label[which.max(chm_re_df$re_mean)],
                     max(chm_re_df$re_mean)))

dtm_site_re <- ranef(fit_dtm, summary = TRUE,
                     probs = c(0.025, 0.5, 0.975))[[dtm_site_grp]]
dtm_re_df <- data.table(
  site_label = dimnames(dtm_site_re)[[1]],
  re_mean   = dtm_site_re[, "Estimate",  "Intercept"],
  re_q025   = dtm_site_re[, "Q2.5",      "Intercept"],
  re_q975   = dtm_site_re[, "Q97.5",     "Intercept"]
)
dtm_re_df[, sort_key := suppressWarnings(as.integer(site_label))]
setorder(dtm_re_df, sort_key, na.last = TRUE)
dtm_re_df[, sort_key := NULL]

log_progress(sprintf("  DTM 18-site RE posterior-mean range:  [%+.3f, %+.3f] m  (context)",
                     min(dtm_re_df$re_mean), max(dtm_re_df$re_mean)))

# ---------------------------------------------------------------------
# S9b: ranef summaries across ALL grouping levels for both products
# ---------------------------------------------------------------------

log_subsection("S9b: per-level ranef summaries across grouping factors")

extract_all_levels <- function(fit, product_label) {
  re_list <- ranef(fit, summary = TRUE, probs = c(0.025, 0.5, 0.975))
  out <- rbindlist(lapply(names(re_list), function(grp) {
    arr <- re_list[[grp]]
    data.table(
      product       = product_label,
      grouping      = grp,
      level         = dimnames(arr)[[1]],
      ranef_mean    = arr[, "Estimate",  "Intercept"],
      ranef_sd      = arr[, "Est.Error", "Intercept"],
      ranef_q025    = arr[, "Q2.5",      "Intercept"],
      ranef_q50     = arr[, "Q50",       "Intercept"],
      ranef_q975    = arr[, "Q97.5",     "Intercept"]
    )
  }))
  out
}

s9b_chm <- extract_all_levels(fit_chm, "CHM_16_site")
s9b_dtm <- extract_all_levels(fit_dtm, "DTM_18_site")

fwrite(s9b_chm, file.path(manuscript_tables_dir,
                          "section_L_s9b_chm_grouping_levels.csv"))
fwrite(s9b_dtm, file.path(manuscript_tables_dir,
                          "section_L_s9b_dtm_grouping_levels.csv"))

log_progress(sprintf("  S9b CHM rows: %d (across %d grouping factors: %s)",
                     nrow(s9b_chm),
                     length(unique(s9b_chm$grouping)),
                     paste(unique(s9b_chm$grouping), collapse = ", ")))
log_progress(sprintf("  S9b DTM rows: %d (across %d grouping factors: %s)",
                     nrow(s9b_dtm),
                     length(unique(s9b_dtm$grouping)),
                     paste(unique(s9b_dtm$grouping), collapse = ", ")))

extract_groupwise_sigma <- function(fit, product_label) {
  vc <- VarCorr(fit, summary = TRUE, probs = c(0.025, 0.5, 0.975))
  out <- rbindlist(lapply(names(vc), function(grp) {
    sd_arr <- vc[[grp]]$sd
    if (is.null(sd_arr)) return(NULL)
    data.table(
      product       = product_label,
      grouping      = grp,
      parameter     = rownames(sd_arr),
      sigma_mean    = sd_arr[, "Estimate"],
      sigma_q025    = sd_arr[, "Q2.5"],
      sigma_q50     = sd_arr[, "Q50"],
      sigma_q975    = sd_arr[, "Q97.5"]
    )
  }))
  out
}

sigma_chm <- extract_groupwise_sigma(fit_chm, "CHM_16_site")
sigma_dtm <- extract_groupwise_sigma(fit_dtm, "DTM_18_site")
sigma_all <- rbind(sigma_chm, sigma_dtm)
fwrite(sigma_all, file.path(manuscript_tables_dir,
                            "section_L_s9b_groupwise_sigma.csv"))

log_progress("  Groupwise sigma posterior medians (for cross-check vs story-lock):")
for (i in seq_len(nrow(sigma_all))) {
  log_progress(sprintf("    %s | %-14s | %-12s | sigma = %.3f [%.3f, %.3f]",
                       sigma_all$product[i],
                       sigma_all$grouping[i],
                       sigma_all$parameter[i],
                       sigma_all$sigma_q50[i],
                       sigma_all$sigma_q025[i],
                       sigma_all$sigma_q975[i]))
}

# ---------------------------------------------------------------------
# S21: CHM site RE vs metadata bivariates
# ---------------------------------------------------------------------

log_subsection("S21: CHM 16-site RE vs metadata covariates")

site_cov_path <- file.path(manuscript_tables_dir,
                           "coupling_site_covariates.csv")
site_cov <- NULL
if (file.exists(site_cov_path)) {
  log_progress(sprintf("  Reading per-site metadata covariates from %s",
                       basename(site_cov_path)))
  site_cov <- fread(site_cov_path)
  log_progress(sprintf("  CSV columns: %s",
                       paste(names(site_cov), collapse = ", ")))
} else {
  log_progress("  per-site covariates CSV not found")
  log_progress("  Falling back to aggregation from CHM data ingest checkpoint")
  if (checkpoint_exists("01_data_ingest")) {
    di <- load_checkpoint("01_data_ingest")
    chm_df <- as.data.table(di$chm_df)
    log_progress(sprintf("  Aggregating from %d footprints", nrow(chm_df)))
    meta_cols <- intersect(c("meta_fwd_ratio", "meta_rev_ratio",
                             "meta_off_nadir_avg", "meta_sun_elev_avg",
                             "meta_abs_geoacc_avg", "meta_rel_geoacc_avg",
                             "meta_leafon", "meta_stereo_pairs",
                             "meta_az_concurrency"),
                           names(chm_df))
    log_progress(sprintf("  Metadata cols available: %s",
                         paste(meta_cols, collapse = ", ")))
    site_key <- intersect(c("site_id", "manuscript_site_id", "site"),
                          names(chm_df))[1]
    site_cov <- chm_df[, lapply(.SD, mean, na.rm = TRUE),
                       by = c(site_key),
                       .SDcols = meta_cols]
    setnames(site_cov, site_key, "site_label")
  } else {
    stop("Cannot find per-site metadata covariates and no data checkpoint.")
  }
}

key_candidates <- c("site_label", "manuscript_site_id", "site_id",
                    "tracker_site_id", "site")
hit <- intersect(key_candidates, names(site_cov))
if (length(hit) == 0) {
  stop(sprintf("site_cov has no recognized site-key column. Columns: %s",
               paste(names(site_cov), collapse = ", ")))
}
if (hit[1] != "site_label") {
  log_progress(sprintf("  Renaming site key '%s' -> 'site_label'", hit[1]))
  setnames(site_cov, hit[1], "site_label")
}

chm_re_df[, site_label := as.character(site_label)]
site_cov[, site_label := as.character(site_label)]

chm_labels <- sort(unique(chm_re_df$site_label))
cov_labels <- sort(unique(site_cov$site_label))
log_progress(sprintf("  CHM ranef site labels: %s",
                     paste(chm_labels, collapse = ", ")))
log_progress(sprintf("  Site-cov labels:        %s",
                     paste(cov_labels, collapse = ", ")))

s21_merged <- merge(chm_re_df[, .(site_label, re_mean)],
                    site_cov, by = "site_label", all.x = TRUE)

fwrite(s21_merged, file.path(manuscript_tables_dir,
                             "section_L_s21_chm_site_re_metadata.csv"))

meta_cols <- setdiff(names(s21_merged), c("site_label", "re_mean"))
biv <- rbindlist(lapply(meta_cols, function(cv) {
  x <- s21_merged[[cv]]
  if (!is.numeric(x)) return(NULL)
  y <- s21_merged$re_mean
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 4) return(NULL)
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "pearson"))
  data.table(
    covariate = cv,
    n_sites   = sum(ok),
    pearson_r = unname(ct$estimate),
    p_value   = ct$p.value,
    r_squared = unname(ct$estimate)^2
  )
}))
if (nrow(biv) == 0) {
  log_progress("  WARNING: no numeric metadata columns yielded valid bivariate fits")
} else {
  biv <- biv[order(-r_squared)]
  fwrite(biv, file.path(manuscript_tables_dir,
                        "section_L_s21_chm_site_re_bivariate_r2.csv"))
  log_progress("  Top metadata covariates by site-RE bivariate r^2:")
  for (i in seq_len(min(8, nrow(biv)))) {
    log_progress(sprintf("    %-30s | n = %2d | r = %+.3f | p = %.4f | r^2 = %.3f",
                         biv$covariate[i],
                         biv$n_sites[i],
                         biv$pearson_r[i],
                         biv$p_value[i],
                         biv$r_squared[i]))
  }
}

# ---------------------------------------------------------------------
# Console summary
# ---------------------------------------------------------------------

log_progress("================================================================")
log_progress("SUMMARY")
log_progress("================================================================")
log_progress("")
log_progress("[S17 RE range — CHM 16-site]")
log_progress(sprintf("  posterior-mean range: [%+.2f, %+.2f] m",
                     s17_range[1], s17_range[2]))
log_progress(sprintf("  95%% CI envelope:      [%+.2f, %+.2f] m",
                     s17_q_range[1], s17_q_range[2]))
log_progress("")
log_progress("[S9b sigmas — for cross-check]")
for (i in seq_len(nrow(sigma_all))) {
  log_progress(sprintf("  %s | %-14s | sigma = %.3f [%.3f, %.3f]",
                       sigma_all$product[i],
                       sigma_all$grouping[i],
                       sigma_all$sigma_q50[i],
                       sigma_all$sigma_q025[i],
                       sigma_all$sigma_q975[i]))
}
log_progress("")
if (exists("biv") && nrow(biv) > 0) {
  log_progress("[S21 top covariates by r^2 — for caption rewrite]")
  for (i in seq_len(min(5, nrow(biv)))) {
    log_progress(sprintf("  %-30s | r = %+.3f | r^2 = %.3f",
                         biv$covariate[i],
                         biv$pearson_r[i],
                         biv$r_squared[i]))
  }
}
log_progress("")
log_progress("Outputs written:")
log_progress(sprintf("  %s/section_L_s17_chm_site_re_summary.csv", manuscript_tables_dir))
log_progress(sprintf("  %s/section_L_s9b_chm_grouping_levels.csv", manuscript_tables_dir))
log_progress(sprintf("  %s/section_L_s9b_dtm_grouping_levels.csv", manuscript_tables_dir))
log_progress(sprintf("  %s/section_L_s9b_groupwise_sigma.csv", manuscript_tables_dir))
log_progress(sprintf("  %s/section_L_s21_chm_site_re_metadata.csv", manuscript_tables_dir))
log_progress(sprintf("  %s/section_L_s21_chm_site_re_bivariate_r2.csv", manuscript_tables_dir))
log_progress("================================================================")
log_progress("Pre-flight complete.")
log_progress("================================================================")
