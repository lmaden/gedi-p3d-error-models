#!/usr/bin/env Rscript
# =====================================================================
# holdout_diagnostic.R
#
# Single-pass diagnostic: probes EVERYTHING needed to write working
# the downstream holdout scripts. No computation, just metadata extraction
# and column-name verification. Runs in under a minute.
#
# Run once. Output is the column / metadata reference for all
# the downstream holdout scripts.
# =====================================================================

suppressPackageStartupMessages({
  library(brms)
  library(data.table)
})

if (!exists("log_progress"))  source("analysis_utils.R")
if (!exists("PROJECT_ROOT"))  source("analysis_config.R")

cat("\n================================================================\n")
cat("Section L diagnostic — probing all metadata in one pass\n")
cat("================================================================\n\n")

# ---------------------------------------------------------------------
# 1. CHM brmsfit metadata
# ---------------------------------------------------------------------

cat("--- (1) CHM 16-site brmsfit metadata ---\n")
cp_chm <- load_checkpoint("10b_chm_sensitivity")
fit_chm <- cp_chm$fit_chm_16

cat(sprintf("  Class:                %s\n",
            paste(class(fit_chm), collapse = "/")))
cat(sprintf("  N observations:       %d\n", nrow(fit_chm$data)))
cat(sprintf("  brms formula$resp (Stan name, may differ from data column):\n"))
cat(sprintf("    '%s'\n", fit_chm$formula$resp))
cat(sprintf("  Full formula:\n"))
cat(sprintf("    %s\n",
            paste(deparse(fit_chm$formula$formula), collapse = " ")))
if (!is.null(fit_chm$formula$pforms)) {
  cat(sprintf("  Auxiliary parameter formulas (sigma, etc.):\n"))
  for (nm in names(fit_chm$formula$pforms)) {
    cat(sprintf("    %s: %s\n", nm,
                paste(deparse(fit_chm$formula$pforms[[nm]]),
                      collapse = " ")))
  }
}
cat(sprintf("  Columns in fit_chm$data (%d total):\n", ncol(fit_chm$data)))
fit_cols <- names(fit_chm$data)
for (i in seq(1, length(fit_cols), by = 6)) {
  end <- min(i + 5, length(fit_cols))
  cat(sprintf("    %s\n", paste(fit_cols[i:end], collapse = ", ")))
}
cat(sprintf("\n  Identifying response column in fit_chm$data...\n"))
# brms strips underscores in $formula$resp; reverse-engineer
resp_stripped <- fit_chm$formula$resp
# Try matching by stripping underscores from each column name
fit_cols_stripped <- gsub("_", "", fit_cols)
resp_match_idx <- which(fit_cols_stripped == resp_stripped)
if (length(resp_match_idx) >= 1) {
  response_col <- fit_cols[resp_match_idx[1]]
  cat(sprintf("    Resolved: brms resp '%s' -> data column '%s'\n",
              resp_stripped, response_col))
} else {
  cat(sprintf("    WARNING: no column in fit_chm$data matches '%s' after underscore strip\n",
              resp_stripped))
  response_col <- NA_character_
}

# Confirm via value range
if (!is.na(response_col)) {
  rng <- range(fit_chm$data[[response_col]], na.rm = TRUE)
  cat(sprintf("    Response value range: [%+.2f, %+.2f]\n", rng[1], rng[2]))
}

# Grouping factors
cat(sprintf("\n  Grouping factors (random effects):\n"))
ranef_list <- ranef(fit_chm, summary = FALSE)
for (grp in names(ranef_list)) {
  cat(sprintf("    %s: %d levels: %s\n",
              grp,
              dim(ranef_list[[grp]])[2],
              paste(dimnames(ranef_list[[grp]])[[2]], collapse = ", ")))
}

# ---------------------------------------------------------------------
# 2. DTM brmsfit metadata
# ---------------------------------------------------------------------

cat("\n--- (2) DTM 18-site brmsfit metadata ---\n")
cp_dtm <- load_checkpoint("10_models_stage2")
fit_dtm <- cp_dtm$fit_dtm_s2

cat(sprintf("  N observations:    %d\n", nrow(fit_dtm$data)))
cat(sprintf("  brms formula$resp: '%s'\n", fit_dtm$formula$resp))
dtm_fit_cols <- names(fit_dtm$data)
dtm_resp_match_idx <- which(gsub("_", "", dtm_fit_cols) == fit_dtm$formula$resp)
if (length(dtm_resp_match_idx) >= 1) {
  dtm_response_col <- dtm_fit_cols[dtm_resp_match_idx[1]]
  cat(sprintf("  Resolved DTM response data column: '%s'\n", dtm_response_col))
}

# ---------------------------------------------------------------------
# 3. CHM data ingest column inventory
# ---------------------------------------------------------------------

cat("\n--- (3) CHM data ingest column inventory ---\n")
di <- load_checkpoint("01_data_ingest")
chm_df <- as.data.table(di$chm_df)
cat(sprintf("  N rows:  %d\n", nrow(chm_df)))
cat(sprintf("  N cols:  %d\n", ncol(chm_df)))
cat(sprintf("  All columns:\n"))
all_cols <- sort(names(chm_df))
for (i in seq(1, length(all_cols), by = 4)) {
  end <- min(i + 3, length(all_cols))
  cat(sprintf("    %s\n", paste(all_cols[i:end], collapse = ", ")))
}

# Spot-check: which columns have the response data?
cat(sprintf("\n  Columns containing 'error' in name:\n"))
err_cols <- grep("error", names(chm_df), ignore.case = TRUE, value = TRUE)
for (cv in err_cols) {
  v <- chm_df[[cv]]
  n_fin <- sum(is.finite(v))
  if (n_fin > 0) {
    rng <- range(v, na.rm = TRUE)
    cat(sprintf("    %-25s: n_finite = %d, range [%+.2f, %+.2f]\n",
                cv, n_fin, rng[1], rng[2]))
  } else {
    cat(sprintf("    %-25s: all NA\n", cv))
  }
}

# Coordinate columns
cat(sprintf("\n  Columns that look like coordinates:\n"))
coord_candidates <- c("x", "y", "lon", "lat", "longitude", "latitude",
                     "easting", "northing", "x_alb", "y_alb")
present_coords <- intersect(coord_candidates, names(chm_df))
for (cv in present_coords) {
  v <- chm_df[[cv]]
  n_fin <- sum(is.finite(v))
  rng <- if (n_fin > 0) range(v, na.rm = TRUE) else c(NA, NA)
  cat(sprintf("    %-12s: n_finite = %d, range [%+.3f, %+.3f]\n",
              cv, n_fin, rng[1], rng[2]))
}

# Site column and its values
cat(sprintf("\n  Site grouping column:\n"))
if ("site" %in% names(chm_df)) {
  sites <- sort(unique(as.character(chm_df$site)))
  cat(sprintf("    site column present; %d unique values: %s\n",
              length(sites), paste(sites, collapse = ", ")))
  # Per-site counts (top 5 / bottom 5)
  site_counts <- chm_df[, .N, by = site][order(-N)]
  cat(sprintf("    Top 3 sites by N: %s\n",
              paste(sprintf("Site %s (n=%d)",
                            site_counts$site[1:3],
                            site_counts$N[1:3]),
                    collapse = "; ")))
  cat(sprintf("    Bottom 3 sites by N: %s\n",
              paste(sprintf("Site %s (n=%d)",
                            site_counts$site[(nrow(site_counts)-2):nrow(site_counts)],
                            site_counts$N[(nrow(site_counts)-2):nrow(site_counts)]),
                    collapse = "; ")))
}

# ---------------------------------------------------------------------
# 4. Predictor mapping: brmsfit predictors <-> chm_df columns
# ---------------------------------------------------------------------

cat("\n--- (4) Predictor column mapping (fit -> chm_df) ---\n")
fit_predictors <- setdiff(fit_cols,
                        c(response_col, "site", "ecoregion", "lc_l1_code"))
cat(sprintf("  %d predictor columns in fit_chm$data:\n", length(fit_predictors)))
for (pv in fit_predictors) {
  in_chm <- pv %in% names(chm_df)
  cat(sprintf("    %-20s: %s\n", pv,
              if (in_chm) "present in chm_df ✓" else "MISSING from chm_df"))
}

# ---------------------------------------------------------------------
# 5. 16-site subset summary (for residuals/variogram setup)
# ---------------------------------------------------------------------

cat("\n--- (5) 16-site frame summary (residuals + variogram setup) ---\n")
fit_sites <- sort(unique(as.character(fit_chm$data$site)))
cat(sprintf("  Sites in fit: %s\n", paste(fit_sites, collapse = ", ")))
chm_df_16 <- chm_df[as.character(site) %in% fit_sites]
cat(sprintf("  N footprints in 16-site frame: %d\n", nrow(chm_df_16)))
if (!is.na(response_col)) {
  finite_resp <- is.finite(chm_df_16[[response_col]])
  cat(sprintf("  Finite response: %d (%.1f%%)\n",
              sum(finite_resp), 100 * mean(finite_resp)))
}
if (length(present_coords) >= 2) {
  cx <- present_coords[1]; cy <- present_coords[2]
  finite_xy <- is.finite(chm_df_16[[cx]]) & is.finite(chm_df_16[[cy]])
  cat(sprintf("  Finite (%s, %s): %d (%.1f%%)\n", cx, cy,
              sum(finite_xy), 100 * mean(finite_xy)))
}

# ---------------------------------------------------------------------
# 6. Existing step 01 / step 02 output files
# ---------------------------------------------------------------------

manuscript_tables_dir <- file.path(PROJECT_ROOT, "manuscript_tables")
cat("\n--- (6) the holdout sequence outputs already on disk ---\n")
section_l_files <- list.files(manuscript_tables_dir,
                              pattern = "^section_L_",
                              full.names = FALSE)
if (length(section_l_files) == 0) {
  cat("  (none)\n")
} else {
  for (f in section_l_files) {
    fp <- file.path(manuscript_tables_dir, f)
    info <- file.info(fp)
    cat(sprintf("  %-50s  %8.1f KB  %s\n",
                f, info$size / 1024,
                format(info$mtime, "%Y-%m-%d %H:%M")))
  }
}

cat("\n================================================================\n")
cat("Diagnostic complete.\n")
cat("================================================================\n")
