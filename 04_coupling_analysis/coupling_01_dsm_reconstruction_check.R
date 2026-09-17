# =====================================================================
# coupling_01_dsm_reconstruction_check.R
#
# Step 01: Sanity check that reconstructing the P3D DSM via
#   p3d_dsm_reconstructed = p3d_chm_mean + p3d_dtm_mean
# gives a plausible surface elevation at GEDI footprints, BEFORE we
# commit to the full chm_alt partition analysis (step 02).
#
# What this script does:
#   1. Loads a handful of per-site enriched CSVs (3 sites by default).
#   2. At 50 random forest footprints per site, computes the
#      reconstructed DSM and sanity-checks it against the 3DEP DTM
#      (which should be close to terrain elevation) and against the
#      site's expected elevation range.
#   3. Writes one CSV (`groundwork_task4_sanity_check.csv`) and prints
#      a concise verdict to the console.
#
# What this script does NOT do:
#   - No full-data processing. Sample size is 50 footprints × 3 sites.
#   - No error metrics. That's step 02.
#   - No figures. Console output + one CSV.
#
# Runtime: ~1–2 minutes.
#
# USAGE (from R in Pane 1, interactively):
#   setwd("/gpfs/data1/vclgp/lmaden/chpt1")  # or wherever the project root is
#   source("coupling_01_dsm_reconstruction_check.R")
#
# After running, inspect:
#   1. The console verdict (PASS/FAIL per check).
#   2. manuscript_tables/groundwork_task4_sanity_check.csv
#
# Expected behavior if the DSM identity is correct:
#   - p3d_dsm_reconstructed should be within a few meters of
#     (dep_dtm_mean + als_chm_mean), i.e. the ALS-derived top-of-canopy
#     elevation, give or take CHM/DTM bias.
#   - p3d_dsm_reconstructed should be numerically plausible (tens to
#     thousands of meters above sea level, matching the site's terrain).
# =====================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# ---- Dependencies: config + utils for paths, log_progress, lookup ----

if (!exists("log_progress")) source("analysis_utils.R")
if (!exists("PROJECT_ROOT")) source("analysis_config.R")

log_progress("=== Step 01: DSM reconstruction sanity check ===")

# ---- Environment sanity: verify we can find everything we need ----
# Script lives in scripts/reviewed/ (relative source() calls above work),
# but all data/outputs live under PROJECT_ROOT. Fail fast if any of the
# expected paths resolve wrong.

cat(sprintf("  PROJECT_ROOT = %s\n",  PROJECT_ROOT))
cat(sprintf("  ENRICHED_DIR = %s\n",  ENRICHED_DIR))
cat(sprintf("  getwd()      = %s\n",  getwd()))

.required_dirs <- c(
  "PROJECT_ROOT" = PROJECT_ROOT,
  "ENRICHED_DIR" = ENRICHED_DIR,
  "manuscript_tables" = file.path(PROJECT_ROOT, "manuscript_tables")
)
.missing <- .required_dirs[!vapply(.required_dirs, dir.exists, logical(1))]
if (length(.missing)) {
  stop("Required directories are missing:\n  ",
       paste(names(.missing), "->", .missing, collapse = "\n  "))
}

# ---- Site selection: pick 3 sites with good coverage ----
#
# We want a mix of: one flagged, one SE-US NEON forest, one grassland.
# These are chosen using TRACKER numbers (which key the enriched CSV
# filenames), then mapped to manuscript numbers for reporting.
#
# Tracker 2 = nasa_howland   -> manuscript Site 2  (FLAGGED)
# Tracker 7 = neon_tall2021  -> manuscript Site 7  (ENF, SE USA)
# Tracker 11 = neon_konz2020 -> manuscript Site 11 (BDF)

sanity_sites_tracker <- c(2L, 7L, 11L)
N_FOOTPRINTS_PER_SITE <- 50L

# ---- Load the site ID lookup ----
# Anchor to PROJECT_ROOT rather than CWD, since this script may be
# source()'d from any working directory.

MS_TABLES_DIR <- file.path(PROJECT_ROOT, "manuscript_tables")
lookup_path <- file.path(MS_TABLES_DIR, "site_id_lookup.csv")
if (!file.exists(lookup_path)) {
  stop("site_id_lookup.csv not found at ", lookup_path,
       "\n  PROJECT_ROOT = ", PROJECT_ROOT,
       "\n  getwd()      = ", getwd(),
       "\n  If PROJECT_ROOT looks wrong, fix the env var and re-source ",
       "analysis_config.R before re-running.")
}
lookup <- data.table::fread(lookup_path)
log_progress(sprintf("Loaded site_id_lookup.csv from %s (%d rows)",
                     lookup_path, nrow(lookup)))

# ---- Read a per-site enriched CSV, keeping only what step 01 needs ----

.read_site_min <- function(site_id_tracker) {
  patterns <- c(
    file.path(ENRICHED_DIR, sprintf("site_%02d_enriched.csv.gz", site_id_tracker)),
    file.path(ENRICHED_DIR, sprintf("site_%d_enriched.csv.gz",   site_id_tracker)),
    file.path(ENRICHED_DIR, sprintf("site_%02d_enriched.csv",    site_id_tracker)),
    file.path(ENRICHED_DIR, sprintf("site_%d_enriched.csv",      site_id_tracker))
  )
  path <- patterns[file.exists(patterns)][1]
  if (is.na(path)) {
    stop("No enriched CSV found for tracker site ", site_id_tracker,
         "\n  Tried: ", paste(patterns, collapse = "\n         "))
  }

  hdr <- names(data.table::fread(path, nrows = 0))

  needed <- c(
    "shot_number", "site",
    "p3d_chm_mean", "p3d_dtm_mean",
    "als_chm_mean", "dep_dtm_mean",
    "als_chm_p90",
    "p3d_chm_valid_frac", "p3d_dtm_valid_frac",
    "als_chm_valid_frac", "dep_dtm_valid_frac"
  )
  have <- intersect(needed, hdr)
  missing <- setdiff(needed, hdr)
  if (length(missing)) {
    log_progress(sprintf(
      "  tracker site %d: CSV is missing columns: %s",
      site_id_tracker, paste(missing, collapse = ", ")
    ))
  }

  DT <- data.table::fread(path, select = have,
                          colClasses = list(character = "shot_number"),
                          nThread = min(4L, CPU_BUDGET),
                          showProgress = FALSE)
  DT[, .tracker_site := site_id_tracker]
  DT[]
}

# ---- Evaluate one site ----

.eval_one <- function(site_id_tracker) {
  ms_row <- lookup[tracker_site == site_id_tracker]
  ms_id  <- if (nrow(ms_row)) ms_row$manuscript_site[1] else NA_integer_
  short  <- if (nrow(ms_row)) ms_row$site_short_name[1] else NA_character_
  flag   <- if (nrow(ms_row)) ms_row$flag_status[1]     else NA_character_

  log_progress(sprintf(
    "--- Tracker site %d (Manuscript Site %s, %s, flag=%s) ---",
    site_id_tracker, as.character(ms_id), short, flag
  ))

  DT <- .read_site_min(site_id_tracker)
  log_progress(sprintf("  read %s rows", format(nrow(DT), big.mark = ",")))

  # Apply the same four-way all-finite + valid-fraction + forest filter
  # that step 02 will use.
  DT_all_finite <- DT[
    is.finite(p3d_chm_mean) &
    is.finite(p3d_dtm_mean) &
    is.finite(als_chm_mean) &
    is.finite(dep_dtm_mean) &
    is.finite(als_chm_p90) &
    is.finite(p3d_chm_valid_frac) & p3d_chm_valid_frac >= 0.5 &
    is.finite(p3d_dtm_valid_frac) & p3d_dtm_valid_frac >= 0.5 &
    is.finite(als_chm_valid_frac) & als_chm_valid_frac >= 0.5 &
    is.finite(dep_dtm_valid_frac) & dep_dtm_valid_frac >= 0.5 &
    als_chm_p90 >= chm_forest_thresh_m
  ]
  log_progress(sprintf(
    "  after all-finite + valid-frac + forest filter (als_chm_p90 >= %g m): %s rows",
    chm_forest_thresh_m, format(nrow(DT_all_finite), big.mark = ",")
  ))

  if (nrow(DT_all_finite) < N_FOOTPRINTS_PER_SITE) {
    log_progress(sprintf(
      "  NOTE: only %d footprints survive filter (< %d requested)",
      nrow(DT_all_finite), N_FOOTPRINTS_PER_SITE
    ))
  }

  n_draw <- min(N_FOOTPRINTS_PER_SITE, nrow(DT_all_finite))
  if (n_draw == 0) {
    log_progress("  WARNING: 0 footprints survive; skipping site")
    return(NULL)
  }

  set.seed(2025L + site_id_tracker)
  S <- DT_all_finite[sample.int(.N, n_draw)]

  # Reconstructed DSM + sanity-comparison surface.
  S[, `:=`(
    p3d_dsm_reconstructed = p3d_chm_mean + p3d_dtm_mean,
    als_top_of_canopy     = dep_dtm_mean + als_chm_mean,
    dsm_minus_als_top     = (p3d_chm_mean + p3d_dtm_mean) -
                            (dep_dtm_mean + als_chm_mean)
  )]

  S[, `:=`(
    tracker_site    = site_id_tracker,
    manuscript_site = ms_id,
    site_short_name = short,
    flag_status     = flag
  )]

  # Per-site console summary.
  dsm   <- S$p3d_dsm_reconstructed
  terr  <- S$dep_dtm_mean
  top   <- S$als_top_of_canopy
  diffs <- S$dsm_minus_als_top

  cat(sprintf(
    paste0(
      "  Summary (n=%d footprints):\n",
      "    p3d_dsm_reconstructed :  min=%7.2f  median=%7.2f  max=%7.2f  (m above sea level)\n",
      "    dep_dtm_mean (terrain):  min=%7.2f  median=%7.2f  max=%7.2f\n",
      "    als_top_of_canopy     :  min=%7.2f  median=%7.2f  max=%7.2f\n",
      "    DSM - ALS_top (diff)  :  min=%7.2f  median=%7.2f  max=%7.2f  sd=%5.2f\n"
    ),
    n_draw,
    min(dsm),  median(dsm),  max(dsm),
    min(terr), median(terr), max(terr),
    min(top),  median(top),  max(top),
    min(diffs), median(diffs), max(diffs), sd(diffs)
  ))

  S[, .(tracker_site, manuscript_site, site_short_name, flag_status,
        shot_number,
        p3d_chm_mean, p3d_dtm_mean,
        als_chm_mean, dep_dtm_mean,
        als_chm_p90,
        p3d_dsm_reconstructed,
        als_top_of_canopy,
        dsm_minus_als_top)]
}

# ---- Run over the three sites ----

all_samples <- rbindlist(lapply(sanity_sites_tracker, .eval_one),
                         fill = TRUE, use.names = TRUE)

# ---- Write the sanity-check CSV ----

dir.create(MS_TABLES_DIR, showWarnings = FALSE, recursive = TRUE)
out_path <- file.path(MS_TABLES_DIR, "groundwork_task4_sanity_check.csv")
data.table::fwrite(all_samples, out_path)
log_progress(sprintf("Wrote %s (%d rows)", out_path, nrow(all_samples)))

# ---- Automated pass/fail checks ----

cat("\n\n====================  VERDICT  ====================\n")

# CHECK 1: reconstructed DSM values are physically plausible.
# (Very loose floor/ceiling: anywhere from sea level to 5000 m.)
chk1_ok <- all(all_samples$p3d_dsm_reconstructed > -50 &
               all_samples$p3d_dsm_reconstructed < 5000)
cat(sprintf("  [%s] 1. Reconstructed DSM in plausible range [-50, 5000] m\n",
            ifelse(chk1_ok, "PASS", "FAIL")))

# CHECK 2: DSM >= DTM at every footprint (surface above ground).
chk2_ok <- all(all_samples$p3d_dsm_reconstructed >= all_samples$dep_dtm_mean - 5)
cat(sprintf("  [%s] 2. Reconstructed DSM >= 3DEP DTM (within 5 m tolerance)\n",
            ifelse(chk2_ok, "PASS", "FAIL")))

# CHECK 3: For each site, most footprints should have
#          |DSM - ALS_top_of_canopy| < 20 m.  Large systematic shifts
#          (tens of meters) would indicate a unit/datum problem.
#          This check is per-site; if ANY site fails, we flag it.
per_site_medabs <- all_samples[, .(
  median_abs_diff = median(abs(dsm_minus_als_top)),
  n               = .N
), by = .(manuscript_site, site_short_name, flag_status)]
chk3_ok <- all(per_site_medabs$median_abs_diff < 20)
cat(sprintf("  [%s] 3. Per-site median |DSM - ALS_top| < 20 m\n",
            ifelse(chk3_ok, "PASS", "FAIL")))
cat("      per-site median |DSM - ALS_top| values:\n")
print(per_site_medabs, row.names = FALSE)

# CHECK 4: Flagged site (nasa_howland) should be the most divergent.
#          This is a positive control on the flag: we EXPECT the
#          flagged site to show a larger |diff| than the non-flagged
#          sites, given the -7.92 m offset documented in Supp Table S13.
if ("FLAGGED" %in% per_site_medabs$flag_status) {
  flag_val <- per_site_medabs[flag_status == "FLAGGED", max(median_abs_diff)]
  ok_val   <- per_site_medabs[flag_status != "FLAGGED", max(median_abs_diff)]
  cat(sprintf(
    "  [INFO] 4. Flagged site median |diff| = %.2f m; worst non-flagged = %.2f m\n",
    flag_val, ok_val
  ))
}

overall_ok <- chk1_ok && chk2_ok && chk3_ok
cat("\n")
if (overall_ok) {
  cat("  ==> Sanity checks PASSED. Ready to proceed to step 02.\n")
} else {
  cat("  ==> Sanity checks FAILED. Inspect CSV and report back before step 02.\n")
}
cat("===================================================\n\n")

log_progress("step 01 complete.")
