# =====================================================================
# figS9_16site_group_values.R
#
# Extracts 16-site CHM sensitivity-fit
# random-effect summaries at the ecoregion and lc_l1_code grouping
# levels, mirroring the 19-site rows already in Supp Table S9a (rows
# R07-R16 of the existing table).
#
# Inputs:
#   ${PROJECT_ROOT}/checkpoints/10b_chm_sensitivity.rds
#     -> $fit_chm_16   (16-site CHM sensitivity refit; brms object)
#
# Outputs:
#   ${PROJECT_ROOT}/manuscript_tables/figS9_16site_group_values.csv
#     - rows ready to drop into Supp Table S9a as
#       "CHM (16-site sensitivity)" rows for Ecoregion and Land cover.
#
# What this produces:
#   For the 16-site CHM fit, summarize the same RE parameters that
#   the 19-site rows summarize:
#     Ecoregion  | Intercept              (mirrors 19-site R07)
#     Ecoregion  | Slope (slope_mean_z)   (mirrors 19-site R08)
#     Land cover | Intercept              (mirrors 19-site R11)
#     Land cover | Slope (wsci_z)         (mirrors 19-site R13)
#     Land cover | sigma Intercept        (mirrors 19-site R12,
#                                          if the sigma submodel has
#                                          a v_lc term in this fit)
#
# Numerics: posterior medians at each level, n_levels, min, max,
# range, SD across levels (matching the existing R01-R16 statistics).
# =====================================================================

suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
})

# ---- 1. Locate inputs ----------------------------------------------

if (!exists("PROJECT_ROOT") || !exists("log_progress")) {
  cfg_candidates <- c(
    file.path(getwd(), "analysis_config.R"),
    "/gpfs/data1/vclgp/lmaden/chpt1/scripts/reviewed/analysis_config.R"
  )
  util_candidates <- c(
    file.path(getwd(), "analysis_utils.R"),
    "/gpfs/data1/vclgp/lmaden/chpt1/scripts/reviewed/analysis_utils.R"
  )
  cfg  <- cfg_candidates[file.exists(cfg_candidates)][1]
  util <- util_candidates[file.exists(util_candidates)][1]
  if (is.na(cfg) || is.na(util)) {
    stop("Could not locate analysis_config.R / analysis_utils.R")
  }
  source(cfg)
  source(util)
}

cp_path <- file.path(CHECKPOINT_DIR, "10b_chm_sensitivity.rds")
out_dir <- file.path(PROJECT_ROOT, "manuscript_tables")
csv_out <- file.path(out_dir, "figS9_16site_group_values.csv")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(cp_path)) {
  stop(sprintf("Checkpoint not found: %s\n  Need fit_chm_16 from 16-site.",
               cp_path))
}

log_progress(sprintf("Loading %s...", cp_path))
cp <- readRDS(cp_path)

# Several wrapping conventions exist depending on how it was saved.
fit_chm_16 <- cp$fit_chm_16
if (is.null(fit_chm_16) && !is.null(cp$data)) {
  fit_chm_16 <- cp$data$fit_chm_16
}
if (is.null(fit_chm_16)) {
  stop(sprintf(
    "fit_chm_16 not found in %s (top-level keys: %s)",
    cp_path, paste(names(cp), collapse = ", ")))
}
log_progress(sprintf("  Loaded fit_chm_16 (n = %d obs, %d sites)",
                     nrow(fit_chm_16$data),
                     length(unique(as.character(fit_chm_16$data$site)))))

# ---- 2. Pull random effects at ecoregion and lc_l1_code levels -----
#
# ranef() returns a list; one entry per grouping factor in the model.
# Each entry is a 3-D array: levels x summary x parameters
#   summary = c("Estimate", "Est.Error", "Q2.5", "Q97.5")  (default)
# We use the posterior median ("Estimate" with default brms is the
# posterior mean, not median; for consistency with the existing S9a
# row construction we explicitly request median here).

re_all <- ranef(fit_chm_16, robust = TRUE,
                probs = c(0.025, 0.975))

log_progress(sprintf("  Random-effect grouping factors in fit: %s",
                     paste(names(re_all), collapse = ", ")))

# ---- 3. Helper: compute n_levels / min / max / range / SD for a
#                 named (group, parameter) cell.
# --------------------------------------------------------------------

summarize_re <- function(re_arr, grouping, parameter, label) {
  if (!grouping %in% names(re_arr)) {
    return(data.frame(
      product = "CHM (16-site sensitivity)",
      grouping_level = label$grouping_label,
      parameter = label$parameter_label,
      n_levels = NA_integer_,
      min_m = NA_real_, max_m = NA_real_,
      range_m = NA_real_, sd_m = NA_real_,
      note = sprintf("grouping factor '%s' not present in fit", grouping),
      stringsAsFactors = FALSE
    ))
  }
  arr <- re_arr[[grouping]]
  if (!parameter %in% dimnames(arr)[[3]]) {
    return(data.frame(
      product = "CHM (16-site sensitivity)",
      grouping_level = label$grouping_label,
      parameter = label$parameter_label,
      n_levels = dim(arr)[1],
      min_m = NA_real_, max_m = NA_real_,
      range_m = NA_real_, sd_m = NA_real_,
      note = sprintf("parameter '%s' not present in %s effects",
                     parameter, grouping),
      stringsAsFactors = FALSE
    ))
  }
  est <- arr[, "Estimate", parameter]
  data.frame(
    product = "CHM (16-site sensitivity)",
    grouping_level = label$grouping_label,
    parameter = label$parameter_label,
    n_levels = length(est),
    min_m = min(est),
    max_m = max(est),
    range_m = max(est) - min(est),
    sd_m = sd(est),
    note = "",
    stringsAsFactors = FALSE
  )
}

# ---- 4. Build the new S9a rows -------------------------------------
#
# Mirror the 19-site rows R07, R08, R11, R12, R13 from Table S9a.
# The 19-site sigma_lc row (R12) is "Land cover | sigma Intercept";
# in brms the corresponding ranef arr has parameter name
# "sigma_Intercept" when there's a v_lc term in the sigma submodel.

rows <- list(
  # Ecoregion
  list(grouping = "ecoregion", parameter = "Intercept",
       label = list(grouping_label = "Ecoregion",
                    parameter_label = "Intercept")),
  list(grouping = "ecoregion", parameter = "slope_mean_z",
       label = list(grouping_label = "Ecoregion",
                    parameter_label = "Slope (slope_mean_z)")),
  # Land cover
  list(grouping = "lc_l1_code", parameter = "Intercept",
       label = list(grouping_label = "Land cover",
                    parameter_label = "Intercept")),
  list(grouping = "lc_l1_code", parameter = "sigma_Intercept",
       label = list(grouping_label = "Land cover",
                    parameter_label = "sigma Intercept")),
  list(grouping = "lc_l1_code", parameter = "wsci_z",
       label = list(grouping_label = "Land cover",
                    parameter_label = "Slope (wsci_z)"))
)

S9a_rows <- do.call(rbind,
  lapply(rows, function(r) {
    summarize_re(re_all, r$grouping, r$parameter, r$label)
  })
)

# Print before saving so the user sees them in the console.
cat("\n", strrep("=", 70), "\n", sep = "")
cat("16-site sensitivity CHM: random-effect summaries (ecoregion + lc)\n")
cat(strrep("=", 70), "\n", sep = "")
print(S9a_rows, digits = 3, row.names = FALSE)

write.csv(S9a_rows, csv_out, row.names = FALSE)
log_progress(sprintf("Wrote %s", csv_out))

# ---- 5. Sanity-check echo against 19-site existing rows ------------
#
# Just so it's easy to eyeball the CSV next to the v02 supplement
# without flipping documents, print the 19-site rows that the new
# rows should be compared against.

cat("\n", strrep("=", 70), "\n", sep = "")
cat("Reference: existing 19-site rows in Supp Table S9a\n",
    "(from supplemental_figures_v02.docx, Table 8)\n", sep = "")
cat(strrep("=", 70), "\n", sep = "")
ref <- data.frame(
  product = c("CHM", "CHM", "CHM", "CHM", "CHM"),
  grouping_level = c("Ecoregion", "Ecoregion", "Land cover", "Land cover",
                     "Land cover"),
  parameter = c("Intercept", "Slope (slope_mean_z)",
                "Intercept", "sigma Intercept", "Slope (wsci_z)"),
  n_levels = c(11, 11, 16, 16, 16),
  min_m = c(-0.24, -0.34, -0.54, -0.29, -0.22),
  max_m = c(+0.31, +0.18, +0.28, +0.14, +0.06),
  range_m = c(0.55, 0.52, 0.82, 0.43, 0.27),
  sd_m = c(0.16, 0.17, 0.21, 0.14, 0.07)
)
print(ref, row.names = FALSE)

cat("\n", strrep("=", 70), "\n", sep = "")
cat("Drop new rows directly under existing 19-site rows in Table S9a.\n", sep = "")
cat("Order them: ecoregion intercept, ecoregion slope,\n",
    "lc intercept, lc sigma intercept, lc slope.\n", sep = "")
cat(strrep("=", 70), "\n\n", sep = "")
log_progress("Extraction complete.")
