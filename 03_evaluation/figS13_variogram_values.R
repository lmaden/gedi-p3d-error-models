# =====================================================================
# figS13_variogram_values.R
#
# Recomputes the five Supp Fig S13
# correlations on the 15-non-flagged frame, alongside the existing
# 18-site frame. Updates the S13 caption and the §4.1 Limitations
# passage that currently cites r = +0.70 (signed-gap vs CHM bias).
#
# Inputs:
#   ${PROJECT_ROOT}/tables/temporal_gap_chm.csv
#     (produced by temporal_mismatch_analysis.R; per-site columns
#      site, gap_years, abs_gap_years, re_intercept, mean_error,
#      mean_abs_error, rmse, etc.)
#
# Outputs:
#   ${PROJECT_ROOT}/manuscript_tables/figS13_variogram_values.csv
#     - long-format table: frame x correlation x (r, p, n)
#
# Notes:
#   The "site" column in temporal_gap_chm.csv is the tracker number,
#   not the manuscript number. Flagged tracker sites = 1, 2, 3.
#   Manuscript sites 1, 2, 3 are the same three sites (Sites 1-15 keep
#   their original numbers; only 16-20 were renumbered after dropping
#   tracker site 16). So filtering tracker {1,2,3} drops the same three
#   sites the manuscript flags.
# =====================================================================

suppressPackageStartupMessages({
  library(dplyr)
})

# ---- 1. Locate input CSV --------------------------------------------

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT",
                           "/gpfs/data1/vclgp/lmaden/chpt1")

csv_in  <- file.path(PROJECT_ROOT, "tables", "temporal_gap_chm.csv")
out_dir <- file.path(PROJECT_ROOT, "manuscript_tables")
csv_out <- file.path(out_dir, "figS13_variogram_values.csv")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(csv_in)) {
  stop(sprintf(paste0(
    "Cannot find %s.\n",
    "If temporal_mismatch_analysis.R was last run with a different\n",
    "out_tables path, set PROJECT_ROOT to the correct value first or\n",
    "edit csv_in directly."
  ), csv_in))
}

d <- read.csv(csv_in, stringsAsFactors = FALSE)
cat(sprintf("Read %s\n  %d rows, %d cols\n", csv_in, nrow(d), ncol(d)))

# ---- 2. Sanity-check columns and frames -----------------------------

required <- c("site", "gap_years", "abs_gap_years",
              "re_intercept", "mean_error", "mean_abs_error", "rmse")
missing  <- setdiff(required, names(d))
if (length(missing) > 0) {
  stop("Missing columns in temporal_gap_chm.csv: ",
       paste(missing, collapse = ", "))
}

flagged_sites <- c(1L, 2L, 3L)  # tracker numbers; see header note

d <- d %>% mutate(site = as.integer(site))
n_total <- nrow(d)
d_15    <- d %>% filter(!site %in% flagged_sites)
n_15    <- nrow(d_15)

cat(sprintf("Sites in CSV:           %d\n", n_total))
cat(sprintf("Flagged sites present:  %s\n",
            paste(intersect(flagged_sites, d$site), collapse = ", ")))
cat(sprintf("Frame after dropping:   %d\n", n_15))

if (n_total != 18L)  warning("Expected 18 sites in the source CSV; got ", n_total)
if (n_15    != 15L)  warning("Expected 15 sites after dropping flagged; got ", n_15)

# ---- 3. Recompute the five correlations on each frame ---------------

# Same five tests as run_correlations() in temporal_mismatch_analysis.R.
# Names match S13/S14 figure-caption usage.

run_one <- function(df) {
  ct <- function(x, y) {
    r <- suppressWarnings(cor.test(x, y))
    list(r = unname(r$estimate), p = r$p.value)
  }

  rows <- list(
    list(test = "abs_gap_vs_re_intercept",
         desc = "|gap| vs CHM site random intercept (panel a)",
         res  = ct(df$abs_gap_years, df$re_intercept)),
    list(test = "abs_gap_vs_mae",
         desc = "|gap| vs site mean |CHM error| (panel b)",
         res  = ct(df$abs_gap_years, df$mean_abs_error)),
    list(test = "abs_gap_vs_bias",
         desc = "|gap| vs site mean CHM bias",
         res  = ct(df$abs_gap_years, df$mean_error)),
    list(test = "signed_gap_vs_bias",
         desc = "Signed gap vs site mean CHM bias (panel c, headline)",
         res  = ct(df$gap_years,     df$mean_error)),
    list(test = "abs_gap_vs_rmse",
         desc = "|gap| vs site CHM RMSE",
         res  = ct(df$abs_gap_years, df$rmse))
  )

  data.frame(
    test = vapply(rows, function(z) z$test, character(1)),
    desc = vapply(rows, function(z) z$desc, character(1)),
    n    = nrow(df),
    r    = vapply(rows, function(z) z$res$r, numeric(1)),
    p    = vapply(rows, function(z) z$res$p, numeric(1)),
    stringsAsFactors = FALSE
  )
}

res18 <- cbind(frame = "18_sites_full",      run_one(d),    stringsAsFactors = FALSE)
res15 <- cbind(frame = "15_non_flagged",     run_one(d_15), stringsAsFactors = FALSE)

result <- rbind(res18, res15)

# ---- 4. Console echo and CSV ----------------------------------------

cat("\n", strrep("=", 70), "\n", sep = "")
cat("S13 RECOMPUTE  ::  18-site full vs 15-non-flagged\n")
cat(strrep("=", 70), "\n", sep = "")

print(result, digits = 3, row.names = FALSE)

write.csv(result, csv_out, row.names = FALSE)
cat(sprintf("\nWrote %s\n", csv_out))

# ---- 5. Summary block: lifted directly to caption + §4.1 ------------

cat("\n", strrep("=", 70), "\n", sep = "")
cat("HEADLINE NUMBERS\n")
cat(strrep("=", 70), "\n", sep = "")

c18 <- result[result$frame == "18_sites_full"      & result$test == "signed_gap_vs_bias", ]
c15 <- result[result$frame == "15_non_flagged"     & result$test == "signed_gap_vs_bias", ]

cat(sprintf("\nPanel (c) signed-gap vs CHM bias:\n"))
cat(sprintf("  18 sites:        r = %+0.3f (p = %.3f)\n", c18$r, c18$p))
cat(sprintf("  15 non-flagged:  r = %+0.3f (p = %.3f)\n", c15$r, c15$p))
cat(sprintf("  delta r:         %+0.3f\n", c15$r - c18$r))

a18 <- result[result$frame == "18_sites_full"  & result$test == "abs_gap_vs_re_intercept", ]
a15 <- result[result$frame == "15_non_flagged" & result$test == "abs_gap_vs_re_intercept", ]
cat(sprintf("\nPanel (a) |gap| vs CHM site RE:\n"))
cat(sprintf("  18 sites:        r = %+0.3f (p = %.3f)\n", a18$r, a18$p))
cat(sprintf("  15 non-flagged:  r = %+0.3f (p = %.3f)\n", a15$r, a15$p))

b18 <- result[result$frame == "18_sites_full"  & result$test == "abs_gap_vs_mae", ]
b15 <- result[result$frame == "15_non_flagged" & result$test == "abs_gap_vs_mae", ]
cat(sprintf("\nPanel (b) |gap| vs MAE:\n"))
cat(sprintf("  18 sites:        r = %+0.3f (p = %.3f)\n", b18$r, b18$p))
cat(sprintf("  15 non-flagged:  r = %+0.3f (p = %.3f)\n", b15$r, b15$p))

cat("\nS13 recompute complete.\n")
