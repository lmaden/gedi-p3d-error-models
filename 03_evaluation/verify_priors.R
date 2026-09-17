# =====================================================================
# verify_priors.R
#
# Phase S2-C verification: does Table S3 of supplemental_figures_v09 give
# the COMPLETE prior specification for Eqs. S1-S4, or is there a gap for
# the variance-model random-effect SDs (v_site, v_lc in Eq. S4)?
#
# This script loads both canonical brms fits and prints:
#   (1) The full prior_summary() — every prior brms is using, set + default
#   (2) The model formula (to confirm v_site/v_lc submodel parameterization)
#   (3) The parameter classes brms exposes (to confirm which SDs exist)
#   (4) A focused diff: u-side vs v-side random-effect SD priors
#
# Run inside Pane 1 (the R console started by tmux-quad-chpt1-v3.sh).
# Wall time: ~30 seconds (load + print only; no refit).
#
# =====================================================================

suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
})

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/gpfs/data1/vclgp/lmaden/chpt1")
ckpt <- function(p) readRDS(file.path(PROJECT_ROOT, "checkpoints", p))

cat("============================================================\n")
cat("Phase S2-C prior-spec verification\n")
cat("PROJECT_ROOT =", PROJECT_ROOT, "\n")
cat("============================================================\n\n")

# ---- Load fits ----
cat("[1/4] Loading checkpoints...\n")
cp_chm <- ckpt("10b_chm_sensitivity.rds")
cp_dtm <- ckpt("10_models_stage2.rds")
fit_chm <- cp_chm$data$fit_chm_16
fit_dtm <- cp_dtm$data$fit_dtm_s2
cat("  fit_chm class:", paste(class(fit_chm), collapse=","), "\n")
cat("  fit_dtm class:", paste(class(fit_dtm), collapse=","), "\n\n")

# ---- Formulas (canonical record of mean + sigma submodels) ----
cat("[2/4] Model formulas — confirms v_site / v_lc are in the sigma submodel:\n\n")
cat("--- CHM 16-site formula ---\n")
print(formula(fit_chm))
cat("\n--- DTM 18-site formula ---\n")
print(formula(fit_dtm))
cat("\n")

# ---- Full prior_summary (the ground truth: every prior, including defaults) ----
cat("[3/4] FULL prior_summary() — every prior in use (user-set OR brms default):\n\n")
cat("--- CHM 16-site priors ---\n")
ps_chm <- prior_summary(fit_chm)
print(ps_chm)
cat("\n--- DTM 18-site priors ---\n")
ps_dtm <- prior_summary(fit_dtm)
print(ps_dtm)
cat("\n")

# ---- Focused diff: u-side vs v-side random-effect SDs ----
# Random-effect SDs in brms have class="sd"; the sigma-submodel ones have dpar="sigma".
# Mean-model ones have dpar="" (empty / NA).
cat("[4/4] Focused split — random-effect SDs by submodel (u-side vs v-side):\n\n")

split_sd_priors <- function(ps, model_label) {
  cat(sprintf("--- %s ---\n", model_label))
  sd_rows <- ps[ps$class == "sd", , drop = FALSE]
  if (nrow(sd_rows) == 0) {
    cat("  (no class='sd' rows found — unexpected)\n")
    return(invisible(NULL))
  }
  # dpar column exists in brms >= 2.x; may be NA or "" for mean-model SDs
  if (!"dpar" %in% names(sd_rows)) {
    cat("  (no 'dpar' column — old brms; cannot split u/v automatically)\n")
    print(sd_rows)
    return(invisible(NULL))
  }
  dpar_norm <- ifelse(is.na(sd_rows$dpar) | sd_rows$dpar == "", "MEAN_MODEL_u", sd_rows$dpar)
  cat("\n  u-side SDs (mean model — Eqs. S3a/b/c — random intercepts/slopes on mu):\n")
  u_rows <- sd_rows[dpar_norm == "MEAN_MODEL_u", , drop = FALSE]
  print(u_rows[, intersect(names(u_rows), c("prior","class","coef","group","resp","dpar","source")), drop=FALSE])
  cat("\n  v-side SDs (sigma submodel — Eq. S4 — v_site, v_lc random effects on log-sigma):\n")
  v_rows <- sd_rows[dpar_norm == "sigma", , drop = FALSE]
  if (nrow(v_rows) == 0) {
    cat("  (NONE — sigma submodel has no random effects; v_site/v_lc may be using\n")
    cat("   a different parameterization, OR the model object differs from Eq. S4)\n")
  } else {
    print(v_rows[, intersect(names(v_rows), c("prior","class","coef","group","resp","dpar","source")), drop=FALSE])
  }
  invisible(NULL)
}

split_sd_priors(ps_chm, "CHM 16-site")
cat("\n")
split_sd_priors(ps_dtm, "DTM 18-site")

cat("\n============================================================\n")
cat("END Phase S2-C verification.\n")
cat("\n")
cat("HOW TO INTERPRET:\n")
cat("Path A — v-side priors are EXPLICITLY SET to the same Student-t(3, 0, 2)\n")
cat("         as the u-side priors (source = 'user' in prior_summary):\n")
cat("         → Table S3's row 4 legitimately covers both. The 'complete prior\n")
cat("           specification' claim in para 93 is accurate; we propose a small\n")
cat("           clarifying parenthetical to row 4's Parameter label.\n\n")
cat("Path B — v-side priors are brms DEFAULT (source = 'default' in prior_summary,\n")
cat("         typically student_t(3, 0, 2.5)):\n")
cat("         → Table S3 understates: either set explicit priors to match what\n")
cat("           the supplement claims, or add a row to Table S3 documenting the\n")
cat("           brms-default used for v-side SDs.\n\n")
cat("Path C — v-side priors are GENUINELY DIFFERENT from u-side (any source):\n")
cat("         → Real gap. Add a separate row to Table S3 for v-side SDs.\n")
cat("============================================================\n")
