#!/usr/bin/env Rscript
# =====================================================================
# loo_decomposition_audit.R  --  read-only
#
# Answers M3's "Required" in one run, and settles which fit is the
# manuscript's counterfactual.
#
# The reviewer (M3) asks for the decomposition stated explicitly:
#   operational coefficient, DTM coefficient, their sum, and the
#   fitted counterfactual -- plus why the sum does not close.
# Manuscript v113 Fig. 4d reports -0.44 -> +0.39; the review says the
# naive sum is +0.55 and the counterfactual CI excludes it.
#
# This script prints the four numbers from their authoritative fits,
# with N and site count attached to each so the "different site set,
# different sample" explanation is evidenced rather than asserted.
#
# It also tests one manuscript claim directly: v113 section 2.3 states the
# counterfactual's specification is "identical to the 16-site CHM
# specification". The refit script pulls chm_formula_s2 from
# 10_models_stage2 (the 19-site fit), so the formulas are compared.
#
# WRITES NOTHING. Loads one checkpoint at a time and frees it, so it is
# safe to run while the M2 DTM refit is using 64 cores on gsapp22.
#
# USAGE (exports first; from any node):
#   export PROJECT_ROOT=/gpfs/data1/vclgp/lmaden/chpt1
#   export R_LIBS=/gpfs/data1/vclgp/lmaden/Rlib
#   nice -n 15 Rscript $PROJECT_ROOT/scripts/reviewed/loo_decomposition_audit.R \
#     2>&1 | tee $PROJECT_ROOT/manuscript_tables/loo_decomposition_audit.txt
# =====================================================================

suppressPackageStartupMessages(library(brms))

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/gpfs/data1/vclgp/lmaden/chpt1")
CKPT <- file.path(PROJECT_ROOT, "checkpoints")

hr  <- function() cat(strrep("=", 78), "\n")
sec <- function(x) { cat("\n"); hr(); cat("## ", x, "\n"); hr() }

unwrap <- function(o) if (is.list(o) && !is.null(o$data)) o$data else o

read_ckpt <- function(name) {
  p <- file.path(CKPT, paste0(name, ".rds"))
  if (!file.exists(p)) stop("checkpoint not found: ", p)
  cat("  reading ", basename(p), " ...\n", sep = "")
  unwrap(readRDS(p))
}

# Pull the RH98 fixed effect out of a brmsfit, whatever the row is called.
rh98_row <- function(fit) {
  fe <- fixef(fit)
  ix <- grep("^rh_?98", rownames(fe))
  if (length(ix) != 1L)
    stop("expected exactly one RH98 row, found ", length(ix), ": ",
         paste(rownames(fe)[ix], collapse = ", "))
  fe[ix, , drop = FALSE]
}

describe <- function(fit, key, label) {
  d  <- fit$data
  ns <- length(unique(as.character(d$site)))
  fe <- rh98_row(fit)
  sd_site <- tryCatch(round(VarCorr(fit)$site$sd["Intercept", "Estimate"], 4),
                      error = function(e) NA_real_)
  cat(sprintf("\n  %-38s  N = %s | sites = %d | sd_site = %s\n",
              label, format(nrow(d), big.mark = ","), ns, format(sd_site)))
  cat(sprintf("    RH98 (%s): %+.4f  [%+.4f, %+.4f]\n",
              rownames(fe)[1], fe[1, "Estimate"], fe[1, "Q2.5"], fe[1, "Q97.5"]))
  list(key = key, label = label, n = nrow(d), n_sites = ns,
       est = fe[1, "Estimate"], lo = fe[1, "Q2.5"], hi = fe[1, "Q97.5"])
}

cat("loo_decomposition_audit.R -- read-only\n")
cat("started: ", format(Sys.time()), "\n", sep = "")
cat("PROJECT_ROOT: ", PROJECT_ROOT, "\n", sep = "")

res <- list()

# ---------------------------------------------------------------------
sec("1 -- OPERATIONAL: the primary 16-site CHM fit (expect -0.44)")
# Per checkpoint v6 sec 6.4 the PRIMARY 16-site CHM fit lives in the file
# named "sensitivity". That naming hazard is why this is stated here.
ck <- read_ckpt("10b_chm_sensitivity")
cat("  objects: ", paste(names(ck), collapse = ", "), "\n", sep = "")
res$chm16 <- describe(ck$fit_chm_16, "chm16", "16-site CHM (operational)")
formula_chm16 <- ck$fit_chm_16$formula
rm(ck); gc(verbose = FALSE)

# ---------------------------------------------------------------------
sec("2 -- DTM: the 18-site terrain fit (expect +0.99), and the 19-site CHM")
ck <- read_ckpt("10_models_stage2")
cat("  objects: ", paste(names(ck), collapse = ", "), "\n", sep = "")
res$dtm18 <- describe(ck$fit_dtm_s2, "dtm18", "18-site DTM")
if (!is.null(ck$fit_chm_s2))
  res$chm19 <- describe(ck$fit_chm_s2, "chm19", "19-site CHM (NOT primary)")
formula_s2 <- ck$chm_formula_s2
rm(ck); gc(verbose = FALSE)

# ---------------------------------------------------------------------
sec("3 -- COUNTERFACTUAL: the err_alt refits (expect +0.39 on 15 sites)")
ck <- read_ckpt("groundwork_task4_refit")
cat("  objects: ", paste(names(ck), collapse = ", "), "\n", sep = "")
res$alt15 <- describe(ck$fit_alt_15, "alt15", "15-site counterfactual (err_alt)")
res$alt18 <- describe(ck$fit_alt_18, "alt18", "18-site counterfactual (err_alt)")
rm(ck); gc(verbose = FALSE)

# ---------------------------------------------------------------------
sec("4 -- THE M3 DECOMPOSITION")

op  <- res$chm16$est
dtm <- res$dtm18$est
cf  <- res$alt15$est

cat(sprintf("\n  operational  (16-site CHM, RH98)      %+.4f  [%+.4f, %+.4f]\n",
            op, res$chm16$lo, res$chm16$hi))
cat(sprintf("  DTM          (18-site DTM,  RH98)      %+.4f  [%+.4f, %+.4f]\n",
            dtm, res$dtm18$lo, res$dtm18$hi))
cat(sprintf("  naive sum    (operational + DTM)       %+.4f\n", op + dtm))
cat(sprintf("  fitted counterfactual (15-site)        %+.4f  [%+.4f, %+.4f]\n",
            cf, res$alt15$lo, res$alt15$hi))
cat(sprintf("\n  gap (naive sum - counterfactual)       %+.4f\n", (op + dtm) - cf))
cat(sprintf("  sum inside the counterfactual CI?      %s\n",
            if ((op + dtm) >= res$alt15$lo && (op + dtm) <= res$alt15$hi)
              "YES" else "NO -- this is the reviewer's point"))

cat("\n  Why the three numbers cannot be expected to close:\n")
cat(sprintf("    operational   N = %-9s sites = %d\n",
            format(res$chm16$n, big.mark = ","), res$chm16$n_sites))
cat(sprintf("    DTM           N = %-9s sites = %d\n",
            format(res$dtm18$n, big.mark = ","), res$dtm18$n_sites))
cat(sprintf("    counterfactual N = %-9s sites = %d\n",
            format(res$alt15$n, big.mark = ","), res$alt15$n_sites))
cat("    Three different site sets and three different samples; the DTM fit\n")
cat("    additionally carries its own land-cover reference level. The identity\n")
cat("    err_alt = err_p3d + err_dtm is exact per footprint, but coefficient\n")
cat("    additivity across three separately-fitted hierarchical models is not.\n")

# ---------------------------------------------------------------------
sec("5 -- SPEC CHECK: is the counterfactual's formula the 16-site one?")
# v113 section 2.3 claims the counterfactual specification is identical to the
# 16-site CHM. counterfactual_chm_refit.R takes chm_formula_s2 from
# 10_models_stage2 (the 19-site fit). Compare them literally.
f_a <- paste(deparse(formula_chm16), collapse = " ")
f_b <- paste(deparse(formula_s2),    collapse = " ")
f_a <- gsub("[[:space:]]+", " ", f_a); f_b <- gsub("[[:space:]]+", " ", f_b)
cat("\n  fit_chm_16$formula:\n    ", f_a, "\n", sep = "")
cat("\n  chm_formula_s2 (used by the counterfactual refit):\n    ", f_b, "\n", sep = "")
cat(sprintf("\n  IDENTICAL: %s\n", identical(f_a, f_b)))
cat("  (If FALSE, v113 section 2.3's \"identical to the 16-site CHM\n")
cat("   specification\" needs rewording -- report the difference, do not refit.)\n")

sec("DONE -- nothing was written except this log")
cat("finished: ", format(Sys.time()), "\n", sep = "")
