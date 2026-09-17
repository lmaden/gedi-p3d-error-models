suppressPackageStartupMessages(library(brms))
CHKPT  <- file.path(Sys.getenv("PROJECT_ROOT"), "checkpoints")
unwrap <- function(p) if (is.list(p) && !is.null(p$data)) p$data else p
cat("loading 10_models_stage2.rds (198 MB), please wait ...\n")
fit <- unwrap(readRDS(file.path(CHKPT, "10_models_stage2.rds")))$fit_dtm_s2
gc(verbose = FALSE)

cat("\n===== 1. FORMULA =====\n"); print(fit$formula)
cat("\n--- mu formula, one line ---\n")
cat(paste(deparse(fit$formula$formula, width.cutoff = 500L), collapse = " "), "\n")
cat("\n--- distributional parts ---\n"); print(fit$formula$pforms)

cat("\n===== 2. FAMILY =====\n");  print(fit$family)
cat("\n===== 3. PRIORS =====\n");  print(fit$prior)

cat("\n===== 4. WEIGHTS =====\n")
cat("'weights' appears in formula:",
    grepl("weights", paste(deparse(fit$formula$formula), collapse = " ")), "\n")
cat("w_dtm retained in fit$data:", "w_dtm" %in% names(fit$data), "\n")

cat("\n===== 5. SIZE =====\n")
cat("nobs:", nobs(fit), " ndraws:", brms::ndraws(fit),
    " nchains:", brms::nchains(fit), "\n")
cat("columns brms kept:\n"); print(names(fit$data))

cat("\n===== 6. PUBLISHED SAMPLER BEHAVIOUR =====\n")
np <- brms::nuts_params(fit); td <- np$Value[np$Parameter == "treedepth__"]
cat("divergences:", sum(np$Value[np$Parameter == "divergent__"]), "\n")
cat("mean treedepth:", round(mean(td), 3), " max:", max(td), "\n")
cat("pct >= 12:", round(100 * mean(td >= 12), 3),
    "%   pct >= 13:", round(100 * mean(td >= 13), 3), "%\n")
cat("stepsizes:", paste(signif(unique(np$Value[np$Parameter == "stepsize__"]), 4),
    collapse = ", "), "\n")
et <- tryCatch(rstan::get_elapsed_time(fit$fit), error = function(e) NULL)
if (!is.null(et)) { cat("\npublished elapsed hours per chain:\n")
  print(round(cbind(et / 3600, total = rowSums(et) / 3600), 2)) }

cat("\n===== 7. FIXED EFFECTS =====\n"); print(rownames(brms::fixef(fit)))
