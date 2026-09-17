# s3_vif_corr_probdir_16site.R  (14 Sep 2026)
# Read-only. Exports, from the published 16-site CHM fit (checkpoints/10b_chm_sensitivity.rds, fit_chm_16):
#   manuscript_tables/s3_corr_matrix_chm16.csv        17 x 17 Pearson correlations of the z-scored predictors (fit data)
#   manuscript_tables/s3_vif_chm16.csv                variance inflation factors of the same 17 predictors
#   manuscript_tables/tableS1_chm16_interactions.csv  posterior mean, median, 95% CI and P(direction) for every
#                                                     interaction coefficient (Table S1 / S4 rebuild)
# Run from the chapter root:  Rscript scripts/reviewed/s3_vif_corr_probdir_16site.R
suppressPackageStartupMessages({ library(brms); library(posterior) })
root <- Sys.getenv("CHPT1_ROOT", "/gpfs/data1/vclgp/lmaden/chpt1")
ck   <- file.path(root, "checkpoints", "10b_chm_sensitivity.rds")
stopifnot(file.exists(ck))
obj <- readRDS(ck)
fit <- if (!is.null(obj$data) && !is.null(obj$data$fit_chm_16)) obj$data$fit_chm_16 else obj$fit_chm_16
stopifnot(inherits(fit, "brmsfit"))
d <- fit$data
preds <- c("slope_mean_z","slope_sd_z","aspect_sin_z","aspect_cos_z","wsci_z","rh_98_z","cover_z",
           "meta_offnad_z","meta_sunel_z","meta_az_conc_z","view_az_sin_z","view_az_cos_z",
           "meta_stereo_z","meta_fwdrev_z","meta_absgeo_z","meta_relgeo_z","meta_leafon_z")
miss <- setdiff(preds, names(d)); if (length(miss)) stop("predictor columns missing from fit$data: ", paste(miss, collapse = ", "))
X <- as.matrix(d[, preds]); X <- X[stats::complete.cases(X), ]
R <- stats::cor(X); vif <- diag(solve(R))
out <- file.path(root, "manuscript_tables"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
utils::write.csv(round(R, 4), file.path(out, "s3_corr_matrix_chm16.csv"))
utils::write.csv(data.frame(predictor = preds, vif = round(vif, 3)), file.path(out, "s3_vif_chm16.csv"), row.names = FALSE)
cat(sprintf("fit data: n = %d footprints, %d sites | max VIF = %.3f (%s) | all VIF < 10: %s\n",
            nrow(X), length(unique(d$site)), max(vif), preds[which.max(vif)], all(vif < 10)))
dr <- posterior::as_draws_matrix(fit, variable = "^b_", regex = TRUE)
cn <- colnames(dr); ix <- grepl(":", cn, fixed = TRUE)
M  <- dr[, ix, drop = FALSE]
tab <- data.frame(term = sub("^b_", "", cn[ix]),
                  mean = colMeans(M), median = apply(M, 2, stats::median),
                  lo95 = apply(M, 2, stats::quantile, 0.025), hi95 = apply(M, 2, stats::quantile, 0.975),
                  p_positive = colMeans(M > 0), row.names = NULL)
utils::write.csv(tab, file.path(out, "tableS1_chm16_interactions.csv"), row.names = FALSE)
cat(sprintf("wrote %d interaction rows; credible (95%% CI excludes 0): %d\n", nrow(tab), sum(tab$lo95 > 0 | tab$hi95 < 0)))
