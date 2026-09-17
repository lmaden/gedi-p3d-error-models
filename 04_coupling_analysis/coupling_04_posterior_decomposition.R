## coupling_04_posterior_decomposition.R
##
## Posterior-based decomposition of the four CHM error fits, plus diagnostics
## needed to write the Ch1 Discussion paragraph.
##
## Inputs (existing checkpoints; no extraction or refit needed):
##   10_models_stage2          -> fit_chm_s2 (Ch1 full, 19 sites, err_p3d)
##   10b_chm_sensitivity       -> fit_chm_16 (16-site,    16 sites, err_p3d)
##   groundwork_task4_refit    -> fit_alt_18, fit_alt_15  (step 04, err_alt)
##
## Outputs (PROJECT_ROOT/manuscript_tables/):
##   coupling_variance_summary.csv   -- four-fit posterior summary
##   coupling_variance_decomposition.csv    -- decomposition with 95% CIs
##   coupling_r2_partition.csv       -- Bayesian R^2 + FE/RE/site shares
##   coupling_site10_diagnostic.csv  -- bounds the 16-vs-15 confound
##   coupling_random_effect_shifts.csv          -- per-site RE Ch1 vs step 04
## Plus: plots/groundwork/coupling_04_posterior_decomposition_summary.pdf
##
## Usage: source("coupling_04_posterior_decomposition.R") from Pane 1.
## Runtime: a few minutes per fit for posterior_predict; ~10-30 min total.
## ============================================================================

suppressPackageStartupMessages({
  library(brms)
  library(posterior)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

source("analysis_config.R")
source("analysis_utils.R")

stopifnot(nzchar(PROJECT_ROOT))
TBL_DIR  <- file.path(PROJECT_ROOT, "manuscript_tables")
PLOT_DIR <- file.path(PROJECT_ROOT, "plots", "groundwork")
dir.create(TBL_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

log_msg <- function(...) {
  message(sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                  paste0(..., collapse = "")))
}
banner <- function(s) {
  log_msg(paste(rep("=", 72), collapse = ""))
  log_msg(s)
  log_msg(paste(rep("=", 72), collapse = ""))
}

banner("Step 04 -- posterior decomposition")
log_msg("PROJECT_ROOT: ", PROJECT_ROOT)

## ----------------------------------------------------------------------------
## Section 1 -- Load checkpoints
## ----------------------------------------------------------------------------
banner("1. Load checkpoints")

read_ckpt <- function(name) {
  p <- file.path(PROJECT_ROOT, "checkpoints", paste0(name, ".rds"))
  if (!file.exists(p)) stop("checkpoint not found: ", p)
  obj <- readRDS(p)
  if (is.list(obj) && !is.null(obj$data)) return(obj$data)
  obj
}

ck_full <- read_ckpt("10_models_stage2")
ck_trkB <- read_ckpt("10b_chm_sensitivity")
ck_p3   <- read_ckpt("groundwork_task4_refit")

fits <- list(
  full   = ck_full$fit_chm_s2,
  trackB = ck_trkB$fit_chm_16,
  alt18  = ck_p3$fit_alt_18,
  alt15  = ck_p3$fit_alt_15
)

fit_meta <- data.table(
  fit_key  = c("full", "trackB", "alt18", "alt15"),
  label    = c("Ch1 full (19, err_p3d)",
               "16-site (16, err_p3d)",
               "step 04 (18, err_alt)",
               "step 04 (15, err_alt)"),
  response = c("err_p3d", "err_p3d", "err_alt", "err_alt"),
  n_sites  = c(19L, 16L, 18L, 15L)
)

for (k in names(fits)) {
  if (is.null(fits[[k]])) stop("fit '", k, "' not found in checkpoint")
  log_msg("  ", k, ": brmsfit, ",
          nrow(fits[[k]]$data), " obs, ",
          length(unique(fits[[k]]$data$site)), " sites")
}

## ----------------------------------------------------------------------------
## Section 2 -- Posterior sigma_site and V_site for each fit
## ----------------------------------------------------------------------------
banner("2. Posterior sigma_site and V_site")

# Extract posterior of sd_site (the random-intercept SD parameter for `site`)
# via as_draws_df, which is more stable across brms versions than VarCorr().
get_sd_site <- function(fit) {
  d <- as_draws_df(fit)
  if ("sd_site__Intercept" %in% names(d)) return(as.numeric(d$sd_site__Intercept))
  nms <- grep("^sd_site", names(d), value = TRUE)
  if (length(nms) == 0L) stop("no sd_site* parameter found in fit")
  log_msg("    using fallback sd parameter: ", nms[1])
  as.numeric(d[[nms[1]]])
}

draws_sd <- lapply(fits, get_sd_site)
draws_V  <- lapply(draws_sd, function(x) x^2)

summarise_post <- function(x) {
  c(mean   = mean(x),
    median = median(x),
    lo     = quantile(x, 0.025, names = FALSE),
    hi     = quantile(x, 0.975, names = FALSE))
}

var_summary <- data.table(
  fit_key  = names(draws_sd),
  label    = fit_meta$label,
  response = fit_meta$response,
  n_sites  = fit_meta$n_sites
)
var_summary[, c("sd_mean", "sd_median", "sd_lo", "sd_hi") :=
              as.list(rep(NA_real_, 4))]
var_summary[, c("V_mean", "V_median", "V_lo", "V_hi") :=
              as.list(rep(NA_real_, 4))]

for (i in seq_along(draws_sd)) {
  s <- summarise_post(draws_sd[[i]])
  v <- summarise_post(draws_V[[i]])
  set(var_summary, i, c("sd_mean","sd_median","sd_lo","sd_hi"),
      as.list(unname(s)))
  set(var_summary, i, c("V_mean","V_median","V_lo","V_hi"),
      as.list(unname(v)))
}

fwrite(var_summary,
       file.path(TBL_DIR, "coupling_variance_summary.csv"))
log_msg("Wrote variance_summary.csv")
print(var_summary[, .(label, n_sites,
                      sd_mean, sd_lo, sd_hi,
                      V_mean,  V_lo,  V_hi)])

## ----------------------------------------------------------------------------
## Section 3 -- Variance decomposition (with posterior CIs)
## ----------------------------------------------------------------------------
banner("3. Variance decomposition")

# All four chains have the same number of post-warmup draws (1500 * 4 = 6000)
# but check before doing draw-wise differences.
n_draws <- sapply(draws_V, length)
if (length(unique(n_draws)) != 1L)
  stop("posterior draw counts differ across fits: ",
       paste(n_draws, collapse = ", "))

# Component contributions, per draw, in m^2:
#   ref_quality_contrib = V_full - V_trackB           (clean, both err_p3d)
#   compensation_upper  = V_trackB - V_alt15          (DTM swap + Site 10)
#   alt_18_drop         = V_full - V_alt18            (DTM swap + Site 10, flagged in)
#   total_reduction     = V_full - V_alt15
#   remainder           = V_alt15
draws_decomp <- data.table(
  ref_quality_contrib = draws_V$full   - draws_V$trackB,
  compensation_upper  = draws_V$trackB - draws_V$alt15,
  alt18_drop          = draws_V$full   - draws_V$alt18,
  total_reduction     = draws_V$full   - draws_V$alt15,
  remainder           = draws_V$alt15,
  V_full              = draws_V$full
)

# Component as fraction of V_full, draw-wise
draws_frac <- draws_decomp[, .(
  ref_quality_frac  = ref_quality_contrib / V_full,
  comp_upper_frac   = compensation_upper  / V_full,
  alt18_drop_frac   = alt18_drop          / V_full,
  total_reduce_frac = total_reduction     / V_full,
  remainder_frac    = remainder           / V_full
)]

# Build one decomposition row with explicit, named columns. Avoids the brittle
# pattern of passing `as.list(named_vec)` positionally to data.table().
decomp_row <- function(component_name, interp_text, V_diff_draws, frac_draws) {
  data.table(
    component   = component_name,
    interp      = interp_text,
    V_mean_m2   = mean(V_diff_draws),
    V_median_m2 = median(V_diff_draws),
    V_lo_m2     = unname(quantile(V_diff_draws, 0.025)),
    V_hi_m2     = unname(quantile(V_diff_draws, 0.975)),
    frac_mean   = mean(frac_draws),
    frac_lo     = unname(quantile(frac_draws, 0.025)),
    frac_hi     = unname(quantile(frac_draws, 0.975))
  )
}

decomp_summary <- rbindlist(list(
  decomp_row("Reference quality (Ch1 full -> 16-site)",
             "Clean: both err_p3d, drops 3 flagged sites",
             draws_decomp$ref_quality_contrib,
             draws_frac$ref_quality_frac),
  decomp_row("CNN-DTM compensation upper bound (16-site -> step 04 15-site)",
             "Confounded: DTM swap + Site 10 removal",
             draws_decomp$compensation_upper,
             draws_frac$comp_upper_frac),
  decomp_row("Alt: DTM swap with flagged in (Ch1 -> step 04 18-site)",
             "Confounded: DTM swap + Site 10 removal, flagged still in",
             draws_decomp$alt18_drop,
             draws_frac$alt18_drop_frac),
  decomp_row("Total reduction (Ch1 full -> step 04 15-site)",
             "Reference + DTM + Site 10",
             draws_decomp$total_reduction,
             draws_frac$total_reduce_frac),
  decomp_row("Remainder (step 04 15-site V_site)",
             "Unexplained, target of Tasks 1/2/3/6",
             draws_decomp$remainder,
             draws_frac$remainder_frac)
))

fwrite(decomp_summary,
       file.path(TBL_DIR, "coupling_variance_decomposition.csv"))
log_msg("Wrote variance_decomp.csv")
print(decomp_summary[, .(component,
                         V_mean_m2, V_lo_m2, V_hi_m2,
                         frac_mean, frac_lo, frac_hi)])

## ----------------------------------------------------------------------------
## Section 4 -- Bayesian R^2 and FE/RE/site-RE partition
## ----------------------------------------------------------------------------
banner("4. Bayesian R^2 and variance partition")

# For each fit:
#   R2_total = bayes_R2(fit)                          (random + fixed)
#   R2_FE    = bayes_R2(fit, re_formula = NA)         (fixed only)
#   R2_RE    = R2_total - R2_FE
#
# Plus, decompose RE: extract sd_* for each grouping and compute % at site.
# Note: this is a posterior-mean partition, not a draw-wise one. Adequate for
# the 9.8/90.2/98.9 -style summary the manuscript uses.

# Use the model's training data to evaluate R^2. brms uses model.frame(fit).

partition_one <- function(fit, key, ndraws = 1000L) {
  log_msg("  computing R^2 for ", key, " (ndraws = ", ndraws, ") ...")
  # Subsample posterior draws to keep runtime tractable on the 200K-obs Ch1 fit.
  r2_total <- bayes_R2(fit, robust = FALSE, re_formula = NULL, ndraws = ndraws)
  r2_fe    <- bayes_R2(fit, robust = FALSE, re_formula = NA,   ndraws = ndraws)
  r2_total_mean <- as.numeric(r2_total[, "Estimate"])
  r2_fe_mean    <- as.numeric(r2_fe[,    "Estimate"])
  r2_re_mean    <- r2_total_mean - r2_fe_mean

  # Within-RE: variance components from VarCorr (posterior means)
  vc <- VarCorr(fit, summary = TRUE)
  vc_sd <- sapply(names(vc), function(g) {
    s <- vc[[g]]$sd
    # `s` is a matrix with rows for each random parameter; take Intercept
    if (is.matrix(s)) {
      ix <- which(rownames(s) == "Intercept")
      if (length(ix)) as.numeric(s[ix, "Estimate"]) else as.numeric(s[1, "Estimate"])
    } else {
      as.numeric(s["Estimate"])
    }
  })
  vc_var <- vc_sd^2  # variances by grouping factor

  total_re_var  <- sum(vc_var[names(vc_var) != "residual__"])
  site_var      <- if ("site" %in% names(vc_var)) vc_var[["site"]] else NA_real_
  pct_at_site   <- if (is.finite(site_var) && total_re_var > 0)
                     100 * site_var / total_re_var else NA_real_

  data.table(
    fit_key             = key,
    R2_total            = r2_total_mean,
    R2_FE               = r2_fe_mean,
    R2_RE               = r2_re_mean,
    pct_FE_of_R2        = 100 * r2_fe_mean / r2_total_mean,
    pct_RE_of_R2        = 100 * r2_re_mean / r2_total_mean,
    sd_site             = if ("site" %in% names(vc_sd)) vc_sd[["site"]] else NA_real_,
    pct_RE_var_at_site  = pct_at_site,
    n_grouping_factors  = length(vc_sd) - ("residual__" %in% names(vc_sd)),
    grouping_factors    = paste(setdiff(names(vc_sd), "residual__"), collapse = ";")
  )
}

r2_partition <- rbindlist(
  lapply(names(fits), function(k) partition_one(fits[[k]], k))
)
r2_partition <- merge(fit_meta[, .(fit_key, label, response, n_sites)],
                      r2_partition, by = "fit_key")
setcolorder(r2_partition, c("fit_key","label","response","n_sites"))

fwrite(r2_partition,
       file.path(TBL_DIR, "coupling_r2_partition.csv"))
log_msg("Wrote r2_partition.csv")
print(r2_partition[, .(label,
                       R2_total, pct_FE_of_R2, pct_RE_of_R2,
                       pct_RE_var_at_site)])

## ----------------------------------------------------------------------------
## Section 5 -- Site 10 diagnostic for the 16-vs-15 confound
## ----------------------------------------------------------------------------
banner("5. Site 10 diagnostic")

# Extract 16-site's posterior random-intercept draws for Site 10. If small,
# the 16-vs-15 sigma_site change is mostly due to the DTM swap. Approximate
# the "16-site without Site 10" by recomputing sd over remaining sites'
# random-intercept draws.

re_trackB <- ranef(fits$trackB, summary = FALSE)$site  # array: draws x sites x parms
if (is.null(re_trackB)) stop("no site ranef in 16-site fit?")

# Locate the Intercept slot
parm_dim <- dim(re_trackB)
parm_names <- dimnames(re_trackB)[[3]]
int_ix <- which(parm_names == "Intercept")
stopifnot(length(int_ix) == 1L)

re_int <- re_trackB[, , int_ix]   # draws x sites, intercepts only
site_lbls <- dimnames(re_trackB)[[2]]

# Find Site 10 (manuscript number "10")
s10_col <- which(site_lbls == "10")
if (!length(s10_col)) {
  log_msg("WARNING: site '10' not in 16-site ranef; columns are: ",
          paste(site_lbls, collapse = ","))
  site10_diag <- data.table(note = "Site 10 not found in 16-site random intercepts")
} else {
  s10_draws <- re_int[, s10_col]

  # Posterior of sd over all 16 sites (matches 16-site sd_site)
  sd_full16 <- apply(re_int, 1, sd)
  # Posterior of sd over the 15 sites excluding Site 10
  sd_minus10 <- apply(re_int[, -s10_col, drop = FALSE], 1, sd)

  site10_diag <- data.table(
    quantity = c(
      "Site 10 RE intercept (16-site)",
      "16-site sample sd over 16 sites' RE intercepts",
      "16-site sample sd over 15 sites' RE intercepts (Site 10 dropped)",
      "Implied sigma_site change attributable to Site 10 alone (mean-of-draws diff)"
    ),
    mean   = c(mean(s10_draws),  mean(sd_full16),  mean(sd_minus10),
               mean(sd_full16 - sd_minus10)),
    lo     = c(quantile(s10_draws,  .025, names = FALSE),
               quantile(sd_full16,  .025, names = FALSE),
               quantile(sd_minus10, .025, names = FALSE),
               quantile(sd_full16 - sd_minus10, .025, names = FALSE)),
    hi     = c(quantile(s10_draws,  .975, names = FALSE),
               quantile(sd_full16,  .975, names = FALSE),
               quantile(sd_minus10, .975, names = FALSE),
               quantile(sd_full16 - sd_minus10, .975, names = FALSE))
  )

  # Important caveat in the printed note:
  site10_diag <- rbind(
    site10_diag,
    data.table(quantity = "NOTE",
               mean = NA_real_, lo = NA_real_, hi = NA_real_),
    fill = TRUE
  )
  site10_diag[is.na(mean) & quantity == "NOTE",
              mean_label := "Sample SD over draws is NOT equal to the model's sd_site hyperparameter; this gives an empirical-Bayes-like proxy for the per-site spread, suitable for bounding the Site 10 effect, not for replacing a refit."]
}

fwrite(site10_diag,
       file.path(TBL_DIR, "coupling_site10_diagnostic.csv"))
log_msg("Wrote site10_diagnostic.csv")
print(site10_diag)

## ----------------------------------------------------------------------------
## Section 6 -- Per-site RE shifts: Ch1 19-site vs step 04 18-site
## ----------------------------------------------------------------------------
banner("6. Per-site RE shifts (Ch1 vs step 04 18-site)")

# Posterior means of site-level random intercepts for the two fits.
re_full <- ranef(fits$full,   summary = TRUE)$site
re_18   <- ranef(fits$alt18,  summary = TRUE)$site

# Both arrays: sites x stats x parms (with summary = TRUE)
get_intercept_table <- function(re_arr) {
  stopifnot(length(dim(re_arr)) == 3L)
  parm_names <- dimnames(re_arr)[[3]]
  ix <- which(parm_names == "Intercept")
  data.table(
    site     = dimnames(re_arr)[[1]],
    estimate = as.numeric(re_arr[, "Estimate",  ix]),
    lo       = as.numeric(re_arr[, "Q2.5",      ix]),
    hi       = as.numeric(re_arr[, "Q97.5",     ix])
  )
}

re_full_t <- get_intercept_table(re_full); setnames(re_full_t, c("site","RE_full","RE_full_lo","RE_full_hi"))
re_18_t   <- get_intercept_table(re_18);   setnames(re_18_t,   c("site","RE_alt18","RE_alt18_lo","RE_alt18_hi"))

re_shifts <- merge(re_full_t, re_18_t, by = "site", all = TRUE)
re_shifts[, delta := RE_alt18 - RE_full]
re_shifts[, abs_delta := abs(delta)]
setorder(re_shifts, -abs_delta)

# Spearman rank correlation, on sites present in both
matched <- re_shifts[!is.na(RE_full) & !is.na(RE_alt18)]
spearman_r <- suppressWarnings(cor(matched$RE_full, matched$RE_alt18,
                                   method = "spearman"))
log_msg(sprintf("Spearman correlation across %d common sites: %.3f",
                nrow(matched), spearman_r))

fwrite(re_shifts,
       file.path(TBL_DIR, "coupling_random_effect_shifts.csv"))
log_msg("Wrote re_shifts.csv")
print(re_shifts)

## ----------------------------------------------------------------------------
## Section 7 -- Visual summary
## ----------------------------------------------------------------------------
banner("7. Visual summary PDF")

# (a) sigma_site posterior densities for all 4 fits
sigma_long <- rbindlist(lapply(seq_along(draws_sd), function(i) {
  data.table(label = factor(fit_meta$label[i], levels = fit_meta$label),
             sigma_site = draws_sd[[i]])
}))

pa <- ggplot(sigma_long, aes(x = sigma_site, fill = label, colour = label)) +
  geom_density(alpha = 0.35, linewidth = 0.4) +
  scale_x_continuous(limits = c(0, NA)) +
  labs(x = expression(sigma[site] ~ "(m)"), y = "posterior density",
       title = "Posterior of site-level random-intercept SD",
       fill = NULL, colour = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

# (b) Variance decomposition stacked-bar (Ch1 full broken into pieces)
decomp_bar <- data.table(
  component = factor(
    c("Reference quality\n(flagged 3 sites)",
      "Compensation +\nSite 10 (upper bound)",
      "Remainder\n(step 04 15-site)"),
    levels = c("Remainder\n(step 04 15-site)",
               "Compensation +\nSite 10 (upper bound)",
               "Reference quality\n(flagged 3 sites)")),
  V_mean = c(mean(draws_decomp$ref_quality_contrib),
             mean(draws_decomp$compensation_upper),
             mean(draws_decomp$remainder)),
  V_lo   = c(quantile(draws_decomp$ref_quality_contrib, .025, names=FALSE),
             quantile(draws_decomp$compensation_upper,  .025, names=FALSE),
             quantile(draws_decomp$remainder,           .025, names=FALSE)),
  V_hi   = c(quantile(draws_decomp$ref_quality_contrib, .975, names=FALSE),
             quantile(draws_decomp$compensation_upper,  .975, names=FALSE),
             quantile(draws_decomp$remainder,           .975, names=FALSE))
)

pb <- ggplot(decomp_bar, aes(x = "", y = V_mean, fill = component)) +
  geom_col(width = 0.6, colour = "white") +
  labs(x = NULL, y = expression(V[site] ~ "(m"^2*")"),
       title = "Variance decomposition (V_site means)",
       fill  = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "right")

# (c) R^2 partition bar chart
r2_long <- melt(r2_partition,
                id.vars = "label",
                measure.vars = c("R2_FE", "R2_RE"),
                variable.name = "component", value.name = "share")
r2_long[, label := factor(label, levels = fit_meta$label)]
r2_long[, component := factor(component,
                              levels = c("R2_RE", "R2_FE"),
                              labels = c("Random effects", "Fixed effects"))]

pc <- ggplot(r2_long, aes(x = label, y = share, fill = component)) +
  geom_col() +
  labs(x = NULL, y = expression(R^2),
       title = "Fixed vs random share of explained variance",
       fill = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        legend.position = "bottom")

# (d) Per-site RE shifts: Ch1 vs step 04 18-site scatter
matched_plot <- matched[, .(site, RE_full, RE_alt18)]
pd <- ggplot(matched_plot, aes(x = RE_full, y = RE_alt18, label = site)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_vline(xintercept = 0, colour = "grey80") +
  geom_point(size = 2, alpha = 0.8) +
  geom_text(nudge_y = 0.4, size = 2.6, alpha = 0.8) +
  labs(x = "Ch1 19-site RE intercept (err_p3d, m)",
       y = "step 04 18-site RE intercept (err_alt, m)",
       title = sprintf("Per-site RE shifts (Spearman = %.2f)", spearman_r)) +
  theme_minimal(base_size = 11)

pdf_path <- file.path(PLOT_DIR, "coupling_04_posterior_decomposition_summary.pdf")
pdf(pdf_path, width = 11, height = 8)
print((pa | pb) / (pc | pd))
invisible(dev.off())
log_msg("Wrote PDF: ", pdf_path)

## ----------------------------------------------------------------------------
## Section 8 -- Interpretive printout
## ----------------------------------------------------------------------------
banner("8. Interpretive cheat-sheet")

cat("\n")
cat("Posterior summary (sigma_site, m):\n")
print(var_summary[, .(label,
                      sd_mean = round(sd_mean, 3),
                      sd_95CI = sprintf("[%.3f, %.3f]", sd_lo, sd_hi))])

cat("\n")
cat("Variance decomposition (m^2, fraction of V_full):\n")
print(decomp_summary[, .(component,
                         V = sprintf("%.3f [%.3f, %.3f]",
                                     V_mean_m2, V_lo_m2, V_hi_m2),
                         frac = sprintf("%.1f%% [%.1f%%, %.1f%%]",
                                        100*frac_mean, 100*frac_lo, 100*frac_hi))])

cat("\n")
cat("R^2 partition (analogue to Ch1's 9.8 / 90.2 / 98.9):\n")
print(r2_partition[, .(label,
                       R2     = sprintf("%.3f", R2_total),
                       pct_FE = sprintf("%.1f%%", pct_FE_of_R2),
                       pct_RE = sprintf("%.1f%%", pct_RE_of_R2),
                       pct_RE_at_site = sprintf("%.1f%%", pct_RE_var_at_site))])

cat("\n")
cat("Site 10 diagnostic (bounds the 16-vs-15 confound):\n")
print(site10_diag)

cat("\n")
cat("Per-site RE shift summary (Ch1 19-site -> step 04 18-site, top 5):\n")
print(head(re_shifts, 5L))

banner("step 04 decomposition COMPLETE")
log_msg("Tables written to: ", TBL_DIR)
log_msg("Plot written to:  ", pdf_path)
log_msg("CSVs written.")
