suppressPackageStartupMessages(library(brms))
CHKPT <- file.path(Sys.getenv("PROJECT_ROOT"), "checkpoints")
unwrap <- function(p) if (is.list(p) && !is.null(p$data)) p$data else p

mp <- unwrap(readRDS(file.path(CHKPT, "08_model_prep.rds")))
d  <- mp$mod_chm_s2
cat("=== mod_chm_s2 ===\n")
cat("rows:", nrow(d), " cols:", ncol(d), "\n")
cat("sites:", length(unique(as.character(d$site))), "->",
    paste(sort(unique(as.character(d$site))), collapse = " "), "\n")
cat("rows per site:\n"); print(sort(table(as.character(d$site)), decreasing = TRUE))
cat("error cols:", paste(grep("error", names(d), value = TRUE), collapse = " "), "\n")
cat("lc_l1_code:", paste(sort(unique(as.character(d$lc_l1_code))), collapse = " "), "\n")
cat("ecoregions:", length(unique(as.character(d$ecoregion))), "\n")
rm(mp); gc(verbose = FALSE)

cat("\n=== fit_chm_16_tuned.rds (the tuned confirmation fit) ===\n")
f <- file.path(Sys.getenv("PROJECT_ROOT"), "models", "confirm", "fit_chm_16_tuned.rds")
cat("exists:", file.exists(f), "\n")
if (file.exists(f)) {
  fit <- readRDS(f)
  cat("nobs:", nobs(fit), " ndraws:", brms::ndraws(fit), "\n")
  cat("sites in fit data:", length(unique(as.character(fit$data$site))), "->",
      paste(sort(unique(as.character(fit$data$site))), collapse = " "), "\n")
  cat("lc column used:", paste(grep("lc", names(fit$data), value = TRUE), collapse = " "), "\n")
  cat("lc levels:", paste(levels(fit$data$lc_grp), collapse = " "), "\n")
  cat("\nformula:\n"); print(fit$formula)
  np <- brms::nuts_params(fit); td <- np$Value[np$Parameter == "treedepth__"]
  cat("\ndivergences:", sum(np$Value[np$Parameter == "divergent__"]),
      " mean treedepth:", round(mean(td), 2), " max:", max(td),
      " pct>=12:", round(100 * mean(td >= 12), 2), "%\n")
}
