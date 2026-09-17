suppressPackageStartupMessages({library(data.table); library(ranger)})
PR <- Sys.getenv("PROJECT_ROOT"); set.seed(2026)
pick <- function(o, nm) if (!is.null(o[[nm]])) o[[nm]] else o$data[[nm]]
rmse <- function(p, y) sqrt(mean((y - p)^2))

harmonise <- function(te, tr, cats) {
  for (cc in cats) {
    a <- as.character(tr[[cc]]); b <- as.character(te[[cc]])
    lv <- names(sort(table(a), decreasing = TRUE))
    b[!b %in% lv] <- lv[1]
    tr[[cc]] <- factor(a, levels = lv); te[[cc]] <- factor(b, levels = lv)
  }
  list(te = te, tr = tr)
}

run <- function(prod, dat, resp) {
  d <- as.data.frame(dat); d$site <- as.integer(as.character(d$site))
  preds <- setdiff(names(d), c(resp, "site"))
  cats  <- preds[sapply(d[preds], function(x) is.factor(x) || is.character(x))]
  eco   <- grep("eco",  cats, value = TRUE, ignore.case = TRUE)
  lcc   <- grep("^lc",  cats, value = TRUE, ignore.case = TRUE)
  rh    <- grep("rh.?98", preds, value = TRUE, ignore.case = TRUE)
  slp   <- grep("slope_mean", preds, value = TRUE, ignore.case = TRUE)
  cov   <- grep("^cover", preds, value = TRUE, ignore.case = TRUE)
  lin_b <- c(rh, slp, cov, lcc)

  cat("  DETECTED  categorical:", paste(cats, collapse = " "), "\n")
  cat("            ecoregion  :", paste(eco, collapse = " "), "\n")
  cat("            linear base:", paste(lin_b, collapse = " "), "\n")
  cat("            n predictors:", length(preds), "\n\n")
  stopifnot(length(eco) == 1, length(lcc) == 1, length(rh) == 1,
            length(slp) == 1, length(cov) == 1)

  out <- list()
  for (s in sort(unique(d$site))) {
    te <- d[d$site == s, ]; tr <- d[d$site != s, ]; y <- te[[resp]]
    r <- data.frame(site = s, n = nrow(te),
                    site_offset = rmse(rep(mean(y), nrow(te)), y),
                    global      = rmse(rep(mean(tr[[resp]]), nrow(te)), y))
    for (tag in c("blind", "aware")) {
      fxl <- if (tag == "aware") c(lin_b, eco) else lin_b
      fxr <- if (tag == "aware") preds         else setdiff(preds, eco)
      zl <- harmonise(te, tr, intersect(fxl, cats))
      r[[paste0("lin_", tag)]] <- rmse(
        predict(lm(reformulate(fxl, resp), data = zl$tr), newdata = zl$te), y)
      zr <- harmonise(te, tr, intersect(fxr, cats))
      rf <- ranger(reformulate(fxr, resp), data = zr$tr,
                   num.trees = 300, num.threads = 16, seed = 2026)
      r[[paste0("rf_", tag)]] <- rmse(predict(rf, data = zr$te)$predictions, y)
    }
    out[[length(out) + 1]] <- r
    cat(sprintf("  %s site %-3d n=%-6d done\n", prod, s, nrow(te))); flush.console()
  }
  rbindlist(out)
}

for (cfg in list(list("chm","10b_chm_sensitivity.rds","fit_chm_16","chm_error_mean"),
                 list("dtm","10_models_stage2.rds","fit_dtm_s2","dtm_error_mean"))) {
  cat("\n\n########", toupper(cfg[[1]]), "########\n")
  o <- readRDS(file.path(PR, "checkpoints", cfg[[2]]))
  res <- run(cfg[[1]], pick(o, cfg[[3]])$data, cfg[[4]]); rm(o); gc()
  fwrite(res, file.path(PR, "manuscript_tables", "loo",
                        paste0("m7_baselines_", cfg[[1]], "_persite.csv")))
  cols <- setdiff(names(res), c("site", "n"))
  cat("\n-- POOLED (n-weighted RMSE) --\n")
  print(data.table(baseline = cols,
    rmse = sapply(cols, function(cc) round(sqrt(sum(res[[cc]]^2 * res$n) / sum(res$n)), 3))))
}
