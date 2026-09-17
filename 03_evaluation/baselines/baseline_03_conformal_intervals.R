suppressPackageStartupMessages({library(data.table); library(ranger)})
PR <- Sys.getenv("PROJECT_ROOT"); set.seed(2026)
A <- commandArgs(trailingOnly = TRUE); SMOKE <- length(A) > 0 && A[1] == "smoke"
NDRAWS <- 1000; CHUNK <- 5000; PP <- (seq_len(NDRAWS) - 0.5) / NDRAWS
DEGEN <- list(chm = c(4,5,15,19,20), dtm = c(5,15,19,20))
BL <- c("site_offset","global","lin_blind","rf_blind","lin_aware","rf_aware")

pick <- function(o, nm) if (!is.null(o[[nm]])) o[[nm]] else o$data[[nm]]
harmonise <- function(te, tr, cats) {
  for (cc in cats) {
    a <- as.character(tr[[cc]]); b <- as.character(te[[cc]])
    lv <- names(sort(table(a), decreasing = TRUE))
    b[!b %in% lv] <- lv[1]
    tr[[cc]] <- factor(a, levels = lv); te[[cc]] <- factor(b, levels = lv)
  }
  list(te = te, tr = tr)
}
crps_cols <- function(draws, obs) {
  m <- nrow(draws)
  vapply(seq_along(obs), function(j) {
    x <- sort(draws[, j]); yy <- obs[j]; i <- seq_len(m)
    (2 / m^2) * sum((x - yy) * (m * (yy < x) - i + 0.5))
  }, numeric(1))
}
score <- function(draws, obs) {
  pm <- colMeans(draws)
  q  <- apply(draws, 2, quantile, probs = c(0.025,0.05,0.25,0.75,0.95,0.975))
  res <- pm - obs
  data.table(n = length(obs), rmse = sqrt(mean(res^2)), bias = mean(res),
             mae = mean(abs(res)), crps = mean(crps_cols(draws, obs)),
             cov50 = mean(obs >= q[3,] & obs <= q[4,]),
             cov90 = mean(obs >= q[2,] & obs <= q[5,]),
             cov95 = mean(obs >= q[1,] & obs <= q[6,]),
             width90 = mean(q[5,] - q[2,]))
}
pool_w <- function(tabs) {
  d <- rbindlist(tabs); N <- sum(d$n)
  data.table(n = N, rmse = sqrt(sum(d$rmse^2*d$n)/N), bias = sum(d$bias*d$n)/N,
             mae = sum(d$mae*d$n)/N, crps = sum(d$crps*d$n)/N,
             cov50 = sum(d$cov50*d$n)/N, cov90 = sum(d$cov90*d$n)/N,
             cov95 = sum(d$cov95*d$n)/N, width90 = sum(d$width90*d$n)/N)
}
r3 <- function(x) { for (j in names(x)) if (is.numeric(x[[j]]) && j != "n")
  set(x, j = j, value = round(x[[j]], 3)); x[] }

for (cfg in list(list("chm","10b_chm_sensitivity.rds","fit_chm_16","chm_error_mean"),
                 list("dtm","10_models_stage2.rds","fit_dtm_s2","dtm_error_mean"))) {
  prod <- cfg[[1]]; resp <- cfg[[4]]
  if (SMOKE && prod == "dtm") next
  cat("\n########", toupper(prod), "########\n")
  o <- readRDS(file.path(PR, "checkpoints", cfg[[2]]))
  d <- as.data.frame(pick(o, cfg[[3]])$data); rm(o); gc()
  d$site <- as.integer(as.character(d$site))
  preds <- setdiff(names(d), c(resp, "site"))
  cats  <- preds[sapply(d[preds], function(x) is.factor(x) || is.character(x))]
  eco <- grep("eco", cats, value=TRUE, ignore.case=TRUE)
  lcc <- grep("^lc", cats, value=TRUE, ignore.case=TRUE)
  rh  <- grep("rh.?98", preds, value=TRUE, ignore.case=TRUE)
  slp <- grep("slope_mean", preds, value=TRUE, ignore.case=TRUE)
  cvv <- grep("^cover", preds, value=TRUE, ignore.case=TRUE)
  lin_b <- c(rh, slp, cvv, lcc)
  cat("  DETECTED eco:", eco, "| lc:", lcc, "| lin:", paste(lin_b, collapse=" "),
      "| n preds:", length(preds), "\n")
  stopifnot(length(eco)==1, length(lcc)==1, length(rh)==1, length(slp)==1, length(cvv)==1)

  sites <- sort(unique(d$site))
  if (SMOKE) sites <- as.integer(names(sort(table(d$site)))[1:3])
  cat("  folds this run:", paste(sites, collapse=" "), "\n")

  P <- rbindlist(lapply(sites, function(s) {
    te <- d[d$site == s, ]; tr <- d[d$site != s, ]; y <- te[[resp]]
    r <- data.table(site = s, obs = y, site_offset = mean(y),
                    global = mean(tr[[resp]]))
    for (tag in c("blind","aware")) {
      fxl <- if (tag=="aware") c(lin_b, eco) else lin_b
      fxr <- if (tag=="aware") preds else setdiff(preds, eco)
      zl <- harmonise(te, tr, intersect(fxl, cats))
      set(r, j = paste0("lin_", tag), value = as.numeric(
        predict(lm(reformulate(fxl, resp), data = zl$tr), newdata = zl$te)))
      zr <- harmonise(te, tr, intersect(fxr, cats))
      rf <- ranger(reformulate(fxr, resp), data = zr$tr, num.trees = 300,
                   num.threads = 16, seed = 2026)
      set(r, j = paste0("rf_", tag),
          value = as.numeric(predict(rf, data = zr$te)$predictions))
    }
    cat(sprintf("  preds %s site %-3d n=%-6d done\n", prod, s, nrow(te)))
    flush.console(); r
  }))

  chk <- P[, lapply(.SD, function(p) sqrt(mean((obs - p)^2))), by = site, .SDcols = BL]
  old <- fread(file.path(PR,"manuscript_tables","loo",
               paste0("m7_baselines_", prod, "_persite.csv")))
  mm  <- merge(chk, old[, c("site", BL), with=FALSE], by="site", suffixes=c("",".old"))
  dif <- sapply(BL, function(b) max(abs(mm[[b]] - mm[[paste0(b,".old")]])))
  cat("\n-- RMSE REPRODUCTION CHECK, max abs diff vs saved --\n"); print(signif(dif, 3))
  cat("  ALL MATCH (< 1e-6):", all(dif < 1e-6), "\n\n")

  res <- list()
  for (s in sites) {
    Ps <- P[site == s]; oth <- P[site != s]
    for (b in BL) {
      ev <- oth$obs - oth[[b]]
      bp <- list(
        nw = as.numeric(quantile(ev, PP, names=FALSE, type=7)),
        eq = as.numeric(quantile(unlist(lapply(split(ev, oth$site),
               function(z) quantile(z, PP, names=FALSE, type=7))),
               PP, names=FALSE, type=7)))
      pools <- list(nw_ctr = bp$nw - mean(bp$nw), nw_raw = bp$nw,
                    eq_ctr = bp$eq - mean(bp$eq), eq_raw = bp$eq)
      for (pn in names(pools)) {
        R <- pools[[pn]]
        p <- Ps[[b]]; y <- Ps$obs; n <- length(y)
        tabs <- lapply(seq(1, n, by=CHUNK), function(k) {
          rng <- k:min(k + CHUNK - 1, n)
          score(outer(R, p[rng], "+"), y[rng])
        })
        res[[length(res)+1]] <- cbind(data.table(site=s, baseline=b, pool=pn),
                                      pool_w(tabs))
      }
    }
    cat(sprintf("  score %s site %-3d done\n", prod, s)); flush.console()
  }
  RES <- rbindlist(res)
  fwrite(RES, file.path(PR,"manuscript_tables","loo", paste0(
    "baseline_03_conformal_intervals_", prod, if (SMOKE) "_SMOKE" else "", "_persite.csv")))

  SC <- c("n","rmse","bias","mae","crps","cov50","cov90","cov95","width90")
  for (pn in c("eq_ctr","eq_raw")) for (basis in c("all","valid")) {
    keep <- if (basis=="valid") setdiff(sites, DEGEN[[prod]]) else sites
    cat(sprintf("\n-- %s | pool=%s | basis=%s | folds=%d --\n",
                toupper(prod), pn, basis, length(keep)))
    print(r3(RES[pool==pn & site %in% keep][, pool_w(list(.SD)),
             by=baseline, .SDcols=SC]))
  }
}
