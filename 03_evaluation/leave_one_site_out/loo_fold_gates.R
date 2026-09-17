#!/usr/bin/env Rscript
# =====================================================================
# loo_fold_gates.R  --  read-only.  Supersedes phase1_gates.R.
#
# Two fixes to v1, both v1's fault:
#   A. The rung labels are "1_oracle_site_known",
#      "2_new_site_eco_is_known", "3_new_site_new_ecoregion". v1 matched
#      rung 2 as "new site AND NOT eco", but rung 2's own label contains
#      "eco", so nothing matched. Now matched on the leading digit.
#   B. mod_chm_s2 holds ALL 19 sites -- the 1/2/3 exclusion is applied
#      downstream at fit time, not in the object. v1 built the CHM
#      ecoregion map from 19 sites and printed the wrong verdict. The
#      analysis site sets are now constructed explicitly.
#
# GATE A -- the 3.362 vs 3.363 DTM rung-2 discrepancy (plan sec 4 item 7:
#   resolve before any table is built). Recomputes the pooled figure
#   every plausible way and names which recipe yields which value.
#
# GATE B -- site-4 ecoregion orphan (checkpoint sec 6.3b, "inferred, not
#   verified"). Rather than check site 4 alone, this DERIVES eco_known
#   from scratch for every site in each product and compares the result
#   against the recorded sets (CHM 4/5/15/19/20, DTM 5/15/19/20). If the
#   derivation reproduces both, the fold classification behind sec 6.2,
#   sec 6.3 and sec 6.3b is verified end to end, not just site 4.
#
# GATE C -- ecoregion inventory on the real analysis bases, to settle the
#   "nine vs ten Level II" inconsistency (review sec 5 item 5) and to
#   identify the level whose label printed blank in the v1 run.
#
# WRITES NOTHING. Safe to run alongside the M2 DTM refit.
#
# USAGE:
#   export PROJECT_ROOT=/gpfs/data1/vclgp/lmaden/chpt1
#   export R_LIBS=/gpfs/data1/vclgp/lmaden/Rlib
#   nice -n 15 Rscript $PROJECT_ROOT/scripts/reviewed/loo_fold_gates.R \
#     2>&1 | tee $PROJECT_ROOT/manuscript_tables/loo_fold_gates.txt
# =====================================================================

suppressPackageStartupMessages(library(data.table))

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/gpfs/data1/vclgp/lmaden/chpt1")
LOO_DIR <- file.path(PROJECT_ROOT, "manuscript_tables", "loo")

hr  <- function() cat(strrep("=", 78), "\n")
sec <- function(x) { cat("\n"); hr(); cat("## ", x, "\n"); hr() }

cat("loo_fold_gates.R -- read-only\n")
cat("started: ", format(Sys.time()), "\n", sep = "")

# =====================================================================
sec("GATE A -- pooled DTM rung-2 RMSE, every plausible recipe")
# =====================================================================

all_csv <- list.files(LOO_DIR, pattern = "^loo_metrics_dtm_site[0-9]+(_ad099)?\\.csv$")
orig  <- grep("^loo_metrics_dtm_site[0-9]+\\.csv$",       all_csv, value = TRUE)
ad099 <- grep("^loo_metrics_dtm_site[0-9]+_ad099\\.csv$", all_csv, value = TRUE)
if (!length(orig)) stop("no DTM metrics files under ", LOO_DIR)

site_of <- function(f) as.integer(sub("^loo_metrics_dtm_site([0-9]+).*$", "\\1", f))
read_set <- function(files) rbindlist(lapply(files, function(f) {
  d <- fread(file.path(LOO_DIR, f)); d[, `:=`(.site = site_of(f), .file = f)]; d
}), fill = TRUE, use.names = TRUE)

preferred <- c(ad099, orig[!site_of(orig) %in% site_of(ad099)])
cat(sprintf("%d folds total | %d substituted by _ad099: %s\n",
            length(preferred), length(ad099),
            paste(sort(site_of(ad099)), collapse = ", ")))

DT_pref <- read_set(preferred)
DT_orig <- read_set(orig)

pick <- function(d, pats, what) {
  for (p in pats) { h <- grep(p, names(d), ignore.case = TRUE, value = TRUE)
                    if (length(h)) return(h[1]) }
  cat("  !! no column for ", what, "; available: ",
      paste(names(d), collapse = ", "), "\n", sep = ""); NA_character_
}
C_rung <- pick(DT_pref, c("^rung$", "rung"),                  "rung")
C_rmse <- pick(DT_pref, c("^rmse$", "rmse"),                  "RMSE")
C_n    <- pick(DT_pref, c("^n$", "^n_", "^nobs$", "^n_obs$"), "n")
C_eco  <- pick(DT_pref, c("eco_known", "eco.*known"),         "eco_known")
cat("columns -> rung: ", C_rung, " | rmse: ", C_rmse, " | n: ", C_n,
    " | eco_known: ", C_eco, "\n", sep = "")
cat("rung levels: ",
    paste(sort(unique(as.character(DT_pref[[C_rung]]))), collapse = " | "), "\n", sep = "")

# FIX: match the leading digit. Labels carry "eco" in BOTH rungs 2 and 3.
rung_is <- function(v, k) grepl(paste0("^\\s*", k), as.character(v))

pooled <- function(d, k, valid_only = FALSE, sqrt_mse = FALSE, unweighted = FALSE) {
  x <- d[rung_is(get(C_rung), k)]
  if (valid_only) { if (is.na(C_eco)) return(NA_real_); x <- x[as.logical(get(C_eco))] }
  if (!nrow(x)) return(NA_real_)
  r <- x[[C_rmse]]; n <- x[[C_n]]
  v <- if (unweighted) mean(r) else if (sqrt_mse) sqrt(sum(n * r^2) / sum(n)) else sum(n * r) / sum(n)
  attr(v, "folds") <- nrow(x); attr(v, "n") <- sum(n); v
}
show <- function(label, v) {
  if (is.na(v)) { cat(sprintf("  %-56s  --\n", label)); return(invisible()) }
  cat(sprintf("  %-56s  %.4f   (folds %2d, n %s)\n", label, v,
              attr(v, "folds"), format(attr(v, "n"), big.mark = ",")))
}

cat("\nWhich recipe gives 3.362 and which gives 3.363?\n\n")
show("1. n-weighted RMSE, all folds, _ad099 preferred",   pooled(DT_pref, 2))
show("2. n-weighted RMSE, eco_known folds only",          pooled(DT_pref, 2, valid_only = TRUE))
show("3. unweighted mean of per-site RMSE",               pooled(DT_pref, 2, unweighted = TRUE))
show("4. n-weighted RMSE, ORIGINALS ONLY (glob trap)",    pooled(DT_orig, 2))
show("5. sqrt(n-weighted mean MSE), _ad099 preferred",    pooled(DT_pref, 2, sqrt_mse = TRUE))
show("6. sqrt(n-weighted mean MSE), ORIGINALS ONLY",      pooled(DT_orig, 2, sqrt_mse = TRUE))

cat("\nCross-check against checkpoint sec 6.2 (expect 3.003 / 3.362 / 3.383 all-fold):\n")
show("rung 1 oracle,        n-weighted, all folds",       pooled(DT_pref, 1))
show("rung 2 new site,      n-weighted, all folds",       pooled(DT_pref, 2))
show("rung 3 new site+eco,  n-weighted, all folds",       pooled(DT_pref, 3))
cat("and the 14 valid folds (expect 3.106 / 3.534 / 3.567):\n")
show("rung 1 oracle,        eco_known only",              pooled(DT_pref, 1, valid_only = TRUE))
show("rung 2 new site,      eco_known only",              pooled(DT_pref, 2, valid_only = TRUE))
show("rung 3 new site+eco,  eco_known only",              pooled(DT_pref, 3, valid_only = TRUE))

cat("\nPer-fold rung-2 rows used in recipe 1:\n")
keep <- c(".site", C_n, C_rmse, if (!is.na(C_eco)) C_eco, ".file")
print(DT_pref[rung_is(get(C_rung), 2), ..keep][order(.site)])

# =====================================================================
sec("GATE B -- eco_known derived from scratch, both products")
# =====================================================================

unwrap <- function(o) if (is.list(o) && !is.null(o$data)) o$data else o
mp <- unwrap(readRDS(file.path(PROJECT_ROOT, "checkpoints", "08_model_prep.rds")))

CHM_EXCLUDE <- c("1", "2", "3")   # applied at fit time, NOT in mod_chm_s2
pairs_of <- function(el, drop = character(0)) {
  d <- as.data.table(mp[[el]])
  d <- d[!as.character(site) %in% drop]
  unique(d[, .(site = as.character(site), eco = as.character(ecoregion))])
}
chm_pairs <- pairs_of("mod_chm_s2", CHM_EXCLUDE)
dtm_pairs <- pairs_of("mod_dtm_s2")

cat("CHM analysis sites (", uniqueN(chm_pairs$site), "): ",
    paste(sort(unique(chm_pairs$site)), collapse = ", "), "\n", sep = "")
cat("DTM analysis sites (", uniqueN(dtm_pairs$site), "): ",
    paste(sort(unique(dtm_pairs$site)), collapse = ", "), "\n", sep = "")

# eco_known(s) == every ecoregion of s also occurs at some OTHER site.
derive <- function(p) {
  sites <- sort(unique(p$site))
  rbindlist(lapply(sites, function(s) {
    mine <- p[site == s, eco]
    elsewhere <- unique(p[site != s, eco])
    orphan <- setdiff(mine, elsewhere)
    data.table(site = s, n_eco = length(mine),
               eco_known = length(orphan) == 0L,
               orphan_eco = paste(orphan, collapse = "; "))
  }))
}
chm_ek <- derive(chm_pairs); dtm_ek <- derive(dtm_pairs)

cat("\nCHM:\n"); print(chm_ek[order(eco_known, as.integer(site))])
cat("\nDTM:\n"); print(dtm_ek[order(eco_known, as.integer(site))])

cmp <- function(ek, expected, label) {
  got <- sort(as.integer(ek[eco_known == FALSE, site]))
  cat(sprintf("\n%s eco_known == FALSE\n  derived : %s\n  recorded: %s\n  MATCH   : %s\n",
              label, paste(got, collapse = ", "),
              paste(expected, collapse = ", "),
              identical(got, as.integer(expected))))
}
cmp(chm_ek, c(4, 5, 15, 19, 20), "CHM")
cmp(dtm_ek, c(5, 15, 19, 20),    "DTM")

cat("\n-- sec 6.3b's specific claim about site 4 --\n")
s4 <- chm_pairs[site == "4", eco]
cat("site 4 ecoregions: ", paste(s4, collapse = " | "), "\n", sep = "")
for (e in s4) {
  chm_o <- setdiff(pairs_of("mod_chm_s2")[eco == e, site], "4")   # 19-site view
  cat("  ", e, "\n    all sites sharing it (19-site view): ",
      paste(sort(chm_o), collapse = ", "), "\n", sep = "")
  cat("    still present in the 16-site CHM: ",
      { k <- setdiff(chm_pairs[eco == e, site], "4")
        if (length(k)) paste(sort(k), collapse = ", ") else "NONE" }, "\n", sep = "")
  cat("    still present in the 18-site DTM: ",
      { k <- setdiff(dtm_pairs[eco == e, site], "4")
        if (length(k)) paste(sort(k), collapse = ", ") else "NONE" }, "\n", sep = "")
}

# =====================================================================
sec("GATE C -- ecoregion inventory (review sec 5 item 5: nine vs ten)")
# =====================================================================

inventory <- function(p, label) {
  tb <- p[, .(n_sites = uniqueN(site)), by = eco][order(-n_sites, eco)]
  cat("\n", label, " -- ", nrow(tb), " distinct ecoregion levels\n", sep = "")
  tb[, blank := is.na(eco) | trimws(eco) == ""]
  print(tb)
  nb <- tb[blank == TRUE]
  if (nrow(nb)) {
    cat("  !! ", nrow(nb), " level(s) with a blank/NA label. Raw value(s):\n", sep = "")
    print(lapply(nb$eco, function(x) list(value = x, is_na = is.na(x),
                                          nchar = if (is.na(x)) NA_integer_ else nchar(x))))
    cat("  sites carrying it: ",
        paste(sort(p[is.na(eco) | trimws(eco) == "", site]), collapse = ", "), "\n", sep = "")
  }
  cat("  NAMED (non-blank) levels: ", nrow(tb) - nrow(nb), "\n", sep = "")
}
inventory(chm_pairs, "16-site CHM")
inventory(dtm_pairs, "18-site DTM")
cat("\nManuscript says nine Level II in four places and ten in sec 4.4.\n")
cat("The NAMED count above is the number to reconcile against.\n")

sec("DONE -- nothing was written except this log")
cat("finished: ", format(Sys.time()), "\n", sep = "")
