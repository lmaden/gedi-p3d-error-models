# find_dtm_fit_standalone.R
# Lists every brmsfit object across the checkpoint .rds files, with its response
# variable and N. Use to locate fit_dtm_s2 (response = dtm_error_mean) if the
# auto-discovery in section_G_fig6_regen.R does not find it. Paste the output back.
suppressPackageStartupMessages(library(brms))
CHPT1_ROOT <- Sys.getenv("CHPT1_ROOT", "/gpfs/data1/vclgp/lmaden/chpt1")
CHKPT_DIR  <- Sys.getenv("CHKPT_DIR",  file.path(CHPT1_ROOT, "checkpoints"))
resp_of <- function(x) tryCatch({ r <- x$formula$resp
  if (is.null(r) || !nzchar(r)) all.vars(x$formula$formula)[1] else r },
  error = function(e) NA_character_)
files <- list.files(CHKPT_DIR, pattern="\\.rds$", full.names=TRUE)
files <- files[!grepl("_cache|01_data_ingest", basename(files))]
files <- files[order(file.info(files)$size)]
cat("Scanning", length(files), "checkpoint files for brmsfit objects...\n\n")
for (fp in files) {
  obj <- tryCatch(readRDS(fp), error=function(e) NULL); if (is.null(obj)) next
  bag <- list(); if (inherits(obj,"brmsfit")) bag <- c(bag, list(`<top-level>`=obj))
  if (is.list(obj)) bag <- c(bag, obj, if (is.list(obj$data)) obj$data else NULL)
  nms <- names(bag)
  for (i in seq_along(bag)) {
    x <- bag[[i]]
    if (inherits(x,"brmsfit"))
      cat(sprintf("  %-32s $ %-22s  resp=%-16s N=%d\n",
                  basename(fp), if (length(nms)>=i && nzchar(nms[i])) nms[i] else "?",
                  resp_of(x), nrow(x$data)))
  }
  rm(obj, bag); gc(verbose=FALSE)
}
cat("\nDone. The DTM model is the brmsfit with resp=dtm_error_mean.\n")
