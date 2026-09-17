CHKPT <- file.path(Sys.getenv("PROJECT_ROOT"), "checkpoints")
unwrap <- function(p) if (is.list(p) && !is.null(p$data)) p$data else p
mp <- unwrap(readRDS(file.path(CHKPT, "08_model_prep.rds")))
cat("elements in 08_model_prep.rds:\n"); print(names(mp))
for (nm in grep("dtm", names(mp), ignore.case = TRUE, value = TRUE)) {
  d <- mp[[nm]]
  cat("\n---", nm, "--- class:", paste(class(d), collapse = ","), "\n")
  if (is.data.frame(d)) {
    cat("rows:", nrow(d), " cols:", ncol(d), "\n")
    cat("sites:", length(unique(as.character(d$site))), "->",
        paste(sort(unique(as.character(d$site))), collapse = " "), "\n")
    cat("error columns:", paste(grep("error", names(d), value = TRUE), collapse = " "), "\n")
    cat("lc_l1_code:", paste(sort(unique(as.character(d$lc_l1_code))), collapse = " "), "\n")
    cat("all columns:\n"); print(names(d))
  }
}
