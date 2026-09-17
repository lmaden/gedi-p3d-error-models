# =====================================================================
# forest_type_accuracy_holdout_full.R
#
# PURPOSE:
#   Compute CHM and DTM predictive accuracy by forest type on FULL
#   holdout data (no subsampling caps) for consistent Tables 13 & 14
#
# OUTPUT:
#   - chm_forest_type_accuracy_holdout_full.csv
#   - dtm_forest_type_accuracy_holdout_full.csv
#
# NOTE: This uses the full holdout set, so posterior_predict will be
#       called on large datasets. Expect 20-40 min depending on
#       cluster load. Can Ctrl+C if needed.
#
# RUN: source("forest_type_accuracy_holdout_full.R")
# =====================================================================

Sys.setenv(DISPLAY = "")
options(device = pdf)

suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
})

cat("\n", strrep("=", 70), "\n")
cat("FOREST-TYPE ACCURACY ON FULL HOLDOUT (NO CAPS)\n")
cat(strrep("=", 70), "\n\n")

# =====================================================================
# CONFIGURATION
# =====================================================================

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/gpfs/data1/vclgp/lmaden/chpt1")
CHECKPOINT_DIR <- file.path(PROJECT_ROOT, "checkpoints")
out_tables <- file.path(PROJECT_ROOT, "tables")
dir.create(out_tables, showWarnings = FALSE, recursive = TRUE)

stage2_frac <- as.numeric(Sys.getenv("STAGE2_FRAC", "0.33"))

forested_classes <- c("BDF", "DNF", "EBF", "ENF")
lc_names <- c(
  "BDF" = "Broadleaf Deciduous",
  "DNF" = "Deciduous Needleleaf",
  "EBF" = "Evergreen Broadleaf",
  "ENF" = "Evergreen Needleleaf"
)

# Posterior draws for predictions (100 is sufficient for mean predictions)
n_draws <- 100

# =====================================================================
# LOAD MODELS
# =====================================================================

cat("Loading models...\n")
model_data <- readRDS(file.path(CHECKPOINT_DIR, "10_models_stage2.rds"))
fit_chm_s2 <- model_data$data$fit_chm_s2
fit_dtm_s2 <- model_data$data$fit_dtm_s2
cat("  Models loaded\n\n")

# =====================================================================
# LOAD DATA AND RECREATE HOLDOUT
# =====================================================================

cat("Loading full dataset...\n")
full_data <- readRDS(file.path(CHECKPOINT_DIR, "01_data_ingest.rds"))

chm_df_full <- full_data$data$chm_df %>%
  filter(is.finite(chm_error_mean), !is.na(lc_l1_code), !is.na(site), !is.na(ecoregion)) %>%
  mutate(lc_l1_code = factor(lc_l1_code),
         ecoregion = factor(ecoregion),
         site = factor(site))

dtm_df_full <- full_data$data$dtm_df %>%
  filter(is.finite(dtm_error_mean), !is.na(lc_l1_code), !is.na(site), !is.na(ecoregion)) %>%
  mutate(lc_l1_code = factor(lc_l1_code),
         ecoregion = factor(ecoregion),
         site = factor(site))

cat(sprintf("  Full data: %s CHM, %s DTM\n",
            format(nrow(chm_df_full), big.mark = ","),
            format(nrow(dtm_df_full), big.mark = ",")))

# Recreate training split (same seed and fraction as section_08)
set.seed(2025)

chm_df_full$original_row <- 1:nrow(chm_df_full)
dtm_df_full$original_row <- 1:nrow(dtm_df_full)

chm_train <- chm_df_full %>% group_by(site) %>% sample_frac(stage2_frac) %>% ungroup()
dtm_train <- dtm_df_full %>% group_by(site) %>% sample_frac(stage2_frac) %>% ungroup()

# Factor levels from training
train_lc_chm <- levels(droplevels(chm_train$lc_l1_code))
train_site_chm <- levels(droplevels(chm_train$site))
train_eco_chm <- levels(droplevels(chm_train$ecoregion))

train_lc_dtm <- levels(droplevels(dtm_train$lc_l1_code))
train_site_dtm <- levels(droplevels(dtm_train$site))
train_eco_dtm <- levels(droplevels(dtm_train$ecoregion))

# Create holdout sets
chm_holdout <- chm_df_full %>%
  filter(!original_row %in% chm_train$original_row) %>%
  select(-original_row) %>%
  filter(lc_l1_code %in% train_lc_chm,
         site %in% train_site_chm,
         ecoregion %in% train_eco_chm) %>%
  mutate(lc_l1_code = factor(lc_l1_code, levels = train_lc_chm),
         site = factor(site, levels = train_site_chm),
         ecoregion = factor(ecoregion, levels = train_eco_chm))

dtm_holdout <- dtm_df_full %>%
  filter(!original_row %in% dtm_train$original_row) %>%
  select(-original_row) %>%
  filter(lc_l1_code %in% train_lc_dtm,
         site %in% train_site_dtm,
         ecoregion %in% train_eco_dtm) %>%
  mutate(lc_l1_code = factor(lc_l1_code, levels = train_lc_dtm),
         site = factor(site, levels = train_site_dtm),
         ecoregion = factor(ecoregion, levels = train_eco_dtm))

# Exclude Site 10 from DTM
dtm_holdout <- dtm_holdout %>% filter(site != "10")

cat(sprintf("  CHM holdout: %s observations\n", format(nrow(chm_holdout), big.mark = ",")))
cat(sprintf("  DTM holdout: %s observations\n", format(nrow(dtm_holdout), big.mark = ",")))

# Filter to forested classes
chm_forest <- chm_holdout %>% filter(lc_l1_code %in% forested_classes)
dtm_forest <- dtm_holdout %>% filter(lc_l1_code %in% forested_classes)

cat(sprintf("  CHM forested holdout: %s observations\n", format(nrow(chm_forest), big.mark = ",")))
cat(sprintf("  DTM forested holdout: %s observations\n", format(nrow(dtm_forest), big.mark = ",")))

# Print class breakdown
cat("\n  CHM class counts:\n")
print(table(chm_forest$lc_l1_code))
cat("\n  DTM class counts:\n")
print(table(dtm_forest$lc_l1_code))

# =====================================================================
# HELPER: Compute accuracy metrics
# =====================================================================

compute_accuracy <- function(observed, predicted) {
  valid <- is.finite(observed) & is.finite(predicted)
  obs <- observed[valid]
  pred <- predicted[valid]
  
  data.frame(
    N = length(obs),
    RMSE = round(sqrt(mean((obs - pred)^2)), 2),
    MAE = round(mean(abs(obs - pred)), 2),
    Bias = round(mean(pred - obs), 2),
    R2 = round(cor(obs, pred)^2, 2),
    stringsAsFactors = FALSE
  )
}

# =====================================================================
# CHM FOREST-TYPE ACCURACY
# =====================================================================

cat("\n", strrep("=", 70), "\n")
cat("CHM FOREST-TYPE ACCURACY (FULL HOLDOUT)\n")
cat(strrep("=", 70), "\n\n")

chm_results <- list()

for (lc in forested_classes) {
  lc_data <- chm_forest %>% filter(lc_l1_code == lc)
  
  if (nrow(lc_data) == 0) {
    cat(sprintf("  %s: No data, skipping\n", lc))
    next
  }
  
  cat(sprintf("  %s (%s): Predicting %s observations...\n",
              lc, lc_names[lc], format(nrow(lc_data), big.mark = ",")))
  
  t_start <- Sys.time()
  lc_pred <- colMeans(posterior_predict(fit_chm_s2, newdata = lc_data,
                                         ndraws = n_draws, allow_new_levels = TRUE))
  t_elapsed <- difftime(Sys.time(), t_start, units = "mins")
  
  acc <- compute_accuracy(lc_data$chm_error_mean, lc_pred)
  acc$ForestType <- lc
  acc$ForestName <- lc_names[lc]
  chm_results[[lc]] <- acc
  
  cat(sprintf("    RMSE=%.2f m, MAE=%.2f m, Bias=%+.2f m, R2=%.2f  (%.1f min)\n",
              acc$RMSE, acc$MAE, acc$Bias, acc$R2, as.numeric(t_elapsed)))
}

# All forest types combined
cat(sprintf("  All forest types: Predicting %s observations...\n",
            format(nrow(chm_forest), big.mark = ",")))
t_start <- Sys.time()
chm_all_pred <- colMeans(posterior_predict(fit_chm_s2, newdata = chm_forest,
                                            ndraws = n_draws, allow_new_levels = TRUE))
t_elapsed <- difftime(Sys.time(), t_start, units = "mins")

acc_all <- compute_accuracy(chm_forest$chm_error_mean, chm_all_pred)
acc_all$ForestType <- "All"
acc_all$ForestName <- "All Forest Types"
chm_results[["All"]] <- acc_all

cat(sprintf("    RMSE=%.2f m, MAE=%.2f m, Bias=%+.2f m, R2=%.2f  (%.1f min)\n",
            acc_all$RMSE, acc_all$MAE, acc_all$Bias, acc_all$R2, as.numeric(t_elapsed)))

# =====================================================================
# DTM FOREST-TYPE ACCURACY
# =====================================================================

cat("\n", strrep("=", 70), "\n")
cat("DTM FOREST-TYPE ACCURACY (FULL HOLDOUT)\n")
cat(strrep("=", 70), "\n\n")

dtm_results <- list()

for (lc in forested_classes) {
  lc_data <- dtm_forest %>% filter(lc_l1_code == lc)
  
  if (nrow(lc_data) == 0) {
    cat(sprintf("  %s: No data, skipping\n", lc))
    next
  }
  
  cat(sprintf("  %s (%s): Predicting %s observations...\n",
              lc, lc_names[lc], format(nrow(lc_data), big.mark = ",")))
  
  t_start <- Sys.time()
  lc_pred <- colMeans(posterior_predict(fit_dtm_s2, newdata = lc_data,
                                         ndraws = n_draws, allow_new_levels = TRUE))
  t_elapsed <- difftime(Sys.time(), t_start, units = "mins")
  
  acc <- compute_accuracy(lc_data$dtm_error_mean, lc_pred)
  acc$ForestType <- lc
  acc$ForestName <- lc_names[lc]
  dtm_results[[lc]] <- acc
  
  cat(sprintf("    RMSE=%.2f m, MAE=%.2f m, Bias=%+.2f m, R2=%.2f  (%.1f min)\n",
              acc$RMSE, acc$MAE, acc$Bias, acc$R2, as.numeric(t_elapsed)))
}

# All forest types combined
cat(sprintf("  All forest types: Predicting %s observations...\n",
            format(nrow(dtm_forest), big.mark = ",")))
t_start <- Sys.time()
dtm_all_pred <- colMeans(posterior_predict(fit_dtm_s2, newdata = dtm_forest,
                                            ndraws = n_draws, allow_new_levels = TRUE))
t_elapsed <- difftime(Sys.time(), t_start, units = "mins")

acc_all_dtm <- compute_accuracy(dtm_forest$dtm_error_mean, dtm_all_pred)
acc_all_dtm$ForestType <- "All"
acc_all_dtm$ForestName <- "All Forest Types"
dtm_results[["All"]] <- acc_all_dtm

cat(sprintf("    RMSE=%.2f m, MAE=%.2f m, Bias=%+.2f m, R2=%.2f  (%.1f min)\n",
            acc_all_dtm$RMSE, acc_all_dtm$MAE, acc_all_dtm$Bias, acc_all_dtm$R2,
            as.numeric(t_elapsed)))

# =====================================================================
# PRINT SUMMARY TABLES
# =====================================================================

cat("\n", strrep("=", 70), "\n")
cat("SUMMARY TABLES (FULL HOLDOUT, NO CAPS)\n")
cat(strrep("=", 70), "\n\n")

cat("--- Table 13: CHM Forest-Type Accuracy (Holdout) ---\n")
cat(sprintf("%-30s %10s %8s %8s %8s %6s\n", "Forest Type", "N", "RMSE", "MAE", "Bias", "R2"))
cat(strrep("-", 72), "\n")
for (lc in c(forested_classes, "All")) {
  r <- chm_results[[lc]]
  cat(sprintf("%-30s %10s %8.2f %8.2f %+8.2f %6.2f\n",
              r$ForestName, format(r$N, big.mark = ","),
              r$RMSE, r$MAE, r$Bias, r$R2))
}

cat("\n--- Table 14: DTM Forest-Type Accuracy (Holdout) ---\n")
cat(sprintf("%-30s %10s %8s %8s %8s %6s\n", "Forest Type", "N", "RMSE", "MAE", "Bias", "R2"))
cat(strrep("-", 72), "\n")
for (lc in c(forested_classes, "All")) {
  r <- dtm_results[[lc]]
  cat(sprintf("%-30s %10s %8.2f %8.2f %+8.2f %6.2f\n",
              r$ForestName, format(r$N, big.mark = ","),
              r$RMSE, r$MAE, r$Bias, r$R2))
}

# =====================================================================
# SAVE CSVs
# =====================================================================

chm_table <- bind_rows(chm_results) %>%
  select(ForestName, N, RMSE, MAE, Bias, R2)

dtm_table <- bind_rows(dtm_results) %>%
  select(ForestName, N, RMSE, MAE, Bias, R2)

write.csv(chm_table, file.path(out_tables, "chm_forest_type_accuracy_holdout_full.csv"),
          row.names = FALSE)
write.csv(dtm_table, file.path(out_tables, "dtm_forest_type_accuracy_holdout_full.csv"),
          row.names = FALSE)

cat("\n  Saved: chm_forest_type_accuracy_holdout_full.csv\n")
cat("  Saved: dtm_forest_type_accuracy_holdout_full.csv\n")

cat("\n", strrep("=", 70), "\n")
cat("SCRIPT COMPLETE\n")
cat(strrep("=", 70), "\n")
