# =====================================================================
# verify_predictive_accuracy.R
#
# PURPOSE:
#   1. Generate DTM forest-type accuracy table (matching CHM Table 13)
#   2. Verify CHM/DTM overall RMSE values (4.09 vs 4.12 discrepancy)
#   3. Clarify holdout vs training accuracy metrics
#
# RUN: source("verify_predictive_accuracy.R")
# =====================================================================

Sys.setenv(DISPLAY = "")
options(device = pdf)

suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
  library(tidyr)
})

cat("\n", strrep("=", 70), "\n")
cat("PREDICTIVE ACCURACY VERIFICATION\n")
cat(strrep("=", 70), "\n\n")

# =====================================================================
# CONFIGURATION
# =====================================================================

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/gpfs/data1/vclgp/lmaden/chpt1")
CHECKPOINT_DIR <- file.path(PROJECT_ROOT, "checkpoints")
out_tables <- file.path(PROJECT_ROOT, "tables")
dir.create(out_tables, showWarnings = FALSE, recursive = TRUE)

stage2_frac <- as.numeric(Sys.getenv("STAGE2_FRAC", "0.33"))
set.seed(2025)

forested_classes <- c("BDF", "DNF", "EBF", "ENF")
lc_names <- c(
  "BDF" = "Broadleaf Deciduous",
  "DNF" = "Deciduous Needleleaf",
  "EBF" = "Evergreen Broadleaf",
  "ENF" = "Evergreen Needleleaf"
)

# =====================================================================
# LOAD MODELS
# =====================================================================

cat("Loading models...\n")
model_data <- readRDS(file.path(CHECKPOINT_DIR, "10_models_stage2.rds"))
fit_chm_s2 <- model_data$data$fit_chm_s2
fit_dtm_s2 <- model_data$data$fit_dtm_s2
cat("  Models loaded\n\n")

# =====================================================================
# LOAD FULL DATA AND RECREATE HOLDOUT
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

# Recreate training/holdout split (same seed/fraction as section_08)
set.seed(2025)

chm_df_full$original_row <- 1:nrow(chm_df_full)
dtm_df_full$original_row <- 1:nrow(dtm_df_full)

chm_train <- chm_df_full %>% group_by(site) %>% sample_frac(stage2_frac) %>% ungroup()
dtm_train <- dtm_df_full %>% group_by(site) %>% sample_frac(stage2_frac) %>% ungroup()

# Get training factor levels
train_lc_chm <- levels(droplevels(chm_train$lc_l1_code))
train_site_chm <- levels(droplevels(chm_train$site))
train_eco_chm <- levels(droplevels(chm_train$ecoregion))

train_lc_dtm <- levels(droplevels(dtm_train$lc_l1_code))
train_site_dtm <- levels(droplevels(dtm_train$site))
train_eco_dtm <- levels(droplevels(dtm_train$ecoregion))

# Holdout sets
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

# Exclude Site 10 from DTM (matches model training)
dtm_holdout <- dtm_holdout %>% filter(site != "10")

cat(sprintf("  Training: %s CHM, %s DTM\n",
            format(nrow(chm_train), big.mark = ","),
            format(nrow(dtm_train), big.mark = ",")))
cat(sprintf("  Holdout:  %s CHM, %s DTM\n",
            format(nrow(chm_holdout), big.mark = ","),
            format(nrow(dtm_holdout), big.mark = ",")))

# =====================================================================
# HELPER: Compute accuracy metrics
# =====================================================================

compute_accuracy <- function(observed, predicted) {
  valid <- is.finite(observed) & is.finite(predicted)
  obs <- observed[valid]
  pred <- predicted[valid]
  
  data.frame(
    N = length(obs),
    RMSE = sqrt(mean((obs - pred)^2)),
    MAE = mean(abs(obs - pred)),
    Bias = mean(pred - obs),
    R2 = cor(obs, pred)^2,
    stringsAsFactors = FALSE
  )
}

# =====================================================================
# PART 1: VERIFY OVERALL RMSE (TRAINING vs HOLDOUT)
# =====================================================================

cat("\n", strrep("=", 70), "\n")
cat("PART 1: OVERALL RMSE VERIFICATION\n")
cat(strrep("=", 70), "\n\n")

# --- Training data predictions ---
cat("Computing training data predictions...\n")
cat("  CHM training predictions...\n")
chm_train_pred <- colMeans(posterior_predict(fit_chm_s2, ndraws = 100))
chm_train_obs <- fit_chm_s2$data$chm_error_mean

cat("  DTM training predictions...\n")
dtm_train_pred <- colMeans(posterior_predict(fit_dtm_s2, ndraws = 100))
dtm_train_obs <- fit_dtm_s2$data$dtm_error_mean

train_chm_acc <- compute_accuracy(chm_train_obs, chm_train_pred)
train_dtm_acc <- compute_accuracy(dtm_train_obs, dtm_train_pred)

cat("\n--- TRAINING DATA ACCURACY ---\n")
cat(sprintf("  CHM: N=%s, RMSE=%.2f m, MAE=%.2f m, Bias=%.2f m, R2=%.3f\n",
            format(train_chm_acc$N, big.mark = ","),
            train_chm_acc$RMSE, train_chm_acc$MAE, train_chm_acc$Bias, train_chm_acc$R2))
cat(sprintf("  DTM: N=%s, RMSE=%.2f m, MAE=%.2f m, Bias=%.2f m, R2=%.3f\n",
            format(train_dtm_acc$N, big.mark = ","),
            train_dtm_acc$RMSE, train_dtm_acc$MAE, train_dtm_acc$Bias, train_dtm_acc$R2))

# --- Holdout data predictions ---
# Subsample holdout for computational feasibility
cat("\nComputing holdout predictions (subsampled for speed)...\n")
n_holdout_sample <- 20000

set.seed(42)
chm_holdout_sub <- chm_holdout %>% slice_sample(n = min(n_holdout_sample, nrow(.)))
dtm_holdout_sub <- dtm_holdout %>% slice_sample(n = min(n_holdout_sample, nrow(.)))

cat(sprintf("  CHM holdout subsample: %d observations\n", nrow(chm_holdout_sub)))
cat(sprintf("  DTM holdout subsample: %d observations\n", nrow(dtm_holdout_sub)))

cat("  CHM holdout predictions...\n")
chm_hold_pred <- colMeans(posterior_predict(fit_chm_s2, newdata = chm_holdout_sub,
                                             ndraws = 100, allow_new_levels = TRUE))
chm_hold_obs <- chm_holdout_sub$chm_error_mean

cat("  DTM holdout predictions...\n")
dtm_hold_pred <- colMeans(posterior_predict(fit_dtm_s2, newdata = dtm_holdout_sub,
                                             ndraws = 100, allow_new_levels = TRUE))
dtm_hold_obs <- dtm_holdout_sub$dtm_error_mean

hold_chm_acc <- compute_accuracy(chm_hold_obs, chm_hold_pred)
hold_dtm_acc <- compute_accuracy(dtm_hold_obs, dtm_hold_pred)

cat("\n--- HOLDOUT DATA ACCURACY ---\n")
cat(sprintf("  CHM: N=%s, RMSE=%.2f m, MAE=%.2f m, Bias=%.2f m, R2=%.3f\n",
            format(hold_chm_acc$N, big.mark = ","),
            hold_chm_acc$RMSE, hold_chm_acc$MAE, hold_chm_acc$Bias, hold_chm_acc$R2))
cat(sprintf("  DTM: N=%s, RMSE=%.2f m, MAE=%.2f m, Bias=%.2f m, R2=%.3f\n",
            format(hold_dtm_acc$N, big.mark = ","),
            hold_dtm_acc$RMSE, hold_dtm_acc$MAE, hold_dtm_acc$Bias, hold_dtm_acc$R2))

cat("\n--- MANUSCRIPT VALUES ---\n")
cat("  Figure caption: CHM RMSE = 4.09 m, R2 = 0.44 | DTM RMSE = 3.03 m, R2 = 0.24\n")
cat("  Table 12:       CHM RMSE = 4.12 m           | DTM RMSE = 3.02 m\n")
cat("\n  Compare with values above to identify which is training vs holdout.\n")

# =====================================================================
# PART 2: DTM FOREST-TYPE ACCURACY
# =====================================================================

cat("\n", strrep("=", 70), "\n")
cat("PART 2: DTM FOREST-TYPE PREDICTIVE ACCURACY\n")
cat(strrep("=", 70), "\n\n")

# Filter holdout to forested classes
dtm_forest_holdout <- dtm_holdout %>%
  filter(lc_l1_code %in% forested_classes)

cat(sprintf("DTM forested holdout: %s observations\n",
            format(nrow(dtm_forest_holdout), big.mark = ",")))

# Sample per forest type
n_per_class <- 20000

dtm_forest_results <- list()

for (lc in forested_classes) {
  lc_data <- dtm_forest_holdout %>% filter(lc_l1_code == lc)
  
  if (nrow(lc_data) == 0) {
    cat(sprintf("  %s: No data, skipping\n", lc))
    next
  }
  
  # Subsample if necessary
  set.seed(42)
  lc_sub <- lc_data %>% slice_sample(n = min(n_per_class, nrow(.)))
  
  cat(sprintf("  %s (%s): Predicting %d observations...\n",
              lc, lc_names[lc], nrow(lc_sub)))
  
  lc_pred <- colMeans(posterior_predict(fit_dtm_s2, newdata = lc_sub,
                                         ndraws = 100, allow_new_levels = TRUE))
  lc_obs <- lc_sub$dtm_error_mean
  
  acc <- compute_accuracy(lc_obs, lc_pred)
  acc$ForestType <- lc
  acc$ForestName <- lc_names[lc]
  
  dtm_forest_results[[lc]] <- acc
  
  cat(sprintf("    RMSE=%.2f m, MAE=%.2f m, Bias=%+.2f m, R2=%.2f\n",
              acc$RMSE, acc$MAE, acc$Bias, acc$R2))
}

# All forest types combined
cat("  All forest types combined...\n")
set.seed(42)
dtm_all_forest_sub <- dtm_forest_holdout %>%
  slice_sample(n = min(n_per_class * 2, nrow(.)))

dtm_all_pred <- colMeans(posterior_predict(fit_dtm_s2, newdata = dtm_all_forest_sub,
                                            ndraws = 100, allow_new_levels = TRUE))
dtm_all_obs <- dtm_all_forest_sub$dtm_error_mean

acc_all <- compute_accuracy(dtm_all_obs, dtm_all_pred)
acc_all$ForestType <- "All"
acc_all$ForestName <- "All Forest Types"
dtm_forest_results[["All"]] <- acc_all

cat(sprintf("    RMSE=%.2f m, MAE=%.2f m, Bias=%+.2f m, R2=%.2f\n",
            acc_all$RMSE, acc_all$MAE, acc_all$Bias, acc_all$R2))

# =====================================================================
# ALSO RECOMPUTE CHM FOREST-TYPE FOR CONSISTENCY
# =====================================================================

cat("\n", strrep("=", 70), "\n")
cat("PART 3: CHM FOREST-TYPE PREDICTIVE ACCURACY (VERIFICATION)\n")
cat(strrep("=", 70), "\n\n")

chm_forest_holdout <- chm_holdout %>%
  filter(lc_l1_code %in% forested_classes)

cat(sprintf("CHM forested holdout: %s observations\n",
            format(nrow(chm_forest_holdout), big.mark = ",")))

chm_forest_results <- list()

for (lc in forested_classes) {
  lc_data <- chm_forest_holdout %>% filter(lc_l1_code == lc)
  
  if (nrow(lc_data) == 0) {
    cat(sprintf("  %s: No data, skipping\n", lc))
    next
  }
  
  set.seed(42)
  lc_sub <- lc_data %>% slice_sample(n = min(n_per_class, nrow(.)))
  
  cat(sprintf("  %s (%s): Predicting %d observations...\n",
              lc, lc_names[lc], nrow(lc_sub)))
  
  lc_pred <- colMeans(posterior_predict(fit_chm_s2, newdata = lc_sub,
                                         ndraws = 100, allow_new_levels = TRUE))
  lc_obs <- lc_sub$chm_error_mean
  
  acc <- compute_accuracy(lc_obs, lc_pred)
  acc$ForestType <- lc
  acc$ForestName <- lc_names[lc]
  
  chm_forest_results[[lc]] <- acc
  
  cat(sprintf("    RMSE=%.2f m, MAE=%.2f m, Bias=%+.2f m, R2=%.2f\n",
              acc$RMSE, acc$MAE, acc$Bias, acc$R2))
}

# All forest types combined
cat("  All forest types combined...\n")
set.seed(42)
chm_all_forest_sub <- chm_forest_holdout %>%
  slice_sample(n = min(n_per_class * 2, nrow(.)))

chm_all_pred <- colMeans(posterior_predict(fit_chm_s2, newdata = chm_all_forest_sub,
                                            ndraws = 100, allow_new_levels = TRUE))
chm_all_obs <- chm_all_forest_sub$chm_error_mean

acc_all_chm <- compute_accuracy(chm_all_obs, chm_all_pred)
acc_all_chm$ForestType <- "All"
acc_all_chm$ForestName <- "All Forest Types"
chm_forest_results[["All"]] <- acc_all_chm

cat(sprintf("    RMSE=%.2f m, MAE=%.2f m, Bias=%+.2f m, R2=%.2f\n",
            acc_all_chm$RMSE, acc_all_chm$MAE, acc_all_chm$Bias, acc_all_chm$R2))

# =====================================================================
# PRINT COMPARISON TABLES
# =====================================================================

cat("\n", strrep("=", 70), "\n")
cat("SUMMARY TABLES\n")
cat(strrep("=", 70), "\n\n")

# CHM table
cat("--- CHM Forest-Type Accuracy (Holdout) ---\n")
cat(sprintf("%-30s %8s %8s %8s %8s %6s\n", "Forest Type", "N", "RMSE", "MAE", "Bias", "R2"))
cat(strrep("-", 70), "\n")
for (lc in c(forested_classes, "All")) {
  r <- chm_forest_results[[lc]]
  cat(sprintf("%-30s %8s %8.2f %8.2f %+8.2f %6.2f\n",
              r$ForestName, format(r$N, big.mark = ","),
              r$RMSE, r$MAE, r$Bias, r$R2))
}

cat("\n--- CHM Manuscript Values (for comparison) ---\n")
cat(sprintf("%-30s %8s %8s %8s %8s %6s\n", "Forest Type", "N", "RMSE", "MAE", "Bias", "R2"))
cat(strrep("-", 70), "\n")
cat(sprintf("%-30s %8s %8.2f %8.2f %+8.2f %6.2f\n", "Broadleaf Deciduous", "55,744", 3.77, 2.81, -0.08, 0.57))
cat(sprintf("%-30s %8s %8.2f %8.2f %+8.2f %6.2f\n", "Deciduous Needleleaf", "14,575", 3.52, 2.67, 0.00, 0.56))
cat(sprintf("%-30s %8s %8.2f %8.2f %+8.2f %6.2f\n", "Evergreen Broadleaf", "2,912", 3.66, 2.75, -0.11, 0.50))
cat(sprintf("%-30s %8s %8.2f %8.2f %+8.2f %6.2f\n", "Evergreen Needleleaf", "81,333", 4.25, 2.96, 0.11, 0.28))
cat(sprintf("%-30s %8s %8.2f %8.2f %+8.2f %6.2f\n", "All Forest Types", "154,564", 4.01, 2.87, 0.03, 0.47))

# DTM table
cat("\n--- DTM Forest-Type Accuracy (Holdout) ---\n")
cat(sprintf("%-30s %8s %8s %8s %8s %6s\n", "Forest Type", "N", "RMSE", "MAE", "Bias", "R2"))
cat(strrep("-", 70), "\n")
for (lc in c(forested_classes, "All")) {
  r <- dtm_forest_results[[lc]]
  cat(sprintf("%-30s %8s %8.2f %8.2f %+8.2f %6.2f\n",
              r$ForestName, format(r$N, big.mark = ","),
              r$RMSE, r$MAE, r$Bias, r$R2))
}

# =====================================================================
# SAVE CSVs
# =====================================================================

chm_table <- bind_rows(chm_forest_results) %>%
  select(ForestName, N, RMSE, MAE, Bias, R2) %>%
  mutate(across(c(RMSE, MAE, Bias, R2), ~round(., 2)))

dtm_table <- bind_rows(dtm_forest_results) %>%
  select(ForestName, N, RMSE, MAE, Bias, R2) %>%
  mutate(across(c(RMSE, MAE, Bias, R2), ~round(., 2)))

write.csv(chm_table, file.path(out_tables, "chm_forest_type_accuracy_holdout.csv"),
          row.names = FALSE)
write.csv(dtm_table, file.path(out_tables, "dtm_forest_type_accuracy_holdout.csv"),
          row.names = FALSE)

cat("\n  Saved: chm_forest_type_accuracy_holdout.csv\n")
cat("  Saved: dtm_forest_type_accuracy_holdout.csv\n")

# =====================================================================
# OVERALL ACCURACY COMPARISON
# =====================================================================

cat("\n", strrep("=", 70), "\n")
cat("OVERALL ACCURACY: TRAINING vs HOLDOUT\n")
cat(strrep("=", 70), "\n\n")

cat(sprintf("%-12s %10s %10s %10s %10s\n", "", "RMSE(trn)", "RMSE(hld)", "R2(trn)", "R2(hld)"))
cat(strrep("-", 54), "\n")
cat(sprintf("%-12s %10.2f %10.2f %10.3f %10.3f\n",
            "CHM", train_chm_acc$RMSE, hold_chm_acc$RMSE,
            train_chm_acc$R2, hold_chm_acc$R2))
cat(sprintf("%-12s %10.2f %10.2f %10.3f %10.3f\n",
            "DTM", train_dtm_acc$RMSE, hold_dtm_acc$RMSE,
            train_dtm_acc$R2, hold_dtm_acc$R2))

cat("\nManuscript reports: CHM RMSE=4.09/4.12, DTM RMSE=3.03/3.02\n")
cat("Compare with training/holdout values above to resolve discrepancy.\n")

cat("\n", strrep("=", 70), "\n")
cat("SCRIPT COMPLETE\n")
cat(strrep("=", 70), "\n")
