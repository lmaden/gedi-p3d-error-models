#!/usr/bin/env Rscript
################################################################################
# Section 14: Multi-Site DTM Error Diagnostics
# Purpose: Identify patterns in DTM errors across all sites
# 
# This script:
# 1. Loads all site_XX_enriched.csv files
# 2. Computes site-level error statistics for CHM and DTM
# 3. Identifies problematic sites and patterns
# 4. Creates comparative visualizations
# 5. Tests for relationships with site characteristics
################################################################################

# Load tidyverse components individually
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(purrr)
library(stringr)
library(tibble)
library(scales)

# Configuration
DATA_DIR <- "/gpfs/data1/vclgp/lmaden/chpt1/data/enriched_by_site"
OUTPUT_DIR <- "/gpfs/data1/vclgp/lmaden/chpt1/diagnostics"

# Create output directory if it doesn't exist
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
setwd(OUTPUT_DIR)

cat("Multi-Site DTM Error Diagnostics\n")
cat("================================\n\n")
cat("Data directory:", DATA_DIR, "\n")
cat("Output directory:", OUTPUT_DIR, "\n\n")

################################################################################
# 1. LOAD ALL SITE DATA
################################################################################

cat("Loading site data...\n")

# Find all site CSV files
site_files <- list.files(DATA_DIR, pattern = "site_\\d+_enriched\\.csv(\\.gz)?$", 
                         full.names = TRUE)

cat(sprintf("Found %d site files\n", length(site_files)))

# Extract site IDs
site_ids <- as.integer(str_extract(basename(site_files), "\\d+"))

# Function to safely load a site and add site ID
load_site <- function(filepath, site_id) {
  tryCatch({
    df <- read_csv(filepath, show_col_types = FALSE)
    df$site <- site_id
    cat(sprintf("  Site %02d: %d observations\n", site_id, nrow(df)))
    return(df)
  }, error = function(e) {
    cat(sprintf("  ERROR loading Site %02d: %s\n", site_id, e$message))
    return(NULL)
  })
}

# Load all sites
all_data <- map2_dfr(site_files, site_ids, load_site)

cat(sprintf("\nTotal observations across all sites: %d\n", nrow(all_data)))
cat(sprintf("Number of sites loaded: %d\n", n_distinct(all_data$site)))

################################################################################
# 2. COMPUTE SITE-LEVEL STATISTICS
################################################################################

cat("\n" , rep("=", 80), "\n", sep = "")
cat("Computing Site-Level Statistics\n")
cat(rep("=", 80), "\n", sep = "")

# Helper function for MAD
mad_val <- function(x) {
  median(abs(x - median(x, na.rm = TRUE)), na.rm = TRUE)
}

# Compute comprehensive statistics by site
site_stats <- all_data %>%
  group_by(site) %>%
  summarise(
    # Sample size
    n_obs = n(),
    n_chm_valid = sum(!is.na(error_mean)),
    n_dtm_valid = sum(!is.na(dtm_error_mean)),
    n_complete_pairs = sum(complete.cases(error_mean, dtm_error_mean)),
    n_als_chm_missing = sum(is.na(als_chm_mean)),
    frac_als_chm_missing = n_als_chm_missing / n_obs,
    
    # CHM error statistics
    chm_mean = mean(error_mean, na.rm = TRUE),
    chm_median = median(error_mean, na.rm = TRUE),
    chm_sd = sd(error_mean, na.rm = TRUE),
    chm_mad = mad_val(error_mean),
    chm_iqr = IQR(error_mean, na.rm = TRUE),
    chm_range = max(error_mean, na.rm = TRUE) - min(error_mean, na.rm = TRUE),
    chm_min = min(error_mean, na.rm = TRUE),
    chm_max = max(error_mean, na.rm = TRUE),
    chm_skew = moments::skewness(error_mean, na.rm = TRUE),
    chm_kurtosis = moments::kurtosis(error_mean, na.rm = TRUE),
    
    # DTM error statistics  
    dtm_mean = mean(dtm_error_mean, na.rm = TRUE),
    dtm_median = median(dtm_error_mean, na.rm = TRUE),
    dtm_sd = sd(dtm_error_mean, na.rm = TRUE),
    dtm_mad = mad_val(dtm_error_mean),
    dtm_iqr = IQR(dtm_error_mean, na.rm = TRUE),
    dtm_range = max(dtm_error_mean, na.rm = TRUE) - min(dtm_error_mean, na.rm = TRUE),
    dtm_min = min(dtm_error_mean, na.rm = TRUE),
    dtm_max = max(dtm_error_mean, na.rm = TRUE),
    dtm_skew = moments::skewness(dtm_error_mean, na.rm = TRUE),
    dtm_kurtosis = moments::kurtosis(dtm_error_mean, na.rm = TRUE),
    
    # Count extreme errors (beyond 3*IQR)
    chm_extreme_low = sum(error_mean < (quantile(error_mean, 0.25, na.rm = TRUE) - 
                          3 * IQR(error_mean, na.rm = TRUE)), na.rm = TRUE),
    chm_extreme_high = sum(error_mean > (quantile(error_mean, 0.75, na.rm = TRUE) + 
                           3 * IQR(error_mean, na.rm = TRUE)), na.rm = TRUE),
    chm_extreme_total = chm_extreme_low + chm_extreme_high,
    chm_extreme_frac = chm_extreme_total / n_chm_valid,
    
    dtm_extreme_low = sum(dtm_error_mean < (quantile(dtm_error_mean, 0.25, na.rm = TRUE) - 
                          3 * IQR(dtm_error_mean, na.rm = TRUE)), na.rm = TRUE),
    dtm_extreme_high = sum(dtm_error_mean > (quantile(dtm_error_mean, 0.75, na.rm = TRUE) + 
                           3 * IQR(dtm_error_mean, na.rm = TRUE)), na.rm = TRUE),
    dtm_extreme_total = dtm_extreme_low + dtm_extreme_high,
    dtm_extreme_frac = dtm_extreme_total / n_dtm_valid,
    
    # Count catastrophic errors (>20m)
    dtm_catastrophic = sum(abs(dtm_error_mean) > 20, na.rm = TRUE),
    dtm_catastrophic_frac = dtm_catastrophic / n_dtm_valid,
    
    # Error correlation (safe version that handles no complete pairs)
    error_correlation = {
      complete_pairs <- complete.cases(error_mean, dtm_error_mean)
      if (sum(complete_pairs) >= 2) {
        cor(error_mean[complete_pairs], dtm_error_mean[complete_pairs])
      } else {
        NA_real_
      }
    },
    
    # Mean landscape characteristics
    mean_als_chm = mean(als_chm_mean, na.rm = TRUE),
    mean_slope = mean(slope_mean, na.rm = TRUE),
    mean_cover = mean(cover, na.rm = TRUE),
    
    # Get ecoregion (mode)
    ecoregion = names(sort(table(ecoregion), decreasing = TRUE))[1],
    
    .groups = 'drop'
  )

# Print summary table
cat("\nSite-Level Summary Statistics:\n")
print(site_stats %>% 
        select(site, n_obs, dtm_sd, dtm_kurtosis, dtm_catastrophic, 
               frac_als_chm_missing, ecoregion) %>%
        arrange(site), n = Inf)

# Write detailed statistics
write_csv(site_stats, "multisite_statistics.csv")
cat("\nSaved: multisite_statistics.csv\n")

################################################################################
# 3. IDENTIFY PROBLEMATIC SITES
################################################################################

cat("\n" , rep("=", 80), "\n", sep = "")
cat("Identifying Problematic Sites\n")
cat(rep("=", 80), "\n", sep = "")

# Flag problematic characteristics
problematic <- site_stats %>%
  mutate(
    # Flag 1: High DTM kurtosis (heavy tails)
    flag_high_kurtosis = dtm_kurtosis > quantile(dtm_kurtosis, 0.75, na.rm = TRUE),
    
    # Flag 2: High DTM variance
    flag_high_variance = dtm_sd > quantile(dtm_sd, 0.75, na.rm = TRUE),
    
    # Flag 3: Catastrophic errors present
    flag_catastrophic = dtm_catastrophic > 0,
    
    # Flag 4: High extreme outlier fraction
    flag_high_outliers = dtm_extreme_frac > quantile(dtm_extreme_frac, 0.75, na.rm = TRUE),
    
    # Flag 5: High missing reference data
    flag_missing_data = frac_als_chm_missing > 0.2,
    
    # Flag 6: DTM variance > CHM variance (unexpected)
    flag_dtm_worse = dtm_sd > chm_sd,
    
    # Total flags
    total_flags = flag_high_kurtosis + flag_high_variance + flag_catastrophic +
                  flag_high_outliers + flag_missing_data + flag_dtm_worse,
    
    # Classification
    severity = case_when(
      total_flags >= 5 ~ "Severe",
      total_flags >= 3 ~ "Moderate",
      total_flags >= 1 ~ "Minor",
      TRUE ~ "None"
    )
  )

# Summary of problematic sites
cat("\nSite Classifications:\n")
table(problematic$severity) %>% print()

cat("\nSites with issues (1+ flags):\n")
problematic %>%
  filter(total_flags > 0) %>%
  select(site, severity, total_flags, dtm_kurtosis, dtm_catastrophic, 
         frac_als_chm_missing, ecoregion) %>%
  arrange(desc(total_flags)) %>%
  print(n = Inf)

# Write flags
write_csv(problematic, "multisite_flags.csv")
cat("\nSaved: multisite_flags.csv\n")

################################################################################
# 4. ANALYZE PATTERNS
################################################################################

cat("\n" , rep("=", 80), "\n", sep = "")
cat("Analyzing Cross-Site Patterns\n")
cat(rep("=", 80), "\n", sep = "")

# Test for relationships between site characteristics and DTM problems

# 1. Missing data vs DTM error variance
cat("\n1. Missing Reference Data vs DTM Error Variance:\n")
cor_missing_var <- cor(problematic$frac_als_chm_missing, 
                        problematic$dtm_sd, 
                        use = "complete.obs")
cat(sprintf("   Correlation: %.3f\n", cor_missing_var))

# 2. Kurtosis distribution
cat("\n2. DTM Error Kurtosis Distribution:\n")
cat(sprintf("   Median: %.1f\n", median(problematic$dtm_kurtosis, na.rm = TRUE)))
cat(sprintf("   Q75: %.1f\n", quantile(problematic$dtm_kurtosis, 0.75, na.rm = TRUE)))
cat(sprintf("   Max: %.1f\n", max(problematic$dtm_kurtosis, na.rm = TRUE)))
cat(sprintf("   Sites with kurtosis > 10: %d\n", 
            sum(problematic$dtm_kurtosis > 10, na.rm = TRUE)))

# 3. Catastrophic errors summary
cat("\n3. Catastrophic Errors (|error| > 20m):\n")
cat(sprintf("   Total across all sites: %d\n", sum(problematic$dtm_catastrophic)))
cat(sprintf("   Sites with catastrophic errors: %d / %d\n",
            sum(problematic$dtm_catastrophic > 0),
            nrow(problematic)))
cat(sprintf("   Maximum at single site: %d\n", max(problematic$dtm_catastrophic)))

# 4. Ecoregion patterns
cat("\n4. Issues by Ecoregion:\n")
ecoregion_summary <- problematic %>%
  group_by(ecoregion) %>%
  summarise(
    n_sites = n(),
    n_problematic = sum(total_flags >= 3),
    mean_kurtosis = mean(dtm_kurtosis, na.rm = TRUE),
    mean_catastrophic = mean(dtm_catastrophic),
    .groups = 'drop'
  ) %>%
  arrange(desc(n_problematic))
print(ecoregion_summary, n = Inf)

################################################################################
# 5. VISUALIZATIONS
################################################################################

cat("\n" , rep("=", 80), "\n", sep = "")
cat("Creating Diagnostic Plots\n")
cat(rep("=", 80), "\n", sep = "")

# Plot 1: DTM vs CHM error characteristics
p1 <- ggplot(site_stats, aes(x = chm_sd, y = dtm_sd)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  geom_point(aes(size = n_obs, color = dtm_kurtosis), alpha = 0.7) +
  scale_color_viridis_c(option = "plasma", trans = "log10") +
  scale_size_continuous(range = c(2, 10)) +
  geom_text(aes(label = site), size = 3, nudge_x = 0.1, nudge_y = 0.1) +
  labs(
    title = "DTM vs CHM Error Variability Across Sites",
    subtitle = "Points above diagonal indicate higher DTM variance",
    x = "CHM Error SD (m)",
    y = "DTM Error SD (m)",
    color = "DTM Kurtosis\n(log scale)",
    size = "N obs"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

ggsave("multisite_01_dtm_vs_chm_variance.pdf", p1, width = 12, height = 8)
cat("Saved: multisite_01_dtm_vs_chm_variance.pdf\n")

# Plot 2: Kurtosis by site
p2 <- problematic %>%
  arrange(dtm_kurtosis) %>%
  mutate(site = factor(site, levels = site)) %>%
  ggplot(aes(x = site, y = dtm_kurtosis, fill = severity)) +
  geom_col(alpha = 0.8) +
  geom_hline(yintercept = 3, linetype = "dashed", color = "blue", linewidth = 0.8) +
  scale_fill_manual(values = c("None" = "#4daf4a", "Minor" = "#ffff33", 
                                "Moderate" = "#ff7f00", "Severe" = "#e41a1c")) +
  labs(
    title = "DTM Error Kurtosis by Site",
    subtitle = "Dashed line = Normal distribution (kurtosis = 3)",
    x = "Site",
    y = "Kurtosis",
    fill = "Severity"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("multisite_02_kurtosis_by_site.pdf", p2, width = 12, height = 6)
cat("Saved: multisite_02_kurtosis_by_site.pdf\n")

# Plot 3: Catastrophic errors
p3 <- problematic %>%
  filter(dtm_catastrophic > 0) %>%
  ggplot(aes(x = reorder(factor(site), dtm_catastrophic), 
             y = dtm_catastrophic)) +
  geom_col(aes(fill = frac_als_chm_missing), alpha = 0.8) +
  scale_fill_gradient(low = "#fee5d9", high = "#a50f15", 
                      labels = percent_format()) +
  coord_flip() +
  labs(
    title = "Catastrophic DTM Errors (|error| > 20m) by Site",
    x = "Site",
    y = "Count of Catastrophic Errors",
    fill = "Fraction Missing\nALS CHM"
  ) +
  theme_minimal()

ggsave("multisite_03_catastrophic_errors.pdf", p3, width = 10, height = 8)
cat("Saved: multisite_03_catastrophic_errors.pdf\n")

# Plot 4: Missing data vs error variance
p4 <- ggplot(site_stats, aes(x = frac_als_chm_missing, y = dtm_sd)) +
  geom_point(aes(size = n_obs, color = dtm_catastrophic), alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", linewidth = 1) +
  scale_x_continuous(labels = percent_format()) +
  scale_color_gradient(low = "#fee5d9", high = "#a50f15") +
  geom_text(aes(label = site), size = 3, nudge_y = 0.2) +
  labs(
    title = "Missing Reference Data vs DTM Error Variance",
    subtitle = sprintf("Correlation: r = %.3f", cor_missing_var),
    x = "Fraction of Observations with Missing ALS CHM",
    y = "DTM Error SD (m)",
    color = "Catastrophic\nErrors",
    size = "N obs"
  ) +
  theme_minimal()

ggsave("multisite_04_missing_data_vs_variance.pdf", p4, width = 12, height = 8)
cat("Saved: multisite_04_missing_data_vs_variance.pdf\n")

# Plot 5: Error distributions comparison
p5 <- site_stats %>%
  select(site, chm_sd, dtm_sd, dtm_kurtosis) %>%
  pivot_longer(cols = c(chm_sd, dtm_sd), 
               names_to = "product", 
               values_to = "error_sd") %>%
  mutate(product = ifelse(product == "chm_sd", "CHM", "DTM")) %>%
  ggplot(aes(x = error_sd, fill = product)) +
  geom_histogram(alpha = 0.6, bins = 20, position = "identity") +
  scale_fill_manual(values = c("CHM" = "#2166ac", "DTM" = "#b2182b")) +
  labs(
    title = "Distribution of Error SD Across Sites",
    x = "Error SD (m)",
    y = "Number of Sites",
    fill = "Product"
  ) +
  theme_minimal()

ggsave("multisite_05_error_sd_distributions.pdf", p5, width = 10, height = 6)
cat("Saved: multisite_05_error_sd_distributions.pdf\n")

################################################################################
# 6. DETAILED ANALYSIS OF WORST SITES
################################################################################

cat("\n" , rep("=", 80), "\n", sep = "")
cat("Analyzing Worst Sites in Detail\n")
cat(rep("=", 80), "\n", sep = "")

# Get top 5 worst sites by multiple criteria
worst_kurtosis <- problematic %>% 
  arrange(desc(dtm_kurtosis)) %>% 
  slice(1:5) %>% 
  pull(site)

worst_catastrophic <- problematic %>%
  filter(dtm_catastrophic > 0) %>%
  arrange(desc(dtm_catastrophic)) %>%
  slice(1:5) %>%
  pull(site)

worst_variance <- problematic %>%
  arrange(desc(dtm_sd)) %>%
  slice(1:5) %>%
  pull(site)

worst_sites <- unique(c(worst_kurtosis, worst_catastrophic, worst_variance))

cat(sprintf("\nIdentified %d worst sites for detailed analysis:\n", length(worst_sites)))
cat("Sites:", paste(worst_sites, collapse = ", "), "\n")

# Extract and save worst observations from these sites
worst_obs_list <- list()

for (site_id in worst_sites) {
  site_data <- all_data %>% filter(site == site_id)
  
  # Get top 10 worst DTM errors
  worst_neg <- site_data %>%
    filter(!is.na(dtm_error_mean)) %>%
    arrange(dtm_error_mean) %>%
    slice(1:10) %>%
    select(site, shot_number, als_chm_mean, p3d_chm_mean, error_mean,
           dep_dtm, mhrsi_dtm, dtm_error_mean, slope_mean, cover)
  
  worst_pos <- site_data %>%
    filter(!is.na(dtm_error_mean)) %>%
    arrange(desc(dtm_error_mean)) %>%
    slice(1:10) %>%
    select(site, shot_number, als_chm_mean, p3d_chm_mean, error_mean,
           dep_dtm, mhrsi_dtm, dtm_error_mean, slope_mean, cover)
  
  worst_obs_list[[as.character(site_id)]] <- bind_rows(worst_neg, worst_pos)
}

# Combine all worst observations
all_worst_obs <- bind_rows(worst_obs_list)
write_csv(all_worst_obs, "multisite_worst_observations.csv")
cat(sprintf("\nSaved worst observations from %d sites: multisite_worst_observations.csv\n", 
            length(worst_sites)))

################################################################################
# 7. SUMMARY REPORT
################################################################################

cat("\n" , rep("=", 80), "\n", sep = "")
cat("MULTI-SITE DIAGNOSTIC SUMMARY\n")
cat(rep("=", 80), "\n", sep = "")

cat(sprintf("\nTotal sites analyzed: %d\n", nrow(site_stats)))
cat(sprintf("Total observations: %d\n", sum(site_stats$n_obs)))

cat("\n--- Error Variance ---\n")
cat(sprintf("CHM error SD: %.3f ± %.3f m (mean ± SD across sites)\n",
            mean(site_stats$chm_sd, na.rm = TRUE),
            sd(site_stats$chm_sd, na.rm = TRUE)))
cat(sprintf("DTM error SD: %.3f ± %.3f m (mean ± SD across sites)\n",
            mean(site_stats$dtm_sd, na.rm = TRUE),
            sd(site_stats$dtm_sd, na.rm = TRUE)))
cat(sprintf("Sites where DTM SD > CHM SD: %d / %d (%.1f%%)\n",
            sum(site_stats$dtm_sd > site_stats$chm_sd, na.rm = TRUE),
            nrow(site_stats),
            100 * mean(site_stats$dtm_sd > site_stats$chm_sd, na.rm = TRUE)))

cat("\n--- Heavy Tails ---\n")
cat(sprintf("Median DTM kurtosis: %.1f\n", median(site_stats$dtm_kurtosis, na.rm = TRUE)))
cat(sprintf("Sites with DTM kurtosis > 10: %d / %d (%.1f%%)\n",
            sum(site_stats$dtm_kurtosis > 10, na.rm = TRUE),
            nrow(site_stats),
            100 * mean(site_stats$dtm_kurtosis > 10, na.rm = TRUE)))
cat(sprintf("Sites with DTM kurtosis > 20: %d / %d (%.1f%%)\n",
            sum(site_stats$dtm_kurtosis > 20, na.rm = TRUE),
            nrow(site_stats),
            100 * mean(site_stats$dtm_kurtosis > 20, na.rm = TRUE)))

cat("\n--- Catastrophic Errors ---\n")
cat(sprintf("Total catastrophic DTM errors (|error| > 20m): %d\n",
            sum(site_stats$dtm_catastrophic)))
cat(sprintf("Sites with catastrophic errors: %d / %d (%.1f%%)\n",
            sum(site_stats$dtm_catastrophic > 0),
            nrow(site_stats),
            100 * mean(site_stats$dtm_catastrophic > 0)))

cat("\n--- Missing Reference Data ---\n")
cat(sprintf("Mean fraction missing ALS CHM: %.1f%%\n",
            100 * mean(site_stats$frac_als_chm_missing)))
cat(sprintf("Correlation (missing data vs DTM SD): r = %.3f\n",
            cor_missing_var))

cat("\n--- Site Classifications ---\n")
severity_counts <- table(problematic$severity)
for (sev in names(severity_counts)) {
  cat(sprintf("%s: %d sites (%.1f%%)\n", 
              sev, 
              severity_counts[sev],
              100 * severity_counts[sev] / nrow(problematic)))
}

cat("\n--- Key Findings ---\n")

if (mean(site_stats$dtm_kurtosis > 10, na.rm = TRUE) > 0.5) {
  cat("⚠ CRITICAL: Majority of sites show extreme heavy tails (kurtosis > 10)\n")
  cat("  → This confirms that ν ≈ 2 pattern is widespread, not site-specific\n")
}

if (sum(site_stats$dtm_catastrophic) > 0) {
  cat(sprintf("⚠ WARNING: %d catastrophic errors detected across %d sites\n",
              sum(site_stats$dtm_catastrophic),
              sum(site_stats$dtm_catastrophic > 0)))
  cat("  → These extreme outliers dominate DTM error unpredictability\n")
}

if (abs(cor_missing_var) > 0.3) {
  cat(sprintf("⚠ PATTERN DETECTED: Missing reference data correlates with DTM variance (r=%.2f)\n",
              cor_missing_var))
  cat("  → Suggests reference data quality issues contribute to DTM unpredictability\n")
}

if (mean(site_stats$dtm_sd > site_stats$chm_sd, na.rm = TRUE) > 0.5) {
  cat("⚠ UNUSUAL: DTM error variance exceeds CHM at majority of sites\n")
  cat("  → Unexpected pattern suggesting potential reference data issues\n")
}

cat("\n" , rep("=", 80), "\n", sep = "")
cat("\nOutput Files Generated:\n")
cat("  - multisite_statistics.csv\n")
cat("  - multisite_flags.csv\n")
cat("  - multisite_worst_observations.csv\n")
cat("  - multisite_01_dtm_vs_chm_variance.pdf\n")
cat("  - multisite_02_kurtosis_by_site.pdf\n")
cat("  - multisite_03_catastrophic_errors.pdf\n")
cat("  - multisite_04_missing_data_vs_variance.pdf\n")
cat("  - multisite_05_error_sd_distributions.pdf\n")

cat("\n✓ Multi-site diagnostic analysis complete!\n")
cat(sprintf("All outputs saved to: %s\n", OUTPUT_DIR))
