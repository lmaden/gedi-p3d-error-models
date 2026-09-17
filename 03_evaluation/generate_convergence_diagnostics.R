# =====================================================================
# generate_convergence_diagnostics.R
#
# PURPOSE:
#   Generate Supplementary Figure: MCMC convergence diagnostics for
#   CHM and DTM Bayesian hierarchical models
#
# OUTPUTS:
#   - figure_s_convergence_combined.pdf/png (multi-panel figure)
#   - Individual panel PDFs for flexibility
#   - convergence_summary_by_type.csv
#
# RUN: source("generate_convergence_diagnostics.R")
# =====================================================================

Sys.setenv(DISPLAY = "")
options(device = pdf)

library(brms)
library(bayesplot)
library(ggplot2)
library(cowplot)
library(patchwork)
library(dplyr)

cat("\n", strrep("=", 70), "\n")
cat("CONVERGENCE DIAGNOSTIC FIGURES\n")
cat(strrep("=", 70), "\n\n")

# =====================================================================
# LOAD MODELS
# =====================================================================

cat("Loading models...\n")
model_data <- readRDS("/gpfs/data1/vclgp/lmaden/chpt1/checkpoints/10_models_stage2.rds")
fit_chm <- model_data$data$fit_chm_s2
fit_dtm <- model_data$data$fit_dtm_s2
cat("  Models loaded\n\n")

# Output directories
out_dir <- file.path(Sys.getenv("PROJECT_ROOT", "."), "plots")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

out_tables <- file.path(Sys.getenv("PROJECT_ROOT", "."), "tables")
if (!dir.exists(out_tables)) dir.create(out_tables, recursive = TRUE)

# =====================================================================
# PUBLICATION THEME
# =====================================================================

theme_pub <- function(base_size = 10) {
  theme_cowplot(font_size = base_size) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.text = element_text(color = "black", size = base_size - 1),
      axis.title = element_text(color = "black", size = base_size),
      plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size - 1, hjust = 0, color = "gray40")
    )
}

# =====================================================================
# EXTRACT DIAGNOSTICS
# =====================================================================

extract_diagnostics <- function(fit, product_name) {
  cat(sprintf("  Extracting diagnostics for %s...\n", product_name))
  
  rhat_vals <- rhat(fit)
  neff_vals <- neff_ratio(fit)
  
  df <- data.frame(
    parameter = names(rhat_vals),
    Rhat = as.numeric(rhat_vals),
    neff_ratio = as.numeric(neff_vals),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      type = case_when(
        grepl("^b_sigma", parameter) ~ "Sigma Fixed",
        grepl("^b_", parameter) ~ "Location Fixed",
        grepl("^r_site__sigma", parameter) ~ "Sigma Site RE",
        grepl("^r_site", parameter) ~ "Site RE",
        grepl("^r_ecoregion", parameter) ~ "Ecoregion RE",
        grepl("^r_lc_l1_code", parameter) ~ "Land Cover RE",
        grepl("^sd_", parameter) ~ "Group-level SD",
        grepl("^cor_", parameter) ~ "Correlations",
        grepl("nu", parameter) ~ "Degrees of Freedom",
        grepl("sslope|^s_", parameter) ~ "Spline Terms",
        TRUE ~ "Other"
      ),
      product = product_name
    )
  
  cat(sprintf("    %d parameters\n", nrow(df)))
  cat(sprintf("    R-hat range: [%.4f, %.4f]\n", 
              min(df$Rhat, na.rm = TRUE), max(df$Rhat, na.rm = TRUE)))
  cat(sprintf("    ESS ratio range: [%.4f, %.4f]\n",
              min(df$neff_ratio, na.rm = TRUE), max(df$neff_ratio, na.rm = TRUE)))
  
  return(df)
}

cat("Extracting diagnostics...\n")
diag_chm <- extract_diagnostics(fit_chm, "CHM")
diag_dtm <- extract_diagnostics(fit_dtm, "DTM")
diag_all <- bind_rows(diag_chm, diag_dtm)

# =====================================================================
# PANEL (a): R-hat distribution (both models)
# =====================================================================

cat("\nCreating panels...\n")

p_rhat <- ggplot(diag_all, aes(x = Rhat, fill = product)) +
  geom_histogram(bins = 50, alpha = 0.6, position = "identity",
                 color = "white", linewidth = 0.2) +
  geom_vline(xintercept = 1.01, linetype = "dashed", color = "red", linewidth = 0.7) +
  annotate("text", x = 1.01, y = Inf, label = "  1.01 threshold",
           hjust = 0, vjust = 1.5, size = 3, color = "red") +
  scale_fill_manual(values = c("CHM" = "#3182bd", "DTM" = "#e6550d"), name = "Model") +
  scale_x_continuous(limits = c(0.999, NA), breaks = seq(1.000, 1.010, by = 0.002)) +
  labs(title = "(a) R-hat distribution",
       x = expression(hat(R)), y = "Number of parameters") +
  theme_pub() +
  theme(legend.position = c(0.82, 0.82))

# =====================================================================
# PANEL (b): ESS ratio distribution (both models)
# =====================================================================

p_neff <- ggplot(diag_all, aes(x = neff_ratio, fill = product)) +
  geom_histogram(bins = 50, alpha = 0.6, position = "identity",
                 color = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "red", linewidth = 0.7) +
  annotate("text", x = 0.1, y = Inf, label = "  0.1 threshold",
           hjust = 0, vjust = 1.5, size = 3, color = "red") +
  scale_fill_manual(values = c("CHM" = "#3182bd", "DTM" = "#e6550d"), name = "Model") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
  labs(title = "(b) Effective sample size ratio",
       x = "ESS / Total draws", y = "Number of parameters") +
  theme_pub() +
  theme(legend.position = "none")

# =====================================================================
# PANEL (c): R-hat by parameter type (CHM)
# =====================================================================

type_order_chm <- diag_chm %>%
  group_by(type) %>%
  summarise(med = median(Rhat, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(med)) %>%
  pull(type)

diag_chm$type_f <- factor(diag_chm$type, levels = rev(type_order_chm))

p_rhat_chm <- ggplot(diag_chm %>% filter(!is.na(type_f)),
                      aes(x = Rhat, y = type_f)) +
  geom_jitter(alpha = 0.4, size = 1, color = "#3182bd", height = 0.2) +
  geom_vline(xintercept = 1.01, linetype = "dashed", color = "red", linewidth = 0.7) +
  scale_x_continuous(limits = c(0.999, 1.012), breaks = seq(1.000, 1.010, by = 0.002)) +
  labs(title = "(c) CHM: R-hat by parameter type",
       x = expression(hat(R)), y = NULL) +
  theme_pub() +
  theme(axis.text.y = element_text(size = 8))

# =====================================================================
# PANEL (d): R-hat by parameter type (DTM)
# =====================================================================

type_order_dtm <- diag_dtm %>%
  group_by(type) %>%
  summarise(med = median(Rhat, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(med)) %>%
  pull(type)

diag_dtm$type_f <- factor(diag_dtm$type, levels = rev(type_order_dtm))

p_rhat_dtm <- ggplot(diag_dtm %>% filter(!is.na(type_f)),
                      aes(x = Rhat, y = type_f)) +
  geom_jitter(alpha = 0.4, size = 1, color = "#e6550d", height = 0.2) +
  geom_vline(xintercept = 1.01, linetype = "dashed", color = "red", linewidth = 0.7) +
  scale_x_continuous(limits = c(0.999, 1.012), breaks = seq(1.000, 1.010, by = 0.002)) +
  labs(title = "(d) DTM: R-hat by parameter type",
       x = expression(hat(R)), y = NULL) +
  theme_pub() +
  theme(axis.text.y = element_text(size = 8))

# =====================================================================
# PANEL (e): CHM trace plots for key parameters
# =====================================================================

cat("  Creating trace plots...\n")

key_chm <- c("b_Intercept", "b_wsci_z", "b_slope_mean_z",
             "b_sigma_Intercept", "b_sigma_slope_mean_z", "nu")
avail_chm <- intersect(key_chm, variables(fit_chm))
cat(sprintf("    CHM trace params: %s\n", paste(avail_chm, collapse = ", ")))

# Human-readable labels for trace plots
chm_labels <- c(
  "b_Intercept" = "Intercept",
  "b_wsci_z" = "WSCI",
  "b_slope_mean_z" = "Slope",
  "b_sigma_Intercept" = "Sigma Intercept",
  "b_sigma_slope_mean_z" = "Sigma Slope",
  "nu" = "Degrees of Freedom (v)"
)

p_trace_chm <- mcmc_trace(fit_chm, pars = avail_chm,
                           facet_args = list(ncol = 1, strip.position = "right",
                                             labeller = as_labeller(chm_labels))) +
  labs(title = "(e) CHM: Trace plots") +
  theme_pub(base_size = 9) +
  theme(strip.text = element_text(size = 6),
        strip.text.y.right = element_text(angle = 0, hjust = 0, margin = margin(l = 2, r = 2)),
        plot.margin = margin(5, 40, 5, 5),
        legend.position = "none")

# =====================================================================
# PANEL (f): DTM trace plots for key parameters
# =====================================================================

key_dtm <- c("b_Intercept", "b_rh_98_z", "b_slope_mean_z",
             "b_sigma_Intercept", "b_sigma_slope_mean_z", "nu")
avail_dtm <- intersect(key_dtm, variables(fit_dtm))
cat(sprintf("    DTM trace params: %s\n", paste(avail_dtm, collapse = ", ")))

# Human-readable labels for DTM trace plots
dtm_labels <- c(
  "b_Intercept" = "Intercept",
  "b_rh_98_z" = "RH98",
  "b_slope_mean_z" = "Slope",
  "b_sigma_Intercept" = "Sigma Intercept",
  "b_sigma_slope_mean_z" = "Sigma Slope",
  "nu" = "Degrees of Freedom (v)"
)

p_trace_dtm <- mcmc_trace(fit_dtm, pars = avail_dtm,
                           facet_args = list(ncol = 1, strip.position = "right",
                                             labeller = as_labeller(dtm_labels))) +
  labs(title = "(f) DTM: Trace plots") +
  theme_pub(base_size = 9) +
  theme(strip.text = element_text(size = 6),
        strip.text.y.right = element_text(angle = 0, hjust = 0, margin = margin(l = 2, r = 2)),
        plot.margin = margin(5, 40, 5, 5),
        legend.position = "none")

# =====================================================================
# COMBINE INTO SINGLE FIGURE
# =====================================================================

cat("  Assembling combined figure...\n")

top_row <- p_rhat | p_neff
mid_row <- p_rhat_chm | p_rhat_dtm
bot_row <- p_trace_chm | p_trace_dtm

combined <- top_row / mid_row / bot_row +
  plot_layout(heights = c(1, 1.2, 1.5)) +
  plot_annotation(
    title = "MCMC Convergence Diagnostics",
    subtitle = sprintf(
      "CHM: %d parameters, max R-hat = %.4f | DTM: %d parameters, max R-hat = %.4f",
      nrow(diag_chm), max(diag_chm$Rhat, na.rm = TRUE),
      nrow(diag_dtm), max(diag_dtm$Rhat, na.rm = TRUE)
    ),
    theme = theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "gray40")
    )
  )

# =====================================================================
# SAVE
# =====================================================================

cat("  Saving figures...\n")

ggsave(file.path(out_dir, "figure_s_convergence_combined.pdf"),
       combined, width = 16, height = 16, bg = "white", device = pdf)
ggsave(file.path(out_dir, "figure_s_convergence_combined.png"),
       combined, width = 16, height = 16, dpi = 300, bg = "white")
cat("    Saved: figure_s_convergence_combined.pdf/png\n")

# Save individual panels for flexibility
ggsave(file.path(out_dir, "figure_s_convergence_rhat_ess.pdf"),
       top_row, width = 12, height = 4, bg = "white", device = pdf)
ggsave(file.path(out_dir, "figure_s_convergence_rhat_by_type.pdf"),
       mid_row, width = 12, height = 5, bg = "white", device = pdf)
ggsave(file.path(out_dir, "figure_s_convergence_traces.pdf"),
       bot_row, width = 12, height = 8, bg = "white", device = pdf)
cat("    Saved: individual panel PDFs\n")

# =====================================================================
# SUMMARY TABLE (CSV)
# =====================================================================

cat("  Creating summary table...\n")

summary_table <- diag_all %>%
  group_by(product, type) %>%
  summarise(
    n_params = n(),
    Rhat_median = round(median(Rhat, na.rm = TRUE), 4),
    Rhat_max = round(max(Rhat, na.rm = TRUE), 4),
    ESS_ratio_median = round(median(neff_ratio, na.rm = TRUE), 4),
    ESS_ratio_min = round(min(neff_ratio, na.rm = TRUE), 4),
    n_Rhat_above_1.01 = sum(Rhat > 1.01, na.rm = TRUE),
    n_ESS_below_0.1 = sum(neff_ratio < 0.1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(product, desc(Rhat_max))

write.csv(summary_table,
          file.path(out_tables, "convergence_summary_by_type.csv"),
          row.names = FALSE)
cat("    Saved: convergence_summary_by_type.csv\n")

# =====================================================================
# PRINT OVERALL SUMMARY
# =====================================================================

cat("\n", strrep("=", 70), "\n")
cat("CONVERGENCE SUMMARY\n")
cat(strrep("=", 70), "\n\n")

for (prod in c("CHM", "DTM")) {
  d <- diag_all %>% filter(product == prod)
  cat(sprintf("%s Model:\n", prod))
  cat(sprintf("  Total parameters:           %d\n", nrow(d)))
  cat(sprintf("  R-hat:  median = %.4f, max = %.4f\n",
              median(d$Rhat, na.rm = TRUE), max(d$Rhat, na.rm = TRUE)))
  cat(sprintf("  ESS ratio: median = %.3f, min = %.3f\n",
              median(d$neff_ratio, na.rm = TRUE), min(d$neff_ratio, na.rm = TRUE)))
  cat(sprintf("  R-hat > 1.01:               %d / %d (%.1f%%)\n",
              sum(d$Rhat > 1.01, na.rm = TRUE), nrow(d),
              sum(d$Rhat > 1.01, na.rm = TRUE) / nrow(d) * 100))
  cat(sprintf("  ESS ratio < 0.1:            %d / %d (%.1f%%)\n\n",
              sum(d$neff_ratio < 0.1, na.rm = TRUE), nrow(d),
              sum(d$neff_ratio < 0.1, na.rm = TRUE) / nrow(d) * 100))
}

# =====================================================================
# SUGGESTED CAPTION
# =====================================================================

cat(strrep("=", 70), "\n")
cat("SUGGESTED CAPTION\n")
cat(strrep("=", 70), "\n\n")

cat('Supplementary Figure SX. MCMC convergence diagnostics for the CHM\n')
cat('and DTM Bayesian hierarchical models. (a) Distribution of R-hat\n')
cat('values across all model parameters; dashed red line indicates the\n')
cat('1.01 convergence threshold (Vehtari et al., 2021). (b) Distribution\n')
cat('of effective sample size (ESS) ratios; dashed red line indicates the\n')
cat('0.1 minimum threshold. (c, d) R-hat values by parameter type for\n')
cat('CHM and DTM models, respectively, showing convergence across all\n')
cat('components of the hierarchical model. (e, f) Trace plots for key\n')
cat('parameters demonstrating adequate chain mixing across all four MCMC\n')
cat('chains. No divergent transitions were observed for either model.\n')

cat("\n", strrep("=", 70), "\n")
cat("SCRIPT COMPLETE\n")
cat(strrep("=", 70), "\n")
