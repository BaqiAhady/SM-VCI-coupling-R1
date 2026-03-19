# ============================================================
# Title: Identifying soil moisture thresholds for agricultural drought sensitivity 
#        through coupling modeled and satellite data in a humid temperate catchment
# Purpose: SM-VCI Threshold optimization and validation
# Authors: Abdul Baqi Ahady, Stefanie Wolf, Elena-Maria Klopries
# Date: March 2026
# ============================================================

# Load required libraries
library(dplyr)
library(ggplot2)
library(sf)
library(patchwork)

# ----------------------------------------------------------------------------
# STEP 1: LOAD AND PREPARE DATA
# ----------------------------------------------------------------------------

# Load data
data <- read.csv("Data/SM_with_VCI_00-24.csv")

# Clean and prepare data
data_clean <- data %>%
  mutate(
    Date = as.Date(Date),
    Month = as.integer(Month),
    SoilMoisture = as.numeric(gsub(" % nFK,", "", SoilMoisture))
  ) %>%
  mutate(
    SoilMoisture = na_if(SoilMoisture, -9999),
    VCI = na_if(VCI, -9999)
  ) %>%
  filter(Month >= 4 & Month <= 9, !is.na(VCI))

# Add Year column
data_clean <- data_clean %>%
  mutate(Year = as.numeric(format(Date, "%Y")))

# ----------------------------------------------------------------------------
# STEP 2: BASELINE CORRELATIONS
# ----------------------------------------------------------------------------

# Calculate baseline correlations
original_corr_pearson <- cor(data_clean$SoilMoisture, data_clean$VCI, 
                             use = "complete.obs", method = "pearson")
original_corr_spearman <- cor(data_clean$SoilMoisture, data_clean$VCI, 
                              use = "complete.obs", method = "spearman")

cat("\n=== BASELINE CORRELATIONS ===\n")
cat("Baseline Pearson r =", round(original_corr_pearson, 3), "\n")
cat("Baseline Spearman ρ =", round(original_corr_spearman, 3), "\n")

# ----------------------------------------------------------------------------
# STEP 3: OPTIMAL THRESHOLD DETERMINATION
# ----------------------------------------------------------------------------

# Test multiple soil moisture thresholds
percentiles <- seq(10, 90, by = 5)
threshold_results_pearson <- data.frame()
threshold_results_spearman <- data.frame()

for(p in percentiles) {
  threshold <- quantile(data_clean$SoilMoisture, p/100, na.rm = TRUE)
  data_filtered <- data_clean %>% filter(SoilMoisture <= threshold)
  
  if(nrow(data_filtered) > 100) {
    # Pearson correlation
    cor_val_pearson <- cor(data_filtered$SoilMoisture, data_filtered$VCI, 
                           use = "complete.obs", method = "pearson")
    
    threshold_results_pearson <- rbind(threshold_results_pearson,
                                       data.frame(percentile = p,
                                                  threshold = threshold,
                                                  n_points = nrow(data_filtered),
                                                  proportion = nrow(data_filtered)/nrow(data_clean)*100,
                                                  correlation = cor_val_pearson,
                                                  method = "Pearson"))
    
    # Spearman correlation
    cor_val_spearman <- cor(data_filtered$SoilMoisture, data_filtered$VCI, 
                            use = "complete.obs", method = "spearman")
    
    threshold_results_spearman <- rbind(threshold_results_spearman,
                                        data.frame(percentile = p,
                                                   threshold = threshold,
                                                   n_points = nrow(data_filtered),
                                                   proportion = nrow(data_filtered)/nrow(data_clean)*100,
                                                   correlation = cor_val_spearman,
                                                   method = "Spearman"))
  }
}

# Identify optimal thresholds
optimal_pearson <- threshold_results_pearson[which.max(threshold_results_pearson$correlation), ]
optimal_spearman <- threshold_results_spearman[which.max(threshold_results_spearman$correlation), ]

cat("\n=== OPTIMAL THRESHOLDS ===\n")
cat("Optimal Pearson threshold:", round(optimal_pearson$threshold, 1), 
    "% nFK (r =", round(optimal_pearson$correlation, 3), ")\n")
cat("Optimal Spearman threshold:", round(optimal_spearman$threshold, 1), 
    "% nFK (ρ =", round(optimal_spearman$correlation, 3), ")\n")

# ----------------------------------------------------------------------------
# STEP 4: LEAVE-ONE-YEAR-OUT CROSS-VALIDATION
# ----------------------------------------------------------------------------

cat("\n=== LEAVE-ONE-YEAR-OUT CROSS-VALIDATION RESULTS ===\n")

years <- unique(data_clean$Year)
loo_results <- data.frame()
all_loo_curves <- data.frame()

for(y in years) {
  # Train on all years except y
  train_data <- data_clean %>% filter(Year != y)
  
  # Test threshold optimization on training data
  loo_percentiles <- seq(10, 90, by = 5)
  loo_threshold_results <- data.frame()
  
  for(p in loo_percentiles) {
    threshold <- quantile(train_data$SoilMoisture, p/100, na.rm = TRUE)
    data_filtered <- train_data %>% filter(SoilMoisture <= threshold)
    
    if(nrow(data_filtered) > 100) {
      cor_val <- cor(data_filtered$SoilMoisture, data_filtered$VCI, 
                     use = "complete.obs", method = "pearson")
      
      loo_threshold_results <- rbind(loo_threshold_results,
                                     data.frame(percentile = p,
                                                threshold = threshold,
                                                correlation = cor_val))
    }
  }
  
  # Store full curve for this iteration
  loo_curve <- loo_threshold_results %>%
    mutate(left_out_year = y)
  all_loo_curves <- rbind(all_loo_curves, loo_curve)
  
  # Find optimal threshold for training set
  loo_optimal <- loo_threshold_results[which.max(loo_threshold_results$correlation), ]
  
  # Test on left-out year
  test_data <- data_clean %>% filter(Year == y)
  test_filtered <- test_data %>% filter(SoilMoisture <= loo_optimal$threshold)
  
  if(nrow(test_filtered) > 10) {
    test_cor <- cor(test_filtered$SoilMoisture, test_filtered$VCI, 
                    use = "complete.obs", method = "pearson")
  } else {
    test_cor <- NA
  }
  
  loo_results <- rbind(loo_results,
                       data.frame(left_out_year = y,
                                  train_threshold = loo_optimal$threshold,
                                  train_correlation = loo_optimal$correlation,
                                  test_correlation = test_cor,
                                  n_test_retained = nrow(test_filtered)))
  
  # Print individual year results
  cat(sprintf("Year %d: threshold = %.1f %% nFK, train r = %.3f, test r = %.3f (n = %d)\n", 
              y, loo_optimal$threshold, loo_optimal$correlation, 
              ifelse(is.na(test_cor), 0, test_cor), nrow(test_filtered)))
}

# LOOCV summary statistics
loo_mean <- mean(loo_results$train_threshold, na.rm = TRUE)
loo_sd <- sd(loo_results$train_threshold, na.rm = TRUE)
loo_min <- min(loo_results$train_threshold, na.rm = TRUE)
loo_max <- max(loo_results$train_threshold, na.rm = TRUE)
mean_test_cor <- mean(loo_results$test_correlation, na.rm = TRUE)
sd_test_cor <- sd(loo_results$test_correlation, na.rm = TRUE)

cat("\n=== LOOCV SUMMARY ===\n")
cat("Training thresholds:\n")
cat(sprintf("  Mean: %.1f %% nFK (SD = %.1f)\n", loo_mean, loo_sd))
cat(sprintf("  Range: %.1f - %.1f %% nFK\n", loo_min, loo_max))
cat(sprintf("  CV: %.1f%%\n", (loo_sd/loo_mean)*100))
cat("\nTest correlations:\n")
cat(sprintf("  Mean r = %.3f (SD = %.3f)\n", mean_test_cor, sd_test_cor))
cat(sprintf("  Range: %.3f - %.3f\n", 
            min(loo_results$test_correlation, na.rm = TRUE),
            max(loo_results$test_correlation, na.rm = TRUE)))

# Compare with full dataset optimum
cat("\n=== COMPARISON WITH FULL DATASET ===\n")
cat(sprintf("Full dataset optimum: %.1f %% nFK\n", optimal_pearson$threshold))
cat(sprintf("LOOCV mean: %.1f %% nFK\n", loo_mean))
cat(sprintf("Difference: %.1f %% nFK\n", abs(optimal_pearson$threshold - loo_mean)))
cat(sprintf("Full dataset optimum within 1 SD of LOOCV mean: %s\n", 
            ifelse(abs(optimal_pearson$threshold - loo_mean) <= loo_sd, "YES", "NO")))

# ----------------------------------------------------------------------------
# STEP 5: GENERATE PLOTS FOR PUBLICATION
# ----------------------------------------------------------------------------

# Shared arrow style for both plots
open_arrow <- arrow(angle = 30, length = unit(0.15, "inches"), 
                    ends = "last", type = "open")

# ----------------------------------------------------------------------------
# PLOT 1: LOO THRESHOLD DISTRIBUTION
# ----------------------------------------------------------------------------

p1 <- ggplot(loo_results, aes(x = train_threshold)) +
  
  geom_histogram(binwidth = 1, fill = "steelblue", color = "white", alpha = 0.7) +
  
  geom_vline(aes(xintercept = optimal_pearson$threshold, color = "Full dataset"),
             linetype = "dashed", linewidth = 0.7) +
  
  geom_vline(aes(xintercept = loo_mean, color = "LOO mean"),
             linetype = "solid", linewidth = 0.7) +
  
  scale_color_manual(
    name = "Threshold estimate",
    values = c("Full dataset" = "red", "LOO mean" = "darkblue"),
    labels = c(
      paste0("Full dataset (", round(optimal_pearson$threshold, 0), "% nFK)"),
      paste0("LOO mean (", round(loo_mean, 1), "% nFK)")
    )
  ) +
  
  labs(x = "Optimal Threshold (% nFK)",
       y = "Frequency (No. of LOO iterations)") +
  
  scale_x_continuous(
    limits = c(60, 66),
    breaks = seq(59, 65, by = 2),
    expand = c(0.02, 0)
  ) +
  
  scale_y_continuous(
    limits = c(0, 11),
    breaks = seq(0, 10, by = 2),
    expand = c(0.02, 0)
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    axis.line.x = element_line(linewidth = 0.5, color = "black", arrow = open_arrow),
    axis.line.y = element_line(linewidth = 0.5, color = "black", arrow = open_arrow),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    legend.position = c(0.50, 0.90),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = scales::alpha("white", 0.85),
                                     color = "gray70"),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10)
  )

# ----------------------------------------------------------------------------
# PLOT 2: THRESHOLD STABILITY WITH PEARSON CURVE
# ----------------------------------------------------------------------------

# Calculate confidence bands from LOOCV
loo_mean_threshold <- mean(loo_results$train_threshold, na.rm = TRUE)
loo_sd_threshold <- sd(loo_results$train_threshold, na.rm = TRUE)
loo_lower <- loo_mean_threshold - loo_sd_threshold
loo_upper <- loo_mean_threshold + loo_sd_threshold

p2 <- ggplot() +
  
  # LOOCV confidence band
  annotate("rect", 
           xmin = loo_lower, xmax = loo_upper, 
           ymin = -Inf, ymax = Inf,
           fill = "gray60", alpha = 0.3) +
  
  # Pearson correlation curve
  geom_line(data = threshold_results_pearson, 
            aes(x = threshold, y = correlation),
            color = "black", linewidth = 0.8) +
  
  geom_point(data = threshold_results_pearson, 
             aes(x = threshold, y = correlation),
             color = "black", size = 2) +
  
  # Optimal point
  geom_point(data = optimal_pearson, 
             aes(x = threshold, y = correlation), 
             size = 3.5, color = "red", shape = 19) +
  
  # Reference lines
  geom_hline(yintercept = original_corr_pearson, 
             linetype = "dashed", color = "gray50", linewidth = 0.6) +
  
  geom_vline(xintercept = optimal_pearson$threshold, 
             color = "red", linetype = "dotted", linewidth = 0.8) +
  
  geom_vline(xintercept = loo_mean_threshold, 
             color = "blue", linetype = "dotted", linewidth = 0.8) +
  
  # Baseline annotation
  annotate("text",
           x = min(threshold_results_pearson$threshold) + 2,
           y = original_corr_pearson,
           label = paste("Baseline r =", round(original_corr_pearson, 3)),
           color = "gray30", size = 4, vjust = -0.5, hjust = 0) +
  
  labs(x = expression("Soil Moisture Threshold (% nFK)"),
       y = expression("Pearson Correlation Coefficient (r)")) +
  
  scale_x_continuous(
    limits = range(threshold_results_pearson$threshold),
    breaks = seq(floor(min(threshold_results_pearson$threshold)/10)*10,
                 ceiling(max(threshold_results_pearson$threshold)/10)*10,
                 by = 5)
  ) +
  
  scale_y_continuous(
    limits = c(0, max(threshold_results_pearson$correlation) * 1.05),
    breaks = scales::pretty_breaks(n = 6)
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    axis.line.x = element_line(linewidth = 0.5, color = "black", arrow = open_arrow),
    axis.line.y = element_line(linewidth = 0.5, color = "black", arrow = open_arrow),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, face = "bold", color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(5, 5, 10, 5), "mm")
  )

# ----------------------------------------------------------------------------
# STEP 6: COMBINE PLOTS
# ----------------------------------------------------------------------------

combined_plot <- p2 + p1 + 
  plot_layout(ncol = 2, widths = c(0.6, 0.4)) +
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(",
                  tag_suffix = ")") &
  theme(plot.tag = element_text(size = 16, face = "bold", 
                                hjust = 0.5, vjust = 1),
        plot.tag.position = c(0.02, 0.98))

# Display combined plot
print(combined_plot)

# Save plot (optional - adjust dimensions as needed)
# ggsave("Figure_X_Threshold_Analysis.tiff", combined_plot, 
#        width = 10, height = 5, dpi = 300)