# ============================================================
# Title: Revealing agricultural drought sensitivity through a soil moisture
#        threshold-based framework integrating modeled and satellite data
# Purpose: SM threshold-based framework definition
# Authors: Abdul Baqi Ahady, Stefanie Wolf, Elena-Maria Klopries
# Date: Jan 2026
# ============================================================


# ------------------------------ Required Packages -----------------------------
library(dplyr)
library(ggplot2)
library(patchwork)
library(scales)

# ------------------------------ Loading the data ------------------------------

data <- read.csv("Data/SM_with_VCI_00-24.csv")

# the data is cleaned to make sure the -9999 
# values are excluded and is filtered for the growing season
data_clean <- data %>%
  mutate(
    Date = as.Date(Date),
    Month = as.integer(Month),
    SoilMoisture = as.numeric(gsub(" % nFK,", "", SoilMoisture)),
    SoilMoisture = na_if(SoilMoisture, -9999),
    VCI = na_if(VCI, -9999)
  ) %>%
  filter(Month >= 4 & Month <= 9, !is.na(VCI))


# -------------------- Optimal threshold determination -------------------------
# We calculate the baseline correlation between SM and VCI
original_corr <- cor(data_clean$SoilMoisture, data_clean$VCI, use = "complete.obs")
cat("Baseline SM–VCI correlation:", round(original_corr, 3), "\n")

# We introduce the percentiles with 10 step
percentiles <- seq(10, 90, by = 10)
threshold_results <- data.frame()

for (p in percentiles) {
  
  threshold <- quantile(data_clean$SoilMoisture, p / 100, na.rm = TRUE)
  data_tmp <- data_clean %>% filter(SoilMoisture <= threshold)
  
  if (nrow(data_tmp) > 100) {
    threshold_results <- rbind(
      threshold_results,
      data.frame(
        percentile  = p,
        threshold   = threshold,
        n_points    = nrow(data_tmp),
        proportion  = nrow(data_tmp) / nrow(data_clean) * 100,
        correlation = cor(data_tmp$SoilMoisture, data_tmp$VCI,
                          use = "complete.obs", method = "pearson")
      )
    )
  }
}
# Here we chose the threshold value aligns with the maximum correlation. 
optimal <- threshold_results[which.max(threshold_results$correlation), ]

# the results are printed here
cat("Optimal threshold:", round(optimal$threshold, 1), "% nFK\n")
cat("Correlation improvement:",
    round(original_corr, 3), "→", round(optimal$correlation, 3), "\n")


# -------------------- Appliying optimal threshold -----------------------------
# We filter out the data based on the optimal threshold
data_filtered <- data_clean %>%
  filter(SoilMoisture <= optimal$threshold)

# and see here how many observations are retained
cat("Total Observations:", nrow(data_clean), "\n")
cat("Retained observations: ", round(optimal$proportion, 1),
          "% (", optimal$n_points, " observations)\n")


# -------------------- Threshold optimization visualization --------------------
# setting for axes of the plots
open_arrow <- arrow(angle = 30, length = unit(0.15, "inches"),
                    ends = "last", type = "open")

# --- Panel (a): Correlation vs Threshold ---
p_sm_threshold <- ggplot(threshold_results,
                         aes(threshold, correlation)) +
  geom_line(linewidth = 0.8, color = "black") +
  geom_point(size = 2, color = "darkgray") +
  geom_point(data = optimal, color = "red", size = 3.5) +
  geom_hline(yintercept = original_corr,
             linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = optimal$threshold,
             linetype = "dotted", color = "red") +
  labs(x = "Soil Moisture Threshold (% nFK)",
       y = expression("Pearson correlation ("*italic("r")*")")) +
  theme_classic(base_size = 12) +
  theme(axis.line = element_line(arrow = open_arrow))


# --- Panel (b): Before–After Comparison ---
p_comparison <- ggplot() +
  geom_point(data = data_clean,
             aes(SoilMoisture, VCI, color = "Removed"),
             alpha = 0.05, size = 1.5) +
  geom_point(data = data_filtered,
             aes(SoilMoisture, VCI, color = "Retained"),
             alpha = 0.3, size = 1.5) +
  geom_smooth(data = data_clean,
              aes(SoilMoisture, VCI),
              method = "lm", se = FALSE,
              linetype = "dashed", color = "black") +
  geom_smooth(data = data_filtered,
              aes(SoilMoisture, VCI),
              method = "lm", se = TRUE,
              color = "red", fill = "pink") +
  geom_vline(xintercept = optimal$threshold,
             linetype = "dashed", color = "red") +
  scale_color_manual(values = c("Retained" = "steelblue",
                                "Removed"  = "gray60")) +
  labs(x = "Soil Moisture (% nFK)",
       y = "Vegetation Condition Index (VCI)",
       color = NULL) +
  theme_classic(base_size = 12) +
  theme(axis.line = element_line(arrow = open_arrow),
        legend.position = c(0.25, 0.12),
        legend.background = element_rect(fill = "white", color = "black"))


# --- Combined plot ---
final_plot <- p_sm_threshold + p_comparison +
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(size = 16, face = "bold"))

print(final_plot)

#ggsave(
  "Figures/Figure_9.png",
  final_plot, width = 12, height = 5.5, dpi = 600, bg = "white"
)
# ----------------------------- End of the script --------------------------------
