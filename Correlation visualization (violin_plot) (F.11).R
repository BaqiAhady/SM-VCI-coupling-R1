# ============================================================
# Title: Revealing agricultural drought sensitivity through a soil moisture
#        threshold-based framework integrating modeled and satellite data
# Purpose: Violin plots of correlation distribution under water-limited conditions
# Authors: Abdul Baqi Ahady, Stefanie Wolf, Elena-Maria Klopries
# Date: Jan 2026
# ============================================================

# -------------------- Required Packages -----------------------------
library(dplyr)
library(ggplot2)
library(sf)
library(patchwork)
library(gstat)
library(raster)

# -------------------- calling  the data ----------------------------------
# the threshold sm-vci data unter water-limited condition
opt_threshold_data <- read.csv("Data/SM_VCI_filtered_optimal.csv")

# Month names (growing season)
month_names <- c("April", "May", "June", "July", "August", "September")


# -------- Monthly calculation correlation for growing season -----------

# Here we calculate the Pearson coefficient for those points where at least 20% of the complete study period (5 observations) are available
calculate_monthly_correlations <- function(data, min_obs = 5) {
  monthly_correlations <- list()
  
  for(month_num in 4:9) {
    month_name <- month_names[month_num - 3]
    month_data <- data %>% filter(Month == month_num)
    
    # here we calculate correlation for each location
    location_corrs <- month_data %>%
      group_by(SDO_ID, Long, Lat) %>%
      summarize(
        n_observations = n(),
        correlation = ifelse(n() >= min_obs,
                             cor(SoilMoisture, VCI, use = "complete.obs"),
                             NA_real_),
        mean_sm = mean(SoilMoisture, na.rm = TRUE),
        mean_vci = mean(VCI, na.rm = TRUE),
        month = month_num,
        month_name = month_name,
        .groups = 'drop'
      ) %>%
      filter(!is.na(correlation), n_observations >= min_obs)
    
    monthly_correlations[[as.character(month_num)]] <- location_corrs
    
  }
  
  return(bind_rows(monthly_correlations))
}

# Below we check, in general, for how many location-months the r has been computed.
monthly_correlations <- calculate_monthly_correlations(opt_threshold_data, min_obs = 5)
cat("Monthly correlations calculated for", nrow(monthly_correlations), 
    "location-month combinations\n")


# -------- Correlation visualization with violin plots -------------------------
# Calculating N (number of observations) for each month (Ensuring statistical transparency)
sample_sizes <- monthly_correlations %>%
  group_by(month_name) %>%
  summarise(n = n(), mean_val = mean(correlation))
# Printing the results
cat("Sample sizes - ", 
    paste(sample_sizes$month_name, ":", sample_sizes$n, collapse = ", "),
    "\nTotal pairs:", nrow(monthly_correlations), "\n")

# Define Arrows for axes
open_arrow <- arrow(angle = 30, length = unit(0.12, "inches"), ends = "last", type = "open")

# Fixing month ordering
month_order <- c("May", "June", "July", "August", "September")
monthly_correlations$month_name <- factor(monthly_correlations$month_name, levels = month_order, ordered = TRUE)
sample_sizes$month_name <- factor(sample_sizes$month_name, levels = month_order, ordered = TRUE)

# Plot configuration
violin_plot <- ggplot(monthly_correlations, aes(x = month_name, y = correlation)) +
  geom_jitter(width = 0.18, alpha = 0.2, color = "steelblue", size = 1.2) +
  geom_violin(fill = "steelblue", alpha = 0.15, color = "steelblue",
              scale = "width", trim = TRUE, linewidth = 0.6) +
   geom_boxplot(width = 0.07, fill = "white", alpha = 0.8,
               outlier.shape = NA, color = "gray30", linewidth = 0.5) +
  
  # Connecting Mean Line between all the plots
  stat_summary(fun = mean, geom = "line", aes(group = 1),
               color = "darkred", linewidth = 0.5) +
  
  # Mean points 
  stat_summary(fun = mean, geom = "point", shape = 16, size = 3, color = "darkred") +
  
  # Adding reference line
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.7) +
  # adding lables to the right part
  geom_text(data = sample_sizes,
            aes(x = month_name, y = mean_val, label = sprintf("%.3f", mean_val)),
            nudge_x = 0.22, hjust = 0, size = 3.8, color = "darkred", fontface = "plain") +
  
  # Adding sample size at the bottom for statistical context
  geom_text(data = sample_sizes, aes(x = month_name, y = 1.15, label = paste0("n=", n)),
            size = 3.5, color = "gray30", vjust = 0) +
  
  # Labels and Scale
  labs(x = "Month", y = expression("Pearson Correlation Coefficient ("*italic("r")*")")) +
  scale_y_continuous(limits = c(-1.1, 1.2), breaks = seq(-1, 1, 0.5), expand = c(0,0)) +
  
  # Theme
  theme_classic(base_size = 14) +
  theme(
    axis.line.x = element_line(linewidth = 0.6, color = "black", arrow = open_arrow),
    axis.line.y = element_line(linewidth = 0.6, color = "black", arrow = open_arrow),
    axis.title = element_text(face = "plain"),
    axis.text = element_text(color = "black", size = 11),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  ) +
  coord_cartesian(clip = "off")

# Printing the violin griaph
print(violin_plot)

# Export the figure
#ggsave("Figures/Figure_11.png",
    violin_plot, width = 11, height = 6, dpi = 300, bg = "white")

# ----------------------------- End of the script --------------------------------




