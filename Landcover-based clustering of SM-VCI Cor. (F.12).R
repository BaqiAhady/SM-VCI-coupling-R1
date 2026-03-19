# ============================================================
# Title: Revealing agricultural drought sensitivity through a soil moisture
#        threshold-based framework integrating modeled and satellite data
# Purpose: Landcover-based clustering of SM-VCI correlation  
# Authors: Abdul Baqi Ahady, Stefanie Wolf, Elena-Maria Klopries
# Date: Jan 2026
# ============================================================

# -------------------- Required Packages -----------------------------
library(dplyr)
library(ggplot2)

# -------------------- calling the data ----------------------------------
# We load the threshold sm-vci data with correlation 
monthly_correlations <- read.csv("Data/Monthly_Cor_filtered_optimal.csv")
# Landcover with points location
LC_data <- read.csv("Data/landcover_data_SDOID.csv")

# Month names (growing season)
month_names <- c("April", "May", "June", "July", "August", "September")

# -------------------- Clustering correlations by landcover type ----------------

# Here we define custom colors for each landcover type 
landcover_colors <- c(
  "Agricultural Land" = "yellow",
  "Grassland" = "lightgreen",
  "Coniferous Forest" = "darkgreen",
  "Deciduous Forest" = "green",
  "Mixed Forest" = "cyan4",
  "Other Semi-Natural Vegetation" = "yellow4",
  "Urban and Industrial" = "red",
  "Mining" = "brown",
  "Water" = "blue"
)

# Joining correlation data with landcover data (using station ID, SDO_ID)
correlations_with_landcover <- monthly_correlations %>%
  left_join(LC_data %>% dplyr::select(SDO_ID, Landcover), by = "SDO_ID") %>%
  filter(!is.na(Landcover),  # Here we remove locations without landcover data
         month_name %in% c("May", "June", "July", "August", "September"))  

# Summary statistics by landcover and month, 
# here we explore the tabulated summary of the correlations 
landcover_summary <- correlations_with_landcover %>%
  group_by(Landcover, month_name) %>%
  summarize(
    n_locations = n(),
    mean_correlation = mean(correlation, na.rm = TRUE),
    median_correlation = median(correlation, na.rm = TRUE),
    sd_correlation = sd(correlation, na.rm = TRUE),
    min_correlation = min(correlation, na.rm = TRUE),
    max_correlation = max(correlation, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(month_name, desc(mean_correlation))

print(head(landcover_summary, 20))

# Creating a filtered dataset for plotting (exclude unrelated landcovers)
correlations_for_plot <- correlations_with_landcover %>%
  filter(!Landcover %in% c("Water", "Urban and Industrial", "Mining")) %>%
  mutate(month_name = factor(month_name, 
                             levels = c("May", "June", "July", "August", "September"),
                             ordered = TRUE))

# ----------------------------- Ploting settings--------------------------------

p_landcover <- ggplot(correlations_for_plot, 
                            aes(x = Landcover, y = correlation, fill = Landcover)) +
  geom_boxplot(alpha = 0.8, outlier.size = 0.5, outlier.alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.5) +
  
  facet_grid(. ~ month_name, scales = "free_x", space = "free_x") +
  
  scale_fill_manual(values = landcover_colors) +
  scale_y_continuous(
    limits = c(-1, 1),
    breaks = c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1)
  ) +
  
  labs(
    x = NULL,
    y = "Pearson Correlation Coefficient (r)",
    fill = "Landcover Type"
  ) +
  
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_rect(
      fill = "gray90", 
      color = "black", 
      size = 0.5,
      linetype = "solid"
    ),
    strip.text = element_text(
      face = "bold", 
      size = 10, 
      color = "black",
      margin = margin(5, 0, 5, 0)
    ),
    panel.spacing = unit(0.5, "lines"),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10),
    legend.box = "horizontal"
  ) +
  guides(fill = guide_legend(nrow = 2))

print(p_landcover)

# Saving the plot
#ggsave("Figures/Figure_12.png",
       p_landcover, width = 10, height = 6, dpi = 600, bg = "white")
# ----------------------------- End of the script --------------------------------
