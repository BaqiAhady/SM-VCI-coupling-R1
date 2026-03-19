# ============================================================
# Title: Revealing agricultural drought sensitivity through a soil moisture
#        threshold-based framework integrating modeled and satellite data
# Purpose: Drought characterization and identification
# Authors: Abdul Baqi Ahady, Stefanie Wolf, Elena-Maria Klopries
# Date: Jan 2026
# ============================================================

# Required packages
library(dplyr)
library(lubridate)
library(ggplot2)
library(viridis)
library(patchwork)
library(cowplot)

# -------------------- Soil Moisture (SM) -------------------------

# Access to SM data and preparing the data
sm_data <- readRDS("Data/Rur_SM_2000-2024.rds") %>%
  mutate(
    Date  = as.Date(Date),
    Year  = year(Date),
    Month = month(Date)
  )

sm_monthly <- sm_data %>%
  group_by(Year, Month) %>%
  summarise(SM = mean(SoilMoisture, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    Month = factor(sprintf("%02d", Month),
                   levels = sprintf("%02d", 1:12),
                   labels = month.abb),
    Year  = factor(Year)
  )

# Calculation of annual mean SM
sm_annual <- sm_monthly %>%
  group_by(Year) %>%
  summarise(SM = mean(SM, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    SM_z = scale(SM)[, 1],
    category = case_when(
      SM_z < -1 ~ "below",
      SM_z >  1 ~ "above",
      TRUE      ~ "normal"
    )
  )

# Monthly SM heatmap plot preparation
p_sm_heat <- ggplot(sm_monthly, aes(Month, Year, fill = SM)) +
  geom_tile(colour = "white") +
  scale_fill_viridis_c(option = "H", direction = -1, name = "SM (% nFK)") +
  scale_x_discrete(position = "top") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank())
# Monthly SM Boxplot preparation 
p_sm_box <- ggplot(sm_monthly, aes(Month, SM)) +
  geom_boxplot(fill = "lightblue", colour = "darkblue") +
  labs(x = "Month", y = "Soil Moisture (% nFK)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank())
# Annual Z-score calculation and plotting for SM
p_sm_z <- ggplot(sm_annual, aes(SM_z, Year, fill = category)) +
  geom_col(colour = "darkblue") +
  geom_vline(xintercept = c(-1, 0, 1), linetype = "dashed") +
  scale_fill_manual(values = c(below = "#FF7043",
                               normal = "#80DEEA",
                               above = "#388E3C")) +
  labs(x = "Z-score (SM)", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")
# We add legend to the heatmap
legend_sm <- get_legend(p_sm_heat)

sm_panel <- (p_sm_heat + theme(legend.position = "none")) /
  p_sm_box |
  (p_sm_z / wrap_elements(legend_sm)) +
  plot_layout(widths = c(0.8, 0.2), heights = c(1.4, 1))


# -------------------- Vegetation Condition Index (VCI) -------------------------

# Access and preparing the data
vci_data <- read.csv("Data/VCI_Rur_Monthly_Mean.csv") %>%
  mutate(
    Date  = as.Date(Date, format = "%d-%b-%y"),
    Year  = year(Date),
    Month = month(Date)
  )

vci_monthly <- vci_data %>%
  group_by(Year, Month) %>%
  summarise(VCI = mean(VCI, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    Month = factor(sprintf("%02d", Month),
                   levels = sprintf("%02d", 1:12),
                   labels = month.abb),
    Year  = factor(Year)
  )

vci_annual <- vci_data %>%
  filter(Month %in% 4:9) %>%
  group_by(Year) %>%
  summarise(VCI = mean(VCI, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    VCI_z = scale(VCI)[, 1],
    category = case_when(
      VCI_z < -1 ~ "below",
      VCI_z >  1 ~ "above",
      TRUE       ~ "normal"
    )
  )

# Monthly heatmap preparation for VCI
p_vci_heat <- ggplot(vci_monthly, aes(Month, Year, fill = VCI)) +
  geom_tile(colour = "white") +
  scale_fill_gradientn(
    colours = c("chocolate4", "gold", "yellowgreen", "darkgreen"),
    name = "VCI (%)"
  ) +
  scale_x_discrete(position = "top") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank())

# Monthly boxplot for VCI
p_vci_box <- ggplot(vci_monthly, aes(Month, VCI)) +
  geom_boxplot(fill = "lightgreen", colour = "darkgreen") +
  labs(x = "Month", y = "VCI (%)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank())

# Annual Z-score anomalies for VCI
p_vci_z <- ggplot(vci_annual, aes(VCI_z, factor(Year), fill = category)) +
  geom_col(colour = "darkblue") +
  geom_vline(xintercept = c(-1, 0, 1), linetype = "dashed") +
  scale_fill_manual(values = c(below = "#FF7043",
                               normal = "#80DEEA",
                               above = "#388E3C")) +
  labs(x = "Z-score (VCI)", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")

legend_vci <- get_legend(p_vci_heat)

vci_panel <- (p_vci_heat + theme(legend.position = "none")) /
  p_vci_box |
  (p_vci_z / wrap_elements(legend_vci)) +
  plot_layout(widths = c(0.8, 0.2), heights = c(1.4, 1))


# -------------------- Final labeled figure (EXPORT ONLY) -------------------------

combined <- wrap_plots(
  sm_panel,
  vci_panel,
  ncol = 2,
  widths = c(1, 1)
)


labeled_plot <- ggdraw(combined) +
  draw_plot_label(
    label = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"),
    x = c(0.02, 0.40, 0.02, 0.52, 0.90, 0.52),
    y = c(0.98, 0.98, 0.48, 0.98, 0.98, 0.48),
    size = 14,
    fontface = "bold"
  )

print(labeled_plot)
# ggsave("Figures/Figure_7.png",
#        labeled_plot, width = 14, height = 10, dpi = 600)
# ----------------------------- End of the script --------------------------------