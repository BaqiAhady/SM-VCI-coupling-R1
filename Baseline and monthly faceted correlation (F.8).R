# ============================================================
# Title: Revealing agricultural drought sensitivity through a soil moisture
#        threshold-based framework integrating modeled and satellite data
# Purpose: Baseline and Monthly faceted SM–VCI correlation
# Authors: Abdul Baqi Ahady, Stefanie Wolf, Elena-Maria Klopries
# Date: Jan 2026
# ============================================================

# -------------------- Required packages -------------------------
library(dplyr)
library(ggplot2)
library(ggpubr)
library(scales)
library(patchwork)

# -------------------- Loading data and baseline calculation -----------------
data <- read.csv("Data/SM_with_VCI_00-24.csv")

# Growing-season filtering and NA handling

data_gs <- data %>%
  filter(Month %in% 4:9) %>%
  mutate(
    SoilMoisture = na_if(SoilMoisture, -9999),
    VCI          = na_if(VCI, -9999),
    Month_Name   = factor(month.abb[Month], levels = month.abb[4:9])
  )

# Density calculation (monthly baseline)
data_monthly <- data_gs %>%
  mutate(
    SM_Bin  = cut(SoilMoisture, breaks = 40, include.lowest = TRUE),
    VCI_Bin = cut(VCI,          breaks = 40, include.lowest = TRUE)
  ) %>%
  group_by(SM_Bin, VCI_Bin) %>%
  mutate(Density_Count = n()) %>%
  ungroup()


# Density calculation (per-month, for facets)
data_gs_dense <- data_gs %>%
  group_by(Month) %>%
  mutate(
    SM_Bin  = cut(SoilMoisture, breaks = 40, include.lowest = TRUE),
    VCI_Bin = cut(VCI,          breaks = 40, include.lowest = TRUE)
  ) %>%
  group_by(Month, SM_Bin, VCI_Bin) %>%
  mutate(Density_Count = n()) %>%
  ungroup()

# -------------------- Plotting and visualization -----------------
# Color palette and styling
smooth_density_palette <- colorRampPalette(c(
  "#D1E5F0", "royalblue3", "navyblue",
  "darkorchid4", "firebrick2", "orange", "gold"
))(100)

# Axes array configuration
open_arrow <- arrow(
  angle  = 30,
  length = unit(0.1, "inches"),
  ends   = "last",
  type   = "open"
)
# Legend dimensions
legend_height <- unit(1.5, "cm")
legend_width  <- unit(0.3, "cm")


# Plot (a): Baseline monthly correlation

p1 <- ggplot(data_monthly, aes(SoilMoisture, VCI)) +
  geom_point(aes(color = Density_Count), size = 0.7, alpha = 0.7) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "red", linewidth = 1) +
  scale_color_gradientn(
    colors = smooth_density_palette,
    limits = c(0, 2500),
    oob    = squish,
    name   = "Density",
    values = c(0, 0.01, seq(0.02, 1, length.out = 98))
  ) +
  stat_cor(
    method = "pearson",
    cor.coef.name = "r",
    aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
    label.x.npc = "left",
    label.y.npc = "bottom",
    r.accuracy  = 0.001,
    p.accuracy  = 0.001,
    size = 4.5
  ) +
  labs(x = "Soil Moisture (% nFK)", y = "VCI (%)") +
  theme_classic(base_size = 14) +
  theme(
    axis.line = element_line(linewidth = 0.6, arrow = open_arrow),
    legend.key.height = legend_height,
    legend.key.width  = legend_width
  )


# Plot (b): Monthly faceted correlations
p2 <- ggplot(data_gs_dense, aes(SoilMoisture, VCI)) +
  geom_point(aes(color = Density_Count), size = 0.5, alpha = 0.7) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "red", linewidth = 1) +
  scale_color_gradientn(
    colors = smooth_density_palette,
    limits = c(0, 400),
    oob    = squish,
    name   = "Density",
    values = c(0, 0.01, seq(0.02, 1, length.out = 98))
  ) +
  facet_wrap(~ Month_Name, ncol = 3, scales = "free") +
  stat_cor(
    method = "pearson",
    cor.coef.name = "r",
    aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
    label.x.npc = "left",
    label.y.npc = "bottom",
    r.accuracy  = 0.001,
    p.accuracy  = 0.001,
    size = 4
  ) +
  labs(x = "Soil Moisture (% nFK)", y = "VCI (%)") +
  theme_classic(base_size = 14) +
  theme(
    axis.line = element_line(linewidth = 0.5, arrow = open_arrow),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.key.height = legend_height,
    legend.key.width  = legend_width,
    panel.spacing = unit(1.5, "lines")
  )


# Final figure
final_figure <- (p1 + p2) +
  plot_layout(widths = c(1, 1.8)) +
  plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')') &
  theme(plot.tag = element_text(face = "bold", size = 18))

print(final_figure)

# ----------------------------
# Save high-resolution output
# ----------------------------
#ggsave("Figures/Figure_8.png",
  final_figure,
  width = 14, height = 7, dpi = 600, bg = "white"
)
# ----------------------------- End of the script --------------------------------





