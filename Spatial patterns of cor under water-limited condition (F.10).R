# ============================================================
# Title: Revealing agricultural drought sensitivity through a soil moisture
#        threshold-based framework integrating modeled and satellite data
# Purpose: Spatial patterns of correlation under water-limited condition
# Authors: Abdul Baqi Ahady, Stefanie Wolf, Elena-Maria Klopries
# Date: Jan 2026
# ============================================================

# ---------------------------- Required Packages -------------------------------
library(dplyr)
library(ggplot2)
library(sf)
library(patchwork)
library(gstat)
library(raster)

# -------------------- calling  the data ---------------------------------------
# the threshold sm-vci data unter water-limited condition
opt_threshold_data <- read.csv("Data/SM_VCI_filtered_optimal.csv")

# Load shapefile for the Rur catchment
shapefile_data <- st_read("Data/Rur_Shapefile/Rur.shp") %>%
  st_transform(crs = 4326)

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

# Here we export the monthly correlation values for future use:

#write.csv(monthly_correlations, "Results/Monthly_Cor_filtered_optimal.csv", row.names = FALSE)

# ---- Spatial maps interpolation and visualization under water-limited conditions-----------

# Setting resolution 
res_val <- 0.01

# to apply the spatial patterns, we create template raster
shapefile_sp <- as_Spatial(shapefile_data)
template_raster <- raster(extent(shapefile_sp), resolution = res_val, crs = crs(shapefile_sp))
values(template_raster) <- 1

# Legend configuration
legend_plot <- ggplot(data.frame(x=1, y=1, r=1)) +
  geom_point(aes(x, y, color=r)) +
  scale_color_gradientn(
    colours = c("darkred", "red", "gray95", "blue", "darkblue"),
    values = scales::rescale(c(-1, -0.5, 0, 0.5, 1)),
    limits = c(-1, 1),
    name = "Correlation (r)",
    breaks = c(-1, -0.5, 0, 0.5, 1)
  ) +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 11, face = "plain", hjust = 0.5),
    legend.text = element_text(size = 9)
  ) +
  guides(color = guide_colorbar(
    title.position = "top", title.hjust = 0.5,
    barwidth = 0.5, barheight = 10,
    frame.colour = "black", ticks.colour = "black"
  ))

legend <- ggpubr::get_legend(legend_plot)

# Here we set the map Generation Function 
create_interpolated_map <- function(month_data, month_name, shapefile_data, limit_to_convex = FALSE, show_x = TRUE, show_y = FALSE) {
  
  pts_sf <- st_as_sf(month_data, coords = c("Long", "Lat"), crs = st_crs(shapefile_data))
  local_area_mask <- if(limit_to_convex) st_convex_hull(st_union(pts_sf)) else NULL
  
  # IDW Interpolation is functioned with with the power parameter set to 2
  points_sp <- as_Spatial(pts_sf)
  idw_model <- gstat(formula = correlation ~ 1, locations = points_sp, set = list(idp = 2))
  idw_raster <- interpolate(template_raster, idw_model)
  
  # Applying Masks to restrict the resulted interpolation into the study area
  idw_masked <- mask(idw_raster, shapefile_sp)
  if(!is.null(local_area_mask)) idw_masked <- mask(idw_masked, as_Spatial(local_area_mask))
  
  # Converting interpolated points to a raster data frame
  raster_df <- as.data.frame(idw_masked, xy = TRUE) %>%
    rename(correlation = var1.pred) %>% filter(!is.na(correlation))
  
  # Below we set the plot configuration 
  p <- ggplot() +
    geom_sf(data = shapefile_data, fill = "gray98", color = "black", linewidth = 0.2) +
    geom_raster(data = raster_df, aes(x = x, y = y, fill = correlation)) +
    geom_sf(data = shapefile_data, fill = NA, color = "black", linewidth = 0.4) +
    scale_fill_gradientn(
      colours = c("darkred", "red", "gray95", "blue", "darkblue"),
      values = scales::rescale(c(-1, -0.5, 0, 0.5, 1)),
      limits = c(-1, 1)
    ) +
    labs(title = month_name, face = "plain"
         #     subtitle = sprintf("Mean r = %.3f | n = %d",
         #                       mean(month_data$correlation, na.rm=T), nrow(month_data))
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "plain"),
          plot.subtitle = element_text(hjust = 0.5, size = 9),
          legend.position = "none",
          axis.title = element_blank(),
          # Control Axis Text visibility
          axis.text.x = if(show_x) element_text(size = 10, color = "gray40") else element_blank(),
          axis.text.y = if(show_y) element_text(size = 10, color = "gray40") else element_blank(),
          axis.ticks = element_line(linewidth = 0.2),
          panel.grid = element_line(color = "gray70", linewidth = 0.2),
    ) +
    coord_sf()
  
  return(p)
}
# Pre-filtering the data
monthly_correlations_sf <- st_as_sf(monthly_correlations, coords = c("Long", "Lat"), crs = st_crs(shapefile_data))
monthly_correlations_filtered <- st_intersection(monthly_correlations_sf, shapefile_data)
monthly_spatial_data <- monthly_correlations_filtered %>%
  mutate(Long = st_coordinates(.)[,1], Lat = st_coordinates(.)[,2]) %>% st_drop_geometry()

monthly_maps <- list()

# May: Localized data is visualized here for the month of May with only 2.5% of area coverage
may_data <- monthly_spatial_data %>% filter(month_name == "May")
monthly_maps[["May"]] <- create_interpolated_map(may_data, "May", shapefile_data,
                                                 limit_to_convex = TRUE, show_x = TRUE, show_y = TRUE)

# June-September: Full Catchment area is visualized here + Show X (E) + Hide Y (N)
for(m in c("June", "July", "August", "September")) {
  m_data <- monthly_spatial_data %>% filter(month_name == m)
  monthly_maps[[m]] <- create_interpolated_map(m_data, m, shapefile_data,
                                               limit_to_convex = FALSE, show_x = TRUE, show_y = FALSE)
}

# Assembling all the maps for the months with valid data
all_months_grid <- wrap_plots(monthly_maps[c("May", "June", "July", "August", "September")], ncol = 5)
final_spatial_map <- (all_months_grid | legend) + plot_layout(widths = c(12, 1))


# Printing the plot
print(final_spatial_map)

# Export the combined figure
#ggsave(filename = "Figures/Figure_10.png",
       plot = final_spatial_map,
       width = 13, height = 4.5, dpi = 600, bg = 'white')

# ----------------------------- End of the script --------------------------------