# ============================================================
# Title: Revealing agricultural drought sensitivity through a soil moisture
#        threshold-based framework integrating modeled and satellite data
# Purpose: Preprocessing and VCI extraction at SM locations
# Authors: Abdul Baqi Ahady, Stefanie Wolf, Elena-Maria Klopries
# Date: Jan 2026
# ============================================================

# -------------------- Configuration -------------------------
sm_dir  <- "PATH_TO_OUTPUT/New_SM_1X1"
vci_dir <- "PATH_TO_VCI_DATA/Monthly_VCI_Rur"
out_csv <- "PATH_TO_OUTPUT/SM_with_VCI.csv"

start_date <- as.Date("2000-01-01")
end_date <- as.Date("2024-12-31")
growing_months <- 4:9

# -------------------- Libraries -----------------------------
library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)
library(raster)
library(sp)

# -------------------- STEP 1: SM processing -----------------
# Load and merge SM files
sm_files <- list.files(sm_dir, pattern = "\\.csv$", full.names = TRUE)
sm_raw <- bind_rows(lapply(sm_files, read_csv))

# Cleaning and standardization
sm_clean <- sm_raw %>%
  distinct(SDO_ID, Zeitstempel, .keep_all = TRUE) %>%
  filter(Wert != -9999) %>%
  separate(SDO_GEOM, into = c("Geom_Type", "Lat", "Long"), sep = " ") %>%
  rename(
    Date = Zeitstempel,
    SoilMoisture = Wert,
    Unit = Einheit
  ) %>%
  mutate(
    Date = as.Date(Date),
    Lat  = as.numeric(gsub("\\(", "", Lat)),
    Long = as.numeric(gsub("\\)", "", Long))
  ) %>%
  filter(Date >= start_date) %>%
  select(SDO_ID, Date, SoilMoisture, Unit, Lat, Long)

# -------------------- STEP 2: VCI extraction ----------------
# Temporal attributes and growing season filter
sm_gs <- sm_clean %>%
  mutate(
    Year  = year(Date),
    Month = month(Date)
  ) %>%
  filter(Month %in% growing_months) %>%
  mutate(
    VCI_layer = paste0("VCI_", Year, "_", str_pad(Month, 2, pad = "0"))
  )

# Load VCI raster stack
vci_files <- list.files(vci_dir, pattern = "VCI_\\d{4}_\\d{2}\\.tif$", full.names = TRUE)
vci_dates <- str_extract(basename(vci_files), "\\d{4}_\\d{2}")
vci_stack <- stack(vci_files[order(vci_dates)])
names(vci_stack) <- paste0("VCI_", vci_dates)

# Spatial alignment
vci_ext <- extent(vci_stack[[1]])
sm_gs <- sm_gs %>%
  filter(
    Long >= xmin(vci_ext) & Long <= xmax(vci_ext),
    Lat  >= ymin(vci_ext) & Lat  <= ymax(vci_ext)
  )

coordinates(sm_gs) <- ~Long + Lat
proj4string(sm_gs) <- CRS("+proj=longlat +datum=WGS84")
sm_gs <- spTransform(sm_gs, crs(vci_stack))

# Layer-wise VCI extraction
vci_values <- numeric(nrow(sm_gs))
for (lyr in unique(sm_gs$VCI_layer)) {
  idx <- which(sm_gs$VCI_layer == lyr)
  if (lyr %in% names(vci_stack)) {
    vci_values[idx] <- extract(vci_stack[[lyr]], sm_gs[idx, ])
  } else {
    vci_values[idx] <- NA
  }
}

# Final dataset
sm_final <- as.data.frame(sm_gs) %>%
  mutate(VCI = vci_values) %>%
  rename(Long = coords.x1, Lat = coords.x2)

# -------------------- Output -------------------------------
# write.csv(sm_final, out_csv, row.names = FALSE)

# ----------------------------- End of the script --------------------------------