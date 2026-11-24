
library(sf)
library(terra)
library(dplyr)
library(tidyr)
library(rnaturalearth)
library(marmap) # Convenient for downloading NOAA ETOPO bathymetry

# ==============================================================================
# 1. DOWNLOAD BATHYMETRY & DEFINE SHELF DOMAIN (< 200 m DEPTH)
# ==============================================================================
# Download NOAA ETOPO bathymetry for the Gulf of Mexico
bathy_raw <- getNOAA.bathy(
  lon1 = -98, lon2 = -80,
  lat1 = 24,  lat2 = 31,
  resolution = 4 # ~4 arc-minutes (~7 km); use 1 or 2 for finer resolution
)

# Convert to terra SpatRaster (note: ocean depths in marmap are negative)
bathy_rast <- marmap::as.raster(bathy_raw) |> rast()

# Keep only Gulf continental shelf waters (between 5 m and 200 m deep)
shelf_rast <- bathy_rast
shelf_rast[shelf_rast > -5 | shelf_rast < -200] <- NA

# Project raster to UTM Zone 16N in meters, then convert to 5x5 km grid
shelf_utm <- project(shelf_rast, "EPSG:32616", res = 5000) # 5000 m = 5 km

# ==============================================================================
# 2. CONVERT RASTER TO DATA FRAME & CLIP LAND POLYGONS
# ==============================================================================
base_grid <- as.data.frame(shelf_utm, xy = TRUE, na.rm = TRUE) |>
  rename(X_m = x, Y_m = y, depth_neg = layer) |>
  mutate(
    X         = X_m / 1000,         # Convert meters to km to match sdmTMB model
    Y         = Y_m / 1000,
    depth_m   = abs(depth_neg),     # Positive depth in meters
    log_depth = log(depth_m),
    area_km2  = 25                  # 5 km x 5 km = 25 km^2 per grid cell
  ) |>
  select(X, Y, depth_m, log_depth, area_km2)

# Double-check that no grid cells fall inside coastline polygons
coast <- ne_countries(scale = "medium", returnclass = "sf") |>
  st_transform(crs = 32616)
st_geometry(coast) <- st_geometry(coast) / 1000 # Convert to km

grid_sf <- st_as_sf(base_grid, coords = c("X", "Y"), crs = NA, remove = FALSE)
on_land <- lengths(st_intersects(grid_sf, coast)) > 0
base_grid <- base_grid[!on_land, ]

# ==============================================================================
# 3. EXPAND ACROSS YEARS/SEASONS & ATTACH DYNAMIC COVARIATES
# ==============================================================================
# Replicate the spatial grid for each year and season you want to predict
gulf_grid <- expand_grid(
  base_grid,
  year   = sort(unique(kmack_data$year)),
  season = factor(c("Spring", "Summer", "Fall", "Winter"),
                  levels = levels(kmack_data$season))
)

# Attach dynamic SST and Chlorophyll-a for each cell/year/season.
# Option A (Real data): Extract from satellite SpatRasters (e.g., MODIS/VIIRS/OISST)
#   using `terra::extract(sst_stack, vect(gulf_grid_ll))`
# Option B (Quick test / climatology): Populate placeholder values to test predict()
gulf_grid <- gulf_grid |>
  mutate(
    sst     = 26.5,       # Replace with extracted seasonal SST (°C)
    chl_a   = 0.45,       # Replace with extracted seasonal Chl-a (mg/m^3)
    log_chl = log(chl_a + 0.01)
  )

head(gulf_grid)