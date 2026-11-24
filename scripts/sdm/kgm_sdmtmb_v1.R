
library(sdmTMB)
library(fmesher)
library(sf)
library(rnaturalearth)
library(dplyr)
library(ggplot2)
# remotes::install_github("pbs-assess/sdmTMBextra", dependencies = TRUE)
library(sdmTMBextra)

# ==============================================================================
# 1. PREPARE DATA & PROJECT COORDINATES TO KM
# ==============================================================================
# Assume `kmack_data` has columns:
# lon, lat, year, season, catch_kg, effort_hrs, sst, chl_a, depth_m
load('c:/Users/brendan.turley/Documents/R_projects/KGM-environ-sensitivity/data/kmk_smdat.RData')
colnames(kmk_pos) <- tolower(colnames(kmk_pos))
colnames(kmk_neg) <- tolower(colnames(kmk_neg))

kmk_neg <- subset(kmk_neg, 
  depth <= max(kmk_pos$depth, na.rm = TRUE) & 
  depth >= min(kmk_pos$depth, na.rm = TRUE))
kmk_pos <- subset(kmk_pos, lon<=(-88))
kmk_neg <- subset(kmk_neg, lon<=(-88))
nrow(kmk_neg)/ nrow(kmk_pos)

kmack_data <- rbind(kmk_pos, kmk_neg) |> 
  # select(lon, lat, year, season, catch_kg, effort_hrs, sst, chl_a, depth_m)
  select(lon, lat, year, month, catch_kg = select_bgs, effort_hrs = hrs_fish, sst = tempsurf, chl_a = chlorsurf, depth_m = depth) |>
  type.convert(as.is = TRUE)

### convert months to season
kmack_data <- kmack_data |>
  mutate(
    season = case_when(
      month %in% c(12, 1, 2) ~ "winter",
      month %in% c(3, 4, 5)  ~ "spring",
      month %in% c(6, 7, 8)  ~ "summer",
      month %in% c(9, 10, 11) ~ "fall",
      TRUE ~ NA_character_
    )
  )

# Convert lat/lon to UTM Zone 16N (EPSG:32616) in kilometers
kmack_data <- add_utm_columns(
  kmack_data,
  ll_names = c("lon", "lat"),
  units = "km",
  utm_crs = 32616
)

# Transform skewed environmental covariates & scale predictors
kmack_data <- kmack_data |>
  mutate(
    log_depth = log(depth_m),
    log_chl   = log(chl_a + 0.01),
    year_f    = as.factor(year),
    season    = as.factor(season)
  )
kmack_data <- kmack_data |> na.omit()

# ==============================================================================
# 2. BUILD SPDE MESH WITH GULF OF MEXICO COASTLINE BARRIER
# ==============================================================================
# Base mesh (cutoff = 20 km minimum triangle edge; adjust finer if data allows)
base_mesh <- make_mesh(
  kmack_data,
  xy_cols = c("X", "Y"),
  cutoff = 20
)

# Download coastline polygon and project to matching UTM (km)
coast <- ne_countries(scale = "medium", returnclass = "sf") |>
  st_make_valid() #|>
  # st_transform(crs = 32616)

# Crop coastline to Gulf bounding box + buffer
gulf_bbox <- st_bbox(
  st_transform(
    st_as_sfc(st_bbox(c(xmin = -99, xmax = -80, ymin = 24, ymax = 31), crs = 4326)),
    crs = st_crs(coast)
  )
)
coast_gulf <- st_crop(coast, gulf_bbox) |>
  st_make_valid() |>
  st_transform(crs = 32616)

# Convert sf geometry from meters to kilometers to match `add_utm_columns(..., units = "km")`
st_geometry(coast_gulf) <- st_geometry(coast_gulf) / 1000

# Attach barrier structure (range_fraction = 0.1 means correlation decays 10x faster on land)
barrier_mesh <- sdmTMBextra::add_barrier_mesh(
  spde_obj = base_mesh,
  barrier_sf = coast_gulf,
  range_fraction = 0.1,
  proj_scaling = 1
)

# Visual check: blue vertices = water, green = land barrier
plot(barrier_mesh)

# ==============================================================================
# 3. FIT DELTA-GAMMA HABITAT SUITABILITY & ABUNDANCE MODEL
# ==============================================================================
# Component 1 (Binomial): Probability of occurrence (Habitat Suitability)
# Component 2 (Gamma):    Positive catch rate when present
fit_kmack <- sdmTMB(
  formula = catch_kg ~ s(sst, k = 4) + s(log_depth, k = 4) + s(log_chl, k = 4) + season,
  data = kmack_data,
  mesh = barrier_mesh,
  offset = log(kmack_data$effort_hrs),
  family = delta_gamma(link1 = "logit", link2 = "log"),
  spatial = list("off", "on"),
  spatiotemporal = list("iid", "iid"),
  # spatiotemporal = "ar1", # Captures interannual persistence in anomalies
  time = "year"
)

# Verify model convergence (all checks should return TRUE)
sanity(fit_kmack)
summary(fit_kmack)

# ==============================================================================
# 4. DIAGNOSTICS & MARGINAL HABITAT RESPONSE CURVES
# ==============================================================================
# Randomized quantile residuals
kmack_data$resid <- residuals(fit_kmack, type = "mle-mvn")
qqnorm(kmack_data$resid); abline(0, 1)

# Plot thermal preference (Model 1 = Encounter Probability / Suitability)
visreg::visreg(fit_kmack, xvar = "sst", model = 1, scale = "response",
               ylab = "Encounter Probability (Suitability)", xlab = "SST (°C)")

# ==============================================================================
# 5. PREDICT HABITAT SUITABILITY ONTO GULF SHELF GRID
# ==============================================================================
# Assume `gulf_grid` is a 5x5 km data frame of the Gulf shelf (< 200m depth)
# containing columns: X, Y, sst, log_depth, log_chl, season, year, area_km2

pred_grid <- predict(
  fit_kmack,
  newdata = gulf_grid,
  type = "response",
  return_tmb_object = TRUE # Needed if running get_index() or get_cog() next
)

# Extract predictions data frame:
# `est1` = predicted probability of presence (Habitat Suitability Index: 0 to 1)
# `est2` = predicted positive catch rate
# `est`  = combined expected CPUE (est1 * est2)
pred_df <- pred_grid$data

# Map Summer 2025 Habitat Suitability
pred_df |>
  filter(year == 2025, season == "Summer") |>
  ggplot(aes(X, Y, fill = est1)) +
  geom_raster() +
  scale_fill_viridis_c(name = "Habitat\nSuitability", limits = c(0, 1)) +
  coord_fixed() +
  labs(
    title = "King Mackerel Habitat Suitability (Summer 2025)",
    x = "Easting (km)", y = "Northing (km)"
  ) +
  theme_minimal()

# ==============================================================================
# 6. OPTIONAL: TRACK DISTRIBUTION SHIFTS (CENTER OF GRAVITY)
# ==============================================================================
cog <- get_cog(pred_grid, format = "wide")
head(cog) # Shows annual Easting/Northing centroid + SEs




#### convergence failures; procedure to step through model simplification:
# 1. Start with spatial = "off" and spatiotemporal = "iid"
# 2. If convergence issues persist, try removing spatiotemporal random effects
# 3. If still failing, simplify the family or link functions
# 4. Gradually reintroduce complexity once the model converges

summary(kmack_data[c("catch_kg", "effort_hrs", "sst", "log_depth", "log_chl")])
table(kmack_data$year)
table(kmack_data$season, useNA = "ifany")

fit_simple <- sdmTMB(
  formula = catch_kg ~ sst + log_depth + log_chl + season,
  data = kmack_data,
  mesh = barrier_mesh,
  family = delta_gamma(link1 = "logit", link2 = "log"),
  offset = log(kmack_data$effort_hrs),
  time = "year",
  ### converges
  # spatial = list("off", "off"),
  # spatiotemporal = list("off", "off")
  ### converges
  # spatial = list("off", "on"),
  # spatiotemporal = list("off", "off")
  spatial = list("off", "on"),
  # spatiotemporal = list("off", "iid") # does not converge
  spatiotemporal = list("off", "ar1"), # The model may not have converged. Maximum final gradient: 15.287682570989. Non-linear minimizer did not converge: do not trust this model!
  extra_time = c(2005, 2020, 2023) # for use with AR1 spatiotemporal model
)
### latest errors:
# ✖ `ln_tau_O` standard error may be large
# ℹ `ln_tau_O` is an internal parameter affecting `sigma_O`
# ℹ `sigma_O` is the spatial standard deviation
# ℹ Try simplifying the model, adjusting the mesh, or adding priors
# ✖ `sigma_O` is smaller than 0.01
# ℹ Consider omitting this part of the model

sanity(fit_simple)
summary(fit_simple)
kmack_data$resid <- residuals(fit_simple, type = "mle-mvn")
qqnorm(kmack_data$resid); abline(0, 1)

# Plot thermal preference (Model 1 = Encounter Probability / Suitability)
sst_grid <- data.frame(
    X = median(kmack_data$X, na.rm = TRUE),
    Y = median(kmack_data$Y, na.rm = TRUE),
  sst = seq(
    min(kmack_data$sst, na.rm = TRUE),
    max(kmack_data$sst, na.rm = TRUE),
    length.out = 100
  ),
  log_depth = median(kmack_data$log_depth, na.rm = TRUE),
  log_chl = median(kmack_data$log_chl, na.rm = TRUE),
  season = factor(
    levels(kmack_data$season)[1],
    levels = levels(kmack_data$season)
  ),
  year = kmack_data$year[1]
)

pred_sst <- predict(
  fit_simple,
  newdata = sst_grid,
  type = "response",
  model = 1,
  offset = rep(
    log(median(kmack_data$effort_hrs, na.rm = TRUE)),
    nrow(sst_grid)
  )
)

ggplot(pred_sst, aes(x = sst, y = est)) +
  geom_line() +
  labs(
    x = "SST (°C)",
    y = "Encounter probability"
  ) +
  theme_minimal()

# If that converges, add complexity in this order:
# spatial = list("off", "on")
# spatiotemporal = list("off", "iid")
# then add the smooth terms one at a time. This is especially relevant if you want
# to omit `sigma_O`: the first `"off"` removes the occurrence component’s spatial field.


### or simple prescence/absence model

# Convert catch to binary presence (1) / absence (0)
kmack_data$present <- ifelse(kmack_data$catch_kg > 0, 1, 0)

# Fit a fast, stable presence/absence habitat suitability model
# (Note on offset in a binomial model: Using offset = log(effort_hrs) works well when encounter probabilities are low,
#  like your ~5% rate, because under a Poisson process $P(\text{catch} > 0) = 1 - e^{-\lambda \cdot \text{effort}} \approx \lambda \cdot \text{effort}$ when $\lambda$ is small.
#  Alternatively, using cloglog link—binomial(link = "cloglog")—makes offset = log(effort_hrs) mathematically exact for presence/absence data!)

fit_hsi <- sdmTMB(
  present ~ s(sst, k = 4) + s(log_depth, k = 4) + s(log_chl, k = 4) + season,
  data = kmack_data,
  mesh = barrier_mesh,
  family = binomial(link = "logit"),
  offset = log(kmack_data$effort_hrs), # Accounts for varying tow/soak times
  time = "year",
  spatial = "on",
  spatiotemporal = "ar1"
)

sanity(fit_hsi)
summary(fit_hsi)

kmack_data$resid_hsi <- residuals(fit_hsi, type = "mle-mvn")
qqnorm(kmack_data$resid_hsi); abline(0, 1)

### convergence issues
# 1. Change spatiotemporal = "ar1" to "iid" (estimating the AR1 correlation parameter $\rho$ is often noisy 
# when presences are patchy).
# 2. Coarsen the mesh slightly (e.g., increase cutoff = 20 to cutoff = 30 km so fewer mesh vertices fall in 
# data-sparse areas).
# 3. Reduce spline wiggliness (k = 3 instead of k = 4 or 5).
# 4. Replace fixed year with random (1 | year_f) if you included year in the fixed formula and any year had 
# very few positives.