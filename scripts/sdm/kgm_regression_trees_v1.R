### untested Gulf of Mexico spatial XGBoost model for KGM data

library(sf)
library(dplyr)
library(blockCV)
library(xgboost)
library(pROC)
library(ggplot2)

# ==============================================================================
# 1. PREPARE SPATIAL DATA & FEATURES
# ==============================================================================
# Assume `kmack_data` has: lon, lat, catch_kg, effort_hrs, sst, depth_m, chl_a, season

load('c:/Users/brendan.turley/Documents/R_projects/KGM-environ-sensitivity/data/kmk_smdat.RData')
colnames(kmk_pos) <- tolower(colnames(kmk_pos))
colnames(kmk_neg) <- tolower(colnames(kmk_neg))

kmk_neg <- subset(kmk_neg, 
  depth <= max(kmk_pos$depth, na.rm = TRUE) & 
  depth >= min(kmk_pos$depth, na.rm = TRUE))
kmk_pos <- subset(kmk_pos, lon<=(-88))
kmk_neg <- subset(kmk_neg, lon<=(-88))

kmack_data <- rbind(kmk_pos, kmk_neg) |> 
  # select(lon, lat, year, season, catch_kg, effort_hrs, sst, chl_a, depth_m)
  select(lon, lat, year, month, catch_kg = select_bgs, effort_hrs = hrs_fish,
    sst = tempsurf, chl_a = chlorsurf, depth_m = depth) |>
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

kmack_df <- kmack_data |>
  mutate(
    present   = ifelse(catch_kg > 0, 1, 0),
    log_depth = log(depth_m),
    log_chl   = log(chl_a + 0.01),
    # Convert season to numeric/dummy encoding for xgboost
    season_n  = as.numeric(factor(season))
  )

# Convert to sf object in projected UTM coordinates (meters) for accurate block sizing
kmack_sf <- st_as_sf(kmack_df, coords = c("lon", "lat"), crs = 4326) |>
  st_transform(crs = 32616) # UTM Zone 16N (Gulf of Mexico)

# ==============================================================================
# 2. CREATE SPATIAL BLOCKS WITH blockCV
# ==============================================================================
# Choose block size larger than the spatial autocorrelation range of your data
# (e.g., 150,000 m = 150 km blocks across the Gulf shelf)
set.seed(42)
sb <- cv_spatial(
  x = kmack_sf,
  column = "present",      # Stratifies folds so each fold gets ~equal presences
  k = 5,                   # 5-fold spatial cross-validation
  size = 150000,           # 150 km x 150 km spatial blocks
  selection = "random",    # Randomly assign blocks to the 5 folds
  iteration = 50,          # Finds balanced presence/absence distribution across folds
  progress = FALSE
)

# Inspect fold map to confirm spatial separation across the Gulf
cv_plot(cv = sb, x = kmack_sf)

# ==============================================================================
# 3. RUN SPATIAL CROSS-VALIDATION WITH XGBOOST
# ==============================================================================
# Select environmental predictors (DO NOT include lat/lon or X/Y!)
predictors <- c("sst", "log_depth", "log_chl", "season_n")
X_mat <- as.matrix(kmack_df[, predictors])
y_vec <- kmack_df$present

# Account for variable survey effort via log(effort_hrs) base margin (offset)
effort_offset <- log(kmack_df$effort_hrs)

# Store out-of-fold predictions and AUC scores
oof_preds <- rep(NA, nrow(kmack_df))
fold_aucs <- numeric(5)

# Conservative tree hyperparameters to prevent overfitting on 5% prevalence data
xgb_params <- list(
  objective        = "binary:logistic",
  eval_metric      = "auc",
  max_depth        = 3,       # Shallow trees capture 2-3 way interactions safely
  eta              = 0.03,    # Slow learning rate
  subsample        = 0.75,    # Row subsampling per tree
  colsample_bytree = 0.80,    # Feature subsampling per tree
  min_child_weight = 10       # Prevents splits on tiny clusters of rare positives
)

best_ntrees <- numeric(5)

for (i in seq_len(5)) {
  train_idx <- sb$folds_list[[i]][[1]]
  test_idx  <- sb$folds_list[[i]][[2]]

  # Build xgb.DMatrix with effort offset attached
  dtrain <- xgb.DMatrix(data = X_mat[train_idx, ], label = y_vec[train_idx])
  setinfo(dtrain, "base_margin", effort_offset[train_idx])

  dtest <- xgb.DMatrix(data = X_mat[test_idx, ], label = y_vec[test_idx])
  setinfo(dtest, "base_margin", effort_offset[test_idx])

  # Train with early stopping on the spatially held-out fold
  bst <- xgb.train(
    params = xgb_params,
    data = dtrain,
    nrounds = 1500,
    evals = list(test = dtest),
    early_stopping_rounds = 50,
    verbose = 0
  )

  best_iter <- xgb.attributes(bst)$best_iteration
  if (length(best_iter) == 0L || is.na(best_iter)) {
    best_iter <- nrow(bst$evaluation_log) - 1L
  }
  # XGBoost iteration numbers are zero-based
  best_ntrees[i] <- as.integer(best_iter) + 1L

  # Predict standardized to 1 unit of effort (base_margin = log(1) = 0)
  dtest_std <- xgb.DMatrix(data = X_mat[test_idx, ])
  setinfo(dtest_std, "base_margin", rep(0, length(test_idx)))

  oof_preds[test_idx] <- predict(bst, dtest_std)
  fold_aucs[i] <- as.numeric(auc(y_vec[test_idx], oof_preds[test_idx], quiet = TRUE))
}

# Report spatially cross-validated performance
cat("Spatial Block Fold AUCs:", round(fold_aucs, 3), "\n")
cat("Mean Spatial CV AUC:    ", round(mean(fold_aucs), 3), "\n")

# ==============================================================================
# 4. FIT FINAL MODEL & PREDICT HABITAT SUITABILITY ONTO GULF GRID
# ==============================================================================
# Train final model on all data using median optimal trees from spatial CV
dfull <- xgb.DMatrix(data = X_mat, label = y_vec)
setinfo(dfull, "base_margin", effort_offset)

final_xgb <- xgb.train(
  params  = xgb_params,
  data    = dfull,
  nrounds = round(median(best_ntrees)),
  verbose = 0
)

# Variable importance (Gain = relative contribution to reducing loss)
imp <- xgb.importance(feature_names = predictors, model = final_xgb)
print(imp)

season_levels <- levels(factor(kmack_df$season))
gulf_grid <- expand.grid(
  lon = seq(-98, -80, length.out = 270),
  lat = seq(24, 31, length.out = 105)
) |>
  mutate(
    sst = median(kmack_df$sst, na.rm = TRUE),
    log_depth = median(kmack_df$log_depth, na.rm = TRUE),
    log_chl = median(kmack_df$log_chl, na.rm = TRUE),
    season = factor("summer", levels = season_levels)
  )

# Predict onto `gulf_grid` (standardized to 1 hour of effort: base_margin = 0)
gulf_grid_mat <- gulf_grid |>
  mutate(season_n = as.numeric(factor(season, levels = levels(factor(kmack_df$season))))) |>
  select(all_of(predictors)) |>
  as.matrix()

dgrid <- xgb.DMatrix(data = gulf_grid_mat)
setinfo(dgrid, "base_margin", rep(0, nrow(gulf_grid)))

gulf_grid$hsi_xgb <- predict(final_xgb, dgrid)


library(ggplot2)
library(dplyr)
library(tidyr)
library(sf)
library(rnaturalearth)
library(pROC)
library(patchwork) # For combining plots cleanly

# ==============================================================================
# 1. GULF OF MEXICO HABITAT SUITABILITY MAP (WITH COASTLINE)
# ==============================================================================
# Prepare coastline in matching UTM Zone 16N (km)
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
  st_make_valid() #|>
  # st_transform(crs = 32616)
# st_geometry(coast_gulf) <- st_geometry(coast_gulf) / 1000

# Filter prediction grid to a target year/season (or average across years by season)
map_df <- gulf_grid |>
  group_by(lon, lat, season) |>
  summarise(mean_hsi = mean(hsi_xgb, na.rm = TRUE), .groups = "drop")

p_map <- ggplot() +
  geom_raster(data = map_df, aes(x = lon, y = lat, fill = mean_hsi)) +
  geom_sf(data =  coast_gulf, fill = "grey85", color = "grey40", linewidth = 0.3) +
  facet_wrap(~season, ncol = 2) +
  scale_fill_viridis_c(
    option = "turbo",
    name = "Encounter\nProbability",
    limits = c(0, max(map_df$mean_hsi, na.rm = TRUE))
  ) +
  coord_sf(
    xlim = range(map_df$lon),
    ylim = range(map_df$lat),
    expand = FALSE
  ) +
  labs(
    title = "Seasonal King Mackerel Habitat Suitability (XGBoost)",
    subtitle = "Standardized to 1 unit of survey effort across the Gulf shelf (< 200 m)",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank())

print(p_map)

# ==============================================================================
# 2. VARIABLE IMPORTANCE PLOT (GAIN)
# ==============================================================================
imp_df <- xgb.importance(feature_names = predictors, model = final_xgb) |>
  mutate(
    Feature = recode(
      Feature,
      "sst"       = "Sea Surface Temp (°C)",
      "log_depth" = "Log Depth (m)",
      "log_chl"   = "Log Chlorophyll-a",
      "season_n"  = "Season"
    )
  )

p_imp <- ggplot(imp_df, aes(x = reorder(Feature, Gain), y = Gain)) +
  geom_col(fill = "#1f78b4", width = 0.65) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Environmental Driver Importance",
    subtitle = "Relative contribution (Gain) to splits across all boosted trees",
    x = NULL,
    y = "Relative Importance (Gain)"
  ) +
  theme_minimal(base_size = 12)

print(p_imp)

# ==============================================================================
# 3. PARTIAL DEPENDENCE CURVES (ECOLOGICAL RESPONSE SHAPES)
# ==============================================================================
# Helper function to compute partial dependence on probability scale (1 hr effort)
calc_pdp <- function(var_name, data_mat, model, grid_len = 50) {
  vals <- seq(
    quantile(data_mat[, var_name], 0.02, na.rm = TRUE),
    quantile(data_mat[, var_name], 0.98, na.rm = TRUE),
    length.out = grid_len
  )
  
  # Subsample background rows for speed if dataset is very large
  bg <- data_mat[sample(nrow(data_mat), min(1000, nrow(data_mat))), , drop = FALSE]
  
  preds <- sapply(vals, function(v) {
    tmp <- bg
    tmp[, var_name] <- v
    dmat <- xgb.DMatrix(data = tmp)
    setinfo(dmat, "base_margin", rep(0, nrow(tmp))) # Standardized effort
    mean(predict(model, dmat))
  })
  
  data.frame(variable = var_name, x_val = vals, yhat = preds)
}

pdp_df <- bind_rows(
  calc_pdp("sst", X_mat, final_xgb),
  calc_pdp("log_depth", X_mat, final_xgb),
  calc_pdp("log_chl", X_mat, final_xgb)
) |>
  # Back-transform log variables for intuitive axis labels
  mutate(
    x_plot = case_when(
      variable == "log_depth" ~ exp(x_val),
      variable == "log_chl"   ~ pmax(exp(x_val) - 0.01, 0),
      TRUE                    ~ x_val
    ),
    label = recode(
      variable,
      "sst"       = "Sea Surface Temperature (°C)",
      "log_depth" = "Bottom Depth (m)",
      "log_chl"   = "Chlorophyll-a (mg/m³)"
    )
  )

p_pdp <- ggplot(pdp_df, aes(x = x_plot, y = yhat)) +
  geom_line(color = "#d95f02", linewidth = 1.2) +
  facet_wrap(~label, scales = "free_x", ncol = 3) +
  labs(
    title = "Partial Dependence Habitat Response Curves",
    subtitle = "Marginal effect of each variable holding others at observed distributions",
    x = NULL,
    y = "Predicted Encounter Probability"
  ) +
  theme_minimal(base_size = 12)

print(p_pdp)

# ==============================================================================
# 4. SPATIAL CV ROC CURVE & PROBABILITY CALIBRATION PLOT
# ==============================================================================
# A. ROC Curve from spatially held-out predictions (`oof_preds`)
roc_obj <- roc(y_vec, oof_preds, quiet = TRUE)
auc_val <- round(as.numeric(auc(roc_obj)), 3)

roc_df <- data.frame(
  fpr = 1 - roc_obj$specificities,
  tpr = roc_obj$sensitivities
)

p_roc <- ggplot(roc_df, aes(x = fpr, y = tpr)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(color = "#1b9e77", linewidth = 1.2) +
  annotate("text", x = 0.65, y = 0.25, label = paste("Spatial CV AUC =", auc_val),
           size = 4.5, fontface = "bold") +
  labs(
    title = "Spatial Block CV ROC Curve",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  theme_minimal(base_size = 12)

# B. Calibration Plot (Crucial for 5% prevalence rare-event data!)
# Checks if predicted 5% probability actually corresponds to 5% observed catches
calib_df <- data.frame(obs = y_vec, pred = oof_preds) |>
  mutate(bin = ntile(pred, 10)) |>
  group_by(bin) |>
  summarise(
    mean_pred = mean(pred, na.rm = TRUE),
    obs_rate  = mean(obs, na.rm = TRUE),
    se        = sqrt(obs_rate * (1 - obs_rate) / n()),
    .groups   = "drop"
  )

p_cal <- ggplot(calib_df, aes(x = mean_pred, y = obs_rate)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(ymin = pmax(0, obs_rate - 1.96 * se),
                    ymax = pmin(1, obs_rate + 1.96 * se)),
                width = 0.003, color = "#7570b3") +
  geom_point(size = 3, color = "#7570b3") +
  labs(
    title = "Probability Calibration (Deciles)",
    x = "Mean Predicted Suitability",
    y = "Observed Occurrence Rate"
  ) +
  theme_minimal(base_size = 12)

# Combine ROC and Calibration side-by-side
print(p_roc + p_cal)
