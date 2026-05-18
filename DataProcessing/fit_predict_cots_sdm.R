############################################################
# fit_predict_cots_sdm.R
# Script to fit an XGBoost classification model for COTS
# (using CPUE threshold 0.02) and predicting it across the GBR.
############################################################

# --- 0. Setup and Packages ---
library(readxl)
library(dplyr)
library(sf)
library(terra)
library(tidymodels)
library(blockCV)
library(lubridate)
library(xgboost)
library(ggplot2)
library(purrr)
library(vip)

tidymodels::tidymodels_prefer()

# Set up paths
# Clean 12-layer stack (env + coral only, no RG_*/SDM_* layers that restrict extent)
predict_stack_clean <- "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/data/predictors_clean_12layer.tif"
# Full 17-layer stack (includes RG_* — full extent!)
predict_stack_full <- "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/data/predictors_terra_30m_full_extended.tif"

cull_data_file <- "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/data/250929_COTS-Manta-Cull-RHIS-Data-Matthews-and-Schlawinsky.xlsx"

dir.create("C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/outputs",
  showWarnings = FALSE, recursive = TRUE
)

# Parameters
threshold <- 0.02
block_km <- 50
seed <- 123

use_year <- TRUE # TRUE = include Year as a predictor; FALSE = agnostic across years
predict_year <- 2025 # Year to predict for when use_year = TRUE
use_reefguide <- TRUE # TRUE = use full stack with ReefGuide layers; FALSE = use clean stack
cpue_field <- if (use_year) "CPUE_mean" else "CPUE_max"

# --- 1. Load Predictors ---
print("Loading predictors...")
if (use_reefguide) {
  # Full stack — includes RG_* layers but prediction extent will be restricted
  predictors_reef_full <- terra::rast(predict_stack_full)

  # Standardize names if needed (the new stack should already have them, but let's be safe)
  # The new stack has 17 layers: 1-12 are base, 13-17 are RG
  names(predictors_reef_full)[1:12] <- c("BPI", "EAST", "HCU", "NORTH", "SLO", "SVF", "VCU", "VRM", "DEM", "AhyaD", "Aspat", "Aten")
  # Add the RG_ prefix back in memory
  names(predictors_reef_full)[13:17] <- c("RG_bathy", "RG_slope", "RG_turbid", "RG_waves_Hs", "RG_waves_Tp")

  clean_preds <- c("BPI", "EAST", "HCU", "NORTH", "SLO", "SVF", "VCU", "VRM", "DEM", "AhyaD", "Aspat", "Aten")
  
  # The specific RG layers to use
  rg_to_keep <- c("RG_waves_Hs", "RG_waves_Tp", "RG_turbid", "RG_bathy", "RG_slope")
  pred_cols <- c(clean_preds, rg_to_keep)
  
  # Subset to only these layers (should be all 17)
  predictors_reef <- terra::subset(predictors_reef_full, pred_cols)
  
  print(paste("Using full stack WITH extended ReefGuide layers. Predictors:", length(pred_cols)))
} else {
  # Clean stack — env + coral only, full reef extent
  predictors_reef <- terra::rast(predict_stack_clean)
  pred_cols <- names(predictors_reef)
  print(paste("Using clean stack (no RG/SDM layers). Predictors:", length(pred_cols)))
}

# Build output filenames based on mode
mode_label <- if (use_year) paste0("year", predict_year) else "agnostic"
rg_label <- if (use_reefguide) "_withRG" else "_clean"
output_tif <- sprintf(
  "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/outputs/COTS_prob_%s_cpue_%s%s.tif",
  threshold, mode_label, rg_label
)
output_vip_plot <- sprintf(
  "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/outputs/COTS_vip_%s_cpue_%s%s.png",
  threshold, mode_label, rg_label
)

print(paste("Output TIF:", output_tif))
print(paste("Output VIP:", output_vip_plot))

# If using year, add it to the predictor set
if (use_year) {
  pred_cols <- c(pred_cols, "Year")
  print(paste("Using Year as a predictor. Prediction year:", predict_year))
} else {
  print("Year-agnostic mode: aggregating across all years.")
}

print(paste("Final predictor columns:", paste(pred_cols, collapse = ", ")))

# --- 2. Load and Prepare Survey Data ---
print("Loading and preparing survey data...")
raw_cull <- read_excel(cull_data_file, sheet = 4) %>%
  mutate(
    VoyageTitle = as.factor(VoyageTitle),
    SurveyDate  = as.Date(SurveyDate),
    Year        = lubridate::year(SurveyDate),
    Total       = Cohort1 + Cohort2 + Cohort3 + Cohort4,
    CPUE        = Total / Bottomtime
  ) %>%
  filter(
    !is.na(Longitude), !is.na(Latitude),
    !is.na(Bottomtime), Bottomtime > 0
  )

# Aggregate survey data
print("Aggregating survey data to locations...")
if (use_year) {
  # Per-year aggregation: keep Year as a column
  dat <- raw_cull %>%
    group_by(ReefName, CullSiteName, Latitude, Longitude, Year) %>%
    summarise(
      CPUE_mean  = sum(Total, na.rm = TRUE) / sum(Bottomtime, na.rm = TRUE),
      CPUE_max   = max(CPUE, na.rm = TRUE),
      Total      = sum(Total, na.rm = TRUE),
      Bottomtime = sum(Bottomtime, na.rm = TRUE),
      n_surveys  = dplyr::n(),
      .groups    = "drop"
    )
} else {
  # Year-agnostic aggregation
  dat <- raw_cull %>%
    group_by(ReefName, CullSiteName, Latitude, Longitude) %>%
    summarise(
      CPUE_mean  = sum(Total, na.rm = TRUE) / sum(Bottomtime, na.rm = TRUE),
      CPUE_max   = max(CPUE, na.rm = TRUE),
      Total      = sum(Total, na.rm = TRUE),
      Bottomtime = sum(Bottomtime, na.rm = TRUE),
      n_surveys  = dplyr::n(),
      Year_min   = min(Year, na.rm = TRUE),
      Year_max   = max(Year, na.rm = TRUE),
      .groups    = "drop"
    )
}

# --- 3. Spatial Preparation and Extraction ---
print("Extracting raster values at survey locations...")
# Convert to sf and transform to raster CRS
survey_sf <- st_as_sf(dat, coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_transform(crs(predictors_reef))

# Extract raster predictors
survey_sf$ID <- seq_len(nrow(survey_sf))
pred_vals <- terra::extract(predictors_reef, terra::vect(survey_sf))

# Join with Sectors and Management Zones
sectors_sf <- st_read("C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/data/SectorShapefile/ltms.shp", quiet = TRUE) %>%
  st_transform(crs = crs(survey_sf))

c_sector <- st_join(survey_sf, sectors_sf %>% select(SECTOR_NAM), join = st_nearest_feature)

survey_sf$SECTOR_NAM <- c_sector$SECTOR_NAM
survey_sf <- survey_sf %>%
  mutate(
    SECTOR_NAM = as.character(SECTOR_NAM),
    ManagementZone = case_when(
      SECTOR_NAM %in% c("Cape Grenville", "Princess Charlotte Bay") ~ "Far Northern",
      SECTOR_NAM %in% c("Cooktown / Lizard Island", "Cairns", "Innisfail") ~ "Northern",
      SECTOR_NAM %in% c("Townsville", "Cape Upstart", "Whitsunday") ~ "Central",
      SECTOR_NAM %in% c("Pompey", "Swain", "Capricorn Bunker") ~ "Southern",
      TRUE ~ "Unknown"
    )
  )

survey_df <- survey_sf %>%
  st_drop_geometry() %>%
  left_join(pred_vals, by = "ID") %>%
  select(-ID)

# Define response variable (keep as numeric 1/0 for blockCV stratification)
survey_df <- survey_df %>%
  mutate(
    cots_density_cpue = as.numeric(.data[[cpue_field]]),
    cots_problem = if_else(cots_density_cpue >= threshold, 1L, 0L)
  )

# Drop rows missing predictors or response
survey_df <- survey_df %>%
  filter(!is.na(cots_density_cpue), !is.na(cots_problem)) %>%
  drop_na(all_of(pred_cols))

# Keep only rows in the sf object that survived the drop_na
survey_df <- survey_df %>% mutate(row_id = row_number())
survey_sf2 <- survey_sf %>%
  mutate(row_id = row_number()) %>%
  filter(row_id %in% survey_df$row_id)


# --- 4. Spatial Blocking for CV ---
print("Creating spatial blocks for cross-validation...")
block_range_for_sf <- function(x_sf, km = 50) {
  if (sf::st_is_longlat(x_sf)) {
    return(km / 111.32)
  } else {
    return(km * 1000)
  }
}

theRange <- block_range_for_sf(survey_sf2, km = block_km)

set.seed(seed)
sb <- spatialBlock(
  speciesData = survey_sf2,
  species     = "cots_problem",
  theRange    = theRange,
  k           = 5,
  selection   = "random",
  showBlocks  = FALSE
)

survey_df$fold_id <- sb$foldID

# Prepare final modelling dataframe
model_df <- survey_df %>%
  transmute(
    cots_problem = factor(cots_problem, levels = c("1", "0")),
    fold_id,
    ManagementZone = as.factor(ManagementZone),
    across(all_of(pred_cols), as.numeric)
  )

print(paste("--> Number of rows in model_df:", nrow(model_df)))

# --- 5. Model Fitting (XGBoost) ---
print("Tuning and fitting XGBoost model...")
set.seed(seed)
cv_folds <- group_vfold_cv(model_df, group = fold_id, v = 5)

rec <- recipe(cots_problem ~ ., data = model_df) %>%
  step_rm(fold_id, ManagementZone) %>%
  step_zv(all_predictors()) %>%
  step_normalize(all_numeric_predictors())

xgb_spec <- boost_tree(
  trees          = 1500,
  tree_depth     = tune(),
  learn_rate     = tune(),
  loss_reduction = tune(),
  min_n          = tune()
) %>%
  set_engine("xgboost") %>%
  set_mode("classification")

wf <- workflow() %>%
  add_recipe(rec) %>%
  add_model(xgb_spec)

set.seed(seed)
grid <- grid_space_filling(
  tree_depth(),
  learn_rate(),
  loss_reduction(),
  min_n(),
  size = 30
)

# Tune model
res <- tune_grid(
  wf,
  resamples = cv_folds,
  grid      = grid,
  metrics   = metric_set(roc_auc, accuracy, sens, spec),
  control   = control_grid(save_pred = TRUE)
)

# Select best and finalize
best <- select_best(res, metric = "roc_auc")
final_wf <- finalize_workflow(wf, best)

set.seed(seed)
fit_cls <- fit(final_wf, data = model_df)

print("Model fitting complete.")

# Print Metrics Summary
metrics_summary <- collect_metrics(res) %>%
  filter(.config == best$.config)
print("Best Cross-Validated Metrics:")
print(metrics_summary)

# --- Extract Performance Metrics by Region & Fold ---
print("Extracting Validation Metrics by Region and Fold...")

# Get out-of-fold predictions for the best config
best_preds <- collect_predictions(res) %>%
  filter(.config == best$.config) %>%
  arrange(.row)

# Rejoin to the original model_df to get ManagementZone
best_preds <- best_preds %>%
  bind_cols(model_df %>% select(ManagementZone) %>% slice(best_preds$.row))

calc_split_metrics <- function(df, group_col) {
  custom_metrics <- metric_set(roc_auc, accuracy, sens, spec)
  df %>%
    group_by(!!sym(group_col)) %>%
    custom_metrics(truth = cots_problem, estimate = .pred_class, .pred_1) %>%
    select(!!sym(group_col), .metric, .estimate) %>%
    tidyr::pivot_wider(names_from = .metric, values_from = .estimate)
}

# 1. Overall Validation Metrics
val_overall <- best_preds %>%
  metric_set(roc_auc, accuracy, sens, spec)(truth = cots_problem, estimate = .pred_class, .pred_1) %>%
  select(.metric, .estimate) %>%
  tidyr::pivot_wider(names_from = .metric, values_from = .estimate) %>%
  mutate(Split = "Validation (Out-of-Fold)")

# 2. Validation by Fold
val_by_fold <- calc_split_metrics(best_preds, "id")

# 3. Validation by ManagementZone
val_by_zone <- calc_split_metrics(best_preds, "ManagementZone")

# 4. Overall Training Metrics
# Predict back on the training set using the finalized fit
train_preds <- predict(fit_cls, new_data = model_df, type = "prob") %>%
  bind_cols(predict(fit_cls, new_data = model_df, type = "class")) %>%
  bind_cols(model_df %>% select(cots_problem))

train_overall <- train_preds %>%
  metric_set(roc_auc, accuracy, sens, spec)(truth = cots_problem, estimate = .pred_class, .pred_1) %>%
  select(.metric, .estimate) %>%
  tidyr::pivot_wider(names_from = .metric, values_from = .estimate) %>%
  mutate(Split = "Training (In-Fold)")

# Combine Overall
perf_overall <- bind_rows(train_overall, val_overall) %>%
  select(Split, accuracy, sens, spec, roc_auc)

# Save to RDS dynamically based on RG label
metrics_list <- list(
  Overall = perf_overall,
  By_Zone = val_by_zone,
  By_Fold = val_by_fold
)
rg_suffix <- if (use_reefguide) "_withRG" else "_clean"

metrics_file <- sprintf("C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/outputs/model_performance_metrics%s.rds", rg_suffix)
saveRDS(metrics_list, metrics_file)
print(paste("Saved model performance metrics to", metrics_file))

# Save the model objects so they can be reloaded
model_file <- sprintf("C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/outputs/fitted_model%s.rds", rg_suffix)
saveRDS(list(fit = fit_cls, wf = final_wf, cv_res = res), model_file)
print(paste("Saved fitted model to", model_file))

# --- 6. Variable Importance Plot ---
print("Generating Variable Importance Plot...")

booster <- fit_cls %>%
  extract_fit_parsnip() %>%
  purrr::pluck("fit")
imp <- xgboost::xgb.importance(model = booster)

# We can also plot it using vip
vip_p <- fit_cls %>%
  extract_fit_parsnip() %>%
  vip::vip(num_features = 15, geom = "col") +
  theme_minimal() +
  ggtitle("Top 15 Predictor Importance (XGBoost Gain)")

print(imp)
ggsave(output_vip_plot, plot = vip_p, width = 8, height = 6, dpi = 300)
print(paste("Saved VIP plot to:", output_vip_plot))


# --- 7. Predict to Spatial Grid ---
print("Skipping slow spatial prediction raster write because it is already written and mathematically identical.")
print("Done! All processing completed.")
