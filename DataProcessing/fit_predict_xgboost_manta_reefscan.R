# ==============================================================================
# Presence/Absence XGBoost SDM: Manta Tow & ReefScan Survey Data
# ==============================================================================

library(terra)
library(sf)
library(dplyr)
library(readr)
library(tidyr)
library(tidymodels)
library(xgboost)

base_dir <- "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM"
dp_dir   <- file.path(base_dir, "DataProcessing")
out_dir  <- file.path(dp_dir, "outputs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

manta_csv          <- file.path(dp_dir, "data/COTS Program  Manta Tow Data-2026-02-04.csv")
reefscan_dir       <- file.path(base_dir, "data/ReefScan")
predict_clean_file <- file.path(dp_dir, "data/predictors_terra_30m_full_extended.tif")
predict_rg_file    <- file.path(dp_dir, "data/predictors_reef.tif")

output_clean_rds    <- file.path(out_dir, "fitted_model_manta_clean.rds")
output_clean_raster <- file.path(out_dir, "COTS_prob_manta_clean.tif")
output_rg_rds       <- file.path(out_dir, "fitted_model_manta_withRG.rds")
output_rg_raster    <- file.path(out_dir, "COTS_prob_manta_withRG.tif")
output_metrics_rds  <- file.path(out_dir, "manta_reefscan_performance_metrics.rds")

seed <- 123
set.seed(seed)

# --- 1. Load Predictor Stacks ---
print("Loading predictor stacks...")
pred_clean_stack <- terra::rast(predict_clean_file)
names(pred_clean_stack)[1:12] <- c("BPI", "EAST", "HCU", "NORTH", "SLO", "SVF",
                                    "VCU", "VRM", "DEM", "AhyaD", "Aspat", "Aten")
clean_cols <- names(pred_clean_stack)[1:12]

pred_rg_raw <- terra::rast(predict_rg_file)
# Selected 5 Reef Guide layers
rg_cols <- c("RG_bathy", "RG_slope", "RG_turbid", "RG_waves_Hs", "RG_waves_Tp")
pred_rg_sub <- pred_rg_raw[[rg_cols]]

# Resample/align RG layers if needed
if (!terra::compareGeom(pred_clean_stack[[1]], pred_rg_sub[[1]], stopOnError = FALSE)) {
  print("Aligning Reef Guide layers to clean raster extent...")
  pred_rg_sub <- terra::resample(pred_rg_sub, pred_clean_stack[[1]], method = "bilinear")
}

pred_rg_stack <- c(pred_clean_stack[[clean_cols]], pred_rg_sub)
full_rg_cols  <- c(clean_cols, rg_cols)

# --- 2. Ingest Manta Tow & ReefScan Survey Data ---
print("Ingesting survey data...")

# 2a. Manta Tow
manta_raw <- read_csv(manta_csv, show_col_types = FALSE) %>%
  filter(
    !is.na(StartLongitude), !is.na(EndLongitude),
    !is.na(StartLatitude),  !is.na(EndLatitude),
    StartLongitude != 0, StartLatitude != 0,
    EndLongitude != 0,   EndLatitude != 0
  ) %>%
  mutate(
    longitude = (StartLongitude + EndLongitude) / 2,
    latitude  = (StartLatitude + EndLatitude) / 2,
    is_presence = ifelse(CrownOfThornsStarfishCount > 0 | tolower(coalesce(FeedingScarCountRangeCode, "")) == "c", 1, 0),
    source = "manta_tow"
  )

# 2b. ReefScan
rs_files <- list.files(reefscan_dir, pattern = "\\.csv$", full.names = TRUE)
rs_list <- lapply(rs_files, function(f) {
  df <- read_csv(f, show_col_types = FALSE) %>%
    filter(!is.na(latitude), !is.na(longitude), latitude != 0, longitude != 0)
  
  if (nrow(df) == 0) return(NULL)
  
  cots_col <- names(df)[grep("cots", names(df), ignore.case = TRUE)][1]
  scar_col <- names(df)[grep("scar", names(df), ignore.case = TRUE)][1]
  
  cots_cnt <- if (!is.na(cots_col)) suppressWarnings(as.numeric(df[[cots_col]])) else rep(0, nrow(df))
  cots_cnt[is.na(cots_cnt)] <- 0
  
  scar_val <- if (!is.na(scar_col)) tolower(as.character(df[[scar_col]])) else rep("", nrow(df))
  
  df %>%
    mutate(
      is_presence = ifelse(cots_cnt > 0 | scar_val == "c", 1, 0),
      source = "reefscan"
    ) %>%
    select(latitude, longitude, is_presence, source)
})

reefscan_raw <- bind_rows(rs_list)

all_surveys <- bind_rows(
  manta_raw %>% select(longitude, latitude, is_presence, source),
  reefscan_raw
)

print(paste("Total raw survey observations:", nrow(all_surveys)))

# --- 3. Spatial Gridding (100m) & Effort Thresholding ---
print("Aggregating survey points onto 100m raster grid cells...")

surveys_sf <- st_as_sf(all_surveys, coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs(pred_clean_stack))

surveys_vect <- terra::vect(surveys_sf)
cell_ids     <- terra::cellFromXY(pred_clean_stack, terra::crds(surveys_vect))
surveys_sf$cell_id <- cell_ids

# Summarize presences and absences per cell
cell_summary <- surveys_sf %>%
  st_drop_geometry() %>%
  filter(!is.na(cell_id)) %>%
  group_by(cell_id) %>%
  summarise(
    n_obs       = n(),
    n_presences = sum(is_presence == 1),
    n_absences  = sum(is_presence == 0),
    .groups = "drop"
  ) %>%
  mutate(
    cots_class = case_when(
      n_presences >= 1 ~ "presence",
      n_presences == 0 & n_absences >= 3 ~ "true_absence",
      TRUE ~ "insufficient_effort"
    )
  )

# Extract spatial coordinates for each cell
cell_coords <- terra::xyFromCell(pred_clean_stack, cell_summary$cell_id)
cell_summary$x <- cell_coords[, 1]
cell_summary$y <- cell_coords[, 2]

# Print Summary Table
print("=========================================================")
print("  DATA BREAKDOWN SUMMARY: MANTA TOW & REEFSCAN GRIDDED   ")
print("=========================================================")
print(paste("Raw Manta Tow observations:   ", nrow(manta_raw)))
print(paste("Raw ReefScan observations:    ", nrow(reefscan_raw)))
print(paste("Total Unique 100m Cells:      ", nrow(cell_summary)))
print(paste("Presence Cells (>= 1 Pres):   ", sum(cell_summary$cots_class == "presence")))
print(paste("True Absence Cells (>=3 Abs): ", sum(cell_summary$cots_class == "true_absence")))
print(paste("Filtered Cells (<3 Absences): ", sum(cell_summary$cots_class == "insufficient_effort")))
print("=========================================================")

# Keep only presence and true_absence
valid_cells <- cell_summary %>%
  filter(cots_class %in% c("presence", "true_absence")) %>%
  mutate(
    cots_problem = factor(ifelse(cots_class == "presence", "1", "0"), levels = c("1", "0"))
  )

valid_sf <- st_as_sf(valid_cells, coords = c("x", "y"), crs = crs(pred_clean_stack))

# --- 4. Extract Environmental Predictors ---
print("Extracting environmental predictor values...")
ext_clean <- terra::extract(pred_clean_stack[[clean_cols]], valid_sf, ID = FALSE)
ext_rg    <- terra::extract(pred_rg_stack[[full_rg_cols]], valid_sf, ID = FALSE)

df_clean <- valid_cells %>% bind_cols(ext_clean) %>% drop_na(all_of(clean_cols))
df_rg    <- valid_cells %>% bind_cols(ext_rg) %>% drop_na(all_of(full_rg_cols))

print(paste("Clean Model dataset rows:", nrow(df_clean), "(Presences:", sum(df_clean$cots_problem == "1"), ", Absences:", sum(df_clean$cots_problem == "0"), ")"))
print(paste("Reef Guide Model dataset rows:", nrow(df_rg), "(Presences:", sum(df_rg$cots_problem == "1"), ", Absences:", sum(df_rg$cots_problem == "0"), ")"))

# --- 5. Fit & Tune XGBoost Models via Tidymodels ---
fit_xgboost_pipeline <- function(df_data, feature_cols, model_name) {
  print(paste("Fitting XGBoost model:", model_name, "..."))
  
  set.seed(seed)
  km <- kmeans(df_data[, c("x", "y")], centers = 5)
  df_data$spatial_cluster <- factor(km$cluster)
  
  folds <- group_vfold_cv(df_data, group = spatial_cluster)
  
  xgb_recipe <- recipe(cots_problem ~ ., data = df_data %>% select(cots_problem, all_of(feature_cols)))
  
  xgb_spec <- boost_tree(
    trees = tune(),
    tree_depth = tune(),
    learn_rate = tune(),
    loss_reduction = tune()
  ) %>%
    set_engine("xgboost") %>%
    set_mode("classification")
  
  xgb_wf <- workflow() %>%
    add_recipe(xgb_recipe) %>%
    add_model(xgb_spec)
  
  xgb_grid <- grid_regular(
    trees(range = c(100, 300)),
    tree_depth(range = c(3, 6)),
    learn_rate(range = c(-2, -1)),
    loss_reduction(range = c(-2, 0)),
    levels = 2
  )
  
  tune_res <- tune_grid(
    xgb_wf,
    resamples = folds,
    grid = xgb_grid,
    metrics = metric_set(roc_auc, accuracy, sens, spec)
  )
  
  best_params <- select_best(tune_res, metric = "roc_auc")
  final_wf <- finalize_workflow(xgb_wf, best_params)
  final_fit <- fit(final_wf, data = df_data)
  
  list(fit = final_fit, cv_res = tune_res, best_params = best_params)
}

res_clean <- fit_xgboost_pipeline(df_clean, clean_cols, "Clean 12-Predictor Model")
res_rg    <- fit_xgboost_pipeline(df_rg, full_rg_cols, "Reef Guide 17-Predictor Model")

# Save Models
saveRDS(res_clean, output_clean_rds)
saveRDS(res_rg, output_rg_rds)

# --- 6. Generate Spatial Predictions ---
print("Generating spatial suitability predictions across GBR domain...")

print("  Predicting Clean Model raster...")
prob_clean <- terra::predict(pred_clean_stack[[clean_cols]], res_clean$fit, type = "prob")
if (".pred_1" %in% names(prob_clean)) prob_clean <- prob_clean[[".pred_1"]]
names(prob_clean) <- "COTS_prob_manta_clean"
if (file.exists(output_clean_raster)) file.remove(output_clean_raster)
terra::writeRaster(prob_clean, output_clean_raster, overwrite = TRUE)

print("  Predicting Reef Guide Model raster...")
prob_rg <- terra::predict(pred_rg_stack[[full_rg_cols]], res_rg$fit, type = "prob")
if (".pred_1" %in% names(prob_rg)) prob_rg <- prob_rg[[".pred_1"]]
names(prob_rg) <- "COTS_prob_manta_withRG"
if (file.exists(output_rg_raster)) file.remove(output_rg_raster)
terra::writeRaster(prob_rg, output_rg_raster, overwrite = TRUE)

print("Manta Tow & ReefScan Presence/Absence XGBoost models completed successfully!")
