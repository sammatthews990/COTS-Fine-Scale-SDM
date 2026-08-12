############################################################
# fit_predict_cots_pseudoreps.R
# Script to fit XGBoost classification model for COTS
# using Ecological Centroids or Weighted Pseudo-Replicate Points.
# Preserves existing models and baseline files.
############################################################

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

# Setup Parameters
if (!exists("point_strategy")) point_strategy <- "ecoCentroid" # "ecoCentroid" or "pseudorepsN5"
if (!exists("use_year")) use_year <- TRUE
if (!exists("predict_year")) predict_year <- 2025
if (!exists("use_reefguide")) use_reefguide <- FALSE # FALSE = clean 12-layer stack, TRUE = full stack with RG
if (!exists("cpue_metric")) cpue_metric <- NULL

threshold <- 0.02
seed <- 123

cat(sprintf("=== Running COTS SDM: Strategy = %s | ReefGuide = %s | Year = %s ===\n",
            point_strategy, use_reefguide, ifelse(use_year, predict_year, "Agnostic")))

# Paths
predict_stack_clean <- "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/data/predictors_clean_12layer.tif"
predict_stack_full  <- "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/data/predictors_terra_30m_full_extended.tif"
cull_data_file       <- "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/data/250929_COTS-Manta-Cull-RHIS-Data-Matthews-and-Schlawinsky.xlsx"

eco_gpkg_path    <- "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/data/cull_ecological_centroids.gpkg"
pseudo_gpkg_path <- "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/data/cull_pseudo_replicates_N5.gpkg"

dir.create("C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/outputs",
           showWarnings = FALSE, recursive = TRUE)

cpue_field <- if (!is.null(cpue_metric)) cpue_metric else if (use_year) "CPUE_mean" else "CPUE_max"
print(paste("Using CPUE field:", cpue_field))

# --- 1. Load Predictors ---
print("Loading predictors...")
if (use_reefguide) {
  predictors_reef_full <- terra::rast(predict_stack_full)
  names(predictors_reef_full)[1:12] <- c("BPI", "EAST", "HCU", "NORTH", "SLO", "SVF", "VCU", "VRM", "DEM", "AhyaD", "Aspat", "Aten")
  names(predictors_reef_full)[13:17] <- c("RG_bathy", "RG_slope", "RG_turbid", "RG_waves_Hs", "RG_waves_Tp")
  clean_preds <- c("BPI", "EAST", "HCU", "NORTH", "SLO", "SVF", "VCU", "VRM", "DEM", "AhyaD", "Aspat", "Aten")
  rg_to_keep  <- c("RG_waves_Hs", "RG_waves_Tp", "RG_turbid", "RG_bathy", "RG_slope")
  pred_cols   <- c(clean_preds, rg_to_keep)
  predictors_reef <- terra::subset(predictors_reef_full, pred_cols)
} else {
  predictors_reef <- terra::rast(predict_stack_clean)
  pred_cols <- names(predictors_reef)
}

metric_suffix <- if (cpue_field == "CPUE_max" && use_year) "_maxCPUE" else ""
mode_label    <- if (use_year) paste0("year", predict_year, metric_suffix) else "agnostic"
rg_label      <- if (use_reefguide) "_withRG" else "_clean"
strat_label   <- paste0("_", point_strategy)

output_tif <- sprintf(
  "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/outputs/COTS_prob_%s_cpue_%s%s%s.tif",
  threshold, mode_label, rg_label, strat_label
)
output_vip_plot <- sprintf(
  "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/outputs/COTS_vip_%s_cpue_%s%s%s.png",
  threshold, mode_label, rg_label, strat_label
)

print(paste("Output TIF:", output_tif))
print(paste("Output VIP:", output_vip_plot))

if (use_year) {
  pred_cols <- c(pred_cols, "Year")
}

# --- 2. Load Survey Data & Join Spatial Points ---
print("Loading survey data...")
raw_cull <- read_excel(cull_data_file, sheet = 4) %>%
  mutate(
    VoyageTitle = as.factor(VoyageTitle),
    SurveyDate  = as.Date(SurveyDate),
    Year        = lubridate::year(SurveyDate),
    Total       = Cohort1 + Cohort2 + Cohort3 + Cohort4,
    CPUE        = Total / Bottomtime
  ) %>%
  filter(
    !is.na(Bottomtime), Bottomtime > 0
  )

# Aggregate survey data to CullSiteName (and Year if use_year=TRUE)
if (use_year) {
  dat_agg <- raw_cull %>%
    group_by(ReefName, CullSiteName, Year) %>%
    summarise(
      CPUE_mean  = sum(Total, na.rm = TRUE) / sum(Bottomtime, na.rm = TRUE),
      CPUE_max   = max(CPUE, na.rm = TRUE),
      Total      = sum(Total, na.rm = TRUE),
      Bottomtime = sum(Bottomtime, na.rm = TRUE),
      n_surveys  = dplyr::n(),
      .groups    = "drop"
    )
} else {
  dat_agg <- raw_cull %>%
    group_by(ReefName, CullSiteName) %>%
    summarise(
      CPUE_mean  = sum(Total, na.rm = TRUE) / sum(Bottomtime, na.rm = TRUE),
      CPUE_max   = max(CPUE, na.rm = TRUE),
      Total      = sum(Total, na.rm = TRUE),
      Bottomtime = sum(Bottomtime, na.rm = TRUE),
      n_surveys  = dplyr::n(),
      .groups    = "drop"
    )
}

cat("Loading spatial points dataset:", ifelse(point_strategy == "ecoCentroid", eco_gpkg_path, pseudo_gpkg_path), "\n")
sp_points <- st_read(ifelse(point_strategy == "ecoCentroid", eco_gpkg_path, pseudo_gpkg_path), quiet = TRUE) %>%
  select(any_of(c("CullSiteName", "ReefName", "weight", "point_type", "point_id")))

if (st_crs(sp_points) != crs(predictors_reef)) {
  sp_points <- st_transform(sp_points, crs(predictors_reef))
}

# Join aggregated survey data with spatial point locations
survey_sf <- sp_points %>%
  inner_join(dat_agg, by = c("CullSiteName", "ReefName"))

cat(sprintf("Joined survey data with spatial points: %d rows.\n", nrow(survey_sf)))

# --- 3. Spatial Predictor Extraction ---
print("Extracting raster values at survey locations...")
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

survey_df <- survey_df %>%
  mutate(
    cots_density_cpue = as.numeric(.data[[cpue_field]]),
    cots_problem = if_else(cots_density_cpue >= threshold, 1L, 0L)
  ) %>%
  filter(!is.na(cots_density_cpue), !is.na(cots_problem)) %>%
  drop_na(all_of(pred_cols))

if (!("weight" %in% names(survey_df))) {
  survey_df$weight <- 1.0
}

# --- 4. Spatial Grouping & Block CV ---
survey_df$cots_problem_fac <- factor(survey_df$cots_problem, levels = c("1", "0"))

model_df <- survey_df %>%
  transmute(
    cots_problem = cots_problem_fac,
    CullSiteName = as.factor(CullSiteName),
    ManagementZone = as.factor(ManagementZone),
    sample_weight = importance_weights(weight),
    across(all_of(pred_cols), as.numeric)
  )

print(paste("--> Final rows in model_df:", nrow(model_df)))

print("Creating grouped cross-validation folds by CullSiteName...")
set.seed(seed)
cv_folds <- group_vfold_cv(model_df, group = CullSiteName, v = 5)

# --- 5. Model Fitting (XGBoost) ---
print("Tuning and fitting XGBoost model...")

form_rec <- as.formula(paste("cots_problem ~ sample_weight +", paste(pred_cols, collapse = " + ")))

rec <- recipe(form_rec, data = model_df) %>%
  step_zv(all_predictors()) %>%
  step_normalize(all_numeric_predictors())

if (requireNamespace("doParallel", quietly = TRUE)) {
  doParallel::registerDoParallel(cores = min(4, parallel::detectCores()))
}

xgb_spec <- boost_tree(
  trees          = 1000,
  tree_depth     = tune(),
  learn_rate     = tune(),
  loss_reduction = tune(),
  min_n          = tune()
) %>%
  set_engine("xgboost", nthread = 4) %>%
  set_mode("classification")

wf <- workflow() %>%
  add_recipe(rec) %>%
  add_model(xgb_spec) %>%
  add_case_weights(sample_weight)

set.seed(seed)
grid <- grid_space_filling(
  tree_depth(range = c(3, 10)),
  learn_rate(range = c(-2.5, -0.5)),
  loss_reduction(),
  min_n(range = c(2, 20)),
  size = 15
)

res <- tune_grid(
  wf,
  resamples = cv_folds,
  grid      = grid,
  metrics   = metric_set(roc_auc, accuracy, sens, spec),
  control   = control_grid(save_pred = TRUE)
)

best <- select_best(res, metric = "roc_auc")
final_wf <- finalize_workflow(wf, best)

set.seed(seed)
fit_cls <- fit(final_wf, data = model_df)

print("Model fitting complete.")

# Collect Metrics
metrics_summary <- collect_metrics(res) %>%
  filter(.config == best$.config)
print("Best Cross-Validated Metrics:")
print(metrics_summary)

best_preds <- collect_predictions(res) %>%
  filter(.config == best$.config) %>%
  arrange(.row)

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

val_overall <- best_preds %>%
  metric_set(roc_auc, accuracy, sens, spec)(truth = cots_problem, estimate = .pred_class, .pred_1) %>%
  select(.metric, .estimate) %>%
  tidyr::pivot_wider(names_from = .metric, values_from = .estimate) %>%
  mutate(Split = "Validation (Out-of-Fold)")

val_by_fold <- calc_split_metrics(best_preds, "id")
val_by_zone <- calc_split_metrics(best_preds, "ManagementZone")

train_preds <- predict(fit_cls, new_data = model_df, type = "prob") %>%
  bind_cols(predict(fit_cls, new_data = model_df, type = "class")) %>%
  bind_cols(model_df %>% select(cots_problem))

train_overall <- train_preds %>%
  metric_set(roc_auc, accuracy, sens, spec)(truth = cots_problem, estimate = .pred_class, .pred_1) %>%
  select(.metric, .estimate) %>%
  tidyr::pivot_wider(names_from = .metric, values_from = .estimate) %>%
  mutate(Split = "Training (In-Fold)")

perf_overall <- bind_rows(train_overall, val_overall) %>%
  select(Split, accuracy, sens, spec, roc_auc)

metrics_list <- list(
  Overall = perf_overall,
  By_Zone = val_by_zone,
  By_Fold = val_by_fold
)

metrics_file <- sprintf("C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/outputs/model_performance_metrics%s%s.rds", rg_label, strat_label)
saveRDS(metrics_list, metrics_file)
print(paste("Saved model performance metrics to", metrics_file))

model_file <- sprintf("C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/outputs/fitted_model%s%s.rds", rg_label, strat_label)
saveRDS(list(fit = fit_cls, wf = final_wf, cv_res = res), model_file)
print(paste("Saved fitted model to", model_file))

# --- 6. Variable Importance Plot ---
print("Generating Variable Importance Plot...")
vip_p <- fit_cls %>%
  extract_fit_parsnip() %>%
  vip::vip(num_features = 15, geom = "col") +
  theme_minimal() +
  ggtitle(sprintf("Predictor Importance (%s - %s)", point_strategy, ifelse(use_reefguide, "With RG", "Clean")))

ggsave(output_vip_plot, plot = vip_p, width = 8, height = 6, dpi = 300)
print(paste("Saved VIP plot to:", output_vip_plot))

# --- 7. Spatial Grid Prediction ---
print(paste("Predicting spatial raster grid to:", output_tif))

pred_fun <- function(model, newdata) {
  newdata <- as.data.frame(newdata)
  if (use_year) {
    newdata$Year <- predict_year
  }
  preds <- predict(model, new_data = newdata, type = "prob")
  return(as.numeric(preds$.pred_1))
}

prob_raster <- terra::predict(
  predictors_reef,
  fit_cls,
  fun = pred_fun,
  na.rm = TRUE
)

names(prob_raster) <- "cots_prob"
terra::writeRaster(prob_raster, output_tif, overwrite = TRUE, gdal = c("COMPRESS=LZW", "TILED=YES"))
print(paste("Successfully wrote spatial prediction raster to:", output_tif))
print("Done! Model fitting and spatial prediction completed.")
