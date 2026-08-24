############################################################
# validate_model_comprehensive.R
# Multi-level validation framework for COTS SDMs.
# Validates any fitted model RDS against:
#   Level 1: Internal CV metrics (from saved model)
#   Level 2: Out-of-time (train ≤2025, test 2026)
#   Level 3: Leave-one-sector-out spatial transferability
#   Level 4: Operational metrics (ghost culls, financial savings)
#   Level 5: Calibration & partial dependence plots
############################################################

library(readxl)
library(dplyr)
library(sf)
library(terra)
library(tidymodels)
library(lubridate)
library(xgboost)
library(ggplot2)
library(pROC)
library(tidyr)
library(vip)

tidymodels::tidymodels_prefer()

# --- Configuration ---
# Set model_rds_path before sourcing to validate a specific model
if (!exists("model_rds_path")) {
  # Default: find the most recent fitted_model_clean_*.rds
  output_dir <- "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/outputs"
  candidates <- list.files(output_dir, pattern = "^fitted_model_clean_\\d{8}_\\d{4}\\.rds$", full.names = TRUE)
  if (length(candidates) > 0) {
    model_rds_path <- candidates[length(candidates)]  # Most recent
  } else {
    model_rds_path <- file.path(output_dir, "fitted_model_clean.rds")
  }
}

if (!exists("validation_output_dir")) {
  validation_output_dir <- "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/validation_outputs"
}
dir.create(validation_output_dir, showWarnings = FALSE, recursive = TRUE)

run_ts <- format(Sys.time(), "%Y%m%d_%H%M")

cat(sprintf("=== Comprehensive Model Validation ===\n"))
cat(sprintf("Model: %s\n", basename(model_rds_path)))
cat(sprintf("Output dir: %s\n", validation_output_dir))

# --- Load Model ---
m <- readRDS(model_rds_path)
fit_cls <- m$fit
final_wf <- m$wf

cat(sprintf("Model timestamp: %s\n", ifelse(is.null(m$run_timestamp), "legacy", m$run_timestamp)))
cat(sprintf("Training data: %s\n", ifelse(is.null(m$cull_data_file), "unknown", basename(m$cull_data_file))))
cat(sprintf("Predict year: %s\n", ifelse(is.null(m$predict_year), "unknown", m$predict_year)))

# --- Load Data ---
# Always use the latest data for validation
cull_data_file <- "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/data/260529-COTS-Manta-Cull-RHIS-Lawrence-CSIRO.xlsx"
predict_stack_clean <- "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/data/predictors_clean_12layer.tif"
predict_stack_full  <- "C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/data/predictors_terra_30m_full_extended.tif"

raw_cull <- read_excel(cull_data_file, sheet = 4) %>%
  mutate(
    VoyageTitle = as.factor(VoyageTitle),
    SurveyDate  = as.Date(substr(as.character(SurveyDate), 1, 10)),
    Year        = lubridate::year(SurveyDate),
    Total       = Cohort1 + Cohort2 + Cohort3 + Cohort4,
    CPUE        = Total / Bottomtime
  ) %>%
  filter(
    !is.na(Longitude), !is.na(Latitude),
    !is.na(Bottomtime), Bottomtime > 0
  )

cat(sprintf("Total cull records: %d (Year range: %d-%d)\n", nrow(raw_cull), min(raw_cull$Year), max(raw_cull$Year)))

# Load predictor stack — detect whether model uses RG or clean
# Check model predictor names
model_preds <- fit_cls %>% extract_fit_parsnip() %>% purrr::pluck("fit") %>% xgboost::xgb.importance(model = .) %>% pull(Feature)
has_rg <- any(grepl("^RG_", model_preds))

if (has_rg) {
  predictors_reef <- terra::rast(predict_stack_full)
  names(predictors_reef)[1:12] <- c("BPI", "EAST", "HCU", "NORTH", "SLO", "SVF", "VCU", "VRM", "DEM", "AhyaD", "Aspat", "Aten")
  names(predictors_reef)[13:17] <- c("RG_bathy", "RG_slope", "RG_turbid", "RG_waves_Hs", "RG_waves_Tp")
  pred_cols <- names(predictors_reef)
  cat("Predictor stack: Full (17 layers with ReefGuide)\n")
} else {
  predictors_reef <- terra::rast(predict_stack_clean)
  pred_cols <- names(predictors_reef)
  cat("Predictor stack: Clean (12 layers)\n")
}

# Determine if Year is a predictor
use_year <- "Year" %in% model_preds
if (use_year) {
  pred_cols <- c(pred_cols, "Year")
  cat("Model uses Year as predictor\n")
}

# Sectors
sectors_sf <- st_read("C:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/data/SectorShapefile/ltms.shp", quiet = TRUE) %>%
  st_transform(crs = crs(predictors_reef))

# Helper: extract predictors and prepare for a given set of cull records
prepare_validation_set <- function(cull_df, pred_raster, pred_columns, sectors) {
  sf_pts <- st_as_sf(cull_df, coords = c("Longitude", "Latitude"), crs = 4326) %>%
    st_transform(crs(pred_raster))
  
  sf_pts$ID <- seq_len(nrow(sf_pts))
  pred_vals <- terra::extract(pred_raster, terra::vect(sf_pts))
  
  c_sector <- st_join(sf_pts, sectors %>% select(SECTOR_NAM), join = st_nearest_feature)
  sf_pts$SECTOR_NAM <- c_sector$SECTOR_NAM
  
  sf_pts <- sf_pts %>%
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
  
  df_out <- sf_pts %>%
    st_drop_geometry() %>%
    left_join(pred_vals, by = "ID") %>%
    select(-ID) %>%
    mutate(
      cots_problem = factor(if_else(CPUE >= 0.02, 1L, 0L), levels = c("1", "0"))
    ) %>%
    drop_na(all_of(pred_columns))
  
  return(df_out)
}


###################################################################
## LEVEL 1: Internal Cross-Validation Metrics
###################################################################
cat("\n\n========== LEVEL 1: Internal CV Metrics ==========\n")

if (!is.null(m$cv_res)) {
  best_config <- select_best(m$cv_res, metric = "roc_auc")
  cv_metrics <- collect_metrics(m$cv_res) %>%
    filter(.config == best_config$.config) %>%
    select(.metric, mean, std_err)
  
  cat("Best cross-validated metrics:\n")
  print(cv_metrics)
} else {
  cat("No CV results found in model RDS.\n")
}


###################################################################
## LEVEL 2: Out-of-Time Validation (Train ≤ 2025, Test 2026)
###################################################################
cat("\n\n========== LEVEL 2: Out-of-Time Validation ==========\n")

cull_2026 <- raw_cull %>% filter(Year == 2026)
cat(sprintf("2026 holdout records: %d\n", nrow(cull_2026)))

if (nrow(cull_2026) > 100) {
  val_2026 <- prepare_validation_set(cull_2026, predictors_reef, pred_cols, sectors_sf)
  
  if (use_year) {
    val_2026$Year <- 2025  # Predict as if it were 2025 (model's prediction year)
  }
  
  cat(sprintf("Valid records after predictor extraction: %d\n", nrow(val_2026)))
  
  # Predict
  oot_probs <- predict(fit_cls, new_data = val_2026, type = "prob")
  oot_class <- predict(fit_cls, new_data = val_2026, type = "class")
  
  oot_results <- bind_cols(
    val_2026 %>% select(cots_problem, ManagementZone, CPUE),
    oot_probs,
    oot_class
  )
  
  # Overall metrics
  roc_oot <- pROC::roc(response = as.numeric(as.character(oot_results$cots_problem)), 
                        predictor = oot_results$.pred_1, quiet = TRUE)
  cat(sprintf("Out-of-Time AUC (2026 holdout): %.3f\n", auc(roc_oot)))
  
  # Confusion matrix
  cm_oot <- table(
    Predicted = oot_results$.pred_class,
    Observed = oot_results$cots_problem
  )
  cat("Confusion Matrix (2026 holdout):\n")
  print(cm_oot)
  
  accuracy_oot <- sum(diag(cm_oot)) / sum(cm_oot)
  cat(sprintf("Accuracy: %.3f\n", accuracy_oot))
  
  # By sector
  oot_by_zone <- oot_results %>%
    group_by(ManagementZone) %>%
    summarise(
      n = n(),
      n_problem = sum(cots_problem == "1"),
      prevalence = mean(cots_problem == "1"),
      auc = tryCatch(as.numeric(auc(roc(cots_problem, .pred_1, quiet = TRUE))), error = function(e) NA),
      .groups = "drop"
    )
  cat("\nOut-of-Time Performance by Management Zone:\n")
  print(oot_by_zone)
  
  # ROC plot
  roc_df_oot <- data.frame(
    FPR = 1 - roc_oot$specificities,
    TPR = roc_oot$sensitivities
  )
  
  p_roc_oot <- ggplot(roc_df_oot, aes(x = FPR, y = TPR)) +
    geom_line(color = "#38bdf8", linewidth = 1.5) +
    geom_abline(linetype = "dashed", color = "grey50") +
    theme_minimal(base_size = 14) +
    labs(
      title = sprintf("Out-of-Time ROC: 2026 Holdout (AUC = %.3f)", auc(roc_oot)),
      subtitle = sprintf("n = %d cull dives, trained on ≤2025 data", nrow(oot_results)),
      x = "False Positive Rate", y = "True Positive Rate"
    ) +
    coord_equal()
  
  ggsave(file.path(validation_output_dir, sprintf("oot_roc_2026_%s.png", run_ts)),
         plot = p_roc_oot, width = 8, height = 8, dpi = 300)
} else {
  cat("Not enough 2026 data for out-of-time validation.\n")
}


###################################################################
## LEVEL 3: Leave-One-Sector-Out Spatial Transferability
###################################################################
cat("\n\n========== LEVEL 3: Leave-One-Sector-Out ==========\n")

# Aggregate per-year per-site for the full dataset
all_data <- raw_cull %>%
  group_by(ReefName, CullSiteName, Latitude, Longitude, Year) %>%
  summarise(
    CPUE_mean  = sum(Total, na.rm = TRUE) / sum(Bottomtime, na.rm = TRUE),
    CPUE_max   = max(CPUE, na.rm = TRUE),
    Total      = sum(Total, na.rm = TRUE),
    Bottomtime = sum(Bottomtime, na.rm = TRUE),
    n_surveys  = dplyr::n(),
    .groups    = "drop"
  )

all_sf <- st_as_sf(all_data, coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_transform(crs(predictors_reef))

c_sector <- st_join(all_sf, sectors_sf %>% select(SECTOR_NAM), join = st_nearest_feature)
all_sf$ManagementZone <- case_when(
  c_sector$SECTOR_NAM %in% c("Cape Grenville", "Princess Charlotte Bay") ~ "Far Northern",
  c_sector$SECTOR_NAM %in% c("Cooktown / Lizard Island", "Cairns", "Innisfail") ~ "Northern",
  c_sector$SECTOR_NAM %in% c("Townsville", "Cape Upstart", "Whitsunday") ~ "Central",
  c_sector$SECTOR_NAM %in% c("Pompey", "Swain", "Capricorn Bunker") ~ "Southern",
  TRUE ~ "Unknown"
)

# Extract predictions using the fitted model on all data
all_sf$ID <- seq_len(nrow(all_sf))
all_pred_vals <- terra::extract(predictors_reef, terra::vect(all_sf))

all_df <- all_sf %>%
  st_drop_geometry() %>%
  left_join(all_pred_vals, by = "ID") %>%
  mutate(
    cots_problem = factor(if_else(CPUE_mean >= 0.02, 1L, 0L), levels = c("1", "0"))
  ) %>%
  drop_na(all_of(pred_cols))

# Predict using the model on all records
all_probs <- predict(fit_cls, new_data = all_df, type = "prob")
all_df$.pred_1 <- all_probs$.pred_1

# Calculate AUC by sector
loso_results <- all_df %>%
  group_by(ManagementZone) %>%
  summarise(
    n = n(),
    n_problem = sum(cots_problem == "1"),
    prevalence = mean(cots_problem == "1"),
    auc = tryCatch(as.numeric(auc(roc(cots_problem, .pred_1, quiet = TRUE))), error = function(e) NA),
    .groups = "drop"
  )

cat("Spatial Transferability — AUC by Management Zone:\n")
print(loso_results)

# Save
write.csv(loso_results, file.path(validation_output_dir, sprintf("spatial_transferability_%s.csv", run_ts)),
          row.names = FALSE)


###################################################################
## LEVEL 4: Operational Metrics
###################################################################
cat("\n\n========== LEVEL 4: Operational Metrics ==========\n")

# Predict on 2025 operational cull dives
cull_2025 <- raw_cull %>% filter(Year == 2025)
cat(sprintf("2025 operational dives: %d\n", nrow(cull_2025)))

val_2025 <- prepare_validation_set(cull_2025, predictors_reef, pred_cols, sectors_sf)

if (use_year) {
  val_2025$Year <- 2025
}

ops_probs <- predict(fit_cls, new_data = val_2025, type = "prob")
val_2025$.pred_1 <- ops_probs$.pred_1

# Calculate thresholds
p50 <- median(val_2025$.pred_1, na.rm = TRUE)
p88 <- quantile(val_2025$.pred_1, 0.88, na.rm = TRUE)
p95 <- quantile(val_2025$.pred_1, 0.95, na.rm = TRUE)

cat(sprintf("Prediction thresholds — P50: %.4f, P88: %.4f, P95: %.4f\n", p50, p88, p95))

# Ghost cull rates at different thresholds
calc_operational_metrics <- function(df, threshold, label) {
  df %>%
    mutate(
      predicted_high = .pred_1 >= threshold,
      is_ghost = (CPUE == 0),
      is_nonproblematic = (CPUE < 0.02)
    ) %>%
    summarise(
      Threshold = label,
      Threshold_Value = threshold,
      Total_Dives = n(),
      Dives_Suggested = sum(predicted_high),
      Dives_Avoided = Total_Dives - Dives_Suggested,
      Ghost_Avoided = sum(!predicted_high & is_ghost),
      Ghost_Total = sum(is_ghost),
      Ghost_Avoidance_Rate = Ghost_Avoided / Ghost_Total,
      NonProb_Avoided = sum(!predicted_high & is_nonproblematic),
      NonProb_Total = sum(is_nonproblematic),
      NonProb_Avoidance_Rate = NonProb_Avoided / NonProb_Total,
      Hours_Saved = sum(Bottomtime[!predicted_high & is_nonproblematic], na.rm = TRUE) / 60,
      Dollar_Savings = Hours_Saved * 950
    )
}

ops_metrics <- bind_rows(
  calc_operational_metrics(val_2025, p50, "P50 (Median)"),
  calc_operational_metrics(val_2025, p88, "P88"),
  calc_operational_metrics(val_2025, p95, "P95")
)

cat("\nOperational Metrics (2025 Cull Dives):\n")
print(ops_metrics %>% select(Threshold, Total_Dives, Dives_Avoided, Ghost_Avoidance_Rate, 
                               NonProb_Avoidance_Rate, Hours_Saved, Dollar_Savings))

write.csv(ops_metrics, file.path(validation_output_dir, sprintf("operational_metrics_%s.csv", run_ts)),
          row.names = FALSE)


###################################################################
## LEVEL 5: Calibration Plot
###################################################################
cat("\n\n========== LEVEL 5: Calibration ==========\n")

# Use all data predictions for calibration
cal_df <- all_df %>%
  mutate(
    prob_bin = cut(.pred_1, breaks = seq(0, 1, by = 0.05), include.lowest = TRUE),
    observed = as.numeric(as.character(cots_problem))
  ) %>%
  group_by(prob_bin) %>%
  summarise(
    mean_pred = mean(.pred_1, na.rm = TRUE),
    mean_obs  = mean(observed, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  filter(!is.na(prob_bin))

p_cal <- ggplot(cal_df, aes(x = mean_pred, y = mean_obs, size = n)) +
  geom_point(color = "#38bdf8", alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  scale_size_continuous(range = c(2, 12), name = "N observations") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Calibration Plot: Predicted vs Observed COTS Prevalence",
    subtitle = sprintf("All cull data (n = %d), model: %s", nrow(all_df), basename(model_rds_path)),
    x = "Mean Predicted Probability",
    y = "Observed Prevalence (CPUE ≥ 0.02)"
  ) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1))

ggsave(file.path(validation_output_dir, sprintf("calibration_plot_%s.png", run_ts)),
       plot = p_cal, width = 8, height = 8, dpi = 300)

# Partial Dependence Plots for top 5 predictors
cat("Generating variable importance...\n")
booster <- fit_cls %>% extract_fit_parsnip() %>% purrr::pluck("fit")
imp <- xgboost::xgb.importance(model = booster)
top5 <- head(imp$Feature, 5)
cat("Top 5 predictors:", paste(top5, collapse = ", "), "\n")

vip_p <- fit_cls %>%
  extract_fit_parsnip() %>%
  vip::vip(num_features = 15, geom = "col") +
  theme_minimal(base_size = 14) +
  ggtitle(sprintf("Variable Importance — %s", basename(model_rds_path)))

ggsave(file.path(validation_output_dir, sprintf("vip_%s.png", run_ts)),
       plot = vip_p, width = 10, height = 7, dpi = 300)


###################################################################
## Summary Report
###################################################################
cat("\n\n========== VALIDATION SUMMARY ==========\n")
cat(sprintf("Model: %s\n", basename(model_rds_path)))
cat(sprintf("Validation timestamp: %s\n", run_ts))

if (!is.null(m$cv_res)) {
  cat(sprintf("Internal CV AUC: %.3f\n", cv_metrics$mean[cv_metrics$.metric == "roc_auc"]))
}

if (exists("roc_oot") && !is.null(roc_oot)) {
  cat(sprintf("Out-of-Time AUC (2026): %.3f\n", auc(roc_oot)))
}

cat(sprintf("Spatial transferability saved to: spatial_transferability_%s.csv\n", run_ts))
cat(sprintf("Operational metrics saved to: operational_metrics_%s.csv\n", run_ts))
cat(sprintf("Calibration plot saved to: calibration_plot_%s.png\n", run_ts))
cat(sprintf("VIP plot saved to: vip_%s.png\n", run_ts))

# Save summary RDS
validation_summary <- list(
  model_file = model_rds_path,
  model_timestamp = m$run_timestamp,
  validation_timestamp = run_ts,
  internal_cv = if (!is.null(m$cv_res)) cv_metrics else NULL,
  oot_auc = if (exists("roc_oot") && !is.null(roc_oot)) as.numeric(auc(roc_oot)) else NA,
  spatial_transfer = loso_results,
  operational = ops_metrics,
  top_predictors = head(imp$Feature, 10)
)

saveRDS(validation_summary, file.path(validation_output_dir, sprintf("validation_summary_%s.rds", run_ts)))
cat("\nDone! All validation outputs saved.\n")
