############################################################
# fit_predict_ensemble_stack.R
# Script to build, cross-validate, and predict Ensemble / Stacked
# Species Distribution Models (SDMs) for COTS across the GBR.
############################################################

library(terra)
library(sf)
library(dplyr)
library(readxl)
library(pROC)
library(yardstick)
library(tidyr)
library(xgboost)

out_dir <- "c:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/outputs"
data_dir <- "c:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/data"
cull_xlsx <- file.path(data_dir, "250929_COTS-Manta-Cull-RHIS-Data-Matthews-and-Schlawinsky.xlsx")

cat("=== 1. Loading Base Model Rasters ===\n")
r_cull_clean  <- terra::rast(file.path(out_dir, "COTS_prob_0.02_cpue_year2025_clean.tif"))
r_cull_rg     <- terra::rast(file.path(out_dir, "COTS_prob_0.02_cpue_year2025_withRG.tif"))
r_manta_clean <- terra::rast(file.path(out_dir, "COTS_prob_manta_clean.tif"))
r_manta_rg    <- terra::rast(file.path(out_dir, "COTS_prob_manta_withRG.tif"))
r_max_clean   <- terra::rast(file.path(out_dir, "COTS_maxent_suitability_clean.tif"))
r_max_rg      <- terra::rast(file.path(out_dir, "COTS_maxent_suitability_reefguide.tif"))

base_rasters <- list(
  cull_rg     = r_cull_rg,
  cull_clean  = r_cull_clean,
  max_rg      = r_max_rg,
  manta_clean = r_manta_clean,
  max_clean   = r_max_clean,
  manta_rg    = r_manta_rg
)

cat("=== 2. Extracting Base Predictions on Cull Validation Dataset ===\n")
cull_raw <- read_excel(cull_xlsx, sheet = 4) %>%
  filter(!is.na(Longitude), !is.na(Latitude), !is.na(Bottomtime), Bottomtime > 0, Longitude != 0, Latitude != 0) %>%
  filter(as.Date(SurveyDate) > as.Date("2018-11-30")) %>%
  mutate(Total = Cohort1 + Cohort2 + Cohort3 + Cohort4, CPUE = Total / Bottomtime, ReefName = as.character(ReefName))

reef_select_eval <- cull_raw %>% 
  filter(as.Date(SurveyDate) >= as.Date("2025-01-01")) %>%
  group_by(ReefName) %>%
  summarise(Total_Hours = sum(Bottomtime, na.rm=TRUE) / 60) %>%
  arrange(desc(Total_Hours)) %>% slice_head(n = 30) %>% pull(ReefName)

cull_sub <- cull_raw %>% filter(ReefName %in% reef_select_eval, as.Date(SurveyDate) >= as.Date("2025-01-01"))
cull_pts <- st_as_sf(cull_sub, coords = c("Longitude", "Latitude"), crs = 4326)
cull_pts_proj <- st_transform(cull_pts, crs = crs(r_cull_clean))
v_pts <- terra::vect(cull_pts_proj)

df_val <- tibble(
  truth       = as.numeric(pull(cull_sub, CPUE) >= 0.02),
  ReefName    = pull(cull_sub, ReefName),
  cull_rg     = as.numeric(terra::extract(r_cull_rg, v_pts)[[2]]),
  cull_clean  = as.numeric(terra::extract(r_cull_clean, v_pts)[[2]]),
  max_rg      = as.numeric(terra::extract(r_max_rg, v_pts)[[2]]),
  manta_clean = as.numeric(terra::extract(r_manta_clean, v_pts)[[2]]),
  max_clean   = as.numeric(terra::extract(r_max_clean, v_pts)[[2]]),
  manta_rg    = as.numeric(terra::extract(r_manta_rg, v_pts)[[2]])
) %>% drop_na()

y_val <- pull(df_val, truth)
cat("Validation sample size:", nrow(df_val), "cull dives across top 30 reefs\n\n")

# --- Strategy A: Performance-Weighted Average ---
auc_vals <- c(
  cull_rg     = 0.865,
  cull_clean  = 0.660,
  max_rg      = 0.559,
  manta_clean = 0.535,
  max_clean   = 0.534,
  manta_rg    = 0.499
)

w_auc <- (pmax(0, auc_vals - 0.5))^2
w_auc <- w_auc / sum(w_auc)
names(w_auc) <- names(auc_vals)

cat("=== Performance Ensemble Weights ===\n")
print(round(w_auc, 4))
cat("\n")

m_matrix <- as.matrix(df_val %>% select(all_of(names(w_auc))))
pred_weighted_ens <- as.numeric(m_matrix %*% w_auc)

auc_weighted <- as.numeric(pROC::auc(pROC::roc(response = y_val, predictor = pred_weighted_ens, quiet = TRUE)))
cat("Performance-Weighted Ensemble AUC:", round(auc_weighted, 3), "\n\n")

# --- Strategy B: Spatially Cross-Validated Logistic / Ridge Meta-Learner ---
set.seed(123)
# Group folds by ReefName to evaluate true Spatial CV
unique_reefs <- unique(df_val$ReefName)
reef_folds   <- sample(rep(1:5, length.out = length(unique_reefs)))
names(reef_folds) <- unique_reefs

df_val$fold <- reef_folds[df_val$ReefName]

oof_logistic <- numeric(nrow(df_val))
oof_xgb      <- numeric(nrow(df_val))

features <- names(w_auc)

for (k in 1:5) {
  train_df <- df_val %>% filter(fold != k)
  test_df  <- df_val %>% filter(fold == k)
  test_idx <- which(df_val$fold == k)
  
  # Fit Logistic Meta-Learner
  meta_log <- glm(truth ~ cull_rg + cull_clean + max_rg + manta_clean + max_clean + manta_rg,
                  data = train_df, family = binomial())
  oof_logistic[test_idx] <- predict(meta_log, newdata = test_df, type = "response")
  
  # Fit XGBoost Meta-Learner (shallow tree, heavy regularization to prevent overfitting)
  dtrain <- xgb.DMatrix(data = as.matrix(train_df %>% select(all_of(features))), label = train_df$truth)
  dtest  <- xgb.DMatrix(data = as.matrix(test_df %>% select(all_of(features))), label = test_df$truth)
  
  params <- list(
    booster = "gbtree",
    objective = "binary:logistic",
    eval_metric = "auc",
    max_depth = 2,
    eta = 0.05,
    subsample = 0.8,
    colsample_bytree = 0.8
  )
  
  meta_xgb <- xgb.train(params = params, data = dtrain, nrounds = 40, verbose = 0)
  oof_xgb[test_idx] <- predict(meta_xgb, dtest)
}

auc_logistic_oof <- as.numeric(pROC::auc(pROC::roc(response = y_val, predictor = oof_logistic, quiet = TRUE)))
auc_xgb_oof      <- as.numeric(pROC::auc(pROC::roc(response = y_val, predictor = oof_xgb, quiet = TRUE)))

cat("=== Spatial Out-of-Fold (OOF) Ensemble Meta-Learner Performance ===\n")
cat("Spatial OOF Stacked Logistic Meta-Model AUC:", round(auc_logistic_oof, 3), "\n")
cat("Spatial OOF Stacked XGBoost Meta-Model AUC: ", round(auc_xgb_oof, 3), "\n\n")

cat("=== 3. Generating Spatial GBR Ensemble Raster (`COTS_prob_ensemble.tif`) ===\n")

# Combine rasters using the performance-weighted strategy for stable, leak-free spatial extrapolation
r_ensemble <- r_cull_rg * w_auc["cull_rg"]
for (m_name in setdiff(names(w_auc), "cull_rg")) {
  if (w_auc[m_name] > 0) {
    # Resample to match r_cull_rg geometry if necessary
    r_sub <- base_rasters[[m_name]]
    if (!terra::compareGeom(r_ensemble, r_sub, stopOnError = FALSE)) {
      r_sub <- terra::resample(r_sub, r_ensemble, method = "bilinear")
    }
    r_ensemble <- r_ensemble + (r_sub * w_auc[m_name])
  }
}

names(r_ensemble) <- "ensemble_prob"
ensemble_tif_path <- file.path(out_dir, "COTS_prob_ensemble.tif")
terra::writeRaster(r_ensemble, ensemble_tif_path, overwrite = TRUE)
cat("Ensemble raster successfully written to:", ensemble_tif_path, "\n\n")

# Save summary metrics
ensemble_summary <- list(
  weights = w_auc,
  auc_weighted = auc_weighted,
  auc_logistic_oof = auc_logistic_oof,
  auc_xgb_oof = auc_xgb_oof
)

saveRDS(ensemble_summary, file.path(out_dir, "fitted_ensemble_model.rds"))
cat("Ensemble summary saved to outputs/fitted_ensemble_model.rds\n")
