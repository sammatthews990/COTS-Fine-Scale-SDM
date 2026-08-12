############################################################
# fit_predict_ensemble_stack.R
# Script to build, cross-validate, and predict Ensemble / Stacked
# Species Distribution Models (SDMs) for COTS across the GBR
# using ALL historical cull dive data (2018-2025) and time-agnostic models.
############################################################

library(terra)
library(sf)
library(dplyr)
library(readxl)
library(pROC)
library(yardstick)
library(tidyr)
library(xgboost)

out_dir   <- "c:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/outputs"
data_dir  <- "c:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM/DataProcessing/data"
cull_xlsx <- file.path(data_dir, "250929_COTS-Manta-Cull-RHIS-Data-Matthews-and-Schlawinsky.xlsx")

cat("=== 1. Loading Base Model Rasters ===\n")
r_cull_2025_clean <- terra::rast(file.path(out_dir, "COTS_prob_0.02_cpue_year2025_clean.tif"))
r_cull_2025_rg    <- terra::rast(file.path(out_dir, "COTS_prob_0.02_cpue_year2025_withRG.tif"))
r_cull_agn_clean  <- terra::rast(file.path(out_dir, "COTS_prob_0.02_cpue_agnostic_clean.tif"))
r_cull_agn_rg     <- terra::rast(file.path(out_dir, "COTS_prob_0.02_cpue_agnostic_withRG.tif"))
r_manta_clean     <- terra::rast(file.path(out_dir, "COTS_prob_manta_clean.tif"))
r_manta_rg        <- terra::rast(file.path(out_dir, "COTS_prob_manta_withRG.tif"))
r_max_clean       <- terra::rast(file.path(out_dir, "COTS_maxent_suitability_clean.tif"))
r_max_rg          <- terra::rast(file.path(out_dir, "COTS_maxent_suitability_reefguide.tif"))

base_rasters <- list(
  cull_2025_rg    = r_cull_2025_rg,
  cull_2025_clean = r_cull_2025_clean,
  cull_agn_rg     = r_cull_agn_rg,
  cull_agn_clean  = r_cull_agn_clean,
  max_rg          = r_max_rg,
  max_clean       = r_max_clean,
  manta_rg        = r_manta_rg,
  manta_clean     = r_manta_clean
)

cat("=== 2. Extracting Base Predictions on ALL Historical Cull Dives (2018-2025) ===\n")
cull_raw <- read_excel(cull_xlsx, sheet = 4) %>%
  filter(!is.na(Longitude), !is.na(Latitude), !is.na(Bottomtime), Bottomtime > 0, Longitude != 0, Latitude != 0) %>%
  filter(as.Date(SurveyDate) > as.Date("2018-11-30")) %>%
  mutate(Total = Cohort1 + Cohort2 + Cohort3 + Cohort4, CPUE = Total / Bottomtime, ReefName = as.character(ReefName))

cat("Total historical cull dives loaded:", nrow(cull_raw), "\n")

cull_pts <- st_as_sf(cull_raw, coords = c("Longitude", "Latitude"), crs = 4326)
cull_pts_proj <- st_transform(cull_pts, crs = crs(r_cull_2025_clean))
v_pts <- terra::vect(cull_pts_proj)

df_all <- tibble(
  truth           = as.numeric(pull(cull_raw, CPUE) >= 0.02),
  CPUE            = pull(cull_raw, CPUE),
  ReefName        = pull(cull_raw, ReefName),
  SurveyDate      = pull(cull_raw, SurveyDate),
  cull_2025_rg    = as.numeric(terra::extract(r_cull_2025_rg, v_pts)[[2]]),
  cull_2025_clean = as.numeric(terra::extract(r_cull_2025_clean, v_pts)[[2]]),
  cull_agn_rg     = as.numeric(terra::extract(r_cull_agn_rg, v_pts)[[2]]),
  cull_agn_clean  = as.numeric(terra::extract(r_cull_agn_clean, v_pts)[[2]]),
  max_rg          = as.numeric(terra::extract(r_max_rg, v_pts)[[2]]),
  max_clean       = as.numeric(terra::extract(r_max_clean, v_pts)[[2]]),
  manta_rg        = as.numeric(terra::extract(r_manta_rg, v_pts)[[2]]),
  manta_clean     = as.numeric(terra::extract(r_manta_clean, v_pts)[[2]])
) %>% drop_na()

cat("Valid historical cull dives with predictions:", nrow(df_all), "\n\n")

y_all <- pull(df_all, truth)

# --- Calculate Base Model AUCs across All Historical Dives ---
model_names <- c("cull_2025_rg", "cull_2025_clean", "cull_agn_rg", "cull_agn_clean", 
                 "max_rg", "max_clean", "manta_rg", "manta_clean")

auc_all <- numeric(length(model_names))
names(auc_all) <- model_names

for (m in model_names) {
  auc_all[m] <- as.numeric(pROC::auc(pROC::roc(response = y_all, predictor = df_all[[m]], quiet = TRUE)))
}

cat("=== AUC across ALL Historical Cull Dives (2018-2025) ===\n")
print(round(auc_all, 3))
cat("\n")

# Calculate Performance Ensemble Weights (All Dives)
w_auc_all <- (pmax(0, auc_all - 0.5))^2
names(w_auc_all) <- model_names
w_auc_all <- w_auc_all / sum(w_auc_all)

cat("=== Performance Ensemble Weights (All Dives) ===\n")
print(round(w_auc_all, 4))
cat("\n")

m_matrix_all <- as.matrix(df_all %>% select(all_of(model_names)))
pred_ens_all <- as.numeric(m_matrix_all %*% w_auc_all)

auc_weighted_all <- as.numeric(pROC::auc(pROC::roc(response = y_all, predictor = pred_ens_all, quiet = TRUE)))
cat("All-Historical Performance-Weighted Ensemble AUC:", round(auc_weighted_all, 3), "\n\n")

# --- Strategy B: Spatial Out-of-Fold (OOF) Stacking ---
set.seed(123)
unique_reefs <- unique(df_all$ReefName)
reef_folds   <- sample(rep(1:5, length.out = length(unique_reefs)))
names(reef_folds) <- unique_reefs

df_all$fold <- reef_folds[df_all$ReefName]

oof_logistic <- numeric(nrow(df_all))
oof_xgb      <- numeric(nrow(df_all))

for (k in 1:5) {
  train_df <- df_all %>% filter(fold != k)
  test_df  <- df_all %>% filter(fold == k)
  test_idx <- which(df_all$fold == k)
  
  meta_log <- glm(truth ~ cull_2025_rg + cull_2025_clean + cull_agn_rg + cull_agn_clean + max_rg + max_clean + manta_rg + manta_clean,
                  data = train_df, family = binomial())
  oof_logistic[test_idx] <- predict(meta_log, newdata = test_df, type = "response")
  
  dtrain <- xgb.DMatrix(data = as.matrix(train_df %>% select(all_of(model_names))), label = train_df$truth)
  dtest  <- xgb.DMatrix(data = as.matrix(test_df %>% select(all_of(model_names))), label = test_df$truth)
  
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

auc_logistic_oof <- as.numeric(pROC::auc(pROC::roc(response = y_all, predictor = oof_logistic, quiet = TRUE)))
auc_xgb_oof      <- as.numeric(pROC::auc(pROC::roc(response = y_all, predictor = oof_xgb, quiet = TRUE)))

cat("=== Spatial Out-of-Fold (OOF) Ensemble Performance (All Historical Dives) ===\n")
cat("Spatial OOF Stacked Logistic Meta-Model AUC:", round(auc_logistic_oof, 3), "\n")
cat("Spatial OOF Stacked XGBoost Meta-Model AUC: ", round(auc_xgb_oof, 3), "\n\n")

cat("=== 3. Writing Ensemble Prediction Rasters ===\n")

# 1. Operational 2025 Ensemble Raster (`COTS_prob_ensemble_2025.tif`)
r_ens_2025 <- r_cull_2025_rg * unname(w_auc_all["cull_2025_rg"]) +
              r_cull_2025_clean * unname(w_auc_all["cull_2025_clean"])

for (m_name in c("max_rg", "max_clean", "manta_rg", "manta_clean")) {
  w_val <- unname(w_auc_all[m_name])
  if (w_val > 0) {
    r_sub <- base_rasters[[m_name]]
    if (!terra::compareGeom(r_ens_2025, r_sub, stopOnError = FALSE)) {
      r_sub <- terra::resample(r_sub, r_ens_2025, method = "bilinear")
    }
    r_ens_2025 <- r_ens_2025 + (r_sub * w_val)
  }
}
names(r_ens_2025) <- "ensemble_prob_2025"
terra::writeRaster(r_ens_2025, file.path(out_dir, "COTS_prob_ensemble_2025.tif"), overwrite = TRUE)
terra::writeRaster(r_ens_2025, file.path(out_dir, "COTS_prob_ensemble.tif"), overwrite = TRUE) # Default alias
cat("2025 Ensemble raster saved to outputs/COTS_prob_ensemble_2025.tif\n")

# 2. Outbreak Persistence Time-Agnostic Ensemble Raster (`COTS_prob_ensemble_agnostic.tif`)
w_agnostic <- c(
  cull_agn_rg    = unname(w_auc_all["cull_agn_rg"]) + unname(w_auc_all["cull_2025_rg"])*0.5,
  cull_agn_clean = unname(w_auc_all["cull_agn_clean"]) + unname(w_auc_all["cull_2025_clean"])*0.5,
  max_rg         = unname(w_auc_all["max_rg"]),
  max_clean      = unname(w_auc_all["max_clean"]),
  manta_rg       = unname(w_auc_all["manta_rg"]),
  manta_clean    = unname(w_auc_all["manta_clean"])
)
w_agnostic <- w_agnostic / sum(w_agnostic)

r_ens_agnostic <- r_cull_agn_rg * unname(w_agnostic["cull_agn_rg"]) +
                  r_cull_agn_clean * unname(w_agnostic["cull_agn_clean"])

for (m_name in c("max_rg", "max_clean", "manta_rg", "manta_clean")) {
  w_val <- unname(w_agnostic[m_name])
  if (w_val > 0) {
    r_sub <- base_rasters[[m_name]]
    if (!terra::compareGeom(r_ens_agnostic, r_sub, stopOnError = FALSE)) {
      r_sub <- terra::resample(r_sub, r_ens_agnostic, method = "bilinear")
    }
    r_ens_agnostic <- r_ens_agnostic + (r_sub * w_val)
  }
}
names(r_ens_agnostic) <- "ensemble_prob_agnostic"
terra::writeRaster(r_ens_agnostic, file.path(out_dir, "COTS_prob_ensemble_agnostic.tif"), overwrite = TRUE)
cat("Agnostic Persistence Ensemble raster saved to outputs/COTS_prob_ensemble_agnostic.tif\n\n")

# Save summary metrics
ensemble_summary <- list(
  auc_all = auc_all,
  weights_all = w_auc_all,
  auc_weighted_all = auc_weighted_all,
  auc_logistic_oof = auc_logistic_oof,
  auc_xgb_oof = auc_xgb_oof,
  sample_size_all = nrow(df_all)
)

saveRDS(ensemble_summary, file.path(out_dir, "fitted_ensemble_model.rds"))
cat("Ensemble summary saved to outputs/fitted_ensemble_model.rds\n")
