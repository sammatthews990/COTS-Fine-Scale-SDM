# DataProcessing/run_prediction_only.R
# Predict COTS suitability across full domain using saved 12-layer MaxEnt model

library(terra)
library(maxnet)

project_dir <- "c:/Users/smatthew/Documents/GitKraken/COTS Fine Scale SDM"
data_dir    <- file.path(project_dir, "DataProcessing/data")
out_dir     <- file.path(project_dir, "DataProcessing/outputs")

predict_stack_full <- file.path(data_dir, "predictors_terra_30m_full_extended.tif")
output_model       <- file.path(out_dir, "maxent_model_cots.rds")
output_tif         <- file.path(out_dir, "COTS_maxent_suitability.tif")
output_tif_clean   <- file.path(out_dir, "COTS_maxent_suitability_clean.tif")

saved_obj  <- readRDS(output_model)
best_model <- saved_obj$model

print("Loading 12-layer clean predictor stack...")
predictors_all <- terra::rast(predict_stack_full)
names(predictors_all)[1:12] <- c("BPI", "EAST", "HCU", "NORTH", "SLO", "SVF",
                                  "VCU", "VRM", "DEM", "AhyaD", "Aspat", "Aten")
predictors <- predictors_all[[1:12]]
pred_names <- names(predictors)
print(paste("Predictors:", paste(pred_names, collapse = ", ")))

pred_fun <- function(model, newdata) {
  newdata <- as.data.frame(newdata)
  names(newdata) <- pred_names
  p <- predict(model, newdata, type = "cloglog")
  as.numeric(p)
}

print("Predicting suitability raster stack across full GBR & Torres Strait domain...")
suitability <- terra::predict(
  predictors,
  best_model,
  fun   = pred_fun,
  na.rm = TRUE,
  cores = 1
)
names(suitability) <- "suitability"

print("Writing COTS_maxent_suitability_clean.tif...")
if (file.exists(output_tif_clean)) file.remove(output_tif_clean)
writeRaster(
  suitability,
  filename  = output_tif_clean,
  overwrite = TRUE,
  gdal      = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES")
)

print("Writing COTS_maxent_suitability.tif...")
if (file.exists(output_tif)) file.remove(output_tif)
writeRaster(
  suitability,
  filename  = output_tif,
  overwrite = TRUE,
  gdal      = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES")
)

print("Raster prediction completed successfully!")
suit_vals <- terra::global(suitability, fun = c("mean", "min", "max", "sd"), na.rm = TRUE)
print(suit_vals)
