library(terra)
library(tidymodels)
library(xgboost)

fitted_obj <- readRDS("DataProcessing/outputs/fitted_model_clean_pseudorepsN5.rds")
fit_cls <- fitted_obj$fit

predictors_reef <- terra::rast("DataProcessing/data/predictors_clean_12layer.tif")
names(predictors_reef) <- c("BPI", "EAST", "HCU", "NORTH", "SLO", "SVF", "VCU", "VRM", "DEM", "AhyaD", "Aspat", "Aten")

pred_fun <- function(model, newdata) {
  newdata <- as.data.frame(newdata)
  newdata$Year <- 2025
  preds <- predict(model, new_data = newdata, type = "prob")
  return(as.numeric(preds$.pred_1))
}

cat("Predicting spatial raster grid...\n")
prob_raster <- terra::predict(predictors_reef, fit_cls, fun = pred_fun, na.rm = TRUE)
names(prob_raster) <- "cots_prob"

output_tif <- "DataProcessing/outputs/COTS_prob_0.02_cpue_year2025_clean_pseudorepsN5.tif"
terra::writeRaster(prob_raster, output_tif, overwrite = TRUE, gdal = c("COMPRESS=LZW", "TILED=YES"))
cat("Successfully generated:", output_tif, "\n")
