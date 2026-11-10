# Layer 1: legacy clean + year (operational incumbent), 2026.
#
# The model and raster are frozen products of DataProcessing/fit_predict_cots_sdm.R
# (run 20260820_1041). This script does not refit: it verifies both checksums,
# regenerates a random sample of cells from the saved model to prove the raster
# came from it, and copies the raster byte-for-byte into the deliverable folder.

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "..", "R", "common.R"))
cfg <- initialise_final()
suppressPackageStartupMessages({
  library(tidymodels)
  library(xgboost)
})

legacy <- readRDS(checked_input(cfg, "legacy_clean_model"))
raster_path <- checked_input(cfg, "legacy_clean_raster")
stopifnot(legacy$predict_year == cfg$legacy$predict_year)

prediction <- terra::rast(raster_path)
predictors <- terra::rast(checked_input(cfg, "predictor_stack_clean"))
set.seed(cfg$hmsc$seed)
pts <- terra::spatSample(prediction, size = 2000, method = "random", na.rm = TRUE,
                         as.points = TRUE, values = FALSE, exhaust = TRUE)
environment <- terra::extract(predictors, pts, ID = FALSE)
environment$Year <- legacy$predict_year
regenerated <- predict(legacy$fit, new_data = environment, type = "prob")$.pred_1
stored <- terra::extract(prediction, pts, ID = FALSE)[[1]]
check <- tibble::tibble(
  layer_id = "legacy_clean_year2026",
  sampled_cells = length(stored),
  max_absolute_error = max(abs(regenerated - stored)),
  correlation = stats::cor(regenerated, stored)
)
write_table(check, cfg, "legacy_clean_year2026_reproduction_check.csv")
if (check$max_absolute_error > 1e-6) stop("Frozen clean raster does not match its model.")

out <- file.path(cfg$dirs$layers, "cots_legacy_clean_year2026_30m.tif")
invisible(file.copy(raster_path, out, overwrite = TRUE))
register_layer(
  cfg, "legacy_clean_year2026", out,
  description = "P(CPUE >= 0.02) for survey year 2026; XGBoost on 12-layer clean stack + Year",
  model = "legacy clean + year (fitted_model_clean_20260820_1041)",
  units = "probability",
  extra = list(source_sha256 = cfg$inputs$legacy_clean_raster$sha256)
)
print(check)
