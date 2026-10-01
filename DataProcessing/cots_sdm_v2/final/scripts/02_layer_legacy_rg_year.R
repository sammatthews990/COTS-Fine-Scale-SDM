# Layer 2: legacy RG + year, 2026.
#
# Predicts the frozen ReefGuide-stack model (run 20260820_1217; its raster was
# never written) for survey year 2026. RG layers do not cover every reef, so
# two products are written:
#   - cots_legacy_rg_year2026_30m.tif: RG model only, NA outside RG coverage;
#   - cots_legacy_rg_year2026_clean_fallback_30m.tif: RG where available, else
#     the clean + year layer mapped onto the RG probability scale by a
#     probit-linear fit over the overlap (monotone, so clean ranks are kept),
#     with a companion source layer (1 = RG, 2 = clean fallback).

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "..", "R", "common.R"))
cfg <- initialise_final()
suppressPackageStartupMessages({
  library(tidymodels)
  library(xgboost)
})

legacy <- readRDS(checked_input(cfg, "legacy_rg_model"))
stopifnot(legacy$predict_year == cfg$legacy$predict_year)
stack <- terra::rast(checked_input(cfg, "predictor_stack_rg"))
names(stack) <- cfg$legacy$rg_layer_names

rg_path <- file.path(cfg$dirs$layers, "cots_legacy_rg_year2026_30m.tif")
if (!file.exists(rg_path) || identical(script_args()$force, "1")) {
  log_message("Predicting RG + year for ", legacy$predict_year, " (single core; ~1 h)")
  predict_legacy_workflow(stack, legacy$fit, legacy$predict_year, rg_path)
}
rg <- terra::rast(rg_path)

# Spot-check the written raster against direct model predictions.
set.seed(cfg$hmsc$seed)
pts <- terra::spatSample(rg, size = 2000, method = "random", na.rm = TRUE,
                         as.points = TRUE, values = FALSE, exhaust = TRUE)
environment <- terra::extract(stack, pts, ID = FALSE)
environment$Year <- legacy$predict_year
regenerated <- predict(legacy$fit, new_data = environment, type = "prob")$.pred_1
stored <- terra::extract(rg, pts, ID = FALSE)[[1]]
check <- tibble::tibble(
  layer_id = "legacy_rg_year2026", sampled_cells = length(stored),
  max_absolute_error = max(abs(regenerated - stored)),
  correlation = stats::cor(regenerated, stored)
)
write_table(check, cfg, "legacy_rg_year2026_reproduction_check.csv")
if (check$max_absolute_error > 1e-5) stop("RG raster does not match its model.")

register_layer(
  cfg, "legacy_rg_year2026", rg_path,
  description = "P(CPUE >= 0.02) for survey year 2026; XGBoost on 12 clean + 5 ReefGuide layers + Year; NA outside RG coverage",
  model = "legacy RG + year (fitted_model_withRG_20260820_1217)",
  units = "probability"
)

# Clean-fallback composite.
clean <- terra::rast(file.path(cfg$dirs$layers, "cots_legacy_clean_year2026_30m.tif"))
if (!terra::compareGeom(clean, rg, stopOnError = FALSE)) stop("Clean and RG grids differ.")
overlap <- paired_sample(rg, clean)
calibration <- fit_probit_calibration(primary = overlap$a, fallback = overlap$b)
write_table(tibble::as_tibble(calibration), cfg, "legacy_rg_clean_fallback_calibration.csv")

calibrated_clean <- terra::app(clean, function(p) apply_probit_calibration(p, calibration))
composite <- terra::cover(rg, calibrated_clean)
names(composite) <- "cots_prob"
source_layer <- terra::ifel(!is.na(rg), 1L, terra::ifel(!is.na(clean), 2L, NA))
names(source_layer) <- "prediction_source"

composite_path <- file.path(cfg$dirs$layers, "cots_legacy_rg_year2026_clean_fallback_30m.tif")
source_path <- file.path(cfg$dirs$layers, "cots_legacy_rg_year2026_clean_fallback_source_30m.tif")
write_layer(composite, composite_path)
write_layer(source_layer, source_path, datatype = "INT1U")

coverage <- terra::freq(source_layer)
coverage$source <- c("1" = "rg", "2" = "clean_fallback")[as.character(coverage$value)]
coverage$fraction <- coverage$count / sum(coverage$count)
write_table(coverage, cfg, "legacy_rg_clean_fallback_coverage.csv")

register_layer(
  cfg, "legacy_rg_year2026_clean_fallback", composite_path,
  description = "RG + year where RG layers exist, else clean + year calibrated to the RG scale",
  model = "legacy RG + year with clean + year fallback",
  units = "probability",
  extra = list(rg_fraction = coverage$fraction[coverage$value == 1])
)
register_layer(
  cfg, "legacy_rg_year2026_clean_fallback_source", source_path,
  description = "Source of each composite cell: 1 = RG model, 2 = clean fallback",
  model = "legacy RG + year with clean + year fallback",
  units = "category"
)
print(check)
print(coverage)
