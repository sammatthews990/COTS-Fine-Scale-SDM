args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v2_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v2_dir, "R", "bootstrap.R"))
ctx <- initialise_v2(script_file)
cfg <- ctx$config

suppressPackageStartupMessages({
  library(tidymodels)
  library(xgboost)
})

model_path <- file.path(
  cfg$paths$repo_root, "DataProcessing", "outputs",
  "fitted_model_clean_20260820_1041.rds"
)
raster_path <- cfg$paths$legacy_2026_raster
input_paths <- c(
  model = model_path,
  prediction = raster_path,
  predictors = cfg$paths$predictor_stack,
  observations = cfg$paths$cull_workbook
)
if (any(!file.exists(input_paths))) {
  stop("Missing legacy baseline input(s): ", paste(input_paths[!file.exists(input_paths)], collapse = "; "))
}

legacy <- readRDS(model_path)
parsnip_fit <- workflows::extract_fit_parsnip(legacy$fit)
args_fit <- parsnip_fit$spec$args
predictors <- terra::rast(cfg$paths$predictor_stack)
prediction <- terra::rast(raster_path)

set.seed(20260831)
sample_points <- terra::spatSample(
  prediction,
  size = 2000,
  method = "random",
  na.rm = TRUE,
  as.points = TRUE,
  values = FALSE,
  exhaust = TRUE
)
sample_environment <- terra::extract(predictors, sample_points, ID = FALSE)
sample_environment$Year <- legacy$predict_year
regenerated <- predict(
  legacy$fit,
  new_data = as.data.frame(sample_environment),
  type = "prob"
)$.pred_1
stored <- terra::extract(prediction, sample_points, ID = FALSE)[[1]]
valid <- is.finite(regenerated) & is.finite(stored)
verification <- tibble::tibble(
  sampled_cells = sum(valid),
  correlation = stats::cor(regenerated[valid], stored[valid]),
  mean_absolute_error = mean(abs(regenerated[valid] - stored[valid])),
  maximum_absolute_error = max(abs(regenerated[valid] - stored[valid])),
  exact_matches = sum(regenerated[valid] == stored[valid]),
  exact_match_fraction = mean(regenerated[valid] == stored[valid])
)

sha256 <- function(path) {
  if (!requireNamespace("digest", quietly = TRUE)) return(NA_character_)
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
file_manifest <- tibble::tibble(
  role = names(input_paths),
  path = normalizePath(input_paths, winslash = "/"),
  bytes = file.info(input_paths)$size,
  sha256 = vapply(input_paths, sha256, character(1))
)

hyperparameters <- tibble::tibble(
  trees = as.integer(rlang::eval_tidy(args_fit$trees)),
  tree_depth = as.integer(rlang::eval_tidy(args_fit$tree_depth)),
  learn_rate = as.numeric(rlang::eval_tidy(args_fit$learn_rate)),
  loss_reduction = as.numeric(rlang::eval_tidy(args_fit$loss_reduction)),
  min_n = as.numeric(rlang::eval_tidy(args_fit$min_n)),
  predict_year = legacy$predict_year,
  threshold = cfg$response$cpue_threshold
)

manifest <- list(
  baseline_name = "legacy_clean_year_20260820_1041",
  target_definition = paste(
    "annual site aggregation: sum(Cohort1:4) / sum(Bottomtime) >=",
    cfg$response$cpue_threshold
  ),
  aggregation_keys = c(
    "ReefName", "CullSiteName", "Latitude", "Longitude", "Year"
  ),
  predictors = c(names(predictors), "Year"),
  prediction_year = legacy$predict_year,
  hyperparameters = hyperparameters,
  files = file_manifest,
  raster_geometry = list(
    dimensions = dim(prediction),
    resolution = terra::res(prediction),
    extent = as.vector(terra::ext(prediction)),
    crs = terra::crs(prediction)
  ),
  verification = verification,
  package_versions = vapply(
    c("R", "terra", "tidymodels", "xgboost", "recipes", "workflows"),
    function(x) if (x == "R") as.character(getRversion()) else as.character(utils::packageVersion(x)),
    character(1)
  )
)

saveRDS(
  manifest,
  file.path(cfg$paths$artifacts_dir, "legacy_clean_year_manifest.rds")
)
utils::write.csv(
  file_manifest,
  file.path(cfg$paths$outputs_dir, "legacy_clean_year_file_manifest.csv"),
  row.names = FALSE
)
utils::write.csv(
  hyperparameters,
  file.path(cfg$paths$outputs_dir, "legacy_clean_year_hyperparameters.csv"),
  row.names = FALSE
)
utils::write.csv(
  verification,
  file.path(cfg$paths$outputs_dir, "legacy_clean_year_reproduction_check.csv"),
  row.names = FALSE
)

if (verification$maximum_absolute_error != 0) {
  stop("Legacy raster reproduction was not exact.")
}
message("Legacy clean-year baseline frozen and reproduced exactly.")
print(verification)
