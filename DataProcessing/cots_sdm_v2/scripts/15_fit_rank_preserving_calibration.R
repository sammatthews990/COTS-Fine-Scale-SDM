args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v2_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v2_dir, "R", "bootstrap.R"))
ctx <- initialise_v2(script_file)
cfg <- ctx$config

site_year_path <- file.path(
  cfg$paths$artifacts_dir, "site_year_first_visit.rds"
)
site_habitat_path <- file.path(
  cfg$paths$artifacts_dir, "site_habitat_all_sites.rds"
)
prediction_path <- file.path(
  cfg$paths$outputs_dir, "annual_rank_predictions.rds"
)
model_path <- file.path(
  cfg$paths$models_dir, "annual_rank_fold_models.rds"
)
required <- c(
  site_year_path, site_habitat_path, prediction_path, model_path
)
if (any(!file.exists(required))) {
  stop(
    "Missing rank-preserving calibration inputs: ",
    paste(required[!file.exists(required)], collapse = "; ")
  )
}

clean_features <- names(terra::rast(cfg$paths$predictor_stack))
rg_features <- c(
  clean_features,
  "RG_bathy", "RG_slope", "RG_turbid",
  "RG_waves_Hs", "RG_waves_Tp"
)
feature_sets <- list(
  annual_clean_static = clean_features,
  annual_rg_static = rg_features
)

site_covariates <- readRDS(site_habitat_path) |>
  dplyr::filter(.data$evidence_source == "cull_site_history") |>
  dplyr::select(
    "site_id", "site_fold", "clean_complete", "rg_complete",
    "RG_bathy", "RG_slope", "RG_turbid",
    "RG_waves_Hs", "RG_waves_Tp"
  )
site_year <- readRDS(site_year_path) |>
  dplyr::left_join(site_covariates, by = "site_id") |>
  dplyr::group_by(.data$site_id) |>
  dplyr::mutate(
    history_years = dplyr::n_distinct(.data$year),
    observation_weight = 1 / .data$history_years
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    reef_year_id = paste(.data$reef_name, .data$year, sep = "||")
  ) |>
  dplyr::filter(
    .data$clean_complete == 1L,
    .data$rg_complete == 1L
  )

base_predictions <- readRDS(prediction_path) |>
  dplyr::filter(.data$model %in% names(feature_sets))
fold_models <- readRDS(model_path)
calibrated <- list()
calibration_audit <- list()
counter <- 0L

for (model_name in names(feature_sets)) {
  features <- feature_sets[[model_name]]
  for (fold in sort(unique(site_year$site_fold))) {
    counter <- counter + 1L
    train <- dplyr::filter(site_year, .data$site_fold != fold)
    test_prediction <- base_predictions |>
      dplyr::filter(
        .data$model == .env$model_name,
        .data$fold == .env$fold
      )
    model_key <- paste0(model_name, "__fold_", fold)
    model <- fold_models[[model_key]]
    if (is.null(model)) stop("Missing saved fold model: ", model_key)

    train_probability <- predict_xgb_binary(model, train, features)
    calibrator <- fit_group_intercept_calibrator(
      probability = train_probability,
      truth = train$high_density,
      group = train$year,
      weight = train$observation_weight,
      prior_weight = 20
    )
    calibrated[[counter]] <- test_prediction |>
      dplyr::mutate(
        model = paste0(model_name, "_year_intercept"),
        probability = predict_group_intercept_calibrator(
          calibrator, .data$probability, .data$year
        )
      )
    calibration_audit[[counter]] <- calibrator$group_table |>
      dplyr::mutate(
        model = model_name,
        fold = fold,
        global_intercept_shift = calibrator$global_intercept_shift,
        prior_weight = calibrator$prior_weight
      )
  }
}

predictions <- dplyr::bind_rows(
  base_predictions,
  dplyr::bind_rows(calibrated)
) |>
  dplyr::mutate(
    reef_year_fold_id = paste(
      .data$reef_year_id, .data$fold, sep = "||fold||"
    )
  )
metrics <- suppressWarnings(summarise_benchmark_predictions(predictions))
decomposition <- suppressWarnings(
  prediction_scale_decomposition(predictions, group = "reef_year_id")
)
fold_conditional_decomposition <- suppressWarnings(
  prediction_scale_decomposition(
    predictions, group = "reef_year_fold_id"
  )
)
capture <- within_group_capture_metrics(
  predictions, group = "reef_year_id"
)
contrast <- within_group_percentile_contrast(
  predictions, group = "reef_year_id"
)

tables <- list(
  rank_preserving_calibration_metrics = metrics,
  rank_preserving_calibration_decomposition = decomposition,
  rank_preserving_calibration_fold_conditional_decomposition =
    fold_conditional_decomposition,
  rank_preserving_calibration_capture = capture,
  rank_preserving_calibration_percentile_contrast = contrast,
  rank_preserving_calibration_audit =
    dplyr::bind_rows(calibration_audit)
)
for (table_name in names(tables)) {
  utils::write.csv(
    tables[[table_name]],
    file.path(cfg$paths$outputs_dir, paste0(table_name, ".csv")),
    row.names = FALSE
  )
}
saveRDS(
  predictions,
  file.path(
    cfg$paths$outputs_dir,
    "rank_preserving_calibration_predictions.rds"
  )
)

message("Rank-preserving year calibration complete.")
print(dplyr::left_join(
  metrics,
  dplyr::select(decomposition, "model", "within_reef_roc_auc"),
  by = "model"
))
