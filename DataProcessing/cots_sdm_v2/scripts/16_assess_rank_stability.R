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
original_prediction_path <- file.path(
  cfg$paths$outputs_dir, "annual_rank_predictions.rds"
)
required <- c(
  site_year_path, site_habitat_path, original_prediction_path
)
if (any(!file.exists(required))) {
  stop(
    "Missing rank-stability inputs: ",
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
    "site_id", "reef_name", "validation_group_id",
    "clean_complete", "rg_complete",
    "RG_bathy", "RG_slope", "RG_turbid",
    "RG_waves_Hs", "RG_waves_Tp"
  )
base_site_year <- readRDS(site_year_path) |>
  dplyr::select(-dplyr::any_of("site_fold")) |>
  dplyr::left_join(
    dplyr::select(site_covariates, -"reef_name"),
    by = "site_id"
  ) |>
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

original_seed <- cfg$validation$spatial_seed
additional_seeds <- original_seed + c(1L, 2L)
prediction_sets <- list(
  readRDS(original_prediction_path) |>
    dplyr::filter(.data$model %in% names(feature_sets)) |>
    dplyr::mutate(resample_seed = original_seed)
)
model_sets <- list()

for (seed in additional_seeds) {
  message("Running annual rank stability seed ", seed, "...")
  reassigned <- assign_within_reef_site_folds(
    dplyr::select(
      site_covariates,
      "site_id", "reef_name", "validation_group_id"
    ),
    k = cfg$site_habitat$within_reef_folds,
    seed = seed
  ) |>
    dplyr::select("site_id", "site_fold")
  site_year <- base_site_year |>
    dplyr::left_join(reassigned, by = "site_id")
  model_cfg <- cfg$model
  model_cfg$seed <- seed
  result <- run_site_year_rank_benchmark(
    site_year, feature_sets, model_cfg
  )
  prediction_sets[[length(prediction_sets) + 1L]] <-
    result$predictions |>
    dplyr::mutate(resample_seed = seed)
  model_sets[[as.character(seed)]] <- result$models
}

predictions <- dplyr::bind_rows(prediction_sets)
seeds <- sort(unique(predictions$resample_seed))
metric_tables <- list()
decomposition_tables <- list()
capture_tables <- list()
contrast_tables <- list()
for (seed in seeds) {
  x <- dplyr::filter(
    predictions, .data$resample_seed == .env$seed
  )
  metric_tables[[as.character(seed)]] <- suppressWarnings(
    summarise_benchmark_predictions(x)
  ) |>
    dplyr::mutate(resample_seed = seed)
  decomposition_tables[[as.character(seed)]] <- suppressWarnings(
    prediction_scale_decomposition(x, group = "reef_year_id")
  ) |>
    dplyr::mutate(resample_seed = seed)
  capture_tables[[as.character(seed)]] <-
    within_group_capture_metrics(x, group = "reef_year_id") |>
    dplyr::mutate(resample_seed = seed)
  contrast_tables[[as.character(seed)]] <-
    within_group_percentile_contrast(
      x, group = "reef_year_id"
    ) |>
    dplyr::mutate(resample_seed = seed)
}

metrics <- dplyr::bind_rows(metric_tables)
decomposition <- dplyr::bind_rows(decomposition_tables)
capture <- dplyr::bind_rows(capture_tables)
contrast <- dplyr::bind_rows(contrast_tables)
stability_summary <- decomposition |>
  dplyr::left_join(
    dplyr::filter(
      capture, abs(.data$target_fraction - 0.20) < 1e-9
    ) |>
      dplyr::select("model", "resample_seed", "capture", "lift"),
    by = c("model", "resample_seed")
  ) |>
  dplyr::left_join(
    dplyr::select(
      contrast,
      "model", "resample_seed", "positive_negative_percentile_gap"
    ),
    by = c("model", "resample_seed")
  ) |>
  dplyr::group_by(.data$model) |>
  dplyr::summarise(
    seeds = dplyr::n(),
    within_reef_year_auc_mean = mean(.data$within_reef_roc_auc),
    within_reef_year_auc_sd = stats::sd(.data$within_reef_roc_auc),
    within_reef_year_auc_min = min(.data$within_reef_roc_auc),
    within_reef_year_auc_max = max(.data$within_reef_roc_auc),
    top20_capture_mean = mean(.data$capture),
    top20_capture_sd = stats::sd(.data$capture),
    top20_lift_mean = mean(.data$lift),
    percentile_gap_mean = mean(
      .data$positive_negative_percentile_gap
    ),
    .groups = "drop"
  )

tables <- list(
  annual_rank_stability_metrics = metrics,
  annual_rank_stability_decomposition = decomposition,
  annual_rank_stability_capture = capture,
  annual_rank_stability_percentile_contrast = contrast,
  annual_rank_stability_summary = stability_summary
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
  file.path(cfg$paths$outputs_dir, "annual_rank_stability_predictions.rds")
)
saveRDS(
  model_sets,
  file.path(cfg$paths$models_dir, "annual_rank_stability_models.rds")
)

message("Annual rank stability assessment complete.")
print(stability_summary)
