args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v2_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v2_dir, "R", "bootstrap.R"))
ctx <- initialise_v2(script_file)
cfg <- ctx$config

site_year_path <- file.path(cfg$paths$artifacts_dir, "site_year_all.rds")
site_habitat_path <- file.path(cfg$paths$artifacts_dir, "site_habitat_all_sites.rds")
legacy_model_path <- file.path(
  cfg$paths$repo_root, "DataProcessing", "outputs",
  "fitted_model_clean_20260820_1041.rds"
)
required <- c(site_year_path, site_habitat_path, legacy_model_path, cfg$paths$predictor_stack)
if (any(!file.exists(required))) stop("Missing multiscale input(s).")

predictors <- terra::rast(cfg$paths$predictor_stack)
base_names <- names(predictors)
site_year <- readRDS(site_year_path)
site_folds <- readRDS(site_habitat_path) |>
  dplyr::filter(.data$evidence_source == "cull_site_history") |>
  dplyr::select("site_id", "site_fold", "clean_complete")

unique_sites <- site_year |>
  dplyr::distinct(.data$site_id, .data$longitude, .data$latitude)
message("Extracting aligned 30/60/120 m block means for ", nrow(unique_sites), " sites.")
site_scales <- attach_multiscale_block_means(unique_sites, predictors)
site_year_scales <- site_year |>
  dplyr::select(-dplyr::any_of(base_names)) |>
  dplyr::left_join(site_scales, by = c("site_id", "longitude", "latitude")) |>
  dplyr::left_join(site_folds, by = "site_id") |>
  dplyr::mutate(reef_year_id = paste(.data$reef_name, .data$year, sep = "||")) |>
  dplyr::group_by(.data$site_id) |>
  dplyr::mutate(observation_weight = 1 / dplyr::n_distinct(.data$year)) |>
  dplyr::ungroup()

feature_sets <- multiscale_feature_sets(base_names)
complete_features <- unique(unlist(feature_sets))
complete_features <- setdiff(complete_features, "year")
site_year_scales <- site_year_scales |>
  dplyr::filter(
    .data$clean_complete == 1L,
    stats::complete.cases(dplyr::pick(dplyr::all_of(complete_features)))
  )

legacy_saved <- readRDS(legacy_model_path)
legacy_parsnip <- workflows::extract_fit_parsnip(legacy_saved$fit)
legacy_args <- legacy_parsnip$spec$args
legacy_cfg <- list(
  seed = cfg$model$seed,
  nrounds = as.integer(rlang::eval_tidy(legacy_args$trees)),
  eta = as.numeric(rlang::eval_tidy(legacy_args$learn_rate)),
  max_depth = as.integer(rlang::eval_tidy(legacy_args$tree_depth)),
  min_child_weight = as.numeric(rlang::eval_tidy(legacy_args$min_n)),
  gamma = as.numeric(rlang::eval_tidy(legacy_args$loss_reduction)),
  subsample = 1,
  colsample_bytree = 1,
  lambda = 1,
  alpha = 0,
  nthread = cfg$model$nthread
)

message("Fitting multiscale point-footprint benchmark on ", nrow(site_year_scales), " site-years.")
point_cache <- file.path(cfg$paths$artifacts_dir, "multiscale_point_benchmark_cache.rds")
if (file.exists(point_cache)) {
  point_benchmark <- readRDS(point_cache)
} else {
  point_benchmark <- run_site_year_rank_benchmark(
    site_year_scales, feature_sets, legacy_cfg
  )
  saveRDS(point_benchmark, point_cache)
}

cell_predictions <- list()
cell_models <- list()
cell_audit <- list()
for (resolution in c(30, 60, 120)) {
  cell_data <- aggregate_site_year_to_resolution(site_year_scales, resolution)
  scale_features <- c(paste0(base_names, "_s", resolution, "m"), "year")
  cell_data <- cell_data |>
    dplyr::filter(stats::complete.cases(dplyr::pick(dplyr::all_of(scale_features))))
  benchmark <- run_site_year_rank_benchmark(
    cell_data,
    setNames(list(scale_features), paste0("cell_pooled_", resolution, "m_year")),
    legacy_cfg
  )
  cell_predictions[[as.character(resolution)]] <- benchmark$predictions
  cell_models <- c(cell_models, benchmark$models)
  cell_audit[[as.character(resolution)]] <- tibble::tibble(
    resolution_m = resolution,
    rows = nrow(cell_data),
    cells = dplyr::n_distinct(cell_data$site_id),
    reefs = dplyr::n_distinct(cell_data$reef_name),
    positives = sum(cell_data$high_density),
    prevalence = mean(cell_data$high_density),
    mean_source_sites_per_cell_year = mean(cell_data$n_source_sites),
    maximum_source_sites_per_cell_year = max(cell_data$n_source_sites)
  )
}

predictions <- dplyr::bind_rows(
  point_benchmark$predictions,
  dplyr::bind_rows(cell_predictions)
)
metrics <- suppressWarnings(summarise_benchmark_predictions(predictions))
decomposition <- suppressWarnings(
  prediction_scale_decomposition(predictions, group = "reef_year_id")
)
capture <- within_group_capture_metrics(predictions, group = "reef_year_id")
contrast <- within_group_percentile_contrast(predictions, group = "reef_year_id")
resolution_audit <- dplyr::bind_rows(cell_audit)

importance <- dplyr::bind_rows(lapply(names(point_benchmark$models), function(key) {
  value <- xgboost::xgb.importance(model = point_benchmark$models[[key]])
  if (nrow(value) == 0L) return(NULL)
  tibble::as_tibble(value) |>
    dplyr::mutate(
      model = sub("__fold_[0-9]+$", "", key),
      fold = as.integer(sub("^.*__fold_", "", key))
    )
})) |>
  dplyr::group_by(.data$model, .data$Feature) |>
  dplyr::summarise(mean_gain = mean(.data$Gain), folds = dplyr::n(), .groups = "drop") |>
  dplyr::arrange(.data$model, dplyr::desc(.data$mean_gain))

tables <- list(
  multiscale_resolution_metrics = metrics,
  multiscale_resolution_decomposition = decomposition,
  multiscale_resolution_capture = capture,
  multiscale_resolution_percentile_contrast = contrast,
  multiscale_resolution_cell_audit = resolution_audit,
  multiscale_point_importance = importance
)
for (name in names(tables)) {
  utils::write.csv(tables[[name]], file.path(cfg$paths$outputs_dir, paste0(name, ".csv")), row.names = FALSE)
}
saveRDS(site_scales, file.path(cfg$paths$artifacts_dir, "site_multiscale_30_60_120m.rds"))
saveRDS(predictions, file.path(cfg$paths$outputs_dir, "multiscale_resolution_predictions.rds"))
saveRDS(c(point_benchmark$models, cell_models), file.path(cfg$paths$models_dir, "multiscale_resolution_fold_models.rds"))

message("Multiscale and resolution benchmark complete.")
print(dplyr::left_join(metrics, dplyr::select(decomposition, "model", "within_reef_roc_auc"), by = "model"))
