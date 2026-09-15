args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v2_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v2_dir, "R", "bootstrap.R"))
ctx <- initialise_v2(script_file)
cfg <- ctx$config

site_year_path <- file.path(cfg$paths$artifacts_dir, "site_year_all.rds")
site_habitat_path <- file.path(
  cfg$paths$artifacts_dir, "site_habitat_all_sites.rds"
)
legacy_model_path <- file.path(
  cfg$paths$repo_root,
  "DataProcessing", "outputs", "fitted_model_clean_20260820_1041.rds"
)
required <- c(site_year_path, site_habitat_path, legacy_model_path)
if (any(!file.exists(required))) {
  stop(
    "Missing legacy-signal validation inputs: ",
    paste(required[!file.exists(required)], collapse = "; ")
  )
}

clean_features <- names(terra::rast(cfg$paths$predictor_stack))
site_folds <- readRDS(site_habitat_path) |>
  dplyr::filter(.data$evidence_source == "cull_site_history") |>
  dplyr::select("site_id", "site_fold", "clean_complete")

site_year <- readRDS(site_year_path) |>
  dplyr::left_join(site_folds, by = "site_id") |>
  dplyr::filter(.data$clean_complete == 1L) |>
  dplyr::mutate(
    reef_year_id = paste(.data$reef_name, .data$year, sep = "||"),
    observation_weight = 1
  )

legacy_saved <- readRDS(legacy_model_path)
legacy_fit <- workflows::extract_fit_parsnip(legacy_saved$fit)
legacy_args <- legacy_fit$spec$args
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

feature_sets <- list(
  legacy_fixed_clean_static = clean_features,
  legacy_fixed_clean_year = c(clean_features, "year")
)

message(
  "Validating fixed legacy year settings on ", nrow(site_year),
  " annual pooled site records at ",
  dplyr::n_distinct(site_year$site_id), " held-out sites."
)
benchmark <- run_site_year_rank_benchmark(
  site_year, feature_sets, legacy_cfg
)
metrics <- suppressWarnings(
  summarise_benchmark_predictions(benchmark$predictions)
)
decomposition <- suppressWarnings(
  prediction_scale_decomposition(
    benchmark$predictions, group = "reef_year_id"
  )
)
capture <- within_group_capture_metrics(
  benchmark$predictions, group = "reef_year_id"
)
contrast <- within_group_percentile_contrast(
  benchmark$predictions, group = "reef_year_id"
)
percentile_decomposition <- suppressWarnings(
  prediction_scale_decomposition(
    within_group_percentile_predictions(
      benchmark$predictions, group = "reef_year_id"
    ),
    group = "reef_year_id"
  )
)
rank_invariance <- dplyr::left_join(
  dplyr::select(
    decomposition,
    "model", raw_within_group_auc = "within_reef_roc_auc"
  ),
  dplyr::select(
    percentile_decomposition,
    "model", percentile_within_group_auc = "within_reef_roc_auc"
  ),
  by = "model"
) |>
  dplyr::mutate(
    auc_change_after_percentile_rescale =
      .data$percentile_within_group_auc - .data$raw_within_group_auc
  )

tables <- list(
  legacy_fixed_rank_metrics = metrics,
  legacy_fixed_rank_decomposition = decomposition,
  legacy_fixed_within_reef_year_capture = capture,
  legacy_fixed_percentile_contrast = contrast,
  legacy_fixed_rank_invariance = rank_invariance
)
for (table_name in names(tables)) {
  utils::write.csv(
    tables[[table_name]],
    file.path(cfg$paths$outputs_dir, paste0(table_name, ".csv")),
    row.names = FALSE
  )
}
saveRDS(
  benchmark$predictions,
  file.path(cfg$paths$outputs_dir, "legacy_fixed_rank_predictions.rds")
)
saveRDS(
  benchmark$models,
  file.path(cfg$paths$models_dir, "legacy_fixed_rank_fold_models.rds")
)

message("Fixed legacy year-signal validation complete.")
print(dplyr::left_join(
  metrics,
  dplyr::select(decomposition, "model", "within_reef_roc_auc"),
  by = "model"
))
print(dplyr::filter(capture, abs(.data$target_fraction - 0.20) < 1e-9))
