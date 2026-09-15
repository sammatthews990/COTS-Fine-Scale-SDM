args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v2_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v2_dir, "R", "bootstrap.R"))
ctx <- initialise_v2(script_file)
cfg <- ctx$config

site_year <- readRDS(file.path(cfg$paths$artifacts_dir, "site_year_all.rds"))
site_scales <- readRDS(file.path(cfg$paths$artifacts_dir, "site_multiscale_30_60_120m.rds"))
site_folds <- readRDS(file.path(cfg$paths$artifacts_dir, "site_habitat_all_sites.rds")) |>
  dplyr::filter(.data$evidence_source == "cull_site_history") |>
  dplyr::select("site_id", "site_fold", "clean_complete")
base_names <- names(terra::rast(cfg$paths$predictor_stack))

data <- site_year |>
  dplyr::select(-dplyr::any_of(base_names)) |>
  dplyr::left_join(site_scales, by = c("site_id", "longitude", "latitude")) |>
  dplyr::left_join(site_folds, by = "site_id") |>
  dplyr::mutate(
    reef_year_id = paste(.data$reef_name, .data$year, sep = "||"),
    observation_weight = 1
  )

all_sets <- multiscale_feature_sets(base_names)
feature_sets <- all_sets[c(
  "clean_30m_year", "clean_60m_year", "clean_120m_year",
  "multiscale_topography_year", "multiscale_topography_static"
)]
required_features <- unique(setdiff(unlist(feature_sets), "year"))
data <- data |>
  dplyr::filter(
    .data$clean_complete == 1L,
    stats::complete.cases(dplyr::pick(dplyr::all_of(required_features)))
  )

legacy <- readRDS(file.path(
  cfg$paths$repo_root, "DataProcessing", "outputs",
  "fitted_model_clean_20260820_1041.rds"
))
parsnip <- workflows::extract_fit_parsnip(legacy$fit)
model_args <- parsnip$spec$args
legacy_cfg <- list(
  seed = cfg$model$seed,
  nrounds = as.integer(rlang::eval_tidy(model_args$trees)),
  eta = as.numeric(rlang::eval_tidy(model_args$learn_rate)),
  max_depth = as.integer(rlang::eval_tidy(model_args$tree_depth)),
  min_child_weight = as.numeric(rlang::eval_tidy(model_args$min_n)),
  gamma = as.numeric(rlang::eval_tidy(model_args$loss_reduction)),
  subsample = 1, colsample_bytree = 1, lambda = 1, alpha = 0,
  nthread = cfg$model$nthread
)

benchmark <- run_site_year_rank_benchmark(data, feature_sets, legacy_cfg)
predictions <- benchmark$predictions |>
  dplyr::mutate(model = paste0(.data$model, "__legacy_unweighted"))
metrics <- suppressWarnings(summarise_benchmark_predictions(predictions))
decomposition <- suppressWarnings(
  prediction_scale_decomposition(predictions, group = "reef_year_id")
)
capture <- within_group_capture_metrics(predictions, group = "reef_year_id")

utils::write.csv(metrics, file.path(cfg$paths$outputs_dir, "legacy_unweighted_scale_metrics.csv"), row.names = FALSE)
utils::write.csv(decomposition, file.path(cfg$paths$outputs_dir, "legacy_unweighted_scale_decomposition.csv"), row.names = FALSE)
utils::write.csv(capture, file.path(cfg$paths$outputs_dir, "legacy_unweighted_scale_capture.csv"), row.names = FALSE)
saveRDS(predictions, file.path(cfg$paths$outputs_dir, "legacy_unweighted_scale_predictions.rds"))
saveRDS(benchmark$models, file.path(cfg$paths$models_dir, "legacy_unweighted_scale_fold_models.rds"))

message("Legacy-unweighted scale sensitivity complete.")
print(dplyr::left_join(metrics, dplyr::select(decomposition, "model", "within_reef_roc_auc"), by = "model"))
