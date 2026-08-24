args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v2_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v2_dir, "R", "bootstrap.R"))
ctx <- initialise_v2(script_file)
cfg <- ctx$config

metadata_path <- file.path(cfg$paths$artifacts_dir, "preparation_metadata.rds")
if (!file.exists(metadata_path)) {
  stop("Run scripts/01_prepare_data.R before fitting the benchmark.")
}
metadata <- readRDS(metadata_path)
target_name <- cfg$response$default_target
target_path <- file.path(cfg$paths$artifacts_dir, paste0(target_name, ".rds"))
if (!file.exists(target_path)) stop("Prepared target not found: ", target_path)

data <- readRDS(target_path)
feature_sets <- default_feature_sets(metadata$predictor_names)
all_features <- unique(unlist(feature_sets, use.names = FALSE))
assert_required_columns(
  data,
  c(
    "site_year_id", "site_id", "reef_name", "year", "high_density",
    "spatial_fold", all_features
  ),
  "Prepared benchmark data"
)

# Use identical rows for every feature ablation. Current-year outbreak
# probability is intentionally excluded; only prior-year values enter.
model_data <- data |>
  dplyr::filter(
    .data$predictor_complete == 1L,
    .data$outbreak_history_complete == 1L
  ) |>
  dplyr::filter(stats::complete.cases(dplyr::pick(
    dplyr::all_of(c("high_density", all_features))
  )))

test_year <- cfg$validation$forward_test_year
spatial_data <- model_data |>
  dplyr::filter(.data$year < test_year)
if (sum(model_data$year == test_year) == 0L) {
  stop("No rows are available for forward test year ", test_year)
}

cat(sprintf(
  "Benchmark target: %s; common rows: %d; forward-%d rows: %d.\n",
  target_name,
  nrow(model_data),
  test_year,
  sum(model_data$year == test_year)
))

cat("Running unseen-reef-cluster spatial benchmark...\n")
spatial <- run_spatial_benchmark(
  spatial_data,
  feature_sets,
  cfg$model
)

cat("Running genuine forward-year benchmark...\n")
forward <- run_forward_benchmark(
  model_data,
  feature_sets,
  cfg$model,
  test_year
)

predictions <- dplyr::bind_rows(
  spatial$predictions,
  forward$predictions
)
metrics <- summarise_benchmark_predictions(predictions)
metrics_by_fold <- predictions |>
  dplyr::filter(.data$validation == "spatial_unseen_reef_cluster") |>
  dplyr::group_by(.data$model, .data$validation, .data$fold) |>
  dplyr::group_modify(~ binary_probability_metrics(
    .x$truth, .x$probability
  )) |>
  dplyr::ungroup()

saveRDS(predictions, file.path(cfg$paths$outputs_dir, "benchmark_predictions.rds"))
saveRDS(
  list(spatial = spatial$models, forward = forward$models),
  file.path(cfg$paths$models_dir, "benchmark_models.rds")
)
saveRDS(
  list(
    target = target_name,
    common_rows = nrow(model_data),
    feature_sets = feature_sets,
    forward_split = forward$split,
    config = cfg
  ),
  file.path(cfg$paths$outputs_dir, "benchmark_manifest.rds")
)
utils::write.csv(
  predictions,
  file.path(cfg$paths$outputs_dir, "benchmark_predictions.csv"),
  row.names = FALSE
)
utils::write.csv(
  metrics,
  file.path(cfg$paths$outputs_dir, "benchmark_metrics.csv"),
  row.names = FALSE
)
utils::write.csv(
  metrics_by_fold,
  file.path(cfg$paths$outputs_dir, "benchmark_metrics_by_fold.csv"),
  row.names = FALSE
)

cat("Benchmark complete.\n")
print(metrics)

