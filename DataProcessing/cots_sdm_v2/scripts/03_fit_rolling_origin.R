args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v2_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v2_dir, "R", "bootstrap.R"))
ctx <- initialise_v2(script_file)
cfg <- ctx$config

metadata_path <- file.path(cfg$paths$artifacts_dir, "preparation_metadata.rds")
if (!file.exists(metadata_path)) {
  stop("Run scripts/01_prepare_data.R before rolling-origin fitting.")
}
metadata <- readRDS(metadata_path)
target_name <- cfg$response$default_target
data <- readRDS(file.path(
  cfg$paths$artifacts_dir, paste0(target_name, ".rds")
))

feature_sets <- default_feature_sets(metadata$predictor_names)
all_features <- unique(unlist(feature_sets, use.names = FALSE))
model_data <- data |>
  dplyr::filter(
    .data$predictor_complete == 1L,
    .data$outbreak_history_complete == 1L
  ) |>
  dplyr::filter(stats::complete.cases(dplyr::pick(
    dplyr::all_of(c("high_density", all_features))
  )))

test_years <- cfg$validation$rolling_origin_years
missing_years <- setdiff(test_years, unique(model_data$year))
if (length(missing_years) > 0L) {
  stop("Rolling-origin test years missing from prepared data: ",
       paste(missing_years, collapse = ", "))
}

cat(sprintf(
  "Running rolling-origin benchmark for %s on %d common rows.\n",
  paste(test_years, collapse = ", "), nrow(model_data)
))
rolling <- run_rolling_origin_benchmark(
  model_data,
  feature_sets,
  cfg$model,
  test_years
)
predictions <- rolling$predictions

metrics_by_year <- predictions |>
  dplyr::group_by(.data$model, .data$test_year) |>
  dplyr::group_modify(~ binary_probability_metrics(
    .x$truth, .x$probability
  )) |>
  dplyr::ungroup()

metrics_overall <- predictions |>
  dplyr::group_by(.data$model) |>
  dplyr::group_modify(~ binary_probability_metrics(
    .x$truth, .x$probability
  )) |>
  dplyr::ungroup()

metrics_by_reef_history <- predictions |>
  dplyr::mutate(
    reef_history = dplyr::if_else(
      .data$reef_seen_before, "seen_reef", "new_reef"
    )
  ) |>
  dplyr::group_by(.data$model, .data$test_year, .data$reef_history) |>
  dplyr::group_modify(~ binary_probability_metrics(
    .x$truth, .x$probability
  )) |>
  dplyr::ungroup()

saveRDS(
  predictions,
  file.path(cfg$paths$outputs_dir, "rolling_origin_predictions.rds")
)
saveRDS(
  rolling$models,
  file.path(cfg$paths$models_dir, "rolling_origin_models.rds")
)
utils::write.csv(
  predictions,
  file.path(cfg$paths$outputs_dir, "rolling_origin_predictions.csv"),
  row.names = FALSE
)
utils::write.csv(
  metrics_by_year,
  file.path(cfg$paths$outputs_dir, "rolling_origin_metrics_by_year.csv"),
  row.names = FALSE
)
utils::write.csv(
  metrics_overall,
  file.path(cfg$paths$outputs_dir, "rolling_origin_metrics_overall.csv"),
  row.names = FALSE
)
utils::write.csv(
  metrics_by_reef_history,
  file.path(
    cfg$paths$outputs_dir,
    "rolling_origin_metrics_by_reef_history.csv"
  ),
  row.names = FALSE
)

cat("Rolling-origin benchmark complete.\n")
print(metrics_by_year)
