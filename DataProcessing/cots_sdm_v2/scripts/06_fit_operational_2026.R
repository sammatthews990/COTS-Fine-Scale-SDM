args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v2_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v2_dir, "R", "bootstrap.R"))
ctx <- initialise_v2(script_file)
cfg <- ctx$config

metadata_path <- file.path(cfg$paths$artifacts_dir, "preparation_metadata.rds")
if (!file.exists(metadata_path)) {
  stop("Run scripts/01_prepare_data.R before operational fitting.")
}
metadata <- readRDS(metadata_path)
target_name <- cfg$response$default_target
data <- readRDS(file.path(
  cfg$paths$artifacts_dir, paste0(target_name, ".rds")
)) |>
  dplyr::mutate(
    reef_outbreak_probability_operational = dplyr::if_else(
      .data$year == cfg$outbreak$forecast_year,
      .data$reef_outbreak_probability_forecast,
      .data$reef_outbreak_probability_current
    ),
    outbreak_value_source = dplyr::if_else(
      .data$year == cfg$outbreak$forecast_year,
      "adjusted_forecast",
      "adjusted_hindcast"
    )
  )

static <- metadata$predictor_names
feature_sets <- list(
  operational_static = static,
  operational_outbreak_only = "reef_outbreak_probability_operational",
  operational_static_outbreak = c(
    static, "reef_outbreak_probability_operational"
  ),
  operational_static_year_outbreak = c(
    static, "year", "reef_outbreak_probability_operational"
  )
)
all_features <- unique(unlist(feature_sets, use.names = FALSE))
model_data <- data |>
  dplyr::filter(.data$predictor_complete == 1L) |>
  dplyr::filter(stats::complete.cases(dplyr::pick(
    dplyr::all_of(c("high_density", all_features))
  )))

test_year <- cfg$outbreak$forecast_year
cat(sprintf(
  paste0(
    "Running provisional operational-%d benchmark on %d common rows; ",
    "%d forecast rows.\n"
  ),
  test_year, nrow(model_data), sum(model_data$year == test_year)
))

operational <- run_forward_benchmark(
  model_data,
  feature_sets,
  cfg$model,
  test_year
)
test <- model_data[operational$split$test, , drop = FALSE]
raw_forecast <- tibble::tibble(
  model = "raw_reef_forecast",
  validation = paste0("provisional_operational_", test_year),
  fold = NA_integer_,
  site_year_id = test$site_year_id,
  reef_name = test$reef_name,
  site_id = test$site_id,
  year = test$year,
  truth = test$high_density,
  probability = test$reef_outbreak_probability_operational
)

predictions <- operational$predictions |>
  dplyr::mutate(
    validation = paste0("provisional_operational_", test_year)
  ) |>
  dplyr::bind_rows(raw_forecast)
metrics <- summarise_benchmark_predictions(predictions)

manifest <- list(
  generated_at = Sys.time(),
  target = target_name,
  forecast_year = test_year,
  common_rows = nrow(model_data),
  forecast_rows = nrow(test),
  feature_sets = feature_sets,
  split = operational$split,
  hindcast_provenance = cfg$outbreak$hindcast_provenance,
  forecast_provenance = cfg$outbreak$forecast_provenance,
  status = "provisional_pending_provenance_audit"
)

saveRDS(
  predictions,
  file.path(cfg$paths$outputs_dir, "operational_2026_predictions.rds")
)
saveRDS(
  operational$models,
  file.path(cfg$paths$models_dir, "operational_2026_models.rds")
)
saveRDS(
  manifest,
  file.path(cfg$paths$outputs_dir, "operational_2026_manifest.rds")
)
utils::write.csv(
  predictions,
  file.path(cfg$paths$outputs_dir, "operational_2026_predictions.csv"),
  row.names = FALSE
)
utils::write.csv(
  metrics,
  file.path(cfg$paths$outputs_dir, "operational_2026_metrics.csv"),
  row.names = FALSE
)

cat("Provisional operational benchmark complete.\n")
cat("Do not treat these metrics as unbiased until provenance is confirmed.\n")
print(metrics)
