args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v2_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v2_dir, "R", "bootstrap.R"))
ctx <- initialise_v2(script_file)
cfg <- ctx$config

rolling_path <- file.path(
  cfg$paths$outputs_dir, "rolling_origin_predictions.rds"
)
operational_path <- file.path(
  cfg$paths$outputs_dir, "operational_2026_predictions.rds"
)
site_year_path <- file.path(
  cfg$paths$artifacts_dir, "site_year_first_visit.rds"
)
required <- c(rolling_path, operational_path, site_year_path)
if (any(!file.exists(required))) {
  stop(
    "Run scripts 01, 03, and 06 first. Missing: ",
    paste(required[!file.exists(required)], collapse = "; ")
  )
}

rolling <- readRDS(rolling_path)
site_year <- readRDS(site_year_path)
training_years <- sort(setdiff(
  unique(rolling$test_year),
  cfg$outbreak$forecast_year
))

bridge_training <- rolling |>
  dplyr::filter(
    .data$model == "static",
    .data$test_year %in% training_years
  ) |>
  dplyr::select(
    dplyr::all_of(c("site_year_id", "test_year", "truth")),
    static_probability = "probability"
  ) |>
  dplyr::left_join(
    dplyr::select(
      site_year,
      dplyr::all_of(c(
        "site_year_id", "reef_outbreak_probability_current"
      ))
    ),
    by = "site_year_id"
  ) |>
  dplyr::filter(stats::complete.cases(dplyr::pick(
    dplyr::all_of(c(
      "truth", "static_probability",
      "reef_outbreak_probability_current"
    ))
  )))

bridge <- fit_additive_operational_bridge(
  bridge_training$truth,
  bridge_training$static_probability,
  bridge_training$reef_outbreak_probability_current
)
bridge$training_years <- training_years
bridge$static_score_source <- "rolling_origin_out_of_time_predictions"
bridge$outbreak_source <- "adjusted_same_year_hindcasts"
bridge$status <- "provisional_pending_outbreak_provenance_audit"

operational <- readRDS(operational_path)
wide <- operational |>
  dplyr::filter(.data$model %in% c(
    "operational_static", "raw_reef_forecast"
  )) |>
  dplyr::select(
    dplyr::all_of(c(
      "model", "validation", "site_year_id", "reef_name",
      "site_id", "year", "truth", "probability"
    ))
  ) |>
  tidyr::pivot_wider(
    names_from = "model",
    values_from = "probability"
  )

additive <- wide |>
  dplyr::transmute(
    model = "operational_additive_bridge",
    validation = .data$validation,
    fold = NA_integer_,
    site_year_id = .data$site_year_id,
    reef_name = .data$reef_name,
    site_id = .data$site_id,
    year = .data$year,
    truth = .data$truth,
    probability = predict_additive_operational_bridge(
      bridge,
      .data$operational_static,
      .data$raw_reef_forecast
    )
  )

comparison <- dplyr::bind_rows(operational, additive)
metrics <- suppressWarnings(summarise_benchmark_predictions(comparison))
decomposition <- suppressWarnings(
  prediction_scale_decomposition(comparison)
)
coefficients <- tibble::tibble(
  term = names(bridge$coefficients),
  estimate = unname(bridge$coefficients),
  training_rows = bridge$training_rows,
  training_year_min = min(training_years),
  training_year_max = max(training_years),
  rank_preserving_within_reef =
    bridge$rank_preserving_within_reef
)

saveRDS(
  bridge,
  file.path(cfg$paths$models_dir, "additive_operational_bridge.rds")
)
saveRDS(
  comparison,
  file.path(
    cfg$paths$outputs_dir,
    "operational_2026_decomposed_predictions.rds"
  )
)
utils::write.csv(
  comparison,
  file.path(
    cfg$paths$outputs_dir,
    "operational_2026_decomposed_predictions.csv"
  ),
  row.names = FALSE
)
utils::write.csv(
  metrics,
  file.path(
    cfg$paths$outputs_dir,
    "operational_2026_decomposed_metrics.csv"
  ),
  row.names = FALSE
)
utils::write.csv(
  decomposition,
  file.path(
    cfg$paths$outputs_dir,
    "operational_2026_scale_decomposition.csv"
  ),
  row.names = FALSE
)
utils::write.csv(
  coefficients,
  file.path(
    cfg$paths$outputs_dir,
    "additive_operational_bridge_coefficients.csv"
  ),
  row.names = FALSE
)

message("Additive rank-preserving bridge fitted.")
message("Training years: ", paste(training_years, collapse = ", "))
print(coefficients)
print(dplyr::filter(
  metrics,
  .data$model %in% c(
    "operational_static",
    "operational_static_outbreak",
    "operational_additive_bridge"
  )
))
print(dplyr::filter(
  decomposition,
  .data$model %in% c(
    "operational_static",
    "operational_static_outbreak",
    "operational_additive_bridge",
    "raw_reef_forecast"
  )
))
