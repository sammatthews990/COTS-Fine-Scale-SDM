args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v2_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v2_dir, "R", "bootstrap.R"))
ctx <- initialise_v2(script_file)
cfg <- ctx$config

cat("Reading and validating cull observations...\n")
observations <- read_cull_observations(cfg$paths$cull_workbook)
threshold <- cfg$response$cpue_threshold

cat("Building effort-aware site-year targets...\n")
site_year_all <- summarise_site_year(
  observations,
  threshold = threshold,
  first_visit_only = FALSE
)
site_year_first_visit <- summarise_site_year(
  observations,
  threshold = threshold,
  first_visit_only = TRUE
)
site_history <- summarise_site_history_diagnostics(
  site_year_all,
  observations,
  threshold
)

cat("Loading reef-year outbreak probabilities and constructing lag features...\n")
hindcast <- read_outbreak_probabilities(
  cfg$paths$outbreak_hindcast_file, "adjusted_hindcast"
)
forecast <- read_outbreak_probabilities(
  cfg$paths$outbreak_forecast_file, "adjusted_forecast"
)
site_year_all <- add_outbreak_history_features(site_year_all, hindcast)
site_year_first_visit <- add_outbreak_history_features(
  site_year_first_visit, hindcast
)
site_year_all <- attach_outbreak_forecast(site_year_all, forecast)
site_year_first_visit <- attach_outbreak_forecast(
  site_year_first_visit, forecast
)

cat("Extracting static predictors once per immutable site ID...\n")
extraction <- extract_static_site_predictors(
  site_year_all,
  cfg$paths$predictor_stack
)
site_year_all <- extraction$data
site_year_first_visit <- attach_existing_site_predictors(
  site_year_first_visit,
  extraction
)

cat("Assigning whole reefs to deterministic spatial clusters...\n")
folds <- assign_spatial_reef_folds(
  site_year_all,
  k = cfg$validation$spatial_folds,
  seed = cfg$validation$spatial_seed
)
site_year_all <- folds$data
site_year_first_visit <- site_year_first_visit |>
  dplyr::left_join(
    dplyr::select(
      folds$reef_folds,
      dplyr::all_of(c("reef_name", "spatial_fold"))
    ),
    by = "reef_name"
  )

if (any(is.na(site_year_first_visit$spatial_fold))) {
  stop("First-visit data contain reefs without a spatial fold.")
}

coverage <- site_year_all |>
  dplyr::group_by(.data$year) |>
  dplyr::summarise(
    rows = dplyr::n(),
    predictor_complete = sum(.data$predictor_complete == 1L),
    outbreak_lag1_available = sum(!is.na(
      .data$reef_outbreak_probability_lag1
    )),
    outbreak_history_complete = sum(.data$outbreak_history_complete == 1L),
    outbreak_forecast_available = sum(!is.na(
      .data$reef_outbreak_probability_forecast
    )),
    prevalence = mean(.data$high_density),
    .groups = "drop"
  )

fold_balance <- fold_balance_summary(site_year_all)
metadata <- list(
  generated_at = Sys.time(),
  config = cfg,
  predictor_names = extraction$predictor_names,
  raster_crs = extraction$raster_crs,
  raster_resolution = extraction$raster_resolution,
  counts = list(
    raw_observations = nrow(observations),
    sites = dplyr::n_distinct(observations$site_id),
    site_year_all = nrow(site_year_all),
    site_year_first_visit = nrow(site_year_first_visit),
    outbreak_hindcast_reef_years = nrow(hindcast),
    outbreak_forecast_reef_years = nrow(forecast)
  )
)

saveRDS(observations, file.path(cfg$paths$artifacts_dir, "cull_observations_clean.rds"))
saveRDS(
  hindcast,
  file.path(cfg$paths$artifacts_dir, "reef_outbreak_hindcasts.rds")
)
saveRDS(
  forecast,
  file.path(cfg$paths$artifacts_dir, "reef_outbreak_forecast_2026.rds")
)
saveRDS(site_year_all, file.path(cfg$paths$artifacts_dir, "site_year_all.rds"))
saveRDS(
  site_year_first_visit,
  file.path(cfg$paths$artifacts_dir, "site_year_first_visit.rds")
)
saveRDS(site_history, file.path(cfg$paths$artifacts_dir, "site_history_diagnostics.rds"))
saveRDS(folds$reef_folds, file.path(cfg$paths$artifacts_dir, "reef_spatial_folds.rds"))
saveRDS(metadata, file.path(cfg$paths$artifacts_dir, "preparation_metadata.rds"))

utils::write.csv(
  coverage,
  file.path(cfg$paths$outputs_dir, "data_coverage_by_year.csv"),
  row.names = FALSE
)
utils::write.csv(
  fold_balance,
  file.path(cfg$paths$outputs_dir, "spatial_fold_balance.csv"),
  row.names = FALSE
)
utils::write.csv(
  site_history,
  file.path(cfg$paths$outputs_dir, "site_history_diagnostics.csv"),
  row.names = FALSE
)

cat("Preparation complete.\n")
cat(sprintf(
  "Rows: %d all-dives site-years; %d first-visit site-years; %d static predictors.\n",
  nrow(site_year_all), nrow(site_year_first_visit),
  length(extraction$predictor_names)
))
cat(sprintf(
  "Predictor-complete rows: %d. Complete lag-1-to-3 outbreak histories: %d.\n",
  sum(site_year_all$predictor_complete == 1L),
  sum(site_year_all$outbreak_history_complete == 1L)
))
print(fold_balance)
