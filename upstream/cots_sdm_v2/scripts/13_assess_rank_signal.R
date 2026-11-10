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
operational_path <- file.path(
  cfg$paths$outputs_dir, "operational_2026_predictions.rds"
)
required <- c(
  site_year_path, site_habitat_path, operational_path,
  cfg$paths$legacy_2026_raster,
  cfg$paths$legacy_agnostic_raster
)
if (any(!file.exists(required))) {
  stop("Missing rank-assessment inputs: ",
       paste(required[!file.exists(required)], collapse = "; "))
}

clean_features <- names(terra::rast(cfg$paths$predictor_stack))
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

rg_features <- c(
  clean_features,
  "RG_bathy", "RG_slope", "RG_turbid",
  "RG_waves_Hs", "RG_waves_Tp"
)
feature_sets <- list(
  annual_clean_static = clean_features,
  annual_clean_year = c(clean_features, "year"),
  annual_rg_static = rg_features,
  annual_rg_year = c(rg_features, "year")
)

message(
  "Fitting annual-response rank benchmark on ", nrow(site_year),
  " site-years at ", dplyr::n_distinct(site_year$site_id), " held-out sites."
)
annual <- run_site_year_rank_benchmark(
  site_year,
  feature_sets,
  cfg$model
)
annual_metrics <- suppressWarnings(
  summarise_benchmark_predictions(annual$predictions)
)
annual_decomposition <- suppressWarnings(
  prediction_scale_decomposition(
    annual$predictions,
    group = "reef_year_id"
  )
)
annual_capture <- within_group_capture_metrics(
  annual$predictions,
  group = "reef_year_id"
)
annual_resolution <- prediction_scale_summary(annual$predictions)
annual_contrast <- within_group_percentile_contrast(
  annual$predictions,
  group = "reef_year_id"
)
annual_percentile_decomposition <- suppressWarnings(
  prediction_scale_decomposition(
    within_group_percentile_predictions(
      annual$predictions,
      group = "reef_year_id"
    ),
    group = "reef_year_id"
  )
)
annual_rank_invariance <- dplyr::left_join(
  dplyr::select(
    annual_decomposition,
    "model",
    raw_within_group_auc = "within_reef_roc_auc"
  ),
  dplyr::select(
    annual_percentile_decomposition,
    "model",
    percentile_within_group_auc = "within_reef_roc_auc"
  ),
  by = "model"
) |>
  dplyr::mutate(
    auc_change_after_percentile_rescale =
      .data$percentile_within_group_auc -
      .data$raw_within_group_auc
  )

operational <- readRDS(operational_path)
operational_static <- operational |>
  dplyr::filter(.data$model == "operational_static")
test_sites <- operational_static |>
  dplyr::select(
    "site_year_id", "reef_name", "site_id", "year", "truth"
  ) |>
  dplyr::left_join(
    dplyr::select(
      readRDS(site_year_path),
      "site_year_id", "longitude", "latitude"
    ),
    by = "site_year_id"
  )

extract_legacy <- function(path, sites) {
  raster <- terra::rast(path)
  points <- terra::vect(
    sites,
    geom = c("longitude", "latitude"),
    crs = "EPSG:4326"
  ) |>
    terra::project(terra::crs(raster))
  as.numeric(terra::extract(raster, points, ID = FALSE)[[1]])
}
legacy_predictions <- dplyr::bind_rows(
  dplyr::mutate(
    test_sites,
    model = "legacy_year_raster_descriptive",
    validation = "descriptive_2026_raster_extraction",
    probability = extract_legacy(
      cfg$paths$legacy_2026_raster, test_sites
    )
  ),
  dplyr::mutate(
    test_sites,
    model = "legacy_agnostic_raster_descriptive",
    validation = "descriptive_2026_raster_extraction",
    probability = extract_legacy(
      cfg$paths$legacy_agnostic_raster, test_sites
    )
  )
) |>
  dplyr::filter(is.finite(.data$probability))

comparison_2026 <- dplyr::bind_rows(
  operational |>
    dplyr::filter(.data$model %in% c(
      "operational_static",
      "operational_static_outbreak",
      "operational_static_year_outbreak"
    )),
  legacy_predictions
)
comparison_2026_metrics <- suppressWarnings(
  summarise_benchmark_predictions(comparison_2026)
)
comparison_2026_decomposition <- suppressWarnings(
  prediction_scale_decomposition(comparison_2026)
)
comparison_2026_capture <- within_group_capture_metrics(
  comparison_2026,
  group = "reef_name"
)
comparison_2026_resolution <- prediction_scale_summary(comparison_2026)
comparison_2026_contrast <- within_group_percentile_contrast(
  comparison_2026,
  group = "reef_name"
)
comparison_2026_percentile_decomposition <- suppressWarnings(
  prediction_scale_decomposition(
    within_group_percentile_predictions(
      comparison_2026,
      group = "reef_name"
    ),
    group = "reef_name"
  )
)
comparison_2026_rank_invariance <- dplyr::left_join(
  dplyr::select(
    comparison_2026_decomposition,
    "model",
    raw_within_group_auc = "within_reef_roc_auc"
  ),
  dplyr::select(
    comparison_2026_percentile_decomposition,
    "model",
    percentile_within_group_auc = "within_reef_roc_auc"
  ),
  by = "model"
) |>
  dplyr::mutate(
    auc_change_after_percentile_rescale =
      .data$percentile_within_group_auc -
      .data$raw_within_group_auc
  )

tables <- list(
  annual_rank_metrics = annual_metrics,
  annual_rank_decomposition = annual_decomposition,
  annual_within_reef_year_capture = annual_capture,
  annual_within_reef_year_percentile_contrast = annual_contrast,
  annual_rank_invariance = annual_rank_invariance,
  annual_signal_resolution_diagnostic = annual_resolution,
  comparison_2026_rank_metrics = comparison_2026_metrics,
  comparison_2026_rank_decomposition = comparison_2026_decomposition,
  comparison_2026_within_reef_capture = comparison_2026_capture,
  comparison_2026_percentile_contrast = comparison_2026_contrast,
  comparison_2026_rank_invariance = comparison_2026_rank_invariance,
  comparison_2026_signal_resolution_diagnostic =
    comparison_2026_resolution
)
for (table_name in names(tables)) {
  utils::write.csv(
    tables[[table_name]],
    file.path(cfg$paths$outputs_dir, paste0(table_name, ".csv")),
    row.names = FALSE
  )
}
saveRDS(
  annual$predictions,
  file.path(cfg$paths$outputs_dir, "annual_rank_predictions.rds")
)
saveRDS(
  annual$models,
  file.path(cfg$paths$models_dir, "annual_rank_fold_models.rds")
)
saveRDS(
  comparison_2026,
  file.path(cfg$paths$outputs_dir, "comparison_2026_rank_predictions.rds")
)

message("Continuous rank-signal assessment complete.")
print(dplyr::left_join(
  annual_metrics,
  dplyr::select(
    annual_decomposition,
    "model", "within_reef_roc_auc"
  ),
  by = "model"
))
print(dplyr::filter(
  annual_capture,
  abs(.data$target_fraction - 0.20) < 1e-9
))
