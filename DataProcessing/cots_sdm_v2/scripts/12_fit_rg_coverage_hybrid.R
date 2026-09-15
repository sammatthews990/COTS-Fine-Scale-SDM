args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v2_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v2_dir, "R", "bootstrap.R"))
ctx <- initialise_v2(script_file)
cfg <- ctx$config

sites_path <- file.path(
  cfg$paths$artifacts_dir, "site_habitat_all_sites.rds"
)
if (!file.exists(sites_path)) {
  stop("Run scripts/09_fit_site_habitat.R first.")
}

all_sites <- readRDS(sites_path)
clean_sites <- all_sites |>
  dplyr::filter(
    .data$evidence_source == "cull_site_history",
    .data$clean_complete == 1L
  )
rg_sites <- clean_sites |>
  dplyr::filter(.data$rg_complete == 1L)
clean_features <- names(terra::rast(cfg$paths$predictor_stack))
rg_features <- c(
  clean_features,
  "RG_bathy", "RG_slope", "RG_turbid",
  "RG_waves_Hs", "RG_waves_Tp"
)

message(
  "Running full-clean-domain coverage benchmark: ", nrow(clean_sites),
  " clean sites; ", nrow(rg_sites), " with RG coverage."
)
benchmark <- run_rg_coverage_hybrid_benchmark(
  clean_sites,
  rg_sites,
  clean_features,
  rg_features,
  cfg$model
)
metrics <- suppressWarnings(
  summarise_benchmark_predictions(benchmark$predictions)
)
decomposition <- suppressWarnings(
  prediction_scale_decomposition(benchmark$predictions)
)
coverage <- benchmark$predictions |>
  dplyr::filter(.data$model == "rg_with_clean_fallback") |>
  dplyr::group_by(.data$prediction_source) |>
  dplyr::summarise(
    sites = dplyr::n(),
    reefs = dplyr::n_distinct(.data$reef_name),
    positives = sum(.data$truth),
    prevalence = mean(.data$truth),
    .groups = "drop"
  )

utils::write.csv(
  metrics,
  file.path(
    cfg$paths$outputs_dir, "rg_coverage_hybrid_metrics.csv"
  ),
  row.names = FALSE
)
utils::write.csv(
  decomposition,
  file.path(
    cfg$paths$outputs_dir, "rg_coverage_hybrid_decomposition.csv"
  ),
  row.names = FALSE
)
utils::write.csv(
  coverage,
  file.path(
    cfg$paths$outputs_dir, "rg_coverage_hybrid_coverage.csv"
  ),
  row.names = FALSE
)
saveRDS(
  benchmark$predictions,
  file.path(
    cfg$paths$outputs_dir, "rg_coverage_hybrid_predictions.rds"
  )
)
saveRDS(
  benchmark$models,
  file.path(
    cfg$paths$models_dir, "rg_coverage_hybrid_fold_models.rds"
  )
)

print(dplyr::left_join(metrics, decomposition, by = "model"))
print(coverage)
