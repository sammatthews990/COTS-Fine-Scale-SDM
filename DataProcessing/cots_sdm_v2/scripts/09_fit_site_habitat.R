args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v2_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v2_dir, "R", "bootstrap.R"))
ctx <- initialise_v2(script_file)
cfg <- ctx$config

site_year_path <- file.path(
  cfg$paths$artifacts_dir, "site_year_first_visit.rds"
)
manta_path <- file.path(
  cfg$paths$artifacts_dir, "manta_observations_clean.rds"
)
required <- c(
  site_year_path,
  manta_path,
  cfg$paths$predictor_stack,
  cfg$paths$predictor_stack_rg,
  cfg$paths$geomorphic_raster
)
if (any(!file.exists(required))) {
  stop(
    "Missing site-habitat inputs: ",
    paste(required[!file.exists(required)], collapse = "; ")
  )
}

message("Preparing lifetime cull-site response from first-visit site-years.")
site_year <- readRDS(site_year_path)
cull_sites <- summarise_cull_site_history(site_year)
template <- terra::rast(cfg$paths$predictor_stack)[[1]]

message("Aggregating conservative manta zeros to site-scale cells.")
manta <- readRDS(manta_path)
manta_summary <- summarise_manta_sites(
  manta,
  template,
  grid_m = cfg$site_habitat$manta_site_grid_m,
  minimum_tows = cfg$site_habitat$manta_minimum_tows,
  minimum_years = cfg$site_habitat$manta_minimum_years
)
conflict_check <- exclude_manta_local_conflicts(
  manta_summary,
  cull_sites,
  template,
  distance_m =
    cfg$site_habitat$manta_local_conflict_distance_m
)
manta_sites <- finalise_manta_zero_sites(
  conflict_check$retained,
  template,
  cfg$site_habitat$manta_site_grid_m
)

message(
  "Cull sites: ", nrow(cull_sites),
  "; manta candidates: ", nrow(manta_summary$candidates),
  "; locally uncontradicted manta zeros: ", nrow(manta_sites)
)

site_columns <- c(
  "site_id", "reef_name", "site_name", "longitude", "latitude",
  "high_density", "n_surveyed_years", "n_high_density_years",
  "high_density_year_fraction", "first_year", "last_year",
  "total_bottom_time", "total_cots", "evidence_source",
  "observation_weight"
)
all_sites <- dplyr::bind_rows(
  dplyr::select(cull_sites, dplyr::all_of(site_columns)),
  dplyr::select(manta_sites, dplyr::all_of(site_columns))
)

message("Extracting clean, RG, and categorical geomorphic predictors.")
all_sites <- extract_site_habitat_predictors(
  all_sites,
  cfg$paths$predictor_stack,
  cfg$paths$predictor_stack_rg,
  cfg$paths$geomorphic_raster
)
all_sites <- add_site_validation_group(
  all_sites,
  template,
  cfg$site_habitat$manta_site_grid_m
) |>
  assign_within_reef_site_folds(
    k = cfg$site_habitat$within_reef_folds,
    seed = cfg$validation$spatial_seed
  )

clean_names <- names(terra::rast(cfg$paths$predictor_stack))
feature_sets <- site_habitat_feature_sets(clean_names)
common_sites <- all_sites |>
  dplyr::filter(
    .data$clean_complete == 1L,
    .data$rg_complete == 1L
  )
cull_model_sites <- common_sites |>
  dplyr::filter(.data$evidence_source == "cull_site_history")
manta_model_sites <- common_sites |>
  dplyr::filter(
    .data$evidence_source == "manta_repeated_site_zero"
  )

message(
  "Common comparison domain: ", nrow(cull_model_sites),
  " cull sites and ", nrow(manta_model_sites), " manta-zero sites."
)
message("Running controlled within-reef site cross-validation.")
benchmark <- run_site_habitat_benchmark(
  cull_model_sites,
  manta_model_sites,
  feature_sets,
  cfg$model,
  maximum_manta_weight_ratio =
    cfg$site_habitat$manta_maximum_negative_weight_ratio
)

predictions <- benchmark$predictions
metrics <- suppressWarnings(
  summarise_benchmark_predictions(predictions)
)
decomposition <- suppressWarnings(
  prediction_scale_decomposition(predictions)
)
importance <- site_habitat_importance(
  benchmark$models, feature_sets
)

effort_metrics <- predictions |>
  dplyr::mutate(
    effort_subset = dplyr::if_else(
      .data$truth == 1L | .data$n_surveyed_years >= 2L,
      "positive_or_repeated_negative",
      "single_year_negative_excluded"
    )
  ) |>
  dplyr::filter(
    .data$effort_subset == "positive_or_repeated_negative"
  ) |>
  dplyr::group_by(.data$model, .data$effort_subset) |>
  dplyr::group_modify(~ binary_probability_metrics(
    .x$truth, .x$probability
  )) |>
  dplyr::ungroup()

manta_audit <- tibble::tibble(
  site_grid_m = cfg$site_habitat$manta_site_grid_m,
  minimum_tows = cfg$site_habitat$manta_minimum_tows,
  minimum_years = cfg$site_habitat$manta_minimum_years,
  local_conflict_distance_m =
    cfg$site_habitat$manta_local_conflict_distance_m,
  candidate_sites = nrow(manta_summary$candidates),
  positive_cull_conflicts = sum(
    conflict_check$audit$conflict_positive_cull
  ),
  positive_manta_conflicts = sum(
    conflict_check$audit$conflict_positive_manta
  ),
  retained_sites_before_predictor_coverage = nrow(manta_sites),
  retained_common_domain_sites = nrow(manta_model_sites),
  retained_common_domain_reefs =
    dplyr::n_distinct(manta_model_sites$reef_name)
)
domain_summary <- common_sites |>
  dplyr::group_by(.data$evidence_source) |>
  dplyr::summarise(
    sites = dplyr::n(),
    reefs = dplyr::n_distinct(.data$reef_name),
    positives = sum(.data$high_density),
    prevalence = mean(.data$high_density),
    geomorphic_unknown = sum(.data$geomorphic_class == "Unknown"),
    .groups = "drop"
  )

saveRDS(
  all_sites,
  file.path(cfg$paths$artifacts_dir, "site_habitat_all_sites.rds")
)
saveRDS(
  conflict_check$audit,
  file.path(
    cfg$paths$artifacts_dir,
    "manta_site_zero_conflict_audit.rds"
  )
)
saveRDS(
  predictions,
  file.path(
    cfg$paths$outputs_dir, "site_habitat_predictions.rds"
  )
)
saveRDS(
  benchmark$models,
  file.path(
    cfg$paths$models_dir, "site_habitat_fold_models.rds"
  )
)

output_tables <- list(
  site_habitat_metrics = metrics,
  site_habitat_scale_decomposition = decomposition,
  site_habitat_effort_sensitivity = effort_metrics,
  site_habitat_training_audit = benchmark$training_audit,
  site_habitat_feature_importance = importance,
  site_habitat_manta_zero_audit = manta_audit,
  site_habitat_domain_summary = domain_summary,
  site_habitat_predictions = predictions
)
for (table_name in names(output_tables)) {
  utils::write.csv(
    output_tables[[table_name]],
    file.path(
      cfg$paths$outputs_dir, paste0(table_name, ".csv")
    ),
    row.names = FALSE
  )
}

message("Site-habitat benchmark complete.")
print(manta_audit)
print(domain_summary)
print(dplyr::left_join(
  dplyr::select(
    metrics,
    dplyr::all_of(c(
      "model", "roc_auc", "pr_auc", "brier", "capture_top_20pct"
    ))
  ),
  dplyr::select(
    decomposition,
    dplyr::all_of(c(
      "model", "between_reef_only_roc_auc",
      "within_reef_roc_auc"
    ))
  ),
  by = "model"
))
