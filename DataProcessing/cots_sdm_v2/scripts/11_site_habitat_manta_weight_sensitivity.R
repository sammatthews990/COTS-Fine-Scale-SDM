args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v2_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v2_dir, "R", "bootstrap.R"))
ctx <- initialise_v2(script_file)
cfg <- ctx$config

sites_path <- file.path(
  cfg$paths$artifacts_dir, "site_habitat_all_sites.rds"
)
baseline_predictions_path <- file.path(
  cfg$paths$outputs_dir, "site_habitat_predictions.rds"
)
required <- c(sites_path, baseline_predictions_path)
if (any(!file.exists(required))) {
  stop("Run scripts/09_fit_site_habitat.R first.")
}

sites <- readRDS(sites_path) |>
  dplyr::filter(
    .data$clean_complete == 1L,
    .data$rg_complete == 1L
  )
cull_sites <- sites |>
  dplyr::filter(.data$evidence_source == "cull_site_history")
manta_sites <- sites |>
  dplyr::filter(
    .data$evidence_source == "manta_repeated_site_zero"
  )
clean_names <- names(terra::rast(cfg$paths$predictor_stack))
feature_sets <- site_habitat_feature_sets(clean_names)

baseline <- readRDS(baseline_predictions_path) |>
  dplyr::filter(.data$manta_augmented) |>
  dplyr::mutate(manta_weight_ratio = 1.0)

predictions <- list(baseline)
training_audit <- list()
ratios <- cfg$site_habitat$manta_weight_sensitivity_ratios
for (i in seq_along(ratios)) {
  ratio <- ratios[[i]]
  message("Fitting manta-weight sensitivity ratio ", ratio, ".")
  benchmark <- run_site_habitat_benchmark(
    cull_sites,
    manta_sites,
    feature_sets,
    cfg$model,
    maximum_manta_weight_ratio = ratio,
    augmentation_modes = TRUE
  )
  predictions[[i + 1L]] <- benchmark$predictions |>
    dplyr::mutate(manta_weight_ratio = ratio)
  training_audit[[i]] <- benchmark$training_audit |>
    dplyr::mutate(manta_weight_ratio = ratio)
}

predictions <- dplyr::bind_rows(predictions) |>
  dplyr::mutate(
    base_model = .data$model,
    model = paste0(
      .data$model,
      "__manta_ratio_",
      format(.data$manta_weight_ratio, trim = TRUE)
    )
  )
metrics <- suppressWarnings(summarise_benchmark_predictions(predictions))
decomposition <- suppressWarnings(
  prediction_scale_decomposition(predictions)
)

model_lookup <- predictions |>
  dplyr::distinct(
    .data$model, .data$base_model, .data$feature_family,
    .data$manta_weight_ratio
  )
summary <- metrics |>
  dplyr::left_join(decomposition, by = "model") |>
  dplyr::left_join(model_lookup, by = "model") |>
  dplyr::arrange(.data$feature_family, .data$manta_weight_ratio)

utils::write.csv(
  summary,
  file.path(
    cfg$paths$outputs_dir,
    "site_habitat_manta_weight_sensitivity.csv"
  ),
  row.names = FALSE
)
utils::write.csv(
  dplyr::bind_rows(training_audit),
  file.path(
    cfg$paths$outputs_dir,
    "site_habitat_manta_weight_sensitivity_training.csv"
  ),
  row.names = FALSE
)
saveRDS(
  predictions,
  file.path(
    cfg$paths$outputs_dir,
    "site_habitat_manta_weight_sensitivity_predictions.rds"
  )
)

message("Manta-weight sensitivity complete.")
print(dplyr::select(
  summary,
  dplyr::all_of(c(
    "base_model", "manta_weight_ratio", "roc_auc", "pr_auc",
    "brier", "between_reef_only_roc_auc", "within_reef_roc_auc"
  ))
))
