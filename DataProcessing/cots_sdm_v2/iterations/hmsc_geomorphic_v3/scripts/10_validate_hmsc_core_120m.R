args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v3_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v3_dir, "R", "bootstrap.R"))
ctx <- initialise_v3(script_file)
cfg <- ctx$config

model <- readRDS(file.path(cfg$paths$models_dir, "hmsc_core_120m.rds"))
artifact <- readRDS(file.path(cfg$paths$artifacts_dir, "hmsc_core_120m_fit_data.rds"))
data <- artifact$data
target <- which(!is.na(data$cull_high_density))

make_reef_folds <- function(data, k = 5L) {
  reefs <- data |>
    dplyr::group_by(.data$reef_key) |>
    dplyr::summarise(
      x = mean(.data$projected_x), y = mean(.data$projected_y), .groups = "drop"
    )
  set.seed(20260908L)
  reefs$fold <- stats::kmeans(scale(reefs[, c("x", "y")]), k, nstart = 50)$cluster
  reefs$fold[match(data$reef_key, reefs$reef_key)]
}

validation_model <- model
validation_model$samples <- 250L
validation_model$thin <- 5L
validation_model$transient <- 750L

run_partition <- function(partition, validation, assess_fold = NULL) {
  message("Hmsc cross-validation: ", validation)
  draws <- Hmsc::computePredictedValues(
    validation_model, partition = partition, expected = TRUE,
    nChains = 4L, nParallel = 4L, verbose = 0
  )
  response_index <- match("cull_high_density", artifact$responses)
  probability <- apply(draws, c(1, 2), mean)[, response_index]
  keep <- target
  if (!is.null(assess_fold)) keep <- keep[partition[keep] %in% assess_fold]
  tibble::tibble(
    model = "hmsc_geomorphic_core_120m", validation = validation,
    fold = partition[keep], cell_id = data$cell_id[keep],
    reef_key = data$reef_key[keep],
    reef_year_id = paste(data$reef_key[keep], data$year[keep], sep = "||"),
    year = data$year[keep], truth = data$cull_high_density[keep],
    probability = probability[keep]
  )
}

held_partition <- 1L + abs(
  floor(data$projected_x / 3000) + 3L * floor(data$projected_y / 3000)
) %% 5L
leave_partition <- make_reef_folds(data, 5L)
forward_partition <- ifelse(data$year <= 2023, 1L, 2L)

predictions <- dplyr::bind_rows(
  run_partition(held_partition, "held_cell_5fold"),
  run_partition(leave_partition, "leave_reef_5fold"),
  run_partition(forward_partition, "forward_train_through_2023", assess_fold = 2L)
)
saveRDS(predictions, file.path(
  cfg$paths$outputs_dir, "hmsc_core_120m_blocked_predictions.rds"
))
metrics <- predictions |>
  dplyr::group_by(.data$model, .data$validation) |>
  dplyr::group_modify(~ {
    overall <- binary_probability_metrics(.x$truth, .x$probability)
    within <- within_group_auc(.x, group = "reef_year_id")$summary
    tibble::tibble(
      rows = nrow(.x), positives = sum(.x$truth),
      overall_auc = overall$roc_auc, pr_auc = overall$pr_auc,
      brier = overall$brier, calibration_slope = overall$calibration_slope,
      within_reef_year_auc = within$within_group_roc_auc,
      within_pairs = within$positive_negative_pairs,
      probability_sd = stats::sd(.x$probability)
    )
  }) |>
  dplyr::ungroup()
utils::write.csv(metrics, file.path(
  cfg$paths$outputs_dir, "hmsc_core_120m_blocked_metrics.csv"
), row.names = FALSE)
print(metrics)
