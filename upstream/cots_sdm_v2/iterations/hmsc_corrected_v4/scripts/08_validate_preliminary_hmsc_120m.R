args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
runtime_args <- commandArgs(trailingOnly = TRUE)
if (length(runtime_args) < 2L) {
  stop("Usage: Rscript 08_validate_preliminary_hmsc_120m.R MODEL VALIDATION")
}
model_name <- runtime_args[[1]]
validation <- runtime_args[[2]]
if (!model_name %in% c("enhanced_cots", "enhanced_joint", "clean_joint")) {
  stop("Unknown model: ", model_name)
}
if (!validation %in% c("held_cell", "leave_reef", "forward")) {
  stop("Unknown validation: ", validation)
}
v4_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v4_dir, "R", "bootstrap.R"))
ctx <- initialise_v4(script_file)
cfg <- ctx$config

model <- readRDS(file.path(cfg$paths$models_dir, paste0(model_name, "_120m.rds")))
artifact <- readRDS(file.path(cfg$paths$artifacts_dir, paste0(model_name, "_fit_data_120m.rds")))
data <- artifact$data
target_rows <- which(!is.na(data$cull_high_density))

make_reef_folds <- function(data, k = 5L) {
  reefs <- data |>
    dplyr::group_by(.data$reef_key) |>
    dplyr::summarise(
      x = mean(.data$projected_x), y = mean(.data$projected_y), .groups = "drop"
    )
  set.seed(20260909L)
  reefs$fold <- stats::kmeans(scale(reefs[, c("x", "y")]), k, nstart = 50)$cluster
  reefs$fold[match(data$reef_key, reefs$reef_key)]
}

if (validation == "held_cell") {
  partition <- make_spatial_fold_v4(
    data$projected_x, data$projected_y,
    cfg$validation$spatial_block_m, cfg$validation$folds
  )
  assess <- sort(unique(partition))
  label <- "held_3km_cell_5fold"
} else if (validation == "leave_reef") {
  partition <- make_reef_folds(data, cfg$validation$folds)
  assess <- sort(unique(partition))
  label <- "leave_reef_5fold"
} else {
  partition <- ifelse(data$year <= 2023, 1L, 2L)
  assess <- 2L
  label <- "forward_train_through_2023"
}

validation_model <- model
validation_model$samples <- 100L
validation_model$thin <- 5L
validation_model$transient <- 300L
message("Running ", label, " for ", model_name)
started <- Sys.time()
draws <- Hmsc::computePredictedValues(
  validation_model, partition = partition, expected = TRUE,
  nChains = 4L, nParallel = 4L, verbose = 0
)
elapsed <- as.numeric(difftime(Sys.time(), started, units = "mins"))
response_index <- match("cull_high_density", artifact$responses)
probability <- apply(draws, c(1, 2), mean)[, response_index]
keep <- target_rows[partition[target_rows] %in% assess]
predictions <- tibble::tibble(
  model = model_name, validation = label, fold = partition[keep],
  cell_id = data$cell_id[keep], reef_key = data$reef_key[keep],
  reef_year_id = data$reef_year_id[keep], year = data$year[keep],
  truth = data$cull_high_density[keep], probability = probability[keep]
)
metric <- tibble::tibble(
  model = model_name, validation = label, rows = nrow(predictions),
  positives = sum(predictions$truth),
  overall_auc = safe_auc_v4(predictions$truth, predictions$probability),
  within_reef_year_auc = within_reef_year_auc_v4(predictions),
  probability_sd = stats::sd(predictions$probability), elapsed_minutes = elapsed,
  cv_chain_samples = validation_model$samples,
  cv_transient = validation_model$transient
)
stem <- paste(model_name, validation, "120m", sep = "_")
saveRDS(predictions, file.path(cfg$paths$outputs_dir, paste0(stem, "_predictions.rds")))
utils::write.csv(metric, file.path(cfg$paths$outputs_dir, paste0(stem, "_metrics.csv")), row.names = FALSE)
print(metric)
