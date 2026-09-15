args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v4_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v4_dir, "R", "bootstrap.R"))
ctx <- initialise_v4(script_file)
cfg <- ctx$config
name <- "enhanced_cots_manta_strong"
model <- readRDS(file.path(cfg$paths$models_dir, paste0(name, "_120m.rds")))
artifact <- readRDS(file.path(cfg$paths$artifacts_dir, paste0(name, "_fit_data_120m.rds")))
data <- artifact$data

coda_beta <- Hmsc::convertToCodaObject(model)$Beta
ess <- coda::effectiveSize(coda_beta)
psrf <- tryCatch(
  coda::gelman.diag(coda_beta, multivariate = FALSE)$psrf,
  error = function(e) matrix(NA_real_, nrow = length(ess), ncol = 2)
)
prediction_draws <- Hmsc::computePredictedValues(model, expected = TRUE)
prediction_mean <- if (length(dim(prediction_draws)) == 3L) {
  apply(prediction_draws, c(1, 2), mean)
} else as.matrix(prediction_draws)
colnames(prediction_mean) <- artifact$responses
target <- which(!is.na(data$cull_high_density))
prediction <- tibble::tibble(
  model = name, validation = "apparent_fit", cell_id = data$cell_id[target],
  reef_key = data$reef_key[target], reef_year_id = data$reef_year_id[target],
  year = data$year[target], truth = data$cull_high_density[target],
  probability = prediction_mean[target, "cull_high_density"]
)
is_manta <- data$observation_support == "manta_tow_footprint"
metric <- tibble::tibble(
  model = name, responses = length(artifact$responses), rows = nrow(data),
  localized_rows = sum(!is_manta), manta_rows = sum(is_manta),
  manta_positives = sum(data$manta_strong[is_manta] == 1L, na.rm = TRUE),
  target_rows = length(target), positives = sum(prediction$truth),
  overall_auc = safe_auc_v4(prediction$truth, prediction$probability),
  within_reef_year_auc = within_reef_year_auc_v4(prediction),
  probability_sd = stats::sd(prediction$probability),
  elapsed_minutes = artifact$elapsed_minutes,
  minimum_beta_ess = min(ess, na.rm = TRUE),
  maximum_psrf_upper = max(psrf[, 2], na.rm = TRUE),
  convergence_pass = min(ess, na.rm = TRUE) >= 400 & max(psrf[, 2], na.rm = TRUE) <= 1.1
)
saveRDS(prediction, file.path(cfg$paths$outputs_dir, paste0(name, "_apparent_predictions_120m.rds")))
utils::write.csv(metric, file.path(cfg$paths$outputs_dir, paste0(name, "_metrics_120m.csv")), row.names = FALSE)
print(metric)
