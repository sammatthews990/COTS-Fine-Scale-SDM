args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v2_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v2_dir, "R", "bootstrap.R"))
ctx <- initialise_v2(script_file)
cfg <- ctx$config

if (!requireNamespace("Hmsc", quietly = TRUE)) stop("Hmsc is required.")
model_path <- file.path(
  cfg$paths$models_dir, "hmsc_v2_geomorphic_g360_core_model.rds"
)
data_path <- file.path(
  cfg$paths$artifacts_dir, "hmsc_v2_geomorphic_g360_fit_data.rds"
)
if (!all(file.exists(c(model_path, data_path)))) {
  stop("Run script 42 before script 43.")
}
model <- readRDS(model_path)
fit_data <- readRDS(data_path)
target <- "cull_high_density"
target_index <- match(target, model$spNames)
target_rows <- which(!is.na(fit_data[[target]]))

coda_object <- Hmsc::convertToCodaObject(model)$Beta
ess <- coda::effectiveSize(coda_object)
psrf <- tryCatch(
  coda::gelman.diag(coda_object, multivariate = FALSE)$psrf,
  error = function(e) matrix(NA_real_, nrow = length(ess), ncol = 2)
)
convergence <- tibble::tibble(
  parameter = names(ess), effective_sample_size = as.numeric(ess),
  psrf_point = as.numeric(psrf[, 1]), psrf_upper = as.numeric(psrf[, 2])
)

beta <- Hmsc::getPostEstimate(model, parName = "Beta")
beta_table <- expand.grid(
  covariate = model$covNames, response = model$spNames,
  stringsAsFactors = FALSE
) |>
  tibble::as_tibble() |>
  dplyr::mutate(
    posterior_mean = as.numeric(beta$mean),
    posterior_support_positive = as.numeric(beta$support)
  )
predictor_sd <- stats::setNames(
  as.numeric(model$XScalePar[2, seq_along(model$covNames)]), model$covNames
)
beta_table <- beta_table |>
  dplyr::mutate(
    predictor_sd = unname(predictor_sd[.data$covariate]),
    standardised_effect = .data$posterior_mean * .data$predictor_sd
  )

message(
  "Streaming target predictions for ", length(target_rows),
  " observed Cull rows across ", sum(lengths(model$postList)), " draws."
)
X_target <- model$X[target_rows, , drop = FALSE]
probability_sum <- numeric(length(target_rows))
draw_count <- 0L
for (chain in model$postList) {
  for (draw in chain) {
    probability_sum <- probability_sum + stats::pnorm(as.numeric(
      X_target %*% draw$Beta[, target_index]
    ))
    draw_count <- draw_count + 1L
  }
}
probability <- probability_sum / draw_count
target_predictions <- tibble::tibble(
  model = "hmsc_v2_geomorphic_g360",
  validation = "apparent_fit_not_promotion_evidence",
  site_year_id = paste(
    fit_data$cell_id[target_rows], fit_data$year[target_rows], sep = "||"
  ),
  site_id = fit_data$cell_id[target_rows], reef_name = fit_data$reef_name[target_rows],
  reef_year_id = fit_data$reef_year_id[target_rows], year = fit_data$year[target_rows],
  truth = fit_data[[target]][target_rows], probability = probability,
  longitude = fit_data$longitude[target_rows], latitude = fit_data$latitude[target_rows]
)
metrics <- binary_probability_metrics(
  target_predictions$truth, target_predictions$probability
)
decomposition <- prediction_scale_decomposition(
  target_predictions, group = "reef_year_id"
)
minimum_ess <- min(convergence$effective_sample_size, na.rm = TRUE)
maximum_psrf <- max(convergence$psrf_upper, na.rm = TRUE)
feature_names <- names(model$XData)
audit <- tibble::tibble(
  rows = nrow(fit_data), cells = dplyr::n_distinct(fit_data$cell_id),
  reefs = dplyr::n_distinct(fit_data$reef_key), years = dplyr::n_distinct(fit_data$year),
  responses = model$ns, predictors = length(feature_names), selected_support_m = 360L,
  incumbent_rank_used = FALSE, reef_identity_in_X = FALSE,
  cell_identity_in_X = FALSE, realised_effort_in_X = FALSE,
  chains = length(model$postList), transient = model$transient,
  thin = model$thin, retained_samples_per_chain = model$samples,
  retained_draws = draw_count, minimum_beta_ess = minimum_ess,
  maximum_beta_psrf_upper = maximum_psrf,
  production_convergence_gate = minimum_ess >= 400 & maximum_psrf <= 1.1,
  apparent_prediction_mode = "streamed_target_only"
)

saveRDS(
  target_predictions,
  file.path(cfg$paths$outputs_dir, "hmsc_v2_geomorphic_g360_apparent_predictions.rds")
)
tables <- list(
  hmsc_v2_geomorphic_g360_audit = audit,
  hmsc_v2_geomorphic_g360_convergence = convergence,
  hmsc_v2_geomorphic_g360_beta = beta_table,
  hmsc_v2_geomorphic_g360_apparent_metrics = metrics,
  hmsc_v2_geomorphic_g360_apparent_decomposition = decomposition,
  hmsc_v2_geomorphic_g360_feature_manifest = tibble::tibble(feature = feature_names)
)
for (name in names(tables)) {
  utils::write.csv(
    tables[[name]], file.path(cfg$paths$outputs_dir, paste0(name, ".csv")),
    row.names = FALSE
  )
}
message("Long reduced geomorphic g360 Hmsc summary complete.")
print(audit)
print(metrics)
print(decomposition)
