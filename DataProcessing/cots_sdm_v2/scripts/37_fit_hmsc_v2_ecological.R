args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v2_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v2_dir, "R", "bootstrap.R"))
ctx <- initialise_v2(script_file)
cfg <- ctx$config

if (!requireNamespace("Hmsc", quietly = TRUE)) stop("Hmsc is required.")
data_path <- file.path(cfg$paths$artifacts_dir, "hmsc_v2_model_data.rds")
if (!file.exists(data_path)) stop("Run script 36 before fitting Hmsc v2.")
data <- readRDS(data_path)

binary_responses <- c(
  "cull_direct_presence", "cull_high_density",
  "manta_direct_presence", "manta_scar_presence",
  "rhis_adult_presence", "rhis_juvenile_presence",
  "rhis_branching_scar_presence", "rhis_table_scar_presence",
  "salad_direct_presence", "salad_scar_presence"
)
normal_responses <- c(
  "rhis_branching_live_asin", "rhis_table_live_asin",
  "salad_log_density_incl_scars"
)
mortality_responses <- c(
  "rhis_branching_recent_dead_asin", "rhis_table_recent_dead_asin"
)
core_responses <- c(binary_responses, normal_responses)

gaussian_features <- grep(
  "_(g30m|g60m|g120m|g240m|g360m)$", names(data), value = TRUE
)
geomorphic_features <- grep("^geom_", names(data), value = TRUE)
geomorphic_features <- setdiff(
  geomorphic_features, "geom_unknown"
)
feature_names <- c(
  gaussian_features,
  "coral_suitability_mean", "RG_turbid", "RG_waves_Hs", "RG_waves_Tp",
  geomorphic_features
)
feature_names <- unique(feature_names[feature_names %in% names(data)])

finite_feature <- vapply(data[feature_names], function(x) {
  all(is.na(x) | is.finite(x)) && stats::sd(x, na.rm = TRUE) > 1e-10
}, logical(1))
feature_names <- feature_names[finite_feature]
complete_predictors <- stats::complete.cases(data[, feature_names, drop = FALSE])
response_observed <- rowSums(!is.na(data[, core_responses, drop = FALSE])) > 0L
fit_data <- data[complete_predictors & response_observed, , drop = FALSE]
if (nrow(fit_data) < 1000L) stop("Too few predictor-complete Hmsc v2 rows.")

trait_table <- function(responses) {
  method <- sub("_.*$", "", responses)
  evidence <- dplyr::case_when(
    grepl("high_density", responses) ~ "high_density",
    grepl("scar", responses) ~ "scar",
    grepl("live_asin", responses) ~ "coral_resource",
    grepl("recent_dead", responses) ~ "mortality_context",
    grepl("density", responses) ~ "density",
    grepl("juvenile", responses) ~ "juvenile_direct",
    TRUE ~ "adult_or_unsized_direct"
  )
  result <- data.frame(
    method = factor(method),
    evidence = factor(evidence),
    scar_uncertain = as.integer(grepl("scar", responses)),
    coral_resource = as.integer(grepl("live_asin", responses)),
    mortality_nonspecific = as.integer(grepl("recent_dead", responses)),
    stringsAsFactors = TRUE
  )
  rownames(result) <- responses
  result
}

fit_model <- function(responses, label, chains, samples, thin, transient) {
  Y <- as.matrix(fit_data[, responses, drop = FALSE])
  storage.mode(Y) <- "double"
  XData <- as.data.frame(fit_data[, feature_names, drop = FALSE])
  distribution <- ifelse(responses %in% binary_responses, "probit", "normal")
  tr_formula <- if (any(grepl("recent_dead", responses))) {
    ~ method + evidence + scar_uncertain + coral_resource +
      mortality_nonspecific
  } else {
    ~ method + evidence + scar_uncertain + coral_resource
  }
  model <- Hmsc::Hmsc(
    Y = Y,
    XData = XData,
    XFormula = stats::as.formula(paste("~", paste(feature_names, collapse = " + "))),
    TrData = trait_table(responses),
    TrFormula = tr_formula,
    distr = distribution,
    XScale = TRUE,
    YScale = TRUE
  )
  started <- Sys.time()
  model <- Hmsc::sampleMcmc(
    model, thin = thin, samples = samples, transient = transient,
    nChains = chains, nParallel = chains, verbose = 0
  )
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "mins"))
  saveRDS(model, file.path(cfg$paths$models_dir, paste0(label, ".rds")))
  list(model = model, elapsed_minutes = elapsed, responses = responses)
}

message(
  "Fitting fixed-effect ecological Hmsc on ", nrow(fit_data), " rows, ",
  length(feature_names), " predictors and ", length(core_responses), " responses."
)
core <- fit_model(
  core_responses, "hmsc_v2_ecological_core_model",
  chains = 4L, samples = 250L, thin = 5L, transient = 1000L
)
message("Fitting nonspecific coral-mortality sensitivity model.")
sensitivity <- fit_model(
  c(core_responses, mortality_responses),
  "hmsc_v2_ecological_mortality_sensitivity_model",
  chains = 2L, samples = 100L, thin = 5L, transient = 500L
)

summarise_model <- function(fit, label) {
  coda_object <- Hmsc::convertToCodaObject(fit$model)$Beta
  ess <- coda::effectiveSize(coda_object)
  psrf <- tryCatch(
    coda::gelman.diag(coda_object, multivariate = FALSE)$psrf,
    error = function(e) matrix(NA_real_, nrow = length(ess), ncol = 2)
  )
  beta <- Hmsc::getPostEstimate(fit$model, parName = "Beta")
  beta_table <- expand.grid(
    covariate = c("(Intercept)", feature_names),
    response = fit$responses,
    stringsAsFactors = FALSE
  ) |>
    tibble::as_tibble() |>
    dplyr::mutate(
      posterior_mean = as.numeric(beta$mean),
      posterior_support_positive = as.numeric(beta$support),
      model = label
    )
  convergence <- tibble::tibble(
    parameter = names(ess),
    effective_sample_size = as.numeric(ess),
    psrf_point = as.numeric(psrf[, 1]),
    psrf_upper = as.numeric(psrf[, 2]),
    model = label
  )
  list(beta = beta_table, convergence = convergence)
}
core_summary <- summarise_model(core, "core")
sensitivity_summary <- summarise_model(sensitivity, "mortality_sensitivity")

predicted <- Hmsc::computePredictedValues(core$model, expected = TRUE)
prediction_mean <- if (length(dim(predicted)) == 3L) {
  apply(predicted, c(1, 2), mean)
} else {
  as.matrix(predicted)
}
colnames(prediction_mean) <- core_responses
target_index <- which(!is.na(fit_data$cull_high_density))
target_predictions <- tibble::tibble(
  model = "hmsc_v2_ecological_core_apparent_fit",
  validation = "apparent_fit_not_promotion_evidence",
  site_year_id = paste(
    fit_data$cell_id[target_index], fit_data$year[target_index], sep = "||"
  ),
  site_id = fit_data$cell_id[target_index],
  reef_name = fit_data$reef_name[target_index],
  reef_year_id = fit_data$reef_year_id[target_index],
  year = fit_data$year[target_index],
  truth = fit_data$cull_high_density[target_index],
  probability = prediction_mean[target_index, "cull_high_density"],
  longitude = fit_data$longitude[target_index],
  latitude = fit_data$latitude[target_index]
)
apparent_metrics <- binary_probability_metrics(
  target_predictions$truth, target_predictions$probability
) |>
  dplyr::mutate(
    model = "hmsc_v2_ecological_core",
    evaluation = "apparent_fit_not_promotion_evidence"
  )
apparent_decomposition <- suppressWarnings(
  prediction_scale_decomposition(target_predictions, group = "reef_year_id")
)

core_response_count <- length(core_responses)
sensitivity_response_count <- length(c(core_responses, mortality_responses))
core_minimum_beta_ess <- min(
  core_summary$convergence$effective_sample_size, na.rm = TRUE
)
core_maximum_beta_psrf_upper <- max(
  core_summary$convergence$psrf_upper, na.rm = TRUE
)
audit <- tibble::tibble(
  rows = nrow(fit_data),
  cells = dplyr::n_distinct(fit_data$cell_id),
  reefs = dplyr::n_distinct(fit_data$reef_key),
  years = dplyr::n_distinct(fit_data$year),
  core_responses = core_response_count,
  mortality_sensitivity_responses = sensitivity_response_count,
  gaussian_predictors = length(gaussian_features),
  fitted_predictors = length(feature_names),
  reef_identity_in_X = FALSE,
  cell_identity_in_X = FALSE,
  realised_effort_in_X = FALSE,
  connectivity_available = FALSE,
  core_mcmc_minutes = core$elapsed_minutes,
  sensitivity_mcmc_minutes = sensitivity$elapsed_minutes,
  core_minimum_beta_ess = core_minimum_beta_ess,
  core_maximum_beta_psrf_upper = core_maximum_beta_psrf_upper,
  production_convergence_gate = core_minimum_beta_ess >= 400 &
    core_maximum_beta_psrf_upper <= 1.1
)

saveRDS(fit_data, file.path(
  cfg$paths$artifacts_dir, "hmsc_v2_ecological_fit_data.rds"
))
saveRDS(target_predictions, file.path(
  cfg$paths$outputs_dir, "hmsc_v2_apparent_target_predictions.rds"
))
tables <- list(
  hmsc_v2_ecological_audit = audit,
  hmsc_v2_ecological_beta = dplyr::bind_rows(
    core_summary$beta, sensitivity_summary$beta
  ),
  hmsc_v2_ecological_convergence = dplyr::bind_rows(
    core_summary$convergence, sensitivity_summary$convergence
  ),
  hmsc_v2_apparent_metrics = apparent_metrics,
  hmsc_v2_apparent_decomposition = apparent_decomposition,
  hmsc_v2_feature_manifest = tibble::tibble(feature = feature_names)
)
for (name in names(tables)) {
  utils::write.csv(
    tables[[name]], file.path(cfg$paths$outputs_dir, paste0(name, ".csv")),
    row.names = FALSE
  )
}
message("Hmsc v2 ecological and mortality-sensitivity fits complete.")
print(audit)
print(apparent_metrics)
print(apparent_decomposition)
