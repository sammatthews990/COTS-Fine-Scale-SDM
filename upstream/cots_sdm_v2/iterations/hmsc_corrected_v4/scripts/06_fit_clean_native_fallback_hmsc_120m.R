args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v4_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v4_dir, "R", "bootstrap.R"))
ctx <- initialise_v4(script_file)
cfg <- ctx$config

data <- readRDS(file.path(cfg$paths$artifacts_dir, "corrected_cell_units_120m.rds"))
binary_responses <- c(
  "cull_high_density", "cull_direct", "salad_high_density", "salad_direct",
  "salad_scar", "rhis_adult", "rhis_juvenile", "rhis_branching_scar",
  "rhis_table_scar"
)
responses <- c(
  binary_responses, "salad_log_density",
  "rhis_branching_live_asin", "rhis_table_live_asin"
)
clean_native_terms <- c(
  "BPI", "EAST", "HCU", "NORTH", "SLO", "SVF", "VCU", "VRM", "DEM",
  "AhyaD", "Aspat", "Aten"
)
data <- data[
  stats::complete.cases(data[, clean_native_terms, drop = FALSE]) &
    rowSums(!is.na(data[, responses, drop = FALSE])) > 0L,
  , drop = FALSE
]
XData <- dplyr::select(data, dplyr::all_of(clean_native_terms))
keep <- vapply(XData, function(x) stats::sd(x) > 1e-10, logical(1))
XData <- XData[, keep, drop = FALSE]
assert_deployable_habitat_terms(names(XData), cfg$model$forbidden_habitat_terms)
XFormula <- stats::as.formula(paste("~", paste(names(XData), collapse = " + ")))
Y <- as.matrix(data[, responses, drop = FALSE])
storage.mode(Y) <- "double"
traits <- data.frame(
  method = factor(dplyr::case_when(
    grepl("^cull", responses) ~ "cull",
    grepl("^salad", responses) ~ "salad",
    grepl("^rhis", responses) ~ "rhis",
    TRUE ~ "other"
  )),
  response_block = factor(ifelse(
    responses %in% c("rhis_branching_live_asin", "rhis_table_live_asin"),
    "coral_resource", "cots_evidence"
  )),
  target_response = as.integer(responses == "cull_high_density"),
  row.names = responses
)
model <- Hmsc::Hmsc(
  Y = Y, XData = as.data.frame(XData), XFormula = XFormula,
  TrData = traits, TrFormula = ~ method + response_block + target_response,
  distr = ifelse(responses %in% binary_responses, "probit", "normal"),
  XScale = TRUE, YScale = TRUE
)
started <- Sys.time()
model <- Hmsc::sampleMcmc(
  model, samples = cfg$model$preliminary_samples,
  thin = cfg$model$preliminary_thin,
  transient = cfg$model$preliminary_transient,
  nChains = cfg$model$preliminary_chains,
  nParallel = cfg$model$preliminary_chains,
  verbose = 0
)
elapsed <- as.numeric(difftime(Sys.time(), started, units = "mins"))
name <- "clean_native_joint"
saveRDS(model, file.path(cfg$paths$models_dir, paste0(name, "_120m.rds")))
saveRDS(list(
  name = name, domain = "clean_native", data = data, XData = XData,
  XFormula = XFormula, pca = list(), responses = responses,
  binary_responses = binary_responses,
  distribution = ifelse(responses %in% binary_responses, "probit", "normal"),
  support_m = 120L, elapsed_minutes = elapsed,
  active_reef_year_conditioning = TRUE,
  fallback_priority = 2L
), file.path(cfg$paths$artifacts_dir, paste0(name, "_fit_data_120m.rds")))
coda_beta <- Hmsc::convertToCodaObject(model)$Beta
ess <- coda::effectiveSize(coda_beta)
psrf <- tryCatch(
  coda::gelman.diag(coda_beta, multivariate = FALSE)$psrf,
  error = function(e) matrix(NA_real_, nrow = length(ess), ncol = 2)
)
convergence <- tibble::tibble(
  model = name, parameter = names(ess), effective_sample_size = as.numeric(ess),
  psrf_point = as.numeric(psrf[, 1]), psrf_upper = as.numeric(psrf[, 2])
)
predicted <- Hmsc::computePredictedValues(model, expected = TRUE)
prediction_mean <- if (length(dim(predicted)) == 3L) apply(predicted, c(1, 2), mean) else as.matrix(predicted)
colnames(prediction_mean) <- responses
target <- which(!is.na(data$cull_high_density))
prediction <- tibble::tibble(
  model = name, validation = "apparent_fit", cell_id = data$cell_id[target],
  reef_key = data$reef_key[target], reef_year_id = data$reef_year_id[target],
  year = data$year[target], truth = data$cull_high_density[target],
  probability = prediction_mean[target, "cull_high_density"]
)
metric <- tibble::tibble(
  model = name, domain = "clean_native", responses = length(responses), rows = nrow(data),
  target_rows = length(target), positives = sum(prediction$truth),
  overall_auc = safe_auc_v4(prediction$truth, prediction$probability),
  within_reef_year_auc = within_reef_year_auc_v4(prediction),
  probability_sd = stats::sd(prediction$probability), elapsed_minutes = elapsed,
  minimum_beta_ess = min(ess, na.rm = TRUE),
  maximum_psrf_upper = max(psrf[, 2], na.rm = TRUE),
  convergence_pass = min(ess, na.rm = TRUE) >= 400 & max(psrf[, 2], na.rm = TRUE) <= 1.1
)
saveRDS(prediction, file.path(cfg$paths$outputs_dir, paste0(name, "_apparent_predictions_120m.rds")))
utils::write.csv(convergence, file.path(cfg$paths$outputs_dir, paste0(name, "_convergence_120m.csv")), row.names = FALSE)
utils::write.csv(metric, file.path(cfg$paths$outputs_dir, paste0(name, "_apparent_metrics_120m.csv")), row.names = FALSE)
print(metric)
