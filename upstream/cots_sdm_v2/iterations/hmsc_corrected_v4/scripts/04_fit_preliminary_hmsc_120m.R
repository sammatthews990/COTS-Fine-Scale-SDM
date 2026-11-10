args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v4_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v4_dir, "R", "bootstrap.R"))
ctx <- initialise_v4(script_file)
cfg <- ctx$config

source_data <- readRDS(file.path(cfg$paths$artifacts_dir, "corrected_cell_units_120m.rds"))
binary_responses <- c(
  "cull_high_density", "cull_direct",
  "salad_high_density", "salad_direct", "salad_scar",
  "rhis_adult", "rhis_juvenile", "rhis_branching_scar", "rhis_table_scar"
)
continuous_responses <- c(
  "salad_log_density", "rhis_branching_live_asin", "rhis_table_live_asin"
)
cots_responses <- c(binary_responses, "salad_log_density")
joint_responses <- c(binary_responses, continuous_responses)
families <- terrain_families(names(source_data))
raw_terms <- unlist(families, use.names = FALSE)

variants <- list(
  enhanced_cots = list(
    domain = "enhanced", responses = cots_responses,
    base_terms = c("geomorphic_class", "AhyaD", "Aspat", "Aten", "turbid", "waves_Hs", "waves_Tp")
  ),
  enhanced_joint = list(
    domain = "enhanced", responses = joint_responses,
    base_terms = c("geomorphic_class", "AhyaD", "Aspat", "Aten", "turbid", "waves_Hs", "waves_Tp")
  ),
  clean_joint = list(
    domain = "clean", responses = joint_responses,
    base_terms = c("AhyaD", "Aspat", "Aten")
  )
)

trait_table <- function(responses) {
  method <- dplyr::case_when(
    grepl("^cull", responses) ~ "cull",
    grepl("^salad", responses) ~ "salad",
    grepl("^rhis", responses) ~ "rhis",
    TRUE ~ "other"
  )
  data.frame(
    method = factor(method),
    response_block = factor(ifelse(
      responses %in% c("rhis_branching_live_asin", "rhis_table_live_asin"),
      "coral_resource", "cots_evidence"
    )),
    target_response = as.integer(responses == "cull_high_density"),
    row.names = responses
  )
}

fit_variant <- function(name, specification) {
  message("Fitting preliminary Hmsc variant: ", name)
  required <- unique(c(specification$base_terms, raw_terms))
  data <- source_data[
    stats::complete.cases(source_data[, required, drop = FALSE]) &
      rowSums(!is.na(source_data[, specification$responses, drop = FALSE])) > 0L,
    , drop = FALSE
  ]
  pca <- fit_family_pca(data, variance = 0, max_components = 1L)
  pcs <- apply_family_pca(data, pca)
  XData <- dplyr::bind_cols(
    dplyr::select(data, dplyr::all_of(specification$base_terms)), pcs
  )
  if ("geomorphic_class" %in% names(XData)) {
    XData$geomorphic_class <- factor(XData$geomorphic_class)
  }
  continuous <- setdiff(names(XData), "geomorphic_class")
  valid <- vapply(XData[continuous], function(x) stats::sd(x) > 1e-10, logical(1))
  continuous <- continuous[valid]
  XData <- dplyr::select(
    XData, dplyr::any_of("geomorphic_class"), dplyr::all_of(continuous)
  )
  assert_deployable_habitat_terms(names(XData), cfg$model$forbidden_habitat_terms)
  XFormula <- stats::as.formula(paste("~", paste(names(XData), collapse = " + ")))

  responses <- specification$responses
  Y <- as.matrix(data[, responses, drop = FALSE])
  storage.mode(Y) <- "double"
  distribution <- ifelse(responses %in% binary_responses, "probit", "normal")
  traits <- trait_table(responses)
  trait_formula <- if (length(unique(traits$response_block)) > 1L) {
    ~ method + response_block + target_response
  } else {
    ~ method + target_response
  }
  model <- Hmsc::Hmsc(
    Y = Y, XData = as.data.frame(XData), XFormula = XFormula,
    TrData = traits, TrFormula = trait_formula,
    distr = distribution, XScale = TRUE, YScale = TRUE
  )
  started <- Sys.time()
  model <- Hmsc::sampleMcmc(
    model,
    samples = cfg$model$preliminary_samples,
    thin = cfg$model$preliminary_thin,
    transient = cfg$model$preliminary_transient,
    nChains = cfg$model$preliminary_chains,
    nParallel = cfg$model$preliminary_chains,
    verbose = 0
  )
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "mins"))
  saveRDS(model, file.path(cfg$paths$models_dir, paste0(name, "_120m.rds")))
  artifact <- list(
    name = name, domain = specification$domain, data = data, XData = XData,
    XFormula = XFormula, pca = pca, responses = responses,
    binary_responses = binary_responses, distribution = distribution,
    support_m = 120L, elapsed_minutes = elapsed,
    active_reef_year_conditioning = TRUE,
    manta_active_definition = "direct COTS OR >=2 counted scars OR legacy common scars"
  )
  saveRDS(artifact, file.path(cfg$paths$artifacts_dir, paste0(name, "_fit_data_120m.rds")))

  coda_beta <- Hmsc::convertToCodaObject(model)$Beta
  ess <- coda::effectiveSize(coda_beta)
  psrf <- tryCatch(
    coda::gelman.diag(coda_beta, multivariate = FALSE)$psrf,
    error = function(e) matrix(NA_real_, nrow = length(ess), ncol = 2)
  )
  convergence <- tibble::tibble(
    model = name, parameter = names(ess),
    effective_sample_size = as.numeric(ess),
    psrf_point = as.numeric(psrf[, 1]), psrf_upper = as.numeric(psrf[, 2])
  )
  predicted <- Hmsc::computePredictedValues(model, expected = TRUE)
  prediction_mean <- if (length(dim(predicted)) == 3L) {
    apply(predicted, c(1, 2), mean)
  } else as.matrix(predicted)
  colnames(prediction_mean) <- responses
  target <- which(!is.na(data$cull_high_density))
  prediction <- tibble::tibble(
    model = name, validation = "apparent_fit", cell_id = data$cell_id[target],
    reef_key = data$reef_key[target], reef_year_id = data$reef_year_id[target],
    year = data$year[target], truth = data$cull_high_density[target],
    probability = prediction_mean[target, "cull_high_density"]
  )
  metric <- tibble::tibble(
    model = name, domain = specification$domain, responses = length(responses),
    rows = nrow(data), target_rows = length(target), positives = sum(prediction$truth),
    overall_auc = safe_auc_v4(prediction$truth, prediction$probability),
    within_reef_year_auc = within_reef_year_auc_v4(prediction),
    probability_sd = stats::sd(prediction$probability),
    elapsed_minutes = elapsed,
    minimum_beta_ess = min(ess, na.rm = TRUE),
    maximum_psrf_upper = max(psrf[, 2], na.rm = TRUE),
    convergence_pass = min(ess, na.rm = TRUE) >= 400 & max(psrf[, 2], na.rm = TRUE) <= 1.1
  )
  saveRDS(prediction, file.path(cfg$paths$outputs_dir, paste0(name, "_apparent_predictions_120m.rds")))
  utils::write.csv(convergence, file.path(cfg$paths$outputs_dir, paste0(name, "_convergence_120m.csv")), row.names = FALSE)
  list(metric = metric, convergence = convergence)
}

results <- lapply(names(variants), function(name) fit_variant(name, variants[[name]]))
metrics <- dplyr::bind_rows(lapply(results, `[[`, "metric"))
utils::write.csv(metrics, file.path(cfg$paths$outputs_dir, "preliminary_hmsc_metrics_120m.csv"), row.names = FALSE)
print(metrics)
