args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v2_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v2_dir, "R", "bootstrap.R"))
ctx <- initialise_v2(script_file)
cfg <- ctx$config

if (!requireNamespace("Hmsc", quietly = TRUE)) stop("Hmsc is required.")
data_path <- file.path(cfg$paths$artifacts_dir, "hmsc_v2_model_data.rds")
if (!file.exists(data_path)) stop("Run script 36 before script 42.")
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
responses <- c(binary_responses, normal_responses)
feature_names <- c(
  setdiff(grep("^geom_", names(data), value = TRUE), "geom_unknown"),
  "coral_suitability_mean", "RG_turbid", "RG_waves_Hs", "RG_waves_Tp",
  grep("_g360m$", names(data), value = TRUE)
)
feature_names <- unique(feature_names[feature_names %in% names(data)])
finite_feature <- vapply(data[feature_names], function(x) {
  all(is.na(x) | is.finite(x)) && stats::sd(x, na.rm = TRUE) > 1e-10
}, logical(1))
feature_names <- feature_names[finite_feature]
fit_data <- data[
  stats::complete.cases(data[, feature_names, drop = FALSE]) &
    rowSums(!is.na(data[, responses, drop = FALSE])) > 0L,
  , drop = FALSE
]

trait_table <- function(responses) {
  method <- sub("_.*$", "", responses)
  evidence <- dplyr::case_when(
    grepl("high_density", responses) ~ "high_density",
    grepl("scar", responses) ~ "scar",
    grepl("live_asin", responses) ~ "coral_resource",
    grepl("density", responses) ~ "density",
    grepl("juvenile", responses) ~ "juvenile_direct",
    TRUE ~ "adult_or_unsized_direct"
  )
  result <- data.frame(
    method = factor(method), evidence = factor(evidence),
    scar_uncertain = as.integer(grepl("scar", responses)),
    coral_resource = as.integer(grepl("live_asin", responses)),
    stringsAsFactors = TRUE
  )
  rownames(result) <- responses
  result
}

Y <- as.matrix(fit_data[, responses, drop = FALSE])
storage.mode(Y) <- "double"
distribution <- ifelse(responses %in% binary_responses, "probit", "normal")
model <- Hmsc::Hmsc(
  Y = Y,
  XData = as.data.frame(fit_data[, feature_names, drop = FALSE]),
  XFormula = stats::as.formula(paste("~", paste(feature_names, collapse = " + "))),
  TrData = trait_table(responses),
  TrFormula = ~ method + evidence + scar_uncertain + coral_resource,
  distr = distribution, XScale = TRUE, YScale = TRUE
)
message(
  "Fitting reduced geomorphic g360 Hmsc on ", nrow(fit_data), " rows, ",
  length(feature_names), " predictors and ", length(responses), " responses."
)
started <- Sys.time()
model <- Hmsc::sampleMcmc(
  model, thin = 5L, samples = 2500L, transient = 2000L,
  nChains = 4L, nParallel = 4L, verbose = 0
)
elapsed_minutes <- as.numeric(difftime(Sys.time(), started, units = "mins"))
saveRDS(
  model,
  file.path(cfg$paths$models_dir, "hmsc_v2_geomorphic_g360_core_model.rds")
)
saveRDS(
  fit_data,
  file.path(cfg$paths$artifacts_dir, "hmsc_v2_geomorphic_g360_fit_data.rds")
)

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
  covariate = c("(Intercept)", feature_names), response = responses,
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

target_rows <- which(!is.na(fit_data$cull_high_density))
target_index <- match("cull_high_density", model$spNames)
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
  site_id = fit_data$cell_id[target_rows],
  reef_name = fit_data$reef_name[target_rows],
  reef_year_id = fit_data$reef_year_id[target_rows],
  year = fit_data$year[target_rows],
  truth = fit_data$cull_high_density[target_rows],
  probability = probability,
  longitude = fit_data$longitude[target_rows],
  latitude = fit_data$latitude[target_rows]
)
metrics <- binary_probability_metrics(
  target_predictions$truth, target_predictions$probability
)
decomposition <- prediction_scale_decomposition(
  target_predictions, group = "reef_year_id"
)
minimum_ess <- min(convergence$effective_sample_size, na.rm = TRUE)
maximum_psrf <- max(convergence$psrf_upper, na.rm = TRUE)
audit <- tibble::tibble(
  rows = nrow(fit_data), cells = dplyr::n_distinct(fit_data$cell_id),
  reefs = dplyr::n_distinct(fit_data$reef_key), years = dplyr::n_distinct(fit_data$year),
  responses = length(responses), predictors = length(feature_names),
  selected_support_m = 360L, incumbent_rank_used = FALSE,
  reef_identity_in_X = FALSE, cell_identity_in_X = FALSE,
  realised_effort_in_X = FALSE, elapsed_minutes = elapsed_minutes,
  retained_draws = draw_count,
  apparent_prediction_mode = "streamed_target_only",
  minimum_beta_ess = minimum_ess, maximum_beta_psrf_upper = maximum_psrf,
  production_convergence_gate = minimum_ess >= 400 & maximum_psrf <= 1.1
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
message("Reduced geomorphic g360 Hmsc fit complete.")
print(audit)
print(metrics)
print(decomposition)
