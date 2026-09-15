args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v3_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v3_dir, "R", "bootstrap.R"))
ctx <- initialise_v3(script_file)
cfg <- ctx$config
support_m <- 120L

responses <- c(
  "cull_high_density", "cull_direct", "salad_high_density",
  "salad_direct", "rhis_adult", "rhis_juvenile"
)
units <- readRDS(file.path(
  cfg$paths$artifacts_dir, paste0("hmsc_core_units_", support_m, "m.rds")
))
predictors <- readRDS(file.path(
  cfg$paths$artifacts_dir, paste0("support_predictors_", support_m, "m.rds")
))
data <- dplyr::left_join(units, predictors, by = c(
  "cell_id", "projected_x", "projected_y"
))
raw_terms <- unlist(terrain_families(names(data)), use.names = FALSE)
base_terms <- c("AhyaD", "Aspat", "Aten", "turbid", "waves_Hs", "waves_Tp")
required <- c("geomorphic_class", base_terms, raw_terms)
data <- data[
  stats::complete.cases(data[, required, drop = FALSE]) &
    rowSums(!is.na(data[, responses, drop = FALSE])) > 0,
  , drop = FALSE
]

pca <- fit_family_pca(data, variance = 0, max_components = 1L)
pcs <- apply_family_pca(data, pca)
XData <- dplyr::bind_cols(
  dplyr::select(data, dplyr::all_of(c("geomorphic_class", base_terms))), pcs
) |>
  dplyr::mutate(geomorphic_class = factor(.data$geomorphic_class))
continuous <- setdiff(names(XData), "geomorphic_class")
valid <- vapply(XData[continuous], function(x) stats::sd(x) > 1e-10, logical(1))
continuous <- continuous[valid]
XData <- dplyr::select(XData, "geomorphic_class", dplyr::all_of(continuous))
assert_deployable_habitat_terms(names(XData), cfg$model$forbidden_habitat_terms)
XFormula <- stats::as.formula(paste(
  "~ geomorphic_class +", paste(continuous, collapse = " + ")
))

Y <- as.matrix(data[, responses, drop = FALSE])
storage.mode(Y) <- "double"
traits <- data.frame(
  method = factor(c("cull", "cull", "salad", "salad", "rhis", "rhis")),
  evidence = factor(c(
    "problem_density", "direct", "density_support", "direct", "adult", "juvenile"
  )),
  target_response = as.integer(responses == "cull_high_density"),
  row.names = responses
)

model <- Hmsc::Hmsc(
  Y = Y, XData = as.data.frame(XData), XFormula = XFormula,
  TrData = traits, TrFormula = ~ method + evidence + target_response,
  distr = rep("probit", length(responses)), XScale = TRUE, YScale = TRUE
)
started <- Sys.time()
model <- Hmsc::sampleMcmc(
  model, samples = 2000L, thin = 5L, transient = 2000L,
  nChains = 4L, nParallel = 4L, verbose = 0
)
elapsed <- as.numeric(difftime(Sys.time(), started, units = "mins"))
saveRDS(model, file.path(cfg$paths$models_dir, "hmsc_core_120m.rds"))
saveRDS(list(
  data = data, XData = XData, pca = pca, responses = responses,
  XFormula = XFormula, support_m = support_m
), file.path(cfg$paths$artifacts_dir, "hmsc_core_120m_fit_data.rds"))

coda_beta <- Hmsc::convertToCodaObject(model)$Beta
ess <- coda::effectiveSize(coda_beta)
psrf <- tryCatch(
  coda::gelman.diag(coda_beta, multivariate = FALSE)$psrf,
  error = function(e) matrix(NA_real_, nrow = length(ess), ncol = 2)
)
convergence <- tibble::tibble(
  parameter = names(ess), effective_sample_size = as.numeric(ess),
  psrf_point = as.numeric(psrf[, 1]), psrf_upper = as.numeric(psrf[, 2])
)
predicted <- Hmsc::computePredictedValues(model, expected = TRUE)
prediction_mean <- if (length(dim(predicted)) == 3L) {
  apply(predicted, c(1, 2), mean)
} else as.matrix(predicted)
colnames(prediction_mean) <- responses
target_rows <- which(!is.na(data$cull_high_density))
target <- tibble::tibble(
  model = "hmsc_geomorphic_core_120m", validation = "apparent_fit",
  cell_id = data$cell_id[target_rows], reef_key = data$reef_key[target_rows],
  reef_year_id = paste(data$reef_key[target_rows], data$year[target_rows], sep = "||"),
  year = data$year[target_rows], truth = data$cull_high_density[target_rows],
  probability = prediction_mean[target_rows, "cull_high_density"]
)
metrics <- binary_probability_metrics(target$truth, target$probability) |>
  dplyr::bind_cols(within_group_auc(target, group = "reef_year_id")$summary) |>
  dplyr::mutate(
    support_m = support_m, elapsed_minutes = elapsed,
    minimum_beta_ess = min(ess, na.rm = TRUE),
    maximum_psrf_upper = max(psrf[, 2], na.rm = TRUE),
    convergence_pass = .data$minimum_beta_ess >= 400 & .data$maximum_psrf_upper <= 1.1
  )
saveRDS(target, file.path(cfg$paths$outputs_dir, "hmsc_core_120m_apparent_predictions.rds"))
utils::write.csv(convergence, file.path(
  cfg$paths$outputs_dir, "hmsc_core_120m_convergence.csv"
), row.names = FALSE)
utils::write.csv(metrics, file.path(
  cfg$paths$outputs_dir, "hmsc_core_120m_apparent_metrics.csv"
), row.names = FALSE)
print(metrics)
