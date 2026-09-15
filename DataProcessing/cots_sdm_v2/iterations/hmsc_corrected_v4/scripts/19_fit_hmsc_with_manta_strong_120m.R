args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v4_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v4_dir, "R", "bootstrap.R"))
ctx <- initialise_v4(script_file)
cfg <- ctx$config
set.seed(20260909L)

localized <- readRDS(file.path(cfg$paths$artifacts_dir, "corrected_cell_units_120m.rds"))
manta <- readRDS(file.path(cfg$paths$artifacts_dir, "manta_footprints_120m.rds")) |>
  dplyr::filter(.data$active_primary, .data$enhanced_complete) |>
  dplyr::mutate(
    manta_strong = dplyr::case_when(
      .data$manta_strong_positive ~ 1L,
      .data$manta_strict_negative ~ 0L,
      TRUE ~ NA_integer_
    )
  ) |>
  dplyr::filter(!is.na(.data$manta_strong)) |>
  dplyr::group_by(.data$reef_year_id) |>
  dplyr::mutate(n_positive = sum(.data$manta_strong == 1L)) |>
  dplyr::arrange(.data$reef_year_id, dplyr::desc(.data$manta_strong), .data$visit_id) |>
  dplyr::filter(
    .data$manta_strong == 1L |
      dplyr::row_number() <= .data$n_positive + dplyr::if_else(
        .data$n_positive > 0L,
        cfg$manta$negative_to_positive_ratio * .data$n_positive,
        10L
      )
  ) |>
  dplyr::ungroup()

binary_responses <- c(
  "cull_high_density", "cull_direct", "salad_high_density", "salad_direct",
  "salad_scar", "rhis_adult", "rhis_juvenile", "rhis_branching_scar",
  "rhis_table_scar", "manta_strong"
)
continuous_responses <- "salad_log_density"
responses <- c(binary_responses, continuous_responses)
base_terms <- c(
  "geomorphic_class", "AhyaD", "Aspat", "Aten",
  "turbid", "waves_Hs", "waves_Tp"
)
families <- terrain_families(names(localized))
raw_terms <- unlist(families, use.names = FALSE)
required <- unique(c(base_terms, raw_terms))

localized <- localized[
  stats::complete.cases(localized[, required, drop = FALSE]) &
    rowSums(!is.na(localized[, setdiff(responses, "manta_strong"), drop = FALSE])) > 0L,
  , drop = FALSE
]
localized$manta_strong <- NA_integer_

manta_rows <- manta |>
  dplyr::transmute(
    cell_id = paste0("manta||", .data$visit_id),
    year, reef_key, reef_name, projected_x, projected_y, reef_year_id,
    dplyr::across(dplyr::all_of(required)),
    cull_high_density = NA_integer_, cull_direct = NA_integer_,
    salad_high_density = NA_integer_, salad_direct = NA_integer_,
    salad_scar = NA_integer_, rhis_adult = NA_integer_,
    rhis_juvenile = NA_integer_, rhis_branching_scar = NA_integer_,
    rhis_table_scar = NA_integer_, manta_strong = as.integer(.data$manta_strong),
    salad_log_density = NA_real_, observation_support = "manta_tow_footprint"
  )
localized_rows <- localized |>
  dplyr::select(
    cell_id, year, reef_key, reef_name, projected_x, projected_y, reef_year_id,
    dplyr::all_of(required), dplyr::all_of(responses), observation_support
  )
data <- dplyr::bind_rows(localized_rows, manta_rows)

pca <- fit_family_pca(data, variance = 0, max_components = 1L)
pcs <- apply_family_pca(data, pca)
XData <- dplyr::bind_cols(dplyr::select(data, dplyr::all_of(base_terms)), pcs)
XData$geomorphic_class <- factor(XData$geomorphic_class)
continuous <- setdiff(names(XData), "geomorphic_class")
continuous <- continuous[vapply(XData[continuous], function(x) stats::sd(x) > 1e-10, logical(1))]
XData <- dplyr::select(XData, "geomorphic_class", dplyr::all_of(continuous))
assert_deployable_habitat_terms(names(XData), cfg$model$forbidden_habitat_terms)
XFormula <- stats::as.formula(paste("~", paste(names(XData), collapse = " + ")))

Y <- as.matrix(data[, responses, drop = FALSE])
storage.mode(Y) <- "double"
distribution <- ifelse(responses %in% binary_responses, "probit", "normal")
method <- dplyr::case_when(
  grepl("^cull", responses) ~ "cull",
  grepl("^salad", responses) ~ "salad",
  grepl("^rhis", responses) ~ "rhis",
  grepl("^manta", responses) ~ "manta",
  TRUE ~ "other"
)
traits <- data.frame(
  method = factor(method),
  target_response = as.integer(responses == "cull_high_density"),
  row.names = responses
)
model <- Hmsc::Hmsc(
  Y = Y, XData = as.data.frame(XData), XFormula = XFormula,
  TrData = traits, TrFormula = ~ method + target_response,
  distr = distribution, XScale = TRUE, YScale = TRUE
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
name <- "enhanced_cots_manta_strong"
saveRDS(model, file.path(cfg$paths$models_dir, paste0(name, "_120m.rds")))
artifact <- list(
  name = name, domain = "enhanced", data = data, XData = XData,
  XFormula = XFormula, pca = pca, responses = responses,
  binary_responses = binary_responses, distribution = distribution,
  support_m = 120L, elapsed_minutes = elapsed,
  active_reef_year_conditioning = TRUE,
  manta_active_definition = "direct COTS OR >=2 counted scars OR legacy common scars",
  manta_response_support = "whole-tow length-weighted footprint",
  manta_weak_evidence_in_primary_fit = FALSE,
  manta_negative_sampling = "all strong positives plus max 3 strict negatives per positive within reef-year; max 10 negatives when no positive"
)
saveRDS(artifact, file.path(cfg$paths$artifacts_dir, paste0(name, "_fit_data_120m.rds")))

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
colnames(prediction_mean) <- responses
target <- which(!is.na(data$cull_high_density))
prediction <- tibble::tibble(
  model = name, validation = "apparent_fit", cell_id = data$cell_id[target],
  reef_key = data$reef_key[target], reef_year_id = data$reef_year_id[target],
  year = data$year[target], truth = data$cull_high_density[target],
  probability = prediction_mean[target, "cull_high_density"]
)
metric <- tibble::tibble(
  model = name, responses = length(responses), rows = nrow(data),
  localized_rows = nrow(localized_rows), manta_rows = nrow(manta_rows),
  manta_positives = sum(manta_rows$manta_strong == 1L),
  target_rows = length(target), positives = sum(prediction$truth),
  overall_auc = safe_auc_v4(prediction$truth, prediction$probability),
  within_reef_year_auc = within_reef_year_auc_v4(prediction),
  probability_sd = stats::sd(prediction$probability), elapsed_minutes = elapsed,
  minimum_beta_ess = min(ess, na.rm = TRUE),
  maximum_psrf_upper = max(psrf[, 2], na.rm = TRUE),
  convergence_pass = min(ess, na.rm = TRUE) >= 400 & max(psrf[, 2], na.rm = TRUE) <= 1.1
)
saveRDS(prediction, file.path(cfg$paths$outputs_dir, paste0(name, "_apparent_predictions_120m.rds")))
utils::write.csv(metric, file.path(cfg$paths$outputs_dir, paste0(name, "_metrics_120m.csv")), row.names = FALSE)
print(metric)
