args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v2_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v2_dir, "R", "bootstrap.R"))
ctx <- initialise_v2(script_file)
cfg <- ctx$config

if (!requireNamespace("Hmsc", quietly = TRUE)) {
  stop("Hmsc is required for this preliminary fit.")
}

site_year <- readRDS(file.path(cfg$paths$artifacts_dir, "site_year_all.rds"))
site_scales <- readRDS(file.path(
  cfg$paths$artifacts_dir, "site_multiscale_30_60_120m.rds"
))

data <- site_year |>
  dplyr::select(
    "site_year_id", "site_id", "reef_name", "year", "total_cots",
    "high_density", "longitude", "latitude"
  ) |>
  dplyr::left_join(site_scales, by = c("site_id", "longitude", "latitude")) |>
  dplyr::transmute(
    site_year_id, site_id, reef_name, year,
    cots_presence = as.integer(.data$total_cots > 0),
    high_density = as.integer(.data$high_density),
    year_z = as.numeric(scale(.data$year)),
    BPI_s30m = .data$BPI_s30m,
    SLO_s60m = .data$SLO_s60m,
    VRM_s60m = .data$VRM_s60m,
    DEM_s30m = .data$DEM_s30m,
    SVF_s120m = .data$SVF_s120m,
    coral_suitability_mean = rowMeans(
      dplyr::pick("AhyaD_s30m", "Aspat_s30m", "Aten_s30m"),
      na.rm = TRUE
    )
  ) |>
  dplyr::filter(stats::complete.cases(dplyr::pick(dplyr::everything())))

set.seed(20260831)
if (nrow(data) > 3000L) {
  sample_rows <- sample.int(nrow(data), 3000L, replace = FALSE)
  data <- data[sample_rows, , drop = FALSE]
}

Y <- as.matrix(data[, c("cots_presence", "high_density")])
XData <- data.frame(
  year_z = data$year_z,
  BPI_s30m = data$BPI_s30m,
  SLO_s60m = data$SLO_s60m,
  VRM_s60m = data$VRM_s60m,
  DEM_s30m = data$DEM_s30m,
  SVF_s120m = data$SVF_s120m,
  coral_suitability_mean = data$coral_suitability_mean
)
studyDesign <- data.frame(reef = factor(data$reef_name))
reef_level <- Hmsc::HmscRandomLevel(units = studyDesign$reef)

model <- Hmsc::Hmsc(
  Y = Y,
  XData = XData,
  XFormula = ~ year_z + BPI_s30m + SLO_s60m + VRM_s60m +
    DEM_s30m + SVF_s120m + coral_suitability_mean,
  studyDesign = studyDesign,
  ranLevels = list(reef = reef_level),
  distr = c("probit", "probit")
)

started <- Sys.time()
model <- Hmsc::sampleMcmc(
  model,
  thin = 2,
  samples = 100,
  transient = 200,
  nChains = 2,
  nParallel = 2,
  verbose = 0
)
elapsed_minutes <- as.numeric(difftime(Sys.time(), started, units = "mins"))

coda_objects <- Hmsc::convertToCodaObject(model)
beta_ess <- coda::effectiveSize(coda_objects$Beta)
beta_psrf <- tryCatch(
  coda::gelman.diag(coda_objects$Beta, multivariate = FALSE)$psrf,
  error = function(e) matrix(NA_real_, nrow = length(beta_ess), ncol = 2)
)
convergence <- tibble::tibble(
  parameter = names(beta_ess),
  effective_sample_size = as.numeric(beta_ess),
  psrf_point = as.numeric(beta_psrf[, 1]),
  psrf_upper = as.numeric(beta_psrf[, 2])
)

predicted <- Hmsc::computePredictedValues(model, expected = TRUE, nParallel = 2)
predicted_mean <- if (length(dim(predicted)) == 3L) {
  apply(predicted, c(1, 2), mean)
} else {
  as.matrix(predicted)
}
colnames(predicted_mean) <- colnames(Y)
performance <- dplyr::bind_rows(lapply(seq_len(ncol(Y)), function(i) {
  binary_probability_metrics(Y[, i], predicted_mean[, i]) |>
    dplyr::mutate(response = colnames(Y)[i], evaluation = "explanatory_not_cross_validated")
})) |>
  dplyr::select("response", "evaluation", dplyr::everything())

beta <- Hmsc::getPostEstimate(model, parName = "Beta")
beta_mean <- as.matrix(beta$mean)
beta_support <- as.matrix(beta$support)
beta_table <- expand.grid(
  covariate = c("(Intercept)", colnames(XData)),
  response = colnames(Y),
  stringsAsFactors = FALSE
) |>
  tibble::as_tibble() |>
  dplyr::mutate(
    posterior_mean = as.numeric(beta_mean),
    posterior_support_positive = as.numeric(beta_support)
  )

audit <- tibble::tibble(
  package_version = as.character(utils::packageVersion("Hmsc")),
  rows = nrow(data),
  reefs = dplyr::n_distinct(data$reef_name),
  site_years = dplyr::n_distinct(data$site_year_id),
  responses = ncol(Y),
  covariates = ncol(XData),
  chains = 2L,
  posterior_samples_per_chain = 100L,
  thin = 2L,
  transient = 200L,
  elapsed_minutes = elapsed_minutes,
  minimum_beta_ess = min(convergence$effective_sample_size, na.rm = TRUE),
  maximum_beta_psrf_upper = max(convergence$psrf_upper, na.rm = TRUE),
  preliminary_only = TRUE
)

utils::write.csv(convergence, file.path(cfg$paths$outputs_dir, "hmsc_preliminary_convergence.csv"), row.names = FALSE)
utils::write.csv(performance, file.path(cfg$paths$outputs_dir, "hmsc_preliminary_performance.csv"), row.names = FALSE)
utils::write.csv(beta_table, file.path(cfg$paths$outputs_dir, "hmsc_preliminary_beta.csv"), row.names = FALSE)
utils::write.csv(audit, file.path(cfg$paths$outputs_dir, "hmsc_preliminary_audit.csv"), row.names = FALSE)
saveRDS(model, file.path(cfg$paths$models_dir, "hmsc_preliminary_model.rds"))
saveRDS(data, file.path(cfg$paths$artifacts_dir, "hmsc_preliminary_data.rds"))

message("Preliminary Hmsc fit complete.")
print(audit)
print(performance)
