args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v2_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v2_dir, "R", "bootstrap.R"))
ctx <- initialise_v2(script_file)
cfg <- ctx$config

if (!requireNamespace("Hmsc", quietly = TRUE)) stop("Hmsc is required.")

grid_m <- 200
data_cache <- file.path(cfg$paths$artifacts_dir, "hmsc_joint_200m_cell_year.rds")
if (file.exists(data_cache)) {
  data <- readRDS(data_cache)
} else {
  template <- terra::rast(cfg$paths$predictor_stack)
  cull <- readRDS(file.path(cfg$paths$artifacts_dir, "site_year_all.rds"))
  cull_membership <- aligned_grid_membership(cull, template, 1L)
  cull$x <- cull_membership$projected_x
  cull$y <- cull_membership$projected_y
  cull_cells <- cull |>
    dplyr::filter(is.finite(.data$x), is.finite(.data$y)) |>
    dplyr::mutate(
      x_bin = floor(.data$x / grid_m),
      y_bin = floor(.data$y / grid_m),
      cell_id = paste(.data$x_bin, .data$y_bin, sep = "||")
    ) |>
    dplyr::group_by(.data$cell_id, .data$x_bin, .data$y_bin, .data$year) |>
    dplyr::summarise(
      reef_name = dplyr::first(.data$reef_name),
      cull_total_cots = sum(.data$total_cots),
      cull_bottom_time = sum(.data$bottom_time),
      cull_n_dives = sum(.data$n_dives),
      cull_n_sites = dplyr::n_distinct(.data$site_id),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      cull_direct_presence = as.integer(.data$cull_total_cots > 0),
      cull_cpue = .data$cull_total_cots / .data$cull_bottom_time,
      cull_high_density = as.integer(.data$cull_cpue >= cfg$response$cpue_threshold)
    )

  manta <- read_manta_observations(cfg$paths$manta_history_csv) |>
    dplyr::filter(
      .data$valid_coordinates, .data$informative,
      dplyr::between(.data$longitude, 140, 160),
      dplyr::between(.data$latitude, -30, -5)
    )
  manta_membership <- aligned_grid_membership(manta, template, 1L)
  manta$x <- manta_membership$projected_x
  manta$y <- manta_membership$projected_y
  manta_cells <- manta |>
    dplyr::filter(is.finite(.data$x), is.finite(.data$y)) |>
    dplyr::mutate(
      x_bin = floor(.data$x / grid_m),
      y_bin = floor(.data$y / grid_m),
      cell_id = paste(.data$x_bin, .data$y_bin, sep = "||"),
      direct_presence = dplyr::coalesce(.data$cots_count > 0, FALSE),
      scar_presence = .data$scar_code %in% c("p", "c") |
        dplyr::coalesce(.data$scar_count > 0, FALSE)
    ) |>
    dplyr::group_by(.data$cell_id, .data$x_bin, .data$y_bin, .data$year) |>
    dplyr::summarise(
      manta_direct_presence = as.integer(any(.data$direct_presence)),
      manta_scar_presence = as.integer(any(.data$scar_presence)),
      manta_n_tows = dplyr::n(),
      manta_n_dates = dplyr::n_distinct(.data$survey_date),
      manta_total_distance_m = sum(.data$tow_distance_m, na.rm = TRUE),
      manta_strict_negative_fraction = mean(.data$strict_negative),
      manta_hard_coral_fraction = mean(.data$hard_coral_observed),
      .groups = "drop"
    )

  joint <- dplyr::inner_join(
    cull_cells, manta_cells,
    by = c("cell_id", "x_bin", "y_bin", "year")
  )
  cell_coordinates <- joint |>
    dplyr::distinct(.data$cell_id, .data$x_bin, .data$y_bin) |>
    dplyr::mutate(
      projected_x = (.data$x_bin + 0.5) * grid_m,
      projected_y = (.data$y_bin + 0.5) * grid_m
    )
  cell_points <- sf::st_as_sf(
    cell_coordinates,
    coords = c("projected_x", "projected_y"),
    crs = terra::crs(template),
    remove = FALSE
  ) |>
    sf::st_transform(4326)
  lonlat <- sf::st_coordinates(cell_points)
  cell_coordinates$longitude <- lonlat[, 1]
  cell_coordinates$latitude <- lonlat[, 2]
  cell_scales <- attach_multiscale_block_means(
    dplyr::select(cell_coordinates, "cell_id", "longitude", "latitude"),
    template
  )

  data <- joint |>
    dplyr::left_join(cell_scales, by = "cell_id") |>
    dplyr::mutate(
      year_z = as.numeric(scale(.data$year)),
      coral_suitability_mean = rowMeans(
        dplyr::pick("AhyaD_s30m", "Aspat_s30m", "Aten_s30m"),
        na.rm = TRUE
      ),
      log_cull_bottom_time = log1p(.data$cull_bottom_time),
      log_cull_dives = log1p(.data$cull_n_dives),
      log_manta_tows = log1p(.data$manta_n_tows),
      log_manta_distance = log1p(.data$manta_total_distance_m)
    )
  saveRDS(data, data_cache)
}

response_names <- c(
  "cull_direct_presence", "cull_high_density",
  "manta_direct_presence", "manta_scar_presence"
)
covariate_names <- c(
  "year_z", "BPI_s30m", "SLO_s60m", "VRM_s60m", "DEM_s30m",
  "SVF_s120m", "coral_suitability_mean",
  "log_cull_bottom_time", "log_cull_dives",
  "log_manta_tows", "log_manta_distance",
  "manta_hard_coral_fraction"
)
data <- data |>
  dplyr::filter(stats::complete.cases(dplyr::pick(
    dplyr::all_of(c(response_names, covariate_names, "reef_name", "cell_id"))
  ))) |>
  dplyr::arrange(.data$cell_id, .data$year)

Y <- as.matrix(data[, response_names])
XData <- as.data.frame(data[, covariate_names])
cell_levels <- unique(data$cell_id)
reef_levels <- unique(data$reef_name)
studyDesign <- data.frame(
  cell = factor(data$cell_id, levels = cell_levels),
  reef = factor(data$reef_name, levels = reef_levels)
)
cell_level <- Hmsc::HmscRandomLevel(units = factor(cell_levels, levels = cell_levels))
reef_level <- Hmsc::HmscRandomLevel(units = factor(reef_levels, levels = reef_levels))

formula <- stats::as.formula(paste("~", paste(covariate_names, collapse = " + ")))
model <- Hmsc::Hmsc(
  Y = Y,
  XData = XData,
  XFormula = formula,
  studyDesign = studyDesign,
  ranLevels = list(cell = cell_level, reef = reef_level),
  distr = rep("probit", length(response_names))
)

mcmc <- list(chains = 4L, samples = 250L, thin = 5L, transient = 500L)
started <- Sys.time()
model <- Hmsc::sampleMcmc(
  model,
  thin = mcmc$thin,
  samples = mcmc$samples,
  transient = mcmc$transient,
  nChains = mcmc$chains,
  nParallel = mcmc$chains,
  verbose = 0
)
mcmc_minutes <- as.numeric(difftime(Sys.time(), started, units = "mins"))
model_path <- file.path(cfg$paths$models_dir, "hmsc_joint_200m_production_model.rds")
saveRDS(model, model_path)

coda_objects <- Hmsc::convertToCodaObject(model)
summarise_coda <- function(object, family) {
  ess <- coda::effectiveSize(object)
  psrf <- tryCatch(
    coda::gelman.diag(object, multivariate = FALSE)$psrf,
    error = function(e) matrix(NA_real_, nrow = length(ess), ncol = 2)
  )
  tibble::tibble(
    family = family,
    parameter = names(ess),
    effective_sample_size = as.numeric(ess),
    psrf_point = as.numeric(psrf[, 1]),
    psrf_upper = as.numeric(psrf[, 2])
  )
}
convergence <- dplyr::bind_rows(
  summarise_coda(coda_objects$Beta, "Beta"),
  dplyr::bind_rows(lapply(seq_along(coda_objects$Omega), function(i) {
    summarise_coda(coda_objects$Omega[[i]], paste0("Omega_random_level_", i))
  }))
)

partition_path <- file.path(cfg$paths$artifacts_dir, "hmsc_joint_200m_cell_partition.rds")
set.seed(20260901)
partition <- Hmsc::createPartition(model, nfolds = 5, column = "cell")
saveRDS(partition, partition_path)
cv_started <- Sys.time()
cv_predicted <- Hmsc::computePredictedValues(
  model,
  partition = partition,
  expected = TRUE,
  nParallel = mcmc$chains
)
cv_minutes <- as.numeric(difftime(Sys.time(), cv_started, units = "mins"))
saveRDS(cv_predicted, file.path(cfg$paths$outputs_dir, "hmsc_joint_200m_cv_predictions_raw.rds"))

prediction_mean <- if (length(dim(cv_predicted)) == 3L) {
  apply(cv_predicted, c(1, 2), mean)
} else {
  as.matrix(cv_predicted)
}
colnames(prediction_mean) <- response_names
performance <- dplyr::bind_rows(lapply(seq_along(response_names), function(i) {
  binary_probability_metrics(Y[, i], prediction_mean[, i]) |>
    dplyr::mutate(
      response = response_names[i],
      evaluation = "five_fold_held_out_200m_cell"
    )
})) |>
  dplyr::select("response", "evaluation", dplyr::everything())

high_density_predictions <- tibble::tibble(
  model = "hmsc_joint_200m",
  validation = "five_fold_held_out_200m_cell",
  fold = as.integer(partition),
  site_year_id = paste(data$cell_id, data$year, sep = "||"),
  site_id = data$cell_id,
  reef_name = data$reef_name,
  reef_year_id = paste(data$reef_name, data$year, sep = "||"),
  year = data$year,
  truth = data$cull_high_density,
  probability = prediction_mean[, "cull_high_density"],
  longitude = data$longitude,
  latitude = data$latitude
)
decomposition <- suppressWarnings(
  prediction_scale_decomposition(high_density_predictions, group = "reef_year_id")
)
capture <- within_group_capture_metrics(
  high_density_predictions, group = "reef_year_id"
)

beta <- Hmsc::getPostEstimate(model, parName = "Beta")
beta_mean <- as.matrix(beta$mean)
beta_support <- as.matrix(beta$support)
beta_table <- expand.grid(
  covariate = c("(Intercept)", covariate_names),
  response = response_names,
  stringsAsFactors = FALSE
) |>
  tibble::as_tibble() |>
  dplyr::mutate(
    posterior_mean = as.numeric(beta_mean),
    posterior_support_positive = as.numeric(beta_support)
  )

audit <- tibble::tibble(
  package_version = as.character(utils::packageVersion("Hmsc")),
  grid_m = grid_m,
  rows = nrow(data),
  cells = dplyr::n_distinct(data$cell_id),
  reefs = dplyr::n_distinct(data$reef_name),
  years = dplyr::n_distinct(data$year),
  response_count = ncol(Y),
  covariate_count = ncol(XData),
  chains = mcmc$chains,
  posterior_samples_per_chain = mcmc$samples,
  thin = mcmc$thin,
  transient = mcmc$transient,
  mcmc_minutes = mcmc_minutes,
  cross_validation_minutes = cv_minutes,
  minimum_ess = min(convergence$effective_sample_size, na.rm = TRUE),
  maximum_psrf_upper = max(convergence$psrf_upper, na.rm = TRUE),
  site_blocked_folds = 5L
)

tables <- list(
  hmsc_joint_200m_audit = audit,
  hmsc_joint_200m_response_prevalence = tibble::tibble(
    response = response_names,
    positives = colSums(Y),
    prevalence = colMeans(Y)
  ),
  hmsc_joint_200m_convergence = convergence,
  hmsc_joint_200m_beta = beta_table,
  hmsc_joint_200m_cv_performance = performance,
  hmsc_joint_200m_high_density_decomposition = decomposition,
  hmsc_joint_200m_high_density_capture = capture
)
for (name in names(tables)) {
  utils::write.csv(
    tables[[name]],
    file.path(cfg$paths$outputs_dir, paste0(name, ".csv")),
    row.names = FALSE
  )
}
saveRDS(high_density_predictions, file.path(
  cfg$paths$outputs_dir, "hmsc_joint_200m_high_density_predictions.rds"
))

message("Joint production-style Hmsc run complete.")
print(audit)
print(performance)
print(decomposition)
