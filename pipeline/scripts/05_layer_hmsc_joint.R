# Layer 3: Hmsc joint habitat model, 120 m, full domain.
#
# Static habitat score: given an active COTS reef-year, the probability that a
# 120 m cell holds problematic-density COTS (Cull CPUE >= 0.02), from a joint
# model of Cull, SALAD and RHIS evidence. Reef, year, outbreak probability,
# effort and method are excluded by design.
#
# Each cell takes the first component that can predict it:
#   1 enhanced_joint      geomorphic class + ReefGuide exposure + Acropora + multiscale terrain
#   2 clean_joint         Acropora + multiscale terrain
#   3 clean_native_joint  clean 30 m stack (block-mean) + Acropora
# Fallbacks are calibrated onto the enhanced probability scale (probit-linear,
# fitted on training cells all components can predict) and never replace an
# available enhanced prediction. Calibration is applied to every posterior
# draw, so the uncertainty bands are on the same scale.

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "..", "R", "common.R"))
cfg <- initialise_final()
target_name <- cfg$hmsc$target_response

components <- stats::setNames(lapply(cfg$hmsc$components, function(name) {
  model_path <- file.path(cfg$dirs$models, paste0("hmsc_", name, "_production.rds"))
  summary_path <- file.path(cfg$dirs$tables, paste0("hmsc_", name, "_production_summary.csv"))
  if (!file.exists(model_path)) stop("Run 03_fit_hmsc_production.R --component=", name, " first.")
  summary <- utils::read.csv(summary_path)
  if (!isTRUE(summary$convergence_pass)) {
    warning(name, " does not pass the convergence gate (min ESS ", round(summary$min_ess),
            ", max upper PSRF ", round(summary$max_psrf_upper, 3), ").")
  }
  model <- readRDS(model_path)
  artifact <- readRDS(checked_input(cfg, paste0("hmsc_", name)))
  target <- match(target_name, artifact$responses)
  list(
    model = model, artifact = artifact,
    beta = target_beta_draws(model, target, cfg$hmsc$posterior_draws_for_maps, cfg$hmsc$seed)
  )
}), cfg$hmsc$components)

posterior_mean <- function(data, component) {
  design <- component_design(data, component$artifact, component$model)
  out <- rep(NA_real_, nrow(data))
  if (length(design$rows)) {
    out[design$rows] <- summarise_target_posterior(design$X, component$beta)[, "mean"]
  }
  out
}

# Seam calibration on the training reference cells.
reference <- readRDS(checked_input(cfg, "hmsc_reference_cells"))
reference_pred <- lapply(components, posterior_mean, data = reference)
calibrations <- list(
  clean_joint = fit_probit_calibration(reference_pred$enhanced_joint, reference_pred$clean_joint),
  clean_native_joint = fit_probit_calibration(reference_pred$enhanced_joint, reference_pred$clean_native_joint)
)
write_table(
  dplyr::bind_rows(lapply(names(calibrations), function(n) {
    tibble::as_tibble(calibrations[[n]]) |> dplyr::mutate(fallback = n, .before = 1)
  })),
  cfg, "hmsc_fallback_calibration.csv"
)

# Full-domain prediction.
grid <- readRDS(file.path(cfg$dirs$cache, "support_grid_120m.rds"))
log_message("Predicting ", format(nrow(grid), big.mark = ","), " cells")
result <- matrix(NA_real_, nrow(grid), 4, dimnames = list(NULL, c("mean", "sd", "q05", "q95")))
source_code <- integer(nrow(grid))
for (k in seq_along(components)) {
  name <- names(components)[k]
  remaining <- which(source_code == 0L)
  if (!length(remaining)) break
  design <- component_design(grid[remaining, , drop = FALSE], components[[name]]$artifact, components[[name]]$model)
  if (!length(design$rows)) next
  transform <- if (name %in% names(calibrations)) {
    calibration <- calibrations[[name]]
    function(p) apply_probit_calibration(p, calibration)
  } else identity
  rows <- remaining[design$rows]
  result[rows, ] <- summarise_target_posterior(design$X, components[[name]]$beta, transform = transform)
  source_code[rows] <- k
  log_message(name, ": ", format(length(rows), big.mark = ","), " cells")
}

# Rasterise onto the 120 m support grid.
template <- terra::rast(file.path(cfg$dirs$cache, "base_support_120m.tif"))[[1]]
to_raster <- function(values, name) {
  x <- terra::rast(template)
  x[grid$cell] <- values
  names(x) <- name
  x
}
probability <- to_raster(result[, "mean"], "cots_prob")
uncertainty <- c(
  to_raster(result[, "sd"], "posterior_sd"),
  to_raster(result[, "q05"], "posterior_q05"),
  to_raster(result[, "q95"], "posterior_q95")
)
source_layer <- to_raster(source_code, "prediction_source")
probability <- terra::trim(probability)
uncertainty <- terra::crop(uncertainty, probability)
source_layer <- terra::crop(source_layer, probability)

layer_path <- file.path(cfg$dirs$layers, "cots_hmsc_joint_120m.tif")
uncertainty_path <- file.path(cfg$dirs$layers, "cots_hmsc_joint_uncertainty_120m.tif")
source_path <- file.path(cfg$dirs$layers, "cots_hmsc_joint_source_120m.tif")
write_layer(probability, layer_path)
write_layer(uncertainty, uncertainty_path)
write_layer(source_layer, source_path, datatype = "INT1U")

cells_by_source <- vapply(c(1:3, 0L), function(k) sum(source_code == k), numeric(1))
coverage <- tibble::tibble(
  source_code = c(1:3, 0L),
  component = c(cfg$hmsc$components, "uncovered"),
  cells = cells_by_source
) |>
  dplyr::mutate(fraction = .data$cells / sum(.data$cells))
write_table(coverage, cfg, "hmsc_joint_coverage.csv")

convergence <- dplyr::bind_rows(lapply(cfg$hmsc$components, function(name) {
  utils::read.csv(file.path(cfg$dirs$tables, paste0("hmsc_", name, "_production_summary.csv")))
}))
register_layer(
  cfg, "hmsc_joint", layer_path,
  description = "Posterior mean P(problematic-density COTS | active reef-year); static habitat, 120 m support",
  model = "Hmsc v4 joint (enhanced -> clean multiscale -> clean native fallback), production chains",
  units = "probability",
  extra = list(
    enhanced_fraction = coverage$fraction[1],
    all_components_converged = all(convergence$convergence_pass)
  )
)
register_layer(
  cfg, "hmsc_joint_uncertainty", uncertainty_path,
  description = "Bands: posterior SD, 5th and 95th percentiles of P (400 posterior draws)",
  model = "Hmsc v4 joint", units = "probability"
)
register_layer(
  cfg, "hmsc_joint_source", source_path,
  description = "Component used per cell: 1 enhanced_joint, 2 clean_joint, 3 clean_native_joint",
  model = "Hmsc v4 joint", units = "category"
)
print(coverage)
