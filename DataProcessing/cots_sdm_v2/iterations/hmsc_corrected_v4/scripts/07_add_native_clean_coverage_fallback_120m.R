args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v4_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
v2_dir <- normalizePath(file.path(v4_dir, "..", ".."), winslash = "/")
source(file.path(v4_dir, "R", "bootstrap.R"))
source(file.path(v2_dir, "R", "mapping.R"))
ctx <- initialise_v4(script_file)
cfg <- ctx$config
support_m <- 120L
base <- terra::rast(cfg$paths$predictor_stack_rg)
legacy <- terra::rast(cfg$paths$legacy_2026_raster)
reference <- readRDS(file.path(cfg$paths$artifacts_dir, "corrected_cell_units_120m.rds"))

load_component <- function(name) {
  list(
    model = readRDS(file.path(cfg$paths$models_dir, paste0(name, "_120m.rds"))),
    artifact = readRDS(file.path(cfg$paths$artifacts_dir, paste0(name, "_fit_data_120m.rds")))
  )
}
enhanced <- load_component("enhanced_joint")
native <- load_component("clean_native_joint")

predict_component <- function(data, component) {
  model <- component$model
  artifact <- component$artifact
  pc_names <- grep("_pc[0-9]+$", names(artifact$XData), value = TRUE)
  base_names <- setdiff(names(artifact$XData), pc_names)
  raw_names <- unlist(lapply(artifact$pca, `[[`, "columns"), use.names = FALSE)
  required <- unique(c(base_names, raw_names))
  complete <- stats::complete.cases(data[, required, drop = FALSE])
  result <- rep(NA_real_, nrow(data))
  if (!any(complete)) return(result)
  z <- data[complete, , drop = FALSE]
  XData <- dplyr::select(z, dplyr::all_of(base_names))
  if (length(artifact$pca)) {
    XData <- dplyr::bind_cols(XData, apply_family_pca(z, artifact$pca))
  }
  if ("geomorphic_class" %in% names(XData)) {
    XData$geomorphic_class <- factor(
      XData$geomorphic_class, levels = levels(artifact$XData$geomorphic_class)
    )
  }
  valid <- stats::complete.cases(XData)
  X <- stats::model.matrix(artifact$XFormula, XData[valid, , drop = FALSE])
  missing <- setdiff(model$covNames, colnames(X))
  if (length(missing)) X <- cbind(X, matrix(
    0, nrow(X), length(missing), dimnames = list(NULL, missing)
  ))
  X <- X[, model$covNames, drop = FALSE]
  beta <- Hmsc::getPostEstimate(model, parName = "Beta")$mean
  target <- match("cull_high_density", artifact$responses)
  result[which(complete)[valid]] <- stats::pnorm(as.numeric(X %*% beta[, target]))
  result
}

clip <- function(x) pmin(pmax(x, 1e-5), 1 - 1e-5)
enhanced_reference <- predict_component(reference, enhanced)
native_reference <- predict_component(reference, native)
overlap <- is.finite(enhanced_reference) & is.finite(native_reference)
calibration <- stats::lm(y ~ x, data = data.frame(
  y = stats::qnorm(clip(enhanced_reference[overlap])),
  x = stats::qnorm(clip(native_reference[overlap]))
))
calibrate_native <- function(probability) {
  result <- rep(NA_real_, length(probability))
  keep <- is.finite(probability)
  result[keep] <- stats::pnorm(stats::predict(
    calibration, newdata = data.frame(x = stats::qnorm(clip(probability[keep])))
  ))
  result
}

extract_base_support <- function(xy) {
  offsets <- seq(-support_m / 2 + 15, support_m / 2 - 15, by = 30)
  grid <- expand.grid(dx = offsets, dy = offsets)
  n <- nrow(grid)
  sample_xy <- cbind(
    rep(xy[, 1], each = n) + rep(grid$dx, times = nrow(xy)),
    rep(xy[, 2], each = n) + rep(grid$dy, times = nrow(xy))
  )
  extracted <- terra::extract(base, sample_xy)
  values <- as.data.frame(vapply(extracted, function(column) {
    z <- matrix(column, nrow = nrow(xy), byrow = TRUE)
    value <- rowMeans(z, na.rm = TRUE)
    value[!is.finite(value)] <- NA_real_
    value
  }, numeric(nrow(xy))))
  names(values) <- names(extracted)
  values
}

reef_names <- c(
  "Eyrie Reef (14-118)", "Batt Reef (16-029)",
  "(Big) Broadhurst Reef (No 1) (18-100a)", "Bait Reef (19-137)"
)
slugs <- c("eyrie", "batt", "broadhurst", "bait")
map_dir <- file.path(cfg$paths$outputs_dir, "four_reef_maps_120m")
figure_dir <- file.path(cfg$paths$reports_dir, "figures")
summary <- list()

for (i in seq_along(slugs)) {
  slug <- slugs[[i]]
  hybrid <- terra::rast(file.path(map_dir, paste0(slug, "_v4_hybrid_120m.tif")))
  enhanced_raster <- terra::rast(file.path(map_dir, paste0(slug, "_v4_enhanced_joint_120m.tif")))
  source_raster <- terra::rast(file.path(map_dir, paste0(slug, "_v4_prediction_source_120m.tif")))
  xy <- terra::xyFromCell(hybrid, seq_len(terra::ncell(hybrid)))
  cells <- extract_base_support(xy)
  cells$projected_x <- xy[, 1]
  cells$projected_y <- xy[, 2]
  native_probability <- predict_component(cells, native)
  native_calibrated <- calibrate_native(native_probability)
  hybrid_value <- terra::values(hybrid, mat = FALSE)
  source_value <- terra::values(source_raster, mat = FALSE)
  dem_domain <- is.finite(cells$DEM)
  fill <- dem_domain & !is.finite(hybrid_value) & is.finite(native_calibrated)
  hybrid_value[fill] <- native_calibrated[fill]
  source_value[fill] <- 3L
  final <- terra::setValues(hybrid, hybrid_value)
  names(final) <- "v4_enhanced_clean_multiscale_native_hybrid"
  source_final <- terra::setValues(source_raster, source_value)
  names(source_final) <- "prediction_source"
  terra::writeRaster(final, file.path(map_dir, paste0(slug, "_v4_final_hybrid_120m.tif")), overwrite = TRUE)
  terra::writeRaster(source_final, file.path(map_dir, paste0(slug, "_v4_final_source_120m.tif")), overwrite = TRUE)

  legacy_aligned <- align_raster_to_template(legacy, final)
  v2_file <- file.path(
    cfg$paths$v2_outputs_dir, "hmsc_v2_reef_maps",
    sprintf("%02d_%s_hmsc_v2_static_habitat.tif", i, slug)
  )
  v2_aligned <- align_raster_to_template(terra::rast(v2_file), final)
  grDevices::png(
    file.path(figure_dir, paste0("v4_final_coverage_comparison_", slug, ".png")),
    width = 3000, height = 1050, res = 180
  )
  graphics::par(mfrow = c(1, 5), mar = c(1.5, 1.5, 3, 3))
  terra::plot(rank_raster_v4(legacy_aligned), col = viridisLite::viridis(100), range = c(0, 1), main = "Legacy clean + year", axes = FALSE)
  terra::plot(rank_raster_v4(v2_aligned), col = viridisLite::viridis(100), range = c(0, 1), main = "Hmsc v2 static joint", axes = FALSE)
  terra::plot(rank_raster_v4(enhanced_raster), col = viridisLite::viridis(100), range = c(0, 1), main = "v4 enhanced joint", axes = FALSE)
  terra::plot(rank_raster_v4(final), col = viridisLite::viridis(100), range = c(0, 1), main = "v4 final hybrid", axes = FALSE)
  terra::plot(source_final, col = c("grey20", "#E69F00", "#56B4E9", "#009E73"), range = c(0, 3), main = "Source: gap / multi / enhanced / native", axes = FALSE)
  grDevices::dev.off()

  summary[[i]] <- tibble::tibble(
    reef_name = reef_names[[i]], dem_domain_cells = sum(dem_domain),
    enhanced_cells = sum(source_value == 2L, na.rm = TRUE),
    clean_multiscale_cells = sum(source_value == 1L, na.rm = TRUE),
    clean_native_cells = sum(source_value == 3L, na.rm = TRUE),
    uncovered_cells = sum(dem_domain & !is.finite(hybrid_value)),
    coverage_fraction = mean(is.finite(hybrid_value[dem_domain])),
    figure = file.path("figures", paste0("v4_final_coverage_comparison_", slug, ".png")),
    expert_review_status = "pending"
  )
}

utils::write.csv(dplyr::bind_rows(summary), file.path(
  cfg$paths$outputs_dir, "four_reef_final_coverage_manifest_120m.csv"
), row.names = FALSE)
utils::write.csv(data.frame(
  intercept = stats::coef(calibration)[[1]], slope = stats::coef(calibration)[[2]],
  overlap_rows = sum(overlap),
  spearman = stats::cor(enhanced_reference[overlap], native_reference[overlap], method = "spearman")
), file.path(cfg$paths$outputs_dir, "native_fallback_calibration_120m.csv"), row.names = FALSE)
print(dplyr::bind_rows(summary))
