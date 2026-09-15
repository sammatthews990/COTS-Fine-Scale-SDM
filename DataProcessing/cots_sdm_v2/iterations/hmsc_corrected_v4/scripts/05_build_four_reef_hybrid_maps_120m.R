args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v4_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v4_dir, "R", "bootstrap.R"))
ctx <- initialise_v4(script_file)
cfg <- ctx$config
support_m <- 120L
base <- terra::rast(cfg$paths$predictor_stack_rg)
geomorphic <- terra::rast(cfg$paths$geomorphic_raster)
legacy <- terra::rast(cfg$paths$legacy_2026_raster)
reference_data <- readRDS(file.path(cfg$paths$artifacts_dir, "corrected_cell_units_120m.rds"))

model_names <- c("enhanced_cots", "enhanced_joint", "clean_joint")
models <- stats::setNames(lapply(model_names, function(name) readRDS(file.path(
  cfg$paths$models_dir, paste0(name, "_120m.rds")
))), model_names)
artifacts <- stats::setNames(lapply(model_names, function(name) readRDS(file.path(
  cfg$paths$artifacts_dir, paste0(name, "_fit_data_120m.rds")
))), model_names)
betas <- lapply(models, function(model) Hmsc::getPostEstimate(model, parName = "Beta")$mean)

predict_fixed <- function(data, name) {
  model <- models[[name]]
  artifact <- artifacts[[name]]
  base_names <- setdiff(names(artifact$XData), grep("_pc[0-9]+$", names(artifact$XData), value = TRUE))
  raw_names <- unlist(lapply(artifact$pca, `[[`, "columns"), use.names = FALSE)
  required <- unique(c(base_names, raw_names))
  complete <- stats::complete.cases(data[, required, drop = FALSE])
  result <- rep(NA_real_, nrow(data))
  if (!any(complete)) return(result)
  z <- data[complete, , drop = FALSE]
  pcs <- apply_family_pca(z, artifact$pca)
  XData <- dplyr::bind_cols(dplyr::select(z, dplyr::all_of(base_names)), pcs)
  if ("geomorphic_class" %in% names(XData)) {
    XData$geomorphic_class <- factor(
      XData$geomorphic_class,
      levels = levels(artifact$XData$geomorphic_class)
    )
  }
  valid_level <- stats::complete.cases(XData)
  if (!any(valid_level)) return(result)
  X <- stats::model.matrix(artifact$XFormula, XData[valid_level, , drop = FALSE])
  missing <- setdiff(model$covNames, colnames(X))
  if (length(missing)) {
    X <- cbind(X, matrix(0, nrow(X), length(missing), dimnames = list(NULL, missing)))
  }
  X <- X[, model$covNames, drop = FALSE]
  target <- match("cull_high_density", artifact$responses)
  position <- which(complete)[valid_level]
  result[position] <- stats::pnorm(as.numeric(X %*% betas[[name]][, target]))
  result
}

# Calibrate only the clean model's probability scale to the enhanced model on
# their common training domain. This is a seam correction, not a legacy-rank
# constraint and it does not change the clean model's ordering.
calibration_data <- reference_data
enhanced_overlap <- predict_fixed(calibration_data, "enhanced_joint")
clean_overlap <- predict_fixed(calibration_data, "clean_joint")
calibration_keep <- is.finite(enhanced_overlap) & is.finite(clean_overlap)
clip <- function(x) pmin(pmax(x, 1e-5), 1 - 1e-5)
fallback_calibration <- stats::lm(
  stats::qnorm(clip(enhanced_overlap[calibration_keep])) ~
    stats::qnorm(clip(clean_overlap[calibration_keep]))
)
calibrate_clean <- function(probability) {
  result <- rep(NA_real_, length(probability))
  keep <- is.finite(probability)
  result[keep] <- stats::pnorm(stats::predict(
    fallback_calibration,
    newdata = data.frame(
      `stats::qnorm(clip(clean_overlap[calibration_keep]))` = stats::qnorm(clip(probability[keep])),
      check.names = FALSE
    )
  ))
  result
}

extract_base_support <- function(cells) {
  offsets <- seq(-support_m / 2 + 15, support_m / 2 - 15, by = 30)
  grid <- expand.grid(dx = offsets, dy = offsets)
  n <- nrow(grid)
  xy <- cbind(
    rep(cells$projected_x, each = n) + rep(grid$dx, times = nrow(cells)),
    rep(cells$projected_y, each = n) + rep(grid$dy, times = nrow(cells))
  )
  values <- terra::extract(base, xy)[, -1, drop = FALSE]
  as.data.frame(vapply(values, function(column) {
    z <- matrix(column, nrow = nrow(cells), byrow = TRUE)
    value <- rowMeans(z, na.rm = TRUE)
    value[!is.finite(value)] <- NA_real_
    value
  }, numeric(nrow(cells))))
}

reef_names <- c(
  "Eyrie Reef (14-118)", "Batt Reef (16-029)",
  "(Big) Broadhurst Reef (No 1) (18-100a)", "Bait Reef (19-137)"
)
reef_slugs <- c("eyrie", "batt", "broadhurst", "bait")
figure_dir <- file.path(cfg$paths$reports_dir, "figures")
map_dir <- file.path(cfg$paths$outputs_dir, "four_reef_maps_120m")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(map_dir, recursive = TRUE, showWarnings = FALSE)
summaries <- list()

for (i in seq_along(reef_names)) {
  reef <- reef_names[[i]]
  slug <- reef_slugs[[i]]
  message("Mapping corrected hybrid at 120 m: ", reef)
  observed <- reference_data |>
    dplyr::filter(.data$reef_name == reef) |>
    dplyr::distinct(.data$projected_x, .data$projected_y)
  if (!nrow(observed)) next
  margin <- 2500
  x_bins <- seq(
    floor((min(observed$projected_x) - margin) / support_m),
    floor((max(observed$projected_x) + margin) / support_m)
  )
  y_bins <- seq(
    floor((min(observed$projected_y) - margin) / support_m),
    floor((max(observed$projected_y) + margin) / support_m)
  )
  cells <- expand.grid(x_bin = x_bins, y_bin = y_bins) |>
    dplyr::mutate(
      cell_id = paste0("map120||", .data$x_bin, "||", .data$y_bin),
      projected_x = (.data$x_bin + 0.5) * support_m,
      projected_y = (.data$y_bin + 0.5) * support_m
    )
  base_values <- extract_base_support(cells)
  cells <- dplyr::bind_cols(cells, tibble::as_tibble(base_values)) |>
    dplyr::filter(is.finite(.data$DEM))
  points <- terra::vect(cells, geom = c("projected_x", "projected_y"), crs = terra::crs(base))
  lonlat_points <- terra::project(points, terra::crs(geomorphic))
  geom_crop <- terra::crop(geomorphic, terra::ext(lonlat_points) + 0.005, snap = "out")
  cells$geomorphic_class <- as.character(terra::extract(geom_crop, lonlat_points)[[2]])
  terrain <- extract_gaussian_terrain_points(
    cells, base, tile_m = 5000, buffer_m = 900,
    cache_dir = file.path(cfg$paths$v3_dir, "tmp", "map_terrain_120m", slug),
    id_field = "cell_id"
  ) |>
    dplyr::select(-dplyr::any_of("terrain_row_id"))
  cells <- dplyr::left_join(cells, terrain, by = "cell_id")

  cells$enhanced_cots <- predict_fixed(cells, "enhanced_cots")
  cells$enhanced_joint <- predict_fixed(cells, "enhanced_joint")
  cells$clean_joint <- predict_fixed(cells, "clean_joint")
  cells$clean_calibrated <- calibrate_clean(cells$clean_joint)
  cells$hybrid_joint <- cells$enhanced_joint
  fallback <- !is.finite(cells$hybrid_joint) & is.finite(cells$clean_calibrated)
  cells$hybrid_joint[fallback] <- cells$clean_calibrated[fallback]
  cells$prediction_source <- dplyr::case_when(
    is.finite(cells$enhanced_joint) ~ "enhanced",
    fallback ~ "clean_fallback",
    TRUE ~ "uncovered"
  )

  extent <- terra::ext(
    min(x_bins) * support_m, (max(x_bins) + 1) * support_m,
    min(y_bins) * support_m, (max(y_bins) + 1) * support_m
  )
  template <- terra::rast(extent, resolution = support_m, crs = terra::crs(base))
  make_raster <- function(values, name) {
    x <- template
    terra::values(x) <- NA_real_
    x[terra::cellFromXY(x, as.matrix(cells[, c("projected_x", "projected_y")]))] <- values
    names(x) <- name
    x
  }
  enhanced_cots <- make_raster(cells$enhanced_cots, "enhanced_cots")
  enhanced_joint <- make_raster(cells$enhanced_joint, "enhanced_joint")
  hybrid <- make_raster(cells$hybrid_joint, "hybrid_joint")
  source_value <- match(cells$prediction_source, c("uncovered", "clean_fallback", "enhanced")) - 1L
  source_raster <- make_raster(source_value, "prediction_source")
  terra::writeRaster(hybrid, file.path(map_dir, paste0(slug, "_v4_hybrid_120m.tif")), overwrite = TRUE)
  terra::writeRaster(enhanced_joint, file.path(map_dir, paste0(slug, "_v4_enhanced_joint_120m.tif")), overwrite = TRUE)
  terra::writeRaster(source_raster, file.path(map_dir, paste0(slug, "_v4_prediction_source_120m.tif")), overwrite = TRUE)

  legacy_aligned <- align_raster_to_template(legacy, template)
  v2_file <- file.path(
    cfg$paths$v2_outputs_dir, "hmsc_v2_reef_maps",
    sprintf("%02d_%s_hmsc_v2_static_habitat.tif", i, slug)
  )
  v2_aligned <- align_raster_to_template(terra::rast(v2_file), template)
  grDevices::png(
    file.path(figure_dir, paste0("v4_hybrid_comparison_", slug, ".png")),
    width = 3200, height = 1800, res = 180
  )
  graphics::par(mfrow = c(2, 3), mar = c(1.5, 1.5, 3, 3))
  terra::plot(rank_raster_v4(legacy_aligned), col = viridisLite::viridis(100), range = c(0, 1), main = "Legacy clean + year percentile", axes = FALSE)
  terra::plot(rank_raster_v4(v2_aligned), col = viridisLite::viridis(100), range = c(0, 1), main = "Hmsc v2 static-joint percentile", axes = FALSE)
  terra::plot(rank_raster_v4(enhanced_cots), col = viridisLite::viridis(100), range = c(0, 1), main = "v4 enhanced COTS-only percentile", axes = FALSE)
  terra::plot(rank_raster_v4(enhanced_joint), col = viridisLite::viridis(100), range = c(0, 1), main = "v4 enhanced joint percentile", axes = FALSE)
  terra::plot(rank_raster_v4(hybrid), col = viridisLite::viridis(100), range = c(0, 1), main = "v4 enhanced + clean fallback", axes = FALSE)
  terra::plot(source_raster, col = c("grey20", "#E69F00", "#56B4E9"), range = c(0, 2), main = "Coverage: uncovered / clean / enhanced", axes = FALSE)
  graphics::mtext(reef, outer = TRUE, line = -1, cex = 1.25, font = 2)
  grDevices::dev.off()

  summaries[[length(summaries) + 1L]] <- tibble::tibble(
    reef_name = reef, support_m = support_m, dem_domain_cells = nrow(cells),
    enhanced_cells = sum(cells$prediction_source == "enhanced"),
    clean_fallback_cells = sum(cells$prediction_source == "clean_fallback"),
    uncovered_cells = sum(cells$prediction_source == "uncovered"),
    coverage_fraction = mean(cells$prediction_source != "uncovered"),
    hybrid_q05 = stats::quantile(cells$hybrid_joint, 0.05, na.rm = TRUE),
    hybrid_median = stats::median(cells$hybrid_joint, na.rm = TRUE),
    hybrid_q95 = stats::quantile(cells$hybrid_joint, 0.95, na.rm = TRUE),
    figure = file.path("figures", paste0("v4_hybrid_comparison_", slug, ".png")),
    expert_review_status = "pending"
  )
}

utils::write.csv(dplyr::bind_rows(summaries), file.path(
  cfg$paths$outputs_dir, "four_reef_hybrid_map_manifest_120m.csv"
), row.names = FALSE)
utils::write.csv(data.frame(
  intercept = stats::coef(fallback_calibration)[[1]],
  slope = stats::coef(fallback_calibration)[[2]],
  overlap_rows = sum(calibration_keep),
  spearman = stats::cor(
    enhanced_overlap[calibration_keep], clean_overlap[calibration_keep], method = "spearman"
  )
), file.path(cfg$paths$outputs_dir, "clean_fallback_calibration_120m.csv"), row.names = FALSE)
print(dplyr::bind_rows(summaries))
