args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v2_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v2_dir, "R", "bootstrap.R"))
ctx <- initialise_v2(script_file)
cfg <- ctx$config

if (!requireNamespace("Hmsc", quietly = TRUE)) stop("Hmsc is required.")
if (!requireNamespace("viridisLite", quietly = TRUE)) {
  stop("viridisLite is required.")
}
candidate_arg <- grep("^--candidate=", commandArgs(trailingOnly = TRUE), value = TRUE)
candidate <- if (length(candidate_arg)) {
  sub("^--candidate=", "", candidate_arg[[1]])
} else {
  "full_multiscale"
}
reduced_candidate <- identical(candidate, "geomorphic_g360")
artifact_prefix <- if (reduced_candidate) {
  "hmsc_v2_geomorphic_g360"
} else {
  "hmsc_v2"
}
candidate_title <- if (reduced_candidate) {
  "Hmsc v2 geomorphic 360 m"
} else {
  "Hmsc v2 static joint habitat"
}
required <- c(
  cfg$paths$predictor_stack, cfg$paths$predictor_stack_rg,
  cfg$paths$geomorphic_raster, cfg$paths$legacy_2026_raster,
  file.path(
    cfg$paths$models_dir,
    if (reduced_candidate) "hmsc_v2_geomorphic_g360_core_model.rds" else
      "hmsc_v2_ecological_core_model.rds"
  ),
  file.path(
    cfg$paths$artifacts_dir,
    if (reduced_candidate) "hmsc_v2_geomorphic_g360_fit_data.rds" else
      "hmsc_v2_model_data.rds"
  )
)
if (any(!file.exists(required))) stop("Missing Hmsc v2 map input(s).")

template <- terra::rast(cfg$paths$predictor_stack)
legacy_2026 <- terra::rast(cfg$paths$legacy_2026_raster)
model <- readRDS(required[[5]])
model_data <- readRDS(required[[6]])
target_index <- match("cull_high_density", model$spNames)
if (is.na(target_index)) stop("Cull high-density response is absent.")
beta <- Hmsc::getPostEstimate(model, parName = "Beta")$mean

reef_names <- c(
  "Eyrie Reef (14-118)",
  "Batt Reef (16-029)",
  "(Big) Broadhurst Reef (No 1) (18-100a)",
  "Bait Reef (19-137)"
)
reef_slugs <- c("eyrie", "batt", "broadhurst", "bait")
fig_dir <- file.path(v2_dir, "reports", "figures")
map_dir <- file.path(
  cfg$paths$outputs_dir, paste0(artifact_prefix, "_reef_maps")
)
cache_dir <- file.path(v2_dir, "tmp", "hmsc_v2_map_gaussian_tiles")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(map_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

rank_raster <- function(x) {
  value <- terra::values(x, mat = FALSE)
  valid <- is.finite(value)
  ranked <- rep(NA_real_, length(value))
  ranked[valid] <- (rank(value[valid], ties.method = "average") - 1) /
    max(1, sum(valid) - 1)
  terra::setValues(x, ranked)
}

predict_fixed_habitat <- function(data) {
  data <- data |>
    dplyr::mutate(
      coral_suitability_mean = rowMeans(
        dplyr::pick("AhyaD", "Aspat", "Aten"), na.rm = TRUE
      )
    )
  feature_names <- names(model$XData)
  complete <- stats::complete.cases(data[, feature_names, drop = FALSE])
  result <- rep(NA_real_, nrow(data))
  if (!any(complete)) return(result)
  X <- stats::model.matrix(
    model$XFormula, data = data[complete, , drop = FALSE]
  )
  X <- X[, model$covNames, drop = FALSE]
  eta <- as.numeric(X %*% beta[, target_index])
  result[complete] <- stats::pnorm(eta)
  result
}

summaries <- list()
for (i in seq_along(reef_names)) {
  reef <- reef_names[[i]]
  reef_sites <- model_data |>
    dplyr::filter(.data$reef_name == reef) |>
    dplyr::distinct(.data$cell_id, .data$longitude, .data$latitude)
  if (nrow(reef_sites) == 0L) {
    warning("No model-domain observations found for ", reef)
    next
  }
  bbox_wgs <- terra::ext(
    min(reef_sites$longitude) - 0.025,
    max(reef_sites$longitude) + 0.025,
    min(reef_sites$latitude) - 0.025,
    max(reef_sites$latitude) + 0.025
  )
  bbox_native <- terra::project(
    bbox_wgs, from = "EPSG:4326", to = terra::crs(template)
  )
  native <- terra::crop(template[[1]], bbox_native, snap = "out")
  cells <- terra::as.data.frame(
    native, cells = TRUE, xy = TRUE, na.rm = TRUE
  ) |>
    dplyr::transmute(
      raster_cell = .data$cell,
      cell_id = paste0("map_", reef_slugs[[i]], "_", .data$cell),
      projected_x = .data$x,
      projected_y = .data$y
    )
  point_vector <- terra::vect(
    cells, geom = c("projected_x", "projected_y"), crs = terra::crs(template)
  ) |>
    terra::project("EPSG:4326")
  lonlat <- terra::crds(point_vector)
  cells$longitude <- lonlat[, 1]
  cells$latitude <- lonlat[, 2]

  site_predictors <- extract_site_habitat_predictors(
    dplyr::select(cells, "cell_id", "longitude", "latitude"),
    cfg$paths$predictor_stack, cfg$paths$predictor_stack_rg,
    cfg$paths$geomorphic_raster
  ) |>
    dplyr::filter(.data$clean_complete == 1L, .data$rg_complete == 1L)
  map_cells <- cells |>
    dplyr::semi_join(site_predictors, by = "cell_id")
  terrain <- extract_gaussian_terrain_points(
    map_cells, template, tile_m = 5000, buffer_m = 900,
    cache_dir = file.path(cache_dir, reef_slugs[[i]]), id_field = "cell_id"
  ) |>
    dplyr::select(-"terrain_row_id")
  prediction_data <- map_cells |>
    dplyr::inner_join(site_predictors, by = "cell_id") |>
    dplyr::inner_join(terrain, by = "cell_id")
  prediction_data$probability <- predict_fixed_habitat(prediction_data)

  probability <- native
  terra::values(probability) <- NA_real_
  cell_position <- match(prediction_data$raster_cell, cells$raster_cell)
  native_position <- cells$raster_cell[cell_position]
  probability[native_position] <- prediction_data$probability
  names(probability) <- "hmsc_v2_static_habitat_probability"
  probability_file <- file.path(
    map_dir, sprintf("%02d_%s_%s_static_habitat.tif", i, reef_slugs[[i]], artifact_prefix)
  )
  probability <- terra::writeRaster(
    probability, probability_file, overwrite = TRUE,
    wopt = list(datatype = "FLT4S", gdal = c("COMPRESS=DEFLATE"))
  )

  legacy <- align_raster_to_template(legacy_2026, probability)
  incumbent_rank <- rank_raster(legacy)
  candidate_rank <- rank_raster(probability)
  difference <- candidate_rank - incumbent_rank
  names(incumbent_rank) <- "incumbent_percentile"
  names(candidate_rank) <- "hmsc_v2_percentile"
  names(difference) <- "hmsc_v2_minus_incumbent"

  site_points <- terra::project(
    terra::vect(
      reef_sites, geom = c("longitude", "latitude"), crs = "EPSG:4326"
    ), terra::crs(template)
  )
  png_name <- sprintf(
    "%s_vs_incumbent_%02d_%s.png", artifact_prefix, i, reef_slugs[[i]]
  )
  grDevices::png(
    file.path(fig_dir, png_name), width = 2550, height = 950, res = 180
  )
  graphics::par(
    mfrow = c(1, 3), mar = c(2.2, 2.2, 4.0, 4.2), oma = c(0, 0, 3, 0)
  )
  terra::plot(
    incumbent_rank, col = viridisLite::viridis(100), range = c(0, 1),
    main = "Current best: clean + year (2026)", axes = FALSE,
    plg = list(title = "Within-map percentile", cex = 0.72)
  )
  terra::points(site_points, pch = 3, cex = 0.30, col = "white")
  graphics::box()
  terra::plot(
    candidate_rank, col = viridisLite::viridis(100), range = c(0, 1),
    main = candidate_title, axes = FALSE,
    plg = list(title = "Within-map percentile", cex = 0.72)
  )
  terra::points(site_points, pch = 3, cex = 0.30, col = "white")
  graphics::box()
  terra::plot(
    difference,
    col = grDevices::hcl.colors(101, "Blue-Red 3", rev = TRUE),
    range = c(-1, 1), main = "Hmsc v2 minus incumbent rank", axes = FALSE,
    plg = list(title = "Percentile difference", cex = 0.72)
  )
  terra::points(site_points, pch = 3, cex = 0.30, col = "black")
  graphics::box()
  graphics::mtext(reef, outer = TRUE, cex = 1.28, font = 2)
  grDevices::dev.off()

  incumbent_value <- terra::values(incumbent_rank, mat = FALSE)
  candidate_value <- terra::values(candidate_rank, mat = FALSE)
  probability_value <- terra::values(probability, mat = FALSE)
  valid <- is.finite(incumbent_value) & is.finite(candidate_value)
  probability_valid <- probability_value[is.finite(probability_value)]
  summaries[[length(summaries) + 1L]] <- tibble::tibble(
    reef_name = reef,
    figure = file.path("figures", png_name),
    valid_cells = sum(valid),
    probability_q05 = unname(stats::quantile(probability_valid, 0.05)),
    probability_median = stats::median(probability_valid),
    probability_q95 = unname(stats::quantile(probability_valid, 0.95)),
    spearman_rank_correlation = stats::cor(
      incumbent_value[valid], candidate_value[valid], method = "spearman"
    ),
    mean_absolute_percentile_difference = mean(abs(
      candidate_value[valid] - incumbent_value[valid]
    )),
    review_status = "pending_expert_visual_acceptance"
  )
}

manifest <- dplyr::bind_rows(summaries)
utils::write.csv(
  manifest,
  file.path(
    cfg$paths$outputs_dir, paste0(artifact_prefix, "_visual_map_manifest.csv")
  ),
  row.names = FALSE
)
message("Hmsc v2 reef comparison maps complete.")
print(manifest)
