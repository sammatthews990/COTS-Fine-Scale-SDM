args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v3_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v3_dir, "R", "bootstrap.R"))
ctx <- initialise_v3(script_file)
cfg <- ctx$config
support_m <- 120L
base <- terra::rast(cfg$paths$predictor_stack_rg)
geomorphic <- terra::rast(cfg$paths$geomorphic_raster)
legacy <- terra::rast(cfg$paths$legacy_2026_raster)
model <- readRDS(file.path(cfg$paths$models_dir, "hmsc_core_120m.rds"))
fit_artifact <- readRDS(file.path(cfg$paths$artifacts_dir, "hmsc_core_120m_fit_data.rds"))
reference_data <- readRDS(file.path(
  cfg$paths$v2_artifacts_dir, "hmsc_v2_model_data.rds"
))
target_index <- match("cull_high_density", fit_artifact$responses)
beta <- Hmsc::getPostEstimate(model, parName = "Beta")$mean

reef_names <- c(
  "Eyrie Reef (14-118)", "Batt Reef (16-029)",
  "(Big) Broadhurst Reef (No 1) (18-100a)", "Bait Reef (19-137)"
)
reef_slugs <- c("eyrie", "batt", "broadhurst", "bait")
figure_dir <- file.path(cfg$paths$reports_dir, "figures")
map_dir <- file.path(cfg$paths$outputs_dir, "four_reef_maps_120m")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(map_dir, recursive = TRUE, showWarnings = FALSE)

rank_raster <- function(x) {
  value <- terra::values(x, mat = FALSE)
  keep <- is.finite(value)
  out <- rep(NA_real_, length(value))
  out[keep] <- (rank(value[keep], ties.method = "average") - 1) / max(sum(keep) - 1, 1)
  terra::setValues(x, out)
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

predict_hmsc <- function(data) {
  pcs <- apply_family_pca(data, fit_artifact$pca)
  XData <- dplyr::bind_cols(
    dplyr::select(data, "geomorphic_class", "AhyaD", "Aspat", "Aten",
                  "turbid", "waves_Hs", "waves_Tp"), pcs
  )
  XData$geomorphic_class <- factor(
    XData$geomorphic_class,
    levels = levels(fit_artifact$XData$geomorphic_class)
  )
  X <- stats::model.matrix(model$XFormula, XData)
  missing <- setdiff(model$covNames, colnames(X))
  if (length(missing)) X <- cbind(X, matrix(
    0, nrow(X), length(missing), dimnames = list(NULL, missing)
  ))
  X <- X[, model$covNames, drop = FALSE]
  stats::pnorm(as.numeric(X %*% beta[, target_index]))
}

summaries <- list()
for (i in seq_along(reef_names)) {
  reef <- reef_names[[i]]
  message("Mapping ", reef, " at fitted 120 m support.")
  observed <- reference_data |>
    dplyr::filter(.data$reef_name == reef) |>
    dplyr::distinct(.data$longitude, .data$latitude)
  if (!nrow(observed)) next
  observed_native <- terra::project(
    terra::vect(observed, geom = c("longitude", "latitude"), crs = "EPSG:4326"),
    terra::crs(base)
  )
  xy_observed <- terra::crds(observed_native)
  margin <- 2500
  x_bins <- seq(
    floor((min(xy_observed[, 1]) - margin) / support_m),
    floor((max(xy_observed[, 1]) + margin) / support_m)
  )
  y_bins <- seq(
    floor((min(xy_observed[, 2]) - margin) / support_m),
    floor((max(xy_observed[, 2]) + margin) / support_m)
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
  points <- terra::vect(
    cells, geom = c("projected_x", "projected_y"), crs = terra::crs(base)
  )
  lonlat_points <- terra::project(points, terra::crs(geomorphic))
  geom_crop <- terra::crop(geomorphic, terra::ext(lonlat_points) + 0.005, snap = "out")
  cells$geomorphic_class <- as.character(
    terra::extract(geom_crop, lonlat_points)[[2]]
  )
  terrain <- extract_gaussian_terrain_points(
    cells, base, tile_m = 5000, buffer_m = 900,
    cache_dir = file.path(cfg$paths$tmp_dir, "map_terrain_120m", reef_slugs[[i]]),
    id_field = "cell_id"
  ) |>
    dplyr::select(-dplyr::any_of("terrain_row_id"))
  cells <- dplyr::left_join(cells, terrain, by = "cell_id")
  required <- c(
    "geomorphic_class", "AhyaD", "Aspat", "Aten", "turbid", "waves_Hs", "waves_Tp",
    unlist(terrain_families(names(cells)), use.names = FALSE)
  )
  cells <- cells[stats::complete.cases(cells[, required, drop = FALSE]), , drop = FALSE]
  cells$probability <- predict_hmsc(cells)

  extent <- terra::ext(
    min(x_bins) * support_m, (max(x_bins) + 1) * support_m,
    min(y_bins) * support_m, (max(y_bins) + 1) * support_m
  )
  candidate <- terra::rast(extent, resolution = support_m, crs = terra::crs(base))
  terra::values(candidate) <- NA_real_
  candidate[terra::cellFromXY(candidate, as.matrix(cells[, c("projected_x", "projected_y")]))] <-
    cells$probability
  names(candidate) <- "hmsc_geomorphic_probability"
  candidate <- terra::writeRaster(
    candidate, file.path(map_dir, paste0(reef_slugs[[i]], "_hmsc_120m.tif")),
    overwrite = TRUE, wopt = list(datatype = "FLT4S", gdal = "COMPRESS=DEFLATE")
  )
  legacy_aligned <- align_raster_to_template(legacy, candidate)
  candidate_rank <- rank_raster(candidate)
  legacy_rank <- rank_raster(legacy_aligned)
  geom_codes <- match(cells$geomorphic_class, levels(fit_artifact$XData$geomorphic_class))
  geom_map <- candidate
  terra::values(geom_map) <- NA_real_
  geom_map[terra::cellFromXY(geom_map, as.matrix(cells[, c("projected_x", "projected_y")]))] <-
    geom_codes

  png_file <- file.path(figure_dir, paste0("v3_hmsc_120m_", reef_slugs[[i]], ".png"))
  grDevices::png(png_file, width = 3000, height = 950, res = 180)
  graphics::par(mfrow = c(1, 4), mar = c(2, 2, 3.5, 3.5), oma = c(0, 0, 2.5, 0))
  terra::plot(legacy_rank, col = viridisLite::viridis(100), range = c(0, 1),
              main = "Legacy 2026 percentile", axes = FALSE)
  terra::plot(candidate_rank, col = viridisLite::viridis(100), range = c(0, 1),
              main = "Hmsc 120 m percentile", axes = FALSE)
  terra::plot(candidate, col = viridisLite::magma(100),
              main = "Hmsc raw habitat probability", axes = FALSE)
  terra::plot(geom_map, col = grDevices::hcl.colors(11, "Set 3"),
              main = "Geomorphic class", axes = FALSE, plg = FALSE)
  graphics::mtext(reef, outer = TRUE, font = 2, cex = 1.2)
  grDevices::dev.off()

  valid <- is.finite(terra::values(candidate_rank, mat = FALSE)) &
    is.finite(terra::values(legacy_rank, mat = FALSE))
  summaries[[length(summaries) + 1L]] <- tibble::tibble(
    reef = reef, support_m = support_m, mapped_cells = nrow(cells),
    probability_q05 = unname(stats::quantile(cells$probability, 0.05)),
    probability_median = stats::median(cells$probability),
    probability_q95 = unname(stats::quantile(cells$probability, 0.95)),
    legacy_rank_correlation_context_only = stats::cor(
      terra::values(candidate_rank, mat = FALSE)[valid],
      terra::values(legacy_rank, mat = FALSE)[valid], method = "spearman"
    ),
    figure = file.path("figures", basename(png_file)),
    expert_review_status = "pending"
  )
}
manifest <- dplyr::bind_rows(summaries)
utils::write.csv(manifest, file.path(
  cfg$paths$outputs_dir, "four_reef_map_manifest.csv"
), row.names = FALSE)
print(manifest)
