args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v4_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
v2_dir <- normalizePath(file.path(v4_dir, "..", ".."), winslash = "/")
source(file.path(v4_dir, "R", "bootstrap.R"))
source(file.path(v2_dir, "R", "gaussian_terrain.R"))
source(file.path(v2_dir, "R", "mapping.R"))
ctx <- initialise_v4(script_file)
cfg <- ctx$config
support_m <- 200L

model_name <- "enhanced_cots_200support_pca2"
model <- readRDS(file.path(cfg$paths$models_dir, paste0(model_name, "_200m.rds")))
artifact <- readRDS(file.path(cfg$paths$artifacts_dir, paste0(model_name, "_fit_data_200m.rds")))
beta <- Hmsc::getPostEstimate(model, parName = "Beta")$mean
base <- terra::rast(cfg$paths$predictor_stack_rg)
geomorphic <- terra::rast(cfg$paths$geomorphic_raster)
legacy <- terra::rast(cfg$paths$legacy_2026_raster)
reference <- readRDS(file.path(cfg$paths$artifacts_dir, "corrected_cell_units_120m.rds"))

predict_fixed <- function(data) {
  base_names <- setdiff(names(artifact$XData), grep("_pc[0-9]+$", names(artifact$XData), value = TRUE))
  raw_names <- unlist(lapply(artifact$pca, `[[`, "columns"), use.names = FALSE)
  required <- unique(c(base_names, raw_names))
  complete <- stats::complete.cases(data[, required, drop = FALSE])
  result <- rep(NA_real_, nrow(data))
  if (!any(complete)) return(result)
  z <- data[complete, , drop = FALSE]
  pcs <- apply_family_pca(z, artifact$pca)
  XData <- dplyr::bind_cols(dplyr::select(z, dplyr::all_of(base_names)), pcs)
  XData$geomorphic_class <- factor(
    XData$geomorphic_class, levels = levels(artifact$XData$geomorphic_class)
  )
  valid <- stats::complete.cases(XData)
  if (!any(valid)) return(result)
  X <- stats::model.matrix(artifact$XFormula, XData[valid, , drop = FALSE])
  missing <- setdiff(model$covNames, colnames(X))
  if (length(missing)) {
    X <- cbind(X, matrix(0, nrow(X), length(missing), dimnames = list(NULL, missing)))
  }
  X <- X[, model$covNames, drop = FALSE]
  target <- match("cull_high_density", artifact$responses)
  result[which(complete)[valid]] <- stats::pnorm(as.numeric(X %*% beta[, target]))
  result
}

cell_polygons <- function(data, crs) {
  half <- support_m / 2
  geometry <- sf::st_sfc(lapply(seq_len(nrow(data)), function(i) {
    x <- data$projected_x[[i]]
    y <- data$projected_y[[i]]
    sf::st_polygon(list(matrix(c(
      x - half, y - half, x + half, y - half,
      x + half, y + half, x - half, y + half,
      x - half, y - half
    ), ncol = 2, byrow = TRUE)))
  }), crs = crs)
  sf::st_sf(cell_id = data$cell_id, geometry = geometry)
}

extract_base <- function(cells) {
  chunks <- split(seq_len(nrow(cells)), ceiling(seq_len(nrow(cells)) / 500L))
  dplyr::bind_rows(lapply(chunks, function(rows) {
    polygons <- terra::vect(cell_polygons(cells[rows, , drop = FALSE], terra::crs(base)))
    value <- terra::extract(base, polygons, fun = mean, na.rm = TRUE)
    value <- value[, setdiff(names(value), "ID"), drop = FALSE]
    dplyr::bind_cols(cells[rows, , drop = FALSE], tibble::as_tibble(value))
  }))
}

cor_common <- function(a, b) {
  av <- terra::values(a, mat = FALSE)
  bv <- terra::values(b, mat = FALSE)
  keep <- is.finite(av) & is.finite(bv)
  if (sum(keep) < 3L) return(NA_real_)
  stats::cor(av[keep], bv[keep], method = "spearman")
}

reef_names <- c(
  "Eyrie Reef (14-118)", "Batt Reef (16-029)",
  "(Big) Broadhurst Reef (No 1) (18-100a)", "Bait Reef (19-137)"
)
slugs <- c("eyrie", "batt", "broadhurst", "bait")
figure_dir <- file.path(cfg$paths$reports_dir, "figures")
map_dir <- file.path(cfg$paths$outputs_dir, "four_reef_maps_200m")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(map_dir, recursive = TRUE, showWarnings = FALSE)
summaries <- list()

for (i in seq_along(reef_names)) {
  reef <- reef_names[[i]]
  slug <- slugs[[i]]
  message("Mapping 200 m support candidate: ", reef)
  observed <- reference |>
    dplyr::filter(.data$reef_name == reef) |>
    dplyr::distinct(.data$projected_x, .data$projected_y)
  margin <- 2500
  x_bins <- seq(floor((min(observed$projected_x) - margin) / support_m), floor((max(observed$projected_x) + margin) / support_m))
  y_bins <- seq(floor((min(observed$projected_y) - margin) / support_m), floor((max(observed$projected_y) + margin) / support_m))
  cells <- expand.grid(x_bin = x_bins, y_bin = y_bins) |>
    dplyr::mutate(
      cell_id = paste0("map200||", .data$x_bin, "||", .data$y_bin),
      projected_x = (.data$x_bin + 0.5) * support_m,
      projected_y = (.data$y_bin + 0.5) * support_m
    )
  cells <- extract_base(cells) |>
    dplyr::filter(is.finite(.data$DEM))
  points <- terra::vect(cells, geom = c("projected_x", "projected_y"), crs = terra::crs(base))
  geom_points <- terra::project(points, terra::crs(geomorphic))
  cells$geomorphic_class <- as.character(terra::extract(geomorphic, geom_points)[[2]])
  terrain <- extract_gaussian_terrain_points(
    cells, base, tile_m = 20000, buffer_m = 900,
    cache_dir = file.path(cfg$paths$tmp_dir, "map_terrain_200m", slug),
    id_field = "cell_id"
  ) |>
    dplyr::select(-dplyr::any_of("terrain_row_id"))
  cells <- dplyr::left_join(cells, terrain, by = "cell_id")
  cells$probability <- predict_fixed(cells)

  extent <- terra::ext(
    min(x_bins) * support_m, (max(x_bins) + 1) * support_m,
    min(y_bins) * support_m, (max(y_bins) + 1) * support_m
  )
  template <- terra::rast(extent, resolution = support_m, crs = terra::crs(base))
  make_raster <- function(value, name) {
    result <- template
    terra::values(result) <- NA_real_
    result[terra::cellFromXY(result, as.matrix(cells[, c("projected_x", "projected_y")]))] <- value
    names(result) <- name
    result
  }
  candidate <- make_raster(cells$probability, "enhanced_cots_200m_pca2")
  coverage <- make_raster(as.integer(is.finite(cells$probability)), "enhanced_coverage")
  terra::writeRaster(candidate, file.path(map_dir, paste0(slug, "_enhanced_cots_pca2_200m.tif")), overwrite = TRUE)
  terra::writeRaster(coverage, file.path(map_dir, paste0(slug, "_enhanced_coverage_200m.tif")), overwrite = TRUE)

  legacy_aligned <- align_raster_to_template(legacy, candidate)
  v2_file <- file.path(cfg$paths$v2_outputs_dir, "hmsc_v2_reef_maps", sprintf("%02d_%s_hmsc_v2_static_habitat.tif", i, slug))
  v2_aligned <- align_raster_to_template(terra::rast(v2_file), candidate)
  v4_120_file <- file.path(cfg$paths$outputs_dir, "four_reef_maps_120m", paste0(slug, "_v4_final_hybrid_120m.tif"))
  v4_120_aligned <- align_raster_to_template(terra::rast(v4_120_file), candidate)

  grDevices::png(
    file.path(figure_dir, paste0("v4_support_comparison_", slug, "_120m_200m.png")),
    width = 3200, height = 1800, res = 180
  )
  graphics::par(mfrow = c(2, 3), mar = c(1.5, 1.5, 3, 3))
  terra::plot(rank_raster_v4(legacy_aligned), col = viridisLite::viridis(100), range = c(0, 1), main = "Legacy clean + year percentile", axes = FALSE)
  terra::plot(rank_raster_v4(v2_aligned), col = viridisLite::viridis(100), range = c(0, 1), main = "Hmsc v2 static-joint percentile", axes = FALSE)
  terra::plot(rank_raster_v4(v4_120_aligned), col = viridisLite::viridis(100), range = c(0, 1), main = "v4 final hybrid 120 m", axes = FALSE)
  terra::plot(rank_raster_v4(candidate), col = viridisLite::viridis(100), range = c(0, 1), main = "v4 enhanced COTS PCA2 200 m", axes = FALSE)
  terra::plot(candidate, col = viridisLite::magma(100), main = "200 m raw probability", axes = FALSE)
  terra::plot(coverage, col = c("grey20", "#56B4E9"), range = c(0, 1), main = "200 m enhanced coverage", axes = FALSE)
  graphics::mtext(reef, outer = TRUE, line = -1, cex = 1.25, font = 2)
  grDevices::dev.off()

  summaries[[length(summaries) + 1L]] <- tibble::tibble(
    reef_name = reef, support_m = support_m,
    dem_domain_cells = nrow(cells), enhanced_cells = sum(is.finite(cells$probability)),
    coverage_fraction = mean(is.finite(cells$probability)),
    probability_sd = stats::sd(cells$probability, na.rm = TRUE),
    spearman_legacy = cor_common(candidate, legacy_aligned),
    spearman_hmsc_v2 = cor_common(candidate, v2_aligned),
    spearman_v4_120m_hybrid = cor_common(candidate, v4_120_aligned),
    expert_review_status = "pending",
    figure = file.path("figures", paste0("v4_support_comparison_", slug, "_120m_200m.png"))
  )
}

manifest <- dplyr::bind_rows(summaries)
utils::write.csv(
  manifest,
  file.path(cfg$paths$outputs_dir, "four_reef_support_comparison_manifest_200m.csv"),
  row.names = FALSE
)
print(manifest)
