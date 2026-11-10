args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v4_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
v2_dir <- normalizePath(file.path(v4_dir, "..", ".."), winslash = "/")
source(file.path(v4_dir, "R", "bootstrap.R"))
source(file.path(v2_dir, "R", "gaussian_terrain.R"))
source(file.path(v2_dir, "R", "mapping.R"))
ctx <- initialise_v4(script_file)
cfg <- ctx$config
support_m <- 120L
base <- terra::rast(cfg$paths$predictor_stack_rg)
geomorphic <- terra::rast(cfg$paths$geomorphic_raster)
legacy <- terra::rast(cfg$paths$legacy_2026_raster)

model_names <- c("enhanced_cots", "enhanced_cots_manta_strong")
models <- stats::setNames(lapply(model_names, function(name) {
  readRDS(file.path(cfg$paths$models_dir, paste0(name, "_120m.rds")))
}), model_names)
artifacts <- stats::setNames(lapply(model_names, function(name) {
  readRDS(file.path(cfg$paths$artifacts_dir, paste0(name, "_fit_data_120m.rds")))
}), model_names)
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
  XData$geomorphic_class <- factor(
    XData$geomorphic_class, levels = levels(artifact$XData$geomorphic_class)
  )
  valid <- stats::complete.cases(XData)
  if (!any(valid)) return(result)
  X <- stats::model.matrix(artifact$XFormula, XData[valid, , drop = FALSE])
  missing <- setdiff(model$covNames, colnames(X))
  if (length(missing)) {
    X <- cbind(X, matrix(0, nrow(X), length(missing), dimnames = list(NULL, missing))))))
  }
  X <- X[, model$covNames, drop = FALSE]
  target <- match("cull_high_density", artifact$responses)
  result[which(complete)[valid]] <- stats::pnorm(as.numeric(X %*% betas[[name]][, target]))
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
  values <- as.data.frame(terra::extract(base, xy))
  as.data.frame(vapply(values, function(column) {
    matrix_value <- matrix(column, nrow = nrow(cells), byrow = TRUE)
    result <- rowMeans(matrix_value, na.rm = TRUE)
    result[!is.finite(result)] <- NA_real_
    result
  }, numeric(nrow(cells))))
}

cor_common <- function(a, b) {
  keep <- is.finite(a) & is.finite(b)
  if (sum(keep) < 3L) return(NA_real_)
  stats::cor(a[keep], b[keep], method = "spearman")
}

reef_names <- c(
  "Eyrie Reef (14-118)", "Batt Reef (16-029)",
  "(Big) Broadhurst Reef (No 1) (18-100a)", "Bait Reef (19-137)"
)
slugs <- c("eyrie", "batt", "broadhurst", "bait")
figure_dir <- file.path(cfg$paths$reports_dir, "figures")
map_dir <- file.path(cfg$paths$outputs_dir, "four_reef_maps_120m")
summaries <- list()

for (i in seq_along(slugs)) {
  slug <- slugs[[i]]
  template <- terra::rast(file.path(map_dir, paste0(slug, "_v4_final_hybrid_120m.tif")))
  xy <- terra::xyFromCell(template, seq_len(terra::ncell(template)))
  x_bin <- floor(xy[, 1] / support_m)
  y_bin <- floor(xy[, 2] / support_m)
  cells <- tibble::tibble(
    cell_id = paste0("map120||", x_bin, "||", y_bin),
    projected_x = xy[, 1], projected_y = xy[, 2]
  )
  cells <- dplyr::bind_cols(cells, tibble::as_tibble(extract_base_support(cells))) |>
    dplyr::filter(is.finite(.data$DEM))
  points <- terra::vect(cells, geom = c("projected_x", "projected_y"), crs = terra::crs(base))
  geom_points <- terra::project(points, terra::crs(geomorphic))
  cells$geomorphic_class <- as.character(terra::extract(geomorphic, geom_points)[[2]])
  terrain <- extract_gaussian_terrain_points(
    cells, base, tile_m = 5000, buffer_m = 900,
    cache_dir = file.path(cfg$paths$v3_dir, "tmp", "map_terrain_120m", slug),
    id_field = "cell_id"
  ) |>
    dplyr::select(-dplyr::any_of("terrain_row_id"))
  cells <- dplyr::left_join(cells, terrain, by = "cell_id")
  cells$cots_only <- predict_fixed(cells, "enhanced_cots")
  cells$manta_augmented <- predict_fixed(cells, "enhanced_cots_manta_strong")

  make_raster <- function(value, name) {
    result <- template
    terra::values(result) <- NA_real_
    target <- terra::cellFromXY(result, as.matrix(cells[, c("projected_x", "projected_y")]))
    result[target] <- value
    names(result) <- name
    result
  }
  cots_raster <- make_raster(cells$cots_only, "enhanced_cots_only")
  manta_raster <- make_raster(cells$manta_augmented, "enhanced_cots_manta_strong")
  terra::writeRaster(manta_raster, file.path(map_dir, paste0(slug, "_v4_enhanced_cots_manta_strong_120m.tif")), overwrite = TRUE)

  legacy_aligned <- align_raster_to_template(legacy, template)
  v2_file <- file.path(cfg$paths$v2_outputs_dir, "hmsc_v2_reef_maps", sprintf("%02d_%s_hmsc_v2_static_habitat.tif", i, slug))
  v2_aligned <- align_raster_to_template(terra::rast(v2_file), template)
  final_hybrid <- terra::rast(file.path(map_dir, paste0(slug, "_v4_final_hybrid_120m.tif")))
  cots_rank <- terra::values(rank_raster_v4(cots_raster), mat = FALSE)
  manta_rank <- terra::values(rank_raster_v4(manta_raster), mat = FALSE)
  rank_change <- make_raster(manta_rank[terra::cellFromXY(template, as.matrix(cells[, c("projected_x", "projected_y")]))] - cots_rank[terra::cellFromXY(template, as.matrix(cells[, c("projected_x", "projected_y")]))], "rank_change")

  grDevices::png(
    file.path(figure_dir, paste0("v4_manta_hmsc_comparison_", slug, ".png")),
    width = 3200, height = 1800, res = 180
  )
  graphics::par(mfrow = c(2, 3), mar = c(1.5, 1.5, 3, 3))
  terra::plot(rank_raster_v4(legacy_aligned), col = viridisLite::viridis(100), range = c(0, 1), main = "Legacy clean + year percentile", axes = FALSE)
  terra::plot(rank_raster_v4(v2_aligned), col = viridisLite::viridis(100), range = c(0, 1), main = "Hmsc v2 static-joint percentile", axes = FALSE)
  terra::plot(rank_raster_v4(final_hybrid), col = viridisLite::viridis(100), range = c(0, 1), main = "v4 final hybrid + fallback", axes = FALSE)
  terra::plot(rank_raster_v4(cots_raster), col = viridisLite::viridis(100), range = c(0, 1), main = "v4 COTS-only Hmsc", axes = FALSE)
  terra::plot(rank_raster_v4(manta_raster), col = viridisLite::viridis(100), range = c(0, 1), main = "v4 + strong Manta response", axes = FALSE)
  terra::plot(rank_change, col = grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(101), range = c(-0.1, 0.1), main = "Manta rank change (-0.1 to 0.1)", axes = FALSE)
  graphics::mtext(reef_names[[i]], outer = TRUE, line = -1, cex = 1.25, font = 2)
  grDevices::dev.off()

  summaries[[i]] <- tibble::tibble(
    reef_name = reef_names[[i]],
    common_cells = sum(is.finite(cells$cots_only) & is.finite(cells$manta_augmented)),
    cots_manta_spearman = cor_common(cells$cots_only, cells$manta_augmented),
    mean_absolute_probability_change = mean(abs(cells$cots_only - cells$manta_augmented), na.rm = TRUE),
    mean_absolute_rank_change = mean(abs(terra::values(rank_change, mat = FALSE)), na.rm = TRUE),
    expert_review_status = "pending",
    figure = file.path("figures", paste0("v4_manta_hmsc_comparison_", slug, ".png"))
  )
}

manifest <- dplyr::bind_rows(summaries)
utils::write.csv(manifest, file.path(cfg$paths$outputs_dir, "four_reef_manta_hmsc_comparison_manifest_120m.csv"), row.names = FALSE)
print(manifest)
