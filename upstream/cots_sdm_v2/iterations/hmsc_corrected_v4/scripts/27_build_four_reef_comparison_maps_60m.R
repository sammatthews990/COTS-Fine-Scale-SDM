args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v4_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
v2_dir <- normalizePath(file.path(v4_dir, "..", ".."), winslash = "/")
source(file.path(v4_dir, "R", "bootstrap.R"))
source(file.path(v2_dir, "R", "gaussian_terrain.R"))
source(file.path(v2_dir, "R", "mapping.R"))
ctx <- initialise_v4(script_file)
cfg <- ctx$config
support_m <- 60L
model_name <- "enhanced_cots_60support_pca1"
model <- readRDS(file.path(cfg$paths$models_dir, paste0(model_name, "_60m.rds")))
artifact <- readRDS(file.path(cfg$paths$artifacts_dir, paste0(model_name, "_fit_data_60m.rds")))
beta <- Hmsc::getPostEstimate(model, parName = "Beta")$mean
base <- terra::rast(cfg$paths$predictor_stack_rg)
geomorphic <- terra::rast(cfg$paths$geomorphic_raster)
legacy <- terra::rast(cfg$paths$legacy_2026_raster)

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
  XData$geomorphic_class <- factor(XData$geomorphic_class, levels = levels(artifact$XData$geomorphic_class))
  valid <- stats::complete.cases(XData)
  if (!any(valid)) return(result)
  X <- stats::model.matrix(artifact$XFormula, XData[valid, , drop = FALSE])
  missing <- setdiff(model$covNames, colnames(X))
  if (length(missing)) X <- cbind(X, matrix(0, nrow(X), length(missing), dimnames = list(NULL, missing)))
  X <- X[, model$covNames, drop = FALSE]
  target <- match("cull_high_density", artifact$responses)
  result[which(complete)[valid]] <- stats::pnorm(as.numeric(X %*% beta[, target]))
  result
}

extract_base_support <- function(cells) {
  offsets <- c(-15, 15)
  grid <- expand.grid(dx = offsets, dy = offsets)
  n <- nrow(grid)
  xy <- cbind(
    rep(cells$projected_x, each = n) + rep(grid$dx, times = nrow(cells)),
    rep(cells$projected_y, each = n) + rep(grid$dy, times = nrow(cells))
  )
  values <- as.data.frame(terra::extract(base, xy))
  as.data.frame(vapply(values, function(column) {
    z <- matrix(column, nrow = nrow(cells), byrow = TRUE)
    value <- rowMeans(z, na.rm = TRUE)
    value[!is.finite(value)] <- NA_real_
    value
  }, numeric(nrow(cells))))
}

cor_raster <- function(a, b) {
  av <- terra::values(a, mat = FALSE); bv <- terra::values(b, mat = FALSE)
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
map_dir <- file.path(cfg$paths$outputs_dir, "four_reef_maps_60m")
dir.create(map_dir, recursive = TRUE, showWarnings = FALSE)
summaries <- list()

for (i in seq_along(slugs)) {
  slug <- slugs[[i]]
  baseline <- terra::rast(file.path(cfg$paths$outputs_dir, "four_reef_maps_120m", paste0(slug, "_v4_final_hybrid_120m.tif")))
  extent <- terra::ext(baseline)
  template <- terra::rast(extent, resolution = support_m, crs = terra::crs(base))
  xy <- terra::xyFromCell(template, seq_len(terra::ncell(template)))
  x_bin <- floor(xy[, 1] / support_m); y_bin <- floor(xy[, 2] / support_m)
  cells <- tibble::tibble(
    cell_id = paste0("map60||", x_bin, "||", y_bin),
    projected_x = xy[, 1], projected_y = xy[, 2]
  )
  cells <- dplyr::bind_cols(cells, tibble::as_tibble(extract_base_support(cells))) |>
    dplyr::filter(is.finite(.data$DEM))
  points <- terra::vect(cells, geom = c("projected_x", "projected_y"), crs = terra::crs(base))
  geom_points <- terra::project(points, terra::crs(geomorphic))
  cells$geomorphic_class <- as.character(terra::extract(geomorphic, geom_points)[[2]])
  terrain <- extract_gaussian_terrain_points(
    cells, base, tile_m = 20000, buffer_m = 900,
    cache_dir = file.path(cfg$paths$tmp_dir, "map_terrain_60m", slug),
    id_field = "cell_id"
  ) |>
    dplyr::select(-dplyr::any_of("terrain_row_id"))
  cells <- dplyr::left_join(cells, terrain, by = "cell_id")
  cells$probability <- predict_fixed(cells)
  candidate <- template; terra::values(candidate) <- NA_real_
  target <- terra::cellFromXY(candidate, as.matrix(cells[, c("projected_x", "projected_y")]))
  candidate[target] <- cells$probability
  names(candidate) <- "enhanced_cots_60m_pca1"
  terra::writeRaster(candidate, file.path(map_dir, paste0(slug, "_enhanced_cots_pca1_60m.tif")), overwrite = TRUE)

  legacy_aligned <- align_raster_to_template(legacy, candidate)
  v2_file <- file.path(cfg$paths$v2_outputs_dir, "hmsc_v2_reef_maps", sprintf("%02d_%s_hmsc_v2_static_habitat.tif", i, slug))
  v2_aligned <- align_raster_to_template(terra::rast(v2_file), candidate)
  baseline_aligned <- align_raster_to_template(baseline, candidate)
  candidate_rank <- rank_raster_v4(candidate)
  baseline_rank <- rank_raster_v4(baseline_aligned)
  rank_delta <- candidate_rank - baseline_rank

  grDevices::png(
    file.path(figure_dir, paste0("v4_support_comparison_", slug, "_60m_120m.png")),
    width = 3200, height = 1800, res = 180
  )
  graphics::par(mfrow = c(2, 3), mar = c(1.5, 1.5, 3, 3))
  terra::plot(rank_raster_v4(legacy_aligned), col = viridisLite::viridis(100), range = c(0, 1), main = "Legacy clean + year percentile", axes = FALSE)
  terra::plot(rank_raster_v4(v2_aligned), col = viridisLite::viridis(100), range = c(0, 1), main = "Hmsc v2 static-joint percentile", axes = FALSE)
  terra::plot(baseline_rank, col = viridisLite::viridis(100), range = c(0, 1), main = "v4 final hybrid 120 m", axes = FALSE)
  terra::plot(candidate_rank, col = viridisLite::viridis(100), range = c(0, 1), main = "v4 enhanced COTS 60 m", axes = FALSE)
  terra::plot(candidate, col = viridisLite::magma(100), main = "60 m raw probability", axes = FALSE)
  terra::plot(rank_delta, col = grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(101), range = c(-0.5, 0.5), main = "60 m minus 120 m rank", axes = FALSE)
  graphics::mtext(reef_names[[i]], outer = TRUE, line = -1, cex = 1.25, font = 2)
  grDevices::dev.off()

  summaries[[i]] <- tibble::tibble(
    reef_name = reef_names[[i]], support_m = support_m,
    dem_domain_cells = nrow(cells), enhanced_cells = sum(is.finite(cells$probability)),
    coverage_fraction = mean(is.finite(cells$probability)),
    probability_sd = stats::sd(cells$probability, na.rm = TRUE),
    spearman_legacy = cor_raster(candidate, legacy_aligned),
    spearman_hmsc_v2 = cor_raster(candidate, v2_aligned),
    spearman_v4_120m_hybrid = cor_raster(candidate, baseline_aligned),
    expert_review_status = "pending",
    figure = file.path("figures", paste0("v4_support_comparison_", slug, "_60m_120m.png"))
  )
}
manifest <- dplyr::bind_rows(summaries)
utils::write.csv(manifest, file.path(cfg$paths$outputs_dir, "four_reef_support_comparison_manifest_60m.csv"), row.names = FALSE)
print(manifest)
