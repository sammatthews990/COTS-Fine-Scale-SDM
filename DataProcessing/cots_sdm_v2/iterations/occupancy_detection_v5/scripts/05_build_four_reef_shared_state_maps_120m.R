args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v5_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
repo_root <- normalizePath(file.path(v5_dir, "..", "..", "..", ".."), winslash = "/")
v2_dir <- file.path(repo_root, "DataProcessing", "cots_sdm_v2")
v4_dir <- file.path(v2_dir, "iterations", "hmsc_corrected_v4")
source(file.path(v4_dir, "R", "bootstrap.R"))
source(file.path(v2_dir, "R", "gaussian_terrain.R"))
source(file.path(v2_dir, "R", "mapping.R"))
ctx <- initialise_v4(file.path(v4_dir, "scripts", "01a_prepare_corrected_cells_120m.R"))
cfg <- ctx$config
support_m <- 120L
base <- terra::rast(cfg$paths$predictor_stack_rg)
geomorphic <- terra::rast(cfg$paths$geomorphic_raster)
legacy <- terra::rast(cfg$paths$legacy_2026_raster)
artifact <- readRDS(file.path(v5_dir, "artifacts", "shared_state_data_120m.rds"))
model <- readRDS(file.path(v5_dir, "models", "shared_state_track_refined_120m.rds"))
beta <- model$parameter_table |>
  dplyr::filter(.data$component == "state") |>
  dplyr::pull(.data$estimate)

predict_state <- function(data) {
  raw_names <- unlist(lapply(artifact$pca, `[[`, "columns"), use.names = FALSE)
  required <- unique(c(artifact$state_base_terms, raw_names))
  complete <- stats::complete.cases(data[, required, drop = FALSE])
  result <- rep(NA_real_, nrow(data))
  if (!any(complete)) return(result)
  z <- data[complete, , drop = FALSE]
  pcs <- apply_family_pca(z, artifact$pca)
  state_data <- dplyr::bind_cols(dplyr::select(z, dplyr::all_of(artifact$state_base_terms)), pcs)
  state_data$geomorphic_class <- factor(
    state_data$geomorphic_class,
    levels = levels(factor(artifact$sites$geomorphic_class))
  )
  X <- stats::model.matrix(
    ~ geomorphic_class + AhyaD + Aspat + Aten + turbid + waves_Hs + waves_Tp + .,
    data = dplyr::select(state_data, -"geomorphic_class") |>
      dplyr::mutate(geomorphic_class = state_data$geomorphic_class, .before = 1)
  )
  missing <- setdiff(artifact$design_columns, colnames(X))
  if (length(missing)) X <- cbind(X, matrix(0, nrow(X), length(missing), dimnames = list(NULL, missing)))
  X <- X[, artifact$design_columns, drop = FALSE]
  X <- sweep(X, 2, artifact$design_center, "-")
  X <- sweep(X, 2, artifact$design_scale, "/")
  valid <- stats::complete.cases(X)
  result[which(complete)[valid]] <- stats::plogis(as.numeric(X[valid, , drop = FALSE] %*% beta))
  result
}

extract_base_support <- function(cells) {
  offsets <- c(-45, -15, 15, 45)
  grid <- expand.grid(dx = offsets, dy = offsets); n <- nrow(grid)
  xy <- cbind(
    rep(cells$projected_x, each = n) + rep(grid$dx, times = nrow(cells)),
    rep(cells$projected_y, each = n) + rep(grid$dy, times = nrow(cells))
  )
  values <- as.data.frame(terra::extract(base, xy))
  as.data.frame(vapply(values, function(column) {
    matrix_value <- matrix(column, nrow = nrow(cells), byrow = TRUE)
    value <- rowMeans(matrix_value, na.rm = TRUE); value[!is.finite(value)] <- NA_real_; value
  }, numeric(nrow(cells))))
}

cor_raster <- function(a, b) {
  av <- terra::values(a, mat = FALSE); bv <- terra::values(b, mat = FALSE)
  keep <- is.finite(av) & is.finite(bv)
  if (sum(keep) < 3L) return(NA_real_)
  stats::cor(av[keep], bv[keep], method = "spearman")
}

reef_names <- c("Eyrie Reef (14-118)", "Batt Reef (16-029)", "(Big) Broadhurst Reef (No 1) (18-100a)", "Bait Reef (19-137)")
slugs <- c("eyrie", "batt", "broadhurst", "bait")
map_dir <- file.path(v5_dir, "outputs", "four_reef_maps_120m")
figure_dir <- file.path(v5_dir, "reports", "figures")
dir.create(map_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
summaries <- list()

for (i in seq_along(slugs)) {
  slug <- slugs[[i]]
  baseline <- terra::rast(file.path(v4_dir, "outputs", "four_reef_maps_120m", paste0(slug, "_v4_final_hybrid_120m.tif")))
  xy <- terra::xyFromCell(baseline, seq_len(terra::ncell(baseline)))
  cells <- tibble::tibble(
    cell_id = paste0("map120||", floor(xy[, 1] / support_m), "||", floor(xy[, 2] / support_m)),
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
  cells$probability <- predict_state(cells)
  candidate <- baseline; terra::values(candidate) <- NA_real_
  target <- terra::cellFromXY(candidate, as.matrix(cells[, c("projected_x", "projected_y")]))
  candidate[target] <- cells$probability; names(candidate) <- "shared_state_probability"
  terra::writeRaster(candidate, file.path(map_dir, paste0(slug, "_shared_state_120m.tif")), overwrite = TRUE)
  legacy_aligned <- align_raster_to_template(legacy, candidate)
  v2_file <- file.path(cfg$paths$v2_outputs_dir, "hmsc_v2_reef_maps", sprintf("%02d_%s_hmsc_v2_static_habitat.tif", i, slug))
  v2_aligned <- align_raster_to_template(terra::rast(v2_file), candidate)
  candidate_rank <- rank_raster_v4(candidate); baseline_rank <- rank_raster_v4(baseline)
  rank_delta <- candidate_rank - baseline_rank

  grDevices::png(file.path(figure_dir, paste0("shared_state_comparison_", slug, ".png")), width = 3200, height = 1800, res = 180)
  graphics::par(mfrow = c(2, 3), mar = c(1.5, 1.5, 3, 3))
  terra::plot(rank_raster_v4(legacy_aligned), col = viridisLite::viridis(100), range = c(0, 1), main = "Legacy clean + year percentile", axes = FALSE)
  terra::plot(rank_raster_v4(v2_aligned), col = viridisLite::viridis(100), range = c(0, 1), main = "Hmsc v2 static-joint percentile", axes = FALSE)
  terra::plot(baseline_rank, col = viridisLite::viridis(100), range = c(0, 1), main = "v4 final hybrid 120 m", axes = FALSE)
  terra::plot(candidate_rank, col = viridisLite::viridis(100), range = c(0, 1), main = "v5 shared-state 120 m", axes = FALSE)
  terra::plot(candidate, col = viridisLite::magma(100), main = "Shared-state raw probability", axes = FALSE)
  terra::plot(rank_delta, col = grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(101), range = c(-0.5, 0.5), main = "Shared-state minus v4 rank", axes = FALSE)
  graphics::mtext(reef_names[[i]], outer = TRUE, line = -1, cex = 1.25, font = 2)
  grDevices::dev.off()
  summaries[[i]] <- tibble::tibble(
    reef_name = reef_names[[i]], enhanced_cells = sum(is.finite(cells$probability)),
    coverage_fraction = mean(is.finite(cells$probability)), probability_sd = stats::sd(cells$probability, na.rm = TRUE),
    spearman_legacy = cor_raster(candidate, legacy_aligned), spearman_hmsc_v2 = cor_raster(candidate, v2_aligned),
    spearman_v4_120m_hybrid = cor_raster(candidate, baseline), expert_review_status = "pending",
    figure = file.path("figures", paste0("shared_state_comparison_", slug, ".png"))
  )
}
manifest <- dplyr::bind_rows(summaries)
utils::write.csv(manifest, file.path(v5_dir, "outputs", "four_reef_shared_state_manifest_120m.csv"), row.names = FALSE)
print(manifest)
