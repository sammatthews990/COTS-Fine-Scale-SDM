args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v6_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
v2_dir <- normalizePath(file.path(v6_dir, "..", ".."), winslash = "/")
v4_dir <- normalizePath(file.path(v6_dir, "..", "hmsc_corrected_v4"), winslash = "/")
source(file.path(v6_dir, "R", "bootstrap.R"))
source(file.path(v2_dir, "R", "gaussian_terrain.R"))
source(file.path(v2_dir, "R", "mapping.R"))
ctx <- initialise_density_v6(script_file)
cfg <- ctx$cfg

support_m <- 60L
artifact <- readRDS(file.path(ctx$paths$artifacts, "density_data_60m.rds"))
models <- list(
  cull_native = readRDS(file.path(ctx$paths$models, "cull_native_60m.rds")),
  salad_native = readRDS(file.path(ctx$paths$models, "salad_native_60m.rds")),
  joint_balanced = readRDS(file.path(ctx$paths$models, "joint_rosetta_balanced_60m.rds"))
)
base <- terra::rast(cfg$paths$predictor_stack_rg)
geomorphic <- terra::rast(cfg$paths$geomorphic_raster)
legacy <- terra::rast(cfg$paths$legacy_2026_raster)

build_density_X <- function(data) {
  design <- artifact$design
  complete <- stats::complete.cases(data[, unique(c(
    design$base_names, design$raw_names
  )), drop = FALSE])
  result <- matrix(NA_real_, nrow(data), length(design$columns),
                   dimnames = list(NULL, design$columns))
  if (!any(complete)) return(result)
  z <- data[complete, , drop = FALSE]
  pcs <- apply_family_pca(z, design$fit_artifact$pca)
  xdata <- dplyr::bind_cols(
    dplyr::select(z, dplyr::all_of(design$base_names)), pcs
  )
  xdata$geomorphic_class <- factor(
    xdata$geomorphic_class,
    levels = levels(design$fit_artifact$XData$geomorphic_class)
  )
  valid <- stats::complete.cases(xdata)
  if (!any(valid)) return(result)
  X <- stats::model.matrix(
    design$fit_artifact$XFormula, xdata[valid, , drop = FALSE]
  )
  missing <- setdiff(design$columns, colnames(X))
  if (length(missing)) {
    X <- cbind(X, matrix(0, nrow(X), length(missing),
                         dimnames = list(NULL, missing)))
  }
  X <- X[, design$columns, drop = FALSE]
  X <- sweep(X, 2, design$centre, "-")
  X <- sweep(X, 2, design$scale, "/")
  result[which(complete)[valid], ] <- X
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
figure_dir <- file.path(v6_dir, "reports", "figures")
map_dir <- file.path(ctx$paths$outputs, "four_reef_density_maps_60m")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(map_dir, recursive = TRUE, showWarnings = FALSE)
summaries <- list()

for (i in seq_along(slugs)) {
  slug <- slugs[[i]]
  baseline_file <- file.path(
    v4_dir, "outputs", "four_reef_maps_120m",
    paste0(slug, "_v4_final_hybrid_120m.tif")
  )
  baseline <- terra::rast(baseline_file)
  template <- terra::rast(
    terra::ext(baseline), resolution = support_m, crs = terra::crs(base)
  )
  xy <- terra::xyFromCell(template, seq_len(terra::ncell(template)))
  cells <- tibble::tibble(
    cell_id = paste0("density60||", floor(xy[, 1] / 60), "||", floor(xy[, 2] / 60)),
    projected_x = xy[, 1], projected_y = xy[, 2]
  )
  cells <- dplyr::bind_cols(
    cells, tibble::as_tibble(extract_base_support(cells))
  ) |> dplyr::filter(is.finite(.data$DEM))
  points <- terra::vect(
    cells, geom = c("projected_x", "projected_y"), crs = terra::crs(base)
  )
  geom_points <- terra::project(points, terra::crs(geomorphic))
  cells$geomorphic_class <- as.character(
    terra::extract(geomorphic, geom_points)[[2]]
  )
  terrain <- extract_gaussian_terrain_points(
    cells, base, tile_m = 20000, buffer_m = 900,
    cache_dir = file.path(v6_dir, "tmp", "map_terrain_60m", slug),
    id_field = "cell_id"
  ) |> dplyr::select(-dplyr::any_of("terrain_row_id"))
  cells <- dplyr::left_join(cells, terrain, by = "cell_id")
  X <- build_density_X(cells)

  rasters <- list()
  for (model_name in names(models)) {
    values <- rep(NA_real_, nrow(cells))
    valid <- stats::complete.cases(X)
    values[valid] <- exp(as.numeric(X[valid, , drop = FALSE] %*% models[[model_name]]$beta))
    r <- template; terra::values(r) <- NA_real_
    target <- terra::cellFromXY(r, as.matrix(cells[, c("projected_x", "projected_y")]))
    r[target] <- values
    names(r) <- model_name
    rasters[[model_name]] <- r
    terra::writeRaster(
      r, file.path(map_dir, paste0(slug, "_", model_name, "_60m.tif")),
      overwrite = TRUE
    )
  }

  legacy_aligned <- align_raster_to_template(legacy, rasters$salad_native)
  v2_file <- file.path(
    cfg$paths$v2_outputs_dir, "hmsc_v2_reef_maps",
    sprintf("%02d_%s_hmsc_v2_static_habitat.tif", i, slug)
  )
  v2_aligned <- align_raster_to_template(terra::rast(v2_file), rasters$salad_native)
  baseline_aligned <- align_raster_to_template(baseline, rasters$salad_native)

  grDevices::png(
    file.path(figure_dir, paste0("density_model_comparison_", slug, ".png")),
    width = 3200, height = 1800, res = 180
  )
  graphics::par(mfrow = c(2, 3), mar = c(1.5, 1.5, 3, 3))
  pal <- viridisLite::viridis(100)
  terra::plot(rank_raster_v4(legacy_aligned), col = pal, range = c(0, 1),
              main = "Legacy clean + year\nwithin-reef percentile", axes = FALSE)
  terra::plot(rank_raster_v4(v2_aligned), col = pal, range = c(0, 1),
              main = "v2 static joint habitat\nwithin-reef percentile", axes = FALSE)
  terra::plot(rank_raster_v4(baseline_aligned), col = pal, range = c(0, 1),
              main = "v4 hybrid 120 m\nwithin-reef percentile", axes = FALSE)
  terra::plot(rank_raster_v4(rasters$cull_native), col = pal, range = c(0, 1),
              main = "Cull-native density 60 m\nwithin-reef percentile", axes = FALSE)
  terra::plot(rank_raster_v4(rasters$salad_native), col = pal, range = c(0, 1),
              main = "SALAD-native density 60 m\nwithin-reef percentile", axes = FALSE)
  terra::plot(rank_raster_v4(rasters$joint_balanced), col = pal, range = c(0, 1),
              main = "Balanced Rosetta joint 60 m\nwithin-reef percentile", axes = FALSE)
  graphics::mtext(reef_names[[i]], outer = TRUE, line = -1, cex = 1.25, font = 2)
  grDevices::dev.off()

  summaries[[i]] <- tibble::tibble(
    reef_name = reef_names[[i]], support_m = support_m,
    dem_domain_cells = nrow(cells),
    enhanced_coverage = mean(is.finite(terra::values(rasters$salad_native))),
    salad_density_median = stats::median(terra::values(rasters$salad_native), na.rm = TRUE),
    salad_density_p95 = stats::quantile(terra::values(rasters$salad_native), 0.95, na.rm = TRUE),
    salad_vs_legacy_spearman = cor_raster(rasters$salad_native, legacy_aligned),
    salad_vs_v2_spearman = cor_raster(rasters$salad_native, v2_aligned),
    salad_vs_v4_spearman = cor_raster(rasters$salad_native, baseline_aligned),
    cull_vs_salad_spearman = cor_raster(rasters$cull_native, rasters$salad_native),
    balanced_vs_salad_spearman = cor_raster(rasters$joint_balanced, rasters$salad_native),
    expert_review_status = "pending",
    figure = file.path("figures", paste0("density_model_comparison_", slug, ".png"))
  )
}

manifest <- dplyr::bind_rows(summaries)
utils::write.csv(
  manifest, file.path(ctx$paths$outputs, "four_reef_density_map_manifest.csv"),
  row.names = FALSE
)
print(manifest)
