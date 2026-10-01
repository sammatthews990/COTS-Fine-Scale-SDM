# Build the full-domain 120 m predictor table for the Hmsc layer.
#
#   Rscript DataProcessing/cots_sdm_v2/final/scripts/04_build_support_grid_120m.R [--workers=18]
#
# Stages are cached under final/cache/ and resume if interrupted. Before the
# full run, every predictor definition is checked against the frozen v4
# training cells; the script stops if they do not reproduce.

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "..", "R", "common.R"))
cfg <- initialise_final()
args <- script_args()
workers <- as.integer(if (is.null(args$workers)) 18L else args$workers)

stack_path <- checked_input(cfg, "predictor_stack_rg")
geomorphic_path <- checked_input(cfg, "geomorphic")
stack <- terra::rast(stack_path)
cache <- cfg$dirs$cache

# --- Stage A: 4x4 block-mean base layers on the 0-anchored 120 m grid -------
base_path <- file.path(cache, "base_support_120m.tif")
if (!file.exists(base_path)) {
  alignment <- quadrature_alignment(stack)
  log_message("Aggregating ", length(hmsc_base_layers), " base layers to 120 m (shift ",
              round(alignment$dx, 2), ", ", round(alignment$dy, 2), " m)")
  cropped <- terra::crop(stack[[hmsc_base_layers]], alignment$crop_ext, snap = "near")
  stopifnot(isTRUE(all.equal(as.vector(terra::ext(cropped)), as.vector(alignment$crop_ext))))
  aggregated <- terra::aggregate(cropped, fact = 4, fun = "mean", na.rm = TRUE)
  aggregated <- terra::shift(aggregated, dx = alignment$dx, dy = alignment$dy)
  terra::writeRaster(aggregated, base_path, overwrite = TRUE, datatype = "FLT4S",
                     gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3", "TILED=YES", "BIGTIFF=IF_SAFER"))
}
base <- terra::rast(base_path)
stopifnot(all(terra::ext(base)[c(1, 3)] %% 120 == 0))

# --- Verification against frozen training cells -----------------------------
reference <- readRDS(checked_input(cfg, "hmsc_reference_cells")) |>
  dplyr::distinct(.data$projected_x, .data$projected_y, .keep_all = TRUE)
check_path <- file.path(cfg$dirs$tables, "support_grid_reference_check.csv")
if (!file.exists(check_path)) {
  log_message("Verifying predictor definitions on ", nrow(reference), " training cells")
  xy <- cbind(reference$projected_x, reference$projected_y)
  base_ref <- terra::extract(base, xy)
  base_check <- dplyr::bind_rows(lapply(hmsc_base_layers, function(v) {
    a <- base_ref[[v]]
    b <- reference[[v]]
    both <- is.finite(a) & is.finite(b)
    tibble::tibble(
      predictor = v, family = "base", compared = sum(both),
      missing_in_grid = sum(!is.finite(a) & is.finite(b)),
      extra_in_grid = sum(is.finite(a) & !is.finite(b)),
      max_relative_error = max(abs(a[both] - b[both]) / pmax(abs(b[both]), 1e-3))
    )
  }))

  geom <- extract_geomorphic_class(xy[, 1], xy[, 2], terra::crs(stack), geomorphic_path)
  both <- !is.na(geom) & !is.na(reference$geomorphic_class)
  geom_check <- tibble::tibble(
    predictor = "geomorphic_class", family = "geomorphic", compared = sum(both),
    missing_in_grid = sum(is.na(geom) & !is.na(reference$geomorphic_class)),
    extra_in_grid = sum(!is.na(geom) & is.na(reference$geomorphic_class)),
    max_relative_error = mean(geom[both] != reference$geomorphic_class[both])
  )

  # As in v3, terrain is only computed where the block-mean DEM is finite.
  points <- data.frame(cell = seq_len(nrow(reference)), x = xy[, 1], y = xy[, 2])
  points <- points[is.finite(base_ref$DEM), ]
  tiles <- split(points, paste(floor(points$x / 5000), floor(points$y / 5000)))
  started <- Sys.time()
  cl <- parallel::makeCluster(workers)
  parallel::clusterExport(cl, c("point_gaussian_terrain", "fast_gaussian_terrain", "tile_gaussian_terrain",
                                "gaussian_scales_m", "terrain_column_names"))
  terrain_ref <- dplyr::bind_rows(parallel::parLapplyLB(cl, tiles, tile_gaussian_terrain, dem_path = stack_path))
  parallel::stopCluster(cl)
  seconds_per_point <- as.numeric(difftime(Sys.time(), started, units = "secs")) * workers / nrow(points)
  terrain_ref <- dplyr::left_join(data.frame(cell = seq_len(nrow(reference))), terrain_ref, by = "cell")
  terrain_check <- dplyr::bind_rows(lapply(terrain_column_names, function(v) {
    a <- terrain_ref[[v]]
    b <- reference[[v]]
    both <- is.finite(a) & is.finite(b)
    tibble::tibble(
      predictor = v, family = "terrain", compared = sum(both),
      missing_in_grid = sum(!is.finite(a) & is.finite(b)),
      extra_in_grid = sum(is.finite(a) & !is.finite(b)),
      max_relative_error = max(abs(a[both] - b[both]) / pmax(abs(b[both]), 1e-3))
    )
  }))
  check <- dplyr::bind_rows(base_check, geom_check, terrain_check)
  check$cpu_seconds_per_terrain_point <- seconds_per_point
  write_table(check, cfg, basename(check_path))
  # The grid may have values where v3 recorded NA (those training rows were
  # dropped as incomplete); it must never miss or disagree with a value that
  # training used. Base layers are float32, so allow single-precision rounding.
  failed <- check$missing_in_grid > 0 | check$max_relative_error > 1e-4
  if (any(failed)) {
    print(check[failed, ])
    file.rename(check_path, sub("\\.csv$", "_FAILED.csv", check_path))
    stop("Support grid does not reproduce the training predictors.")
  }
}

# --- Stage B: domain cells with base values and geomorphic class -------------
cells_path <- file.path(cache, "support_cells_120m.rds")
if (!file.exists(cells_path)) {
  # Every component needs DEM and the three Acropora layers; cells lacking
  # them cannot be predicted and are left uncovered.
  core <- c("DEM", "AhyaD", "Aspat", "Aten")
  ok <- terra::noNA(base[[core]])
  cell <- which(terra::values(ok, mat = FALSE) == 1L)
  log_message("Domain: ", format(length(cell), big.mark = ","), " cells at 120 m")
  xy <- terra::xyFromCell(base, cell)
  values <- terra::extract(base, cell)
  cells <- data.frame(cell = cell, x = xy[, 1], y = xy[, 2], values)
  log_message("Extracting geomorphic class")
  cells$geomorphic_class <- extract_geomorphic_class(cells$x, cells$y, terra::crs(base), geomorphic_path)
  saveRDS(cells, cells_path)
}
cells <- readRDS(cells_path)

# --- Stage C: multiscale Gaussian terrain at every domain cell ---------------
grid_path <- file.path(cache, "support_grid_120m.rds")
if (!file.exists(grid_path)) {
  tile_dir <- file.path(cache, "terrain_120m_tiles")
  dir.create(tile_dir, showWarnings = FALSE)
  tile_key <- paste(floor(cells$x / 5000), floor(cells$y / 5000), sep = "_")
  tiles <- split(cells[, c("cell", "x", "y")], tile_key)
  todo <- names(tiles)[!file.exists(file.path(tile_dir, paste0(names(tiles), ".rds")))]
  log_message("Terrain: ", length(tiles), " tiles, ", length(todo), " to compute on ", workers, " workers")
  if (length(todo)) {
    cl <- parallel::makeCluster(workers)
    parallel::clusterExport(cl, c("point_gaussian_terrain", "fast_gaussian_terrain", "tile_gaussian_terrain",
                                  "gaussian_scales_m", "terrain_column_names"))
    jobs <- lapply(todo, function(key) list(key = key, points = tiles[[key]]))
    parallel::parLapplyLB(cl, jobs, function(job, dem_path, tile_dir) {
      result <- tile_gaussian_terrain(job$points, dem_path)
      saveRDS(result, file.path(tile_dir, paste0(job$key, ".rds")))
      NULL
    }, dem_path = stack_path, tile_dir = tile_dir)
    parallel::stopCluster(cl)
  }
  terrain <- dplyr::bind_rows(lapply(file.path(tile_dir, paste0(names(tiles), ".rds")), readRDS))
  grid <- dplyr::left_join(cells, terrain, by = "cell")
  stopifnot(nrow(grid) == nrow(cells))
  saveRDS(grid, grid_path)
}
grid <- readRDS(grid_path)

summary <- tibble::tibble(
  domain_cells = nrow(grid),
  enhanced_complete = sum(stats::complete.cases(grid[, c(hmsc_base_layers, terrain_column_names)]) &
                            !is.na(grid$geomorphic_class)),
  clean_multiscale_complete = sum(stats::complete.cases(grid[, c("AhyaD", "Aspat", "Aten", terrain_column_names)])),
  clean_native_complete = sum(stats::complete.cases(grid[, setdiff(hmsc_base_layers, c("turbid", "waves_Hs", "waves_Tp"))]))
)
write_table(summary, cfg, "support_grid_120m_summary.csv")
print(summary)
