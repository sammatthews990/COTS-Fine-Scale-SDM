args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v3_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v3_dir, "R", "bootstrap.R"))
ctx <- initialise_v3(script_file)
cfg <- ctx$config
base <- terra::rast(cfg$paths$predictor_stack_rg)
geomorphic <- terra::rast(cfg$paths$geomorphic_raster)
requested_support <- suppressWarnings(as.integer(commandArgs(trailingOnly = TRUE)))
supports <- if (length(requested_support) && all(is.finite(requested_support))) {
  intersect(requested_support, cfg$validation$cell_support_m)
} else cfg$validation$cell_support_m

cell_polygons <- function(cells, support_m, crs) {
  half <- support_m / 2
  geometry <- sf::st_sfc(lapply(seq_len(nrow(cells)), function(i) {
    x <- cells$projected_x[[i]]
    y <- cells$projected_y[[i]]
    sf::st_polygon(list(matrix(c(
      x - half, y - half,
      x + half, y - half,
      x + half, y + half,
      x - half, y + half,
      x - half, y - half
    ), ncol = 2, byrow = TRUE)))
  }), crs = crs)
  sf::st_sf(cell_id = cells$cell_id, geometry = geometry)
}

extract_support <- function(cell_size_m) {
  message("Extracting support-matched predictors at ", cell_size_m, " m.")
  cells <- readRDS(file.path(
    cfg$paths$artifacts_dir, paste0("predictor_cells_", cell_size_m, "m.rds")
  ))
  offsets <- seq(-cell_size_m / 2 + 15, cell_size_m / 2 - 15, by = 30)
  offset_grid <- expand.grid(dx = offsets, dy = offsets)
  base_cache <- file.path(cfg$paths$tmp_dir, paste0("base_support_", cell_size_m, "m.rds"))
  if (file.exists(base_cache)) {
    base_values <- readRDS(base_cache)
  } else {
    extraction_chunks <- split(
      seq_len(nrow(cells)), ceiling(seq_len(nrow(cells)) / 2000L)
    )
    base_values <- dplyr::bind_rows(lapply(seq_along(extraction_chunks), function(i) {
      message("  support quadrature chunk ", i, " of ", length(extraction_chunks))
      rows <- extraction_chunks[[i]]
      x <- cells$projected_x[rows]
      y <- cells$projected_y[rows]
      samples_per_cell <- nrow(offset_grid)
      xy <- cbind(
        rep(x, each = samples_per_cell) + rep(offset_grid$dx, times = length(rows)),
        rep(y, each = samples_per_cell) + rep(offset_grid$dy, times = length(rows))
      )
      values <- terra::extract(base, xy)[, -1, drop = FALSE]
      means <- vapply(values, function(column) {
        matrix_values <- matrix(column, nrow = length(rows), byrow = TRUE)
        value <- rowMeans(matrix_values, na.rm = TRUE)
        value[!is.finite(value)] <- NA_real_
        value
      }, numeric(length(rows)))
      tibble::as_tibble(means)
    }))
    saveRDS(base_values, base_cache)
  }

  centres <- terra::vect(
    as.data.frame(cells), geom = c("projected_x", "projected_y"),
    crs = terra::crs(base)
  )
  centres_lonlat <- terra::project(centres, terra::crs(geomorphic))
  lonlat <- terra::crds(centres_lonlat)
  geom_cache <- file.path(cfg$paths$tmp_dir, paste0("geomorphic_", cell_size_m, "m.rds"))
  if (file.exists(geom_cache)) {
    geom_labels <- readRDS(geom_cache)
  } else {
    tile_id <- paste(floor(lonlat[, 1] / 0.1), floor(lonlat[, 2] / 0.1), sep = "||")
    groups <- split(seq_len(nrow(cells)), tile_id)
    geom_labels <- rep(NA_character_, nrow(cells))
    for (i in seq_along(groups)) {
      if (i %% 25L == 1L || i == length(groups)) {
        message("  geomorphic tile ", i, " of ", length(groups))
      }
      rows <- groups[[i]]
      extent <- terra::ext(
        min(lonlat[rows, 1]) - 0.002, max(lonlat[rows, 1]) + 0.002,
        min(lonlat[rows, 2]) - 0.002, max(lonlat[rows, 2]) + 0.002
      )
      raster_extent <- terra::ext(geomorphic)
      overlaps <- terra::xmin(extent) <= terra::xmax(raster_extent) &&
        terra::xmax(extent) >= terra::xmin(raster_extent) &&
        terra::ymin(extent) <= terra::ymax(raster_extent) &&
        terra::ymax(extent) >= terra::ymin(raster_extent)
      if (!overlaps) next
      tile <- terra::crop(geomorphic, extent, snap = "out")
      extracted <- terra::extract(tile, lonlat[rows, , drop = FALSE])
      if (ncol(extracted) > 0L) {
        geom_labels[rows] <- as.character(extracted[[ncol(extracted)]])
      }
    }
    saveRDS(geom_labels, geom_cache)
  }

  terrain_cells <- cells[is.finite(base_values$DEM), , drop = FALSE]
  terrain <- extract_gaussian_terrain_points(
    terrain_cells, base, tile_m = 5000, buffer_m = 900,
    cache_dir = file.path(cfg$paths$tmp_dir, paste0("terrain_", cell_size_m, "m")),
    id_field = "cell_id"
  ) |>
    dplyr::select(-dplyr::any_of("terrain_row_id"))

  result <- dplyr::bind_cols(
    cells,
    tibble::as_tibble(base_values),
    tibble::tibble(geomorphic_class = geom_labels)
  ) |>
    dplyr::left_join(terrain, by = "cell_id")
  saveRDS(result, file.path(
    cfg$paths$artifacts_dir, paste0("support_predictors_", cell_size_m, "m.rds")
  ))
  tibble::tibble(
    support_m = cell_size_m, cells = nrow(result),
    support_extraction = paste0(nrow(offset_grid), "_point_30m_subcell_mean"),
    base_complete = sum(stats::complete.cases(base_values)),
    geomorphic_observed = sum(!is.na(result$geomorphic_class)),
    gaussian_complete = sum(stats::complete.cases(
      result[, grep("_g(30|60|120|240|360)m$", names(result)), drop = FALSE]
    ))
  )
}

audit <- dplyr::bind_rows(lapply(supports, extract_support))
audit_path <- file.path(cfg$paths$outputs_dir, "support_predictor_audit.csv")
if (file.exists(audit_path)) {
  prior <- suppressWarnings(readr::read_csv(audit_path, show_col_types = FALSE))
  audit <- dplyr::bind_rows(prior, audit) |>
    dplyr::arrange(.data$support_m) |>
    dplyr::distinct(.data$support_m, .keep_all = TRUE)
}
utils::write.csv(audit, audit_path, row.names = FALSE)
print(audit)
