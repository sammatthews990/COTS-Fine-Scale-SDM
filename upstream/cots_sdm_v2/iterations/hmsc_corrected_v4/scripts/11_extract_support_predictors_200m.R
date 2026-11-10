args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v4_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
v2_dir <- normalizePath(file.path(v4_dir, "..", ".."), winslash = "/")
source(file.path(v4_dir, "R", "bootstrap.R"))
source(file.path(v2_dir, "R", "gaussian_terrain.R"))
ctx <- initialise_v4(script_file)
cfg <- ctx$config
support_m <- 200L
base <- terra::rast(cfg$paths$predictor_stack_rg)
geomorphic <- terra::rast(cfg$paths$geomorphic_raster)
visits <- readRDS(file.path(cfg$paths$artifacts_dir, "cell_visits_200m_window_3d.rds"))
cells <- visits |>
  dplyr::filter(.data$method != "manta") |>
  dplyr::distinct(.data$cell_id, .data$projected_x, .data$projected_y)

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

message("Extracting exact area-weighted base predictors for 200 m cells.")
chunks <- split(seq_len(nrow(cells)), ceiling(seq_len(nrow(cells)) / 500L))
base_rows <- lapply(seq_along(chunks), function(i) {
  if (i == 1L || i == length(chunks)) message("Base polygon chunk ", i, " of ", length(chunks))
  rows <- chunks[[i]]
  polygons <- terra::vect(cell_polygons(cells[rows, , drop = FALSE], terra::crs(base)))
  value <- terra::extract(base, polygons, fun = mean, na.rm = TRUE)
  value <- value[, setdiff(names(value), "ID"), drop = FALSE]
  dplyr::bind_cols(cells[rows, , drop = FALSE], tibble::as_tibble(value))
})
base_values <- dplyr::bind_rows(base_rows)

message("Extracting modal geomorphic class for 200 m cells.")
centres <- terra::vect(cells, geom = c("projected_x", "projected_y"), crs = terra::crs(base))
centres_lonlat <- terra::project(centres, terra::crs(geomorphic))
geomorphic_class <- as.character(terra::extract(geomorphic, centres_lonlat)[[2]])
base_values$geomorphic_class <- geomorphic_class

message("Extracting Gaussian terrain at 200 m cell centres.")
terrain <- extract_gaussian_terrain_points(
  cells, base, tile_m = 5000, buffer_m = 900,
  cache_dir = file.path(cfg$paths$tmp_dir, "gaussian_terrain_200m"),
  id_field = "cell_id"
) |>
  dplyr::select(-dplyr::any_of("terrain_row_id"))
result <- dplyr::left_join(base_values, terrain, by = "cell_id")
saveRDS(result, file.path(cfg$paths$artifacts_dir, "support_predictors_200m.rds"))
raw_terms <- unlist(terrain_families(names(result)), use.names = FALSE)
audit <- tibble::tibble(
  support_m = support_m, cells = nrow(result),
  base_complete = sum(stats::complete.cases(result[, names(base), drop = FALSE])),
  native_bpi_present = "BPI" %in% names(result),
  geomorphic_observed = sum(!is.na(result$geomorphic_class)),
  gaussian_complete = sum(stats::complete.cases(result[, raw_terms, drop = FALSE]))
)
utils::write.csv(audit, file.path(cfg$paths$outputs_dir, "support_predictor_audit_200m.csv"), row.names = FALSE)
print(audit)
