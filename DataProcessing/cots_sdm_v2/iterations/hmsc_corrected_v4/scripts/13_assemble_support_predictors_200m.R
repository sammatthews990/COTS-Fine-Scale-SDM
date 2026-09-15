args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v4_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
v2_dir <- normalizePath(file.path(v4_dir, "..", ".."), winslash = "/")
source(file.path(v4_dir, "R", "bootstrap.R"))
source(file.path(v2_dir, "R", "gaussian_terrain.R"))
ctx <- initialise_v4(script_file)
cfg <- ctx$config

base <- readRDS(file.path(cfg$paths$artifacts_dir, "base_support_predictors_200m.rds"))
template <- terra::rast(cfg$paths$predictor_stack_rg)
geomorphic <- terra::rast(cfg$paths$geomorphic_raster)
terrain_cache <- readRDS(file.path(cfg$paths$artifacts_dir, "gaussian_cache_v2_union.rds"))
terrain_names <- grep(
  "^(dem|slope|eastness|northness|hcu|vcu|bpi|sd)_g(30|60|120|240|360)m$",
  names(terrain_cache), value = TRUE
)
terrain_cache <- terrain_cache |>
  dplyr::select(.data$cell_id, dplyr::all_of(terrain_names)) |>
  dplyr::distinct(.data$cell_id, .keep_all = TRUE) |>
  dplyr::rename(cache_cell_id = .data$cell_id)

base <- base |>
  dplyr::mutate(cache_cell_id = sub("^200m\\|\\|", "", .data$cell_id)) |>
  dplyr::left_join(terrain_cache, by = "cache_cell_id")

message("Extracting geomorphic class at valid 200 m cell centres.")
centres <- terra::vect(
  base, geom = c("projected_x", "projected_y"),
  crs = terra::crs(template)
)
centres_geomorphic <- terra::project(centres, terra::crs(geomorphic))
geomorphic_value <- terra::extract(geomorphic, centres_geomorphic)
base$geomorphic_class <- as.character(geomorphic_value[[2]])

missing <- base |>
  dplyr::filter(!is.finite(.data$dem_g30m), is.finite(.data$DEM))
dem_extent <- terra::ext(template[["DEM"]])
missing <- missing |>
  dplyr::filter(
    .data$projected_x >= dem_extent$xmin,
    .data$projected_x <= dem_extent$xmax,
    .data$projected_y >= dem_extent$ymin,
    .data$projected_y <= dem_extent$ymax
  ) |>
  dplyr::select(.data$cell_id, .data$projected_x, .data$projected_y)

if (nrow(missing)) {
  message("Computing Gaussian terrain for ", nrow(missing), " uncached valid cells.")
  computed <- extract_gaussian_terrain_points(
    missing, template, tile_m = 20000, buffer_m = 900,
    cache_dir = file.path(cfg$paths$tmp_dir, "gaussian_terrain_200m_missing"),
    id_field = "cell_id"
  ) |>
    dplyr::select(-dplyr::any_of("terrain_row_id"))
  computed <- computed |>
    dplyr::rename_with(~ paste0(.x, "_computed"), dplyr::all_of(terrain_names))
  base <- base |>
    dplyr::left_join(computed, by = "cell_id")
  for (name in terrain_names) {
    computed_name <- paste0(name, "_computed")
    base[[name]] <- dplyr::coalesce(base[[name]], base[[computed_name]])
    base[[computed_name]] <- NULL
  }
}

result <- base |>
  dplyr::select(-.data$cache_cell_id)
saveRDS(result, file.path(cfg$paths$artifacts_dir, "support_predictors_200m.rds"))

base_names <- names(template)
audit <- tibble::tibble(
  support_m = 200L,
  cells = nrow(result),
  base_complete = sum(stats::complete.cases(result[, base_names, drop = FALSE])),
  native_bpi_present = "BPI" %in% names(result),
  geomorphic_observed = sum(!is.na(result$geomorphic_class)),
  gaussian_30m_observed = sum(is.finite(result$dem_g30m)),
  gaussian_all_complete = sum(stats::complete.cases(result[, terrain_names, drop = FALSE]))
)
utils::write.csv(
  audit,
  file.path(cfg$paths$outputs_dir, "support_predictor_audit_200m.csv"),
  row.names = FALSE
)
print(audit)
