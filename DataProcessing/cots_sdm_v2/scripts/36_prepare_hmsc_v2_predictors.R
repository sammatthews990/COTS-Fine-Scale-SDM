args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v2_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v2_dir, "R", "bootstrap.R"))
ctx <- initialise_v2(script_file)
cfg <- ctx$config

unit_path <- file.path(cfg$paths$artifacts_dir, "hmsc_v2_model_units.rds")
required <- c(
  unit_path, cfg$paths$predictor_stack, cfg$paths$predictor_stack_rg,
  cfg$paths$geomorphic_raster
)
if (any(!file.exists(required))) stop("Missing Hmsc v2 predictor input(s).")

units <- readRDS(unit_path)
binary_or_na <- function(x) {
  if (all(is.na(x))) NA_integer_ else as.integer(max(x, na.rm = TRUE))
}
sum_or_na <- function(x) {
  if (all(is.na(x))) NA_real_ else sum(x, na.rm = TRUE)
}
mean_or_na <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}
weighted_or_na <- function(x, weight) {
  valid <- is.finite(x) & is.finite(weight) & weight > 0
  if (!any(valid)) NA_real_ else stats::weighted.mean(x[valid], weight[valid])
}

cell_year <- units |>
  dplyr::group_by(.data$cell_id, .data$year) |>
  dplyr::summarise(
    reef_label = dplyr::first(stats::na.omit(.data$reef_label), default = NA_character_),
    reef_name = dplyr::first(stats::na.omit(.data$reef_name), default = NA_character_),
    reef_key = dplyr::first(stats::na.omit(.data$reef_key), default = NA_character_),
    longitude = dplyr::first(.data$longitude),
    latitude = dplyr::first(.data$latitude),
    projected_x = dplyr::first(.data$projected_x),
    projected_y = dplyr::first(.data$projected_y),
    cull_direct_presence = binary_or_na(.data$cull_direct_presence),
    cull_high_density = binary_or_na(.data$cull_high_density),
    manta_direct_presence_raw = binary_or_na(.data$manta_direct_presence),
    manta_scar_presence_raw = binary_or_na(.data$manta_scar_presence),
    manta_n_tows = sum(.data$manta_n_tows, na.rm = TRUE),
    manta_total_distance_m = sum(.data$manta_total_distance_m, na.rm = TRUE),
    rhis_adult_presence = binary_or_na(.data$rhis_adult_presence),
    rhis_juvenile_presence = binary_or_na(.data$rhis_juvenile_presence),
    rhis_branching_scar_presence = binary_or_na(
      .data$rhis_branching_scar_presence
    ),
    rhis_table_scar_presence = binary_or_na(.data$rhis_table_scar_presence),
    rhis_branching_live_pct = mean_or_na(.data$rhis_branching_live_pct),
    rhis_table_live_pct = mean_or_na(.data$rhis_table_live_pct),
    rhis_branching_recent_dead_pct = mean_or_na(
      .data$rhis_branching_recent_dead_pct
    ),
    rhis_table_recent_dead_pct = mean_or_na(.data$rhis_table_recent_dead_pct),
    salad_density_incl_scars_ha = weighted_or_na(
      .data$salad_density_incl_scars_ha, .data$salad_total_area_m2
    ),
    salad_direct_presence = binary_or_na(.data$salad_direct_presence),
    salad_scar_presence = binary_or_na(.data$salad_scar_presence),
    salad_total_area_m2 = sum(.data$salad_total_area_m2, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    manta_positive = dplyr::coalesce(.data$manta_direct_presence_raw == 1L, FALSE) |
      dplyr::coalesce(.data$manta_scar_presence_raw == 1L, FALSE),
    manta_repeated_negative = !.data$manta_positive &
      .data$manta_n_tows >= 3L & .data$manta_total_distance_m >= 1000,
    manta_ecological_observed = .data$manta_positive |
      .data$manta_repeated_negative,
    manta_direct_presence = dplyr::if_else(
      .data$manta_ecological_observed,
      .data$manta_direct_presence_raw, NA_integer_
    ),
    manta_scar_presence = dplyr::if_else(
      .data$manta_ecological_observed,
      .data$manta_scar_presence_raw, NA_integer_
    ),
    any_positive_evidence = rowSums(dplyr::pick(
      "cull_direct_presence", "manta_direct_presence", "manta_scar_presence",
      "rhis_adult_presence", "rhis_juvenile_presence",
      "rhis_branching_scar_presence", "rhis_table_scar_presence",
      "salad_direct_presence", "salad_scar_presence"
    ) == 1L, na.rm = TRUE) > 0L |
      dplyr::coalesce(.data$salad_density_incl_scars_ha > 0, FALSE),
    reef_year_id = paste(.data$reef_key, .data$year, sep = "||")
  ) |>
  dplyr::group_by(.data$reef_year_id) |>
  dplyr::mutate(outbreak_active = any(.data$any_positive_evidence)) |>
  dplyr::ungroup() |>
  dplyr::filter(
    .data$outbreak_active,
    nzchar(.data$reef_key),
    rowSums(!is.na(dplyr::pick(
      "cull_high_density", "manta_direct_presence", "manta_scar_presence",
      "rhis_adult_presence", "rhis_juvenile_presence",
      "rhis_branching_scar_presence", "rhis_table_scar_presence",
      "salad_density_incl_scars_ha"
    ))) > 0L
  ) |>
  dplyr::mutate(
    rhis_branching_live_asin = asin(sqrt(pmin(pmax(
      .data$rhis_branching_live_pct, 0
    ), 100) / 100)),
    rhis_table_live_asin = asin(sqrt(pmin(pmax(
      .data$rhis_table_live_pct, 0
    ), 100) / 100)),
    rhis_branching_recent_dead_asin = asin(sqrt(pmin(pmax(
      .data$rhis_branching_recent_dead_pct, 0
    ), 100) / 100)),
    rhis_table_recent_dead_asin = asin(sqrt(pmin(pmax(
      .data$rhis_table_recent_dead_pct, 0
    ), 100) / 100)),
    salad_log_density_incl_scars = log1p(.data$salad_density_incl_scars_ha)
  )

response_audit <- tibble::tibble(
  rows = nrow(cell_year),
  cells = dplyr::n_distinct(cell_year$cell_id),
  reefs = dplyr::n_distinct(cell_year$reef_key),
  reef_years = dplyr::n_distinct(cell_year$reef_year_id),
  years = dplyr::n_distinct(cell_year$year),
  manta_positive_rows = sum(cell_year$manta_positive),
  manta_repeated_negative_rows = sum(cell_year$manta_repeated_negative),
  cull_target_rows = sum(!is.na(cell_year$cull_high_density)),
  rhis_cots_rows = sum(!is.na(cell_year$rhis_adult_presence)),
  salad_density_rows = sum(!is.na(cell_year$salad_density_incl_scars_ha))
)
message("Conditional cell-year response table:")
print(response_audit)
saveRDS(cell_year, file.path(
  cfg$paths$artifacts_dir, "hmsc_v2_cell_year_responses_unfiltered.rds"
))

cells <- cell_year |>
  dplyr::distinct(
    .data$cell_id, .data$longitude, .data$latitude,
    .data$projected_x, .data$projected_y
  )
site_predictors <- extract_site_habitat_predictors(
  dplyr::select(cells, "cell_id", "longitude", "latitude"),
  cfg$paths$predictor_stack,
  cfg$paths$predictor_stack_rg,
  cfg$paths$geomorphic_raster
)
site_predictors <- site_predictors |>
  dplyr::select(
    "cell_id", "AhyaD", "Aspat", "Aten",
    "RG_bathy", "RG_slope", "RG_turbid", "RG_waves_Hs", "RG_waves_Tp",
    "geomorphic_class", dplyr::starts_with("geom_"),
    "clean_complete", "rg_complete"
  ) |>
  dplyr::mutate(
    coral_suitability_mean = rowMeans(
      dplyr::pick("AhyaD", "Aspat", "Aten"), na.rm = TRUE
    )
  )
site_predictors <- site_predictors |>
  dplyr::filter(.data$clean_complete == 1L, .data$rg_complete == 1L)
cells <- cells |>
  dplyr::semi_join(site_predictors, by = "cell_id")
cell_year <- cell_year |>
  dplyr::semi_join(cells, by = "cell_id")
response_audit$predictor_complete_rows <- nrow(cell_year)
response_audit$predictor_complete_cells <- dplyr::n_distinct(cell_year$cell_id)
message("Predictor-complete conditional domain:")
print(response_audit)

terrain_cache <- file.path(cfg$paths$artifacts_dir, "hmsc_v2_gaussian_terrain.rds")
if (file.exists(terrain_cache)) {
  terrain <- readRDS(terrain_cache)
  if (!setequal(terrain$cell_id, cells$cell_id)) {
    stop("Gaussian terrain cache does not match the current response domain.")
  }
} else {
  template <- terra::rast(cfg$paths$predictor_stack)
  terrain <- extract_gaussian_terrain_points(
    cells, template, tile_m = 5000, buffer_m = 900,
    cache_dir = file.path(v2_dir, "tmp", "hmsc_v2_gaussian_point_tiles"),
    id_field = "cell_id"
  ) |>
    dplyr::select(-"terrain_row_id")
  saveRDS(terrain, terrain_cache)
}

model_data <- cell_year |>
  dplyr::inner_join(site_predictors, by = "cell_id") |>
  dplyr::inner_join(terrain, by = "cell_id")

saveRDS(cell_year, file.path(
  cfg$paths$artifacts_dir, "hmsc_v2_cell_year_responses.rds"
))
saveRDS(model_data, file.path(
  cfg$paths$artifacts_dir, "hmsc_v2_model_data.rds"
))
utils::write.csv(response_audit, file.path(
  cfg$paths$outputs_dir, "hmsc_v2_response_domain_audit.csv"
), row.names = FALSE)
utils::write.csv(
  tibble::tibble(
    requested_predictor = c(
      "Gaussian terrain 30-360 m", "Geomorphic class", "Wave exposure",
      "Larval/current connectivity"
    ),
    available = c(TRUE, TRUE, TRUE, FALSE),
    implementation = c(
      "DEM-first Gaussian smoothing with re-derived gradients, curvature, BPI and local SD",
      "GBR10 categorical one-hot encoding",
      "Reef Guidance significant wave height and peak period",
      "No compatible connectivity layer found locally; omitted rather than proxied by reef identity"
    )
  ),
  file.path(cfg$paths$outputs_dir, "hmsc_v2_predictor_availability.csv"),
  row.names = FALSE
)
message("Hmsc v2 response and Gaussian predictor preparation complete.")
