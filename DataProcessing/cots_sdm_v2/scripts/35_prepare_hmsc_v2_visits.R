args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v2_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v2_dir, "R", "bootstrap.R"))
ctx <- initialise_v2(script_file)
cfg <- ctx$config

required <- c(
  cfg$paths$cull_workbook, cfg$paths$manta_current_workbook,
  cfg$paths$salad_workbook, cfg$paths$predictor_stack
)
if (any(!file.exists(required))) stop("Missing Hmsc v2 visit input(s).")

template <- terra::rast(cfg$paths$predictor_stack)
grid_m <- 200

cull_raw <- read_cull_observations(cfg$paths$cull_workbook)
cull <- dplyr::transmute(
  cull_raw,
  visit_id = paste0(
    "cull||", .data$CrownOfThornsStarfishCullDiveId, "||",
    .data$source_record_id
  ),
  method = "cull",
  survey_date = .data$survey_date,
  year = .data$year,
  reef_label = as.character(.data$ReefLabel),
  reef_name = .data$reef_name,
  reef_key = normalise_reef_key(.data$reef_name),
  longitude = .data$longitude,
  latitude = .data$latitude,
  footprint_length_m = NA_real_,
  cohort1 = as.numeric(.data$Cohort1),
  cohort2 = as.numeric(.data$Cohort2),
  cohort3 = as.numeric(.data$Cohort3),
  cohort4 = as.numeric(.data$Cohort4),
  direct_cots_count = .data$total_cots,
  bottom_time = .data$bottom_time,
  cull_dive_type = as.character(.data$CullDiveTypeName)
)

manta_raw <- read_manta_observations(
  cfg$paths$manta_current_workbook, sheet = "Manta"
)
manta <- dplyr::transmute(
  manta_raw,
  visit_id = .data$manta_record_id,
  method = "manta",
  survey_date = .data$survey_date,
  year = .data$year,
  reef_label = NA_character_,
  reef_name = .data$reef_name,
  reef_key = normalise_reef_key(.data$reef_name),
  longitude = .data$longitude,
  latitude = .data$latitude,
  footprint_length_m = .data$tow_distance_m,
  direct_cots_count = .data$cots_count,
  scar_count = .data$scar_count,
  scar_code = .data$scar_code,
  tow_distance_m = .data$tow_distance_m,
  hard_coral_observed = .data$hard_coral_observed,
  strict_negative = .data$strict_negative,
  informative = .data$informative,
  valid_coordinates = .data$valid_coordinates
) |>
  dplyr::filter(.data$valid_coordinates)

rhis <- read_rhis_observations(cfg$paths$cull_workbook)
salad <- read_salad_tracks(cfg$paths$salad_workbook)
salad_events <- interpolate_salad_records(cfg$paths$salad_workbook, salad)

method_tables <- list(cull = cull, manta = manta, rhis = rhis, salad = salad)
method_tables <- lapply(method_tables, assign_visit_cells, template = template, grid_m = grid_m)

common <- c(
  "visit_id", "method", "survey_date", "year", "reef_label", "reef_name",
  "reef_key", "longitude", "latitude", "footprint_length_m", "projected_x",
  "projected_y", "x_bin", "y_bin", "cell_id"
)
visits_index <- dplyr::bind_rows(lapply(method_tables, function(x) {
  dplyr::select(x, dplyr::all_of(common))
})) |>
  annotate_local_cull_history() |>
  cluster_visit_dates(maximum_gap_days = 14L)

add_annotations <- function(data) {
  data |>
    dplyr::left_join(
      dplyr::select(
        visits_index, "visit_id", "local_first_cull_date",
        "days_from_local_first_cull", "cull_timing", "ecological_eligible",
        "visit_cluster_id"
      ),
      by = "visit_id"
    )
}
method_tables <- lapply(method_tables, add_annotations)
cull <- method_tables$cull
manta <- method_tables$manta
rhis <- method_tables$rhis
salad <- method_tables$salad

cull_units <- cull |>
  dplyr::filter(.data$ecological_eligible) |>
  dplyr::group_by(.data$visit_cluster_id, .data$cell_id) |>
  dplyr::summarise(
    cull_total_cots = sum(.data$direct_cots_count, na.rm = TRUE),
    cull_bottom_time = sum(.data$bottom_time, na.rm = TRUE),
    cull_n_dives = dplyr::n(),
    cull_direct_presence = as.integer(.data$cull_total_cots > 0),
    cull_cpue = .data$cull_total_cots / .data$cull_bottom_time,
    cull_high_density = as.integer(.data$cull_cpue >= cfg$response$cpue_threshold),
    .groups = "drop"
  )

manta_units <- manta |>
  dplyr::filter(.data$ecological_eligible) |>
  dplyr::group_by(.data$visit_cluster_id, .data$cell_id) |>
  dplyr::summarise(
    manta_direct_presence = as.integer(any(
      dplyr::coalesce(.data$direct_cots_count > 0, FALSE)
    )),
    manta_scar_presence = as.integer(any(
      dplyr::coalesce(.data$scar_count > 0, FALSE) |
        dplyr::coalesce(.data$scar_code %in% c("p", "c"), FALSE)
    )),
    manta_n_tows = dplyr::n(),
    manta_total_distance_m = sum(.data$tow_distance_m, na.rm = TRUE),
    .groups = "drop"
  )

rhis_units <- rhis |>
  dplyr::filter(.data$ecological_eligible) |>
  dplyr::group_by(.data$visit_cluster_id, .data$cell_id) |>
  dplyr::summarise(
    rhis_adult_presence = if (any(.data$rhis_cots_protocol_observed)) {
      as.integer(any(.data$adult_count > 0, na.rm = TRUE))
    } else NA_integer_,
    rhis_juvenile_presence = if (any(.data$rhis_cots_protocol_observed)) {
      as.integer(any(.data$juvenile_count > 0, na.rm = TRUE))
    } else NA_integer_,
    rhis_branching_scar_presence = if (any(.data$rhis_cots_protocol_observed)) {
      as.integer(any(.data$branching_scar_count > 0, na.rm = TRUE))
    } else NA_integer_,
    rhis_table_scar_presence = if (any(.data$rhis_cots_protocol_observed)) {
      as.integer(any(.data$table_scar_count > 0, na.rm = TRUE))
    } else NA_integer_,
    rhis_branching_live_pct = mean(.data$branching_live_pct, na.rm = TRUE),
    rhis_table_live_pct = mean(.data$table_live_pct, na.rm = TRUE),
    rhis_branching_recent_dead_pct = mean(
      .data$branching_recent_dead_pct, na.rm = TRUE
    ),
    rhis_table_recent_dead_pct = mean(.data$table_recent_dead_pct, na.rm = TRUE),
    rhis_n_surveys = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::mutate(dplyr::across(
    dplyr::starts_with("rhis_"),
    ~ dplyr::if_else(is.nan(.x), NA_real_, as.numeric(.x))
  ))

salad_units <- salad |>
  dplyr::filter(.data$ecological_eligible) |>
  dplyr::group_by(.data$visit_cluster_id, .data$cell_id) |>
  dplyr::summarise(
    salad_density_incl_scars_ha = stats::weighted.mean(
      .data$density_incl_scars_ha, w = .data$footprint_area_m2,
      na.rm = TRUE
    ),
    salad_direct_presence = as.integer(any(.data$direct_cots_count > 0)),
    salad_scar_presence = as.integer(any(.data$scar_count > 0)),
    salad_n_tracks = dplyr::n(),
    salad_total_area_m2 = sum(.data$footprint_area_m2, na.rm = TRUE),
    .groups = "drop"
  )

unit_index <- visits_index |>
  dplyr::filter(.data$ecological_eligible) |>
  dplyr::group_by(.data$visit_cluster_id, .data$cell_id) |>
  dplyr::summarise(
    survey_date_start = min(.data$survey_date),
    survey_date_end = max(.data$survey_date),
    year = lubridate::year(stats::median(.data$survey_date)),
    reef_label = dplyr::first(stats::na.omit(.data$reef_label), default = NA_character_),
    reef_name = dplyr::first(stats::na.omit(.data$reef_name), default = NA_character_),
    reef_key = dplyr::first(stats::na.omit(.data$reef_key), default = NA_character_),
    x_bin = dplyr::first(.data$x_bin),
    y_bin = dplyr::first(.data$y_bin),
    projected_x = (.data$x_bin + 0.5) * grid_m,
    projected_y = (.data$y_bin + 0.5) * grid_m,
    methods = paste(sort(unique(.data$method)), collapse = ";"),
    method_count = dplyr::n_distinct(.data$method),
    .groups = "drop"
  )

model_units <- unit_index |>
  dplyr::left_join(cull_units, by = c("visit_cluster_id", "cell_id")) |>
  dplyr::left_join(manta_units, by = c("visit_cluster_id", "cell_id")) |>
  dplyr::left_join(rhis_units, by = c("visit_cluster_id", "cell_id")) |>
  dplyr::left_join(salad_units, by = c("visit_cluster_id", "cell_id"))

cell_points <- terra::vect(
  as.data.frame(model_units), geom = c("projected_x", "projected_y"),
  crs = terra::crs(template)
)
cell_points <- terra::project(cell_points, "EPSG:4326")
lonlat <- terra::crds(cell_points)
model_units$longitude <- lonlat[, 1]
model_units$latitude <- lonlat[, 2]

saveRDS(visits_index, file.path(
  cfg$paths$artifacts_dir, "hmsc_v2_visit_index.rds"
))
saveRDS(model_units, file.path(
  cfg$paths$artifacts_dir, "hmsc_v2_model_units.rds"
))
saveRDS(method_tables, file.path(
  cfg$paths$artifacts_dir, "hmsc_v2_annotated_method_tables.rds"
))
saveRDS(salad_events, file.path(
  cfg$paths$artifacts_dir, "salad_cots_recorded_interpolation_diagnostic.rds"
))

alignment_audit <- model_units |>
  dplyr::count(.data$methods, .data$method_count, name = "model_units") |>
  dplyr::arrange(dplyr::desc(.data$model_units))
pre_cull_audit <- visits_index |>
  dplyr::count(.data$method, .data$cull_timing, .data$ecological_eligible)
salad_audit <- tibble::tibble(
  tracks_total = nrow(suppressWarnings(readxl::read_excel(
    cfg$paths$salad_workbook, sheet = "Tracks"
  ))),
  tracks_with_valid_midpoint = nrow(salad),
  median_track_length_m = stats::median(salad$footprint_length_m, na.rm = TRUE),
  recorded_events = nrow(salad_events),
  events_matched_to_track = sum(!is.na(salad_events$visit_id)),
  events_with_valid_time_interpolation = sum(salad_events$interpolation_valid)
)
utils::write.csv(alignment_audit, file.path(
  cfg$paths$outputs_dir, "hmsc_v2_alignment_audit.csv"
), row.names = FALSE)
utils::write.csv(pre_cull_audit, file.path(
  cfg$paths$outputs_dir, "hmsc_v2_pre_cull_audit.csv"
), row.names = FALSE)
utils::write.csv(salad_audit, file.path(
  cfg$paths$outputs_dir, "hmsc_v2_salad_audit.csv"
), row.names = FALSE)

message("Hmsc v2 visit alignment complete.")
print(alignment_audit)
print(pre_cull_audit)
print(salad_audit)
