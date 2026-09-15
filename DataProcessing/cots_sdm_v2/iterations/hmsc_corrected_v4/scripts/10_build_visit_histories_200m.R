args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v4_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
v2_dir <- normalizePath(file.path(v4_dir, "..", ".."), winslash = "/")
source(file.path(v4_dir, "R", "bootstrap.R"))
source(file.path(v2_dir, "R", "joint_visits.R"))
ctx <- initialise_v4(script_file)
cfg <- ctx$config
support_m <- 200L
target_crs <- terra::crs(terra::rast(cfg$paths$predictor_stack))
links <- readRDS(file.path(cfg$paths$artifacts_dir, "footprint_links_200m.rds"))
salad_events <- readRDS(file.path(
  cfg$paths$v2_artifacts_dir, "salad_cots_recorded_interpolation_diagnostic.rds"
)) |>
  dplyr::filter(.data$interpolation_valid)

manta <- links$manta |>
  dplyr::transmute(
    visit_id, method = "manta", survey_date, year, reef_name, reef_key,
    cell_id, cell_support_m, projected_x, projected_y,
    effort_m = .data$weighted_length_m, effort_minutes = NA_real_,
    direct_detection = dplyr::if_else(.data$strict_negative, 0L, NA_integer_),
    scar_detection = dplyr::if_else(.data$strict_negative, 0L, NA_integer_),
    high_density = NA_integer_, direct_count = .data$direct_cots_count,
    scar_count, direct_localized = .data$strict_negative,
    strict_negative, geometry_fallback, geometry_quality_reason,
    source_fraction = .data$fraction
  )
salad_event_cells <- point_cell_links(
  dplyr::transmute(
    salad_events, visit_id, longitude, latitude, survey_date, evidence
  ), target_crs, support_m
) |>
  dplyr::group_by(.data$visit_id, .data$cell_id) |>
  dplyr::summarise(
    event_direct = as.integer(any(.data$evidence == "cots")),
    event_scar = as.integer(any(.data$evidence == "scar")), .groups = "drop"
  )
salad <- links$salad |>
  dplyr::left_join(salad_event_cells, by = c("visit_id", "cell_id")) |>
  dplyr::mutate(
    track_negative = dplyr::coalesce(.data$direct_cots_count, 0) == 0 &
      dplyr::coalesce(.data$scar_count, 0) == 0
  ) |>
  dplyr::transmute(
    visit_id, method = "salad", survey_date, year, reef_name, reef_key,
    cell_id, cell_support_m, projected_x, projected_y,
    effort_m = .data$weighted_length_m,
    effort_minutes = .data$fraction * .data$survey_minutes,
    direct_detection = dplyr::case_when(
      .data$event_direct == 1L ~ 1L, .data$track_negative ~ 0L, TRUE ~ NA_integer_
    ),
    scar_detection = dplyr::case_when(
      .data$event_scar == 1L ~ 1L, .data$track_negative ~ 0L, TRUE ~ NA_integer_
    ),
    high_density = as.integer(.data$density_incl_scars_ha > 0.02),
    direct_count = .data$direct_cots_count, scar_count,
    direct_localized = .data$track_negative | .data$event_direct == 1L,
    strict_negative = .data$track_negative,
    geometry_fallback, geometry_quality_reason,
    source_fraction = .data$fraction
  )
cull <- links$cull |>
  dplyr::transmute(
    visit_id, method = "cull", survey_date, year, reef_name, reef_key,
    cell_id, cell_support_m, projected_x, projected_y,
    effort_m = NA_real_, effort_minutes = .data$bottom_time,
    direct_detection = as.integer(.data$direct_cots_count > 0),
    scar_detection = NA_integer_,
    high_density = as.integer(
      .data$bottom_time > 0 & .data$direct_cots_count / .data$bottom_time >= 0.02
    ),
    direct_count = .data$direct_cots_count, scar_count = NA_real_,
    direct_localized = TRUE, strict_negative = .data$direct_cots_count == 0,
    geometry_fallback, geometry_quality_reason = "point_observation",
    source_fraction = 1
  )
rhis <- links$rhis |>
  dplyr::transmute(
    visit_id, method = "rhis", survey_date, year, reef_name, reef_key,
    cell_id, cell_support_m, projected_x, projected_y,
    effort_m = NA_real_, effort_minutes = NA_real_,
    direct_detection = dplyr::if_else(
      .data$rhis_cots_protocol_observed,
      as.integer(dplyr::coalesce(.data$adult_count, 0) > 0 |
        dplyr::coalesce(.data$juvenile_count, 0) > 0), NA_integer_
    ),
    scar_detection = dplyr::if_else(
      .data$rhis_cots_protocol_observed,
      as.integer(dplyr::coalesce(.data$branching_scar_count, 0) > 0 |
        dplyr::coalesce(.data$table_scar_count, 0) > 0), NA_integer_
    ),
    high_density = NA_integer_,
    direct_count = dplyr::coalesce(.data$adult_count, 0) + dplyr::coalesce(.data$juvenile_count, 0),
    scar_count = dplyr::coalesce(.data$branching_scar_count, 0) + dplyr::coalesce(.data$table_scar_count, 0),
    direct_localized = TRUE,
    strict_negative = .data$direct_detection == 0L & .data$scar_detection == 0L,
    geometry_fallback, geometry_quality_reason = "point_observation",
    source_fraction = 1, branching_live_pct, table_live_pct,
    branching_recent_dead_pct, table_recent_dead_pct, visibility, depth_m
  )

visits <- dplyr::bind_rows(manta, salad, cull, rhis) |>
  dplyr::filter(!is.na(.data$survey_date), !is.na(.data$cell_id)) |>
  annotate_support_cull_history()
aligned <- cluster_alignment_windows(visits, cfg$validation$alignment_window_days)
saveRDS(visits, file.path(cfg$paths$artifacts_dir, "cell_visits_200m.rds"))
saveRDS(aligned, file.path(cfg$paths$artifacts_dir, "cell_visits_200m_window_3d.rds"))
audit <- aligned |>
  dplyr::group_by(.data$alignment_id) |>
  dplyr::summarise(method_count = dplyr::n_distinct(.data$method), .groups = "drop") |>
  dplyr::summarise(
    support_m = support_m, alignment_units = dplyr::n(),
    multi_method_units = sum(.data$method_count > 1L),
    four_method_units = sum(.data$method_count == 4L)
  )
utils::write.csv(audit, file.path(cfg$paths$outputs_dir, "alignment_window_audit_200m.csv"), row.names = FALSE)
print(audit)
