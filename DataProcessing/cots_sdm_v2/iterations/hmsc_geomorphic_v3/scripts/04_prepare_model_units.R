args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v3_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v3_dir, "R", "bootstrap.R"))
ctx <- initialise_v3(script_file)
cfg <- ctx$config

observed_any <- function(x) {
  if (all(is.na(x))) return(NA_integer_)
  as.integer(any(x == 1L, na.rm = TRUE))
}

prepare_support <- function(cell_size_m, window_days = 3L, zero_ratio = 5L) {
  message("Preparing model units at ", cell_size_m, " m.")
  visits <- readRDS(file.path(
    cfg$paths$artifacts_dir,
    paste0("cell_visits_", cell_size_m, "m_window_", window_days, "d.rds")
  ))
  links <- readRDS(file.path(
    cfg$paths$artifacts_dir, paste0("footprint_links_", cell_size_m, "m.rds")
  ))

  rhis_detail <- links$rhis |>
    dplyr::transmute(
      visit_id, cell_id,
      adult_detection = dplyr::if_else(
        .data$rhis_cots_protocol_observed,
        as.integer(dplyr::coalesce(.data$adult_count, 0) > 0), NA_integer_
      ),
      juvenile_detection = dplyr::if_else(
        .data$rhis_cots_protocol_observed,
        as.integer(dplyr::coalesce(.data$juvenile_count, 0) > 0), NA_integer_
      ),
      branching_scar_detection = dplyr::if_else(
        .data$rhis_cots_protocol_observed,
        as.integer(dplyr::coalesce(.data$branching_scar_count, 0) > 0), NA_integer_
      ),
      table_scar_detection = dplyr::if_else(
        .data$rhis_cots_protocol_observed,
        as.integer(dplyr::coalesce(.data$table_scar_count, 0) > 0), NA_integer_
      )
    )
  visits <- dplyr::left_join(visits, rhis_detail, by = c("visit_id", "cell_id"))

  event_visits <- visits |>
    dplyr::filter(.data$ecological_eligible, !is.na(.data$direct_detection)) |>
    dplyr::group_by(
      .data$cell_id, .data$year, .data$survey_date, .data$method
    ) |>
    dplyr::summarise(
      reef_key = dplyr::first(.data$reef_key),
      reef_name = dplyr::first(.data$reef_name),
      projected_x = dplyr::first(.data$projected_x),
      projected_y = dplyr::first(.data$projected_y),
      direct_detection = observed_any(.data$direct_detection),
      scar_detection = observed_any(.data$scar_detection),
      effort_m = sum(.data$effort_m, na.rm = TRUE),
      effort_minutes = sum(.data$effort_minutes, na.rm = TRUE),
      geometry_fallback = any(.data$geometry_fallback),
      .groups = "drop"
    ) |>
    dplyr::mutate(site_id = paste(.data$cell_id, .data$year, sep = "||"))

  site_summary <- event_visits |>
    dplyr::group_by(.data$site_id) |>
    dplyr::summarise(
      cell_id = dplyr::first(.data$cell_id), year = dplyr::first(.data$year),
      reef_key = dplyr::first(.data$reef_key), reef_name = dplyr::first(.data$reef_name),
      projected_x = dplyr::first(.data$projected_x),
      projected_y = dplyr::first(.data$projected_y),
      n_visits = dplyr::n(), positive = any(.data$direct_detection == 1L),
      methods = dplyr::n_distinct(.data$method), .groups = "drop"
    ) |>
    dplyr::filter(!is.na(.data$reef_key), .data$reef_key != "")
  positive_counts <- site_summary |>
    dplyr::group_by(.data$reef_key, .data$year) |>
    dplyr::summarise(n_positive = sum(.data$positive), .groups = "drop") |>
    dplyr::filter(.data$n_positive > 0)
  eligible_sites <- dplyr::inner_join(
    site_summary, positive_counts, by = c("reef_key", "year")
  )
  positive_sites <- dplyr::filter(eligible_sites, .data$positive)
  zero_sites <- eligible_sites |>
    dplyr::filter(!.data$positive) |>
    dplyr::arrange(.data$reef_key, .data$year, .data$cell_id) |>
    dplyr::group_by(.data$reef_key, .data$year) |>
    dplyr::filter(dplyr::row_number() <= max(
      1L, zero_ratio * dplyr::first(.data$n_positive)
    )) |>
    dplyr::ungroup()
  selected_sites <- dplyr::bind_rows(positive_sites, zero_sites) |>
    dplyr::distinct(.data$site_id, .keep_all = TRUE)
  occupancy_visits <- dplyr::semi_join(
    event_visits, selected_sites, by = "site_id"
  ) |>
    dplyr::arrange(.data$site_id, .data$survey_date, .data$method)

  aligned <- visits |>
    dplyr::filter(.data$ecological_eligible) |>
    dplyr::group_by(.data$alignment_id) |>
    dplyr::summarise(
      cell_id = dplyr::first(.data$cell_id), year = dplyr::first(.data$year),
      reef_key = dplyr::first(.data$reef_key), reef_name = dplyr::first(.data$reef_name),
      projected_x = dplyr::first(.data$projected_x),
      projected_y = dplyr::first(.data$projected_y),
      cull_high_density = observed_any(.data$high_density[.data$method == "cull"]),
      cull_direct = observed_any(.data$direct_detection[.data$method == "cull"]),
      salad_high_density = observed_any(.data$high_density[.data$method == "salad"]),
      salad_direct = observed_any(.data$direct_detection[.data$method == "salad"]),
      rhis_adult = observed_any(.data$adult_detection),
      rhis_juvenile = observed_any(.data$juvenile_detection),
      rhis_direct = observed_any(.data$direct_detection[.data$method == "rhis"]),
      rhis_branching_scar = observed_any(.data$branching_scar_detection),
      rhis_table_scar = observed_any(.data$table_scar_detection),
      rhis_branching_live_pct = mean(.data$branching_live_pct, na.rm = TRUE),
      rhis_table_live_pct = mean(.data$table_live_pct, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(dplyr::across(
      c("rhis_branching_live_pct", "rhis_table_live_pct"),
      ~ dplyr::if_else(is.nan(.x), NA_real_, .x)
    )) |>
    dplyr::filter(dplyr::if_any(
      c("cull_high_density", "cull_direct", "salad_high_density", "salad_direct",
        "rhis_adult", "rhis_juvenile", "rhis_direct"), ~ !is.na(.x)
    ))

  saveRDS(occupancy_visits, file.path(
    cfg$paths$artifacts_dir, paste0("occupancy_visits_", cell_size_m, "m.rds")
  ))
  saveRDS(selected_sites, file.path(
    cfg$paths$artifacts_dir, paste0("occupancy_sites_", cell_size_m, "m.rds")
  ))
  saveRDS(aligned, file.path(
    cfg$paths$artifacts_dir, paste0("hmsc_core_units_", cell_size_m, "m.rds")
  ))

  all_cells <- dplyr::bind_rows(
    dplyr::select(selected_sites, "cell_id", "projected_x", "projected_y"),
    dplyr::select(aligned, "cell_id", "projected_x", "projected_y")
  ) |>
    dplyr::distinct(.data$cell_id, .keep_all = TRUE) |>
    dplyr::mutate(cell_support_m = cell_size_m)
  saveRDS(all_cells, file.path(
    cfg$paths$artifacts_dir, paste0("predictor_cells_", cell_size_m, "m.rds")
  ))

  tibble::tibble(
    support_m = cell_size_m,
    occupancy_sites = nrow(selected_sites),
    occupancy_positive_sites = sum(selected_sites$positive),
    occupancy_visit_rows = nrow(occupancy_visits),
    repeated_occupancy_sites = sum(selected_sites$n_visits > 1),
    hmsc_units = nrow(aligned), predictor_cells = nrow(all_cells)
  )
}

audit <- dplyr::bind_rows(lapply(cfg$validation$cell_support_m, prepare_support))
utils::write.csv(audit, file.path(
  cfg$paths$outputs_dir, "model_unit_audit.csv"
), row.names = FALSE)
print(audit)
