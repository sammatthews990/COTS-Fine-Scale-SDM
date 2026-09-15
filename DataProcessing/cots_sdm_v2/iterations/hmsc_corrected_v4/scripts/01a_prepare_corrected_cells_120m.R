args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v4_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v4_dir, "R", "bootstrap.R"))
ctx <- initialise_v4(script_file)
cfg <- ctx$config

links <- readRDS(file.path(cfg$paths$v3_artifacts_dir, "footprint_links_120m.rds"))
visits <- readRDS(file.path(
  cfg$paths$v3_artifacts_dir, "cell_visits_120m_window_3d.rds"
))
predictors <- readRDS(file.path(
  cfg$paths$v3_artifacts_dir, "support_predictors_120m.rds"
))

manta_visit <- links$manta |>
  dplyr::distinct(.data$visit_id, .keep_all = TRUE) |>
  classify_manta_evidence(cfg)
local_activity <- visits |>
  dplyr::filter(.data$ecological_eligible, .data$method != "manta") |>
  dplyr::mutate(
    positive = dplyr::coalesce(.data$direct_detection == 1L, FALSE) |
      dplyr::coalesce(.data$scar_detection == 1L, FALSE) |
      dplyr::coalesce(.data$high_density == 1L, FALSE)
  ) |>
  dplyr::group_by(.data$reef_key, .data$year) |>
  dplyr::summarise(
    localized_positive = any(.data$positive),
    localized_positive_visits = dplyr::n_distinct(.data$visit_id[.data$positive]),
    .groups = "drop"
  )
manta_activity <- manta_visit |>
  dplyr::group_by(.data$reef_key, .data$year) |>
  dplyr::summarise(
    manta_strong_positive = any(.data$manta_strong_positive),
    manta_weak_scar = any(.data$manta_weak_scar),
    manta_strong_positive_tows = sum(.data$manta_strong_positive),
    manta_weak_scar_tows = sum(.data$manta_weak_scar),
    .groups = "drop"
  )
activity <- dplyr::full_join(
  local_activity, manta_activity, by = c("reef_key", "year")
) |>
  dplyr::mutate(
    dplyr::across(
      c("localized_positive", "manta_strong_positive", "manta_weak_scar"),
      ~ dplyr::coalesce(.x, FALSE)
    ),
    dplyr::across(
      c("localized_positive_visits", "manta_strong_positive_tows", "manta_weak_scar_tows"),
      ~ dplyr::coalesce(.x, 0L)
    ),
    active_primary = .data$localized_positive | .data$manta_strong_positive,
    active_weak_sensitivity = .data$active_primary | .data$manta_weak_scar
  )

rhis_detail <- links$rhis |>
  dplyr::transmute(
    visit_id, cell_id,
    rhis_adult = dplyr::if_else(
      .data$rhis_cots_protocol_observed,
      as.integer(dplyr::coalesce(.data$adult_count, 0) > 0), NA_integer_
    ),
    rhis_juvenile = dplyr::if_else(
      .data$rhis_cots_protocol_observed,
      as.integer(dplyr::coalesce(.data$juvenile_count, 0) > 0), NA_integer_
    ),
    rhis_branching_scar = dplyr::if_else(
      .data$rhis_cots_protocol_observed,
      as.integer(dplyr::coalesce(.data$branching_scar_count, 0) > 0), NA_integer_
    ),
    rhis_table_scar = dplyr::if_else(
      .data$rhis_cots_protocol_observed,
      as.integer(dplyr::coalesce(.data$table_scar_count, 0) > 0), NA_integer_
    )
  )
salad_detail <- links$salad |>
  dplyr::transmute(
    visit_id, cell_id,
    salad_density_incl_scars = as.numeric(.data$density_incl_scars_ha)
  )
enriched <- visits |>
  dplyr::left_join(rhis_detail, by = c("visit_id", "cell_id")) |>
  dplyr::left_join(salad_detail, by = c("visit_id", "cell_id"))

cell_units_all <- enriched |>
  dplyr::filter(.data$ecological_eligible, .data$method != "manta") |>
  dplyr::group_by(.data$alignment_id) |>
  dplyr::summarise(
    cell_id = dplyr::first(.data$cell_id), year = dplyr::first(.data$year),
    reef_key = dplyr::first(.data$reef_key), reef_name = dplyr::first(.data$reef_name),
    projected_x = dplyr::first(.data$projected_x),
    projected_y = dplyr::first(.data$projected_y),
    cull_high_density = observed_any_v4(.data$high_density[.data$method == "cull"]),
    cull_direct = observed_any_v4(.data$direct_detection[.data$method == "cull"]),
    salad_high_density = observed_any_v4(.data$high_density[.data$method == "salad"]),
    salad_direct = observed_any_v4(.data$direct_detection[.data$method == "salad"]),
    salad_scar = observed_any_v4(.data$scar_detection[.data$method == "salad"]),
    salad_log_density = {
      z <- .data$salad_density_incl_scars[.data$method == "salad"]
      if (!length(z) || all(is.na(z))) NA_real_ else log1p(mean(z, na.rm = TRUE))
    },
    rhis_adult = observed_any_v4(.data$rhis_adult),
    rhis_juvenile = observed_any_v4(.data$rhis_juvenile),
    rhis_branching_scar = observed_any_v4(.data$rhis_branching_scar),
    rhis_table_scar = observed_any_v4(.data$rhis_table_scar),
    rhis_branching_live_asin = {
      z <- .data$branching_live_pct[.data$method == "rhis"]
      if (!length(z) || all(is.na(z))) NA_real_ else asin(sqrt(pmin(pmax(mean(z, na.rm = TRUE), 0), 100) / 100))
    },
    rhis_table_live_asin = {
      z <- .data$table_live_pct[.data$method == "rhis"]
      if (!length(z) || all(is.na(z))) NA_real_ else asin(sqrt(pmin(pmax(mean(z, na.rm = TRUE), 0), 100) / 100))
    },
    methods = paste(sort(unique(.data$method)), collapse = ";"),
    .groups = "drop"
  ) |>
  dplyr::filter(!is.na(.data$reef_key), .data$reef_key != "") |>
  dplyr::left_join(activity, by = c("reef_key", "year")) |>
  dplyr::mutate(
    active_primary = dplyr::coalesce(.data$active_primary, FALSE),
    active_weak_sensitivity = dplyr::coalesce(.data$active_weak_sensitivity, FALSE),
    reef_year_id = paste(.data$reef_key, .data$year, sep = "||"),
    observation_support = "localized_120m_cell"
  )
response_names <- c(
  "cull_high_density", "cull_direct", "salad_high_density", "salad_direct",
  "salad_scar", "salad_log_density", "rhis_adult", "rhis_juvenile",
  "rhis_branching_scar", "rhis_table_scar", "rhis_branching_live_asin",
  "rhis_table_live_asin"
)
cell_units_all <- cell_units_all |>
  dplyr::filter(rowSums(!is.na(dplyr::pick(dplyr::all_of(response_names)))) > 0L)
cell_units <- cell_units_all |>
  dplyr::filter(.data$active_primary) |>
  dplyr::left_join(predictors, by = c("cell_id", "projected_x", "projected_y"))
raw_terms <- unlist(terrain_families(names(cell_units)), use.names = FALSE)
clean_terms <- c("AhyaD", "Aspat", "Aten", raw_terms)
enhanced_terms <- c("geomorphic_class", clean_terms, "turbid", "waves_Hs", "waves_Tp")
cell_units <- cell_units |>
  dplyr::mutate(
    clean_complete = stats::complete.cases(dplyr::pick(dplyr::all_of(clean_terms))),
    enhanced_complete = stats::complete.cases(dplyr::pick(dplyr::all_of(enhanced_terms)))
  )

saveRDS(activity, file.path(cfg$paths$artifacts_dir, "reef_year_activity_120m.rds"))
saveRDS(cell_units_all, file.path(cfg$paths$artifacts_dir, "corrected_cell_units_all_120m.rds"))
saveRDS(cell_units, file.path(cfg$paths$artifacts_dir, "corrected_cell_units_120m.rds"))
audit <- tibble::tibble(
  support_m = 120L, all_localized_units = nrow(cell_units_all),
  primary_active_units = sum(cell_units_all$active_primary),
  excluded_inactive_units = sum(!cell_units_all$active_primary),
  primary_active_reef_years = dplyr::n_distinct(cell_units_all$reef_year_id[cell_units_all$active_primary]),
  inactive_reef_years = dplyr::n_distinct(cell_units_all$reef_year_id[!cell_units_all$active_primary]),
  fitted_cell_units = nrow(cell_units),
  clean_complete_cell_units = sum(cell_units$clean_complete),
  enhanced_complete_cell_units = sum(cell_units$enhanced_complete)
)
utils::write.csv(activity, file.path(cfg$paths$outputs_dir, "reef_year_activity_audit_120m.csv"), row.names = FALSE)
utils::write.csv(audit, file.path(cfg$paths$outputs_dir, "corrected_unit_audit_120m.csv"), row.names = FALSE)
print(audit)
