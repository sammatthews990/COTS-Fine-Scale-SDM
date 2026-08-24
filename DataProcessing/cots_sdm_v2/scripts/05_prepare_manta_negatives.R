args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v2_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v2_dir, "R", "bootstrap.R"))
ctx <- initialise_v2(script_file)
cfg <- ctx$config

cat("Reading long-history manta observations...\n")
manta_history <- read_manta_observations(cfg$paths$manta_history_csv)
cat("Reading current manta observations for contradiction checks...\n")
manta_current <- read_manta_observations(
  cfg$paths$manta_current_workbook,
  sheet = "Manta"
)

reef_evidence <- summarise_manta_reef_evidence(manta_history)
classified <- classify_manta_negative_reefs(
  reef_evidence,
  manta_current,
  minimum_tows = cfg$manta$minimum_tows,
  minimum_survey_years = cfg$manta$minimum_survey_years
)
strong_reefs <- classified |>
  dplyr::filter(.data$strong_negative_evidence)

cull <- if (file.exists(file.path(
  cfg$paths$artifacts_dir, "cull_observations_clean.rds"
))) {
  readRDS(file.path(
    cfg$paths$artifacts_dir, "cull_observations_clean.rds"
  ))
} else {
  read_cull_observations(cfg$paths$cull_workbook)
}
cull_check <- cull |>
  dplyr::group_by(.data$reef_name) |>
  dplyr::summarise(
    cull_dives = dplyr::n(),
    cull_positive_dives = sum(.data$total_cots > 0),
    cull_total_cots = sum(.data$total_cots),
    first_cull_year = min(.data$year),
    last_cull_year = max(.data$year),
    .groups = "drop"
  )
strong_reefs <- strong_reefs |>
  dplyr::left_join(cull_check, by = "reef_name") |>
  dplyr::mutate(
    cull_dives = dplyr::coalesce(.data$cull_dives, 0L),
    cull_positive_dives = dplyr::coalesce(
      .data$cull_positive_dives, 0L
    ),
    cull_total_cots = dplyr::coalesce(.data$cull_total_cots, 0),
    cross_method_status = dplyr::case_when(
      .data$cull_positive_dives > 0L ~ "manta_cull_conflict",
      .data$cull_dives > 0L ~ "cross_method_zero_supported",
      TRUE ~ "manta_only_uncontradicted"
    ),
    static_prior_negative_candidate =
      .data$cull_positive_dives == 0L
  )

forecast <- read_outbreak_probabilities(
  cfg$paths$outbreak_forecast_file, "adjusted_forecast"
) |>
  dplyr::select(dplyr::all_of(c(
    "reef_name", "reef_longitude", "reef_latitude",
    "reef_outbreak_probability", "outbreak_probability_cutoff"
  ))) |>
  dplyr::rename(
    reef_outbreak_probability_2026 = "reef_outbreak_probability"
  )
strong_reefs <- strong_reefs |>
  dplyr::left_join(forecast, by = "reef_name")

strong_tows <- manta_history |>
  dplyr::filter(
    .data$reef_name %in% strong_reefs$reef_name,
    .data$strict_negative
  ) |>
  dplyr::left_join(
    strong_reefs |>
      dplyr::select(dplyr::all_of(c(
        "reef_name", "cross_method_status",
        "static_prior_negative_candidate"
      ))),
    by = "reef_name"
  ) |>
  dplyr::mutate(
    adequate_tow_distance = is.na(.data$tow_distance_m) |
      .data$tow_distance_m >= cfg$manta$minimum_tow_distance_m,
    fine_scale_eligible = .data$valid_coordinates &
      .data$hard_coral_observed & .data$adequate_tow_distance,
    evidence_role = "presence_stage_strong_negative"
  )

reef_year_evidence <- strong_tows |>
  dplyr::group_by(.data$reef_name, .data$year) |>
  dplyr::summarise(
    negative_tows = dplyr::n(),
    total_tow_distance_m = sum(.data$tow_distance_m, na.rm = TRUE),
    tow_distance_recorded = sum(!is.na(.data$tow_distance_m)),
    coordinate_tows = sum(.data$valid_coordinates),
    fine_scale_eligible_tows = sum(.data$fine_scale_eligible),
    cross_method_status = dplyr::first(.data$cross_method_status),
    static_prior_negative_candidate = dplyr::first(
      .data$static_prior_negative_candidate
    ),
    evidence_role = "presence_stage_strong_negative",
    .groups = "drop"
  )

sensitivity <- manta_negative_sensitivity(reef_evidence)

summary <- tibble::tibble(
  manta_rows = nrow(manta_history),
  manta_reefs = dplyr::n_distinct(manta_history$reef_name),
  minimum_tows = cfg$manta$minimum_tows,
  minimum_survey_years = cfg$manta$minimum_survey_years,
  strong_negative_reefs = nrow(strong_reefs),
  strong_negative_tows = nrow(strong_tows),
  fine_scale_eligible_tows = sum(strong_tows$fine_scale_eligible),
  manta_only_uncontradicted_reefs = sum(
    strong_reefs$cross_method_status == "manta_only_uncontradicted"
  ),
  cross_method_zero_supported_reefs = sum(
    strong_reefs$cross_method_status == "cross_method_zero_supported"
  ),
  manta_cull_conflict_reefs = sum(
    strong_reefs$cross_method_status == "manta_cull_conflict"
  ),
  static_prior_negative_candidate_reefs = sum(
    strong_reefs$static_prior_negative_candidate
  ),
  contradicted_by_current_workbook = sum(
    classified$candidate_before_current_check &
      classified$contradicted_by_current_data
  )
)

saveRDS(
  manta_history,
  file.path(cfg$paths$artifacts_dir, "manta_observations_clean.rds")
)
saveRDS(
  reef_year_evidence,
  file.path(cfg$paths$artifacts_dir, "manta_strong_negative_reef_years.rds")
)
saveRDS(
  strong_tows,
  file.path(cfg$paths$artifacts_dir, "manta_strong_negative_tows.rds")
)
utils::write.csv(
  classified,
  file.path(cfg$paths$outputs_dir, "manta_reef_evidence.csv"),
  row.names = FALSE
)
utils::write.csv(
  strong_reefs,
  file.path(cfg$paths$outputs_dir, "manta_strong_negative_reefs.csv"),
  row.names = FALSE
)
utils::write.csv(
  reef_year_evidence,
  file.path(cfg$paths$outputs_dir, "manta_strong_negative_reef_years.csv"),
  row.names = FALSE
)
utils::write.csv(
  strong_tows,
  file.path(cfg$paths$outputs_dir, "manta_strong_negative_tows.csv"),
  row.names = FALSE
)
utils::write.csv(
  sensitivity,
  file.path(cfg$paths$outputs_dir, "manta_negative_sensitivity.csv"),
  row.names = FALSE
)
utils::write.csv(
  summary,
  file.path(cfg$paths$outputs_dir, "manta_negative_summary.csv"),
  row.names = FALSE
)

cat("Manta negative-evidence preparation complete.\n")
print(summary)
print(sensitivity)
