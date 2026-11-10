summarise_manta_reef_evidence <- function(manta) {
  manta |>
    dplyr::group_by(.data$reef_name) |>
    dplyr::summarise(
      manta_tows = dplyr::n(),
      informative_tows = sum(.data$informative),
      strict_negative_tows = sum(.data$strict_negative),
      positive_evidence_tows = sum(.data$positive_evidence),
      survey_years = dplyr::n_distinct(.data$year),
      first_survey_year = min(.data$year),
      last_survey_year = max(.data$year),
      temporal_span_years = .data$last_survey_year - .data$first_survey_year,
      coordinate_tows = sum(.data$valid_coordinates),
      hard_coral_tows = sum(.data$hard_coral_observed),
      all_tows_informative = all(.data$informative),
      never_positive = !any(.data$positive_evidence),
      all_tows_strict_negative = all(.data$strict_negative),
      .groups = "drop"
    )
}

classify_manta_negative_reefs <- function(
    reef_evidence, current_manta, minimum_tows, minimum_survey_years) {
  current_check <- current_manta |>
    dplyr::group_by(.data$reef_name) |>
    dplyr::summarise(
      current_dataset_tows = dplyr::n(),
      current_dataset_last_year = max(.data$year),
      current_dataset_positive_tows = sum(.data$positive_evidence),
      .groups = "drop"
    )

  reef_evidence |>
    dplyr::left_join(current_check, by = "reef_name") |>
    dplyr::mutate(
      current_dataset_tows = dplyr::coalesce(
        .data$current_dataset_tows, 0L
      ),
      current_dataset_positive_tows = dplyr::coalesce(
        .data$current_dataset_positive_tows, 0L
      ),
      contradicted_by_current_data =
        .data$current_dataset_positive_tows > 0L,
      candidate_before_current_check =
        .data$all_tows_informative &
        .data$all_tows_strict_negative &
        .data$manta_tows >= minimum_tows &
        .data$survey_years >= minimum_survey_years,
      strong_negative_evidence =
        .data$candidate_before_current_check &
        !.data$contradicted_by_current_data,
      evidence_role = dplyr::if_else(
        .data$strong_negative_evidence,
        "presence_stage_strong_negative",
        "not_selected"
      )
    )
}

manta_negative_sensitivity <- function(reef_evidence) {
  thresholds <- tidyr::expand_grid(
    minimum_tows = c(1L, 5L, 10L, 20L, 50L),
    minimum_survey_years = c(1L, 2L, 3L, 5L, 10L)
  ) |>
    dplyr::filter(
      (.data$minimum_tows == 1L & .data$minimum_survey_years == 1L) |
      (.data$minimum_tows == 5L & .data$minimum_survey_years == 2L) |
      (.data$minimum_tows == 10L & .data$minimum_survey_years == 3L) |
      (.data$minimum_tows == 20L & .data$minimum_survey_years == 5L) |
      (.data$minimum_tows == 50L & .data$minimum_survey_years == 10L)
    )

  thresholds |>
    dplyr::rowwise() |>
    dplyr::mutate(
      reefs = sum(
        reef_evidence$all_tows_informative &
          reef_evidence$all_tows_strict_negative &
          reef_evidence$manta_tows >= .data$minimum_tows &
          reef_evidence$survey_years >= .data$minimum_survey_years
      ),
      tows = sum(reef_evidence$manta_tows[
        reef_evidence$all_tows_informative &
          reef_evidence$all_tows_strict_negative &
          reef_evidence$manta_tows >= .data$minimum_tows &
          reef_evidence$survey_years >= .data$minimum_survey_years
      ])
    ) |>
    dplyr::ungroup()
}
