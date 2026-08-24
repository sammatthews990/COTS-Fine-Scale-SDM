safe_max <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  max(x, na.rm = TRUE)
}

summarise_site_year <- function(observations, threshold, first_visit_only = FALSE) {
  working <- observations
  if (first_visit_only) {
    working <- working |>
      dplyr::group_by(.data$site_id, .data$year) |>
      dplyr::filter(.data$survey_date == min(.data$survey_date)) |>
      dplyr::ungroup()
  }

  working |>
    dplyr::group_by(
      .data$site_id, .data$reef_name, .data$site_name,
      .data$latitude, .data$longitude, .data$year
    ) |>
    dplyr::summarise(
      financial_year = dplyr::first(.data$financial_year),
      first_survey_date = min(.data$survey_date),
      last_survey_date = max(.data$survey_date),
      total_cots = sum(.data$total_cots, na.rm = TRUE),
      bottom_time = sum(.data$bottom_time, na.rm = TRUE),
      n_dives = dplyr::n(),
      n_positive_dives = sum(.data$total_cots > 0, na.rm = TRUE),
      n_high_density_dives = sum(.data$dive_cpue >= threshold, na.rm = TRUE),
      cpue_mean = .data$total_cots / .data$bottom_time,
      cpue_max_diagnostic = safe_max(.data$dive_cpue),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      site_year_id = paste(.data$site_id, .data$year, sep = "||"),
      presence = as.integer(.data$total_cots > 0),
      high_density = as.integer(.data$cpue_mean >= threshold),
      conditional_high_density = dplyr::if_else(
        .data$presence == 1L, .data$high_density, NA_integer_
      ),
      positive_dive_fraction = .data$n_positive_dives / .data$n_dives,
      high_density_dive_fraction = .data$n_high_density_dives / .data$n_dives,
      target_definition = dplyr::if_else(
        first_visit_only,
        "first surveyed date within site-year",
        "all dives within site-year"
      )
    ) |>
    dplyr::arrange(.data$reef_name, .data$site_name, .data$year)
}

summarise_site_history_diagnostics <- function(site_year, observations, threshold) {
  dive_history <- observations |>
    dplyr::group_by(.data$site_id) |>
    dplyr::summarise(
      lifetime_dives = dplyr::n(),
      lifetime_cpue_max = safe_max(.data$dive_cpue),
      ever_high_by_dive_max = as.integer(.data$lifetime_cpue_max >= threshold),
      .groups = "drop"
    )

  site_year |>
    dplyr::group_by(
      .data$site_id, .data$reef_name, .data$site_name,
      .data$latitude, .data$longitude
    ) |>
    dplyr::summarise(
      first_year = min(.data$year),
      last_year = max(.data$year),
      n_surveyed_years = dplyr::n_distinct(.data$year),
      n_high_density_years = sum(.data$high_density),
      high_density_year_fraction = mean(.data$high_density),
      any_high_density_year = as.integer(any(.data$high_density == 1L)),
      both_annual_states = any(.data$high_density == 1L) &&
        any(.data$high_density == 0L),
      .groups = "drop"
    ) |>
    dplyr::left_join(dive_history, by = "site_id")
}

add_outbreak_history_features <- function(site_year, outbreak) {
  current <- outbreak |>
    dplyr::transmute(
      .data$reef_name, .data$year,
      reef_outbreak_probability_current = .data$reef_outbreak_probability
    )

  lag_table <- site_year |>
    dplyr::select(dplyr::all_of(c("reef_name", "year"))) |>
    dplyr::distinct()

  for (lag_n in 1:3) {
    lag_values <- outbreak |>
      dplyr::transmute(
        .data$reef_name,
        year = .data$year + lag_n,
        value = .data$reef_outbreak_probability
      )
    names(lag_values)[names(lag_values) == "value"] <- paste0(
      "reef_outbreak_probability_lag", lag_n
    )
    lag_table <- dplyr::left_join(
      lag_table, lag_values, by = c("reef_name", "year")
    )
  }

  lag_table <- lag_table |>
    dplyr::mutate(
      reef_outbreak_probability_mean_lag1_3 = rowMeans(
        dplyr::pick(
          dplyr::all_of(c(
            "reef_outbreak_probability_lag1",
            "reef_outbreak_probability_lag2",
            "reef_outbreak_probability_lag3"
          ))
        ),
        na.rm = TRUE
      ),
      reef_outbreak_probability_max_lag1_3 = pmax(
        .data$reef_outbreak_probability_lag1,
        .data$reef_outbreak_probability_lag2,
        .data$reef_outbreak_probability_lag3,
        na.rm = TRUE
      ),
      reef_outbreak_probability_trend_lag1_3 =
        .data$reef_outbreak_probability_lag1 -
        .data$reef_outbreak_probability_lag3,
      outbreak_history_complete = as.integer(stats::complete.cases(dplyr::pick(
        dplyr::starts_with("reef_outbreak_probability_lag")
      )))
    )

  lag_table$reef_outbreak_probability_mean_lag1_3[
    !is.finite(lag_table$reef_outbreak_probability_mean_lag1_3)
  ] <- NA_real_
  lag_table$reef_outbreak_probability_max_lag1_3[
    !is.finite(lag_table$reef_outbreak_probability_max_lag1_3)
  ] <- NA_real_

  site_year |>
    dplyr::left_join(lag_table, by = c("reef_name", "year")) |>
    dplyr::left_join(current, by = c("reef_name", "year"))
}

attach_outbreak_forecast <- function(site_year, forecast) {
  forecast_features <- forecast |>
    dplyr::transmute(
      .data$reef_name,
      .data$year,
      reef_outbreak_probability_forecast =
        .data$reef_outbreak_probability,
      reef_outbreak_probability_cutoff =
        .data$outbreak_probability_cutoff
    )

  site_year |>
    dplyr::left_join(
      forecast_features,
      by = c("reef_name", "year")
    )
}
