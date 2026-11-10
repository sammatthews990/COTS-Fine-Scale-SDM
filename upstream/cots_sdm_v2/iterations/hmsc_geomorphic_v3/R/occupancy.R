make_unmarked_frame <- function(histories, site_covariates, observation_covariates) {
  if (!requireNamespace("unmarked", quietly = TRUE)) stop("Package unmarked is required.")
  unmarked::unmarkedFrameOccu(
    y = histories,
    siteCovs = site_covariates,
    obsCovs = observation_covariates
  )
}

wide_observation_matrix <- function(data, value, sites, maximum_visits = 12L) {
  x <- data |>
    dplyr::arrange(.data$site_id, .data$survey_date, .data$method) |>
    dplyr::group_by(.data$site_id) |>
    dplyr::slice_head(n = maximum_visits) |>
    dplyr::mutate(visit_number = dplyr::row_number()) |>
    dplyr::ungroup() |>
    dplyr::select("site_id", "visit_number", dplyr::all_of(value)) |>
    tidyr::pivot_wider(
      names_from = "visit_number", values_from = dplyr::all_of(value),
      names_prefix = "visit_"
    )
  result <- dplyr::left_join(tibble::tibble(site_id = sites), x, by = "site_id")
  as.matrix(result[, setdiff(names(result), "site_id"), drop = FALSE])
}
