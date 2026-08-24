assign_spatial_reef_folds <- function(data, k = 5L, seed = 1L) {
  required <- c("reef_name", "longitude", "latitude")
  assert_required_columns(data, required, "Spatial fold input")

  reefs <- data |>
    dplyr::group_by(.data$reef_name) |>
    dplyr::summarise(
      longitude = mean(.data$longitude, na.rm = TRUE),
      latitude = mean(.data$latitude, na.rm = TRUE),
      .groups = "drop"
    )

  if (nrow(reefs) < k) stop("Fewer reefs than requested spatial folds.")

  reef_sf <- sf::st_as_sf(
    reefs,
    coords = c("longitude", "latitude"),
    crs = 4326,
    remove = FALSE
  ) |>
    sf::st_transform(3577)
  xy <- sf::st_coordinates(reef_sf)
  xy_scaled <- scale(xy)

  set.seed(seed)
  clustering <- stats::kmeans(xy_scaled, centers = k, nstart = 100)
  reefs$spatial_fold <- as.integer(clustering$cluster)

  joined <- data |>
    dplyr::left_join(
      dplyr::select(reefs, "reef_name", "spatial_fold"),
      by = "reef_name"
    )

  if (any(is.na(joined$spatial_fold))) {
    stop("At least one observation did not receive a spatial fold.")
  }
  if (any((joined |>
    dplyr::select(dplyr::all_of(c("reef_name", "spatial_fold"))) |>
    dplyr::distinct() |>
    dplyr::count(.data$reef_name) |>
    dplyr::pull(.data$n)) != 1L)) {
    stop("A reef was assigned to more than one spatial fold.")
  }

  list(data = joined, reef_folds = reefs)
}

make_forward_year_split <- function(data, test_year) {
  train <- which(data$year < test_year)
  test <- which(data$year == test_year)
  if (length(train) == 0L || length(test) == 0L) {
    stop("Forward split has no training or testing rows for year ", test_year)
  }
  list(
    split = paste0("forward_", test_year),
    train = train,
    test = test,
    train_year_max = max(data$year[train]),
    test_year = test_year
  )
}

fold_balance_summary <- function(data, fold_col = "spatial_fold") {
  data |>
    dplyr::group_by(.data[[fold_col]]) |>
    dplyr::summarise(
      rows = dplyr::n(),
      reefs = dplyr::n_distinct(.data$reef_name),
      sites = dplyr::n_distinct(.data$site_id),
      positives = sum(.data$high_density),
      prevalence = mean(.data$high_density),
      min_year = min(.data$year),
      max_year = max(.data$year),
      .groups = "drop"
    )
}
