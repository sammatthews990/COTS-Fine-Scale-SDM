extract_static_site_predictors <- function(site_year, predictor_path) {
  if (!file.exists(predictor_path)) {
    stop("Predictor stack not found: ", predictor_path)
  }

  predictors <- terra::rast(predictor_path)
  predictor_names <- names(predictors)

  sites <- site_year |>
    dplyr::select(dplyr::all_of(c(
      "site_id", "reef_name", "site_name", "longitude", "latitude"
    ))) |>
    dplyr::distinct() |>
    dplyr::mutate(.extract_id = dplyr::row_number())

  site_sf <- sf::st_as_sf(
    sites,
    coords = c("longitude", "latitude"),
    crs = 4326,
    remove = FALSE
  ) |>
    sf::st_transform(terra::crs(predictors))

  extracted <- terra::extract(
    predictors,
    terra::vect(site_sf),
    ID = TRUE
  )
  names(extracted)[1] <- ".extract_id"

  site_predictors <- sites |>
    dplyr::left_join(extracted, by = ".extract_id") |>
    dplyr::mutate(
      predictor_complete = as.integer(stats::complete.cases(dplyr::pick(
        dplyr::all_of(predictor_names)
      )))
    )

  joined <- site_year |>
    dplyr::left_join(
      dplyr::select(
        site_predictors,
        dplyr::all_of(c("site_id", ".extract_id")),
        dplyr::all_of(predictor_names),
        "predictor_complete"
      ),
      by = "site_id"
    )

  list(
    data = joined,
    sites = site_predictors,
    predictor_names = predictor_names,
    raster_crs = terra::crs(predictors),
    raster_resolution = terra::res(predictors)
  )
}

attach_existing_site_predictors <- function(site_year, extraction) {
  predictor_names <- extraction$predictor_names
  site_year |>
    dplyr::left_join(
      dplyr::select(
        extraction$sites,
        dplyr::all_of(c("site_id", ".extract_id")),
        dplyr::all_of(predictor_names),
        "predictor_complete"
      ),
      by = "site_id"
    )
}
