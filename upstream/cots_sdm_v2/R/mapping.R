v2_visual_regions <- function() {
  list(
    list(
      name = "lizard", label = "Lizard Island to Cooktown",
      south = -15.75, north = -14.50, west = 144.80, east = 146.10
    ),
    list(
      name = "cairns", label = "Cairns Reef Cluster",
      south = -17.30, north = -16.30, west = 145.60, east = 146.60
    ),
    list(
      name = "townsville", label = "Townsville Reef Cluster",
      south = -19.30, north = -18.20, west = 146.40, east = 147.70
    ),
    list(
      name = "whitsunday", label = "Whitsunday Reef Cluster",
      south = -20.50, north = -19.50, west = 148.50, east = 150.00
    )
  )
}

predict_xgb_raster <- function(
    x, model, filename, maximum_in_memory_cells = 5e6) {
  if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop("xgboost is required for raster prediction.")
  }
  feature_names <- model$feature_names
  if (!identical(names(x), feature_names)) {
    stop(
      "Raster features do not match model features. Expected: ",
      paste(feature_names, collapse = ", "), "; received: ",
      paste(names(x), collapse = ", ")
    )
  }
  if (terra::ncell(x) > maximum_in_memory_cells) {
    stop(
      "Visual prediction grid has ", format(terra::ncell(x), big.mark = ","),
      " cells; aggregate it below ", format(maximum_in_memory_cells, big.mark = ","),
      " cells before prediction."
    )
  }

  values <- terra::values(x, mat = TRUE)
  colnames(values) <- feature_names
  out <- rep(NA_real_, nrow(values))
  complete <- stats::complete.cases(values)
  if (any(complete)) {
    out[complete] <- as.numeric(stats::predict(
      model,
      xgboost::xgb.DMatrix(values[complete, , drop = FALSE])
    ))
  }
  result <- terra::setValues(x[[1]], out)
  names(result) <- "prediction"
  terra::writeRaster(
    result,
    filename,
    overwrite = TRUE,
    wopt = list(datatype = "FLT4S", gdal = c("COMPRESS=DEFLATE"))
  )
}

predict_additive_bridge_raster <- function(
    static_raster, outbreak_raster, bridge, filename) {
  if (!terra::compareGeom(
    static_raster, outbreak_raster, stopOnError = FALSE
  )) {
    stop("Static and outbreak rasters are not aligned.")
  }
  static_probability <- terra::values(static_raster, mat = FALSE)
  outbreak_probability <- terra::values(outbreak_raster, mat = FALSE)
  probability <- rep(NA_real_, length(static_probability))
  complete <- is.finite(static_probability) &
    is.finite(outbreak_probability)
  probability[complete] <- predict_additive_operational_bridge(
    bridge,
    static_probability[complete],
    outbreak_probability[complete]
  )
  result <- terra::setValues(static_raster, probability)
  names(result) <- "prediction"
  terra::writeRaster(
    result,
    filename,
    overwrite = TRUE,
    wopt = list(datatype = "FLT4S", gdal = c("COMPRESS=DEFLATE"))
  )
}

make_forecast_voronoi_raster <- function(
    forecast, template, region, buffer_degrees = 1.0) {
  local <- forecast |>
    dplyr::filter(
      .data$reef_longitude >= region$west - buffer_degrees,
      .data$reef_longitude <= region$east + buffer_degrees,
      .data$reef_latitude >= region$south - buffer_degrees,
      .data$reef_latitude <= region$north + buffer_degrees
    ) |>
    dplyr::transmute(
      longitude = .data$reef_longitude,
      latitude = .data$reef_latitude,
      reef_outbreak_probability_operational =
        .data$reef_outbreak_probability
    )
  if (nrow(local) < 3L) {
    stop("Fewer than three forecast reef centroids surround region ", region$name)
  }

  points <- terra::vect(
    local,
    geom = c("longitude", "latitude"),
    crs = "EPSG:4326"
  ) |>
    terra::project(terra::crs(template))
  polygons <- terra::voronoi(points) |>
    terra::crop(terra::ext(template))
  outbreak <- terra::rasterize(
    polygons,
    template,
    field = "reef_outbreak_probability_operational",
    touches = TRUE
  )
  outbreak <- terra::mask(outbreak, template)
  names(outbreak) <- "reef_outbreak_probability_operational"
  outbreak
}

align_raster_to_template <- function(raster, template) {
  cropped <- terra::crop(raster, terra::ext(template), snap = "out")
  if (!terra::compareGeom(cropped, template, stopOnError = FALSE)) {
    cropped <- terra::resample(cropped, template, method = "bilinear")
  }
  terra::mask(cropped, template)
}

project_layers_to_wgs84 <- function(layers, native_template) {
  wgs_template <- terra::project(native_template, "EPSG:4326")
  lapply(layers, function(x) {
    if (terra::compareGeom(x, native_template, stopOnError = FALSE)) {
      terra::project(x, wgs_template, method = "bilinear")
    } else {
      terra::project(
        terra::resample(x, native_template, method = "bilinear"),
        wgs_template,
        method = "bilinear"
      )
    }
  })
}

raster_to_gray_png <- function(raster, filename) {
  values <- terra::values(raster, mat = FALSE)
  valid <- is.finite(values)
  rows <- nrow(raster)
  cols <- ncol(raster)
  image <- array(0, dim = c(rows, cols, 2L))
  grayscale <- rep(0, length(values))
  grayscale[valid] <- pmin(pmax(values[valid], 0), 1)
  image[, , 1L] <- matrix(grayscale, nrow = rows, byrow = TRUE)
  alpha <- rep(0, length(values))
  alpha[valid] <- 1
  image[, , 2L] <- matrix(alpha, nrow = rows, byrow = TRUE)
  png::writePNG(image, filename)
  invisible(filename)
}

raster_difference_to_gray_png <- function(raster, filename, limit) {
  scaled <- (raster + limit) / (2 * limit)
  raster_to_gray_png(terra::clamp(scaled, 0, 1), filename)
}

point_feature_collection <- function(
    data, longitude = "longitude", latitude = "latitude", properties) {
  valid <- is.finite(data[[longitude]]) & is.finite(data[[latitude]])
  data <- data[valid, , drop = FALSE]
  features <- lapply(seq_len(nrow(data)), function(i) {
    property_values <- as.list(data[i, properties, drop = FALSE])
    property_values <- lapply(property_values, function(x) {
      if (length(x) == 0L || is.na(x)) NULL else unname(x)
    })
    list(
      type = "Feature",
      properties = property_values,
      geometry = list(
        type = "Point",
        coordinates = c(data[[longitude]][i], data[[latitude]][i])
      )
    )
  })
  list(type = "FeatureCollection", features = features)
}

raster_layer_summary <- function(raster, region, layer) {
  quantiles <- as.numeric(terra::global(
    raster,
    fun = function(x) stats::quantile(
      x, probs = c(0, 0.05, 0.5, 0.95, 1), na.rm = TRUE
    )
  ))
  tibble::tibble(
    region = region,
    layer = layer,
    valid_cells = as.numeric(terra::global(
      !is.na(raster), "sum", na.rm = TRUE
    )[1, 1]),
    minimum = quantiles[1],
    q05 = quantiles[2],
    median = quantiles[3],
    q95 = quantiles[4],
    maximum = quantiles[5]
  )
}
