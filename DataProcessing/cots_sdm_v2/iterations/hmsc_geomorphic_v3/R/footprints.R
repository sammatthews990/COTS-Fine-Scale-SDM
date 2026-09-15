valid_lonlat <- function(longitude, latitude) {
  is.finite(longitude) & is.finite(latitude) &
    dplyr::between(longitude, 140, 160) & dplyr::between(latitude, -30, -5)
}

track_geometry_audit <- function(data, method) {
  required <- c("start_longitude", "start_latitude", "end_longitude", "end_latitude")
  if (!all(required %in% names(data))) stop(method, " track endpoints are absent.")
  good_start <- valid_lonlat(data$start_longitude, data$start_latitude)
  good_end <- valid_lonlat(data$end_longitude, data$end_latitude)
  tibble::tibble(
    method = method,
    records = nrow(data),
    valid_start = sum(good_start),
    valid_end = sum(good_end),
    valid_line = sum(good_start & good_end),
    point_fallback = sum(good_start & !good_end),
    unusable_geometry = sum(!good_start),
    line_fraction = mean(good_start & good_end)
  )
}

make_track_sf <- function(data, crs, id = "visit_id") {
  good_start <- valid_lonlat(data$start_longitude, data$start_latitude)
  good_end <- valid_lonlat(data$end_longitude, data$end_latitude)
  keep <- good_start & good_end
  geometry <- sf::st_sfc(lapply(which(keep), function(i) sf::st_linestring(matrix(c(
    data$start_longitude[[i]], data$start_latitude[[i]],
    data$end_longitude[[i]], data$end_latitude[[i]]
  ), ncol = 2, byrow = TRUE))), crs = 4326)
  sf::st_transform(sf::st_sf(
    visit_id = as.character(data[[id]][keep]), source_row = which(keep),
    geometry = geometry
  ), crs)
}

compile_grid_traversal <- function() {
  if (!requireNamespace("Rcpp", quietly = TRUE)) stop("Package Rcpp is required.")
  if (exists("line_grid_links_cpp", mode = "function", inherits = TRUE)) return(invisible(TRUE))
  Rcpp::cppFunction(code = '
    Rcpp::DataFrame line_grid_links_cpp(
      Rcpp::NumericVector x0, Rcpp::NumericVector y0,
      Rcpp::NumericVector x1, Rcpp::NumericVector y1, double cell_size) {
      std::vector<int> row_id, x_bin, y_bin;
      std::vector<double> fraction, geometric_length;
      for (int i = 0; i < x0.size(); ++i) {
        if (!R_finite(x0[i]) || !R_finite(y0[i])) continue;
        if (!R_finite(x1[i]) || !R_finite(y1[i])) {
          row_id.push_back(i + 1);
          x_bin.push_back((int)std::floor(x0[i] / cell_size));
          y_bin.push_back((int)std::floor(y0[i] / cell_size));
          fraction.push_back(1.0);
          geometric_length.push_back(0.0);
          continue;
        }
        double dx = x1[i] - x0[i], dy = y1[i] - y0[i];
        double total_length = std::sqrt(dx * dx + dy * dy);
        if (total_length < 1e-8) {
          row_id.push_back(i + 1);
          x_bin.push_back((int)std::floor(x0[i] / cell_size));
          y_bin.push_back((int)std::floor(y0[i] / cell_size));
          fraction.push_back(1.0);
          geometric_length.push_back(0.0);
          continue;
        }
        std::vector<double> cuts; cuts.push_back(0.0); cuts.push_back(1.0);
        if (std::abs(dx) > 1e-12) {
          int lo = (int)std::floor(std::min(x0[i], x1[i]) / cell_size) + 1;
          int hi = (int)std::floor(std::max(x0[i], x1[i]) / cell_size);
          for (int k = lo; k <= hi; ++k) {
            double t = (k * cell_size - x0[i]) / dx;
            if (t > 0.0 && t < 1.0) cuts.push_back(t);
          }
        }
        if (std::abs(dy) > 1e-12) {
          int lo = (int)std::floor(std::min(y0[i], y1[i]) / cell_size) + 1;
          int hi = (int)std::floor(std::max(y0[i], y1[i]) / cell_size);
          for (int k = lo; k <= hi; ++k) {
            double t = (k * cell_size - y0[i]) / dy;
            if (t > 0.0 && t < 1.0) cuts.push_back(t);
          }
        }
        std::sort(cuts.begin(), cuts.end());
        cuts.erase(std::unique(cuts.begin(), cuts.end(), [](double a, double b) {
          return std::abs(a - b) < 1e-10;
        }), cuts.end());
        for (unsigned int j = 0; j + 1 < cuts.size(); ++j) {
          double dt = cuts[j + 1] - cuts[j];
          if (dt <= 1e-12) continue;
          double tm = (cuts[j + 1] + cuts[j]) / 2.0;
          row_id.push_back(i + 1);
          x_bin.push_back((int)std::floor((x0[i] + tm * dx) / cell_size));
          y_bin.push_back((int)std::floor((y0[i] + tm * dy) / cell_size));
          fraction.push_back(dt);
          geometric_length.push_back(dt * total_length);
        }
      }
      return Rcpp::DataFrame::create(
        Rcpp::_["source_row"] = row_id, Rcpp::_["x_bin"] = x_bin,
        Rcpp::_["y_bin"] = y_bin, Rcpp::_["fraction"] = fraction,
        Rcpp::_["geometric_length_m"] = geometric_length
      );
    }
  ', depends = "Rcpp")
  invisible(TRUE)
}

line_grid_links_r <- function(x0, y0, x1, y1, cell_size, chunk_size = 20000L) {
  one_line <- function(i) {
    if (!is.finite(x0[[i]]) || !is.finite(y0[[i]])) return(NULL)
    if (!is.finite(x1[[i]]) || !is.finite(y1[[i]])) {
      return(data.frame(
        source_row = i, x_bin = floor(x0[[i]] / cell_size),
        y_bin = floor(y0[[i]] / cell_size), fraction = 1,
        geometric_length_m = 0
      ))
    }
    dx <- x1[[i]] - x0[[i]]
    dy <- y1[[i]] - y0[[i]]
    total_length <- sqrt(dx^2 + dy^2)
    if (!is.finite(total_length) || total_length < 1e-8) {
      return(data.frame(
        source_row = i, x_bin = floor(x0[[i]] / cell_size),
        y_bin = floor(y0[[i]] / cell_size), fraction = 1,
        geometric_length_m = 0
      ))
    }
    cuts <- c(0, 1)
    if (abs(dx) > 1e-12) {
      lower <- (floor(min(x0[[i]], x1[[i]]) / cell_size) + 1) * cell_size
      upper <- floor(max(x0[[i]], x1[[i]]) / cell_size) * cell_size
      boundaries <- if (lower <= upper) seq(lower, upper, by = cell_size) else numeric()
      if (length(boundaries)) cuts <- c(cuts, (boundaries - x0[[i]]) / dx)
    }
    if (abs(dy) > 1e-12) {
      lower <- (floor(min(y0[[i]], y1[[i]]) / cell_size) + 1) * cell_size
      upper <- floor(max(y0[[i]], y1[[i]]) / cell_size) * cell_size
      boundaries <- if (lower <= upper) seq(lower, upper, by = cell_size) else numeric()
      if (length(boundaries)) cuts <- c(cuts, (boundaries - y0[[i]]) / dy)
    }
    cuts <- sort(unique(cuts[is.finite(cuts) & cuts >= 0 & cuts <= 1]))
    dt <- diff(cuts)
    keep <- dt > 1e-12
    tm <- (head(cuts, -1) + tail(cuts, -1)) / 2
    data.frame(
      source_row = rep.int(i, sum(keep)),
      x_bin = floor((x0[[i]] + tm[keep] * dx) / cell_size),
      y_bin = floor((y0[[i]] + tm[keep] * dy) / cell_size),
      fraction = dt[keep], geometric_length_m = dt[keep] * total_length
    )
  }
  starts <- seq.int(1L, length(x0), by = chunk_size)
  pieces <- lapply(seq_along(starts), function(k) {
    idx <- starts[[k]]:min(length(x0), starts[[k]] + chunk_size - 1L)
    message("  grid traversal chunk ", k, " of ", length(starts))
    rows <- lapply(idx, one_line)
    if (requireNamespace("data.table", quietly = TRUE)) {
      as.data.frame(data.table::rbindlist(rows))
    } else {
      dplyr::bind_rows(rows)
    }
  })
  dplyr::bind_rows(pieces)
}

project_endpoints <- function(data, target_crs) {
  start <- sf::sf_project(
    "EPSG:4326", target_crs,
    as.matrix(data[, c("start_longitude", "start_latitude")])
  )
  end_lon <- dplyr::if_else(
    valid_lonlat(data$end_longitude, data$end_latitude),
    data$end_longitude, data$start_longitude
  )
  end_lat <- dplyr::if_else(
    valid_lonlat(data$end_longitude, data$end_latitude),
    data$end_latitude, data$start_latitude
  )
  end <- sf::sf_project("EPSG:4326", target_crs, cbind(end_lon, end_lat))
  dplyr::mutate(
    data, start_x = start[, 1], start_y = start[, 2],
    end_x = end[, 1], end_y = end[, 2],
    geometry_fallback = !valid_lonlat(.data$end_longitude, .data$end_latitude)
  )
}

line_cell_links <- function(data, target_crs, cell_size_m, reported_length = NULL) {
  projected <- project_endpoints(data, target_crs)
  projected$geometry_quality_reason <- dplyr::if_else(
    projected$geometry_fallback, "missing_endpoint", "valid_line"
  )
  if (!is.null(reported_length)) {
    reported <- as.numeric(projected[[reported_length]])
    geometric <- sqrt(
      (projected$end_x - projected$start_x)^2 +
        (projected$end_y - projected$start_y)^2
    )
    mismatch <- is.finite(reported) & reported > 0 &
      is.finite(geometric) & geometric / reported > 3
    projected$end_x[mismatch] <- projected$start_x[mismatch]
    projected$end_y[mismatch] <- projected$start_y[mismatch]
    projected$geometry_fallback[mismatch] <- TRUE
    projected$geometry_quality_reason[mismatch] <- "endpoint_length_mismatch"
  }
  links <- line_grid_links_r(
    projected$start_x, projected$start_y, projected$end_x, projected$end_y,
    as.numeric(cell_size_m)
  )
  links <- dplyr::left_join(
    tibble::as_tibble(links),
    dplyr::mutate(projected, source_row = dplyr::row_number()),
    by = "source_row"
  )
  if (is.null(reported_length)) {
    links$weighted_length_m <- links$geometric_length_m
  } else {
    supplied <- as.numeric(links[[reported_length]])
    links$weighted_length_m <- dplyr::if_else(
      is.finite(supplied) & supplied > 0,
      links$fraction * supplied, links$geometric_length_m
    )
  }
  links |>
    dplyr::mutate(
      cell_id = paste0(cell_size_m, "m||", .data$x_bin, "||", .data$y_bin),
      projected_x = (.data$x_bin + 0.5) * cell_size_m,
      projected_y = (.data$y_bin + 0.5) * cell_size_m,
      cell_support_m = as.integer(cell_size_m)
    )
}

point_cell_links <- function(data, target_crs, cell_size_m) {
  xy <- sf::sf_project(
    "EPSG:4326", target_crs,
    as.matrix(data[, c("longitude", "latitude")])
  )
  dplyr::mutate(
    data, x_bin = floor(xy[, 1] / cell_size_m),
    y_bin = floor(xy[, 2] / cell_size_m),
    cell_id = paste0(cell_size_m, "m||", .data$x_bin, "||", .data$y_bin),
    projected_x = (.data$x_bin + 0.5) * cell_size_m,
    projected_y = (.data$y_bin + 0.5) * cell_size_m,
    fraction = 1, weighted_length_m = NA_real_,
    geometry_fallback = TRUE, cell_support_m = as.integer(cell_size_m)
  )
}

cluster_alignment_windows <- function(visits, maximum_gap_days) {
  visits |>
    dplyr::arrange(.data$cell_id, .data$survey_date, .data$method, .data$visit_id) |>
    dplyr::group_by(.data$cell_id) |>
    dplyr::mutate(
      previous_date = dplyr::lag(.data$survey_date),
      date_gap = as.integer(.data$survey_date - .data$previous_date),
      new_cluster = dplyr::row_number() == 1L |
        dplyr::coalesce(.data$date_gap > maximum_gap_days, TRUE),
      cluster_number = cumsum(.data$new_cluster),
      alignment_id = paste0(
        .data$cell_id, "||w", maximum_gap_days, "||", .data$cluster_number
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-"previous_date", -"new_cluster", -"cluster_number")
}

annotate_support_cull_history <- function(visits) {
  first_cull <- visits |>
    dplyr::filter(.data$method == "cull") |>
    dplyr::group_by(.data$cell_id) |>
    dplyr::summarise(local_first_cull_date = min(.data$survey_date), .groups = "drop")
  visits |>
    dplyr::left_join(first_cull, by = "cell_id") |>
    dplyr::mutate(
      days_from_local_first_cull = as.integer(
        .data$survey_date - .data$local_first_cull_date
      ),
      cull_timing = dplyr::case_when(
        is.na(.data$local_first_cull_date) ~ "no_local_cull_record",
        .data$days_from_local_first_cull < 0 ~ "strictly_pre_cull",
        .data$days_from_local_first_cull == 0 ~ "same_day",
        TRUE ~ "post_cull"
      ),
      ecological_eligible = .data$cull_timing %in% c(
        "no_local_cull_record", "strictly_pre_cull"
      ) | (.data$method == "cull" & .data$cull_timing == "same_day")
    )
}
