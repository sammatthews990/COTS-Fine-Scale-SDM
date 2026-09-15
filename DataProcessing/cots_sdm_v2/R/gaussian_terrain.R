# Gaussian terrain functions are intentionally tile-based. The GBR-wide 30 m
# DEM has more than two billion cells, so full in-memory focal processing is
# neither necessary nor reproducible for a point-sampled model.
gaussian_scale_spec <- function() {
  tibble::tibble(
    scale_m = c(30, 60, 120, 240, 360),
    sigma_m = c(0, 30, 60, 120, 180)
  )
}

gaussian_smooth_dem <- function(dem, sigma_m) {
  if (sigma_m <= 0) return(dem)
  weights <- terra::focalMat(dem, sigma_m, type = "Gauss")
  terra::focal(dem, w = weights, fun = "sum", na.rm = TRUE, fillvalue = NA)
}

terrain_derivatives <- function(dem, scale_m) {
  resolution <- terra::res(dem)
  dx <- resolution[[1]]
  dy <- resolution[[2]]
  base_gradient <- matrix(c(-1, 0, 1, -2, 0, 2, -1, 0, 1), 3, 3)
  kx <- base_gradient / (8 * dx)
  ky <- t(base_gradient) / (8 * dy)
  kxx <- matrix(c(0, 0, 0, 1, -2, 1, 0, 0, 0), 3, 3) / dx^2
  kyy <- t(matrix(c(0, 0, 0, 1, -2, 1, 0, 0, 0), 3, 3)) / dy^2
  kxy <- matrix(c(1, 0, -1, 0, 0, 0, -1, 0, 1), 3, 3) / (4 * dx * dy)
  p <- terra::focal(dem, w = kx, fun = "sum", na.rm = TRUE)
  q <- terra::focal(dem, w = ky, fun = "sum", na.rm = TRUE)
  r <- terra::focal(dem, w = kxx, fun = "sum", na.rm = TRUE)
  t <- terra::focal(dem, w = kyy, fun = "sum", na.rm = TRUE)
  s <- terra::focal(dem, w = kxy, fun = "sum", na.rm = TRUE)
  gradient <- sqrt(p^2 + q^2)
  slope <- atan(gradient) * 180 / pi
  gradient_safe <- terra::ifel(gradient < 1e-8, 1e-8, gradient)
  eastness <- p / gradient_safe
  northness <- q / gradient_safe
  denominator_plan_raw <- (p^2 + q^2)^(3 / 2)
  denominator_plan <- terra::ifel(
    denominator_plan_raw < 1e-8, 1e-8, denominator_plan_raw
  )
  denominator_profile_raw <-
    (p^2 + q^2) * (1 + p^2 + q^2)^(3 / 2)
  denominator_profile <- terra::ifel(
    denominator_profile_raw < 1e-8, 1e-8, denominator_profile_raw
  )
  hcu <- (t * p^2 - 2 * s * p * q + r * q^2) / denominator_plan
  vcu <- -(r * p^2 + 2 * s * p * q + t * q^2) / denominator_profile
  outer_radius <- max(60, scale_m * 1.5)
  outer_weights <- terra::focalMat(dem, outer_radius, type = "circle")
  bpi_mean <- terra::focal(
    dem, w = outer_weights, fun = "sum", na.rm = TRUE
  )
  bpi <- dem - bpi_mean
  second_moment <- terra::focal(
    dem^2, w = outer_weights, fun = "sum", na.rm = TRUE
  )
  variance <- second_moment - bpi_mean^2
  local_sd <- sqrt(terra::ifel(variance < 0, 0, variance))
  result <- c(dem, slope, eastness, northness, hcu, vcu, bpi, local_sd)
  names(result) <- paste0(
    c("dem", "slope", "eastness", "northness", "hcu", "vcu", "bpi", "sd"),
    "_g", scale_m, "m"
  )
  result
}

extract_gaussian_terrain_tiled <- function(
    points, template, tile_m = 20000, buffer_m = 900,
    cache_dir = NULL, id_field = "cell_id") {
  if (!all(c("projected_x", "projected_y") %in% names(points))) {
    stop("Projected coordinates are required for Gaussian extraction.")
  }
  dem <- template[["DEM"]]
  points <- points |>
    dplyr::mutate(
      tile_x = floor(.data$projected_x / tile_m),
      tile_y = floor(.data$projected_y / tile_m),
      terrain_row_id = dplyr::row_number()
    )
  groups <- split(points, interaction(points$tile_x, points$tile_y, drop = TRUE))
  spec <- gaussian_scale_spec()
  if (!is.null(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }
  extracted <- lapply(seq_along(groups), function(group_index) {
    group <- groups[[group_index]]
    if (group_index %% 25L == 1L || group_index == length(groups)) {
      message("Gaussian terrain tile ", group_index, " of ", length(groups))
    }
    cache_file <- if (is.null(cache_dir)) NULL else file.path(
      cache_dir,
      paste0("tile_", group$tile_x[[1]], "_", group$tile_y[[1]], ".rds")
    )
    if (!is.null(cache_file) && file.exists(cache_file)) {
      cached <- readRDS(cache_file)
      if (
        id_field %in% names(group) && id_field %in% names(cached) &&
        setequal(as.character(group[[id_field]]), as.character(cached[[id_field]]))
      ) {
        return(dplyr::left_join(
          dplyr::select(group, "terrain_row_id", dplyr::all_of(id_field)),
          cached, by = id_field
        ))
      }
    }
    extent <- terra::ext(
      min(group$projected_x) - buffer_m,
      max(group$projected_x) + buffer_m,
      min(group$projected_y) - buffer_m,
      max(group$projected_y) + buffer_m
    )
    tile_dem <- terra::crop(dem, extent, snap = "out")
    xy <- as.matrix(group[, c("projected_x", "projected_y")])
    values <- lapply(seq_len(nrow(spec)), function(i) {
      smooth <- gaussian_smooth_dem(tile_dem, spec$sigma_m[[i]])
      terrain <- terrain_derivatives(smooth, spec$scale_m[[i]])
      terra::extract(terrain, xy)
    })
    value <- dplyr::bind_cols(
      dplyr::select(group, "terrain_row_id", dplyr::any_of(id_field)),
      dplyr::bind_cols(values)
    )
    if (!is.null(cache_file) && id_field %in% names(value)) {
      saveRDS(dplyr::select(value, -"terrain_row_id"), cache_file)
    }
    value
  })
  dplyr::bind_rows(extracted) |>
    dplyr::arrange(.data$terrain_row_id)
}

point_gaussian_terrain <- function(
    dem_matrix, centre_row, centre_col, resolution, scale_m) {
  sigma_pixels <- if (scale_m == 30) 0.75 else (scale_m / 2) / resolution
  radius <- max(2L, as.integer(ceiling(3 * sigma_pixels)))
  rows <- seq.int(
    max(1L, centre_row - radius),
    min(nrow(dem_matrix), centre_row + radius)
  )
  columns <- seq.int(
    max(1L, centre_col - radius),
    min(ncol(dem_matrix), centre_col + radius)
  )
  value <- as.vector(dem_matrix[rows, columns, drop = FALSE])
  grid <- expand.grid(row = rows, column = columns)
  x <- (grid$column - centre_col) * resolution
  y <- (centre_row - grid$row) * resolution
  sigma_m <- sigma_pixels * resolution
  weight <- exp(-(x^2 + y^2) / (2 * sigma_m^2))
  valid <- is.finite(value) & is.finite(weight)
  if (sum(valid) < 8L) return(rep(NA_real_, 8L))
  value <- value[valid]
  x <- x[valid]
  y <- y[valid]
  weight <- weight[valid]
  weight <- weight / sum(weight)
  design <- cbind(1, x, y, x^2, x * y, y^2)
  penalty <- diag(c(1e-10, rep(1e-8, 5L)))
  coefficient <- tryCatch(
    solve(
      crossprod(design, design * weight) + penalty,
      crossprod(design, weight * value)
    ),
    error = function(e) rep(NA_real_, 6L)
  )
  if (any(!is.finite(coefficient))) return(rep(NA_real_, 8L))
  smoothed_dem <- sum(weight * value)
  p <- coefficient[[2]]
  q <- coefficient[[3]]
  r <- 2 * coefficient[[4]]
  s <- coefficient[[5]]
  t <- 2 * coefficient[[6]]
  gradient <- sqrt(p^2 + q^2)
  slope <- atan(gradient) * 180 / pi
  eastness <- if (gradient < 1e-8) 0 else p / gradient
  northness <- if (gradient < 1e-8) 0 else q / gradient
  plan_denominator <- max((p^2 + q^2)^(3 / 2), 1e-8)
  profile_denominator <- max(
    (p^2 + q^2) * (1 + p^2 + q^2)^(3 / 2), 1e-8
  )
  hcu <- (t * p^2 - 2 * s * p * q + r * q^2) / plan_denominator
  vcu <- -(r * p^2 + 2 * s * p * q + t * q^2) / profile_denominator
  centre_value <- dem_matrix[centre_row, centre_col]
  bpi <- if (is.finite(centre_value)) centre_value - smoothed_dem else NA_real_
  local_sd <- sqrt(sum(weight * (value - smoothed_dem)^2))
  c(smoothed_dem, slope, eastness, northness, hcu, vcu, bpi, local_sd)
}

extract_gaussian_terrain_points <- function(
    points, template, tile_m = 5000, buffer_m = 900,
    cache_dir = NULL, id_field = "cell_id") {
  if (!all(c("projected_x", "projected_y") %in% names(points))) {
    stop("Projected coordinates are required for Gaussian extraction.")
  }
  dem <- template[["DEM"]]
  points <- points |>
    dplyr::mutate(
      tile_x = floor(.data$projected_x / tile_m),
      tile_y = floor(.data$projected_y / tile_m),
      terrain_row_id = dplyr::row_number()
    )
  groups <- split(points, interaction(points$tile_x, points$tile_y, drop = TRUE))
  if (!is.null(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }
  scale_m <- gaussian_scale_spec()$scale_m
  variable <- c("dem", "slope", "eastness", "northness", "hcu", "vcu", "bpi", "sd")
  output_names <- unlist(lapply(scale_m, function(scale) {
    paste0(variable, "_g", scale, "m")
  }))
  pieces <- lapply(seq_along(groups), function(group_index) {
    group <- groups[[group_index]]
    if (group_index %% 50L == 1L || group_index == length(groups)) {
      message("Point Gaussian terrain tile ", group_index, " of ", length(groups))
    }
    cache_file <- if (is.null(cache_dir)) NULL else file.path(
      cache_dir,
      paste0("tile_", group$tile_x[[1]], "_", group$tile_y[[1]], ".rds")
    )
    if (!is.null(cache_file) && file.exists(cache_file)) {
      cached <- readRDS(cache_file)
      if (
        id_field %in% names(cached) &&
        setequal(as.character(group[[id_field]]), as.character(cached[[id_field]]))
      ) {
        return(dplyr::left_join(
          dplyr::select(group, "terrain_row_id", dplyr::all_of(id_field)),
          cached, by = id_field
        ))
      }
    }
    extent <- terra::ext(
      min(group$projected_x) - buffer_m,
      max(group$projected_x) + buffer_m,
      min(group$projected_y) - buffer_m,
      max(group$projected_y) + buffer_m
    )
    tile_dem <- terra::crop(dem, extent, snap = "out")
    dem_matrix <- as.matrix(tile_dem, wide = TRUE)
    cells <- terra::cellFromXY(
      tile_dem,
      as.matrix(group[, c("projected_x", "projected_y")])
    )
    row_column <- terra::rowColFromCell(tile_dem, cells)
    terrain_matrix <- t(vapply(seq_len(nrow(group)), function(i) {
      unlist(lapply(scale_m, function(scale) {
        point_gaussian_terrain(
          dem_matrix, row_column[i, 1], row_column[i, 2],
          terra::res(tile_dem)[[1]], scale
        )
      }))
    }, numeric(length(output_names))))
    colnames(terrain_matrix) <- output_names
    value <- dplyr::bind_cols(
      dplyr::select(group, "terrain_row_id", dplyr::all_of(id_field)),
      tibble::as_tibble(terrain_matrix)
    )
    if (!is.null(cache_file)) {
      saveRDS(dplyr::select(value, -"terrain_row_id"), cache_file)
    }
    value
  })
  dplyr::bind_rows(pieces) |>
    dplyr::arrange(.data$terrain_row_id)
}
