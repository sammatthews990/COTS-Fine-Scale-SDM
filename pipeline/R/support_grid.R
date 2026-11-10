# Full-domain 120 m support grid for the Hmsc layer.
#
# Training cells (v3 05_extract_support_predictors.R) are 120 m bins on a grid
# anchored at 0 m (MGA55). Their predictors are:
#   - base layers: mean of 16 point samples at +/-15 and +/-45 m from the centre,
#     each falling in a distinct 30 m pixel, i.e. a 4x4 block mean;
#   - geomorphic class at the cell centre;
#   - multiscale Gaussian terrain fitted at the cell centre's 30 m pixel.
# The functions below reproduce those definitions for every cell in the domain.

hmsc_base_layers <- c(
  "EAST", "HCU", "NORTH", "SLO", "SVF", "VCU", "VRM", "DEM",
  "AhyaD", "Aspat", "Aten", "turbid", "waves_Hs", "waves_Tp"
)

# Crop window and shift that make aggregate(fact = 4) equal the 16-point
# quadrature mean of 120 m bins anchored at 0.
quadrature_alignment <- function(r, support_m = 120) {
  e <- terra::ext(r)
  res <- terra::res(r)
  stopifnot(all(abs(res - 30) < 1e-9), support_m == 120)
  # First and last bins fully inside the raster.
  bx <- c(ceiling(e$xmin / support_m), floor(e$xmax / support_m) - 1)
  by <- c(ceiling(e$ymin / support_m), floor(e$ymax / support_m) - 1)
  # Column containing the first sample (+15 m) of the first bin; row containing
  # the top sample (+105 m) of the last bin (rows run from the top).
  col0 <- floor((support_m * bx[1] + 15 - e$xmin) / 30)
  row0 <- floor((e$ymax - (support_m * by[2] + 105)) / 30)
  n_bx <- diff(bx) + 1
  n_by <- diff(by) + 1
  crop_ext <- terra::ext(
    e$xmin + 30 * col0, e$xmin + 30 * col0 + support_m * n_bx,
    e$ymax - 30 * row0 - support_m * n_by, e$ymax - 30 * row0
  )
  list(
    crop_ext = crop_ext,
    dx = support_m * bx[1] - crop_ext$xmin,
    dy = support_m * by[1] - crop_ext$ymin,
    bins_x = bx, bins_y = by
  )
}

# Point-wise local quadratic Gaussian terrain (copied unchanged from
# cots_sdm_v2/R/gaussian_terrain.R so this pipeline is self-contained).
gaussian_scales_m <- c(30, 60, 120, 240, 360)
terrain_variables <- c("dem", "slope", "eastness", "northness", "hcu", "vcu", "bpi", "sd")
terrain_column_names <- unlist(lapply(gaussian_scales_m, function(s) paste0(terrain_variables, "_g", s, "m")))

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

# Vectorised point_gaussian_terrain for points whose window is complete (inside
# the matrix and all values finite). For such windows the weighted
# least-squares fit is one fixed linear operator per scale, so all points are
# solved with a single matrix product. Returns NULL rows for points that need
# the exact per-point routine (window clipped or containing NA).
fast_gaussian_terrain <- function(dem_matrix, centre_rows, centre_cols, resolution, scale_m) {
  sigma_pixels <- if (scale_m == 30) 0.75 else (scale_m / 2) / resolution
  radius <- max(2L, as.integer(ceiling(3 * sigma_pixels)))
  offsets <- expand.grid(dr = -radius:radius, dc = -radius:radius)
  x <- offsets$dc * resolution
  y <- -offsets$dr * resolution
  sigma_m <- sigma_pixels * resolution
  weight <- exp(-(x^2 + y^2) / (2 * sigma_m^2))
  weight <- weight / sum(weight)
  design <- cbind(1, x, y, x^2, x * y, y^2)
  penalty <- diag(c(1e-10, rep(1e-8, 5L)))
  operator <- solve(crossprod(design, design * weight) + penalty, t(design * weight))

  inside <- centre_rows - radius >= 1L & centre_rows + radius <= nrow(dem_matrix) &
    centre_cols - radius >= 1L & centre_cols + radius <= ncol(dem_matrix)
  out <- matrix(NA_real_, length(centre_rows), 8L)
  fast <- rep(FALSE, length(centre_rows))
  if (!any(inside)) return(list(values = out, fast = fast))
  idx <- which(inside)
  cell_index <- outer(offsets$dr, centre_rows[idx], `+`) +
    (outer(offsets$dc, centre_cols[idx], `+`) - 1L) * nrow(dem_matrix)
  # as.vector: a 2-column index matrix would be read as (row, col) pairs.
  v <- matrix(dem_matrix[as.vector(cell_index)], nrow = nrow(offsets))
  complete <- colSums(!is.finite(v)) == 0L
  idx <- idx[complete]
  v <- v[, complete, drop = FALSE]
  if (!length(idx)) return(list(values = out, fast = fast))
  coefficient <- operator %*% v
  smoothed <- colSums(weight * v)
  p <- coefficient[2, ]
  q <- coefficient[3, ]
  r <- 2 * coefficient[4, ]
  s <- coefficient[5, ]
  t <- 2 * coefficient[6, ]
  gradient <- sqrt(p^2 + q^2)
  flat <- gradient < 1e-8
  centre <- v[(nrow(offsets) + 1L) / 2L, ]
  out[idx, ] <- cbind(
    smoothed,
    atan(gradient) * 180 / pi,
    ifelse(flat, 0, p / gradient),
    ifelse(flat, 0, q / gradient),
    (t * p^2 - 2 * s * p * q + r * q^2) / pmax((p^2 + q^2)^(3 / 2), 1e-8),
    -(r * p^2 + 2 * s * p * q + t * q^2) / pmax((p^2 + q^2) * (1 + p^2 + q^2)^(3 / 2), 1e-8),
    centre - smoothed,
    sqrt(colSums(weight * (v - rep(smoothed, each = nrow(v)))^2))
  )
  fast[idx] <- TRUE
  list(values = out, fast = fast)
}

# Terrain for all points in one tile. `dem_path` is opened inside the call so
# the function can run on a parallel worker.
tile_gaussian_terrain <- function(points, dem_path, dem_layer = "DEM", buffer_m = 900) {
  dem <- terra::rast(dem_path)[[dem_layer]]
  extent <- terra::ext(
    min(points$x) - buffer_m, max(points$x) + buffer_m,
    min(points$y) - buffer_m, max(points$y) + buffer_m
  )
  if (is.null(terra::intersect(extent, terra::ext(dem)))) {
    values <- matrix(NA_real_, nrow(points), length(terrain_column_names),
                     dimnames = list(NULL, terrain_column_names))
    return(cbind(data.frame(cell = points$cell), values))
  }
  tile_dem <- terra::crop(dem, extent, snap = "out")
  dem_matrix <- terra::as.matrix(tile_dem, wide = TRUE)
  rc <- terra::rowColFromCell(tile_dem, terra::cellFromXY(tile_dem, cbind(points$x, points$y)))
  resolution <- terra::res(tile_dem)[[1]]
  values <- do.call(cbind, lapply(gaussian_scales_m, function(scale) {
    result <- fast_gaussian_terrain(dem_matrix, rc[, 1], rc[, 2], resolution, scale)
    for (i in which(!result$fast)) {
      result$values[i, ] <- point_gaussian_terrain(dem_matrix, rc[i, 1], rc[i, 2], resolution, scale)
    }
    result$values
  }))
  colnames(values) <- terrain_column_names
  cbind(data.frame(cell = points$cell), values)
}

# Geomorphic class at projected cell centres, tiled to keep reads local.
extract_geomorphic_class <- function(x, y, crs, geomorphic_path, tile_deg = 0.5) {
  geomorphic <- terra::rast(geomorphic_path)
  lonlat <- terra::crds(terra::project(terra::vect(cbind(x, y), crs = crs), terra::crs(geomorphic)))
  labels <- rep(NA_character_, length(x))
  tile <- paste(floor(lonlat[, 1] / tile_deg), floor(lonlat[, 2] / tile_deg))
  groups <- split(seq_along(x), tile)
  ge <- terra::ext(geomorphic)
  for (i in seq_along(groups)) {
    rows <- groups[[i]]
    e <- terra::ext(
      min(lonlat[rows, 1]) - 0.002, max(lonlat[rows, 1]) + 0.002,
      min(lonlat[rows, 2]) - 0.002, max(lonlat[rows, 2]) + 0.002
    )
    if (e$xmin > ge$xmax || e$xmax < ge$xmin || e$ymin > ge$ymax || e$ymax < ge$ymin) next
    crop <- terra::crop(geomorphic, e, snap = "out")
    extracted <- terra::extract(crop, lonlat[rows, , drop = FALSE])
    labels[rows] <- as.character(extracted[[ncol(extracted)]])
  }
  labels
}
