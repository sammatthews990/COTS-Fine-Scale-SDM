# Reef-context predictors at the 0.5-1 km scale, on the 120 m support grid.
#
# Two sources:
#   - reef geometry from the GBR reef-feature polygons (no new data);
#   - eReefs GBR1 v2.0 hydrodynamics, whole-period mean ("all-one"), from the
#     AIMS THREDDS server (regular 0.005 deg grid).

ereefs_all_one_url <- paste0(
  "https://thredds.ereefs.aims.gov.au/thredds/dodsC/ereefs/gbr1_2.0/all-one/",
  "EREEFS_AIMS-CSIRO_gbr1_2.0_hydro_all-one.nc"
)
ereefs_grid <- list(
  nlat = 4237L, nlon = 2671L,
  lat0 = -28.589904, lon0 = 142.020371, step = 0.005
)
# zc index (0-based) -> depth: 12 = -9 m, 13 = -5.35 m, 14 = -2.35 m, 15 = -0.5 m
ereefs_layers <- c(`12` = -9, `13` = -5.35, `14` = -2.35, `15` = -0.5)

# Download one depth slice of one variable as an OPeNDAP binary response and
# return it as a lon/lat SpatRaster. Cached in `cache_dir`.
read_ereefs_slice <- function(variable, k, cache_dir) {
  path <- file.path(cache_dir, sprintf("%s_k%d.dods", variable, k))
  if (!file.exists(path)) {
    query <- sprintf("%s.dods?%s[0][%d][0:1:%d][0:1:%d]", ereefs_all_one_url, variable, k,
                     ereefs_grid$nlat - 1L, ereefs_grid$nlon - 1L)
    utils::download.file(query, path, mode = "wb", quiet = TRUE)
  }
  bytes <- readBin(path, "raw", file.info(path)$size)
  marker <- charToRaw("\nData:\n")
  start <- NA_integer_
  for (i in seq_len(min(length(bytes), 4000L) - length(marker))) {
    if (all(bytes[i:(i + length(marker) - 1L)] == marker)) { start <- i + length(marker); break }
  }
  if (is.na(start)) stop("Unexpected OPeNDAP response for ", variable)
  n <- ereefs_grid$nlat * ereefs_grid$nlon
  # XDR array: two 4-byte length words, then big-endian float32 values.
  con <- rawConnection(bytes[(start + 8L):length(bytes)])
  on.exit(close(con))
  values <- readBin(con, "numeric", n = n, size = 4, endian = "big")
  values[!is.finite(values) | abs(values) > 1e10] <- NA
  m <- matrix(values, nrow = ereefs_grid$nlat, ncol = ereefs_grid$nlon, byrow = TRUE)
  m <- m[nrow(m):1, , drop = FALSE] # latitude ascends in the file; rasters run north -> south
  half <- ereefs_grid$step / 2
  terra::rast(
    m, crs = "EPSG:4326",
    extent = terra::ext(
      ereefs_grid$lon0 - half, ereefs_grid$lon0 + (ereefs_grid$nlon - 1) * ereefs_grid$step + half,
      ereefs_grid$lat0 - half, ereefs_grid$lat0 + (ereefs_grid$nlat - 1) * ereefs_grid$step + half
    )
  )
}

# Shallowest-available current fields: the -9 m layer, filled upward where the
# model column is shallower than 9 m.
ereefs_current_layers <- function(cache_dir) {
  stack_for <- function(variable) {
    layers <- lapply(names(ereefs_layers), function(k) read_ereefs_slice(variable, as.integer(k), cache_dir))
    Reduce(terra::cover, layers)
  }
  speed <- stack_for("mean_cur")
  u <- stack_for("u")
  v <- stack_for("v")
  residual <- sqrt(u^2 + v^2)
  steadiness <- terra::clamp(residual / speed, 0, 1)
  out <- c(speed, residual, steadiness)
  names(out) <- c("current_speed", "residual_current", "flow_steadiness")
  out
}

# Fill NA cells from the nearest valued neighbours, up to `iterations` cells.
fill_nearby <- function(x, iterations = 5L) {
  for (i in seq_len(iterations)) {
    x <- terra::focal(x, w = 3, fun = "mean", na.policy = "only", na.rm = TRUE)
  }
  x
}

# Inter-reef gap for each reef-edge cell: distance to the nearest edge cell of
# a different reef. Searches widen k up to `max_k`; edge cells whose own reef
# fills every neighbour get the k-th neighbour distance as a lower bound.
edge_gap_distance <- function(edge_xy, edge_id, max_k = 1024L) {
  out <- rep(NA_real_, nrow(edge_xy))
  todo <- seq_len(nrow(edge_xy))
  k <- 32L
  repeat {
    nn <- FNN::get.knnx(edge_xy, edge_xy[todo, , drop = FALSE], k = min(k, nrow(edge_xy)))
    other <- matrix(edge_id[nn$nn.index], nrow = length(todo)) != edge_id[todo]
    first <- max.col(other, ties.method = "first")
    found <- other[cbind(seq_along(todo), first)]
    out[todo[found]] <- nn$nn.dist[cbind(which(found), first[found])]
    if (k >= max_k || all(found)) {
      out[todo[!found]] <- nn$nn.dist[!found, ncol(nn$nn.dist)]
      break
    }
    todo <- todo[!found]
    k <- k * 4L
  }
  out
}
