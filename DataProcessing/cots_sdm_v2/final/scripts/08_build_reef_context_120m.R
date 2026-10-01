# Reef-context predictors at the 0.5-1 km scale on the 120 m support grid.
#
# Geometry (from GBR reef polygons):
#   dist_to_reef_edge_m     signed: + inside a reef polygon, - outside
#   windward_position       0-1 position along the SE axis within the reef (1 = SE / windward edge)
#   seaward_position        0-1 position along the NE axis within the reef (1 = NE / seaward edge)
#   upwind_open_water_km    distance SE (into the trades) to the next reef, capped at 10 km
#   nearest_other_reef_km   distance to the nearest edge plus the inter-reef gap there (channels / passes)
#   log_reef_area_km2       size of the assigned reef
# Hydrodynamics (eReefs GBR1 v2.0 whole-period mean, -9 m, filled upward):
#   current_speed, residual_current, flow_steadiness (residual / mean speed;
#   low = oscillating tidal flow, a retention proxy)
#
# Cells outside a polygon take the nearest reef within 600 m.

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "..", "R", "common.R"))
cfg <- initialise_final()
out_path <- file.path(cfg$dirs$cache, "reef_context_120m.tif")
template <- terra::rast(file.path(cfg$dirs$cache, "base_support_120m.tif"))[[1]]
grid <- readRDS(file.path(cfg$dirs$cache, "support_grid_120m.rds"))[, c("cell", "x", "y")]
n <- nrow(grid)
log_message("Domain cells: ", format(n, big.mark = ","))

# --- Reef assignment --------------------------------------------------------
reefs <- sf::st_read(checked_input(cfg, "reef_features"), quiet = TRUE) |>
  dplyr::filter(.data$FEAT_NAME == "Reef") |>
  sf::st_transform(terra::crs(template))
reefs$reef_index <- seq_len(nrow(reefs))
reef_area_km2 <- as.numeric(sf::st_area(reefs)) / 1e6
reef_id <- terra::classify(terra::rast(file.path(cfg$dirs$cache, "reef_id_120m.tif")), cbind(0, NA))
expanded_path <- file.path(cfg$dirs$cache, "reef_id_expanded_120m.tif")
if (!file.exists(expanded_path)) {
  expanded <- reef_id
  for (i in 1:5) expanded <- terra::focal(expanded, w = 3, fun = "modal", na.policy = "only", na.rm = TRUE)
  terra::writeRaster(expanded, expanded_path, overwrite = TRUE, datatype = "INT2U")
}
expanded <- terra::rast(expanded_path)
own_inside <- terra::extract(reef_id, grid$cell)[[1]]
own_id <- terra::extract(expanded, grid$cell)[[1]]
log_message("Cells inside a reef polygon: ", round(100 * mean(!is.na(own_inside)), 1),
            "%; within 600 m of one: ", round(100 * mean(!is.na(own_id)), 1), "%")

# --- Edges, gaps and orientation --------------------------------------------
log_message("Reef edge cells")
# A reef cell is an edge if any 4-neighbour is non-reef or another reef.
id_matrix <- terra::as.matrix(reef_id, wide = TRUE)
nr <- nrow(id_matrix)
nc <- ncol(id_matrix)
differs <- function(a, b) is.na(b) | a != b
is_edge <- matrix(FALSE, nr, nc)
is_edge[-1, ] <- is_edge[-1, ] | differs(id_matrix[-1, ], id_matrix[-nr, ])
is_edge[-nr, ] <- is_edge[-nr, ] | differs(id_matrix[-nr, ], id_matrix[-1, ])
is_edge[, -1] <- is_edge[, -1] | differs(id_matrix[, -1], id_matrix[, -nc])
is_edge[, -nc] <- is_edge[, -nc] | differs(id_matrix[, -nc], id_matrix[, -1])
is_edge <- is_edge & !is.na(id_matrix)
edge_rc <- which(is_edge, arr.ind = TRUE)
rm(is_edge)
edge_cells <- terra::cellFromRowCol(template, edge_rc[, 1], edge_rc[, 2])
edge_xy <- terra::xyFromCell(template, edge_cells)
edge_id <- id_matrix[edge_rc]
query_xy <- as.matrix(grid[, c("x", "y")])
log_message("Edge cells: ", format(length(edge_cells), big.mark = ","), "; nearest-edge distances")
nn_edge <- FNN::get.knnx(edge_xy, query_xy, k = 1)
edge_distance <- as.numeric(nn_edge$nn.dist)
dist_to_reef_edge_m <- ifelse(is.na(own_inside), -edge_distance, edge_distance)
log_message("Inter-reef gaps at edge cells")
edge_gap <- edge_gap_distance(edge_xy, edge_id)
# Distance to a channel: reach the nearest edge, then cross the gap there.
nearest_other_reef_km <- (edge_distance + edge_gap[as.integer(nn_edge$nn.index)]) / 1000
log_message("Orientation within reefs")

assigned_cells <- which(!is.na(terra::values(expanded, mat = FALSE)))
assigned_xy <- terra::xyFromCell(template, assigned_cells)
assigned_id <- terra::values(expanded, mat = FALSE)[assigned_cells]
axis_position <- function(projection_all, projection_query) {
  lo <- tapply(projection_all, assigned_id, min)
  hi <- tapply(projection_all, assigned_id, max)
  key <- as.character(own_id)
  span <- hi[key] - lo[key]
  ifelse(is.na(own_id) | span <= 0, NA_real_, (projection_query - lo[key]) / span)
}
# SE unit vector (1, -1)/sqrt(2); NE unit vector (1, 1)/sqrt(2).
windward_position <- axis_position((assigned_xy[, 1] - assigned_xy[, 2]) / sqrt(2),
                                   (grid$x - grid$y) / sqrt(2))
seaward_position <- axis_position((assigned_xy[, 1] + assigned_xy[, 2]) / sqrt(2),
                                  (grid$x + grid$y) / sqrt(2))
log_reef_area_km2 <- log(reef_area_km2[own_id])

# --- Upwind shelter: step SE (row + k, col + k) until a different reef ------
log_message("Ray-casting upwind (SE) shelter")
rc <- terra::rowColFromCell(template, grid$cell)
max_steps <- 60L
first_hit <- rep(NA_integer_, n)
own_for_ray <- ifelse(is.na(own_inside), -1L, own_inside)
for (k in seq_len(max_steps)) {
  open <- is.na(first_hit)
  r <- rc[open, 1] + k
  c <- rc[open, 2] + k
  valid <- r <= nrow(id_matrix) & c <= ncol(id_matrix)
  hit_id <- rep(NA_integer_, sum(open))
  hit_id[valid] <- id_matrix[cbind(r[valid], c[valid])]
  hit <- !is.na(hit_id) & hit_id != own_for_ray[open]
  first_hit[which(open)[hit]] <- k
}
upwind_open_water_km <- ifelse(is.na(first_hit), max_steps, first_hit) * 120 * sqrt(2) / 1000
rm(id_matrix)

# --- eReefs hydrodynamics ------------------------------------------------------
log_message("eReefs GBR1 currents")
currents <- fill_nearby(ereefs_current_layers(file.path(cfg$dirs$cache, "ereefs")), iterations = 3L)
currents_120 <- terra::project(currents, template, method = "bilinear")
current_values <- terra::extract(currents_120, grid$cell)

# --- Assemble --------------------------------------------------------------
context <- data.frame(
  dist_to_reef_edge_m, windward_position, seaward_position, upwind_open_water_km,
  nearest_other_reef_km, log_reef_area_km2, current_values
)
layers <- terra::rast(template, nlyrs = ncol(context))
names(layers) <- names(context)
for (j in seq_along(context)) {
  v <- rep(NA_real_, terra::ncell(template))
  v[grid$cell] <- context[[j]]
  layers[[j]] <- terra::setValues(layers[[j]], v)
}
terra::writeRaster(layers, out_path, overwrite = TRUE, datatype = "FLT4S", gdal = layer_gdal_options)
log_message("Wrote ", out_path)

summary <- tibble::tibble(
  predictor = names(context),
  coverage = vapply(context, function(x) mean(is.finite(x)), numeric(1)),
  q05 = vapply(context, function(x) stats::quantile(x, 0.05, na.rm = TRUE), numeric(1)),
  median = vapply(context, function(x) stats::median(x, na.rm = TRUE), numeric(1)),
  q95 = vapply(context, function(x) stats::quantile(x, 0.95, na.rm = TRUE), numeric(1))
)
write_table(summary, cfg, "reef_context_120m_summary.csv")
print(summary)
