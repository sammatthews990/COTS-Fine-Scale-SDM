library(terra)
library(sf)
library(leaflet)
library(leafem)
library(raster)
library(igraph)
library(units)
library(ggplot2)
library(smoothr)
library(dplyr)
library(gdistance)

coralSDM <- terra::rast("data/maxent_predrast_GBR_AhyaD_015_lq2.tif")
crs(coralSDM)
bb_wgs <- ext(146.6, 146.86, -18.3, -18.15)
bb_wgs_poly <- as.polygons(bb_wgs, crs = "EPSG:4326")
bb_mga_poly <- project(bb_wgs_poly, crs(coralSDM))
coralSDM_crop <- crop(coralSDM, bb_mga_poly)
plot(coralSDM_crop)

r <- coralSDM_crop
summary(values(r))
crs_ <- st_crs(crs(r))
r <- terra::aggregate(coralSDM_crop, fact = 3, fun = mean)
alpha = 1.5
cond <- app(r, function(x) (pmax(0, pmin(1, x)))^alpha)
cond[is.na(cond)] <- 0
plot(cond)
hist(values(cond), main = "Conductance Values", xlab = "Conductance")
cond_r <- raster::raster(cond)
tr <- gdistance::transition(cond_r, function(x) mean(x, na.rm = TRUE), directions = 8)
tr <- gdistance::geoCorrection(tr, type = "c")
snap_to_tr <- function(pt_sp, tr, rRL) {
  stopifnot(inherits(pt_sp, "SpatialPoints"))
  valid_cells <- which(raster::values(rRL) > 0)
  if (!length(valid_cells)) stop("No valid conductance cells in raster.")
  cell0 <- raster::cellFromXY(rRL, sp::coordinates(pt_sp))
  ok <- is.finite(cell0) && (cell0 %in% valid_cells)
  if (!ok) {
    valid_xy <- raster::xyFromCell(rRL, valid_cells)
    xy <- sp::coordinates(pt_sp)
    i <- which.min((valid_xy[, 1] - xy[1, 1])^2 + (valid_xy[, 2] - xy[1, 2])^2)
    pt_sp <- SpatialPoints(matrix(valid_xy[i, ], ncol = 2), proj4string = raster::crs(rRL))
  }
  pt_sp
}

v <- values(r, mat = FALSE, na.rm = TRUE)
v <- ifelse(v > 0.2, v, NA)
o <- sort(v, decreasing = TRUE)
cs <- cumsum(o) / sum(o)
t50 <- o[which(cs >= 0.50)[1]]
top50 <- r >= t50
top50 <- ifel(top50, 1, NA)
poly_top50 <- as.polygons(top50, dissolve = TRUE) |> st_as_sf()
sf_top50 <- st_as_sf(poly_top50)
sf_top50 <- st_transform(sf_top50, crs_r)
drop_small <- function(sf_obj, min_area_m2 = 2000) {
  sf_obj %>% mutate(a = st_area(geometry)) %>% filter(a >= set_units(min_area_m2, m^2)) %>% st_make_valid()
}
sf_top50 <- drop_small(sf_top50, 2000)

reef_line <- sf_top50 |>
  st_cast("MULTILINESTRING") |>
  st_cast("LINESTRING")
pts <- st_line_sample(reef_line, n = 1, type = "regular") |>
  st_cast("POINT") |> as("Spatial")
pts_snapped <- lapply(seq_len(length(pts)), function(i) {
  snap_to_tr(pts[i, ], tr, cond_r)
})
pts_sf <- st_as_sf(do.call(rbind, pts_snapped))
st_crs(pts_sf) <- st_crs(reef_crs)
reef_centroid <- st_centroid(st_union(reef_line))
pts_m <- st_transform(pts_sf, reef_crs)
coords <- st_coordinates(pts_m)
db <- dbscan(coords, eps = 200, minPts = 3)
pts_m$cluster <- db$cluster
vals <- terra::extract(r, terra::vect(pts_m))
prob_col <- setdiff(names(vals), "ID")[1]
pts_m$prob <- vals[[prob_col]]
pts_m$prob[is.na(pts_m$prob)] <- 0
anchors <- pts_m |>
  group_by(cluster) |>
  slice_max(order_by = prob, n = 1, with_ties = FALSE) |>
  ungroup()
A <- st_coordinates(anchors)
CM <- as.matrix(costDistance(tr, A, A))
diag(CM) <- Inf
tsp <- TSP::TSP(CM)
ctr <- st_coordinates(st_transform(reef_centroid, reef_crs))
theta <- atan2(A[, 2] - ctr[2], A[, 1] - ctr[1])
init_order <- order(theta)
init_tour <- TOUR(init_order)
tour <- solve_TSP(tsp, method = "two_opt", control = list(tour = init_tour))
ord <- as.integer(tour)
edges <- cbind(ord, c(ord[-1], ord[1]))
seg_sf <- purrr::map(seq_len(nrow(edges)), \(k) {
  i <- edges[k, 1]; j <- edges[k, 2]
  sp <- shortestPath(tr, A[i, ], A[j, ], output = "SpatialLines")
})
coords_list <- lapply(seg_sf, function(x) {
  if (!is.null(x)) {
    sp::coordinates(x)[[1]][[1]]
  } else {
    NULL
  }
})
coords_list <- Filter(Negate(is.null), coords_list)
all_coords <- do.call(rbind, coords_list)
path_line <- st_linestring(all_coords)
path_line_sf <- st_sf(geometry = st_sfc(path_line), crs = reef_crs)
route <- smooth(path_line_sf, method = "ksmooth", smoothness = 10)

path_50 <- path_line_sf
path_50 <- st_transform(path_50, crs(r))
path_50_buffer <- sf::st_buffer(path_50, dist = 200)
buffer_raster <- terra::rasterize(path_50_buffer, r)
buffer_raster_aligned <- terra::project(buffer_raster, r, method = "near")
r.masked <- terra::mask(r, buffer_raster_aligned, maskvalues = 1, inverse = F)
plot(r.masked)
buffer_raster_aligned <- terra::project(buffer_raster, cond, method = "near")
cond_30_masked <- terra::mask(cond, buffer_raster_aligned, maskvalues = 1, inverse = F)
plot(cond_30_masked)
cond_30_r <- raster::raster(cond_30_masked)
tr2 <- gdistance::transition(cond_30_r, function(x) mean(x, na.rm = TRUE), directions = 8)
tr2 <- gdistance::geoCorrection(tr2, type = "c")
r <- r.masked
crs_r <- st_crs(crs(r))
v <- values(r, mat = FALSE, na.rm = TRUE)
v <- ifelse(v > 0.1, v, NA)
o <- sort(v, decreasing = TRUE)
cs <- cumsum(o) / sum(o)
t50 <- o[which(cs >= 0.50)[1]]
top50 <- r >= t50
top50 <- ifel(top50, 1, NA)
poly_top50 <- as.polygons(top50, dissolve = TRUE) |> st_as_sf()
sf_top50 <- st_as_sf(poly_top50)
sf_top50 <- st_transform(sf_top50, crs_r)
drop_small <- function(sf_obj, min_area_m2 = 2000) {
  sf_obj %>% mutate(a = st_area(geometry)) %>% filter(a >= set_units(min_area_m2, m^2)) %>% st_make_valid()
}
sf_top50 <- drop_small(sf_top50, 2000)

library(sf)
library(dplyr)
library(dbscan)
library(gdistance)
library(TSP)
library(purrr)
library(smoothr)

reef_line.next <- sf_top50 |>
  st_cast("MULTILINESTRING") |>
  st_cast("LINESTRING")
pts.next <- st_line_sample(reef_line.next, n = 1, type = "regular") |>
  st_cast("POINT") |> as("Spatial")
pts_snapped.next <- lapply(seq_len(length(pts.next)), function(i) {
  snap_to_tr(pts.next[i, ], tr2, cond_30_r)
})
pts_sf.next <- st_as_sf(do.call(rbind, pts_snapped.next))
st_crs(pts_sf.next) <- st_crs(reef_crs)
eps <- 1e-6
cond.mask <- terra::ifel(is.na(cond_30_masked), eps, cond_30_masked)
cond.mask <- terra::focal(cond.mask, w = 3, fun = max, na.policy = "omit")
cond_maskr <- raster::raster(cond.mask)
tr2 <- gdistance::transition(cond_maskr, function(x) mean(x), directions = 8)
tr2 <- gdistance::geoCorrection(tr2, type = "c")
if (terra::is.lonlat(r)) {
  r <- terra::project(r, as.character(reef_crs))
  pts_m <- st_transform(pts_sf.next, reef_crs)
} else {
  pts_m <- st_transform(pts_sf.next, crs(r))
}
pts_m <- st_transform(pts_sf.next, reef_crs)
coords <- st_coordinates(pts_m)
db <- dbscan(coords, eps = 200, minPts = 3)
pts_m$cluster <- db$cluster
vals <- terra::extract(r, terra::vect(pts_m))
prob_col <- setdiff(names(vals), "ID")[1]
pts_m$prob <- vals[[prob_col]]
pts_m$prob[is.na(pts_m$prob)] <- 0
anchors2 <- pts_m |>
  group_by(cluster) |>
  slice_max(order_by = prob, n = 1, with_ties = FALSE) |>
  ungroup()
A <- st_coordinates(anchors2)
CM <- as.matrix(costDistance(tr2, A, A))
diag(CM) <- Inf
tsp <- TSP::TSP(CM)
ctr <- st_coordinates(st_transform(reef_centroid, reef_crs))
theta <- atan2(A[, 2] - ctr[2], A[, 1] - ctr[1])
init_order <- order(theta)
init_tour <- TOUR(init_order)
tour <- solve_TSP(tsp, method = "two_opt", control = list(tour = init_tour))
ord <- as.integer(tour)
edges <- cbind(ord, c(ord[-1], ord[1]))
seg_sf <- purrr::map(seq_len(nrow(edges)), \(k) {
  i <- edges[k, 1]; j <- edges[k, 2]
  sp <- shortestPath(tr2, A[i, ], A[j, ], output = "SpatialLines")
})
coords_list2 <- lapply(seg_sf, function(x) {
  if (!is.null(x)) {
    sp::coordinates(x)[[1]][[1]]
  } else {
    NULL
  }
})
coords_list2 <- Filter(Negate(is.null), coords_list2)
all_coords2 <- do.call(rbind, coords_list2)
path_line2 <- st_linestring(all_coords2)
path_line2_sf <- st_sf(geometry = st_sfc(path_line2), crs = reef_crs)
route2 <- smooth(path_line2_sf, method = "ksmooth", smoothness = 10)

path_line_sf <- st_transform(path_line_sf, 4326)
path_smooth <- st_transform(route, 4326)
path_smooth.next <- st_transform(path_smooth.next, 4326)
pts_combined <- do.call(rbind, pts_snapped)
pts_sf <- st_as_sf(pts_combined)
pts_sf <- st_transform(pts_sf, 4326)
pts_sf.next <- st_transform(pts_sf.next, 4326)
pts_leaf <- st_transform(anchors, 4326)
path.next_leaf <- st_transform(route2, 4326)
pts.next_leaf <- st_transform(anchors2[init_order, ], 4326)
m <- leaflet() |>
  addProviderTiles(providers$Esri.WorldImagery, group = "Esri Satellite") |>
  addProviderTiles(providers$Esri.OceanBasemap, group = "Esri Ocean") |>
  addRasterImage(
    coralSDM_crop,
    colors  = pal_coral,
    opacity = 0.7,
    project = TRUE,
    group   = "CoralSDM"
  ) |>
  addRasterImage(
    cond_r,
    opacity = 0.7,
    project = TRUE,
    group   = "Conductance"
  ) |>
  addPolylines(
    data  = path_line_sf,
    color = "orange",
    weight = 3,
    group = "Tow Path"
  ) |>
  addPolylines(
    data  = path_smooth,
    color = "springgreen",
    weight = 3,
    group = "Tow Path - Smooth"
  ) |>
  addPolylines(
    data  = path.next_leaf,
    color = "skyblue",
    weight = 3,
    group = "Tow Path - Smooth 2nd"
  ) |>
  addCircleMarkers(
    data = pts_leaf,
    radius = 12,
    stroke = TRUE,
    color = "black",
    fill = TRUE,
    fillColor = "white",
    fillOpacity = 1,
    label = ~rownames(pts_leaf),
    labelOptions = labelOptions(
      noHide = TRUE, direction = "center", textsize = 25,
      textOnly = TRUE
    ),
    group = "Snapped Points"
  ) |>
  addCircleMarkers(
    data = pts.next_leaf,
    radius = 12,
    stroke = TRUE,
    color = "blue",
    fill = TRUE,
    fillColor = "white",
    fillOpacity = 1,
    label = ~rownames(pts.next_leaf),
    labelOptions = labelOptions(
      noHide = TRUE, direction = "center", textsize = 25,
      textOnly = TRUE
    ),
    group = "Snapped Points 2nd"
  ) |>
  addLegend(
    pal = pal_coral, values = rng, title = "Probability",
    position = "bottomright", opacity = 1
  ) |>
  addLayersControl(
    baseGroups    = c("Esri Satellite", "Esri Ocean"),
    overlayGroups = c(
      "CoralSDM", "Conductance", "Tow Path",
      "Tow Path - Smooth", "Snapped Points",
      "Snapped Points 2nd", "Tow Path - Smooth 2nd"
    ),
    options = layersControlOptions(collapsed = FALSE)
  ) |>
  leafem::addMouseCoordinates()
m
