# Exploration zones: high-ranked habitat outside the designated cull sites and
# away from recent manta tows, so search effort does not only consolidate on
# the same sites and tow paths.
#
# Part 1, check (held-out tows from 17): when tows went over new ground (no
# tow within 300 m in the previous 5 seasons), did the habitat models still
# rank strong-evidence tows above negatives within the reef-year? And is new
# ground more or less productive than repeat ground?
#
# Part 2, map (reefs with designated cull sites):
#   habitat   ensemble of the two layers best at locating COTS evidence: the
#             cull-trained primary layer and the within reef-year manta model
#             (mean of within-reef percentiles, re-ranked within the reef;
#             primary layer alone where the manta layer is missing),
#             smoothed with a 5 x 5 cell (150 m) mean so zones read as ground
#             to search rather than scattered cells
#   eligible  reef cell 0.5-16 m deep (clean DEM), > 100 m from any designated
#             cull site and > 150 m from any manta tow in the last 5 seasons
#   outputs   cots_exploration_2027_30m.tif
#               band 1 exploration score: smoothed habitat percentile on
#               eligible cells
#               band 2 distance to the nearest manta tow, last 5 seasons (m)
#               band 3 distance to the nearest designated cull site (m)
#             exploration zones: all eligible ground in the reef's top 40, 30,
#             20 and 10 % of habitat, as nested shaded bands (outlines closed
#             with a +30 m / -30 m buffer, simplified to 20 m to cut the cell
#             stairs, Chaikin corner-cutting; holes < 0.5 ha filled, pieces
#             < 0.5 ha dropped, clipped to the reef and outside cull sites,
#             simplified to 8 m). Kernel smoothing took ~35 s per reef on
#             these large bands; Chaikin takes a few seconds. Each band lies inside the one
#             below it, so stacked semi-transparent fills darken towards the
#             best habitat. Added as an "Exploration zones" folder to the 2027
#             cull-site KMZs (20, 20b) and written to a GeoPackage and a
#             per-reef CSV.
# Zones replace the earlier top-10 discrete areas per reef: nothing separated
# those areas from the surrounding unsearched habitat, so the whole of it is
# shown, graded by habitat rank.
# The manta layer is extrapolated away from tow paths, which is exactly where
# these zones are; treat them as places to look, not predictions to trust.

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "..", "R", "common.R"))
cfg <- initialise_final()
source(file.path(cfg$final_dir, "R", "cull_sites.R"))
source(file.path(cfg$final_dir, "R", "site_selection.R"))
sr <- cfg$site_ranking
set.seed(cfg$hmsc$seed + 21L)
new_ground_m <- 300
lookback_seasons <- 5L
site_clearance_m <- 100
tow_clearance_m <- 150
depth_range <- c(-16, -0.5)
zone_floors <- c(0.6, 0.7, 0.8, 0.9)   # top 40, 30, 20, 10 % of the reef's habitat
smooth_cells <- 5   # 150 m mean, so zones read as ground rather than scattered cells
hole_ha <- 0.5
min_piece_ha <- 0.5
simplify_m <- 8   # smoothing adds many vertices; 8 m keeps the shape (area change < 1 %)
boot_reps <- 500L
out_dir <- file.path(cfg$dirs$outputs, "site_ranking")

# --- Part 1: do the habitat models find COTS on new ground? -----------------------
tow_scores_path <- file.path(cfg$dirs$cache, "manta_response_tow_scores.rds")
if (!file.exists(tow_scores_path)) stop("Run 17_manta_response_model.R first (needs ", tow_scores_path, ")")
tw <- readRDS(tow_scores_path) |>
  dplyr::filter(is.finite(.data$mx), is.finite(.data$cull_cv), is.finite(.data$within_cv))
tw$prior_dist <- NA_real_
for (T in sort(unique(tw$year))) {
  now <- which(tw$year == T)
  prior <- which(tw$year >= T - lookback_seasons & tw$year < T)
  if (length(prior)) {
    tw$prior_dist[now] <- FNN::get.knnx(cbind(tw$mx[prior], tw$my[prior]), cbind(tw$mx[now], tw$my[now]), k = 1)$nn.dist[, 1]
  }
}
tw <- tw |>
  dplyr::filter(!is.na(.data$prior_dist), !is.na(.data$label)) |>
  dplyr::mutate(ground = ifelse(.data$prior_dist > new_ground_m, "new ground", "repeat ground")) |>
  dplyr::group_by(.data$reef_year) |>
  dplyr::mutate(ensemble = pct_rank(pct_rank(.data$cull_cv) + pct_rank(.data$within_cv))) |>
  dplyr::ungroup()
pieces <- function(d, col) {
  dplyr::bind_rows(lapply(split(seq_len(nrow(d)), d$reef_year), function(i) {
    y <- d$label[i] == 1
    n1 <- sum(y)
    n0 <- sum(!y)
    if (n1 == 0 || n0 == 0) return(NULL)
    r <- rank(d[[col]][i])
    top <- rank(-d[[col]][i], ties.method = "random") <= ceiling(0.2 * length(i))
    data.frame(auc = (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0), pairs = n1 * n0, hits = sum(y & top), pos = n1)
  }))
}
summ <- function(p) c(auc = stats::weighted.mean(p$auc, p$pairs), lift = sum(p$hits) / sum(p$pos) / 0.2)
new_ground_check <- dplyr::bind_rows(lapply(c("new ground", "repeat ground"), function(g) {
  d <- tw[tw$ground == g, ]
  dplyr::bind_rows(lapply(c(cull = "cull_cv", manta_within = "within_cv", ensemble = "ensemble"), function(col) {
    p <- pieces(d, col)
    point <- summ(p)
    boot <- replicate(boot_reps, summ(p[sample.int(nrow(p), replace = TRUE), ]))
    tibble::tibble(ground = g, score = sub("_cv$", "", col), metric = names(point), value = point,
                   low = apply(boot, 1, stats::quantile, 0.025), high = apply(boot, 1, stats::quantile, 0.975),
                   tows = nrow(d), strong_share = mean(d$label == 1), reef_years_with_both = nrow(p))
  }))
}))
new_ground_check$score[new_ground_check$score == "within"] <- "manta_within"
write_table(new_ground_check, cfg, "exploration_new_ground_check.csv")
print(new_ground_check, n = Inf)

# --- Part 2: exploration zones on reefs with designated cull sites ---------------
sites_sf <- read_designated_sites(cfg)
template <- terra::rast(file.path(cfg$dirs$cache, "reef_id_30m.tif"))
reefs_sf <- sf::st_read(checked_input(cfg, "reef_features"), quiet = TRUE) |>
  dplyr::filter(.data$FEAT_NAME == "Reef") |>
  sf::st_transform(terra::crs(template))
reefs_sf$reef_index <- seq_len(nrow(reefs_sf))
site_xy <- sf::st_coordinates(suppressWarnings(sf::st_point_on_surface(sf::st_geometry(sites_sf))))[, 1:2]
site_reef <- terra::extract(template, site_xy)[[1]]
reef_map <- unique(data.frame(reef_name = sites_sf$reef_name, reef_index = site_reef))
reef_map <- reef_map[!is.na(reef_map$reef_index) & reef_map$reef_index >= 1 & reef_map$reef_index <= nrow(reefs_sf), ]
reef_map <- reef_map[!sf::st_is_empty(reefs_sf[reef_map$reef_index, ]), ]
log_message(dplyr::n_distinct(reef_map$reef_name), " cull-site reefs on ", nrow(reef_map), " GBR reef polygons")

cull_pct <- terra::rast(file.path(cfg$dirs$cache, "search_planning_percentile_30m.tif"))
manta_pct <- terra::rast(file.path(cfg$dirs$layers, "cots_manta_within_30m.tif"))[["within_reef_percentile"]]
dem <- terra::rast(checked_input(cfg, "predictor_stack_clean"))[["DEM"]]

geom <- readRDS(checked_input(cfg, "footprint_links"))$manta |>
  dplyr::group_by(.data$visit_id) |>
  dplyr::summarise(sx = dplyr::first(.data$start_x), sy = dplyr::first(.data$start_y),
                   ex = dplyr::first(.data$end_x), ey = dplyr::first(.data$end_y), .groups = "drop")
recent_tows <- readRDS(checked_input(cfg, "manta_tows")) |>
  dplyr::select("visit_id", "survey_date") |>
  dplyr::inner_join(geom, by = "visit_id") |>
  dplyr::filter(is.finite(.data$sx), is.finite(.data$ex),
                .data$survey_date > sr$evidence_end - 365 * lookback_seasons, .data$survey_date <= sr$evidence_end)
tow_lines <- terra::vect(sf::st_sfc(lapply(seq_len(nrow(recent_tows)), function(i)
  sf::st_linestring(matrix(c(recent_tows$sx[i], recent_tows$ex[i], recent_tows$sy[i], recent_tows$ey[i]), 2))),
  crs = sf::st_crs(reefs_sf)))
site_vect <- terra::vect(sites_sf)
log_message(format(nrow(recent_tows), big.mark = ","), " manta tows in the last ", lookback_seasons, " seasons")

# Rounded outline of a set of cells, clipped to the reef and outside cull
# sites, without pieces too small to matter.
zone_outline <- function(raw, reef_geom, site_union) {
  usable <- function(g) length(g) && !all(sf::st_is_empty(g)) && !anyNA(sf::st_coordinates(g))
  g <- sf::st_simplify(sf::st_buffer(sf::st_buffer(sf::st_make_valid(raw), 30), -30), dTolerance = 20)   # close gaps, cut the 30 m stairs
  if (!usable(g)) g <- raw
  smoothed <- tryCatch(smoothr::smooth(g, method = "chaikin", refinements = 3), error = function(e) NULL)
  if (!is.null(smoothed) && usable(smoothed)) g <- sf::st_make_valid(smoothed)
  g <- smoothr::fill_holes(g, threshold = units::set_units(hole_ha, "ha"))
  g <- sf::st_intersection(g, reef_geom)
  if (length(site_union)) g <- sf::st_difference(g, site_union)
  g <- sf::st_simplify(g, dTolerance = simplify_m, preserveTopology = TRUE)
  g <- sf::st_cast(sf::st_collection_extract(sf::st_make_valid(g), "POLYGON"), "POLYGON")
  g <- g[as.numeric(sf::st_area(g)) / 1e4 >= min_piece_ha]
  if (!length(g)) return(NULL)
  sf::st_union(g)
}

cells <- list()
zones <- list()
for (i in seq_len(nrow(reef_map))) {
  idx <- reef_map$reef_index[i]
  poly <- terra::vect(reefs_sf[idx, ])
  frame <- terra::ext(terra::buffer(poly, 500))
  rid <- terra::crop(template, frame)
  in_reef <- terra::values(rid, mat = FALSE) == idx
  in_reef[is.na(in_reef)] <- FALSE
  if (sum(in_reef) < 20) next
  cp <- terra::values(terra::crop(cull_pct, frame), mat = FALSE)
  mp <- terra::values(terra::crop(manta_pct, frame), mat = FALSE)
  z <- terra::values(terra::crop(dem, frame), mat = FALSE)
  hab <- rep(NA_real_, length(cp))
  both <- in_reef & is.finite(cp) & is.finite(mp)
  hab[both] <- (cp[both] + mp[both]) / 2
  only_cull <- in_reef & is.finite(cp) & !is.finite(mp)
  hab[only_cull] <- cp[only_cull]
  ok <- in_reef & is.finite(hab)
  hab[ok] <- pct_rank(hab[ok])
  # Smoothed habitat (smooth_cells-square mean over reef cells), re-ranked within the reef.
  hab_r <- rid
  terra::values(hab_r) <- hab
  hs <- terra::values(terra::focal(hab_r, matrix(1, smooth_cells, smooth_cells), fun = mean, na.rm = TRUE), mat = FALSE)
  hs[!ok] <- NA
  hs[ok] <- pct_rank(hs[ok])
  # Distances to designated sites and recent tows (NA-cells are the targets).
  near_sites <- terra::crop(site_vect, frame)
  d_site <- if (length(near_sites)) terra::values(terra::distance(terra::rasterize(near_sites, rid, touches = TRUE)), mat = FALSE) else rep(Inf, length(cp))
  near_tows <- terra::crop(tow_lines, frame)
  d_tow <- if (length(near_tows)) terra::values(terra::distance(terra::rasterize(near_tows, rid, touches = TRUE)), mat = FALSE) else rep(Inf, length(cp))
  eligible <- ok & is.finite(z) & z >= depth_range[1] & z <= depth_range[2] &
    d_site > site_clearance_m & d_tow > tow_clearance_m
  global <- terra::cellFromXY(template, terra::xyFromCell(rid, which(ok)))
  cells[[i]] <- data.frame(cell = global, score = ifelse(eligible[ok], hs[ok], NA_real_),
                           d_tow = pmin(d_tow[ok], 1e5), d_site = pmin(d_site[ok], 1e5))
  # Nested zones: eligible ground above each habitat floor, each band kept
  # inside the one below it.
  reef_geom <- sf::st_geometry(reefs_sf[idx, ])
  frame_sf <- sf::st_geometry(sf::st_as_sf(terra::as.polygons(frame, crs = terra::crs(rid))))
  site_union <- sf::st_union(sf::st_geometry(sites_sf)[lengths(sf::st_intersects(sites_sf, frame_sf)) > 0])
  outer <- NULL
  reef_zones <- list()
  for (b in zone_floors) {
    m <- rid
    terra::values(m) <- ifelse(eligible & hs >= b, 1, NA)
    if (!any(terra::values(m, mat = FALSE) == 1, na.rm = TRUE)) break
    raw <- sf::st_geometry(sf::st_as_sf(terra::as.polygons(m, dissolve = TRUE)))
    g <- zone_outline(raw, reef_geom, site_union)
    if (!is.null(g) && !is.null(outer)) {
      g <- sf::st_collection_extract(sf::st_make_valid(sf::st_intersection(g, outer)), "POLYGON")
      g <- if (length(g)) sf::st_union(g) else NULL
    }
    if (is.null(g)) break
    outer <- g
    pieces <- sf::st_cast(g, "POLYGON")
    piece_ha <- as.numeric(sf::st_area(pieces)) / 1e4
    zm <- ok & eligible & hs >= b
    reef_zones[[length(reef_zones) + 1L]] <- sf::st_sf(
      reef_name = reef_map$reef_name[i], habitat_floor = b, zone = sprintf("Top %d %% habitat", round(100 * (1 - b))),
      area_ha = sum(piece_ha), pieces = length(pieces), largest_piece_ha = max(piece_ha),
      mean_habitat_percentile = mean(hs[zm]), depth_min_m = -max(z[zm]), depth_max_m = -min(z[zm]),
      largest_piece = sf::st_geometry(pieces)[which.max(piece_ha)], geometry = g, sf_column_name = "geometry")
  }
  if (length(reef_zones)) zones[[i]] <- do.call(rbind, reef_zones)
  if (i %% 50 == 0) log_message(i, " of ", nrow(reef_map), " reef polygons")
}
cells <- dplyr::bind_rows(cells)
cells <- cells[!duplicated(cells$cell), ]
zones <- do.call(rbind, Filter(Negate(is.null), zones))
log_message(format(sum(is.finite(cells$score)), big.mark = ","), " eligible cells; exploration zones on ",
            dplyr::n_distinct(zones$reef_name), " reefs")

# --- Write the 30 m exploration layer --------------------------------------------
band_paths <- file.path(cfg$dirs$cache, paste0("exploration_", c("score", "dist_tow", "dist_site"), "_30m.tif"))
write_cells(template, cells$cell[is.finite(cells$score)], cells$score[is.finite(cells$score)], band_paths[1], "exploration_score")
write_cells(template, cells$cell, cells$d_tow, band_paths[2], "distance_recent_manta_tow_m")
write_cells(template, cells$cell, cells$d_site, band_paths[3], "distance_cull_site_m")
layer <- terra::rast(band_paths)
names(layer) <- c("exploration_score", "distance_recent_manta_tow_m", "distance_cull_site_m")
layer_path <- file.path(cfg$dirs$layers, "cots_exploration_2027_30m.tif")
write_layer(layer, layer_path)
register_layer(
  cfg, "exploration_2027", layer_path,
  description = paste("Exploration score: within-reef habitat percentile (primary layer + manta within reef-year ensemble,",
                      "5 x 5 cell mean) on reef cells 0.5-16 m deep, > 100 m from designated cull sites and > 150 m from",
                      "manta tows in the last 5 seasons; bands 2-3: distance to recent tows and to cull sites.",
                      "Reefs with designated cull sites only."),
  model = "ensemble of primary search-planning layer and manta within reef-year layer (scripts 11, 18, 21)",
  units = "percentile / m"
)

# --- Zones: GeoPackage, per-reef CSV and KMZ folder ---------------------------------
zones_ll <- sf::st_transform(zones, 4326)
centre <- sf::st_coordinates(suppressWarnings(sf::st_point_on_surface(sf::st_transform(zones$largest_piece, 4326))))
zone_table <- sf::st_drop_geometry(zones_ll) |>
  dplyr::select(-"largest_piece") |>
  dplyr::mutate(dplyr::across(c("area_ha", "largest_piece_ha", "depth_min_m", "depth_max_m"), ~ round(.x, 1)),
                mean_habitat_percentile = round(.data$mean_habitat_percentile, 3),
                largest_piece_longitude = round(centre[, 1], 5), largest_piece_latitude = round(centre[, 2], 5))
utils::write.csv(zone_table, file.path(out_dir, "cots_exploration_zones_2027.csv"), row.names = FALSE)
sf::st_write(dplyr::select(zones, "reef_name", "habitat_floor", "zone", "area_ha", "pieces", "mean_habitat_percentile"),
             file.path(out_dir, "cots_exploration_zones_2027.gpkg"), delete_dsn = TRUE, quiet = TRUE)
# The discrete areas these zones replace.
old <- file.path(out_dir, paste0("cots_exploration_areas_2027.", c("csv", "gpkg")))
if (any(file.exists(old))) file.remove(old[file.exists(old)])

xml_escape <- function(x) gsub(">", "&gt;", gsub("<", "&lt;", gsub("&", "&amp;", x)))
ring_text <- function(m) paste(sprintf("%.6f,%.6f,0", m[, 1], m[, 2]), collapse = " ")
polygon_kml <- function(g) {
  polys <- if (inherits(g, "MULTIPOLYGON")) unclass(g) else list(unclass(g))
  paste0(vapply(polys, function(p) {
    inner <- if (length(p) > 1) paste0(vapply(p[-1], function(r) paste0("<innerBoundaryIs><LinearRing><coordinates>", ring_text(r),
                                                                            "</coordinates></LinearRing></innerBoundaryIs>"), character(1)), collapse = "") else ""
    paste0("<Polygon><altitudeMode>clampToGround</altitudeMode><outerBoundaryIs><LinearRing><coordinates>",
           ring_text(p[[1]]), "</coordinates></LinearRing></outerBoundaryIs>", inner, "</Polygon>")
  }, character(1)), collapse = "")
}
style_id <- function(b) paste0("exploration-", round(100 * b))
placemark <- function(j) {
  a <- zone_table[j, ]
  fields <- c("Exploration zone" = a$zone, "Reef" = a$reef_name, "Area (ha)" = a$area_ha,
              "Depth (m)" = sprintf("%.1f-%.1f", a$depth_min_m, a$depth_max_m))
  rows <- paste0("<tr><td class='property-name'>", names(fields), "</td><td class='property-value'>", xml_escape(fields), "</td></tr>", collapse = "")
  data <- paste0("<Data name=\"", names(fields), "\"><value>", xml_escape(fields), "</value></Data>", collapse = "")
  paste0("<Placemark><name>", xml_escape(a$zone), "</name><styleUrl>#", style_id(a$habitat_floor), "</styleUrl>",
         "<description><![CDATA[<table>", rows, "</table><p>Habitat ranked in the reef's ", tolower(a$zone),
         ", outside designated cull sites and more than 150 m from any manta tow of the last 5 seasons. Darker purple is ",
         "better-ranked habitat. A place to look, not a prediction to trust.</p>]]></description>",
         "<ExtendedData>", data, "</ExtendedData><MultiGeometry>", polygon_kml(sf::st_geometry(zones_ll)[[j]]),
         "</MultiGeometry></Placemark>")
}
# Lowest band first, so better bands draw on top and the fills stack.
order_rows <- order(zone_table$reef_name, zone_table$habitat_floor)
folders <- vapply(split(order_rows, zone_table$reef_name[order_rows]), function(js) {
  paste0("<Folder><name>Reef: ", xml_escape(zone_table$reef_name[js[1]]), "</name><visibility>1</visibility>",
         paste0(vapply(js, placemark, character(1)), collapse = ""), "</Folder>")
}, character(1))
# Purple (KML aabbggrr); each stacked band adds about a quarter opacity, and
# the best band gets a thin outline.
styles <- paste0(vapply(zone_floors, function(b) paste0(
  "<Style id=\"", style_id(b), "\"><LineStyle><color>", if (b == max(zone_floors)) "ccffffff" else "00ffffff",
  "</color><width>", if (b == max(zone_floors)) 1 else 0, "</width></LineStyle><PolyStyle><color>40a73a4a</color>",
  "<fill>1</fill><outline>", as.integer(b == max(zone_floors)), "</outline></PolyStyle></Style>"), character(1)), collapse = "\n    ")
exploration_folder <- paste0(
  "<Folder><name>Exploration zones ", sr$target_season,
  " (purple, darker = better habitat): rarely searched ground outside cull sites</name>",
  "<open>0</open><visibility>1</visibility>", paste0(folders, collapse = ""), "</Folder>")

# Add the zones to each cull-site KMZ staged by 20 and 20b.
kmz_targets <- c(site_ranking_kmz = "cots_cull_site_priority_2027.kmz",
                 site_ranking_density_kmz = "cots_cull_site_priority_density_2027.kmz")
stage_kml <- file.path(cfg$dirs$cache, names(kmz_targets), "doc.kml")
if (!file.exists(stage_kml[1])) stop("Run 20_site_ranking_2027.R first (needs ", stage_kml[1], ")")
for (k in which(file.exists(stage_kml))) {
  kml <- paste(readLines(stage_kml[k], encoding = "UTF-8", warn = FALSE), collapse = "\n")
  kml <- sub("<Style id=", paste0(styles, "\n    <Style id="), kml, fixed = TRUE)
  kml <- sub("</Document>", paste0(exploration_folder, "\n  </Document>"), kml, fixed = TRUE)
  combined_dir <- file.path(cfg$dirs$cache, paste0(names(kmz_targets)[k], "_exploration"))
  dir.create(combined_dir, showWarnings = FALSE)
  writeLines(kml, file.path(combined_dir, "doc.kml"), useBytes = TRUE)
  kmz_path <- file.path(out_dir, kmz_targets[[k]])
  if (file.exists(kmz_path)) file.remove(kmz_path)
  utils::zip(kmz_path, file.path(combined_dir, "doc.kml"), flags = "-j9Xq")
  log_message("Wrote ", kmz_path, " (cull-site priority + exploration zones)")
}

print(zone_table |>
        dplyr::group_by(.data$zone) |>
        dplyr::summarise(reefs = dplyr::n(), total_ha = sum(.data$area_ha), median_ha = stats::median(.data$area_ha),
                         median_pieces = stats::median(.data$pieces), .groups = "drop"))
