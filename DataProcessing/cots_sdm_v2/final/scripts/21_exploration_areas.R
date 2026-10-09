# Exploration areas: high-ranked habitat outside the designated cull sites and
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
#             primary layer alone where the manta layer is missing)
#   eligible  reef cell 0.5-16 m deep (clean DEM), > 100 m from any designated
#             cull site and > 150 m from any manta tow in the last 5 seasons
#   outputs   cots_exploration_2027_30m.tif
#               band 1 exploration score: habitat percentile on eligible cells
#               band 2 distance to the nearest manta tow, last 5 seasons (m)
#               band 3 distance to the nearest designated cull site (m)
#             candidate areas: eligible cells in the reef's top 30 % of
#             habitat; one-cell gaps closed (dilate then erode, eligible
#             cells only) so nearby pockets join; patches over ~9 ha split
#             into compact chunks (k-means on cell centres, about one cull
#             site each); outlines rounded (+45 m / -30 m buffer, kernel
#             smoothing, small holes filled, simplified to 8 m) and clipped to the reef and
#             outside cull sites. Areas >= 0.5 ha are ranked by mean habitat
#             percentile and the top 10 per reef, at least 350 m apart, are
#             added as an "Exploration areas" folder to the 2027 cull-site
#             KMZ (20) and written to a CSV and GeoPackage.
# The manta layer is extrapolated away from tow paths, which is exactly where
# these areas are; treat them as places to look, not predictions to trust.

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
habitat_floor <- 0.7
min_area_ha <- 0.5
chunk_ha <- 9
min_spacing_m <- 350
areas_per_reef <- 10L
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

# --- Part 2: exploration surface on reefs with designated cull sites -------------
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

cells <- list()
areas <- list()
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
  # Distances to designated sites and recent tows (NA-cells are the targets).
  near_sites <- terra::crop(site_vect, frame)
  d_site <- if (length(near_sites)) terra::values(terra::distance(terra::rasterize(near_sites, rid, touches = TRUE)), mat = FALSE) else rep(Inf, length(cp))
  near_tows <- terra::crop(tow_lines, frame)
  d_tow <- if (length(near_tows)) terra::values(terra::distance(terra::rasterize(near_tows, rid, touches = TRUE)), mat = FALSE) else rep(Inf, length(cp))
  eligible <- ok & is.finite(z) & z >= depth_range[1] & z <= depth_range[2] &
    d_site > site_clearance_m & d_tow > tow_clearance_m
  global <- terra::cellFromXY(template, terra::xyFromCell(rid, which(ok)))
  cells[[i]] <- data.frame(cell = global, score = ifelse(eligible[ok], hab[ok], NA_real_),
                           d_tow = pmin(d_tow[ok], 1e5), d_site = pmin(d_site[ok], 1e5))
  # Candidate areas: eligible top-habitat cells, one-cell gaps closed
  # through eligible cells so nearby pockets join into one area.
  cand01 <- rid
  terra::values(cand01) <- as.numeric(eligible & hab >= habitat_floor)
  if (!any(terra::values(cand01, mat = FALSE) == 1, na.rm = TRUE)) next
  w3 <- matrix(1, 3, 3)
  closed <- terra::focal(terra::focal(cand01, w3, fun = "max", na.rm = TRUE), w3, fun = "min", na.rm = TRUE)
  closed_v <- terra::values(closed, mat = FALSE) == 1 & eligible
  closed_v[is.na(closed_v)] <- FALSE
  cand <- rid
  terra::values(cand) <- ifelse(closed_v, 1, NA)
  pv <- terra::values(terra::patches(cand, directions = 8), mat = FALSE)
  # Split large patches into compact, roughly site-sized chunks.
  cxy <- terra::xyFromCell(rid, seq_len(terra::ncell(rid)))
  cell_ha <- prod(terra::res(rid)) / 1e4
  group <- rep(NA_real_, length(pv))
  for (pid in unique(stats::na.omit(pv))) {
    members <- which(pv == pid)
    k <- ceiling(length(members) * cell_ha / chunk_ha)
    group[members] <- if (k <= 1) pid * 1000 else {
      set.seed(pid)
      pid * 1000 + stats::kmeans(cxy[members, ], centers = k, nstart = 3, iter.max = 50)$cluster
    }
  }
  keep <- !is.na(group)
  stats_df <- data.frame(patch = group[keep], hab = hab[keep], cp = cp[keep], mp = mp[keep], z = z[keep],
                         d_site = d_site[keep], d_tow = d_tow[keep], x = cxy[keep, 1], y = cxy[keep, 2]) |>
    dplyr::group_by(.data$patch) |>
    dplyr::summarise(cells = dplyr::n(), habitat_percentile = mean(.data$hab), weight = mean(.data$hab),
                     cull_layer_percentile = mean(.data$cp, na.rm = TRUE), manta_layer_percentile = mean(.data$mp, na.rm = TRUE),
                     depth_min_m = -max(.data$z), depth_max_m = -min(.data$z),
                     nearest_cull_site_m = min(.data$d_site), nearest_recent_tow_m = min(.data$d_tow),
                     cx = mean(.data$x), cy = mean(.data$y), .groups = "drop") |>
    dplyr::filter(.data$cells * cell_ha >= min_area_ha) |>
    dplyr::arrange(dplyr::desc(.data$weight), dplyr::desc(.data$cells))
  if (!nrow(stats_df)) next
  # Greedy pick: best chunks first, at least min_spacing_m apart.
  chosen <- integer(0)
  for (j in seq_len(nrow(stats_df))) {
    if (length(chosen) >= areas_per_reef) break
    if (!length(chosen) || all(sqrt((stats_df$cx[chosen] - stats_df$cx[j])^2 + (stats_df$cy[chosen] - stats_df$cy[j])^2) >= min_spacing_m)) {
      chosen <- c(chosen, j)
    }
  }
  stats_df <- stats_df[chosen, ]
  # Outlines: dissolve cells, round the staircase edges, clip to the reef
  # and keep outside the designated sites.
  groups_r <- rid
  terra::values(groups_r) <- ifelse(group %in% stats_df$patch, group, NA)
  polys_sf <- sf::st_as_sf(terra::as.polygons(groups_r, values = TRUE))
  names(polys_sf)[1] <- "patch"
  reef_geom <- sf::st_geometry(reefs_sf[idx, ])
  frame_sf <- sf::st_geometry(sf::st_as_sf(terra::as.polygons(frame, crs = terra::crs(rid))))
  site_union <- sf::st_union(sf::st_geometry(sites_sf)[lengths(sf::st_intersects(sites_sf, frame_sf)) > 0])
  usable <- function(g) length(g) && !all(sf::st_is_empty(g)) && !anyNA(sf::st_coordinates(g))
  geoms <- lapply(seq_len(nrow(polys_sf)), function(r) {
    raw <- sf::st_make_valid(sf::st_geometry(polys_sf)[r])
    g <- sf::st_buffer(sf::st_buffer(raw, 45), -30)
    if (!usable(g)) g <- raw
    smoothed <- tryCatch(smoothr::smooth(g, method = "ksmooth", smoothness = 3), error = function(e) NULL)
    if (!is.null(smoothed) && usable(smoothed)) g <- sf::st_make_valid(smoothed)
    g <- smoothr::fill_holes(g, threshold = units::set_units(1, "ha"))
    g <- sf::st_intersection(g, reef_geom)
    if (length(site_union)) g <- sf::st_difference(g, site_union)
    g <- sf::st_simplify(g, dTolerance = simplify_m, preserveTopology = TRUE)
    g <- sf::st_collection_extract(sf::st_make_valid(g), "POLYGON")
    if (!length(g)) return(sf::st_multipolygon())
    sf::st_union(g)[[1]]
  })
  polys_sf <- sf::st_sf(patch = polys_sf$patch, geometry = sf::st_sfc(geoms, crs = sf::st_crs(reefs_sf)))
  polys_sf <- polys_sf[!sf::st_is_empty(polys_sf), ]
  if (!nrow(polys_sf)) next
  polys_sf$area_ha <- as.numeric(sf::st_area(polys_sf)) / 1e4
  polys_sf <- polys_sf[polys_sf$area_ha >= min_area_ha, ] |>
    dplyr::left_join(dplyr::select(stats_df, -"cx", -"cy", -"cells"), by = "patch") |>
    dplyr::mutate(reef_name = reef_map$reef_name[i])
  areas[[i]] <- polys_sf
}
cells <- dplyr::bind_rows(cells)
cells <- cells[!duplicated(cells$cell), ]
areas <- do.call(rbind, Filter(Negate(is.null), areas)) |>
  dplyr::group_by(.data$reef_name) |>
  dplyr::arrange(dplyr::desc(.data$weight), .by_group = TRUE) |>
  dplyr::mutate(area_rank_on_reef = dplyr::row_number()) |>
  dplyr::ungroup()
log_message(format(sum(is.finite(cells$score)), big.mark = ","), " eligible cells; ", nrow(areas), " candidate areas on ",
            dplyr::n_distinct(areas$reef_name), " reefs")

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
  description = paste("Exploration score: within-reef habitat percentile (primary layer + manta within reef-year ensemble) on reef cells",
                      "0.5-16 m deep, > 100 m from designated cull sites and > 150 m from manta tows in the last 5 seasons;",
                      "bands 2-3: distance to recent tows and to cull sites. Reefs with designated cull sites only."),
  model = "ensemble of primary search-planning layer and manta within reef-year layer (scripts 11, 18, 21)",
  units = "percentile / m"
)

# --- Candidate areas: CSV and KMZ folder -------------------------------------------
areas_ll <- sf::st_transform(areas, 4326)
centre <- sf::st_coordinates(suppressWarnings(sf::st_point_on_surface(sf::st_geometry(areas_ll))))
area_table <- sf::st_drop_geometry(areas_ll) |>
  dplyr::transmute(reef_name = .data$reef_name, area_rank_on_reef = .data$area_rank_on_reef,
                   area_ha = round(.data$area_ha, 1), habitat_percentile = round(.data$habitat_percentile, 3),
                   cull_layer_percentile = round(.data$cull_layer_percentile, 3),
                   manta_layer_percentile = round(.data$manta_layer_percentile, 3),
                   depth_min_m = round(.data$depth_min_m, 1), depth_max_m = round(.data$depth_max_m, 1),
                   nearest_cull_site_m = round(.data$nearest_cull_site_m), nearest_recent_tow_m = round(pmin(.data$nearest_recent_tow_m, 1e5)),
                   centre_longitude = round(centre[, 1], 5), centre_latitude = round(centre[, 2], 5))
utils::write.csv(area_table, file.path(out_dir, "cots_exploration_areas_2027.csv"), row.names = FALSE)
sf::st_write(dplyr::select(areas, "reef_name", "area_rank_on_reef", "area_ha", "habitat_percentile",
                           "nearest_cull_site_m", "nearest_recent_tow_m"),
             file.path(out_dir, "cots_exploration_areas_2027.gpkg"), delete_dsn = TRUE, quiet = TRUE)

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
placemark <- function(j) {
  a <- area_table[j, ]
  fields <- c("Exploration area" = paste0(a$area_rank_on_reef, " of ", sum(area_table$reef_name == a$reef_name), " on reef"),
              "Reef" = a$reef_name, "Area (ha)" = a$area_ha,
              "Habitat percentile" = sprintf("%.0f%%", 100 * a$habitat_percentile),
              "Depth (m)" = sprintf("%.1f-%.1f", a$depth_min_m, a$depth_max_m),
              "Nearest cull site (m)" = a$nearest_cull_site_m,
              "Nearest manta tow, last 5 seasons (m)" = ifelse(a$nearest_recent_tow_m >= 1e5, "none on reef", a$nearest_recent_tow_m))
  rows <- paste0("<tr><td class='property-name'>", names(fields), "</td><td class='property-value'>", xml_escape(fields), "</td></tr>", collapse = "")
  data <- paste0("<Data name=\"", names(fields), "\"><value>", xml_escape(fields), "</value></Data>", collapse = "")
  paste0("<Placemark><name>Explore ", a$area_rank_on_reef, "</name><styleUrl>#exploration-area</styleUrl>",
         "<description><![CDATA[<table>", rows, "</table><p>High-ranked habitat outside designated cull sites and away from ",
         "manta tows of the last 5 seasons. A place to look, not a prediction to trust.</p>]]></description>",
         "<ExtendedData>", data, "</ExtendedData><MultiGeometry>", polygon_kml(sf::st_geometry(areas_ll)[[j]]),
         "</MultiGeometry></Placemark>")
}
folders <- vapply(split(seq_len(nrow(area_table)), area_table$reef_name), function(js) {
  paste0("<Folder><name>Reef: ", xml_escape(area_table$reef_name[js[1]]), "</name><visibility>1</visibility>",
         paste0(vapply(js, placemark, character(1)), collapse = ""), "</Folder>")
}, character(1))
style <- paste0("<Style id=\"exploration-area\"><LineStyle><color>ffffffff</color><width>1.5</width></LineStyle>",
                "<PolyStyle><color>b3a73a4a</color><fill>1</fill><outline>1</outline></PolyStyle>",
                "<LabelStyle><scale>0.8</scale></LabelStyle></Style>")
exploration_folder <- paste0(
  "<Folder><name>Exploration areas ", sr$target_season, " (purple): outside cull sites, away from recent manta tows</name>",
  "<open>0</open><visibility>1</visibility>", paste0(folders, collapse = ""), "</Folder>")

stage_kml <- file.path(cfg$dirs$cache, "site_ranking_kmz", "doc.kml")
if (!file.exists(stage_kml)) stop("Run 20_site_ranking_2027.R first (needs ", stage_kml, ")")
kml <- paste(readLines(stage_kml, encoding = "UTF-8", warn = FALSE), collapse = "\n")
kml <- sub("<Style id=", paste0(style, "\n    <Style id="), kml, fixed = TRUE)
kml <- sub("</Document>", paste0(exploration_folder, "\n  </Document>"), kml, fixed = TRUE)
combined_dir <- file.path(cfg$dirs$cache, "site_ranking_kmz_exploration")
dir.create(combined_dir, showWarnings = FALSE)
writeLines(kml, file.path(combined_dir, "doc.kml"), useBytes = TRUE)
kmz_path <- file.path(out_dir, "cots_cull_site_priority_2027.kmz")
if (file.exists(kmz_path)) file.remove(kmz_path)
utils::zip(kmz_path, file.path(combined_dir, "doc.kml"), flags = "-j9Xq")
log_message("Wrote ", kmz_path, " (cull-site priority + exploration areas)")

print(summary(area_table[, c("area_ha", "habitat_percentile", "nearest_cull_site_m", "nearest_recent_tow_m")]))
