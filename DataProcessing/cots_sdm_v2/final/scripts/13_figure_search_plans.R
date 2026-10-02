# Illustrative search plans on Batt and Eyrie reefs from the primary
# search-planning layer (11), over satellite imagery.
#
#   Manta panel: the largest recent survey with strong COTS evidence; tows the
#     layer would keep (top 60 % within the survey) vs skip, strong tows marked.
#   Dive panel: the same number of dive sites as that reef's most recent full
#     cull year, placed greedily on the highest-ranked diveable cells (DEM
#     -16 to -0.5 m, >= 400 m apart, at most 15) and ordered into a short route
#     (nearest neighbour + 2-opt), with that year's actual sites by outcome.
# These reefs are in the layer's training data, so the panels illustrate the
# planning logic; the held-out savings estimates come from 12.

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "..", "R", "common.R"))
cfg <- initialise_final()
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})
reefs_to_plot <- c(batt = "16-029", eyrie = "14-118")
keep_fraction <- cfg$operations$headline_keep
dive_depth <- c(-16, -0.5)
min_spacing_m <- 400
max_planned_sites <- 15
strong_classes <- c("direct_cots", "multiple_counted_scars", "legacy_common_scars")

layer <- terra::rast(file.path(cfg$dirs$layers, "cots_search_planning_2026_30m.tif"))[["within_reef_percentile"]]
dem <- terra::rast(checked_input(cfg, "predictor_stack_clean"))[["DEM"]]
template <- terra::rast(file.path(cfg$dirs$cache, "base_support_120m.tif"))[[1]]
search120 <- terra::rast(file.path(cfg$dirs$cache, "search_planning_120m.tif"))
features <- sf::st_read(checked_input(cfg, "reef_features"), quiet = TRUE) |> sf::st_transform(terra::crs(layer))
tows_all <- readRDS(checked_input(cfg, "manta_tows"))
links_all <- readRDS(checked_input(cfg, "footprint_links"))$manta
dives_all <- readRDS(file.path(cfg$exploration_dir, "artifacts", "cull_observations_clean.rds"))
tile_dir <- file.path(cfg$dirs$cache, "tiles")
to_mga <- function(lon, lat) sf::st_coordinates(sf::st_transform(sf::st_as_sf(data.frame(lon, lat), coords = c("lon", "lat"), crs = 4326), terra::crs(layer)))

route_order <- function(xy) {
  n <- nrow(xy)
  if (n < 3) return(seq_len(n))
  d <- as.matrix(stats::dist(xy))
  tour <- which.min(xy[, 1]) # start from the landward (western) site
  while (length(tour) < n) {
    rest <- setdiff(seq_len(n), tour)
    tour <- c(tour, rest[which.min(d[tail(tour, 1), rest])])
  }
  improved <- TRUE
  while (improved) {
    improved <- FALSE
    for (i in 2:(n - 2)) for (j in (i + 1):(n - 1)) {
      a <- tour[i - 1]; b <- tour[i]; c <- tour[j]; e <- tour[j + 1]
      if (d[a, c] + d[b, e] < d[a, b] + d[c, e] - 1e-9) {
        tour[i:j] <- rev(tour[i:j]); improved <- TRUE
      }
    }
  }
  tour
}

plan_reef <- function(label) {
  polygon <- features[features$LABEL_ID == label & features$FEAT_NAME == "Reef", ]
  name <- polygon$GBR_NAME[1]
  reef_name <- unique(tows_all$reef_name[grepl(label, tows_all$reef_name, fixed = TRUE)])[1]
  frame <- sf::st_bbox(sf::st_buffer(polygon, 400))
  v <- terra::vect(polygon)
  pct <- terra::mask(terra::crop(layer, v, snap = "out"), v)
  pct_df <- raster_df(pct)

  # Manta survey: largest since 2023 with strong evidence.
  surveys <- tows_all |>
    dplyr::filter(.data$reef_name == !!reef_name, .data$year >= 2023) |>
    dplyr::group_by(.data$survey_date) |>
    dplyr::summarise(tows = dplyr::n(), strong = sum(.data$manta_evidence_class %in% strong_classes), .groups = "drop") |>
    dplyr::filter(.data$strong > 0) |>
    dplyr::arrange(dplyr::desc(.data$tows))
  survey_date <- surveys$survey_date[1]
  tows <- tows_all |> dplyr::filter(.data$reef_name == !!reef_name, .data$survey_date == !!survey_date)
  lk <- links_all |> dplyr::filter(.data$visit_id %in% tows$visit_id)
  lk$value <- terra::extract(search120, terra::cellFromXY(template, cbind((lk$x_bin + 0.5) * 120, (lk$y_bin + 0.5) * 120)))[[1]]
  geometry <- lk |>
    dplyr::group_by(.data$visit_id) |>
    dplyr::summarise(score = stats::weighted.mean(.data$value, pmax(.data$weighted_length_m, 1), na.rm = TRUE),
                     start_lon = dplyr::first(.data$start_longitude), start_lat = dplyr::first(.data$start_latitude),
                     end_lon = dplyr::first(.data$end_longitude), end_lat = dplyr::first(.data$end_latitude), .groups = "drop")
  tows <- dplyr::inner_join(tows, geometry, by = "visit_id") |>
    dplyr::mutate(strong = .data$manta_evidence_class %in% strong_classes,
                  rank = rank(-.data$score, ties.method = "first"),
                  kept = .data$rank <= ceiling(keep_fraction * dplyr::n()))
  s <- to_mga(tows$start_lon, tows$start_lat); e <- to_mga(tows$end_lon, tows$end_lat)
  tows$x <- s[, 1]; tows$y <- s[, 2]; tows$xend <- e[, 1]; tows$yend <- e[, 2]

  # Cull year: most recent complete year with >= 10 sites.
  dives <- dives_all |> dplyr::filter(.data$reef_name == !!reef_name)
  years <- dives |> dplyr::group_by(.data$year) |> dplyr::summarise(sites = dplyr::n_distinct(.data$site_id), .groups = "drop") |>
    dplyr::filter(.data$sites >= 10, .data$year < max(dives_all$year))
  cull_year <- max(years$year)
  sites <- dives |>
    dplyr::filter(.data$year == cull_year) |>
    dplyr::group_by(.data$site_id) |>
    dplyr::summarise(lon = dplyr::first(.data$longitude), lat = dplyr::first(.data$latitude), dives = dplyr::n(),
                     ghost_dives = sum(.data$total_cots == 0), cots = sum(.data$total_cots),
                     cpue = sum(.data$total_cots) / sum(.data$bottom_time), .groups = "drop") |>
    dplyr::mutate(outcome = dplyr::case_when(.data$cots == 0 ~ "0 COTS (ghost)",
                                             .data$cpue < cfg$legacy$cpue_threshold ~ "< 0.02 CPUE",
                                             TRUE ~ "≥ 0.02 CPUE"))
  sxy <- to_mga(sites$lon, sites$lat); sites$x <- sxy[, 1]; sites$y <- sxy[, 2]
  sites$percentile <- terra::extract(layer, cbind(sites$x, sites$y))[[1]]

  # Planned dive sites on top-ranked diveable cells.
  depth <- terra::mask(terra::crop(dem, v, snap = "out"), v)
  cand <- data.frame(terra::xyFromCell(pct, seq_len(terra::ncell(pct))),
                     p = terra::values(pct, mat = FALSE), z = terra::values(depth, mat = FALSE))
  cand <- cand[is.finite(cand$p) & is.finite(cand$z) & cand$z >= dive_depth[1] & cand$z <= dive_depth[2], ]
  cand <- cand[order(-cand$p), ]
  n_sites <- min(nrow(sites), max_planned_sites)
  picked <- integer(0)
  for (i in seq_len(nrow(cand))) {
    if (!length(picked) || all((cand$x[picked] - cand$x[i])^2 + (cand$y[picked] - cand$y[i])^2 >= min_spacing_m^2)) picked <- c(picked, i)
    if (length(picked) == n_sites) break
  }
  plan <- cand[picked, ]
  plan <- plan[route_order(as.matrix(plan[, c("x", "y")])), ]
  plan$stop <- seq_len(nrow(plan))
  route_km <- sum(sqrt(diff(plan$x)^2 + diff(plan$y)^2)) / 1000

  imagery <- maptiles::get_tiles(sf::st_as_sfc(frame), provider = "Esri.WorldImagery", zoom = 14,
                                 crop = TRUE, project = TRUE, cachedir = tile_dir)
  summary <- tibble::tibble(
    reef = name, label = label, survey_date = as.character(survey_date), n_tows = nrow(tows),
    tows_kept = sum(tows$kept), strong_tows = sum(tows$strong), strong_kept = sum(tows$strong & tows$kept),
    cull_year = cull_year, actual_sites = nrow(sites), actual_dives = sum(sites$dives),
    actual_ghost_dives = sum(sites$ghost_dives),
    actual_sites_in_top30 = sum(sites$percentile >= 0.7, na.rm = TRUE),
    problematic_sites = sum(sites$outcome == "≥ 0.02 CPUE"),
    problematic_sites_in_top30 = sum(sites$outcome == "≥ 0.02 CPUE" & sites$percentile >= 0.7, na.rm = TRUE),
    ghost_only_sites = sum(sites$outcome == "0 COTS (ghost)"),
    ghost_only_sites_in_top30 = sum(sites$outcome == "0 COTS (ghost)" & sites$percentile >= 0.7, na.rm = TRUE),
    planned_sites = nrow(plan), planned_route_km = route_km
  )
  list(name = name, label = label, polygon = polygon, frame = frame, pct = pct_df, imagery = imagery,
       tows = tows, sites = sites, plan = plan, summary = summary)
}
plans <- lapply(reefs_to_plot, plan_reef)
write_table(dplyr::bind_rows(lapply(plans, `[[`, "summary")), cfg, "search_plan_examples.csv")

# --- Figure ---------------------------------------------------------------------
base_map <- function(p, title) {
  d <- p$pct
  d$alpha <- ifelse(d$value >= 0.7, 0.75, 0.25)
  ggplot() +
    tidyterra::geom_spatraster_rgb(data = p$imagery, maxcell = 2e6) +
    geom_raster(data = d, aes(x, y, fill = value, alpha = alpha)) +
    geom_sf(data = p$polygon, fill = NA, colour = "white", linewidth = 0.25) +
    scale_fill_viridis_c(option = "viridis", limits = c(0.7, 1), oob = scales::squish,
                         breaks = c(0.7, 0.8, 0.9, 1), labels = c("≤ 70th", "80th", "90th", "100th"),
                         name = "Within-reef percentile") +
    scale_alpha_identity() +
    labs(title = title) +
    theme_map(8)
}
coords <- function(p) coord_sf(xlim = p$frame[c(1, 3)], ylim = p$frame[c(2, 4)], expand = FALSE, datum = NA)
manta_panel <- function(p) {
  s <- p$summary
  base_map(p, sprintf("%s: manta survey %s", sub(" \\(.*", "", p$name), format(as.Date(s$survey_date), "%d %b %Y"))) +
    geom_segment(data = p$tows, aes(x, y, xend = xend, yend = yend, linetype = ifelse(kept, "Kept (top 60%)", "Skipped")),
                 colour = "white", linewidth = 0.55) +
    geom_point(data = dplyr::filter(p$tows, strong), aes((x + xend) / 2, (y + yend) / 2, shape = "Strong COTS evidence"),
               colour = "#e34948", fill = "#e34948", size = 1.1) +
    scale_linetype_manual(values = c(`Kept (top 60%)` = "solid", Skipped = "22"), name = "Manta tows",
                          guide = guide_legend(override.aes = list(colour = ink$primary))) +
    scale_shape_manual(values = c(`Strong COTS evidence` = 21), name = "Manta evidence") +
    labs(subtitle = sprintf("Keeping %d of %d tows retains
%d of %d strong-evidence tows",
                            s$tows_kept, s$n_tows, s$strong_kept, s$strong_tows)) +
    coords(p)
}
cull_panel <- function(p) {
  s <- p$summary
  base_map(p, sprintf("%s: dive plan vs %d cull sites", sub(" \\(.*", "", p$name), s$cull_year)) +
    geom_point(data = p$sites, aes(x, y, colour = outcome), size = 1.4, shape = 16) +
    geom_path(data = p$plan, aes(x, y), colour = "white", linewidth = 0.4, linetype = "solid") +
    geom_point(data = p$plan, aes(x, y), shape = 21, fill = "white", colour = ink$primary, size = 2.5, stroke = 0.3) +
    geom_text(data = p$plan, aes(x, y, label = stop), size = 1.5, colour = ink$primary, family = report_font) +
    scale_colour_manual(values = c(`0 COTS (ghost)` = "#c3c2b7", `< 0.02 CPUE` = "#eda100",
                                   `≥ 0.02 CPUE` = "#e34948"),
                        name = sprintf("Actual cull sites")) +
    labs(subtitle = sprintf("%d planned sites, %.1f km of straight legs. In the top 30%%:
%d of %d problematic sites, %d of %d ghost-only sites",
                            s$planned_sites, s$planned_route_km, s$problematic_sites_in_top30, s$problematic_sites,
                            s$ghost_only_sites_in_top30, s$ghost_only_sites)) +
    coords(p)
}
fig <- wrap_plots(
  manta_panel(plans$batt), cull_panel(plans$batt),
  manta_panel(plans$eyrie), cull_panel(plans$eyrie),
  ncol = 2, guides = "collect", heights = c(0.62, 1)
) +
  plot_annotation(
    title = "Planning a reef visit with the primary search-planning layer",
    subtitle = paste(
      "Colour shows each reef's top 30% (RG + year, clean fallback). Left: tows the layer would keep in a real survey.",
      "Right: numbered dive sites on the top-ranked diveable cells (16-0.5 m), ordered into a short route.",
      sep = "\n"
    ),
    caption = paste(
      "Illustrative: these reefs are in the layer's training data; held-out savings estimates are in fig11. Route legs are straight lines, not navigation paths.",
      paste("Imagery:", maptiles::get_credit("Esri.WorldImagery")),
      sep = "\n"
    ),
    theme = theme_report()
  ) &
  theme(legend.position = "bottom", legend.box = "horizontal", legend.title.position = "top",
        legend.key.width = grid::unit(0.7, "cm"), legend.key.height = grid::unit(0.18, "cm"),
        legend.text = element_text(size = 6.5), legend.title = element_text(size = 7),
        legend.spacing.x = grid::unit(0.15, "cm"), legend.box.spacing = grid::unit(0.1, "cm"),
        plot.subtitle = element_text(size = 6.5, colour = ink$secondary))
save_figure(fig, cfg, "fig12_search_plans_batt_eyrie.png", width = 7.5, height = 8.6)
print(dplyr::bind_rows(lapply(plans, `[[`, "summary")))
