# Manuscript Figure 1 (Frontiers): the operational habitat layer (RG + year
# with the clean + year fallback) across the GBR with a 30 m inset of one reef, beside the 2027 workflow on four reefs
# (the report's fig22 without satellite imagery, because the article is CC-BY
# and Esri imagery may not be relicensed). Reef and land outlines come from the
# GBR features layer. Built by 07_build_report_figures.R
# (--manuscript_map_only=1 rebuilds only this figure).

# Background imagery: one Sentinel-2 L2A true-colour scene per reef, from the
# Earth Search catalogue (Element 84, public COGs on AWS). Copernicus Sentinel
# data are free and open with attribution, so they can be published under
# CC BY. Mosaics such as EOX Sentinel-2 cloudless mask the open ocean and show
# no offshore reefs. The first candidate (lowest scene cloud cover) whose reef
# window is complete and nearly free of bright cloud pixels is cached in
# cache/s2_<reef>.tif; the scenes used are written to
# manuscript_figure1_imagery.csv.
s2_background <- function(cfg, slug, frame, map_crs) {
  out <- file.path(cfg$dirs$cache, paste0("s2_", slug, ".tif"))
  if (file.exists(out)) return(terra::rast(out))
  window <- sf::st_buffer(sf::st_as_sfc(frame), 200)
  body <- list(collections = list("sentinel-2-l2a"), bbox = as.numeric(sf::st_bbox(sf::st_transform(window, 4326))),
               limit = 15, query = list(`eo:cloud_cover` = list(lt = 5)),
               sortby = list(list(field = "properties.eo:cloud_cover", direction = "asc")))
  res <- httr::POST("https://earth-search.aws.element84.com/v1/search", body = body, encode = "json")
  items <- jsonlite::fromJSON(httr::content(res, "text", encoding = "UTF-8"), simplifyVector = FALSE)$features
  for (f in items) {
    r <- terra::rast(paste0("/vsicurl/", f$assets$visual$href))
    x <- tryCatch(terra::crop(r, terra::ext(terra::vect(sf::st_transform(window, terra::crs(r))))), error = function(e) NULL)
    if (is.null(x)) next
    v <- terra::values(x)
    valid <- stats::complete.cases(v) & rowSums(v, na.rm = TRUE) > 0
    if (mean(valid) < 0.99 || mean(rowSums(v[valid, , drop = FALSE] > 180) == 3) > 0.02) next
    y <- terra::project(x, map_crs, method = "bilinear", res = 10)
    names(y) <- c("red", "green", "blue")
    terra::writeRaster(y, out, datatype = "INT1U", overwrite = TRUE)
    utils::write.csv(data.frame(reef = slug, scene = f$id, date = substr(f$properties$datetime, 1, 10)),
                     sub("\\.tif$", ".csv", out), row.names = FALSE)
    return(terra::rast(out))
  }
  stop("No clear Sentinel-2 scene for ", slug)
}

build_manuscript_map_figure <- function(cfg) {
  suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
  })
  ranking_dir <- file.path(cfg$dirs$outputs, "site_ranking")
  # The operational habitat layer: RG + year, with the calibrated clean + year
  # fallback where ReefGuidance layers are missing (script 05).
  overview_raster <- terra::rast(file.path(cfg$dirs$cache, "search_planning_120m.tif"))
  native <- terra::rast(file.path(cfg$dirs$layers, "cots_search_planning_2026_30m.tif"))[["within_reef_percentile"]]
  map_crs <- terra::crs(overview_raster)
  features <- sf::st_read(checked_input(cfg, "reef_features"), quiet = TRUE) |>
    sf::st_transform(map_crs)
  land <- features[features$FEAT_NAME %in% c("Mainland", "Island"), ]
  land_coarse <- sf::st_simplify(land, dTolerance = 300)
  reefs <- features[!features$FEAT_NAME %in% c("Mainland", "Island"), ]

  reef_ids <- c(batt = "16-029", eyrie = "14-118", broadhurst = "18-100a", lady_musgrave = "23-082a")
  inset_reef <- "eyrie"
  reef_frame <- function(slug, polygon) {
    window <- cfg$figures$four_reef_windows[[slug]]
    if (is.null(window)) return(sf::st_bbox(sf::st_buffer(polygon, 300)))
    sf::st_bbox(window, crs = sf::st_crs(polygon))
  }
  geo <- lapply(stats::setNames(nm = names(reef_ids)), function(slug) {
    polygon <- features[features$LABEL_ID == reef_ids[[slug]], ]
    frame <- reef_frame(slug, polygon)
    list(name = polygon$GBR_NAME[1], polygon = polygon, frame = frame, image = s2_background(cfg, slug, frame, map_crs))
  })
  scenes <- do.call(rbind, lapply(names(reef_ids), function(slug)
    utils::read.csv(file.path(cfg$dirs$cache, paste0("s2_", slug, ".csv")), stringsAsFactors = FALSE)))
  write_table(scenes, cfg, "manuscript_figure1_imagery.csv")
  short_name <- function(slug) sub(" \\(.*", "", sub("^\\((Big|Little)\\) ", "\\1 ", geo[[slug]]$name))
  frame_sf <- function(frame) sf::st_as_sfc(frame)
  in_frame <- function(x, frame) x[lengths(sf::st_intersects(x, frame_sf(frame))) > 0, ]
  text_pt <- 8

  # --- GBR-wide clean + year, ~2.9 km cell means ----------------------------------
  overview <- terra::aggregate(overview_raster, fact = 24, fun = "mean", na.rm = TRUE)
  bbox <- sf::st_bbox(terra::ext(terra::trim(overview)) |> terra::vect(crs = map_crs) |> sf::st_as_sf()) +
    c(-20000, -30000, 20000, 20000)
  d <- raster_df(overview)
  limits <- stats::quantile(d$value, c(0.02, 0.98))
  markers <- dplyr::bind_rows(lapply(names(reef_ids), function(slug) {
    xy <- sf::st_coordinates(sf::st_centroid(sf::st_union(geo[[slug]]$polygon)))
    data.frame(slug = slug, label = short_name(slug), x = xy[1], y = xy[2])
  }))
  gbr <- ggplot() +
    geom_sf(data = land_coarse, fill = ink$land, colour = NA) +
    geom_raster(data = d, aes(x, y, fill = value)) +
    geom_point(data = markers, aes(x, y), shape = 22, size = 2.2, stroke = 0.6, colour = ink$primary, fill = NA) +
    geom_text(data = markers, aes(x, y, label = label), hjust = 1, nudge_x = -25000, size = text_pt / ggplot2::.pt,
              colour = ink$primary, family = report_font) +
    scale_fill_viridis_c(option = "viridis", limits = limits, oob = scales::squish,
                         name = "RG + year (clean + year fallback)\nP(problematic site-year), 2026",
                         labels = scales::label_number(accuracy = 0.01)) +
    coord_sf(xlim = bbox[c(1, 3)], ylim = bbox[c(2, 4)], expand = FALSE, datum = NA) +
    theme_map(text_pt) +
    theme(legend.position = "inside", legend.position.inside = c(0.03, 0.03), legend.justification = c(0, 0),
          legend.direction = "horizontal", legend.title.position = "top",
          legend.key.width = grid::unit(0.7, "cm"), legend.key.height = grid::unit(0.18, "cm"),
          legend.background = element_rect(fill = scales::alpha(ink$surface, 0.85), colour = NA),
          legend.title = element_text(size = text_pt), legend.text = element_text(size = text_pt - 0.5))

  # --- inset: the same layer at 30 m on one reef, as its deployed within-reef percentile
  g <- geo[[inset_reef]]
  v <- terra::vect(g$polygon)
  cells <- raster_df(terra::mask(terra::crop(native, v, snap = "out"), v))
  inset <- ggplot() +
    geom_sf(data = in_frame(land, g$frame), fill = ink$land, colour = NA) +
    geom_raster(data = cells, aes(x, y, fill = value)) +
    geom_sf(data = g$polygon, fill = NA, colour = ink$secondary, linewidth = 0.25) +
    scale_fill_viridis_c(option = "viridis", limits = c(0, 1), guide = "none") +
    coord_sf(xlim = g$frame[c(1, 3)], ylim = g$frame[c(2, 4)], expand = FALSE, datum = NA) +
    labs(title = short_name(inset_reef), subtitle = "30 m, rank on reef") +
    theme_map(text_pt) +
    theme(plot.title = element_text(size = text_pt, face = "bold"),
          plot.subtitle = element_text(size = text_pt - 0.5, colour = ink$secondary),
          plot.background = element_rect(fill = ink$surface, colour = ink$secondary, linewidth = 0.4),
          plot.margin = margin(3, 3, 3, 3))
  left <- gbr + inset_element(inset, left = 0.56, bottom = 0.56, right = 1, top = 1, align_to = "panel")

  # --- four reefs: 2027 density-aware v2 classes, tow zone and exploration zones -----
  explore <- sf::st_read(file.path(ranking_dir, "cots_exploration_zones_2027.gpkg"), quiet = TRUE)
  zone_levels <- unique(explore$zone[order(explore$habitat_floor)])
  explore$zone <- factor(explore$zone, zone_levels)
  zone_names <- stats::setNames(sub(" habitat$", "", sub(" %", "%", zone_levels)), zone_levels)
  zone_alpha <- 0.28
  cull_sites <- sf::st_read(file.path(cfg$dirs$cache, "cull_sites_designated.gpkg"), quiet = TRUE)
  cull_sites <- cull_sites[!duplicated(paste(cull_sites$reef_name, cull_sites$site_name)), ]
  density_sites <- utils::read.csv(file.path(ranking_dir, "cots_cull_site_priority_density_v2_2027.csv"), stringsAsFactors = FALSE)
  density_reefs <- utils::read.csv(file.path(ranking_dir, "cots_cull_site_priority_density_v2_2027_reefs.csv"), stringsAsFactors = FALSE)
  cull_sites$priority_density <- density_sites$priority_density_v2[match(paste(cull_sites$reef_name, cull_sites$site_name),
                                                                         paste(density_sites$reef_name, density_sites$site_name))]
  strong_classes <- c("direct_cots", "multiple_counted_scars", "legacy_common_scars")
  tow_geom <- readRDS(checked_input(cfg, "footprint_links"))$manta |>
    dplyr::group_by(.data$visit_id) |>
    dplyr::summarise(sx = dplyr::first(.data$start_x), sy = dplyr::first(.data$start_y),
                     ex = dplyr::first(.data$end_x), ey = dplyr::first(.data$end_y), .groups = "drop")
  tows <- readRDS(checked_input(cfg, "manta_tows")) |>
    dplyr::select("visit_id", "survey_date", "manta_evidence_class") |>
    dplyr::inner_join(tow_geom, by = "visit_id") |>
    dplyr::filter(is.finite(.data$sx), .data$survey_date > cfg$site_ranking$evidence_end - 365,
                  .data$survey_date <= cfg$site_ranking$evidence_end) |>
    dplyr::mutate(evidence = ifelse(.data$manta_evidence_class %in% strong_classes, "Strong evidence", "Weak or none"))
  tows <- tows[order(tows$evidence == "Strong evidence"), ]
  priority_cols <- c(A = "#e34948", B = "#eb6834", C = "#eda100", D = "#f2d98d", E = "#898781")
  zone_label <- "Follow-up tow zone: A–C sites + 100 m"
  reef_panel <- function(slug) {
    g <- geo[[slug]]
    frame <- g$frame
    label <- reef_ids[[slug]]
    in_reef <- function(x) grepl(paste0("(", label, ")"), x, fixed = TRUE)
    sites_r <- sf::st_transform(cull_sites[in_reef(cull_sites$reef_name), ], sf::st_crs(g$polygon))
    zone_c <- sf::st_sf(geometry = sf::st_union(sf::st_buffer(sites_r[sites_r$priority_density %in% c("A", "B", "C"), ], 100)))
    explore_r <- sf::st_transform(explore[in_reef(explore$reef_name), ], sf::st_crs(g$polygon))
    tows_r <- tows[tows$sx > frame[1] & tows$sx < frame[3] & tows$sy > frame[2] & tows$sy < frame[4], ]
    rd <- density_reefs[in_reef(density_reefs$reef_name), ]
    ggplot() +
      tidyterra::geom_spatraster_rgb(data = g$image, maxcell = 2e6) +
      geom_sf(data = g$polygon, fill = NA, colour = "white", linewidth = 0.25) +
      geom_sf(data = explore_r[order(explore_r$habitat_floor), ], aes(alpha = zone), fill = "#7b68d9", colour = NA) +
      geom_sf(data = sites_r, aes(fill = priority_density), colour = NA, alpha = 0.85) +
      geom_segment(data = tows_r, aes(x = sx, y = sy, xend = ex, yend = ey, colour = evidence, linewidth = evidence)) +
      geom_sf(data = zone_c, aes(linetype = zone_label), fill = NA, colour = "white", linewidth = 0.6) +
      scale_fill_manual(values = priority_cols, limits = names(priority_cols), name = "Site class (density-aware)",
                        na.translate = FALSE) +
      scale_colour_manual(values = c("Strong evidence" = "#2ee6ff", "Weak or none" = "#cde2fb"),
                          limits = c("Strong evidence", "Weak or none"), name = "Manta tows, past 12 months") +
      scale_linewidth_manual(values = c("Strong evidence" = 0.8, "Weak or none" = 0.25),
                             limits = c("Strong evidence", "Weak or none"), guide = "none") +
      scale_linetype_manual(values = stats::setNames("22", zone_label), limits = zone_label, name = NULL) +
      scale_alpha_manual(values = stats::setNames(rep(zone_alpha, length(zone_levels)), zone_levels),
                         labels = zone_names, name = "Exploration zones", drop = FALSE) +
      coord_sf(xlim = frame[c(1, 3)], ylim = frame[c(2, 4)], expand = FALSE, datum = NA) +
      labs(title = short_name(slug),
           subtitle = sprintf("Expected problematic %.0f%%\nA–C %d, A–D %d of %d sites", 100 * rd$reef_density_v2,
                              rd$priority_a + rd$priority_b + rd$priority_c, rd$sites - rd$priority_e, rd$sites)) +
      theme_map(text_pt) +
      theme(plot.title = element_text(size = text_pt + 0.5, face = "bold"),
            plot.subtitle = element_text(size = text_pt - 0.5, colour = ink$secondary, lineheight = 0.95),
            plot.margin = margin(2, 3, 2, 3), legend.position = "none")
  }
  # Legends drawn once, from Eyrie (every key present), in two blocks.
  legend_grob <- function(p) {
    g <- ggplotGrob(p)
    i <- grep("^guide-box", g$layout$name)
    i <- i[!vapply(g$grobs[i], inherits, logical(1), "zeroGrob")]
    g$grobs[[i[1]]]
  }
  legend_theme <- theme(legend.position = "bottom", legend.box = "horizontal", legend.spacing.x = grid::unit(0.2, "cm"),
                        legend.direction = "horizontal", legend.title = element_text(size = text_pt),
                        legend.text = element_text(size = text_pt - 0.5), legend.spacing.y = grid::unit(0.05, "cm"),
                        legend.margin = margin(0, 0, 0, 0), legend.key.size = grid::unit(0.3, "cm"))
  eyrie <- reef_panel("eyrie")
  legend_a <- legend_grob(eyrie + legend_theme +
    guides(fill = guide_legend(order = 1, nrow = 1, override.aes = list(alpha = 0.9)),
           linetype = guide_legend(order = 2, override.aes = list(colour = ink$secondary)),
           colour = "none", alpha = "none"))
  legend_b <- legend_grob(eyrie + legend_theme +
    guides(fill = "none", linetype = "none",
           colour = guide_legend(order = 1, override.aes = list(linewidth = c(0.8, 0.4))),
           alpha = guide_legend(order = 2, nrow = 1, override.aes = list(
             fill = "#7b68d9", alpha = 1 - (1 - zone_alpha)^seq_along(zone_levels)))))
  # Columns pair a wide reef with a tall one: Batt over Lady Musgrave, Eyrie
  # over Big Broadhurst.
  right <- (reef_panel("batt") | eyrie) / (reef_panel("lady_musgrave") | reef_panel("broadhurst")) +
    plot_layout(widths = c(1.35, 0.85))
  legends <- wrap_elements(full = legend_a) / wrap_elements(full = legend_b)
  fig <- (wrap_plots(left, right, widths = c(1, 1.45)) / legends) + plot_layout(heights = c(1, 0.12), guides = "keep")
  save_figure(fig, cfg, "fig35_manuscript_gbr_and_workflow.png", width = 7.1, height = 6.1)
}
