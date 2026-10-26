# Portrait report versions of concentration figures lc2, lc4 and lc6.
# Use completed concentration outputs only; do not imply source attribution.
build_larval_concentration_figures <- function(cfg) {
  read_lc <- function(n) read.csv(file.path(cfg$dirs$tables, paste0("larval_concentration_", n, ".csv")))
  z <- read_lc("reef_zones")
  s <- read_lc("reef_summary")
  r <- read_lc("releases")
  geo <- readRDS(file.path(cfg$dirs$cache, "larval_concentration_report.rds"))
  zoom <- c(xmin = 145.28, xmax = 145.74, ymin = -14.95, ymax = -14.45)
  box <- sf::st_as_sfc(sf::st_bbox(zoom, crs = 4326))
  # A small set of named reefs keeps labels legible at final report width.
  labels <- z |> dplyr::mutate(label = sub(" Reef.*$| Reefs$| \\(.*$", "", .data$name)) |>
    dplyr::filter(!grepl("^U/N", .data$name)) |>
    dplyr::group_by(.data$label) |> dplyr::slice_max(.data$area, n = 1, with_ties = FALSE) |> dplyr::ungroup()
  panel <- function(xlim, ylim, selected, domain = FALSE) {
    p <- ggplot2::ggplot() + ggplot2::geom_sf(data = geo$land, fill = ink$land, colour = NA) +
      ggplot2::geom_sf(data = geo$sites, ggplot2::aes(fill = .data$top20_share, colour = .data$top20_share), linewidth = .15) +
      ggplot2::scale_fill_gradientn(colours = ramp_orange, limits = c(0, 1), breaks = c(0, .5, 1),
        labels = scales::label_percent(), name = "Releases in\ntop fifth", aesthetics = c("fill", "colour")) +
      ggrepel::geom_text_repel(data = labels |> dplyr::filter(.data$label %in% selected),
        ggplot2::aes(.data$lon, .data$lat, label = .data$label), size = 2.85, family = report_font,
        seed = 23, min.segment.length = 0, max.overlaps = Inf, box.padding = .3) +
      theme_map(9)
    if (domain) p <- p + ggplot2::geom_sf(data = box, fill = NA, colour = ink$secondary, linetype = 2, linewidth = .35)
    p + ggplot2::coord_sf(xlim = xlim, ylim = ylim, expand = FALSE, datum = NA)
  }
  a <- panel(c(145.1, 145.9), c(-15.8, -14.27), c("Lizard Island", "Jewell", "Mackay", "Williamson", "Osterlund"), TRUE) +
    ggplot2::labs(title = "A  Model domain")
  b <- panel(zoom[c("xmin", "xmax")], zoom[c("ymin", "ymax")],
             c("Lizard Island", "MacGillivray", "Eyrie", "Martin", "North Direction Island", "Ribbon No 10", "Yonge", "Carter")) +
    ggplot2::labs(title = "B  Lizard Island group")
  main <- (a | b) + patchwork::plot_layout(widths = c(1, 1.3), guides = "collect") &
    ggplot2::theme(legend.position = "bottom", legend.direction = "horizontal", legend.text = ggplot2::element_text(size = 8),
      legend.title = ggplot2::element_text(size = 8), legend.key.width = grid::unit(1.2, "cm"))
  save_figure(main, cfg, "fig28_larval_hotspot_persistence.png", width = 6.26, height = 6.6)

  zones <- c("Inner shelf", "Lizard group", "Other mid-shelf", "Ribbon Reefs", "Outer barrier N of Lizard")
  cols <- setNames(c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"), zones)
  r$release <- as.Date(r$release)
  r$zone <- factor(r$zone, rev(zones))
  r$age <- factor(paste("Day", r$age), paste("Day", c(5, 9, 10)))
  r$display <- ifelse(r$rel_conc < .1, "<0.1", sprintf("%.1f", r$rel_conc))
  heat <- ggplot2::ggplot(r, ggplot2::aes(.data$release, .data$zone, fill = .data$log2_rel)) +
    ggplot2::geom_tile(colour = ink$surface, linewidth = .6) +
    ggplot2::geom_text(ggplot2::aes(label = .data$display, colour = abs(.data$log2_rel) > 1.6), size = 2.85) +
    ggplot2::scale_colour_manual(values = c("FALSE" = ink$primary, "TRUE" = "white"), guide = "none") +
    ggplot2::scale_fill_gradientn(colours = rev(diverging_blue_red), limits = c(-3, 3), oob = scales::squish,
      breaks = c(-3, 0, 3), labels = c("1/8", "1", "8"), name = "Zone mean / domain mean") +
    ggplot2::scale_x_date(breaks = sort(unique(r$release)), date_labels = "%d\n%b") +
    ggplot2::scale_y_discrete(labels = function(x) gsub("Outer barrier N of Lizard", "Outer barrier\nnorth of Lizard", x)) +
    ggplot2::facet_wrap(~age, ncol = 1) + ggplot2::labs(x = "Release date", y = NULL) + theme_report(10) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "bottom", axis.text = ggplot2::element_text(size = 8))
  save_figure(heat, cfg, "fig29_larval_release_comparison.png", width = 6.26, height = 6.8)

  bands <- c("South: Cooktown to Ribbon 4", "Central: Ribbon 5 to 9", "Lizard: Ribbon 10 to Lizard", "North of Lizard")
  s <- s |> dplyr::left_join(z |> dplyr::select(.data$GBRMPA_ID, .data$dist_coast_km), by = "GBRMPA_ID") |>
    dplyr::mutate(band = factor(cut(.data$lat, c(-15.8, -15.35, -14.9, -14.6, -14.25), labels = bands), rev(bands)),
                  age = factor(paste("Day", .data$age), paste("Day", c(5, 9, 10))))
  cross <- ggplot2::ggplot(s, ggplot2::aes(.data$dist_coast_km, .data$mean_pct, colour = .data$zone)) +
    ggplot2::geom_hline(yintercept = .5, colour = ink$axis, linewidth = .3) + ggplot2::geom_point(size = 1.4) +
    ggplot2::scale_colour_manual(values = cols, name = NULL) +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = c(0, .5, 1), labels = c("0", "0.5", "1")) +
    ggplot2::scale_x_continuous(breaks = c(10, 30, 50)) +
    ggplot2::facet_grid(band ~ age, labeller = ggplot2::labeller(band = ggplot2::label_wrap_gen(13))) +
    ggplot2::labs(x = "Distance from mainland (km)", y = "Mean within-release concentration rank") +
    theme_report(10) + ggplot2::theme(legend.position = "bottom", legend.text = ggplot2::element_text(size = 8),
      axis.text = ggplot2::element_text(size = 8), strip.text.y = ggplot2::element_text(angle = 0, size = 8)) +
    ggplot2::guides(colour = ggplot2::guide_legend(ncol = 2))
  save_figure(cross, cfg, "fig30_larval_cross_shelf.png", width = 6.26, height = 7.0)
}
