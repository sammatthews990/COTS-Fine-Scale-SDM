# Conceptual diagram of the COTS-coral metapopulation model (report Methods,
# beside the mechanisms table).
#
# Question: what does the Lizard Island prototype represent, and which model
# component does what? Method: (a) the real sub-reef domain (2,895 sites on 113
# reefs, the same GeoPackage as fig29 panel d); (b) a schematic of the COTS life
# cycle and its coral feedbacks, following the mechanisms table and Appendix E;
# (c) the annual time-step order from the ADRIA COTS submodel documentation
# (docs/src/concepts/cots_submodel.md, "Simulation Timestep Order of Events").
# Validation: none; this is a description of model structure, not a result.
# Larval connections are drawn schematically because the population model's
# connectivity matrices are not pinned inputs of this pipeline.

concept_star <- function(x0, y0, r, n = 14, inner = 0.55, id = 1) {
  a <- seq(0, 2 * pi, length.out = 2 * n + 1)[-1] + pi / 2
  rr <- rep(c(r, r * inner), n)
  data.frame(x = x0 + rr * cos(a), y = y0 + rr * sin(a), id = id)
}

concept_ellipse <- function(x0, y0, w, h, from = 0, to = pi, n = 60, id = 1) {
  a <- seq(from, to, length.out = n)
  data.frame(x = x0 + w * cos(a), y = y0 + h * sin(a), id = id)
}

# Tabular and corymbose (fast, preferred) and massive (slow) coral glyphs.
concept_fast_coral <- function(x0, y0, s = 1) {
  plate <- rbind(concept_ellipse(x0 - 1.2 * s, y0 + 3.2 * s, 4.2 * s, 1.1 * s, 0, pi, id = "plate"),
                 concept_ellipse(x0 - 1.2 * s, y0 + 3.2 * s, 4.2 * s, 0.35 * s, pi, 2 * pi, id = "plate"))
  stalk <- data.frame(x = x0 - 1.2 * s + c(-0.45, 0.45, 0.7, -0.7) * s, y = y0 + c(3.0, 3.0, 0, 0) * s, id = "stalk")
  bush <- do.call(rbind, lapply(seq_len(7), function(i) {
    cx <- x0 + 4.2 * s + c(-1.4, -0.7, 0, 0.7, 1.4, -0.35, 0.35)[i] * s
    cy <- y0 + c(1.1, 1.6, 1.9, 1.6, 1.1, 2.5, 2.5)[i] * s
    concept_ellipse(cx, cy, 0.42 * s, 0.75 * s, 0, 2 * pi, n = 30, id = paste0("bush", i))
  }))
  rbind(plate, stalk, bush)
}

concept_slow_coral <- function(x0, y0, s = 1) {
  dome <- concept_ellipse(x0, y0, 3.6 * s, 2.8 * s, 0, pi, id = "dome")
  bumps <- do.call(rbind, lapply(seq_len(5), function(i) {
    a <- c(0.35, 0.95, 1.57, 2.2, 2.8)[i]
    concept_ellipse(x0 + 3.3 * s * cos(a), y0 + 2.55 * s * sin(a), 0.55 * s, 0.45 * s, 0, 2 * pi, n = 24, id = paste0("b", i))
  }))
  rbind(dome, bumps)
}

build_metapop_concept_figure <- function(cfg) {
  txt <- 2.45   # about 7 pt
  nm <- 2.75    # node names
  cots_ink <- "#3a3a3a"
  fast_fill <- "#e7a17a"
  slow_fill <- "#c9b48a"
  arrow_col <- "#4d4d4d"
  cull_col <- key_colours[["culled"]]
  arr <- ggplot2::arrow(length = grid::unit(1.6, "mm"), type = "closed")

  # (a) Spatial domain --------------------------------------------------------
  s2_old <- suppressMessages(sf::sf_use_s2(FALSE))
  on.exit(suppressMessages(sf::sf_use_s2(s2_old)), add = TRUE)
  path <-checked_input(cfg, "reef_features")
  layer <- sf::st_layers(path)$name[1]
  land <- sf::st_read(path, query = paste0('SELECT * FROM "', layer, '" WHERE FEAT_NAME IN (\'Mainland\',\'Island\')'), quiet = TRUE)
  land <- sf::st_make_valid(sf::st_transform(land, 4326))
  domain <- sf::st_read(checked_input(cfg, "report_metapop_domain"), quiet = TRUE)
  stopifnot(nrow(domain) == 2895, length(unique(domain$GBRMPA_ID)) == 113)
  domain <- sf::st_transform(domain, 4326)
  site_fill <- okabe_ito[["blue"]]
  land_dom <- suppressWarnings(sf::st_crop(land, c(xmin = 144.9, xmax = 146.1, ymin = -16, ymax = -14.1)))
  zoom <- c(xmin = 145.42, xmax = 145.49, ymin = -14.715, ymax = -14.643)
  dom_map <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = land_dom, fill = ink$land, colour = NA) +
    ggplot2::geom_sf(data = domain, fill = site_fill, colour = site_fill, linewidth = 0.05) +
    ggplot2::annotate("rect", xmin = zoom[["xmin"]], xmax = zoom[["xmax"]], ymin = zoom[["ymin"]], ymax = zoom[["ymax"]],
                      fill = NA, colour = ink$primary, linewidth = 0.4) +
    ggplot2::annotate("text", x = 145.08, y = -15.12, label = "Cape York\nPeninsula", size = txt, colour = ink$secondary, hjust = 0) +
    ggplot2::coord_sf(xlim = c(145.05, 145.9), ylim = c(-15.82, -14.25), expand = FALSE) +
    ggplot2::scale_x_continuous(breaks = c(145.2, 145.5, 145.8)) +
    ggplot2::scale_y_continuous(breaks = c(-14.5, -15, -15.5)) +
    theme_report(7) +
    ggplot2::theme(axis.title = ggplot2::element_blank(), axis.text = ggplot2::element_text(size = 6))
  land_loc <- suppressWarnings(sf::st_crop(land, c(xmin = 141, xmax = 154, ymin = -26, ymax = -10)))
  locator <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = land_loc, fill = ink$land, colour = NA) +
    ggplot2::annotate("rect", xmin = 145.05, xmax = 145.9, ymin = -15.82, ymax = -14.25, fill = NA, colour = ink$primary, linewidth = 0.4) +
    ggplot2::coord_sf(xlim = c(141.5, 153.5), ylim = c(-25.5, -10.5), expand = FALSE) +
    ggplot2::theme_void() +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = ink$surface, colour = ink$primary, linewidth = 0.3))
  land_zoom <- suppressWarnings(sf::st_crop(land, zoom))
  dom_zoom <- suppressWarnings(sf::st_crop(sf::st_make_valid(domain), zoom))
  zoom_map <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = land_zoom, fill = ink$land, colour = NA) +
    ggplot2::geom_sf(data = dom_zoom, fill = site_fill, colour = ink$surface, alpha = 0.8, linewidth = 0.25) +
    ggplot2::annotate("text", x = 145.4215, y = -14.6445, label = "Lizard Island", size = txt, colour = ink$primary, fontface = "plain", hjust = 0, vjust = 1) +
    ggplot2::coord_sf(xlim = zoom[c("xmin", "xmax")], ylim = zoom[c("ymin", "ymax")], expand = FALSE) +
    ggplot2::theme_void() +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = ink$surface, colour = ink$primary, linewidth = 0.4))
  panel_a <- dom_map +
    patchwork::inset_element(locator, left = 0.70, bottom = 0.79, right = 0.99, top = 0.99, align_to = "panel") +
    patchwork::inset_element(zoom_map, left = 0.01, bottom = 0.01, right = 0.37, top = 0.215, align_to = "panel")

  # (b) Mechanisms -----------------------------------------------------------
  st <- data.frame(stage = c("Recruit", "Juvenile", "Adult"), x = c(29, 49, 72), y = 50, r = c(1.6, 2.7, 4.6))
  stars <- do.call(rbind, lapply(seq_len(nrow(st)), function(i) concept_star(st$x[i], st$y[i], st$r[i], id = st$stage[i])))
  larvae <- data.frame(x = c(7.2, 8.8, 7.9, 9.6, 6.6, 8.4), y = c(51.5, 52.1, 49.6, 50.4, 49.9, 48.6))
  larvae_out <- transform(larvae, x = x + 86, y = y + 12)
  fast <- concept_fast_coral(70, 9, 1.25)
  slow <- concept_slow_coral(92, 9, 1.2)
  band <- data.frame(xmin = c(0, 0), xmax = c(100, 100), ymin = c(31, 0), ymax = c(77, 28.5), fill = c("#f3f3f3", "#f8f4ee"))
  seg <- function(x, y, xend, yend, linewidth = 0.45)
    ggplot2::annotate("segment", x = x, y = y, xend = xend, yend = yend, arrow = arr, colour = arrow_col, linewidth = linewidth)
  crv <- function(x, y, xend, yend, curvature, colour = arrow_col, linetype = "solid", lw = 0.45)
    ggplot2::annotate("curve", x = x, y = y, xend = xend, yend = yend, curvature = curvature, arrow = arr, colour = colour, linetype = linetype, linewidth = lw)
  lab <- function(x, y, label, hjust = 0.5, colour = ink$primary, size = txt, fontface = "plain", vjust = 0.5, parse = FALSE)
    ggplot2::annotate("text", x = x, y = y, label = label, hjust = hjust, vjust = vjust, colour = colour, size = size, fontface = fontface, lineheight = 0.9, parse = parse)
  loop <- function(x, y) crv(x + 1.4, y, x - 1.4, y, curvature = -2.4)
  panel_b <- ggplot2::ggplot() +
    ggplot2::geom_rect(data = band, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = band$fill) +
    lab(1.2, 75.6, "COTS on each site", hjust = 0, vjust = 1, fontface = "bold", colour = ink$secondary) +
    lab(1.2, 27.0, "Coral on each site", hjust = 0, vjust = 1, fontface = "bold", colour = ink$secondary) +
    # life cycle
    ggplot2::geom_point(data = larvae, ggplot2::aes(x, y), size = 0.7, colour = cots_ink) +
    ggplot2::geom_point(data = larvae_out, ggplot2::aes(x, y), size = 0.7, colour = cots_ink) +
    ggplot2::geom_polygon(data = stars, ggplot2::aes(x, y, group = id), fill = cots_ink, colour = NA) +
    lab(6.6, 45.4, "Larvae", fontface = "bold", size = nm) +
    lab(st$x, 50 + st$r + 2.3, st$stage, fontface = "bold", size = nm) +
    lab(94.4, 57.6, "Larvae", fontface = "bold", size = nm) +
    seg(11.5, 50, 26.4, 50) + lab(19, 52.6, "Settlement") +
    lab(20.0, 46.6, "lower on reefs\nexposed to swell", colour = ink$secondary) +
    seg(31.4, 50, 45.6, 50) + seg(52.4, 50, 66.6, 50) + lab(59.5, 52.4, "Maturation") +
    crv(76.8, 52.5, 90.6, 61.6, curvature = -0.3) +
    lab(99, 53.0, "Reproduction", hjust = 1) +
    lab(99, 48.0, "fewer eggs in poor\ncondition; fertilisation\nfails when sparse", hjust = 1, colour = ink$secondary) +
    # dispersal back to larvae (schematic)
    crv(93, 66.4, 10.2, 55.6, curvature = 0.16) +
    lab(55, 74.2, "Dispersal among sites and reefs; larval survival varies with chlorophyll") +
    crv(2.6, 63.0, 6.6, 54.6, curvature = -0.25) +
    lab(1.2, 66.6, "Supply from outside\nthe domain", hjust = 0, colour = ink$secondary) +
    # survival
    loop(30, 45.2) + loop(49, 44.1) + loop(79.6, 45.6) +
    lab(30, 38.8, "Survival") +
    lab(46, 37.4, "Predator pit:\npredation saturates") +
    lab(91.5, 40.2, "Starvation lags\npreferred-coral loss") +
    # culling
    crv(51.5, 47.0, 58.6, 40.8, curvature = 0.25, colour = cull_col, lw = 0.6) +
    crv(68.4, 46.6, 61.4, 40.8, curvature = -0.25, colour = cull_col, lw = 0.6) +
    lab(60, 38.4, "Culling", colour = cull_col, fontface = "bold") +
    # coral
    ggplot2::geom_polygon(data = fast, ggplot2::aes(x, y, group = id), fill = fast_fill, colour = "#b8724f", linewidth = 0.2) +
    ggplot2::geom_polygon(data = slow, ggplot2::aes(x, y, group = id), fill = slow_fill, colour = "#8f7d58", linewidth = 0.2) +
    lab(65.5, 5.6, "Fast-growing, preferred") +
    lab(65.5, 2.6, '"(tabular, corymbose "*italic(Acropora)*")"', parse = TRUE) +
    lab(92, 5.6, "Slow-growing") + lab(92, 2.6, "(massive, other)") +
    # grazing and feedback
    seg(72.6, 45, 72.6, 16.6, linewidth = 0.6) +
    lab(73.8, 21.5, "Grazing,\npreferred first", hjust = 0) +
    crv(75.6, 45.6, 88.4, 16.8, curvature = -0.08, linetype = "22") +
    lab(94.2, 21.0, "when preferred\ncoral is scarce", colour = ink$secondary) +
    crv(66.0, 16.4, 70.2, 45.2, curvature = -0.18, linetype = "22") +
    lab(52.6, 25.2, "Preferred cover sets\nbody condition and survival", colour = ink$secondary) +
    # coral processes
    crv(15.5, 17.0, 63.6, 12.6, curvature = -0.12) +
    lab(1.2, 17.2, "Growth and\nrecruitment", hjust = 0) +
    crv(15.5, 7.4, 63.6, 9.4, curvature = 0.12) +
    lab(1.2, 7.0, "Bleaching and\ncyclones", hjust = 0) +
    ggplot2::coord_fixed(xlim = c(0, 100), ylim = c(0, 77), expand = FALSE, clip = "off") +
    ggplot2::theme_void()

  # (c) Model coupling and annual order --------------------------------------
  steps <- data.frame(
    k = 1:6,
    label = c("Coral spawning\nand recruitment", "Bleaching", "Cyclones", "COTS grazing,\nsurvival, reproduction\nand culling",
              "Larval dispersal\namong sites", "Outside larval\nsupply"),
    owner = c("CoralBlox.jl", "ADRIA.jl", "ADRIA.jl", "COTSMod.jl", "COTSMod.jl", "ADRIA.jl"))
  steps$x <- c(8, 23, 36, 54, 73, 90)
  steps$w <- c(13, 10, 10, 18, 15, 13)
  owner_fill <- c(`ADRIA.jl` = "#e6e6e6", `CoralBlox.jl` = "#f6dccd", `COTSMod.jl` = "#cfcfcf")
  boxes <- transform(steps, xmin = x - w / 2, xmax = x + w / 2, ymin = 6, ymax = 18)
  panel_c <- ggplot2::ggplot() +
    ggplot2::geom_rect(data = boxes, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = owner), colour = NA) +
    ggplot2::scale_fill_manual(values = owner_fill, guide = "none") +
    ggplot2::geom_text(data = boxes, ggplot2::aes(x, 12, label = label), size = txt, lineheight = 0.9) +
    ggplot2::geom_text(data = boxes, ggplot2::aes(x, 3.6, label = owner), size = txt, fontface = "italic", colour = ink$secondary) +
    ggplot2::geom_segment(data = boxes[-nrow(boxes), ], ggplot2::aes(x = xmax + 0.3, xend = boxes$xmin[-1] - 0.3, y = 12, yend = 12),
                          arrow = arr, colour = arrow_col, linewidth = 0.45) +
    lab(0, 23.5, "Each annual step in ADRIA.jl, which holds the 2,895-site domain, disturbance forcing, connectivity and scenarios", hjust = 0, fontface = "bold", colour = ink$secondary) +
    ggplot2::coord_cartesian(xlim = c(0, 100), ylim = c(1.5, 25), expand = FALSE, clip = "off") +
    ggplot2::theme_void()

  top <- patchwork::wrap_elements(full = panel_a) + panel_b + patchwork::plot_layout(widths = c(1, 2.05))
  p <- top / panel_c + patchwork::plot_layout(heights = c(3.1, 1)) +
    patchwork::plot_annotation(tag_levels = list(c("a", "b", "c")))
  p <- p & ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 9))
  out <- file.path(cfg$dirs$figures, "fig36_metapop_concept.png")
  ggplot2::ggsave(out, p, width = 7.5, height = 5.7, dpi = 300, bg = ink$surface)
  invisible(out)
}
