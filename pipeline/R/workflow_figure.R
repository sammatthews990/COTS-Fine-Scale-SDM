# Proposed on-water workflow for a reef visit, drawn as a staged decision tree
# in the style of the IPM operational strategy (Fletcher et al. 2020).
#
# Question: where does the ranking enter a reef visit? Method: the decision
# steps of the chosen protocol (scripts 35-36; density-aware classes, A-D dive
# rule, A-C follow-up tow zone with season memory). Validation: none; this is a
# description of the proposed rules, back-tested in script 36 (Appendix D).

build_workflow_figure <- function(cfg) {
  txt <- 2.55
  new_fill <- "#d4e6f3"
  new_line <- okabe_ito[["blue"]]
  old_fill <- ink$surface
  old_line <- ink$muted
  arr <- grid::arrow(length = grid::unit(1.5, "mm"), type = "closed")

  # Boxes: id, centre, size, new step?, dashed?, label.
  box <- tibble::tribble(
    ~id,       ~x,   ~y,   ~w,  ~h,  ~new,  ~dash, ~label,
    "plan",    55,   93.5, 62,  8,   TRUE,  FALSE, "Refresh the ranking to the departure date;\nissue A–E site classes, tow zones and exploration zones",
    "open",    55,   80.5, 62,  6.5, FALSE, FALSE, "Dive open sites first (last dive above 0.04)\nand cull them below 0.04",
    "full",    29,   57,   30,  7,   FALSE, FALSE, "Tow the whole\nreef perimeter",
    "zone",    75,   57,   40,  8.5, TRUE,  FALSE, "Tow only the A–C zone, open sites\nand sites near strong tows\nearlier in the season",
    "explore", 80.5, 46,   29,  5.5, TRUE,  TRUE,  "Use tows saved to explore\nshaded habitat",
    "dive1",   80,   31,   24,  6,   FALSE, FALSE, "Dive",
    "dive2",   80,   15,   30,  6,   FALSE, FALSE, "Dive at crew discretion",
    "defer",   40,   3.8,  38,  5.2,   TRUE,  FALSE, "Defer to the next visit\nor the annual check")
  dia <- tibble::tribble(
    ~id,      ~x, ~y,   ~hw, ~hh, ~new,  ~label,
    "first",  55, 68.5, 17,  5.8, FALSE, "First visit to this\nreef this season?",
    "strong", 40, 31,   19,  6.5, FALSE, "Strong tow within 100 m?\nCOTS, two or more scars\nor scars common",
    "rank",   40, 15,   15,  5.4, TRUE,  "Site ranked A–D?\n(density-aware\nclasses)")
  dia_poly <- do.call(rbind, lapply(seq_len(nrow(dia)), function(i) with(dia[i, ],
    data.frame(id = id, new = new, x = c(x, x + hw, x, x - hw), y = c(y + hh, y, y - hh, y)))))

  # Connectors: each path ends in an arrow; labels sit beside the branch.
  path <- function(id, xs, ys, dash = FALSE) data.frame(id = id, x = xs, y = ys, dash = dash)
  paths <- rbind(
    path("p1", c(55, 55), c(89.5, 83.75)),
    path("p2", c(55, 55), c(77.25, 74.3)),
    path("p3", c(38, 29, 29), c(68.5, 68.5, 60.5)),
    path("p4", c(72, 75, 75), c(68.5, 68.5, 61.25)),
    path("p5", c(80.5, 80.5), c(52.75, 48.75), dash = TRUE),
    path("p6", c(29, 29, 40, 40), c(53.5, 42, 42, 37.5)),
    path("p7", c(58, 58, 40), c(52.75, 42, 42)),
    path("p8", c(59, 68), c(31, 31)),
    path("p9", c(40, 40), c(24.5, 20.4)),
    path("p10", c(55, 65), c(15, 15)),
    path("p11", c(40, 40), c(9.6, 6.4)))
  ends <- c("p1", "p2", "p3", "p4", "p5", "p6", "p8", "p9", "p10", "p11")
  seglab <- tibble::tribble(
    ~x,   ~y,   ~label,
    33,   70.4, "yes",
    75.5, 70.4, "no",
    63.5, 33,   "yes",
    42.6, 22.4, "no",
    60,   17,   "yes",
    42.6, 8,    "no")

  bands <- tibble::tribble(
    ~ymin, ~ymax, ~label,               ~fill,
    88,    99,    "Before the\nvoyage",  "#f2f2f2",
    41,    87,    "On arrival",          "#f8f8f8",
    0,     40,    "For each\nremaining\nsite", "#f2f2f2")

  key <- tibble::tibble(x = c(68, 82), y = 4, new = c(TRUE, FALSE), label = c("New step", "Current practice"))

  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(data = bands, ggplot2::aes(xmin = 0, xmax = 100, ymin = ymin, ymax = ymax), fill = bands$fill) +
    ggplot2::geom_text(data = bands, ggplot2::aes(x = 1.2, y = (ymin + ymax) / 2, label = label), hjust = 0,
                       size = txt, fontface = "bold", colour = ink$secondary, lineheight = 0.9) +
    ggplot2::geom_path(data = paths[!paths$id %in% ends | TRUE, ], ggplot2::aes(x, y, group = id, linetype = dash),
                       colour = ink$secondary, linewidth = 0.45) +
    ggplot2::geom_path(data = do.call(rbind, lapply(ends, function(i) utils::tail(paths[paths$id == i, ], 2))),
                       ggplot2::aes(x, y, group = id, linetype = dash), colour = ink$secondary, linewidth = 0.45, arrow = arr) +
    ggplot2::scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "22"), guide = "none") +
    ggplot2::geom_text(data = seglab, ggplot2::aes(x, y, label = label), size = txt, colour = ink$secondary) +
    ggplot2::geom_rect(data = box, ggplot2::aes(xmin = x - w / 2, xmax = x + w / 2, ymin = y - h / 2, ymax = y + h / 2),
                       fill = ifelse(box$new, new_fill, old_fill), colour = ifelse(box$new, new_line, old_line),
                       linetype = ifelse(box$dash, "22", "solid"), linewidth = 0.5) +
    ggplot2::geom_polygon(data = dia_poly, ggplot2::aes(x, y, group = id),
                          fill = ifelse(dia_poly$new, new_fill, old_fill), colour = ifelse(dia_poly$new, new_line, old_line), linewidth = 0.5) +
    ggplot2::geom_text(data = box, ggplot2::aes(x, y, label = label), size = txt, lineheight = 0.92) +
    ggplot2::geom_text(data = dia, ggplot2::aes(x, y, label = label), size = txt, lineheight = 0.92) +
    ggplot2::geom_tile(data = key, ggplot2::aes(x, y), width = 2.6, height = 2.6,
                       fill = ifelse(key$new, new_fill, old_fill), colour = ifelse(key$new, new_line, old_line), linewidth = 0.4) +
    ggplot2::geom_text(data = key, ggplot2::aes(x + 2.2, y, label = label), hjust = 0, size = txt, colour = ink$secondary) +
    ggplot2::coord_cartesian(xlim = c(0, 100), ylim = c(0, 99), expand = FALSE, clip = "off") +
    ggplot2::theme_void() +
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = ink$surface, colour = NA))
  out <- file.path(cfg$dirs$figures, "fig23_workflow_decision.png")
  ggplot2::ggsave(out, p, width = 7.5, height = 6.4, dpi = 300, bg = ink$surface)
  invisible(out)
}
