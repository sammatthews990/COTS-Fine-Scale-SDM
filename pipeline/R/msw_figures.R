# Figure helpers for marine-science papers (marine-science-figures skill).
# Style: Matthews and McWilliam figures (references/figure_style.md).
# Copied from the marine-science-figures skill (13 November 2026).
# Requires ggplot2 and patchwork.

# Okabe-Ito (Okabe & Ito 2008; colour-blind safe) for classes; viridis
# (ggplot2::scale_*_viridis_c) for gradients. Colour the series that carry the
# finding; grey for context.
okabe_ito <- c(orange = "#E69F00", skyblue = "#56B4E9", green = "#009E73", yellow = "#F0E442",
               blue = "#0072B2", vermillion = "#D55E00", purple = "#CC79A7", grey = "#999999",
               black = "#000000")
msw_pal <- c(okabe_ito, light = "#D9D9D9", dark = "#4D4D4D")

# Discrete scales from a named key (names = levels) or Okabe-Ito in order.
scale_colour_okabe <- function(values = NULL, ...) {
  ggplot2::scale_colour_manual(values = if (is.null(values)) unname(okabe_ito) else values, ...)
}
scale_fill_okabe <- function(values = NULL, ...) {
  ggplot2::scale_fill_manual(values = if (is.null(values)) unname(okabe_ito) else values, ...)
}

# Classic axes, no gridlines, small type. base_size 8 suits figures saved at
# their final printed width (save_msw()).
theme_msw <- function(base_size = 8, base_family = "sans") {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(colour = "black", linewidth = 0.3),
      axis.ticks = ggplot2::element_line(colour = "black", linewidth = 0.3),
      axis.text = ggplot2::element_text(colour = "black"),
      plot.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1), hjust = 0),
      plot.title.position = "panel",
      plot.tag = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.25)),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.title = ggplot2::element_blank(),
      legend.key.size = ggplot2::unit(3, "mm"),
      legend.background = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
}

# Bold lower-case panel letters (a, b, c) for a patchwork composite.
msw_tags <- function(levels = "a") patchwork::plot_annotation(tag_levels = levels)

# Save at the final printed size. width: "single" (89 mm), "double" (183 mm),
# or millimetres; height in mm.
save_msw <- function(plot, path, width = "double", height = 120, dpi = 300) {
  w <- if (is.character(width)) c(single = 89, double = 183)[[width]] else width
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, plot, width = w, height = height, units = "mm", dpi = dpi, bg = "white")
  invisible(path)
}

# Report the common clutter in a ggplot or patchwork object (see the
# declutter checklist in SKILL.md). Returns the issues invisibly.
msw_lint <- function(p) {
  issues <- character()
  plots <- if (inherits(p, "patchwork")) c(list(p), p$patches$plots) else list(p)
  ann <- if (inherits(p, "patchwork")) p$patches$annotation else list()
  for (f in c("title", "subtitle", "caption")) {
    if (!is.null(ann[[f]]) && !inherits(ann[[f]], "waiver")) issues <- c(issues, paste0("figure-level ", f, " in the image"))
  }
  for (i in seq_along(plots)) {
    q <- plots[[i]]
    if (!inherits(q, "ggplot")) next
    lab <- q$labels
    if (!is.null(lab$subtitle)) issues <- c(issues, sprintf("panel %d: subtitle", i))
    if (!is.null(lab$caption)) issues <- c(issues, sprintf("panel %d: caption", i))
    if (!is.null(lab$title) && nchar(gsub("^[a-zA-Z]\\s+", "", lab$title)) > 35) {
      issues <- c(issues, sprintf("panel %d: long panel title (%d chars)", i, nchar(lab$title)))
    }
    txt <- unlist(lab[c("title", "x", "y")])
    if (any(grepl(">=|<=|_[a-z]", txt))) issues <- c(issues, sprintf("panel %d: code-style text in labels", i))
    n_text <- sum(vapply(q$layers, \(l) inherits(l$geom, c("GeomText", "GeomLabel")), logical(1)))
    if (n_text > 1) issues <- c(issues, sprintf("panel %d: %d text layers", i, n_text))
    b <- ggplot2::ggplot_build(q)
    cols <- unique(unlist(lapply(b$data, \(d) c(d$colour, d$fill))))
    cols <- setdiff(toupper(cols[!is.na(cols)]), c("BLACK", "WHITE", "GREY50", "#000000", "#FFFFFF"))
    if (length(cols) > 8) issues <- c(issues, sprintf("panel %d: %d colours (more than Okabe-Ito can separate)", i, length(cols)))
    fv <- setdiff(names(b$layout$layout), c("PANEL", "ROW", "COL", "SCALE_X", "SCALE_Y"))
    strips <- unlist(lapply(b$layout$layout[fv], as.character))
    if (length(strips) && any(grepl("[0-9]\\.[0-9]|\\(|×", as.character(strips)))) {
      issues <- c(issues, sprintf("panel %d: statistics in facet strips", i))
    }
  }
  if (length(issues)) message(paste("-", unique(issues), collapse = "\n")) else message("msw_lint: no issues")
  invisible(unique(issues))
}
