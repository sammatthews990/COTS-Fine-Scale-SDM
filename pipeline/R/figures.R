# Figure styling for the report and manuscript (static PNG, white surface).
#
# Journal style (marine-science-figures skill, 13 November 2026): classic
# axes, no gridlines, no figure title, subtitle or caption in the image (the
# legend text carries them), bold lower-case panel letters, Okabe-Ito for
# classes and viridis for gradients. Helpers from msw_figures.R.

layer_labels <- c(
  legacy_clean_year2026 = "Clean + year",
  legacy_rg_year2026 = "RG + year",
  hmsc_joint = "HMSC joint"
)

# Project colour key: one meaning per colour in every figure.
#   orange       habitat alone (cull-trained habitat layer)
#   blue         habitat + recent evidence (the ranking; within-reef classes)
#   purple       density-aware classes, first version
#   bluish green density-aware classes, second version (recommended)
#   grey         random selection or other baseline; black systematic thinning
#   vermillion / sky blue  metacommunity runs with / without recorded culling
key_colours <- c(habitat = okabe_ito[["orange"]], combined = okabe_ito[["blue"]],
                 density_v1 = okabe_ito[["purple"]], density_v2 = okabe_ito[["green"]],
                 random = okabe_ito[["grey"]], systematic = "#000000", evidence = "#666666",
                 culled = okabe_ito[["vermillion"]], unculled = okabe_ito[["skyblue"]])
layer_colours <- c(
  legacy_clean_year2026 = okabe_ito[["orange"]],
  legacy_rg_year2026 = okabe_ito[["vermillion"]],
  hmsc_joint = okabe_ito[["green"]]
)

ink <- list(
  surface = "#ffffff", primary = "#000000", secondary = "#4d4d4d",
  muted = "#808080", grid = "#e1e0d9", axis = "#c3c2b7", land = "#d9d8d2"
)

# Single-hue sequential ramps, light -> dark (kept for legacy panels; new
# gradients use viridis).
ramp_blue <- c("#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b")
ramp_orange <- c("#fde3d6", "#f8bfa3", "#f39670", "#eb6834", "#c24f20", "#933a16", "#62260d")
# Diverging blue <-> vermillion with a neutral grey midpoint.
diverging_blue_red <- c("#0072B2", "#56B4E9", "#a9d6f2", "#f0efec", "#f2b38c", "#E69F00", "#D55E00")
# Ordinal priority (primary -> fallback) for coverage bars.
ordinal_blue <- c("#184f95", "#3987e5", "#86b6ef")

report_font <- "Segoe UI"

theme_report <- function(base_size = 9) {
  ggplot2::theme_classic(base_size = base_size, base_family = report_font) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = ink$surface, colour = NA),
      panel.background = ggplot2::element_rect(fill = ink$surface, colour = NA),
      axis.line = ggplot2::element_line(colour = "black", linewidth = 0.3),
      axis.ticks = ggplot2::element_line(colour = "black", linewidth = 0.3),
      axis.text = ggplot2::element_text(colour = "black"),
      axis.title = ggplot2::element_text(colour = "black"),
      plot.title = ggplot2::element_text(colour = ink$primary, face = "bold", size = base_size, hjust = 0),
      plot.subtitle = ggplot2::element_text(colour = ink$secondary),
      plot.caption = ggplot2::element_text(colour = ink$muted, hjust = 0),
      plot.tag = ggplot2::element_text(face = "bold", size = base_size + 2),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(colour = ink$primary, face = "bold", hjust = 0),
      legend.title = ggplot2::element_text(colour = "black"),
      legend.text = ggplot2::element_text(colour = "black"),
      legend.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      plot.title.position = "panel",
      plot.caption.position = "plot"
    )
}

theme_map <- function(base_size = 9) {
  theme_report(base_size) +
    ggplot2::theme(
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank()
    )
}

# Drop figure-level text: the manuscript and report legends carry the title,
# what is plotted and the caveats. Panel headings and panel letters stay.
strip_figure_text <- function(plot) {
  if (inherits(plot, "patchwork")) {
    for (f in c("title", "subtitle", "caption")) plot$patches$annotation[f] <- list(NULL)
  } else if (inherits(plot, "ggplot")) {
    for (f in c("title", "subtitle", "caption")) plot$labels[f] <- list(NULL)
  }
  plot
}

save_figure <- function(plot, cfg, name, width, height) {
  path <- file.path(cfg$dirs$figures, name)
  ggplot2::ggsave(path, strip_figure_text(plot), width = width, height = height, dpi = 300,
                  device = ragg::agg_png, bg = ink$surface)
  log_message("Wrote ", path)
  invisible(path)
}

# Raster cells to a data frame for geom_raster.
raster_df <- function(x, name = "value") {
  d <- terra::as.data.frame(x, xy = TRUE, na.rm = TRUE)
  names(d)[3] <- name
  d
}

# Pooled within-group AUC (Mann-Whitney), weighting each group by its pairs.
within_group_auc <- function(score, truth, group) {
  pieces <- lapply(split(seq_along(score), group), function(i) {
    y <- truth[i]
    n1 <- sum(y)
    n0 <- sum(!y)
    if (n1 == 0 || n0 == 0) return(NULL)
    r <- rank(score[i])
    c(auc = (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0), pairs = n1 * n0)
  })
  pieces <- do.call(rbind, pieces)
  if (is.null(pieces)) return(NA_real_)
  stats::weighted.mean(pieces[, "auc"], pieces[, "pairs"])
}

# Share of positives falling in each group's top 20 % of scores, as a lift
# over the 0.2 expected by chance.
within_group_top20_lift <- function(score, truth, group) {
  hits <- 0
  total <- 0
  for (i in split(seq_along(score), group)) {
    y <- truth[i]
    if (sum(y) == 0 || sum(!y) == 0) next
    top <- rank(-score[i], ties.method = "random") <= ceiling(0.2 * length(i))
    hits <- hits + sum(y & top)
    total <- total + sum(y)
  }
  (hits / total) / 0.2
}
