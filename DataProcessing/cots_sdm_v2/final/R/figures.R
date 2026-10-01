# Figure styling for the technical report (static PNG, light surface).
#
# Colours follow the dataviz reference palette. Each delivered layer keeps one
# categorical hue in every figure (validated all-pairs for CVD: worst deutan
# Delta E 9.2). Aqua is below 3:1 contrast on the surface, so every chart that
# uses it also names the layer in text (axis label or direct label).

layer_labels <- c(
  legacy_clean_year2026 = "Clean + year",
  legacy_rg_year2026 = "RG + year",
  hmsc_joint = "HMSC joint"
)
layer_colours <- c(
  legacy_clean_year2026 = "#2a78d6",
  legacy_rg_year2026 = "#eb6834",
  hmsc_joint = "#1baf7a"
)

ink <- list(
  surface = "#fcfcfb", primary = "#0b0b0b", secondary = "#52514e",
  muted = "#898781", grid = "#e1e0d9", axis = "#c3c2b7", land = "#d9d8d2"
)

# Single-hue sequential ramps, light -> dark.
ramp_blue <- c("#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b")
ramp_orange <- c("#fde3d6", "#f8bfa3", "#f39670", "#eb6834", "#c24f20", "#933a16", "#62260d")
# Diverging blue <-> red with a neutral grey midpoint.
diverging_blue_red <- c("#184f95", "#3987e5", "#9ec5f4", "#f0efec", "#f2a9a8", "#e34948", "#9e2a29")
# Ordinal priority (primary -> fallback) for coverage bars.
ordinal_blue <- c("#184f95", "#3987e5", "#86b6ef")

report_font <- "Segoe UI"

theme_report <- function(base_size = 9) {
  ggplot2::theme_minimal(base_size = base_size, base_family = report_font) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = ink$surface, colour = NA),
      panel.background = ggplot2::element_rect(fill = ink$surface, colour = NA),
      panel.grid.major = ggplot2::element_line(colour = ink$grid, linewidth = 0.25),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(colour = ink$secondary),
      axis.title = ggplot2::element_text(colour = ink$secondary),
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(colour = ink$primary, face = "bold", size = base_size + 2),
      plot.subtitle = ggplot2::element_text(colour = ink$secondary),
      plot.caption = ggplot2::element_text(colour = ink$muted, hjust = 0),
      strip.text = ggplot2::element_text(colour = ink$primary, face = "bold", hjust = 0),
      legend.title = ggplot2::element_text(colour = ink$secondary),
      legend.text = ggplot2::element_text(colour = ink$secondary),
      plot.title.position = "plot",
      plot.caption.position = "plot"
    )
}

theme_map <- function(base_size = 9) {
  theme_report(base_size) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(colour = ink$grid, fill = NA, linewidth = 0.3)
    )
}

save_figure <- function(plot, cfg, name, width, height) {
  path <- file.path(cfg$dirs$figures, name)
  ggplot2::ggsave(path, plot, width = width, height = height, dpi = 300,
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
