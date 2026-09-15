args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v2_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v2_dir, "R", "bootstrap.R"))
ctx <- initialise_v2(script_file)
cfg <- ctx$config

required <- file.path(cfg$paths$outputs_dir, c(
  "two_part_static_repeat_decomposition.csv",
  "two_part_static_repeat_capture.csv",
  "two_part_component_repeat_decomposition.csv",
  "two_part_component_repeat_metrics.csv",
  "two_part_component_coefficients.csv"
))
if (any(!file.exists(required))) {
  stop("Run script 28 before building two-part diagnostic figures.")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("ggplot2 is required for diagnostic figures.")
}

fig_dir <- file.path(v2_dir, "reports", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

static_decomposition <- utils::read.csv(required[[1]])
static_capture <- utils::read.csv(required[[2]]) |>
  dplyr::filter(.data$target_fraction == 0.20)
component_decomposition <- utils::read.csv(required[[3]])
component_metrics <- utils::read.csv(required[[4]])
coefficients <- utils::read.csv(required[[5]])

static_labels <- c(
  clean_30m_annual = "Clean 30 m",
  clean_120m_annual = "Clean 120 m",
  multiscale_topography_annual = "Multiscale annual",
  multiscale_topography_propensity = "Multiscale site propensity"
)
static_plot <- static_decomposition |>
  dplyr::left_join(
    dplyr::select(
      static_capture,
      "model", "validation_repeat", "lift"
    ),
    by = c("model", "validation_repeat")
  ) |>
  dplyr::mutate(
    model_label = factor(
      unname(static_labels[.data$model]),
      levels = unname(static_labels)
    )
  )

p_static <- ggplot2::ggplot(
  static_plot,
  ggplot2::aes(
    .data$model_label, .data$within_reef_roc_auc,
    group = .data$validation_repeat
  )
) +
  ggplot2::geom_hline(
    yintercept = 0.568344817887779,
    linetype = 2, colour = "#555555"
  ) +
  ggplot2::geom_line(
    colour = "#9ecae1", linewidth = 0.55
  ) +
  ggplot2::geom_point(
    ggplot2::aes(colour = factor(.data$validation_repeat)),
    size = 2.6
  ) +
  ggplot2::stat_summary(
    ggplot2::aes(group = 1),
    fun = mean, geom = "point", shape = 18, size = 4,
    colour = "#b2182b"
  ) +
  ggplot2::scale_colour_brewer(
    palette = "Blues", name = "Blocked repeat"
  ) +
  ggplot2::labs(
    x = NULL, y = "Within-reef-year ROC-AUC",
    title = "Static habitat candidates across repeated site-blocked folds",
    subtitle = "Dashed line: legacy clean-year benchmark; diamond: repeat mean"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 20, hjust = 1),
    panel.grid.major.x = ggplot2::element_blank(),
    legend.position = "top"
  )
ggplot2::ggsave(
  file.path(fig_dir, "two_part_static_repeated_auc.png"),
  p_static, width = 9, height = 5.4, dpi = 180
)

component_labels <- c(
  habitat_static = "Static habitat",
  state_bridge = "Habitat + outbreak state",
  fixed_effort_pooled = "Fixed-effort pooled",
  realised_effort_diagnostic = "Realised effort (diagnostic)",
  fixed_protocol_dive_detection = "Fixed-protocol dive detection"
)
component_plot <- component_decomposition |>
  dplyr::mutate(
    component = sub("^.*__", "", .data$model),
    component_label = factor(
      unname(component_labels[.data$component]),
      levels = unname(component_labels)
    )
  ) |>
  dplyr::select(
    "component_label", "validation_repeat",
    Overall = "overall_roc_auc",
    `Within reef-year` = "within_reef_roc_auc"
  ) |>
  tidyr::pivot_longer(
    c("Overall", "Within reef-year"),
    names_to = "evaluation_scale", values_to = "auc"
  )

p_component <- ggplot2::ggplot(
  component_plot,
  ggplot2::aes(
    .data$component_label, .data$auc,
    colour = .data$evaluation_scale
  )
) +
  ggplot2::geom_hline(
    yintercept = 0.5, linetype = 3, colour = "grey55"
  ) +
  ggplot2::geom_point(
    ggplot2::aes(shape = factor(.data$validation_repeat)),
    position = ggplot2::position_dodge(width = 0.45),
    size = 2.5
  ) +
  ggplot2::stat_summary(
    ggplot2::aes(group = .data$evaluation_scale),
    fun = mean, geom = "line", linewidth = 0.8
  ) +
  ggplot2::scale_colour_manual(
    values = c("Overall" = "#2166ac", "Within reef-year" = "#b2182b")
  ) +
  ggplot2::labs(
    x = NULL, y = "Held-out ROC-AUC",
    colour = NULL, shape = "Blocked repeat",
    title = "Separated model components under nested validation"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 22, hjust = 1),
    panel.grid.major.x = ggplot2::element_blank(),
    legend.position = "top"
  )
ggplot2::ggsave(
  file.path(fig_dir, "two_part_component_auc.png"),
  p_component, width = 10, height = 5.7, dpi = 180
)

effort_terms <- coefficients |>
  dplyr::filter(
    .data$component == "effort_standardiser",
    .data$term %in% c("log_dives", "log_bottom_time")
  ) |>
  dplyr::mutate(
    term = dplyr::recode(
      .data$term,
      log_dives = "Number of dives",
      log_bottom_time = "Bottom time"
    )
  )
p_effort <- ggplot2::ggplot(
  effort_terms,
  ggplot2::aes(.data$term, .data$estimate, fill = .data$term)
) +
  ggplot2::geom_hline(
    yintercept = 0, linetype = 2, colour = "grey50"
  ) +
  ggplot2::geom_boxplot(width = 0.55, outlier.shape = NA) +
  ggplot2::geom_jitter(width = 0.08, alpha = 0.65, size = 1.6) +
  ggplot2::scale_fill_manual(
    values = c(
      "Number of dives" = "#67a9cf",
      "Bottom time" = "#ef8a62"
    ),
    guide = "none"
  ) +
  ggplot2::labs(
    x = NULL, y = "Nested-fold log-odds coefficient",
    title = "Survey-effort effects estimated then held fixed for mapping"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())
ggplot2::ggsave(
  file.path(fig_dir, "two_part_effort_coefficients.png"),
  p_effort, width = 7.5, height = 4.8, dpi = 180
)

component_calibration <- component_metrics |>
  dplyr::mutate(
    component = sub("^.*__", "", .data$model),
    component_label = unname(component_labels[.data$component])
  )
utils::write.csv(
  component_calibration,
  file.path(
    cfg$paths$outputs_dir,
    "two_part_component_calibration_summary.csv"
  ),
  row.names = FALSE
)

message("Two-part diagnostic figures complete.")
