args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v2_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v2_dir, "R", "bootstrap.R"))
ctx <- initialise_v2(script_file)
cfg <- ctx$config

if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")
fig_dir <- file.path(v2_dir, "reports", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

metrics <- utils::read.csv(file.path(
  cfg$paths$outputs_dir, "legacy_unweighted_scale_metrics.csv"
))
decomp <- utils::read.csv(file.path(
  cfg$paths$outputs_dir, "legacy_unweighted_scale_decomposition.csv"
))
capture <- utils::read.csv(file.path(
  cfg$paths$outputs_dir, "legacy_unweighted_scale_capture.csv"
)) |>
  dplyr::filter(abs(.data$target_fraction - 0.2) < 1e-8) |>
  dplyr::mutate(
    eligible_prevalence = .data$eligible_positives / .data$eligible_rows,
    top20_lift = .data$selected_prevalence / .data$eligible_prevalence
  )

labels <- c(
  clean_30m_year__legacy_unweighted = "30 m + year",
  clean_60m_year__legacy_unweighted = "60 m + year",
  clean_120m_year__legacy_unweighted = "120 m + year",
  multiscale_topography_year__legacy_unweighted = "Multiscale + year",
  multiscale_topography_static__legacy_unweighted = "Multiscale static"
)
auc <- metrics |>
  dplyr::select("model", overall_auc = "roc_auc") |>
  dplyr::left_join(
    dplyr::select(decomp, "model", within_auc = "within_reef_roc_auc"),
    by = "model"
  ) |>
  tidyr::pivot_longer(
    c("overall_auc", "within_auc"), names_to = "scale", values_to = "auc"
  ) |>
  dplyr::mutate(
    model_label = unname(labels[.data$model]),
    scale = dplyr::recode(
      .data$scale, overall_auc = "Overall", within_auc = "Within reef-year"
    ),
    model_label = factor(.data$model_label, levels = unname(labels))
  )

p <- ggplot2::ggplot(auc, ggplot2::aes(.data$model_label, .data$auc, fill = .data$scale)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.78), width = 0.7) +
  ggplot2::geom_hline(yintercept = 0.5, linetype = 2, colour = "grey45") +
  ggplot2::coord_cartesian(ylim = c(0.48, 0.69)) +
  ggplot2::scale_fill_manual(values = c("Overall" = "#2166ac", "Within reef-year" = "#b2182b")) +
  ggplot2::labs(
    x = NULL, y = "Held-out ROC-AUC", fill = NULL,
    title = "Overall discrimination and fine-scale within-reef ranking are different tasks",
    subtitle = "All candidates use the legacy hyperparameters and unweighted annual Cull target"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 25, hjust = 1),
    panel.grid.major.x = ggplot2::element_blank(),
    legend.position = "top"
  )
ggplot2::ggsave(
  file.path(fig_dir, "multiscale_auc_comparison.png"), p,
  width = 9.2, height = 5.5, dpi = 180
)

lift <- capture |>
  dplyr::mutate(
    model_label = factor(unname(labels[.data$model]), levels = unname(labels))
  )
p <- ggplot2::ggplot(lift, ggplot2::aes(.data$model_label, .data$top20_lift)) +
  ggplot2::geom_col(width = 0.65, fill = "#1b7837") +
  ggplot2::geom_hline(yintercept = 1, linetype = 2, colour = "grey45") +
  ggplot2::labs(
    x = NULL, y = "High-density lift in top 20%", title = "Operational within-reef capture"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 25, hjust = 1),
    panel.grid.major.x = ggplot2::element_blank()
  )
ggplot2::ggsave(
  file.path(fig_dir, "multiscale_top20_lift.png"), p,
  width = 8.5, height = 4.8, dpi = 180
)

maxent <- utils::read.csv(file.path(
  cfg$paths$outputs_dir, "maxent_multiscale_metrics.csv"
)) |>
  dplyr::mutate(
    evidence = dplyr::case_when(
      grepl("high_density", .data$model) ~ "Cull >0.02",
      grepl("direct_cots", .data$model) ~ "Direct COTS",
      TRUE ~ "COTS or scars"
    ),
    scale = ifelse(grepl("multiscale", .data$model), "Selected multiscale", "Native 30 m")
  )
p <- ggplot2::ggplot(maxent, ggplot2::aes(.data$evidence, .data$roc_auc, fill = .data$scale)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75), width = 0.68) +
  ggplot2::geom_hline(yintercept = 0.5, linetype = 2, colour = "grey45") +
  ggplot2::coord_cartesian(ylim = c(0.49, 0.54)) +
  ggplot2::labs(
    x = NULL, y = "Descriptive AUC on annual Cull rows", fill = NULL,
    title = "Presence-only screening does not recover the high-density target",
    subtitle = "Full-fit diagnostic only; values are not held-out estimates"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(), legend.position = "top")
ggplot2::ggsave(
  file.path(fig_dir, "maxent_descriptive_comparison.png"), p,
  width = 8.2, height = 4.8, dpi = 180
)

hmsc_performance_path <- file.path(
  cfg$paths$outputs_dir, "hmsc_joint_200m_cv_performance.csv"
)
if (file.exists(hmsc_performance_path)) {
  hmsc <- utils::read.csv(hmsc_performance_path) |>
    dplyr::mutate(
      response_label = dplyr::recode(
        .data$response,
        cull_direct_presence = "Cull direct",
        cull_high_density = "Cull >0.02",
        manta_direct_presence = "Manta direct",
        manta_scar_presence = "Manta scars"
      )
    )
  p <- ggplot2::ggplot(hmsc, ggplot2::aes(.data$response_label, .data$roc_auc)) +
    ggplot2::geom_col(width = 0.65, fill = "#762a83") +
    ggplot2::geom_hline(yintercept = 0.5, linetype = 2, colour = "grey45") +
    ggplot2::coord_cartesian(ylim = c(0.45, max(0.7, hmsc$roc_auc + 0.03))) +
    ggplot2::labs(
      x = NULL, y = "Five-fold cell-blocked ROC-AUC",
      title = "Joint Hmsc held-out performance by evidence type"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())
  ggplot2::ggsave(
    file.path(fig_dir, "hmsc_joint_cv_auc.png"), p,
    width = 7.7, height = 4.8, dpi = 180
  )
}

message("Report figures complete.")
