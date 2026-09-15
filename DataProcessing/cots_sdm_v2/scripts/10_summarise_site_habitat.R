args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v2_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v2_dir, "R", "bootstrap.R"))
ctx <- initialise_v2(script_file)
cfg <- ctx$config

metrics_path <- file.path(
  cfg$paths$outputs_dir, "site_habitat_metrics.csv"
)
decomposition_path <- file.path(
  cfg$paths$outputs_dir,
  "site_habitat_scale_decomposition.csv"
)
sites_path <- file.path(
  cfg$paths$artifacts_dir, "site_habitat_all_sites.rds"
)
required <- c(metrics_path, decomposition_path, sites_path)
if (any(!file.exists(required))) {
  stop("Run scripts/09_fit_site_habitat.R first.")
}

metrics <- utils::read.csv(metrics_path, stringsAsFactors = FALSE)
decomposition <- utils::read.csv(
  decomposition_path, stringsAsFactors = FALSE
)
comparison <- dplyr::left_join(
  metrics,
  dplyr::select(
    decomposition,
    dplyr::all_of(c(
      "model", "between_reef_only_roc_auc",
      "within_reef_roc_auc"
    ))
  ),
  by = "model"
) |>
  dplyr::mutate(
    feature_family = dplyr::case_when(
      grepl("^clean_geomorphic", .data$model) ~ "Clean + geomorphic",
      grepl("^clean", .data$model) ~ "Clean",
      grepl("^rg_geomorphic", .data$model) ~ "RG + geomorphic",
      TRUE ~ "RG"
    ),
    training_evidence = dplyr::if_else(
      grepl("plus_manta$", .data$model),
      "Cull + manta zeros",
      "Cull only"
    )
  )

sites <- readRDS(sites_path) |>
  dplyr::filter(
    .data$clean_complete == 1L,
    .data$rg_complete == 1L
  )
geomorphic_summary <- sites |>
  dplyr::group_by(
    .data$geomorphic_class, .data$evidence_source
  ) |>
  dplyr::summarise(
    sites = dplyr::n(),
    reefs = dplyr::n_distinct(.data$reef_name),
    high_density_sites = sum(.data$high_density),
    high_density_prevalence = mean(.data$high_density),
    .groups = "drop"
  )

utils::write.csv(
  comparison,
  file.path(
    cfg$paths$outputs_dir,
    "site_habitat_model_comparison.csv"
  ),
  row.names = FALSE
)
utils::write.csv(
  geomorphic_summary,
  file.path(
    cfg$paths$outputs_dir,
    "site_habitat_geomorphic_summary.csv"
  ),
  row.names = FALSE
)

plot_data <- comparison |>
  dplyr::select(
    feature_family,
    training_evidence,
    overall = roc_auc,
    within_reef = within_reef_roc_auc
  ) |>
  tidyr::pivot_longer(
    cols = c("overall", "within_reef"),
    names_to = "validation_scale",
    values_to = "roc_auc"
  ) |>
  dplyr::mutate(
    validation_scale = factor(
      .data$validation_scale,
      levels = c("overall", "within_reef"),
      labels = c("Overall", "Within reef")
    ),
    feature_family = factor(
      .data$feature_family,
      levels = c(
        "Clean", "Clean + geomorphic",
        "RG", "RG + geomorphic"
      )
    )
  )

comparison_plot <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(
    x = .data$feature_family,
    y = .data$roc_auc,
    fill = .data$validation_scale
  )
) +
  ggplot2::geom_col(
    position = ggplot2::position_dodge(width = 0.78),
    width = 0.7
  ) +
  ggplot2::geom_hline(
    yintercept = 0.5,
    linetype = "dashed",
    colour = "#64748b"
  ) +
  ggplot2::facet_wrap(~training_evidence) +
  ggplot2::scale_fill_manual(
    values = c("Overall" = "#2563eb", "Within reef" = "#f59e0b")
  ) +
  ggplot2::coord_cartesian(ylim = c(0.5, 0.7)) +
  ggplot2::labs(
    title = "Site-scale COTS habitat model comparison",
    subtitle = "Held-out cull sites only; no reef-level predictors",
    x = NULL,
    y = "ROC-AUC",
    fill = "Validation scale"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 20, hjust = 1
    ),
    legend.position = "bottom"
  )

ggplot2::ggsave(
  file.path(
    cfg$paths$outputs_dir,
    "site_habitat_auc_comparison.png"
  ),
  comparison_plot,
  width = 10,
  height = 5.8,
  dpi = 180
)
message("Site-habitat summaries and comparison plot written.")
