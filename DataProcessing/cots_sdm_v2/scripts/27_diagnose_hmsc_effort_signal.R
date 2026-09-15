args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v2_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v2_dir, "R", "bootstrap.R"))
ctx <- initialise_v2(script_file)
cfg <- ctx$config

response_names <- c(
  "cull_direct_presence", "cull_high_density",
  "manta_direct_presence", "manta_scar_presence"
)
habitat_names <- c(
  "BPI_s30m", "SLO_s60m", "VRM_s60m", "DEM_s30m",
  "SVF_s120m", "coral_suitability_mean"
)
effort_names <- c(
  "log_cull_bottom_time", "log_cull_dives",
  "log_manta_tows", "log_manta_distance",
  "manta_hard_coral_fraction"
)
covariate_names <- c("year_z", habitat_names, effort_names)

data <- readRDS(file.path(
  cfg$paths$artifacts_dir, "hmsc_joint_200m_cell_year.rds"
)) |>
  dplyr::filter(stats::complete.cases(dplyr::pick(
    dplyr::all_of(c(response_names, covariate_names, "reef_name", "cell_id"))
  ))) |>
  dplyr::arrange(.data$cell_id, .data$year) |>
  dplyr::mutate(
    reef_year_id = paste(.data$reef_name, .data$year, sep = "||"),
    site_year_id = paste(.data$cell_id, .data$year, sep = "||")
  )
partition <- as.integer(readRDS(file.path(
  cfg$paths$artifacts_dir, "hmsc_joint_200m_cell_partition.rds"
)))
hmsc <- readRDS(file.path(
  cfg$paths$outputs_dir, "hmsc_joint_200m_high_density_predictions.rds"
))

stopifnot(nrow(data) == length(partition), nrow(data) == nrow(hmsc))
match_index <- match(data$site_year_id, hmsc$site_year_id)
if (anyNA(match_index)) stop("Hmsc predictions do not match the model data.")
hmsc_probability <- hmsc$probability[match_index]

feature_sets <- list(
  glm_static_habitat = habitat_names,
  glm_habitat_year = c("year_z", habitat_names),
  glm_effort_year = c("year_z", effort_names),
  glm_habitat_effort_year = c("year_z", habitat_names, effort_names)
)

fit_fold_glm <- function(features, name) {
  predictions <- rep(NA_real_, nrow(data))
  for (fold in sort(unique(partition))) {
    train <- partition != fold
    test <- partition == fold
    formula <- stats::reformulate(features, response = "cull_high_density")
    model <- suppressWarnings(stats::glm(
      formula, data = data[train, , drop = FALSE], family = stats::binomial()
    ))
    predictions[test] <- stats::predict(
      model, newdata = data[test, , drop = FALSE], type = "response"
    )
  }
  tibble::tibble(
    model = name,
    validation = "five_fold_held_out_200m_cell",
    fold = partition,
    site_year_id = data$site_year_id,
    site_id = data$cell_id,
    reef_name = data$reef_name,
    reef_year_id = data$reef_year_id,
    year = data$year,
    truth = data$cull_high_density,
    probability = predictions
  )
}

ablation_predictions <- dplyr::bind_rows(lapply(names(feature_sets), function(name) {
  fit_fold_glm(feature_sets[[name]], name)
}))
ablation_predictions <- dplyr::bind_rows(
  ablation_predictions,
  dplyr::transmute(
    data,
    model = "hmsc_joint_detection_conditioned",
    validation = "five_fold_held_out_200m_cell",
    fold = partition,
    site_year_id = .data$site_year_id,
    site_id = .data$cell_id,
    reef_name = .data$reef_name,
    reef_year_id = .data$reef_year_id,
    year = .data$year,
    truth = .data$cull_high_density,
    probability = hmsc_probability
  )
)

metrics <- suppressWarnings(summarise_benchmark_predictions(ablation_predictions))
decomposition <- suppressWarnings(prediction_scale_decomposition(
  ablation_predictions, group = "reef_year_id"
))
capture <- within_group_capture_metrics(
  ablation_predictions, group = "reef_year_id"
)

within_demean <- function(x, group) {
  x - ave(x, group, FUN = function(z) mean(z, na.rm = TRUE))
}
correlations <- dplyr::bind_rows(lapply(c("year_z", habitat_names, effort_names), function(name) {
  tibble::tibble(
    covariate = name,
    spearman_overall = stats::cor(
      hmsc_probability, data[[name]], method = "spearman", use = "complete.obs"
    ),
    pearson_within_reef_year = stats::cor(
      within_demean(hmsc_probability, data$reef_year_id),
      within_demean(data[[name]], data$reef_year_id),
      method = "pearson", use = "complete.obs"
    )
  )
}))

tables <- list(
  hmsc_effort_ablation_metrics = metrics,
  hmsc_effort_ablation_decomposition = decomposition,
  hmsc_effort_ablation_capture = capture,
  hmsc_prediction_covariate_correlations = correlations
)
for (name in names(tables)) {
  utils::write.csv(
    tables[[name]], file.path(cfg$paths$outputs_dir, paste0(name, ".csv")),
    row.names = FALSE
  )
}
saveRDS(ablation_predictions, file.path(
  cfg$paths$outputs_dir, "hmsc_effort_ablation_predictions.rds"
))

if (requireNamespace("ggplot2", quietly = TRUE)) {
  fig_dir <- file.path(v2_dir, "reports", "figures")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  labels <- c(
    glm_static_habitat = "Static habitat",
    glm_habitat_year = "Habitat + year",
    glm_effort_year = "Effort + year",
    glm_habitat_effort_year = "Habitat + effort + year",
    hmsc_joint_detection_conditioned = "Joint Hmsc"
  )
  plot_data <- dplyr::select(metrics, "model", overall = "roc_auc") |>
    dplyr::left_join(
      dplyr::select(decomposition, "model", within = "within_reef_roc_auc"),
      by = "model"
    ) |>
    tidyr::pivot_longer(
      c("overall", "within"), names_to = "evaluation_scale", values_to = "auc"
    ) |>
    dplyr::mutate(
      model_label = factor(unname(labels[.data$model]), levels = unname(labels)),
      evaluation_scale = dplyr::recode(
        .data$evaluation_scale,
        overall = "Overall", within = "Within reef-year"
      )
    )
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(.data$model_label, .data$auc, fill = .data$evaluation_scale)
  ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = 0.78), width = 0.7
    ) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = 2, colour = "grey45") +
    ggplot2::scale_fill_manual(values = c(
      "Overall" = "#2166ac", "Within reef-year" = "#b2182b"
    )) +
    ggplot2::coord_cartesian(ylim = c(0.48, 0.88)) +
    ggplot2::labs(
      x = NULL, y = "Five-fold cell-blocked ROC-AUC", fill = NULL,
      title = "Realised survey effort, not habitat, drives the Hmsc ranking"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 25, hjust = 1),
      panel.grid.major.x = ggplot2::element_blank(),
      legend.position = "top"
    )
  ggplot2::ggsave(
    file.path(fig_dir, "hmsc_effort_ablation_auc.png"), p,
    width = 9, height = 5.2, dpi = 180
  )
}

message("Hmsc effort-signal diagnostic complete.")
print(dplyr::left_join(
  metrics,
  dplyr::select(decomposition, "model", "within_reef_roc_auc"),
  by = "model"
))
print(correlations)
