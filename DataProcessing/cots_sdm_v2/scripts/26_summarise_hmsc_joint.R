args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v2_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v2_dir, "R", "bootstrap.R"))
ctx <- initialise_v2(script_file)
cfg <- ctx$config

model <- readRDS(file.path(
  cfg$paths$models_dir, "hmsc_joint_200m_production_model.rds"
))
response_names <- colnames(model$Y)
level_names <- names(model$ranLevels)
partition <- readRDS(file.path(
  cfg$paths$artifacts_dir, "hmsc_joint_200m_cell_partition.rds"
))

coda_objects <- Hmsc::convertToCodaObject(model)
summarise_coda <- function(object, family) {
  ess <- coda::effectiveSize(object)
  psrf <- tryCatch(
    coda::gelman.diag(object, multivariate = FALSE)$psrf,
    error = function(e) matrix(NA_real_, nrow = length(ess), ncol = 2)
  )
  tibble::tibble(
    family = family,
    parameter = names(ess),
    effective_sample_size = as.numeric(ess),
    psrf_point = as.numeric(psrf[, 1]),
    psrf_upper = as.numeric(psrf[, 2])
  )
}
convergence <- dplyr::bind_rows(
  summarise_coda(coda_objects$Beta, "Beta"),
  dplyr::bind_rows(lapply(seq_along(coda_objects$Omega), function(i) {
    summarise_coda(coda_objects$Omega[[i]], paste0("Omega_random_level_", i))
  }))
)

fold_balance <- dplyr::bind_cols(
  tibble::tibble(
    fold = as.integer(partition),
    cell = as.character(model$studyDesign$cell)
  ),
  as.data.frame(model$Y)
) |>
  dplyr::group_by(.data$fold) |>
  dplyr::summarise(
    rows = dplyr::n(),
    cells = dplyr::n_distinct(.data$cell),
    dplyr::across(
      dplyr::all_of(response_names),
      list(positives = sum, prevalence = mean),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

omega_tables <- lapply(seq_along(level_names), function(r) {
  omega <- Hmsc::getPostEstimate(model, parName = "OmegaCor", r = r)
  grid <- expand.grid(
    response_1 = response_names,
    response_2 = response_names,
    stringsAsFactors = FALSE
  )
  grid |>
    tibble::as_tibble() |>
    dplyr::mutate(
      random_level = level_names[r],
      residual_correlation = as.numeric(omega$mean),
      posterior_support_positive = as.numeric(omega$support),
      posterior_support_negative = as.numeric(omega$supportNeg),
      response_1_index = match(.data$response_1, response_names),
      response_2_index = match(.data$response_2, response_names)
    ) |>
    dplyr::filter(.data$response_1_index < .data$response_2_index) |>
    dplyr::select(-dplyr::ends_with("_index"))
}) |>
  dplyr::bind_rows()

group <- c(
  1L,
  2L,
  rep(3L, 6L),
  rep(4L, 2L),
  rep(5L, 3L)
)
group_names <- c(
  "Intercept", "Year trend", "Multiscale habitat",
  "Cull effort", "Manta effort and coral visibility"
)
variance <- Hmsc::computeVariancePartitioning(
  model, group = group, groupnames = group_names
)
variance_matrix <- as.matrix(variance$vals)
variance_table <- expand.grid(
  component = rownames(variance_matrix),
  response = colnames(variance_matrix),
  stringsAsFactors = FALSE
) |>
  tibble::as_tibble() |>
  dplyr::mutate(variance_fraction = as.numeric(variance_matrix))

utils::write.csv(
  omega_tables,
  file.path(cfg$paths$outputs_dir, "hmsc_joint_200m_residual_correlations.csv"),
  row.names = FALSE
)
utils::write.csv(
  variance_table,
  file.path(cfg$paths$outputs_dir, "hmsc_joint_200m_variance_partition.csv"),
  row.names = FALSE
)
utils::write.csv(
  fold_balance,
  file.path(cfg$paths$outputs_dir, "hmsc_joint_200m_fold_balance.csv"),
  row.names = FALSE
)
utils::write.csv(
  convergence,
  file.path(cfg$paths$outputs_dir, "hmsc_joint_200m_convergence_live.csv"),
  row.names = FALSE
)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  fig_dir <- file.path(v2_dir, "reports", "figures")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  plot_data <- variance_table |>
    dplyr::mutate(
      response = dplyr::recode(
        .data$response,
        cull_direct_presence = "Cull direct",
        cull_high_density = "Cull >0.02",
        manta_direct_presence = "Manta direct",
        manta_scar_presence = "Manta scars"
      )
    )
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(.data$response, .data$variance_fraction, fill = .data$component)
  ) +
    ggplot2::geom_col(width = 0.72) +
    ggplot2::scale_y_continuous(labels = scales::label_percent()) +
    ggplot2::labs(
      x = NULL, y = "Modelled variance fraction", fill = NULL,
      title = "Hmsc variance partition by evidence type"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      legend.position = "right"
    )
  ggplot2::ggsave(
    file.path(fig_dir, "hmsc_joint_variance_partition.png"), p,
    width = 9, height = 5.2, dpi = 180
  )
}

message("Joint Hmsc posterior summaries complete.")
print(omega_tables)
print(variance_table)
print(fold_balance)
print(dplyr::summarise(
  dplyr::group_by(convergence, .data$family),
  minimum_ess = min(.data$effective_sample_size, na.rm = TRUE),
  maximum_psrf_upper = max(.data$psrf_upper, na.rm = TRUE)
))
