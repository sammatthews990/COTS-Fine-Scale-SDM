args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v2_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v2_dir, "R", "bootstrap.R"))
ctx <- initialise_v2(script_file)
cfg <- ctx$config

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(Hmsc)
  library(patchwork)
  library(png)
  library(tidymodels)
  library(xgboost)
})

repo_root <- cfg$paths$repo_root
fig_dir <- file.path(v2_dir, "reports", "figures", "final_synthesis")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- cfg$paths$outputs_dir

v4_dir <- file.path(v2_dir, "iterations", "hmsc_corrected_v4")
v5_dir <- file.path(v2_dir, "iterations", "occupancy_detection_v5")
v6_dir <- file.path(v2_dir, "iterations", "density_rosetta_v6")

required_files <- c(
  legacy_model = file.path(repo_root, "DataProcessing", "outputs", "fitted_model_clean_20260820_1041.rds"),
  legacy_predictions = file.path(out_dir, "legacy_fixed_rank_predictions.rds"),
  legacy_surface = cfg$paths$legacy_2026_raster,
  cull_observations = file.path(cfg$paths$artifacts_dir, "cull_observations_clean.rds"),
  predictors = cfg$paths$predictor_stack,
  hmsc_model = file.path(v4_dir, "models", "enhanced_joint_120m.rds"),
  hmsc_data = file.path(v4_dir, "artifacts", "enhanced_joint_fit_data_120m.rds"),
  density_model = file.path(v6_dir, "models", "salad_native_60m.rds"),
  density_data = file.path(v6_dir, "artifacts", "density_data_60m.rds"),
  density_metrics = file.path(v6_dir, "outputs", "density_model_metrics.csv")
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing final synthesis input(s): ", paste(missing_files, collapse = "; "))
}

theme_final <- function() {
  theme_minimal(base_size = 11) +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, colour = "grey30"),
      plot.caption = element_text(size = 8.5, colour = "grey35", hjust = 0),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

read_one <- function(path) {
  if (!file.exists(path)) stop("Required metric file is missing: ", path)
  utils::read.csv(path, check.names = FALSE)
}

# -------------------------------------------------------------------------
# Figure 1: classification validation.
# -------------------------------------------------------------------------

legacy_predictions <- readRDS(required_files[["legacy_predictions"]])
legacy_year <- legacy_predictions |>
  filter(.data$model == "legacy_fixed_clean_year")
legacy_decomp <- prediction_scale_decomposition(legacy_year, group = "reef_year_id")

v2_full <- read_one(file.path(out_dir, "hmsc_v2_blocked_decomposition.csv"))
v2_reduced <- read_one(file.path(out_dir, "hmsc_v2_geomorphic_g360_blocked_decomposition.csv"))

metric_row <- function(model, validation, auc, status) {
  tibble(model = model, validation = validation, within_auc = auc, status = status)
}

v4_metric <- function(stem, validation, model, status) {
  z <- read_one(file.path(v4_dir, "outputs", paste0(stem, "_metrics.csv")))
  metric_row(model, validation, z$within_reef_year_auc[[1]], status)
}

v5_metric <- function(stem, validation) {
  z <- read_one(file.path(v5_dir, "outputs", paste0(stem, "_120m_metrics.csv")))
  metric_row("Shared-state v5 120 m", validation, z$within_reef_year_auc[[1]], "Incomplete")
}

classification_validation <- bind_rows(
  metric_row(
    "Legacy clean + year", "Historical held site",
    legacy_decomp$within_reef_roc_auc[[1]], "Accepted incumbent"
  ),
  v2_full |>
    transmute(
      model = "Hmsc v2 full 30-360 m",
      validation = recode(
        .data$validation,
        held_cell_3fold = "Held cell",
        leave_reef_3fold = "Leave reef",
        forward_year_train_through_2023 = "Forward year"
      ),
      within_auc = .data$within_reef_roc_auc,
      status = "Convergence failed"
    ),
  v2_reduced |>
    transmute(
      model = "Hmsc v2 reduced 360 m",
      validation = recode(
        .data$validation,
        held_cell_3fold = "Held cell",
        leave_reef_3fold = "Leave reef",
        forward_year_train_through_2023 = "Forward year"
      ),
      within_auc = .data$within_reef_roc_auc,
      status = "Converged"
    ),
  v4_metric(
    "enhanced_joint_held_cell_120m", "Held cell",
    "Hmsc v4 enhanced 120 m", "Convergence failed"
  ),
  v4_metric(
    "enhanced_joint_leave_reef_120m", "Leave reef",
    "Hmsc v4 enhanced 120 m", "Convergence failed"
  ),
  v4_metric(
    "enhanced_joint_forward_120m", "Forward year",
    "Hmsc v4 enhanced 120 m", "Convergence failed"
  ),
  v4_metric(
    "enhanced_cots_60support_pca1_held_cell_60m", "Held cell",
    "Hmsc v4 enhanced 60 m", "Convergence failed"
  ),
  v4_metric(
    "enhanced_cots_60support_pca1_leave_reef_60m", "Leave reef",
    "Hmsc v4 enhanced 60 m", "Convergence failed"
  ),
  v4_metric(
    "enhanced_cots_60support_pca1_forward_60m", "Forward year",
    "Hmsc v4 enhanced 60 m", "Convergence failed"
  ),
  v5_metric("shared_state_held_cell", "Held cell"),
  v5_metric("shared_state_leave_reef", "Leave reef"),
  v5_metric("shared_state_forward", "Forward year")
)

model_levels <- c(
  "Legacy clean + year",
  "Hmsc v2 full 30-360 m",
  "Hmsc v2 reduced 360 m",
  "Hmsc v4 enhanced 120 m",
  "Hmsc v4 enhanced 60 m",
  "Shared-state v5 120 m"
)
classification_validation <- classification_validation |>
  mutate(
    model = factor(.data$model, levels = rev(model_levels)),
    validation = factor(
      .data$validation,
      levels = c("Historical held site", "Held cell", "Leave reef", "Forward year")
    )
  )

p_class <- ggplot(
  classification_validation,
  aes(x = .data$within_auc, y = .data$model, colour = .data$validation, shape = .data$status)
) +
  geom_vline(xintercept = 0.5, colour = "grey55", linewidth = 0.5) +
  geom_point(size = 3, stroke = 1.1) +
  geom_text(
    aes(label = sprintf("%.3f", .data$within_auc)),
    nudge_x = 0.006, hjust = 0, size = 3, show.legend = FALSE
  ) +
  scale_colour_manual(values = c(
    "Historical held site" = "#7A7A7A",
    "Held cell" = "#0072B2",
    "Leave reef" = "#D55E00",
    "Forward year" = "#009E73"
  )) +
  scale_shape_manual(values = c(
    "Accepted incumbent" = 18,
    "Converged" = 16,
    "Convergence failed" = 1,
    "Incomplete" = 4
  )) +
  coord_cartesian(xlim = c(0.49, 0.63), clip = "off") +
  labs(
    title = "Blocked within-reef discrimination of classification candidates",
    subtitle = "AUC values separate site-scale ordering from reef-level prevalence",
    x = "Within-reef-year ROC AUC",
    y = NULL,
    colour = "Validation",
    shape = "Model status",
    caption = paste(
      "Legacy used its historical held-site protocol and is not directly matched to the modern folds.",
      "Hollow/crossed symbols identify candidates without a promotion-eligible converged fit."
    )
  ) +
  theme_final() +
  theme(plot.margin = margin(5.5, 35, 5.5, 5.5))

fig1 <- file.path(fig_dir, "final_fig01_classification_validation.png")
ggsave(fig1, p_class, width = 11.5, height = 6.8, dpi = 300, bg = "white")
write.csv(
  classification_validation |>
    mutate(model = as.character(.data$model), validation = as.character(.data$validation)),
  file.path(out_dir, "final_synthesis_classification_validation.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# Figure 2: density validation. Keep rank correlation separate from AUC.
# -------------------------------------------------------------------------

density_validation <- read_one(required_files[["density_metrics"]]) |>
  filter(.data$support_m == 60, .data$validation != "apparent") |>
  mutate(
    series = case_when(
      .data$model == "cull_native" ~ "Cull-native: Cull CPUE",
      .data$model == "salad_native" ~ "SALAD-native: SALAD density",
      .data$model == "joint_rosetta" & .data$method == "cull_cpue" ~
        "Rosetta joint: Cull CPUE",
      .data$model == "joint_rosetta" & .data$method == "salad_density_ha" ~
        "Rosetta joint: SALAD density",
      TRUE ~ paste(.data$model, .data$method)
    ),
    validation_label = recode(
      .data$validation,
      held_cell = "Held cell",
      leave_reef = "Leave reef",
      forward_year = "Forward year"
    ),
    validation_label = factor(
      .data$validation_label,
      levels = c("Held cell", "Leave reef", "Forward year")
    )
  )

p_density <- ggplot(
  density_validation,
  aes(
    x = .data$validation_label,
    y = .data$within_reef_year_spearman,
    group = .data$series,
    colour = .data$series,
    shape = .data$series
  )
) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.5) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.8) +
  geom_text(
    aes(label = sprintf("%.2f", .data$within_reef_year_spearman)),
    nudge_y = 0.018, size = 2.8, show.legend = FALSE
  ) +
  scale_colour_manual(values = c(
    "Cull-native: Cull CPUE" = "#0072B2",
    "SALAD-native: SALAD density" = "#009E73",
    "Rosetta joint: Cull CPUE" = "#E69F00",
    "Rosetta joint: SALAD density" = "#CC79A7"
  )) +
  scale_shape_manual(values = c(16, 17, 15, 18)) +
  coord_cartesian(ylim = c(-0.14, 0.25), clip = "off") +
  labs(
    title = "Continuous-density models do not yet transfer stable within-reef rank",
    subtitle = "Spearman correlation is shown separately from classification AUC",
    x = NULL,
    y = "Within-reef-year Spearman correlation",
    colour = NULL,
    shape = NULL,
    caption = paste(
      "SALAD-native has a positive forward result but does not reproduce it under held-cell validation.",
      "The joint calibration degrades the SALAD signal."
    )
  ) +
  theme_final()

fig2 <- file.path(fig_dir, "final_fig02_density_validation.png")
ggsave(fig2, p_density, width = 10.5, height = 6.5, dpi = 300, bg = "white")
write.csv(
  density_validation,
  file.path(out_dir, "final_synthesis_density_validation.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# Figure 3: model-specific variable importance/effect summaries.
# -------------------------------------------------------------------------

legacy <- readRDS(required_files[["legacy_model"]])
legacy_engine <- workflows::extract_fit_engine(legacy$fit)
legacy_recipe <- workflows::extract_recipe(legacy$fit, estimated = TRUE)
legacy_importance <- xgboost::xgb.importance(model = legacy_engine) |>
  as_tibble() |>
  transmute(
    model = "Legacy clean + year",
    term = .data$Feature,
    value = .data$Gain,
    direction = "Unsigned gain",
    measure = "XGBoost gain"
  )

hmsc_model <- readRDS(required_files[["hmsc_model"]])
hmsc_artifact <- readRDS(required_files[["hmsc_data"]])
target_index <- match("cull_high_density", hmsc_model$spNames)
hmsc_beta <- Hmsc::getPostEstimate(hmsc_model, parName = "Beta")
hmsc_effect <- as.numeric(hmsc_beta$mean[, target_index]) *
  as.numeric(hmsc_model$XScalePar[2, ])
hmsc_importance <- tibble(
  model = "Hmsc v4 enhanced 120 m",
  term = hmsc_model$covNames,
  value = hmsc_effect
) |>
  filter(.data$term != "(Intercept)") |>
  slice_max(abs(.data$value), n = 12, with_ties = FALSE) |>
  mutate(
    direction = if_else(.data$value >= 0, "Positive association", "Negative association"),
    measure = "Standardised probit coefficient"
  )

density_model <- readRDS(required_files[["density_model"]])
density_artifact <- readRDS(required_files[["density_data"]])
density_importance <- tibble(
  model = "SALAD-native density 60 m",
  term = density_artifact$design$columns,
  value = as.numeric(density_model$beta)
) |>
  filter(.data$term != "(Intercept)") |>
  slice_max(abs(.data$value), n = 12, with_ties = FALSE) |>
  mutate(
    direction = if_else(.data$value >= 0, "Positive association", "Negative association"),
    measure = "Standardised log-density coefficient"
  )

pretty_term <- function(x) {
  x <- sub("^geomorphic_class", "Geomorph: ", x)
  x <- sub("_pc1$", " multiscale PC1", x)
  x <- gsub("_", " ", x)
  recode(
    x,
    "Year" = "Survey year",
    "DEM" = "Depth (DEM)",
    "SLO" = "Slope",
    "SVF" = "Sky-view factor",
    "VRM" = "Vector ruggedness",
    "BPI" = "Bathymetric position",
    "EAST" = "Eastness",
    "NORTH" = "Northness",
    "HCU" = "Horizontal curvature",
    "VCU" = "Vertical curvature",
    "AhyaD" = "Acropora hyacinthus",
    "Aspat" = "Acropora spathulata",
    "Aten" = "Acropora tenuis",
    "waves Hs" = "Significant wave height",
    "waves Tp" = "Peak wave period"
  )
}

legacy_importance$label <- pretty_term(legacy_importance$term)
hmsc_importance$label <- pretty_term(hmsc_importance$term)
density_importance$label <- pretty_term(density_importance$term)

importance_panel <- function(data, title, x_title, unsigned = FALSE) {
  data <- data |>
    arrange(.data$value) |>
    mutate(label = factor(.data$label, levels = .data$label))
  p <- ggplot(data, aes(x = .data$value, y = .data$label))
  if (unsigned) {
    p <- p +
      geom_col(fill = "#0072B2", width = 0.72)
  } else {
    p <- p +
      geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.45) +
      geom_col(aes(fill = .data$direction), width = 0.72) +
      scale_fill_manual(values = c(
        "Negative association" = "#D55E00",
        "Positive association" = "#0072B2"
      ), guide = "none")
  }
  p +
    labs(title = title, x = x_title, y = NULL) +
    theme_final() +
    theme(
      plot.title = element_text(size = 11.5),
      axis.text.y = element_text(size = 8.5),
      legend.position = "none"
    )
}

p_imp_legacy <- importance_panel(
  legacy_importance,
  "A. Legacy clean + year",
  "XGBoost gain",
  unsigned = TRUE
)
p_imp_hmsc <- importance_panel(
  hmsc_importance,
  "B. Hmsc v4 target response",
  "Standardised probit coefficient"
)
p_imp_density <- importance_panel(
  density_importance,
  "C. SALAD-native density",
  "Standardised log-density coefficient"
)

p_importance <- p_imp_legacy + p_imp_hmsc + p_imp_density +
  plot_layout(widths = c(0.9, 1.15, 1.15)) +
  plot_annotation(
    title = "Variable influence differs among the operational, habitat and density models",
    subtitle = paste(
      "Legacy gain is unsigned; Hmsc and density panels show signed associations.",
      "Hmsc effects are diagnostic because the v4 posterior did not pass convergence."
    ),
    caption = paste(
      "Importance measures are model-specific and should not be compared by bar length across panels.",
      "Geomorphic classes are contrasts against the fitted reference class."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10),
      plot.caption = element_text(size = 8.5, colour = "grey35", hjust = 0)
    )
  )

fig3 <- file.path(fig_dir, "final_fig03_variable_importance.png")
ggsave(fig3, p_importance, width = 15, height = 8.2, dpi = 300, bg = "white")
write.csv(
  bind_rows(legacy_importance, hmsc_importance, density_importance),
  file.path(out_dir, "final_synthesis_variable_importance.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# Figure 4: marginal partial dependence, centred within each model.
# -------------------------------------------------------------------------

set.seed(20260914)
predictor_stack <- terra::rast(required_files[["predictors"]])
legacy_surface <- terra::rast(required_files[["legacy_surface"]])
legacy_sample_locations <- terra::spatSample(
  legacy_surface,
  size = 1800,
  method = "random",
  na.rm = TRUE,
  as.points = TRUE,
  values = FALSE,
  exhaust = TRUE
)
legacy_points <- terra::extract(
  predictor_stack, legacy_sample_locations, ID = FALSE
) |>
  as.data.frame()
legacy_points$Year <- legacy$predict_year
legacy_features <- legacy_engine$feature_names
legacy_points <- legacy_points[
  complete.cases(legacy_points[, legacy_features, drop = FALSE]),
  legacy_features,
  drop = FALSE
]
if (nrow(legacy_points) > 800) {
  legacy_points <- legacy_points[sample(seq_len(nrow(legacy_points)), 800), , drop = FALSE]
}

legacy_predict <- function(new_data) {
  baked <- recipes::bake(
    legacy_recipe,
    new_data = as.data.frame(new_data),
    recipes::all_predictors()
  )
  stats::predict(
    legacy_engine,
    newdata = as.matrix(baked[, legacy_features, drop = FALSE])
  )
}
legacy_baseline <- mean(legacy_predict(legacy_points))
legacy_year_range <- range(
  readRDS(required_files[["cull_observations"]])$year,
  na.rm = TRUE
)
legacy_pdp_vars <- c("Year", "BPI", "DEM", "AhyaD")
legacy_labels <- c(
  Year = "Legacy: survey year",
  BPI = "Legacy: bathymetric position",
  DEM = "Legacy: depth (DEM)",
  AhyaD = "Legacy: A. hyacinthus suitability"
)

legacy_pdp <- bind_rows(lapply(legacy_pdp_vars, function(variable) {
  grid <- if (variable == "Year") {
    seq(legacy_year_range[[1]], min(legacy_year_range[[2]], legacy$predict_year))
  } else {
    as.numeric(quantile(
      legacy_points[[variable]], probs = seq(0.02, 0.98, length.out = 35),
      na.rm = TRUE, names = FALSE
    ))
  }
  grid <- unique(grid[is.finite(grid)])
  groups <- rep(seq_along(grid), each = nrow(legacy_points))
  changed <- legacy_points[
    rep(seq_len(nrow(legacy_points)), times = length(grid)),
    , drop = FALSE
  ]
  changed[[variable]] <- rep(grid, each = nrow(legacy_points))
  values <- as.numeric(rowsum(legacy_predict(changed), groups) / nrow(legacy_points))
  tibble(
    model = "Legacy clean + year",
    variable = variable,
    panel = legacy_labels[[variable]],
    x = grid,
    probability = values,
    lower = values,
    upper = values,
    baseline = legacy_baseline,
    posterior_interval = FALSE
  )
}))

hmsc_xdata <- as.data.frame(hmsc_artifact$XData)
if (nrow(hmsc_xdata) > 500) {
  hmsc_xdata <- hmsc_xdata[sample(seq_len(nrow(hmsc_xdata)), 500), , drop = FALSE]
}
hmsc_draws <- unlist(hmsc_model$postList, recursive = FALSE)
draw_index <- unique(round(seq(1, length(hmsc_draws), length.out = min(100, length(hmsc_draws)))))
hmsc_beta_draws <- vapply(
  hmsc_draws[draw_index],
  function(draw) draw$Beta[, target_index],
  numeric(length(hmsc_model$covNames))
)
hmsc_beta_mean <- rowMeans(hmsc_beta_draws)

hmsc_matrix <- function(data) {
  X <- model.matrix(hmsc_artifact$XFormula, data = data)
  missing <- setdiff(hmsc_model$covNames, colnames(X))
  if (length(missing)) {
    X <- cbind(X, matrix(
      0, nrow = nrow(X), ncol = length(missing),
      dimnames = list(NULL, missing)
    ))
  }
  X[, hmsc_model$covNames, drop = FALSE]
}

hmsc_baseline_matrix <- hmsc_matrix(hmsc_xdata)
hmsc_baseline <- mean(pnorm(as.numeric(hmsc_baseline_matrix %*% hmsc_beta_mean)))
hmsc_pdp_vars <- c("sd_pc1", "vcu_pc1", "AhyaD", "waves_Hs")
hmsc_labels <- c(
  sd_pc1 = "Hmsc: local variability PC1",
  vcu_pc1 = "Hmsc: vertical curvature PC1",
  AhyaD = "Hmsc: A. hyacinthus suitability",
  waves_Hs = "Hmsc: significant wave height"
)

hmsc_pdp <- bind_rows(lapply(hmsc_pdp_vars, function(variable) {
  grid <- as.numeric(quantile(
    hmsc_xdata[[variable]], probs = seq(0.02, 0.98, length.out = 35),
    na.rm = TRUE, names = FALSE
  ))
  grid <- unique(grid[is.finite(grid)])
  groups <- rep(seq_along(grid), each = nrow(hmsc_xdata))
  changed <- hmsc_xdata[
    rep(seq_len(nrow(hmsc_xdata)), times = length(grid)),
    , drop = FALSE
  ]
  changed[[variable]] <- rep(grid, each = nrow(hmsc_xdata))
  probability_matrix <- pnorm(hmsc_matrix(changed) %*% hmsc_beta_draws)
  grid_draw_probability <- rowsum(probability_matrix, groups) / nrow(hmsc_xdata)
  rm(probability_matrix)
  tibble(
    model = "Hmsc v4 enhanced 120 m",
    variable = variable,
    panel = hmsc_labels[[variable]],
    x = grid,
    probability = rowMeans(grid_draw_probability),
    lower = apply(grid_draw_probability, 1, quantile, 0.1, names = FALSE),
    upper = apply(grid_draw_probability, 1, quantile, 0.9, names = FALSE),
    baseline = hmsc_baseline,
    posterior_interval = TRUE
  )
}))

partial_dependence <- bind_rows(legacy_pdp, hmsc_pdp) |>
  mutate(
    delta_pp = 100 * (.data$probability - .data$baseline),
    lower_delta_pp = 100 * (.data$lower - .data$baseline),
    upper_delta_pp = 100 * (.data$upper - .data$baseline),
    panel = factor(.data$panel, levels = c(unname(legacy_labels), unname(hmsc_labels)))
  )

p_pdp <- ggplot(
  partial_dependence,
  aes(x = .data$x, y = .data$delta_pp, colour = .data$model, fill = .data$model)
) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.45) +
  geom_ribbon(
    data = filter(partial_dependence, .data$posterior_interval),
    aes(ymin = .data$lower_delta_pp, ymax = .data$upper_delta_pp),
    alpha = 0.18, colour = NA
  ) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ panel, scales = "free_x", ncol = 4) +
  scale_colour_manual(values = c(
    "Legacy clean + year" = "#0072B2",
    "Hmsc v4 enhanced 120 m" = "#D55E00"
  ), guide = "none") +
  scale_fill_manual(values = c(
    "Legacy clean + year" = "#0072B2",
    "Hmsc v4 enhanced 120 m" = "#D55E00"
  ), guide = "none") +
  labs(
    title = "Partial dependence reveals different sources of spatial and temporal contrast",
    subtitle = "Curves show marginal probability change relative to each model's baseline",
    x = "Predictor value",
    y = "Change in predicted probability (percentage points)",
    caption = paste(
      "Legacy curves use the 2026 GBR predictor distribution; its year curve varies year while holding habitat fixed.",
      "Hmsc ribbons are 10th-90th posterior intervals but remain diagnostic because v4 convergence failed.",
      "Multiscale PC axes are unitless and their sign is specific to the fitted PCA."
    )
  ) +
  theme_final() +
  theme(
    strip.text = element_text(size = 8.5),
    axis.text = element_text(size = 8),
    legend.position = "none"
  )

fig4 <- file.path(fig_dir, "final_fig04_partial_dependence.png")
ggsave(fig4, p_pdp, width = 14, height = 7.8, dpi = 300, bg = "white")
write.csv(
  partial_dependence,
  file.path(out_dir, "final_synthesis_partial_dependence.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# Figure 5: four-reef expert map synthesis.
# -------------------------------------------------------------------------

reef_slugs <- c("eyrie", "batt", "broadhurst", "bait")
map_files <- file.path(
  v4_dir, "reports", "figures",
  paste0("v4_support_comparison_", reef_slugs, "_60m_120m.png")
)
if (any(!file.exists(map_files))) {
  stop("Missing four-reef synthesis map(s): ", paste(map_files[!file.exists(map_files)], collapse = "; "))
}

fig5 <- file.path(fig_dir, "final_fig05_four_reef_maps.png")
grDevices::png(fig5, width = 4200, height = 3200, res = 260, bg = "white")
grid::grid.newpage()
grid::pushViewport(grid::viewport(
  layout = grid::grid.layout(
    nrow = 4, ncol = 4,
    widths = grid::unit(rep(1, 4), "null"),
    heights = grid::unit(rep(1, 4), "null")
  )
))
for (i in seq_along(map_files)) {
  map_image <- png::readPNG(map_files[[i]])
  image_height <- dim(map_image)[[1]]
  image_width <- dim(map_image)[[2]]
  row_break <- floor(image_height / 2)
  column_breaks <- round(seq(1, image_width + 1, length.out = 4))
  panels <- list(
    map_image[seq_len(row_break), column_breaks[[1]]:(column_breaks[[2]] - 1), , drop = FALSE],
    map_image[seq_len(row_break), column_breaks[[2]]:(column_breaks[[3]] - 1), , drop = FALSE],
    map_image[seq_len(row_break), column_breaks[[3]]:(column_breaks[[4]] - 1), , drop = FALSE],
    map_image[(row_break + 1):image_height, column_breaks[[1]]:(column_breaks[[2]] - 1), , drop = FALSE]
  )
  for (column in seq_along(panels)) {
    grid::grid.raster(
      panels[[column]],
      interpolate = TRUE,
      vp = grid::viewport(layout.pos.row = i, layout.pos.col = column)
    )
  }
}
grid::popViewport()
grDevices::dev.off()

figure_manifest <- tibble(
  figure = 1:5,
  topic = c(
    "Classification validation",
    "Continuous-density validation",
    "Variable importance and signed associations",
    "Partial dependence",
    "Four-reef expert map comparison"
  ),
  file = file.path("figures", "final_synthesis", basename(c(fig1, fig2, fig3, fig4, fig5))),
  primary_models = c(
    "Legacy, Hmsc v2/v4, shared-state v5",
    "Cull-native, SALAD-native, Rosetta joint",
    "Legacy, Hmsc v4 enhanced 120 m, SALAD-native 60 m",
    "Legacy and Hmsc v4 enhanced 120 m",
    "Legacy, Hmsc v2 static joint, Hmsc v4 60 m and 120 m"
  ),
  acceptance_use = c(
    "Statistical context; protocols labelled",
    "Density branch validation; not comparable to AUC",
    "Interpretation only; model-specific scales",
    "Interpretation only; not causal",
    "Required named-zone expert review"
  )
)
write.csv(
  figure_manifest,
  file.path(out_dir, "final_synthesis_figure_manifest.csv"),
  row.names = FALSE
)

message("Final synthesis figures complete.")
print(figure_manifest)
