args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v2_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v2_dir, "R", "bootstrap.R"))
ctx <- initialise_v2(script_file)
cfg <- ctx$config

required <- c(
  file.path(cfg$paths$artifacts_dir, "site_year_all.rds"),
  file.path(cfg$paths$artifacts_dir, "site_multiscale_30_60_120m.rds"),
  file.path(cfg$paths$artifacts_dir, "site_habitat_all_sites.rds"),
  file.path(
    cfg$paths$repo_root, "DataProcessing", "outputs",
    "fitted_model_clean_20260820_1041.rds"
  )
)
if (any(!file.exists(required))) stop("Missing two-part model input(s).")

site_year <- readRDS(required[[1]])
site_scales <- readRDS(required[[2]])
site_reference <- readRDS(required[[3]]) |>
  dplyr::filter(.data$evidence_source == "cull_site_history") |>
  dplyr::select(
    "site_id", "validation_group_id", "clean_complete"
  )
base_names <- names(terra::rast(cfg$paths$predictor_stack))
all_sets <- multiscale_feature_sets(base_names)

candidates <- list(
  clean_30m_annual = list(
    features = setdiff(all_sets$clean_30m_year, "year"),
    target_mode = "annual"
  ),
  clean_120m_annual = list(
    features = setdiff(all_sets$clean_120m_year, "year"),
    target_mode = "annual"
  ),
  multiscale_topography_annual = list(
    features = all_sets$multiscale_topography_static,
    target_mode = "annual"
  ),
  multiscale_topography_propensity = list(
    features = all_sets$multiscale_topography_static,
    target_mode = "site_propensity"
  )
)

data <- site_year |>
  dplyr::select(-dplyr::any_of(base_names)) |>
  dplyr::left_join(
    site_scales,
    by = c("site_id", "longitude", "latitude")
  ) |>
  dplyr::left_join(site_reference, by = "site_id") |>
  dplyr::mutate(
    reef_year_id = paste(.data$reef_name, .data$year, sep = "||"),
    outbreak_probability = .data$reef_outbreak_probability_current
  )

required_features <- unique(unlist(lapply(candidates, `[[`, "features")))
data <- data |>
  dplyr::filter(
    .data$clean_complete == 1L,
    is.finite(.data$outbreak_probability),
    .data$n_dives > 0,
    .data$bottom_time > 0,
    stats::complete.cases(dplyr::pick(
      dplyr::all_of(required_features)
    ))
  )

legacy <- readRDS(required[[4]])
parsnip <- workflows::extract_fit_parsnip(legacy$fit)
model_args <- parsnip$spec$args
model_cfg <- list(
  seed = cfg$model$seed,
  nrounds = as.integer(rlang::eval_tidy(model_args$trees)),
  eta = as.numeric(rlang::eval_tidy(model_args$learn_rate)),
  max_depth = as.integer(rlang::eval_tidy(model_args$tree_depth)),
  min_child_weight = as.numeric(rlang::eval_tidy(model_args$min_n)),
  gamma = as.numeric(rlang::eval_tidy(model_args$loss_reduction)),
  subsample = 1,
  colsample_bytree = 1,
  lambda = 1,
  alpha = 0,
  nthread = cfg$model$nthread
)

validation_seeds <- c(20260831L, 20260901L, 20260902L)
k <- 5L
site_frame <- data |>
  dplyr::distinct(
    .data$site_id, .data$reef_name, .data$validation_group_id
  )

fold_maps <- lapply(seq_along(validation_seeds), function(repeat_id) {
  assigned <- assign_within_reef_site_folds(
    site_frame, k = k, seed = validation_seeds[[repeat_id]]
  )
  dplyr::transmute(
    assigned,
    site_id = .data$site_id,
    validation_repeat = repeat_id,
    validation_seed = validation_seeds[[repeat_id]],
    fold = .data$site_fold
  )
})
fold_map <- dplyr::bind_rows(fold_maps)

evaluate_repeats <- function(predictions) {
  pieces <- split(predictions, interaction(
    predictions$model, predictions$validation_repeat, drop = TRUE
  ))
  metrics <- dplyr::bind_rows(lapply(pieces, function(x) {
    dplyr::bind_cols(
      dplyr::distinct(
        x, .data$model, .data$validation_repeat,
        .data$validation_seed
      ),
      binary_probability_metrics(x$truth, x$probability)
    )
  }))
  decomposition <- dplyr::bind_rows(lapply(pieces, function(x) {
    dplyr::bind_cols(
      dplyr::distinct(
        x, .data$validation_repeat, .data$validation_seed
      ),
      prediction_scale_decomposition(x, group = "reef_year_id")
    )
  }))
  capture <- dplyr::bind_rows(lapply(pieces, function(x) {
    dplyr::bind_cols(
      dplyr::distinct(
        x, .data$validation_repeat, .data$validation_seed
      ),
      within_group_capture_metrics(x, group = "reef_year_id")
    )
  }))
  list(
    metrics = metrics,
    decomposition = decomposition,
    capture = capture
  )
}

static_checkpoint_path <- file.path(
  cfg$paths$artifacts_dir,
  "two_part_static_screen_checkpoint.rds"
)
if (file.exists(static_checkpoint_path)) {
  checkpoint <- readRDS(static_checkpoint_path)
  static_predictions <- checkpoint$predictions
  static_evaluation <- checkpoint$evaluation
  static_summary <- checkpoint$summary
  selected_name <- checkpoint$selected_name
  message(
    "Loaded completed static screen checkpoint: ", selected_name
  )
} else {
message(
  "Running repeated blocked static-habitat screening on ",
  nrow(data), " site-years."
)
static_predictions <- list()
counter <- 0L
for (repeat_id in seq_along(validation_seeds)) {
  repeat_data <- data |>
    dplyr::left_join(
      dplyr::filter(
        fold_map, .data$validation_repeat == repeat_id
      ),
      by = "site_id"
    )
  for (candidate_name in names(candidates)) {
    candidate <- candidates[[candidate_name]]
    for (fold_id in seq_len(k)) {
      counter <- counter + 1L
      train <- dplyr::filter(repeat_data, .data$fold != fold_id)
      test <- dplyr::filter(repeat_data, .data$fold == fold_id)
      fit <- fit_static_habitat_candidate(
        train, candidate$features, candidate$target_mode,
        model_cfg, seed_offset = 60000L + counter
      )
      probability <- predict_xgb_binary(
        fit, test, candidate$features
      )
      static_predictions[[counter]] <- dplyr::transmute(
        test,
        model = candidate_name,
        validation = "repeated_held_out_site_within_reef",
        validation_repeat = repeat_id,
        validation_seed = validation_seeds[[repeat_id]],
        fold = fold_id,
        site_year_id = .data$site_year_id,
        site_id = .data$site_id,
        reef_name = .data$reef_name,
        reef_year_id = .data$reef_year_id,
        year = .data$year,
        truth = .data$high_density,
        probability = probability
      )
    }
  }
}
static_predictions <- dplyr::bind_rows(static_predictions)

static_evaluation <- suppressWarnings(
  evaluate_repeats(static_predictions)
)
static_summary <- static_evaluation$decomposition |>
  dplyr::select(
    "model", "validation_repeat", "within_reef_roc_auc",
    "overall_roc_auc"
  ) |>
  dplyr::left_join(
    static_evaluation$capture |>
      dplyr::filter(.data$target_fraction == 0.20) |>
      dplyr::select(
        "model", "validation_repeat", "capture", "lift"
      ),
    by = c("model", "validation_repeat")
  ) |>
  dplyr::group_by(.data$model) |>
  dplyr::summarise(
    repeats = dplyr::n(),
    mean_overall_auc = mean(.data$overall_roc_auc),
    mean_within_reef_year_auc = mean(.data$within_reef_roc_auc),
    sd_within_reef_year_auc = stats::sd(.data$within_reef_roc_auc),
    minimum_within_reef_year_auc = min(.data$within_reef_roc_auc),
    mean_top20_capture = mean(.data$capture),
    mean_top20_lift = mean(.data$lift),
    .groups = "drop"
  ) |>
  dplyr::arrange(
    dplyr::desc(.data$mean_within_reef_year_auc),
    dplyr::desc(.data$mean_top20_lift)
  )
selected_name <- static_summary$model[[1]]
selected <- candidates[[selected_name]]
message(
  "Selected static candidate for nested composition: ",
  selected_name
)
saveRDS(
  list(
    predictions = static_predictions,
    evaluation = static_evaluation,
    summary = static_summary,
    selected_name = selected_name
  ),
  file.path(
    cfg$paths$artifacts_dir,
    "two_part_static_screen_checkpoint.rds"
  )
)
utils::write.csv(
  static_summary,
  file.path(
    cfg$paths$outputs_dir,
    "two_part_static_candidate_summary.csv"
  ),
  row.names = FALSE
)
}
selected <- candidates[[selected_name]]

fixed_protocol <- list(n_dives = 2, bottom_time = 600)
nested_checkpoint_path <- file.path(
  cfg$paths$artifacts_dir,
  "two_part_nested_component_checkpoint.rds"
)
if (file.exists(nested_checkpoint_path)) {
  nested_checkpoint <- readRDS(nested_checkpoint_path)
  component_predictions <- nested_checkpoint$predictions
  component_coefficients <- nested_checkpoint$coefficients
  completed_outer_folds <- nested_checkpoint$completed_outer_folds
  message(
    "Loaded nested checkpoint with ",
    length(completed_outer_folds), " completed outer folds."
  )
} else {
  component_predictions <- list()
  component_coefficients <- list()
  completed_outer_folds <- character()
}
component_counter <- length(component_predictions)
for (repeat_id in seq_along(validation_seeds)) {
  repeat_data <- data |>
    dplyr::left_join(
      dplyr::filter(
        fold_map, .data$validation_repeat == repeat_id
      ),
      by = "site_id"
    )
  for (outer_fold in seq_len(k)) {
    outer_key <- paste(repeat_id, outer_fold, sep = "__")
    if (outer_key %in% completed_outer_folds) next
    outer_train <- dplyr::filter(
      repeat_data, .data$fold != outer_fold
    )
    outer_test <- dplyr::filter(
      repeat_data, .data$fold == outer_fold
    )
    outer_static <- static_predictions |>
      dplyr::filter(
        .data$model == selected_name,
        .data$validation_repeat == repeat_id,
        .data$fold == outer_fold
      ) |>
      dplyr::select(
        "site_year_id", habitat_probability = "probability"
      )
    outer_test <- dplyr::left_join(
      outer_test, outer_static, by = "site_year_id"
    )

    inner_predictions <- list()
    for (inner_fold in setdiff(seq_len(k), outer_fold)) {
      inner_fit_data <- dplyr::filter(
        outer_train, .data$fold != inner_fold
      )
      inner_test <- dplyr::filter(
        outer_train, .data$fold == inner_fold
      )
      inner_fit <- fit_static_habitat_candidate(
        inner_fit_data, selected$features, selected$target_mode,
        model_cfg,
        seed_offset = 70000L + repeat_id * 100L +
          outer_fold * 10L + inner_fold
      )
      inner_test$habitat_probability <- predict_xgb_binary(
        inner_fit, inner_test, selected$features
      )
      inner_predictions[[length(inner_predictions) + 1L]] <-
        inner_test
    }
    calibration <- dplyr::bind_rows(inner_predictions)
    bridge <- fit_rank_locked_state_bridge(
      calibration$high_density,
      calibration$habitat_probability,
      calibration$outbreak_probability
    )
    calibration$state_probability <-
      predict_rank_locked_state_bridge(
        bridge,
        calibration$habitat_probability,
        calibration$outbreak_probability
      )
    effort_fit <- fit_effort_standardiser(
      calibration$high_density,
      calibration$state_probability,
      calibration$n_dives,
      calibration$bottom_time
    )
    dive_fit <- fit_repeated_dive_detection(
      calibration$n_high_density_dives,
      calibration$n_dives,
      calibration$state_probability,
      calibration$bottom_time
    )

    outer_test$state_probability <-
      predict_rank_locked_state_bridge(
        bridge,
        outer_test$habitat_probability,
        outer_test$outbreak_probability
      )
    outer_test$fixed_effort_probability <-
      predict_effort_standardised(
        effort_fit, outer_test$state_probability,
        fixed_protocol$n_dives, fixed_protocol$bottom_time
      )
    outer_test$realised_effort_probability <-
      predict_effort_standardised(
        effort_fit, outer_test$state_probability,
        outer_test$n_dives, outer_test$bottom_time
      )
    outer_test$fixed_dive_detection_probability <-
      predict_fixed_protocol_detection(
        dive_fit, outer_test$state_probability,
        fixed_protocol$n_dives, fixed_protocol$bottom_time
      )

    common <- dplyr::transmute(
      outer_test,
      validation = "nested_repeated_held_out_site",
      validation_repeat = repeat_id,
      validation_seed = validation_seeds[[repeat_id]],
      fold = outer_fold,
      site_year_id = .data$site_year_id,
      site_id = .data$site_id,
      reef_name = .data$reef_name,
      reef_year_id = .data$reef_year_id,
      year = .data$year,
      truth = .data$high_density,
      habitat_probability = .data$habitat_probability,
      outbreak_probability = .data$outbreak_probability,
      n_dives = .data$n_dives,
      bottom_time = .data$bottom_time
    )
    values <- list(
      habitat_static = outer_test$habitat_probability,
      state_bridge = outer_test$state_probability,
      fixed_effort_pooled = outer_test$fixed_effort_probability,
      realised_effort_diagnostic =
        outer_test$realised_effort_probability,
      fixed_protocol_dive_detection =
        outer_test$fixed_dive_detection_probability
    )
    for (component in names(values)) {
      component_counter <- component_counter + 1L
      component_predictions[[component_counter]] <- dplyr::mutate(
        common,
        model = paste(
          selected_name, .env$component, sep = "__"
        ),
        component = .env$component,
        probability = .env$values[[.env$component]]
      )
    }
    component_coefficients[[
      length(component_coefficients) + 1L
    ]] <- tibble::tibble(
      validation_repeat = repeat_id,
      fold = outer_fold,
      component = c(
        rep("state_bridge", length(bridge$coefficients)),
        rep(
          "effort_standardiser",
          length(stats::coef(effort_fit))
        ),
        rep(
          "repeated_dive_detection",
          length(stats::coef(dive_fit))
        )
      ),
      term = c(
        names(bridge$coefficients),
        names(stats::coef(effort_fit)),
        names(stats::coef(dive_fit))
      ),
      estimate = c(
        unname(bridge$coefficients),
        unname(stats::coef(effort_fit)),
        unname(stats::coef(dive_fit))
      )
    )
    completed_outer_folds <- c(completed_outer_folds, outer_key)
    saveRDS(
      list(
        predictions = component_predictions,
        coefficients = component_coefficients,
        completed_outer_folds = completed_outer_folds
      ),
      nested_checkpoint_path
    )
    message(
      "Checkpointed nested outer fold ", outer_key,
      " (", length(completed_outer_folds), "/", k * length(validation_seeds),
      ")."
    )
  }
}
component_predictions <- dplyr::bind_rows(component_predictions)
rank_locked_components <- c(
  "habitat_static", "state_bridge", "fixed_effort_pooled",
  "fixed_protocol_dive_detection"
)
component_predictions <- component_predictions |>
  dplyr::group_by(
    .data$model, .data$validation_repeat, .data$reef_year_id
  ) |>
  dplyr::mutate(
    raw_component_probability = .data$probability,
    probability = dplyr::if_else(
      .data$component %in% rank_locked_components,
      rank_align_probability(
        .data$habitat_probability, .data$probability
      ),
      .data$probability
    )
  ) |>
  dplyr::ungroup()
component_evaluation <- suppressWarnings(
  evaluate_repeats(component_predictions)
)
component_summary <- component_evaluation$decomposition |>
  dplyr::select(
    "model", "validation_repeat", "within_reef_roc_auc",
    "overall_roc_auc"
  ) |>
  dplyr::left_join(
    component_evaluation$metrics |>
      dplyr::select(
        "model", "validation_repeat", "pr_auc", "brier",
        "log_loss"
      ),
    by = c("model", "validation_repeat")
  ) |>
  dplyr::left_join(
    component_evaluation$capture |>
      dplyr::filter(.data$target_fraction == 0.20) |>
      dplyr::select(
        "model", "validation_repeat", "capture", "lift"
      ),
    by = c("model", "validation_repeat")
  ) |>
  dplyr::group_by(.data$model) |>
  dplyr::summarise(
    repeats = dplyr::n(),
    mean_overall_auc = mean(.data$overall_roc_auc),
    mean_within_reef_year_auc = mean(
      .data$within_reef_roc_auc
    ),
    sd_within_reef_year_auc = stats::sd(
      .data$within_reef_roc_auc
    ),
    mean_pr_auc = mean(.data$pr_auc),
    mean_brier = mean(.data$brier),
    mean_log_loss = mean(.data$log_loss),
    mean_top20_capture = mean(.data$capture),
    mean_top20_lift = mean(.data$lift),
    .groups = "drop"
  )

rank_groups <- split(
  component_predictions,
  interaction(
    component_predictions$model,
    component_predictions$validation_repeat,
    drop = TRUE
  )
)
rank_checks <- rank_groups |>
  lapply(function(x) {
    dplyr::bind_cols(
      dplyr::distinct(
        x, .data$model, .data$validation_repeat
      ),
      rank_lock_check(x, group = "reef_year_id")
    )
  }) |>
  dplyr::bind_rows()

full_habitat_model <- fit_static_habitat_candidate(
  data, selected$features, selected$target_mode,
  model_cfg, seed_offset = 90000L
)
deployment_habitat <- component_predictions |>
  dplyr::filter(.data$component == "habitat_static")
deployment_bridge <- fit_rank_locked_state_bridge(
  deployment_habitat$truth,
  deployment_habitat$probability,
  deployment_habitat$outbreak_probability
)
deployment_rows <- component_predictions |>
  dplyr::filter(.data$component == "state_bridge")
deployment_effort <- fit_effort_standardiser(
  deployment_rows$truth,
  deployment_rows$probability,
  deployment_rows$n_dives,
  deployment_rows$bottom_time
)

model_bundle <- list(
  iteration = "two_part_fixed_effort_01",
  status = paste(
    "candidate_not_promoted_pending_outbreak_provenance",
    "and_expert_map_review",
    sep = "_"
  ),
  selected_static_candidate = selected_name,
  static_features = selected$features,
  static_target_mode = selected$target_mode,
  static_model = full_habitat_model,
  state_bridge = deployment_bridge,
  effort_standardiser = deployment_effort,
  fixed_protocol = fixed_protocol,
  outbreak_source = basename(cfg$paths$outbreak_hindcast_file),
  rank_locked_within_reef_year = TRUE,
  prediction_warning = paste(
    "Use fixed protocol values only; realised survey effort is",
    "diagnostic and must not be supplied to operational habitat maps."
  )
)

tables <- list(
  two_part_static_repeat_metrics = static_evaluation$metrics,
  two_part_static_repeat_decomposition =
    static_evaluation$decomposition,
  two_part_static_repeat_capture = static_evaluation$capture,
  two_part_static_candidate_summary = static_summary,
  two_part_component_repeat_metrics = component_evaluation$metrics,
  two_part_component_repeat_decomposition =
    component_evaluation$decomposition,
  two_part_component_repeat_capture = component_evaluation$capture,
  two_part_component_summary = component_summary,
  two_part_component_coefficients =
    dplyr::bind_rows(component_coefficients),
  two_part_rank_lock_checks = rank_checks,
  two_part_validation_fold_map = fold_map
)
for (name in names(tables)) {
  utils::write.csv(
    tables[[name]],
    file.path(cfg$paths$outputs_dir, paste0(name, ".csv")),
    row.names = FALSE
  )
}
saveRDS(
  static_predictions,
  file.path(
    cfg$paths$outputs_dir, "two_part_static_predictions.rds"
  )
)
saveRDS(
  component_predictions,
  file.path(
    cfg$paths$outputs_dir, "two_part_component_predictions.rds"
  )
)
saveRDS(
  model_bundle,
  file.path(
    cfg$paths$models_dir,
    "two_part_fixed_effort_iteration_01.rds"
  )
)
saveRDS(
  list(
    rows = nrow(data),
    sites = dplyr::n_distinct(data$site_id),
    reefs = dplyr::n_distinct(data$reef_name),
    years = range(data$year),
    validation_seeds = validation_seeds,
    fixed_protocol = fixed_protocol,
    selected_candidate = selected_name,
    model_cfg = model_cfg
  ),
  file.path(
    cfg$paths$artifacts_dir, "two_part_iteration_01_audit.rds"
  )
)

message("Two-part fixed-effort iteration complete.")
print(static_summary)
print(component_summary)
