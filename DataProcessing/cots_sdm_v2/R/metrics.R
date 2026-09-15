safe_probability <- function(probability, epsilon = 1e-6) {
  pmin(pmax(as.numeric(probability), epsilon), 1 - epsilon)
}

capture_at_fraction <- function(truth, probability, fraction) {
  n_select <- max(1L, ceiling(length(truth) * fraction))
  selected <- order(probability, decreasing = TRUE)[seq_len(n_select)]
  total_positive <- sum(truth == 1L)
  if (total_positive == 0L) return(NA_real_)
  sum(truth[selected] == 1L) / total_positive
}

binary_probability_metrics <- function(truth, probability) {
  truth <- as.integer(truth)
  probability <- safe_probability(probability)
  truth_factor <- factor(as.character(truth), levels = c("1", "0"))

  calibration <- tryCatch(
    stats::coef(stats::glm(
      truth ~ stats::qlogis(probability),
      family = stats::binomial()
    )),
    error = function(e) c(`(Intercept)` = NA_real_, slope = NA_real_)
  )

  tibble::tibble(
    n = length(truth),
    positives = sum(truth == 1L),
    prevalence = mean(truth),
    roc_auc = tryCatch(
      yardstick::roc_auc_vec(
        truth_factor, probability, event_level = "first"
      ),
      error = function(e) NA_real_
    ),
    pr_auc = tryCatch(
      yardstick::pr_auc_vec(
        truth_factor, probability, event_level = "first"
      ),
      error = function(e) NA_real_
    ),
    brier = mean((truth - probability)^2),
    log_loss = -mean(
      truth * log(probability) + (1 - truth) * log(1 - probability)
    ),
    calibration_intercept = unname(calibration[1]),
    calibration_slope = unname(calibration[2]),
    capture_top_10pct = capture_at_fraction(truth, probability, 0.10),
    capture_top_20pct = capture_at_fraction(truth, probability, 0.20)
  )
}

within_group_auc <- function(
    data, group = "reef_name", truth = "truth",
    probability = "probability") {
  groups <- split(data, data[[group]], drop = TRUE)
  per_group <- lapply(names(groups), function(group_name) {
    x <- groups[[group_name]]
    positives <- sum(x[[truth]] == 1L)
    negatives <- sum(x[[truth]] == 0L)
    pairs <- positives * negatives
    auc <- if (pairs > 0L) {
      binary_probability_metrics(
        x[[truth]], x[[probability]]
      )$roc_auc
    } else {
      NA_real_
    }
    tibble::tibble(
      group = group_name,
      rows = nrow(x),
      positives = positives,
      negatives = negatives,
      positive_negative_pairs = pairs,
      roc_auc = auc
    )
  })
  per_group <- dplyr::bind_rows(per_group)
  mixed <- dplyr::filter(
    per_group,
    .data$positive_negative_pairs > 0L,
    is.finite(.data$roc_auc)
  )
  pooled_auc <- if (nrow(mixed) > 0L) {
    stats::weighted.mean(
      mixed$roc_auc,
      mixed$positive_negative_pairs
    )
  } else {
    NA_real_
  }
  list(
    summary = tibble::tibble(
      groups = nrow(per_group),
      mixed_groups = nrow(mixed),
      positive_negative_pairs = sum(
        mixed$positive_negative_pairs
      ),
      within_group_roc_auc = pooled_auc
    ),
    by_group = per_group
  )
}

prediction_scale_decomposition <- function(
    predictions, group = "reef_name") {
  models <- split(predictions, predictions$model, drop = TRUE)
  dplyr::bind_rows(lapply(names(models), function(model_name) {
    x <- models[[model_name]]
    overall <- binary_probability_metrics(x$truth, x$probability)
    within <- within_group_auc(x, group = group)$summary
    group_mean_probability <- ave(
      x$probability, x[[group]], FUN = mean
    )
    between_only <- binary_probability_metrics(
      x$truth, group_mean_probability
    )
    tibble::tibble(
      model = model_name,
      rows = nrow(x),
      reefs = dplyr::n_distinct(x[[group]]),
      mixed_reefs = within$mixed_groups,
      overall_roc_auc = overall$roc_auc,
      between_reef_only_roc_auc = between_only$roc_auc,
      within_reef_roc_auc = within$within_group_roc_auc,
      within_positive_negative_pairs =
        within$positive_negative_pairs,
      mean_within_reef_sd = mean(
        ave(x$probability, x[[group]], FUN = stats::sd),
        na.rm = TRUE
      )
    )
  }))
}

within_group_capture_metrics <- function(
    predictions, group = "reef_name",
    fractions = c(0.10, 0.20, 0.30),
    minimum_group_rows = 3L) {
  models <- split(predictions, predictions$model, drop = TRUE)
  dplyr::bind_rows(lapply(names(models), function(model_name) {
    x <- models[[model_name]]
    groups <- split(x, x[[group]], drop = TRUE)
    eligible <- groups[vapply(groups, function(z) {
      nrow(z) >= minimum_group_rows &&
        any(z$truth == 1L) &&
        any(z$truth == 0L)
    }, logical(1))]

    dplyr::bind_rows(lapply(fractions, function(fraction) {
      selected <- lapply(eligible, function(z) {
        n_select <- max(1L, ceiling(nrow(z) * fraction))
        index <- order(z$probability, decreasing = TRUE)[
          seq_len(n_select)
        ]
        tibble::tibble(
          rows = nrow(z),
          positives = sum(z$truth == 1L),
          selected = n_select,
          selected_positives = sum(z$truth[index] == 1L)
        )
      }) |>
        dplyr::bind_rows()
      if (nrow(selected) == 0L) {
        return(tibble::tibble(
          model = model_name,
          group = group,
          target_fraction = fraction,
          eligible_groups = 0L,
          eligible_rows = 0L,
          eligible_positives = 0L,
          selected_rows = 0L,
          captured_positives = 0L,
          capture = NA_real_,
          selected_prevalence = NA_real_,
          baseline_prevalence = NA_real_,
          lift = NA_real_
        ))
      }
      positives <- sum(selected$positives)
      selected_rows <- sum(selected$selected)
      captured <- sum(selected$selected_positives)
      baseline <- positives / sum(selected$rows)
      selected_prevalence <- captured / selected_rows
      tibble::tibble(
        model = model_name,
        group = group,
        target_fraction = fraction,
        eligible_groups = nrow(selected),
        eligible_rows = sum(selected$rows),
        eligible_positives = positives,
        selected_rows = selected_rows,
        captured_positives = captured,
        capture = captured / positives,
        selected_prevalence = selected_prevalence,
        baseline_prevalence = baseline,
        lift = selected_prevalence / baseline
      )
    }))
  }))
}

prediction_scale_summary <- function(predictions) {
  predictions |>
    dplyr::group_by(.data$model, .data$validation) |>
    dplyr::summarise(
      rows = dplyr::n(),
      unique_values = dplyr::n_distinct(.data$probability),
      minimum = min(.data$probability, na.rm = TRUE),
      q01 = stats::quantile(.data$probability, 0.01, na.rm = TRUE),
      q05 = stats::quantile(.data$probability, 0.05, na.rm = TRUE),
      median = stats::median(.data$probability, na.rm = TRUE),
      q95 = stats::quantile(.data$probability, 0.95, na.rm = TRUE),
      q99 = stats::quantile(.data$probability, 0.99, na.rm = TRUE),
      maximum = max(.data$probability, na.rm = TRUE),
      standard_deviation = stats::sd(.data$probability, na.rm = TRUE),
      q05_q95_span = .data$q95 - .data$q05,
      .groups = "drop"
    )
}

within_group_percentile_predictions <- function(
    predictions, group = "reef_name") {
  predictions |>
    dplyr::group_by(.data$model, .data[[group]]) |>
    dplyr::mutate(
      group_rows = dplyr::n(),
      probability = dplyr::if_else(
        .data$group_rows > 1L,
        (rank(.data$probability, ties.method = "average") - 1) /
          (.data$group_rows - 1),
        0.5
      )
    ) |>
    dplyr::ungroup()
}

within_group_percentile_contrast <- function(
    predictions, group = "reef_name",
    minimum_group_rows = 3L) {
  ranked <- within_group_percentile_predictions(
    predictions, group = group
  )
  eligible_groups <- ranked |>
    dplyr::group_by(.data$model, .data[[group]]) |>
    dplyr::summarise(
      rows = dplyr::n(),
      positives = sum(.data$truth == 1L),
      negatives = sum(.data$truth == 0L),
      .groups = "drop"
    ) |>
    dplyr::filter(
      .data$rows >= minimum_group_rows,
      .data$positives > 0L,
      .data$negatives > 0L
    ) |>
    dplyr::mutate(eligible = TRUE)

  ranked |>
    dplyr::left_join(
      dplyr::select(
        eligible_groups,
        "model", dplyr::all_of(group), "eligible"
      ),
      by = c("model", group)
    ) |>
    dplyr::filter(.data$eligible %in% TRUE) |>
    dplyr::group_by(.data$model) |>
    dplyr::summarise(
      group = group,
      eligible_groups = dplyr::n_distinct(.data[[group]]),
      rows = dplyr::n(),
      positives = sum(.data$truth == 1L),
      mean_positive_percentile = mean(
        .data$probability[.data$truth == 1L]
      ),
      mean_negative_percentile = mean(
        .data$probability[.data$truth == 0L]
      ),
      positive_negative_percentile_gap =
        .data$mean_positive_percentile -
        .data$mean_negative_percentile,
      .groups = "drop"
    )
}
