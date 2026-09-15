safe_component_logit <- function(probability, epsilon = 1e-5) {
  stats::qlogis(pmin(pmax(as.numeric(probability), epsilon), 1 - epsilon))
}

prepare_static_habitat_training <- function(
    data, features, target_mode = c("annual", "site_propensity")) {
  target_mode <- match.arg(target_mode)
  if (target_mode == "annual") {
    result <- data
    result$observation_weight <- 1
    return(result)
  }

  data |>
    dplyr::group_by(.data$site_id) |>
    dplyr::summarise(
      high_density = mean(.data$high_density),
      observation_weight = sqrt(dplyr::n()),
      dplyr::across(dplyr::all_of(features), dplyr::first),
      .groups = "drop"
    )
}

fit_static_habitat_candidate <- function(
    train, features, target_mode, model_cfg, seed_offset = 0L) {
  training <- prepare_static_habitat_training(
    train, features, target_mode = target_mode
  )
  fit_xgb_binary(
    training, features, model_cfg,
    seed_offset = seed_offset,
    weight_col = "observation_weight"
  )
}

fit_rank_locked_state_bridge <- function(
    truth, habitat_probability, outbreak_probability) {
  training <- data.frame(
    truth = as.integer(truth),
    habitat_logit = safe_component_logit(habitat_probability),
    outbreak_logit = safe_component_logit(outbreak_probability)
  )
  fit <- suppressWarnings(stats::glm(
    truth ~ habitat_logit + outbreak_logit,
    data = training,
    family = stats::binomial()
  ))
  coefficients <- stats::coef(fit)
  valid <- isTRUE(fit$converged) && all(is.finite(coefficients)) &&
    coefficients[["habitat_logit"]] > 0
  if (!valid) {
    stop("Rank-locked state bridge failed or reversed habitat rank.")
  }
  structure(list(
    coefficients = coefficients,
    training_rows = nrow(training),
    rank_locked_within_reef_year = TRUE
  ), class = "cots_rank_locked_state_bridge")
}

predict_rank_locked_state_bridge <- function(
    bridge, habitat_probability, outbreak_probability) {
  b <- bridge$coefficients
  stats::plogis(
    b[["(Intercept)"]] +
      b[["habitat_logit"]] * safe_component_logit(habitat_probability) +
      b[["outbreak_logit"]] * safe_component_logit(outbreak_probability)
  )
}

fit_effort_standardiser <- function(
    truth, state_probability, n_dives, bottom_time) {
  training <- data.frame(
    truth = as.integer(truth),
    state_logit = safe_component_logit(state_probability),
    log_dives = log1p(n_dives),
    log_bottom_time = log1p(bottom_time)
  )
  fit <- suppressWarnings(stats::glm(
    truth ~ offset(state_logit) + log_dives + log_bottom_time,
    data = training,
    family = stats::binomial()
  ))
  if (!isTRUE(fit$converged) || any(!is.finite(stats::coef(fit)))) {
    stop("Effort standardisation model did not converge.")
  }
  fit
}

predict_effort_standardised <- function(
    fit, state_probability, n_dives, bottom_time) {
  new_data <- data.frame(
    state_logit = safe_component_logit(state_probability),
    log_dives = log1p(n_dives),
    log_bottom_time = log1p(bottom_time)
  )
  as.numeric(stats::predict(fit, newdata = new_data, type = "response"))
}

fit_repeated_dive_detection <- function(
    n_high_density_dives, n_dives, state_probability,
    bottom_time) {
  training <- data.frame(
    successes = as.integer(n_high_density_dives),
    failures = as.integer(pmax(n_dives - n_high_density_dives, 0)),
    state_logit = safe_component_logit(state_probability),
    log_time_per_dive = log1p(bottom_time / pmax(n_dives, 1))
  )
  fit <- suppressWarnings(stats::glm(
    cbind(successes, failures) ~
      offset(state_logit) + log_time_per_dive,
    data = training,
    family = stats::quasibinomial()
  ))
  if (!isTRUE(fit$converged) || any(!is.finite(stats::coef(fit)))) {
    stop("Repeated-dive detection model did not converge.")
  }
  fit
}

predict_fixed_protocol_detection <- function(
    fit, state_probability, n_dives, bottom_time) {
  n_dives <- rep_len(as.numeric(n_dives), length(state_probability))
  bottom_time <- rep_len(as.numeric(bottom_time), length(state_probability))
  new_data <- data.frame(
    state_logit = safe_component_logit(state_probability),
    log_time_per_dive = log1p(bottom_time / pmax(n_dives, 1))
  )
  per_dive <- as.numeric(stats::predict(
    fit, newdata = new_data, type = "response"
  ))
  1 - (1 - per_dive)^n_dives
}

rank_align_probability <- function(reference, probability) {
  valid <- is.finite(reference) & is.finite(probability)
  result <- rep(NA_real_, length(reference))
  if (!any(valid)) return(result)
  reference_valid <- reference[valid]
  probability_valid <- probability[valid]
  order_reference <- order(reference_valid)
  sorted_probability <- sort(probability_valid)
  reference_ordered <- reference_valid[order_reference]
  tie_groups <- split(
    seq_along(reference_ordered), reference_ordered, drop = TRUE
  )
  aligned_ordered <- sorted_probability
  for (index in tie_groups) {
    aligned_ordered[index] <- mean(sorted_probability[index])
  }
  aligned <- numeric(length(reference_valid))
  aligned[order_reference] <- aligned_ordered
  result[valid] <- aligned
  result
}

rank_lock_check <- function(
    predictions, group = "reef_year_id", tolerance = 1e-12) {
  groups <- split(predictions, predictions[[group]], drop = TRUE)
  checks <- lapply(groups, function(x) {
    if (nrow(x) < 2L || stats::sd(x$habitat_probability) == 0) return(NULL)
    association <- suppressWarnings(stats::cor(
      x$habitat_probability, x$probability,
      method = "spearman", use = "complete.obs"
    ))
    data.frame(rows = nrow(x), association = association)
  })
  checks <- dplyr::bind_rows(checks)
  tibble::tibble(
    groups_checked = nrow(checks),
    minimum_spearman = if (nrow(checks)) min(checks$association) else NA_real_,
    violations = if (nrow(checks)) {
      sum(!is.finite(checks$association) | checks$association < 1 - tolerance)
    } else {
      0L
    }
  )
}
