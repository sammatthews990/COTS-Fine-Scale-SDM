solve_intercept_shift <- function(
    probability, truth, weight = rep(1, length(truth))) {
  probability <- safe_probability(probability)
  truth <- as.numeric(truth)
  weight <- as.numeric(weight)
  valid <- is.finite(probability) & is.finite(truth) &
    is.finite(weight) & weight > 0
  probability <- probability[valid]
  truth <- truth[valid]
  weight <- weight[valid]
  if (length(truth) == 0L) return(NA_real_)

  observed <- stats::weighted.mean(truth, weight)
  observed <- safe_probability(observed)
  score <- stats::qlogis(probability)
  objective <- function(shift) {
    stats::weighted.mean(stats::plogis(score + shift), weight) -
      observed
  }
  stats::uniroot(objective, interval = c(-20, 20))$root
}

fit_group_intercept_calibrator <- function(
    probability, truth, group,
    weight = rep(1, length(truth)), prior_weight = 20) {
  data <- tibble::tibble(
    probability = safe_probability(probability),
    truth = as.numeric(truth),
    group = as.character(group),
    weight = as.numeric(weight)
  ) |>
    dplyr::filter(
      is.finite(.data$probability),
      is.finite(.data$truth),
      is.finite(.data$weight),
      .data$weight > 0,
      !is.na(.data$group)
    )
  if (nrow(data) == 0L) stop("No valid rows for intercept calibration.")

  global_shift <- solve_intercept_shift(
    data$probability, data$truth, data$weight
  )
  by_group <- split(data, data$group, drop = TRUE)
  table <- dplyr::bind_rows(lapply(names(by_group), function(group_name) {
    x <- by_group[[group_name]]
    raw_shift <- solve_intercept_shift(
      x$probability, x$truth, x$weight
    )
    effective_weight <- sum(x$weight)
    shrinkage <- effective_weight / (effective_weight + prior_weight)
    shift <- global_shift + shrinkage * (raw_shift - global_shift)
    tibble::tibble(
      group = group_name,
      rows = nrow(x),
      effective_weight = effective_weight,
      observed_prevalence = stats::weighted.mean(x$truth, x$weight),
      uncalibrated_prevalence = stats::weighted.mean(
        x$probability, x$weight
      ),
      raw_intercept_shift = raw_shift,
      shrinkage = shrinkage,
      intercept_shift = shift
    )
  }))

  structure(
    list(
      global_intercept_shift = global_shift,
      prior_weight = prior_weight,
      group_table = table
    ),
    class = "cots_group_intercept_calibrator"
  )
}

predict_group_intercept_calibrator <- function(
    calibrator, probability, group) {
  group <- as.character(group)
  shift <- calibrator$group_table$intercept_shift[
    match(group, calibrator$group_table$group)
  ]
  shift[!is.finite(shift)] <- calibrator$global_intercept_shift
  stats::plogis(stats::qlogis(safe_probability(probability)) + shift)
}
