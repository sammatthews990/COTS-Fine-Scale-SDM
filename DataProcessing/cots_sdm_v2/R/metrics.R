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

