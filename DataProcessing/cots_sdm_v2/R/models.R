default_feature_sets <- function(static_predictors) {
  outbreak_features <- c(
    "reef_outbreak_probability_lag1",
    "reef_outbreak_probability_lag2",
    "reef_outbreak_probability_lag3",
    "reef_outbreak_probability_mean_lag1_3",
    "reef_outbreak_probability_max_lag1_3",
    "reef_outbreak_probability_trend_lag1_3"
  )

  list(
    static = static_predictors,
    static_year = c(static_predictors, "year"),
    static_outbreak_history = c(static_predictors, outbreak_features),
    static_year_outbreak_history = c(
      static_predictors, "year", outbreak_features
    )
  )
}

fit_xgb_binary <- function(
    train, features, model_cfg, seed_offset = 0L,
    weight_col = NULL) {
  assert_required_columns(train, c("high_density", features), "Model training data")
  x <- data.matrix(train[, features, drop = FALSE])
  y <- as.numeric(train$high_density)
  fit_seed <- as.integer(model_cfg$seed + seed_offset)
  set.seed(fit_seed)

  weights <- if (is.null(weight_col)) {
    rep(1, nrow(train))
  } else {
    as.numeric(train[[weight_col]])
  }
  if (any(!is.finite(weights)) || any(weights <= 0)) {
    stop("Training weights must be finite and positive.")
  }
  dtrain <- xgboost::xgb.DMatrix(
    data = x, label = y, weight = weights
  )
  params <- list(
    objective = "binary:logistic",
    eval_metric = "logloss",
    eta = model_cfg$eta,
    max_depth = model_cfg$max_depth,
    min_child_weight = model_cfg$min_child_weight,
    gamma = if (is.null(model_cfg$gamma)) 0 else model_cfg$gamma,
    subsample = model_cfg$subsample,
    colsample_bytree = model_cfg$colsample_bytree,
    lambda = model_cfg$lambda,
    alpha = model_cfg$alpha,
    nthread = model_cfg$nthread
  )

  xgboost::xgb.train(
    params = params,
    data = dtrain,
    nrounds = model_cfg$nrounds,
    verbose = 0
  )
}

predict_xgb_binary <- function(model, new_data, features) {
  requireNamespace("xgboost", quietly = TRUE)
  x <- data.matrix(new_data[, features, drop = FALSE])
  as.numeric(stats::predict(model, xgboost::xgb.DMatrix(x)))
}

run_spatial_benchmark <- function(data, feature_sets, model_cfg) {
  folds <- sort(unique(data$spatial_fold))
  predictions <- list()
  fitted_models <- list()
  counter <- 0L

  for (model_name in names(feature_sets)) {
    features <- feature_sets[[model_name]]
    for (fold in folds) {
      counter <- counter + 1L
      train <- data[data$spatial_fold != fold, , drop = FALSE]
      test <- data[data$spatial_fold == fold, , drop = FALSE]

      model <- fit_xgb_binary(
        train,
        features,
        model_cfg,
        seed_offset = counter
      )
      probability <- predict_xgb_binary(model, test, features)

      predictions[[counter]] <- tibble::tibble(
        model = model_name,
        validation = "spatial_unseen_reef_cluster",
        fold = fold,
        site_year_id = test$site_year_id,
        reef_name = test$reef_name,
        site_id = test$site_id,
        year = test$year,
        truth = test$high_density,
        probability = probability
      )
      fitted_models[[paste(model_name, fold, sep = "__fold_")]] <- model
    }
  }

  list(
    predictions = dplyr::bind_rows(predictions),
    models = fitted_models
  )
}

run_forward_benchmark <- function(data, feature_sets, model_cfg, test_year) {
  split <- make_forward_year_split(data, test_year)
  train <- data[split$train, , drop = FALSE]
  test <- data[split$test, , drop = FALSE]

  predictions <- list()
  fitted_models <- list()
  for (i in seq_along(feature_sets)) {
    model_name <- names(feature_sets)[i]
    features <- feature_sets[[i]]
    model <- fit_xgb_binary(train, features, model_cfg, seed_offset = 1000L + i)
    probability <- predict_xgb_binary(model, test, features)

    predictions[[i]] <- tibble::tibble(
      model = model_name,
      validation = paste0("forward_year_", test_year),
      fold = NA_integer_,
      site_year_id = test$site_year_id,
      reef_name = test$reef_name,
      site_id = test$site_id,
      year = test$year,
      truth = test$high_density,
      probability = probability
    )
    fitted_models[[paste0(model_name, "__train_through_", test_year - 1L)]] <- model
  }

  list(
    predictions = dplyr::bind_rows(predictions),
    models = fitted_models,
    split = split
  )
}

run_rolling_origin_benchmark <- function(
    data, feature_sets, model_cfg, test_years) {
  predictions <- list()
  fitted_models <- list()
  counter <- 0L

  for (test_year in sort(unique(test_years))) {
    split <- make_forward_year_split(data, test_year)
    train <- data[split$train, , drop = FALSE]
    test <- data[split$test, , drop = FALSE]
    test$reef_seen_before <- test$reef_name %in% unique(train$reef_name)
    test$site_seen_before <- test$site_id %in% unique(train$site_id)

    for (i in seq_along(feature_sets)) {
      counter <- counter + 1L
      model_name <- names(feature_sets)[i]
      features <- feature_sets[[i]]
      model <- fit_xgb_binary(
        train,
        features,
        model_cfg,
        seed_offset = 10000L + test_year * 10L + i
      )
      probability <- predict_xgb_binary(model, test, features)

      predictions[[counter]] <- tibble::tibble(
        model = model_name,
        validation = "rolling_origin",
        test_year = test_year,
        train_year_max = split$train_year_max,
        n_train = nrow(train),
        site_year_id = test$site_year_id,
        reef_name = test$reef_name,
        site_id = test$site_id,
        reef_seen_before = test$reef_seen_before,
        site_seen_before = test$site_seen_before,
        truth = test$high_density,
        probability = probability
      )
      model_key <- paste0(
        model_name, "__rolling_train_through_", split$train_year_max
      )
      fitted_models[[model_key]] <- model
    }
  }

  list(
    predictions = dplyr::bind_rows(predictions),
    models = fitted_models
  )
}

summarise_benchmark_predictions <- function(predictions) {
  predictions |>
    dplyr::group_by(.data$model, .data$validation) |>
    dplyr::group_modify(~ binary_probability_metrics(
      .x$truth, .x$probability
    )) |>
    dplyr::ungroup()
}

fit_additive_operational_bridge <- function(
    truth, static_probability, outbreak_probability) {
  training <- data.frame(
    truth = as.integer(truth),
    static_logit = stats::qlogis(
      safe_probability(static_probability)
    ),
    outbreak_logit = stats::qlogis(
      safe_probability(outbreak_probability)
    )
  )
  training <- training[stats::complete.cases(training), , drop = FALSE]
  fit <- stats::glm(
    truth ~ static_logit + outbreak_logit,
    data = training,
    family = stats::binomial(),
    model = FALSE,
    x = FALSE,
    y = FALSE
  )
  coefficients <- stats::coef(fit)
  if (!isTRUE(fit$converged) || any(!is.finite(coefficients))) {
    stop("The additive operational bridge did not converge.")
  }
  if (coefficients[["static_logit"]] <= 0) {
    stop(
      "The fitted static coefficient is not positive; the bridge would ",
      "not preserve the fine-scale habitat ranking."
    )
  }
  structure(
    list(
      coefficients = coefficients,
      training_rows = nrow(training),
      converged = fit$converged,
      rank_preserving_within_reef = TRUE
    ),
    class = "cots_additive_operational_bridge"
  )
}

predict_additive_operational_bridge <- function(
    bridge, static_probability, outbreak_probability) {
  coefficients <- bridge$coefficients
  linear_predictor <- coefficients[["(Intercept)"]] +
    coefficients[["static_logit"]] * stats::qlogis(
      safe_probability(static_probability)
    ) +
    coefficients[["outbreak_logit"]] * stats::qlogis(
      safe_probability(outbreak_probability)
    )
  stats::plogis(linear_predictor)
}
