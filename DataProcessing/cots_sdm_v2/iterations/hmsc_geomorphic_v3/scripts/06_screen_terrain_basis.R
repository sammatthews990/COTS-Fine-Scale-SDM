args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v3_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v3_dir, "R", "bootstrap.R"))
ctx <- initialise_v3(script_file)
cfg <- ctx$config

base_terms <- c("AhyaD", "Aspat", "Aten", "turbid", "waves_Hs", "waves_Tp")
assert_deployable_habitat_terms(base_terms, cfg$model$forbidden_habitat_terms)

make_reef_folds <- function(data, k = 5L) {
  reefs <- data |>
    dplyr::group_by(.data$reef_key) |>
    dplyr::summarise(
      x = mean(.data$projected_x), y = mean(.data$projected_y), .groups = "drop"
    )
  set.seed(20260908L)
  reefs$fold <- stats::kmeans(scale(reefs[, c("x", "y")]), k, nstart = 50)$cluster
  reefs$fold[match(data$reef_key, reefs$reef_key)]
}

build_features <- function(train, test, representation) {
  all_levels <- sort(unique(c(train$geomorphic_class, test$geomorphic_class)))
  train$geomorphic_class <- factor(train$geomorphic_class, levels = all_levels)
  test$geomorphic_class <- factor(test$geomorphic_class, levels = all_levels)
  base_train <- stats::model.matrix(
    ~ geomorphic_class + AhyaD + Aspat + Aten + turbid + waves_Hs + waves_Tp - 1,
    data = train
  )
  base_test <- stats::model.matrix(
    ~ geomorphic_class + AhyaD + Aspat + Aten + turbid + waves_Hs + waves_Tp - 1,
    data = test
  )
  missing <- setdiff(colnames(base_train), colnames(base_test))
  if (length(missing)) {
    base_test <- cbind(base_test, matrix(
      0, nrow(test), length(missing), dimnames = list(NULL, missing)
    ))
  }
  base_test <- base_test[, colnames(base_train), drop = FALSE]
  transform <- NULL
  if (representation == "geomorphic_base") {
    return(list(train = base_train, test = base_test, transform = transform))
  }
  terrain_cols <- unlist(terrain_families(names(train)), use.names = FALSE)
  if (representation == "raw_multiscale_shrinkage") {
    return(list(
      train = cbind(base_train, as.matrix(train[, terrain_cols, drop = FALSE])),
      test = cbind(base_test, as.matrix(test[, terrain_cols, drop = FALSE])),
      transform = "raw"
    ))
  }
  max_components <- if (representation == "family_pca_1") 1L else 2L
  variance <- if (representation == "family_pca_1") 0 else 0.80
  transform <- fit_family_pca(
    train, variance = variance, max_components = max_components
  )
  pca_train <- apply_family_pca(train, transform)
  pca_test <- apply_family_pca(test, transform)
  list(
    train = cbind(base_train, as.matrix(pca_train)),
    test = cbind(base_test, as.matrix(pca_test)), transform = transform
  )
}

fit_one <- function(data, train_rows, test_rows, representation, validation, fold) {
  train <- data[train_rows, , drop = FALSE]
  test <- data[test_rows, , drop = FALSE]
  feature <- build_features(train, test, representation)
  terms <- colnames(feature$train)
  assert_deployable_habitat_terms(terms, cfg$model$forbidden_habitat_terms)
  inner_fold <- 1L + abs(
    floor(train$projected_x / 3000) + 3L * floor(train$projected_y / 3000)
  ) %% 4L
  alpha <- if (representation == "raw_multiscale_shrinkage") 1 else 0.25
  set.seed(20260908L + fold)
  cv <- glmnet::cv.glmnet(
    x = feature$train, y = as.integer(train$positive),
    family = "binomial", alpha = alpha, foldid = inner_fold,
    type.measure = "deviance", nlambda = 50, standardize = TRUE
  )
  probability <- as.numeric(stats::predict(
    cv, newx = feature$test, s = "lambda.1se", type = "response"
  ))
  list(
    predictions = tibble::tibble(
      support_m = unique(data$cell_support_m), model = representation,
      validation = validation, fold = fold, site_id = test$site_id,
      cell_id = test$cell_id, reef_key = test$reef_key,
      reef_year_id = paste(test$reef_key, test$year, sep = "||"),
      year = test$year, truth = as.integer(test$positive), probability = probability
    ),
    fit = list(
      cv = cv, representation = representation, validation = validation,
      fold = fold, transform = feature$transform, terms = terms
    )
  )
}

screen_support <- function(cell_size_m) {
  message("Screening terrain bases at ", cell_size_m, " m.")
  sites <- readRDS(file.path(
    cfg$paths$artifacts_dir, paste0("occupancy_sites_", cell_size_m, "m.rds")
  ))
  predictors <- readRDS(file.path(
    cfg$paths$artifacts_dir, paste0("support_predictors_", cell_size_m, "m.rds")
  ))
  data <- dplyr::left_join(sites, predictors, by = c(
    "cell_id", "projected_x", "projected_y"
  ))
  raw_terms <- unlist(terrain_families(names(data)), use.names = FALSE)
  required <- c(base_terms, "geomorphic_class", raw_terms)
  data <- data[stats::complete.cases(data[, required, drop = FALSE]), , drop = FALSE]
  data$cell_support_m <- cell_size_m
  data$held_cell_fold <- 1L + abs(
    floor(data$projected_x / 3000) + 3L * floor(data$projected_y / 3000)
  ) %% 5L
  data$leave_reef_fold <- make_reef_folds(data, 5L)
  representations <- c(
    "geomorphic_base", "family_pca_1", "family_pca_80", "raw_multiscale_shrinkage"
  )
  fits <- list()
  index <- 0L
  for (representation in representations) {
    for (validation in c("held_cell", "leave_reef")) {
      fold_column <- paste0(validation, "_fold")
      for (fold in sort(unique(data[[fold_column]]))) {
        index <- index + 1L
        message("  ", representation, " ", validation, " fold ", fold)
        fits[[index]] <- fit_one(
          data, data[[fold_column]] != fold, data[[fold_column]] == fold,
          representation, validation, fold
        )
      }
    }
    train <- data$year <= 2023
    test <- data$year >= 2024
    if (sum(train) > 100 && sum(test) > 20 && length(unique(data$positive[test])) == 2L) {
      index <- index + 1L
      fits[[index]] <- fit_one(
        data, train, test, representation, "forward_year_2024_2026", 1L
      )
    }
  }
  saveRDS(lapply(fits, `[[`, "fit"), file.path(
    cfg$paths$models_dir, paste0("terrain_screen_models_", cell_size_m, "m.rds")
  ))
  dplyr::bind_rows(lapply(fits, `[[`, "predictions"))
}

predictions <- dplyr::bind_rows(lapply(cfg$validation$cell_support_m, screen_support))
saveRDS(predictions, file.path(cfg$paths$outputs_dir, "terrain_basis_predictions.rds"))
metrics <- predictions |>
  dplyr::group_by(.data$support_m, .data$model, .data$validation) |>
  dplyr::group_modify(~ {
    overall <- binary_probability_metrics(.x$truth, .x$probability)
    within <- within_group_auc(.x, group = "reef_year_id")$summary
    tibble::tibble(
      rows = nrow(.x), positives = sum(.x$truth),
      overall_auc = overall$roc_auc, pr_auc = overall$pr_auc,
      brier = overall$brier,
      within_reef_year_auc = within$within_group_roc_auc,
      within_pairs = within$positive_negative_pairs
    )
  }) |>
  dplyr::ungroup()
utils::write.csv(metrics, file.path(
  cfg$paths$outputs_dir, "terrain_basis_blocked_metrics.csv"
), row.names = FALSE)
print(metrics)
