args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v4_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v4_dir, "R", "bootstrap.R"))
ctx <- initialise_v4(script_file)
cfg <- ctx$config
set.seed(20260909L)

source_data <- readRDS(file.path(cfg$paths$artifacts_dir, "manta_footprints_120m.rds")) |>
  dplyr::filter(.data$active_primary)
families <- terrain_families(names(source_data))
raw_terms <- unlist(families, use.names = FALSE)

make_labels <- function(data, definition) {
  label <- if (definition == "strong_primary") {
    dplyr::case_when(
      data$manta_strong_positive ~ 1L,
      data$manta_strict_negative ~ 0L,
      TRUE ~ NA_integer_
    )
  } else {
    dplyr::case_when(
      data$manta_strong_positive | data$manta_weak_scar ~ 1L,
      data$manta_strict_negative ~ 0L,
      TRUE ~ NA_integer_
    )
  }
  dplyr::mutate(data, truth = label, evidence_definition = definition) |>
    dplyr::filter(!is.na(.data$truth)) |>
    dplyr::group_by(.data$reef_year_id) |>
    dplyr::mutate(n_positive = sum(.data$truth == 1L)) |>
    dplyr::arrange(.data$reef_year_id, dplyr::desc(.data$truth), .data$visit_id) |>
    dplyr::filter(
      .data$truth == 1L |
        dplyr::row_number() <= .data$n_positive + dplyr::if_else(
          .data$n_positive > 0L,
          cfg$manta$negative_to_positive_ratio * .data$n_positive,
          10L
        )
    ) |>
    dplyr::ungroup()
}

make_leave_folds <- function(data, k = 5L) {
  reefs <- data |>
    dplyr::group_by(.data$reef_key) |>
    dplyr::summarise(
      x = mean(.data$projected_x), y = mean(.data$projected_y), .groups = "drop"
    )
  set.seed(20260909L)
  reefs$fold <- stats::kmeans(scale(reefs[, c("x", "y")]), k, nstart = 50)$cluster
  reefs$fold[match(data$reef_key, reefs$reef_key)]
}

fit_predict <- function(train, test, base_terms) {
  pca <- fit_family_pca(train, variance = 0, max_components = 1L)
  train_pc <- apply_family_pca(train, pca)
  test_pc <- apply_family_pca(test, pca)
  train <- dplyr::bind_cols(train, train_pc)
  test <- dplyr::bind_cols(test, test_pc)
  pair <- model_matrix_pair_v4(train, test, c(base_terms, names(train_pc)))
  keep <- apply(pair$train, 2, stats::sd) > 1e-10
  fit <- glmnet::glmnet(
    pair$train[, keep, drop = FALSE], train$truth,
    family = "binomial", alpha = 0.25, lambda = 0.01, standardize = TRUE
  )
  as.numeric(stats::predict(
    fit, newx = pair$test[, keep, drop = FALSE], type = "response", s = 0.01
  ))
}

domains <- list(
  clean = c("AhyaD", "Aspat", "Aten"),
  enhanced = c("geomorphic_class", "AhyaD", "Aspat", "Aten", "turbid", "waves_Hs", "waves_Tp")
)
predictions <- list()
append_partition <- function(data, train_index, test_index, definition, validation, fold_label) {
  for (domain in names(domains)) {
    needed <- c(domains[[domain]], raw_terms)
    train <- data[train_index & stats::complete.cases(data[, needed, drop = FALSE]), , drop = FALSE]
    test <- data[test_index & stats::complete.cases(data[, needed, drop = FALSE]), , drop = FALSE]
    if (nrow(train) < 200L || nrow(test) < 20L || length(unique(test$truth)) < 2L) next
    probability <- fit_predict(train, test, domains[[domain]])
    predictions[[length(predictions) + 1L]] <<- tibble::tibble(
      evidence_definition = definition, domain = domain,
      validation = validation, fold = as.character(fold_label),
      visit_id = test$visit_id, reef_key = test$reef_key,
      reef_year_id = test$reef_year_id, year = test$year,
      truth = test$truth, probability = probability
    )
  }
}

for (definition in c("strong_primary", "any_scar_sensitivity")) {
  data <- make_labels(source_data, definition)
  held <- make_spatial_fold_v4(
    data$projected_x, data$projected_y,
    cfg$validation$spatial_block_m, cfg$validation$folds
  )
  for (f in sort(unique(held))) {
    append_partition(data, held != f, held == f, definition, "held_3km_footprint", f)
  }
  leave <- make_leave_folds(data, cfg$validation$folds)
  for (f in sort(unique(leave))) {
    append_partition(data, leave != f, leave == f, definition, "leave_reef", f)
  }
  append_partition(
    data, data$year <= 2023, data$year >= 2024,
    definition, "forward_2024_2026", "future"
  )
}

predictions <- dplyr::bind_rows(predictions)
metrics <- predictions |>
  dplyr::group_by(.data$evidence_definition, .data$domain, .data$validation) |>
  dplyr::summarise(
    rows = dplyr::n(), positives = sum(.data$truth),
    overall_auc = safe_auc_v4(.data$truth, .data$probability),
    within_reef_year_auc = within_reef_year_auc_v4(dplyr::pick(dplyr::everything())),
    probability_sd = stats::sd(.data$probability),
    .groups = "drop"
  )
saveRDS(predictions, file.path(cfg$paths$outputs_dir, "manta_footprint_predictions_120m.rds"))
utils::write.csv(metrics, file.path(cfg$paths$outputs_dir, "manta_footprint_metrics_120m.csv"), row.names = FALSE)
print(metrics)
