args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v4_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v4_dir, "R", "bootstrap.R"))
ctx <- initialise_v4(script_file)
cfg <- ctx$config
set.seed(20260909L)

data <- readRDS(file.path(cfg$paths$artifacts_dir, "corrected_cell_units_120m.rds")) |>
  dplyr::filter(!is.na(.data$cull_high_density))
families <- terrain_families(names(data))
raw_terms <- unlist(families, use.names = FALSE)
data <- data[stats::complete.cases(data[, raw_terms, drop = FALSE]), , drop = FALSE]

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

fit_ridge <- function(train, test, terms, lambda = 0.01) {
  pair <- model_matrix_pair_v4(train, test, terms)
  keep <- apply(pair$train, 2, stats::sd) > 1e-10
  x_train <- pair$train[, keep, drop = FALSE]
  x_test <- pair$test[, keep, drop = FALSE]
  fit <- glmnet::glmnet(
    x_train, train$cull_high_density, family = "binomial", alpha = 0.25,
    lambda = lambda, standardize = TRUE
  )
  as.numeric(stats::predict(fit, newx = x_test, type = "response", s = lambda))
}

inner_score <- function(train, terms) {
  fold <- make_spatial_fold_v4(
    train$projected_x, train$projected_y,
    cfg$validation$spatial_block_m, cfg$validation$folds
  )
  probability <- rep(NA_real_, nrow(train))
  for (f in sort(unique(fold))) {
    probability[fold == f] <- fit_ridge(
      train[fold != f, , drop = FALSE], train[fold == f, , drop = FALSE], terms
    )
  }
  within_reef_year_auc_v4(tibble::tibble(
    reef_year_id = train$reef_year_id,
    truth = train$cull_high_density,
    probability = probability
  ))
}

select_scales <- function(train, base_terms) {
  dplyr::bind_rows(lapply(names(families), function(family) {
    candidates <- families[[family]]
    scores <- vapply(candidates, function(term) {
      inner_score(train, c(base_terms, term))
    }, numeric(1))
    tibble::tibble(
      family = family, term = candidates, inner_within_auc = scores,
      selected = candidates == candidates[which.max(scores)]
    )
  }))
}

make_representation <- function(train, test, representation, base_terms, selected = NULL) {
  if (representation == "base_only") {
    return(list(train = train, test = test, terms = base_terms))
  }
  if (representation == "all_scale_shrinkage") {
    return(list(train = train, test = test, terms = c(base_terms, raw_terms)))
  }
  if (representation == "selected_actual_scale") {
    return(list(
      train = train, test = test,
      terms = c(base_terms, selected$term[selected$selected])
    ))
  }
  max_components <- if (representation == "family_pca_1") 1L else 2L
  variance <- if (representation == "family_pca_1") 0 else 0.80
  pca <- fit_family_pca(train, variance = variance, max_components = max_components)
  train_pc <- apply_family_pca(train, pca)
  test_pc <- apply_family_pca(test, pca)
  pc_names <- names(train_pc)
  list(
    train = dplyr::bind_cols(train, train_pc),
    test = dplyr::bind_cols(test, test_pc),
    terms = c(base_terms, pc_names)
  )
}

domains <- list(
  clean = c("AhyaD", "Aspat", "Aten"),
  enhanced = c("geomorphic_class", "AhyaD", "Aspat", "Aten", "turbid", "waves_Hs", "waves_Tp")
)
representations <- c(
  "base_only", "selected_actual_scale", "family_pca_1",
  "family_pca_2", "all_scale_shrinkage"
)

predictions <- list()
selection_rows <- list()
append_partition <- function(train_index, test_index, validation, fold_label) {
  train_all <- data[train_index, , drop = FALSE]
  test_all <- data[test_index, , drop = FALSE]
  for (domain in names(domains)) {
    base_terms <- domains[[domain]]
    needed <- unique(c(base_terms, raw_terms))
    train <- train_all[stats::complete.cases(train_all[, needed, drop = FALSE]), , drop = FALSE]
    test <- test_all[stats::complete.cases(test_all[, needed, drop = FALSE]), , drop = FALSE]
    if (nrow(train) < 200L || nrow(test) < 20L || length(unique(test$cull_high_density)) < 2L) next
    selected <- select_scales(train, base_terms)
    selection_rows[[length(selection_rows) + 1L]] <<- dplyr::mutate(
      selected, domain = domain, validation = validation, fold = as.character(fold_label)
    )
    for (representation in representations) {
      represented <- make_representation(
        train, test, representation, base_terms, selected
      )
      probability <- fit_ridge(
        represented$train, represented$test, represented$terms
      )
      predictions[[length(predictions) + 1L]] <<- tibble::tibble(
        domain = domain, representation = representation,
        validation = validation, fold = as.character(fold_label),
        cell_id = test$cell_id, reef_key = test$reef_key,
        reef_year_id = test$reef_year_id, year = test$year,
        truth = test$cull_high_density, probability = probability
      )
    }
  }
  invisible(NULL)
}

held <- make_spatial_fold_v4(
  data$projected_x, data$projected_y,
  cfg$validation$spatial_block_m, cfg$validation$folds
)
for (f in sort(unique(held))) {
  append_partition(held != f, held == f, "held_3km_cell", f)
}
leave <- make_leave_folds(data, cfg$validation$folds)
for (f in sort(unique(leave))) {
  append_partition(leave != f, leave == f, "leave_reef", f)
}
append_partition(data$year <= 2023, data$year >= 2024, "forward_2024_2026", "future")

predictions <- dplyr::bind_rows(predictions)
selection <- dplyr::bind_rows(selection_rows)
metrics <- predictions |>
  dplyr::group_by(.data$domain, .data$representation, .data$validation) |>
  dplyr::summarise(
    rows = dplyr::n(), positives = sum(.data$truth),
    overall_auc = safe_auc_v4(.data$truth, .data$probability),
    within_reef_year_auc = within_reef_year_auc_v4(dplyr::pick(dplyr::everything())),
    probability_sd = stats::sd(.data$probability),
    .groups = "drop"
  )
selection_summary <- selection |>
  dplyr::filter(.data$selected) |>
  dplyr::count(.data$domain, .data$family, .data$term, name = "outer_selections") |>
  dplyr::group_by(.data$domain, .data$family) |>
  dplyr::mutate(selection_fraction = .data$outer_selections / sum(.data$outer_selections)) |>
  dplyr::ungroup() |>
  dplyr::arrange(.data$domain, .data$family, dplyr::desc(.data$outer_selections))

saveRDS(predictions, file.path(cfg$paths$outputs_dir, "terrain_domain_screen_predictions_120m.rds"))
utils::write.csv(metrics, file.path(cfg$paths$outputs_dir, "terrain_domain_screen_metrics_120m.csv"), row.names = FALSE)
utils::write.csv(selection, file.path(cfg$paths$outputs_dir, "nested_scale_selection_120m.csv"), row.names = FALSE)
utils::write.csv(selection_summary, file.path(cfg$paths$outputs_dir, "nested_scale_selection_summary_120m.csv"), row.names = FALSE)
print(metrics)
print(dplyr::filter(selection_summary, .data$selection_fraction == max(.data$selection_fraction), .by = c("domain", "family")))
