args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v2_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v2_dir, "R", "bootstrap.R"))
ctx <- initialise_v2(script_file)
cfg <- ctx$config

data_path <- file.path(
  cfg$paths$artifacts_dir, "hmsc_v2_ecological_fit_data.rds"
)
if (!file.exists(data_path)) stop("Run script 37 before script 41.")
data <- readRDS(data_path) |>
  dplyr::filter(!is.na(.data$cull_high_density))

base_features <- c(
  setdiff(grep("^geom_", names(data), value = TRUE), "geom_unknown"),
  "coral_suitability_mean", "RG_turbid", "RG_waves_Hs", "RG_waves_Tp"
)
scale_names <- c(30L, 60L, 120L, 240L, 360L)
feature_sets <- c(
  list(geomorph_only = base_features),
  stats::setNames(lapply(scale_names, function(scale) {
    c(base_features, grep(paste0("_g", scale, "m$"), names(data), value = TRUE))
  }), paste0("geomorph_g", scale_names, "m"))
)

make_partition <- function(group, folds = 3L, seed = 1L) {
  set.seed(seed)
  group_value <- sample(unique(as.character(group)))
  lookup <- stats::setNames(
    rep(seq_len(folds), length.out = length(group_value)), group_value
  )
  unname(lookup[as.character(group)])
}
cell_fold <- make_partition(data$cell_id, seed = 24082026L)
reef_fold <- make_partition(data$reef_key, seed = 25082026L)

fit_predict_glm <- function(train, test, features) {
  train_x <- as.matrix(train[, features, drop = FALSE])
  test_x <- as.matrix(test[, features, drop = FALSE])
  centre <- colMeans(train_x)
  spread <- apply(train_x, 2L, stats::sd)
  spread[!is.finite(spread) | spread < 1e-10] <- 1
  train_x <- sweep(sweep(train_x, 2L, centre, "-"), 2L, spread, "/")
  test_x <- sweep(sweep(test_x, 2L, centre, "-"), 2L, spread, "/")
  fit <- suppressWarnings(stats::glm.fit(
    x = cbind("(Intercept)" = 1, train_x),
    y = train$cull_high_density, family = stats::binomial()
  ))
  coefficient <- fit$coefficients
  coefficient[!is.finite(coefficient)] <- 0
  stats::plogis(as.numeric(cbind(1, test_x) %*% coefficient))
}

prediction_frame <- function(test, probability, candidate, validation) {
  tibble::tibble(
    model = candidate,
    validation = validation,
    site_year_id = paste(test$cell_id, test$year, sep = "||"),
    site_id = test$cell_id,
    reef_name = test$reef_name,
    reef_year_id = test$reef_year_id,
    year = test$year,
    truth = test$cull_high_density,
    probability = probability,
    longitude = test$longitude,
    latitude = test$latitude
  )
}

run_partition <- function(features, candidate, fold, validation) {
  pieces <- lapply(sort(unique(fold)), function(k) {
    train <- data[fold != k, , drop = FALSE]
    test <- data[fold == k, , drop = FALSE]
    prediction_frame(
      test, fit_predict_glm(train, test, features), candidate, validation
    )
  })
  dplyr::bind_rows(pieces)
}

all_predictions <- lapply(names(feature_sets), function(candidate) {
  features <- feature_sets[[candidate]]
  held <- run_partition(
    features, candidate, cell_fold, "held_cell_3fold"
  )
  reef <- run_partition(
    features, candidate, reef_fold, "leave_reef_3fold"
  )
  train <- data$year <= 2023L
  forward <- prediction_frame(
    data[!train, , drop = FALSE],
    fit_predict_glm(
      data[train, , drop = FALSE], data[!train, , drop = FALSE], features
    ),
    candidate, "forward_year_train_through_2023"
  )
  dplyr::bind_rows(held, reef, forward)
}) |>
  dplyr::bind_rows()

metrics <- all_predictions |>
  dplyr::group_by(.data$model, .data$validation) |>
  dplyr::group_modify(~ binary_probability_metrics(
    .x$truth, .x$probability
  )) |>
  dplyr::ungroup()
decomposition <- all_predictions |>
  dplyr::group_split(.data$model, .data$validation) |>
  lapply(function(x) {
    prediction_scale_decomposition(x, group = "reef_year_id") |>
      dplyr::mutate(
        validation = x$validation[[1]],
        model = x$model[[1]], .before = 1L
      )
  }) |>
  dplyr::bind_rows()

selection <- metrics |>
  dplyr::select("model", "validation", overall_roc_auc = "roc_auc") |>
  dplyr::left_join(
    dplyr::select(
      decomposition, "model", "validation", "within_reef_roc_auc"
    ), by = c("model", "validation")
  ) |>
  dplyr::group_by(.data$model) |>
  dplyr::summarise(
    mean_within_reef_auc = mean(.data$within_reef_roc_auc, na.rm = TRUE),
    minimum_within_reef_auc = min(.data$within_reef_roc_auc, na.rm = TRUE),
    mean_overall_auc = mean(.data$overall_roc_auc, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(
    dplyr::desc(.data$mean_within_reef_auc),
    dplyr::desc(.data$minimum_within_reef_auc)
  ) |>
  dplyr::mutate(
    selected_for_reduced_hmsc = dplyr::row_number() == 1L,
    incumbent_rank_used = FALSE
  )

saveRDS(
  all_predictions,
  file.path(cfg$paths$outputs_dir, "hmsc_v2_geomorphic_scale_screen_predictions.rds")
)
tables <- list(
  hmsc_v2_geomorphic_scale_screen_metrics = metrics,
  hmsc_v2_geomorphic_scale_screen_decomposition = decomposition,
  hmsc_v2_geomorphic_scale_selection = selection
)
for (name in names(tables)) {
  utils::write.csv(
    tables[[name]], file.path(cfg$paths$outputs_dir, paste0(name, ".csv")),
    row.names = FALSE
  )
}
message("Independent geomorphic scale screen complete.")
print(selection)
print(decomposition |>
  dplyr::select("model", "validation", "within_reef_roc_auc"))
