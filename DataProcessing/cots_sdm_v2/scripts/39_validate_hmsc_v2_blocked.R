args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v2_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v2_dir, "R", "bootstrap.R"))
ctx <- initialise_v2(script_file)
cfg <- ctx$config

if (!requireNamespace("Hmsc", quietly = TRUE)) stop("Hmsc is required.")
candidate_arg <- grep("^--candidate=", commandArgs(trailingOnly = TRUE), value = TRUE)
candidate <- if (length(candidate_arg)) {
  sub("^--candidate=", "", candidate_arg[[1]])
} else {
  "full_multiscale"
}
reduced_candidate <- identical(candidate, "geomorphic_g360")
artifact_prefix <- if (reduced_candidate) {
  "hmsc_v2_geomorphic_g360"
} else {
  "hmsc_v2"
}
model_label <- if (reduced_candidate) {
  "hmsc_v2_geomorphic_g360"
} else {
  "hmsc_v2_ecological_core"
}
model_path <- file.path(
  cfg$paths$models_dir,
  if (reduced_candidate) {
    "hmsc_v2_geomorphic_g360_core_model.rds"
  } else {
    "hmsc_v2_ecological_core_model.rds"
  }
)
data_path <- file.path(
  cfg$paths$artifacts_dir,
  if (reduced_candidate) {
    "hmsc_v2_geomorphic_g360_fit_data.rds"
  } else {
    "hmsc_v2_ecological_fit_data.rds"
  }
)
if (!all(file.exists(c(model_path, data_path)))) {
  stop("Run script 37 before script 39.")
}
core <- readRDS(model_path)
data <- readRDS(data_path)
target <- "cull_high_density"
target_index <- match(target, core$spNames)
if (is.na(target_index)) stop("Cull high-density response is absent.")

make_group_partition <- function(group, folds = 3L, seed = 24082026L) {
  set.seed(seed)
  groups <- sample(unique(as.character(group)))
  fold_by_group <- stats::setNames(
    rep(seq_len(folds), length.out = length(groups)), groups
  )
  unname(fold_by_group[as.character(group)])
}

cell_partition <- make_group_partition(data$cell_id, seed = 24082026L)
reef_partition <- make_group_partition(data$reef_key, seed = 25082026L)
partition_manifest <- dplyr::bind_rows(
  tibble::tibble(
    validation = "held_cell_3fold", group = data$cell_id,
    fold = cell_partition
  ) |> dplyr::distinct(),
  tibble::tibble(
    validation = "leave_reef_3fold", group = data$reef_key,
    fold = reef_partition
  ) |> dplyr::distinct()
)
utils::write.csv(
  partition_manifest,
  file.path(cfg$paths$outputs_dir, paste0(artifact_prefix, "_partition_manifest.csv")),
  row.names = FALSE
)

validation_model <- core
validation_model$samples <- 40L
validation_model$thin <- 3L
validation_model$transient <- 250L
validation_model$adaptNf <- rep(250L, validation_model$nr)
validation_model$postList <- lapply(
  core$postList[seq_len(2L)], function(chain) utils::tail(chain, 40L)
)
empty_conditioning <- matrix(
  NA_real_, nrow = validation_model$ny, ncol = validation_model$ns
)

target_prediction_frame <- function(probability, validation) {
  index <- which(!is.na(data[[target]]))
  tibble::tibble(
    model = model_label,
    validation = validation,
    site_year_id = paste(data$cell_id[index], data$year[index], sep = "||"),
    site_id = data$cell_id[index],
    reef_name = data$reef_name[index],
    reef_year_id = data$reef_year_id[index],
    year = data$year[index],
    truth = data[[target]][index],
    probability = probability[index],
    longitude = data$longitude[index],
    latitude = data$latitude[index]
  )
}

run_blocked <- function(partition, label) {
  message("Running ", label, " Hmsc validation.")
  predicted <- Hmsc::computePredictedValues(
    validation_model, partition = partition, Yc = empty_conditioning,
    expected = TRUE, nChains = 2L, nParallel = 2L, verbose = 0
  )
  probability <- apply(
    predicted[, target_index, , drop = FALSE], 1L, mean
  )
  rm(predicted)
  gc()
  target_prediction_frame(probability, label)
}

held_cell_checkpoint <- file.path(
  cfg$paths$outputs_dir,
  paste0(artifact_prefix, "_held_cell_predictions_checkpoint.rds")
)
leave_reef_checkpoint <- file.path(
  cfg$paths$outputs_dir,
  paste0(artifact_prefix, "_leave_reef_predictions_checkpoint.rds")
)
if (file.exists(held_cell_checkpoint)) {
  held_cell <- readRDS(held_cell_checkpoint)
} else {
  held_cell <- run_blocked(cell_partition, "held_cell_3fold")
  saveRDS(held_cell, held_cell_checkpoint)
}
if (file.exists(leave_reef_checkpoint)) {
  leave_reef <- readRDS(leave_reef_checkpoint)
} else {
  leave_reef <- run_blocked(reef_partition, "leave_reef_3fold")
  saveRDS(leave_reef, leave_reef_checkpoint)
}

forward_cutoff <- 2023L
train <- data$year <= forward_cutoff
test <- data$year > forward_cutoff
if (sum(test & !is.na(data[[target]])) < 30L) {
  stop("Too few post-2023 Cull target observations for forward validation.")
}
forward <- Hmsc::Hmsc(
  Y = core$Y[train, , drop = FALSE],
  XData = data[train, names(core$XData), drop = FALSE],
  XFormula = core$XFormula,
  TrData = core$TrData,
  TrFormula = core$TrFormula,
  distr = core$distr,
  XScale = TRUE,
  YScale = TRUE
)
message("Fitting forward-year Hmsc through ", forward_cutoff, ".")
forward <- Hmsc::sampleMcmc(
  forward, samples = 40L, thin = 3L, transient = 250L,
  nChains = 2L, nParallel = 2L, verbose = 0
)
X_test <- stats::model.matrix(
  forward$XFormula,
  data = data[test, names(core$XData), drop = FALSE]
)
X_test <- X_test[, forward$covNames, drop = FALSE]
forward_draws <- unlist(forward$postList, recursive = FALSE)
forward_probability <- rowMeans(vapply(
  forward_draws,
  function(draw) stats::pnorm(
    as.numeric(X_test %*% draw$Beta[, target_index])
  ),
  numeric(sum(test))
))
test_index <- which(test & !is.na(data[[target]]))
test_position <- match(test_index, which(test))
forward_year <- tibble::tibble(
  model = model_label,
  validation = "forward_year_train_through_2023",
  site_year_id = paste(data$cell_id[test_index], data$year[test_index], sep = "||"),
  site_id = data$cell_id[test_index],
  reef_name = data$reef_name[test_index],
  reef_year_id = data$reef_year_id[test_index],
  year = data$year[test_index],
  truth = data[[target]][test_index],
  probability = forward_probability[test_position],
  longitude = data$longitude[test_index],
  latitude = data$latitude[test_index]
)

predictions <- dplyr::bind_rows(held_cell, leave_reef, forward_year)
metrics <- predictions |>
  dplyr::group_by(.data$model, .data$validation) |>
  dplyr::group_modify(~ binary_probability_metrics(
    .x$truth, .x$probability
  )) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    mcmc_status = "short_blocked_validation_refit",
    promotion_eligible = FALSE
  )
decomposition <- predictions |>
  dplyr::group_split(.data$validation) |>
  lapply(function(x) {
    prediction_scale_decomposition(x, group = "reef_year_id") |>
      dplyr::mutate(validation = x$validation[[1]], .before = 1L)
  }) |>
  dplyr::bind_rows()

saveRDS(
  predictions,
  file.path(
    cfg$paths$outputs_dir, paste0(artifact_prefix, "_blocked_predictions.rds")
  )
)
saveRDS(
  forward,
  file.path(
    cfg$paths$models_dir, paste0(artifact_prefix, "_forward_2023_model.rds")
  )
)
utils::write.csv(
  metrics,
  file.path(
    cfg$paths$outputs_dir, paste0(artifact_prefix, "_blocked_metrics.csv")
  ),
  row.names = FALSE
)
utils::write.csv(
  decomposition,
  file.path(
    cfg$paths$outputs_dir,
    paste0(artifact_prefix, "_blocked_decomposition.csv")
  ),
  row.names = FALSE
)
message("Hmsc v2 blocked validations complete.")
print(metrics)
print(decomposition)
