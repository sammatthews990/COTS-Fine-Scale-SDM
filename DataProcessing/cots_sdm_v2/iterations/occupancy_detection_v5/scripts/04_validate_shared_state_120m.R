args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
runtime_args <- commandArgs(trailingOnly = TRUE)
if (length(runtime_args) != 1L || !runtime_args[[1]] %in% c("held_cell", "leave_reef", "forward")) {
  stop("Usage: Rscript 04_validate_shared_state_120m.R held_cell|leave_reef|forward")
}
validation <- runtime_args[[1]]
v5_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
repo_root <- normalizePath(file.path(v5_dir, "..", "..", "..", ".."), winslash = "/")
v4_dir <- file.path(repo_root, "DataProcessing", "cots_sdm_v2", "iterations", "hmsc_corrected_v4")
source(file.path(v4_dir, "R", "bootstrap.R"))
artifact <- readRDS(file.path(v5_dir, "artifacts", "shared_state_data_120m.rds"))
sites <- artifact$sites

make_reef_folds <- function(data, k = 5L) {
  reefs <- data |>
    dplyr::group_by(.data$reef_key) |>
    dplyr::summarise(x = mean(.data$projected_x), y = mean(.data$projected_y), .groups = "drop")
  set.seed(20260910L)
  reefs$fold <- stats::kmeans(scale(reefs[, c("x", "y")]), k, nstart = 50)$cluster
  reefs$fold[match(data$reef_key, reefs$reef_key)]
}

if (validation == "held_cell") {
  partition <- make_spatial_fold_v4(sites$projected_x, sites$projected_y, 3000L, 5L)
  folds <- sort(unique(partition)); label <- "held_3km_cell_5fold"
} else if (validation == "leave_reef") {
  partition <- make_reef_folds(sites, 5L)
  folds <- sort(unique(partition)); label <- "leave_reef_5fold"
} else {
  partition <- ifelse(sites$year <= 2023, 1L, 2L)
  folds <- 2L; label <- "forward_train_through_2023"
}

truth <- artifact$point |>
  dplyr::filter(.data$method == "cull", .data$post_cull == 0L, !is.na(.data$direct_detection)) |>
  dplyr::arrange(.data$site_id, .data$survey_date) |>
  dplyr::distinct(.data$site_id, .keep_all = TRUE) |>
  dplyr::transmute(
    site_id, truth = as.integer(.data$direct_detection), reef_key,
    reef_year_id = paste(.data$reef_key, .data$year, sep = "||"), year
  )
truth$site <- match(truth$site_id, sites$site_id) - 1L

parameters <- list(
  beta_state = rep(0, ncol(artifact$data$X)),
  alpha_direct = c(-1.5, -1.5, -1.5), effort_direct = rep(0.2, 3),
  depletion_post = -0.5, depletion_log_days = -0.1,
  alpha_scar = c(-1.5, -1.5), effort_scar = rep(0, 2),
  logit_fp_scar = c(-4, -4), scar_prior_effect = 0.5,
  scar_decay = -0.1, scar_post_cull_effect = 0,
  alpha_manta_direct = -2, effort_manta_direct = 0.2,
  alpha_manta_scar = -2, effort_manta_scar = 0.2,
  logit_fp_manta_scar = -4
)

subset_data <- function(train_site) {
  value <- artifact$data
  keep_direct <- train_site[value$direct_site + 1L]
  for (name in c("direct_site", "direct_method", "direct_y", "direct_log_effort", "direct_post_cull", "direct_log_days_after_cull")) {
    value[[name]] <- value[[name]][keep_direct]
  }
  keep_scar <- train_site[value$scar_site + 1L]
  for (name in c("scar_site", "scar_method", "scar_y", "scar_prior_direct_known", "scar_log_days_since_direct", "scar_post_cull")) {
    value[[name]] <- value[[name]][keep_scar]
  }
  point_sites <- unique(c(value$direct_site, value$scar_site))
  value$has_point_observation <- integer(nrow(value$X))
  value$has_point_observation[point_sites + 1L] <- 1L

  links <- artifact$manta_links
  group_ok <- links |>
    dplyr::group_by(.data$group_index) |>
    dplyr::summarise(keep = all(train_site[.data$site + 1L]), .groups = "drop") |>
    dplyr::filter(.data$keep)
  links <- dplyr::semi_join(links, group_ok, by = "group_index") |>
    dplyr::arrange(.data$group_index)
  groups <- artifact$manta_groups |>
    dplyr::mutate(old_group_index = dplyr::row_number()) |>
    dplyr::semi_join(group_ok, by = c("old_group_index" = "group_index")) |>
    dplyr::arrange(.data$old_group_index)
  remap <- stats::setNames(seq_len(nrow(groups)), groups$old_group_index)
  links$new_group <- unname(remap[as.character(links$group_index)])
  links <- dplyr::filter(links, !is.na(.data$new_group)) |>
    dplyr::arrange(.data$new_group)
  counts <- tabulate(links$new_group, nbins = nrow(groups))
  value$manta_group_start <- as.integer(c(0L, cumsum(counts)))
  value$manta_group_type <- as.integer(groups$manta_group_type)
  value$manta_link_site <- as.integer(links$site)
  value$manta_link_log_effort <- as.numeric(links$log_effort)
  value$manta_group_weight <- rep(unique(artifact$data$manta_group_weight)[1], nrow(groups))
  value
}

source_dir <- file.path(v5_dir, "src")
old <- setwd(source_dir); on.exit(setwd(old), add = TRUE)
dll <- TMB::dynlib("shared_state_track")
dyn.load(dll); on.exit(try(dyn.unload(dll), silent = TRUE), add = TRUE)
predictions <- list(); fit_audit <- list()
for (fold in folds) {
  test_site <- partition == fold
  train_site <- !test_site
  data_fold <- subset_data(train_site)
  objective <- TMB::MakeADFun(data_fold, parameters, DLL = "shared_state_track", silent = TRUE)
  fit <- stats::nlminb(
    objective$par, objective$fn, objective$gr,
    control = list(iter.max = 1500, eval.max = 2500, rel.tol = 1e-10)
  )
  state_par <- fit$par[grepl("^beta_state", names(fit$par))]
  psi <- stats::plogis(as.numeric(artifact$data$X %*% state_par))
  use <- truth$site >= 0L & test_site[truth$site + 1L]
  predictions[[length(predictions) + 1L]] <- truth[use, , drop = FALSE] |>
    dplyr::mutate(
      validation = label, fold = as.character(fold),
      probability = psi[.data$site + 1L]
    )
  fit_audit[[length(fit_audit) + 1L]] <- tibble::tibble(
    validation = label, fold = as.character(fold), objective = fit$objective,
    convergence_code = fit$convergence,
    max_absolute_gradient = max(abs(objective$gr(fit$par))),
    training_point_observations = sum(data_fold$has_point_observation),
    training_manta_groups = length(data_fold$manta_group_type)
  )
}
predictions <- dplyr::bind_rows(predictions)
metrics <- tibble::tibble(
  model = "shared_state_track_120m", validation = label,
  rows = nrow(predictions), positives = sum(predictions$truth),
  overall_auc = safe_auc_v4(predictions$truth, predictions$probability),
  within_reef_year_auc = within_reef_year_auc_v4(predictions),
  probability_sd = stats::sd(predictions$probability),
  converged_folds = sum(vapply(fit_audit, function(x) x$convergence_code[[1]] == 0L, logical(1))),
  maximum_fold_gradient = max(vapply(fit_audit, function(x) x$max_absolute_gradient[[1]], numeric(1)))
)
stem <- paste0("shared_state_", validation, "_120m")
saveRDS(predictions, file.path(v5_dir, "outputs", paste0(stem, "_predictions.rds")))
utils::write.csv(metrics, file.path(v5_dir, "outputs", paste0(stem, "_metrics.csv")), row.names = FALSE)
utils::write.csv(dplyr::bind_rows(fit_audit), file.path(v5_dir, "outputs", paste0(stem, "_fit_audit.csv")), row.names = FALSE)
print(metrics)
print(dplyr::bind_rows(fit_audit))
