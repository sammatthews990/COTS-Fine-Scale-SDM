args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v6_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v6_dir, "R", "bootstrap.R"))
ctx <- initialise_density_v6(script_file)

old <- setwd(ctx$paths$src); on.exit(setwd(old), add = TRUE)
if (!file.exists(TMB::dynlib("density_joint"))) TMB::compile("density_joint.cpp")
dyn.load(TMB::dynlib("density_joint")); on.exit(try(dyn.unload(TMB::dynlib("density_joint")), silent = TRUE), add = TRUE)

model_types <- c("cull_native", "salad_native", "joint_rosetta")
all_metrics <- list(); all_audits <- list(); all_predictions <- list()

for (support_m in c(60L, 120L)) {
  artifact <- readRDS(file.path(ctx$paths$artifacts, sprintf("density_data_%dm.rds", support_m)))
  for (model_type in model_types) {
    message("Full fit: ", model_type, " at ", support_m, " m")
    started <- Sys.time()
    full <- fit_density_tmb_v6(artifact, model_type, "density_joint")
    full$artifact_file <- sprintf("density_data_%dm.rds", support_m)
    saveRDS(full, file.path(ctx$paths$models, sprintf("%s_%dm.rds", model_type, support_m)))
    all_audits[[length(all_audits) + 1L]] <- tibble::tibble(
      support_m = support_m, model = model_type, validation = "apparent",
      convergence_code = full$fit$convergence, convergence_message = full$fit$message,
      objective = full$fit$objective, max_gradient = full$max_gradient,
      positive_definite_hessian = full$pd_hessian,
      elapsed_minutes = as.numeric(difftime(Sys.time(), started, units = "mins"))
    )
    apparent <- predict_observations_v6(
      full, artifact, artifact$cull$observation_id, artifact$salad$observation_id
    )
    if (nrow(apparent)) {
      m <- density_metrics_v6(apparent) |>
        dplyr::mutate(support_m = support_m, model = model_type, validation = "apparent", .before = 1)
      all_metrics[[length(all_metrics) + 1L]] <- m
    }

    schemes <- list(
      held_cell = lapply(1:5, function(k) list(
        train_cull = artifact$cull$observation_id[artifact$cull$held_fold != k],
        test_cull = artifact$cull$observation_id[artifact$cull$held_fold == k],
        train_salad = artifact$salad$observation_id[artifact$salad$held_fold != k],
        test_salad = artifact$salad$observation_id[artifact$salad$held_fold == k]
      )),
      leave_reef = lapply(1:5, function(k) list(
        train_cull = artifact$cull$observation_id[artifact$cull$reef_fold != k],
        test_cull = artifact$cull$observation_id[artifact$cull$reef_fold == k],
        train_salad = artifact$salad$observation_id[artifact$salad$reef_fold != k],
        test_salad = artifact$salad$observation_id[artifact$salad$reef_fold == k]
      )),
      forward_year = list(list(
        train_cull = artifact$cull$observation_id[artifact$cull$year <= 2023],
        test_cull = artifact$cull$observation_id[artifact$cull$year >= 2024],
        train_salad = artifact$salad$observation_id[artifact$salad$year <= 2023],
        test_salad = artifact$salad$observation_id[artifact$salad$year >= 2024]
      ))
    )

    for (scheme_name in names(schemes)) {
      scheme_predictions <- list()
      for (fold in seq_along(schemes[[scheme_name]])) {
        split <- schemes[[scheme_name]][[fold]]
        message("  ", scheme_name, " fold ", fold)
        fit <- fit_density_tmb_v6(
          artifact, model_type, "density_joint",
          cull_ids = split$train_cull, salad_ids = split$train_salad
        )
        pred <- predict_observations_v6(
          fit, artifact, split$test_cull, split$test_salad
        ) |>
          dplyr::mutate(fold = fold)
        scheme_predictions[[fold]] <- pred
        all_audits[[length(all_audits) + 1L]] <- tibble::tibble(
          support_m = support_m, model = model_type,
          validation = paste0(scheme_name, "_fold_", fold),
          convergence_code = fit$fit$convergence,
          convergence_message = fit$fit$message, objective = fit$fit$objective,
          max_gradient = fit$max_gradient, positive_definite_hessian = fit$pd_hessian,
          elapsed_minutes = NA_real_
        )
      }
      predictions <- dplyr::bind_rows(scheme_predictions)
      if (nrow(predictions)) {
        all_predictions[[paste(support_m, model_type, scheme_name, sep = "_")]] <- predictions
        metrics <- density_metrics_v6(predictions) |>
          dplyr::mutate(support_m = support_m, model = model_type, validation = scheme_name, .before = 1)
        all_metrics[[length(all_metrics) + 1L]] <- metrics
      }
    }
  }
}

metrics <- dplyr::bind_rows(all_metrics)
audits <- dplyr::bind_rows(all_audits)
utils::write.csv(metrics, file.path(ctx$paths$outputs, "density_model_metrics.csv"), row.names = FALSE)
utils::write.csv(audits, file.path(ctx$paths$outputs, "density_model_fit_audit.csv"), row.names = FALSE)
saveRDS(all_predictions, file.path(ctx$paths$outputs, "density_blocked_predictions.rds"), compress = "xz")
print(metrics)
print(dplyr::filter(audits, .data$validation == "apparent"))
