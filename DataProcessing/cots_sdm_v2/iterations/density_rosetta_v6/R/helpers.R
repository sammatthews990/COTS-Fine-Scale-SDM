rosetta_summary_v6 <- function(lib_dir) {
  .libPaths(c(lib_dir, .libPaths()))
  sample_file <- system.file(
    "extdata", "samples_cullFromSALAD.RDS", package = "RosettaCOTS"
  )
  if (!nzchar(sample_file)) stop("RosettaCOTS calibration samples not found.")
  samples <- as.matrix(readRDS(sample_file))
  list(
    samples = samples,
    median = stats::setNames(apply(samples, 2, stats::median), colnames(samples)),
    intervals = apply(samples, 2, stats::quantile, c(0.025, 0.5, 0.975))
  )
}

predictor_design_v6 <- function(support_m, paths) {
  fit_name <- if (support_m == 60L) {
    "enhanced_cots_60support_pca1_fit_data_60m.rds"
  } else {
    "enhanced_cots_fit_data_120m.rds"
  }
  fit_artifact <- readRDS(file.path(paths$v4, "artifacts", fit_name))
  predictors <- readRDS(file.path(
    paths$v3, "artifacts", sprintf("support_predictors_%dm.rds", support_m)
  )) |>
    dplyr::distinct(.data$cell_id, .keep_all = TRUE)
  base_names <- setdiff(
    names(fit_artifact$XData),
    grep("_pc[0-9]+$", names(fit_artifact$XData), value = TRUE)
  )
  raw_names <- unique(unlist(lapply(fit_artifact$pca, `[[`, "columns"), use.names = FALSE))
  required <- unique(c(base_names, raw_names))
  complete <- stats::complete.cases(predictors[, required, drop = FALSE])
  predictors <- predictors[complete, , drop = FALSE]
  pcs <- apply_family_pca(predictors, fit_artifact$pca)
  xdata <- dplyr::bind_cols(
    dplyr::select(predictors, dplyr::all_of(base_names)), pcs
  )
  xdata$geomorphic_class <- factor(
    xdata$geomorphic_class,
    levels = levels(fit_artifact$XData$geomorphic_class)
  )
  keep <- stats::complete.cases(xdata)
  predictors <- predictors[keep, , drop = FALSE]
  xdata <- xdata[keep, , drop = FALSE]
  X <- stats::model.matrix(fit_artifact$XFormula, xdata)
  continuous <- setdiff(
    colnames(X), c("(Intercept)", grep("^geomorphic_class", colnames(X), value = TRUE))
  )
  centre <- stats::setNames(rep(0, ncol(X)), colnames(X))
  scale_value <- stats::setNames(rep(1, ncol(X)), colnames(X))
  centre[continuous] <- colMeans(X[, continuous, drop = FALSE])
  scale_value[continuous] <- apply(X[, continuous, drop = FALSE], 2, stats::sd)
  scale_value[!is.finite(scale_value) | scale_value <= 1e-10] <- 1
  X <- sweep(X, 2, centre, "-")
  X <- sweep(X, 2, scale_value, "/")
  predictors$design_row <- seq_len(nrow(predictors))
  list(
    cells = predictors, X = unname(X), columns = colnames(X),
    centre = centre, scale = scale_value, fit_artifact = fit_artifact,
    base_names = base_names, raw_names = raw_names
  )
}

prepare_density_support_v6 <- function(support_m, ctx, minimum_track_coverage = 0.95) {
  paths <- ctx$paths
  design <- predictor_design_v6(support_m, paths)
  activity <- readRDS(file.path(
    paths$v4, "artifacts", sprintf("reef_year_activity_%dm.rds", support_m)
  )) |>
    dplyr::select("reef_key", "year", "active_primary")
  visits <- readRDS(file.path(
    paths$v3, "artifacts", sprintf("cell_visits_%dm_window_3d.rds", support_m)
  ))
  links <- readRDS(file.path(
    paths$v3, "artifacts", sprintf("footprint_links_%dm.rds", support_m)
  ))

  cull <- visits |>
    dplyr::filter(
      .data$method == "cull", .data$ecological_eligible,
      is.finite(.data$direct_count), .data$direct_count >= 0,
      is.finite(.data$effort_minutes), .data$effort_minutes > 0
    ) |>
    dplyr::inner_join(activity, by = c("reef_key", "year")) |>
    dplyr::filter(.data$active_primary) |>
    dplyr::inner_join(
      dplyr::select(design$cells, "cell_id", "design_row"), by = "cell_id"
    ) |>
    dplyr::transmute(
      observation_id = .data$visit_id, .data$visit_id, .data$survey_date,
      .data$year, .data$reef_key, .data$reef_name, .data$cell_id,
      .data$projected_x, .data$projected_y, .data$design_row,
      count = as.integer(round(.data$direct_count)),
      effort_minutes = as.numeric(.data$effort_minutes),
      observed_cpue = .data$count / .data$effort_minutes
    )

  salad_raw <- suppressWarnings(readxl::read_excel(
    ctx$cfg$paths$salad_workbook, sheet = "Tracks", guess_max = 10000
  )) |>
    dplyr::transmute(
      visit_id = paste0("salad||", .data$Count),
      curated_total = as.integer(round(as.numeric(.data$Total))),
      reported_area_m2 = as.numeric(.data$`AREA(5MWIDE)`),
      reported_density_ha = as.numeric(.data$Density_inclScars)
    )

  salad_all_links <- links$salad |>
    dplyr::inner_join(activity, by = c("reef_key", "year")) |>
    dplyr::filter(.data$active_primary) |>
    dplyr::inner_join(salad_raw, by = "visit_id")
  salad_complete_links <- salad_all_links |>
    dplyr::inner_join(
      dplyr::select(design$cells, "cell_id", "design_row"), by = "cell_id"
    )
  coverage <- salad_all_links |>
    dplyr::group_by(.data$visit_id) |>
    dplyr::summarise(linked_area_m2 = sum(.data$weighted_area_m2), .groups = "drop") |>
    dplyr::left_join(
      salad_complete_links |>
        dplyr::group_by(.data$visit_id) |>
        dplyr::summarise(complete_area_m2 = sum(.data$weighted_area_m2), .groups = "drop"),
      by = "visit_id"
    ) |>
    dplyr::mutate(
      complete_area_m2 = dplyr::coalesce(.data$complete_area_m2, 0),
      predictor_coverage = .data$complete_area_m2 / .data$linked_area_m2
    )
  retained_ids <- coverage$visit_id[coverage$predictor_coverage >= minimum_track_coverage]
  salad_links <- salad_complete_links |>
    dplyr::filter(.data$visit_id %in% retained_ids) |>
    dplyr::group_by(.data$visit_id) |>
    dplyr::mutate(
      model_area_ha = (.data$weighted_area_m2 / sum(.data$weighted_area_m2)) *
        dplyr::first(.data$reported_area_m2) / 10000
    ) |>
    dplyr::ungroup() |>
    dplyr::select("visit_id", "cell_id", "design_row", "model_area_ha")
  salad <- salad_complete_links |>
    dplyr::filter(.data$visit_id %in% retained_ids) |>
    dplyr::group_by(.data$visit_id) |>
    dplyr::summarise(
      observation_id = dplyr::first(.data$visit_id),
      survey_date = dplyr::first(.data$survey_date), year = dplyr::first(.data$year),
      reef_key = dplyr::first(.data$reef_key), reef_name = dplyr::first(.data$reef_name),
      projected_x = stats::weighted.mean(.data$projected_x, .data$weighted_area_m2),
      projected_y = stats::weighted.mean(.data$projected_y, .data$weighted_area_m2),
      count = dplyr::first(.data$curated_total),
      effort_area_ha = dplyr::first(.data$reported_area_m2) / 10000,
      observed_density_ha = dplyr::first(.data$reported_density_ha),
      .groups = "drop"
    )

  reef_levels <- sort(unique(c(cull$reef_key, salad$reef_key)))
  add_folds <- function(x) dplyr::mutate(
    x,
    held_fold = make_spatial_fold_v4(.data$projected_x, .data$projected_y, 3000L, 5L),
    reef_fold = 1L + (match(.data$reef_key, reef_levels) - 1L) %% 5L,
    reef_year_id = paste(.data$reef_key, .data$year, sep = "||")
  )
  cull <- add_folds(cull)
  salad <- add_folds(salad)
  calibration <- rosetta_summary_v6(paths$lib)
  artifact <- list(
    support_m = support_m, design = design, cull = cull, salad = salad,
    salad_links = salad_links, coverage = coverage,
    calibration = calibration, minimum_track_coverage = minimum_track_coverage,
    target = "active-reef-year static habitat; SALAD-equivalent COTS-plus-scars density per hectare"
  )
  artifact
}

subset_tmb_data_v6 <- function(artifact, model_type, cull_ids = NULL, salad_ids = NULL,
                               calibration_draw = NULL) {
  if (is.null(cull_ids)) cull_ids <- artifact$cull$observation_id
  if (is.null(salad_ids)) salad_ids <- artifact$salad$observation_id
  use_cull <- model_type %in% c("cull_native", "joint_rosetta")
  use_salad <- model_type %in% c("salad_native", "joint_rosetta")
  cull <- artifact$cull |>
    dplyr::filter(.data$observation_id %in% cull_ids)
  salad <- artifact$salad |>
    dplyr::filter(.data$observation_id %in% salad_ids)
  if (!use_cull) cull <- cull[0, , drop = FALSE]
  if (!use_salad) salad <- salad[0, , drop = FALSE]
  salad_links <- artifact$salad_links |>
    dplyr::filter(.data$visit_id %in% salad$observation_id)
  salad <- salad[match(unique(salad_links$visit_id), salad$observation_id), , drop = FALSE]
  salad_links$track_index <- match(salad_links$visit_id, salad$observation_id) - 1L
  p <- ncol(artifact$design$X)
  Xc <- if (nrow(cull)) artifact$design$X[cull$design_row, , drop = FALSE] else matrix(numeric(), 0, p)
  Xs <- if (nrow(salad_links)) artifact$design$X[salad_links$design_row, , drop = FALSE] else matrix(numeric(), 0, p)

  if (model_type == "joint_rosetta") {
    cal <- if (is.null(calibration_draw)) artifact$calibration$median else calibration_draw
    cal_alpha <- unname(cal[["alpha"]]); cal_beta <- unname(cal[["beta"]])
    cal_offset <- -log(10000)
  } else {
    cal_alpha <- 0; cal_beta <- 1; cal_offset <- 0
  }
  list(
    data = list(
      X_cull = unname(Xc), cull_y = as.integer(cull$count),
      cull_effort = as.numeric(cull$effort_minutes),
      X_salad_link = unname(Xs),
      salad_link_track = as.integer(salad_links$track_index),
      salad_link_area_ha = as.numeric(salad_links$model_area_ha),
      salad_y = as.integer(salad$count),
      calibration_alpha = cal_alpha, calibration_beta = cal_beta,
      calibration_log_density_offset = cal_offset
    ),
    cull = cull, salad = salad, salad_links = salad_links
  )
}

fit_density_tmb_v6 <- function(artifact, model_type, dll_name,
                               cull_ids = NULL, salad_ids = NULL,
                               calibration_draw = NULL, start = NULL) {
  input <- subset_tmb_data_v6(
    artifact, model_type, cull_ids, salad_ids, calibration_draw
  )
  p <- ncol(artifact$design$X)
  if (is.null(start)) {
    beta <- rep(0, p)
    beta[[1]] <- if (model_type == "cull_native") {
      log(max(stats::median(input$cull$count / input$cull$effort_minutes), 0.005))
    } else {
      log(max(stats::median(input$salad$count / input$salad$effort_area_ha), 1))
    }
    start <- list(beta = beta, log_phi_cull = log(5), log_phi_salad = log(5))
  }
  objective <- TMB::MakeADFun(input$data, start, DLL = dll_name, silent = TRUE)
  fit <- stats::nlminb(
    objective$par, objective$fn, objective$gr,
    control = list(iter.max = 1500, eval.max = 2500, rel.tol = 1e-10)
  )
  gradient <- objective$gr(fit$par)
  hessian <- try(stats::optimHess(fit$par, objective$fn, objective$gr), silent = TRUE)
  pd_hessian <- !inherits(hessian, "try-error") && all(is.finite(hessian)) &&
    all(eigen(hessian, symmetric = TRUE, only.values = TRUE)$values > 0)
  list(
    model_type = model_type, fit = fit, beta = fit$par[grepl("^beta", names(fit$par))],
    max_gradient = max(abs(gradient)), pd_hessian = pd_hessian,
    calibration_draw = calibration_draw, input_counts = c(cull = nrow(input$cull), salad = nrow(input$salad))
  )
}

cell_target_v6 <- function(model, artifact, design_rows) {
  exp(as.numeric(artifact$design$X[design_rows, , drop = FALSE] %*% model$beta))
}

predict_observations_v6 <- function(model, artifact, cull_ids, salad_ids) {
  input <- subset_tmb_data_v6(
    artifact, model$model_type, cull_ids, salad_ids, model$calibration_draw
  )
  cal <- input$data
  cull <- input$cull
  if (nrow(cull)) {
    target <- cell_target_v6(model, artifact, cull$design_row)
    cull$observed <- cull$count / cull$effort_minutes
    cull$predicted <- if (model$model_type == "joint_rosetta") {
      exp(cal$calibration_alpha + cal$calibration_beta *
        (log(target) + cal$calibration_log_density_offset))
    } else target
    cull$method <- "cull_cpue"
  }
  salad <- input$salad
  if (nrow(salad)) {
    links <- input$salad_links
    links$density <- cell_target_v6(model, artifact, links$design_row)
    pred <- links |>
      dplyr::group_by(.data$visit_id) |>
      dplyr::summarise(
        predicted = sum(.data$model_area_ha * .data$density) / sum(.data$model_area_ha),
        .groups = "drop"
      )
    salad <- salad |>
      dplyr::left_join(pred, by = c("observation_id" = "visit_id")) |>
      dplyr::mutate(observed = .data$count / .data$effort_area_ha, method = "salad_density_ha")
  }
  dplyr::bind_rows(
    if (nrow(cull)) dplyr::select(cull, "observation_id", "year", "reef_key", "reef_year_id", "method", "observed", "predicted") else NULL,
    if (nrow(salad)) dplyr::select(salad, "observation_id", "year", "reef_key", "reef_year_id", "method", "observed", "predicted") else NULL
  )
}

within_group_spearman_v6 <- function(data) {
  z <- data |>
    dplyr::group_by(.data$reef_year_id) |>
    dplyr::summarise(
      n = dplyr::n(),
      rho = if (dplyr::n() >= 3L && stats::sd(.data$observed) > 0 && stats::sd(.data$predicted) > 0) {
        stats::cor(.data$observed, .data$predicted, method = "spearman")
      } else NA_real_,
      .groups = "drop"
    ) |>
    dplyr::filter(is.finite(.data$rho))
  if (!nrow(z)) return(NA_real_)
  stats::weighted.mean(z$rho, z$n)
}

density_metrics_v6 <- function(data) {
  data |>
    dplyr::group_by(.data$method) |>
    dplyr::group_modify(~ {
      z <- .x |>
        dplyr::filter(is.finite(.data$observed), is.finite(.data$predicted))
      if (!nrow(z)) return(tibble::tibble())
      threshold <- stats::quantile(z$predicted, 0.8, na.rm = TRUE)
      overall_mean <- mean(z$observed)
      tibble::tibble(
        rows = nrow(z),
        spearman = if (stats::sd(z$observed) > 0 && stats::sd(z$predicted) > 0) stats::cor(z$observed, z$predicted, method = "spearman") else NA_real_,
        within_reef_year_spearman = within_group_spearman_v6(z),
        log1p_rmse = sqrt(mean((log1p(z$observed) - log1p(z$predicted))^2)),
        log1p_mae = mean(abs(log1p(z$observed) - log1p(z$predicted))),
        top20_density_lift = if (overall_mean > 0) mean(z$observed[z$predicted >= threshold]) / overall_mean else NA_real_
      )
    }) |>
    dplyr::ungroup()
}
