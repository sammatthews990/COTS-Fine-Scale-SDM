args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v3_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v3_dir, "R", "bootstrap.R"))
ctx <- initialise_v3(script_file)
cfg <- ctx$config

fit_support <- function(cell_size_m) {
  message("Fitting formal occupancy/detection model at ", cell_size_m, " m.")
  observed_any <- function(x) {
    if (all(is.na(x))) return(NA_integer_)
    as.integer(any(x == 1L, na.rm = TRUE))
  }
  visits <- readRDS(file.path(
    cfg$paths$artifacts_dir,
    paste0("cell_visits_", cell_size_m, "m_window_3d.rds")
  )) |>
    dplyr::filter(
      .data$ecological_eligible,
      .data$method %in% c("cull", "rhis", "salad"),
      !is.na(.data$direct_detection)
    ) |>
    dplyr::group_by(.data$cell_id, .data$year, .data$survey_date, .data$method) |>
    dplyr::summarise(
      reef_key = dplyr::first(.data$reef_key),
      reef_name = dplyr::first(.data$reef_name),
      projected_x = dplyr::first(.data$projected_x),
      projected_y = dplyr::first(.data$projected_y),
      direct_detection = observed_any(.data$direct_detection),
      effort_m = sum(.data$effort_m, na.rm = TRUE),
      effort_minutes = sum(.data$effort_minutes, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(site_id = paste(.data$cell_id, .data$year, sep = "||"))
  sites <- visits |>
    dplyr::group_by(.data$site_id) |>
    dplyr::summarise(
      cell_id = dplyr::first(.data$cell_id), year = dplyr::first(.data$year),
      reef_key = dplyr::first(.data$reef_key), reef_name = dplyr::first(.data$reef_name),
      projected_x = dplyr::first(.data$projected_x),
      projected_y = dplyr::first(.data$projected_y),
      n_visits = dplyr::n(), positive = any(.data$direct_detection == 1L),
      .groups = "drop"
    ) |>
    dplyr::filter(!is.na(.data$reef_key), .data$reef_key != "")
  positive_counts <- sites |>
    dplyr::group_by(.data$reef_key, .data$year) |>
    dplyr::summarise(n_positive = sum(.data$positive), .groups = "drop") |>
    dplyr::filter(.data$n_positive > 0)
  eligible <- dplyr::inner_join(sites, positive_counts, by = c("reef_key", "year"))
  positive_sites <- dplyr::filter(eligible, .data$positive)
  zero_sites <- eligible |>
    dplyr::filter(!.data$positive) |>
    dplyr::arrange(.data$reef_key, .data$year, .data$cell_id) |>
    dplyr::group_by(.data$reef_key, .data$year) |>
    dplyr::filter(dplyr::row_number() <= 5L * dplyr::first(.data$n_positive)) |>
    dplyr::ungroup()
  sites <- dplyr::bind_rows(positive_sites, zero_sites) |>
    dplyr::distinct(.data$site_id, .keep_all = TRUE)
  visits <- dplyr::semi_join(visits, sites, by = "site_id")
  predictors <- readRDS(file.path(
    cfg$paths$artifacts_dir, paste0("support_predictors_", cell_size_m, "m.rds")
  ))
  sites <- sites |>
    dplyr::left_join(predictors, by = c("cell_id", "projected_x", "projected_y"))

  raw_terms <- unlist(terrain_families(names(sites)), use.names = FALSE)
  required <- c("geomorphic_class", "AhyaD", "Aspat", "Aten", "turbid",
                "waves_Hs", "waves_Tp", raw_terms)
  sites <- sites[stats::complete.cases(sites[, required, drop = FALSE]), , drop = FALSE]
  visits <- dplyr::semi_join(visits, sites, by = "site_id")
  pca <- fit_family_pca(sites, variance = 0.80, max_components = 2L)
  pcs <- apply_family_pca(sites, pca)
  site_design <- stats::model.matrix(
    ~ geomorphic_class + AhyaD + Aspat + Aten + turbid + waves_Hs + waves_Tp - 1,
    data = sites
  )
  site_design <- cbind(site_design, as.matrix(pcs))
  colnames(site_design) <- make.names(colnames(site_design), unique = TRUE)
  design_center <- colMeans(site_design)
  design_scale <- apply(site_design, 2, stats::sd)
  keep <- is.finite(design_scale) & design_scale > 0
  site_design <- sweep(site_design[, keep, drop = FALSE], 2, design_center[keep], "-")
  site_design <- sweep(site_design, 2, design_scale[keep], "/")
  site_covs <- as.data.frame(site_design)

  site_order <- sites$site_id
  y <- wide_observation_matrix(visits, "direct_detection", site_order)
  visits <- visits |>
    dplyr::mutate(
      is_rhis = as.integer(.data$method == "rhis"),
      is_salad = as.integer(.data$method == "salad"),
      log_effort = log1p(dplyr::coalesce(.data$effort_minutes, 0) +
                           dplyr::coalesce(.data$effort_m, 0) / 100)
    )
  obs_covs <- list(
    is_rhis = wide_observation_matrix(visits, "is_rhis", site_order),
    is_salad = wide_observation_matrix(visits, "is_salad", site_order),
    log_effort = wide_observation_matrix(visits, "log_effort", site_order)
  )
  umf <- make_unmarked_frame(y, site_covs, obs_covs)
  occupancy_formula <- paste(colnames(site_design), collapse = " + ")
  formula <- stats::as.formula(paste(
    "~ is_rhis + is_salad + log_effort ~", occupancy_formula
  ))
  fit <- unmarked::occu(formula, data = umf, method = "BFGS")
  saveRDS(list(
    fit = fit, pca = pca, site_design_center = design_center[keep],
    site_design_scale = design_scale[keep], sites = sites,
    site_terms = colnames(site_design), formula = formula,
    response_policy = list(
      included = c("Cull point", "RHIS point", "SALAD time-localized event/negative track"),
      excluded = "Manta positive lines: along-track detection location unknown"
    )
  ), file.path(
    cfg$paths$models_dir, paste0("multimethod_occupancy_", cell_size_m, "m.rds")
  ))

  coefficients <- unlist(lapply(names(fit@estimates@estimates), function(component) {
    value <- fit@estimates@estimates[[component]]@estimates
    stats::setNames(value, paste(component, names(value), sep = "::"))
  }))
  coef_table <- tibble::tibble(
    support_m = cell_size_m,
    parameter = names(coefficients), coefficient = as.numeric(coefficients)
  )
  fitted_psi <- unmarked::predict(fit, type = "state")
  audit <- tibble::tibble(
    support_m = cell_size_m, total_sites = nrow(sites), total_visits = nrow(visits),
    positive_sites = sum(sites$positive), repeated_sites = sum(sites$n_visits > 1),
    converged = isTRUE(fit@opt$convergence == 0),
    negative_log_likelihood = fit@opt$value,
    psi_min = min(fitted_psi$Predicted), psi_median = stats::median(fitted_psi$Predicted),
    psi_max = max(fitted_psi$Predicted)
  )
  list(coef = coef_table, audit = audit)
}

result <- lapply(cfg$validation$cell_support_m, fit_support)
coefs <- dplyr::bind_rows(lapply(result, `[[`, "coef"))
audit <- dplyr::bind_rows(lapply(result, `[[`, "audit"))
utils::write.csv(coefs, file.path(
  cfg$paths$outputs_dir, "multimethod_occupancy_coefficients.csv"
), row.names = FALSE)
utils::write.csv(audit, file.path(
  cfg$paths$outputs_dir, "multimethod_occupancy_audit.csv"
), row.names = FALSE)
print(audit)
print(dplyr::filter(coefs, grepl("is_rhis|is_salad|log_effort", .data$parameter)))
