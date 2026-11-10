args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v3_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v3_dir, "R", "bootstrap.R"))
ctx <- initialise_v3(script_file)
cfg <- ctx$config

visits <- readRDS(file.path(
  cfg$paths$artifacts_dir, "cell_visits_120m_window_3d.rds"
))

direct_dates <- visits |>
  dplyr::filter(
    .data$direct_localized, .data$direct_detection == 1L,
    .data$method %in% c("cull", "rhis", "salad")
  ) |>
  dplyr::distinct(.data$cell_id, .data$survey_date) |>
  dplyr::arrange(.data$cell_id, .data$survey_date)

scar <- visits |>
  dplyr::filter(
    .data$method %in% c("rhis", "salad"), !is.na(.data$scar_detection),
    !is.na(.data$reef_key), .data$reef_key != ""
  ) |>
  dplyr::group_by(.data$cell_id, .data$survey_date, .data$method) |>
  dplyr::summarise(
    scar_detection = as.integer(any(.data$scar_detection == 1L)),
    reef_key = dplyr::first(.data$reef_key),
    cull_timing = dplyr::first(.data$cull_timing), .groups = "drop"
  )

date_lookup <- split(direct_dates$survey_date, direct_dates$cell_id)
scar$days_since_prior_direct <- vapply(seq_len(nrow(scar)), function(i) {
  dates <- date_lookup[[scar$cell_id[[i]]]]
  if (is.null(dates)) return(NA_real_)
  prior <- dates[dates < scar$survey_date[[i]]]
  if (!length(prior)) return(NA_real_)
  as.numeric(scar$survey_date[[i]] - max(prior))
}, numeric(1))
scar <- scar |>
  dplyr::mutate(
    prior_direct_known = as.integer(is.finite(.data$days_since_prior_direct)),
    log_days_since_direct = log1p(dplyr::coalesce(.data$days_since_prior_direct, 0)),
    post_cull = as.integer(.data$cull_timing == "post_cull")
  )
scar_fit <- glmmTMB::glmmTMB(
  scar_detection ~ method + prior_direct_known + log_days_since_direct +
    post_cull + (1 | reef_key),
  data = scar, family = stats::binomial()
)
saveRDS(scar_fit, file.path(cfg$paths$models_dir, "scar_persistence_120m.rds"))
scar_coef <- as.data.frame(summary(scar_fit)$coefficients$cond) |>
  tibble::rownames_to_column("parameter") |>
  dplyr::mutate(process = "scar_persistence")

cull <- visits |>
  dplyr::filter(
    .data$method == "cull", is.finite(.data$direct_count),
    is.finite(.data$effort_minutes), .data$effort_minutes > 0,
    !is.na(.data$reef_key), .data$reef_key != ""
  ) |>
  dplyr::distinct(.data$visit_id, .keep_all = TRUE) |>
  dplyr::mutate(
    post_first_cull = as.integer(.data$days_from_local_first_cull > 0),
    days_after_first = pmax(dplyr::coalesce(.data$days_from_local_first_cull, 0), 0),
    log_days_after_first = log1p(.data$days_after_first),
    log_effort = log(pmax(.data$effort_minutes, 1))
  )
cull_fit <- glmmTMB::glmmTMB(
  direct_count ~ post_first_cull + log_days_after_first +
    (1 | reef_key) + (1 | cell_id),
  offset = log_effort, data = cull, family = glmmTMB::nbinom2()
)
saveRDS(cull_fit, file.path(cfg$paths$models_dir, "cull_depletion_120m.rds"))
cull_coef <- as.data.frame(summary(cull_fit)$coefficients$cond) |>
  tibble::rownames_to_column("parameter") |>
  dplyr::mutate(process = "cull_depletion")

coefficients <- dplyr::bind_rows(scar_coef, cull_coef)
utils::write.csv(coefficients, file.path(
  cfg$paths$outputs_dir, "observation_process_coefficients.csv"
), row.names = FALSE)
audit <- tibble::tibble(
  process = c("scar_persistence", "cull_depletion"),
  rows = c(nrow(scar), nrow(cull)),
  reefs = c(dplyr::n_distinct(scar$reef_key), dplyr::n_distinct(cull$reef_key)),
  cells = c(dplyr::n_distinct(scar$cell_id), dplyr::n_distinct(cull$cell_id)),
  response_events = c(sum(scar$scar_detection), sum(cull$direct_count)),
  converged = c(
    isTRUE(scar_fit$fit$convergence == 0), isTRUE(cull_fit$fit$convergence == 0)
  )
)
utils::write.csv(audit, file.path(
  cfg$paths$outputs_dir, "observation_process_audit.csv"
), row.names = FALSE)
print(audit)
print(dplyr::select(coefficients, "process", "parameter", "Estimate", "Std. Error", "Pr(>|z|)"))
