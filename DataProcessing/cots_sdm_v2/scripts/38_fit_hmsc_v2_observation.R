args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v2_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v2_dir, "R", "bootstrap.R"))
ctx <- initialise_v2(script_file)
cfg <- ctx$config

method_path <- file.path(
  cfg$paths$artifacts_dir, "hmsc_v2_annotated_method_tables.rds"
)
if (!file.exists(method_path)) stop("Run script 35 before script 38.")
method_tables <- readRDS(method_path)

to_binary <- function(x) as.integer(dplyr::coalesce(x, 0) > 0)
common_fields <- c(
  "visit_id", "method", "survey_date", "year", "reef_key", "reef_name",
  "cell_id", "visit_cluster_id", "cull_timing", "local_first_cull_date",
  "days_from_local_first_cull"
)

cull_records <- method_tables$cull |>
  dplyr::transmute(
    dplyr::across(dplyr::all_of(common_fields)),
    detected_direct = to_binary(.data$direct_cots_count),
    detected_scar = NA_integer_,
    scar_observed = FALSE,
    effort = as.numeric(.data$bottom_time),
    count = as.numeric(.data$direct_cots_count),
    cull_dive_type = dplyr::coalesce(.data$cull_dive_type, "missing")
  )
manta_records <- method_tables$manta |>
  dplyr::transmute(
    dplyr::across(dplyr::all_of(common_fields)),
    detected_direct = to_binary(.data$direct_cots_count),
    detected_scar = as.integer(
      dplyr::coalesce(.data$scar_count > 0, FALSE) |
        dplyr::coalesce(.data$scar_code %in% c("p", "c"), FALSE)
    ),
    scar_observed = TRUE,
    effort = as.numeric(.data$tow_distance_m),
    count = as.numeric(.data$direct_cots_count),
    cull_dive_type = NA_character_
  )
rhis_records <- method_tables$rhis |>
  dplyr::filter(.data$rhis_cots_protocol_observed) |>
  dplyr::transmute(
    dplyr::across(dplyr::all_of(common_fields)),
    detected_direct = as.integer(
      dplyr::coalesce(.data$adult_count > 0, FALSE) |
        dplyr::coalesce(.data$juvenile_count > 0, FALSE)
    ),
    detected_scar = as.integer(
      dplyr::coalesce(.data$branching_scar_count > 0, FALSE) |
        dplyr::coalesce(.data$table_scar_count > 0, FALSE)
    ),
    scar_observed = TRUE,
    effort = NA_real_,
    count = rowSums(cbind(
      dplyr::coalesce(.data$adult_count, 0),
      dplyr::coalesce(.data$juvenile_count, 0)
    )),
    cull_dive_type = NA_character_
  )
salad_records <- method_tables$salad |>
  dplyr::transmute(
    dplyr::across(dplyr::all_of(common_fields)),
    detected_direct = to_binary(.data$direct_cots_count),
    detected_scar = to_binary(.data$scar_count),
    scar_observed = TRUE,
    effort = as.numeric(.data$footprint_area_m2),
    count = as.numeric(.data$direct_cots_count),
    cull_dive_type = NA_character_
  )

records <- dplyr::bind_rows(
  cull_records, manta_records, rhis_records, salad_records
) |>
  dplyr::filter(!is.na(.data$cell_id), !is.na(.data$visit_cluster_id)) |>
  dplyr::mutate(
    cluster_cell_id = paste(.data$visit_cluster_id, .data$cell_id, sep = "||"),
    effort_positive = is.finite(.data$effort) & .data$effort > 0,
    log_effort = dplyr::if_else(
      .data$effort_positive, log1p(.data$effort), NA_real_
    )
  )

method_units <- records |>
  dplyr::group_by(
    .data$cluster_cell_id, .data$method, .data$year, .data$reef_key,
    .data$cull_timing
  ) |>
  dplyr::summarise(
    detected_direct = as.integer(any(.data$detected_direct == 1L)),
    detected_scar = if (any(.data$scar_observed)) {
      as.integer(any(.data$detected_scar == 1L, na.rm = TRUE))
    } else NA_integer_,
    effort = sum(.data$effort, na.rm = TRUE),
    effort_observed = any(.data$effort_positive),
    records = dplyr::n(),
    days_from_local_first_cull = stats::median(
      .data$days_from_local_first_cull, na.rm = TRUE
    ),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    days_from_local_first_cull = dplyr::if_else(
      is.nan(.data$days_from_local_first_cull), NA_real_,
      .data$days_from_local_first_cull
    )
  )

group_evidence <- method_units |>
  dplyr::group_by(.data$cluster_cell_id) |>
  dplyr::summarise(
    methods_present = dplyr::n_distinct(.data$method),
    any_direct = any(.data$detected_direct == 1L),
    any_scar = any(.data$detected_scar == 1L, na.rm = TRUE),
    .groups = "drop"
  )
cross_method <- method_units |>
  dplyr::left_join(group_evidence, by = "cluster_cell_id") |>
  dplyr::group_by(.data$cluster_cell_id) |>
  dplyr::mutate(
    other_method_direct = vapply(seq_len(dplyr::n()), function(i) {
      any(.data$detected_direct[-i] == 1L & .data$method[-i] != .data$method[i])
    }, logical(1))
  ) |>
  dplyr::ungroup()

detection_summary <- cross_method |>
  dplyr::filter(.data$methods_present > 1L, .data$other_method_direct) |>
  dplyr::group_by(.data$method) |>
  dplyr::summarise(
    confirmed_footprints = dplyr::n(),
    direct_detections = sum(.data$detected_direct),
    apparent_detection = mean(.data$detected_direct),
    .groups = "drop"
  )

scar_summary <- cross_method |>
  dplyr::filter(!is.na(.data$detected_scar)) |>
  dplyr::group_by(.data$method) |>
  dplyr::summarise(
    scar_observed_footprints = dplyr::n(),
    scar_positive = sum(.data$detected_scar),
    scar_with_any_direct = sum(.data$detected_scar == 1L & .data$any_direct),
    scar_without_same_cluster_direct = sum(
      .data$detected_scar == 1L & !.data$any_direct
    ),
    scar_direct_discordance = .data$scar_without_same_cluster_direct /
      pmax(.data$scar_positive, 1L),
    .groups = "drop"
  )

depletion_data <- method_units |>
  dplyr::filter(
    .data$cull_timing %in% c("strictly_pre_cull", "post_cull")
  ) |>
  dplyr::mutate(
    post_cull = as.integer(.data$cull_timing == "post_cull"),
    log_effort = dplyr::if_else(
      .data$effort_observed & .data$effort > 0, log1p(.data$effort), 0
    ),
    effort_observed = as.integer(.data$effort_observed),
    method = stats::relevel(factor(.data$method), ref = "manta")
  )

depletion_model <- stats::glm(
  detected_direct ~ method * post_cull + effort_observed + log_effort,
  family = stats::binomial(), data = depletion_data
)

cull_depletion_data <- records |>
  dplyr::filter(
    .data$method == "cull", .data$cull_timing %in% c("same_day_ambiguous", "post_cull"),
    is.finite(.data$count), .data$effort_positive
  ) |>
  dplyr::mutate(
    post_cull = as.integer(.data$cull_timing == "post_cull"),
    log_effort = log(.data$effort),
    log_days_after_first = log1p(pmax(.data$days_from_local_first_cull, 0)),
    cull_dive_type = factor(.data$cull_dive_type)
  )
cull_count_model <- stats::glm(
  count ~ log_effort + post_cull + log_days_after_first + cull_dive_type,
  family = stats::quasipoisson(), data = cull_depletion_data
)

tidy_glm <- function(model, name) {
  coefficient_table <- as.data.frame(summary(model)$coefficients)
  coefficient_table$term <- rownames(coefficient_table)
  rownames(coefficient_table) <- NULL
  names(coefficient_table)[seq_len(4)] <- c(
    "estimate", "standard_error", "statistic", "p_value"
  )
  tibble::as_tibble(coefficient_table) |>
    dplyr::mutate(model = name, .before = 1L)
}
coefficient_table <- dplyr::bind_rows(
  tidy_glm(depletion_model, "method_detection_pre_vs_post_cull"),
  tidy_glm(cull_count_model, "cull_count_depletion")
)

record_count <- nrow(records)
cluster_cell_count <- dplyr::n_distinct(records$cluster_cell_id)
cross_method_cluster_count <- sum(group_evidence$methods_present > 1L)
cross_method_confirmed_count <- sum(cross_method$other_method_direct)
same_day_ambiguous_count <- sum(
  records$cull_timing == "same_day_ambiguous"
)
pre_cull_count <- sum(records$cull_timing == "strictly_pre_cull")
post_cull_count <- sum(records$cull_timing == "post_cull")
no_local_cull_count <- sum(records$cull_timing == "no_local_cull_record")
audit <- tibble::tibble(
  records = record_count,
  cluster_cells = cluster_cell_count,
  cross_method_cluster_cells = cross_method_cluster_count,
  cross_method_confirmed_records = cross_method_confirmed_count,
  same_day_ambiguous_records = same_day_ambiguous_count,
  pre_cull_records = pre_cull_count,
  post_cull_records = post_cull_count,
  no_local_cull_records = no_local_cull_count,
  deployable_habitat_predictors_used = FALSE,
  interpretation = paste(
    "Observation/depletion diagnostic only; cross-method confirmation is an",
    "imperfect reference and scar-direct discordance is not a false-positive rate."
  )
)

saveRDS(
  list(
    method_units = method_units, cross_method = cross_method,
    depletion_model = depletion_model, cull_count_model = cull_count_model
  ),
  file.path(cfg$paths$models_dir, "hmsc_v2_observation_depletion_models.rds")
)
tables <- list(
  hmsc_v2_observation_audit = audit,
  hmsc_v2_cross_method_detection = detection_summary,
  hmsc_v2_scar_direct_discordance = scar_summary,
  hmsc_v2_observation_coefficients = coefficient_table
)
for (name in names(tables)) {
  utils::write.csv(
    tables[[name]], file.path(cfg$paths$outputs_dir, paste0(name, ".csv")),
    row.names = FALSE
  )
}
message("Hmsc v2 observation and depletion diagnostics complete.")
print(audit)
print(detection_summary)
print(scar_summary)
print(coefficient_table)
