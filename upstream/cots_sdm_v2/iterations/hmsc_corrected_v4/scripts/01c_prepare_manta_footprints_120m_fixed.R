args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v4_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v4_dir, "R", "bootstrap.R"))
ctx <- initialise_v4(script_file)
cfg <- ctx$config

links <- readRDS(file.path(cfg$paths$v3_artifacts_dir, "footprint_links_120m.rds"))
predictors <- readRDS(file.path(cfg$paths$v3_artifacts_dir, "support_predictors_120m.rds"))
activity <- readRDS(file.path(cfg$paths$artifacts_dir, "reef_year_activity_120m.rds")) |>
  dplyr::select("reef_key", "year", "active_primary", "active_weak_sensitivity")
manta <- links$manta |>
  classify_manta_evidence(cfg) |>
  dplyr::left_join(predictors, by = c("cell_id", "projected_x", "projected_y")) |>
  dplyr::mutate(
    predictor_weight = normalise_manta_weights_v4(.data$weighted_length_m, .data$fraction)
  )

numeric_predictors <- names(predictors)[vapply(predictors, is.numeric, logical(1))]
numeric_predictors <- setdiff(numeric_predictors, c("projected_x", "projected_y", "cell_support_m"))
value_names <- c("projected_x", "projected_y", numeric_predictors)
group_levels <- unique(manta$visit_id)
group_index <- match(manta$visit_id, group_levels)
weight <- manta$predictor_weight
value_matrix <- as.matrix(manta[, value_names, drop = FALSE])
finite <- is.finite(value_matrix)
value_matrix[!finite] <- 0
numerator <- rowsum(value_matrix * weight, group_index, reorder = FALSE)
denominator <- rowsum(finite * weight, group_index, reorder = FALSE)
weighted_values <- numerator / denominator
weighted_values[denominator <= 0] <- NA_real_
colnames(weighted_values) <- value_names

first_row <- match(group_levels, manta$visit_id)
metadata <- manta[first_row, c(
  "visit_id", "survey_date", "year", "reef_key", "reef_name",
  "manta_strong_positive", "manta_weak_scar", "manta_strict_negative",
  "manta_evidence_class"
), drop = FALSE]
dominant_geomorph <- manta |>
  dplyr::filter(!is.na(.data$geomorphic_class)) |>
  dplyr::arrange(.data$visit_id, dplyr::desc(.data$predictor_weight)) |>
  dplyr::distinct(.data$visit_id, .keep_all = TRUE) |>
  dplyr::select("visit_id", "geomorphic_class")
footprint_length <- as.numeric(rowsum(weight, group_index, reorder = FALSE))

manta_footprints <- dplyr::bind_cols(
  tibble::as_tibble(metadata), tibble::as_tibble(weighted_values)
) |>
  dplyr::mutate(footprint_length_m = footprint_length) |>
  dplyr::left_join(dominant_geomorph, by = "visit_id") |>
  dplyr::left_join(activity, by = c("reef_key", "year")) |>
  dplyr::mutate(
    active_primary = dplyr::coalesce(.data$active_primary, FALSE),
    active_weak_sensitivity = dplyr::coalesce(.data$active_weak_sensitivity, FALSE),
    reef_year_id = paste(.data$reef_key, .data$year, sep = "||"),
    manta_binary = dplyr::case_when(
      .data$manta_strong_positive ~ 1L,
      .data$manta_strict_negative ~ 0L,
      TRUE ~ NA_integer_
    ),
    observation_support = "manta_tow_footprint"
  )
raw_terms <- unlist(terrain_families(names(manta_footprints)), use.names = FALSE)
clean_terms <- c("AhyaD", "Aspat", "Aten", raw_terms)
enhanced_terms <- c("geomorphic_class", clean_terms, "turbid", "waves_Hs", "waves_Tp")
manta_footprints <- manta_footprints |>
  dplyr::mutate(
    clean_complete = stats::complete.cases(dplyr::pick(dplyr::all_of(clean_terms))),
    enhanced_complete = stats::complete.cases(dplyr::pick(dplyr::all_of(enhanced_terms)))
  )
saveRDS(manta_footprints, file.path(cfg$paths$artifacts_dir, "manta_footprints_120m.rds"))
evidence_audit <- manta_footprints |>
  dplyr::count(.data$manta_evidence_class, name = "tows") |>
  dplyr::arrange(dplyr::desc(.data$tows))
coverage_audit <- tibble::tibble(
  footprints = nrow(manta_footprints),
  primary_active_footprints = sum(manta_footprints$active_primary),
  strong_positive = sum(manta_footprints$manta_strong_positive),
  weak_scar = sum(manta_footprints$manta_weak_scar),
  strict_negative = sum(manta_footprints$manta_strict_negative),
  primary_binary_eligible = sum(!is.na(manta_footprints$manta_binary) & manta_footprints$active_primary),
  clean_complete = sum(manta_footprints$clean_complete),
  enhanced_complete = sum(manta_footprints$enhanced_complete)
)
utils::write.csv(evidence_audit, file.path(cfg$paths$outputs_dir, "manta_evidence_audit.csv"), row.names = FALSE)
utils::write.csv(coverage_audit, file.path(cfg$paths$outputs_dir, "manta_footprint_audit_120m.csv"), row.names = FALSE)
print(evidence_audit)
print(coverage_audit)
