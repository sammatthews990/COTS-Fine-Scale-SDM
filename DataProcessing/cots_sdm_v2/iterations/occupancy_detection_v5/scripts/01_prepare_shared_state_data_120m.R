args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v5_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
repo_root <- (function(p) { while (!dir.exists(file.path(p, ".git"))) { if (dirname(p) == p) stop("Repository root not found"); p <- dirname(p) }; p })(v5_dir)
v4_dir <- normalizePath(file.path(v5_dir, "..", "hmsc_corrected_v4"), winslash = "/")
v3_dir <- normalizePath(file.path(v5_dir, "..", "hmsc_geomorphic_v3"), winslash = "/")
source(file.path(v4_dir, "R", "bootstrap.R"))
ctx <- initialise_v4(file.path(v4_dir, "scripts", "01a_prepare_corrected_cells_120m.R"))
cfg <- ctx$config
set.seed(20260910L)

artifact_dir <- file.path(v5_dir, "artifacts")
output_dir <- file.path(v5_dir, "outputs")
dir.create(artifact_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

visits <- readRDS(file.path(v3_dir, "artifacts", "cell_visits_120m_window_3d.rds"))
links <- readRDS(file.path(v3_dir, "artifacts", "footprint_links_120m.rds"))
predictors <- readRDS(file.path(v3_dir, "artifacts", "support_predictors_120m.rds"))
activity <- readRDS(file.path(v4_dir, "artifacts", "reef_year_activity_120m.rds")) |>
  dplyr::select("reef_key", "year", "active_primary")

raw_terms <- unlist(terrain_families(names(predictors)), use.names = FALSE)
state_base_terms <- c(
  "geomorphic_class", "AhyaD", "Aspat", "Aten",
  "turbid", "waves_Hs", "waves_Tp"
)
required <- unique(c(state_base_terms, raw_terms))

point <- visits |>
  dplyr::filter(
    .data$ecological_eligible,
    .data$method %in% c("cull", "rhis", "salad"),
    !is.na(.data$reef_key), .data$reef_key != ""
  ) |>
  dplyr::left_join(activity, by = c("reef_key", "year")) |>
  dplyr::filter(dplyr::coalesce(.data$active_primary, FALSE)) |>
  dplyr::group_by(
    .data$cell_id, .data$year, .data$survey_date, .data$method,
    .data$reef_key, .data$reef_name
  ) |>
  dplyr::summarise(
    projected_x = dplyr::first(.data$projected_x),
    projected_y = dplyr::first(.data$projected_y),
    direct_detection = observed_any_v4(.data$direct_detection),
    scar_detection = observed_any_v4(.data$scar_detection),
    effort_m = sum(.data$effort_m, na.rm = TRUE),
    effort_minutes = sum(.data$effort_minutes, na.rm = TRUE),
    post_cull = as.integer(any(.data$cull_timing == "post_cull", na.rm = TRUE)),
    days_after_cull = max(pmax(dplyr::coalesce(.data$days_from_local_first_cull, 0), 0), na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(site_id = paste(.data$cell_id, .data$year, sep = "||"))

manta_visit <- links$manta |>
  dplyr::distinct(.data$visit_id, .keep_all = TRUE) |>
  classify_manta_evidence(cfg) |>
  dplyr::left_join(activity, by = c("reef_key", "year")) |>
  dplyr::filter(
    dplyr::coalesce(.data$active_primary, FALSE),
    .data$manta_strong_positive | .data$manta_strict_negative
  ) |>
  dplyr::mutate(
    manta_group_type = dplyr::case_when(
      .data$manta_direct_positive ~ 1L,
      .data$manta_multiple_scar_positive | .data$manta_common_scar_positive ~ 2L,
      TRUE ~ 0L
    )
  ) |>
  dplyr::group_by(.data$reef_key, .data$year) |>
  dplyr::mutate(n_positive = sum(.data$manta_group_type > 0L)) |>
  dplyr::arrange(.data$reef_key, .data$year, dplyr::desc(.data$manta_group_type), .data$visit_id) |>
  dplyr::filter(
    .data$manta_group_type > 0L |
      dplyr::row_number() <= .data$n_positive + dplyr::if_else(.data$n_positive > 0L, .data$n_positive, 5L)
  ) |>
  dplyr::ungroup()

manta_links <- links$manta |>
  dplyr::semi_join(manta_visit, by = "visit_id") |>
  dplyr::left_join(
    dplyr::select(manta_visit, "visit_id", "manta_group_type"),
    by = "visit_id"
  ) |>
  dplyr::mutate(site_id = paste(.data$cell_id, .data$year, sep = "||"))

site_seed <- dplyr::bind_rows(
  point |>
    dplyr::select("site_id", "cell_id", "year", "reef_key", "reef_name", "projected_x", "projected_y"),
  manta_links |>
    dplyr::select("site_id", "cell_id", "year", "reef_key", "reef_name", "projected_x", "projected_y")
) |>
  dplyr::distinct(.data$site_id, .keep_all = TRUE) |>
  dplyr::left_join(predictors, by = c("cell_id", "projected_x", "projected_y"))
sites <- site_seed[stats::complete.cases(site_seed[, required, drop = FALSE]), , drop = FALSE]
point <- dplyr::semi_join(point, sites, by = "site_id")
manta_links <- dplyr::semi_join(manta_links, sites, by = "site_id")
manta_visit <- dplyr::semi_join(manta_visit, manta_links, by = "visit_id")

pca <- fit_family_pca(sites, variance = 0, max_components = 1L)
pcs <- apply_family_pca(sites, pca)
state_data <- dplyr::bind_cols(
  dplyr::select(sites, dplyr::all_of(state_base_terms)), pcs
)
state_data$geomorphic_class <- factor(state_data$geomorphic_class)
X <- stats::model.matrix(
  ~ geomorphic_class + AhyaD + Aspat + Aten + turbid + waves_Hs + waves_Tp + .,
  data = dplyr::select(state_data, -"geomorphic_class") |>
    dplyr::mutate(geomorphic_class = state_data$geomorphic_class, .before = 1)
)
centre <- rep(0, ncol(X)); names(centre) <- colnames(X)
scale_value <- rep(1, ncol(X)); names(scale_value) <- colnames(X)
continuous <- setdiff(colnames(X), c("(Intercept)", grep("^geomorphic_class", colnames(X), value = TRUE)))
centre[continuous] <- colMeans(X[, continuous, drop = FALSE])
scale_value[continuous] <- apply(X[, continuous, drop = FALSE], 2, stats::sd)
valid_columns <- is.finite(scale_value) & scale_value > 1e-10
X <- X[, valid_columns, drop = FALSE]
centre <- centre[valid_columns]
scale_value <- scale_value[valid_columns]
X <- sweep(X, 2, centre, "-")
X <- sweep(X, 2, scale_value, "/")
site_index <- stats::setNames(seq_len(nrow(sites)) - 1L, sites$site_id)

direct <- point |>
  dplyr::filter(!is.na(.data$direct_detection)) |>
  dplyr::mutate(
    site = unname(site_index[.data$site_id]),
    method_index = match(.data$method, c("cull", "rhis", "salad")) - 1L,
    log_effort = log1p(.data$effort_minutes + .data$effort_m / 100),
    log_days_after_cull = log1p(.data$days_after_cull)
  )

positive_dates <- direct |>
  dplyr::filter(.data$direct_detection == 1L) |>
  dplyr::distinct(.data$cell_id, .data$survey_date) |>
  dplyr::arrange(.data$cell_id, .data$survey_date)
date_lookup <- split(positive_dates$survey_date, positive_dates$cell_id)
scar <- point |>
  dplyr::filter(.data$method %in% c("rhis", "salad"), !is.na(.data$scar_detection)) |>
  dplyr::mutate(
    site = unname(site_index[.data$site_id]),
    method_index = match(.data$method, c("rhis", "salad")) - 1L
  )
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
    log_days_since_direct = log1p(dplyr::coalesce(.data$days_since_prior_direct, 0))
  )

manta_links <- manta_links |>
  dplyr::mutate(
    site = unname(site_index[.data$site_id]),
    log_effort = log1p(pmax(.data$weighted_length_m, 1) / 100)
  ) |>
  dplyr::arrange(.data$visit_id, .data$cell_id)
group_table <- manta_visit |>
  dplyr::semi_join(manta_links, by = "visit_id") |>
  dplyr::arrange(match(.data$visit_id, unique(manta_links$visit_id))) |>
  dplyr::distinct(.data$visit_id, .keep_all = TRUE)
manta_links$group_index <- match(manta_links$visit_id, group_table$visit_id)
manta_links <- manta_links |>
  dplyr::filter(!is.na(.data$group_index)) |>
  dplyr::arrange(.data$group_index)
group_counts <- tabulate(manta_links$group_index, nbins = nrow(group_table))
group_start <- c(0L, cumsum(group_counts))
point_sites <- unique(c(direct$site, scar$site))
group_weight <- rep(min(1, length(point_sites) / max(nrow(group_table), 1)), nrow(group_table))

has_point <- integer(nrow(sites))
has_point[point_sites + 1L] <- 1L
tmb_data <- list(
  X = unname(X), has_point_observation = has_point,
  direct_site = as.integer(direct$site),
  direct_method = as.integer(direct$method_index),
  direct_y = as.integer(direct$direct_detection),
  direct_log_effort = as.numeric(direct$log_effort),
  direct_post_cull = as.numeric(direct$post_cull),
  direct_log_days_after_cull = as.numeric(direct$log_days_after_cull),
  scar_site = as.integer(scar$site),
  scar_method = as.integer(scar$method_index),
  scar_y = as.integer(scar$scar_detection),
  scar_prior_direct_known = as.numeric(scar$prior_direct_known),
  scar_log_days_since_direct = as.numeric(scar$log_days_since_direct),
  scar_post_cull = as.numeric(scar$post_cull),
  manta_group_start = as.integer(group_start),
  manta_group_type = as.integer(group_table$manta_group_type),
  manta_link_site = as.integer(manta_links$site),
  manta_link_log_effort = as.numeric(manta_links$log_effort),
  manta_group_weight = as.numeric(group_weight)
)

artifact <- list(
  data = tmb_data, sites = sites, point = point, direct = direct, scar = scar,
  manta_groups = group_table, manta_links = manta_links,
  pca = pca, state_base_terms = state_base_terms,
  design_columns = colnames(X), design_center = centre, design_scale = scale_value,
  support_m = 120L,
  likelihood_note = "Exact cell-year occupancy likelihood for point repeats plus a weighted marginal track-level composite likelihood for Manta OR detections."
)
saveRDS(artifact, file.path(artifact_dir, "shared_state_data_120m.rds"), compress = "xz")

audit <- tibble::tibble(
  support_m = 120L, state_sites = nrow(sites), point_site_years = length(point_sites),
  repeated_point_site_years = sum(table(c(direct$site, scar$site)) > 1L),
  direct_observations = nrow(direct), direct_positives = sum(direct$direct_detection),
  scar_observations = nrow(scar), scar_positives = sum(scar$scar_detection),
  manta_tows = nrow(group_table), manta_direct_positive = sum(group_table$manta_group_type == 1L),
  manta_scar_positive = sum(group_table$manta_group_type == 2L),
  manta_strict_negative = sum(group_table$manta_group_type == 0L),
  manta_cell_links = nrow(manta_links), manta_composite_weight = unique(group_weight)[1]
)
utils::write.csv(audit, file.path(output_dir, "shared_state_data_audit_120m.csv"), row.names = FALSE)
print(audit)
