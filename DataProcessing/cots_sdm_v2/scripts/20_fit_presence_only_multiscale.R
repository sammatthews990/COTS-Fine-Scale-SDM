args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v2_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v2_dir, "R", "bootstrap.R"))
ctx <- initialise_v2(script_file)
cfg <- ctx$config

required <- c(
  file.path(cfg$paths$artifacts_dir, "site_year_all.rds"),
  file.path(cfg$paths$artifacts_dir, "site_habitat_all_sites.rds"),
  file.path(cfg$paths$artifacts_dir, "site_multiscale_30_60_120m.rds"),
  file.path(cfg$paths$repo_root, "DataProcessing", "outputs", "maxent_bg_points_10k.gpkg")
)
if (any(!file.exists(required))) stop("Missing presence-only input(s).")

predictors <- terra::rast(cfg$paths$predictor_stack)
base_names <- names(predictors)
features <- maxent_feature_sets(base_names)
all_features <- unique(unlist(features))

site_year <- readRDS(required[1])
site_folds <- readRDS(required[2]) |>
  dplyr::filter(.data$evidence_source == "cull_site_history") |>
  dplyr::select("site_id", "site_fold")
site_scales <- readRDS(required[3]) |>
  dplyr::left_join(site_folds, by = "site_id")

cull_history <- site_year |>
  dplyr::group_by(.data$site_id) |>
  dplyr::summarise(
    ever_high_density = any(.data$high_density == 1L),
    ever_direct_cots = any(.data$total_cots > 0),
    .groups = "drop"
  ) |>
  dplyr::left_join(site_scales, by = "site_id") |>
  dplyr::mutate(source = "cull")

cache_path <- file.path(cfg$paths$artifacts_dir, "external_presence_multiscale.rds")
if (file.exists(cache_path)) {
  external <- readRDS(cache_path)
} else {
  message("Preparing Manta and ReefScan presence points.")
  manta <- read_manta_observations(cfg$paths$manta_history_csv) |>
    dplyr::filter(.data$valid_coordinates, .data$informative) |>
    dplyr::transmute(
      longitude, latitude,
      direct_cots = dplyr::coalesce(.data$cots_count > 0, FALSE),
      cots_or_scars = .data$positive_evidence,
      source = "manta"
    ) |>
    dplyr::filter(.data$cots_or_scars)

  reefscan_files <- list.files(
    file.path(cfg$paths$repo_root, "data", "ReefScan"),
    pattern = "[.]csv$", full.names = TRUE
  )
  reefscan <- dplyr::bind_rows(lapply(
    reefscan_files,
    function(path) readr::read_csv(path, show_col_types = FALSE) |>
      dplyr::transmute(
        longitude = as.numeric(.data$longitude),
        latitude = as.numeric(.data$latitude),
        direct_cots = TRUE,
        cots_or_scars = TRUE,
        source = "reefscan"
      )
  ))
  external_raw <- dplyr::bind_rows(manta, reefscan)
  external_raw <- external_raw |>
    dplyr::filter(
      dplyr::between(.data$longitude, 140, 160),
      dplyr::between(.data$latitude, -30, -5)
    )
  membership <- aligned_grid_membership(external_raw, predictors, 1L)
  external_raw$cell_s30m <- membership$coarse_cell
  external_raw <- external_raw |>
    dplyr::filter(is.finite(.data$cell_s30m)) |>
    dplyr::group_by(.data$cell_s30m) |>
    dplyr::summarise(
      longitude = mean(.data$longitude),
      latitude = mean(.data$latitude),
      direct_cots = any(.data$direct_cots),
      cots_or_scars = any(.data$cots_or_scars),
      source = paste(sort(unique(.data$source)), collapse = "+"),
      .groups = "drop"
    )
  external <- attach_multiscale_block_means(
    dplyr::select(external_raw, -"cell_s30m"), predictors
  ) |>
    dplyr::mutate(
      direct_cots = external_raw$direct_cots,
      cots_or_scars = external_raw$cots_or_scars,
      source = external_raw$source
    )
  saveRDS(external, cache_path)
}

background_sf <- sf::st_read(required[4], quiet = TRUE) |>
  sf::st_transform(4326)
background_xy <- sf::st_coordinates(background_sf)
background_points <- tibble::tibble(
  longitude = background_xy[, 1], latitude = background_xy[, 2]
)
background_cache <- file.path(cfg$paths$artifacts_dir, "maxent_background_multiscale.rds")
if (file.exists(background_cache)) {
  background <- readRDS(background_cache)
} else {
  background <- attach_multiscale_block_means(background_points, predictors)
  saveRDS(background, background_cache)
}

thin_external <- function(data, grid_m = 500) {
  data |>
    dplyr::mutate(
      external_x_bin = floor(.data$projected_x / grid_m),
      external_y_bin = floor(.data$projected_y / grid_m)
    ) |>
    dplyr::arrange(.data$source, .data$cell_s30m) |>
    dplyr::distinct(.data$external_x_bin, .data$external_y_bin, .keep_all = TRUE)
}
external_direct_thinned <- external |>
  dplyr::filter(.data$direct_cots) |>
  thin_external()
external_any_thinned <- thin_external(external)

presence_sets <- list(
  high_density_history = dplyr::filter(cull_history, .data$ever_high_density),
  direct_cots_history = dplyr::bind_rows(
    dplyr::filter(cull_history, .data$ever_direct_cots),
    external_direct_thinned
  ),
  cots_or_scars_history = dplyr::bind_rows(
    dplyr::filter(cull_history, .data$ever_direct_cots),
    external_any_thinned
  )
)
presence_sets <- lapply(presence_sets, deduplicate_presence_cells)

evaluation <- site_year |>
  dplyr::select(
    "site_year_id", "site_id", "reef_name", "year", "high_density",
    "longitude", "latitude"
  ) |>
  dplyr::left_join(site_scales, by = c("site_id", "longitude", "latitude")) |>
  dplyr::mutate(reef_year_id = paste(.data$reef_name, .data$year, sep = "||")) |>
  dplyr::filter(stats::complete.cases(dplyr::pick(dplyr::all_of(all_features))))

predictions <- list()
models <- list()
audit <- list()
counter <- 0L
for (definition in names(presence_sets)) {
  presence <- presence_sets[[definition]]
  for (scale_name in names(features)) {
    feature_names <- features[[scale_name]]
    for (fold in 0L) {
      counter <- counter + 1L
      test <- evaluation
      train_presence <- presence |>
        dplyr::filter(
          stats::complete.cases(dplyr::pick(dplyr::all_of(feature_names)))
        )
      train_background <- background |>
        dplyr::filter(
          !.data$cell_s30m %in% train_presence$cell_s30m,
          stats::complete.cases(dplyr::pick(dplyr::all_of(feature_names)))
        )
      fit <- fit_maxnet_lq(
        train_presence, train_background, feature_names,
        regmult = 2, classes = "l"
      )
      probability <- predict_cots_maxnet(fit, test)
      model_name <- paste("maxent_linear", definition, scale_name, sep = "__")
      predictions[[counter]] <- tibble::tibble(
        model = model_name,
        validation = "descriptive_full_fit_not_held_out",
        fold = fold,
        site_year_id = test$site_year_id,
        site_id = test$site_id,
        reef_name = test$reef_name,
        reef_year_id = test$reef_year_id,
        year = test$year,
        truth = test$high_density,
        probability = probability,
        longitude = test$longitude,
        latitude = test$latitude
      )
      models[[paste0(model_name, "__fold_", fold)]] <- fit
      audit[[counter]] <- tibble::tibble(
        model = model_name, fold = fold,
        training_presences = nrow(train_presence),
        training_background = nrow(train_background),
        cull_presences = sum(train_presence$source == "cull"),
        external_presences = sum(train_presence$source != "cull")
      )
    }
  }
}

predictions <- dplyr::bind_rows(predictions)
metrics <- suppressWarnings(summarise_benchmark_predictions(predictions))
decomposition <- suppressWarnings(
  prediction_scale_decomposition(predictions, group = "reef_year_id")
)
capture <- within_group_capture_metrics(predictions, group = "reef_year_id")
tables <- list(
  maxent_multiscale_metrics = metrics,
  maxent_multiscale_decomposition = decomposition,
  maxent_multiscale_capture = capture,
  maxent_multiscale_training_audit = dplyr::bind_rows(audit)
)
for (name in names(tables)) {
  utils::write.csv(tables[[name]], file.path(cfg$paths$outputs_dir, paste0(name, ".csv")), row.names = FALSE)
}
saveRDS(predictions, file.path(cfg$paths$outputs_dir, "maxent_multiscale_predictions.rds"))
saveRDS(models, file.path(cfg$paths$models_dir, "maxent_multiscale_fold_models.rds"))
message("Presence-only multiscale benchmark complete.")
print(dplyr::left_join(metrics, dplyr::select(decomposition, "model", "within_reef_roc_auc"), by = "model"))
