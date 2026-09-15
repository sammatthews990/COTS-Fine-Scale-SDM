args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v2_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v2_dir, "R", "bootstrap.R"))
ctx <- initialise_v2(script_file)
cfg <- ctx$config

required_packages <- c("terra", "xgboost", "png", "jsonlite", "base64enc")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages) > 0L) {
  stop("Missing blind-review packages: ", paste(missing_packages, collapse = ", "))
}

visual_dir <- file.path(cfg$paths$outputs_dir, "visual_comparison")
template_path <- file.path(
  cfg$paths$v2_dir, "templates", "expert_blind_review.html"
)
required <- c(
  cfg$paths$predictor_stack,
  cfg$paths$predictor_stack_rg,
  cfg$paths$legacy_2026_raster,
  cfg$paths$legacy_agnostic_raster,
  file.path(cfg$paths$artifacts_dir, "site_year_first_visit.rds"),
  file.path(cfg$paths$artifacts_dir, "site_habitat_all_sites.rds"),
  template_path
)
if (any(!file.exists(required))) {
  stop("Missing blind-review inputs: ", paste(
    required[!file.exists(required)], collapse = "; "
  ))
}

predictor_stack <- terra::rast(cfg$paths$predictor_stack)
rg_predictor_stack <- terra::rast(cfg$paths$predictor_stack_rg)
names(rg_predictor_stack)[13:17] <- c(
  "RG_bathy", "RG_slope", "RG_turbid",
  "RG_waves_Hs", "RG_waves_Tp"
)
clean_features <- names(predictor_stack)
rg_features <- names(rg_predictor_stack)

site_covariates <- readRDS(file.path(
  cfg$paths$artifacts_dir, "site_habitat_all_sites.rds"
)) |>
  dplyr::filter(.data$evidence_source == "cull_site_history") |>
  dplyr::select(
    "site_id", "clean_complete", "rg_complete",
    "RG_bathy", "RG_slope", "RG_turbid",
    "RG_waves_Hs", "RG_waves_Tp"
  )
annual <- readRDS(file.path(
  cfg$paths$artifacts_dir, "site_year_first_visit.rds"
)) |>
  dplyr::left_join(site_covariates, by = "site_id") |>
  dplyr::group_by(.data$site_id) |>
  dplyr::mutate(
    history_years = dplyr::n_distinct(.data$year),
    observation_weight = 1 / .data$history_years
  ) |>
  dplyr::ungroup()
clean_train <- dplyr::filter(annual, .data$clean_complete == 1L)
rg_train <- dplyr::filter(
  annual, .data$clean_complete == 1L, .data$rg_complete == 1L
)

message("Fitting full annual clean and RG models for blind visual review...")
annual_models <- list(
  clean = fit_xgb_binary(
    clean_train, clean_features, cfg$model,
    seed_offset = 61001L, weight_col = "observation_weight"
  ),
  rg = fit_xgb_binary(
    rg_train, rg_features, cfg$model,
    seed_offset = 61002L, weight_col = "observation_weight"
  )
)
saveRDS(
  annual_models,
  file.path(cfg$paths$models_dir, "annual_static_full_models.rds")
)

candidate_labels <- c(
  legacy_year = "Legacy year",
  legacy_agnostic = "Legacy agnostic",
  lifetime_rg_fallback = "Lifetime RG plus clean fallback",
  annual_rg_fallback = "Annual-response RG plus clean fallback",
  v2_static = "v2 operational static"
)
set.seed(20260827L)
randomised_keys <- sample(names(candidate_labels))
surface_codes <- paste("Surface", LETTERS[seq_along(randomised_keys)])
review_key <- tibble::tibble(
  surface = surface_codes,
  model_key = randomised_keys,
  model_label = unname(candidate_labels[randomised_keys])
)
utils::write.csv(
  review_key,
  file.path(cfg$paths$outputs_dir, "expert_blind_review_key.csv"),
  row.names = FALSE
)

percentile_raster <- function(x) {
  values <- terra::values(x, mat = FALSE)
  valid <- is.finite(values)
  ranked <- rep(NA_real_, length(values))
  if (sum(valid) == 1L) ranked[valid] <- 0.5
  if (sum(valid) > 1L) {
    ranked[valid] <- (rank(values[valid], ties.method = "average") - 1) /
      (sum(valid) - 1)
  }
  terra::setValues(x, ranked)
}
base64_png <- function(path) {
  paste0("data:image/png;base64,", base64enc::base64encode(path))
}

legacy_year <- terra::rast(cfg$paths$legacy_2026_raster)
legacy_agnostic <- terra::rast(cfg$paths$legacy_agnostic_raster)
regions <- v2_visual_regions()
map_data <- list()
for (region in regions) {
  message("Preparing blind review region: ", region$label)
  region_dir <- file.path(visual_dir, region$name)
  review_dir <- file.path(region_dir, "blind_review")
  dir.create(review_dir, recursive = TRUE, showWarnings = FALSE)

  extent_wgs84 <- terra::ext(
    region$west, region$east, region$south, region$north
  )
  extent_native <- terra::project(
    extent_wgs84, from = "EPSG:4326", to = terra::crs(predictor_stack)
  )
  clean_x <- terra::crop(predictor_stack, extent_native)
  rg_x <- terra::crop(rg_predictor_stack, extent_native)
  aggregation_factor <- cfg$visual$raster_aggregation_factor
  if (!is.null(aggregation_factor) && aggregation_factor > 1L) {
    clean_x <- terra::aggregate(
      clean_x, fact = aggregation_factor, fun = "mean", na.rm = TRUE
    )
    rg_x <- terra::aggregate(
      rg_x, fact = aggregation_factor, fun = "mean", na.rm = TRUE
    )
  }
  native_template <- clean_x[[1]]
  annual_clean <- predict_xgb_raster(
    clean_x, annual_models$clean,
    file.path(review_dir, "annual_clean.tif")
  )
  annual_rg <- predict_xgb_raster(
    rg_x, annual_models$rg,
    file.path(review_dir, "annual_rg.tif")
  )
  annual_clean_fallback <- terra::resample(
    annual_clean, annual_rg, method = "bilinear"
  )
  annual_rg_fallback <- terra::cover(annual_rg, annual_clean_fallback)

  lifetime_path <- file.path(
    region_dir, "site_habitat_rg_with_clean_fallback.tif"
  )
  static_path <- file.path(region_dir, "v2_static.tif")
  if (!file.exists(lifetime_path) || !file.exists(static_path)) {
    stop("Run 07_build_visual_comparison.R before blind review.")
  }
  native_layers <- list(
    legacy_year = align_raster_to_template(legacy_year, native_template),
    legacy_agnostic = align_raster_to_template(
      legacy_agnostic, native_template
    ),
    lifetime_rg_fallback = align_raster_to_template(
      terra::rast(lifetime_path), native_template
    ),
    annual_rg_fallback = annual_rg_fallback,
    v2_static = align_raster_to_template(
      terra::rast(static_path), native_template
    )
  )
  wgs_layers <- project_layers_to_wgs84(native_layers, native_template)
  images <- list()
  for (i in seq_len(nrow(review_key))) {
    surface <- review_key$surface[i]
    model_key <- review_key$model_key[i]
    image_path <- file.path(
      review_dir,
      paste0(tolower(gsub(" ", "_", surface)), ".png")
    )
    raster_to_gray_png(
      percentile_raster(wgs_layers[[model_key]]), image_path
    )
    images[[surface]] <- base64_png(image_path)
  }
  map_data[[region$name]] <- list(
    name = region$name,
    label = region$label,
    bounds = list(
      south = region$south, north = region$north,
      west = region$west, east = region$east
    ),
    images = images
  )
}

map_json <- jsonlite::toJSON(map_data, auto_unbox = TRUE, null = "null")
html <- paste(readLines(template_path, warn = FALSE), collapse = "\n")
html <- sub("__MAP_DATA__", map_json, html, fixed = TRUE)
html_path <- file.path(cfg$paths$outputs_dir, "expert_blind_review.html")
writeLines(html, html_path, useBytes = TRUE)

score_template <- tidyr::expand_grid(
  reviewer = "",
  region = names(map_data),
  surface = review_key$surface
) |>
  dplyr::mutate(
    known_rank = NA_integer_,
    fine_scale = NA_integer_,
    few_false_hotspots = NA_integer_,
    operational_use = NA_integer_,
    confidence = "",
    notes = ""
  )
utils::write.csv(
  score_template,
  file.path(cfg$paths$outputs_dir, "expert_review_scores_template.csv"),
  row.names = FALSE, na = ""
)

message("Blind expert review written to: ", html_path)
message("Keep the model key hidden until scoring is complete.")
