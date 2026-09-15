args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v2_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v2_dir, "R", "bootstrap.R"))
ctx <- initialise_v2(script_file)
cfg <- ctx$config

required <- c("terra", "xgboost", "viridisLite")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))

out_dir <- file.path(cfg$paths$outputs_dir, "multiscale_reef_maps")
fig_dir <- file.path(v2_dir, "reports", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

predictors <- terra::rast(cfg$paths$predictor_stack)
base_names <- names(predictors)
site_year <- readRDS(file.path(cfg$paths$artifacts_dir, "site_year_all.rds"))
site_scales <- readRDS(file.path(cfg$paths$artifacts_dir, "site_multiscale_30_60_120m.rds"))
site_folds <- readRDS(file.path(cfg$paths$artifacts_dir, "site_habitat_all_sites.rds")) |>
  dplyr::filter(.data$evidence_source == "cull_site_history") |>
  dplyr::select("site_id", "clean_complete")

training <- site_year |>
  dplyr::select(-dplyr::any_of(base_names)) |>
  dplyr::left_join(site_scales, by = c("site_id", "longitude", "latitude")) |>
  dplyr::left_join(site_folds, by = "site_id") |>
  dplyr::mutate(observation_weight = 1)

sets <- multiscale_feature_sets(base_names)
sets <- sets[c(
  "clean_120m_year", "multiscale_topography_year",
  "multiscale_topography_static"
)]
needed <- unique(setdiff(unlist(sets), "year"))
training <- training |>
  dplyr::filter(
    .data$clean_complete == 1L,
    stats::complete.cases(dplyr::pick(dplyr::all_of(needed)))
  )

legacy_workflow <- readRDS(file.path(
  cfg$paths$repo_root, "DataProcessing", "outputs",
  "fitted_model_clean_20260820_1041.rds"
))
legacy_fit <- workflows::extract_fit_parsnip(legacy_workflow$fit)
legacy_args <- legacy_fit$spec$args
legacy_cfg <- list(
  seed = cfg$model$seed,
  nrounds = as.integer(rlang::eval_tidy(legacy_args$trees)),
  eta = as.numeric(rlang::eval_tidy(legacy_args$learn_rate)),
  max_depth = as.integer(rlang::eval_tidy(legacy_args$tree_depth)),
  min_child_weight = as.numeric(rlang::eval_tidy(legacy_args$min_n)),
  gamma = as.numeric(rlang::eval_tidy(legacy_args$loss_reduction)),
  subsample = 1, colsample_bytree = 1, lambda = 1, alpha = 0,
  nthread = cfg$model$nthread
)

model_path <- file.path(cfg$paths$models_dir, "multiscale_full_visual_models.rds")
if (file.exists(model_path)) {
  models <- readRDS(model_path)
} else {
  models <- lapply(seq_along(sets), function(i) {
    fit_xgb_binary(training, sets[[i]], legacy_cfg, seed_offset = 24000L + i)
  })
  names(models) <- names(sets)
  saveRDS(models, model_path)
}

reef_names <- c(
  "Eyrie Reef (14-118)",
  "Batt Reef (16-029)",
  "(Big) Broadhurst Reef (No 1) (18-100a)",
  "Bait Reef (19-137)"
)

aligned_crop <- function(x, extent_native, factor = 4L) {
  corners <- rbind(
    c(terra::xmin(extent_native), terra::ymax(extent_native)),
    c(terra::xmax(extent_native), terra::ymin(extent_native))
  )
  rc <- terra::rowColFromCell(x, terra::cellFromXY(x, corners))
  r1 <- max(1L, ((min(rc[, 1]) - 1L) %/% factor) * factor + 1L)
  r2 <- min(terra::nrow(x), ceiling(max(rc[, 1]) / factor) * factor)
  c1 <- max(1L, ((min(rc[, 2]) - 1L) %/% factor) * factor + 1L)
  c2 <- min(terra::ncol(x), ceiling(max(rc[, 2]) / factor) * factor)
  corner_cells <- c(
    terra::cellFromRowCol(x, r1, c1),
    terra::cellFromRowCol(x, r2, c2)
  )
  terra::crop(x, terra::ext(x, cells = corner_cells), snap = "out")
}

make_multiscale_raster <- function(native) {
  at_scale <- function(factor, metres) {
    coarse <- if (factor == 1L) native else terra::aggregate(
      native, fact = factor, fun = "mean", na.rm = TRUE
    )
    if (factor > 1L) coarse <- terra::resample(coarse, native, method = "near")
    names(coarse) <- paste0(base_names, "_s", metres, "m")
    coarse
  }
  c(at_scale(1L, 30L), at_scale(2L, 60L), at_scale(4L, 120L))
}

rank_raster <- function(x) {
  v <- terra::values(x, mat = FALSE)
  ok <- is.finite(v)
  out <- rep(NA_real_, length(v))
  out[ok] <- rank(v[ok], ties.method = "average") / sum(ok)
  terra::setValues(x, out)
}

raw_summary <- function(x, reef, model, legacy_rank) {
  v <- terra::values(x, mat = FALSE)
  lr <- terra::values(legacy_rank, mat = FALSE)
  ok <- is.finite(v)
  q <- stats::quantile(v[ok], c(0.05, 0.5, 0.95), na.rm = TRUE)
  xr <- rank_raster(x)
  rr <- terra::values(xr, mat = FALSE)
  pair <- is.finite(rr) & is.finite(lr)
  tibble::tibble(
    reef_name = reef, model = model, valid_cells = sum(ok),
    q05 = unname(q[1]), median = unname(q[2]), q95 = unname(q[3]),
    q95_minus_q05 = unname(q[3] - q[1]),
    spearman_rank_vs_legacy = if (sum(pair) > 2L) stats::cor(rr[pair], lr[pair]) else NA_real_
  )
}

summaries <- list()
for (i in seq_along(reef_names)) {
  reef <- reef_names[i]
  sites <- site_year |>
    dplyr::filter(.data$reef_name == reef) |>
    dplyr::distinct(.data$site_id, .data$longitude, .data$latitude)
  if (!nrow(sites)) next

  bbox_wgs <- terra::ext(
    min(sites$longitude) - 0.025, max(sites$longitude) + 0.025,
    min(sites$latitude) - 0.025, max(sites$latitude) + 0.025
  )
  bbox_native <- terra::project(
    bbox_wgs, from = "EPSG:4326", to = terra::crs(predictors)
  )
  native <- aligned_crop(predictors, bbox_native, factor = 4L)
  multiscale <- make_multiscale_raster(native)

  reef_dir <- file.path(out_dir, sprintf("%02d_%s", i, gsub("[^A-Za-z0-9]+", "_", reef)))
  dir.create(reef_dir, recursive = TRUE, showWarnings = FALSE)
  layers <- list()
  layers$legacy_clean_year_2026 <- align_raster_to_template(
    terra::rast(cfg$paths$legacy_2026_raster), native[[1]]
  )
  for (name in names(models)) {
    features <- sets[[name]]
    x <- multiscale[[intersect(features, names(multiscale))]]
    if ("year" %in% features) {
      yr <- native[[1]]
      terra::values(yr) <- 2026
      names(yr) <- "year"
      x <- c(x, yr)
    }
    x <- x[[features]]
    layers[[name]] <- predict_xgb_raster(
      x, models[[name]], file.path(reef_dir, paste0(name, ".tif"))
    )
  }

  legacy_rank <- rank_raster(layers$legacy_clean_year_2026)
  for (name in names(layers)) {
    summaries[[length(summaries) + 1L]] <- raw_summary(
      layers[[name]], reef, name, legacy_rank
    )
  }

  site_points <- terra::project(
    terra::vect(as.data.frame(sites), geom = c("longitude", "latitude"), crs = "EPSG:4326"),
    terra::crs(native)
  )
  png_file <- file.path(fig_dir, sprintf("multiscale_rank_map_%02d.png", i))
  grDevices::png(png_file, width = 2000, height = 1650, res = 180)
  graphics::par(mfrow = c(2, 2), mar = c(2.2, 2.2, 3.5, 4.2), oma = c(0, 0, 3, 0))
  titles <- c(
    legacy_clean_year_2026 = "Legacy clean + year (2026)",
    clean_120m_year = "All variables averaged to 120 m + year",
    multiscale_topography_year = "30/60/120 m topography + year",
    multiscale_topography_static = "30/60/120 m static habitat"
  )
  for (name in names(titles)) {
    terra::plot(
      rank_raster(layers[[name]]), col = viridisLite::viridis(100),
      range = c(0, 1), main = titles[[name]], axes = FALSE,
      plg = list(title = "Within-map percentile", cex = 0.75)
    )
    terra::points(site_points, pch = 3, cex = 0.35, col = grDevices::adjustcolor("white", 0.65))
    graphics::box()
  }
  graphics::mtext(reef, outer = TRUE, cex = 1.25, font = 2)
  grDevices::dev.off()
}

summary_table <- dplyr::bind_rows(summaries)
utils::write.csv(
  summary_table,
  file.path(cfg$paths$outputs_dir, "multiscale_reef_map_summary.csv"),
  row.names = FALSE
)
message("Multiscale reef maps complete.")
print(summary_table)
