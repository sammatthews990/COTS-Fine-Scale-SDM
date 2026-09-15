args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v2_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v2_dir, "R", "bootstrap.R"))
ctx <- initialise_v2(script_file)
cfg <- ctx$config

required <- c(
  cfg$paths$predictor_stack,
  cfg$paths$legacy_2026_raster,
  file.path(
    cfg$paths$models_dir,
    "two_part_fixed_effort_iteration_01.rds"
  ),
  file.path(cfg$paths$artifacts_dir, "site_year_all.rds")
)
if (any(!file.exists(required))) {
  stop("Missing candidate-versus-incumbent map input(s).")
}
if (!requireNamespace("viridisLite", quietly = TRUE)) {
  stop("viridisLite is required for reef comparison maps.")
}

bundle <- readRDS(required[[3]])
predictors <- terra::rast(required[[1]])
legacy_2026 <- terra::rast(required[[2]])
site_year <- readRDS(required[[4]])
base_names <- names(predictors)

reef_names <- c(
  "Eyrie Reef (14-118)",
  "Batt Reef (16-029)",
  "(Big) Broadhurst Reef (No 1) (18-100a)",
  "Bait Reef (19-137)"
)
reef_slugs <- c("eyrie", "batt", "broadhurst", "bait")
fig_dir <- file.path(v2_dir, "reports", "figures")
map_dir <- file.path(
  cfg$paths$outputs_dir,
  "two_part_iteration_01_candidate_maps"
)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(map_dir, recursive = TRUE, showWarnings = FALSE)

aligned_crop <- function(x, extent_native, factor = 4L) {
  corners <- rbind(
    c(terra::xmin(extent_native), terra::ymax(extent_native)),
    c(terra::xmax(extent_native), terra::ymin(extent_native))
  )
  rc <- terra::rowColFromCell(
    x, terra::cellFromXY(x, corners)
  )
  r1 <- max(
    1L, ((min(rc[, 1]) - 1L) %/% factor) * factor + 1L
  )
  r2 <- min(
    terra::nrow(x), ceiling(max(rc[, 1]) / factor) * factor
  )
  c1 <- max(
    1L, ((min(rc[, 2]) - 1L) %/% factor) * factor + 1L
  )
  c2 <- min(
    terra::ncol(x), ceiling(max(rc[, 2]) / factor) * factor
  )
  cells <- c(
    terra::cellFromRowCol(x, r1, c1),
    terra::cellFromRowCol(x, r2, c2)
  )
  terra::crop(x, terra::ext(x, cells = cells), snap = "out")
}

make_multiscale_raster <- function(native) {
  at_scale <- function(factor, metres) {
    coarse <- if (factor == 1L) {
      native
    } else {
      terra::aggregate(
        native, fact = factor, fun = "mean", na.rm = TRUE
      )
    }
    if (factor > 1L) {
      coarse <- terra::resample(
        coarse, native, method = "near"
      )
    }
    names(coarse) <- paste0(
      base_names, "_s", metres, "m"
    )
    coarse
  }
  c(
    at_scale(1L, 30L),
    at_scale(2L, 60L),
    at_scale(4L, 120L)
  )
}

rank_raster <- function(x) {
  value <- terra::values(x, mat = FALSE)
  valid <- is.finite(value)
  ranked <- rep(NA_real_, length(value))
  ranked[valid] <- (
    rank(value[valid], ties.method = "average") - 1
  ) / max(1, sum(valid) - 1)
  terra::setValues(x, ranked)
}

summaries <- list()
for (i in seq_along(reef_names)) {
  reef <- reef_names[[i]]
  sites <- site_year |>
    dplyr::filter(.data$reef_name == reef) |>
    dplyr::distinct(
      .data$site_id, .data$longitude, .data$latitude
    )
  if (nrow(sites) == 0L) next

  bbox_wgs <- terra::ext(
    min(sites$longitude) - 0.025,
    max(sites$longitude) + 0.025,
    min(sites$latitude) - 0.025,
    max(sites$latitude) + 0.025
  )
  bbox_native <- terra::project(
    bbox_wgs, from = "EPSG:4326", to = terra::crs(predictors)
  )
  native <- aligned_crop(
    predictors, bbox_native, factor = 4L
  )
  multiscale <- make_multiscale_raster(native)
  candidate_inputs <- multiscale[[
    bundle$static_features
  ]]
  candidate_file <- file.path(
    map_dir,
    sprintf("%02d_%s_candidate_static.tif", i, reef_slugs[[i]])
  )
  candidate <- predict_xgb_raster(
    candidate_inputs,
    bundle$static_model,
    candidate_file
  )
  legacy <- align_raster_to_template(
    legacy_2026, candidate
  )
  legacy_rank <- rank_raster(legacy)
  candidate_rank <- rank_raster(candidate)
  rank_difference <- candidate_rank - legacy_rank
  names(legacy_rank) <- "incumbent_percentile"
  names(candidate_rank) <- "candidate_percentile"
  names(rank_difference) <- "candidate_minus_incumbent"

  site_points <- terra::project(
    terra::vect(
      as.data.frame(sites),
      geom = c("longitude", "latitude"),
      crs = "EPSG:4326"
    ),
    terra::crs(native)
  )
  png_name <- sprintf(
    "candidate_vs_incumbent_%02d_%s.png",
    i, reef_slugs[[i]]
  )
  png_file <- file.path(fig_dir, png_name)
  grDevices::png(
    png_file, width = 2550, height = 950, res = 180
  )
  graphics::par(
    mfrow = c(1, 3), mar = c(2.2, 2.2, 4.0, 4.2),
    oma = c(0, 0, 3, 0)
  )
  terra::plot(
    legacy_rank,
    col = viridisLite::viridis(100),
    range = c(0, 1),
    main = "Current best: clean + year (2026)",
    axes = FALSE,
    plg = list(title = "Within-reef percentile", cex = 0.72)
  )
  terra::points(
    site_points, pch = 3, cex = 0.30,
    col = grDevices::adjustcolor("white", 0.65)
  )
  graphics::box()
  terra::plot(
    candidate_rank,
    col = viridisLite::viridis(100),
    range = c(0, 1),
    main = "Candidate: multiscale static habitat",
    axes = FALSE,
    plg = list(title = "Within-reef percentile", cex = 0.72)
  )
  terra::points(
    site_points, pch = 3, cex = 0.30,
    col = grDevices::adjustcolor("white", 0.65)
  )
  graphics::box()
  terra::plot(
    rank_difference,
    col = grDevices::hcl.colors(
      101, palette = "Blue-Red 3", rev = TRUE
    ),
    range = c(-1, 1),
    main = "Candidate minus current-best rank",
    axes = FALSE,
    plg = list(title = "Percentile difference", cex = 0.72)
  )
  terra::points(
    site_points, pch = 3, cex = 0.30,
    col = grDevices::adjustcolor("black", 0.55)
  )
  graphics::box()
  graphics::mtext(
    reef, outer = TRUE, cex = 1.28, font = 2
  )
  grDevices::dev.off()

  legacy_value <- terra::values(
    legacy_rank, mat = FALSE
  )
  candidate_value <- terra::values(
    candidate_rank, mat = FALSE
  )
  valid <- is.finite(legacy_value) &
    is.finite(candidate_value)
  summaries[[length(summaries) + 1L]] <- tibble::tibble(
    reef_name = reef,
    figure = file.path("figures", png_name),
    candidate_model = bundle$iteration,
    incumbent_model = "legacy_clean_year_20260820_1041",
    valid_cells = sum(valid),
    spearman_rank_correlation = stats::cor(
      legacy_value[valid], candidate_value[valid],
      method = "spearman"
    ),
    mean_absolute_percentile_difference = mean(abs(
      candidate_value[valid] - legacy_value[valid]
    )),
    review_status = "pending_visual_acceptance"
  )
}

manifest <- dplyr::bind_rows(summaries)
utils::write.csv(
  manifest,
  file.path(
    cfg$paths$outputs_dir,
    "two_part_iteration_01_visual_map_manifest.csv"
  ),
  row.names = FALSE
)
review_template <- dplyr::transmute(
  manifest,
  reef_name = .data$reef_name,
  candidate_model = .data$candidate_model,
  reviewer = "",
  spatial_pattern_score_1_to_5 = NA_integer_,
  known_high_density_alignment_1_to_5 = NA_integer_,
  misleading_hotspots_1_to_5 = NA_integer_,
  operational_usefulness_1_to_5 = NA_integer_,
  accept_candidate = NA,
  notes = ""
)
utils::write.csv(
  review_template,
  file.path(
    v2_dir, "templates",
    "two_part_iteration_01_visual_review.csv"
  ),
  row.names = FALSE, na = ""
)

message("Candidate-versus-incumbent reef maps complete.")
print(manifest)
