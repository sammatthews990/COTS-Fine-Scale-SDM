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
    "hmsc_joint_200m_production_model.rds"
  ),
  file.path(
    cfg$paths$outputs_dir, "hmsc_joint_200m_beta.csv"
  ),
  file.path(
    cfg$paths$artifacts_dir,
    "hmsc_joint_200m_cell_year.rds"
  ),
  file.path(cfg$paths$artifacts_dir, "site_year_all.rds")
)
if (any(!file.exists(required))) {
  stop("Missing Hmsc reef-map input(s).")
}
if (!requireNamespace("viridisLite", quietly = TRUE)) {
  stop("viridisLite is required for Hmsc reef maps.")
}

predictors <- terra::rast(required[[1]])
legacy_2026 <- terra::rast(required[[2]])
beta <- utils::read.csv(required[[4]]) |>
  dplyr::filter(.data$response == "cull_high_density")
beta_value <- stats::setNames(
  beta$posterior_mean, beta$covariate
)
hmsc_data <- readRDS(required[[5]])
site_year <- readRDS(required[[6]])
base_names <- names(predictors)

year_transform <- stats::coef(stats::lm(
  year_z ~ year, data = hmsc_data
))
fixed_protocol <- list(
  year = 2026,
  year_z = unname(
    year_transform[[1]] + year_transform[[2]] * 2026
  ),
  cull_bottom_time = stats::median(
    hmsc_data$cull_bottom_time, na.rm = TRUE
  ),
  cull_n_dives = stats::median(
    hmsc_data$cull_n_dives, na.rm = TRUE
  ),
  manta_n_tows = stats::median(
    hmsc_data$manta_n_tows, na.rm = TRUE
  ),
  manta_total_distance_m = stats::median(
    hmsc_data$manta_total_distance_m, na.rm = TRUE
  ),
  manta_hard_coral_fraction = stats::median(
    hmsc_data$manta_hard_coral_fraction, na.rm = TRUE
  )
)

reef_names <- c(
  "Eyrie Reef (14-118)",
  "Batt Reef (16-029)",
  "(Big) Broadhurst Reef (No 1) (18-100a)",
  "Bait Reef (19-137)"
)
reef_slugs <- c("eyrie", "batt", "broadhurst", "bait")
fig_dir <- file.path(v2_dir, "reports", "figures")
map_dir <- file.path(
  cfg$paths$outputs_dir, "hmsc_fixed_protocol_reef_maps"
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

predict_hmsc_fixed_protocol <- function(multiscale) {
  coral <- terra::app(
    multiscale[[
      c("AhyaD_s30m", "Aspat_s30m", "Aten_s30m")
    ]],
    fun = mean, na.rm = TRUE
  )
  names(coral) <- "coral_suitability_mean"
  spatial <- c(
    multiscale[[
      c(
        "BPI_s30m", "SLO_s60m", "VRM_s60m",
        "DEM_s30m", "SVF_s120m"
      )
    ]],
    coral
  )
  eta <- spatial[[1]] * 0 + beta_value[["(Intercept)"]]
  for (name in names(spatial)) {
    eta <- eta + spatial[[name]] * beta_value[[name]]
  }
  constants <- c(
    year_z = fixed_protocol$year_z,
    log_cull_bottom_time = log1p(
      fixed_protocol$cull_bottom_time
    ),
    log_cull_dives = log1p(fixed_protocol$cull_n_dives),
    log_manta_tows = log1p(fixed_protocol$manta_n_tows),
    log_manta_distance = log1p(
      fixed_protocol$manta_total_distance_m
    ),
    manta_hard_coral_fraction =
      fixed_protocol$manta_hard_coral_fraction
  )
  for (name in names(constants)) {
    eta <- eta + beta_value[[name]] * constants[[name]]
  }
  probability <- terra::app(
    eta, fun = function(x) stats::pnorm(x)
  )
  names(probability) <- "hmsc_fixed_protocol_probability"
  probability
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
  hmsc_probability <- predict_hmsc_fixed_protocol(multiscale)
  hmsc_file <- file.path(
    map_dir,
    sprintf("%02d_%s_hmsc_fixed_protocol.tif", i, reef_slugs[[i]])
  )
  hmsc_probability <- terra::writeRaster(
    hmsc_probability, hmsc_file, overwrite = TRUE,
    wopt = list(
      datatype = "FLT4S",
      gdal = c("COMPRESS=DEFLATE")
    )
  )
  legacy <- align_raster_to_template(
    legacy_2026, hmsc_probability
  )
  legacy_rank <- rank_raster(legacy)
  hmsc_rank <- rank_raster(hmsc_probability)
  rank_difference <- hmsc_rank - legacy_rank
  names(legacy_rank) <- "incumbent_percentile"
  names(hmsc_rank) <- "hmsc_percentile"
  names(rank_difference) <- "hmsc_minus_incumbent"

  site_points <- terra::project(
    terra::vect(
      as.data.frame(sites),
      geom = c("longitude", "latitude"),
      crs = "EPSG:4326"
    ),
    terra::crs(native)
  )
  png_name <- sprintf(
    "hmsc_vs_incumbent_%02d_%s.png",
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
    col = viridisLite::viridis(100), range = c(0, 1),
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
    hmsc_rank,
    col = viridisLite::viridis(100), range = c(0, 1),
    main = "Hmsc fixed effects at fixed protocol",
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
    main = "Hmsc minus current-best rank",
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

  legacy_value <- terra::values(legacy_rank, mat = FALSE)
  hmsc_value <- terra::values(hmsc_rank, mat = FALSE)
  raw_value <- terra::values(
    hmsc_probability, mat = FALSE
  )
  valid <- is.finite(legacy_value) & is.finite(hmsc_value)
  raw_valid <- raw_value[is.finite(raw_value)]
  summaries[[length(summaries) + 1L]] <- tibble::tibble(
    reef_name = reef,
    figure = file.path("figures", png_name),
    hmsc_surface = "posterior_mean_fixed_effects_fixed_protocol",
    incumbent_model = "legacy_clean_year_20260820_1041",
    valid_cells = sum(valid),
    raw_probability_q05 = unname(stats::quantile(
      raw_valid, 0.05
    )),
    raw_probability_median = stats::median(raw_valid),
    raw_probability_q95 = unname(stats::quantile(
      raw_valid, 0.95
    )),
    spearman_rank_correlation = stats::cor(
      legacy_value[valid], hmsc_value[valid],
      method = "spearman"
    ),
    mean_absolute_percentile_difference = mean(abs(
      hmsc_value[valid] - legacy_value[valid]
    )),
    review_status = "pending_visual_acceptance"
  )
}

manifest <- dplyr::bind_rows(summaries)
utils::write.csv(
  manifest,
  file.path(
    cfg$paths$outputs_dir,
    "hmsc_fixed_protocol_visual_map_manifest.csv"
  ),
  row.names = FALSE
)
utils::write.csv(
  tibble::as_tibble(fixed_protocol),
  file.path(
    cfg$paths$outputs_dir,
    "hmsc_fixed_protocol_map_settings.csv"
  ),
  row.names = FALSE
)

message("Hmsc fixed-protocol reef maps complete.")
print(manifest)
