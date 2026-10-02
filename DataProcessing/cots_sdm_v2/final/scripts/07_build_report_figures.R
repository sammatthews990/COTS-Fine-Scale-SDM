# Report figures for the three delivered layers. Run after 01-06.
#
#   fig01  GBR-wide overview of each layer
#   fig02  Four expert-review reefs: within-reef percentiles + Hmsc uncertainty
#   fig03  Independent Manta check of within-reef ranking
#   fig04  Hmsc target-response habitat effects
#   fig05  Hmsc joint coefficients across all evidence streams
#   fig06  Legacy XGBoost variable importance, clean vs RG
#   fig07  Within-reef agreement between layers
#   fig08  Coverage by prediction source

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "..", "R", "common.R"))
cfg <- initialise_final()
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})
layers_dir <- cfg$dirs$layers
tables_dir <- cfg$dirs$tables
layer_ids <- names(layer_labels)
read_table <- function(name) utils::read.csv(file.path(tables_dir, name), stringsAsFactors = FALSE)

support <- c(
  terra::rast(file.path(cfg$dirs$cache, "legacy_clean_year2026_120m.tif")),
  terra::rast(file.path(cfg$dirs$cache, "legacy_rg_year2026_120m.tif")),
  terra::crop(
    terra::extend(
      terra::rast(file.path(layers_dir, "cots_hmsc_joint_120m.tif")),
      terra::rast(file.path(cfg$dirs$cache, "legacy_clean_year2026_120m.tif"))
    ),
    terra::rast(file.path(cfg$dirs$cache, "legacy_clean_year2026_120m.tif"))
  )
)
names(support) <- layer_ids
features <- sf::st_read(checked_input(cfg, "reef_features"), quiet = TRUE) |>
  sf::st_transform(terra::crs(support))
land <- features |>
  dplyr::filter(.data$FEAT_NAME %in% c("Mainland", "Island")) |>
  sf::st_simplify(dTolerance = 300)

value_titles <- c(
  legacy_clean_year2026 = "P(CPUE ≥ 0.02), 2026",
  legacy_rg_year2026 = "P(CPUE ≥ 0.02), 2026",
  hmsc_joint = "P(problematic density | active reef-year)"
)

# --- fig01: GBR overview ------------------------------------------------------
overview <- terra::aggregate(support, fact = 24, fun = "mean", na.rm = TRUE)
bbox <- sf::st_bbox(terra::ext(terra::trim(overview[[1]])) |> terra::vect(crs = terra::crs(support)) |> sf::st_as_sf())
overview_panel <- function(id) {
  d <- raster_df(overview[[id]])
  limits <- stats::quantile(d$value, c(0.02, 0.98))
  ggplot() +
    geom_sf(data = land, fill = ink$land, colour = NA) +
    geom_raster(data = d, aes(x, y, fill = value)) +
    scale_fill_gradientn(
      colours = ramp_blue, limits = limits, oob = scales::squish,
      name = value_titles[[id]], labels = scales::label_number(accuracy = 0.01)
    ) +
    coord_sf(xlim = bbox[c(1, 3)], ylim = bbox[c(2, 4)], expand = FALSE, datum = NA) +
    labs(title = layer_labels[[id]]) +
    theme_map() +
    theme(
      legend.position = "bottom", legend.key.width = grid::unit(0.9, "cm"),
      legend.key.height = grid::unit(0.18, "cm"), legend.title.position = "top"
    )
}
fig01 <- wrap_plots(lapply(layer_ids, overview_panel), nrow = 1) +
  plot_annotation(
    title = "Three delivered COTS layers across the GBR",
    subtitle = "Cell means at ~2.9 km for display; colour stretch is each layer's 2nd-98th percentile.",
    caption = "Gaps in RG + year are reefs without ReefGuide exposure layers (see fig08 for the clean fallback).",
    theme = theme_report()
  )
save_figure(fig01, cfg, "fig01_gbr_overview.png", width = 7.5, height = 5.2)

# --- fig02: four expert-review reefs ----------------------------------------
native <- list(
  legacy_clean_year2026 = terra::rast(file.path(layers_dir, "cots_legacy_clean_year2026_30m.tif")),
  legacy_rg_year2026 = terra::rast(file.path(layers_dir, "cots_legacy_rg_year2026_30m.tif")),
  hmsc_joint = terra::rast(file.path(layers_dir, "cots_hmsc_joint_120m.tif"))
)
# The history-updated layer (script 10) is shown when it exists.
combined_file <- file.path(layers_dir, "cots_hmsc_history_combined_120m.tif")
fig02_layers <- c(
  legacy_clean_year2026 = "Clean + year", legacy_rg_year2026 = "RG + year", hmsc_joint = "HMSC habitat"
)
if (file.exists(combined_file)) {
  native$hmsc_history_combined <- terra::rast(combined_file)
  fig02_layers <- c(fig02_layers, hmsc_history_combined = "HMSC + history")
}
uncertainty <- terra::rast(file.path(layers_dir, "cots_hmsc_joint_uncertainty_120m.tif"))
interval_width <- uncertainty[["posterior_q95"]] - uncertainty[["posterior_q05"]]
reef_ids <- cfg$figures$four_reefs

percentile_in_polygon <- function(r, polygon) {
  v <- terra::vect(polygon)
  x <- terra::mask(terra::crop(r, v, snap = "out"), v, touches = terra::res(r)[1] > 60)
  values <- terra::values(x, mat = FALSE)
  ok <- is.finite(values)
  values[ok] <- (rank(values[ok]) - 1) / max(sum(ok) - 1, 1)
  terra::setValues(x, values)
}
width_max <- 0
reef_panels <- lapply(names(reef_ids), function(slug) {
  polygon <- features[features$LABEL_ID == reef_ids[[slug]], ]
  frame <- sf::st_bbox(sf::st_buffer(polygon, 300))
  cells <- lapply(names(fig02_layers), function(id) raster_df(percentile_in_polygon(native[[id]], polygon)))
  names(cells) <- names(fig02_layers)
  v <- terra::vect(polygon)
  width_df <- raster_df(terra::mask(terra::crop(interval_width, v, snap = "out"), v, touches = TRUE))
  width_max <<- max(width_max, width_df$value)
  list(
    slug = slug, name = polygon$GBR_NAME[1], label = reef_ids[[slug]], polygon = polygon,
    frame = frame, cells = cells, width = width_df
  )
})
map_panel <- function(d, polygon, frame, fill_scale, title = NULL, ylab = NULL) {
  ggplot() +
    geom_raster(data = d, aes(x, y, fill = value)) +
    geom_sf(data = polygon, fill = NA, colour = ink$secondary, linewidth = 0.2) +
    fill_scale +
    coord_sf(xlim = frame[c(1, 3)], ylim = frame[c(2, 4)], expand = FALSE, datum = NA) +
    labs(title = title, y = ylab) +
    theme_map(8) +
    theme(axis.title.y = element_text(angle = 90, colour = ink$primary, face = "bold", size = 8))
}
percentile_scale <- scale_fill_gradientn(
  colours = ramp_blue, limits = c(0, 1), name = "Within-reef percentile",
  breaks = c(0, 0.5, 1), labels = c("lowest", "median", "highest")
)
width_scale <- scale_fill_gradientn(
  colours = ramp_orange, limits = c(0, width_max),
  name = "HMSC 90% interval width", labels = scales::label_number(accuracy = 0.01)
)
# Top-rank emphasis for fig02: everything below the 70th within-reef
# percentile shares the darkest viridis tone, so the full ramp separates the
# top 30% of each reef.
top_rank_floor <- 0.7
top_rank_scale <- scale_fill_viridis_c(
  option = "viridis", limits = c(top_rank_floor, 1), oob = scales::squish,
  breaks = c(0.7, 0.8, 0.9, 1), labels = c("\u2264 70th", "80th", "90th", "100th"),
  name = "Within-reef percentile (top 30% emphasised)"
)
column_titles <- c(unname(fig02_layers), "HMSC uncertainty")
n_cols <- length(column_titles)
fig02_panels <- unlist(lapply(seq_along(reef_panels), function(i) {
  p <- reef_panels[[i]]
  titles <- if (i == 1) column_titles else rep(list(NULL), n_cols)
  c(
    lapply(seq_along(fig02_layers), function(j) {
      map_panel(
        p$cells[[names(fig02_layers)[j]]], p$polygon, p$frame, top_rank_scale, titles[[j]],
        if (j == 1) paste0(sub(" \\(.*", "", p$name), "\n", p$label) else NULL
      )
    }),
    list(map_panel(p$width, p$polygon, p$frame, width_scale, titles[[n_cols]]))
  )
}), recursive = FALSE)
fig02 <- wrap_plots(fig02_panels, ncol = n_cols, guides = "collect") +
  plot_annotation(
    title = "Four expert-review reefs: where each layer ranks cells highest",
    subtitle = paste(
      "Each layer is ranked within the reef polygon at its delivered resolution (legacy 30 m, HMSC 120 m);",
      "colour resolves only the top 30% of each reef.",
      sep = "\n"
    ),
    caption = paste(
      "HMSC + history is the 2027 forecast: habitat updated by prior Cull records (equals habitat where no Cull history lies nearby).",
      "Uncertainty is the posterior 5th-95th percentile width of the HMSC habitat probability.",
      sep = "\n"
    ),
    theme = theme_report()
  ) &
  theme(
    legend.position = "bottom", legend.key.width = grid::unit(1.1, "cm"),
    legend.key.height = grid::unit(0.18, "cm"), legend.title.position = "top",
    legend.spacing.x = grid::unit(1.2, "cm")
  )
save_figure(fig02, cfg, "fig02_four_reef_comparison.png", width = 7.5, height = 7.4)

# --- fig02b: same panels over satellite imagery -------------------------------
# Esri World Imagery tiles, cached under cache/tiles. Cells below the top-rank
# floor are drawn nearly transparent so the reef shows through; the top 30%
# stays more opaque.
imagery_provider <- "Esri.WorldImagery"
tile_dir <- file.path(cfg$dirs$cache, "tiles")
dir.create(tile_dir, recursive = TRUE, showWarnings = FALSE)
reef_imagery <- lapply(reef_panels, function(p) {
  frame_sf <- sf::st_as_sfc(p$frame)
  maptiles::get_tiles(frame_sf,
    provider = imagery_provider, zoom = 14, crop = TRUE,
    project = TRUE, cachedir = tile_dir
  )
})
imagery_panel <- function(d, image, polygon, frame, fill_scale, title = NULL, ylab = NULL, rank_alpha = TRUE) {
  d$alpha <- if (rank_alpha) ifelse(d$value >= top_rank_floor, 0.8, 0.4) else 0.6
  ggplot() +
    tidyterra::geom_spatraster_rgb(data = image, maxcell = 2e6) +
    geom_raster(data = d, aes(x, y, fill = value, alpha = alpha)) +
    geom_sf(data = polygon, fill = NA, colour = "white", linewidth = 0.25) +
    fill_scale +
    scale_alpha_identity() +
    coord_sf(xlim = frame[c(1, 3)], ylim = frame[c(2, 4)], expand = FALSE, datum = NA) +
    labs(title = title, y = ylab) +
    theme_map(8) +
    theme(axis.title.y = element_text(angle = 90, colour = ink$primary, face = "bold", size = 8))
}
fig02b_panels <- unlist(lapply(seq_along(reef_panels), function(i) {
  p <- reef_panels[[i]]
  image <- reef_imagery[[i]]
  titles <- if (i == 1) column_titles else rep(list(NULL), n_cols)
  c(
    lapply(seq_along(fig02_layers), function(j) {
      imagery_panel(
        p$cells[[names(fig02_layers)[j]]], image, p$polygon, p$frame, top_rank_scale, titles[[j]],
        if (j == 1) paste0(sub(" \\(.*", "", p$name), "\n", p$label) else NULL
      )
    }),
    list(imagery_panel(p$width, image, p$polygon, p$frame, width_scale, titles[[n_cols]], rank_alpha = FALSE))
  )
}), recursive = FALSE)
fig02b <- wrap_plots(fig02b_panels, ncol = n_cols, guides = "collect") +
  plot_annotation(
    title = "Four expert-review reefs over satellite imagery",
    subtitle = paste(
      "Same rankings as fig02. Cells below each reef's 70th percentile are near-transparent;",
      "the top 30% are drawn more opaque so they stand out against the reef.",
      sep = "\n"
    ),
    caption = paste(
      "HMSC + history is the 2027 forecast: habitat updated by prior Cull records. Uncertainty: posterior 5th-95th percentile width.",
      paste("Imagery:", maptiles::get_credit(imagery_provider)),
      sep = "\n"
    ),
    theme = theme_report()
  ) &
  theme(
    legend.position = "bottom", legend.key.width = grid::unit(1.1, "cm"),
    legend.key.height = grid::unit(0.18, "cm"), legend.title.position = "top",
    legend.spacing.x = grid::unit(1.2, "cm")
  )
save_figure(fig02b, cfg, "fig02b_four_reef_comparison_imagery.png", width = 7.5, height = 7.4)

# --- fig03: independent Manta check -----------------------------------------
manta <- read_table("manta_independent_check.csv") |>
  dplyr::mutate(layer = factor(.data$layer, levels = rev(unname(layer_labels))))
metric_panel <- function(estimate, low, high, reference, xlab, accuracy) {
  ggplot(manta, aes(y = layer, colour = layer_id)) +
    geom_vline(xintercept = reference, colour = ink$muted, linetype = "22", linewidth = 0.4) +
    geom_segment(aes(x = .data[[low]], xend = .data[[high]], yend = layer), linewidth = 2, lineend = "round", alpha = 0.35) +
    geom_point(aes(x = .data[[estimate]]), size = 3) +
    geom_text(aes(x = .data[[high]], label = scales::number(.data[[estimate]], accuracy = accuracy)),
      hjust = -0.3, colour = ink$primary, size = 3, family = report_font
    ) +
    scale_colour_manual(values = layer_colours, guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.25))) +
    labs(x = xlab, y = NULL) +
    theme_report() +
    theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(colour = ink$primary, size = 9))
}
fig03 <- (metric_panel(
  "within_reef_year_auc", "auc_low", "auc_high", 0.5,
  "Within-reef-year AUC (0.5 = chance)", 0.001
) |
  (metric_panel(
    "top20_lift", "lift_low", "lift_high", 1,
    "Top-20% lift within reef-year (1 = chance)", 0.01
  ) +
    theme(axis.text.y = element_blank()))) +
  plot_annotation(
    title = "Independent check: do Manta tows with strong COTS evidence score higher?",
    subtitle = sprintf(
      "%s tows (%s strong evidence) in %s active reef-years across %s reefs; every layer scores the same tows.",
      format(manta$tows[1], big.mark = ","), format(manta$positives[1], big.mark = ","),
      format(manta$reef_years[1], big.mark = ","), manta$reefs[1]
    ),
    caption = paste(
      "Strong evidence: direct COTS, >= 2 counted scars or legacy 'common' scars; compared with strict-negative tows.",
      "Bars: 95% bootstrap interval over reef-years. Manta is not a response in any model (HMSC uses it only to flag active reef-years).",
      sep = "\n"
    ),
    theme = theme_report()
  )
save_figure(fig03, cfg, "fig03_manta_independent_check.png", width = 7.5, height = 2.8)

# --- fig04/05: Hmsc coefficients ---------------------------------------------
response_labels <- c(
  cull_high_density = "Cull: CPUE ≥ 0.02 (mapped target)", cull_direct = "Cull: COTS present",
  salad_high_density = "SALAD: problematic density", salad_direct = "SALAD: COTS present",
  salad_scar = "SALAD: feeding scars", salad_log_density = "SALAD: log density",
  rhis_adult = "RHIS: adult COTS", rhis_juvenile = "RHIS: juvenile COTS",
  rhis_branching_scar = "RHIS: branching-coral scars", rhis_table_scar = "RHIS: table-coral scars",
  rhis_branching_live_asin = "RHIS: live branching coral", rhis_table_live_asin = "RHIS: live table coral"
)
covariate_labels <- c(
  AhyaD = "Acropora hyacinthus suitability", Aspat = "A. spathulata suitability",
  Aten = "A. tenuis suitability", turbid = "Turbidity", waves_Hs = "Wave height (Hs)",
  waves_Tp = "Wave period (Tp)", dem_pc1 = "Depth (multiscale)", slope_pc1 = "Slope (multiscale)",
  eastness_pc1 = "Eastness (multiscale)", northness_pc1 = "Northness (multiscale)",
  hcu_pc1 = "Plan curvature (multiscale)", vcu_pc1 = "Profile curvature (multiscale)",
  bpi_pc1 = "Bathymetric position (multiscale)", sd_pc1 = "Local relief (multiscale)"
)
hmsc_model <- readRDS(file.path(cfg$dirs$models, "hmsc_enhanced_joint_production.rds"))
hmsc_artifact <- readRDS(checked_input(cfg, "hmsc_enhanced_joint"))
pooled <- Hmsc::poolMcmcChains(hmsc_model$postList)
beta <- simplify2array(lapply(pooled, `[[`, "Beta"))
dimnames(beta)[1:2] <- list(hmsc_model$covNames, hmsc_artifact$responses)
scale_factor <- vapply(hmsc_model$covNames, function(cn) {
  if (cn %in% names(hmsc_artifact$XData) && is.numeric(hmsc_artifact$XData[[cn]])) stats::sd(hmsc_artifact$XData[[cn]]) else 1
}, numeric(1))
coefficients <- tidyr::expand_grid(covariate = hmsc_model$covNames, response = hmsc_artifact$responses) |>
  dplyr::filter(.data$covariate != "(Intercept)") |>
  dplyr::rowwise() |>
  dplyr::mutate(
    draws = list(beta[.data$covariate, .data$response, ] * scale_factor[[.data$covariate]]),
    mean = mean(.data$draws), q05 = stats::quantile(.data$draws, 0.05), q95 = stats::quantile(.data$draws, 0.95),
    support = mean(.data$draws > 0)
  ) |>
  dplyr::ungroup() |>
  dplyr::select(-"draws") |>
  dplyr::mutate(
    standardised = .data$covariate %in% names(scale_factor)[scale_factor != 1],
    covariate_label = dplyr::if_else(
      grepl("^geomorphic_class", .data$covariate),
      paste0("Geomorphic: ", sub("^geomorphic_class", "", .data$covariate)),
      unname(covariate_labels[.data$covariate])
    ),
    credible = .data$q05 > 0 | .data$q95 < 0
  )
write_table(coefficients, cfg, "hmsc_enhanced_joint_coefficients.csv")

target <- coefficients |>
  dplyr::filter(.data$response == cfg$hmsc$target_response) |>
  dplyr::mutate(covariate_label = stats::reorder(.data$covariate_label, .data$mean))
fig04 <- ggplot(target, aes(y = covariate_label)) +
  geom_vline(xintercept = 0, colour = ink$axis, linewidth = 0.4) +
  geom_segment(aes(x = q05, xend = q95, yend = covariate_label),
    colour = layer_colours[["hmsc_joint"]],
    linewidth = 1.6, lineend = "round", alpha = 0.4
  ) +
  geom_point(aes(x = mean, shape = credible), colour = layer_colours[["hmsc_joint"]], size = 2.4, fill = ink$surface) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), labels = c(`TRUE` = "90% CI excludes 0", `FALSE` = "90% CI includes 0"), name = NULL) +
  labs(
    title = "HMSC joint: habitat effects on the mapped target",
    subtitle = "Posterior mean and 90% credible interval, enhanced component (probit scale).",
    x = "Effect on probit P(Cull CPUE ≥ 0.02) per SD (continuous) or vs back reef slope (geomorphic)",
    y = NULL,
    caption = "Multiscale terrain enters as the first principal component of each family across 30-360 m Gaussian scales; PC signs are arbitrary."
  ) +
  theme_report() +
  theme(
    legend.position = "top", legend.justification = "left", panel.grid.major.y = element_blank(),
    axis.text.y = element_text(colour = ink$primary)
  )
save_figure(fig04, cfg, "fig04_hmsc_target_effects.png", width = 7.5, height = 6)

covariate_order <- target$covariate_label[order(target$mean)]
heat <- coefficients |>
  dplyr::mutate(
    response_label = factor(response_labels[.data$response], levels = rev(response_labels)),
    covariate_label = factor(.data$covariate_label, levels = levels(covariate_order))
  )
limit <- max(abs(heat$mean))
fig05 <- ggplot(heat, aes(x = covariate_label, y = response_label)) +
  geom_tile(aes(fill = mean), colour = ink$surface, linewidth = 0.6) +
  geom_point(data = dplyr::filter(heat, credible), size = 0.9, colour = ink$primary) +
  scale_fill_gradientn(colours = diverging_blue_red, limits = c(-limit, limit), name = "Posterior mean effect") +
  labs(
    title = "HMSC joint: how each evidence stream responds to the same habitat variables",
    subtitle = "Rows are responses fitted jointly; columns ordered by effect on the mapped Cull target. Dots: 90% CI excludes 0.",
    x = NULL, y = NULL,
    caption = "Continuous covariates per SD; geomorphic classes vs back reef slope. Probit scale for binary, identity for continuous responses."
  ) +
  theme_report(8) +
  theme(
    axis.text.x = element_text(angle = 50, hjust = 1, colour = ink$primary),
    axis.text.y = element_text(colour = ink$primary), panel.grid.major = element_blank(),
    legend.position = "bottom", legend.key.width = grid::unit(1.2, "cm"),
    legend.key.height = grid::unit(0.18, "cm"), legend.title.position = "top"
  )
save_figure(fig05, cfg, "fig05_hmsc_joint_coefficients.png", width = 7.5, height = 6.2)

# --- fig06: legacy variable importance --------------------------------------
importance <- dplyr::bind_rows(lapply(c(legacy_clean_year2026 = "legacy_clean_model", legacy_rg_year2026 = "legacy_rg_model"), function(input) {
  legacy <- readRDS(checked_input(cfg, input))
  booster <- workflows::extract_fit_parsnip(legacy$fit)$fit
  imp <- xgboost::xgb.importance(model = booster)
  tibble::tibble(feature = imp$Feature, gain_share = imp$Gain / sum(imp$Gain))
}), .id = "layer_id")
write_table(importance, cfg, "legacy_variable_importance.csv")
predictor_labels <- c(
  Year = "Survey year", BPI = "Bathymetric position", EAST = "Eastness", HCU = "Plan curvature",
  NORTH = "Northness", SLO = "Slope", SVF = "Sky-view factor", VCU = "Profile curvature",
  VRM = "Ruggedness (VRM)", DEM = "Depth", AhyaD = "Acropora hyacinthus suitability",
  Aspat = "A. spathulata suitability", Aten = "A. tenuis suitability",
  RG_bathy = "RG bathymetry", RG_slope = "RG slope", RG_turbid = "RG turbidity",
  RG_waves_Hs = "RG wave height (Hs)", RG_waves_Tp = "RG wave period (Tp)"
)
importance_plot <- importance |>
  # A predictor absent from one model gets a zero-length bar so both keep their slot.
  tidyr::complete(layer_id, feature, fill = list(gain_share = 0)) |>
  dplyr::mutate(feature_label = dplyr::coalesce(predictor_labels[.data$feature], .data$feature)) |>
  dplyr::group_by(.data$feature_label) |>
  dplyr::mutate(order = max(.data$gain_share)) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    feature_label = stats::reorder(.data$feature_label, .data$order),
    layer = factor(layer_labels[.data$layer_id], levels = layer_labels[c("legacy_rg_year2026", "legacy_clean_year2026")])
  )
fig06 <- ggplot(importance_plot, aes(x = gain_share, y = feature_label, fill = layer_id)) +
  geom_col(aes(group = layer), position = position_dodge(width = 0.82), width = 0.74, colour = ink$surface, linewidth = 0.4) +
  scale_fill_manual(
    values = layer_colours[c("legacy_clean_year2026", "legacy_rg_year2026")],
    labels = layer_labels[c("legacy_clean_year2026", "legacy_rg_year2026")], name = NULL,
    breaks = c("legacy_clean_year2026", "legacy_rg_year2026")
  ) +
  scale_x_continuous(labels = scales::label_percent(accuracy = 1), expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Legacy XGBoost models: what drives each layer",
    subtitle = "Share of total gain. Survey year is the largest single driver of both; RG wave and turbidity layers take much of its share.",
    x = "Share of XGBoost gain", y = NULL,
    caption = "Gain is unsigned and model-specific; it shows reliance, not direction of effect."
  ) +
  theme_report() +
  theme(
    legend.position = "top", legend.justification = "left", panel.grid.major.y = element_blank(),
    axis.text.y = element_text(colour = ink$primary), legend.key.size = grid::unit(0.3, "cm")
  )
save_figure(fig06, cfg, "fig06_legacy_variable_importance.png", width = 7.5, height = 5)

# --- fig07: within-reef agreement -------------------------------------------
agreement <- read_table("layer_agreement_within_reef.csv") |>
  dplyr::filter(is.finite(.data$spearman)) |>
  dplyr::mutate(pair = factor(
    paste(layer_labels[.data$layer_a], "vs", layer_labels[.data$layer_b]),
    levels = rev(c("Clean + year vs RG + year", "Clean + year vs HMSC joint", "RG + year vs HMSC joint"))
  ))
agreement_summary <- agreement |>
  dplyr::group_by(.data$pair) |>
  dplyr::summarise(
    reefs = dplyr::n(), median = stats::median(.data$spearman),
    q25 = stats::quantile(.data$spearman, 0.25), q75 = stats::quantile(.data$spearman, 0.75), .groups = "drop"
  )
write_table(agreement_summary, cfg, "layer_agreement_summary.csv")
fig07 <- ggplot(agreement, aes(x = spearman, y = pair)) +
  geom_vline(xintercept = 0, colour = ink$axis, linewidth = 0.4) +
  geom_violin(fill = ink$grid, colour = NA, scale = "width", width = 0.7) +
  geom_segment(data = agreement_summary, aes(x = q25, xend = q75, yend = pair), colour = ink$secondary, linewidth = 1.6, lineend = "round") +
  geom_point(data = agreement_summary, aes(x = median), colour = ink$primary, size = 2.4) +
  geom_text(
    data = agreement_summary, aes(x = median, label = sprintf("median %.2f\n%s reefs", median, format(reefs, big.mark = ","))),
    vjust = -0.9, size = 2.7, colour = ink$primary, family = report_font, lineheight = 0.9
  ) +
  scale_x_continuous(limits = c(-1, 1)) +
  labs(
    title = "How similarly do the layers rank cells within the same reef?",
    subtitle = "Spearman correlation between layers across 120 m cells inside each reef polygon (reefs with ≥ 30 shared cells).",
    x = "Within-reef Spearman correlation", y = NULL,
    caption = "Legacy layers are averaged to the HMSC 120 m support for this comparison. Bar: interquartile range."
  ) +
  theme_report() +
  theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(colour = ink$primary, size = 9))
save_figure(fig07, cfg, "fig07_layer_agreement.png", width = 7.5, height = 3.6)

# --- fig08: coverage by prediction source -----------------------------------
rg_cov <- read_table("legacy_rg_clean_fallback_coverage.csv")
hmsc_cov <- read_table("hmsc_joint_coverage.csv")
coverage <- dplyr::bind_rows(
  tibble::tibble(
    product = "RG + year with clean fallback (30 m)",
    source = c("RG model", "Clean + year (calibrated)")[match(rg_cov$value, 1:2)],
    rank = rg_cov$value, fraction = rg_cov$fraction
  ),
  hmsc_cov |>
    dplyr::filter(.data$source_code > 0) |>
    dplyr::transmute(
      product = "HMSC joint (120 m)",
      source = c("Enhanced", "Clean multiscale", "Clean native")[.data$source_code],
      rank = .data$source_code, fraction = .data$fraction
    )
) |>
  dplyr::mutate(
    source = factor(.data$source, levels = unique(.data$source[order(.data$product, -.data$rank)])),
    product = factor(.data$product, levels = rev(unique(.data$product)))
  )
fig08 <- ggplot(coverage, aes(x = fraction, y = product, group = rank)) +
  geom_col(aes(fill = factor(rank)), width = 0.55, colour = ink$surface, linewidth = 0.6, position = position_stack(reverse = TRUE)) +
  geom_text(aes(label = ifelse(fraction > 0.06, paste0(source, "\n", scales::percent(fraction, accuracy = 1)), "")),
    position = position_stack(vjust = 0.5, reverse = TRUE), size = 2.6, colour = "white",
    family = report_font, lineheight = 0.9
  ) +
  scale_fill_manual(values = ordinal_blue, labels = c("1" = "Primary", "2" = "First fallback", "3" = "Second fallback"), name = NULL) +
  scale_x_continuous(labels = scales::label_percent(), expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "Which model fills each cell",
    subtitle = "Fallbacks are calibrated onto the primary model's probability scale and never replace an available primary prediction.",
    x = "Share of predicted cells", y = NULL,
    caption = "HMSC enhanced = geomorphic class + ReefGuide exposure + Acropora + multiscale terrain; clean native (unlabelled, 4%) = block-mean clean stack."
  ) +
  theme_report() +
  theme(
    legend.position = "top", legend.justification = "left", panel.grid.major.y = element_blank(),
    axis.text.y = element_text(colour = ink$primary, size = 9)
  )
save_figure(fig08, cfg, "fig08_coverage_by_source.png", width = 7.5, height = 2.6)

# --- fig09: forward validation of habitat, context, history and combined ----
if (file.exists(file.path(tables_dir, "forward_validation_summary.csv"))) {
  fv <- read_table("forward_validation_summary.csv")
  all_rows <- fv |>
    dplyr::filter(.data$subset == "all") |>
    dplyr::mutate(label = stats::reorder(.data$label, .data$within_reef_year_auc))
  panel_a <- ggplot(all_rows, aes(y = label)) +
    geom_vline(xintercept = 0.5, colour = ink$muted, linetype = "22", linewidth = 0.4) +
    geom_segment(aes(x = auc_low, xend = auc_high, yend = label),
      colour = ink$secondary,
      linewidth = 1.6, lineend = "round", alpha = 0.35
    ) +
    geom_point(aes(x = within_reef_year_auc), colour = ink$primary, size = 2.6) +
    geom_text(aes(x = auc_high, label = sprintf("%.3f", within_reef_year_auc)),
      hjust = -0.3,
      size = 2.8, colour = ink$primary, family = report_font
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.2))) +
    labs(
      x = "Within-reef-year AUC, test years 2022-2026 pooled", y = NULL,
      subtitle = sprintf(
        "All %s test site-years; bars: 95%% bootstrap over reef-years",
        format(all_rows$site_years[1], big.mark = ",")
      )
    ) +
    theme_report() +
    theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(colour = ink$primary))
  key_methods <- c("habitat", "effort", "history_fixed", "combined_fixed")
  by_subset <- fv |>
    dplyr::filter(.data$method %in% key_methods, .data$subset != "all") |>
    dplyr::mutate(
      label = factor(.data$label, levels = rev(unique(fv$label[match(key_methods, fv$method)]))),
      subset = factor(
        sprintf("%s (n = %s)", .data$subset, trimws(format(.data$site_years, big.mark = ","))),
        levels = unique(sprintf("%s (n = %s)", .data$subset, trimws(format(.data$site_years, big.mark = ","))))
      )
    )
  panel_b <- ggplot(by_subset, aes(x = within_reef_year_auc, y = label)) +
    geom_vline(xintercept = 0.5, colour = ink$muted, linetype = "22", linewidth = 0.4) +
    geom_point(colour = ink$primary, size = 2.2) +
    geom_text(aes(label = sprintf("%.3f", within_reef_year_auc)),
      hjust = -0.35, size = 2.6,
      colour = ink$primary, family = report_font
    ) +
    facet_wrap(~subset, ncol = 1) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.3))) +
    labs(x = "Within-reef-year AUC", y = NULL, subtitle = "By how much Cull history lies nearby") +
    theme_report() +
    theme(
      panel.grid.major.y = element_blank(), axis.text.y = element_text(colour = ink$primary),
      strip.clip = "off"
    )
  fig09 <- (panel_a | panel_b) + plot_layout(widths = c(1.15, 1)) +
    plot_annotation(
      title = "Forecasting next season: habitat, reef context and Cull history",
      caption = paste(
        "Every score for test year T uses only data from years before T. 'Nested' settings are re-tuned on T - 1; 'fixed' (500 m, no recency decay, kappa 2) was set in advance.",
        "History and effort partly encode where operators chose to dive, so they reward revisiting known sites.",
        sep = "\n"
      ),
      theme = theme_report()
    )
  save_figure(fig09, cfg, "fig09_forward_validation.png", width = 7.5, height = 5.2)
}

# --- fig10: habitat vs history vs combined on the four reefs -----------------
history_path <- file.path(layers_dir, "cots_cull_history_kernel_120m.tif")
combined_path <- file.path(layers_dir, "cots_hmsc_history_combined_120m.tif")
if (file.exists(history_path) && file.exists(combined_path)) {
  history <- terra::rast(history_path)
  combined <- terra::rast(combined_path)
  weight_max <- 0
  rows10 <- lapply(names(reef_ids), function(slug) {
    polygon <- features[features$LABEL_ID == reef_ids[[slug]], ]
    v <- terra::vect(polygon)
    crop_mask <- function(r) terra::mask(terra::crop(r, v, snap = "out"), v, touches = TRUE)
    w <- raster_df(crop_mask(history[["history_weight"]]))
    weight_max <<- max(weight_max, w$value)
    list(
      name = polygon$GBR_NAME[1], label = reef_ids[[slug]], polygon = polygon,
      frame = sf::st_bbox(sf::st_buffer(polygon, 300)),
      habitat = raster_df(percentile_in_polygon(native$hmsc_joint, polygon)),
      history = raster_df(crop_mask(history[["history_rate"]])),
      combined = raster_df(percentile_in_polygon(combined, polygon)),
      weight = w
    )
  })
  rate_scale <- scale_fill_gradientn(
    colours = ramp_blue, limits = c(0, 1), name = "Prior high-density rate",
    labels = scales::label_percent()
  )
  weight_scale <- scale_fill_gradientn(
    colours = ramp_orange, limits = c(0, weight_max),
    name = "History weight W", labels = scales::label_number(accuracy = 1)
  )
  titles10 <- c("HMSC habitat", "Cull history", "Habitat + history", "History weight")
  panels10 <- unlist(lapply(seq_along(rows10), function(i) {
    p <- rows10[[i]]
    t <- if (i == 1) titles10 else rep(list(NULL), 4)
    list(
      map_panel(
        p$habitat, p$polygon, p$frame, percentile_scale, t[[1]],
        paste0(sub(" \\(.*", "", p$name), "\n", p$label)
      ),
      map_panel(p$history, p$polygon, p$frame, rate_scale, t[[2]]),
      map_panel(p$combined, p$polygon, p$frame, percentile_scale, t[[3]]),
      map_panel(p$weight, p$polygon, p$frame, weight_scale, t[[4]])
    )
  }), recursive = FALSE)
  fig10 <- wrap_plots(panels10, ncol = 4, guides = "collect") +
    plot_annotation(
      title = "Updating the habitat layer with where COTS have been found before",
      subtitle = "Habitat and combined layers ranked within each reef; history is the kernel-weighted rate in all prior years.",
      caption = "Blank history cells have no Cull site-years within reach of the kernel. The combined layer equals habitat where W = 0.",
      theme = theme_report()
    ) &
    theme(
      legend.position = "bottom", legend.key.width = grid::unit(0.9, "cm"),
      legend.key.height = grid::unit(0.18, "cm"), legend.title.position = "top",
      legend.spacing.x = grid::unit(0.8, "cm")
    )
  save_figure(fig10, cfg, "fig10_history_combined_four_reefs.png", width = 7.5, height = 8.7)
}

# --- fig11: operational efficiency of the primary search-planning layer -------
if (file.exists(file.path(tables_dir, "operations_cull_scenarios.csv"))) {
  sel_cols <- c(`Primary layer` = layer_colours[["legacy_rg_year2026"]],
                `History + habitat` = layer_colours[["hmsc_joint"]], Random = ink$muted)
  sel_lty <- c(`Primary layer` = "solid", `History + habitat` = "solid", Random = "22")
  mo <- read_table("operations_manta_scenarios.csv") |> dplyr::filter(.data$subset == "surveys 2022-2026")
  manta_long <- dplyr::bind_rows(
    tibble::tibble(metric = "reef_years", keep = mo$keep, selection = "Primary layer", value = mo$reef_years_detected),
    tibble::tibble(metric = "reef_years", keep = mo$keep, selection = "Random", value = mo$reef_years_detected_random),
    tibble::tibble(metric = "strong_tows", keep = mo$keep, selection = "Primary layer", value = mo$strong_tows_kept),
    tibble::tibble(metric = "strong_tows", keep = mo$keep, selection = "Random", value = mo$strong_tows_kept_random)
  )
  co <- read_table("operations_cull_scenarios.csv")
  cull_long <- dplyr::bind_rows(
    co |> dplyr::filter(.data$model == "annual_rg_year_2022_2026", .data$selection == "layer") |> dplyr::mutate(selection = "Primary layer"),
    co |> dplyr::filter(.data$model == "history_combined_forward", .data$selection == "layer") |> dplyr::mutate(selection = "History + habitat"),
    co |> dplyr::filter(.data$model == "annual_rg_year_2022_2026", .data$selection == "random") |> dplyr::mutate(selection = "Random")
  )
  headline_k <- cfg$operations$headline_keep
  curve_panel <- function(d, y, xlab, ylab, title) {
    ggplot(d, aes(x = keep, y = .data[[y]], colour = selection, linetype = selection)) +
      geom_vline(xintercept = headline_k, colour = ink$grid, linewidth = 0.5) +
      geom_line(linewidth = 0.7) +
      geom_point(size = 1.4) +
      scale_colour_manual(values = sel_cols, name = NULL, limits = names(sel_cols), drop = FALSE) +
      scale_linetype_manual(values = sel_lty, name = NULL, limits = names(sel_cols), drop = FALSE) +
      scale_x_reverse(labels = scales::label_percent()) +
      scale_y_continuous(labels = scales::label_percent(), limits = c(0, 1)) +
      labs(x = xlab, y = ylab, title = title) +
      theme_report(8) +
      theme(plot.title = element_text(size = 9))
  }
  p1 <- curve_panel(dplyr::filter(manta_long, .data$metric == "reef_years"), "value",
                    "Share of manta tows kept in each survey", "Outbreak reef-years still detected", "Manta: reef-level detection")
  p2 <- curve_panel(dplyr::filter(manta_long, .data$metric == "strong_tows"), "value",
                    "Share of manta tows kept in each survey", "Strong-evidence tows kept", "Manta: tows with COTS evidence")
  p3 <- curve_panel(cull_long, "cots_kept", "Share of cull sites kept per reef-year",
                    "COTS at the kept sites", "Cull: COTS still reached")
  p4 <- curve_panel(cull_long, "ghost_avoided", "Share of cull sites kept per reef-year",
                    "Ghost dives (0 COTS) avoided", "Cull: ghost dives avoided")
  # The cull panels carry all three series, so their legend serves the figure.
  p1 <- p1 + theme(legend.position = "none")
  p2 <- p2 + theme(legend.position = "none")
  fig11 <- (p1 | p2) / (p3 | p4) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = "What the primary layer saves compared with cutting effort at random (2022-2026)",
      subtitle = sprintf(paste(
        "Effort is cut within each survey or reef-year by skipping the lowest-ranked tows or sites; grey line: %d%% kept.",
        "Primary-layer cull scores are held out (each site excluded from the model that scored it); history + habitat is forward-only.",
        sep = "
"), round(100 * headline_k)),
      caption = "Manta is not used to train the layer. Strong evidence: direct COTS, >= 2 counted scars or legacy 'common' scars.",
      theme = theme_report()
    ) &
    theme(legend.position = "bottom")
  fig11[[1]] <- fig11[[1]] & theme(legend.position = "none")
  save_figure(fig11, cfg, "fig11_operational_efficiency.png", width = 7.5, height = 6.6)
}

# --- fig13: manta thinning vs aggregation detection (leave-reef-out CV) -------
if (file.exists(file.path(tables_dir, "manta_aggregation_cv_summary.csv"))) {
  agg <- read_table("manta_aggregation_cv_summary.csv") |>
    dplyr::filter(.data$subset == "aggregations (>= 2 strong tows), all years")
  strat_labels <- c(layer_cv = "Primary layer, leave-reef-out", layer_deployed = "Primary layer, deployed (confounded)",
                    systematic = "Systematic (every n-th tow)", random = "Random")
  strat_cols <- c(layer_cv = layer_colours[["legacy_rg_year2026"]], layer_deployed = "#f39670",
                  systematic = ink$secondary, random = ink$muted)
  strat_lty <- c(layer_cv = "solid", layer_deployed = "42", systematic = "solid", random = "22")
  strat_shape <- c(layer_cv = 16, layer_deployed = 1, systematic = 15, random = 17)
  long <- dplyr::bind_rows(lapply(names(strat_labels), function(st) tibble::tibble(
    keep = agg$keep, strategy = st, value = agg[[st]], low = agg[[paste0(st, "_low")]], high = agg[[paste0(st, "_high")]])))
  long$strategy <- factor(long$strategy, levels = names(strat_labels))
  pa <- ggplot(long, aes(keep, value, colour = strategy, linetype = strategy)) +
    geom_hline(yintercept = 0.95, colour = ink$axis, linewidth = 0.4) +
    annotate("text", x = 0.3, y = 0.95, label = "95% of aggregations", hjust = 1, vjust = -0.5, size = 2.6,
             colour = ink$secondary, family = report_font) +
    geom_ribbon(data = dplyr::filter(long, strategy == "layer_cv"), aes(ymin = low, ymax = high),
                fill = strat_cols[["layer_cv"]], alpha = 0.15, colour = NA, show.legend = FALSE) +
    geom_line(linewidth = 0.7) +
    geom_point(aes(shape = strategy), size = 1.6) +
    scale_colour_manual(values = strat_cols, labels = strat_labels, name = NULL) +
    scale_linetype_manual(values = strat_lty, labels = strat_labels, name = NULL) +
    scale_shape_manual(values = strat_shape, labels = strat_labels, name = NULL) +
    scale_x_reverse(labels = scales::label_percent()) +
    scale_y_continuous(labels = scales::label_percent(), limits = c(0.4, 1)) +
    labs(x = "Share of manta tows kept in each survey", y = "Aggregations found (>= 1 tow kept)",
         title = sprintf("Aggregations of >= 2 strong-evidence tows (n = %s, %s surveys)",
                         format(agg$aggregations[1], big.mark = ","), format(agg$surveys[1], big.mark = ","))) +
    theme_report(8) +
    theme(legend.position = "bottom", plot.title = element_text(size = 9)) +
    guides(colour = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))
  sizes <- read_table("manta_aggregation_cv_by_size.csv") |>
    dplyr::filter(.data$keep == cfg$operations$headline_keep) |>
    tidyr::pivot_longer(dplyr::all_of(names(strat_labels)), names_to = "strategy", values_to = "value") |>
    dplyr::mutate(strategy = factor(.data$strategy, levels = names(strat_labels)),
                  size_class = sprintf("%s\n(n = %s)", .data$size_class, trimws(format(.data$aggregations, big.mark = ","))),
                  size_class = factor(.data$size_class, levels = unique(.data$size_class)))
  pb <- ggplot(sizes, aes(value, size_class, colour = strategy, shape = strategy)) +
    geom_point(size = 2.2, position = position_dodge(width = 0.6)) +
    scale_colour_manual(values = strat_cols, labels = strat_labels, guide = "none") +
    scale_shape_manual(values = strat_shape, labels = strat_labels, guide = "none") +
    scale_x_continuous(labels = scales::label_percent(), limits = c(NA, 1)) +
    labs(x = sprintf("Found with %d%% of tows kept", round(100 * cfg$operations$headline_keep)), y = NULL,
         title = "By aggregation size") +
    theme_report(8) +
    theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(colour = ink$primary),
          plot.title = element_text(size = 9))
  fig13 <- (pa | pb) + plot_layout(widths = c(1.5, 1)) +
    plot_annotation(
      title = "Can manta tows be cut without missing aggregations?",
      subtitle = paste(
        "Each reef's tows are scored by a model refitted without that reef's Cull data (5-fold leave-reef-out),",
        "removing the loop where Cull dives follow manta detections. Band: 95% bootstrap over surveys.",
        sep = "\n"),
      caption = paste(
        "Aggregation: strong-evidence tows (direct COTS, >= 2 counted scars, legacy 'common' scars) within 500 m in one survey; found if any of its tows is kept.",
        "Systematic keeps evenly spaced tows around the reef; random is the exact expectation.",
        sep = "\n"),
      theme = theme_report()
    )
  save_figure(fig13, cfg, "fig13_manta_aggregation_cv.png", width = 7.5, height = 4.6)
}

manifest <- tibble::tibble(file = list.files(cfg$dirs$figures, pattern = "^fig.*\\.png$"))
write_table(manifest, cfg, "figure_manifest.csv")
