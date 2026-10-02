# Put the three delivered layers on a common 120 m support and compare them.
#
#   1. Legacy 30 m layers -> 120 m by the same 4x4 block mean used for Hmsc
#      predictors, on the same 0-anchored grid.
#   2. Within-reef agreement: Spearman correlation between layers inside each
#      GBR reef polygon.
#   3. Independent Manta check: within reef-year, do strong-evidence tows
#      (direct COTS, >= 2 counted scars, or legacy "common" scars) score higher
#      than strict-negative tows? Manta was not a response in any of the three
#      models, so this is a method-independent check of within-reef ranking.
#      Tow scores are length-weighted means over the 120 m cells each tow crosses.

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "..", "R", "common.R"))
cfg <- initialise_final()
layer_ids <- names(layer_labels)
template <- terra::rast(file.path(cfg$dirs$cache, "base_support_120m.tif"))[[1]]

# --- 1. Common 120 m layers -------------------------------------------------
to_support_120 <- function(path, out) aggregate_to_support_120(path, out, template)
hmsc <- terra::rast(file.path(cfg$dirs$layers, "cots_hmsc_joint_120m.tif"))
hmsc <- terra::crop(terra::extend(hmsc, template), template)
support_layers <- c(
  to_support_120(file.path(cfg$dirs$layers, "cots_legacy_clean_year2026_30m.tif"),
                 file.path(cfg$dirs$cache, "legacy_clean_year2026_120m.tif")),
  to_support_120(file.path(cfg$dirs$layers, "cots_legacy_rg_year2026_30m.tif"),
                 file.path(cfg$dirs$cache, "legacy_rg_year2026_120m.tif")),
  hmsc
)
names(support_layers) <- layer_ids

# --- 2. Within-reef agreement -----------------------------------------------
reefs <- sf::st_read(checked_input(cfg, "reef_features"), quiet = TRUE) |>
  dplyr::filter(.data$FEAT_NAME == "Reef") |>
  sf::st_transform(terra::crs(template))
reefs$reef_index <- seq_len(nrow(reefs))
reef_id_path <- file.path(cfg$dirs$cache, "reef_id_120m.tif")
if (!file.exists(reef_id_path)) {
  terra::rasterize(terra::vect(reefs), template, field = "reef_index",
                   filename = reef_id_path, overwrite = TRUE, wopt = list(datatype = "INT2U"))
}
# reef_index starts at 1, so 0 can only be an unflagged background value.
reef_id <- terra::classify(terra::rast(reef_id_path), cbind(0, NA))
reef_cells <- which(!is.na(terra::values(reef_id, mat = FALSE)))
cells <- data.frame(
  reef_index = terra::values(reef_id, mat = FALSE)[reef_cells],
  terra::extract(support_layers, reef_cells)
)
pairs <- utils::combn(layer_ids, 2, simplify = FALSE)
agreement <- dplyr::bind_rows(lapply(split(cells, cells$reef_index), function(d) {
  dplyr::bind_rows(lapply(pairs, function(p) {
    ok <- is.finite(d[[p[1]]]) & is.finite(d[[p[2]]])
    tibble::tibble(
      reef_index = d$reef_index[1], layer_a = p[1], layer_b = p[2], cells = sum(ok),
      spearman = if (sum(ok) >= 30) stats::cor(d[[p[1]]][ok], d[[p[2]]][ok], method = "spearman") else NA_real_
    )
  }))
})) |>
  dplyr::left_join(sf::st_drop_geometry(reefs)[, c("reef_index", "GBR_NAME", "LABEL_ID")], by = "reef_index")
write_table(agreement, cfg, "layer_agreement_within_reef.csv")

coverage_on_reefs <- tibble::tibble(
  layer_id = layer_ids,
  reef_cells = length(reef_cells),
  covered_cells = vapply(layer_ids, function(l) sum(is.finite(cells[[l]])), numeric(1))
) |> dplyr::mutate(fraction = .data$covered_cells / .data$reef_cells)
write_table(coverage_on_reefs, cfg, "layer_coverage_on_reef_polygons.csv")

# --- 3. Independent Manta check ---------------------------------------------
tows <- readRDS(checked_input(cfg, "manta_tows")) |>
  dplyr::select("visit_id", "year", "reef_key", "reef_name", "reef_year_id",
                "manta_evidence_class", "active_primary")
links <- readRDS(checked_input(cfg, "footprint_links"))$manta |>
  dplyr::select("visit_id", "x_bin", "y_bin", "weighted_length_m", "fraction") |>
  dplyr::mutate(
    weight = dplyr::case_when(
      is.finite(.data$weighted_length_m) & .data$weighted_length_m > 0 ~ .data$weighted_length_m,
      is.finite(.data$fraction) & .data$fraction > 0 ~ .data$fraction,
      TRUE ~ 1
    ),
    cell = terra::cellFromXY(template, cbind((.data$x_bin + 0.5) * 120, (.data$y_bin + 0.5) * 120))
  )
link_values <- terra::extract(support_layers, links$cell)
tow_scores <- dplyr::bind_rows(lapply(layer_ids, function(l) {
  v <- link_values[[l]]
  ok <- is.finite(v)
  tibble::tibble(visit_id = links$visit_id, w = links$weight, v = v, ok = ok) |>
    dplyr::group_by(.data$visit_id) |>
    dplyr::summarise(
      coverage = sum(.data$w[.data$ok]) / sum(.data$w),
      score = if (any(.data$ok)) stats::weighted.mean(.data$v[.data$ok], .data$w[.data$ok]) else NA_real_,
      .groups = "drop"
    ) |>
    dplyr::mutate(layer_id = l)
}))
tow_scores <- tow_scores |>
  dplyr::mutate(score = dplyr::if_else(.data$coverage >= 0.5, .data$score, NA_real_)) |>
  dplyr::select("visit_id", "layer_id", "score") |>
  tidyr::pivot_wider(names_from = "layer_id", values_from = "score")

strong_classes <- c("direct_cots", "multiple_counted_scars", "legacy_common_scars")
evaluation <- tows |>
  dplyr::filter(.data$active_primary,
                .data$manta_evidence_class %in% c(strong_classes, "strict_negative")) |>
  dplyr::mutate(positive = .data$manta_evidence_class %in% strong_classes) |>
  dplyr::inner_join(tow_scores, by = "visit_id")
matched <- evaluation[stats::complete.cases(evaluation[, layer_ids]), ]
# Keep reef-years that contain both outcomes.
mixed <- matched |>
  dplyr::group_by(.data$reef_year_id) |>
  dplyr::filter(any(.data$positive), any(!.data$positive)) |>
  dplyr::ungroup()

set.seed(cfg$hmsc$seed)
groups <- split(seq_len(nrow(mixed)), mixed$reef_year_id)
boot_index <- replicate(500, sample(length(groups), replace = TRUE), simplify = FALSE)
manta_metrics <- dplyr::bind_rows(lapply(layer_ids, function(l) {
  score <- mixed[[l]]
  boot <- vapply(boot_index, function(g) {
    i <- unlist(groups[g], use.names = FALSE)
    grp <- rep(seq_along(g), lengths(groups[g]))
    c(within_group_auc(score[i], mixed$positive[i], grp),
      within_group_top20_lift(score[i], mixed$positive[i], grp))
  }, numeric(2))
  tibble::tibble(
    layer_id = l, layer = layer_labels[[l]],
    tows = nrow(mixed), positives = sum(mixed$positive),
    reef_years = length(groups), reefs = dplyr::n_distinct(mixed$reef_key),
    within_reef_year_auc = within_group_auc(score, mixed$positive, mixed$reef_year_id),
    auc_low = stats::quantile(boot[1, ], 0.025), auc_high = stats::quantile(boot[1, ], 0.975),
    top20_lift = within_group_top20_lift(score, mixed$positive, mixed$reef_year_id),
    lift_low = stats::quantile(boot[2, ], 0.025), lift_high = stats::quantile(boot[2, ], 0.975)
  )
}))
write_table(manta_metrics, cfg, "manta_independent_check.csv")
print(manta_metrics)
