# Primary search-planning layer: RG + year with clean + year fallback, 30 m.
#
# Built from 02's composite. Adds the operational forms used in the field:
#   cots_search_planning_2026_30m.tif           band 1 probability, band 2 within-reef percentile
#   cots_search_planning_priority_2026_30m.tif  1 = top 10 % of the reef, 2 = 10-20 %, 3 = 20-30 %, 0 = rest
# Percentiles are computed within each GBR reef polygon, so a reef's top cells
# are flagged whatever that reef's overall probability. The composite mixes
# two models on reefs that straddle the RG coverage boundary; the per-reef
# source mix is written to tables/search_planning_reef_sources.csv.

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "..", "R", "common.R"))
cfg <- initialise_final()

composite <- terra::rast(file.path(cfg$dirs$layers, "cots_legacy_rg_year2026_clean_fallback_30m.tif"))
source_layer <- terra::rast(file.path(cfg$dirs$layers, "cots_legacy_rg_year2026_clean_fallback_source_30m.tif"))

reefs <- sf::st_read(checked_input(cfg, "reef_features"), quiet = TRUE) |>
  dplyr::filter(.data$FEAT_NAME == "Reef") |>
  sf::st_transform(terra::crs(composite))
reefs$reef_index <- seq_len(nrow(reefs))
reef_id_path <- file.path(cfg$dirs$cache, "reef_id_30m.tif")
if (!file.exists(reef_id_path)) {
  log_message("Rasterising ", nrow(reefs), " reef polygons at 30 m")
  terra::rasterize(terra::vect(reefs), composite, field = "reef_index", filename = reef_id_path,
                   overwrite = TRUE, wopt = list(datatype = "INT2U", NAflag = 65535))
}
reef_id <- terra::rast(reef_id_path)

# Pass 1: collect (cell, reef, probability, source) for every reef cell.
log_message("Collecting reef cells")
blocks <- terra::blocks(composite, n = 8)
terra::readStart(composite); terra::readStart(reef_id); terra::readStart(source_layer)
pieces <- vector("list", blocks$n)
for (b in seq_len(blocks$n)) {
  id <- terra::readValues(reef_id, blocks$row[b], blocks$nrows[b])
  keep <- which(!is.na(id) & id > 0 & id < 65535)
  if (!length(keep)) next
  p <- terra::readValues(composite, blocks$row[b], blocks$nrows[b])[keep]
  s <- terra::readValues(source_layer, blocks$row[b], blocks$nrows[b])[keep]
  ok <- is.finite(p)
  pieces[[b]] <- data.frame(
    cell = (blocks$row[b] - 1) * terra::ncol(composite) + keep[ok],
    reef = id[keep][ok], p = p[ok], src = s[ok]
  )
}
terra::readStop(composite); terra::readStop(reef_id); terra::readStop(source_layer)
cells <- dplyr::bind_rows(pieces)
log_message(format(nrow(cells), big.mark = ","), " reef cells with a prediction")

cells <- cells |>
  dplyr::group_by(.data$reef) |>
  dplyr::mutate(
    n = dplyr::n(),
    percentile = if (dplyr::n() > 1) (rank(.data$p, ties.method = "average") - 1) / (dplyr::n() - 1) else NA_real_
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(priority = dplyr::case_when(
    .data$percentile >= 0.9 ~ 1L, .data$percentile >= 0.8 ~ 2L, .data$percentile >= 0.7 ~ 3L, TRUE ~ 0L
  ))

# Pass 2: write percentile and priority rasters block by block.
write_blocks <- function(template, values_by_cell, path, datatype, name) {
  out <- terra::rast(template)
  names(out) <- name
  options <- if (datatype == "INT1U") c("COMPRESS=DEFLATE", "TILED=YES", "BIGTIFF=IF_SAFER") else layer_gdal_options
  b <- terra::writeStart(out, path, overwrite = TRUE, datatype = datatype, gdal = options,
                         NAflag = if (datatype == "INT1U") 255 else NA)
  nc <- terra::ncol(template)
  for (i in seq_len(b$n)) {
    first <- (b$row[i] - 1) * nc + 1
    last <- (b$row[i] + b$nrows[i] - 1) * nc
    v <- rep(NA_real_, last - first + 1)
    sel <- which(values_by_cell$cell >= first & values_by_cell$cell <= last)
    v[values_by_cell$cell[sel] - first + 1] <- values_by_cell$value[sel]
    terra::writeValues(out, v, b$row[i], b$nrows[i])
  }
  terra::writeStop(out)
}
cells <- cells[order(cells$cell), ]
percentile_path <- file.path(cfg$dirs$cache, "search_planning_percentile_30m.tif")
priority_path <- file.path(cfg$dirs$layers, "cots_search_planning_priority_2026_30m.tif")
log_message("Writing within-reef percentile")
write_blocks(composite, data.frame(cell = cells$cell, value = cells$percentile), percentile_path, "FLT4S", "within_reef_percentile")
log_message("Writing priority classes")
write_blocks(composite, data.frame(cell = cells$cell, value = cells$priority), priority_path, "INT1U", "priority_class")

layer_path <- file.path(cfg$dirs$layers, "cots_search_planning_2026_30m.tif")
stack <- c(composite, terra::rast(percentile_path))
names(stack) <- c("cots_prob", "within_reef_percentile")
write_layer(stack, layer_path)

reef_sources <- cells |>
  dplyr::group_by(.data$reef) |>
  dplyr::summarise(cells = dplyr::n(), rg_fraction = mean(.data$src == 1), .groups = "drop") |>
  dplyr::left_join(sf::st_drop_geometry(reefs)[, c("reef_index", "GBR_NAME", "LABEL_ID")], by = c("reef" = "reef_index")) |>
  dplyr::mutate(source = dplyr::case_when(.data$rg_fraction == 1 ~ "RG only", .data$rg_fraction == 0 ~ "clean fallback only", TRUE ~ "mixed"))
write_table(reef_sources, cfg, "search_planning_reef_sources.csv")

register_layer(
  cfg, "search_planning_2026", layer_path,
  description = "PRIMARY search-planning layer. Band 1: P(CPUE >= 0.02) in 2026 (RG + year, clean + year fallback); band 2: percentile within the GBR reef polygon",
  model = "legacy RG + year with clean + year fallback",
  units = "probability / percentile",
  extra = list(reefs = nrow(reef_sources),
               reefs_rg_only = sum(reef_sources$source == "RG only"),
               reefs_mixed = sum(reef_sources$source == "mixed"))
)
register_layer(
  cfg, "search_planning_priority_2026", priority_path,
  description = "Within-reef priority: 1 = top 10 %, 2 = 10-20 %, 3 = 20-30 %, 0 = remainder of the reef",
  model = "legacy RG + year with clean + year fallback",
  units = "category"
)
print(table(reef_sources$source))
