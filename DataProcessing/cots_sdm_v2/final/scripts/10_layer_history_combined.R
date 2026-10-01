# Layer 4: Hmsc habitat updated by Cull history (forecast for the next season).
#
#   history rate H and evidence weight W from every active-reef-year Cull
#   site-year up to the latest survey year, weighted by distance and recency;
#   combined = (kappa * p_habitat + W * H) / (kappa + W), so cells with no
#   nearby history keep the habitat value and well-surveyed areas move toward
#   their observed record.
# Kernel settings are the fixed configuration (cfg$history), which the forward
# validation (09b) preferred over per-year nested selection. The history
# component partly encodes where operators have chosen to dive; report it as such.

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "..", "R", "common.R"))
cfg <- initialise_final()
source(file.path(cfg$final_dir, "R", "forward.R"))
setting <- list(sigma = cfg$history$sigma_m, tau = cfg$history$tau_years, kappa = cfg$history$kappa)
log_message("Deployed setting: sigma ", setting$sigma, " m, tau ", setting$tau, " y, kappa ", setting$kappa)

sy <- readRDS(file.path(cfg$exploration_dir, "artifacts", "site_year_first_visit.rds")) |>
  dplyr::filter(!is.na(.data$high_density)) |>
  dplyr::group_by(.data$reef_name, .data$year) |>
  dplyr::filter(any(.data$high_density == 1)) |>
  dplyr::ungroup()
forecast_year <- max(sy$year) + 1L
xy <- sf::st_coordinates(sf::st_transform(
  sf::st_as_sf(sy, coords = c("longitude", "latitude"), crs = 4326), 7855))

habitat <- terra::rast(file.path(cfg$dirs$layers, "cots_hmsc_joint_120m.tif"))
cells <- which(!is.na(terra::values(habitat, mat = FALSE)))
cell_xy <- terra::xyFromCell(habitat, cells)
p_habitat <- terra::values(habitat, mat = FALSE)[cells]

# Kernel only where some site lies within 3 sigma; elsewhere W = 0.
near <- as.numeric(FNN::get.knnx(xy, cell_xy, k = 1)$nn.dist) <= 3 * setting$sigma
H <- rep(NA_real_, length(cells))
W <- rep(0, length(cells))
idx <- which(near)
for (chunk in split(idx, ceiling(seq_along(idx) / 100000L))) {
  h <- history_kernel(xy, sy$high_density, (forecast_year - 1L) - sy$year, cell_xy[chunk, , drop = FALSE],
                      setting$sigma, setting$tau)
  H[chunk] <- h$H
  W[chunk] <- h$W
}
combined <- combine_history(p_habitat, H, W, setting$kappa)
log_message(sprintf("History within reach of %.1f%% of cells; W > 1 at %.1f%%",
                    100 * mean(W > 0), 100 * mean(W > 1)))

to_raster <- function(values, name) {
  x <- terra::rast(habitat)
  x[cells] <- values
  names(x) <- name
  x
}
history_path <- file.path(cfg$dirs$layers, "cots_cull_history_kernel_120m.tif")
combined_path <- file.path(cfg$dirs$layers, "cots_hmsc_history_combined_120m.tif")
write_layer(c(to_raster(H, "history_rate"), to_raster(W, "history_weight")), history_path)
write_layer(to_raster(combined, "cots_prob"), combined_path)
register_layer(
  cfg, "cull_history_kernel", history_path,
  description = sprintf("Band 1: kernel-weighted prior high-density rate (active reef-years <= %d); band 2: evidence weight W", forecast_year - 1L),
  model = sprintf("Gaussian kernel sigma %g m, recency tau %g y", setting$sigma, setting$tau),
  units = "rate / weight"
)
register_layer(
  cfg, "hmsc_history_combined", combined_path,
  description = sprintf("HMSC habitat updated by Cull history, forecast for %d: (kappa p + W H) / (kappa + W)", forecast_year),
  model = sprintf("HMSC joint + history kernel (kappa %g); partly encodes operator targeting", setting$kappa),
  units = "probability",
  extra = list(forecast_year = forecast_year)
)
write_table(tibble::tibble(
  forecast_year = forecast_year, sigma_m = setting$sigma, tau_years = setting$tau, kappa = setting$kappa,
  history_site_years = nrow(sy), cells = length(cells),
  fraction_with_history = mean(W > 0), fraction_w_gt_1 = mean(W > 1),
  fraction_history_dominant = mean(W > setting$kappa)
), cfg, "history_combined_layer_summary.csv")

# --- ReefScan presence-only check -------------------------------------------
# Each detection's within-reef percentile (0 = lowest cell on its reef,
# 1 = highest); 0.5 is chance. ReefScan is not used by any layer.
reefscan_files <- list.files(file.path(cfg$repo_root, "data", "ReefScan"), pattern = "\\.csv$", full.names = TRUE)
if (length(reefscan_files)) {
  rs <- dplyr::bind_rows(lapply(reefscan_files, utils::read.csv)) |>
    dplyr::distinct(.data$latitude, .data$longitude, .keep_all = TRUE)
  rs_xy <- sf::st_coordinates(sf::st_transform(sf::st_as_sf(rs, coords = c("longitude", "latitude"), crs = 4326), 7855))
  reef_id <- terra::rast(file.path(cfg$dirs$cache, "reef_id_expanded_120m.tif"))
  support <- c(
    terra::rast(file.path(cfg$dirs$cache, "legacy_clean_year2026_120m.tif")),
    terra::rast(file.path(cfg$dirs$cache, "legacy_rg_year2026_120m.tif")),
    terra::crop(terra::extend(habitat, reef_id), reef_id),
    terra::crop(terra::extend(terra::rast(combined_path), reef_id), reef_id)
  )
  names(support) <- c("legacy_clean_year2026", "legacy_rg_year2026", "hmsc_joint", "hmsc_history_combined")
  rs_cell <- terra::cellFromXY(reef_id, rs_xy)
  id_values <- terra::values(reef_id, mat = FALSE)
  rs_reef <- id_values[rs_cell]
  in_detected_reefs <- which(id_values %in% unique(stats::na.omit(rs_reef)))
  reef_cells <- split(in_detected_reefs, id_values[in_detected_reefs])
  rm(id_values)
  check <- dplyr::bind_rows(lapply(names(support), function(layer) {
    values <- terra::values(support[[layer]], mat = FALSE)
    pct <- vapply(seq_along(rs_cell), function(i) {
      if (is.na(rs_reef[i])) return(NA_real_)
      v <- values[reef_cells[[as.character(rs_reef[i])]]]
      x <- values[rs_cell[i]]
      if (!is.finite(x) || sum(is.finite(v)) < 10) return(NA_real_)
      mean(v[is.finite(v)] < x) + 0.5 * mean(v[is.finite(v)] == x)
    }, numeric(1))
    ok <- is.finite(pct)
    tibble::tibble(layer = layer, detections = sum(ok), reefs = dplyr::n_distinct(rs_reef[ok]),
                   mean_within_reef_percentile = mean(pct[ok]),
                   wilcoxon_p_vs_0.5 = if (sum(ok) > 5) stats::wilcox.test(pct[ok], mu = 0.5)$p.value else NA_real_)
  }))
  write_table(check, cfg, "reefscan_presence_check.csv")
  print(check)
}
