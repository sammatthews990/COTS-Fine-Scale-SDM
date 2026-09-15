args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v6_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v6_dir, "R", "bootstrap.R"))
ctx <- initialise_density_v6(script_file)

target <- file.path(ctx$paths$artifacts, "density_data_120m.rds")
if (!file.exists(target)) {
  artifact <- prepare_density_support_v6(120L, ctx)
  saveRDS(artifact, target, compress = "xz")
}

audits <- lapply(c(60L, 120L), function(support_m) {
  artifact <- readRDS(file.path(ctx$paths$artifacts, sprintf("density_data_%dm.rds", support_m)))
  anomalies <- artifact$salad |>
    dplyr::filter(abs(.data$observed_density_ha - .data$count / .data$effort_area_ha) > 1e-8)
  tibble::tibble(
    support_m = support_m,
    enhanced_predictor_cells = nrow(artifact$design$cells),
    eligible_active_cull_visits = nrow(artifact$cull),
    cull_zero_fraction = mean(artifact$cull$count == 0),
    cull_year_min = min(artifact$cull$year), cull_year_max = max(artifact$cull$year),
    retained_salad_tracks = nrow(artifact$salad),
    excluded_salad_tracks_for_coverage = sum(artifact$coverage$predictor_coverage < artifact$minimum_track_coverage),
    salad_zero_fraction = mean(artifact$salad$count == 0),
    salad_year_min = min(artifact$salad$year), salad_year_max = max(artifact$salad$year),
    salad_density_reconciliation_failures = nrow(anomalies),
    median_salad_predictor_coverage = stats::median(artifact$coverage$predictor_coverage)
  )
})
audit <- dplyr::bind_rows(audits)
utils::write.csv(audit, file.path(ctx$paths$outputs, "density_input_audit.csv"), row.names = FALSE)
print(audit)
