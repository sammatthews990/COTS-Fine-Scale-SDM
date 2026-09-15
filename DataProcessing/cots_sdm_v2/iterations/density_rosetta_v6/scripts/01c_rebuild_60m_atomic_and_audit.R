args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v6_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v6_dir, "R", "bootstrap.R"))
ctx <- initialise_density_v6(script_file)

artifact <- prepare_density_support_v6(60L, ctx)
target <- file.path(ctx$paths$artifacts, "density_data_60m.rds")
temporary <- tempfile("density_data_60m_", tmpdir = ctx$paths$artifacts, fileext = ".rds")
saveRDS(artifact, temporary, compress = "gzip")
check <- readRDS(temporary)
stopifnot(check$support_m == 60L, nrow(check$cull) > 0L, nrow(check$salad) > 0L)
if (!file.rename(temporary, target)) stop("Could not atomically replace 60 m artifact.")

audits <- lapply(c(60L, 120L), function(support_m) {
  x <- readRDS(file.path(ctx$paths$artifacts, sprintf("density_data_%dm.rds", support_m)))
  tibble::tibble(
    support_m = support_m,
    enhanced_predictor_cells = nrow(x$design$cells),
    eligible_active_cull_visits = nrow(x$cull),
    cull_zero_fraction = mean(x$cull$count == 0),
    cull_year_min = min(x$cull$year), cull_year_max = max(x$cull$year),
    retained_salad_tracks = nrow(x$salad),
    excluded_salad_tracks_for_coverage = sum(x$coverage$predictor_coverage < x$minimum_track_coverage),
    salad_zero_fraction = mean(x$salad$count == 0),
    salad_year_min = min(x$salad$year), salad_year_max = max(x$salad$year),
    salad_density_reconciliation_failures = sum(abs(x$salad$observed_density_ha - x$salad$count / x$salad$effort_area_ha) > 1e-8),
    median_salad_predictor_coverage = stats::median(x$coverage$predictor_coverage)
  )
})
audit <- dplyr::bind_rows(audits)
utils::write.csv(audit, file.path(ctx$paths$outputs, "density_input_audit.csv"), row.names = FALSE)
print(audit)
