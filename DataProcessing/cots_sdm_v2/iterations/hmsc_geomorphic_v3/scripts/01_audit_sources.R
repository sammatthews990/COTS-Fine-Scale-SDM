args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v3_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v3_dir, "R", "bootstrap.R"))
ctx <- initialise_v3(script_file)
cfg <- ctx$config

manta <- read_manta_observations(cfg$paths$manta_current_workbook, sheet = "Manta") |>
  dplyr::mutate(visit_id = .data$manta_record_id)
salad <- read_salad_tracks(cfg$paths$salad_workbook)

geometry_audit <- dplyr::bind_rows(
  track_geometry_audit(manta, "manta"),
  track_geometry_audit(salad, "salad")
)
utils::write.csv(geometry_audit, file.path(
  cfg$paths$outputs_dir, "source_track_geometry_audit.csv"
), row.names = FALSE)

events <- interpolate_salad_records(cfg$paths$salad_workbook, salad)
event_audit <- tibble::tibble(
  total_events = nrow(events),
  matched_track = sum(!is.na(events$visit_id)),
  valid_time_interpolation = sum(events$interpolation_valid),
  direct_events = sum(events$evidence == "cots", na.rm = TRUE),
  scar_events = sum(events$evidence == "scar", na.rm = TRUE)
)
utils::write.csv(event_audit, file.path(
  cfg$paths$outputs_dir, "salad_event_position_audit.csv"
), row.names = FALSE)

read_outbreak <- function(path, layer) {
  x <- readr::read_csv(path, show_col_types = FALSE)
  tibble::tibble(
    layer = layer, rows = nrow(x), reefs = dplyr::n_distinct(x$reefName),
    year_min = min(x$year), year_max = max(x$year),
    missing_probability = sum(!is.finite(x$outbrProb)),
    probability_min = min(x$outbrProb, na.rm = TRUE),
    probability_max = max(x$outbrProb, na.rm = TRUE),
    duplicate_reef_year = sum(duplicated(paste(x$reefName, x$year)))
  )
}
outbreak_audit <- dplyr::bind_rows(
  read_outbreak(cfg$paths$outbreak_hindcast_file, "hindcast"),
  read_outbreak(cfg$paths$outbreak_forecast_file, "forecast")
)
utils::write.csv(outbreak_audit, file.path(
  cfg$paths$outputs_dir, "outbreak_layer_schema_audit.csv"
), row.names = FALSE)

packages <- c("Hmsc", "unmarked", "glmmTMB", "glmnet", "sf", "terra", "exactextractr")
package_audit <- tibble::tibble(
  package = packages,
  installed = vapply(packages, requireNamespace, logical(1), quietly = TRUE)
)
utils::write.csv(package_audit, file.path(
  cfg$paths$outputs_dir, "modelling_package_audit.csv"
), row.names = FALSE)

saveRDS(list(
  created = Sys.time(), geometry = geometry_audit, salad_events = event_audit,
  outbreak = outbreak_audit, packages = package_audit,
  policy = list(
    legacy_role = "benchmark_only",
    legacy_used_in_training = FALSE,
    outbreak_role = "separate_activation_layer",
    outbreak_used_in_habitat = FALSE
  )
), file.path(cfg$paths$artifacts_dir, "v3_source_audit.rds"))

print(geometry_audit)
print(event_audit)
print(outbreak_audit)
