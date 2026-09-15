args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v4_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
v2_dir <- normalizePath(file.path(v4_dir, "..", ".."), winslash = "/")
source(file.path(v4_dir, "R", "bootstrap.R"))
source(file.path(v2_dir, "R", "io.R"))
source(file.path(v2_dir, "R", "joint_visits.R"))
ctx <- initialise_v4(script_file)
cfg <- ctx$config
support_m <- 200L
template <- terra::rast(cfg$paths$predictor_stack)
target_crs <- terra::crs(template)

manta <- read_manta_observations(cfg$paths$manta_current_workbook, sheet = "Manta") |>
  dplyr::transmute(
    visit_id = .data$manta_record_id, method = "manta",
    survey_date, year, reef_name, reef_key = normalise_reef_key(.data$reef_name),
    start_longitude, start_latitude, end_longitude, end_latitude,
    tow_distance_m, direct_cots_count = .data$cots_count,
    scar_count, scar_code, strict_negative, informative, hard_coral_observed
  ) |>
  dplyr::filter(valid_lonlat(.data$start_longitude, .data$start_latitude))
salad <- read_salad_tracks(cfg$paths$salad_workbook)
cull_raw <- read_cull_observations(cfg$paths$cull_workbook)
cull <- dplyr::transmute(
  cull_raw,
  visit_id = paste0("cull||", .data$CrownOfThornsStarfishCullDiveId, "||", .data$source_record_id),
  method = "cull", survey_date, year,
  reef_label = as.character(.data$ReefLabel), reef_name,
  reef_key = normalise_reef_key(.data$reef_name), longitude, latitude,
  direct_cots_count = .data$total_cots, bottom_time = .data$bottom_time,
  cohort1 = as.numeric(.data$Cohort1), cohort2 = as.numeric(.data$Cohort2),
  cohort3 = as.numeric(.data$Cohort3), cohort4 = as.numeric(.data$Cohort4),
  cull_dive_type = as.character(.data$CullDiveTypeName)
)
rhis <- read_rhis_observations(cfg$paths$cull_workbook)

message("Building exact manta/SALAD and point Cull/RHIS links at 200 m.")
links <- list(
  manta = line_cell_links(manta, target_crs, support_m, reported_length = "tow_distance_m"),
  salad = line_cell_links(salad, target_crs, support_m, reported_length = "footprint_length_m") |>
    dplyr::mutate(weighted_area_m2 = .data$fraction * .data$footprint_area_m2),
  cull = point_cell_links(cull, target_crs, support_m),
  rhis = point_cell_links(rhis, target_crs, support_m)
)
saveRDS(links, file.path(cfg$paths$artifacts_dir, "footprint_links_200m.rds"), compress = "xz")
audit <- tibble::tibble(
  support_m = support_m, method = names(links),
  source_visits = vapply(links, function(x) dplyr::n_distinct(x$visit_id), integer(1)),
  footprint_cell_links = vapply(links, nrow, integer(1)),
  cells = vapply(links, function(x) dplyr::n_distinct(x$cell_id), integer(1)),
  fallback_links = vapply(links, function(x) sum(x$geometry_fallback), integer(1))
)
utils::write.csv(audit, file.path(cfg$paths$outputs_dir, "footprint_link_audit_200m.csv"), row.names = FALSE)
print(audit)
