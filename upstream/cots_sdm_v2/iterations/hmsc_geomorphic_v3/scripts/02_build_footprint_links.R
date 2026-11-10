args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v3_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v3_dir, "R", "bootstrap.R"))
ctx <- initialise_v3(script_file)
cfg <- ctx$config
template <- terra::rast(cfg$paths$predictor_stack)
target_crs <- terra::crs(template)

manta <- read_manta_observations(cfg$paths$manta_current_workbook, sheet = "Manta") |>
  dplyr::transmute(
    visit_id = .data$manta_record_id, method = "manta",
    survey_date, year, reef_name, reef_key = normalise_reef_key(.data$reef_name),
    start_longitude, start_latitude, end_longitude, end_latitude,
    tow_distance_m, direct_cots_count = .data$cots_count,
    scar_count, scar_code, strict_negative, informative,
    hard_coral_observed
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

build_support <- function(cell_size_m) {
  message("Building exact footprint links at ", cell_size_m, " m.")
  manta_links <- line_cell_links(
    manta, target_crs, cell_size_m, reported_length = "tow_distance_m"
  )
  salad_links <- line_cell_links(
    salad, target_crs, cell_size_m, reported_length = "footprint_length_m"
  ) |>
    dplyr::mutate(weighted_area_m2 = .data$fraction * .data$footprint_area_m2)
  cull_links <- point_cell_links(cull, target_crs, cell_size_m)
  rhis_links <- point_cell_links(rhis, target_crs, cell_size_m)
  result <- list(manta = manta_links, salad = salad_links, cull = cull_links, rhis = rhis_links)
  saveRDS(result, file.path(
    cfg$paths$artifacts_dir, paste0("footprint_links_", cell_size_m, "m.rds")
  ), compress = "xz")
  tibble::tibble(
    support_m = cell_size_m,
    method = names(result),
    source_visits = vapply(result, function(x) dplyr::n_distinct(x$visit_id), integer(1)),
    footprint_cell_links = vapply(result, nrow, integer(1)),
    cells = vapply(result, function(x) dplyr::n_distinct(x$cell_id), integer(1)),
    fallback_links = vapply(result, function(x) sum(x$geometry_fallback), integer(1)),
    endpoint_mismatch_visits = vapply(result, function(x) {
      if (!"geometry_quality_reason" %in% names(x)) return(0L)
      dplyr::n_distinct(x$visit_id[x$geometry_quality_reason == "endpoint_length_mismatch"])
    }, integer(1))
  )
}

audit <- dplyr::bind_rows(lapply(cfg$validation$cell_support_m, build_support))
utils::write.csv(audit, file.path(
  cfg$paths$outputs_dir, "footprint_link_audit.csv"
), row.names = FALSE)
print(audit)
