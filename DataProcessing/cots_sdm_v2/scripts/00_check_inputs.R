args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v2_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v2_dir, "R", "bootstrap.R"))
ctx <- initialise_v2(script_file)
cfg <- ctx$config

required_packages <- c(
  "readxl", "dplyr", "lubridate", "sf", "terra", "tibble",
  "tidyr", "readr", "yardstick", "xgboost"
)
package_status <- data.frame(
  package = required_packages,
  available = vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
)
if (!all(package_status$available)) {
  stop(
    "Missing required packages: ",
    paste(package_status$package[!package_status$available], collapse = ", ")
  )
}

input_paths <- c(
  cull_workbook = cfg$paths$cull_workbook,
  predictor_stack = cfg$paths$predictor_stack,
  outbreak_hindcast = cfg$paths$outbreak_hindcast_file,
  outbreak_forecast = cfg$paths$outbreak_forecast_file,
  manta_history = cfg$paths$manta_history_csv,
  manta_current = cfg$paths$manta_current_workbook
)
input_status <- data.frame(
  input = names(input_paths),
  path = unname(input_paths),
  exists = file.exists(input_paths),
  size_bytes = ifelse(
    file.exists(input_paths),
    file.info(input_paths)$size,
    NA_real_
  )
)
if (!all(input_status$exists)) {
  stop(
    "Missing configured inputs: ",
    paste(input_status$path[!input_status$exists], collapse = ", ")
  )
}

observations <- read_cull_observations(cfg$paths$cull_workbook)
hindcast <- read_outbreak_probabilities(
  cfg$paths$outbreak_hindcast_file, "adjusted_hindcast"
)
forecast <- read_outbreak_probabilities(
  cfg$paths$outbreak_forecast_file, "adjusted_forecast"
)

if (!all(forecast$year == cfg$outbreak$forecast_year)) {
  stop("Forecast input contains a year other than the configured forecast year.")
}
if (!setequal(hindcast$reef_name, forecast$reef_name)) {
  stop("Adjusted hindcast and forecast reef sets do not match.")
}

cull_reefs <- sort(unique(observations$reef_name))
outbreak_reefs <- sort(unique(hindcast$reef_name))
unmatched <- setdiff(cull_reefs, outbreak_reefs)

audit <- list(
  generated_at = Sys.time(),
  config_version = cfg$version,
  packages = package_status,
  inputs = input_status,
  cull = list(
    rows = nrow(observations),
    sites = dplyr::n_distinct(observations$site_id),
    reefs = length(cull_reefs),
    min_year = min(observations$year),
    max_year = max(observations$year)
  ),
  outbreak_hindcast = list(
    rows = nrow(hindcast),
    reefs = length(outbreak_reefs),
    min_year = min(hindcast$year),
    max_year = max(hindcast$year),
    duplicate_reef_years = sum(duplicated(
      hindcast[c("reef_name", "year")]
    ))
  ),
  outbreak_forecast = list(
    rows = nrow(forecast),
    reefs = dplyr::n_distinct(forecast$reef_name),
    year = unique(forecast$year),
    cutoff = unique(stats::na.omit(
      forecast$outbreak_probability_cutoff
    )),
    duplicate_reef_years = sum(duplicated(
      forecast[c("reef_name", "year")]
    ))
  ),
  join = list(
    matched_cull_reefs = sum(cull_reefs %in% outbreak_reefs),
    match_rate = mean(cull_reefs %in% outbreak_reefs),
    unmatched_cull_reefs = unmatched
  )
)

saveRDS(audit, file.path(cfg$paths$artifacts_dir, "input_audit.rds"))
utils::write.csv(
  input_status,
  file.path(cfg$paths$artifacts_dir, "input_status.csv"),
  row.names = FALSE
)

cat("COTS SDM v2 input audit passed.\n")
cat(sprintf(
  "Cull: %d records, %d sites, %d reefs, %d-%d.\n",
  audit$cull$rows, audit$cull$sites, audit$cull$reefs,
  audit$cull$min_year, audit$cull$max_year
))
cat(sprintf(
  "Adjusted hindcast: %d reef-years, %d reefs, %d-%d.\n",
  audit$outbreak_hindcast$rows, audit$outbreak_hindcast$reefs,
  audit$outbreak_hindcast$min_year, audit$outbreak_hindcast$max_year
))
cat(sprintf(
  "Adjusted forecast: %d reefs for %d; cutoff %s.\n",
  audit$outbreak_forecast$reefs, audit$outbreak_forecast$year,
  paste(audit$outbreak_forecast$cutoff, collapse = ", ")
))
cat(sprintf(
  "Exact reef-name match: %d/%d (%.1f%%).\n",
  audit$join$matched_cull_reefs, audit$cull$reefs,
  100 * audit$join$match_rate
))
if (length(unmatched) > 0L) {
  cat("Unmatched reefs:", paste(unmatched, collapse = "; "), "\n")
}
cat(
  "Forecast probabilities are audit/operational fields and are excluded ",
  "from hindcast model features.\n",
  sep = ""
)
