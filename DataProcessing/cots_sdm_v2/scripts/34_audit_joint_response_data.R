args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v2_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v2_dir, "R", "bootstrap.R"))
ctx <- initialise_v2(script_file)
cfg <- ctx$config

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("readxl is required for the joint-response data audit.")
}
source_file <- cfg$paths$cull_workbook
if (!file.exists(source_file)) stop("Combined Cull-Manta-RHIS workbook missing.")

read_sheet <- function(sheet) {
  suppressMessages(readxl::read_excel(source_file, sheet = sheet))
}
cull <- read_sheet("Cull")
manta <- read_sheet("Manta")
rhis <- read_sheet("RHIS")

as_day <- function(x) as.Date(x)
valid_text <- function(x) !is.na(x) & nzchar(trimws(as.character(x)))
nonzero_rows <- function(data, fields) {
  value <- data[, intersect(fields, names(data)), drop = FALSE]
  if (ncol(value) == 0L) return(NA_integer_)
  numeric_value <- lapply(value, function(x) suppressWarnings(as.numeric(x)))
  numeric_value <- as.data.frame(numeric_value)
  sum(rowSums(numeric_value, na.rm = TRUE) > 0)
}
date_range <- function(x) {
  x <- as_day(x)
  x <- x[!is.na(x)]
  if (length(x) == 0L) return(c(NA_character_, NA_character_))
  as.character(range(x))
}

sheet_summary <- function(
  data, sheet, date_field, reef_field, site_fields,
  response_fields, effort_fields, ecological_fields
) {
  dates <- date_range(data[[date_field]])
  site_fields_present <- intersect(site_fields, names(data))
  response_fields_present <- intersect(response_fields, names(data))
  response_nonzero <- nonzero_rows(data, response_fields_present)
  named_sites <- unique(unlist(lapply(
    data[site_fields_present], function(x) as.character(x[valid_text(x)])
  )))
  tibble::tibble(
    sheet = sheet,
    rows = nrow(data),
    date_min = dates[[1]],
    date_max = dates[[2]],
    reefs = dplyr::n_distinct(data[[reef_field]], na.rm = TRUE),
    named_sites = length(named_sites),
    response_fields = paste(
      response_fields_present, collapse = "; "
    ),
    response_nonzero_rows = response_nonzero,
    effort_detection_fields = paste(
      intersect(effort_fields, names(data)), collapse = "; "
    ),
    ecological_coresponse_fields = paste(
      intersect(ecological_fields, names(data)), collapse = "; "
    )
  )
}

cull_response <- paste0("Cohort", 1:4)
manta_response <- c(
  "CrownOfThornsStarfishCount", "ScarsCount"
)
rhis_cots <- c(
  "PredatorCrownOfThornsStarfishJuvenileIndividualCount",
  "PredatorCrownOfThornsStarfishAdultIndividualCount"
)
rhis_scars <- grep(
  "CrownOfThornsStarfishScarCount$", names(rhis), value = TRUE
)
rhis_coral <- c(
  grep("^Coral.*LiveProportion$", names(rhis), value = TRUE),
  "BenthosLiveCoralPercentage", "BenthosRecentlyDeadCoralPercentage"
)

availability <- dplyr::bind_rows(
  sheet_summary(
    cull, "Cull", "SurveyDate", "ReefLabel", "CullSiteName",
    cull_response,
    c("Bottomtime", "Depth", "CullDiveTypeName", "VesselName_ID"),
    character()
  ),
  sheet_summary(
    manta, "Manta", "SurveyTime", "ReefLabel",
    c("CullSiteName", "SecondaryCullSiteName"),
    manta_response,
    c("TowDistance", "FeedingScarCountRangeCode", "VesselName_ID"),
    c(
      "HardCoralProportionRangeCode", "SoftCoralProportionRangeCode",
      "RecentlyDeadCoralProportionRangeCode"
    )
  ),
  sheet_summary(
    rhis, "RHIS", "SurveyTime", "ReefLabel",
    c("CullSite", "DiveSite"),
    c(rhis_cots, rhis_scars),
    c(
      "VisibilityDistanceRange", "Depth", "SwimType", "SurveyHabitat",
      "VesselName_ID"
    ),
    rhis_coral
  )
)

repeat_summary <- function(data, method, date_field, reef_field, site_field) {
  x <- data |>
    dplyr::transmute(
      reef = as.character(.data[[reef_field]]),
      site = as.character(.data[[site_field]]),
      day = as_day(.data[[date_field]])
    ) |>
    dplyr::filter(valid_text(.data$reef), valid_text(.data$site), !is.na(.data$day)) |>
    dplyr::distinct(.data$reef, .data$site, .data$day) |>
    dplyr::count(.data$reef, .data$site, name = "visit_days")
  tibble::tibble(
    comparison = paste0(method, " named site"),
    units = nrow(x),
    units_with_repeat_days = sum(x$visit_days > 1L),
    maximum_visit_days = max(x$visit_days, 0L),
    shared_reef_days = NA_integer_
  )
}
reef_days <- function(data, date_field) {
  data |>
    dplyr::transmute(
      reef = as.character(.data$ReefLabel),
      day = as_day(.data[[date_field]])
    ) |>
    dplyr::filter(valid_text(.data$reef), !is.na(.data$day)) |>
    dplyr::distinct(.data$reef, .data$day)
}
shared_days <- function(left, right, label) {
  tibble::tibble(
    comparison = label,
    units = NA_integer_,
    units_with_repeat_days = NA_integer_,
    maximum_visit_days = NA_integer_,
    shared_reef_days = nrow(dplyr::inner_join(
      left, right, by = c("reef", "day")
    ))
  )
}

cull_days <- reef_days(cull, "SurveyDate")
manta_days <- reef_days(manta, "SurveyTime")
rhis_days <- reef_days(rhis, "SurveyTime")
alignment <- dplyr::bind_rows(
  repeat_summary(cull, "Cull", "SurveyDate", "ReefLabel", "CullSiteName"),
  repeat_summary(manta, "Manta", "SurveyTime", "ReefLabel", "CullSiteName"),
  repeat_summary(rhis, "RHIS", "SurveyTime", "ReefLabel", "CullSite"),
  shared_days(cull_days, manta_days, "Cull-Manta same reef/day"),
  shared_days(cull_days, rhis_days, "Cull-RHIS same reef/day"),
  shared_days(manta_days, rhis_days, "Manta-RHIS same reef/day")
)

utils::write.csv(
  availability,
  file.path(cfg$paths$outputs_dir, "joint_response_data_inventory.csv"),
  row.names = FALSE
)
utils::write.csv(
  alignment,
  file.path(cfg$paths$outputs_dir, "joint_response_alignment_audit.csv"),
  row.names = FALSE
)
message("Joint-response workbook audit complete.")
print(availability)
print(alignment)
