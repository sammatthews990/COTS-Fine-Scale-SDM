ensure_output_directories <- function(cfg) {
  dirs <- unlist(cfg$paths[c(
    "artifacts_dir", "models_dir", "outputs_dir", "logs_dir"
  )], use.names = FALSE)
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
}

assert_required_columns <- function(data, required, data_name) {
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop(
      data_name, " is missing required columns: ",
      paste(missing, collapse = ", ")
    )
  }
  invisible(TRUE)
}

read_cull_observations <- function(path) {
  if (!file.exists(path)) stop("Cull workbook not found: ", path)

  required <- c(
    "CrownOfThornsStarfishCullDiveId", "SurveyDate", "VoyageTitle",
    "ReefName", "CullSiteName", "Latitude", "Longitude", "Bottomtime",
    "Cohort1", "Cohort2", "Cohort3", "Cohort4"
  )

  raw <- suppressWarnings(readxl::read_excel(path, sheet = "Cull"))
  assert_required_columns(raw, required, "Cull sheet")

  raw |>
    dplyr::mutate(
      source_record_id = dplyr::row_number(),
      survey_date = as.Date(substr(as.character(.data$SurveyDate), 1, 10)),
      year = lubridate::year(.data$survey_date),
      month = lubridate::month(.data$survey_date),
      financial_year_start = dplyr::if_else(
        .data$month >= 7L, .data$year, .data$year - 1L
      ),
      financial_year = sprintf(
        "FY%02d/%02d",
        .data$financial_year_start %% 100,
        (.data$financial_year_start + 1L) %% 100
      ),
      cohort_non_missing = rowSums(!is.na(dplyr::pick(
        Cohort1, Cohort2, Cohort3, Cohort4
      ))),
      total_cots = rowSums(dplyr::pick(
        Cohort1, Cohort2, Cohort3, Cohort4
      ), na.rm = TRUE),
      total_cots = dplyr::if_else(
        .data$cohort_non_missing == 0L, NA_real_, .data$total_cots
      ),
      bottom_time = as.numeric(.data$Bottomtime),
      dive_cpue = .data$total_cots / .data$bottom_time,
      reef_name = trimws(as.character(.data$ReefName)),
      site_name = trimws(as.character(.data$CullSiteName)),
      longitude = as.numeric(.data$Longitude),
      latitude = as.numeric(.data$Latitude),
      site_id = paste(
        .data$reef_name,
        .data$site_name,
        sprintf("%.6f", .data$longitude),
        sprintf("%.6f", .data$latitude),
        sep = "||"
      )
    ) |>
    dplyr::filter(
      !is.na(.data$survey_date),
      !is.na(.data$longitude), !is.na(.data$latitude),
      is.finite(.data$longitude), is.finite(.data$latitude),
      !is.na(.data$bottom_time), .data$bottom_time > 0,
      !is.na(.data$total_cots)
    )
}

read_outbreak_probabilities <- function(paths, product_type) {
  existing <- paths[file.exists(paths)]
  if (length(existing) == 0L) {
    stop("No configured reef outbreak probability files were found.")
  }

  out <- dplyr::bind_rows(lapply(existing, function(path) {
    x <- utils::read.csv(path, stringsAsFactors = FALSE)
    assert_required_columns(
      x, c("year", "reefName", "outbrProb", "x", "y"), basename(path)
    )
    probability_cutoff <- if ("probCutoff" %in% names(x)) {
      as.numeric(x$probCutoff)
    } else {
      rep(NA_real_, nrow(x))
    }
    dplyr::transmute(
      x,
      year = as.integer(.data$year),
      reef_name = trimws(as.character(.data$reefName)),
      reef_outbreak_probability = as.numeric(.data$outbrProb),
      reef_longitude = as.numeric(.data$x),
      reef_latitude = as.numeric(.data$y),
      outbreak_probability_cutoff = .env$probability_cutoff,
      outbreak_product = product_type,
      source_file = basename(path)
    )
  }))

  duplicate_keys <- out |>
    dplyr::count(.data$reef_name, .data$year) |>
    dplyr::filter(.data$n > 1L)
  if (nrow(duplicate_keys) > 0L) {
    stop("Outbreak probabilities contain duplicate reef-year keys.")
  }
  if (any(!is.finite(out$reef_outbreak_probability)) ||
      any(out$reef_outbreak_probability < 0 | out$reef_outbreak_probability > 1)) {
    stop("Outbreak probabilities must be finite values between zero and one.")
  }

  dplyr::arrange(out, .data$reef_name, .data$year)
}

read_manta_observations <- function(path, sheet = NULL) {
  if (!file.exists(path)) stop("Manta input not found: ", path)
  raw <- if (is.null(sheet)) {
    suppressWarnings(readr::read_csv(path, show_col_types = FALSE))
  } else {
    suppressWarnings(readxl::read_excel(path, sheet = sheet))
  }
  required <- c(
    "SurveyTime", "ReefName", "StartLatitude", "StartLongitude",
    "EndLatitude", "EndLongitude", "TowDistance",
    "CrownOfThornsStarfishCount", "FeedingScarCountRangeCode",
    "ScarsCount", "HardCoralProportionRangeCode", "Source"
  )
  assert_required_columns(raw, required, basename(path))

  raw |>
    dplyr::transmute(
      manta_record_id = paste0(basename(path), "||", dplyr::row_number()),
      survey_date = as.Date(substr(as.character(.data$SurveyTime), 1, 10)),
      year = lubridate::year(.data$survey_date),
      reef_name = trimws(as.character(.data$ReefName)),
      start_latitude = as.numeric(.data$StartLatitude),
      start_longitude = as.numeric(.data$StartLongitude),
      end_latitude = as.numeric(.data$EndLatitude),
      end_longitude = as.numeric(.data$EndLongitude),
      tow_distance_m = as.numeric(.data$TowDistance),
      cots_count = as.numeric(.data$CrownOfThornsStarfishCount),
      scar_code = tolower(trimws(as.character(
        .data$FeedingScarCountRangeCode
      ))),
      scar_count = as.numeric(.data$ScarsCount),
      hard_coral_code = trimws(as.character(
        .data$HardCoralProportionRangeCode
      )),
      source = as.character(.data$Source),
      source_file = basename(path)
    ) |>
    dplyr::mutate(
      informative = !is.na(.data$cots_count) &
        .data$scar_code %in% c("a", "p", "c"),
      positive_evidence = dplyr::coalesce(.data$cots_count > 0, FALSE) |
        .data$scar_code %in% c("p", "c") |
        dplyr::coalesce(.data$scar_count > 0, FALSE),
      strict_negative = .data$informative & .data$cots_count == 0 &
        .data$scar_code == "a" &
        (is.na(.data$scar_count) | .data$scar_count == 0),
      hard_coral_observed = !is.na(.data$hard_coral_code) &
        .data$hard_coral_code != "" & .data$hard_coral_code != "0",
      valid_start_coordinates = is.finite(.data$start_longitude) &
        is.finite(.data$start_latitude),
      valid_end_coordinates = is.finite(.data$end_longitude) &
        is.finite(.data$end_latitude),
      longitude = dplyr::if_else(
        .data$valid_start_coordinates & .data$valid_end_coordinates,
        (.data$start_longitude + .data$end_longitude) / 2,
        .data$start_longitude
      ),
      latitude = dplyr::if_else(
        .data$valid_start_coordinates & .data$valid_end_coordinates,
        (.data$start_latitude + .data$end_latitude) / 2,
        .data$start_latitude
      ),
      valid_coordinates = is.finite(.data$longitude) &
        is.finite(.data$latitude)
    ) |>
    dplyr::filter(
      !is.na(.data$survey_date),
      !is.na(.data$reef_name), .data$reef_name != ""
    )
}
