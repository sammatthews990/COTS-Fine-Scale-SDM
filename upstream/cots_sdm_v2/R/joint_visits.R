normalise_key <- function(x) {
  x <- toupper(gsub("[^A-Z0-9]+", "", trimws(as.character(x))))
  vapply(x, function(value) {
    token <- regmatches(value, gregexpr("[A-Z]+|[0-9]+", value))[[1]]
    token <- vapply(token, function(part) {
      if (grepl("^[0-9]+$", part)) as.character(as.integer(part)) else part
    }, character(1))
    paste(token, collapse = "")
  }, character(1))
}

normalise_reef_key <- function(x) {
  tolower(gsub("[^a-z0-9]+", "", trimws(as.character(x))))
}

parse_salad_coordinate <- function(x, coordinate = c("latitude", "longitude")) {
  coordinate <- match.arg(coordinate)
  text <- trimws(as.character(x))
  decimal <- suppressWarnings(as.numeric(text))
  dm_match <- regexec(
    "^(-?[0-9]+)[oO]([0-9]+(?:[.][0-9]+)?)$", text, perl = TRUE
  )
  dm_parts <- regmatches(text, dm_match)
  for (i in seq_along(dm_parts)) {
    if (length(dm_parts[[i]]) == 3L) {
      degree <- as.numeric(dm_parts[[i]][[2]])
      minute <- as.numeric(dm_parts[[i]][[3]])
      if (is.finite(minute) && minute < 60) {
        decimal[[i]] <- sign(ifelse(degree == 0, 1, degree)) *
          (abs(degree) + minute / 60)
      }
    }
  }
  if (coordinate == "latitude") {
    decimal[is.finite(decimal)] <- -abs(decimal[is.finite(decimal)])
    decimal[!dplyr::between(decimal, -30, -5)] <- NA_real_
  } else {
    decimal[!dplyr::between(decimal, 140, 160)] <- NA_real_
  }
  decimal
}

read_rhis_observations <- function(path) {
  raw <- suppressWarnings(readxl::read_excel(path, sheet = "RHIS"))
  required <- c(
    "SurveyId", "SurveyTime", "ReefName", "Latitude", "Longitude",
    "PredatorCrownOfThornsStarfishJuvenileIndividualCount",
    "PredatorCrownOfThornsStarfishAdultIndividualCount",
    "PredatorBranchingCrownOfThornsStarfishScarCount",
    "PredatorPlateTableCrownOfThornsStarfishScarCount",
    "CoralBranchingLiveProportion", "CoralPlateTableLiveProportion",
    "CoralBranchingRecentlyDeadProportion",
    "CoralPlateTableRecentlyDeadProportion"
  )
  assert_required_columns(raw, required, "RHIS sheet")
  reef_label <- if ("ReefLabel" %in% names(raw)) {
    as.character(raw$ReefLabel)
  } else {
    rep(NA_character_, nrow(raw))
  }
  dplyr::transmute(
    raw,
    visit_id = paste0("rhis||", .data$SurveyId, "||", dplyr::row_number()),
    method = "rhis",
    survey_date = as.Date(substr(as.character(.data$SurveyTime), 1, 10)),
    year = lubridate::year(.data$survey_date),
    reef_label = trimws(.env$reef_label),
    reef_name = trimws(as.character(.data$ReefName)),
    reef_key = normalise_reef_key(.data$reef_name),
    longitude = as.numeric(.data$Longitude),
    latitude = as.numeric(.data$Latitude),
    footprint_length_m = NA_real_,
    adult_count = as.numeric(
      .data$PredatorCrownOfThornsStarfishAdultIndividualCount
    ),
    juvenile_count = as.numeric(
      .data$PredatorCrownOfThornsStarfishJuvenileIndividualCount
    ),
    branching_scar_count = as.numeric(
      .data$PredatorBranchingCrownOfThornsStarfishScarCount
    ),
    table_scar_count = as.numeric(
      .data$PredatorPlateTableCrownOfThornsStarfishScarCount
    ),
    branching_live_pct = as.numeric(.data$CoralBranchingLiveProportion),
    table_live_pct = as.numeric(.data$CoralPlateTableLiveProportion),
    branching_recent_dead_pct = as.numeric(
      .data$CoralBranchingRecentlyDeadProportion
    ),
    table_recent_dead_pct = as.numeric(
      .data$CoralPlateTableRecentlyDeadProportion
    ),
    visibility = as.character(.data$VisibilityDistanceRange),
    depth_m = as.numeric(.data$Depth),
    survey_habitat = as.character(.data$SurveyHabitat)
  ) |>
    dplyr::mutate(
      rhis_cots_protocol_observed = rowSums(!is.na(dplyr::pick(
        "adult_count", "juvenile_count", "branching_scar_count",
        "table_scar_count"
      ))) > 0L
    ) |>
    dplyr::filter(
      !is.na(.data$survey_date),
      is.finite(.data$longitude), is.finite(.data$latitude),
      dplyr::between(.data$longitude, 140, 160),
      dplyr::between(.data$latitude, -30, -5)
    )
}

read_salad_tracks <- function(path) {
  raw <- suppressWarnings(readxl::read_excel(path, sheet = "Tracks"))
  required <- c(
    "Year", "Date", "Track", "Diver", "Region", "Reef", "Location",
    "Start_Time", "End_Time", "Time (Min)", "No. COTS", "No. scars",
    "Distance (m)", "AREA(5MWIDE)", "Density_inclScars"
  )
  assert_required_columns(raw, required, "SALAD Tracks sheet")
  coordinate_names <- grep("^(Start|End)_Coord", names(raw), value = TRUE)
  if (length(coordinate_names) != 4L) {
    stop("SALAD Tracks must contain four start/end coordinate columns.")
  }
  number <- function(x) suppressWarnings(as.numeric(as.character(x)))
  start_latitude <- parse_salad_coordinate(
    raw[[coordinate_names[[1]]]], "latitude"
  )
  start_longitude <- parse_salad_coordinate(
    raw[[coordinate_names[[2]]]], "longitude"
  )
  end_latitude <- parse_salad_coordinate(
    raw[[coordinate_names[[3]]]], "latitude"
  )
  end_longitude <- parse_salad_coordinate(
    raw[[coordinate_names[[4]]]], "longitude"
  )
  loop_to_start <- grepl(
    "back|loop", as.character(raw[[coordinate_names[[3]]]]), ignore.case = TRUE
  ) | grepl(
    "back|loop", as.character(raw[[coordinate_names[[4]]]]), ignore.case = TRUE
  )
  end_latitude[loop_to_start] <- start_latitude[loop_to_start]
  end_longitude[loop_to_start] <- start_longitude[loop_to_start]
  start_latitude[is.na(start_latitude)] <- end_latitude[is.na(start_latitude)]
  start_longitude[is.na(start_longitude)] <- end_longitude[is.na(start_longitude)]
  end_latitude[is.na(end_latitude)] <- start_latitude[is.na(end_latitude)]
  end_longitude[is.na(end_longitude)] <- start_longitude[is.na(end_longitude)]
  dplyr::transmute(
    raw,
    visit_id = paste0("salad||", .data$Count),
    method = "salad",
    survey_date = as.Date(.data$Date),
    year = as.integer(.data$Year),
    track = as.character(.data$Track),
    track_key = normalise_key(.data$Track),
    diver = as.character(.data$Diver),
    reef_label = NA_character_,
    reef_name = trimws(as.character(.data$Reef)),
    reef_key = normalise_reef_key(.data$Reef),
    site_name = trimws(as.character(.data$Location)),
    start_time_fraction = number(.data$Start_Time),
    end_time_fraction = number(.data$End_Time),
    start_latitude = .env$start_latitude,
    start_longitude = .env$start_longitude,
    end_latitude = .env$end_latitude,
    end_longitude = .env$end_longitude,
    longitude = (.data$start_longitude + .data$end_longitude) / 2,
    latitude = (.data$start_latitude + .data$end_latitude) / 2,
    footprint_length_m = number(.data$`Distance (m)`),
    footprint_area_m2 = number(.data$`AREA(5MWIDE)`),
    survey_minutes = number(.data$`Time (Min)`),
    direct_cots_count = number(.data$`No. COTS`),
    scar_count = number(.data$`No. scars`),
    density_incl_scars_ha = number(.data$Density_inclScars),
    detectability = number(.data$Detectability)
  ) |>
    dplyr::mutate(
      coordinate_method = "track_endpoint_midpoint",
      valid_coordinates = is.finite(.data$longitude) &
        is.finite(.data$latitude) &
        dplyr::between(.data$longitude, 140, 160) &
        dplyr::between(.data$latitude, -30, -5)
    ) |>
    dplyr::filter(!is.na(.data$survey_date), .data$valid_coordinates)
}

interpolate_salad_records <- function(path, tracks) {
  raw <- suppressWarnings(readxl::read_excel(path, sheet = "COTS recorded"))
  assert_required_columns(
    raw, c("#", "Date", "Track", "Time", "CoTS/ Scar"),
    "SALAD COTS recorded sheet"
  )
  events <- dplyr::transmute(
    raw,
    event_id = paste0("salad_event||", dplyr::row_number()),
    survey_date = as.Date(.data$Date),
    track_key = normalise_key(.data$Track),
    event_time_fraction = suppressWarnings(as.numeric(as.character(.data$Time))),
    evidence = tolower(trimws(as.character(.data$`CoTS/ Scar`))),
    size_mm = as.numeric(.data$`Size(mm)`),
    depth_m = as.numeric(.data$`Depth(m)`)
  )
  track_lookup <- tracks |>
    dplyr::arrange(.data$survey_date, .data$track_key, .data$visit_id) |>
    dplyr::distinct(.data$survey_date, .data$track_key, .keep_all = TRUE) |>
    dplyr::select(
      "survey_date", "track_key", "visit_id",
      "start_time_fraction", "end_time_fraction",
      "start_longitude", "start_latitude", "end_longitude", "end_latitude"
    )
  events |>
    dplyr::left_join(track_lookup, by = c("survey_date", "track_key")) |>
    dplyr::mutate(
      along_track_fraction = (
        .data$event_time_fraction - .data$start_time_fraction
      ) / (.data$end_time_fraction - .data$start_time_fraction),
      interpolation_valid = is.finite(.data$along_track_fraction) &
        dplyr::between(.data$along_track_fraction, 0, 1),
      longitude = dplyr::if_else(
        .data$interpolation_valid,
        .data$start_longitude + .data$along_track_fraction *
          (.data$end_longitude - .data$start_longitude),
        NA_real_
      ),
      latitude = dplyr::if_else(
        .data$interpolation_valid,
        .data$start_latitude + .data$along_track_fraction *
          (.data$end_latitude - .data$start_latitude),
        NA_real_
      )
    )
}

assign_visit_cells <- function(data, template, grid_m = 200) {
  xy <- sf::sf_project(
    from = "EPSG:4326", to = terra::crs(template),
    pts = as.matrix(data[, c("longitude", "latitude")])
  )
  data |>
    dplyr::mutate(
      projected_x = xy[, 1],
      projected_y = xy[, 2],
      x_bin = floor(.data$projected_x / grid_m),
      y_bin = floor(.data$projected_y / grid_m),
      cell_id = paste(.data$x_bin, .data$y_bin, sep = "||")
    )
}

cluster_visit_dates <- function(visits, maximum_gap_days = 14L) {
  visits |>
    dplyr::arrange(.data$cell_id, .data$survey_date, .data$method, .data$visit_id) |>
    dplyr::group_by(.data$cell_id) |>
    dplyr::mutate(
      date_gap = as.integer(.data$survey_date - dplyr::lag(.data$survey_date)),
      new_cluster = dplyr::row_number() == 1L |
        dplyr::coalesce(.data$date_gap > maximum_gap_days, TRUE),
      visit_cluster = cumsum(.data$new_cluster),
      visit_cluster_id = paste(.data$cell_id, .data$visit_cluster, sep = "||")
    ) |>
    dplyr::ungroup()
}

annotate_local_cull_history <- function(visits) {
  cull_history <- visits |>
    dplyr::filter(.data$method == "cull") |>
    dplyr::group_by(.data$cell_id) |>
    dplyr::summarise(
      local_first_cull_date = min(.data$survey_date),
      .groups = "drop"
    )
  visits |>
    dplyr::left_join(cull_history, by = "cell_id") |>
    dplyr::mutate(
      days_from_local_first_cull = as.integer(
        .data$survey_date - .data$local_first_cull_date
      ),
      cull_timing = dplyr::case_when(
        is.na(.data$local_first_cull_date) ~ "no_local_cull_record",
        .data$days_from_local_first_cull < 0L ~ "strictly_pre_cull",
        .data$days_from_local_first_cull == 0L ~ "same_day_ambiguous",
        TRUE ~ "post_cull"
      ),
      ecological_eligible = .data$cull_timing %in% c(
        "no_local_cull_record", "strictly_pre_cull"
      ) | (
        .data$method == "cull" & .data$cull_timing == "same_day_ambiguous"
      )
    )
}
