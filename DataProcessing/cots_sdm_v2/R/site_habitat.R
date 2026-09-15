geomorphic_class_levels <- function() {
  c(
    "Deep", "Shallow Lagoon", "Deep Lagoon", "Inner Reef Flat",
    "Outer Reef Flat", "Reef Crest", "Sheltered Reef Slope",
    "Reef Slope", "Plateau", "Back Reef Slope", "Patch Reef",
    "Unknown"
  )
}

geomorphic_feature_names <- function() {
  paste0(
    "geom_",
    gsub("[^a-z0-9]+", "_", tolower(geomorphic_class_levels()))
  )
}

summarise_cull_site_history <- function(site_year) {
  site_year |>
    dplyr::group_by(
      .data$site_id, .data$reef_name, .data$site_name,
      .data$longitude, .data$latitude
    ) |>
    dplyr::summarise(
      high_density = as.integer(any(.data$high_density == 1L)),
      n_surveyed_years = dplyr::n_distinct(.data$year),
      n_high_density_years = sum(.data$high_density),
      high_density_year_fraction = mean(.data$high_density),
      first_year = min(.data$year),
      last_year = max(.data$year),
      total_bottom_time = sum(.data$bottom_time),
      total_cots = sum(.data$total_cots),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      evidence_source = "cull_site_history",
      observation_weight = 1
    )
}

project_manta_midpoints <- function(manta, template) {
  eligible <- manta |>
    dplyr::filter(
      .data$valid_coordinates,
      .data$informative,
      is.na(.data$tow_distance_m) | .data$tow_distance_m >= 100
    )
  points <- sf::st_as_sf(
    eligible,
    coords = c("longitude", "latitude"),
    crs = 4326,
    remove = FALSE
  ) |>
    sf::st_transform(terra::crs(template))
  xy <- sf::st_coordinates(points)
  eligible$x <- xy[, 1]
  eligible$y <- xy[, 2]
  eligible |>
    dplyr::filter(is.finite(.data$x), is.finite(.data$y))
}

summarise_manta_sites <- function(
    manta, template, grid_m = 200, minimum_tows = 5L,
    minimum_years = 2L) {
  projected <- project_manta_midpoints(manta, template) |>
    dplyr::mutate(
      x_bin = floor((.data$x - terra::xmin(template)) / grid_m),
      y_bin = floor((.data$y - terra::ymin(template)) / grid_m)
    )

  grouped <- projected |>
    dplyr::group_by(.data$reef_name, .data$x_bin, .data$y_bin) |>
    dplyr::summarise(
      n_tows = dplyr::n(),
      n_survey_dates = dplyr::n_distinct(.data$survey_date),
      n_surveyed_years = dplyr::n_distinct(.data$year),
      first_year = min(.data$year),
      last_year = max(.data$year),
      any_positive_evidence = any(.data$positive_evidence),
      all_strict_negative = all(.data$strict_negative),
      all_hard_coral_observed = all(.data$hard_coral_observed),
      total_tow_distance_m = sum(.data$tow_distance_m, na.rm = TRUE),
      x = mean(.data$x),
      y = mean(.data$y),
      .groups = "drop"
    )

  candidates <- grouped |>
    dplyr::filter(
      !.data$any_positive_evidence,
      .data$all_strict_negative,
      .data$all_hard_coral_observed,
      .data$n_tows >= minimum_tows,
      .data$n_surveyed_years >= minimum_years
    )

  list(
    candidates = candidates,
    positive_sites = dplyr::filter(
      grouped, .data$any_positive_evidence
    ),
    all_sites = grouped
  )
}

has_nearby_point <- function(query, reference, distance_m, crs) {
  if (nrow(reference) == 0L) return(rep(FALSE, nrow(query)))
  query_sf <- sf::st_as_sf(
    query, coords = c("x", "y"), crs = crs, remove = FALSE
  )
  reference_sf <- sf::st_as_sf(
    reference, coords = c("x", "y"), crs = crs, remove = FALSE
  )
  lengths(sf::st_is_within_distance(
    query_sf, reference_sf, dist = distance_m
  )) > 0L
}

exclude_manta_local_conflicts <- function(
    manta_summary, cull_sites, template, distance_m = 250) {
  cull_sf <- sf::st_as_sf(
    cull_sites,
    coords = c("longitude", "latitude"),
    crs = 4326,
    remove = FALSE
  ) |>
    sf::st_transform(terra::crs(template))
  cull_xy <- sf::st_coordinates(cull_sf)
  cull_projected <- cull_sites
  cull_projected$x <- cull_xy[, 1]
  cull_projected$y <- cull_xy[, 2]
  positive_cull <- dplyr::filter(
    cull_projected, .data$high_density == 1L
  )

  candidates <- manta_summary$candidates |>
    dplyr::mutate(
      conflict_positive_cull = has_nearby_point(
        dplyr::pick(dplyr::everything()),
        positive_cull,
        distance_m,
        terra::crs(template)
      ),
      conflict_positive_manta = has_nearby_point(
        dplyr::pick(dplyr::everything()),
        manta_summary$positive_sites,
        distance_m,
        terra::crs(template)
      ),
      local_conflict =
        .data$conflict_positive_cull |
        .data$conflict_positive_manta
    )

  list(
    retained = dplyr::filter(candidates, !.data$local_conflict),
    audit = candidates
  )
}

finalise_manta_zero_sites <- function(sites, template, grid_m) {
  points <- sf::st_as_sf(
    sites,
    coords = c("x", "y"),
    crs = terra::crs(template),
    remove = FALSE
  ) |>
    sf::st_transform(4326)
  xy <- sf::st_coordinates(points)

  sites |>
    dplyr::mutate(
      longitude = xy[, 1],
      latitude = xy[, 2],
      site_id = paste(
        "manta", paste0(grid_m, "m"), .data$reef_name,
        .data$x_bin, .data$y_bin, sep = "||"
      ),
      site_name = paste0(
        "manta_", grid_m, "m_", .data$x_bin, "_", .data$y_bin
      ),
      high_density = 0L,
      n_high_density_years = 0L,
      high_density_year_fraction = 0,
      total_bottom_time = NA_real_,
      total_cots = 0,
      evidence_source = "manta_repeated_site_zero",
      observation_weight = 1
    )
}

extract_site_habitat_predictors <- function(
    sites, clean_path, rg_path, geomorphic_path) {
  clean <- terra::rast(clean_path)
  rg_stack <- terra::rast(rg_path)
  names(rg_stack)[13:17] <- c(
    "RG_bathy", "RG_slope", "RG_turbid",
    "RG_waves_Hs", "RG_waves_Tp"
  )
  rg_names <- names(rg_stack)[13:17]
  geomorphic <- terra::rast(geomorphic_path)

  points_wgs <- terra::vect(
    sites,
    geom = c("longitude", "latitude"),
    crs = "EPSG:4326"
  )
  points_native <- terra::project(points_wgs, terra::crs(clean))
  clean_values <- terra::extract(clean, points_native, ID = FALSE)
  rg_values <- terra::extract(
    rg_stack[[rg_names]], points_native, ID = FALSE
  )
  geomorphic_values <- terra::extract(
    geomorphic, points_wgs, ID = FALSE
  )

  geomorphic_class <- as.character(geomorphic_values[[1]])
  geomorphic_class[
    is.na(geomorphic_class) |
      !geomorphic_class %in% geomorphic_class_levels()
  ] <- "Unknown"

  result <- dplyr::bind_cols(
    sites,
    clean_values,
    rg_values,
    tibble::tibble(geomorphic_class = geomorphic_class)
  )
  for (i in seq_along(geomorphic_class_levels())) {
    result[[geomorphic_feature_names()[i]]] <- as.integer(
      result$geomorphic_class == geomorphic_class_levels()[i]
    )
  }
  result |>
    dplyr::mutate(
      clean_complete = as.integer(stats::complete.cases(
        dplyr::pick(dplyr::all_of(names(clean)))
      )),
      rg_complete = as.integer(stats::complete.cases(
        dplyr::pick(dplyr::all_of(rg_names))
      ))
    )
}

add_site_validation_group <- function(sites, template, grid_m = 200) {
  points <- sf::st_as_sf(
    sites,
    coords = c("longitude", "latitude"),
    crs = 4326,
    remove = FALSE
  ) |>
    sf::st_transform(terra::crs(template))
  xy <- sf::st_coordinates(points)
  sites |>
    dplyr::mutate(
      validation_x_bin = floor(
        (xy[, 1] - terra::xmin(template)) / grid_m
      ),
      validation_y_bin = floor(
        (xy[, 2] - terra::ymin(template)) / grid_m
      ),
      validation_group_id = paste(
        .data$reef_name,
        .data$validation_x_bin,
        .data$validation_y_bin,
        sep = "||"
      )
    )
}

assign_within_reef_site_folds <- function(
    sites, k = 5L, seed = 1L) {
  groups <- sites |>
    dplyr::distinct(.data$reef_name, .data$validation_group_id)
  reefs <- split(seq_len(nrow(groups)), groups$reef_name)
  fold <- rep(NA_integer_, nrow(groups))
  set.seed(seed)
  for (index in reefs) {
    shuffled <- if (length(index) > 1L) {
      sample(index, length(index), replace = FALSE)
    } else {
      index
    }
    fold[shuffled] <- rep(seq_len(k), length.out = length(index))
  }
  groups$site_fold <- fold
  sites |>
    dplyr::left_join(
      dplyr::select(
        groups,
        dplyr::all_of(c("validation_group_id", "site_fold"))
      ),
      by = "validation_group_id"
    )
}

site_habitat_feature_sets <- function(clean_names) {
  rg_names <- c(
    "RG_bathy", "RG_slope", "RG_turbid",
    "RG_waves_Hs", "RG_waves_Tp"
  )
  list(
    clean = clean_names,
    clean_geomorphic = c(clean_names, geomorphic_feature_names()),
    rg = c(clean_names, rg_names),
    rg_geomorphic = c(
      clean_names, rg_names, geomorphic_feature_names()
    )
  )
}

run_site_habitat_benchmark <- function(
    cull_sites, manta_sites, feature_sets, model_cfg,
    maximum_manta_weight_ratio = 1,
    augmentation_modes = c(FALSE, TRUE)) {
  folds <- sort(unique(cull_sites$site_fold))
  predictions <- list()
  fitted_models <- list()
  training_audit <- list()
  counter <- 0L

  for (family in names(feature_sets)) {
    features <- feature_sets[[family]]
    for (use_manta in augmentation_modes) {
      model_name <- paste0(
        family,
        if (use_manta) "_plus_manta" else "_cull_only"
      )
      for (fold in folds) {
        counter <- counter + 1L
        cull_train <- dplyr::filter(
          cull_sites, .data$site_fold != fold
        )
        test <- dplyr::filter(
          cull_sites, .data$site_fold == fold
        )
        cull_train$observation_weight <- 1
        train <- cull_train
        manta_weight <- NA_real_

        if (use_manta) {
          manta_train <- dplyr::filter(
            manta_sites, .data$site_fold != fold
          )
          manta_weight <- min(
            1,
            maximum_manta_weight_ratio *
              sum(cull_train$high_density == 0L) /
              max(1, nrow(manta_train))
          )
          manta_train$observation_weight <- manta_weight
          train <- dplyr::bind_rows(cull_train, manta_train)
        }

        model <- fit_xgb_binary(
          train,
          features,
          model_cfg,
          seed_offset = 20000L + counter,
          weight_col = "observation_weight"
        )
        probability <- predict_xgb_binary(model, test, features)
        predictions[[counter]] <- tibble::tibble(
          model = model_name,
          feature_family = family,
          manta_augmented = use_manta,
          validation = "within_reef_site_cv",
          fold = fold,
          site_id = test$site_id,
          reef_name = test$reef_name,
          truth = test$high_density,
          probability = probability,
          n_surveyed_years = test$n_surveyed_years,
          n_high_density_years = test$n_high_density_years,
          longitude = test$longitude,
          latitude = test$latitude
        )
        fitted_models[[paste(model_name, fold, sep = "__fold_")]] <- model
        training_audit[[counter]] <- tibble::tibble(
          model = model_name,
          fold = fold,
          cull_training_sites = nrow(cull_train),
          manta_training_sites = if (use_manta) nrow(manta_train) else 0L,
          manta_site_weight = manta_weight,
          effective_manta_weight = if (use_manta) {
            nrow(manta_train) * manta_weight
          } else {
            0
          },
          test_cull_sites = nrow(test)
        )
      }
    }
  }

  list(
    predictions = dplyr::bind_rows(predictions),
    models = fitted_models,
    training_audit = dplyr::bind_rows(training_audit)
  )
}

run_rg_coverage_hybrid_benchmark <- function(
    clean_sites, rg_sites, clean_features, rg_features, model_cfg) {
  folds <- sort(unique(clean_sites$site_fold))
  predictions <- list()
  models <- list()

  for (fold in folds) {
    clean_train <- dplyr::filter(clean_sites, .data$site_fold != fold)
    rg_train <- dplyr::filter(rg_sites, .data$site_fold != fold)
    test <- dplyr::filter(clean_sites, .data$site_fold == fold)
    clean_train$observation_weight <- 1
    rg_train$observation_weight <- 1

    clean_model <- fit_xgb_binary(
      clean_train,
      clean_features,
      model_cfg,
      seed_offset = 40000L + fold,
      weight_col = "observation_weight"
    )
    rg_model <- fit_xgb_binary(
      rg_train,
      rg_features,
      model_cfg,
      seed_offset = 41000L + fold,
      weight_col = "observation_weight"
    )
    clean_probability <- predict_xgb_binary(
      clean_model, test, clean_features
    )
    hybrid_probability <- clean_probability
    rg_available <- test$rg_complete == 1L
    hybrid_probability[rg_available] <- predict_xgb_binary(
      rg_model,
      test[rg_available, , drop = FALSE],
      rg_features
    )

    common <- tibble::tibble(
      validation = "within_reef_site_cv_full_clean_domain",
      fold = fold,
      site_id = test$site_id,
      reef_name = test$reef_name,
      truth = test$high_density,
      n_surveyed_years = test$n_surveyed_years,
      longitude = test$longitude,
      latitude = test$latitude,
      rg_available = rg_available
    )
    predictions[[length(predictions) + 1L]] <- dplyr::mutate(
      common,
      model = "clean_full_domain_cull_only",
      probability = clean_probability,
      prediction_source = "clean"
    )
    predictions[[length(predictions) + 1L]] <- dplyr::mutate(
      common,
      model = "rg_with_clean_fallback",
      probability = hybrid_probability,
      prediction_source = dplyr::if_else(
        .data$rg_available, "rg", "clean_fallback"
      )
    )
    models[[paste0("clean__fold_", fold)]] <- clean_model
    models[[paste0("rg__fold_", fold)]] <- rg_model
  }

  list(
    predictions = dplyr::bind_rows(predictions),
    models = models
  )
}

run_site_year_rank_benchmark <- function(
    site_year, feature_sets, model_cfg) {
  folds <- sort(unique(site_year$site_fold))
  predictions <- list()
  models <- list()
  counter <- 0L

  for (family in names(feature_sets)) {
    features <- feature_sets[[family]]
    for (fold in folds) {
      counter <- counter + 1L
      train <- dplyr::filter(site_year, .data$site_fold != fold)
      test <- dplyr::filter(site_year, .data$site_fold == fold)
      model <- fit_xgb_binary(
        train,
        features,
        model_cfg,
        seed_offset = 50000L + counter,
        weight_col = "observation_weight"
      )
      probability <- predict_xgb_binary(model, test, features)
      predictions[[counter]] <- tibble::tibble(
        model = family,
        validation = "held_out_site_within_reef_annual_response",
        fold = fold,
        site_year_id = test$site_year_id,
        site_id = test$site_id,
        reef_name = test$reef_name,
        reef_year_id = test$reef_year_id,
        year = test$year,
        truth = test$high_density,
        probability = probability,
        longitude = test$longitude,
        latitude = test$latitude,
        site_history_weight = test$observation_weight
      )
      models[[paste0(family, "__fold_", fold)]] <- model
    }
  }

  list(
    predictions = dplyr::bind_rows(predictions),
    models = models
  )
}

site_habitat_importance <- function(models, feature_sets) {
  dplyr::bind_rows(lapply(names(feature_sets), function(family) {
    dplyr::bind_rows(lapply(
      c("cull_only", "plus_manta"),
      function(augmentation) {
        prefix <- paste0(family, "_", augmentation, "__fold_")
        keys <- grep(
          paste0("^", prefix), names(models), value = TRUE
        )
        fold_importance <- dplyr::bind_rows(lapply(keys, function(key) {
          value <- xgboost::xgb.importance(model = models[[key]])
          if (nrow(value) == 0L) return(NULL)
          tibble::as_tibble(value) |>
            dplyr::mutate(model_key = key)
        }))
        if (nrow(fold_importance) == 0L) return(NULL)
        fold_importance |>
          dplyr::group_by(.data$Feature) |>
          dplyr::summarise(
            mean_gain = mean(.data$Gain),
            folds_used = dplyr::n_distinct(.data$model_key),
            .groups = "drop"
          ) |>
          dplyr::mutate(
            model = paste(family, augmentation, sep = "_")
          )
      }
    ))
  }))
}
