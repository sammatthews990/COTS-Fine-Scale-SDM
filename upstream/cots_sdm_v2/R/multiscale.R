scale_suffix <- function(factor, native_resolution_m = 30) {
  paste0("s", as.integer(factor * native_resolution_m), "m")
}

project_points_to_template <- function(data, template) {
  points <- terra::vect(
    as.data.frame(data),
    geom = c("longitude", "latitude"),
    crs = "EPSG:4326",
    keepgeom = FALSE
  )
  terra::project(points, terra::crs(template))
}

aligned_grid_membership <- function(data, template, factor = 1L) {
  factor <- as.integer(factor)
  if (factor < 1L) stop("Aggregation factor must be at least one.")
  projected <- project_points_to_template(data, template)
  xy <- terra::crds(projected)
  native_cell <- terra::cellFromXY(template, xy)
  row_col <- terra::rowColFromCell(template, native_cell)
  coarse_ncol <- ceiling(terra::ncol(template) / factor)
  coarse_row <- (row_col[, 1] - 1L) %/% factor + 1L
  coarse_col <- (row_col[, 2] - 1L) %/% factor + 1L
  coarse_cell <- (coarse_row - 1L) * coarse_ncol + coarse_col

  tibble::tibble(
    native_cell = native_cell,
    native_row = row_col[, 1],
    native_col = row_col[, 2],
    coarse_row = coarse_row,
    coarse_col = coarse_col,
    coarse_cell = coarse_cell,
    projected_x = xy[, 1],
    projected_y = xy[, 2]
  )
}

coarse_block_native_cells <- function(template, membership, factor) {
  blocks <- unique(membership[, c("coarse_cell", "coarse_row", "coarse_col")])
  pieces <- lapply(seq_len(nrow(blocks)), function(i) {
    row_start <- (blocks$coarse_row[i] - 1L) * factor + 1L
    col_start <- (blocks$coarse_col[i] - 1L) * factor + 1L
    rows <- row_start:min(row_start + factor - 1L, terra::nrow(template))
    cols <- col_start:min(col_start + factor - 1L, terra::ncol(template))
    grid <- expand.grid(row = rows, col = cols)
    tibble::tibble(
      coarse_cell = blocks$coarse_cell[i],
      native_cell = terra::cellFromRowCol(
        template, grid$row, grid$col
      )
    )
  })
  dplyr::bind_rows(pieces)
}

extract_aligned_block_means <- function(
    predictors, data, factor = 1L, native_resolution_m = 30) {
  membership <- aligned_grid_membership(data, predictors, factor)
  valid <- is.finite(membership$coarse_cell)
  cell_map <- coarse_block_native_cells(
    predictors, membership[valid, , drop = FALSE], factor
  )
  values <- as.matrix(predictors[cell_map$native_cell])
  values[!is.finite(values)] <- NA_real_
  group <- as.character(cell_map$coarse_cell)
  sums <- rowsum(replace(values, is.na(values), 0), group, reorder = FALSE)
  counts <- rowsum(1L * !is.na(values), group, reorder = FALSE)
  means <- sums / counts
  means[counts == 0] <- NA_real_
  means <- as.data.frame(means, check.names = FALSE)
  means$coarse_cell <- as.numeric(rownames(means))
  rownames(means) <- NULL

  suffix <- scale_suffix(factor, native_resolution_m)
  names(means)[names(means) != "coarse_cell"] <- paste0(
    names(means)[names(means) != "coarse_cell"], "_", suffix
  )
  dplyr::bind_cols(
    data,
    membership,
    dplyr::select(
      dplyr::left_join(
        tibble::tibble(coarse_cell = membership$coarse_cell),
        means,
        by = "coarse_cell"
      ),
      -"coarse_cell"
    )
  )
}

attach_multiscale_block_means <- function(
    data, predictors, factors = c(1L, 2L, 4L),
    native_resolution_m = 30) {
  result <- data
  for (factor in factors) {
    extracted <- extract_aligned_block_means(
      predictors,
      dplyr::select(data, "longitude", "latitude"),
      factor = factor,
      native_resolution_m = native_resolution_m
    )
    suffix <- scale_suffix(factor, native_resolution_m)
    value_names <- grep(paste0("_", suffix, "$"), names(extracted), value = TRUE)
    if (factor == 1L) {
      membership_names <- c(
        "native_cell", "native_row", "native_col",
        "projected_x", "projected_y"
      )
      result <- dplyr::bind_cols(
        result,
        dplyr::select(extracted, dplyr::all_of(membership_names))
      )
    }
    result[[paste0("cell_", suffix)]] <- extracted$coarse_cell
    result <- dplyr::bind_cols(
      result,
      dplyr::select(extracted, dplyr::all_of(value_names))
    )
  }
  result
}

multiscale_feature_sets <- function(base_names) {
  topographic <- intersect(
    c("BPI", "EAST", "HCU", "NORTH", "SLO", "SVF", "VCU", "VRM", "DEM"),
    base_names
  )
  biological <- setdiff(base_names, topographic)
  at_scale <- function(names, scale) paste0(names, "_s", scale, "m")
  list(
    clean_30m_year = c(at_scale(base_names, 30), "year"),
    clean_60m_year = c(at_scale(base_names, 60), "year"),
    clean_120m_year = c(at_scale(base_names, 120), "year"),
    multiscale_topography_year = c(
      unlist(lapply(c(30, 60, 120), function(x) at_scale(topographic, x))),
      at_scale(biological, 30),
      "year"
    ),
    multiscale_all_year = c(
      unlist(lapply(c(30, 60, 120), function(x) at_scale(base_names, x))),
      "year"
    ),
    multiscale_topography_static = c(
      unlist(lapply(c(30, 60, 120), function(x) at_scale(topographic, x))),
      at_scale(biological, 30)
    )
  )
}

aggregate_site_year_to_resolution <- function(
    data, resolution_m, native_resolution_m = 30) {
  cell_name <- paste0("cell_s", resolution_m, "m")
  if (!cell_name %in% names(data)) stop("Missing ", cell_name)
  scale_names <- grep(paste0("_s", resolution_m, "m$"), names(data), value = TRUE)
  scale_names <- setdiff(scale_names, cell_name)
  result <- data |>
    dplyr::group_by(.data$reef_name, .data$year, .data[[cell_name]]) |>
    dplyr::summarise(
      longitude = mean(.data$longitude),
      latitude = mean(.data$latitude),
      total_cots = sum(.data$total_cots),
      bottom_time = sum(.data$bottom_time),
      n_source_sites = dplyr::n_distinct(.data$site_id),
      dplyr::across(dplyr::all_of(scale_names), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      cpue_mean = .data$total_cots / .data$bottom_time,
      high_density = as.integer(.data$cpue_mean >= 0.02),
      site_id = paste0("cell_", resolution_m, "m||", .data[[cell_name]]),
      site_year_id = paste(.data$site_id, .data$year, sep = "||"),
      reef_year_id = paste(.data$reef_name, .data$year, sep = "||"),
      validation_group_id = paste(.data$reef_name, .data$site_id, sep = "||")
    )
  result <- assign_within_reef_site_folds(result, k = 5L, seed = 20260831L)
  result |>
    dplyr::group_by(.data$site_id) |>
    dplyr::mutate(observation_weight = 1 / dplyr::n_distinct(.data$year)) |>
    dplyr::ungroup()
}
