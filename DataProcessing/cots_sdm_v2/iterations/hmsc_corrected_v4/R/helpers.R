classify_manta_evidence <- function(data, cfg) {
  code <- tolower(trimws(as.character(data$scar_code)))
  direct <- dplyr::coalesce(as.numeric(data$direct_cots_count) > 0, FALSE)
  multiple_scars <- dplyr::coalesce(
    as.numeric(data$scar_count) >= cfg$manta$strong_scar_count, FALSE
  )
  common_scars <- code %in% cfg$manta$strong_scar_codes
  weak_scar <- !direct & !multiple_scars & !common_scars & (
    dplyr::coalesce(as.numeric(data$scar_count) == 1, FALSE) |
      code %in% cfg$manta$weak_scar_codes
  )
  strict_negative <- dplyr::coalesce(as.logical(data$strict_negative), FALSE)
  dplyr::mutate(
    data,
    manta_direct_positive = direct,
    manta_multiple_scar_positive = multiple_scars,
    manta_common_scar_positive = common_scars,
    manta_strong_positive = direct | multiple_scars | common_scars,
    manta_weak_scar = weak_scar,
    manta_strict_negative = strict_negative,
    manta_evidence_class = dplyr::case_when(
      direct ~ "direct_cots",
      multiple_scars ~ "multiple_counted_scars",
      common_scars ~ "legacy_common_scars",
      weak_scar ~ "weak_single_or_present_scar",
      strict_negative ~ "strict_negative",
      TRUE ~ "uninformative"
    )
  )
}

observed_any_v4 <- function(x) {
  if (length(x) == 0L || all(is.na(x))) return(NA_integer_)
  as.integer(any(x == 1L, na.rm = TRUE))
}

mean_or_na_v4 <- function(x) {
  if (!length(x) || all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

weighted_mean_or_na_v4 <- function(x, w) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  if (!any(keep)) return(NA_real_)
  stats::weighted.mean(x[keep], w[keep])
}

make_spatial_fold_v4 <- function(x, y, block_m = 3000L, k = 5L) {
  1L + abs(floor(x / block_m) + 3L * floor(y / block_m)) %% k
}

normalise_manta_weights_v4 <- function(weighted_length_m, fraction) {
  dplyr::if_else(
    is.finite(weighted_length_m) & weighted_length_m > 0,
    weighted_length_m,
    dplyr::if_else(is.finite(fraction) & fraction > 0, fraction, 1)
  )
}

model_matrix_pair_v4 <- function(train, test, terms) {
  combined <- dplyr::bind_rows(
    dplyr::mutate(train, .v4_split = "train"),
    dplyr::mutate(test, .v4_split = "test")
  )
  if ("geomorphic_class" %in% terms) {
    combined$geomorphic_class <- factor(combined$geomorphic_class)
  }
  formula <- stats::as.formula(paste("~", paste(terms, collapse = " + "), "- 1"))
  matrix <- stats::model.matrix(formula, data = combined)
  list(
    train = matrix[combined$.v4_split == "train", , drop = FALSE],
    test = matrix[combined$.v4_split == "test", , drop = FALSE]
  )
}

safe_auc_v4 <- function(truth, probability) {
  keep <- is.finite(truth) & is.finite(probability)
  if (sum(keep) < 2L || length(unique(truth[keep])) < 2L) return(NA_real_)
  as.numeric(yardstick::roc_auc_vec(
    factor(truth[keep], levels = c(0, 1)), probability[keep], event_level = "second"
  ))
}

within_reef_year_auc_v4 <- function(data) {
  groups <- split(data, data$reef_year_id)
  rows <- lapply(groups, function(x) {
    positives <- sum(x$truth == 1L, na.rm = TRUE)
    negatives <- sum(x$truth == 0L, na.rm = TRUE)
    if (!positives || !negatives) return(NULL)
    data.frame(
      auc = safe_auc_v4(x$truth, x$probability),
      pairs = positives * negatives
    )
  })
  rows <- dplyr::bind_rows(rows)
  if (!nrow(rows)) return(NA_real_)
  stats::weighted.mean(rows$auc, rows$pairs, na.rm = TRUE)
}

rank_raster_v4 <- function(x) {
  value <- terra::values(x, mat = FALSE)
  valid <- is.finite(value)
  result <- rep(NA_real_, length(value))
  if (sum(valid) > 1L) {
    result[valid] <- (rank(value[valid], ties.method = "average") - 1) /
      (sum(valid) - 1)
  }
  terra::setValues(x, result)
}
