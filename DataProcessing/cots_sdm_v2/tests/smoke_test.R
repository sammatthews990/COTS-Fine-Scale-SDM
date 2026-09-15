args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v2_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(v2_dir, "R", "bootstrap.R"))
ctx <- initialise_v2(script_file)

synthetic_outbreak <- expand.grid(
  reef_name = c("reef_a", "reef_b"),
  year = 2000:2005,
  stringsAsFactors = FALSE
)
synthetic_outbreak$reef_outbreak_probability <- seq(
  0.1, 0.9, length.out = nrow(synthetic_outbreak)
)

synthetic_site_year <- expand.grid(
  reef_name = c("reef_a", "reef_b", "reef_c", "reef_d"),
  year = 2003:2006,
  stringsAsFactors = FALSE
)
synthetic_site_year$site_id <- paste0(synthetic_site_year$reef_name, "_site")
synthetic_site_year$site_year_id <- paste(
  synthetic_site_year$site_id,
  synthetic_site_year$year,
  sep = "||"
)
synthetic_site_year$longitude <- unname(c(
  reef_a = 145, reef_b = 147, reef_c = 149, reef_d = 151
)[synthetic_site_year$reef_name])
synthetic_site_year$latitude <- unname(c(
  reef_a = -15, reef_b = -17, reef_c = -19, reef_d = -21
)[synthetic_site_year$reef_name])
synthetic_site_year$high_density <- rep(c(0L, 1L), 8)

lagged <- add_outbreak_history_features(synthetic_site_year, synthetic_outbreak)
expected <- synthetic_outbreak$reef_outbreak_probability[
  synthetic_outbreak$reef_name == "reef_a" & synthetic_outbreak$year == 2002
]
observed <- lagged$reef_outbreak_probability_lag1[
  lagged$reef_name == "reef_a" & lagged$year == 2003
]
stopifnot(length(observed) == 1L, isTRUE(all.equal(observed, expected)))

synthetic_forecast <- tibble::tibble(
  reef_name = c("reef_a", "reef_b"),
  year = 2006L,
  reef_outbreak_probability = c(0.8, 0.2),
  outbreak_probability_cutoff = 0.45
)
with_forecast <- attach_outbreak_forecast(lagged, synthetic_forecast)
stopifnot(
  with_forecast$reef_outbreak_probability_forecast[
    with_forecast$reef_name == "reef_a" & with_forecast$year == 2006L
  ] == 0.8,
  is.na(with_forecast$reef_outbreak_probability_forecast[
    with_forecast$reef_name == "reef_a" & with_forecast$year == 2005L
  ])
)

synthetic_reef_evidence <- tibble::tibble(
  reef_name = c("reef_a", "reef_b"),
  manta_tows = c(12L, 12L),
  survey_years = c(3L, 3L),
  all_tows_informative = TRUE,
  all_tows_strict_negative = TRUE
)
synthetic_current_manta <- tibble::tibble(
  reef_name = c("reef_a", "reef_b"),
  year = 2006L,
  positive_evidence = c(FALSE, TRUE)
)
classified <- classify_manta_negative_reefs(
  synthetic_reef_evidence,
  synthetic_current_manta,
  minimum_tows = 10L,
  minimum_survey_years = 3L
)
stopifnot(
  classified$strong_negative_evidence[classified$reef_name == "reef_a"],
  !classified$strong_negative_evidence[classified$reef_name == "reef_b"]
)

folded <- assign_spatial_reef_folds(
  synthetic_site_year,
  k = 2L,
  seed = 1L
)$data
reef_fold_counts <- folded |>
  dplyr::select(dplyr::all_of(c("reef_name", "spatial_fold"))) |>
  dplyr::distinct() |>
  dplyr::count(.data$reef_name)
stopifnot(all(reef_fold_counts$n == 1L))

split <- make_forward_year_split(synthetic_site_year, 2006L)
stopifnot(
  max(synthetic_site_year$year[split$train]) == 2005L,
  all(synthetic_site_year$year[split$test] == 2006L)
)

metric_check <- suppressWarnings(binary_probability_metrics(
  c(0L, 0L, 1L, 1L),
  c(0.1, 0.3, 0.7, 0.9)
))
stopifnot(
  metric_check$roc_auc == 1,
  metric_check$brier < 0.1
)

scale_check <- prediction_scale_decomposition(tibble::tibble(
  model = "synthetic",
  reef_name = rep(c("reef_a", "reef_b"), each = 2L),
  truth = c(0L, 1L, 0L, 1L),
  probability = c(0.1, 0.9, 0.8, 0.2)
))
stopifnot(
  scale_check$mixed_reefs == 2L,
  isTRUE(all.equal(scale_check$within_reef_roc_auc, 0.5))
)

rank_predictions <- tibble::tibble(
  model = "synthetic_rank",
  validation = "smoke_test",
  reef_year_id = rep(c("reef_a||2006", "reef_b||2006"), each = 4L),
  truth = rep(c(0L, 0L, 1L, 1L), 2L),
  probability = c(
    0.490, 0.500, 0.520, 0.530,
    0.491, 0.501, 0.521, 0.531
  )
)
raw_rank_auc <- within_group_auc(
  rank_predictions, group = "reef_year_id"
)$summary$within_group_roc_auc
percentile_rank_auc <- within_group_auc(
  within_group_percentile_predictions(
    rank_predictions, group = "reef_year_id"
  ),
  group = "reef_year_id"
)$summary$within_group_roc_auc
capture_check <- within_group_capture_metrics(
  rank_predictions,
  group = "reef_year_id",
  fractions = 0.5
)
stopifnot(
  raw_rank_auc == 1,
  percentile_rank_auc == raw_rank_auc,
  capture_check$capture == 1,
  capture_check$lift == 2
)

calibration_probability <- c(0.10, 0.20, 0.30, 0.40)
calibration_truth <- c(0L, 1L, 1L, 1L)
calibration_group <- rep("2006", 4L)
calibrator_check <- fit_group_intercept_calibrator(
  calibration_probability,
  calibration_truth,
  calibration_group,
  prior_weight = 0
)
calibrated_probability <- predict_group_intercept_calibrator(
  calibrator_check,
  calibration_probability,
  calibration_group
)
stopifnot(
  isTRUE(all.equal(mean(calibrated_probability), 0.75, tolerance = 1e-6)),
  identical(order(calibrated_probability), order(calibration_probability))
)

bridge_check <- fit_additive_operational_bridge(
  truth = rep(
    c(0L, 0L, 1L, 0L, 1L, 0L, 1L, 1L, 0L, 1L),
    8L
  ),
  static_probability = rep(
    c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.55, 0.45),
    8L
  ),
  outbreak_probability = rep(c(0.3, 0.7), each = 40L)
)
bridge_probability <- predict_additive_operational_bridge(
  bridge_check,
  static_probability = c(0.25, 0.75),
  outbreak_probability = c(0.5, 0.5)
)
stopifnot(
  bridge_check$rank_preserving_within_reef,
  bridge_probability[2] > bridge_probability[1]
)

site_fold_check <- assign_within_reef_site_folds(
  tibble::tibble(
    reef_name = c("reef_a", "reef_a", "reef_b", "reef_c"),
    validation_group_id = c("a1", "a1", "b1", "c1")
  ),
  k = 3L,
  seed = 1L
)
stopifnot(
  !anyNA(site_fold_check$site_fold),
  dplyr::n_distinct(
    site_fold_check$site_fold[
      site_fold_check$validation_group_id == "a1"
    ]
  ) == 1L
)

multiscale_raster <- terra::rast(
  nrows = 4, ncols = 4, xmin = 0, xmax = 4, ymin = 0, ymax = 4,
  crs = "EPSG:4326"
)
terra::values(multiscale_raster) <- 1:16
names(multiscale_raster) <- "terrain"
multiscale_check <- extract_aligned_block_means(
  multiscale_raster,
  data.frame(longitude = c(0.5, 2.5), latitude = c(3.5, 1.5)),
  factor = 2L,
  native_resolution_m = 1
)
stopifnot(isTRUE(all.equal(
  multiscale_check$terrain_s2m, c(3.5, 13.5)
)))

two_part_habitat <- rep(
  rep(c(0.20, 0.35, 0.45, 0.65, 0.75, 0.85), each = 5L),
  2L
)
two_part_truth_block <- unlist(lapply(
  c(1L, 1L, 2L, 3L, 4L, 4L),
  function(n) c(rep(1L, n), rep(0L, 5L - n))
))
two_part_truth <- rep(two_part_truth_block, 2L)
two_part_outbreak <- rep(c(0.15, 0.70), each = 30L)
two_part_bridge <- fit_rank_locked_state_bridge(
  two_part_truth, two_part_habitat, two_part_outbreak
)
two_part_state <- predict_rank_locked_state_bridge(
  two_part_bridge, two_part_habitat, two_part_outbreak
)
two_part_effort <- fit_effort_standardiser(
  two_part_truth,
  two_part_state,
  n_dives = rep(c(1, 2, 3), 20L),
  bottom_time = rep(c(300, 600, 900), 20L)
)
two_part_fixed <- predict_effort_standardised(
  two_part_effort, two_part_state,
  n_dives = 2, bottom_time = 600
)
two_part_check <- rank_lock_check(tibble::tibble(
  reef_year_id = rep(
    c("reef_a||2025", "reef_b||2025"), each = 30L
  ),
  habitat_probability = two_part_habitat,
  probability = two_part_fixed
))
stopifnot(
  two_part_bridge$rank_locked_within_reef_year,
  two_part_check$violations == 0L
)
rank_aligned <- rank_align_probability(
  reference = c(0.1, 0.2, 0.2, 0.4),
  probability = c(0.8, 0.1, 0.9, 0.3)
)
stopifnot(
  isTRUE(all.equal(
    sort(rank_aligned), c(0.1, 0.55, 0.55, 0.9)
  )),
  rank_aligned[[1]] <= rank_aligned[[2]],
  rank_aligned[[2]] == rank_aligned[[3]],
  rank_aligned[[3]] <= rank_aligned[[4]]
)

cat("All COTS SDM v2 smoke tests passed.\n")
