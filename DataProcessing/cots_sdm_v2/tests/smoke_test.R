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

cat("All COTS SDM v2 smoke tests passed.\n")
