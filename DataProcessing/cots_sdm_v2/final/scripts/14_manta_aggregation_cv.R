# Can manta surveillance effort be cut without missing COTS aggregations?
#
# Confound: the primary layer is trained on Cull dives, and Cull crews mostly
# dive where manta tows found COTS or scars. Scoring manta tows with the
# deployed layer therefore partly scores the layer against its own source.
# Here every reef's tows are scored by a model that never saw that reef's
# Cull data (leave-reef-out cross-validation, 5 folds of reefs, refitting the
# RG + year XGBoost with the legacy hyperparameters).
#
# Aggregations: within a survey (reef x date), strong-evidence tows (direct
# COTS, >= 2 counted scars, legacy "common" scars) whose midpoints lie within
# 500 m of each other are linked into one aggregation. An aggregation is found
# if at least one of its tows is kept (that tow triggers a closer search).
#
# Thinning strategies, each keeping the same number of tows per survey:
#   layer       top-ranked tows by the cross-validated score
#   random      exact hypergeometric expectation
#   systematic  every n-th tow around the reef (tows ordered by bearing from
#               the survey centre), averaged over 20 starting offsets
# Uncertainty: 500 bootstrap resamples of surveys.

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "..", "R", "common.R"))
cfg <- initialise_final()
set.seed(cfg$hmsc$seed)
keep_fractions <- c(1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3)
link_distance_m <- 500
n_folds <- 5L
boot_reps <- 500L
strong_classes <- c("direct_cots", "multiple_counted_scars", "legacy_common_scars")
legacy_params <- list(objective = "binary:logistic", eta = 0.00137382379588326, max_depth = 5,
                      min_child_weight = 20, gamma = 1e-10, nthread = 8)
legacy_rounds <- 1500L

stack <- terra::rast(checked_input(cfg, "predictor_stack_rg"))
names(stack) <- cfg$legacy$rg_layer_names
features <- c(cfg$legacy$rg_layer_names, "Year")

# --- Cull training data (site-years) -----------------------------------------
sites <- readRDS(file.path(cfg$exploration_dir, "artifacts", "site_year_all.rds")) |>
  dplyr::filter(!is.na(.data$high_density), is.finite(.data$longitude), is.finite(.data$latitude))
sxy <- sf::st_coordinates(sf::st_transform(sf::st_as_sf(sites, coords = c("longitude", "latitude"), crs = 4326), terra::crs(stack)))
site_x <- terra::extract(stack, sxy)
site_x$Year <- sites$year

# --- Manta tows: geometry, evidence, points along each tow ------------------
tow_cache <- file.path(cfg$dirs$cache, "manta_tow_predictors.rds")
tows <- readRDS(checked_input(cfg, "manta_tows")) |>
  dplyr::select("visit_id", "survey_date", "year", "reef_name", "reef_year_id", "manta_evidence_class")
geom <- readRDS(checked_input(cfg, "footprint_links"))$manta |>
  dplyr::group_by(.data$visit_id) |>
  dplyr::summarise(sx = dplyr::first(.data$start_x), sy = dplyr::first(.data$start_y),
                   ex = dplyr::first(.data$end_x), ey = dplyr::first(.data$end_y), .groups = "drop")
tows <- dplyr::inner_join(tows, geom, by = "visit_id") |>
  dplyr::filter(is.finite(.data$sx), is.finite(.data$ex)) |>
  dplyr::mutate(strong = .data$manta_evidence_class %in% strong_classes,
                survey = paste(.data$reef_name, .data$survey_date),
                mx = (.data$sx + .data$ex) / 2, my = (.data$sy + .data$ey) / 2,
                length = sqrt((.data$ex - .data$sx)^2 + (.data$ey - .data$sy)^2))
if (!file.exists(tow_cache)) {
  log_message("Extracting predictors along ", format(nrow(tows), big.mark = ","), " tows")
  n_pts <- pmin(15L, pmax(2L, ceiling(tows$length / 50)))
  idx <- rep(seq_len(nrow(tows)), n_pts)
  f <- unlist(lapply(n_pts, function(n) seq(0, 1, length.out = n)))
  pts <- cbind(tows$sx[idx] + f * (tows$ex[idx] - tows$sx[idx]), tows$sy[idx] + f * (tows$ey[idx] - tows$sy[idx]))
  saveRDS(list(tow = idx, values = terra::extract(stack, pts)), tow_cache)
}
tow_points <- readRDS(tow_cache)

# --- Leave-reef-out scores ---------------------------------------------------------
reefs <- unique(c(sites$reef_name, tows$reef_name))
fold_of <- stats::setNames(sample(rep(seq_len(n_folds), length.out = length(reefs))), reefs)
point_year <- tows$year[tow_points$tow]
tows$cv_score <- NA_real_
for (k in seq_len(n_folds)) {
  train <- fold_of[sites$reef_name] != k
  dtrain <- xgboost::xgb.DMatrix(as.matrix(site_x[train, features]), label = sites$high_density[train], missing = NA)
  model <- xgboost::xgb.train(params = legacy_params, data = dtrain, nrounds = legacy_rounds, verbose = 0)
  test_points <- which(fold_of[tows$reef_name[tow_points$tow]] == k)
  x <- tow_points$values[test_points, cfg$legacy$rg_layer_names]
  x$Year <- point_year[test_points]
  p <- stats::predict(model, xgboost::xgb.DMatrix(as.matrix(x[, features]), missing = NA))
  s <- tapply(p, tow_points$tow[test_points], mean, na.rm = TRUE)
  tows$cv_score[as.integer(names(s))] <- s
  log_message("Fold ", k, ": trained on ", sum(train), " site-years; scored ", length(s), " tows")
}
# For comparison: the deployed layer (trained on every reef's Cull data).
search120 <- terra::rast(file.path(cfg$dirs$cache, "search_planning_120m.tif"))
tows$deployed_score <- terra::extract(search120, cbind(tows$mx, tows$my))[[1]]

# --- Surveys, aggregations and tow order ------------------------------------------
tows <- tows |>
  dplyr::filter(is.finite(.data$cv_score), is.finite(.data$deployed_score)) |>
  dplyr::group_by(.data$survey) |>
  dplyr::filter(dplyr::n() >= 5, any(.data$strong)) |>
  dplyr::mutate(bearing = atan2(.data$my - mean(.data$my), .data$mx - mean(.data$mx)),
                order = rank(.data$bearing, ties.method = "first")) |>
  dplyr::ungroup()
cluster_strong <- function(d) {
  out <- rep(NA_integer_, nrow(d))
  s <- which(d$strong)
  out[s] <- if (length(s) == 1) 1L else stats::cutree(stats::hclust(stats::dist(cbind(d$mx[s], d$my[s])), "single"), h = link_distance_m)
  out
}
tows <- tows |>
  dplyr::group_by(.data$survey) |>
  dplyr::group_modify(~ dplyr::mutate(.x, aggregation = cluster_strong(.x))) |>
  dplyr::ungroup()
aggs <- tows |>
  dplyr::filter(.data$strong) |>
  dplyr::count(.data$survey, .data$aggregation, name = "size") |>
  dplyr::mutate(size_class = cut(.data$size, c(0, 1, 2, 5, Inf), labels = c("1 tow (isolated)", "2 tows", "3-5 tows", "6+ tows")))
log_message(format(dplyr::n_distinct(tows$survey), big.mark = ","), " surveys with strong evidence; ",
            nrow(aggs), " aggregations (", sum(aggs$size >= 2), " with >= 2 strong tows)")

# Per survey and keep fraction: P(each aggregation found) under each strategy.
survey_outcomes <- function(d, k) {
  n <- nrow(d)
  kept_n <- ceiling(k * n)
  agg_ids <- sort(unique(stats::na.omit(d$aggregation)))
  found <- function(kept) vapply(agg_ids, function(a) as.numeric(any(kept & d$aggregation %in% a)), numeric(1))
  layer_cv <- found(rank(-d$cv_score, ties.method = "first") <= kept_n)
  layer_deployed <- found(rank(-d$deployed_score, ties.method = "first") <= kept_n)
  sizes <- vapply(agg_ids, function(a) sum(d$aggregation %in% a), numeric(1))
  random <- 1 - exp(lchoose(n - sizes, kept_n) - lchoose(n, kept_n))
  random[kept_n > n - sizes] <- 1
  phases <- seq(0, 1, length.out = 21)[-21]
  systematic <- rowMeans(matrix(vapply(phases, function(phase) {
    pos <- (d$order - 1 + phase * n / kept_n) %% (n / kept_n)
    found(pos < 1)
  }, numeric(length(agg_ids))), nrow = length(agg_ids)))
  data.frame(survey = d$survey[1], year = d$year[1], keep = k, tows = n, kept = kept_n,
             aggregation = agg_ids, size = sizes,
             layer_cv = layer_cv, layer_deployed = layer_deployed, random = random, systematic = systematic)
}
surveys <- split(tows, tows$survey)
outcomes <- dplyr::bind_rows(lapply(keep_fractions, function(k) dplyr::bind_rows(lapply(surveys, survey_outcomes, k = k))))
outcomes$size_class <- cut(outcomes$size, c(0, 1, 2, 5, Inf), labels = levels(aggs$size_class))
saveRDS(outcomes, file.path(cfg$dirs$cache, "manta_aggregation_outcomes.rds"))

# --- Bootstrap over surveys ------------------------------------------------------
strategies <- c("layer_cv", "layer_deployed", "random", "systematic")
# Pre-aggregate per survey and keep fraction, then resample surveys with
# weights (equivalent to resampling the rows survey by survey).
bootstrap_summary <- function(o, label) {
  per_survey <- o |>
    dplyr::group_by(.data$survey, .data$keep) |>
    dplyr::summarise(n = dplyr::n(), dplyr::across(dplyr::all_of(strategies), sum), .groups = "drop")
  ids <- sort(unique(per_survey$survey))
  ks <- sort(unique(per_survey$keep), decreasing = TRUE)
  cube <- function(column) {
    m <- matrix(0, length(ids), length(ks))
    m[cbind(match(per_survey$survey, ids), match(per_survey$keep, ks))] <- per_survey[[column]]
    m
  }
  n_m <- cube("n")
  s_m <- lapply(stats::setNames(strategies, strategies), cube)
  rate <- function(w) vapply(s_m, function(m) colSums(m * w) / colSums(n_m * w), numeric(length(ks)))
  point <- rate(rep(1, length(ids)))
  boot <- replicate(boot_reps, rate(tabulate(sample.int(length(ids), replace = TRUE), length(ids))), simplify = "array")
  q <- function(x, p) apply(x, 1, stats::quantile, p)
  gain_r <- boot[, "layer_cv", ] - boot[, "random", ]
  gain_s <- boot[, "layer_cv", ] - boot[, "systematic", ]
  out <- tibble::tibble(subset = label, aggregations = sum(n_m[, 1]), surveys = length(ids), keep = ks)
  for (st in strategies) {
    out[[st]] <- point[, st]
    out[[paste0(st, "_low")]] <- q(boot[, st, , drop = FALSE][, 1, ], 0.025)
    out[[paste0(st, "_high")]] <- q(boot[, st, , drop = FALSE][, 1, ], 0.975)
  }
  out$gain_vs_random <- point[, "layer_cv"] - point[, "random"]
  out$gain_vs_random_low <- q(gain_r, 0.025)
  out$gain_vs_random_high <- q(gain_r, 0.975)
  out$gain_vs_systematic <- point[, "layer_cv"] - point[, "systematic"]
  out$gain_vs_systematic_low <- q(gain_s, 0.025)
  out$gain_vs_systematic_high <- q(gain_s, 0.975)
  out
}
subsets <- list(
  `aggregations (>= 2 strong tows), all years` = outcomes$size >= 2,
  `aggregations (>= 2 strong tows), 2022-2026` = outcomes$size >= 2 & outcomes$year >= 2022,
  `isolated strong tows, all years` = outcomes$size == 1
)
summary <- dplyr::bind_rows(lapply(names(subsets), function(s) bootstrap_summary(outcomes[subsets[[s]], ], s)))
write_table(summary, cfg, "manta_aggregation_cv_summary.csv")

by_size <- outcomes |>
  dplyr::group_by(.data$keep, .data$size_class) |>
  dplyr::summarise(aggregations = dplyr::n(), dplyr::across(dplyr::all_of(strategies), mean), .groups = "drop")
write_table(by_size, cfg, "manta_aggregation_cv_by_size.csv")

# Largest cut that still finds >= 95 % of aggregations (point estimate and
# lower 95 % bootstrap bound).
target <- 0.95
largest_cut <- function(keep, ok) if (any(ok)) 1 - min(keep[ok]) else 0
max_cut <- summary |>
  dplyr::group_by(.data$subset) |>
  dplyr::summarise(
    layer_cv_point = largest_cut(.data$keep, .data$layer_cv >= target),
    layer_cv_assured = largest_cut(.data$keep, .data$layer_cv_low >= target),
    random_point = largest_cut(.data$keep, .data$random >= target),
    systematic_point = largest_cut(.data$keep, .data$systematic >= target),
    systematic_assured = largest_cut(.data$keep, .data$systematic_low >= target),
    .groups = "drop"
  )
write_table(max_cut, cfg, "manta_aggregation_cv_max_cut.csv")
print(summary[, c("subset", "keep", "aggregations", "layer_cv", "layer_cv_low", "layer_cv_high", "layer_deployed", "random", "systematic")], n = Inf)
print(max_cut)
