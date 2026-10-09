# Does a habitat model trained on manta tows (strong COTS evidence vs strict
# negatives) serve the programme goal better than the cull-trained layer?
#
# Goal: predict where on a reef COTS are most likely, so that crews can
#   (a) manta tow less without missing aggregations, and
#   (b) search fewer cull sites without missing problematic dives
#       (CPUE >= 0.02), i.e. avoid low-density and ghost dives.
#
# Manta tows cover far more of each reef than Cull dives, including places
# never dived. Both models use the same predictors (RG stack + survey year) and
# the legacy XGBoost hyperparameters (hist tree method), and are refitted on
# the same 5 folds of reefs, so every reef is scored by models that never saw
# its own tows or dives:
#   cull   response: Cull site-year CPUE >= 0.02 (the deployed layer's target)
#   manta  response: tow with strong evidence (direct COTS, >= 2 counted scars,
#          legacy "common" scars) = 1, strict negative = 0; weak single scars
#          are dropped. Tow predictors: mean along the tow.
#   manta within reef-year: same response, but the reef-year's strong-tow
#          rate is fixed as an offset (xgboost base_margin) and survey year is
#          dropped, so the trees can only learn which locations within a
#          reef-year carry more COTS (the ranking the programme uses).
# Validation:
#   tows   within reef-year AUC and top-20 % lift, strong vs negative; share of
#          aggregations found when thinning tows (as in 14)
#   cull   designated-site ranking as in 15 (problematic dives kept, low and
#          ghost dives avoided), and with last season's evidence as in 16,
#          compared at matched retention of problematic dives. Ensembles of
#          the cull-trained and within reef-year manta ranks are tested two
#          ways: their average rank, and both ranks as separate terms in the
#          evidence model.
# Note: the manta model's own test is in-kind (held-out tows), so the cull
# test is the independent one for it, and vice versa.

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "..", "R", "common.R"))
cfg <- initialise_final()
source(file.path(cfg$final_dir, "R", "cull_sites.R"))
set.seed(cfg$hmsc$seed + 17L)
n_folds <- 5L
boot_reps <- 500L
keep_fractions <- c(1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3)
site_fractions <- round(seq(0.05, 1, 0.05), 2)
retention_targets <- c(0.8, 0.85, 0.9, 0.95)
link_distance_m <- 500
strong_classes <- c("direct_cots", "multiple_counted_scars", "legacy_common_scars")
params <- list(objective = "binary:logistic", eta = 0.00137382379588326, max_depth = 5,
               min_child_weight = 20, gamma = 1e-10, nthread = 8, tree_method = "hist")
rounds <- 1500L

stack <- terra::rast(checked_input(cfg, "predictor_stack_rg"))
names(stack) <- cfg$legacy$rg_layer_names
rg <- cfg$legacy$rg_layer_names
features <- c(rg, "Year")

# --- Manta tows with mean predictors along each tow (cache from 14) ----------
tow_cache <- file.path(cfg$dirs$cache, "manta_tow_predictors.rds")
if (!file.exists(tow_cache)) stop("Run 14_manta_aggregation_cv.R first (needs ", tow_cache, ")")
tows <- readRDS(checked_input(cfg, "manta_tows")) |>
  dplyr::select("visit_id", "survey_date", "year", "reef_name", "manta_evidence_class")
geom <- readRDS(checked_input(cfg, "footprint_links"))$manta |>
  dplyr::group_by(.data$visit_id) |>
  dplyr::summarise(sx = dplyr::first(.data$start_x), sy = dplyr::first(.data$start_y),
                   ex = dplyr::first(.data$end_x), ey = dplyr::first(.data$end_y), .groups = "drop")
tows <- dplyr::inner_join(tows, geom, by = "visit_id") |>
  dplyr::filter(is.finite(.data$sx), is.finite(.data$ex)) |>
  dplyr::mutate(strong = .data$manta_evidence_class %in% strong_classes,
                label = dplyr::case_when(.data$strong ~ 1, .data$manta_evidence_class == "strict_negative" ~ 0),
                survey = paste(.data$reef_name, .data$survey_date),
                reef_year = paste(.data$reef_name, .data$year),
                mx = (.data$sx + .data$ex) / 2, my = (.data$sy + .data$ey) / 2) |>
  dplyr::group_by(.data$reef_year) |>
  dplyr::mutate(ry_margin = stats::qlogis((sum(.data$label %in% 1) + 0.5) / (sum(!is.na(.data$label)) + 1))) |>
  dplyr::ungroup()
tow_points <- readRDS(tow_cache)
v <- as.matrix(tow_points$values[, rg])
ok <- !is.na(v)
v[!ok] <- 0
sums <- rowsum(v, tow_points$tow)
counts <- rowsum(ok * 1, tow_points$tow)
tow_x <- matrix(NA_real_, nrow(tows), length(rg), dimnames = list(NULL, rg))
tow_x[as.integer(rownames(sums)), ] <- ifelse(counts > 0, sums / counts, NA)
tow_x <- data.frame(tow_x, Year = tows$year)
log_message(format(sum(!is.na(tows$label)), big.mark = ","), " labelled tows (",
            format(sum(tows$label %in% 1), big.mark = ","), " strong evidence)")

# --- Cull training data (site-years), designated sites, dives ------------------
train_sites <- readRDS(file.path(cfg$exploration_dir, "artifacts", "site_year_all.rds")) |>
  dplyr::filter(!is.na(.data$high_density), is.finite(.data$longitude), is.finite(.data$latitude))
sxy <- sf::st_coordinates(sf::st_transform(sf::st_as_sf(train_sites, coords = c("longitude", "latitude"), crs = 4326), terra::crs(stack)))
train_x <- terra::extract(stack, sxy)[, rg]
train_x$Year <- train_sites$year

sites_sf <- read_designated_sites(cfg)
site_pts <- readRDS(file.path(cfg$dirs$cache, "cull_sites_points.rds"))
pts_by_site <- split(seq_along(site_pts$site), site_pts$site)
site_years <- readRDS(file.path(cfg$dirs$cache, "cull_site_ranks.rds"))[, c("reef_name", "year", "site_row")]
dives <- matched_cull_dives(cfg, sites_sf)
dives <- dives[!is.na(dives$site_row), ]

# --- Leave-reef-out fits: both responses, same folds --------------------------------
reefs <- unique(c(train_sites$reef_name, tows$reef_name, site_years$reef_name))
fold_of <- stats::setNames(sample(rep(seq_len(n_folds), length.out = length(reefs))), reefs)
fit <- function(x, y) xgboost::xgb.train(params = params, nrounds = rounds, verbose = 0,
                                         data = xgboost::xgb.DMatrix(as.matrix(x[, features]), label = y, missing = NA))
score <- function(model, x) stats::predict(model, xgboost::xgb.DMatrix(as.matrix(x[, features]), missing = NA))
# Within reef-year model: habitat layers only, reef-year offset in training.
# The offset is constant within a reef-year, so it is omitted when ranking.
fit_within <- function(x, y, margin) xgboost::xgb.train(params = params, nrounds = rounds, verbose = 0,
  data = xgboost::xgb.DMatrix(as.matrix(x[, rg]), label = y, base_margin = margin, missing = NA))
score_within <- function(model, x) stats::predict(model, xgboost::xgb.DMatrix(as.matrix(x[, rg]), missing = NA), outputmargin = TRUE)

tows$cull_cv <- NA_real_
tows$manta_cv <- NA_real_
tows$within_cv <- NA_real_
site_years$cull_cv <- NA_real_
site_years$manta_cv <- NA_real_
site_years$within_cv <- NA_real_
importance <- list()
for (k in seq_len(n_folds)) {
  cull_train <- fold_of[train_sites$reef_name] != k
  manta_train <- !is.na(tows$label) & fold_of[tows$reef_name] != k
  cull_model <- fit(train_x[cull_train, ], train_sites$high_density[cull_train])
  manta_model <- fit(tow_x[manta_train, ], tows$label[manta_train])
  within_model <- fit_within(tow_x[manta_train, ], tows$label[manta_train], tows$ry_margin[manta_train])
  importance[[k]] <- dplyr::bind_rows(
    dplyr::mutate(as.data.frame(xgboost::xgb.importance(model = cull_model)), model = "cull", fold = k),
    dplyr::mutate(as.data.frame(xgboost::xgb.importance(model = manta_model)), model = "manta", fold = k),
    dplyr::mutate(as.data.frame(xgboost::xgb.importance(model = within_model)), model = "within", fold = k))

  test_tows <- which(fold_of[tows$reef_name] == k)
  tows$cull_cv[test_tows] <- score(cull_model, tow_x[test_tows, ])
  tows$manta_cv[test_tows] <- score(manta_model, tow_x[test_tows, ])
  tows$within_cv[test_tows] <- score_within(within_model, tow_x[test_tows, ])

  test_sy <- which(fold_of[site_years$reef_name] == k)
  p_idx <- pts_by_site[as.character(site_years$site_row[test_sy])]
  n_pts <- lengths(p_idx)
  x <- site_pts$values[unlist(p_idx, use.names = FALSE), rg]
  x$Year <- rep(site_years$year[test_sy], n_pts)
  owner <- rep(seq_along(test_sy), n_pts)
  for (m in c("cull", "manta", "within")) {
    p <- switch(m, cull = score(cull_model, x), manta = score(manta_model, x), within = score_within(within_model, x))
    s <- tapply(p, owner, mean, na.rm = TRUE)
    site_years[[paste0(m, "_cv")]][test_sy[as.integer(names(s))]] <- s
  }
  log_message(sprintf("Fold %d: cull model on %s site-years, manta model on %s tows", k,
                      format(sum(cull_train), big.mark = ","), format(sum(manta_train), big.mark = ",")))
}
# Held-out tow scores, for the new-ground check in 21.
saveRDS(tows[, c("visit_id", "survey_date", "year", "reef_name", "reef_year", "survey", "label", "strong",
                 "mx", "my", "cull_cv", "manta_cv", "within_cv")],
        file.path(cfg$dirs$cache, "manta_response_tow_scores.rds"))
importance <- dplyr::bind_rows(importance) |>
  dplyr::group_by(.data$model, .data$Feature) |>
  dplyr::summarise(gain = mean(.data$Gain), .groups = "drop") |>
  dplyr::arrange(.data$model, dplyr::desc(.data$gain))
write_table(importance, cfg, "manta_response_importance.csv")

# --- (a) Held-out tows: strong vs negative within reef-years -----------------------
group_pieces <- function(score, y, group) {
  dplyr::bind_rows(lapply(split(seq_along(score), group), function(i) {
    yi <- y[i]
    n1 <- sum(yi)
    n0 <- sum(!yi)
    if (n1 == 0 || n0 == 0) return(NULL)
    r <- rank(score[i])
    top <- rank(-score[i], ties.method = "random") <= ceiling(0.2 * length(i))
    data.frame(auc = (sum(r[yi]) - n1 * (n1 + 1) / 2) / (n1 * n0), pairs = n1 * n0,
               hits = sum(yi & top), positives = n1)
  }))
}
summarise_pieces <- function(p) c(auc = stats::weighted.mean(p$auc, p$pairs), lift = sum(p$hits) / sum(p$positives) / 0.2)
labelled <- tows[!is.na(tows$label) & is.finite(tows$cull_cv) & is.finite(tows$manta_cv) & is.finite(tows$within_cv), ]
tow_check <- dplyr::bind_rows(lapply(c(cull = "cull_cv", manta = "manta_cv", within = "within_cv"), function(col) {
  pieces <- group_pieces(labelled[[col]], labelled$label == 1, labelled$reef_year)
  point <- summarise_pieces(pieces)
  boot <- replicate(boot_reps, summarise_pieces(pieces[sample.int(nrow(pieces), replace = TRUE), ]))
  tibble::tibble(model = sub("_cv$", "", col), metric = names(point), value = point,
                 low = apply(boot, 1, stats::quantile, 0.025), high = apply(boot, 1, stats::quantile, 0.975),
                 tows = nrow(labelled), strong_tows = sum(labelled$label == 1), reef_years = nrow(pieces))
}))
write_table(tow_check, cfg, "manta_response_tow_check.csv")

# --- (a) Thinning tows: aggregations found (as in 14) -------------------------------
agg_tows <- tows |>
  dplyr::filter(is.finite(.data$cull_cv), is.finite(.data$manta_cv), is.finite(.data$within_cv)) |>
  dplyr::group_by(.data$survey) |>
  dplyr::filter(dplyr::n() >= 5, any(.data$strong)) |>
  dplyr::mutate(order = rank(atan2(.data$my - mean(.data$my), .data$mx - mean(.data$mx)), ties.method = "first")) |>
  dplyr::ungroup()
cluster_strong <- function(d) {
  out <- rep(NA_integer_, nrow(d))
  s <- which(d$strong)
  out[s] <- if (length(s) == 1) 1L else stats::cutree(stats::hclust(stats::dist(cbind(d$mx[s], d$my[s])), "single"), h = link_distance_m)
  out
}
agg_tows <- agg_tows |>
  dplyr::group_by(.data$survey) |>
  dplyr::group_modify(~ dplyr::mutate(.x, aggregation = cluster_strong(.x))) |>
  dplyr::ungroup()
strategies <- c("cull_cv", "manta_cv", "within_cv", "random", "systematic")
survey_outcomes <- function(d, k) {
  n <- nrow(d)
  kept_n <- ceiling(k * n)
  agg_ids <- sort(unique(stats::na.omit(d$aggregation)))
  sizes <- vapply(agg_ids, function(a) sum(d$aggregation %in% a), numeric(1))
  found <- function(kept) vapply(agg_ids, function(a) as.numeric(any(kept & d$aggregation %in% a)), numeric(1))
  random <- 1 - exp(lchoose(n - sizes, kept_n) - lchoose(n, kept_n))
  random[kept_n > n - sizes] <- 1
  phases <- seq(0, 1, length.out = 21)[-21]
  systematic <- rowMeans(matrix(vapply(phases, function(phase) {
    found(((d$order - 1 + phase * n / kept_n) %% (n / kept_n)) < 1)
  }, numeric(length(agg_ids))), nrow = length(agg_ids)))
  data.frame(survey = d$survey[1], year = d$year[1], keep = k, size = sizes,
             cull_cv = found(rank(-d$cull_cv, ties.method = "first") <= kept_n),
             manta_cv = found(rank(-d$manta_cv, ties.method = "first") <= kept_n),
             within_cv = found(rank(-d$within_cv, ties.method = "first") <= kept_n),
             random = random, systematic = systematic)
}
surveys <- split(agg_tows, agg_tows$survey)
agg_out <- dplyr::bind_rows(lapply(keep_fractions, function(k) dplyr::bind_rows(lapply(surveys, survey_outcomes, k = k))))
agg_out <- agg_out[agg_out$size >= 2, ]
per_survey <- agg_out |>
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
aggregations <- dplyr::bind_rows(lapply(strategies, function(st) tibble::tibble(
  strategy = st, keep = ks, found = point[, st],
  low = apply(boot[, st, ], 1, stats::quantile, 0.025), high = apply(boot[, st, ], 1, stats::quantile, 0.975),
  aggregations = sum(n_m[, 1]), surveys = length(ids))))
gain <- boot[, "within_cv", ] - boot[, "cull_cv", ]
aggregations <- dplyr::left_join(aggregations, tibble::tibble(
  keep = ks, within_minus_cull = point[, "within_cv"] - point[, "cull_cv"],
  within_minus_cull_low = apply(gain, 1, stats::quantile, 0.025),
  within_minus_cull_high = apply(gain, 1, stats::quantile, 0.975)), by = "keep")
write_table(aggregations, cfg, "manta_response_aggregations.csv")

# --- (b) Designated cull sites: ranks, with and without last season's evidence ------
pct <- function(v) {
  ok <- is.finite(v)
  out <- rep(NA_real_, length(v))
  out[ok] <- (rank(v[ok], ties.method = "random") - 1) / max(sum(ok) - 1, 1)
  out
}
site_years <- site_years |>
  dplyr::group_by(.data$reef_name, .data$year) |>
  dplyr::mutate(rank_cull = pct(.data$cull_cv), rank_manta = pct(.data$manta_cv), rank_within = pct(.data$within_cv),
                rank_ensemble = pct((.data$rank_cull + .data$rank_within) / 2)) |>
  dplyr::ungroup()
saveRDS(site_years, file.path(cfg$dirs$cache, "manta_response_site_ranks.rds"))

evidence_path <- file.path(cfg$dirs$cache, "cull_site_evidence_scores.rds")
if (!file.exists(evidence_path)) stop("Run 16_cull_search_area_evidence.R first (needs ", evidence_path, ")")
ev <- readRDS(evidence_path) |>
  dplyr::filter(.data$window == "last_season") |>
  dplyr::select("reef_name", "year", "site_row", "site_problematic_prev", "site_quiet_prev",
                "cull_rate", "cull_effort", "manta_strong", "manta_effort", "rank_evidence") |>
  dplyr::inner_join(site_years[, c("reef_name", "year", "site_row", "rank_cull", "rank_manta", "rank_within", "rank_ensemble")],
                    by = c("reef_name", "year", "site_row"))
terms_for <- function(f, habitat) {
  data.frame(habitat = ifelse(is.finite(f[[habitat]]), f[[habitat]], 0.5),
             habitat_manta = ifelse(is.finite(f$rank_within), f$rank_within, 0.5),
             site_problematic_prev = f$site_problematic_prev, site_quiet_prev = f$site_quiet_prev,
             cull_rate = f$cull_rate, log_cull_effort = log1p(f$cull_effort),
             log_manta_strong = log1p(f$manta_strong), log_manta_effort = log1p(f$manta_effort))
}
combined_formula <- problematic ~ habitat + site_problematic_prev + site_quiet_prev + cull_rate +
  log_cull_effort + log_manta_strong + log_manta_effort
key <- c("reef_name", "year", "site_row")
dive_cols <- c(key, "outcome", "total_cots", "bottom_time")
habitat_ranks <- c(cull = "rank_cull", manta = "rank_manta", within = "rank_within",
                   ensemble = "rank_ensemble", both = "rank_cull")
both_formula <- stats::update(combined_formula, . ~ . + habitat_manta)
for (h in names(habitat_ranks)) {
  terms <- terms_for(ev, habitat_ranks[[h]])
  formula <- if (h == "both") both_formula else combined_formula
  train_all <- dives[, dive_cols] |>
    dplyr::inner_join(dplyr::bind_cols(ev[, key], terms), by = key) |>
    dplyr::mutate(problematic = as.numeric(.data$outcome == "problematic"))
  lp <- rep(NA_real_, nrow(ev))
  for (k in seq_len(n_folds)) {
    model <- stats::glm(formula, family = stats::binomial(), data = train_all[fold_of[train_all$reef_name] != k, ])
    test <- fold_of[ev$reef_name] == k
    lp[test] <- stats::predict(model, terms[test, ])
  }
  ev[[paste0("lp_", h)]] <- lp
}
ev <- ev |>
  dplyr::group_by(.data$reef_name, .data$year) |>
  dplyr::mutate(rank_cull_evidence = pct(.data$lp_cull), rank_manta_evidence = pct(.data$lp_manta),
                rank_within_evidence = pct(.data$lp_within), rank_ensemble_evidence = pct(.data$lp_ensemble),
                rank_both_evidence = pct(.data$lp_both)) |>
  dplyr::ungroup()

score_cols <- c(cull = "rank_cull", manta = "rank_manta", within = "rank_within", evidence = "rank_evidence",
                cull_evidence = "rank_cull_evidence", manta_evidence = "rank_manta_evidence",
                within_evidence = "rank_within_evidence", ensemble = "rank_ensemble",
                ensemble_evidence = "rank_ensemble_evidence", both_evidence = "rank_both_evidence")
scored <- dives[, dive_cols] |>
  dplyr::inner_join(ev, by = key) |>
  dplyr::filter(dplyr::if_all(dplyr::all_of(unname(score_cols)), is.finite)) |>
  dplyr::mutate(reef_year = paste(.data$reef_name, .data$year), hours = .data$bottom_time / 60)

metrics_for <- function(d, col, a) {
  kept <- d[[col]] >= 1 - a - 1e-9
  c(problematic_dives_kept = sum(kept & d$outcome == "problematic") / sum(d$outcome == "problematic"),
    low_dives_avoided = sum(!kept & d$outcome == "low") / sum(d$outcome == "low"),
    ghost_dives_avoided = sum(!kept & d$outcome == "ghost") / sum(d$outcome == "ghost"),
    non_problematic_hours_avoided = sum(d$hours[!kept & d$outcome != "problematic"]) / sum(d$hours[d$outcome != "problematic"]),
    dives_kept = mean(kept))
}
metric_names <- names(metrics_for(scored[1:2, ], "rank_cull", 1))
grid <- expand.grid(metric = metric_names, sites_searched = site_fractions, score = names(score_cols), stringsAsFactors = FALSE)
evaluate_all <- function(d) unlist(lapply(score_cols, function(col) lapply(site_fractions, function(a) metrics_for(d, col, a))))
matched_metrics <- c("sites_searched", "low_dives_avoided", "ghost_dives_avoided", "non_problematic_hours_avoided")
matched_grid <- expand.grid(metric = matched_metrics, retention = retention_targets, score = names(score_cols), stringsAsFactors = FALSE)
matched_at <- function(v) {
  unlist(lapply(names(score_cols), function(sc) {
    m <- matrix(v[grid$score == sc], nrow = length(metric_names), dimnames = list(metric_names, NULL))
    pk <- c(0, m["problematic_dives_kept", ])
    curves <- list(sites_searched = c(0, site_fractions), low_dives_avoided = c(1, m["low_dives_avoided", ]),
                   ghost_dives_avoided = c(1, m["ghost_dives_avoided", ]),
                   non_problematic_hours_avoided = c(1, m["non_problematic_hours_avoided", ]))
    lapply(retention_targets, function(t) {
      vapply(matched_metrics, function(k) stats::approx(pk, curves[[k]], xout = t, ties = mean)$y, numeric(1))
    })
  }))
}

fixed <- list()
matched <- list()
cull_auc <- list()
for (period in c("2018-2026", "2022-2026")) {
  d <- as.data.frame(scored[scored$year >= as.integer(substr(period, 1, 4)),
                            c("reef_year", "outcome", "hours", unname(score_cols))])
  groups <- split(seq_len(nrow(d)), d$reef_year)
  pt <- evaluate_all(d)
  bt <- vapply(seq_len(boot_reps), function(r) {
    evaluate_all(d[unlist(groups[sample.int(length(groups), replace = TRUE)], use.names = FALSE), ])
  }, numeric(length(pt)))
  fixed[[period]] <- tibble::as_tibble(grid) |>
    dplyr::mutate(period = period, value = pt, low = apply(bt, 1, stats::quantile, 0.025),
                  high = apply(bt, 1, stats::quantile, 0.975), dives = nrow(d), reef_years = length(groups))
  mp <- matched_at(pt)
  mb <- apply(bt, 2, matched_at)
  base <- which(matched_grid$score == "cull")
  matched[[period]] <- dplyr::bind_rows(lapply(names(score_cols), function(sc) {
    idx <- which(matched_grid$score == sc)
    diff <- mb[idx, , drop = FALSE] - mb[base, , drop = FALSE]
    tibble::as_tibble(matched_grid[idx, ]) |>
      dplyr::mutate(period = period, value = mp[idx],
                    low = apply(mb[idx, , drop = FALSE], 1, stats::quantile, 0.025, na.rm = TRUE),
                    high = apply(mb[idx, , drop = FALSE], 1, stats::quantile, 0.975, na.rm = TRUE),
                    gain_over_cull = mp[idx] - mp[base],
                    gain_low = apply(diff, 1, stats::quantile, 0.025, na.rm = TRUE),
                    gain_high = apply(diff, 1, stats::quantile, 0.975, na.rm = TRUE))
  }))
  pl <- d$outcome %in% c("problematic", "low")
  pg <- d$outcome %in% c("problematic", "ghost")
  cull_auc[[period]] <- dplyr::bind_rows(lapply(names(score_cols), function(sc) tibble::tibble(
    period = period, score = sc,
    auc_problematic_vs_low = within_group_auc(d[[score_cols[[sc]]]][pl], d$outcome[pl] == "problematic", d$reef_year[pl]),
    auc_problematic_vs_ghost = within_group_auc(d[[score_cols[[sc]]]][pg], d$outcome[pg] == "problematic", d$reef_year[pg]))))
  log_message(period, ": ", format(nrow(d), big.mark = ","), " dives in ", length(groups), " reef-years")
}
write_table(dplyr::bind_rows(fixed), cfg, "manta_response_cull_cv.csv")
write_table(dplyr::bind_rows(matched), cfg, "manta_response_cull_matched.csv")
write_table(dplyr::bind_rows(cull_auc), cfg, "manta_response_cull_auc.csv")

print(tow_check)
print(aggregations[aggregations$keep %in% c(0.8, 0.6), ], width = Inf)
print(dplyr::bind_rows(cull_auc))
m <- dplyr::bind_rows(matched)
print(m[m$retention %in% c(0.8, 0.9), c("period", "retention", "metric", "score", "value", "gain_over_cull", "gain_low", "gain_high")], n = Inf)
print(importance |> dplyr::group_by(.data$model) |> dplyr::slice_head(n = 6))
