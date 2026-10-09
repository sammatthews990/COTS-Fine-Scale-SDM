# Can the layer plus last season's evidence shrink the area searched on a reef
# without missing problematic dives (CPUE >= 0.02)?
#
# Same question and candidates as 15 (every designated cull site on the reef,
# including never-dived sites), but each site is now ranked using:
#   habitat   the held-out RG + year score from 15 (leave-reef-out XGBoost)
#   evidence  last season's Cull dives and manta tows only (years T-w..T-1;
#             nothing from the season being searched):
#               - dives at the site: any problematic dive, or dived but none
#                 problematic;
#               - Cull dives nearby: Gaussian kernel (sigma 300 m) problematic
#                 rate and dive effort;
#               - manta tows nearby: kernel count of strong-evidence tows and
#                 total tows.
# Scores (each re-ranked among all designated sites on the reef, per year):
#   habitat only              rank from 15
#   evidence only             logistic model of the evidence terms
#   habitat + evidence        logistic model of habitat rank + evidence terms
#   rule                      sites flagged last season first (problematic dive
#                             at the site, or >= 1 strong manta tow nearby),
#                             then by habitat
# The logistic models are fitted on dives from other reefs (5-fold
# leave-reef-out), so no reef is scored using its own dives. Evidence windows:
# last season (w = 1, primary) and last three seasons (w = 3, sensitivity).
# Uncertainty: 500 bootstrap resamples of reef-years, shared across scores so
# differences are paired.
#
# Because evidence moves both problematic and low-density dives up the ranking,
# scores are also compared at matched retention: the share of sites searched,
# and of low-density and ghost dives avoided, when each score keeps the same
# share (80-95 %) of problematic dives (linear interpolation between site
# fractions, inside every bootstrap replicate).

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "..", "R", "common.R"))
cfg <- initialise_final()
source(file.path(cfg$final_dir, "R", "forward.R"))
set.seed(cfg$hmsc$seed + 16L)
n_folds <- 5L
boot_reps <- 500L
site_fractions <- round(seq(0.05, 1, 0.05), 2)
retention_targets <- c(0.8, 0.85, 0.9, 0.95)
kernel_sigma_m <- 300
evidence_windows <- c(last_season = 1L, last_3_seasons = 3L)
strong_classes <- c("direct_cots", "multiple_counted_scars", "legacy_common_scars")
threshold <- cfg$legacy$cpue_threshold
clean_name <- function(x) toupper(gsub("^[[:space:](]+|[[:space:]]+$", "", x))

ranks_path <- file.path(cfg$dirs$cache, "cull_site_ranks.rds")
if (!file.exists(ranks_path)) stop("Run 15_cull_search_area_cv.R first (needs ", ranks_path, ")")
ranks <- readRDS(ranks_path)

# --- Designated sites (same rows as 15) ------------------------------------------
sites_sf <- sf::st_read(file.path(cfg$dirs$cache, "cull_sites_designated.gpkg"), quiet = TRUE)
sites_sf$site_key <- paste(sites_sf$reef_name, clean_name(sites_sf$site_name), sep = "||")
sites_sf <- sites_sf[!duplicated(sites_sf$site_key), ]
crs_mga <- sf::st_crs(sites_sf)
site_xy <- sf::st_coordinates(suppressWarnings(sf::st_point_on_surface(sf::st_geometry(sites_sf))))[, 1:2]
sites <- sf::st_drop_geometry(sites_sf)
sites$row <- seq_len(nrow(sites))

# --- Cull dives, matched to sites as in 15 -------------------------------------
dives <- readRDS(file.path(cfg$exploration_dir, "artifacts", "cull_observations_clean.rds")) |>
  dplyr::filter(is.finite(.data$total_cots), is.finite(.data$bottom_time), .data$bottom_time > 0) |>
  dplyr::mutate(cpue = .data$total_cots / .data$bottom_time,
                outcome = dplyr::case_when(.data$total_cots == 0 ~ "ghost",
                                           .data$cpue < threshold ~ "low",
                                           TRUE ~ "problematic"),
                site_key = paste(.data$reef_name, clean_name(.data$site_name), sep = "||"))
dives$site_row <- match(dives$site_key, sites$site_key)
has_xy <- is.finite(dives$longitude) & is.finite(dives$latitude)
dive_xy <- matrix(NA_real_, nrow(dives), 2)
dive_xy[has_xy, ] <- sf::st_coordinates(sf::st_transform(
  sf::st_as_sf(dives[has_xy, ], coords = c("longitude", "latitude"), crs = 4326), crs_mga))
unmatched <- which(is.na(dives$site_row) & has_xy)
if (length(unmatched)) {
  dpts <- sf::st_as_sf(data.frame(x = dive_xy[unmatched, 1], y = dive_xy[unmatched, 2]), coords = c("x", "y"), crs = crs_mga)
  near <- sf::st_nearest_feature(dpts, sites_sf)
  dist <- as.numeric(sf::st_distance(dpts, sites_sf[near, ], by_element = TRUE))
  ok <- dist <= 100 & sites$reef_name[near] == dives$reef_name[unmatched]
  dives$site_row[unmatched[ok]] <- near[ok]
}
# Evidence location: the dive's own point, else its site.
fill <- !has_xy & !is.na(dives$site_row)
dive_xy[fill, ] <- site_xy[dives$site_row[fill], ]
dives$x <- dive_xy[, 1]
dives$y <- dive_xy[, 2]

# --- Manta tows (midpoints, evidence class) --------------------------------------
tows <- readRDS(checked_input(cfg, "manta_tows")) |>
  dplyr::select("visit_id", "year", "reef_name", "manta_evidence_class")
geom <- readRDS(checked_input(cfg, "footprint_links"))$manta |>
  dplyr::group_by(.data$visit_id) |>
  dplyr::summarise(x = (dplyr::first(.data$start_x) + dplyr::first(.data$end_x)) / 2,
                   y = (dplyr::first(.data$start_y) + dplyr::first(.data$end_y)) / 2, .groups = "drop")
tows <- dplyr::inner_join(tows, geom, by = "visit_id") |>
  dplyr::filter(is.finite(.data$x), is.finite(.data$y)) |>
  dplyr::mutate(strong = .data$manta_evidence_class %in% strong_classes)
log_message(format(nrow(tows), big.mark = ","), " manta tows (", sum(tows$strong), " strong evidence), ",
            min(tows$year), "-", max(tows$year))

# --- Prior-season evidence for every designated site on each dived reef-year --------
evaluated <- ranks |> dplyr::distinct(.data$reef_name, .data$year)
evidence_for <- function(w) {
  dplyr::bind_rows(lapply(sort(unique(evaluated$year)), function(T) {
    rows <- ranks$site_row[ranks$year == T]
    q <- site_xy[rows, , drop = FALSE]
    pd <- dives[dives$year >= T - w & dives$year <= T - 1 & is.finite(dives$x), ]
    pt <- tows[tows$year >= T - w & tows$year <= T - 1, ]
    ck <- history_kernel(cbind(pd$x, pd$y), as.numeric(pd$outcome == "problematic"), rep(0, nrow(pd)),
                         q, kernel_sigma_m, Inf)
    mk <- history_kernel(cbind(pt$x, pt$y), as.numeric(pt$strong), rep(0, nrow(pt)), q, kernel_sigma_m, Inf)
    at_site <- pd[!is.na(pd$site_row), ] |>
      dplyr::group_by(.data$site_row) |>
      dplyr::summarise(site_problematic = any(.data$outcome == "problematic"), .groups = "drop")
    out <- data.frame(year = T, site_row = rows, reef_name = ranks$reef_name[ranks$year == T],
                      cull_rate = ifelse(ck$W > 0, ck$H, 0), cull_effort = ck$W,
                      manta_strong = ifelse(mk$W > 0, mk$H * mk$W, 0), manta_effort = mk$W)
    m <- match(out$site_row, at_site$site_row)
    out$site_problematic_prev <- as.numeric(!is.na(m) & at_site$site_problematic[m] %in% TRUE)
    out$site_quiet_prev <- as.numeric(!is.na(m) & at_site$site_problematic[m] %in% FALSE)
    # Was the reef itself visited in the window (any dives or tows on it)?
    out$reef_cull_prev <- out$reef_name %in% pd$reef_name
    out$reef_manta_prev <- out$reef_name %in% pt$reef_name
    out
  }))
}

pct <- function(v) {
  ok <- is.finite(v)
  out <- rep(NA_real_, length(v))
  out[ok] <- (rank(v[ok], ties.method = "random") - 1) / max(sum(ok) - 1, 1)
  out
}
model_terms <- function(f) {
  data.frame(habitat = ifelse(is.finite(f$rank_cv), f$rank_cv, 0.5),
             site_problematic_prev = f$site_problematic_prev, site_quiet_prev = f$site_quiet_prev,
             cull_rate = f$cull_rate, log_cull_effort = log1p(f$cull_effort),
             log_manta_strong = log1p(f$manta_strong), log_manta_effort = log1p(f$manta_effort))
}
evidence_formula <- problematic ~ site_problematic_prev + site_quiet_prev + cull_rate + log_cull_effort +
  log_manta_strong + log_manta_effort
combined_formula <- stats::update(evidence_formula, . ~ . + habitat)

reefs <- unique(ranks$reef_name)
fold_of <- stats::setNames(sample(rep(seq_len(n_folds), length.out = length(reefs))), reefs)
dive_cols <- c("reef_name", "year", "site_row", "outcome", "total_cots", "bottom_time")

site_scores <- list()
coefficients <- list()
coverage_rows <- list()
for (wname in names(evidence_windows)) {
  w <- evidence_windows[[wname]]
  feat <- ranks |>
    dplyr::select("reef_name", "year", "site_row", "rank_cv") |>
    dplyr::inner_join(evidence_for(w), by = c("reef_name", "year", "site_row"))
  terms <- model_terms(feat)
  train_all <- dives[!is.na(dives$site_row), dive_cols] |>
    dplyr::inner_join(dplyr::bind_cols(feat[, c("reef_name", "year", "site_row")], terms),
                      by = c("reef_name", "year", "site_row")) |>
    dplyr::mutate(problematic = as.numeric(.data$outcome == "problematic"))
  feat$evidence_lp <- NA_real_
  feat$combined_lp <- NA_real_
  for (k in seq_len(n_folds)) {
    train <- train_all[fold_of[train_all$reef_name] != k, ]
    test <- fold_of[feat$reef_name] == k
    ev <- stats::glm(evidence_formula, family = stats::binomial(), data = train)
    co <- stats::glm(combined_formula, family = stats::binomial(), data = train)
    feat$evidence_lp[test] <- stats::predict(ev, terms[test, ])
    feat$combined_lp[test] <- stats::predict(co, terms[test, ])
  }
  full <- stats::glm(combined_formula, family = stats::binomial(), data = train_all)
  coefficients[[wname]] <- tibble::tibble(window = wname, term = names(stats::coef(full)),
                                          estimate = unname(stats::coef(full)),
                                          std_error = unname(sqrt(diag(stats::vcov(full)))))
  feat$flagged <- feat$site_problematic_prev == 1 | feat$manta_strong >= 0.5
  feat$rule <- as.numeric(feat$flagged) + ifelse(is.finite(feat$rank_cv), feat$rank_cv, 0.5) * 0.999
  feat <- feat |>
    dplyr::group_by(.data$reef_name, .data$year) |>
    dplyr::mutate(rank_evidence = pct(.data$evidence_lp), rank_combined = pct(.data$combined_lp),
                  rank_rule = pct(.data$rule)) |>
    dplyr::ungroup()
  feat$window <- wname
  site_scores[[wname]] <- feat
  coverage_rows[[wname]] <- feat |>
    dplyr::distinct(.data$reef_name, .data$year, .data$reef_cull_prev, .data$reef_manta_prev) |>
    dplyr::summarise(window = wname, reef_years = dplyr::n(),
                     share_reef_dived_prev = mean(.data$reef_cull_prev),
                     share_reef_towed_prev = mean(.data$reef_manta_prev),
                     share_reef_any_prev = mean(.data$reef_cull_prev | .data$reef_manta_prev))
  log_message(wname, ": ", format(nrow(feat), big.mark = ","), " site-years scored; flagged share ",
              sprintf("%.1f%%", 100 * mean(feat$flagged)))
}
site_scores <- dplyr::bind_rows(site_scores)
saveRDS(site_scores, file.path(cfg$dirs$cache, "cull_site_evidence_scores.rds"))

# --- Dives scored under every ranking ----------------------------------------------
score_cols <- c(habitat = "rank_cv", evidence = "rank_evidence", combined = "rank_combined", rule = "rank_rule")
scored <- dives[!is.na(dives$site_row), dive_cols] |>
  dplyr::inner_join(site_scores, by = c("reef_name", "year", "site_row"), relationship = "many-to-many") |>
  dplyr::filter(is.finite(.data$rank_cv)) |>
  dplyr::mutate(reef_year = paste(.data$reef_name, .data$year), hours = .data$bottom_time / 60,
                any_prev = .data$reef_cull_prev | .data$reef_manta_prev,
                at_flagged = .data$flagged)

metrics_for <- function(d, rank_col, a) {
  kept <- d[[rank_col]] >= 1 - a - 1e-9
  c(problematic_dives_kept = sum(kept & d$outcome == "problematic") / sum(d$outcome == "problematic"),
    problematic_cots_kept = sum(d$total_cots[kept & d$outcome == "problematic"]) / sum(d$total_cots[d$outcome == "problematic"]),
    low_dives_avoided = sum(!kept & d$outcome == "low") / sum(d$outcome == "low"),
    ghost_dives_avoided = sum(!kept & d$outcome == "ghost") / sum(d$outcome == "ghost"),
    non_problematic_hours_avoided = sum(d$hours[!kept & d$outcome != "problematic"]) / sum(d$hours[d$outcome != "problematic"]),
    problematic_share_of_kept = sum(kept & d$outcome == "problematic") / max(sum(kept), 1),
    dives_kept = mean(kept))
}
period_sets <- function(d) list(`2018-2026` = d$year >= 2018, `2022-2026` = d$year >= 2022)

# Every score x fraction x metric for one set of dives, as one long vector.
evaluate_all <- function(d) {
  unlist(lapply(score_cols, function(col) lapply(site_fractions, function(a) metrics_for(d, col, a))))
}
metric_names <- names(metrics_for(scored[1:2, ], "rank_cv", 1))
grid <- expand.grid(metric = metric_names, sites_searched = site_fractions, score = names(score_cols),
                    stringsAsFactors = FALSE)
# At each retention target, interpolate the other metrics along a score's curve.
matched_metrics <- c("sites_searched", "dives_kept", "low_dives_avoided", "ghost_dives_avoided",
                     "non_problematic_hours_avoided")
matched_grid <- expand.grid(metric = matched_metrics, retention = retention_targets, score = names(score_cols),
                            stringsAsFactors = FALSE)
matched_at <- function(v) {
  unlist(lapply(names(score_cols), function(sc) {
    m <- matrix(v[grid$score == sc], nrow = length(metric_names), dimnames = list(metric_names, NULL))
    pk <- c(0, m["problematic_dives_kept", ])
    curves <- list(sites_searched = c(0, site_fractions), dives_kept = c(0, m["dives_kept", ]),
                   low_dives_avoided = c(1, m["low_dives_avoided", ]), ghost_dives_avoided = c(1, m["ghost_dives_avoided", ]),
                   non_problematic_hours_avoided = c(1, m["non_problematic_hours_avoided", ]))
    lapply(retention_targets, function(t) {
      vapply(matched_metrics, function(k) stats::approx(pk, curves[[k]], xout = t, ties = mean)$y, numeric(1))
    })
  }))
}
results <- list()
paired <- list()
matched <- list()
for (wname in names(evidence_windows)) {
  sw <- scored[scored$window == wname, ]
  sets <- period_sets(sw)
  for (period in names(sets)) {
    d <- as.data.frame(sw[sets[[period]], c("reef_year", "outcome", "total_cots", "hours", score_cols)])
    groups <- split(seq_len(nrow(d)), d$reef_year)
    point <- evaluate_all(d)
    boot <- vapply(seq_len(boot_reps), function(r) {
      evaluate_all(d[unlist(groups[sample.int(length(groups), replace = TRUE)], use.names = FALSE), ])
    }, numeric(length(point)))
    res <- tibble::as_tibble(grid) |>
      dplyr::mutate(window = wname, period = period, value = point,
                    low = apply(boot, 1, stats::quantile, 0.025), high = apply(boot, 1, stats::quantile, 0.975),
                    dives = nrow(d), reef_years = length(groups),
                    problematic_dives = sum(d$outcome == "problematic"), low_dives = sum(d$outcome == "low"),
                    ghost_dives = sum(d$outcome == "ghost"))
    results[[length(results) + 1L]] <- res
    base <- which(grid$score == "habitat")
    for (s in setdiff(names(score_cols), "habitat")) {
      idx <- which(grid$score == s)
      diff <- boot[idx, , drop = FALSE] - boot[base, , drop = FALSE]
      paired[[length(paired) + 1L]] <- tibble::as_tibble(grid[idx, ]) |>
        dplyr::mutate(window = wname, period = period, gain_over_habitat = point[idx] - point[base],
                      low = apply(diff, 1, stats::quantile, 0.025), high = apply(diff, 1, stats::quantile, 0.975))
    }
    mp <- matched_at(point)
    mb <- apply(boot, 2, matched_at)
    hb <- which(matched_grid$score == "habitat")
    matched[[length(matched) + 1L]] <- dplyr::bind_rows(lapply(names(score_cols), function(sc) {
      idx <- which(matched_grid$score == sc)
      diff <- mb[idx, , drop = FALSE] - mb[hb, , drop = FALSE]
      tibble::as_tibble(matched_grid[idx, ]) |>
        dplyr::mutate(window = wname, period = period, value = mp[idx],
                      low = apply(mb[idx, , drop = FALSE], 1, stats::quantile, 0.025, na.rm = TRUE),
                      high = apply(mb[idx, , drop = FALSE], 1, stats::quantile, 0.975, na.rm = TRUE),
                      gain_over_habitat = mp[idx] - mp[hb],
                      gain_low = apply(diff, 1, stats::quantile, 0.025, na.rm = TRUE),
                      gain_high = apply(diff, 1, stats::quantile, 0.975, na.rm = TRUE))
    }))
    log_message(wname, " ", period, ": ", format(nrow(d), big.mark = ","), " dives in ", length(groups), " reef-years")
  }
}
matched <- dplyr::bind_rows(matched)
write_table(matched, cfg, "cull_site_evidence_matched.csv")
results <- dplyr::bind_rows(results)
paired <- dplyr::bind_rows(paired)
write_table(results, cfg, "cull_site_evidence_cv.csv")
write_table(paired, cfg, "cull_site_evidence_gain.csv")

# Within reef-year discrimination, overall and split by whether the reef was
# visited in the evidence window.
auc <- dplyr::bind_rows(lapply(names(evidence_windows), function(wname) {
  sw <- scored[scored$window == wname & scored$year >= 2018, ]
  subsets <- list(all = rep(TRUE, nrow(sw)), `reef visited last season(s)` = sw$any_prev,
                  `reef not visited` = !sw$any_prev)
  dplyr::bind_rows(lapply(names(subsets), function(sub) {
    d <- sw[subsets[[sub]], ]
    pl <- d$outcome %in% c("problematic", "low")
    pg <- d$outcome %in% c("problematic", "ghost")
    dplyr::bind_rows(lapply(names(score_cols), function(s) {
      col <- score_cols[[s]]
      tibble::tibble(window = wname, subset = sub, score = s, dives = nrow(d),
                     reef_years = dplyr::n_distinct(d$reef_year),
                     auc_problematic_vs_low = within_group_auc(d[[col]][pl], d$outcome[pl] == "problematic", d$reef_year[pl]),
                     auc_problematic_vs_ghost = within_group_auc(d[[col]][pg], d$outcome[pg] == "problematic", d$reef_year[pg]))
    }))
  }))
}))
write_table(auc, cfg, "cull_site_evidence_auc.csv")

# How often last season's evidence exists, and what flagged sites look like.
flag_summary <- scored |>
  dplyr::filter(.data$year >= 2018) |>
  dplyr::group_by(.data$window, .data$at_flagged) |>
  dplyr::summarise(dives = dplyr::n(), problematic_share = mean(.data$outcome == "problematic"),
                   low_share = mean(.data$outcome == "low"), ghost_share = mean(.data$outcome == "ghost"),
                   .groups = "drop")
coverage <- dplyr::bind_rows(coverage_rows) |>
  dplyr::left_join(site_scores |> dplyr::filter(.data$year >= 2018) |>
                     dplyr::group_by(.data$window) |>
                     dplyr::summarise(site_years = dplyr::n(), share_sites_flagged = mean(.data$flagged)),
                   by = "window")
write_table(flag_summary, cfg, "cull_site_evidence_flagged.csv")
write_table(coverage, cfg, "cull_site_evidence_coverage.csv")
write_table(dplyr::bind_rows(coefficients), cfg, "cull_site_evidence_coefficients.csv")

print(coverage); print(flag_summary); print(dplyr::bind_rows(coefficients))
print(auc, n = Inf)
print(tidyr::pivot_wider(
  results[results$window == "last_season" & results$period == "2018-2026" &
            results$metric %in% c("problematic_dives_kept", "low_dives_avoided", "ghost_dives_avoided") &
            results$sites_searched %in% c(0.5, 0.6, 0.7, 0.8), c("score", "sites_searched", "metric", "value")],
  names_from = "metric", values_from = "value"), n = Inf)
print(matched[matched$window == "last_season" & matched$retention %in% c(0.8, 0.9) &
                matched$metric %in% c("sites_searched", "low_dives_avoided", "ghost_dives_avoided"), ], n = Inf, width = Inf)
