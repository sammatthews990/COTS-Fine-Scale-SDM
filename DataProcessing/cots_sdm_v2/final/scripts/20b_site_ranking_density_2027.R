# Density-aware priority classes for the 2027 cull-site ranking.
#
# Question: the classes in 20 are shares of each reef's designated sites, so
# every reef gets 20 % A sites however many COTS it holds. If reefs with fewer
# COTS get fewer A-C sites (and dense reefs more), are more > 0.04 dives kept
# for the same total number of sites selected?
#
# Method:
#   reef density  the reef's previous-season evidence: share of problematic
#                 dives and of strong-evidence manta tows (each shrunk by half
#                 a dive or tow), log dive count, indicators for no dives / no
#                 tows, and the GBR-wide reef outbreak forecast (Stewart et al.
#                 hybrid model) for the previous year. A reef with no recent
#                 dives or tows is placed by the forecast alone, and its P is
#                 scaled by the chance such a reef needed culling at all
#                 (cull-site reefs towed with no evidence the season before:
#                 about 1 in 5 had a problematic dive; logistic on the
#                 forecast). Without that, the model would give unknown reefs
#                 the level of reefs crews chose to dive. A reef with no
#                 forecast either gets the median forecast (habitat only: its
#                 within-reef order at a low default density).
#   site          one logistic model of "a dive here is problematic" on the
#                 reef terms plus the site's within-reef percentile in the 20
#                 ranking (19's forward percentiles in the test). The reef terms
#                 are constant within a reef, so the order on a reef is
#                 unchanged; only how many sites reach each class moves.
#   classes       scaled quotas: each reef keeps its own order, and its class
#                 shares shift with its density on the odds scale, so a
#                 typical culled reef (median over the reef-years culled in
#                 cfg$site_ranking$density_calibration_years) keeps A 20, B 50,
#                 C 66, D 80 %; Batt-like density gives about 40 % A, a sparse
#                 reef few A, an unknown low-forecast reef almost none.
#                 Absolute cut-offs on P were tested too: they let density
#                 swamp the within-reef order (dense reefs all A, sparse and
#                 unknown reefs none), leaving crews no order on most reefs.
#
# Validation, strictly forward (test seasons 2020-2026; 19's forward
# percentiles start in 2019): for season T the model is fitted on seasons
# 2019..T-1 and every designated site on each reef dived in T is scored.
#   equal effort    the same total number of sites as within-reef A-C (66 %
#                   of every reef): > 0.04 and problematic dives kept, low and
#                   ghost dives avoided
#   matched         sites selected to keep 80 % and 90 % of > 0.04 dives
#   deployment      the class rule as deployed (absolute: cut-offs calibrated
#                   on the training seasons; scaled: fixed), A-C: sites
#                   selected and dives kept
# Thresholds for equal effort and matching are set per season. Variants:
# absolute P from evidence + forecast, evidence only, forecast only (unknown
# reefs unscaled), "absolute" (evidence + forecast, unknown reefs scaled) and
# "scaled" (the 2027 classes). Unknown reef-years in the test were all dived,
# usually after tows found COTS, so scaling them down is penalised here; in
# use, the first tows of a visit re-rank the reef (33). Baselines: the
# within-reef classes (20) and random selection (exact expectation: keeps the
# share of sites selected). Intervals: bootstrap of reef-years.
# Caveats: forecasts before 2026 are hindcasts, probably fitted on these years.
# Dive counts and past dives encode where crews chose to dive; the test can
# only score reefs that were dived, so it says nothing about reefs never
# culled, where the forecast alone sets the classes.
#
# Outputs (outputs/site_ranking/):
#   cots_cull_site_priority_density_2027.csv        one row per designated site
#   cots_cull_site_priority_density_2027_reefs.csv  per-reef density and classes
#   cots_cull_site_priority_density_2027.kmz        EoTR KMZ coloured by these
#                                                   classes (21 adds exploration)

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "..", "R", "common.R"))
cfg <- initialise_final()
source(file.path(cfg$final_dir, "R", "cull_sites.R"))
source(file.path(cfg$final_dir, "R", "site_selection.R"))
sr <- cfg$site_ranking
set.seed(cfg$hmsc$seed + 200L)
strong_classes <- c("direct_cots", "multiple_counted_scars", "legacy_common_scars")
key <- c("reef_name", "year", "site_row")
test_years <- 2020:2026
boot_reps <- 200L
high_cpue <- 0.04
fractions <- round(seq(0.05, 1, 0.05), 2)
breaks <- sr$priority_breaks
labels <- c(names(breaks), LETTERS[length(breaks) + 1L])
ac_share <- breaks[["C"]]
out_dir <- file.path(cfg$dirs$outputs, "site_ranking")

# --- Inputs ---------------------------------------------------------------------------
reef_id <- function(x) sub(".*[(]([0-9]+-[0-9]+[a-z]?)[)][[:space:]]*$", "\\1", x)
forecast <- dplyr::bind_rows(utils::read.csv(checked_input(cfg, "reef_forecast_2026")),
                             utils::read.csv(checked_input(cfg, "reef_forecast_hindcast"))) |>
  dplyr::transmute(reef_id = reef_id(.data$reefName), year = as.integer(.data$year), forecast = .data$outbrProb) |>
  dplyr::distinct(.data$reef_id, .data$year, .keep_all = TRUE)   # the 2026 forecast wins any overlap

sites_sf <- read_designated_sites(cfg)
site_counts <- dplyr::count(sf::st_drop_geometry(sites_sf), .data$reef_name, name = "sites_on_reef")
dives <- matched_cull_dives(cfg, sites_sf)
dives$high <- dives$cpue > high_cpue
tows <- readRDS(checked_input(cfg, "manta_tows")) |>
  dplyr::select("visit_id", "year", "survey_date", "reef_name", "manta_evidence_class") |>
  dplyr::mutate(strong = .data$manta_evidence_class %in% strong_classes)

# Reef-level evidence from dives `d` and tows `t` (one window), with the forecast.
reef_evidence <- function(reefs, d, t, forecast_year) {
  dv <- d |>
    dplyr::group_by(.data$reef_name) |>
    dplyr::summarise(dives = dplyr::n(), problematic = sum(.data$outcome == "problematic"), .groups = "drop")
  tw <- t |>
    dplyr::group_by(.data$reef_name) |>
    dplyr::summarise(tows = dplyr::n(), strong = sum(.data$strong), .groups = "drop")
  fc <- forecast[forecast$year == forecast_year, ]
  tibble::tibble(reef_name = reefs) |>
    dplyr::left_join(dv, by = "reef_name") |>
    dplyr::left_join(tw, by = "reef_name") |>
    dplyr::mutate(dplyr::across(c("dives", "problematic", "tows", "strong"), ~ dplyr::coalesce(.x, 0L)),
                  forecast_prob = fc$forecast[match(reef_id(.data$reef_name), fc$reef_id)])
}
reef_terms <- function(e, forecast_fill) {
  logit_share <- function(k, n) ifelse(n > 0, stats::qlogis((k + 0.5) / (n + 1)), 0)
  f <- dplyr::coalesce(e$forecast_prob, forecast_fill)
  data.frame(has_dives = as.numeric(e$dives > 0), dive_rate = logit_share(e$problematic, e$dives),
             log_dives = log1p(e$dives), has_tows = as.numeric(e$tows > 0),
             tow_rate = logit_share(e$strong, e$tows), forecast = stats::qlogis(pmin(pmax(f, 1e-3), 1 - 1e-3)))
}
variants <- list(
  evidence_forecast = problematic ~ pct + has_dives + dive_rate + log_dives + has_tows + tow_rate + forecast,
  evidence = problematic ~ pct + has_dives + dive_rate + log_dives + has_tows + tow_rate,
  forecast = problematic ~ pct + forecast)
forecast_median <- stats::median(forecast$forecast[forecast$year >= 2018], na.rm = TRUE)

# Previous-season evidence for every reef-year 2019-2026 (calendar T - 1, as in
# 16 and 19), and the within-reef percentile of every dived site from 19.
fw <- readRDS(file.path(cfg$dirs$cache, "site_ranking_forward.rds")) |>
  dplyr::distinct(dplyr::across(dplyr::all_of(key)), .data$cull_all) |>
  dplyr::rename(pct = "cull_all")
dived <- dives[!is.na(dives$site_row), c(key, "outcome", "high")] |>
  dplyr::inner_join(fw, by = key)
years <- sort(unique(dived$year))
reef_years <- dplyr::bind_rows(lapply(years, function(T) {
  reefs <- unique(dived$reef_name[dived$year == T])
  e <- reef_evidence(reefs, dives[dives$year == T - 1, ], tows[tows$year == T - 1, ], T - sr$density_forecast_lag)
  dplyr::bind_cols(e, reef_terms(e, forecast_median)) |> dplyr::mutate(year = T)
})) |>
  dplyr::left_join(site_counts, by = "reef_name") |>
  dplyr::filter(!is.na(.data$sites_on_reef)) |>
  dplyr::mutate(reef_year = paste(.data$reef_name, .data$year))
dived <- dplyr::inner_join(dived, reef_years, by = c("reef_name", "year")) |>
  dplyr::mutate(problematic = as.numeric(.data$outcome == "problematic"))
log_message(format(nrow(dived), big.mark = ","), " dives on ", nrow(reef_years), " reef-years ", min(years), "-", max(years))

# Every designated site of a reef-year, at its within-reef percentile.
site_grid <- function(ry) {
  n <- ry$sites_on_reef
  idx <- rep(seq_len(nrow(ry)), n)
  out <- ry[idx, ]
  out$pct <- unlist(lapply(n, function(k) (seq_len(k) - 1) / max(k - 1, 1)))
  out
}
# Cut-offs on P giving the cumulative class shares over a set of sites.
class_cutoffs <- function(p, shares) {
  s <- sort(p, decreasing = TRUE)
  vapply(shares, function(a) s[max(1L, ceiling(a * length(s)))], numeric(1))
}

# Unknown reefs (no dives or tows in the previous season). The site model's
# rows for such reefs are reef-years crews dived anyway, usually because tows
# had just found COTS, so it sets their level far too high. Their P is
# scaled by the chance that such a reef needed culling at all: cull-site
# reef-years towed in season T with no dives or tows in T - 1, outcome any
# problematic dive in T, logistic on the logit forecast.
cull_reefs <- unique(sites_sf$reef_name)
surveyed <- dplyr::bind_rows(lapply(sort(unique(tows$year)), function(T) {
  prev <- union(tows$reef_name[tows$year == T - 1], dives$reef_name[dives$year == T - 1])
  r <- setdiff(intersect(unique(tows$reef_name[tows$year == T]), cull_reefs), prev)
  if (!length(r)) return(NULL)
  e <- reef_evidence(r, dives[0, ], tows[0, ], T - sr$density_forecast_lag)
  data.frame(year = T, reef_name = r, forecast = reef_terms(e, forecast_median)$forecast,
             needed = as.numeric(r %in% dives$reef_name[dives$year == T & dives$outcome == "problematic"]))
}))
unknown_fit <- function(max_year) stats::glm(needed ~ forecast, family = stats::binomial(),
                                             data = surveyed[surveyed$year <= max_year, ])
discount <- function(lp, d, um) {
  unknown <- d$has_dives == 0 & d$has_tows == 0
  w <- stats::predict(um, data.frame(forecast = d$forecast[unknown]), type = "response")
  lp[unknown] <- stats::qlogis(w * stats::plogis(lp[unknown]))
  lp
}
log_message(nrow(surveyed), " unknown-reef surveys ", min(surveyed$year), "-", max(surveyed$year), "; ",
            round(100 * mean(surveyed$needed)), " % needed culling")
# Scaled quotas (deployed): absolute cut-offs on P let reef density swamp the
# within-reef order (dense reefs all A, sparse reefs none), so the product
# instead keeps every reef's own order and shifts its class shares with its
# density on the odds scale:
#   share in classes up to k = plogis(qlogis(break_k) + delta),
#   delta = qlogis(reef density) - qlogis(typical culled reef),
# where reef density is the mean P over the reef's sites (evidence + forecast,
# unknown reefs scaled) and the typical reef is the median over culled
# reef-years. A typical reef keeps 20/50/66/80 %. As one score: a site at top
# fraction t of its reef (rank - 0.5) / n scores delta - qlogis(t), and is in
# class k or better when the score >= -qlogis(break_k).
scaled_score <- function(delta, rank, n) delta - stats::qlogis((rank - 0.5) / n)
reef_density <- function(lp, reef) tapply(stats::plogis(lp), reef, mean)
rank_of <- function(pct, n) n - pct * (n - 1)
# Scores compared: each model variant; "absolute" (evidence + forecast,
# unknown reefs scaled, absolute cut-offs on P); "scaled" (deployed).
model_scores <- c(names(variants), "absolute")
score_names <- c(model_scores, "scaled")

# --- Forward test ---------------------------------------------------------------------
scored_sites <- list()
scored_dives <- list()
deploy <- list()
for (T in test_years) {
  tr <- dived[dived$year < T, ]
  te_ry <- reef_years[reef_years$year == T & reef_years$reef_year %in% dived$reef_year, ]
  grid <- site_grid(te_ry)
  dv <- dived[dived$year == T, ]
  tr_grid <- site_grid(reef_years[reef_years$year < T & reef_years$reef_year %in% tr$reef_year, ])
  um <- unknown_fit(T - 1)
  for (vn in model_scores) {
    m <- stats::glm(variants[[if (vn == "absolute") "evidence_forecast" else vn]], family = stats::binomial(), data = tr)
    adjust <- if (vn == "absolute") function(lp, d) discount(lp, d, um) else function(lp, d) lp
    grid[[vn]] <- adjust(stats::predict(m, grid), grid)
    dv[[vn]] <- adjust(stats::predict(m, dv), dv)
    tr_lp <- adjust(stats::predict(m, tr_grid), tr_grid)
    # Deployment: cut-off for A-C calibrated on the training seasons' sites.
    deploy[[length(deploy) + 1L]] <- tibble::tibble(year = T, variant = vn, cutoff = class_cutoffs(tr_lp, ac_share))
  }
  # Scaled quotas from the absolute model's reef densities.
  ref <- stats::median(reef_density(tr_lp, tr_grid$reef_year))
  dens <- reef_density(grid$absolute, grid$reef_year)
  delta <- stats::qlogis(dens) - stats::qlogis(ref)
  grid$scaled <- scaled_score(delta[grid$reef_year], rank_of(grid$pct, grid$sites_on_reef), grid$sites_on_reef)
  dv$scaled <- scaled_score(delta[dv$reef_year], rank_of(dv$pct, dv$sites_on_reef), dv$sites_on_reef)
  deploy[[length(deploy) + 1L]] <- tibble::tibble(year = T, variant = "scaled", cutoff = -stats::qlogis(ac_share))
  scored_sites[[as.character(T)]] <- grid
  scored_dives[[as.character(T)]] <- dv
}
scored_sites <- dplyr::bind_rows(scored_sites)
scored_dives <- dplyr::bind_rows(scored_dives)
deploy <- dplyr::bind_rows(deploy)

kept_metrics <- function(kept, d, prefix) {
  stats::setNames(c(sum(kept & d$high) / sum(d$high),
                    sum(kept & d$outcome == "problematic") / sum(d$outcome == "problematic"),
                    sum(!kept & d$outcome == "low") / sum(d$outcome == "low"),
                    sum(!kept & d$outcome == "ghost") / sum(d$outcome == "ghost")),
                  paste0(prefix, c("high_kept", "problematic_kept", "low_avoided", "ghost_avoided")))
}
# Matched retention: sites selected (share of all designated sites) and low and
# ghost dives avoided where a selection rule keeps a share of > 0.04 dives.
# `kept_at` returns the dives kept when a share `a` of sites is selected.
matched_high <- function(s_share, kept_at, d, targets = c(0.8, 0.9)) {
  curve <- vapply(fractions, function(a) {
    k <- kept_at(a)
    c(s_share(a), sum(d$high & k) / sum(d$high), sum(!k & d$outcome == "low") / sum(d$outcome == "low"),
      sum(!k & d$outcome == "ghost") / sum(d$outcome == "ghost"))
  }, numeric(4))
  unlist(lapply(targets, function(t) {
    v <- vapply(list(c(0, curve[1, ]), c(1, curve[3, ]), c(1, curve[4, ])), function(y)
      stats::approx(c(0, curve[2, ]), y, xout = t, ties = mean)$y, numeric(1))
    stats::setNames(v, paste0(c("sites", "low_avoided", "ghost_avoided"), "_at_", round(100 * t), "_high"))
  }))
}
evaluate <- function(s, d) {
  within_kept <- d$pct >= 1 - ac_share - 1e-9
  budget <- sum(s$pct >= 1 - ac_share - 1e-9)   # within-reef A-C: the equal-effort budget
  out <- c(within = c(kept_metrics(within_kept, d, "ac_"), ac_sites = budget / nrow(s), ac_reef_years_without_sites = 0,
                      matched_high(function(a) mean(s$pct >= 1 - a - 1e-9), function(a) d$pct >= 1 - a - 1e-9, d)))
  # Thresholds are set per season (each season's models differ in calibration,
  # and effort is planned per season): the season's top share `a` of sites.
  ys <- as.character(sort(unique(s$year)))
  by_year <- split(seq_len(nrow(s)), s$year)[ys]
  for (vn in score_names) {
    srt <- lapply(by_year, function(i) sort(s[[vn]][i], decreasing = TRUE))
    thr <- function(a) vapply(srt, function(v) v[max(1L, ceiling(a * length(v)))], numeric(1))
    above <- function(x, yr, t) x >= t[as.character(yr)]
    eq <- vapply(names(by_year), function(y) srt[[y]][sum(s$pct[by_year[[y]]] >= 1 - ac_share - 1e-9)], numeric(1))
    kept <- above(d[[vn]], d$year, eq)
    ry_any <- tapply(above(s[[vn]], s$year, eq), s$reef_year, any)
    cut <- deploy$cutoff[deploy$variant == vn][match(s$year, deploy$year[deploy$variant == vn])]
    dcut <- deploy$cutoff[deploy$variant == vn][match(d$year, deploy$year[deploy$variant == vn])]
    v <- c(kept_metrics(kept, d, "ac_"), ac_sites = mean(above(s[[vn]], s$year, eq)),
           ac_reef_years_without_sites = mean(!ry_any),
           matched_high(function(a) mean(above(s[[vn]], s$year, thr(a))), function(a) above(d[[vn]], d$year, thr(a)), d),
           deployed_ac_sites = mean(s[[vn]] >= cut), kept_metrics(d[[vn]] >= dcut, d, "deployed_ac_"))
    out <- c(out, stats::setNames(v, paste(vn, names(v), sep = ".")))
  }
  out
}
groups_s <- split(seq_len(nrow(scored_sites)), scored_sites$reef_year)
groups_d <- split(seq_len(nrow(scored_dives)), scored_dives$reef_year)
ry_all <- names(groups_d)
point <- evaluate(scored_sites, scored_dives)
boot <- vapply(seq_len(boot_reps), function(r) {
  pick <- sample(ry_all, replace = TRUE)
  evaluate(scored_sites[unlist(groups_s[pick], use.names = FALSE), ], scored_dives[unlist(groups_d[pick], use.names = FALSE), ])
}, numeric(length(point)))
parts <- do.call(rbind, strsplit(names(point), ".", fixed = TRUE))
forward <- tibble::tibble(classes = parts[, 1], metric = parts[, 2], value = point,
                          low = apply(boot, 1, stats::quantile, 0.025, na.rm = TRUE),
                          high = apply(boot, 1, stats::quantile, 0.975, na.rm = TRUE),
                          dives = nrow(scored_dives), reef_years = length(ry_all), seasons = "2020-2026")
# Paired differences from the within-reef classes.
for (vn in score_names) for (mt in c("ac_high_kept", "ac_problematic_kept", "ac_low_avoided", "ac_ghost_avoided",
                                          "sites_at_80_high", "sites_at_90_high", "low_avoided_at_90_high",
                                          "ghost_avoided_at_90_high")) {
  a <- paste(vn, mt, sep = ".")
  b <- paste("within", mt, sep = ".")
  forward <- dplyr::bind_rows(forward, tibble::tibble(
    classes = vn, metric = paste0(mt, "_minus_within"), value = point[[a]] - point[[b]],
    low = stats::quantile(boot[a, ] - boot[b, ], 0.025), high = stats::quantile(boot[a, ] - boot[b, ], 0.975),
    dives = nrow(scored_dives), reef_years = length(ry_all), seasons = "2020-2026"))
}
write_table(forward, cfg, "site_ranking_density_forward.csv")
print(tidyr::pivot_wider(forward[, c("classes", "metric", "value")], names_from = "classes", values_from = "value"), n = Inf)

# Per season (point estimates): does the gain hold every year?
by_season <- dplyr::bind_rows(lapply(test_years, function(T) {
  p <- evaluate(scored_sites[scored_sites$year == T, ], scored_dives[scored_dives$year == T, ])
  parts <- do.call(rbind, strsplit(names(p), ".", fixed = TRUE))
  tibble::tibble(season = T, classes = parts[, 1], metric = parts[, 2], value = p,
                 dives = sum(scored_dives$year == T), high_dives = sum(scored_dives$high[scored_dives$year == T]))
}))
write_table(by_season, cfg, "site_ranking_density_forward_by_season.csv")
print(tidyr::pivot_wider(by_season[by_season$metric %in% c("ac_high_kept", "sites_at_90_high", "ghost_avoided_at_90_high"),
                                   c("season", "metric", "classes", "value")], names_from = "classes", values_from = "value"), n = Inf)

# --- 2027: density-aware classes for every designated site -------------------------
deployed <- "evidence_forecast"
cal <- dived[dived$year >= sr$density_calibration_years[1] & dived$year <= sr$density_calibration_years[2], ]
model <- stats::glm(variants[[deployed]], family = stats::binomial(), data = cal)
um <- unknown_fit(sr$density_calibration_years[2])
saveRDS(list(site = model, unknown_reef = um), file.path(cfg$dirs$models, "site_ranking_density_2027_glm.rds"))
cal_grid <- site_grid(reef_years[reef_years$reef_year %in% cal$reef_year, ])
typical <- stats::median(reef_density(discount(stats::predict(model, cal_grid), cal_grid, um), cal_grid$reef_year))
coef_term <- function(m, prefix) tibble::tibble(term = paste0(prefix, names(stats::coef(m))), estimate = unname(stats::coef(m)),
                                                std_error = unname(sqrt(diag(stats::vcov(m)))))
coef_table <- dplyr::bind_rows(
  coef_term(model, ""), coef_term(um, "unknown_reef_needs_culling:"),
  tibble::tibble(term = "typical_culled_reef_density", estimate = typical)) |>
  dplyr::mutate(training_dives = nrow(cal), training_reef_years = dplyr::n_distinct(cal$reef_year),
                unknown_reef_surveys = sum(surveyed$year <= sr$density_calibration_years[2]),
                training_years = paste(sr$density_calibration_years, collapse = "-"), variant = deployed)
write_table(coef_table, cfg, "site_ranking_density_2027_coefficients.csv")

# Reef evidence in the 12 months used by 20, and the latest forecast (2026).
evidence_start <- sr$evidence_end - sr$evidence_days + 1
ranking <- utils::read.csv(file.path(out_dir, "cots_cull_site_priority_2027.csv"), stringsAsFactors = FALSE)
e27 <- reef_evidence(unique(ranking$reef_name),
                     dives[dives$survey_date >= evidence_start & dives$survey_date <= sr$evidence_end, ],
                     tows[tows$survey_date >= evidence_start & tows$survey_date <= sr$evidence_end, ],
                     sr$target_season - sr$density_forecast_lag)
e27 <- dplyr::bind_cols(e27, reef_terms(e27, forecast_median)) |>
  dplyr::mutate(density_source = dplyr::case_when(
    .data$dives > 0 & .data$tows > 0 ~ "dives and manta tows",
    .data$dives > 0 ~ "dives",
    .data$tows > 0 ~ "manta tows",
    !is.na(.data$forecast_prob) ~ "GBR outbreak forecast only",
    TRUE ~ "none (habitat order at median density)"))
sites <- ranking |>
  dplyr::mutate(pct = (.data$sites_on_reef - .data$rank_on_reef) / pmax(.data$sites_on_reef - 1, 1)) |>
  dplyr::left_join(e27, by = "reef_name")
sites$lp <- discount(stats::predict(model, sites), sites, um)
sites$p_problematic_density <- stats::plogis(sites$lp)
sites <- sites |>
  dplyr::group_by(.data$reef_name) |>
  dplyr::mutate(reef_expected_problematic = mean(.data$p_problematic_density)) |>
  dplyr::ungroup()
sites$score_scaled <- scaled_score(stats::qlogis(sites$reef_expected_problematic) - stats::qlogis(typical),
                                   sites$rank_on_reef, sites$sites_on_reef)
sites$priority_density <- labels[1L + rowSums(outer(sites$score_scaled, -stats::qlogis(breaks), "<"))]

table_out <- sites |>
  dplyr::transmute(reef_name = .data$reef_name, site_name = .data$site_name, cull_site_id = .data$cull_site_id,
                   priority_density = .data$priority_density, priority_within_reef = .data$priority,
                   rank_on_reef = .data$rank_on_reef, sites_on_reef = .data$sites_on_reef,
                   p_problematic_if_dived = round(.data$p_problematic_density, 3),
                   reef_expected_problematic = round(.data$reef_expected_problematic, 3),
                   density_source = .data$density_source, reef_dives_12m = .data$dives,
                   reef_problematic_dives_12m = .data$problematic, reef_tows_12m = .data$tows,
                   reef_strong_tows_12m = .data$strong, reef_outbreak_forecast_2026 = round(.data$forecast_prob, 3),
                   evidence = .data$evidence, score_log_odds_within = .data$score_log_odds)
utils::write.csv(table_out, file.path(out_dir, "cots_cull_site_priority_density_2027.csv"), row.names = FALSE, na = "")
log_message("Wrote ", file.path(out_dir, "cots_cull_site_priority_density_2027.csv"))

count_class <- function(p, cls) sum(p == cls)
reefs <- table_out |>
  dplyr::group_by(.data$reef_name) |>
  dplyr::summarise(sites = dplyr::n(), density_source = dplyr::first(.data$density_source),
                   reef_expected_problematic = dplyr::first(.data$reef_expected_problematic),
                   dives_12m = dplyr::first(.data$reef_dives_12m), problematic_dives_12m = dplyr::first(.data$reef_problematic_dives_12m),
                   tows_12m = dplyr::first(.data$reef_tows_12m), strong_tows_12m = dplyr::first(.data$reef_strong_tows_12m),
                   outbreak_forecast_2026 = dplyr::first(.data$reef_outbreak_forecast_2026),
                   priority_a = count_class(.data$priority_density, "A"), priority_b = count_class(.data$priority_density, "B"),
                   priority_c = count_class(.data$priority_density, "C"), priority_d = count_class(.data$priority_density, "D"),
                   priority_e = count_class(.data$priority_density, "E"),
                   within_reef_a_to_c = sum(.data$priority_within_reef %in% c("A", "B", "C")), .groups = "drop") |>
  dplyr::arrange(dplyr::desc(.data$reef_expected_problematic))
utils::write.csv(reefs, file.path(out_dir, "cots_cull_site_priority_density_2027_reefs.csv"), row.names = FALSE, na = "")

# How the classes moved, by where each reef's density came from.
shift <- table_out |>
  dplyr::group_by(.data$density_source) |>
  dplyr::summarise(reefs = dplyr::n_distinct(.data$reef_name), sites = dplyr::n(),
                   within_a = sum(.data$priority_within_reef == "A"), density_a = sum(.data$priority_density == "A"),
                   within_a_to_c = sum(.data$priority_within_reef %in% c("A", "B", "C")),
                   density_a_to_c = sum(.data$priority_density %in% c("A", "B", "C")),
                   reefs_without_a = sum(tapply(.data$priority_density == "A", .data$reef_name, sum) == 0),
                   reefs_without_a_to_c = sum(tapply(.data$priority_density %in% c("A", "B", "C"), .data$reef_name, sum) == 0),
                   median_reef_expected_problematic = stats::median(tapply(.data$reef_expected_problematic, .data$reef_name, dplyr::first)),
                   .groups = "drop")
shift <- dplyr::bind_rows(shift, shift |>
                            dplyr::summarise(density_source = "all reefs", dplyr::across(-c("density_source", "median_reef_expected_problematic"), sum),
                                             median_reef_expected_problematic = stats::median(.env$reefs$reef_expected_problematic)))
write_table(shift, cfg, "site_ranking_density_2027_classes.csv")

# --- KMZ: the EoTR file coloured by the density-aware classes -------------------------
fmt <- function(x, d) ifelse(is.finite(x), formatC(x, format = "f", digits = d), "")
k <- priority_kml(cfg, table_out, table_out$priority_density,
                  fields = list(`Priority 2027 (density-aware)` = table_out$priority_density,
                                `Within-reef priority` = table_out$priority_within_reef,
                                `Rank on reef` = paste(table_out$rank_on_reef, "of", table_out$sites_on_reef),
                                `P(problematic if dived)` = fmt(table_out$p_problematic_if_dived, 2),
                                `Reef density` = paste0(round(100 * table_out$reef_expected_problematic), " % of dives expected problematic (",
                                                        table_out$density_source, ")"),
                                `Last season` = table_out$evidence),
                  title = paste0("COTS cull site priority ", sr$target_season,
                                 ", density-aware (A red, B orange, C yellow, D pale yellow, E grey)"))
write_staged_kmz(cfg, k$kml, "site_ranking_density_kmz", file.path(out_dir, "cots_cull_site_priority_density_2027.kmz"))

print(coef_table, n = Inf)
print(shift, width = Inf)
print(utils::head(reefs[, c("reef_name", "density_source", "reef_expected_problematic", "sites", "priority_a", "priority_b",
                            "priority_c", "within_reef_a_to_c")], 15), width = Inf)
