# How much surveillance and culling effort could the primary search-planning
# layer save, and at what cost in detections? Each scenario keeps the top
# fraction k of tows (or cull sites) ranked by the layer *within* a survey or
# reef-year, and is compared with keeping the same number chosen at random.
#
# Manta surveillance (unit: one reef on one survey date)
#   score: primary layer (RG + year, clean fallback) block-mean at 120 m,
#          length-weighted along each tow's footprint
#   outcome: strong COTS evidence on the tow (direct COTS, >= 2 counted scars,
#          legacy "common" scars). Manta is not used to train the layer.
#   metrics: tows and km saved; strong-evidence tows kept; outbreak surveys
#          and reef-years still detected (>= 1 strong tow kept)
# Cull dives (unit: one reef in one year; sites ranked within it)
#   score: held-out RG + year prediction for the site (site excluded from the
#          model that scored it), so the estimate is not in-sample
#   metrics: diver-hours saved, ghost dives (0 COTS) avoided, COTS at skipped
#          sites, and the net COTS gain if saved hours are redeployed at the
#          kept sites' catch rate.

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "..", "R", "common.R"))
cfg <- initialise_final()
ops <- cfg$operations
set.seed(cfg$hmsc$seed)
template <- terra::rast(file.path(cfg$dirs$cache, "base_support_120m.tif"))[[1]]
search120 <- aggregate_to_support_120(
  file.path(cfg$dirs$layers, "cots_legacy_rg_year2026_clean_fallback_30m.tif"),
  file.path(cfg$dirs$cache, "search_planning_120m.tif"), template
)
strong_classes <- c("direct_cots", "multiple_counted_scars", "legacy_common_scars")

# --- Manta tows -----------------------------------------------------------------
tows <- readRDS(checked_input(cfg, "manta_tows")) |>
  dplyr::select("visit_id", "survey_date", "year", "reef_key", "reef_name", "reef_year_id", "manta_evidence_class")
links <- readRDS(checked_input(cfg, "footprint_links"))$manta |>
  dplyr::select("visit_id", "x_bin", "y_bin", "weighted_length_m", "fraction", "tow_distance_m") |>
  dplyr::mutate(
    weight = dplyr::case_when(
      is.finite(.data$weighted_length_m) & .data$weighted_length_m > 0 ~ .data$weighted_length_m,
      is.finite(.data$fraction) & .data$fraction > 0 ~ .data$fraction,
      TRUE ~ 1
    ),
    value = terra::extract(search120, terra::cellFromXY(template, cbind((.data$x_bin + 0.5) * 120, (.data$y_bin + 0.5) * 120)))[[1]]
  )
tow_scores <- links |>
  dplyr::group_by(.data$visit_id) |>
  dplyr::summarise(
    coverage = sum(.data$weight[is.finite(.data$value)]) / sum(.data$weight),
    score = if (any(is.finite(.data$value))) stats::weighted.mean(.data$value[is.finite(.data$value)], .data$weight[is.finite(.data$value)]) else NA_real_,
    tow_km = dplyr::first(.data$tow_distance_m) / 1000,
    .groups = "drop"
  )
cull_reefs <- unique(readRDS(file.path(cfg$exploration_dir, "artifacts", "cull_observations_clean.rds"))$reef_name)
manta <- tows |>
  dplyr::inner_join(tow_scores, by = "visit_id") |>
  dplyr::filter(.data$coverage >= 0.5, is.finite(.data$score)) |>
  dplyr::mutate(
    positive = .data$manta_evidence_class %in% strong_classes,
    survey = paste(.data$reef_name, .data$survey_date),
    tow_km = dplyr::if_else(is.finite(.data$tow_km) & .data$tow_km > 0, .data$tow_km, NA_real_),
    never_culled = !(.data$reef_name %in% cull_reefs)
  ) |>
  dplyr::group_by(.data$survey) |>
  dplyr::filter(dplyr::n() >= 3) |>
  dplyr::ungroup()
log_message("Manta: ", format(nrow(manta), big.mark = ","), " tows in ",
            format(dplyr::n_distinct(manta$survey), big.mark = ","), " reef surveys")

# Per-survey outcomes when the top k tows (by layer score) are kept.
manta_by_survey <- function(d, k) {
  d |>
    dplyr::group_by(.data$survey, .data$reef_year_id) |>
    dplyr::mutate(rank = rank(-.data$score, ties.method = "random")) |>
    dplyr::summarise(
      n = dplyr::n(), m = sum(.data$positive), kept = ceiling(k * dplyr::n()),
      kept_pos = sum(.data$positive & .data$rank <= ceiling(k * dplyr::n())),
      km = sum(.data$tow_km, na.rm = TRUE),
      km_kept = sum(.data$tow_km[.data$rank <= ceiling(k * dplyr::n())], na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      detected = .data$kept_pos > 0,
      # Random selection of the same number of tows: hypergeometric chance of
      # keeping at least one strong tow.
      p_detect_random = ifelse(.data$m == 0, 0, 1 - exp(lchoose(.data$n - .data$m, .data$kept) - lchoose(.data$n, .data$kept))),
      kept_pos_random = .data$m * .data$kept / .data$n
    )
}
summarise_manta <- function(s) {
  ry <- s |>
    dplyr::group_by(.data$reef_year_id) |>
    dplyr::summarise(any_pos = any(.data$m > 0), detected = any(.data$detected),
                     p_random = 1 - prod(1 - .data$p_detect_random), .groups = "drop") |>
    dplyr::filter(.data$any_pos)
  pos <- s[s$m > 0, ]
  tibble::tibble(
    tows = sum(s$n), tows_kept_share = sum(s$kept) / sum(s$n),
    km_saved = sum(s$km - s$km_kept),
    strong_tows = sum(s$m),
    strong_tows_kept = sum(s$kept_pos) / sum(s$m),
    strong_tows_kept_random = sum(s$kept_pos_random) / sum(s$m),
    outbreak_surveys = nrow(pos),
    surveys_detected = mean(pos$detected),
    surveys_detected_random = mean(pos$p_detect_random),
    outbreak_reef_years = nrow(ry),
    reef_years_detected = mean(ry$detected),
    reef_years_detected_random = mean(ry$p_random)
  )
}
manta_subsets <- list(
  `all surveys 2012-2026` = rep(TRUE, nrow(manta)),
  `surveys 2022-2026` = manta$year >= 2022,
  `reefs never culled` = manta$never_culled
)
manta_results <- dplyr::bind_rows(lapply(names(manta_subsets), function(sub) {
  d <- manta[manta_subsets[[sub]], ]
  dplyr::bind_rows(lapply(ops$keep_fractions, function(k) {
    s <- manta_by_survey(d, k)
    out <- summarise_manta(s) |> dplyr::mutate(subset = sub, keep = k, .before = 1)
    if (sub == "all surveys 2012-2026" && k == ops$headline_keep) {
      groups <- split(s, s$reef_year_id)
      boot <- replicate(ops$bootstrap_reps, {
        pick <- sample.int(length(groups), replace = TRUE)
        b <- dplyr::bind_rows(groups[pick], .id = "draw")
        b$reef_year_id <- b$draw
        unlist(summarise_manta(b)[, c("strong_tows_kept", "reef_years_detected")])
      })
      out$strong_tows_kept_low <- stats::quantile(boot[1, ], 0.025)
      out$strong_tows_kept_high <- stats::quantile(boot[1, ], 0.975)
      out$reef_years_detected_low <- stats::quantile(boot[2, ], 0.025)
      out$reef_years_detected_high <- stats::quantile(boot[2, ], 0.975)
    }
    out
  }))
}))
write_table(manta_results, cfg, "operations_manta_scenarios.csv")

# --- Cull dives -----------------------------------------------------------------
oof <- readRDS(file.path(cfg$exploration_dir, "outputs", "annual_rank_predictions.rds"))
dives_all <- readRDS(file.path(cfg$exploration_dir, "artifacts", "cull_observations_clean.rds")) |>
  dplyr::filter(is.finite(.data$total_cots), is.finite(.data$bottom_time), .data$bottom_time > 0) |>
  dplyr::mutate(hours = .data$bottom_time / 60, ghost = .data$total_cots == 0,
                low = .data$total_cots / .data$bottom_time < cfg$legacy$cpue_threshold)

# Score sources: held-out legacy predictions (site excluded from the model
# that scored it) and, for 2022-2026, the forward history + habitat score
# from 09b (every input from earlier years only).
oof_scores <- function(model_name) {
  oof |>
    dplyr::filter(.data$model == model_name) |>
    dplyr::distinct(.data$site_id, .data$year, .keep_all = TRUE) |>
    dplyr::select("site_id", "year", score = "probability")
}
forward <- readRDS(file.path(cfg$dirs$cache, "forward_test_scores.rds")) |>
  dplyr::distinct(.data$site_id, .data$year, .keep_all = TRUE)
score_sources <- list(
  annual_rg_year = oof_scores("annual_rg_year"),
  annual_clean_year = oof_scores("annual_clean_year"),
  annual_rg_year_2022_2026 = oof_scores("annual_rg_year") |>
    dplyr::semi_join(forward, by = c("site_id", "year")),
  history_combined_forward = forward |>
    dplyr::semi_join(oof_scores("annual_rg_year"), by = c("site_id", "year")) |>
    dplyr::select("site_id", "year", score = "combined_fixed")
)

cull_scenarios <- function(model_name) {
  scores <- score_sources[[model_name]]
  dives <- dives_all |>
    dplyr::inner_join(scores, by = c("site_id", "year")) |>
    dplyr::mutate(reef_year = paste(.data$reef_name, .data$year))
  sites <- dives |>
    dplyr::group_by(.data$reef_year, .data$site_id) |>
    dplyr::summarise(score = dplyr::first(.data$score), .groups = "drop") |>
    dplyr::group_by(.data$reef_year) |>
    dplyr::filter(dplyr::n() >= 3) |>
    dplyr::ungroup()
  dives <- dplyr::semi_join(dives, sites, by = c("reef_year", "site_id"))
  totals <- list(hours = sum(dives$hours), cots = sum(dives$total_cots),
                 ghost = sum(dives$ghost), ghost_hours = sum(dives$hours[dives$ghost]),
                 low_hours = sum(dives$hours[dives$low]))
  outcome <- function(kept_sites) {
    kept <- paste(dives$reef_year, dives$site_id) %in% kept_sites
    cpue_kept <- sum(dives$total_cots[kept]) / sum(dives$hours[kept])
    hours_saved <- sum(dives$hours[!kept])
    cots_missed <- sum(dives$total_cots[!kept])
    c(hours_kept = sum(dives$hours[kept]) / totals$hours,
      hours_saved = hours_saved,
      cots_kept = sum(dives$total_cots[kept]) / totals$cots,
      ghost_avoided = sum(dives$ghost[!kept]) / totals$ghost,
      ghost_hours_avoided = sum(dives$hours[!kept & dives$ghost]),
      low_hours_avoided = sum(dives$hours[!kept & dives$low]),
      ghost_rate_kept = mean(dives$ghost[kept]),
      cpue_gain = cpue_kept / (totals$cots / totals$hours),
      # Saved hours redeployed at the kept sites' catch rate, minus COTS left behind.
      net_redeploy_gain = (hours_saved * cpue_kept - cots_missed) / totals$cots)
  }
  site_groups <- split(sites, sites$reef_year)
  dplyr::bind_rows(lapply(ops$keep_fractions, function(k) {
    keep_model <- unlist(lapply(site_groups, function(g) {
      n_keep <- ceiling(k * nrow(g))
      paste(g$reef_year, g$site_id)[order(-g$score)][seq_len(n_keep)]
    }), use.names = FALSE)
    model <- outcome(keep_model)
    random <- rowMeans(replicate(ops$random_reps, outcome(unlist(lapply(site_groups, function(g) {
      paste(g$reef_year, g$site_id)[sample.int(nrow(g), ceiling(k * nrow(g)))]
    }), use.names = FALSE))))
    dplyr::bind_rows(
      tibble::as_tibble(as.list(model)) |> dplyr::mutate(selection = "layer"),
      tibble::as_tibble(as.list(random)) |> dplyr::mutate(selection = "random")
    ) |>
      dplyr::mutate(model = model_name, keep = k, dives = nrow(dives), reef_years = length(site_groups),
                    sites = nrow(sites), total_hours = totals$hours, total_cots = totals$cots,
                    ghost_dives = totals$ghost, .before = 1)
  }))
}
cull_results <- dplyr::bind_rows(lapply(names(score_sources), cull_scenarios)) |>
  dplyr::mutate(cost_saved_aud = .data$hours_saved * ops$diver_hour_cost_aud)
write_table(cull_results, cfg, "operations_cull_scenarios.csv")

# --- Headline (current programme era, 2022-2026) ------------------------------
h <- ops$headline_keep
pick_cull <- function(model, selection, column) {
  cull_results[[column]][cull_results$model == model & cull_results$keep == h & cull_results$selection == selection]
}
headline_rows <- function(period, manta_subset, cull_model) {
  mh <- manta_results[manta_results$subset == manta_subset & manta_results$keep == h, ]
  tibble::tribble(
    ~metric, ~primary_layer, ~random,
    "Manta: share of tows kept", mh$tows_kept_share, mh$tows_kept_share,
    "Manta: tow km saved", mh$km_saved, mh$km_saved,
    "Manta: strong-evidence tows kept", mh$strong_tows_kept, mh$strong_tows_kept_random,
    "Manta: outbreak reef-years still detected", mh$reef_years_detected, mh$reef_years_detected_random,
    "Cull: share of diver-hours kept", pick_cull(cull_model, "layer", "hours_kept"), pick_cull(cull_model, "random", "hours_kept"),
    "Cull: share of COTS at kept sites", pick_cull(cull_model, "layer", "cots_kept"), pick_cull(cull_model, "random", "cots_kept"),
    "Cull: share of ghost dives avoided", pick_cull(cull_model, "layer", "ghost_avoided"), pick_cull(cull_model, "random", "ghost_avoided"),
    "Cull: diver-hours saved", pick_cull(cull_model, "layer", "hours_saved"), pick_cull(cull_model, "random", "hours_saved"),
    "Cull: catch-rate ratio, kept sites / all", pick_cull(cull_model, "layer", "cpue_gain"), pick_cull(cull_model, "random", "cpue_gain"),
    "Cull: net COTS gain if saved hours redeployed", pick_cull(cull_model, "layer", "net_redeploy_gain"), pick_cull(cull_model, "random", "net_redeploy_gain")
  ) |> dplyr::mutate(period = period, keep_fraction = h, .before = 1)
}
headline <- dplyr::bind_rows(
  headline_rows("2022-2026", "surveys 2022-2026", "annual_rg_year_2022_2026"),
  headline_rows("all years", "all surveys 2012-2026", "annual_rg_year")
)
history_row <- tibble::tibble(
  period = "2022-2026", keep_fraction = h,
  metric = c("Cull (history + habitat, forward): share of COTS at kept sites",
             "Cull (history + habitat, forward): share of ghost dives avoided"),
  primary_layer = c(pick_cull("history_combined_forward", "layer", "cots_kept"),
                    pick_cull("history_combined_forward", "layer", "ghost_avoided")),
  random = c(pick_cull("history_combined_forward", "random", "cots_kept"),
             pick_cull("history_combined_forward", "random", "ghost_avoided"))
)
headline <- dplyr::bind_rows(headline, history_row)
write_table(headline, cfg, "operations_headline.csv")
print(headline, n = Inf)
