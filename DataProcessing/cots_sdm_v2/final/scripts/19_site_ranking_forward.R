# Year-by-year forward test of the cull-site ranking (habitat + last season's
# evidence), to choose the training window for the 2027 site ranking (20).
#
# For each test season T = 2019..2026, the site-evidence model is fitted only
# on dives from seasons before T (each with its own previous-season evidence),
# then every designated site on each reef dived in T is ranked. Training
# windows: all earlier seasons, the last 5, the last 3. Habitat terms:
#   cull       held-out cull-trained rank (17)
#   ensemble   average of cull-trained and within reef-year manta ranks (17)
#   both       both ranks as separate terms
# plus habitat only (no fitting) as the reference. Habitat ranks are the
# leave-reef-out ranks from 17 (a reef's own dives never inform its rank);
# the evidence model is strictly forward.
# Scored per season and pooled: within reef-year AUC (problematic vs
# low-density, vs ghost) and, at 80 % and 90 % of problematic dives kept, the
# share of designated sites searched and of low-density and ghost dives
# avoided. Intervals: 200 bootstrap resamples of reef-years.

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "..", "R", "common.R"))
cfg <- initialise_final()
source(file.path(cfg$final_dir, "R", "cull_sites.R"))
source(file.path(cfg$final_dir, "R", "site_selection.R"))
set.seed(cfg$hmsc$seed + 19L)
test_years <- 2019:2026
windows <- c(all = Inf, last5 = 5, last3 = 3)
targets <- c(0.8, 0.9)
boot_reps <- 200L
key <- c("reef_name", "year", "site_row")

ev <- readRDS(file.path(cfg$dirs$cache, "cull_site_evidence_scores.rds")) |>
  dplyr::filter(.data$window == "last_season") |>
  dplyr::select(dplyr::all_of(key), "site_problematic_prev", "site_quiet_prev", "cull_rate", "cull_effort",
                "manta_strong", "manta_effort", "reef_cull_prev", "reef_manta_prev") |>
  dplyr::inner_join(readRDS(file.path(cfg$dirs$cache, "manta_response_site_ranks.rds")) |>
                      dplyr::select(dplyr::all_of(key), "rank_cull", "rank_within", "rank_ensemble"), by = key)
sites_sf <- read_designated_sites(cfg)
dives <- matched_cull_dives(cfg, sites_sf)
dives <- dives[!is.na(dives$site_row), c(key, "outcome")]
train_rows <- dplyr::inner_join(dives, ev, by = key) |>
  dplyr::mutate(problematic = as.numeric(.data$outcome == "problematic"))

habitats <- list(cull = list(h = "rank_cull"), ensemble = list(h = "rank_ensemble"),
                 both = list(h = "rank_cull", h2 = "rank_within"))
ranked <- list()
coefs <- list()
for (T in test_years) {
  test <- ev[ev$year == T, ]
  scores <- data.frame(test[, key], habitat_only = test$rank_cull)
  for (hn in names(habitats)) for (wn in names(windows)) {
    hab <- habitats[[hn]]
    w <- windows[[wn]]
    tr <- train_rows[train_rows$year < T & train_rows$year >= T - w, ]
    terms_tr <- site_evidence_terms(tr, hab$h, hab$h2)
    terms_tr$problematic <- tr$problematic
    model <- stats::glm(site_evidence_formula(habitat2 = !is.null(hab$h2)), family = stats::binomial(), data = terms_tr)
    scores[[paste(hn, wn, sep = "_")]] <- stats::predict(model, site_evidence_terms(test, hab$h, hab$h2))
    coefs[[length(coefs) + 1L]] <- tibble::tibble(test_year = T, habitat = hn, window = wn,
                                                  training_dives = nrow(tr), term = names(stats::coef(model)),
                                                  estimate = unname(stats::coef(model)))
  }
  score_names <- setdiff(names(scores), key)
  scores <- scores |>
    dplyr::group_by(.data$reef_name, .data$year) |>
    dplyr::mutate(dplyr::across(dplyr::all_of(score_names), pct_rank)) |>
    dplyr::ungroup()
  ranked[[as.character(T)]] <- dplyr::inner_join(dives[dives$year == T, ], scores, by = key)
  log_message(T, ": ", format(nrow(ranked[[as.character(T)]]), big.mark = ","), " dives scored")
}
ranked <- dplyr::bind_rows(ranked) |> dplyr::mutate(reef_year = paste(.data$reef_name, .data$year))
score_names <- setdiff(names(ranked), c(key, "outcome", "reef_year"))
saveRDS(ranked, file.path(cfg$dirs$cache, "site_ranking_forward.rds"))
write_table(dplyr::bind_rows(coefs), cfg, "site_ranking_forward_coefficients.csv")

# --- Per-season and pooled performance ---------------------------------------------
evaluate <- function(d) {
  unlist(lapply(score_names, function(s) {
    pl <- d$outcome %in% c("problematic", "low")
    pg <- d$outcome %in% c("problematic", "ghost")
    v <- c(auc_problematic_vs_low = within_group_auc(d[[s]][pl], d$outcome[pl] == "problematic", d$reef_year[pl]),
           auc_problematic_vs_ghost = within_group_auc(d[[s]][pg], d$outcome[pg] == "problematic", d$reef_year[pg]),
           matched_retention(d, s, targets))
    stats::setNames(v, paste(s, names(v), sep = "|"))
  }))
}
periods <- c(stats::setNames(as.list(test_years), test_years), list(`2019-2026` = test_years, `2022-2026` = 2022:2026))
results <- dplyr::bind_rows(lapply(names(periods), function(pn) {
  d <- as.data.frame(ranked[ranked$year %in% periods[[pn]], ])
  groups <- split(seq_len(nrow(d)), d$reef_year)
  point <- evaluate(d)
  boot <- vapply(seq_len(boot_reps), function(r) evaluate(d[unlist(groups[sample.int(length(groups), replace = TRUE)], use.names = FALSE), ]),
                 numeric(length(point)))
  parts <- do.call(rbind, strsplit(names(point), "|", fixed = TRUE))
  tibble::tibble(period = pn, score = parts[, 1], metric = parts[, 2], value = point,
                 low = apply(boot, 1, stats::quantile, 0.025, na.rm = TRUE),
                 high = apply(boot, 1, stats::quantile, 0.975, na.rm = TRUE),
                 dives = nrow(d), reef_years = length(groups),
                 problematic_share = mean(d$outcome == "problematic"))
}))
write_table(results, cfg, "site_ranking_forward.csv")

show <- c("auc_problematic_vs_low", "sites_searched_at_90", "low_dives_avoided_at_90", "ghost_dives_avoided_at_90")
print(tidyr::pivot_wider(results[results$period == "2019-2026" & results$metric %in% show, c("score", "metric", "value")],
                         names_from = "metric", values_from = "value"), n = Inf)
print(tidyr::pivot_wider(results[results$score %in% c("habitat_only", "cull_all", "ensemble_all", "both_all") &
                                   results$metric == "low_dives_avoided_at_90" & results$period %in% test_years,
                                 c("period", "score", "value")], names_from = "score", values_from = "value"), n = Inf)
