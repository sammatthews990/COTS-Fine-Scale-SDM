# Deployment-realistic validation: let the models use everything that would
# be known when the decision is made, including the reef's own history.
#
# Question: are the leave-reef-out statistics underselling the ranking? A
# deployed model has seen every reef the program works on, and crews plan a
# voyage knowing everything up to the day it leaves. Leave-reef-out answers
# "how does this work on a reef the model has never seen"; the deployment
# question is "how well does it rank the next visit, given all earlier data".
#
# Design: strictly forward in time. For each test season T = 2019..2026,
# every model is fitted only on data from before T, from ALL reefs (so the
# test reef's own earlier seasons are in the training data). Rankings are
# scored on the same dives, grouped by reef visit (reef x voyage).
#   Habitat (cull-trained RG + year XGBoost, legacy hyperparameters):
#     lro        leave-reef-out rank from 15 (current validation)
#     forward    refitted on site-years before T, all reefs
#     in-sample  the deployed layer, which has seen season T's own dives
#                (a leaked ceiling, reported only to show the leak)
#   Evidence (the site-evidence logistic model, fitted on dives 2013..T-1):
#     previous season   calendar season T-1 (as deployed)
#     to voyage start   the 365 days before the voyage left port: earlier
#                       voyages in the same season count, the voyage's own
#                       manta survey and dives do not
#   Combinations: lro + previous (current product validation), forward +
#   previous, lro + to-voyage, forward + to-voyage.
# Targets: >= 0.02 (training response) and > 0.04 (management target).
# Metrics: within reef-visit AUC, and at 90 % retention the share of
# designated sites selected and of low-density and ghost dives avoided;
# retention of > 0.04 dives at the A-C cutoff (top 80 % of sites).
# Intervals: 300 paired bootstrap resamples of reef-years.
# Manta: the within reef-year manta model (17-18), leave-reef-out vs
# forward (tows before T, all reefs), ranking tows in season T.

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "..", "R", "common.R"))
cfg <- initialise_final()
source(file.path(cfg$final_dir, "R", "forward.R"))
source(file.path(cfg$final_dir, "R", "cull_sites.R"))
source(file.path(cfg$final_dir, "R", "site_selection.R"))
set.seed(cfg$hmsc$seed + 33L)
test_years <- 2019:2026
first_year <- 2013L
boot_reps <- 300L
kernel_sigma_m <- 300
window_days <- 365
fractions <- round(seq(0.02, 1, 0.02), 2)
strong_classes <- c("direct_cots", "multiple_counted_scars", "legacy_common_scars")
params <- list(objective = "binary:logistic", eta = 0.00137382379588326, max_depth = 5,
               min_child_weight = 20, gamma = 1e-10, nthread = 8, tree_method = "hist")
rounds <- 1500L
key <- c("reef_name", "year", "site_row")

stack <- terra::rast(checked_input(cfg, "predictor_stack_rg"))
rg <- cfg$legacy$rg_layer_names
names(stack) <- rg
features <- c(rg, "Year")

# --- Sites, dives and visits -------------------------------------------------------
sites_sf <- read_designated_sites(cfg)
site_xy <- sf::st_coordinates(suppressWarnings(sf::st_point_on_surface(sf::st_geometry(sites_sf))))[, 1:2]
site_reef <- sites_sf$reef_name
site_pts <- readRDS(file.path(cfg$dirs$cache, "cull_sites_points.rds"))
pts_by_site <- split(seq_along(site_pts$site), site_pts$site)
ranks15 <- readRDS(file.path(cfg$dirs$cache, "cull_site_ranks.rds"))[, c(key, "rank_cv", "rank_deployed")]

dives <- matched_cull_dives(cfg, sites_sf) |>
  dplyr::filter(.data$year >= first_year, !is.na(.data$VoyageTitle)) |>
  dplyr::mutate(high = .data$cpue > 0.04, problematic = .data$outcome == "problematic",
                visit = paste(.data$reef_name, .data$VoyageTitle, .data$VesselName_ID, sep = "||"))
visit_info <- dives |>
  dplyr::group_by(.data$visit) |>
  dplyr::summarise(reef_name = dplyr::first(.data$reef_name), year = min(.data$year),
                   cut = min(as.Date(.data$VoyageStartDate), .data$survey_date), .groups = "drop")
dives$year <- visit_info$year[match(dives$visit, visit_info$visit)]   # a visit belongs to one season
has_xy <- is.finite(dives$longitude) & is.finite(dives$latitude)
xy <- matrix(NA_real_, nrow(dives), 2)
xy[has_xy, ] <- sf::st_coordinates(sf::st_transform(sf::st_as_sf(dives[has_xy, ], coords = c("longitude", "latitude"), crs = 4326), sf::st_crs(sites_sf)))
fill <- !has_xy & !is.na(dives$site_row)
xy[fill, ] <- site_xy[dives$site_row[fill], ]
dives$x <- xy[, 1]
dives$y <- xy[, 2]
all_dives <- dives
dives <- dives[!is.na(dives$site_row), ]

tow_geom <- readRDS(checked_input(cfg, "footprint_links"))$manta |>
  dplyr::group_by(.data$visit_id) |>
  dplyr::summarise(x = (dplyr::first(.data$start_x) + dplyr::first(.data$end_x)) / 2,
                   y = (dplyr::first(.data$start_y) + dplyr::first(.data$end_y)) / 2, .groups = "drop")
tows_ev <- readRDS(checked_input(cfg, "manta_tows")) |>
  dplyr::select("visit_id", "survey_date", "manta_evidence_class") |>
  dplyr::inner_join(tow_geom, by = "visit_id") |>
  dplyr::filter(is.finite(.data$x)) |>
  dplyr::mutate(strong = .data$manta_evidence_class %in% strong_classes)
log_message(format(nrow(dives), big.mark = ","), " matched dives in ", dplyr::n_distinct(dives$visit), " reef visits, ",
            first_year, "-", max(dives$year))

# --- Forward habitat ranks: one model per season, trained on earlier seasons ------
train_sites <- readRDS(file.path(cfg$exploration_dir, "artifacts", "site_year_all.rds")) |>
  dplyr::filter(!is.na(.data$high_density), is.finite(.data$longitude), is.finite(.data$latitude))
sxy <- sf::st_coordinates(sf::st_transform(sf::st_as_sf(train_sites, coords = c("longitude", "latitude"), crs = 4326), terra::crs(stack)))
train_x <- terra::extract(stack, sxy)[, rg]
train_x$Year <- train_sites$year
fwd_path <- file.path(cfg$dirs$cache, "deployment_forward_habitat.rds")
if (!file.exists(fwd_path)) {
  fwd <- list()
  for (t in first_year:max(test_years)) {
    tr <- train_sites$year < t
    model <- xgboost::xgb.train(params = params, nrounds = rounds, verbose = 0,
      data = xgboost::xgb.DMatrix(as.matrix(train_x[tr, features]), label = train_sites$high_density[tr], missing = NA))
    sy <- ranks15[ranks15$year == t, ]
    if (!nrow(sy)) next
    p_idx <- pts_by_site[as.character(sy$site_row)]
    x <- site_pts$values[unlist(p_idx, use.names = FALSE), rg]
    x$Year <- t
    p <- stats::predict(model, xgboost::xgb.DMatrix(as.matrix(x[, features]), missing = NA))
    s <- tapply(p, rep(seq_len(nrow(sy)), lengths(p_idx)), mean, na.rm = TRUE)
    sy$habitat_fwd <- NA_real_
    sy$habitat_fwd[as.integer(names(s))] <- s
    fwd[[as.character(t)]] <- sy[, c(key, "habitat_fwd")]
    log_message("Forward habitat ", t, ": trained on ", format(sum(tr), big.mark = ","), " site-years")
  }
  saveRDS(dplyr::bind_rows(fwd), fwd_path)
}
habitat <- readRDS(fwd_path) |>
  dplyr::group_by(.data$reef_name, .data$year) |>
  dplyr::mutate(rank_fwd = pct_rank(.data$habitat_fwd)) |>
  dplyr::ungroup() |>
  dplyr::inner_join(ranks15, by = key)

# --- Evidence: previous season (16) and up to voyage start -----------------------
ev_prev <- readRDS(file.path(cfg$dirs$cache, "cull_site_evidence_scores.rds")) |>
  dplyr::filter(.data$window == "last_season") |>
  dplyr::select(dplyr::all_of(key), "site_problematic_prev", "site_quiet_prev", "cull_rate", "cull_effort",
                "manta_strong", "manta_effort")
roll_path <- file.path(cfg$dirs$cache, "deployment_rolling_evidence.rds")
if (!file.exists(roll_path)) {
  vis <- visit_info[visit_info$reef_name %in% site_reef, ]
  roll <- vector("list", nrow(vis))
  for (i in seq_len(nrow(vis))) {
    rows <- which(site_reef == vis$reef_name[i])
    q <- site_xy[rows, , drop = FALSE]
    box <- c(min(q[, 1]) - 1000, max(q[, 1]) + 1000, min(q[, 2]) - 1000, max(q[, 2]) + 1000)
    start <- vis$cut[i] - window_days
    pd <- all_dives[all_dives$survey_date < vis$cut[i] & all_dives$survey_date >= start & is.finite(all_dives$x) &
                      all_dives$x > box[1] & all_dives$x < box[2] & all_dives$y > box[3] & all_dives$y < box[4], ]
    pt <- tows_ev[tows_ev$survey_date < vis$cut[i] & tows_ev$survey_date >= start &
                    tows_ev$x > box[1] & tows_ev$x < box[2] & tows_ev$y > box[3] & tows_ev$y < box[4], ]
    ck <- history_kernel(cbind(pd$x, pd$y), as.numeric(pd$outcome == "problematic"), rep(0, nrow(pd)), q, kernel_sigma_m, Inf)
    mk <- history_kernel(cbind(pt$x, pt$y), as.numeric(pt$strong), rep(0, nrow(pt)), q, kernel_sigma_m, Inf)
    at <- pd[!is.na(pd$site_row), ] |> dplyr::group_by(.data$site_row) |>
      dplyr::summarise(prob = any(.data$outcome == "problematic"), .groups = "drop")
    m <- match(rows, at$site_row)
    roll[[i]] <- data.frame(visit = vis$visit[i], reef_name = vis$reef_name[i], year = vis$year[i], site_row = rows,
                            site_problematic_prev = as.numeric(!is.na(m) & at$prob[m] %in% TRUE),
                            site_quiet_prev = as.numeric(!is.na(m) & at$prob[m] %in% FALSE),
                            cull_rate = ifelse(ck$W > 0, ck$H, 0), cull_effort = ck$W,
                            manta_strong = ifelse(mk$W > 0, mk$H * mk$W, 0), manta_effort = mk$W)
  }
  saveRDS(dplyr::bind_rows(roll), roll_path)
}
ev_roll <- readRDS(roll_path)
log_message(format(nrow(ev_roll), big.mark = ","), " site-visits with evidence up to voyage start")

# Site-visit table: every designated site on the reef, for every visit.
sv <- ev_roll |>
  dplyr::inner_join(habitat, by = key) |>
  dplyr::left_join(ev_prev, by = key, suffix = c("_roll", "_prev"))
for (v in c("site_problematic_prev", "site_quiet_prev", "cull_rate", "cull_effort", "manta_strong", "manta_effort")) {
  sv[[paste0(v, "_prev")]][is.na(sv[[paste0(v, "_prev")]])] <- 0
}
pick <- function(d, window) {
  out <- d
  for (v in c("site_problematic_prev", "site_quiet_prev", "cull_rate", "cull_effort", "manta_strong", "manta_effort")) {
    out[[v]] <- d[[paste0(v, "_", window)]]
  }
  out
}
dv <- dives[, c(key, "visit", "outcome", "problematic", "high", "cpue")] |>
  dplyr::inner_join(sv[, c("visit", "site_row")], by = c("visit", "site_row"))

# --- Forward evidence models, per test season --------------------------------------
combos <- list(lro_prev = c("rank_cv", "prev"), fwd_prev = c("rank_fwd", "prev"),
               lro_voyage = c("rank_cv", "roll"), fwd_voyage = c("rank_fwd", "roll"))
scored <- list()
coefs <- list()
for (T in test_years) {
  test <- sv[sv$year == T, ]
  out <- test[, c("visit", key)]
  out$habitat_lro <- test$rank_cv
  out$habitat_fwd <- test$rank_fwd
  out$habitat_insample <- test$rank_deployed
  for (cn in names(combos)) {
    h <- combos[[cn]][1]
    w <- combos[[cn]][2]
    tr_sv <- pick(sv[sv$year >= first_year & sv$year < T, ], w)
    tr <- dives[dives$year < T, c("visit", "site_row", "problematic")] |>
      dplyr::inner_join(tr_sv, by = c("visit", "site_row"))
    terms <- site_evidence_terms(tr, h)
    terms$problematic <- as.numeric(tr$problematic)
    model <- stats::glm(site_evidence_formula(), family = stats::binomial(), data = terms)
    out[[cn]] <- stats::predict(model, site_evidence_terms(pick(test, w), h))
    coefs[[length(coefs) + 1L]] <- tibble::tibble(test_year = T, ranking = cn, term = names(stats::coef(model)),
                                                  estimate = unname(stats::coef(model)), training_dives = nrow(tr))
  }
  score_cols <- setdiff(names(out), c("visit", key))
  out <- out |>
    dplyr::group_by(.data$visit) |>
    dplyr::mutate(dplyr::across(dplyr::all_of(score_cols), pct_rank)) |>
    dplyr::ungroup()
  scored[[as.character(T)]] <- dplyr::inner_join(dv[dv$year == T, ], out, by = c("visit", key))
  log_message(T, ": ", format(nrow(scored[[as.character(T)]]), big.mark = ","), " dives in ",
              dplyr::n_distinct(scored[[as.character(T)]]$visit), " visits")
}
scored <- dplyr::bind_rows(scored) |> dplyr::mutate(reef_year = paste(.data$reef_name, .data$year))
score_cols <- c("habitat_lro", "habitat_fwd", "habitat_insample", names(combos))
# Sites without a habitat prediction rank last (as in the deployed product).
for (col in score_cols) scored[[col]][!is.finite(scored[[col]])] <- -1
saveRDS(scored, file.path(cfg$dirs$cache, "deployment_validation_scored.rds"))write_table(dplyr::bind_rows(coefs), cfg, "deployment_validation_coefficients.csv")

# --- Metrics ---------------------------------------------------------------------------
curve <- function(d, col, target) {
  vapply(fractions, function(a) {
    kept <- d[[col]] >= 1 - a - 1e-9
    c(kept = sum(kept & d[[target]]) / sum(d[[target]]),
      low = sum(!kept & d$outcome == "low") / sum(d$outcome == "low"),
      ghost = sum(!kept & d$outcome == "ghost") / sum(d$outcome == "ghost"))
  }, numeric(3))
}
at_retention <- function(m, t) {
  pk <- c(0, m["kept", ])
  if (sum(is.finite(pk)) < 2) return(c(sites = NA_real_, low = NA_real_, ghost = NA_real_))
  c(sites = stats::approx(pk, c(0, fractions), xout = t, ties = mean)$y,
    low = stats::approx(pk, c(1, m["low", ]), xout = t, ties = mean)$y,
    ghost = stats::approx(pk, c(1, m["ghost", ]), xout = t, ties = mean)$y)
}
evaluate <- function(d) {
  unlist(lapply(score_cols, function(col) {
    pl <- d$outcome %in% c("problematic", "low")
    c02 <- curve(d, col, "problematic")
    c04 <- curve(d, col, "high")
    v <- c(auc_problematic_vs_low = within_group_auc(d[[col]][pl], d$problematic[pl], d$visit[pl]),
           auc_high_vs_rest = within_group_auc(d[[col]], d$high, d$visit),
           stats::setNames(at_retention(c02, 0.9), c("sites_at_90_02", "low_avoided_at_90_02", "ghost_avoided_at_90_02")),
           stats::setNames(at_retention(c04, 0.9), c("sites_at_90_04", "low_avoided_at_90_04", "ghost_avoided_at_90_04")),
           high_kept_top80 = sum(d[[col]] >= 0.2 - 1e-9 & d$high) / sum(d$high))
    stats::setNames(v, paste(col, names(v), sep = "|"))
  }))
}
periods <- c(stats::setNames(as.list(test_years), test_years), list(`2019-2026` = test_years, `2022-2026` = 2022:2026))
results <- dplyr::bind_rows(lapply(names(periods), function(pn) {
  d <- as.data.frame(scored[scored$year %in% periods[[pn]], ])
  point <- evaluate(d)
  pooled <- nchar(pn) > 4
  boot <- if (pooled) {
    groups <- split(seq_len(nrow(d)), d$reef_year)
    vapply(seq_len(boot_reps), function(r) evaluate(d[unlist(groups[sample.int(length(groups), replace = TRUE)], use.names = FALSE), ]),
           numeric(length(point)))
  } else NULL
  parts <- do.call(rbind, strsplit(names(point), "|", fixed = TRUE))
  tibble::tibble(period = pn, ranking = parts[, 1], metric = parts[, 2], value = point,
                 low = if (pooled) apply(boot, 1, stats::quantile, 0.025, na.rm = TRUE) else NA_real_,
                 high = if (pooled) apply(boot, 1, stats::quantile, 0.975, na.rm = TRUE) else NA_real_,
                 dives = nrow(d), visits = dplyr::n_distinct(d$visit), reef_years = dplyr::n_distinct(d$reef_year))
}))
write_table(results, cfg, "deployment_validation.csv")

# Paired differences against the current validation (lro + previous season).
d <- as.data.frame(scored)
groups <- split(seq_len(nrow(d)), d$reef_year)
point <- evaluate(d)
boot <- vapply(seq_len(boot_reps), function(r) evaluate(d[unlist(groups[sample.int(length(groups), replace = TRUE)], use.names = FALSE), ]),
               numeric(length(point)))
parts <- do.call(rbind, strsplit(names(point), "|", fixed = TRUE))
base <- match(paste("lro_prev", parts[, 2], sep = "|"), names(point))
diffs <- tibble::tibble(ranking = parts[, 1], metric = parts[, 2], difference = point - point[base],
                        low = apply(boot - boot[base, ], 1, stats::quantile, 0.025, na.rm = TRUE),
                        high = apply(boot - boot[base, ], 1, stats::quantile, 0.975, na.rm = TRUE)) |>
  dplyr::filter(.data$ranking != "lro_prev")
write_table(diffs, cfg, "deployment_validation_differences.csv")

# --- Manta: within reef-year model, leave-reef-out vs forward ----------------------
tow_scores <- readRDS(file.path(cfg$dirs$cache, "manta_response_tow_scores.rds"))
tows <- readRDS(checked_input(cfg, "manta_tows")) |>
  dplyr::select("visit_id", "year", "reef_name", "manta_evidence_class")
geom <- readRDS(checked_input(cfg, "footprint_links"))$manta |>
  dplyr::group_by(.data$visit_id) |>
  dplyr::summarise(sx = dplyr::first(.data$start_x), ex = dplyr::first(.data$end_x), .groups = "drop")
tows <- dplyr::inner_join(tows, geom, by = "visit_id") |>
  dplyr::filter(is.finite(.data$sx), is.finite(.data$ex)) |>
  dplyr::mutate(label = dplyr::case_when(.data$manta_evidence_class %in% strong_classes ~ 1,
                                         .data$manta_evidence_class == "strict_negative" ~ 0),
                reef_year = paste(.data$reef_name, .data$year)) |>
  dplyr::group_by(.data$reef_year) |>
  dplyr::mutate(ry_margin = stats::qlogis((sum(.data$label %in% 1) + 0.5) / (sum(!is.na(.data$label)) + 1))) |>
  dplyr::ungroup()
tp <- readRDS(file.path(cfg$dirs$cache, "manta_tow_predictors.rds"))
v <- as.matrix(tp$values[, rg])
ok <- !is.na(v)
v[!ok] <- 0
sums <- rowsum(v, tp$tow)
counts <- rowsum(ok * 1, tp$tow)
tow_x <- matrix(NA_real_, nrow(tows), length(rg), dimnames = list(NULL, rg))
tow_x[as.integer(rownames(sums)), ] <- ifelse(counts > 0, sums / counts, NA)
tows$within_lro <- tow_scores$within_cv[match(tows$visit_id, tow_scores$visit_id)]
tows$within_fwd <- NA_real_
for (T in test_years) {
  tr <- !is.na(tows$label) & tows$year < T
  model <- xgboost::xgb.train(params = params, nrounds = rounds, verbose = 0,
    data = xgboost::xgb.DMatrix(tow_x[tr, ], label = tows$label[tr], base_margin = tows$ry_margin[tr], missing = NA))
  te <- which(tows$year == T)
  tows$within_fwd[te] <- stats::predict(model, xgboost::xgb.DMatrix(tow_x[te, , drop = FALSE], missing = NA), outputmargin = TRUE)
  log_message("Forward manta model ", T, ": trained on ", format(sum(tr), big.mark = ","), " tows")
}
tt <- tows[tows$year %in% test_years & !is.na(tows$label) & is.finite(tows$within_lro) & is.finite(tows$within_fwd), ]
tow_pieces <- function(d, col) {
  dplyr::bind_rows(lapply(split(seq_len(nrow(d)), d$reef_year), function(i) {
    y <- d$label[i] == 1
    n1 <- sum(y); n0 <- sum(!y)
    if (n1 == 0 || n0 == 0) return(NULL)
    r <- rank(d[[col]][i])
    top <- rank(-d[[col]][i], ties.method = "random") <= ceiling(0.2 * length(i))
    data.frame(reef_year = d$reef_year[i[1]], auc = (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0), pairs = n1 * n0,
               hits = sum(y & top), pos = n1)
  }))
}
pl <- tow_pieces(tt, "within_lro")
pf <- tow_pieces(tt, "within_fwd")
summ <- function(p) c(auc = stats::weighted.mean(p$auc, p$pairs), lift = sum(p$hits) / sum(p$pos) / 0.2)
bt <- replicate(boot_reps, {
  i <- sample.int(nrow(pl), replace = TRUE)
  c(summ(pl[i, ]), summ(pf[i, ]), summ(pf[i, ]) - summ(pl[i, ]))
})
manta <- tibble::tibble(model = rep(c("leave-reef-out", "forward (sees reef's earlier tows)", "forward minus leave-reef-out"), each = 2),
                        metric = rep(c("auc", "lift"), 3),
                        value = c(summ(pl), summ(pf), summ(pf) - summ(pl)),
                        low = apply(bt, 1, stats::quantile, 0.025), high = apply(bt, 1, stats::quantile, 0.975),
                        tows = nrow(tt), reef_years = nrow(pl))
write_table(manta, cfg, "deployment_validation_manta.csv")

show <- c("auc_problematic_vs_low", "auc_high_vs_rest", "sites_at_90_04", "low_avoided_at_90_04", "ghost_avoided_at_90_04", "high_kept_top80")
print(tidyr::pivot_wider(results[results$period == "2019-2026" & results$metric %in% show, c("ranking", "metric", "value")],
                         names_from = "metric", values_from = "value"), width = Inf)
print(diffs[diffs$metric %in% show, ], n = Inf)
print(manta)
