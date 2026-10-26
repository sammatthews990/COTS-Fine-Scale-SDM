# Does manta-tow coral cover help pick problematic cull sites, and does
# refitting the cull-site ranking at CPUE > 0.04 give a different or better
# outcome than the deployed CPUE >= 0.02 ranking?
#
# Coral. The 2027 ranking (20) uses the habitat rank and last season's Cull
# dives and manta tows, but not the coral cover recorded on each tow. Each
# tow's HardCoralProportionRangeCode, and RecentlyDeadCoralProportionRangeCode
# (recorded from Nov 2018), is converted to its bin midpoint (README sheet of
# the manta workbook) and smoothed onto every designated site with the same
# 300 m kernel as the manta evidence (last season's tows only). Variants of
# the site-evidence model (16, 19):
#   base                 habitat + evidence (the deployed formula)
#   coral                + hard coral cover nearby
#   coral_x_manta        + coral, and coral x strong manta tows nearby (COTS
#                        seen where coral is high)
#   coral_x_manta_rdc    + recently dead coral nearby
# Threshold. Every fitted part is refitted with "problematic" defined as
# CPUE >= 0.02 (deployed) and CPUE > 0.04 (as in 23): the leave-reef-out RG +
# year XGBoost habitat rank (17's cull model), the evidence terms built from
# problematic dives at and near the site, and the evidence model.
# Validation: strictly forward as in 19 (fit on seasons before T, rank every
# designated site on each reef dived in T, T = 2019..2026; habitat ranks are
# leave-reef-out). Every ranking is scored against both definitions, so the
# 0.04-trained ranking can be compared with the 0.02-trained one on either
# target. Metrics: within reef-year AUC (problematic vs below-threshold
# positive, vs ghost); at 80/90/95 % of problematic dives kept, the share of
# designated sites searched and of below-threshold and ghost dives avoided
# (random selection: sites searched = share kept). 200 paired bootstrap
# resamples of reef-years. Nothing here changes the 2027 product.

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "..", "R", "common.R"))
cfg <- initialise_final()
source(file.path(cfg$final_dir, "R", "forward.R"))
source(file.path(cfg$final_dir, "R", "cull_sites.R"))
source(file.path(cfg$final_dir, "R", "site_selection.R"))
set.seed(cfg$hmsc$seed + 24L)
thresholds <- list(t02 = function(cpue) cpue >= 0.02, t04 = function(cpue) cpue > 0.04)
threshold_labels <- c(t02 = "CPUE >= 0.02", t04 = "CPUE > 0.04")
test_years <- 2019:2026
targets <- c(0.8, 0.9, 0.95)
fractions <- round(seq(0.05, 1, 0.05), 2)
boot_reps <- 200L
n_folds <- 5L
kernel_sigma_m <- 300
rdc_from <- as.Date("2018-11-01")
strong_classes <- c("direct_cots", "multiple_counted_scars", "legacy_common_scars")
key <- c("reef_name", "year", "site_row")
params <- list(objective = "binary:logistic", eta = 0.00137382379588326, max_depth = 5,
               min_child_weight = 20, gamma = 1e-10, nthread = 8, tree_method = "hist")
rounds <- 1500L
rg <- cfg$legacy$rg_layer_names
features <- c(rg, "Year")
# Cover-class midpoints (%), README sheet of the manta workbook.
cover_mid <- c(`0` = 0, `1-` = 2.5, `1+` = 7.5, `2-` = 15, `2+` = 25, `3-` = 35, `3+` = 45,
               `4-` = 56.25, `4+` = 68.75, `5-` = 81.25, `5+` = 93.75,
               `1` = 5, `2` = 20, `3` = 40, `4` = 62.5, `5` = 87.5)

# --- Designated sites, candidate site-years, dives -----------------------------------
sites_sf <- read_designated_sites(cfg)
site_xy <- sf::st_coordinates(suppressWarnings(sf::st_point_on_surface(sf::st_geometry(sites_sf))))[, 1:2]
site_years <- readRDS(file.path(cfg$dirs$cache, "cull_site_ranks.rds"))[, key]
dives <- matched_cull_dives(cfg, sites_sf)
has_xy <- is.finite(dives$longitude) & is.finite(dives$latitude)
xy <- matrix(NA_real_, nrow(dives), 2)
xy[has_xy, ] <- sf::st_coordinates(sf::st_transform(sf::st_as_sf(dives[has_xy, ], coords = c("longitude", "latitude"), crs = 4326), sf::st_crs(sites_sf)))
fill <- !has_xy & !is.na(dives$site_row)
xy[fill, ] <- site_xy[dives$site_row[fill], ]
dives$x <- xy[, 1]
dives$y <- xy[, 2]
for (t in names(thresholds)) dives[[paste0("problematic_", t)]] <- thresholds[[t]](dives$cpue)

# --- Manta tows with coral cover ------------------------------------------------------
wb_path <- checked_input(cfg, "manta_workbook")
sheet <- readxl::read_excel(wb_path, sheet = "Manta", col_types = "text")
coral <- tibble::tibble(visit_id = paste0(basename(wb_path), "||", seq_len(nrow(sheet))),
                        sheet_date = as.Date(substr(sheet$SurveyTime, 1, 10)),
                        hard_coral = unname(cover_mid[trimws(sheet$HardCoralProportionRangeCode)]),
                        dead_coral = unname(cover_mid[trimws(sheet$RecentlyDeadCoralProportionRangeCode)]))
coral$dead_coral[coral$sheet_date < rdc_from] <- NA
tows <- readRDS(checked_input(cfg, "manta_tows")) |>
  dplyr::select("visit_id", "survey_date", "year", "reef_name", "manta_evidence_class")
geom <- readRDS(checked_input(cfg, "footprint_links"))$manta |>
  dplyr::group_by(.data$visit_id) |>
  dplyr::summarise(x = (dplyr::first(.data$start_x) + dplyr::first(.data$end_x)) / 2,
                   y = (dplyr::first(.data$start_y) + dplyr::first(.data$end_y)) / 2, .groups = "drop")
tows <- dplyr::inner_join(tows, geom, by = "visit_id") |>
  dplyr::left_join(coral, by = "visit_id") |>
  dplyr::filter(is.finite(.data$x), is.finite(.data$y)) |>
  dplyr::mutate(strong = .data$manta_evidence_class %in% strong_classes)
stopifnot(all(tows$sheet_date == tows$survey_date, na.rm = TRUE), !anyNA(tows$sheet_date))
log_message(format(nrow(tows), big.mark = ","), " tows; hard coral coded on ",
            sprintf("%.1f%%", 100 * mean(is.finite(tows$hard_coral))), ", recently dead coral on ",
            sprintf("%.1f%%", 100 * mean(is.finite(tows$dead_coral))))

# --- Leave-reef-out habitat ranks at each threshold (17's cull model) -----------------
stack <- terra::rast(checked_input(cfg, "predictor_stack_rg"))
names(stack) <- rg
train_sites <- readRDS(file.path(cfg$exploration_dir, "artifacts", "site_year_all.rds")) |>
  dplyr::filter(!is.na(.data$high_density), is.finite(.data$longitude), is.finite(.data$latitude))
stopifnot(all(train_sites$high_density == thresholds$t02(train_sites$cpue_mean)))
sxy <- sf::st_coordinates(sf::st_transform(sf::st_as_sf(train_sites, coords = c("longitude", "latitude"), crs = 4326), terra::crs(stack)))
train_x <- terra::extract(stack, sxy)[, rg]
train_x$Year <- train_sites$year
site_pts <- readRDS(file.path(cfg$dirs$cache, "cull_sites_points.rds"))
pts_by_site <- split(seq_along(site_pts$site), site_pts$site)
reefs <- unique(c(train_sites$reef_name, site_years$reef_name))
fold_of <- stats::setNames(sample(rep(seq_len(n_folds), length.out = length(reefs))), reefs)
for (t in names(thresholds)) site_years[[paste0("hab_", t)]] <- NA_real_
for (k in seq_len(n_folds)) {
  train <- fold_of[train_sites$reef_name] != k
  test_sy <- which(fold_of[site_years$reef_name] == k)
  p_idx <- pts_by_site[as.character(site_years$site_row[test_sy])]
  n_pts <- lengths(p_idx)
  x <- site_pts$values[unlist(p_idx, use.names = FALSE), rg]
  x$Year <- rep(site_years$year[test_sy], n_pts)
  owner <- rep(seq_along(test_sy), n_pts)
  for (t in names(thresholds)) {
    model <- xgboost::xgb.train(params = params, nrounds = rounds, verbose = 0,
                                data = xgboost::xgb.DMatrix(as.matrix(train_x[train, features]),
                                                            label = as.numeric(thresholds[[t]](train_sites$cpue_mean[train])), missing = NA))
    p <- stats::predict(model, xgboost::xgb.DMatrix(as.matrix(x[, features]), missing = NA))
    s <- tapply(p, owner, mean, na.rm = TRUE)
    site_years[[paste0("hab_", t)]][test_sy[as.integer(names(s))]] <- s
  }
  log_message("Habitat fold ", k, ": ", format(length(test_sy), big.mark = ","), " site-years scored")
}
site_years <- site_years |>
  dplyr::group_by(.data$reef_name, .data$year) |>
  dplyr::mutate(rank_hab_t02 = pct_rank(.data$hab_t02), rank_hab_t04 = pct_rank(.data$hab_t04)) |>
  dplyr::ungroup()

# --- Last season's evidence for every candidate site-year -----------------------------
evidence <- dplyr::bind_rows(lapply(sort(unique(site_years$year)), function(T) {
  rows <- site_years$site_row[site_years$year == T]
  q <- site_xy[rows, , drop = FALSE]
  pd <- dives[dives$year == T - 1 & is.finite(dives$x), ]
  pt <- tows[tows$year == T - 1, ]
  kern <- function(p, y) history_kernel(cbind(p$x, p$y), y, rep(0, nrow(p)), q, kernel_sigma_m, Inf)
  mk <- kern(pt, as.numeric(pt$strong))
  hc <- kern(pt[is.finite(pt$hard_coral), ], pt$hard_coral[is.finite(pt$hard_coral)])
  rd <- kern(pt[is.finite(pt$dead_coral), ], pt$dead_coral[is.finite(pt$dead_coral)])
  out <- data.frame(year = T, site_row = rows, reef_name = site_years$reef_name[site_years$year == T],
                    manta_strong = ifelse(mk$W > 0, mk$H * mk$W, 0), manta_effort = mk$W,
                    hard_coral = ifelse(hc$W > 0, hc$H, NA), dead_coral = ifelse(rd$W > 0, rd$H, NA))
  for (t in names(thresholds)) {
    pp <- pd[[paste0("problematic_", t)]]
    ck <- kern(pd, as.numeric(pp))
    at_site <- tapply(pp[!is.na(pd$site_row)], pd$site_row[!is.na(pd$site_row)], any)
    m <- at_site[as.character(rows)]
    out[[paste0("cull_rate_", t)]] <- ifelse(ck$W > 0, ck$H, 0)
    out$cull_effort <- ck$W
    out[[paste0("site_problematic_prev_", t)]] <- as.numeric(m %in% TRUE)
    out[[paste0("site_quiet_prev_", t)]] <- as.numeric(m %in% FALSE)
  }
  out
}))
feat <- dplyr::inner_join(site_years, evidence, by = key)
coral_ref <- mean(feat$hard_coral, na.rm = TRUE)
dead_ref <- mean(feat$dead_coral, na.rm = TRUE)

# Model terms at threshold t. Coral is centred and per 10 % cover; recently dead
# coral per 1 % cover; both are 0 (the mean) where no tow was nearby.
terms_for <- function(f, t) {
  data.frame(habitat = ifelse(is.finite(f[[paste0("rank_hab_", t)]]), f[[paste0("rank_hab_", t)]], 0.5),
             site_problematic_prev = f[[paste0("site_problematic_prev_", t)]],
             site_quiet_prev = f[[paste0("site_quiet_prev_", t)]],
             cull_rate = f[[paste0("cull_rate_", t)]], log_cull_effort = log1p(f$cull_effort),
             log_manta_strong = log1p(f$manta_strong), log_manta_effort = log1p(f$manta_effort),
             coral = ifelse(is.finite(f$hard_coral), (f$hard_coral - coral_ref) / 10, 0),
             dead_coral = ifelse(is.finite(f$dead_coral), f$dead_coral - dead_ref, 0))
}
base_formula <- problematic ~ habitat + site_problematic_prev + site_quiet_prev + cull_rate + log_cull_effort +
  log_manta_strong + log_manta_effort
formulas <- list(base = base_formula,
                 coral = stats::update(base_formula, . ~ . + coral),
                 coral_x_manta = stats::update(base_formula, . ~ . + coral + coral:log_manta_strong),
                 coral_x_manta_rdc = stats::update(base_formula, . ~ . + coral + coral:log_manta_strong + dead_coral))

# --- Checks against the deployed pipeline (0.02) -------------------------------------
old_ev <- readRDS(file.path(cfg$dirs$cache, "cull_site_evidence_scores.rds")) |>
  dplyr::filter(.data$window == "last_season") |>
  dplyr::select(dplyr::all_of(key), "cull_rate", "manta_strong", "site_problematic_prev") |>
  dplyr::inner_join(readRDS(file.path(cfg$dirs$cache, "manta_response_site_ranks.rds"))[, c(key, "rank_cull")], by = key) |>
  dplyr::inner_join(feat, by = key, suffix = c("_old", ""))
checks <- tibble::tibble(
  check = c("site-years compared", "max |cull_rate - 16|", "max |manta_strong - 16|",
            "site_problematic_prev agreement with 16", "Spearman habitat rank vs 17 (0.02)",
            "Spearman habitat rank 0.02 vs 0.04", "candidate site-years with nearby coral data",
            "mean nearby hard coral (%)", "mean nearby recently dead coral (%)"),
  value = c(nrow(old_ev), max(abs(old_ev$cull_rate - old_ev$cull_rate_t02)),
            max(abs(old_ev$manta_strong_old - old_ev$manta_strong)),
            mean(old_ev$site_problematic_prev == old_ev$site_problematic_prev_t02),
            stats::cor(old_ev$rank_cull, old_ev$rank_hab_t02, method = "spearman", use = "complete.obs"),
            stats::cor(feat$rank_hab_t02, feat$rank_hab_t04, method = "spearman", use = "complete.obs"),
            mean(is.finite(feat$hard_coral)), coral_ref, dead_ref))
write_table(checks, cfg, "threshold_coral_checks.csv")
print(checks)

# --- Training rows: one per matched dive -----------------------------------------------
rows <- dplyr::inner_join(dives[!is.na(dives$site_row), c(key, "cpue", "total_cots", "problematic_t02", "problematic_t04")],
                          feat, by = key)
fit_glm <- function(tr, t, v) {
  d <- terms_for(tr, t)
  d$problematic <- as.numeric(tr[[paste0("problematic_", t)]])
  stats::glm(formulas[[v]], family = stats::binomial(), data = d)
}

# Coefficients and likelihood-ratio tests on the deployed training seasons
# (2019-2026, as in 20). Reef-year-clustered standard errors; in-sample, so
# the forward test below is what decides.
fit_rows <- rows[rows$year >= 2019 & rows$year <= 2026, ]
fit_rows$reef_year <- paste(fit_rows$reef_name, fit_rows$year)
coef_rows <- list()
lr_rows <- list()
for (t in names(thresholds)) {
  models <- lapply(stats::setNames(names(formulas), names(formulas)), function(v) fit_glm(fit_rows, t, v))
  for (v in names(models)) {
    m <- models[[v]]
    se_cl <- sqrt(diag(sandwich::vcovCL(m, cluster = fit_rows$reef_year)))
    coef_rows[[length(coef_rows) + 1L]] <- tibble::tibble(
      threshold = threshold_labels[[t]], model = v, term = names(stats::coef(m)), estimate = unname(stats::coef(m)),
      std_error = unname(sqrt(diag(stats::vcov(m)))), std_error_reef_year = unname(se_cl),
      z_reef_year = unname(stats::coef(m) / se_cl), p_reef_year = 2 * stats::pnorm(-abs(unname(stats::coef(m) / se_cl))),
      dives = stats::nobs(m))
    if (v != "base") {
      a <- stats::anova(models$base, m, test = "Chisq")
      lr_rows[[length(lr_rows) + 1L]] <- tibble::tibble(threshold = threshold_labels[[t]], model = v,
                                                        deviance_drop = a$Deviance[2], df = a$Df[2], p = a$`Pr(>Chi)`[2],
                                                        aic_change = stats::AIC(m) - stats::AIC(models$base))
    }
  }
}
coefs <- dplyr::bind_rows(coef_rows)
lr <- dplyr::bind_rows(lr_rows)
write_table(coefs, cfg, "threshold_coral_coefficients.csv")
write_table(lr, cfg, "threshold_coral_lr_tests.csv")

# Observed problematic share by nearby coral cover and strong manta tows.
bins <- fit_rows |>
  dplyr::filter(is.finite(.data$hard_coral)) |>
  dplyr::mutate(coral_nearby = cut(.data$hard_coral, c(-Inf, 10, 20, 30, Inf), labels = c("<10%", "10-20%", "20-30%", ">=30%")),
                strong_manta_nearby = .data$manta_strong >= 0.5) |>
  dplyr::group_by(.data$strong_manta_nearby, .data$coral_nearby) |>
  dplyr::summarise(dives = dplyr::n(), reef_years = dplyr::n_distinct(.data$reef_year),
                   share_ge_002 = mean(.data$problematic_t02), share_gt_004 = mean(.data$problematic_t04),
                   share_ghost = mean(.data$total_cots == 0), .groups = "drop")
write_table(bins, cfg, "threshold_coral_bins.csv")
print(bins)

# --- Forward test ------------------------------------------------------------------------
variants <- c("habitat", names(formulas))
ranked <- list()
for (T in test_years) {
  test <- feat[feat$year == T, ]
  tr <- rows[rows$year < T, ]
  scores <- test[, key]
  for (t in names(thresholds)) {
    scores[[paste0("habitat_", t)]] <- test[[paste0("rank_hab_", t)]]
    for (v in names(formulas)) scores[[paste(v, t, sep = "_")]] <- stats::predict(fit_glm(tr, t, v), terms_for(test, t))
  }
  score_names <- setdiff(names(scores), key)
  scores <- scores |>
    dplyr::group_by(.data$reef_name, .data$year) |>
    dplyr::mutate(dplyr::across(dplyr::all_of(score_names), pct_rank)) |>
    dplyr::ungroup()
  ranked[[as.character(T)]] <- dplyr::inner_join(rows[rows$year == T, c(key, "cpue", "total_cots")], scores, by = key)
  log_message(T, ": ", format(nrow(ranked[[as.character(T)]]), big.mark = ","), " dives scored, ",
              format(nrow(tr), big.mark = ","), " training dives")
}
ranked <- dplyr::bind_rows(ranked) |> dplyr::mutate(reef_year = paste(.data$reef_name, .data$year))
score_names <- as.vector(outer(variants, names(thresholds), paste, sep = "_"))
saveRDS(ranked, file.path(cfg$dirs$cache, "threshold_coral_forward.rds"))

# Outcome classes under each evaluation definition. Below-threshold positives
# are split at 0.02 so the 0.04 evaluation shows which low dives are avoided.
outcome_for <- function(d, e) {
  dplyr::case_when(d$total_cots == 0 ~ "ghost", thresholds[[e]](d$cpue) ~ "problematic",
                   d$cpue < 0.02 ~ "low", TRUE ~ "intermediate")
}
for (e in names(thresholds)) ranked[[paste0("outcome_", e)]] <- outcome_for(ranked, e)

evaluate <- function(d, e) {
  outcome <- d[[paste0("outcome_", e)]]
  prob <- outcome == "problematic"
  below <- outcome %in% c("low", "intermediate")
  classes <- list(below = below, low = outcome == "low", intermediate = outcome == "intermediate",
                  ghost = outcome == "ghost")
  classes <- classes[vapply(classes, any, logical(1))]
  pl <- prob | below
  pg <- prob | classes$ghost
  unlist(lapply(score_names, function(s) {
    sc <- d[[s]]
    kept <- outer(sc, 1 - fractions - 1e-9, ">=")
    pk <- c(0, colSums(kept[prob, , drop = FALSE]) / sum(prob))
    v <- c(auc_vs_below = within_group_auc(sc[pl], prob[pl], d$reef_year[pl]),
           auc_vs_ghost = within_group_auc(sc[pg], prob[pg], d$reef_year[pg]))
    for (tg in targets) {
      v[[paste0("sites_searched_at_", 100 * tg)]] <- stats::approx(pk, c(0, fractions), xout = tg, ties = mean)$y
      for (cl in names(classes)) {
        y <- c(1, 1 - colSums(kept[classes[[cl]], , drop = FALSE]) / sum(classes[[cl]]))
        v[[paste0(cl, "_avoided_at_", 100 * tg)]] <- stats::approx(pk, y, xout = tg, ties = mean)$y
      }
    }
    stats::setNames(v, paste(s, names(v), sep = "|"))
  }))
}
d <- as.data.frame(ranked)
groups <- split(seq_len(nrow(d)), d$reef_year)
boot_idx <- lapply(seq_len(boot_reps), function(r) unlist(groups[sample.int(length(groups), replace = TRUE)], use.names = FALSE))
results <- list()
differences <- list()
for (e in names(thresholds)) {
  point <- evaluate(d, e)
  boot <- vapply(boot_idx, function(i) evaluate(d[i, ], e), numeric(length(point)))
  parts <- do.call(rbind, strsplit(names(point), "|", fixed = TRUE))
  info <- tibble::tibble(score = parts[, 1], metric = parts[, 2],
                         variant = sub("_t0[24]$", "", parts[, 1]), trained_on = unname(threshold_labels[sub(".*_(t0[24])$", "\\1", parts[, 1])]))
  results[[e]] <- dplyr::mutate(info, evaluated_on = threshold_labels[[e]], value = point,
                                low = apply(boot, 1, stats::quantile, 0.025, na.rm = TRUE),
                                high = apply(boot, 1, stats::quantile, 0.975, na.rm = TRUE),
                                dives = nrow(d), reef_years = length(groups),
                                problematic_share = mean(d[[paste0("outcome_", e)]] == "problematic"))
  # Paired differences: each variant vs base (same training threshold), and
  # 0.04-trained minus 0.02-trained (same variant).
  pair <- function(a, b, comparison) {
    ia <- which(info$score == a)
    ib <- match(paste(b, info$metric[ia]), paste(info$score, info$metric))
    diff <- boot[ia, , drop = FALSE] - boot[ib, , drop = FALSE]
    tibble::tibble(comparison = comparison, evaluated_on = threshold_labels[[e]], metric = info$metric[ia],
                   difference = point[ia] - point[ib],
                   low = apply(diff, 1, stats::quantile, 0.025, na.rm = TRUE),
                   high = apply(diff, 1, stats::quantile, 0.975, na.rm = TRUE))
  }
  for (t in names(thresholds)) for (v in setdiff(variants, "base")) {
    differences[[length(differences) + 1L]] <- pair(paste(v, t, sep = "_"), paste("base", t, sep = "_"),
                                                    paste0(v, " - base (trained ", threshold_labels[[t]], ")"))
  }
  for (v in variants) {
    differences[[length(differences) + 1L]] <- pair(paste(v, "t04", sep = "_"), paste(v, "t02", sep = "_"),
                                                    paste0(v, ": trained 0.04 - trained 0.02"))
  }
  log_message("Evaluated on ", threshold_labels[[e]])
}
results <- dplyr::bind_rows(results)
differences <- dplyr::bind_rows(differences)
write_table(results, cfg, "threshold_coral_forward.csv")
write_table(differences, cfg, "threshold_coral_differences.csv")

show <- c("auc_vs_below", "auc_vs_ghost", "sites_searched_at_90", "below_avoided_at_90", "ghost_avoided_at_90")
print(tidyr::pivot_wider(results[results$metric %in% show, c("evaluated_on", "trained_on", "variant", "metric", "value")],
                         names_from = "metric", values_from = "value"), n = Inf, width = Inf)
print(differences[differences$metric %in% show, ], n = Inf, width = Inf)
print(lr)
print(coefs[coefs$model == "coral_x_manta_rdc", ], n = Inf, width = Inf)
