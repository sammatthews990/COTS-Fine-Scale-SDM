# If crews dived only the higher-ranked designated cull sites on each reef,
# what share of problematic dives (CPUE >= 0.02) would be missed, and what
# share of low-density dives (COTS present, < 0.02) and ghost dives (0 COTS)
# would be avoided?
#
# Candidates: every designated cull site on the reef (EoTR cull-site KMZ),
# including sites that have never been dived. The layer ranks all of them.
#
# Held-out scoring: the layer is trained on these Cull sites, so each reef is
# scored by an RG + year XGBoost (legacy hyperparameters) refitted without
# that reef's Cull data (5-fold leave-reef-out). A site's score is the mean
# prediction over its polygon (up to 16 points), for the year being assessed.
# Diving the top fraction a of a reef's sites keeps a dive when its site's
# within-reef rank is >= 1 - a. The deployed layer (trained on all reefs) is
# scored the same way for comparison. Uncertainty: 500 bootstrap resamples of
# reef-years.

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "..", "R", "common.R"))
cfg <- initialise_final()
set.seed(cfg$hmsc$seed)
n_folds <- 5L
boot_reps <- 500L
points_per_site <- 16L
site_fractions <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
legacy_params <- list(objective = "binary:logistic", eta = 0.00137382379588326, max_depth = 5,
                      min_child_weight = 20, gamma = 1e-10, nthread = 8)
legacy_rounds <- 1500L
threshold <- cfg$legacy$cpue_threshold

stack <- terra::rast(checked_input(cfg, "predictor_stack_rg"))
names(stack) <- cfg$legacy$rg_layer_names
features <- c(cfg$legacy$rg_layer_names, "Year")
composite <- terra::rast(file.path(cfg$dirs$layers, "cots_legacy_rg_year2026_clean_fallback_30m.tif"))
crs_mga <- terra::crs(stack)
clean_name <- function(x) toupper(gsub("^[[:space:](]+|[[:space:]]+$", "", x))

# --- Designated cull sites from the KMZ -------------------------------------
sites_path <- file.path(cfg$dirs$cache, "cull_sites_designated.gpkg")
if (!file.exists(sites_path)) {
  kml_dir <- file.path(cfg$dirs$cache, "cull_sites_kmz")
  utils::unzip(checked_input(cfg, "cull_sites_kmz"), exdir = kml_dir)
  kml <- list.files(kml_dir, pattern = "\\.kml$", full.names = TRUE, recursive = TRUE)[1]
  # One GDAL pass to a multi-layer GeoPackage; reading KML layer by layer
  # re-parses the whole file each time.
  staged <- file.path(kml_dir, "cull_sites_layers.gpkg")
  sf::gdal_utils("vectortranslate", kml, staged, options = c("-dim", "XY", "-overwrite"))
  layers <- sf::st_layers(staged)$name
  log_message("KMZ: ", length(layers), " reef layers")
  sites_sf <- dplyr::bind_rows(lapply(layers, function(l) {
    x <- sf::st_read(staged, layer = l, quiet = TRUE)
    sf::st_sf(site_name = x$Name, reef_name = sub("^Reef:\\s*", "", l), geometry = sf::st_zm(sf::st_geometry(x)))
  })) |>
    sf::st_make_valid() |>
    sf::st_transform(crs_mga)
  sf::st_write(sites_sf, sites_path, quiet = TRUE, delete_dsn = TRUE)
}
sites_sf <- sf::st_read(sites_path, quiet = TRUE)
sites_sf$site_key <- paste(sites_sf$reef_name, clean_name(sites_sf$site_name), sep = "||")
sites_sf <- sites_sf[!duplicated(sites_sf$site_key), ]
log_message(format(nrow(sites_sf), big.mark = ","), " designated cull sites on ",
            dplyr::n_distinct(sites_sf$reef_name), " reefs")

# Points inside each site polygon (regular sample, centroid fallback).
pts_cache <- file.path(cfg$dirs$cache, "cull_sites_points.rds")
if (!file.exists(pts_cache)) {
  pts <- lapply(seq_len(nrow(sites_sf)), function(i) {
    g <- sf::st_geometry(sites_sf)[i]
    p <- suppressWarnings(tryCatch(sf::st_sample(g, points_per_site, type = "regular"), error = function(e) NULL))
    xy <- if (is.null(p) || !length(p)) sf::st_coordinates(sf::st_point_on_surface(g)) else sf::st_coordinates(p)
    cbind(site = i, xy[, 1:2, drop = FALSE])
  })
  pts <- do.call(rbind, pts)
  saveRDS(list(site = pts[, 1], values = terra::extract(stack, pts[, 2:3]),
               deployed = terra::extract(composite, pts[, 2:3])[[1]]), pts_cache)
}
site_pts <- readRDS(pts_cache)
sites_sf$deployed <- as.numeric(tapply(site_pts$deployed, factor(site_pts$site, levels = seq_len(nrow(sites_sf))), mean, na.rm = TRUE))
sites <- sf::st_drop_geometry(sites_sf)
sites$row <- seq_len(nrow(sites))

# --- Dives matched to designated sites -----------------------------------------
dives <- readRDS(file.path(cfg$exploration_dir, "artifacts", "cull_observations_clean.rds")) |>
  dplyr::filter(is.finite(.data$total_cots), is.finite(.data$bottom_time), .data$bottom_time > 0) |>
  dplyr::mutate(cpue = .data$total_cots / .data$bottom_time,
                outcome = dplyr::case_when(.data$total_cots == 0 ~ "ghost",
                                           .data$cpue < threshold ~ "low",
                                           TRUE ~ "problematic"),
                site_key = paste(.data$reef_name, clean_name(.data$site_name), sep = "||"))
dives$site_row <- match(dives$site_key, sites$site_key)
# Spatial fallback: the dive point inside, or within 100 m of, a site on the same reef.
unmatched <- which(is.na(dives$site_row) & is.finite(dives$longitude))
if (length(unmatched)) {
  dpts <- sf::st_transform(sf::st_as_sf(dives[unmatched, ], coords = c("longitude", "latitude"), crs = 4326), crs_mga)
  near <- sf::st_nearest_feature(dpts, sites_sf)
  dist <- as.numeric(sf::st_distance(dpts, sites_sf[near, ], by_element = TRUE))
  ok <- dist <= 100 & sites$reef_name[near] == dives$reef_name[unmatched]
  dives$site_row[unmatched[ok]] <- near[ok]
}
match_rate <- mean(!is.na(dives$site_row))
log_message(sprintf("Dives matched to a designated site: %.1f%%", 100 * match_rate))
dives <- dives[!is.na(dives$site_row), ]

# --- Training data and folds ---------------------------------------------------
to_mga <- function(lon, lat) sf::st_coordinates(sf::st_transform(sf::st_as_sf(data.frame(lon, lat), coords = c("lon", "lat"), crs = 4326), crs_mga))
train_sites <- readRDS(file.path(cfg$exploration_dir, "artifacts", "site_year_all.rds")) |>
  dplyr::filter(!is.na(.data$high_density), is.finite(.data$longitude), is.finite(.data$latitude))
train_x <- terra::extract(stack, to_mga(train_sites$longitude, train_sites$latitude))
train_x$Year <- train_sites$year
reefs <- unique(c(train_sites$reef_name, dives$reef_name))
fold_of <- stats::setNames(sample(rep(seq_len(n_folds), length.out = length(reefs))), reefs)

# --- Leave-reef-out ranks of every designated site, per dived reef-year ------------
reef_years <- dives |> dplyr::distinct(.data$reef_name, .data$year) |>
  dplyr::mutate(fold = fold_of[.data$reef_name])
pts_by_site <- split(seq_along(site_pts$site), site_pts$site)
ranks <- list()
for (k in seq_len(n_folds)) {
  train <- fold_of[train_sites$reef_name] != k
  model <- xgboost::xgb.train(
    params = legacy_params,
    data = xgboost::xgb.DMatrix(as.matrix(train_x[train, features]), label = train_sites$high_density[train], missing = NA),
    nrounds = legacy_rounds, verbose = 0
  )
  ry <- reef_years[reef_years$fold == k, ]
  for (i in seq_len(nrow(ry))) {
    reef_rows <- sites$row[sites$reef_name == ry$reef_name[i]]
    if (length(reef_rows) < 3) next
    p_idx <- unlist(pts_by_site[as.character(reef_rows)], use.names = FALSE)
    x <- site_pts$values[p_idx, cfg$legacy$rg_layer_names]
    x$Year <- ry$year[i]
    p <- stats::predict(model, xgboost::xgb.DMatrix(as.matrix(x[, features]), missing = NA))
    score <- tapply(p, factor(site_pts$site[p_idx], levels = reef_rows), mean, na.rm = TRUE)
    dep <- sites$deployed[reef_rows]
    pct <- function(v) { ok <- is.finite(v); out <- rep(NA_real_, length(v)); out[ok] <- (rank(v[ok]) - 1) / max(sum(ok) - 1, 1); out }
    ranks[[length(ranks) + 1L]] <- data.frame(
      reef_name = ry$reef_name[i], year = ry$year[i], site_row = reef_rows,
      rank_cv = pct(as.numeric(score)), rank_deployed = pct(dep), reef_sites = length(reef_rows)
    )
  }
  log_message("Fold ", k, ": ", nrow(ry), " reef-years ranked")
}
ranks <- dplyr::bind_rows(ranks)

# Every designated site on a dived reef-year, with whether (and how) it was dived.
site_outcomes <- dives |>
  dplyr::group_by(.data$reef_name, .data$year, .data$site_row) |>
  dplyr::summarise(dives = dplyr::n(), problematic = sum(.data$outcome == "problematic"), .groups = "drop")
ever_dived <- unique(dives$site_row)
site_table <- ranks |>
  dplyr::left_join(site_outcomes, by = c("reef_name", "year", "site_row")) |>
  dplyr::mutate(dived_this_year = !is.na(.data$dives), ever_dived = .data$site_row %in% ever_dived)
saveRDS(site_table, file.path(cfg$dirs$cache, "cull_site_ranks.rds"))

scored <- dives |>
  dplyr::inner_join(ranks, by = c("reef_name", "year", "site_row")) |>
  dplyr::filter(is.finite(.data$rank_cv), is.finite(.data$rank_deployed)) |>
  dplyr::mutate(reef_year = paste(.data$reef_name, .data$year), hours = .data$bottom_time / 60)
log_message(format(nrow(scored), big.mark = ","), " dives scored in ", dplyr::n_distinct(scored$reef_year), " reef-years")

# --- Scenarios: dive only the top fraction of each reef's designated sites --------
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
period_sets <- list(`2018-2026` = scored$year >= 2018, `2022-2026` = scored$year >= 2022)
results <- dplyr::bind_rows(lapply(names(period_sets), function(period) {
  d <- scored[period_sets[[period]], ]
  groups <- split(seq_len(nrow(d)), d$reef_year)
  dplyr::bind_rows(lapply(c("rank_cv", "rank_deployed"), function(rank_col) {
    dplyr::bind_rows(lapply(site_fractions, function(a) {
      point <- metrics_for(d, rank_col, a)
      boot <- replicate(boot_reps, metrics_for(d[unlist(groups[sample.int(length(groups), replace = TRUE)]), ], rank_col, a))
      tibble::tibble(period = period, score = ifelse(rank_col == "rank_cv", "leave-reef-out", "deployed"),
                     sites_searched = a, metric = names(point), value = point,
                     low = apply(boot, 1, stats::quantile, 0.025), high = apply(boot, 1, stats::quantile, 0.975),
                     dives = nrow(d), reef_years = length(groups),
                     problematic_dives = sum(d$outcome == "problematic"), low_dives = sum(d$outcome == "low"),
                     ghost_dives = sum(d$outcome == "ghost"))
    }))
  }))
}))
write_table(results, cfg, "cull_site_selection_cv.csv")

# Within reef-year discrimination among dived sites, and the undived sites the
# layer would put forward.
auc <- dplyr::bind_rows(lapply(names(period_sets), function(period) {
  d <- scored[period_sets[[period]], ]
  dplyr::bind_rows(lapply(c(`leave-reef-out` = "rank_cv", deployed = "rank_deployed"), function(col) {
    pl <- d$outcome %in% c("problematic", "low")
    pg <- d$outcome %in% c("problematic", "ghost")
    tibble::tibble(period = period, score = ifelse(col == "rank_cv", "leave-reef-out", "deployed"),
                   auc_problematic_vs_low = within_group_auc(d[[col]][pl], d$outcome[pl] == "problematic", d$reef_year[pl]),
                   auc_problematic_vs_ghost = within_group_auc(d[[col]][pg], d$outcome[pg] == "problematic", d$reef_year[pg]))
  }))
}))
write_table(auc, cfg, "cull_site_selection_auc.csv")

undived <- site_table |>
  dplyr::filter(is.finite(.data$rank_cv)) |>
  dplyr::group_by(top30 = .data$rank_cv >= 0.7) |>
  dplyr::summarise(site_years = dplyr::n(), share_never_dived = mean(!.data$ever_dived),
                   share_not_dived_that_year = mean(!.data$dived_this_year), .groups = "drop")
coverage <- tibble::tibble(
  designated_sites = nrow(sites), reefs_with_sites = dplyr::n_distinct(sites$reef_name),
  sites_ever_dived = length(ever_dived), dive_match_rate = match_rate,
  reef_years_scored = dplyr::n_distinct(scored$reef_year), dives_scored = nrow(scored)
)
write_table(undived, cfg, "cull_site_selection_undived.csv")
write_table(coverage, cfg, "cull_site_selection_coverage.csv")
print(coverage); print(auc); print(undived)
print(tidyr::pivot_wider(results[results$score == "leave-reef-out" & results$period == "2018-2026", c("sites_searched", "metric", "value")],
                         names_from = "metric", values_from = "value"))
