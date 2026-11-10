# The proposed protocol (35) under the density-aware classes (20b): A-B vs
# A-C tow zone, and how much evidence a dive needs.
#
# Question: if the tow zone and the dive rule use density-aware classes (a
# sparse reef has fewer A-C sites, a dense reef more) instead of the
# within-reef classes (every reef 20 % A, 66 % A-C), how do tows saved, ghost
# and low-density dives avoided, and > 0.04 dives missed change? And how do
# the A-B and A-C tow zones compare, with dive rules of different strictness?
#
# Design: 35's visits, tows, dives and rules (34's skip rules, season memory,
# each tow counted once), test seasons 2020-2026 (20b's forward test starts in
# 2020). Each visit keeps its voyage-start order of sites (33, fwd_voyage);
# only how many sites reach each class changes.
#   within-reef   class k = top break_k of the reef's ranked sites (35)
#   density-aware class k = top plogis(qlogis(break_k) + delta) of the sites,
#                 delta = logit(reef density) - logit(typical culled reef),
#                 reef density from the reef's dives and tows in the 365 days
#                 before the visit's planning date (earlier visits that season
#                 count) and the outbreak forecast for the year before; 20b's
#                 evidence + forecast model and unknown-reef discount, refitted
#                 on seasons before T (strictly forward)
# Tow zone (follow-up visits only; first visits tow the whole reef): A-B (top
# 50 %) or A-C (top 66 %), open sites, and with season memory sites near a
# strong tow earlier that season.
# Dive rule, for dives not within 100 m of a kept strong tow; cut-off A-B, A-C
# or A-D:
#   discretion at A-x   weak, negative-only and no-tow dives only at A-x sites
#                       (35's dive rule)
#   weak evidence at A-x  a weak-scar tow is enough at A-x sites; dives with
#                       only negative tows or no tow nearby are not made
#   weak-scar dives only  only weak-scar dives are filtered by class; negative
#                       and no-tow dives stay at crew discretion
# Baselines: current practice, and for each dive rule alone the same number of
# dives skipped at random within its pool (exact expectation). Intervals: 300
# bootstrap resamples of reef-years, paired for density-aware minus
# within-reef differences.
# Caveats: recorded dives were directed by tows and crew choice, so dives the
# rules skip are only those crews made; the density at voyage start uses
# earlier visits' dives, which encode where crews chose to dive.
#
# Alternative reef densities (10 November), the same replay and statistics.
# 20b's density counts the 12 months evenly and is fitted on dived reef-years,
# so a reef culled down months ago, or towed clean and left, keeps the level
# of reefs crews chose to dive. Each alternative is fitted on earlier visits
# (seasons before T) with evidence at their planning dates:
#   12 m          share of a cull visit's dives that were problematic ~ the
#                 reef's dives and tows in the 365 days before planning +
#                 forecast; unknown reefs discounted as in 20b (control for
#                 the model form)
#   recent        the same, each dive and tow weighted 0.5^(age / 90 days)
#   no forecast   the same without the outbreak forecast
#   two-part      P(a survey occasion needs culling | evidence) x the share
#                 problematic if dived. Occasions are all of a cull-site reef's
#                 dive and tow records split at gaps > 10 days, so tow-only
#                 surveys where crews found nothing count as "not needed"; the
#                 dived-visit model never sees them.
# Each gives delta against the median over training visits. Paired
# differences against the within-reef classes and against 20b's density.
# Extra outputs: density_methods_accuracy.csv (does density rank visits'
# dives, pooled and within season, with a half-life sensitivity that was not
# used to choose 90 days), density_methods_2027_reefs.csv (each method's 2027
# density and A / A-C counts per reef, evidence to cfg evidence_end),
# density_layers_matched.csv (every layer at matched effort, with paired
# differences), density_layers_2027_sites.csv / _totals.csv (every layer's
# 2027 classes; 37 writes the density-aware v2 product from them).

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "..", "R", "common.R"))
cfg <- initialise_final()
source(file.path(cfg$final_dir, "R", "density_classes.R"))
sr <- cfg$site_ranking
set.seed(cfg$hmsc$seed + 36L)
boot_reps <- as.integer(Sys.getenv("COTS_BOOT_REPS", "300"))   # fewer for exploratory runs
min_skip <- 3L
segment_break_m <- 1000
window_days <- 365
test_years <- 2020:2026
zone_cuts <- c(`A-B` = 0.5, `A-C` = 0.66)
rule_cuts <- c(`A-B` = 0.5, `A-C` = 0.66, `A-D` = 0.8)
evidence_rules <- c("discretion at A-x", "weak evidence at A-x", "weak-scar dives only")
weak_class <- "weak_single_or_present_scar"

# --- Inputs: 34's visits and zone ranks, 33's voyage ranking, 20b's training data ---
ctx <- readRDS(file.path(cfg$dirs$cache, "manta_zone_context.rds"))
tows <- ctx$tows
vis <- ctx$vis
guide <- ctx$guide
tow_sites <- ctx$tow_sites
vis$follow_up <- ctx$visit_info$follow_up[match(vis$visit, ctx$visit_info$visit)]
scored <- as.data.frame(readRDS(file.path(cfg$dirs$cache, "deployment_validation_scored.rds")))
dens_in <- readRDS(file.path(cfg$dirs$cache, "site_ranking_density_inputs.rds"))
forecast <- read_reef_forecast(cfg)
sv <- ctx$sv
sv$ranked <- sv$zone_combined >= 0
sv$n_ranked <- stats::ave(as.numeric(sv$ranked), sv$visit, FUN = sum)
# Top fraction (rank - 0.5) / n among the visit's ranked sites; unranked last.
top_fraction <- function(pct, n) ifelse(pct >= 0, (rank_of(pct, n) - 0.5) / n, Inf)
sv$top <- top_fraction(sv$zone_combined, sv$n_ranked)

# --- Reef density at each visit's planning date -------------------------------------
unknown_fit <- function(max_year) stats::glm(needed ~ forecast, family = stats::binomial(),
                                             data = dens_in$surveyed[dens_in$surveyed$year <= max_year, ])
dives_by_reef <- split(dens_in$dives, dens_in$dives$reef_name)
tows_by_reef <- split(dens_in$tows, dens_in$tows$reef_name)
in_window <- function(x, end) if (is.null(x)) x else x[x$survey_date < end & x$survey_date >= end - window_days, ]
vis_t <- vis[vis$year %in% test_years, ]
ev <- dplyr::bind_rows(lapply(seq_len(nrow(vis_t)), function(i) {
  r <- vis_t$reef_name[i]
  d <- in_window(dives_by_reef[[r]], vis_t$plan[i])
  t <- in_window(tows_by_reef[[r]], vis_t$plan[i])
  reef_evidence(r, if (is.null(d)) dens_in$dives[0, ] else d, if (is.null(t)) dens_in$tows[0, ] else t,
                vis_t$year[i] - sr$density_forecast_lag, forecast)
}))
ev <- dplyr::bind_cols(ev, reef_terms(ev, dens_in$forecast_median))
ev$visit <- vis_t$visit
ev$year <- vis_t$year
ev$sites_on_reef <- sv$n_ranked[match(ev$visit, sv$visit)]
ev <- ev[is.finite(ev$sites_on_reef) & ev$sites_on_reef > 0, ]
ev$delta <- NA_real_
ev$density <- NA_real_
for (T in test_years) {
  tr <- dens_in$dived[dens_in$dived$year < T, ]
  model <- stats::glm(density_variants$evidence_forecast, family = stats::binomial(), data = tr)
  um <- unknown_fit(T - 1)
  tr_grid <- site_grid(dens_in$reef_years[dens_in$reef_years$year < T & dens_in$reef_years$reef_year %in% tr$reef_year, ])
  ref <- stats::median(reef_density(discount_unknown(stats::predict(model, tr_grid), tr_grid, um), tr_grid$reef_year))
  i <- which(ev$year == T)
  g <- site_grid(ev[i, ])
  dens <- reef_density(discount_unknown(stats::predict(model, g), g, um), g$visit)
  ev$density[i] <- dens[ev$visit[i]]
  ev$delta[i] <- stats::qlogis(ev$density[i]) - stats::qlogis(ref)
  log_message(T, ": ", length(i), " visits, typical reef density ", round(ref, 3), ", median visit density ",
              round(stats::median(ev$density[i]), 3))
}
sv$delta <- ev$delta[match(sv$visit, ev$visit)]
vis$delta <- ev$delta[match(vis$visit, ev$visit)]

# --- Alternative reef densities (header), at every visit's planning date --------------
density_methods <- list(
  `12 m` = list(half_life = Inf, forecast = TRUE, two_part = FALSE),
  `recent` = list(half_life = 90, forecast = TRUE, two_part = FALSE),
  `recent, no forecast` = list(half_life = 90, forecast = FALSE, two_part = FALSE),
  `recent, two-part` = list(half_life = 90, forecast = TRUE, two_part = TRUE),
  `recent, two-part, no forecast` = list(half_life = 90, forecast = FALSE, two_part = TRUE))
half_life_check <- c(30, 60, 90, 180, Inf)
source(file.path(cfg$final_dir, "R", "cull_sites.R"))
cull_reefs <- unique(read_designated_sites(cfg)$reef_name)
occ <- survey_occasions(dens_in$dives, dens_in$tows, cull_reefs)
occ <- occ[occ$year >= 2013 & occ$year <= max(test_years), ]
vo <- dplyr::summarise(dplyr::group_by(ctx$dives, .data$visit), n = dplyr::n(), k = sum(.data$outcome == "problematic"), .groups = "drop")
vis$n_dives <- dplyr::coalesce(vo$n[match(vis$visit, vo$visit)], 0L)
vis$k_dives <- dplyr::coalesce(vo$k[match(vis$visit, vo$visit)], 0L)
evidence_terms <- function(reef, end, year, h) {
  e <- weighted_evidence(reef, end, year - sr$density_forecast_lag, dives_by_reef, tows_by_reef, h, forecast, window_days)
  dplyr::bind_cols(e[, c("dives", "problematic", "tows", "strong")], reef_terms(e, dens_in$forecast_median))
}
hl_all <- sort(unique(c(vapply(density_methods, `[[`, numeric(1), "half_life"), half_life_check)))
ve <- lapply(stats::setNames(hl_all, hl_all), function(h) evidence_terms(vis$reef_name, vis$plan, vis$year, h))
oe <- lapply(stats::setNames(hl_all, hl_all), function(h) cbind(evidence_terms(occ$reef_name, occ$start, occ$year, h),
                                                                needed = occ$needed))
reef_rhs <- function(fc) paste("has_dives + dive_rate + log_dives + has_tows + tow_rate", if (fc) "+ forecast" else "")
# Density of `new` (evidence terms) for method `m`, fitted on visits and
# occasions before season T; delta against the median training cull visit.
method_density <- function(m, T, new) {
  key <- as.character(m$half_life)
  tr <- cbind(ve[[key]], k = vis$k_dives, n = vis$n_dives)[vis$year < T & vis$year >= 2019 & vis$n_dives > 0, ]
  ms <- stats::glm(stats::as.formula(paste("cbind(k, n - k) ~", reef_rhs(m$forecast))), family = stats::binomial(), data = tr)
  p <- stats::predict(ms, new, type = "response")
  ptr <- stats::predict(ms, tr, type = "response")
  if (m$two_part) {
    mn <- stats::glm(stats::as.formula(paste("needed ~", reef_rhs(m$forecast))), family = stats::binomial(),
                     data = oe[[key]][occ$year < T, ])
    p <- p * stats::predict(mn, new, type = "response")
    ptr <- ptr * stats::predict(mn, tr, type = "response")
  } else {
    um <- stats::glm(if (m$forecast) needed ~ forecast else needed ~ 1, family = stats::binomial(),
                     data = dens_in$surveyed[dens_in$surveyed$year <= T - 1, ])
    disc <- function(p, d) {
      u <- d$has_dives == 0 & d$has_tows == 0
      p[u] <- p[u] * stats::predict(um, d[u, , drop = FALSE], type = "response")
      p
    }
    p <- disc(p, new)
    ptr <- disc(ptr, tr)
  }
  clamp <- function(x) pmin(pmax(x, 1e-4), 1 - 1e-4)
  list(density = p, delta = stats::qlogis(clamp(p)) - stats::qlogis(stats::median(clamp(ptr))))
}
vis_delta <- list(`density-aware` = vis$delta)
vis_density <- list(`density-aware` = ev$density[match(vis$visit, ev$visit)])
for (nm in names(density_methods)) {
  vis_delta[[nm]] <- vis_density[[nm]] <- rep(NA_real_, nrow(vis))
  for (T in test_years) {
    i <- which(vis$year == T)
    r <- method_density(density_methods[[nm]], T, ve[[as.character(density_methods[[nm]]$half_life)]][i, ])
    vis_delta[[nm]][i] <- r$delta
    vis_density[[nm]][i] <- r$density
  }
}
layers <- c("within-reef", "density-aware")

# --- Keeping the sparse parts of dense reefs low (header) ----------------------------
# The scaled quota lifts every site of a dense reef together, so parts of the
# reef that rarely hold COTS rise with it. Three ways to keep the lift off
# them, on the recent two-part density:
#   damped lift          delta > 0 halved (sparse reefs still shrink in full)
#   no lift where clean  sites towed in the last year with no strong tow and
#                        no problematic dive within 300 m keep delta <= 0
#   site model           P(dive problematic) on the reef's logit density, the
#                        site's logit rank on the reef and its own history
#                        (dives within 300 m and linked tows: recent, weighted
#                        by a 90-day half-life within a year, and 1-6 years
#                        before), fitted on earlier seasons' dives
#   site visit model     the same terms plus log dive and tow counts, fitted
#                        on every designated site of earlier visits, outcome
#                        "a problematic dive at this site on this visit": sites
#                        crews towed past and never dived count as zeros, so
#                        a reef's long-unvisited parts can stay low. This
#                        encodes where crews chose to dive as well as COTS.
#   site + visit model   site model + 0.5 x site visit model (log odds): a
#                        site needs a good chance of a problematic dive if
#                        dived, and a record like sites that got problematic
#                        dives counts for half as much. The weight was chosen
#                        on these test seasons (exploratory grid 0.5-5, 11
#                        November): weights above about 2 kept long-clean
#                        flanks low but lost misses and tow-zone performance
#   quota, site + visit order
#                        the recent two-part density's scaled quota (how many
#                        sites per class on each reef), with the reef's sites
#                        ordered by the site + visit score instead of the
#                        voyage rank: between-reef allocation stays with the
#                        reef density, so reefs known only from the forecast
#                        cannot gain sites from the site models
#   density-aware v2     site + visit model on reefs with any dives or tows
#                        in the year before planning (the replay's reefs, where
#                        it is tested), the quota version on reefs without
#                        (known only from the forecast, where the site models
#                        would give about 60 % of sites A-C against 27 %
#                        under the quota). Quota scores are mapped onto the
#                        site + visit scale through the two layers' class
#                        thresholds, so each keeps its own classes.
#   v2, untowed neutral  v2 with one change (10 November): a site with no
#                        linked tow in the past year was scored with a strong-
#                        tow rate of 50 % (logit 0), far above a typical towed
#                        site, so untowed lagoon patches outranked towed ground
#                        with strong scars (Tongue, John Brewer). It now gets the
#                        mean recent strong-tow rate of the visit's towed sites,
#                        as if towed with the reef's average result; its 1-6
#                        year record still counts. Shifting untowed scores to
#                        the towed sites' median was tried first and rejected:
#                        on dense reefs it lifted ground crews skip because it
#                        is usually clean (Batt's NE flank to A).
#   + habitat            the untowed-neutral v2 with the site's habitat
#                        percentile on its reef entered directly in both site
#                        models (centred), plus habitat x untowed so habitat
#                        can count for more where there are no tows (no
#                        untowed main effect, so untowed ground is not lifted
#                        again). Replay: 33's forward cull-habitat layer for
#                        each season (trained on earlier seasons); 2027: the
#                        deployed cull-habitat rank from 20. The manta-trained
#                        layer was fitted on all tows and would leak here.
#                        Extra output: density_untowed_auc.csv (within-visit
#                        AUC among dives at untowed sites).
#   v3, local tow density  the untowed-neutral v2 plus a local tow-density
#                        term in both site models: the site's strong-tow rate
#                        pooled over its reef's designated sites with a 500 m
#                        Gaussian kernel (recent tows in full, 1-6 year tows at
#                        a quarter per year), shrunk to the reef's rate with 2
#                        pseudo-tows, entered as the logit difference from the
#                        reef's rate. One reef-wide density let a dense flank
#                        lift the whole reef (Batt). "+ aspect" also weights
#                        the kernel by ((1 + cos(aspect difference)) / 2)^4,
#                        aspect = the direction the nearest reef edge faces,
#                        so opposite flanks of a narrow reef share little. Width
#                        and terms set on a 2020-2023 check (local tows beat
#                        the site's own 100 m rate, AUC 0.62 vs 0.58; aspect
#                        did not help on average); 2024-2026 confirmed it.
# Every alternative layer (the five densities and these three) is a site
# score: delta - qlogis(top fraction), or the site model's log odds. A site
# is in class k when its score reaches the layer's threshold for break_k, set
# on the training seasons' designated sites so that the same pooled share of
# sites reaches each class as under the within-reef classes. The layers then
# differ only in where that effort goes, not in how much there is.
base_layer <- "recent, two-part"
sites_sf <- read_designated_sites(cfg)
sites_xy <- sf::st_coordinates(suppressWarnings(sf::st_centroid(sf::st_geometry(sites_sf))))
all_dives <- matched_cull_dives(cfg, sites_sf)
site_raw <- site_history(vis, sv, sites_xy, sites_sf$reef_name, all_dives, tows, ctx$tow_sites)
site_feat <- site_terms(site_raw)
# Batt Reef's north-eastern and south-western flanks (either side of its long
# axis), for the diagnostics below.
batt <- which(grepl("16-029", sites_sf$reef_name))
batt_pc <- stats::prcomp(sites_xy[batt, ])
batt_perp <- c(-batt_pc$rotation[2, 1], batt_pc$rotation[1, 1])
if (batt_perp[2] < 0) batt_perp <- -batt_perp
batt_ne <- batt[as.numeric(scale(sites_xy[batt, ], batt_pc$center, FALSE) %*% batt_perp) > 0]
neg_q <- function(t) ifelse(is.finite(t), -stats::qlogis(pmin(t, 1 - 1e-6)), -Inf)
rank_lo <- function(t) ifelse(is.finite(t), stats::qlogis(pmin(t, 1 - 1e-6)), NA_real_)
sv_vis <- match(sv$visit, vis$visit)
sv_key <- paste(sv$visit, sv$site_row)
damp <- function(d) ifelse(d > 0, 0.5 * d, d)
shares <- sort(unique(c(zone_cuts, rule_cuts, sr$priority_breaks)))
# Score functions of (delta, top fraction, site features), and the density
# each takes its delta from.
score_fns <- c(lapply(stats::setNames(nm = names(density_methods)), function(nm) function(delta, top, feat) delta + neg_q(top)),
               list(`two-part, damped lift` = function(delta, top, feat) damp(delta) + neg_q(top),
                    `two-part, no lift where searched clean` = function(delta, top, feat)
                      ifelse(feat$searched_clean, pmin(delta, 0), delta) + neg_q(top)))
score_density <- c(stats::setNames(names(density_methods), names(density_methods)),
                   `two-part, damped lift` = base_layer, `two-part, no lift where searched clean` = base_layer)
visit_weight <- 0.5
score_layers <- c(names(score_fns), "two-part, site model", "two-part, site visit model", "two-part, site + visit model", "two-part quota, site + visit order", "density-aware v2", "density-aware v2, untowed neutral", "density-aware v2, untowed neutral + habitat", "density-aware v3, local tow density", "density-aware v3, local tow density + aspect")
# Top fraction within each group (visit or reef) of a score; -Inf ranks last.
top_within <- function(score, g) {
  r <- stats::ave(-score, g, FUN = function(x) rank(x, ties.method = "first"))
  n <- stats::ave(as.numeric(is.finite(score)), g, FUN = sum)
  ifelse(is.finite(score), (r - 0.5) / n, Inf)
}
sv_tq <- rep(Inf, nrow(sv))
# Untowed sites (no linked tow in the past year) get the mean recent strong-tow
# rate of the towed sites in their group (visit, or reef for 2027).
impute_tows <- function(feat, grp) {
  untowed <- feat$has_tr == 0
  fill <- tapply(feat$rate_tr[!untowed], grp[!untowed], mean)[as.character(grp)]
  i <- untowed & !is.na(fill)
  feat$has_tr[i] <- 1
  feat$rate_tr[i] <- fill[i]
  feat
}
site_feat_n <- impute_tows(site_feat, sv_vis)
# Habitat percentile of each site on its reef in the visit's season (33's
# forward cull layer), centred; missing = reef median (0).
source(file.path(cfg$final_dir, "R", "site_selection.R"))
hab_fwd <- readRDS(file.path(cfg$dirs$cache, "deployment_forward_habitat.rds")) |>
  dplyr::group_by(.data$reef_name, .data$year) |>
  dplyr::mutate(rank_fwd = pct_rank(.data$habitat_fwd)) |>
  dplyr::ungroup()
sv_hab <- hab_fwd$rank_fwd[match(paste(vis$reef_name[sv_vis], vis$year[sv_vis], sv$site_row),
                                 paste(hab_fwd$reef_name, hab_fwd$year, hab_fwd$site_row))] - 0.5
log_message("Forward habitat for ", round(100 * mean(!is.na(sv_hab[vis$year[sv_vis] %in% test_years]))), " % of test site-visits")
sv_hab[is.na(sv_hab)] <- 0
sv_untowed <- as.numeric(site_feat$has_tr == 0)
hab_rhs <- "+ hab + hab:untowed"
with_hab <- function(f, rows) cbind(f, hab = sv_hab[rows], untowed = sv_untowed[rows])
# v3: aspect of each designated site (direction its nearest reef edge faces,
# from the GBR reef outlines) and local tow density relative to the reef.
reef_outlines <- sf::st_read(checked_input(cfg, "reef_features"), quiet = TRUE) |> sf::st_transform(sf::st_crs(sites_sf))
site_lab <- reef_id(sites_sf$reef_name)
site_aspect <- rep(NA_real_, nrow(sites_sf))
for (lab in unique(site_lab)) {
  poly <- reef_outlines[reef_outlines$LABEL_ID == lab, ]
  if (!nrow(poly)) next
  poly <- sf::st_union(sf::st_make_valid(poly))
  i <- which(site_lab == lab)
  pts <- sf::st_sfc(lapply(i, function(k) sf::st_point(sites_xy[k, ])), crs = sf::st_crs(sites_sf))
  nb <- sf::st_nearest_points(pts, sf::st_cast(sf::st_boundary(poly), "MULTILINESTRING"))
  near <- t(vapply(seq_along(nb), function(k) sf::st_coordinates(nb[k])[2, 1:2], numeric(2)))
  v <- near - sites_xy[i, , drop = FALSE]
  outside <- !as.logical(sf::st_intersects(pts, poly, sparse = FALSE)[, 1])
  v[outside, ] <- -v[outside, ]
  site_aspect[i] <- atan2(v[, 2], v[, 1])
}
local_dev <- function(raw, grp, site_rows, k, h = 500, m = 2) {
  out <- rep(0, length(site_rows))
  t_all <- raw$t_recent + 0.25 * raw$t_long / 5
  s_all <- raw$s_recent + 0.25 * raw$s_long / 5
  for (g in split(seq_along(site_rows), grp)) {
    g <- g[!is.na(site_rows[g])]
    if (!length(g)) next
    K <- exp(-as.matrix(stats::dist(sites_xy[site_rows[g], , drop = FALSE]))^2 / (2 * h^2))
    if (k > 0) {
      a <- site_aspect[site_rows[g]]
      a[is.na(a)] <- 0
      K <- K * ((1 + cos(outer(a, a, "-"))) / 2)^k
    }
    reef <- (sum(s_all[g]) + 0.5) / (sum(t_all[g]) + 1)
    out[g] <- stats::qlogis((K %*% s_all[g] + m * reef) / (K %*% t_all[g] + m)) - stats::qlogis(reef)
  }
  out
}
local_layers <- c(k0 = "density-aware v3, local tow density", k4 = "density-aware v3, local tow density + aspect")
local_k <- c(k0 = 0, k4 = 4)
sv_local <- lapply(local_k, function(k) local_dev(site_raw, sv_vis, sv$site_row, k))
with_local <- function(f, rows, v) cbind(f, local_dev = sv_local[[v]][rows])
local_models <- list(k0 = list(), k4 = list())
sv_scores <- lapply(stats::setNames(nm = score_layers), function(l) rep(-Inf, nrow(sv)))
thr_mat <- lapply(stats::setNames(nm = score_layers), function(l)
  matrix(NA_real_, length(test_years), length(shares), dimnames = list(test_years, format(shares))))
calibrate <- function(score) {
  v <- sort(score, decreasing = TRUE)
  vapply(shares, function(a) v[max(1L, ceiling(a * length(v)))], numeric(1))
}
site_rhs <- "reef_lp + rank_lo + has_dr + rate_dr + has_dl + rate_dl + has_tr + rate_tr + has_tl + rate_tl"
site_frame <- function(rows, top, reef_lp, feat = site_feat) cbind(feat[rows, ], rank_lo = rank_lo(top), reef_lp = reef_lp)
visit_rhs <- paste(site_rhs, "+ log_dl + log_tl + log_dr")
prob_key <- unique(paste(ctx$dives$visit, ctx$dives$site_row)[ctx$dives$outcome == "problematic"])
sv_prob <- as.numeric(sv_key %in% prob_key)
site_models <- list()
visit_models <- list()
hab_models <- list()
for (T in test_years) {
  tr_rows <- which(vis$year[sv_vis] < T & vis$year[sv_vis] >= 2019)
  te_rows <- which(vis$year[sv_vis] == T)
  delta_T <- lapply(stats::setNames(nm = names(density_methods)), function(nm)
    method_density(density_methods[[nm]], T, ve[[as.character(density_methods[[nm]]$half_life)]])$delta)
  for (l in names(score_fns)) {
    d <- delta_T[[score_density[[l]]]][sv_vis]
    thr_mat[[l]][as.character(T), ] <- calibrate(score_fns[[l]](d[tr_rows], sv$top[tr_rows], site_feat[tr_rows, ]))
    sv_scores[[l]][te_rows] <- score_fns[[l]](d[te_rows], sv$top[te_rows], site_feat[te_rows, ])
  }
  lp_T <- delta_T[[base_layer]]   # logit density up to a constant per season
  tr_d <- ctx$dives[ctx$dives$visit %in% vis$visit[vis$year < T & vis$year >= 2019], ]
  j <- match(paste(tr_d$visit, tr_d$site_row), sv_key)
  tr_d <- tr_d[!is.na(j) & is.finite(sv$top[j]), ]
  j <- match(paste(tr_d$visit, tr_d$site_row), sv_key)
  fr <- site_frame(j, sv$top[j], lp_T[sv_vis[j]])
  fr$problematic <- as.numeric(tr_d$outcome == "problematic")
  site_models[[as.character(T)]] <- m <- stats::glm(stats::as.formula(paste("problematic ~", site_rhs)),
                                                     family = stats::binomial(), data = fr)
  score_rows <- function(rows) {
    out <- rep(-Inf, length(rows))
    ok <- is.finite(sv$top[rows])
    out[ok] <- stats::predict(m, site_frame(rows[ok], sv$top[rows[ok]], lp_T[sv_vis[rows[ok]]]))
    out
  }
  thr_mat[["two-part, site model"]][as.character(T), ] <- calibrate(score_rows(tr_rows))
  sv_scores[["two-part, site model"]][te_rows] <- score_rows(te_rows)
  attr(site_models[[as.character(T)]], "lp") <- lp_T
  ok_tr <- tr_rows[is.finite(sv$top[tr_rows])]
  fv <- site_frame(ok_tr, sv$top[ok_tr], lp_T[sv_vis[ok_tr]])
  fv$problematic <- sv_prob[ok_tr]
  visit_models[[as.character(T)]] <- mv <- stats::glm(stats::as.formula(paste("problematic ~", visit_rhs)),
                                                       family = stats::binomial(), data = fv)
  attr(visit_models[[as.character(T)]], "lp") <- lp_T
  score_v <- function(rows) {
    out <- rep(-Inf, length(rows))
    ok <- is.finite(sv$top[rows])
    out[ok] <- stats::predict(mv, site_frame(rows[ok], sv$top[rows[ok]], lp_T[sv_vis[rows[ok]]]))
    out
  }
  thr_mat[["two-part, site visit model"]][as.character(T), ] <- calibrate(score_v(tr_rows))
  sv_scores[["two-part, site visit model"]][te_rows] <- score_v(te_rows)
  comb_tr <- score_rows(tr_rows) + visit_weight * score_v(tr_rows)
  comb_te <- score_rows(te_rows) + visit_weight * score_v(te_rows)
  thr_mat[["two-part, site + visit model"]][as.character(T), ] <- calibrate(comb_tr)
  sv_scores[["two-part, site + visit model"]][te_rows] <- comb_te
  thr_mat[["two-part quota, site + visit order"]][as.character(T), ] <- calibrate(delta_T[[base_layer]][sv_vis[tr_rows]] + neg_q(top_within(comb_tr, sv_vis[tr_rows])))
  sv_tq[te_rows] <- top_within(comb_te, sv_vis[te_rows])
  sv_scores[["two-part quota, site + visit order"]][te_rows] <- delta_T[[base_layer]][sv_vis[te_rows]] + neg_q(sv_tq[te_rows])
  comb_n <- function(rows) {
    out <- rep(-Inf, length(rows))
    ok <- is.finite(sv$top[rows])
    f <- site_frame(rows[ok], sv$top[rows[ok]], lp_T[sv_vis[rows[ok]]], site_feat_n)
    out[ok] <- stats::predict(m, f) + visit_weight * stats::predict(mv, f)
    out
  }
  thr_mat[["density-aware v2, untowed neutral"]][as.character(T), ] <- calibrate(comb_n(tr_rows))
  sv_scores[["density-aware v2, untowed neutral"]][te_rows] <- comb_n(te_rows)
  m_h <- stats::glm(stats::as.formula(paste("problematic ~", site_rhs, hab_rhs)), family = stats::binomial(),
                    data = with_hab(fr, j))
  mv_h <- stats::glm(stats::as.formula(paste("problematic ~", visit_rhs, hab_rhs)), family = stats::binomial(),
                     data = with_hab(fv, ok_tr))
  hab_models[[as.character(T)]] <- list(dive = m_h, visit = mv_h)
  comb_h <- function(rows) {
    out <- rep(-Inf, length(rows))
    ok <- is.finite(sv$top[rows])
    f <- with_hab(site_frame(rows[ok], sv$top[rows[ok]], lp_T[sv_vis[rows[ok]]], site_feat_n), rows[ok])
    out[ok] <- stats::predict(m_h, f) + visit_weight * stats::predict(mv_h, f)
    out
  }
  thr_mat[["density-aware v2, untowed neutral + habitat"]][as.character(T), ] <- calibrate(comb_h(tr_rows))
  sv_scores[["density-aware v2, untowed neutral + habitat"]][te_rows] <- comb_h(te_rows)
  for (v in names(local_layers)) {
    m_l <- stats::glm(stats::as.formula(paste("problematic ~", site_rhs, "+ local_dev")), family = stats::binomial(),
                      data = with_local(fr, j, v))
    mv_l <- stats::glm(stats::as.formula(paste("problematic ~", visit_rhs, "+ local_dev")), family = stats::binomial(),
                       data = with_local(fv, ok_tr, v))
    local_models[[v]][[as.character(T)]] <- list(dive = m_l, visit = mv_l)
    comb_l <- function(rows) {
      out <- rep(-Inf, length(rows))
      ok <- is.finite(sv$top[rows])
      f <- with_local(site_frame(rows[ok], sv$top[rows[ok]], lp_T[sv_vis[rows[ok]]], site_feat_n), rows[ok], v)
      out[ok] <- stats::predict(m_l, f) + visit_weight * stats::predict(mv_l, f)
      out
    }
    thr_mat[[local_layers[[v]]]][as.character(T), ] <- calibrate(comb_l(tr_rows))
    sv_scores[[local_layers[[v]]]][te_rows] <- comb_l(te_rows)
  }
}
# density-aware v2: site + visit where the reef has recent evidence, else quota.
ev_key <- as.character(density_methods[[base_layer]]$half_life)
vis_evidence <- ve[[ev_key]]$has_dives > 0 | ve[[ev_key]]$has_tows > 0
# Monotone map of quota scores onto the site + visit scale, exact at the
# class thresholds of season T (linear beyond them).
map_quota <- function(sq, year, target = "two-part, site + visit model") {
  out <- sq
  for (T in unique(year[is.finite(sq)])) {
    i <- which(year == T & is.finite(sq))
    x <- thr_mat[["two-part quota, site + visit order"]][as.character(T), ]
    y <- thr_mat[[target]][as.character(T), ]
    o <- order(x)
    x <- x[o]
    y <- y[o]
    v <- stats::approx(x, y, xout = sq[i], rule = 2)$y
    lo <- sq[i] < x[1]
    hi <- sq[i] > x[length(x)]
    v[lo] <- y[1] + sq[i][lo] - x[1]
    v[hi] <- y[length(y)] + sq[i][hi] - x[length(x)]
    out[i] <- v
  }
  out
}
thr_mat[["density-aware v2"]] <- thr_mat[["two-part, site + visit model"]]
te_all <- which(vis$year[sv_vis] %in% test_years)
sv_scores[["density-aware v2"]][te_all] <- ifelse(vis_evidence[sv_vis[te_all]], sv_scores[["two-part, site + visit model"]][te_all],
                                    map_quota(sv_scores[["two-part quota, site + visit order"]][te_all], vis$year[sv_vis[te_all]]))
sv_scores[["density-aware v2, untowed neutral"]][te_all] <- ifelse(vis_evidence[sv_vis[te_all]], sv_scores[["density-aware v2, untowed neutral"]][te_all],
                                    map_quota(sv_scores[["two-part quota, site + visit order"]][te_all], vis$year[sv_vis[te_all]], "density-aware v2, untowed neutral"))
sv_scores[["density-aware v2, untowed neutral + habitat"]][te_all] <- ifelse(vis_evidence[sv_vis[te_all]], sv_scores[["density-aware v2, untowed neutral + habitat"]][te_all],
                                    map_quota(sv_scores[["two-part quota, site + visit order"]][te_all], vis$year[sv_vis[te_all]], "density-aware v2, untowed neutral + habitat"))
for (L in local_layers) sv_scores[[L]][te_all] <- ifelse(vis_evidence[sv_vis[te_all]], sv_scores[[L]][te_all],
                                                         map_quota(sv_scores[["two-part quota, site + visit order"]][te_all], vis$year[sv_vis[te_all]], L))
site_thr <- lapply(thr_mat, function(tm) function(share, year)
  tm[cbind(match(year, test_years), match(round(share, 6), round(shares, 6)))])
coef_site <- dplyr::bind_rows(lapply(c("dives", "site visits", "dives + habitat", "site visits + habitat", "dives + local tows",
                                         "site visits + local tows", "dives + local tows, aspect", "site visits + local tows, aspect"), function(kind) {
  ms <- switch(kind, dives = site_models, `site visits` = visit_models,
               `dives + habitat` = lapply(hab_models, `[[`, "dive"), `site visits + habitat` = lapply(hab_models, `[[`, "visit"),
               `dives + local tows` = lapply(local_models$k0, `[[`, "dive"), `site visits + local tows` = lapply(local_models$k0, `[[`, "visit"),
               `dives + local tows, aspect` = lapply(local_models$k4, `[[`, "dive"),
               `site visits + local tows, aspect` = lapply(local_models$k4, `[[`, "visit"))
  dplyr::bind_rows(lapply(names(ms), function(T) {
    m <- ms[[T]]
    tibble::tibble(model = kind, fitted_before = as.integer(T), term = names(stats::coef(m)), estimate = unname(stats::coef(m)),
                   std_error = unname(sqrt(diag(stats::vcov(m)))), training_rows = stats::nobs(m))
  }))
}))
write_table(coef_site, cfg, "density_site_model_coefficients.csv")
print(coef_site[coef_site$fitted_before == max(test_years), ], n = Inf)
layers <- c(layers, score_layers)
write_table(dplyr::bind_rows(lapply(names(thr_mat), function(l) tibble::tibble(
  layer = l, fitted_before = test_years, as.data.frame(thr_mat[[l]], check.names = FALSE)))), cfg, "density_layer_thresholds.csv")

# Does density rank dives between visits? Dive-level AUC of "problematic" on
# the visit's density, pooled and within season (test visits with dives).
auc <- function(s, y) {
  ok <- is.finite(s)
  s <- s[ok]
  y <- y[ok]
  if (length(unique(y)) < 2) return(NA_real_)
  (sum(rank(s)[y == 1]) - sum(y == 1) * (sum(y == 1) + 1) / 2) / (sum(y == 1) * sum(y == 0))
}
acc_rows <- function(method, dens) {
  di <- match(ctx$dives$visit, vis$visit)
  sel <- !is.na(di) & vis$year[di] %in% test_years
  y <- as.numeric(ctx$dives$outcome[sel] == "problematic")
  s <- dens[di[sel]]
  yr <- vis$year[di[sel]]
  fu <- ctx$visit_info$follow_up[match(ctx$dives$visit[sel], ctx$visit_info$visit)]
  tibble::tibble(method = method, dives = length(y), auc_pooled = auc(s, y),
                 auc_within_season = mean(vapply(split(seq_along(y), yr), function(j) auc(s[j], y[j]), numeric(1)), na.rm = TRUE),
                 auc_follow_up = auc(s[fu %in% TRUE], y[fu %in% TRUE]), auc_first_visit = auc(s[fu %in% FALSE], y[fu %in% FALSE]),
                 high_auc_pooled = auc(s, as.numeric(ctx$dives$cpue[sel] > 0.04)))
}
accuracy <- dplyr::bind_rows(lapply(names(vis_density), function(nm) acc_rows(nm, vis_density[[nm]])))
for (two in c(FALSE, TRUE)) for (h in half_life_check) {
  m <- list(half_life = h, forecast = TRUE, two_part = two)
  dens <- rep(NA_real_, nrow(vis))
  for (T in test_years) {
    i <- which(vis$year == T)
    dens[i] <- method_density(m, T, ve[[as.character(h)]][i, ])$density
  }
  accuracy <- dplyr::bind_rows(accuracy, acc_rows(paste0("sensitivity: half-life ", h, " d", if (two) ", two-part" else ""), dens))
}
write_table(accuracy, cfg, "density_methods_accuracy.csv")
print(accuracy, width = Inf)

# 2027: every layer's classes for every designated site in 20b's product
# (evidence to cfg evidence_end; models and thresholds fitted on every visit
# and occasion to 2026, as in the replay), and per-reef densities and counts.
p27 <- utils::read.csv(file.path(cfg$dirs$outputs, "site_ranking", "cots_cull_site_priority_density_2027.csv"),
                       stringsAsFactors = FALSE)
T27 <- sr$target_season
r27 <- unique(p27$reef_name)
pv <- data.frame(visit = paste(T27, r27), reef_name = r27, plan = sr$evidence_end + 1, year = T27)
p27$visit <- paste(T27, p27$reef_name)
p27$site_row <- match(paste(p27$reef_name, clean_site_name(p27$site_name), sep = "||"), sites_sf$site_key)
p27$top <- (p27$rank_on_reef - 0.5) / p27$sites_on_reef
raw27 <- site_history(pv, p27[, c("visit", "site_row")], sites_xy, sites_sf$reef_name, all_dives, tows, ctx$tow_sites)
feat27 <- site_terms(raw27)
p_vis <- match(p27$reef_name, r27)
all_rows <- which(vis$year[sv_vis] >= 2019)
class_of <- function(score, thr) {
  b <- thr[match(round(sr$priority_breaks, 6), round(shares, 6))]
  labels27[1L + rowSums(outer(score, b, "<"))]
}
labels27 <- c(names(sr$priority_breaks), LETTERS[length(sr$priority_breaks) + 1L])
cls27 <- list(`within-reef` = p27$priority_within_reef, `density-aware` = p27$priority_density)
dens27 <- list(`density-aware` = p27$reef_expected_problematic[match(r27, p27$reef_name)])
delta_tr27 <- delta_27 <- list()
for (nm in names(density_methods)) {
  m <- density_methods[[nm]]
  delta_tr27[[nm]] <- method_density(m, T27, ve[[as.character(m$half_life)]])$delta
  r <- method_density(m, T27, evidence_terms(r27, pv$plan, pv$year, m$half_life))
  delta_27[[nm]] <- r$delta
  dens27[[nm]] <- r$density
}
for (l in names(score_fns)) {
  d_tr <- delta_tr27[[score_density[[l]]]][sv_vis]
  thr <- calibrate(score_fns[[l]](d_tr[all_rows], sv$top[all_rows], site_feat[all_rows, ]))
  cls27[[l]] <- class_of(score_fns[[l]](delta_27[[score_density[[l]]]][p_vis], p27$top, feat27), thr)
}
tr_d <- ctx$dives[ctx$dives$visit %in% vis$visit[vis$year >= 2019], ]
j <- match(paste(tr_d$visit, tr_d$site_row), sv_key)
ok <- !is.na(j) & is.finite(sv$top[j])
lp_tr27 <- delta_tr27[[base_layer]]
fr <- site_frame(j[ok], sv$top[j[ok]], lp_tr27[sv_vis[j[ok]]])
fr$problematic <- as.numeric(tr_d$outcome[ok] == "problematic")
m27 <- stats::glm(stats::as.formula(paste("problematic ~", site_rhs)), family = stats::binomial(), data = fr)
ok_rows <- all_rows[is.finite(sv$top[all_rows])]
thr <- calibrate(c(stats::predict(m27, site_frame(ok_rows, sv$top[ok_rows], lp_tr27[sv_vis[ok_rows]])),
                   rep(-Inf, length(all_rows) - length(ok_rows))))
cls27[["two-part, site model"]] <- class_of(stats::predict(m27, cbind(feat27, rank_lo = rank_lo(p27$top),
                                                                     reef_lp = delta_27[[base_layer]][p_vis])), thr)
ok_rows <- all_rows[is.finite(sv$top[all_rows])]
fv <- site_frame(ok_rows, sv$top[ok_rows], lp_tr27[sv_vis[ok_rows]])
fv$problematic <- sv_prob[ok_rows]
mv27 <- stats::glm(stats::as.formula(paste("problematic ~", visit_rhs)), family = stats::binomial(), data = fv)
thr <- calibrate(c(stats::predict(mv27, fv), rep(-Inf, length(all_rows) - length(ok_rows))))
new27 <- cbind(feat27, rank_lo = rank_lo(p27$top), reef_lp = delta_27[[base_layer]][p_vis])
cls27[["two-part, site visit model"]] <- class_of(stats::predict(mv27, new27), thr)
comb_all <- rep(-Inf, length(all_rows))
comb_all[match(ok_rows, all_rows)] <- stats::predict(m27, fv) + visit_weight * stats::predict(mv27, fv)
comb27 <- stats::predict(m27, new27) + visit_weight * stats::predict(mv27, new27)
cls27[["two-part, site + visit model"]] <- class_of(comb27, calibrate(comb_all))
thr <- calibrate(delta_tr27[[base_layer]][sv_vis[all_rows]] + neg_q(top_within(comb_all, sv_vis[all_rows])))
cls27[["two-part quota, site + visit order"]] <- class_of(delta_27[[base_layer]][p_vis] + neg_q(top_within(comb27, p27$reef_name)), thr)
ev27 <- evidence_terms(r27, pv$plan, pv$year, density_methods[[base_layer]]$half_life)
ev27 <- (ev27$has_dives > 0 | ev27$has_tows > 0)[p_vis]
cls27[["density-aware v2"]] <- ifelse(ev27, cls27[["two-part, site + visit model"]], cls27[["two-part quota, site + visit order"]])
untowed27 <- feat27$has_tr == 0
new27n <- cbind(impute_tows(feat27, p27$reef_name), rank_lo = rank_lo(p27$top), reef_lp = delta_27[[base_layer]][p_vis])
neut27 <- stats::predict(m27, new27n) + visit_weight * stats::predict(mv27, new27n)
fvn <- site_frame(ok_rows, sv$top[ok_rows], lp_tr27[sv_vis[ok_rows]], site_feat_n)
comb_all_n <- rep(-Inf, length(all_rows))
comb_all_n[match(ok_rows, all_rows)] <- stats::predict(m27, fvn) + visit_weight * stats::predict(mv27, fvn)
cls27[["density-aware v2, untowed neutral"]] <- ifelse(ev27, class_of(neut27, calibrate(comb_all_n)), cls27[["two-part quota, site + visit order"]])
# + habitat: deployed cull-habitat rank from 20 for 2027 sites; forward ranks
# for the 2019-2026 training site-visits.
r20 <- utils::read.csv(file.path(cfg$dirs$outputs, "site_ranking", "cots_cull_site_priority_2027.csv"), stringsAsFactors = FALSE)
hab27 <- r20$habitat_rank_cull[match(paste(p27$reef_name, p27$site_name), paste(r20$reef_name, r20$site_name))] - 0.5
hab27[is.na(hab27)] <- 0
tr_j <- j[ok]
m27h <- stats::glm(stats::as.formula(paste("problematic ~", site_rhs, hab_rhs)), family = stats::binomial(),
                   data = cbind(fr, hab = sv_hab[tr_j], untowed = sv_untowed[tr_j]))
mv27h <- stats::glm(stats::as.formula(paste("problematic ~", visit_rhs, hab_rhs)), family = stats::binomial(),
                    data = with_hab(fv, ok_rows))
fvh <- with_hab(fvn, ok_rows)
comb_all_h <- rep(-Inf, length(all_rows))
comb_all_h[match(ok_rows, all_rows)] <- stats::predict(m27h, fvh) + visit_weight * stats::predict(mv27h, fvh)
new27h <- cbind(new27n, hab = hab27, untowed = as.numeric(untowed27))
hab_score27 <- stats::predict(m27h, new27h) + visit_weight * stats::predict(mv27h, new27h)
cls27[["density-aware v2, untowed neutral + habitat"]] <- ifelse(ev27, class_of(hab_score27, calibrate(comb_all_h)), cls27[["two-part quota, site + visit order"]])
# v3: local tow density for 2027 sites, models fitted on every visit to 2026.
local27 <- lapply(local_k, function(k) local_dev(raw27, p27$reef_name, p27$site_row, k))
local_score27 <- list()
for (v in names(local_layers)) {
  m27l <- stats::glm(stats::as.formula(paste("problematic ~", site_rhs, "+ local_dev")), family = stats::binomial(),
                     data = cbind(fr, local_dev = sv_local[[v]][tr_j]))
  mv27l <- stats::glm(stats::as.formula(paste("problematic ~", visit_rhs, "+ local_dev")), family = stats::binomial(),
                      data = with_local(fv, ok_rows, v))
  fvl <- with_local(fvn, ok_rows, v)
  comb_all_l <- rep(-Inf, length(all_rows))
  comb_all_l[match(ok_rows, all_rows)] <- stats::predict(m27l, fvl) + visit_weight * stats::predict(mv27l, fvl)
  new27l <- cbind(new27n, local_dev = local27[[v]])
  local_score27[[v]] <- stats::predict(m27l, new27l) + visit_weight * stats::predict(mv27l, new27l)
  cls27[[local_layers[[v]]]] <- ifelse(ev27, class_of(local_score27[[v]], calibrate(comb_all_l)), cls27[["two-part quota, site + visit order"]])
}
sites27 <- dplyr::bind_cols(p27[, c("reef_name", "site_name", "rank_on_reef", "sites_on_reef")],
                            searched_clean_12m = feat27$searched_clean, reef_recent_evidence = ev27,
                            site_visit_score = round(comb27, 4), site_visit_rank = stats::ave(-comb27, p27$reef_name,
                                                                                               FUN = function(x) rank(x, ties.method = "first")),
                            untowed_12m = untowed27, site_visit_score_neutral = round(neut27, 4),
                            site_visit_rank_neutral = stats::ave(-neut27, p27$reef_name, FUN = function(x) rank(x, ties.method = "first")),
                            habitat_rank_cull = round(hab27 + 0.5, 3), site_visit_score_habitat = round(hab_score27, 4),
                            site_visit_rank_habitat = stats::ave(-hab_score27, p27$reef_name, FUN = function(x) rank(x, ties.method = "first")),
                            local_tow_dev = round(local27$k0, 3), local_tow_dev_aspect = round(local27$k4, 3),
                            site_visit_score_v3 = round(local_score27$k0, 4), site_visit_score_v3_aspect = round(local_score27$k4, 4),
                            site_visit_rank_v3 = stats::ave(-local_score27$k0, p27$reef_name, FUN = function(x) rank(x, ties.method = "first")),
                            site_visit_rank_v3_aspect = stats::ave(-local_score27$k4, p27$reef_name, FUN = function(x) rank(x, ties.method = "first")),
                            tibble::as_tibble(cls27))
write_table(sites27, cfg, "density_layers_2027_sites.csv")
count_cls <- function(cl, k) as.numeric(tapply(cl %in% k, factor(p27$reef_name, r27), sum))
reef27 <- p27[match(r27, p27$reef_name), c("reef_name", "sites_on_reef", "density_source", "reef_dives_12m",
                                          "reef_problematic_dives_12m", "reef_tows_12m", "reef_strong_tows_12m",
                                          "reef_outbreak_forecast_2026")]
for (nm in names(dens27)) reef27[[paste0(nm, ".density")]] <- dens27[[nm]]
for (l in names(cls27)) {
  reef27[[paste0(l, ".a")]] <- count_cls(cls27[[l]], "A")
  reef27[[paste0(l, ".a_to_c")]] <- count_cls(cls27[[l]], c("A", "B", "C"))
}
write_table(reef27, cfg, "density_methods_2027_reefs.csv")
# Reef-level check: 2027 A and A-C sites by where the reef's density comes from.
totals27 <- dplyr::bind_rows(lapply(names(cls27), function(l) {
  dplyr::summarise(dplyr::group_by(tibble::tibble(source = p27$density_source, reef = p27$reef_name, cl = cls27[[l]]), .data$source),
                   layer = l, reefs = dplyr::n_distinct(.data$reef), sites = dplyr::n(), a = sum(.data$cl == "A"),
                   a_to_c = sum(.data$cl %in% c("A", "B", "C")), share_a_to_c = mean(.data$cl %in% c("A", "B", "C")),
                   reefs_without_a_to_c = sum(tapply(.data$cl %in% c("A", "B", "C"), .data$reef, sum) == 0), .groups = "drop")
}))
write_table(totals27, cfg, "density_layers_2027_totals.csv")
print(as.data.frame(totals27[, c("layer", "source", "reefs", "sites", "a", "a_to_c", "share_a_to_c")]))
show27 <- reef27[grepl("18-100a|18-075|16-029|14-118", reef27$reef_name), ]
print(t(show27[, grepl("density$|a_to_c$", names(show27))]))
# Batt Reef's flanks under each layer (2027), with their tow record since 2019.
b27 <- which(p27$site_row %in% batt)
flank27 <- ifelse(p27$site_row[b27] %in% batt_ne, "NE", "SW")
linked <- as.integer(names(ctx$tow_sites))
since <- tows$survey_date[linked] >= as.Date("2019-01-01")
flank_tows <- function(rows) {
  hit <- since & vapply(ctx$tow_sites, function(j) any(j %in% rows), logical(1))
  c(n = sum(hit), strong = mean(tows$strong[linked[hit]]))
}
ft <- list(NE = flank_tows(batt_ne), SW = flank_tows(setdiff(batt, batt_ne)))
batt_flank <- dplyr::bind_rows(lapply(names(cls27), function(l) dplyr::bind_rows(lapply(c("NE", "SW"), function(f) {
  cl <- cls27[[l]][b27][flank27 == f]
  tibble::tibble(layer = l, flank = f, sites = length(cl), a = sum(cl == "A"), a_to_c = sum(cl %in% c("A", "B", "C")),
                 tows_since_2019 = ft[[f]][["n"]], strong_share = ft[[f]][["strong"]])
}))))
write_table(batt_flank, cfg, "density_batt_flank_2027.csv")
print(as.data.frame(batt_flank))

# A site is in the top `share` of its visit under each layer.
in_class <- function(layer, share, pct, top, delta) {
  if (layer == "within-reef") pct >= 1 - share - 1e-9 else top <= scaled_share(share, delta)
}

# --- Every dive of every test visit: pool, and which tow-zone configurations lose it --
skip_runs <- function(candidate, seg) {
  r <- rle(paste(seg, candidate))
  candidate & rep(r$lengths >= min_skip, r$lengths)
}
sites_near <- function(g) unique(unlist(tow_sites[as.character(g)]))
strong_g <- lapply(guide, function(g) g[tows$strong[g]])
zone_configs <- expand.grid(layer = layers, zone = names(zone_cuts), stringsAsFactors = FALSE)
zone_configs$tag <- paste(zone_configs$layer, zone_configs$zone, sep = "|")
sv_by_visit <- split(seq_len(nrow(sv)), sv$visit)
d_by_visit <- split(seq_len(nrow(ctx$dives)), ctx$dives$visit)
test_idx <- which(vis$year %in% test_years & Reduce(`&`, lapply(vis_delta, is.finite)))
per_visit <- lapply(test_idx, function(i) {
  ev_i <- which(vis$reef_name == vis$reef_name[i] & vis$year == vis$year[i] & vis$first < vis$first[i])
  gs <- unique(unlist(strong_g[ev_i]))
  season_strong <- sites_near(gs[tows$survey_date[gs] < vis$plan[i]])
  g <- guide[[i]]
  s_idx <- sv_by_visit[[vis$visit[i]]]
  s <- sv[s_idx, ]
  d <- ctx$dives[d_by_visit[[vis$visit[i]]], ]
  d <- d[d$site_row %in% s$site_row, ]
  gap <- sqrt((tows$sx[g][-1] - tows$ex[g][-length(g)])^2 + (tows$sy[g][-1] - tows$ey[g][-length(g)])^2)
  seg <- cumsum(c(TRUE, tows$survey_date[g][-1] != tows$survey_date[g][-length(g)] | gap > segment_break_m))
  own <- ctx$own[[i]]
  strong <- tows$strong[g]
  weak <- tows$manta_evidence_class[g] == weak_class
  pool <- dplyr::case_when(d$site_row %in% s$site_row[s$open] ~ "open",
                           d$site_row %in% sites_near(g[strong]) ~ "strong",
                           d$site_row %in% sites_near(g[weak]) ~ "weak",
                           d$site_row %in% sites_near(g) ~ "negative",
                           TRUE ~ "no tow")
  out_d <- data.frame(d[, c("visit", "reef_name", "year", "site_row", "cpue", "outcome", "total_cots", "bottom_time")],
                      pool = pool, follow_up = rep(vis$follow_up[i], nrow(d)))
  out_t <- data.frame(visit = vis$visit[i], reef_name = vis$reef_name[i], year = vis$year[i],
                      follow_up = vis$follow_up[i], tows = sum(own))
  share <- list()
  for (k in seq_len(nrow(zone_configs))) {
    ly <- zone_configs$layer[k]
    zs <- zone_cuts[[zone_configs$zone[k]]]
    zin <- if (ly %in% names(sv_scores)) sv_scores[[ly]][s_idx] >= site_thr[[ly]](zs, vis$year[i]) - 1e-9 else
      in_class(ly, zs, s$zone_combined, s$top, vis_delta[[ly]][i])
    zone_rows <- s$site_row[s$open | zin]
    memory_rows <- union(zone_rows, intersect(season_strong, s$site_row))
    for (mem in c(FALSE, TRUE)) {
      rows <- if (mem) memory_rows else zone_rows
      keep <- vapply(tow_sites[as.character(g)], function(j) any(j %in% rows), logical(1))
      skip <- if (vis$follow_up[i]) skip_runs(!keep & own, seg) else rep(FALSE, length(g))
      tag <- paste(zone_configs$tag[k], if (mem) "memory" else "zone", sep = "|")
      out_d[[tag]] <- pool == "strong" & !(d$site_row %in% sites_near(g[strong & !skip]))
      out_t[[tag]] <- sum(skip)
    }
    share[[zone_configs$tag[k]]] <- mean(s$open | zin)
  }
  list(dives = out_d, tows = out_t, share = data.frame(visit = vis$visit[i], follow_up = vis$follow_up[i],
                                                        delta = vis$delta[i], as.data.frame(share, check.names = FALSE),
                                                        check.names = FALSE))
})
dv <- dplyr::bind_rows(lapply(per_visit, `[[`, "dives"))
tv <- dplyr::bind_rows(lapply(per_visit, `[[`, "tows"))
zone_share <- dplyr::bind_rows(lapply(per_visit, `[[`, "share"))
rank_key <- paste(scored$visit, scored$site_row)
dv$pct <- scored$fwd_voyage[match(paste(dv$visit, dv$site_row), rank_key)]
dv$pct[!is.finite(dv$pct)] <- -1   # unranked sites rank last, as in the product
dv$n_ranked <- sv$n_ranked[match(paste(dv$visit, dv$site_row), paste(sv$visit, sv$site_row))]
dv$top <- top_fraction(dv$pct, dv$n_ranked)
dv_delta <- lapply(vis_delta, function(x) x[match(dv$visit, vis$visit)])
# Site-layer scores of each dive, with the dive's own voyage-start rank.
dv_sv <- match(paste(dv$visit, dv$site_row), sv_key)
dv_feat <- site_feat[dv_sv, ]
dv_feat_n <- site_feat_n[dv_sv, ]
dv_scores <- lapply(stats::setNames(nm = names(score_fns)), function(l) score_fns[[l]](dv_delta[[score_density[[l]]]], dv$top, dv_feat))
dv_scores[["two-part, site model"]] <- dv_scores[["two-part, site visit model"]] <- rep(-Inf, nrow(dv))
dv_comb_n <- dv_comb_h <- rep(-Inf, nrow(dv))
dv_comb_l <- list(k0 = rep(-Inf, nrow(dv)), k4 = rep(-Inf, nrow(dv)))
for (T in test_years) {
  m <- site_models[[as.character(T)]]
  j <- which(dv$year == T & is.finite(dv$top))
  dv_scores[["two-part, site visit model"]][j] <- stats::predict(visit_models[[as.character(T)]],
    cbind(dv_feat[j, ], rank_lo = rank_lo(dv$top[j]), reef_lp = attr(m, "lp")[match(dv$visit[j], vis$visit)]))
  dv_scores[["two-part, site model"]][j] <- stats::predict(m, cbind(dv_feat[j, ], rank_lo = rank_lo(dv$top[j]),
                                                                    reef_lp = attr(m, "lp")[match(dv$visit[j], vis$visit)]))
  fn <- cbind(dv_feat_n[j, ], rank_lo = rank_lo(dv$top[j]), reef_lp = attr(m, "lp")[match(dv$visit[j], vis$visit)])
  dv_comb_n[j] <- stats::predict(m, fn) + visit_weight * stats::predict(visit_models[[as.character(T)]], fn)
  fh <- cbind(fn, hab = sv_hab[dv_sv[j]], untowed = sv_untowed[dv_sv[j]])
  hm <- hab_models[[as.character(T)]]
  dv_comb_h[j] <- stats::predict(hm$dive, fh) + visit_weight * stats::predict(hm$visit, fh)
  for (v in names(local_layers)) {
    lm_ <- local_models[[v]][[as.character(T)]]
    fl <- cbind(fn, local_dev = sv_local[[v]][dv_sv[j]])
    dv_comb_l[[v]][j] <- stats::predict(lm_$dive, fl) + visit_weight * stats::predict(lm_$visit, fl)
  }
}
dv_scores[["two-part, site + visit model"]] <- dv_scores[["two-part, site model"]] + visit_weight * dv_scores[["two-part, site visit model"]]
dv_scores[["two-part quota, site + visit order"]] <- dv_delta[[base_layer]] + neg_q(sv_tq[dv_sv])
dv_ev <- vis_evidence[match(dv$visit, vis$visit)]
dv_scores[["density-aware v2"]] <- ifelse(dv_ev, dv_scores[["two-part, site + visit model"]], map_quota(dv_scores[["two-part quota, site + visit order"]], dv$year))
dv_scores[["density-aware v2, untowed neutral"]] <- ifelse(dv_ev, dv_comb_n,
                            map_quota(dv_scores[["two-part quota, site + visit order"]], dv$year, "density-aware v2, untowed neutral"))
dv_scores[["density-aware v2, untowed neutral + habitat"]] <- ifelse(dv_ev, dv_comb_h, map_quota(dv_scores[["two-part quota, site + visit order"]], dv$year, "density-aware v2, untowed neutral + habitat"))
for (v in names(local_layers)) dv_scores[[local_layers[[v]]]] <- ifelse(dv_ev, dv_comb_l[[v]],
                                                                       map_quota(dv_scores[["two-part quota, site + visit order"]], dv$year, local_layers[[v]]))
# Does habitat help where there are no tows? Within-visit AUC of "problematic"
# (and > 0.04) among dives at sites with no linked tow in the year before the
# visit, for each score.
ut <- which(sv_untowed[dv_sv] == 1 & dv_ev & is.finite(dv_comb_h))
untowed_auc <- dplyr::bind_rows(lapply(c(`habitat percentile alone` = "hab", `voyage-start rank (33)` = "rank",
                                         `v2 untowed neutral` = "n", `v2 untowed neutral + habitat` = "h",
                                         `v3 local tow density` = "l0", `v3 local tow density + aspect` = "l4"), function(k) {
  sc <- switch(k, hab = sv_hab[dv_sv], rank = -dv$top, n = dv_comb_n, h = dv_comb_h, l0 = dv_comb_l$k0, l4 = dv_comb_l$k4)[ut]
  tibble::tibble(dives = length(ut), visits = dplyr::n_distinct(dv$visit[ut]),
                 auc_problematic = within_group_auc(sc, dv$outcome[ut] == "problematic", dv$visit[ut]),
                 auc_high = within_group_auc(sc, dv$cpue[ut] > 0.04, dv$visit[ut]))
}), .id = "score")
write_table(untowed_auc, cfg, "density_untowed_auc.csv")
print(untowed_auc, width = Inf)
dv$high <- dv$cpue > 0.04
dv$inter <- dv$outcome == "problematic" & !dv$high
dv$nonstrong <- dv$pool %in% c("weak", "negative", "no tow")
dv$reef_year <- paste(dv$reef_name, dv$year)
tv$reef_year <- paste(tv$reef_name, tv$year)
log_message(nrow(dv), " dives in ", nrow(tv), " visits (", sum(tv$follow_up), " follow-up), ", min(test_years), "-", max(test_years))

# --- Rules: every combination of layer, tow zone and dive rule -----------------------
dive_rule <- function(layer, cut, evidence) {
  out <- if (layer %in% names(dv_scores)) !(dv_scores[[layer]] >= site_thr[[layer]](rule_cuts[[cut]], dv$year) - 1e-9) else
    !in_class(layer, rule_cuts[[cut]], dv$pct, dv$top, dv_delta[[layer]])
  switch(evidence,
         `discretion at A-x` = dv$nonstrong & out,
         `weak evidence at A-x` = dv$pool %in% c("negative", "no tow") | (dv$pool == "weak" & out),
         `weak-scar dives only` = dv$pool == "weak" & out)
}
rule_pool <- c(`discretion at A-x` = "nonstrong", `weak evidence at A-x` = "nonstrong", `weak-scar dives only` = "weak")
specs <- list(tibble::tibble(layer = "none", zone = "none", memory = FALSE, rule = "none", evidence = "none"))
for (ly in layers) {
  specs[[length(specs) + 1L]] <- tidyr::expand_grid(layer = ly, zone = names(zone_cuts), memory = c(FALSE, TRUE),
                                                    rule = "none", evidence = "none")
  specs[[length(specs) + 1L]] <- tidyr::expand_grid(layer = ly, zone = c("none", names(zone_cuts)), memory = TRUE,
                                                    rule = names(rule_cuts), evidence = evidence_rules)
}
specs <- dplyr::bind_rows(specs)
specs$memory[specs$zone == "none"] <- FALSE
specs <- dplyr::distinct(specs)
specs$id <- with(specs, paste(layer, paste("zone", zone), ifelse(memory, "memory", "no memory"), paste("rule", rule), evidence, sep = " / "))
skips <- lapply(seq_len(nrow(specs)), function(k) {
  sp <- specs[k, ]
  z <- if (sp$zone == "none") rep(FALSE, nrow(dv)) else dv[[paste(sp$layer, sp$zone, if (sp$memory) "memory" else "zone", sep = "|")]]
  r <- if (sp$rule == "none") rep(FALSE, nrow(dv)) else dive_rule(sp$layer, sp$rule, sp$evidence)
  z | r
})
tow_skips <- lapply(seq_len(nrow(specs)), function(k) {
  sp <- specs[k, ]
  if (sp$zone == "none") 0 * tv$tows else tv[[paste(sp$layer, sp$zone, if (sp$memory) "memory" else "zone", sep = "|")]]
})
# Evidence-only rules (no ranking), for reference.
evidence_only <- list(`skip weak-scar dives` = dv$pool == "weak",
                      `skip no-evidence dives (negative-only or no tow)` = dv$pool %in% c("negative", "no tow"),
                      `skip all non-strong dives` = dv$nonstrong)

# --- Metrics, weighted by reef-year resample counts (as in 35) ------------------------
reef_years <- unique(dv$reef_year)
dv_ry <- match(dv$reef_year, reef_years)
tv_ry <- match(tv$reef_year, reef_years)
seasons <- dplyr::n_distinct(dv$year)
measure <- function(skip, dw, tw, sk) {
  hm <- skip & dv$high
  c(dives_avoided = sum(dw * skip), dives_avoided_pct = sum(dw * skip) / sum(dw),
    high_missed = sum(dw * hm), high_missed_pct = sum(dw * hm) / sum(dw * dv$high),
    inter_missed_pct = sum(dw * (skip & dv$inter)) / sum(dw * dv$inter),
    low_avoided = sum(dw * (skip & dv$outcome == "low")),
    low_avoided_pct = sum(dw * (skip & dv$outcome == "low")) / sum(dw * (dv$outcome == "low")),
    ghost_avoided = sum(dw * (skip & dv$outcome == "ghost")),
    ghost_avoided_pct = sum(dw * (skip & dv$outcome == "ghost")) / sum(dw * (dv$outcome == "ghost")),
    ghost_per_high_missed = sum(dw * (skip & dv$outcome == "ghost")) / sum(dw * hm),
    low_ghost_per_high_missed = sum(dw * (skip & dv$outcome %in% c("low", "ghost"))) / sum(dw * hm),
    bottom_time_saved_pct = sum(dw * skip * dv$bottom_time, na.rm = TRUE) / sum(dw * dv$bottom_time, na.rm = TRUE),
    cots_forgone_pct = sum(dw * skip * dv$total_cots, na.rm = TRUE) / sum(dw * dv$total_cots, na.rm = TRUE),
    cpue_missed_high_mean = sum(dw[hm] * dv$cpue[hm]) / sum(dw[hm]),
    tows_saved = sum(tw * sk), tows_saved_pct = sum(tw * sk) / sum(tw * tv$tows),
    tows_saved_follow_up_pct = sum((tw * sk)[tv$follow_up]) / sum((tw * tv$tows)[tv$follow_up]))
}
# Random: the same number of dives skipped at random within the rule's pool.
measure_random <- function(skip, dw, pool) {
  f <- sum(dw * skip) / sum(dw * pool)
  hm <- f * sum(dw * (pool & dv$high))
  c(dives_avoided = sum(dw * skip), dives_avoided_pct = sum(dw * skip) / sum(dw),
    high_missed = hm, high_missed_pct = hm / sum(dw * dv$high),
    low_avoided_pct = f * sum(dw * (pool & dv$outcome == "low")) / sum(dw * (dv$outcome == "low")),
    ghost_avoided_pct = f * sum(dw * (pool & dv$outcome == "ghost")) / sum(dw * (dv$outcome == "ghost")),
    ghost_per_high_missed = sum(dw * (pool & dv$outcome == "ghost")) / sum(dw * (pool & dv$high)))
}
pools <- list(nonstrong = dv$nonstrong, weak = dv$pool == "weak")
all_metrics <- function(w) {
  dw <- w[dv_ry]
  tw <- w[tv_ry]
  out <- list()
  for (scope in c("all visits", "follow-up visits")) {
    sel <- if (scope == "all visits") rep(1, nrow(dv)) else as.numeric(dv$follow_up)
    tsel <- if (scope == "all visits") rep(1, nrow(tv)) else as.numeric(tv$follow_up)
    for (k in seq_len(nrow(specs))) {
      out[[paste(scope, specs$id[k], sep = "#")]] <- measure(skips[[k]], dw * sel, tw * tsel, tow_skips[[k]])
      if (specs$zone[k] == "none" && specs$rule[k] != "none")
        out[[paste(scope, paste("random:", specs$id[k]), sep = "#")]] <-
          measure_random(skips[[k]], dw * sel, pools[[rule_pool[[specs$evidence[k]]]]])
    }
    for (n in names(evidence_only)) out[[paste(scope, paste("evidence only:", n), sep = "#")]] <-
      measure(evidence_only[[n]], dw * sel, tw * tsel, 0 * tv$tows)
  }
  unlist(lapply(names(out), function(n) stats::setNames(out[[n]], paste(n, names(out[[n]]), sep = "#"))))
}
point <- all_metrics(rep(1, length(reef_years)))
boot <- vapply(seq_len(boot_reps), function(r) {
  m <- all_metrics(tabulate(sample.int(length(reef_years), replace = TRUE), length(reef_years)))
  m[names(point)]
}, numeric(length(point)))
rownames(boot) <- names(point)
parts <- do.call(rbind, strsplit(names(point), "#", fixed = TRUE))
spec_of <- specs[match(sub("^random: ", "", parts[, 2]), specs$id), c("layer", "zone", "memory", "rule", "evidence")]
stats_table <- tibble::tibble(scope = parts[, 1], variant = parts[, 2], spec_of, random = startsWith(parts[, 2], "random:"),
                              metric = parts[, 3], value = point,
                              per_season = ifelse(parts[, 3] %in% c("dives_avoided", "high_missed", "low_avoided",
                                                                    "ghost_avoided", "tows_saved"), point / seasons, NA_real_),
                              low = apply(boot, 1, stats::quantile, 0.025, na.rm = TRUE),
                              high = apply(boot, 1, stats::quantile, 0.975, na.rm = TRUE))
write_table(stats_table, cfg, "protocol_density_statistics.csv")

# Paired: each density layer minus the within-reef classes, and each
# alternative density minus 20b's (density-aware), for the same zone and rule.
paired <- function(layer, reference) {
  rows <- which(!stats_table$random & stats_table$layer %in% layer)
  ref_names <- vapply(names(point)[rows], function(n) sub(paste0("#", layer, " / "), paste0("#", reference, " / "), n, fixed = TRUE),
                      character(1))
  stats_table[rows, c("scope", "layer", "zone", "memory", "rule", "evidence", "metric")] |>
    dplyr::mutate(reference = reference, value = point[rows], reference_value = point[ref_names],
                  difference = .data$value - .data$reference_value,
                  low = apply(boot[rows, , drop = FALSE] - boot[ref_names, , drop = FALSE], 1, stats::quantile, 0.025, na.rm = TRUE),
                  high = apply(boot[rows, , drop = FALSE] - boot[ref_names, , drop = FALSE], 1, stats::quantile, 0.975, na.rm = TRUE))
}
differences <- dplyr::bind_rows(lapply(setdiff(layers, "within-reef"), paired, reference = "within-reef"),
                                lapply(score_layers, paired, reference = "density-aware"),
                                lapply(setdiff(score_layers, base_layer), paired, reference = base_layer),
                                lapply("density-aware v2, untowed neutral", paired, reference = "density-aware v2"),
                                lapply(c("density-aware v2, untowed neutral + habitat", "density-aware v3, local tow density", "density-aware v3, local tow density + aspect"), paired, reference = "density-aware v2, untowed neutral"))
write_table(differences, cfg, "protocol_density_differences.csv")

# Matched effort: the layers' cut-offs select different numbers of sites and
# dives, so each layer's trade-off curve (A-D, A-C, A-B, through zero) is
# interpolated at a common point: ghost and low-density dives avoided at 3 %
# of > 0.04 dives missed, > 0.04 missed at 14 % of dives skipped (the
# within-reef A-C rule), and tow-zone misses at 12 % of follow-up tows saved.
# Intervals and paired differences from the same reef-year resamples.
mk <- function(l, z, r, mem, ev, m) paste("all visits", paste(l, paste("zone", z), ifelse(mem, "memory", "no memory"),
                                                             paste("rule", r), ev, sep = " / "), m, sep = "#")
matched_one <- function(l, v) {
  rp <- function(m) vapply(c("A-D", "A-C", "A-B"), function(r) v(mk(l, "none", r, FALSE, "discretion at A-x", m)), 0)
  zp <- function(m) vapply(c("A-C", "A-B"), function(z) v(mk(l, z, "none", TRUE, "none", m)), 0)
  ia <- function(x, y, at) stats::approx(c(0, x), c(0, y), xout = at, ties = mean)$y
  c(ghost_avoided_at_3pct_missed = ia(rp("high_missed_pct"), rp("ghost_avoided_pct"), 0.03),
    low_avoided_at_3pct_missed = ia(rp("high_missed_pct"), rp("low_avoided_pct"), 0.03),
    high_missed_at_14pct_skipped = ia(rp("dives_avoided_pct"), rp("high_missed_pct"), 0.14),
    zone_missed_at_12pct_tows_saved = ia(zp("tows_saved_follow_up_pct"), zp("high_missed"), 0.12))
}
m_point <- lapply(stats::setNames(nm = layers), function(l) matched_one(l, function(n) point[[n]]))
m_boot <- lapply(stats::setNames(nm = layers), function(l)
  vapply(seq_len(ncol(boot)), function(b) matched_one(l, function(n) boot[n, b]), numeric(4)))
q <- function(x, p) stats::quantile(x, p, na.rm = TRUE)
matched <- dplyr::bind_rows(lapply(layers, function(l) tibble::tibble(
  layer = l, metric = names(m_point[[l]]), value = m_point[[l]],
  low = apply(m_boot[[l]], 1, q, 0.025), high = apply(m_boot[[l]], 1, q, 0.975),
  minus_density_aware = m_point[[l]] - m_point[["density-aware"]],
  minus_density_aware_low = apply(m_boot[[l]] - m_boot[["density-aware"]], 1, q, 0.025),
  minus_density_aware_high = apply(m_boot[[l]] - m_boot[["density-aware"]], 1, q, 0.975),
  minus_within_reef = m_point[[l]] - m_point[["within-reef"]],
  minus_within_reef_low = apply(m_boot[[l]] - m_boot[["within-reef"]], 1, q, 0.025),
  minus_within_reef_high = apply(m_boot[[l]] - m_boot[["within-reef"]], 1, q, 0.975))))
write_table(matched, cfg, "density_layers_matched.csv")
print(as.data.frame(matched[, c("layer", "metric", "value", "low", "high", "minus_density_aware", "minus_density_aware_low",
                                "minus_density_aware_high")]))

# Compact comparison of the density layers: tow zone A-B / A-C with memory,
# the A-C dive rule alone, and zone + rule, with intervals.
compare <- stats_table[!stats_table$random & stats_table$scope %in% "all visits" & !is.na(stats_table$layer) &
                         ((stats_table$rule %in% "none" & stats_table$zone %in% names(zone_cuts) & stats_table$memory %in% TRUE) |
                            (stats_table$rule %in% "A-C" & stats_table$evidence %in% "discretion at A-x")) &
                         stats_table$metric %in% c("tows_saved_follow_up_pct", "high_missed", "high_missed_pct", "dives_avoided_pct",
                                                   "ghost_avoided_pct", "low_avoided_pct", "ghost_per_high_missed"),
                       c("layer", "zone", "rule", "metric", "value", "low", "high")]
write_table(compare, cfg, "density_methods_protocol.csv")

# Compact trade-off table: every rule set, all visits.
show <- c("tows_saved_pct", "tows_saved_follow_up_pct", "dives_avoided_pct", "high_missed", "high_missed_pct",
          "inter_missed_pct", "low_avoided_pct", "ghost_avoided_pct", "ghost_per_high_missed", "low_ghost_per_high_missed",
          "bottom_time_saved_pct", "cpue_missed_high_mean")
tradeoff <- stats_table[stats_table$scope == "all visits" & stats_table$metric %in% show,
                        c("variant", "layer", "zone", "memory", "rule", "evidence", "random", "metric", "value")] |>
  tidyr::pivot_wider(names_from = "metric", values_from = "value")
write_table(tradeoff, cfg, "protocol_density_tradeoff.csv")

# How big the zone is: share of a follow-up visit's designated sites in the zone.
zone_size <- tidyr::pivot_longer(zone_share[zone_share$follow_up, ], dplyr::all_of(zone_configs$tag), names_to = "config",
                                 values_to = "share") |>
  tidyr::separate("config", c("layer", "zone"), sep = "[|]") |>
  dplyr::group_by(.data$layer, .data$zone) |>
  dplyr::summarise(follow_up_visits = dplyr::n(), mean_share = mean(.data$share), p10 = stats::quantile(.data$share, 0.1),
                   median = stats::median(.data$share), p90 = stats::quantile(.data$share, 0.9), .groups = "drop")
write_table(zone_size, cfg, "protocol_density_zone_size.csv")
# How each layer treats the sparse parts of dense reefs (test visits on reefs
# denser than typical under the base density): share of designated sites in
# A-C, and what dives there found, for sites searched clean in the last year
# vs the rest; and Batt Reef's north-eastern vs south-western flank.
in_ac <- function(layer, rows) {
  ac <- sr$priority_breaks[["C"]]
  if (layer %in% names(sv_scores)) return(sv_scores[[layer]][rows] >= site_thr[[layer]](ac, vis$year[sv_vis[rows]]) - 1e-9)
  in_class(layer, ac, sv$zone_combined[rows], sv$top[rows], vis_delta[[layer]][sv_vis[rows]])
}
dense_rows <- which(vis$year[sv_vis] %in% test_years & vis$visit[sv_vis] %in% tv$visit &
                      vis_delta[[base_layer]][sv_vis] > 0)
group_of <- function(rows) {
  g <- ifelse(site_feat$searched_clean[rows], "dense reef: searched clean in the last year", "dense reef: other sites")
  b <- sv$site_row[rows] %in% batt
  rbind(data.frame(row = rows, group = g),
        data.frame(row = rows[b], group = ifelse(sv$site_row[rows[b]] %in% batt_ne, "Batt Reef: north-eastern flank", "Batt Reef: south-western flank")))
}
gs <- group_of(which(vis$year[sv_vis] %in% test_years & vis$visit[sv_vis] %in% tv$visit))
gs <- gs[gs$row %in% dense_rows | startsWith(gs$group, "Batt"), ]
dk <- match(sv_key[gs$row], paste(dv$visit, dv$site_row))
suppression <- dplyr::bind_rows(lapply(layers, function(ly) {
  inac <- in_ac(ly, gs$row)
  dplyr::bind_rows(lapply(split(seq_len(nrow(gs)), gs$group), function(j) {
    dd <- dv[dv$visit %in% vis$visit[sv_vis[gs$row[j]]] & paste(dv$visit, dv$site_row) %in% sv_key[gs$row[j]], ]
    tibble::tibble(layer = ly, group = gs$group[j[1]], site_visits = length(j), share_a_to_c = mean(inac[j]),
                   dives = nrow(dd), problematic = mean(dd$outcome == "problematic"), high = mean(dd$high),
                   ghost = mean(dd$outcome == "ghost"))
  }))
}))
write_table(suppression, cfg, "density_site_suppression.csv")
print(as.data.frame(suppression))
# Long-clean ground on dense reefs: sites with >= 5 linked tows and no strong
# tow or dive within 300 m in the 1-6 years before the visit. Does a clean
# record predict the next visit? Share in A-C under each layer, what dives
# there found, and what the A-C dive rule (discretion) missed there.
clean_tow <- site_raw$t_long >= 5 & site_raw$s_long == 0 & site_raw$d_long == 0
ct_rows <- dense_rows[clean_tow[dense_rows]]
ct_dv <- dv_sv %in% ct_rows
clean_tab <- dplyr::bind_rows(lapply(layers, function(ly) {
  skip <- dive_rule(ly, "A-C", "discretion at A-x")
  tibble::tibble(layer = ly, site_visits = length(ct_rows), share_a_to_c = mean(in_ac(ly, ct_rows)), dives = sum(ct_dv),
                 problematic = mean(dv$outcome[ct_dv] == "problematic"), high = sum(dv$high[ct_dv]),
                 high_missed_by_a_to_c_rule = sum(skip & dv$high & ct_dv),
                 ghost_avoided_by_a_to_c_rule = sum(skip & ct_dv & dv$outcome == "ghost"))
}))
write_table(clean_tab, cfg, "density_clean_tow_sites.csv")
print(as.data.frame(clean_tab))

totals <- tibble::tibble(seasons = seasons, visits = nrow(tv), follow_up_visits = sum(tv$follow_up), tows = sum(tv$tows),
                         dives = nrow(dv), high = sum(dv$high), inter = sum(dv$inter), low = sum(dv$outcome == "low"),
                         ghost = sum(dv$outcome == "ghost"), weak_dives = sum(dv$pool == "weak"),
                         no_evidence_dives = sum(dv$pool %in% c("negative", "no tow")), strong_dives = sum(dv$pool == "strong"),
                         open_dives = sum(dv$pool == "open"), median_delta_follow_up = stats::median(vis$delta[vis$visit %in% tv$visit[tv$follow_up]]),
                         median_delta_first = stats::median(vis$delta[vis$visit %in% tv$visit[!tv$follow_up]]))
write_table(totals, cfg, "protocol_density_totals.csv")
print(totals, width = Inf)
print(zone_size, width = Inf)
print(tradeoff[!tradeoff$random & (tradeoff$zone != "none" | tradeoff$rule == "none"),
               c("layer", "zone", "memory", "rule", "evidence", "tows_saved_follow_up_pct", "high_missed", "high_missed_pct",
                 "low_avoided_pct", "ghost_avoided_pct", "ghost_per_high_missed", "cpue_missed_high_mean")], n = Inf, width = Inf)
print(tradeoff[tradeoff$zone %in% "none" & tradeoff$rule != "none" | is.na(tradeoff$layer),
               c("variant", "dives_avoided_pct", "high_missed", "high_missed_pct", "low_avoided_pct", "ghost_avoided_pct",
                 "ghost_per_high_missed")], n = Inf, width = Inf)
