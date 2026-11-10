# Option A: tow only within the priority zone, then let the manta tows guide
# the cull dives, as crews do now.
#
# Question: compared with towing the whole reef, how many tows does a zone
# save, and what does it cost in > 0.04 dives found, and in ghost and
# low-density dives? Is the zone better than skipping the same number of
# tows elsewhere, or than feasible block thinning?
#
# Operational rules, applied identically to every strategy:
#   Open sites   last dive in the 365 days before planning was > 0.04. Dived
#                first under every strategy, so excluded from all outcomes and
#                always inside the towed zone.
#   Guiding tows manta tows on the reef from 30 days before the visit's first
#                dive to its last dive, in recorded order. A new sequence
#                segment starts on a new day or when the next tow starts
#                > 1 km from where the last one ended.
#   Own tows     each tow belongs to the earliest visit whose window holds it.
#                A follow-up visit keeps an earlier visit's tows as evidence
#                (they were towed then) but can skip, and counts as effort,
#                only its own tows.
#   Planning     the zone is set before the visit's own survey: evidence runs
#                to the earlier of voyage departure and its first own tow.
#   Zone         designated sites ranked in the top `cut` of the reef, plus
#                open sites. A tow is in the zone if it passes within 100 m
#                of a zone site. cut = 1 tows only near designated sites.
#   Skipping     only runs of >= 3 consecutive out-of-zone tows (within a
#                segment) are skipped; shorter gaps are towed on the way past.
#   Trigger      crews dive designated sites within 100 m of a kept
#                strong-evidence tow (direct). Neighbour fill: when a kept
#                strong tow borders a skipped run, crews also dive sites within
#                100 m of the adjacent skipped tow.
#   Recent tows  variant: in-zone tows that repeat ground towed in the last
#                365 days (midpoints within 100 m) may also be skipped; a
#                skipped tow whose earlier tow was strong still triggers.
# Rankings: the deployment ranking from 33 (forward habitat + evidence; glm
# coefficients of season T) recomputed at the planning date, and the forward
# habitat rank alone.
# Baselines at the same effort:
#   shifted   the same skip pattern rotated around the visit's tow sequence
#             (20 random offsets): same tows saved, random place
#   blocks    k tows on, m off around the sequence, averaged over phases
#   full      every guiding tow (the reference; crews dived with this)
# Outcomes are actual dives of the visit at non-open sites, so every strategy
# is scored as found or lost relative to full towing. Fill can trigger sites
# with no recorded dive; these are counted as unobserved.
# Intervals: 300 paired bootstrap resamples of reef-years. Seasons 2019-2026.
# Also: absolute dive counts for first vs follow-up visits of the season, and
# where ghost dives come from (open, tow-triggered, other) with the
# deployment ranking's performance in each pool.

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "..", "R", "common.R"))
cfg <- initialise_final()
source(file.path(cfg$final_dir, "R", "forward.R"))
source(file.path(cfg$final_dir, "R", "cull_sites.R"))
source(file.path(cfg$final_dir, "R", "site_selection.R"))
set.seed(cfg$hmsc$seed + 34L)
test_years <- 2019:2026
first_year <- 2013L
boot_reps <- 300L
kernel_sigma_m <- 300
window_days <- 365
lead_days <- 30
site_buffer_m <- 100
min_skip <- 3L
segment_break_m <- 1000
link_distance_m <- 500
zone_cuts <- c(1, 0.9, 0.8, 0.7, 0.66, 0.6, 0.5)   # 0.66 = A-C, 0.5 = A-B
n_shifts <- 20L
blocks <- list(c(9, 3), c(6, 3), c(3, 3), c(3, 6))
strong_classes <- c("direct_cots", "multiple_counted_scars", "legacy_common_scars")
key <- c("reef_name", "year", "site_row")

# --- Sites, dives, visits (as 33) ------------------------------------------------
sites_sf <- read_designated_sites(cfg)
site_xy <- sf::st_coordinates(suppressWarnings(sf::st_point_on_surface(sf::st_geometry(sites_sf))))[, 1:2]
site_reef <- sites_sf$reef_name
dives <- matched_cull_dives(cfg, sites_sf) |>
  dplyr::filter(.data$year >= first_year, !is.na(.data$VoyageTitle)) |>
  dplyr::mutate(visit = paste(.data$reef_name, .data$VoyageTitle, .data$VesselName_ID, sep = "||"))
visit_info <- dives |>
  dplyr::group_by(.data$visit) |>
  dplyr::summarise(reef_name = dplyr::first(.data$reef_name), year = min(.data$year),
                   cut = min(as.Date(.data$VoyageStartDate), .data$survey_date),
                   first = min(.data$survey_date), last = max(.data$survey_date), .groups = "drop")
dives$year <- visit_info$year[match(dives$visit, visit_info$visit)]
has_xy <- is.finite(dives$longitude) & is.finite(dives$latitude)
xy <- matrix(NA_real_, nrow(dives), 2)
xy[has_xy, ] <- sf::st_coordinates(sf::st_transform(sf::st_as_sf(dives[has_xy, ], coords = c("longitude", "latitude"), crs = 4326), sf::st_crs(sites_sf)))
fill <- !has_xy & !is.na(dives$site_row)
xy[fill, ] <- site_xy[dives$site_row[fill], ]
dives$x <- xy[, 1]
dives$y <- xy[, 2]

geom <- readRDS(checked_input(cfg, "footprint_links"))$manta |>
  dplyr::group_by(.data$visit_id) |>
  dplyr::summarise(sx = dplyr::first(.data$start_x), sy = dplyr::first(.data$start_y),
                   ex = dplyr::first(.data$end_x), ey = dplyr::first(.data$end_y), .groups = "drop")
tows <- readRDS(checked_input(cfg, "manta_tows")) |>
  dplyr::select("visit_id", "survey_date", "reef_name", "manta_evidence_class") |>
  dplyr::inner_join(geom, by = "visit_id") |>
  dplyr::filter(is.finite(.data$sx), is.finite(.data$ex)) |>
  dplyr::mutate(strong = .data$manta_evidence_class %in% strong_classes,
                x = (.data$sx + .data$ex) / 2, y = (.data$sy + .data$ey) / 2,
                row = as.integer(sub("^.*[|][|]", "", .data$visit_id)))   # export row = towing order
tows <- tows[order(tows$reef_name, tows$survey_date, tows$row), ]

# --- Visits with guiding tows, and the planning date -------------------------------
vis <- visit_info[visit_info$year %in% test_years & visit_info$reef_name %in% site_reef, ]
tows_by_reef <- split(seq_len(nrow(tows)), tows$reef_name)
guide <- lapply(seq_len(nrow(vis)), function(i) {
  r <- tows_by_reef[[vis$reef_name[i]]]
  r[tows$survey_date[r] >= vis$first[i] - lead_days & tows$survey_date[r] <= vis$last[i]]
})
vis$n_guide <- lengths(guide)
log_message(nrow(vis), " visits ", min(test_years), "-", max(test_years), "; ", sum(vis$n_guide > 0), " with guiding tows")
keep_vis <- vis$n_guide > 0
vis <- vis[keep_vis, ]
guide <- guide[keep_vis]
taken <- logical(nrow(tows))
own <- vector("list", nrow(vis))
for (i in order(vis$first, vis$visit)) {
  own[[i]] <- !taken[guide[[i]]]
  taken[guide[[i]]] <- TRUE
}
log_message(sum(!unlist(own)), " of ", length(unlist(own)), " guiding tows belong to an earlier visit")
vis$plan <- as.Date(vapply(seq_len(nrow(vis)), function(i) {
  r <- guide[[i]][own[[i]]]
  min(as.numeric(vis$cut[i]), if (length(r)) as.numeric(min(tows$survey_date[r])) else Inf)
}, numeric(1)), origin = "1970-01-01")

# --- Evidence, open sites and recent tows at the planning date ---------------------
plan_path <- file.path(cfg$dirs$cache, "manta_zone_plan.rds")
if (!file.exists(plan_path)) {
  site_dives <- dives[!is.na(dives$site_row), ]
  plan <- vector("list", nrow(vis))
  recent <- vector("list", nrow(vis))
  for (i in seq_len(nrow(vis))) {
    rows <- which(site_reef == vis$reef_name[i])
    q <- site_xy[rows, , drop = FALSE]
    box <- c(min(q[, 1]) - 1000, max(q[, 1]) + 1000, min(q[, 2]) - 1000, max(q[, 2]) + 1000)
    start <- vis$plan[i] - window_days
    pd <- dives[dives$survey_date < vis$plan[i] & dives$survey_date >= start & is.finite(dives$x) &
                  dives$x > box[1] & dives$x < box[2] & dives$y > box[3] & dives$y < box[4], ]
    pt <- tows[tows$survey_date < vis$plan[i] & tows$survey_date >= start &
                 tows$x > box[1] & tows$x < box[2] & tows$y > box[3] & tows$y < box[4], ]
    ck <- history_kernel(cbind(pd$x, pd$y), as.numeric(pd$outcome == "problematic"), rep(0, nrow(pd)), q, kernel_sigma_m, Inf)
    mk <- history_kernel(cbind(pt$x, pt$y), as.numeric(pt$strong), rep(0, nrow(pt)), q, kernel_sigma_m, Inf)
    at <- pd[!is.na(pd$site_row), ] |>
      dplyr::arrange(.data$survey_date) |>
      dplyr::group_by(.data$site_row) |>
      dplyr::summarise(prob = any(.data$outcome == "problematic"),
                       open = max(.data$cpue[.data$survey_date == max(.data$survey_date)]) > 0.04, .groups = "drop")
    m <- match(rows, at$site_row)
    plan[[i]] <- data.frame(visit = vis$visit[i], reef_name = vis$reef_name[i], year = vis$year[i], site_row = rows,
                            site_problematic_prev = as.numeric(!is.na(m) & at$prob[m] %in% TRUE),
                            site_quiet_prev = as.numeric(!is.na(m) & at$prob[m] %in% FALSE),
                            cull_rate = ifelse(ck$W > 0, ck$H, 0), cull_effort = ck$W,
                            manta_strong = ifelse(mk$W > 0, mk$H * mk$W, 0), manta_effort = mk$W,
                            open = !is.na(m) & at$open[m] %in% TRUE)
    # Recent ground: an earlier tow (last 365 days) with its midpoint within 100 m.
    g <- guide[[i]]
    covered <- old_strong <- rep(FALSE, length(g))
    if (nrow(pt)) {
      nn <- FNN::get.knnx(cbind(pt$x, pt$y), cbind(tows$x[g], tows$y[g]), k = min(20L, nrow(pt)))
      near <- nn$nn.dist <= site_buffer_m
      covered <- rowSums(near) > 0
      old_strong <- rowSums(near & matrix(pt$strong[nn$nn.index], nrow = length(g))) > 0
    }
    recent[[i]] <- data.frame(visit = vis$visit[i], tow = g, covered = covered, old_strong = old_strong)
  }
  saveRDS(list(plan = dplyr::bind_rows(plan), recent = dplyr::bind_rows(recent)), plan_path)
}
plan_cache <- readRDS(plan_path)

# --- Zone ranks: deployment ranking (33 coefficients) and habitat alone ------------
habitat <- readRDS(file.path(cfg$dirs$cache, "deployment_forward_habitat.rds")) |>
  dplyr::group_by(.data$reef_name, .data$year) |>
  dplyr::mutate(rank_fwd = pct_rank(.data$habitat_fwd)) |>
  dplyr::ungroup()
coefs <- utils::read.csv(file.path(cfg$dirs$tables, "deployment_validation_coefficients.csv")) |>
  dplyr::filter(.data$ranking == "fwd_voyage")
sv <- dplyr::left_join(plan_cache$plan, habitat[, c(key, "rank_fwd")], by = key)
sv$combined <- NA_real_
for (T in test_years) {
  i <- which(sv$year == T)
  b <- coefs[coefs$test_year == T, ]
  beta <- stats::setNames(b$estimate, b$term)
  terms <- site_evidence_terms(sv[i, ], "rank_fwd")
  slopes <- setdiff(names(beta), "(Intercept)")
  sv$combined[i] <- beta[["(Intercept)"]] + as.vector(as.matrix(terms[, slopes]) %*% beta[slopes])
}
sv <- sv |>
  dplyr::group_by(.data$visit) |>
  dplyr::mutate(zone_combined = pct_rank(.data$combined), zone_habitat = pct_rank(.data$rank_fwd)) |>
  dplyr::ungroup()
# Sites without a habitat prediction rank last (as in the deployed product).
for (col in c("zone_combined", "zone_habitat")) sv[[col]][!is.finite(sv[[col]])] <- -1
zone_ranks <- c(combined = "zone_combined", habitat = "zone_habitat")

# --- Tow-site links: designated sites within 100 m of each guiding tow -------------
all_guide <- sort(unique(unlist(guide)))
link_path <- file.path(cfg$dirs$cache, "manta_zone_tow_sites.rds")
if (!file.exists(link_path)) {
  lines <- sf::st_sfc(lapply(all_guide, function(r) sf::st_linestring(rbind(c(tows$sx[r], tows$sy[r]), c(tows$ex[r], tows$ey[r])))),
                      crs = sf::st_crs(sites_sf))
  links <- vector("list", length(all_guide))
  for (reef in unique(tows$reef_name[all_guide])) {
    t_idx <- which(tows$reef_name[all_guide] == reef)
    s_idx <- which(site_reef == reef)
    w <- sf::st_is_within_distance(lines[t_idx], sf::st_geometry(sites_sf)[s_idx], site_buffer_m)
    links[t_idx] <- lapply(w, function(j) s_idx[j])
  }
  saveRDS(stats::setNames(links, all_guide), link_path)
}
tow_sites <- readRDS(link_path)

# --- Per-visit inputs ------------------------------------------------------------
visit_dives <- dives[dives$visit %in% vis$visit & !is.na(dives$site_row), ]
vd_by_visit <- split(seq_len(nrow(visit_dives)), visit_dives$visit)
sv_by_visit <- split(seq_len(nrow(sv)), sv$visit)
rc_by_visit <- split(seq_len(nrow(plan_cache$recent)), plan_cache$recent$visit)
prepare_visit <- function(i) {
  g <- guide[[i]]
  s <- sv[sv_by_visit[[vis$visit[i]]], ]
  local <- function(rows) match(rows, s$site_row)
  adj <- lapply(tow_sites[as.character(g)], local)
  d_all <- visit_dives[vd_by_visit[[vis$visit[i]]], ]
  d_all <- d_all[!is.na(local(d_all$site_row)), ]
  d <- d_all[!s$open[local(d_all$site_row)], ]
  li <- local(d$site_row)
  n_s <- nrow(s)
  count <- function(k) tabulate(li[k %in% TRUE], nbins = n_s)
  gap <- sqrt((tows$sx[g][-1] - tows$ex[g][-length(g)])^2 + (tows$sy[g][-1] - tows$ey[g][-length(g)])^2)
  new_seg <- c(TRUE, tows$survey_date[g][-1] != tows$survey_date[g][-length(g)] | gap > segment_break_m)
  agg <- rep(NA_integer_, length(g))
  st <- which(tows$strong[g])
  if (length(st) == 1) agg[st] <- 1L
  if (length(st) > 1) {
    xy_s <- cbind(tows$x[g][st], tows$y[g][st])
    agg[st] <- stats::cutree(stats::hclust(stats::dist(xy_s), method = "single"), h = link_distance_m)
  }
  rc <- plan_cache$recent[rc_by_visit[[vis$visit[i]]], ]
  list(adj = adj, seg = cumsum(new_seg), strong = tows$strong[g], agg = agg, own = own[[i]],
       covered = rc$covered[match(g, rc$tow)], old_strong = rc$old_strong[match(g, rc$tow)],
       open = s$open, ranks = s[, unname(zone_ranks)], n_sites = n_s,
       high = count(d$cpue > 0.04), inter = count(d$outcome == "problematic" & d$cpue <= 0.04),
       low = count(d$outcome == "low"), ghost = count(d$outcome == "ghost"),
       dived = tabulate(li, nbins = n_s) > 0,
       n_high_open = sum(d_all$cpue > 0.04, na.rm = TRUE) - sum(d$cpue > 0.04, na.rm = TRUE))
}

# Skip only runs of >= min_skip consecutive candidate tows within a segment.
skip_runs <- function(candidate, seg) {
  r <- rle(paste(seg, candidate))
  long <- rep(r$lengths >= min_skip, r$lengths)
  candidate & long
}

count_names <- c("tows", "skipped", "strong", "strong_kept", "aggs", "aggs_found",
                 paste0(rep(c("direct_", "fill_"), each = 6), c("high", "inter", "low", "ghost", "sites", "unobserved")))
evaluate_skip <- function(v, skip, use_old) {
  kept <- !skip
  ks <- kept & v$strong
  trig <- if (use_old) ks | (skip & v$old_strong) else ks
  n <- length(skip)
  same <- v$seg[-1] == v$seg[-n]
  nb <- c(FALSE, ks[-n] & skip[-1] & same) | c(skip[-n] & ks[-1] & same, FALSE)
  direct <- logical(v$n_sites)
  direct[unlist(v$adj[trig])] <- TRUE
  filled <- direct
  filled[unlist(v$adj[nb])] <- TRUE
  outcome <- function(x) {
    x <- x & !v$open
    c(sum(v$high[x]), sum(v$inter[x]), sum(v$low[x]), sum(v$ghost[x]), sum(x), sum(x & !v$dived))
  }
  aggs <- unique(v$agg[!is.na(v$agg)])
  c(sum(v$own), sum(skip), sum(v$strong), sum(ks), length(aggs), length(unique(v$agg[kept & !is.na(v$agg)])),
    outcome(direct), outcome(filled))
}

shift <- function(x, k) x[((seq_along(x) - 1L + k) %% length(x)) + 1L]
run_visit <- function(i) {
  v <- prepare_visit(i)
  n <- length(v$strong)
  o <- which(v$own)
  n_own <- length(o)
  near_site <- lengths(v$adj) > 0 & v$own
  # Baselines move skips among the visit's own tows only.
  shift_own <- function(skip, k) replace(skip, o, shift(skip[o], k))
  offsets <- if (n_own > 1) sample.int(n_own - 1L, min(n_shifts, n_own - 1L)) else 0L
  out <- list(full = evaluate_skip(v, rep(FALSE, n), FALSE))
  for (rk in names(zone_ranks)) {
    rank <- v$ranks[[zone_ranks[[rk]]]]
    for (cut in zone_cuts) {
      in_zone_site <- v$open | rank >= 1 - cut - 1e-9
      in_zone <- vapply(v$adj, function(j) any(in_zone_site[j]), logical(1))
      for (use_old in c(FALSE, TRUE)) {
        cand <- (!in_zone | (use_old & v$covered)) & v$own
        skip <- skip_runs(cand, v$seg)
        tag <- paste(rk, cut, if (use_old) "recent" else "fresh", sep = "|")
        out[[paste("zone", tag, sep = "|")]] <- evaluate_skip(v, skip, use_old)
        out[[paste("shifted", tag, sep = "|")]] <- rowMeans(vapply(offsets, function(k) evaluate_skip(v, shift_own(skip, k), use_old),
                                                                   numeric(length(count_names))))
      }
    }
  }
  for (b in blocks) {
    period <- sum(b)
    pattern <- rep(c(rep(FALSE, b[1]), rep(TRUE, b[2])), length.out = n_own + period)
    out[[paste("block", paste0(b[1], "on", b[2], "off"), NA, "fresh", sep = "|")]] <-
      rowMeans(vapply(seq_len(period) - 1L, function(p) evaluate_skip(v, replace(logical(n), o, pattern[p + seq_len(n_own)]), FALSE), numeric(length(count_names))))
  }
  m <- do.call(rbind, out)
  colnames(m) <- count_names
  parts <- do.call(rbind, strsplit(c("full|all|1|fresh", names(out)[-1]), "|", fixed = TRUE))
  data.frame(visit = vis$visit[i], reef_year = paste(vis$reef_name[i], vis$year[i]), year = vis$year[i],
             strategy = parts[, 1], ranking = parts[, 2], cut = suppressWarnings(as.numeric(parts[, 3])), tows_from = parts[, 4],
             near_site = sum(near_site), high_open = v$n_high_open,
             high_all = sum(v$high), m, row.names = NULL)
}
counts_path <- file.path(cfg$dirs$cache, "manta_zone_counts.rds")
if (!file.exists(counts_path)) {
  res <- vector("list", nrow(vis))
  for (i in seq_len(nrow(vis))) {
    res[[i]] <- run_visit(i)
    if (i %% 250 == 0) log_message("Scored ", i, " of ", nrow(vis), " visits")
  }
  saveRDS(dplyr::bind_rows(res), counts_path)
}
counts <- readRDS(counts_path)
counts$config <- paste(counts$strategy, counts$ranking, counts$cut, counts$tows_from, sep = "|")

# --- Coverage of the trigger rule --------------------------------------------------
full <- counts[counts$strategy == "full", ]
coverage <- tibble::tibble(
  visits = nrow(full), reef_years = dplyr::n_distinct(full$reef_year), guiding_tows = sum(full$tows),
  tows_near_site = sum(full$near_site) / sum(full$tows),
  high_dives = sum(full$high_all) + sum(full$high_open), high_at_open_sites = sum(full$high_open) / (sum(full$high_all) + sum(full$high_open)),
  high_nonopen_triggered_full = sum(full$direct_high) / sum(full$high_all),
  ghost_nonopen_triggered_full = sum(full$direct_ghost),
  strong_tows = sum(full$strong), aggregations = sum(full$aggs))
write_table(coverage, cfg, "manta_zone_coverage.csv")
print(coverage, width = Inf)

# --- Metrics, with paired reef-year bootstrap -----------------------------------------
configs <- unique(counts$config)
reef_years <- unique(counts$reef_year)
vars <- count_names
tot <- lapply(stats::setNames(vars, vars), function(cn) {
  m <- tapply(counts[[cn]], list(factor(counts$reef_year, reef_years), factor(counts$config, configs)), sum)
  m[is.na(m)] <- 0
  m
})
metrics_from <- function(w) {
  s <- lapply(tot, function(m) colSums(m * w))
  f <- which(configs == "full|all|1|fresh")
  unlist(lapply(c("direct", "fill"), function(tr) {
    g <- function(x) s[[paste0(tr, "_", x)]]
    v <- rbind(tows_saved = s$skipped / s$tows,
               strong_kept = s$strong_kept / s$strong,
               aggs_found = s$aggs_found / s$aggs,
               high_found = g("high") / g("high")[f],
               inter_found = g("inter") / g("inter")[f],
               low_incurred = g("low") / g("low")[f],
               ghost_incurred = g("ghost") / g("ghost")[f],
               sites_triggered = g("sites") / g("sites")[f],
               unobserved_sites = g("unobserved") / g("sites")[f])
    stats::setNames(as.vector(v), paste(tr, rep(configs, each = nrow(v)), rownames(v), sep = "#"))
  }))
}
point <- metrics_from(rep(1, length(reef_years)))
boot <- vapply(seq_len(boot_reps), function(r) metrics_from(tabulate(sample.int(length(reef_years), replace = TRUE), length(reef_years))),
               numeric(length(point)))
parts <- do.call(rbind, strsplit(names(point), "#", fixed = TRUE))
cp <- do.call(rbind, strsplit(parts[, 2], "|", fixed = TRUE))
results <- tibble::tibble(trigger = parts[, 1], strategy = cp[, 1], ranking = cp[, 2], cut = suppressWarnings(as.numeric(cp[, 3])),
                          tows_from = cp[, 4], metric = parts[, 3], value = point,
                          low = apply(boot, 1, stats::quantile, 0.025, na.rm = TRUE),
                          high = apply(boot, 1, stats::quantile, 0.975, na.rm = TRUE))
write_table(results, cfg, "manta_zone_strategy.csv")

# Zone minus the same skip pattern placed at random (paired).
zone_rows <- which(results$strategy == "zone")
null_rows <- match(paste(results$trigger, "shifted", results$ranking, results$cut, results$tows_from, results$metric)[zone_rows],
                   paste(results$trigger, results$strategy, results$ranking, results$cut, results$tows_from, results$metric))
diffs <- results[zone_rows, c("trigger", "ranking", "cut", "tows_from", "metric")] |>
  dplyr::mutate(difference = point[zone_rows] - point[null_rows],
                low = apply(boot[zone_rows, ] - boot[null_rows, ], 1, stats::quantile, 0.025, na.rm = TRUE),
                high = apply(boot[zone_rows, ] - boot[null_rows, ], 1, stats::quantile, 0.975, na.rm = TRUE))
write_table(diffs, cfg, "manta_zone_strategy_differences.csv")

# --- Absolute counts, by first visit of the season vs follow-up visit ---------------
# A follow-up visit comes after an earlier visit to the same reef that season
# (the "between sweeps" case on an active-culling reef).
visit_info <- visit_info |>
  dplyr::group_by(.data$reef_name, .data$year) |>
  dplyr::mutate(follow_up = rank(.data$first, ties.method = "min") > 1) |>
  dplyr::ungroup()
counts$visit_type <- ifelse(visit_info$follow_up[match(counts$visit, visit_info$visit)], "follow-up", "first of season")
absolute <- dplyr::bind_rows(lapply(list(all = unique(counts$visit_type), `follow-up` = "follow-up", `first of season` = "first of season"), function(vt) {
  counts[counts$visit_type %in% vt, ] |>
    dplyr::group_by(.data$config, .data$strategy, .data$ranking, .data$cut, .data$tows_from) |>
    dplyr::summarise(visits = dplyr::n(), tows = sum(.data$tows), tows_saved = sum(.data$skipped),
                     high = sum(.data$direct_high), intermediate = sum(.data$direct_inter),
                     low = sum(.data$direct_low), ghost = sum(.data$direct_ghost),
                     aggs_found = sum(.data$aggs_found) / sum(.data$aggs), .groups = "drop")
}), .id = "visit_type")
ref <- absolute[absolute$config == "full|all|1|fresh", ]
r <- match(absolute$visit_type, ref$visit_type)
absolute <- absolute |>
  dplyr::mutate(high_lost = ref$high[r] - .data$high, intermediate_lost = ref$intermediate[r] - .data$intermediate,
                low_avoided = ref$low[r] - .data$low, ghost_avoided = ref$ghost[r] - .data$ghost)
write_table(dplyr::select(absolute, -"config"), cfg, "manta_zone_absolute.csv")

# --- Where ghost dives come from: open, tow-triggered, other ------------------------
# Every dive of the visit, by why it could have been chosen: an open site, a
# site within 100 m of a strong guiding tow, or neither (crew judgement, weak
# or older evidence). Scored with 33's deployment ranking (evidence to voyage
# start), which is what a crew would hold when choosing dives.
scored <- readRDS(file.path(cfg$dirs$cache, "deployment_validation_scored.rds"))
scored <- as.data.frame(scored[scored$visit %in% vis$visit, ])
# The "other" pool is split by the strongest guiding tow within 100 m: a weak
# scar (single or "present"), only negative tows, or no guiding tow at all.
sites_near <- function(g) unique(unlist(tow_sites[as.character(g)]))
weak_class <- "weak_single_or_present_scar"
near_strong <- lapply(guide, function(g) sites_near(g[tows$strong[g]]))
near_weak <- lapply(guide, function(g) sites_near(g[tows$manta_evidence_class[g] == weak_class]))
near_any <- lapply(guide, sites_near)
names(near_strong) <- names(near_weak) <- names(near_any) <- vis$visit
has <- function(lst) mapply(function(v, s) s %in% lst[[v]], scored$visit, scored$site_row)
open_key <- paste(sv$visit, sv$site_row)[sv$open]
scored$pool <- dplyr::case_when(paste(scored$visit, scored$site_row) %in% open_key ~ "open",
                                has(near_strong) ~ "tow-triggered",
                                has(near_weak) ~ "other: weak-scar tow",
                                has(near_any) ~ "other: negative tows only",
                                TRUE ~ "other: no guiding tow")
pool_names <- c("open", "tow-triggered", "other: weak-scar tow", "other: negative tows only", "other: no guiding tow")
pool_metrics <- function(d) {
  d$pool_all_other <- startsWith(d$pool, "other")
  unlist(lapply(c(pool_names, "other: all"), function(p) {
    x <- if (p == "other: all") d[d$pool_all_other, ] else d[d$pool == p, ]
    hg <- x$high | x$outcome == "ghost"
    v <- c(dives = nrow(x), high = sum(x$high), ghost = sum(x$outcome == "ghost"), low = sum(x$outcome == "low"),
           auc_high_vs_ghost = within_group_auc(x$fwd_voyage[hg], x$high[hg], x$visit[hg]),
           auc_high_vs_rest = within_group_auc(x$fwd_voyage, x$high, x$visit))
    for (a in c(0.8, 0.5)) {
      k <- x$fwd_voyage >= 1 - a - 1e-9
      v[paste0(c("high_kept_top", "ghost_avoided_top", "low_avoided_top"), 100 * a)] <-
        c(mean(k[x$high]), mean(!k[x$outcome == "ghost"]), mean(!k[x$outcome == "low"]))
    }
    stats::setNames(v, paste(p, names(v), sep = "#"))
  }))
}
pp <- pool_metrics(scored)
pg <- split(seq_len(nrow(scored)), paste(scored$reef_name, scored$year))
pb <- vapply(seq_len(boot_reps), function(r) pool_metrics(scored[unlist(pg[sample.int(length(pg), replace = TRUE)], use.names = FALSE), ]),
             numeric(length(pp)))
pparts <- do.call(rbind, strsplit(names(pp), "#", fixed = TRUE))
pools <- tibble::tibble(pool = pparts[, 1], metric = pparts[, 2], value = pp,
                        low = apply(pb, 1, stats::quantile, 0.025, na.rm = TRUE),
                        high = apply(pb, 1, stats::quantile, 0.975, na.rm = TRUE))
pools[pools$metric %in% c("dives", "high", "ghost", "low"), c("low", "high")] <- NA
write_table(pools, cfg, "manta_zone_dive_pools.csv")

# Context for the dive-level protocol statistics (35).
saveRDS(list(vis = vis, guide = guide, own = own, tows = tows[, c("survey_date", "sx", "sy", "ex", "ey", "strong", "manta_evidence_class")],
             tow_sites = tow_sites, sv = sv[, c("visit", "site_row", "open", "zone_combined")],
             dives = visit_dives, visit_info = visit_info),
        file.path(cfg$dirs$cache, "manta_zone_context.rds"))
print(tidyr::pivot_wider(pools[, c("pool", "metric", "value")], names_from = "pool", values_from = "value"), n = Inf)
print(absolute[absolute$config %in% c("full|all|1|fresh", "zone|combined|0.5|fresh", "shifted|combined|0.5|fresh", "block|9on3off|NA|fresh"),
               c("visit_type", "strategy", "ranking", "cut", "visits", "tows_saved", "high_lost", "low_avoided", "ghost_avoided", "aggs_found")], n = Inf)

show <- c("tows_saved", "high_found", "inter_found", "low_incurred", "ghost_incurred", "aggs_found", "unobserved_sites")
wide <- tidyr::pivot_wider(results[results$metric %in% show, c("trigger", "strategy", "ranking", "cut", "tows_from", "metric", "value")],
                           names_from = "metric", values_from = "value")
print(wide[wide$tows_from == "fresh" | wide$strategy == "full", ], n = Inf, width = Inf)
print(wide[wide$tows_from == "recent", ], n = Inf, width = Inf)
print(diffs[diffs$metric %in% c("high_found", "ghost_incurred", "low_incurred") & diffs$cut %in% c(1, 0.8, 0.6), ], n = Inf)
