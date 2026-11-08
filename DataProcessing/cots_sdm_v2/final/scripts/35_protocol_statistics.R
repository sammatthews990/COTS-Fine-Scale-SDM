# Statistics for the proposed manta tow and cull dive protocol.
#
# Question: if crews followed the proposed protocol, how many tows, ghost and
# low-density dives would it save, how many > 0.04 dives would it miss, and
# what CPUE did those missed dives have? How much of the gain comes from the
# strong vs weak manta evidence rule alone, without the site ranking? What
# does a single voyage look like, for planning a one-voyage field trial?
#
# Protocol, scored on the recorded dives of 2019-2026 visits (34's rules):
#   1 Open sites are dived first.
#   2 Follow-up visits (after an earlier visit to the reef that season) tow
#     only the A-C zone (top 66 %, open sites, and with season memory sites
#     near a strong tow on an earlier visit that season; 34's skip rules);
#     first visits of the season tow everything.
#   3 Sites within 100 m of a kept strong-evidence tow are dived.
#   4 Any other dive (weak scar, only negative tows, or no guiding tow nearby;
#     crew discretion, e.g. adjacent sites or long-unvisited ground) is made
#     only if the site ranks A-C (top 66 %; or A-D, top 80 %; or A-B, top
#     50 %; classes A 20 %, B 50 %, C 66 %, D 80 %, E the rest) in the
#     deployment ranking at voyage start (33, fwd_voyage).
# A dive whose only strong tows were skipped is counted as missed: crews are
# assumed not to dive there without the tow. Tows an earlier visit already made
# are never skipped again, and effort counts each tow once (34).
# Variants: the site-only rule (skip every dive at a site outside the top 80 %,
# as in 23's A-D cutoff, but on these visits and voyage-start ranks); the
# strong-evidence rule without the ranking (skip weak-scar
# dives; skip all non-strong dives), skipping no-evidence dives, the rank
# filter on weak-scar dives only, the dive rule (rank filter on all non-strong
# dives), the zone alone, and the full protocol. Baselines: the same number
# of dives skipped at random within the same pool (exact expectation).
# Intervals: 300 bootstrap resamples of reef-years.
# Trial planning: per-voyage counts (2024-2026 seasons), and how often one or
# several voyages would show the filter working (flagged non-strong dives with
# a lower > 0.04 rate and a higher ghost rate than kept non-strong dives), by
# resampling historical voyages.

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "..", "R", "common.R"))
cfg <- initialise_final()
set.seed(cfg$hmsc$seed + 35L)
boot_reps <- 300L
trial_reps <- 2000L
min_skip <- 3L
segment_break_m <- 1000
zone_cut <- 0.66   # tow zone A-C (was A-B, 0.5, until 6 November 2026)
filter_keep <- c(`A-D` = 0.8, `A-C` = 0.66, `A-B` = 0.5)
trial_filter <- "A-C"
trial_years <- 2024:2026
pooled_voyages <- c(1L, 2L, 3L, 5L, 10L, 25L, 50L)
weak_class <- "weak_single_or_present_scar"

# --- Inputs: 34's visits, guiding tows and zone ranks; 33's voyage ranking --------
ctx <- readRDS(file.path(cfg$dirs$cache, "manta_zone_context.rds"))
tows <- ctx$tows
vis <- ctx$vis
guide <- ctx$guide
tow_sites <- ctx$tow_sites
vis$follow_up <- ctx$visit_info$follow_up[match(vis$visit, ctx$visit_info$visit)]
scored <- as.data.frame(readRDS(file.path(cfg$dirs$cache, "deployment_validation_scored.rds")))
rank_key <- paste(scored$visit, scored$site_row)

# --- Every dive of every visit: why it was chosen, and whether the zone loses it ---
skip_runs <- function(candidate, seg) {
  r <- rle(paste(seg, candidate))
  candidate & rep(r$lengths >= min_skip, r$lengths)
}
sites_near <- function(g) unique(unlist(tow_sites[as.character(g)]))
# Season memory: sites within 100 m of a strong tow made on an earlier visit
# that season, before this visit's planning date. A variant keeps them in the
# zone whatever their rank.
strong_g <- lapply(guide, function(g) g[tows$strong[g]])
season_strong <- lapply(seq_len(nrow(vis)), function(i) {
  ev <- which(vis$reef_name == vis$reef_name[i] & vis$year == vis$year[i] & vis$first < vis$first[i])
  g <- unique(unlist(strong_g[ev]))
  sites_near(g[tows$survey_date[g] < vis$plan[i]])
})
sv_by_visit <- split(seq_len(nrow(ctx$sv)), ctx$sv$visit)
d_by_visit <- split(seq_len(nrow(ctx$dives)), ctx$dives$visit)
per_visit <- lapply(seq_len(nrow(vis)), function(i) {
  g <- guide[[i]]
  s <- ctx$sv[sv_by_visit[[vis$visit[i]]], ]
  d <- ctx$dives[d_by_visit[[vis$visit[i]]], ]
  d <- d[d$site_row %in% s$site_row, ]
  zone_rows <- s$site_row[s$open | s$zone_combined >= 1 - zone_cut - 1e-9]
  in_zone <- vapply(tow_sites[as.character(g)], function(j) any(j %in% zone_rows), logical(1))
  memory_rows <- union(zone_rows, intersect(season_strong[[i]], s$site_row))
  in_memory <- vapply(tow_sites[as.character(g)], function(j) any(j %in% memory_rows), logical(1))
  gap <- sqrt((tows$sx[g][-1] - tows$ex[g][-length(g)])^2 + (tows$sy[g][-1] - tows$ey[g][-length(g)])^2)
  seg <- cumsum(c(TRUE, tows$survey_date[g][-1] != tows$survey_date[g][-length(g)] | gap > segment_break_m))
  own <- ctx$own[[i]]   # tows an earlier visit already made are kept (34)
  skip <- if (vis$follow_up[i]) skip_runs(!in_zone & own, seg) else rep(FALSE, length(g))
  skip_mem <- if (vis$follow_up[i]) skip_runs(!in_memory & own, seg) else rep(FALSE, length(g))
  strong <- tows$strong[g]
  weak <- tows$manta_evidence_class[g] == weak_class
  pool <- dplyr::case_when(d$site_row %in% s$site_row[s$open] ~ "open",
                           d$site_row %in% sites_near(g[strong]) ~ "strong",
                           d$site_row %in% sites_near(g[weak]) ~ "weak",
                           d$site_row %in% sites_near(g) ~ "negative",
                           TRUE ~ "no tow")
  list(dives = data.frame(d[, c("visit", "reef_name", "year", "site_row", "site_name", "survey_date", "VoyageTitle", "VesselName_ID",
                                "cpue", "outcome", "total_cots", "bottom_time")],
                          pool = pool, follow_up = rep(vis$follow_up[i], nrow(d)),
                          zone_lost = pool == "strong" & !(d$site_row %in% sites_near(g[strong & !skip])),
                          memory_lost = pool == "strong" & !(d$site_row %in% sites_near(g[strong & !skip_mem]))),
       tows = data.frame(visit = vis$visit[i], reef_name = vis$reef_name[i], year = vis$year[i],
                         follow_up = vis$follow_up[i], tows = sum(own), skipped = sum(skip), skipped_mem = sum(skip_mem)))
})
dv <- dplyr::bind_rows(lapply(per_visit, `[[`, "dives"))
tv <- dplyr::bind_rows(lapply(per_visit, `[[`, "tows"))
dv$rank <- scored$fwd_voyage[match(paste(dv$visit, dv$site_row), rank_key)]
dv$rank[!is.finite(dv$rank)] <- -1   # unranked sites rank last, as in the product
dv$high <- dv$cpue > 0.04
dv$inter <- dv$outcome == "problematic" & !dv$high
dv$nonstrong <- dv$pool %in% c("weak", "negative", "no tow")
dv$reef_year <- paste(dv$reef_name, dv$year)
tv$reef_year <- paste(tv$reef_name, tv$year)
dv$voyage <- paste(dv$VoyageTitle, dv$VesselName_ID, sep = "||")
tv$voyage <- sub("^[^|]*[|][|]", "", tv$visit)
log_message(nrow(dv), " dives in ", nrow(tv), " visits; ", sum(dv$rank < 0), " dives at unranked sites")

weak <- dv$pool == "weak"
no_evidence <- dv$pool %in% c("negative", "no tow")
flagged <- lapply(filter_keep, function(a) weak & dv$rank < 1 - a - 1e-9)
flagged_ns <- lapply(filter_keep, function(a) dv$nonstrong & dv$rank < 1 - a - 1e-9)
# The lowest-ranked non-strong dives, as many as the weak-scar rule skips.
matched <- rep(FALSE, nrow(dv))
ns_idx <- which(dv$nonstrong)
matched[ns_idx[order(dv$rank[ns_idx])][seq_len(sum(dv$pool == "weak"))]] <- TRUE
variants <- list(
  `current practice` = rep(FALSE, nrow(dv)),
  `site rule A-D, all dives` = dv$rank < 1 - filter_keep[["A-D"]] - 1e-9,
  `evidence rule: skip weak-scar dives` = dv$pool == "weak",
  `evidence rule: skip all non-strong dives` = dv$nonstrong,
  `ranking filter, as many dives as the weak-scar rule` = matched,
  `skip no-evidence dives` = no_evidence,
  `skip negative-only dives` = dv$pool == "negative",
  `skip no-tow dives` = dv$pool == "no tow",
  `weak-scar filter A-D` = flagged$`A-D`,
  `weak-scar filter A-C` = flagged$`A-C`,
  `weak-scar filter A-B` = flagged$`A-B`,
  `dive rule A-D` = flagged_ns$`A-D`,
  `dive rule A-C` = flagged_ns$`A-C`,
  `dive rule A-B` = flagged_ns$`A-B`,
  `zone only` = dv$zone_lost,
  `zone with season memory` = dv$memory_lost,
  `protocol: zone with season memory + dive rule A-D` = dv$memory_lost | flagged_ns$`A-D`,
  `protocol: zone with season memory + dive rule A-C` = dv$memory_lost | flagged_ns$`A-C`,
  `protocol: zone with season memory + dive rule A-B` = dv$memory_lost | flagged_ns$`A-B`)
uses_zone <- ifelse(grepl("season memory", names(variants)), "memory", ifelse(grepl("zone", names(variants)), "zone", "none"))
random_of <- c(`weak-scar filter A-D` = "random weak-scar filter A-D", `weak-scar filter A-C` = "random weak-scar filter A-C",
               `weak-scar filter A-B` = "random weak-scar filter A-B", `dive rule A-D` = "random dive rule A-D",
               `dive rule A-C` = "random dive rule A-C", `dive rule A-B` = "random dive rule A-B")

# --- Metrics, weighted by reef-year resample counts ---------------------------------
wquant <- function(x, w, p) {
  o <- order(x)
  x <- x[o]
  cw <- cumsum(w[o]) / sum(w)
  vapply(p, function(q) x[which(cw >= q)[1]], numeric(1))
}
reef_years <- unique(dv$reef_year)
dv_ry <- match(dv$reef_year, reef_years)
tv_ry <- match(tv$reef_year, reef_years)
seasons <- dplyr::n_distinct(dv$year)
measure <- function(skip, w, dw, tw, zone) {
  sk <- if (zone == "memory") tv$skipped_mem else if (zone == "zone") tv$skipped else 0 * tv$skipped
  hm <- skip & dv$high
  kh <- !skip & dv$high
  q <- if (sum(dw[hm]) > 0) wquant(dv$cpue[hm], dw[hm], c(0.25, 0.5, 0.75)) else rep(NA_real_, 3)
  c(dives_avoided = sum(dw * skip), dives_avoided_pct = sum(dw * skip) / sum(dw),
    high_missed = sum(dw * hm), high_missed_pct = sum(dw * hm) / sum(dw * dv$high),
    inter_missed = sum(dw * (skip & dv$inter)),
    inter_missed_pct = sum(dw * (skip & dv$inter)) / sum(dw * dv$inter),
    low_avoided = sum(dw * (skip & dv$outcome == "low")),
    low_avoided_pct = sum(dw * (skip & dv$outcome == "low")) / sum(dw * (dv$outcome == "low")),
    ghost_avoided = sum(dw * (skip & dv$outcome == "ghost")),
    ghost_avoided_pct = sum(dw * (skip & dv$outcome == "ghost")) / sum(dw * (dv$outcome == "ghost")),
    ghost_per_high_missed = sum(dw * (skip & dv$outcome == "ghost")) / sum(dw * hm),
    share_deferred_ghost = sum(dw * (skip & dv$outcome == "ghost")) / sum(dw * skip),
    share_deferred_low = sum(dw * (skip & dv$outcome == "low")) / sum(dw * skip),
    share_deferred_intermediate = sum(dw * (skip & dv$inter)) / sum(dw * skip),
    share_deferred_high = sum(dw * hm) / sum(dw * skip),
    bottom_time_saved_pct = sum(dw * skip * dv$bottom_time, na.rm = TRUE) / sum(dw * dv$bottom_time, na.rm = TRUE),
    cots_forgone_pct = sum(dw * skip * dv$total_cots, na.rm = TRUE) / sum(dw * dv$total_cots, na.rm = TRUE),
    cpue_missed_high_mean = sum(dw[hm] * dv$cpue[hm]) / sum(dw[hm]),
    cpue_missed_high_q25 = q[1], cpue_missed_high_median = q[2], cpue_missed_high_q75 = q[3],
    cpue_kept_high_mean = sum(dw[kh] * dv$cpue[kh]) / sum(dw[kh]),
    tows_saved = sum(tw * sk), tows_saved_pct = sum(tw * sk) / sum(tw * tv$tows),
    tows_saved_follow_up_pct = sum((tw * sk)[tv$follow_up]) / sum((tw * tv$tows)[tv$follow_up]))
}
# Random: the same number of dives skipped at random within the same pool
# (weak-scar dives, or all non-strong dives), as an exact expectation.
measure_random <- function(skip, dw, ns) {
  f <- sum(dw * skip) / sum(dw * ns)
  hm <- sum(dw * (ns & dv$high))
  c(dives_avoided = sum(dw * skip), dives_avoided_pct = sum(dw * skip) / sum(dw),
    high_missed = f * hm, high_missed_pct = f * hm / sum(dw * dv$high),
    inter_missed = f * sum(dw * (ns & dv$inter)),
    inter_missed_pct = f * sum(dw * (ns & dv$inter)) / sum(dw * dv$inter),
    low_avoided = f * sum(dw * (ns & dv$outcome == "low")),
    low_avoided_pct = f * sum(dw * (ns & dv$outcome == "low")) / sum(dw * (dv$outcome == "low")),
    ghost_avoided = f * sum(dw * (ns & dv$outcome == "ghost")),
    ghost_avoided_pct = f * sum(dw * (ns & dv$outcome == "ghost")) / sum(dw * (dv$outcome == "ghost")),
    ghost_per_high_missed = sum(dw * (ns & dv$outcome == "ghost")) / hm,
    share_deferred_ghost = sum(dw * (ns & dv$outcome == "ghost")) / sum(dw * ns),
    share_deferred_low = sum(dw * (ns & dv$outcome == "low")) / sum(dw * ns),
    share_deferred_intermediate = sum(dw * (ns & dv$inter)) / sum(dw * ns),
    share_deferred_high = hm / sum(dw * ns),
    cpue_missed_high_mean = sum((dw * dv$cpue)[ns & dv$high]) / hm)
}
all_metrics <- function(w) {
  dw <- w[dv_ry]
  tw <- w[tv_ry]
  out <- list()
  for (scope in c("all visits", "follow-up visits")) {
    sel <- if (scope == "all visits") rep(TRUE, nrow(dv)) else dv$follow_up
    tsel <- if (scope == "all visits") rep(TRUE, nrow(tv)) else tv$follow_up
    for (k in seq_along(variants)) {
      v <- measure(variants[[k]], w, dw * sel, tw * tsel, uses_zone[k])
      out[[paste(scope, names(variants)[k], sep = "#")]] <- v
      if (names(variants)[k] %in% names(random_of))
        out[[paste(scope, random_of[[names(variants)[k]]], sep = "#")]] <- measure_random(variants[[k]], dw * sel, if (startsWith(names(variants)[k], "weak")) weak else dv$nonstrong)
    }
  }
  unlist(lapply(names(out), function(n) stats::setNames(out[[n]], paste(n, names(out[[n]]), sep = "#"))))
}
point <- all_metrics(rep(1, length(reef_years)))
boot <- vapply(seq_len(boot_reps), function(r) {
  m <- all_metrics(tabulate(sample.int(length(reef_years), replace = TRUE), length(reef_years)))
  m[names(point)]
}, numeric(length(point)))
parts <- do.call(rbind, strsplit(names(point), "#", fixed = TRUE))
stats_table <- tibble::tibble(scope = parts[, 1], variant = parts[, 2], metric = parts[, 3], value = point,
                              per_season = ifelse(parts[, 3] %in% c("dives_avoided", "high_missed", "inter_missed", "low_avoided",
                                                                    "ghost_avoided", "tows_saved"), point / seasons, NA_real_),
                              low = apply(boot, 1, stats::quantile, 0.025, na.rm = TRUE),
                              high = apply(boot, 1, stats::quantile, 0.975, na.rm = TRUE))
write_table(stats_table, cfg, "protocol_statistics.csv")

totals <- tibble::tibble(seasons = seasons, visits = nrow(tv), follow_up_visits = sum(tv$follow_up), tows = sum(tv$tows),
                         dives = nrow(dv), high = sum(dv$high), inter = sum(dv$inter), low = sum(dv$outcome == "low"),
                         ghost = sum(dv$outcome == "ghost"), cpue_high_mean = mean(dv$cpue[dv$high]),
                         cpue_high_median = stats::median(dv$cpue[dv$high]), unranked_dives = sum(dv$rank < 0))
write_table(totals, cfg, "protocol_totals.csv")
print(totals, width = Inf)

# Outcome shares among non-strong dives the filter skips vs keeps (the trial's
# decision rule compares these).
filter_rates <- dplyr::bind_rows(lapply(names(flagged_ns), function(f) {
  dv[dv$nonstrong, ] |>
    dplyr::mutate(filter = f, skipped = flagged_ns[[f]][dv$nonstrong]) |>
    dplyr::group_by(.data$filter, .data$skipped) |>
    dplyr::summarise(dives = dplyr::n(), high_share = mean(.data$high), inter_share = mean(.data$inter),
                     low_share = mean(.data$outcome == "low"), ghost_share = mean(.data$outcome == "ghost"), .groups = "drop")
}))
write_table(filter_rates, cfg, "protocol_filter_rates.csv")
print(filter_rates, width = Inf)

# --- Dives above 0.04 the zone would miss: where, when, and whether an earlier
# visit that season had already flagged the site (strong tow or dive > 0.04) ------
strong_sites <- stats::setNames(lapply(guide, function(g) sites_near(g[tows$strong[g]])), vis$visit)
zm <- dv[dv$zone_lost & dv$high, ]
zm$visit_first <- vis$first[match(zm$visit, vis$visit)]
zm$zone_rank <- ctx$sv$zone_combined[match(paste(zm$visit, zm$site_row), paste(ctx$sv$visit, ctx$sv$site_row))]
zm$class <- cut(zm$rank, c(-Inf, 0.2, 0.34, 0.5, 0.8, Inf), c("E", "D", "C", "B", "A"), right = FALSE)
flag <- t(vapply(seq_len(nrow(zm)), function(j) {
  ev <- vis$visit[vis$reef_name == zm$reef_name[j] & vis$year == zm$year[j] & vis$first < zm$visit_first[j]]
  ed <- dv[dv$visit %in% ev & dv$site_row == zm$site_row[j], ]
  c(earlier_visits = length(ev),
    days_since_first_visit = as.numeric(zm$visit_first[j] - min(vis$first[vis$visit %in% c(ev, zm$visit[j])])),
    earlier_strong_tow = any(vapply(strong_sites[ev], function(s) zm$site_row[j] %in% s, logical(1))),
    earlier_dived = nrow(ed) > 0, earlier_max_cpue = if (nrow(ed)) max(ed$cpue) else NA_real_)
}, numeric(5)))
zm <- cbind(zm, flag)
zone_missed <- zm[order(-zm$total_cots), c("reef_name", "year", "voyage", "visit_first", "survey_date", "site_name", "cpue",
                                          "total_cots", "bottom_time", "zone_rank", "class", "earlier_visits",
                                          "days_since_first_visit", "earlier_strong_tow", "earlier_dived", "earlier_max_cpue", "memory_lost")]
write_table(zone_missed, cfg, "protocol_zone_missed_dives.csv")
zone_missed_reefs <- zone_missed |>
  dplyr::group_by(.data$reef_name, .data$year) |>
  dplyr::summarise(dives = dplyr::n(), visits = dplyr::n_distinct(.data$voyage), sites = dplyr::n_distinct(.data$site_name),
                   cots = sum(.data$total_cots), mean_cpue = mean(.data$cpue), max_cpue = max(.data$cpue),
                   earlier_strong_tow = mean(.data$earlier_strong_tow == 1), earlier_dived = mean(.data$earlier_dived == 1),
                   missed_with_season_memory = sum(.data$memory_lost),
                   .groups = "drop") |>
  dplyr::mutate(cots_share = .data$cots / sum(.data$cots)) |>
  dplyr::arrange(dplyr::desc(.data$cots))
write_table(zone_missed_reefs, cfg, "protocol_zone_missed_reefs.csv")
print(zone_missed_reefs, n = 15, width = Inf)
print(tibble::as_tibble(utils::head(zone_missed, 20)), width = Inf)
top <- zone_missed$reef_name == zone_missed_reefs$reef_name[1] & zone_missed$year == zone_missed_reefs$year[1]
log_message("Zone-missed > 0.04: ", nrow(zone_missed), " dives, mean CPUE ", round(mean(zone_missed$cpue), 3),
            "; without the top reef-year ", round(mean(zone_missed$cpue[!top]), 3),
            "; earlier strong tow at the site ", round(mean(zone_missed$earlier_strong_tow == 1), 2),
            "; dived earlier ", round(mean(zone_missed$earlier_dived == 1), 2),
            "; earlier dive > 0.04 ", round(mean(zone_missed$earlier_max_cpue > 0.04, na.rm = TRUE) * mean(zone_missed$earlier_dived == 1), 2))

# --- One voyage: what a crew would see (2024-2026 seasons) -------------------------
fl <- flagged_ns[[trial_filter]]
voy <- dv[dv$year %in% trial_years, ] |>
  dplyr::mutate(flag = fl[dv$year %in% trial_years]) |>
  dplyr::group_by(.data$voyage) |>
  dplyr::summarise(reefs = dplyr::n_distinct(.data$reef_name), dives = dplyr::n(),
                   high = sum(.data$high), ghost = sum(.data$outcome == "ghost"),
                   weak = sum(.data$pool == "weak"), no_evidence = sum(.data$pool %in% c("negative", "no tow")), flagged = sum(.data$flag),
                   flagged_high = sum(.data$flag & .data$high), flagged_ghost = sum(.data$flag & .data$outcome == "ghost"),
                   flagged_low = sum(.data$flag & .data$outcome == "low"),
                   zone_lost_high = sum(.data$zone_lost & .data$high),
                   memory_lost_high = sum(.data$memory_lost & .data$high), .groups = "drop")
voy_tows <- tv[tv$year %in% trial_years, ] |>
  dplyr::group_by(.data$voyage) |>
  dplyr::summarise(follow_up_reefs = sum(.data$follow_up), follow_up_tows = sum(.data$tows[.data$follow_up]),
                   all_tows = sum(.data$tows), tows_saved = sum(.data$skipped),
                   tows_saved_memory = sum(.data$skipped_mem), .groups = "drop")
voy <- dplyr::inner_join(voy, voy_tows, by = "voyage")
voyage_summary <- tidyr::pivot_longer(voy, -"voyage", names_to = "quantity") |>
  dplyr::group_by(.data$quantity) |>
  dplyr::summarise(voyages = dplyr::n(), mean = mean(.data$value), p10 = stats::quantile(.data$value, 0.1),
                   median = stats::median(.data$value), p90 = stats::quantile(.data$value, 0.9), .groups = "drop")
write_table(voyage_summary, cfg, "protocol_voyage_summary.csv")
print(voyage_summary, n = Inf)

# How often k pooled voyages show the trial filter working: flagged non-strong
# dives with a lower > 0.04 rate, and a higher ghost rate, than kept ones.
all_voy <- dv |>
  dplyr::mutate(flag = fl) |>
  dplyr::filter(.data$nonstrong) |>
  dplyr::group_by(.data$voyage) |>
  dplyr::summarise(f_n = sum(.data$flag), f_high = sum(.data$flag & .data$high), f_ghost = sum(.data$flag & .data$outcome == "ghost"),
                   k_n = sum(!.data$flag), k_high = sum(!.data$flag & .data$high), k_ghost = sum(!.data$flag & .data$outcome == "ghost"),
                   .groups = "drop")
trial <- dplyr::bind_rows(lapply(pooled_voyages, function(k) {
  r <- vapply(seq_len(trial_reps), function(b) {
    s <- colSums(all_voy[sample.int(nrow(all_voy), k, replace = TRUE), -1])
    c(none_flagged = s[["f_n"]] == 0,
      lower_high_rate = s[["f_n"]] > 0 && s[["k_n"]] > 0 && s[["f_high"]] / s[["f_n"]] < s[["k_high"]] / s[["k_n"]],
      higher_ghost_rate = s[["f_n"]] > 0 && s[["k_n"]] > 0 && s[["f_ghost"]] / s[["f_n"]] > s[["k_ghost"]] / s[["k_n"]],
      flagged = s[["f_n"]])
  }, numeric(4))
  tibble::tibble(voyages = k, p_none_flagged = mean(r[1, ]), p_lower_high_rate = mean(r[2, ]),
                 p_higher_ghost_rate = mean(r[3, ]), flagged_median = stats::median(r[4, ]),
                 flagged_p10 = stats::quantile(r[4, ], 0.1), flagged_p90 = stats::quantile(r[4, ], 0.9))
}))
write_table(trial, cfg, "protocol_trial_power.csv")
print(trial, width = Inf)

show <- c("tows_saved_pct", "tows_saved_follow_up_pct", "dives_avoided", "dives_avoided_pct", "high_missed", "high_missed_pct",
          "inter_missed", "low_avoided", "low_avoided_pct", "ghost_avoided", "ghost_avoided_pct", "ghost_per_high_missed",
          "bottom_time_saved_pct", "cots_forgone_pct", "cpue_missed_high_mean", "cpue_missed_high_median", "cpue_kept_high_mean")
print(tidyr::pivot_wider(stats_table[stats_table$scope == "all visits" & stats_table$metric %in% show, c("variant", "metric", "value")],
                         names_from = "metric", values_from = "value"), n = Inf, width = Inf)
print(stats_table[stats_table$scope == "all visits" & stats_table$metric %in% c("cpue_missed_high_mean", "high_missed_pct",
                                                                                 "ghost_avoided_pct", "low_avoided_pct"), ], n = Inf)
