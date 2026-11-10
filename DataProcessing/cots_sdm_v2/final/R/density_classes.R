# Density-aware priority classes (20b, 36): reef density from recent dives,
# manta tows and the GBR-wide reef outbreak forecast, and scaled class quotas.

reef_id <- function(x) sub(".*[(]([0-9]+-[0-9]+[a-z]?)[)][[:space:]]*$", "\\1", x)

# Stewart et al. hybrid forecast, one row per reef and target year; the 2026
# file wins any overlap with the hindcasts.
read_reef_forecast <- function(cfg) {
  dplyr::bind_rows(utils::read.csv(checked_input(cfg, "reef_forecast_2026")),
                   utils::read.csv(checked_input(cfg, "reef_forecast_hindcast"))) |>
    dplyr::transmute(reef_id = reef_id(.data$reefName), year = as.integer(.data$year), forecast = .data$outbrProb) |>
    dplyr::distinct(.data$reef_id, .data$year, .keep_all = TRUE)
}

# Reef-level evidence from dives `d` and tows `t` (one window), with the forecast.
reef_evidence <- function(reefs, d, t, forecast_year, forecast) {
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
density_variants <- list(
  evidence_forecast = problematic ~ pct + has_dives + dive_rate + log_dives + has_tows + tow_rate + forecast,
  evidence = problematic ~ pct + has_dives + dive_rate + log_dives + has_tows + tow_rate,
  forecast = problematic ~ pct + forecast)

# Every designated site of a reef-year, at its within-reef percentile.
site_grid <- function(ry) {
  n <- ry$sites_on_reef
  idx <- rep(seq_len(nrow(ry)), n)
  out <- ry[idx, ]
  out$pct <- unlist(lapply(n, function(k) (seq_len(k) - 1) / max(k - 1, 1)))
  out
}

# Unknown reefs (no dives or tows in the window): P scaled by the chance such a
# reef needed culling at all (`um`, logistic of "needed" on the logit forecast).
discount_unknown <- function(lp, d, um) {
  unknown <- d$has_dives == 0 & d$has_tows == 0
  w <- stats::predict(um, data.frame(forecast = d$forecast[unknown]), type = "response")
  lp[unknown] <- stats::qlogis(w * stats::plogis(lp[unknown]))
  lp
}

# Scaled quotas: share in classes up to k = plogis(qlogis(break_k) + delta),
# delta = qlogis(reef density) - qlogis(typical culled reef). As one score: a
# site at top fraction t = (rank - 0.5) / n of its reef scores delta - qlogis(t)
# and is in class k or better when the score >= -qlogis(break_k).
scaled_score <- function(delta, rank, n) delta - stats::qlogis((rank - 0.5) / n)
scaled_share <- function(share, delta) stats::plogis(stats::qlogis(share) + delta)
reef_density <- function(lp, reef) tapply(stats::plogis(lp), reef, mean)
rank_of <- function(pct, n) n - pct * (n - 1)

# --- Alternative reef densities (36) -------------------------------------------------
# Evidence at date `end` (records strictly before it, within `window` days),
# each dive and tow weighted 0.5^(age / half_life); half_life = Inf weights the
# window evenly. Weighted counts feed reef_terms() unchanged.
weighted_evidence <- function(reef, end, year, dives_by_reef, tows_by_reef, half_life, forecast, window = 365) {
  w <- function(x, i) {
    if (is.null(x)) return(numeric(0))
    age <- as.numeric(end[i] - x$survey_date)
    ifelse(age > 0 & age <= window, 0.5^(age / half_life), 0)
  }
  fc <- forecast$forecast[match(paste(reef_id(reef), year), paste(forecast$reef_id, forecast$year))]
  out <- t(vapply(seq_along(reef), function(i) {
    d <- dives_by_reef[[reef[i]]]
    tw <- tows_by_reef[[reef[i]]]
    wd <- w(d, i)
    wt <- w(tw, i)
    c(dives = sum(wd), problematic = sum(wd * (if (is.null(d)) 0 else d$outcome == "problematic")),
      tows = sum(wt), strong = sum(wt * (if (is.null(tw)) 0 else tw$strong)))
  }, numeric(4)))
  data.frame(reef_name = reef, out, forecast_prob = fc)
}

# Survey occasions: a reef's dive and tow dates split where they are more than
# `gap` days apart. `needed`: any problematic dive in the occasion.
survey_occasions <- function(dives, tows, reefs, gap = 10) {
  rec <- rbind(data.frame(reef_name = dives$reef_name, date = dives$survey_date, prob = dives$outcome == "problematic"),
               data.frame(reef_name = tows$reef_name, date = tows$survey_date, prob = FALSE))
  rec <- rec[rec$reef_name %in% reefs, ]
  rec <- rec[order(rec$reef_name, rec$date), ]
  new <- c(TRUE, rec$reef_name[-1] != rec$reef_name[-nrow(rec)] | diff(as.numeric(rec$date)) > gap)
  rec$occ <- cumsum(new)
  dplyr::summarise(dplyr::group_by(rec, .data$occ), reef_name = dplyr::first(.data$reef_name), start = min(.data$date),
                   needed = as.numeric(any(.data$prob)), .groups = "drop") |>
    dplyr::mutate(year = as.integer(format(.data$start, "%Y")))
}

# Site evidence at each visit's planning date (36): dives at designated sites
# within `radius` m of the site and tows linked to it (34's 100 m links).
# Recent: weighted 0.5^(age / half_life) within a year; year: unweighted within
# a year; long: one to `years` years before, unweighted. One row per `sv` row.
site_history <- function(vis, sv, sites_xy, site_reef, dives, tows, tow_sites, radius = 300, half_life = 90, years = 6) {
  cols <- c("d_recent", "p_recent", "d_long", "p_long", "t_recent", "s_recent", "t_long", "s_long",
            "d_year", "p_year", "t_year", "s_year")
  out <- matrix(0, nrow(sv), length(cols), dimnames = list(NULL, cols))
  linked <- as.integer(names(tow_sites))
  tow_reef <- site_reef[vapply(tow_sites, function(j) if (length(j)) j[1] else NA_integer_, integer(1))]
  sv_visit <- match(sv$visit, vis$visit)
  for (r in unique(vis$reef_name)) {
    vi <- which(vis$reef_name == r)
    s_rows <- which(site_reef == r)
    rows <- which(vis$reef_name[sv_visit] == r)
    if (!length(s_rows) || !length(rows)) next
    vpos <- match(sv_visit[rows], vi)
    spos <- match(sv$site_row[rows], s_rows)
    ok <- !is.na(spos)
    xy <- sites_xy[s_rows, , drop = FALSE]
    nb <- Matrix::Matrix(as.matrix(stats::dist(xy)) <= radius, sparse = TRUE) * 1
    weights <- function(dates) {
      age <- outer(as.numeric(vis$plan[vi]), as.numeric(dates), "-")
      list(recent = ifelse(age > 0 & age <= 365, 0.5^(age / half_life), 0),
           year = (age > 0 & age <= 365) * 1, long = (age > 365 & age <= years * 365) * 1)
    }
    take <- function(m) as.matrix(m)[cbind(vpos[ok], spos[ok])]
    dd <- dives[dives$reef_name == r & !is.na(dives$site_row) & dives$site_row %in% s_rows, ]
    if (nrow(dd)) {
      a <- Matrix::sparseMatrix(i = seq_len(nrow(dd)), j = match(dd$site_row, s_rows), x = 1,
                                dims = c(nrow(dd), length(s_rows))) %*% nb
      ap <- a * (dd$outcome == "problematic")
      w <- weights(dd$survey_date)
      out[rows[ok], c("d_recent", "p_recent")] <- cbind(take(w$recent %*% a), take(w$recent %*% ap))
      out[rows[ok], c("d_year", "p_year")] <- cbind(take(w$year %*% a), take(w$year %*% ap))
      out[rows[ok], c("d_long", "p_long")] <- cbind(take(w$long %*% a), take(w$long %*% ap))
    }
    ti <- which(tow_reef == r)
    if (length(ti)) {
      jj <- lapply(tow_sites[ti], function(j) match(j[j %in% s_rows], s_rows))
      a <- Matrix::sparseMatrix(i = rep(seq_along(jj), lengths(jj)), j = unlist(jj), x = 1,
                                dims = c(length(jj), length(s_rows)))
      as_ <- a * tows$strong[linked[ti]]
      w <- weights(tows$survey_date[linked[ti]])
      out[rows[ok], c("t_recent", "s_recent")] <- cbind(take(w$recent %*% a), take(w$recent %*% as_))
      out[rows[ok], c("t_year", "s_year")] <- cbind(take(w$year %*% a), take(w$year %*% as_))
      out[rows[ok], c("t_long", "s_long")] <- cbind(take(w$long %*% a), take(w$long %*% as_))
    }
  }
  as.data.frame(out)
}
site_terms <- function(h) {
  lr <- function(k, n) ifelse(n > 1e-9, stats::qlogis((k + 0.5) / (n + 1)), 0)
  data.frame(has_dr = as.numeric(h$d_recent > 1e-9), rate_dr = lr(h$p_recent, h$d_recent),
             has_dl = as.numeric(h$d_long > 0), rate_dl = lr(h$p_long, h$d_long),
             has_tr = as.numeric(h$t_recent > 1e-9), rate_tr = lr(h$s_recent, h$t_recent),
             has_tl = as.numeric(h$t_long > 0), rate_tl = lr(h$s_long, h$t_long),
             log_dl = log1p(h$d_long), log_tl = log1p(h$t_long), log_dr = log1p(h$d_recent),
             searched_clean = h$t_year > 0 & h$s_year == 0 & h$p_year == 0)
}
