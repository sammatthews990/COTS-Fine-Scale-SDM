# Ranking designated cull sites and scoring the ranking against dive outcomes
# (used by 19 and 20; 15-17 carry their own copies).

# Percentile within a group, ties broken at random (0 = lowest, 1 = highest).
pct_rank <- function(v) {
  ok <- is.finite(v)
  out <- rep(NA_real_, length(v))
  out[ok] <- (rank(v[ok], ties.method = "random") - 1) / max(sum(ok) - 1, 1)
  out
}

# Model terms for the site-evidence model. `habitat` names the habitat rank
# column; `habitat2` optionally adds a second habitat rank.
site_evidence_terms <- function(f, habitat, habitat2 = NULL) {
  fill <- function(x) ifelse(is.finite(x), x, 0.5)
  out <- data.frame(habitat = fill(f[[habitat]]),
                    site_problematic_prev = f$site_problematic_prev, site_quiet_prev = f$site_quiet_prev,
                    cull_rate = f$cull_rate, log_cull_effort = log1p(f$cull_effort),
                    log_manta_strong = log1p(f$manta_strong), log_manta_effort = log1p(f$manta_effort))
  if (!is.null(habitat2)) out$habitat2 <- fill(f[[habitat2]])
  out
}
site_evidence_formula <- function(habitat = TRUE, habitat2 = FALSE) {
  rhs <- c(if (habitat) "habitat", if (habitat2) "habitat2", "site_problematic_prev", "site_quiet_prev",
           "cull_rate", "log_cull_effort", "log_manta_strong", "log_manta_effort")
  stats::as.formula(paste("problematic ~", paste(rhs, collapse = " + ")))
}

# Share of problematic dives kept, and of low-density and ghost dives avoided,
# when only sites ranked in the top fraction `a` of their reef are dived.
site_metrics <- function(d, col, a) {
  kept <- d[[col]] >= 1 - a - 1e-9
  c(problematic_dives_kept = sum(kept & d$outcome == "problematic") / sum(d$outcome == "problematic"),
    low_dives_avoided = sum(!kept & d$outcome == "low") / sum(d$outcome == "low"),
    ghost_dives_avoided = sum(!kept & d$outcome == "ghost") / sum(d$outcome == "ghost"),
    dives_kept = mean(kept))
}

# Sites searched and dives avoided when a ranking keeps a fixed share of
# problematic dives (linear interpolation along its curve over site fractions).
matched_retention <- function(d, col, targets, fractions = round(seq(0.05, 1, 0.05), 2)) {
  m <- vapply(fractions, function(a) site_metrics(d, col, a), numeric(4))
  pk <- c(0, m["problematic_dives_kept", ])
  curves <- list(sites_searched = c(0, fractions), low_dives_avoided = c(1, m["low_dives_avoided", ]),
                 ghost_dives_avoided = c(1, m["ghost_dives_avoided", ]))
  unlist(lapply(targets, function(t) {
    v <- vapply(curves, function(y) stats::approx(pk, y, xout = t, ties = mean)$y, numeric(1))
    stats::setNames(v, paste0(names(curves), "_at_", round(100 * t)))
  }))
}
