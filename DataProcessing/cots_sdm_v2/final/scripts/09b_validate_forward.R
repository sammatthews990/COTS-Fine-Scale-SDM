# Strictly forward validation of habitat, reef-context, history and combined
# scores on Cull site-years. For each test year T every score uses only data
# from years < T:
#   habitat          clean_joint Hmsc refit on rows < T (09a, variant base)
#   habitat+context  same with the 0.5-1 km reef-context terms (variant context)
#   history          kernel-weighted prior high-density rate (active reef-years < T)
#   effort           dives at the same site in years < T (operator-targeting baseline)
#   combined         habitat prior updated by history: (kappa p + sum w y) / (kappa + W)
# Kernel width sigma, recency tau and kappa are chosen for year T using year
# T - 1 only (nested); a fixed pre-registered setting is reported alongside.
# Units: first-visit site-years in active reef-years with both outcomes.

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "..", "R", "common.R"))
cfg <- initialise_final()
source(file.path(cfg$final_dir, "R", "forward.R"))
test_years <- 2022:2026
cutoffs <- 2021:2026
grid_sigma <- c(250, 500, 1000)
grid_tau <- c(1, 3, Inf)
grid_kappa <- c(0.5, 1, 2, 4, 8)
fixed <- list(sigma = cfg$history$sigma_m, tau = cfg$history$tau_years, kappa = cfg$history$kappa)

# --- Site-years ---------------------------------------------------------------
sy <- readRDS(file.path(cfg$exploration_dir, "artifacts", "site_year_first_visit.rds")) |>
  dplyr::filter(!is.na(.data$high_density))
xy <- sf::st_coordinates(sf::st_transform(
  sf::st_as_sf(sy, coords = c("longitude", "latitude"), crs = 4326), 7855))
sy$x <- xy[, 1]
sy$y <- xy[, 2]
sy$reef_year <- paste(sy$reef_name, sy$year)
sy <- sy |>
  dplyr::group_by(.data$reef_year) |>
  dplyr::mutate(active = any(.data$high_density == 1), mixed = active & any(.data$high_density == 0)) |>
  dplyr::ungroup()

template <- terra::rast(file.path(cfg$dirs$cache, "base_support_120m.tif"))[[1]]
sy$cell <- terra::cellFromXY(template, cbind(sy$x, sy$y))
grid <- readRDS(file.path(cfg$dirs$cache, "support_grid_120m.rds"))
site_predictors <- grid[match(sy$cell, grid$cell), setdiff(names(grid), c("cell", "x", "y")), drop = FALSE]
context_raw <- terra::extract(terra::rast(file.path(cfg$dirs$cache, "reef_context_120m.tif")), sy$cell)
rm(grid)

# --- Habitat scores from forward refits -------------------------------------
habitat_score <- function(variant, cutoff, rows) {
  fit <- readRDS(file.path(cfg$dirs$models, "forward", sprintf("clean_joint_%s_before%d.rds", variant, cutoff)))
  d <- site_predictors[rows, , drop = FALSE]
  if (variant == "context") d <- cbind(d, context_terms(context_raw[rows, , drop = FALSE], fit$artifact$context_impute))
  design <- component_design(d, fit$artifact, fit$model)
  out <- rep(NA_real_, length(rows))
  target <- match(cfg$hmsc$target_response, fit$artifact$responses)
  beta <- target_beta_draws(fit$model, target, 200L, cfg$hmsc$seed)
  if (length(design$rows)) out[design$rows] <- summarise_target_posterior(design$X, beta)[, "mean"]
  out
}

# --- Score every year that can serve as a test or a selection year -----------
score_year <- function(T) {
  test <- which(sy$year == T & sy$mixed)
  prior <- which(sy$year < T & sy$active)
  d <- sy[test, c("site_id", "reef_name", "reef_year", "year", "high_density", "x", "y")]
  d$habitat <- habitat_score("base", T, test)
  d$habitat_context <- habitat_score("context", T, test)
  same_site_dives <- tapply(sy$n_dives[sy$year < T], sy$site_id[sy$year < T], sum)
  d$effort <- dplyr::coalesce(as.numeric(same_site_dives[d$site_id]), 0) + stats::runif(nrow(d), 0, 1e-6)
  d$same_site_history <- d$site_id %in% sy$site_id[prior]
  prior_xy <- cbind(sy$x[prior], sy$y[prior])
  base_rate <- mean(sy$high_density[prior])
  kernels <- list()
  for (s in grid_sigma) for (t in grid_tau) {
    h <- history_kernel(prior_xy, sy$high_density[prior], (T - 1) - sy$year[prior], cbind(d$x, d$y), s, t)
    kernels[[paste(s, t)]] <- h
  }
  list(year = T, d = d, kernels = kernels, base_rate = base_rate)
}
set.seed(cfg$hmsc$seed)
scored <- lapply(cutoffs, score_year)
names(scored) <- cutoffs

within_auc <- function(score, truth, group) {
  ok <- is.finite(score)
  within_group_auc(score[ok], truth[ok] == 1, group[ok])
}
history_score <- function(s, sigma, tau) {
  h <- s$kernels[[paste(sigma, tau)]]
  ifelse(h$W > 0, h$H, s$base_rate)
}
combined_score <- function(s, sigma, tau, kappa, habitat = s$d$habitat) {
  h <- s$kernels[[paste(sigma, tau)]]
  combine_history(habitat, h$H, h$W, kappa)
}

# --- Nested selection on T - 1 ---------------------------------------------------
settings <- expand.grid(sigma = grid_sigma, tau = grid_tau, kappa = grid_kappa)
selection <- dplyr::bind_rows(lapply(test_years, function(T) {
  s <- scored[[as.character(T - 1)]]
  hist_auc <- vapply(seq_len(nrow(settings)), function(i) within_auc(
    history_score(s, settings$sigma[i], settings$tau[i]), s$d$high_density, s$d$reef_year), numeric(1))
  comb_auc <- vapply(seq_len(nrow(settings)), function(i) within_auc(
    combined_score(s, settings$sigma[i], settings$tau[i], settings$kappa[i]), s$d$high_density, s$d$reef_year), numeric(1))
  comb_ctx <- vapply(seq_len(nrow(settings)), function(i) within_auc(
    combined_score(s, settings$sigma[i], settings$tau[i], settings$kappa[i], s$d$habitat_context),
    s$d$high_density, s$d$reef_year), numeric(1))
  h <- which.max(hist_auc); c1 <- which.max(comb_auc); c2 <- which.max(comb_ctx)
  tibble::tibble(test_year = T,
                 history_sigma = settings$sigma[h], history_tau = settings$tau[h],
                 combined_sigma = settings$sigma[c1], combined_tau = settings$tau[c1], combined_kappa = settings$kappa[c1],
                 context_sigma = settings$sigma[c2], context_tau = settings$tau[c2], context_kappa = settings$kappa[c2])
}))
write_table(selection, cfg, "forward_nested_selection.csv")

# --- Test-year scores -----------------------------------------------------------
test_scores <- dplyr::bind_rows(lapply(test_years, function(T) {
  s <- scored[[as.character(T)]]
  sel <- selection[selection$test_year == T, ]
  s$d |> dplyr::mutate(
    history = history_score(s, sel$history_sigma, sel$history_tau),
    history_fixed = history_score(s, fixed$sigma, fixed$tau),
    combined = combined_score(s, sel$combined_sigma, sel$combined_tau, sel$combined_kappa),
    combined_fixed = combined_score(s, fixed$sigma, fixed$tau, fixed$kappa),
    combined_context = combined_score(s, sel$context_sigma, sel$context_tau, sel$context_kappa, s$d$habitat_context),
    history_weight = s$kernels[[paste(fixed$sigma, fixed$tau)]]$W
  )
}))
saveRDS(test_scores, file.path(cfg$dirs$cache, "forward_test_scores.rds"))

methods <- c(
  habitat = "Habitat (HMSC, refit < T)",
  habitat_context = "Habitat + reef context",
  effort = "Effort only (prior dives at site)",
  history = "History kernel (nested tuning)",
  history_fixed = "History kernel (fixed)",
  combined = "Habitat + history (nested tuning)",
  combined_fixed = "Habitat + history (fixed, deployed)",
  combined_context = "Habitat + context + history (nested)"
)
# Every method is scored on the same rows: those with a habitat prediction.
eval_rows <- is.finite(test_scores$habitat) & is.finite(test_scores$habitat_context)
ev <- test_scores[eval_rows, ]
subsets <- list(
  all = rep(TRUE, nrow(ev)),
  `with nearby history` = ev$history_weight > 0.5,
  `little or no nearby history` = ev$history_weight <= 0.5,
  `new sites` = !ev$same_site_history  # no Cull record at the same site before T
)
groups <- split(seq_len(nrow(ev)), ev$reef_year)
boot <- replicate(500, sample(length(groups), replace = TRUE), simplify = FALSE)
summary <- dplyr::bind_rows(lapply(names(subsets), function(sub) {
  keep <- subsets[[sub]]
  dplyr::bind_rows(lapply(names(methods), function(m) {
    score <- ev[[m]]
    est <- within_auc(score[keep], ev$high_density[keep], ev$reef_year[keep])
    ci <- if (sub == "all") stats::quantile(vapply(boot, function(g) {
      i <- unlist(groups[g], use.names = FALSE)
      grp <- rep(seq_along(g), lengths(groups[g]))
      within_auc(score[i], ev$high_density[i], grp)
    }, numeric(1)), c(0.025, 0.975), na.rm = TRUE) else c(NA, NA)
    tibble::tibble(
      subset = sub, method = m, label = methods[[m]],
      site_years = sum(keep), reef_years = dplyr::n_distinct(ev$reef_year[keep]),
      within_reef_year_auc = est, auc_low = ci[1], auc_high = ci[2],
      top20_lift = within_group_top20_lift(score[keep], ev$high_density[keep] == 1, ev$reef_year[keep])
    )
  }))
}))
write_table(summary, cfg, "forward_validation_summary.csv")

by_year <- dplyr::bind_rows(lapply(test_years, function(T) {
  e <- ev[ev$year == T, ]
  dplyr::bind_rows(lapply(names(methods), function(m) tibble::tibble(
    test_year = T, method = m, site_years = nrow(e), reef_years = dplyr::n_distinct(e$reef_year),
    within_reef_year_auc = within_auc(e[[m]], e$high_density, e$reef_year))))
}))
write_table(by_year, cfg, "forward_validation_by_year.csv")

forward_fits <- dplyr::bind_rows(lapply(list.files(file.path(cfg$dirs$models, "forward"), full.names = TRUE), function(f) {
  x <- readRDS(f)
  tibble::tibble(variant = x$variant, cutoff = x$cutoff, rows = nrow(x$artifact$data),
                 min_ess = x$min_ess, max_psrf_upper = x$max_psrf_upper)
}))
write_table(forward_fits, cfg, "forward_habitat_fits.csv")
print(summary[summary$subset == "all", c("label", "site_years", "within_reef_year_auc", "auc_low", "auc_high", "top20_lift")])
print(tidyr::pivot_wider(summary[, c("subset", "method", "within_reef_year_auc")], names_from = "subset", values_from = "within_reef_year_auc"))
