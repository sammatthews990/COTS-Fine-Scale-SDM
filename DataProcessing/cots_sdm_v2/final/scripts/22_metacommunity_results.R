# Question: what does the latest archived metacommunity evidence support?
# Method: summarise the five-seed candidate comparison, historical-culling
# diagnostic and single-seed 2026 extension from frozen ADRIA run outputs.
# Validation: calibration and exploratory diagnostics only, not held-out tests.
# Observed catches/targets are not copied; reported aggregates come from audits.
# Run from the repository root, then 07 --metacommunity_only=1 for figures.

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "..", "R", "common.R"))
cfg <- initialise_final()
keys <- grep("^mc_", names(cfg$inputs), value = TRUE)
stopifnot(!anyDuplicated(names(cfg$inputs)))
paths <- stats::setNames(lapply(keys, function(k) checked_input(cfg, k)), keys)
read_mc <- function(k) {
  d <- utils::read.csv(paths[[k]], stringsAsFactors = FALSE)
  for (n in names(d)) {
    if (is.character(d[[n]]) && all(tolower(d[[n]]) %in% c("true", "false")))
      d[[n]] <- tolower(d[[n]]) == "true"
  }
  d
}
summarise_stats <- function(d, groups, values) {
  dplyr::summarise(dplyr::group_by(d, dplyr::across(dplyr::all_of(groups))),
    dplyr::across(dplyr::all_of(values), list(mean = mean, sd = stats::sd)), .groups = "drop")
}

audit <- jsonlite::fromJSON(paths[["mc_candidates_candidates_audit_json"]])
stopifnot(audit$status == "bounded_diagnostic_not_promoted",
          all(unlist(audit$reference_replay_max_rel_delta) == 0))
a <- read_mc("mc_candidates_candidates_audit")
a <- a[a$candidate %in% c("W", "C3a", "C6"), ]
stopifnot(nrow(a) == 15, all(table(a$candidate) == 5), all(a$q == 0.015),
          all(a$matched_v2 == 7), all(a$observed_v2 == 7))
s <- summarise_stats(a, "candidate", c("summed_loss_v2", "coral_J", "liz_peak_adults_1995_2000"))
counts <- dplyr::summarise(dplyr::group_by(a, .data$candidate),
  seeds = dplyr::n(), matched_peaks = min(.data$matched_v2), observed_peaks = min(.data$observed_v2),
  gate10_passes = sum(.data$g2_pass_ahalf10), gate30_passes = sum(.data$g2_pass_ahalf30), .groups = "drop")
s <- dplyr::left_join(s, counts, by = "candidate")
s$run_id <- "20261014T_48_candidates_1991_v2"
s$status <- "replicated_calibration_not_promoted"
write_table(s, cfg, "metacommunity_candidates.csv")
write_table(a, cfg, "metacommunity_candidate_replicates.csv")
p <- read_mc("mc_candidates_candidates_replicates")
p <- p[p$candidate %in% s$candidate, ]
stopifnot(nrow(p) == 60)
write_table(summarise_stats(p, c("candidate", "reef_name"),
  c("first_peak_year", "first_peak_cpue", "second_peak_year", "second_peak_cpue")),
  cfg, "metacommunity_peak_summary.csv")

traj <- read_mc("mc_candidates_trajectories")
traj <- dplyr::inner_join(traj, a[, c("candidate", "seed", "treatment")], by = "treatment")
stopifnot(nrow(traj) == 15 * 4 * 34, all(abs(traj$simulated_cpue - 0.015 * traj$adults_ha) < 1e-8))
write_table(summarise_stats(traj[traj$year >= 1992, ], c("candidate", "reef_name", "year"),
  c("simulated_cpue", "adults_ha", "coral_cover")), cfg, "metacommunity_trajectories.csv")

# Manta observations for the plotted historical comparison. This reads the same
# archived series as the old report figure; do not copy per-survey raw records.
obs <- read_mc("mc_observations")
labels <- c("Lizard Isles" = "Lizard Island Reef", "Macgillivray Reef" = "MacGillivray Reef",
            "North Direction Island" = "North Direction Reef", "Eyrie Reef" = "Eyrie Reef")
obs <- obs[obs$reef_name %in% names(labels) & obs$year >= 1992 & obs$year <= 2024, ]
obs$reef_name <- unname(labels[obs$reef_name])
# Plotting observations remain in the git-ignored cache; only window maxima
# needed for reported comparisons are written to the results tables.
saveRDS(obs[, c("reef_name", "year", "cotsptow")], file.path(cfg$dirs$cache, "metacommunity_observations.rds"))
observed <- dplyr::summarise(dplyr::group_by(obs, .data$reef_name),
  first_wave_max_cpue = max(.data$cotsptow[.data$year <= 2005]),
  second_wave_max_cpue = max(.data$cotsptow[.data$year >= 2010 & .data$year <= 2017]), .groups = "drop")
write_table(observed, cfg, "metacommunity_observed_maxima.csv")

culls <- read_mc("mc_culling_culling_audit")
write_table(culls, cfg, "metacommunity_culling_diagnostic.csv")
rem <- read_mc("mc_culling_culling_removals")
write_table(rem, cfg, "metacommunity_culling_removal_summary.csv")

# V0.4 is a separate single-seed diagnostic. Compare score v3 only within this
# experiment; never compare its numerical value with the v2 candidate loss.
score_rows <- list(); trajectory_rows <- list()
for (dataset in c("original", "corrected")) {
  d <- read_mc(paste0("mc_", dataset, "_score_detectability_scores"))
  d <- d[d$a_half == 0 & d$score_v2_primary, ]
  stopifnot(nrow(d) == 24, all(table(d$treatment) == 4), all(is.finite(d$loss_v3)))
  z <- dplyr::summarise(dplyr::group_by(d, .data$treatment),
    loss_v2 = sum(.data$loss_v2), loss_v3 = sum(.data$loss_v3),
    diver_error = sum(.data$diver_error), .groups = "drop")
  z$dataset <- dataset
  z$coral <- ifelse(grepl("_t60_", z$treatment), "revised thermal settings", "earlier thermal settings")
  z$larval_index <- sub(".*_(chloff|mb15|mb30|pb15|pb30)_cv2_pnd$", "\\1", z$treatment)
  gates <- read_mc(paste0("mc_", dataset, "_gate_gate_g2prime"))
  g <- gates[gates$a_half == 30, c("treatment", "passes_g2prime")]
  z <- dplyr::left_join(z, g, by = "treatment")
  z$seeds <- 1L
  z$status <- "single_seed_diagnostic_not_promoted"
  score_rows[[dataset]] <- z
  t <- read_mc(paste0("mc_", dataset, "_trajectories"))
  t$dataset <- dataset
  trajectory_rows[[dataset]] <- t
}
latest <- dplyr::bind_rows(score_rows)
stopifnot(nrow(latest) == 12, !anyNA(latest$passes_g2prime))
write_table(latest, cfg, "metacommunity_2026_diagnostics.csv")
latest_traj <- dplyr::bind_rows(trajectory_rows)
latest_traj <- latest_traj[latest_traj$year >= 2018 & grepl("_chloff_", latest_traj$treatment), ]
stopifnot(nrow(latest_traj) == 2 * 2 * 4 * 9,
          identical(range(latest_traj$year), c(2018L, 2026L)))
latest_traj$coral <- ifelse(grepl("_t60_", latest_traj$treatment), "revised thermal settings", "earlier thermal settings")
write_table(latest_traj, cfg, "metacommunity_2026_trajectories.csv")

manifest <- data.frame(input = keys,
  path = vapply(paths, function(p) sub("^.*ADRIA[.]jl/", "ADRIA.jl/", gsub("\\\\", "/", p)), ""),
  sha256 = vapply(paths, sha256_file, ""), snapshot = "2026-10-15")
write_table(manifest, cfg, "metacommunity_input_manifest.csv")
