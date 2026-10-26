# Metacommunity model (COTSMod.jl + ADRIA) calibration summary for the
# technical report.
#
# Question: how well do the two current COTSMod configurations reproduce the
# observed COTS outbreak peaks on the four Lizard Island calibration reefs?
# Method: reads two archived calibration runs from the ADRIA.jl repository
# (sandbox/calibration/runs, branch cots-submodel) and the AIMS LTMP manta
# COTS-per-tow series they were scored against, without copying the LTMP data
# here. Writes a figure, a per-reef peak table and a summary table, plus the
# SHA-256 of every input so the numbers can be traced to the exact run files.
# Validation: none here. These are the runs' own audit results; neither
# configuration has passed the calibration gate (status
# "bounded_diagnostic_not_promoted" in both audit JSON files).
#
# Run from the repository root:
#   Rscript DataProcessing/cots_sdm_v2/final/report/external/metacommunity/build_metacommunity_summary.R

suppressPackageStartupMessages(library(ggplot2))

adria <- Sys.getenv("ADRIA_DIR", "C:/Users/smatthew/Documents/GitKraken/ADRIA.jl")
runs <- file.path(adria, "sandbox", "calibration", "runs")
out_dir <- "DataProcessing/cots_sdm_v2/final/report/external/metacommunity"
source("DataProcessing/cots_sdm_v2/final/R/figures.R")

configs <- data.frame(
  configuration = c("Working model", "Leading candidate"),
  run_id = c("20261012T_06_benthic_cyclone_1991_v2", "20261012T_15_food_grace_1991_v1"),
  treatment = c("k05_nocyc", "acrpoc_g0"),
  audit = c("benthic_cyclone_audit", "food_grace_audit"),
  stringsAsFactors = FALSE
)
reefs <- c("Lizard Island Reef" = "Lizard Isles", "MacGillivray Reef" = "Macgillivray Reef",
           "North Direction Reef" = "North Direction Island", "Eyrie Reef" = "Eyrie Reef")
obs_file <- file.path(adria, "sandbox", "data", "reef_cots.csv")

inputs <- c(obs_file, unlist(lapply(seq_len(nrow(configs)), function(i)
  file.path(runs, configs$run_id[i], c("trajectories.csv", paste0(configs$audit[i], c(".csv", ".json")))))))
missing <- inputs[!file.exists(inputs)]
if (length(missing)) stop("Missing inputs: ", paste(missing, collapse = ", "))

# Trajectories and per-reef audit rows for each configuration.
traj <- do.call(rbind, lapply(seq_len(nrow(configs)), function(i) {
  t <- read.csv(file.path(runs, configs$run_id[i], "trajectories.csv"), stringsAsFactors = FALSE)
  t <- t[t$treatment == configs$treatment[i] & t$year >= 1992, c("reef_name", "year", "simulated_cpue")]
  t$configuration <- configs$configuration[i]
  t
}))
keep <- c("reef_name", "first_peak_year", "second_peak_year", "first_peak_cpue", "second_peak_cpue",
          "max_cpue_2010_2017", "observed_max_cpue_2010_2017", "matched_peaks")
peaks <- do.call(rbind, lapply(seq_len(nrow(configs)), function(i) {
  a <- read.csv(file.path(runs, configs$run_id[i], paste0(configs$audit[i], ".csv")), stringsAsFactors = FALSE)
  a <- a[a$treatment == configs$treatment[i], ]
  gap <- if ("mean_gap_adults_ha" %in% names(a)) a$mean_gap_adults_ha else a$mean_gap_adults_2004_2010
  data.frame(configuration = configs$configuration[i], run_id = configs$run_id[i],
             treatment = configs$treatment[i], a[, keep], mean_gap_adults_ha_2004_2010 = gap)
}))
if (nrow(peaks) != 8 || !all(reefs %in% reefs[peaks$reef_name])) stop("Unexpected audit rows")

summary <- do.call(rbind, lapply(split(peaks, peaks$configuration), function(p) {
  j <- jsonlite::fromJSON(file.path(runs, p$run_id[1], paste0(configs$audit[configs$run_id == p$run_id[1]], ".json")))
  data.frame(configuration = p$configuration[1], run_id = p$run_id[1], treatment = p$treatment[1],
             status = j$status, passes_gate = j$treatments[[p$treatment[1]]]$passes_qualitative_gate,
             matched_peaks_of_seven = sum(p$matched_peaks),
             reefs_with_two_peaks = sum(!is.na(p$second_peak_year)),
             mean_first_peak_year = mean(p$first_peak_year),
             lizard_first_peak_cpue = p$first_peak_cpue[p$reef_name == "Lizard Island Reef"],
             mean_gap_adults_ha_2004_2010 = mean(p$mean_gap_adults_ha_2004_2010))
}))

obs <- read.csv(obs_file, stringsAsFactors = FALSE)
obs <- obs[obs$reef_name %in% reefs & obs$year >= 1992 & obs$year <= 2024, c("reef_name", "year", "cotsptow")]
obs$reef_name <- names(reefs)[match(obs$reef_name, reefs)]

# Observed maxima in the two outbreak windows (first wave 1992-2005; second
# wave 2010-2017, the audits' SECOND_WINDOW), as reference values for the table.
wave_max <- function(d, lo, hi) {
  w <- d[d$year >= lo & d$year <= hi, ]
  w[which.max(w$cotsptow), c("year", "cotsptow")]
}
observed <- do.call(rbind, lapply(split(obs, obs$reef_name), function(d) {
  a <- wave_max(d, 1992, 2005); b <- wave_max(d, 2010, 2017)
  data.frame(reef_name = d$reef_name[1], surveys_1992_2024 = nrow(d),
             first_wave_max_year = a$year, first_wave_max_cpue = a$cotsptow,
             second_wave_max_year = b$year, second_wave_max_cpue = b$cotsptow)
}))

write.csv(peaks, file.path(out_dir, "metacommunity_peaks.csv"), row.names = FALSE)
write.csv(observed, file.path(out_dir, "metacommunity_observed.csv"), row.names = FALSE)
write.csv(summary, file.path(out_dir, "metacommunity_summary.csv"), row.names = FALSE)
writeLines(c("input,sha256", paste0(sub(paste0(adria, "/"), "ADRIA.jl/", inputs, fixed = TRUE), ",",
                                    vapply(inputs, function(f) digest::digest(f, algo = "sha256", file = TRUE), ""))),
           file.path(out_dir, "metacommunity_inputs.csv"))

facet_order <- names(reefs)
traj$reef_name <- factor(traj$reef_name, facet_order)
obs$reef_name <- factor(obs$reef_name, facet_order)
cols <- c("Working model" = layer_colours[["legacy_clean_year2026"]],
          "Leading candidate" = layer_colours[["legacy_rg_year2026"]])

p <- ggplot() +
  geom_hline(yintercept = 0.22, colour = ink$muted, linetype = "dotted", linewidth = 0.4) +
  geom_line(data = traj, aes(year, simulated_cpue, colour = configuration), linewidth = 0.7) +
  geom_point(data = obs, aes(year, cotsptow, shape = "LTMP manta survey"), colour = ink$primary, size = 1.4) +
  facet_wrap(~reef_name, ncol = 2, scales = "free_y") +
  scale_colour_manual(values = cols, name = NULL) +
  scale_shape_manual(values = c("LTMP manta survey" = 16), name = NULL) +
  labs(x = "Year", y = "COTS per manta tow",
       title = "Metacommunity model: simulated and observed COTS outbreaks, Lizard Island cluster",
       subtitle = "Simulated adults converted at 0.015 COTS per tow per COTS/ha; dotted line: outbreak level (0.22 COTS per tow)",
       caption = "Neither configuration has passed the calibration gate. Runs: ADRIA.jl sandbox/calibration (branch cots-submodel).") +
  theme_report() +
  theme(legend.position = "bottom")
ggsave(file.path(out_dir, "fig_metacommunity_calibration.png"), p, width = 7.5, height = 5.2, dpi = 300, bg = ink$surface)
message("Wrote ", out_dir)
