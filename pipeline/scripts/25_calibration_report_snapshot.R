# Question: what do the completed 21 October calibration and projection runs show?
# Method: import the dated report snapshot; preserve the October 15 results.
# Validation: verify pinned metadata and all member figures/tables; no refitting.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "..", "R", "common.R"))
cfg <- initialise_final()
cal_path <- checked_input(cfg, "calibration_report_20261021")
proj_path <- checked_input(cfg, "projection_report_20261021")
cal <- jsonlite::fromJSON(cal_path)
proj <- jsonlite::fromJSON(proj_path)
for (p in c(cal_path, proj_path)) {
  meta <- jsonlite::fromJSON(p)
  for (n in names(meta$outputs)) verify_frozen_input(file.path(dirname(p), n), meta$outputs[[n]])
}
write_table(read.csv(file.path(dirname(cal_path), "peak_timing_summary.csv")), cfg,
            "calibration_20261021_peaks.csv")
write_table(read.csv(file.path(dirname(proj_path), "culling_effects_2027_2050.csv")), cfg,
            "projection_20261021_reef_outcomes.csv")
write_table(data.frame(metric = c("base_seeds", "calibration_arms_per_cull_state", "projection_arms", "climate_models"),
                       value = c(length(cal$base_arms), cal$n_arms, proj$arms_used, length(proj$runs))),
            cfg, "calibration_20261021_summary.csv")
# Reefs outside the four-reef calibration score. Martin and Linnet test wave
# timing; the wave-exposure settlement term was designed after false outbreaks
# on Yonge, Carter and No Name, so those three only check the design. Baseline:
# the mean observed series of the four calibration reefs, with no model.
adria <- normalizePath(file.path(dirname(cal_path), "../../../../../.."), winslash = "/")
extra_path <- file.path(adria, gsub("\\\\", "/", cal$run), "trajectories_extra.csv")
obs_path <- file.path(adria, "sandbox/data/reef_cots.csv")
verify_frozen_input(extra_path, cal$inputs[["trajectories_extra.csv"]])
verify_frozen_input(obs_path, cal$inputs[["reef_cots.csv"]])
q_manta <- 0.015
irruption_cpue <- 0.22
extra <- read.csv(extra_path)
extra <- extra[extra$treatment %in% cal$base_arms & extra$year >= 1992 & extra$year <= 2026, ]
model <- aggregate(adults_ha ~ reef_name + year, extra, mean)
model$cpue <- q_manta * model$adults_ha
obs <- aggregate(cotsptow ~ reef_name + year, read.csv(obs_path), mean)
calibration_mean <- aggregate(cotsptow ~ year, obs[obs$reef_name %in% c("Lizard Isles", "Macgillivray Reef", "North Direction Island", "Eyrie Reef"), ], mean)
tss <- function(pred, truth) {
  if (!any(truth) || all(truth)) return(NA_real_)
  mean(pred[truth]) - mean(pred[!truth])
}
fit_row <- function(reef, group) {
  m <- model[model$reef_name == reef, ]
  o <- obs[obs$reef_name == reef & obs$year >= 1992 & obs$year <= 2026, ]
  o$model <- m$cpue[match(o$year, m$year)]
  o$baseline <- calibration_mean$cotsptow[match(o$year, calibration_mean$year)]
  o <- o[is.finite(o$model) & is.finite(o$baseline), ]
  rho <- function(x) if (stats::sd(o$cotsptow) > 0 && stats::sd(x) > 0) stats::cor(x, o$cotsptow, method = "spearman") else NA_real_
  data.frame(reef = reef, group = group, surveys = nrow(o), obs_max_cpue = max(o$cotsptow), model_max_cpue = max(m$cpue),
             obs_years_irruption = sum(o$cotsptow >= irruption_cpue), model_years_irruption = sum(o$model >= irruption_cpue),
             mae = mean(abs(o$model - o$cotsptow)), baseline_mae = mean(abs(o$baseline - o$cotsptow)),
             spearman = rho(o$model), baseline_spearman = rho(o$baseline),
             tss = tss(o$model >= irruption_cpue, o$cotsptow >= irruption_cpue),
             baseline_tss = tss(o$baseline >= irruption_cpue, o$cotsptow >= irruption_cpue))
}
write_table(rbind(fit_row("Martin Reef", "irruption history"), fit_row("Linnet Reef", "irruption history"),
                  fit_row("Yonge Reef", "almost no COTS"), fit_row("Carter Reef", "almost no COTS"),
                  fit_row("No Name Reef", "almost no COTS")),
            cfg, "calibration_20261021_additional_reefs.csv")
log_message("21 October snapshot imported; run 07 --calibration_report_only=1 next.")
