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
log_message("21 October snapshot imported; run 07 --calibration_report_only=1 next.")
