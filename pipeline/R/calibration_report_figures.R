# Dated ADRIA report figures: frozen metadata transitively pins figure members
# and the raw observations/trajectories used to combine source Figures 1 and 4.
build_calibration_report_figures <- function(cfg) {
  cal <- checked_input(cfg, "calibration_report_20261021")
  proj <- checked_input(cfg, "projection_report_20261021")
  adria <- normalizePath(file.path(dirname(cal), "../../../../../.."), winslash = "/")
  status <- system2(Sys.getenv("PYTHON", "python"), shQuote(c(
    file.path(cfg$final_dir, "scripts", "build_calibration_report_figures.py"),
    "--calibration", cal, "--projections", proj, "--adria", adria,
    "--output", cfg$dirs$figures)))
  if (status != 0L) stop("Calibration report figure build failed")
}
