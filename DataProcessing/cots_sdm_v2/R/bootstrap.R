initialise_v2 <- function(script_file) {
  script_file <- normalizePath(script_file, winslash = "/", mustWork = TRUE)
  v2_dir <- normalizePath(
    file.path(dirname(script_file), ".."),
    winslash = "/",
    mustWork = TRUE
  )
  # Walk up to the repository root so this tree also runs from archive/.
  repo_root <- (function(p) { while (!dir.exists(file.path(p, ".git"))) { if (dirname(p) == p) stop("Repository root not found"); p <- dirname(p) }; p })(v2_dir)

  source(file.path(v2_dir, "config", "default_config.R"), local = FALSE)
  cfg <- default_config(repo_root, v2_dir)

  invisible(lapply(
    c(
      "io.R", "targets.R", "manta.R", "predictors.R", "resampling.R",
      "metrics.R", "models.R", "calibration.R", "mapping.R",
      "site_habitat.R", "multiscale.R", "presence_only.R",
      "two_part.R", "joint_visits.R", "gaussian_terrain.R"
    ),
    function(f) source(file.path(v2_dir, "R", f), local = FALSE)
  ))

  ensure_output_directories(cfg)
  list(config = cfg, repo_root = repo_root, v2_dir = v2_dir)
}

current_script_file <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) != 1L) {
    stop("Run this entry point with Rscript so its path can be resolved.")
  }
  sub("^--file=", "", file_arg)
}
