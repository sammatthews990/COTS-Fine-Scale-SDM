initialise_v4 <- function(script_file) {
  script_file <- sub("^--file=", "", script_file)
  v4_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
  v2_dir <- normalizePath(file.path(v4_dir, "..", ".."), winslash = "/")
  v3_dir <- normalizePath(file.path(v4_dir, "..", "hmsc_geomorphic_v3"), winslash = "/")

  source(file.path(v2_dir, "R", "bootstrap.R"))
  source(file.path(v3_dir, "R", "terrain_basis.R"))
  source(file.path(v3_dir, "R", "footprints.R"))
  source(file.path(v4_dir, "config", "default_config.R"))
  source(file.path(v4_dir, "R", "helpers.R"))

  repo_root <- normalizePath(file.path(v2_dir, "..", ".."), winslash = "/")
  cfg <- default_config_v4(repo_root, v4_dir)
  invisible(lapply(c(
    cfg$paths$artifacts_dir, cfg$paths$models_dir, cfg$paths$outputs_dir,
    cfg$paths$reports_dir, cfg$paths$logs_dir, cfg$paths$tmp_dir
  ), dir.create, recursive = TRUE, showWarnings = FALSE))
  list(config = cfg, v2_dir = v2_dir, v3_dir = v3_dir, v4_dir = v4_dir)
}
