initialise_v3 <- function(script_file) {
  script_file <- normalizePath(script_file, winslash = "/", mustWork = TRUE)
  v3_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
  v2_dir <- normalizePath(file.path(v3_dir, "..", ".."), winslash = "/")
  repo_root <- (function(p) { while (!dir.exists(file.path(p, ".git"))) { if (dirname(p) == p) stop("Repository root not found"); p <- dirname(p) }; p })(v2_dir)

  source(file.path(v2_dir, "R", "bootstrap.R"), local = FALSE)
  source(file.path(v3_dir, "config", "default_config.R"), local = FALSE)
  parent <- initialise_v2(file.path(v2_dir, "scripts", "35_prepare_hmsc_v2_visits.R"))
  cfg <- default_config_v3(repo_root, v3_dir)
  invisible(lapply(
    c("footprints.R", "terrain_basis.R", "occupancy.R"),
    function(f) source(file.path(v3_dir, "R", f), local = FALSE)
  ))
  invisible(lapply(
    unlist(cfg$paths[c("artifacts_dir", "models_dir", "outputs_dir",
                       "reports_dir", "logs_dir", "tmp_dir")]),
    dir.create, recursive = TRUE, showWarnings = FALSE
  ))
  list(config = cfg, repo_root = repo_root, v2_dir = v2_dir, v3_dir = v3_dir)
}

