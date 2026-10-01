initialise_density_v6 <- function(script_file) {
  script_file <- sub("^--file=", "", script_file)
  v6_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
  iterations_dir <- normalizePath(file.path(v6_dir, ".."), winslash = "/")
  v2_dir <- normalizePath(file.path(iterations_dir, ".."), winslash = "/")
  v3_dir <- file.path(iterations_dir, "hmsc_geomorphic_v3")
  v4_dir <- file.path(iterations_dir, "hmsc_corrected_v4")
  repo_root <- (function(p) { while (!dir.exists(file.path(p, ".git"))) { if (dirname(p) == p) stop("Repository root not found"); p <- dirname(p) }; p })(v2_dir)

  source(file.path(v2_dir, "R", "bootstrap.R"))
  source(file.path(v3_dir, "R", "terrain_basis.R"))
  source(file.path(v4_dir, "config", "default_config.R"))
  source(file.path(v4_dir, "R", "helpers.R"))
  source(file.path(v6_dir, "R", "helpers.R"))

  cfg <- default_config_v4(repo_root, v4_dir)
  paths <- list(
    root = v6_dir,
    artifacts = file.path(v6_dir, "artifacts"),
    models = file.path(v6_dir, "models"),
    outputs = file.path(v6_dir, "outputs"),
    reports = file.path(v6_dir, "reports"),
    figures = file.path(v6_dir, "reports", "figures"),
    maps = file.path(v6_dir, "outputs", "four_reef_maps"),
    src = file.path(v6_dir, "src"),
    lib = file.path(v6_dir, "lib"),
    v2 = v2_dir, v3 = v3_dir, v4 = v4_dir
  )
  invisible(lapply(paths[c("artifacts", "models", "outputs", "reports", "figures", "maps", "src")],
    dir.create, recursive = TRUE, showWarnings = FALSE))
  .libPaths(c(paths$lib, .libPaths()))
  list(cfg = cfg, paths = paths)
}
