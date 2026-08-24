initialise_v2 <- function(script_file) {
  script_file <- normalizePath(script_file, winslash = "/", mustWork = TRUE)
  v2_dir <- normalizePath(
    file.path(dirname(script_file), ".."),
    winslash = "/",
    mustWork = TRUE
  )
  repo_root <- normalizePath(
    file.path(v2_dir, "..", ".."),
    winslash = "/",
    mustWork = TRUE
  )

  source(file.path(v2_dir, "config", "default_config.R"), local = FALSE)
  cfg <- default_config(repo_root, v2_dir)

  invisible(lapply(
    c(
      "io.R", "targets.R", "manta.R", "predictors.R", "resampling.R",
      "metrics.R", "models.R"
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
