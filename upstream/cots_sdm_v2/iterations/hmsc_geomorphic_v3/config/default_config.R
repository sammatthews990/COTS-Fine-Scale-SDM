default_config_v3 <- function(repo_root, v3_dir) {
  v2_dir <- normalizePath(file.path(v3_dir, "..", ".."), winslash = "/")
  source(file.path(v2_dir, "config", "default_config.R"), local = TRUE)
  cfg <- default_config(repo_root, v2_dir)
  cfg$version <- "0.3.0"
  cfg$paths$v3_dir <- v3_dir
  cfg$paths$artifacts_dir <- file.path(v3_dir, "artifacts")
  cfg$paths$models_dir <- file.path(v3_dir, "models")
  cfg$paths$outputs_dir <- file.path(v3_dir, "outputs")
  cfg$paths$reports_dir <- file.path(v3_dir, "reports")
  cfg$paths$logs_dir <- file.path(v3_dir, "logs")
  cfg$paths$tmp_dir <- file.path(v3_dir, "tmp")
  cfg$paths$v2_artifacts_dir <- file.path(v2_dir, "artifacts")
  cfg$paths$v2_models_dir <- file.path(v2_dir, "models")
  cfg$paths$v2_outputs_dir <- file.path(v2_dir, "outputs")
  cfg$validation$alignment_windows_days <- c(1L, 3L, 7L)
  cfg$validation$cell_support_m <- c(60L, 120L)
  cfg$validation$expert_reefs <- c("Eyrie", "Batt", "Broadhurst", "Bait")
  cfg$model$terrain_basis <- "within_family_pca"
  cfg$model$terrain_pca_variance <- 0.80
  cfg$model$terrain_pca_max_components <- 2L
  cfg$model$forbidden_habitat_terms <- c(
    "year", "reef_key", "reef_name", "reef_label", "reef_year_id",
    "outbrProb", "reef_outbreak_probability", "legacy_probability",
    "legacy_rank", "effort", "method"
  )
  cfg
}

