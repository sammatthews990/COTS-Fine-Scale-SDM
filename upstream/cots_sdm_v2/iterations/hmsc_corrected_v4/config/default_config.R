default_config_v4 <- function(repo_root, v4_dir) {
  v2_dir <- normalizePath(file.path(v4_dir, "..", ".."), winslash = "/")
  v3_dir <- normalizePath(file.path(v4_dir, "..", "hmsc_geomorphic_v3"), winslash = "/")
  source(file.path(v2_dir, "config", "default_config.R"), local = TRUE)
  cfg <- default_config(repo_root, v2_dir)
  cfg$version <- "0.4.0"
  cfg$paths$v2_dir <- v2_dir
  cfg$paths$v3_dir <- v3_dir
  cfg$paths$v4_dir <- v4_dir
  cfg$paths$artifacts_dir <- file.path(v4_dir, "artifacts")
  cfg$paths$models_dir <- file.path(v4_dir, "models")
  cfg$paths$outputs_dir <- file.path(v4_dir, "outputs")
  cfg$paths$reports_dir <- file.path(v4_dir, "reports")
  cfg$paths$logs_dir <- file.path(v4_dir, "logs")
  cfg$paths$tmp_dir <- file.path(v4_dir, "tmp")
  cfg$paths$v2_artifacts_dir <- file.path(v2_dir, "artifacts")
  cfg$paths$v2_models_dir <- file.path(v2_dir, "models")
  cfg$paths$v2_outputs_dir <- file.path(v2_dir, "outputs")
  cfg$paths$v3_artifacts_dir <- file.path(v3_dir, "artifacts")
  cfg$paths$v3_models_dir <- file.path(v3_dir, "models")
  cfg$paths$v3_outputs_dir <- file.path(v3_dir, "outputs")
  cfg$validation$cell_support_m <- c(120L, 200L)
  cfg$validation$alignment_window_days <- 3L
  cfg$validation$spatial_block_m <- 3000L
  cfg$validation$folds <- 5L
  cfg$validation$expert_reefs <- c("Eyrie", "Batt", "Broadhurst", "Bait")
  cfg$manta$strong_scar_count <- 2L
  cfg$manta$strong_scar_codes <- c("c", "common")
  cfg$manta$weak_scar_codes <- c("p", "present")
  cfg$manta$negative_to_positive_ratio <- 3L
  cfg$model$forbidden_habitat_terms <- c(
    "year", "reef_key", "reef_name", "reef_label", "reef_year_id",
    "outbrProb", "reef_outbreak_probability", "legacy_probability",
    "legacy_rank", "effort", "method", "observation_support"
  )
  cfg$model$preliminary_samples <- 500L
  cfg$model$preliminary_thin <- 5L
  cfg$model$preliminary_transient <- 1000L
  cfg$model$preliminary_chains <- 4L
  cfg
}
