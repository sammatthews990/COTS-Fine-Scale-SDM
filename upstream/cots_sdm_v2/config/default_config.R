default_config <- function(repo_root, v2_dir) {
  list(
    version = "0.2.0",
    paths = list(
      repo_root = repo_root,
      v2_dir = v2_dir,
      cull_workbook = file.path(
        repo_root, "DataProcessing", "data",
        "260529-COTS-Manta-Cull-RHIS-Lawrence-CSIRO.xlsx"
      ),
      predictor_stack = file.path(
        repo_root, "DataProcessing", "data",
        "predictors_clean_12layer.tif"
      ),
      predictor_stack_rg = file.path(
        repo_root, "DataProcessing", "data",
        "predictors_terra_30m_full_extended.tif"
      ),
      geomorphic_raster = file.path(
        repo_root, "DataProcessing", "data",
        "GBR10_GBRMP_Geomorphic", "GBR10 GBRMP Geomorphic.tif"
      ),
      outbreak_hindcast_file = file.path(
        repo_root, "DataProcessing", "data", "gbrPredsAdj_20262408.csv"
      ),
      outbreak_forecast_file = file.path(
        repo_root, "DataProcessing", "data", "gbrPreds2026Adj_20262408.csv"
      ),
      manta_history_csv = file.path(
        repo_root, "DataProcessing", "data",
        "COTS Program  Manta Tow Data-2026-02-04.csv"
      ),
      manta_current_workbook = file.path(
        repo_root, "DataProcessing", "data",
        "260529-COTS-Manta-Cull-RHIS-Lawrence-CSIRO.xlsx"
      ),
      salad_workbook = file.path(
        repo_root, "DataProcessing", "data",
        "COTS SALAD MASTER 2021-2026.xlsx"
      ),
      legacy_2026_raster = file.path(
        repo_root, "DataProcessing", "outputs",
        "COTS_prob_0.02_cpue_year2026_clean_20260820_1041.tif"
      ),
      legacy_agnostic_raster = file.path(
        repo_root, "DataProcessing", "outputs",
        "COTS_prob_0.02_cpue_agnostic_clean_20260821_1046.tif"
      ),
      artifacts_dir = file.path(v2_dir, "artifacts"),
      models_dir = file.path(v2_dir, "models"),
      outputs_dir = file.path(v2_dir, "outputs"),
      logs_dir = file.path(v2_dir, "logs")
    ),
    response = list(
      cpue_threshold = 0.02,
      default_target = "site_year_first_visit"
    ),
    validation = list(
      spatial_folds = 5L,
      spatial_seed = 20260824L,
      forward_test_year = 2026L,
      rolling_origin_years = 2022:2026
    ),
    model = list(
      seed = 20260824L,
      nrounds = 400L,
      eta = 0.03,
      max_depth = 4L,
      min_child_weight = 20,
      subsample = 0.8,
      colsample_bytree = 0.8,
      lambda = 1,
      alpha = 0,
      nthread = 4L
    ),
    outbreak = list(
      feature_policy = "lagged_only",
      required_lags = 1:3,
      forecast_year = 2026L,
      hindcast_provenance = "adjusted_hindcast_not_yet_audited",
      forecast_provenance = "adjusted_2026_prediction_not_yet_audited"
    ),
    manta = list(
      minimum_tows = 10L,
      minimum_survey_years = 3L,
      minimum_tow_distance_m = 100,
      require_absent_scar_code = TRUE,
      integration_role = "presence_stage_negative_evidence"
    ),
    site_habitat = list(
      manta_site_grid_m = 200,
      manta_minimum_tows = 5L,
      manta_minimum_years = 2L,
      manta_local_conflict_distance_m = 250,
      manta_maximum_negative_weight_ratio = 1.0,
      manta_weight_sensitivity_ratios = c(0.10, 0.25, 0.50),
      within_reef_folds = 5L
    ),
    visual = list(
      probability_min = 0.20,
      probability_max = 0.80,
      difference_limit = 0.30,
      raster_aggregation_factor = 5L,
      forecast_voronoi_buffer_degrees = 1.0
    )
  )
}
