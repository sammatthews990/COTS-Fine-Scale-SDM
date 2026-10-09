# Configuration for the final COTS fine-scale SDM pipeline.
#
# Every input below is frozen: it was produced by the exploratory work in
# cots_sdm_v2 (see README.md, "Provenance") and is pinned by SHA-256 so that a
# silently regenerated upstream file stops the pipeline rather than changing a
# delivered layer.

final_config <- function(final_dir) {
  repo_root <- find_repo_root(final_dir)
  dp <- file.path(repo_root, "DataProcessing")
  v2_dir <- file.path(dp, "cots_sdm_v2")
  # The exploratory tree lives either at cots_sdm_v2/ or, once archived, at
  # cots_sdm_v2/archive/. Both layouts are supported.
  exploration_dir <- if (dir.exists(file.path(v2_dir, "archive", "iterations"))) {
    file.path(v2_dir, "archive")
  } else {
    v2_dir
  }
  v4_artifacts <- file.path(exploration_dir, "iterations", "hmsc_corrected_v4", "artifacts")

  outputs <- file.path(final_dir, "outputs")
  list(
    repo_root = repo_root,
    final_dir = final_dir,
    exploration_dir = exploration_dir,
    dirs = list(
      outputs = outputs,
      layers = file.path(outputs, "layers"),
      tables = file.path(outputs, "tables"),
      figures = file.path(outputs, "figures"),
      models = file.path(final_dir, "models"),
      cache = file.path(final_dir, "cache"),
      logs = file.path(final_dir, "logs")
    ),

    inputs = list(
      predictor_stack_clean = list(
        path = file.path(dp, "data", "predictors_clean_12layer.tif"),
        sha256 = "dd9e201242f3c367357fcfe911159cf1baeff20782f12924289d687691ba2a8e"
      ),
      predictor_stack_rg = list(
        path = file.path(dp, "data", "predictors_terra_30m_full_extended.tif"),
        sha256 = "3bc1c1b4d4939f45a3043c2a083e1e9a908a9e9e563d8ab3452073a5e41085bb"
      ),
      geomorphic = list(
        path = file.path(dp, "data", "GBR10_GBRMP_Geomorphic", "GBR10 GBRMP Geomorphic.tif"),
        sha256 = "291627f0c15d4ece1910149abd915b0926bd8fdf4d436ba2f0b91435f1f88646"
      ),
      reef_features = list(
        path = file.path(dp, "data", "Great_Barrier_Reef_Features.gpkg"),
        sha256 = "1bc36b741064e106fd46a2a51aeb2a490e3af122f4041ca9824afbae2ac891b4"
      ),
      legacy_clean_model = list(
        path = file.path(dp, "outputs", "fitted_model_clean_20260820_1041.rds"),
        sha256 = "2943dc9c21b15d5a2f7e2871ac5816b06634591d8e49f1d333079f7f58388d93"
      ),
      legacy_clean_raster = list(
        path = file.path(dp, "outputs", "COTS_prob_0.02_cpue_year2026_clean_20260820_1041.tif"),
        sha256 = "cd75e8997ffacaa2b540012464dcc26b1be9680aa07d337d3803a9fbece6ee9f"
      ),
      legacy_clean_metrics = list(
        path = file.path(dp, "outputs", "model_performance_metrics_clean_20260820_1041.rds"),
        sha256 = "547ef4f1f800d7a487d968499d5092b97bc972268b53b6d1b1753dd5973d5aad"
      ),
      legacy_rg_model = list(
        path = file.path(dp, "outputs", "fitted_model_withRG_20260820_1217.rds"),
        sha256 = "a88177cced079c7b5c1d0b8e14dd71841d7ffb0520eb2ea9cb20f5a26c1164e1"
      ),
      legacy_rg_metrics = list(
        path = file.path(dp, "outputs", "model_performance_metrics_withRG_20260820_1217.rds"),
        sha256 = "4f517ee99a5a6b36896d5152feb15c0f753844b1a37f00c84fbad8d3547da9c1"
      ),
      hmsc_enhanced_joint = list(
        path = file.path(v4_artifacts, "enhanced_joint_fit_data_120m.rds"),
        sha256 = "168078eb35d92d481c00540b16ecf0a7f925433b4b6c7c24c8b0deed32e31ec4"
      ),
      hmsc_clean_joint = list(
        path = file.path(v4_artifacts, "clean_joint_fit_data_120m.rds"),
        sha256 = "ad70a28fbc26dfe1cc54124155dc2a6617430fb0e3730555c95543de9a83b876"
      ),
      hmsc_clean_native_joint = list(
        path = file.path(v4_artifacts, "clean_native_joint_fit_data_120m.rds"),
        sha256 = "8e4527272b5e93cd3519e9801bf372bccee789981b0e904177c857026dc6e11e"
      ),
      hmsc_reference_cells = list(
        path = file.path(v4_artifacts, "corrected_cell_units_120m.rds"),
        sha256 = "c66ce09757ef52c52660d5a72e5b408af0d1a49f140c223a058a6422dbb388e4"
      ),
      # Manta tows (evidence class, active reef-year) and their 120 m footprint
      # links; used only for the independent cross-layer check in the figures.
      manta_tows = list(
        path = file.path(v4_artifacts, "manta_footprints_120m.rds"),
        sha256 = "34079b26740473ec3a901547980b06e858d8da8417cd478e002e09d4e4e14154"
      ),
      # Designated cull sites (one polygon per site, one KML layer per reef),
      # including sites that have never been dived.
      cull_sites_kmz = list(
        path = file.path(repo_root, "data", "Eotr_CotsCullSites_2025_11_19_1_58_PM_1.kmz"),
        sha256 = "9ff89578a76a2220009c5250b8965a1eac4640be197f01c7d05d051790d8c3b2"
      ),
      footprint_links = list(
        path = file.path(exploration_dir, "iterations", "hmsc_geomorphic_v3", "artifacts", "footprint_links_120m.rds"),
        sha256 = "8aad8b68a79b5b8d8ab67357037e04ae7cb008e68a3cea0c10be5b7d9a8c61dd"
      )
    ),

    legacy = list(
      predict_year = 2026L,
      cpue_threshold = 0.02,
      # Layer names the RG workflow was trained with (fit_predict_cots_sdm.R).
      rg_layer_names = c(
        "BPI", "EAST", "HCU", "NORTH", "SLO", "SVF", "VCU", "VRM", "DEM",
        "AhyaD", "Aspat", "Aten",
        "RG_bathy", "RG_slope", "RG_turbid", "RG_waves_Hs", "RG_waves_Tp"
      )
    ),

    hmsc = list(
      support_m = 120L,
      target_response = "cull_high_density",
      # Fallback priority: first component with a complete predictor set wins.
      components = c("enhanced_joint", "clean_joint", "clean_native_joint"),
      production = list(samples = 1000L, thin = 30L, transient = 5000L, chains = 4L),
      seed = 20261001L,
      gates = list(min_ess = 400, max_psrf_upper = 1.1),
      posterior_draws_for_maps = 400L
    ),

    # History kernel + habitat blend. Fixed before the forward validation was
    # run (09b); it beat the nested per-year selection overall and where
    # history is sparse, so it is the deployed setting.
    history = list(sigma_m = 500, tau_years = Inf, kappa = 2),

    # Operational efficiency scenarios (12). Cost rate as used in the earlier
    # validation report (SDM_Validation_Report.qmd); hours are reported too.
    operations = list(
      diver_hour_cost_aud = 950,
      keep_fractions = c(1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3),
      headline_keep = 0.6,
      random_reps = 200L,
      bootstrap_reps = 200L
    ),

    # Cull-site ranking for the next season (20). Chosen from the forward test
    # (19): cull-trained habitat rank + previous-season evidence; the manta
    # ensemble and "both ranks" variants performed the same, and training on
    # all earlier seasons matched the last 5 and beat the last 3.
    site_ranking = list(
      target_season = 2027L,
      habitat = "cull",               # cull | ensemble | both
      forward_window = "all",         # column in 19's results that matches
      training_years = c(2019L, 2026L),
      # Calendar 2026 holds only Jan-May, so "last season" for 2027 is the
      # most recent 12 months of dives and tows.
      evidence_end = as.Date("2026-05-31"),
      evidence_days = 365L,
      kernel_sigma_m = 300,
      # Priority classes by share of the reef's designated sites, from the top.
      priority_breaks = c(A = 0.2, B = 0.5, C = 0.8)
    ),

    figures = list(
      # Four-reef comparison panels (fig02, fig02b, fig10): GBR reef-feature
      # LABEL_IDs, north to south.
      four_reefs = c(
        eyrie = "14-118",
        batt = "16-029",
        broadhurst = "18-100a",
        lady_musgrave = "23-082a"
      ),
      # Optional display windows (MGA55 m). Lady Musgrave's reef polygon
      # includes a large western platform; the window shows the ring reef
      # around the island. Percentiles are still ranked over the whole polygon.
      four_reef_windows = list(
        lady_musgrave = c(xmin = 1048300, ymin = 7343800, xmax = 1054100, ymax = 7348700)
      )
    )
  )
}

input_path <- function(cfg, name) {
  entry <- cfg$inputs[[name]]
  if (is.null(entry)) stop("Unknown input: ", name)
  entry$path
}

# Verify a frozen input's checksum (skipped for the multi-GB stacks unless
# full = TRUE; 00_check_inputs.R verifies everything).
checked_input <- function(cfg, name, full = FALSE) {
  entry <- cfg$inputs[[name]]
  large <- file.exists(entry$path) && file.info(entry$path)$size > 1e9
  verify_frozen_input(entry$path, if (large && !full) NA_character_ else entry$sha256)
}
