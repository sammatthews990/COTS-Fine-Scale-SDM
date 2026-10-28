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

  adria_dir <- Sys.getenv("ADRIA_DIR", file.path(dirname(repo_root), "ADRIA.jl"))
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
      ranking_deployed_habitat = list(
        path = file.path(outputs, "layers/cots_legacy_rg_year2026_clean_fallback_30m.tif"),
        sha256 = "876a6b7d1dbc7282651498e64a241063c5ebfdc55ca1214ce239f63aed0343b2"
      ),
      ranking_train_site_years = list(
        path = file.path(exploration_dir, "artifacts/site_year_all.rds"),
        sha256 = "05a3e784c0b486b5ac96f18f0245f4109a6d824b6632eb6e2f611518cf173b9c"
      ),
      ranking_site_sample_points = list(
        path = file.path(final_dir, "cache/cull_sites_points.rds"),
        sha256 = "c2b803450573de551ef1823e85f326ec7f24c19d7a300e907516c504c6ac2218"
      ),
      # Review 09: geometry of the domain used by the 21 October population
      # model; the calibration store now holds the archived ADRIA sandbox data.
      report_metapop_domain = list(
        path = file.path(dirname(repo_root), "COTSMod-Calibration/store/data/Lizard_Historical_v0.4/spatial/lizard_cluster.gpkg"),
        sha256 = "f2ddd55b0ccefebf7aa8003623a20977125cb9ea3a8f7df745c35f0c488098f9"
      ),
      # Review 06: completed concentration analysis, 23 October 2026.
      # Freeze attribute data as well as geometry; concentrations are in DBF.
      larval_report_summary = list(
        path = file.path(dp, "connectivity_lizard/outputs/tables/reef_concentration_summary.csv"),
        sha256 = "4d7472bd822407388a2a57cb0ccbc915edb6283ecd1bde7211405587f8cef2e2"
      ),
      larval_report_zones = list(
        path = file.path(dp, "connectivity_lizard/outputs/tables/reef_zones.csv"),
        sha256 = "48bed3b51d26636c76c3670fcc5aa609218948a03a42c300974f411b523f4d55"
      ),
      larval_report_releases = list(
        path = file.path(dp, "connectivity_lizard/outputs/tables/zone_by_release.csv"),
        sha256 = "20e62dc04e7b022f8e8852bd79332aa69f7c08c78cecafc0a0b4ae6518cdf650"
      ),
      larval_report_latitude = list(
        path = file.path(dp, "connectivity_lizard/outputs/tables/latitude_trend.csv"),
        sha256 = "8e15a3b4723c8248c6f7957c3ec2c6932e2e238701674ec29b9a4005603124f5"
      ),
      larval_report_day10_shp = list(
        path = file.path(dp, "data/fngbr_LarvalConcentrations/fngbr_reef_daily_means_day10.shp"),
        sha256 = "99a8f5a81849fac2b4fac083e6a5af639af25b9bcc1a17233709153040ae0cd9"
      ),
      larval_report_day10_dbf = list(
        path = file.path(dp, "data/fngbr_LarvalConcentrations/fngbr_reef_daily_means_day10.dbf"),
        sha256 = "217e1e637c8eb4a62a9d8a2602dda1ddba16e661fabaad16ccb1747d61cdf5a9"
      ),
      larval_report_day10_shx = list(
        path = file.path(dp, "data/fngbr_LarvalConcentrations/fngbr_reef_daily_means_day10.shx"),
        sha256 = "b5a16f068ce92466d272acd3a322100064897ea24b887317e875720e72d187c9"
      ),
      larval_report_day10_prj = list(
        path = file.path(dp, "data/fngbr_LarvalConcentrations/fngbr_reef_daily_means_day10.prj"),
        sha256 = "a02a27b1d1982c8516d83398e85a3c8b1aef1713c13ef4d84d7bde17430c07c4"
      ),
      # 21 October report snapshot. Each pinned metadata file also freezes its
      # constituent figures/data by SHA-256; script 25 and the figure builder
      # verify those member checksums before using them. Keep October 15 intact.
      calibration_report_20261021 = list(
        path = file.path(adria_dir, "sandbox/data/derived/calibration_report_figures/2026-10-21_wwave_lf04/calibration/metadata.json"),
        sha256 = "b7f8c8dbb82a5b5f77456932182ea6c4f39e30e55f63199a7b7fefb615739917"
      ),
      projection_report_20261021 = list(
        path = file.path(adria_dir, "sandbox/data/derived/calibration_report_figures/2026-10-21_wwave_lf04/projections_2050/metadata.json"),
        sha256 = "ca67abdf433b3cbe998bda0daad4fb2ed47cfb5d75f19b0042ef34e19dd7a4a6"
      ),
      # Review 03: frozen source records and existing ranking snapshots for a
      # threshold-only reassessment. Refresh these pins explicitly after reruns.
      assessment_dives = list(
        path = file.path(exploration_dir, "artifacts/cull_observations_clean.rds"),
        sha256 = "995a566f09559b507cf1f2dd7e0d95a034692e48a87ff9403b4bb1332e718767"
      ),
      assessment_sites = list(
        path = file.path(final_dir, "cache/cull_sites_designated.gpkg"),
        sha256 = "14b0b292bfc163d425d90e87cde39aef4a1256ba16c165802c2cf426880e9176"
      ),
      assessment_evidence = list(
        path = file.path(final_dir, "cache/cull_site_evidence_scores.rds"),
        sha256 = "9774b06af7ab11674bec49d13dc7d0ea70054eb5fb173162fa2c10b97ebe740e"
      ),
      assessment_manta_sites = list(
        path = file.path(final_dir, "cache/manta_response_site_ranks.rds"),
        sha256 = "fe68b13539404cd37c9cd119f1612cd2391bec338014a1844bcd7034de68c2ac"
      ),
      assessment_forward = list(
        path = file.path(final_dir, "cache/site_ranking_forward.rds"),
        sha256 = "2fabc8d6e888eef58e4aef3b94f19bc0119f9a9c7ac3443ecbf8bf6a94056ba5"
      ),
      # Explicit frozen coefficients for replaying all candidate-site ranks,
      # including never-dived sites, without refitting the forward models.
      assessment_forward_coefficients = list(
        path = file.path(final_dir, "outputs/tables/site_ranking_forward_coefficients.csv"),
        sha256 = "a82524e8309b55f14a4acbc2ba303be2fe2e2bd74e47354a509ec9a77ce6924b"
      ),
      # Metacommunity reporting snapshot, 15 October 2026. Archived run files;
      # these pins do not select or promote a biological configuration.
      mc_candidates_candidates_audit = list(
        path = file.path(adria_dir, "sandbox/calibration/runs/20261014T_48_candidates_1991_v2/candidates_audit.csv"),
        sha256 = "fba377b700c63f2526f50e1a886086455b3337cf61ee1b08a75bafc6394b7cf4"
      ),
      mc_candidates_candidates_audit_json = list(
        path = file.path(adria_dir, "sandbox/calibration/runs/20261014T_48_candidates_1991_v2/candidates_audit.json"),
        sha256 = "73fce61126f3b459be215b758100b9503ab6da10817130302ba60d3b9abdb52e"
      ),
      mc_candidates_candidates_replicates = list(
        path = file.path(adria_dir, "sandbox/calibration/runs/20261014T_48_candidates_1991_v2/candidates_replicates.csv"),
        sha256 = "06261c108c8e5ada5ee535e1f429b40393675c245651d5e646b1f43baae1153b"
      ),
      mc_candidates_trajectories = list(
        path = file.path(adria_dir, "sandbox/calibration/runs/20261014T_48_candidates_1991_v2/trajectories.csv"),
        sha256 = "a65e73532586c1142fcd5b643ee977907cc8f77a976127b2f8a7616068cd000b"
      ),
      mc_culling_culling_audit = list(
        path = file.path(adria_dir, "sandbox/calibration/runs/20261014T_42_culling_1991_v1/culling_audit.csv"),
        sha256 = "991320244274b8092391c61e577bf7ffdbfef3e4f59c2cd47ca920a2a9bc59c4"
      ),
      mc_culling_culling_removals = list(
        path = file.path(adria_dir, "sandbox/calibration/runs/20261014T_42_culling_1991_v1/culling_removals.csv"),
        sha256 = "dba0f86f1decfc623b4d3d00bd50547c0f055b447c34e08fc26b90034a025496"
      ),
      mc_original_trajectories = list(
        path = file.path(adria_dir, "sandbox/calibration/runs/20261015T_08_chl_survival_owen_global_2018_2023_1991_v1/trajectories.csv"),
        sha256 = "b0f8c7f334cf5c6514fae8d385b9d64c21f1c70269f70b1cfd2969437c1d798c"
      ),
      mc_original_metadata = list(
        path = file.path(adria_dir, "sandbox/calibration/runs/20261015T_08_chl_survival_owen_global_2018_2023_1991_v1/metadata.toml"),
        sha256 = "ce849387524bd2fab62f27e7a584929f5410065f722e645491a08610cdd15113"
      ),
      mc_corrected_trajectories = list(
        path = file.path(adria_dir, "sandbox/calibration/runs/20261015T_09_chl_survival_owen_global_2018_2023_corrected_1991_v1/trajectories.csv"),
        sha256 = "c00edaaa84a7354c43300a50a96249e0fc1b0574ea36d36e4cc3cd14f67b4883"
      ),
      mc_corrected_metadata = list(
        path = file.path(adria_dir, "sandbox/calibration/runs/20261015T_09_chl_survival_owen_global_2018_2023_corrected_1991_v1/metadata.toml"),
        sha256 = "9e09ba1a6d07bd6dd6f3084aab3d22d8d44f3a19e1a13600fd5eeb119c16e56f"
      ),
      mc_original_score_detectability_scores = list(
        path = file.path(adria_dir, "sandbox/calibration/runs/20261015T_10_rescore_chl_survival_original/detectability_scores.csv"),
        sha256 = "d39b61af63846a3dca2f5c2bc3abdfdd97f474222396f72b213f5a2c4c504892"
      ),
      mc_corrected_score_detectability_scores = list(
        path = file.path(adria_dir, "sandbox/calibration/runs/20261015T_12_rescore_chl_survival_corrected/detectability_scores.csv"),
        sha256 = "8f46320e10f5e6b0692d7c8f7a25a8bddc603f0b9c584b1c36bd7b7bd0e7ae71"
      ),
      mc_original_gate_gate_g2prime = list(
        path = file.path(adria_dir, "sandbox/calibration/runs/20261015T_11_gate_chl_survival_original/gate_g2prime.csv"),
        sha256 = "b3b4a44e29941aa35871019a8856f6bcf132a41baf85be2295e4d2511c517879"
      ),
      mc_corrected_gate_gate_g2prime = list(
        path = file.path(adria_dir, "sandbox/calibration/runs/20261015T_13_gate_chl_survival_corrected/gate_g2prime.csv"),
        sha256 = "5cd7dc7b207ef83bcaf36e1760ad7a12f64c55140d364d5ca204c1748772bbad"
      ),
      mc_observations = list(
        path = file.path(adria_dir, "sandbox/data/reef_cots.csv"),
        sha256 = "0d8871dd90a5e0f6f883fbd22a4c92d0f036d85f7e2dd484197ee7f8ffeeaef5"
      ),
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
      # Manta sheet of the CSIRO workbook. The footprint tows above are its rows
      # (visit_id = file||row); read only for the tow coral cover codes (24).
      manta_workbook = list(
        path = file.path(dp, "data", "260529-COTS-Manta-Cull-RHIS-Lawrence-CSIRO.xlsx"),
        sha256 = "950f5296f6bf206bb95e16fe0e5fa391bfa8f6cf4065988ecc9d022164109b04"
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
    ranking_candidates = list(
      half_life_days = 90,
      habitat_quantile = 0.90,
      bootstrap_reps = 500L,
      evaluation_years = 2019:2026,
      recent_years = 2024:2026
    ),
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
