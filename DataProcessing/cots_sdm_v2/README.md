# COTS fine-scale SDM v2

> **Start here: [`final/`](final/README.md).** The final pipeline delivers the
> three report layers (legacy clean + year, legacy RG + year, Hmsc joint) and
> the report figures from frozen, checksummed inputs. Everything else in this
> directory is the exploratory record that produced those inputs and is kept
> for provenance only.

This directory is an isolated refactor of the COTS modelling workflow. It does
not source, overwrite, or modify the legacy model scripts or their outputs.

The refactor separates three questions that were previously mixed together:

1. **Static habitat prior** — where can high-density COTS be supported?
2. **Dynamic outbreak risk** — which reefs are likely to be occupied or
   outbreaking in a target year?
3. **Conditional high density** — given COTS presence, where is CPUE likely to
   exceed 0.02?

The first implemented milestone is a controlled benchmark. All candidate
models use the same site-year target, rows, spatial folds, and XGBoost settings.
Only the feature set changes:

- static fine-scale habitat predictors;
- static predictors plus calendar year;
- static predictors plus lagged reef outbreak probability;
- static predictors plus both year and lagged outbreak probability.

Lagged adjusted-hindcast features are the leakage-safe default. A separate
provisional operational experiment trains with same-year adjusted hindcasts
through 2025 and deploys the adjusted 2026 forecast. It remains explicitly
provisional until upstream hindcast and forecast provenance is confirmed.

## Run from the repository root

```powershell
Rscript DataProcessing/cots_sdm_v2/scripts/00_check_inputs.R
Rscript DataProcessing/cots_sdm_v2/scripts/01_prepare_data.R
Rscript DataProcessing/cots_sdm_v2/scripts/02_fit_benchmark.R
Rscript DataProcessing/cots_sdm_v2/scripts/03_fit_rolling_origin.R
Rscript DataProcessing/cots_sdm_v2/scripts/04_diagnose_targets.R
Rscript DataProcessing/cots_sdm_v2/scripts/05_prepare_manta_negatives.R
Rscript DataProcessing/cots_sdm_v2/scripts/06_fit_operational_2026.R
Rscript DataProcessing/cots_sdm_v2/scripts/08_fit_additive_bridge.R
Rscript DataProcessing/cots_sdm_v2/scripts/07_build_visual_comparison.R
Rscript DataProcessing/cots_sdm_v2/scripts/09_fit_site_habitat.R
Rscript DataProcessing/cots_sdm_v2/scripts/10_summarise_site_habitat.R
Rscript DataProcessing/cots_sdm_v2/scripts/11_site_habitat_manta_weight_sensitivity.R
Rscript DataProcessing/cots_sdm_v2/scripts/12_fit_rg_coverage_hybrid.R
Rscript DataProcessing/cots_sdm_v2/scripts/13_assess_rank_signal.R
Rscript DataProcessing/cots_sdm_v2/scripts/14_validate_legacy_year_signal.R
Rscript DataProcessing/cots_sdm_v2/scripts/15_fit_rank_preserving_calibration.R
Rscript DataProcessing/cots_sdm_v2/scripts/16_assess_rank_stability.R
Rscript DataProcessing/cots_sdm_v2/scripts/17_build_blind_expert_review.R
Rscript DataProcessing/cots_sdm_v2/scripts/18_freeze_legacy_clean_year.R
Rscript DataProcessing/cots_sdm_v2/scripts/19_fit_multiscale_resolution.R
Rscript DataProcessing/cots_sdm_v2/scripts/20_fit_presence_only_multiscale.R
Rscript DataProcessing/cots_sdm_v2/scripts/21_fit_hmsc_preliminary.R
Rscript DataProcessing/cots_sdm_v2/scripts/22_fit_legacy_unweighted_scale_check.R
Rscript DataProcessing/cots_sdm_v2/scripts/23_fit_hmsc_joint_production.R
Rscript DataProcessing/cots_sdm_v2/scripts/24_build_multiscale_reef_maps.R
Rscript DataProcessing/cots_sdm_v2/scripts/25_build_multiscale_report_figures.R
Rscript DataProcessing/cots_sdm_v2/scripts/26_summarise_hmsc_joint.R
Rscript DataProcessing/cots_sdm_v2/scripts/27_diagnose_hmsc_effort_signal.R
Rscript DataProcessing/cots_sdm_v2/scripts/28_fit_two_part_fixed_effort.R
Rscript DataProcessing/cots_sdm_v2/scripts/29_plot_two_part_diagnostics.R
Rscript DataProcessing/cots_sdm_v2/scripts/31_build_candidate_reef_maps.R
Rscript DataProcessing/cots_sdm_v2/scripts/32_build_hmsc_reef_maps.R
Rscript DataProcessing/cots_sdm_v2/scripts/34_audit_joint_response_data.R
Rscript DataProcessing/cots_sdm_v2/scripts/33_build_hmsc_joint_roadmap.R
Rscript DataProcessing/cots_sdm_v2/scripts/30_update_model_reports.R
Rscript DataProcessing/cots_sdm_v2/scripts/35_prepare_hmsc_v2_visits.R
Rscript DataProcessing/cots_sdm_v2/scripts/36_prepare_hmsc_v2_predictors.R
Rscript DataProcessing/cots_sdm_v2/scripts/37_fit_hmsc_v2_ecological.R
Rscript DataProcessing/cots_sdm_v2/scripts/38_fit_hmsc_v2_observation.R
Rscript DataProcessing/cots_sdm_v2/scripts/39_validate_hmsc_v2_blocked.R --candidate=geomorphic_g360
Rscript DataProcessing/cots_sdm_v2/scripts/40_build_hmsc_v2_reef_maps.R --candidate=geomorphic_g360
Rscript DataProcessing/cots_sdm_v2/scripts/41_screen_geomorphic_scales.R
Rscript DataProcessing/cots_sdm_v2/scripts/42_fit_hmsc_v2_geomorphic_g360.R
Rscript DataProcessing/cots_sdm_v2/scripts/43_summarise_hmsc_v2_geomorphic_g360.R
Rscript DataProcessing/cots_sdm_v2/scripts/91_build_final_synthesis_figures.R
Rscript DataProcessing/cots_sdm_v2/tests/smoke_test.R
```

`01_prepare_data.R` writes derived tables under `artifacts/`. The benchmark
writes fold predictions and metrics under `outputs/`, and fitted fold models
under `models/`. These generated files are ignored by Git.

## Data contracts

- Cull observations are read from the `Cull` sheet of the configured workbook.
- `site_year_all` pools all dives at a site within a calendar year.
- `site_year_first_visit` pools only dives on the first surveyed date at a site
  within a calendar year, reducing within-year depletion bias.
- `high_density` is based on pooled CPUE (`total COTS / bottom time >= 0.02`).
  `CPUE_max` is retained only as a diagnostic and is never the primary target.
- Predictor extraction preserves an immutable `site_id` and extraction ID.
- All observations from one reef are assigned to the same spatial fold.
- Forward-year validation trains only on years before the test year.
- Adjusted hindcasts (1991-2025) and the adjusted 2026 forecast are ingested as
  separate products; the forecast never enters lag-only benchmark training.
- Repeated manta non-detections are prepared as presence-stage evidence with
  method conflicts and effort retained. They are not appended as unweighted
  pseudo-absences to the high-density target.
- The site-habitat ablation excludes all reef-level outbreak predictors. It
  compares the original clean and RG predictor families, categorical
  geomorphic classes, and conservative site-aggregated manta zeros using
  held-out cull sites and within-reef discrimination.

See `IMPLEMENTATION_PLAN.md` for the staged roadmap and decision gates, and
`BENCHMARK_FINDINGS.md` for the first controlled and rolling-origin results.
`VISUAL_CHECK_FINDINGS.md` records the matched-scale map review and the
between-reef versus within-reef score decomposition.
`SITE_HABITAT_FINDINGS.md` records the clean/RG, geomorphic, and manta-zero
ablation focused specifically on within-reef high-density habitat.
`RANK_SIGNAL_FINDINGS.md` explains why compressed probability values are not
discarded, compares annual and lifetime targets, and defines the rank-based
acceptance tests for fine-scale habitat.
`EXPERT_REVIEW_PROTOCOL.md` defines the blinded spatial-plausibility review and
the separation between expert evidence and model tuning.
`reports/MULTISCALE_RESOLUTION_MAXENT_HMSC_REPORT.md` records the exact legacy
reproduction, multiscale and resolution experiments, presence-only review,
joint Hmsc analysis, and percentile-based reef map comparison.
`reports/TWO_PART_FIXED_EFFORT_ITERATION_01.md` records the separated static
habitat, outbreak-state and fixed-effort observation experiment.
`reports/HMSC_VISUAL_JOINT_MODELLING_ROADMAP.md` compares the fitted Hmsc
fixed-protocol surface with the incumbent on four reefs and defines the data
and modelling work needed for a genuinely joint ecological model.
`reports/HMSC_V2_GEOMORPHIC_G360_ITERATION_03.md` records the independently
selected, converged reduced 360 m geomorphic Hmsc benchmark and embeds its four
production-posterior reef comparisons.
`reports/FINAL_MODEL_SYNTHESIS_20260914.md` provides the final five-figure
summary of model validation, density validation, variable importance, partial
dependence and four-reef spatial review.
`reports/CURRENT_BEST_MODEL.md` is the promotion-controlled current model
record; it changes incumbent only when statistical, provenance and expert-map
gates all pass.
Every successor iteration must generate two to four named reef comparisons
against the current-best raster. Both reports embed those maps, and the report
generator stops if their manifest or image files are missing.
