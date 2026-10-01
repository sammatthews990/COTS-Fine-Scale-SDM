# COTS fine-scale SDM: final layers

This folder delivers the COTS habitat layers (three static layers plus a history-updated forecast) for the technical report,
and the figures that describe them. It is self-contained: every input is a
frozen file pinned by SHA-256 in [`config.R`](config.R), and nothing here
refits the exploratory models except for lengthening the Hmsc MCMC chains.

## The three layers

All layers are in GDA2020 / MGA zone 55 and live in `outputs/layers/`, with
checksums in [`layer_manifest.csv`](outputs/layers/layer_manifest.csv).

| Layer | File | Resolution | Answers |
|---|---|---|---|
| **Clean + year** | `cots_legacy_clean_year2026_30m.tif` | 30 m | P(Cull CPUE ≥ 0.02) at a site surveyed in 2026. XGBoost on the 12-layer clean stack + survey year. Full clean-stack coverage. The operational incumbent. |
| **RG + year** | `cots_legacy_rg_year2026_30m.tif` | 30 m | Same target, with five ReefGuide layers added (wave period/height, turbidity, bathymetry, slope). NA where the RG layers are missing. |
| | `cots_legacy_rg_year2026_clean_fallback_30m.tif` (+ `_source`) | 30 m | RG + year where available, otherwise clean + year calibrated onto the RG scale. Source: 1 = RG, 2 = clean fallback. |
| **HMSC joint** | `cots_hmsc_joint_120m.tif` | 120 m | Posterior mean P(problematic-density COTS ∣ active reef-year): a static habitat score from a joint model of Cull, SALAD and RHIS evidence. |
| | `cots_hmsc_joint_uncertainty_120m.tif` | 120 m | Bands: posterior SD, 5th and 95th percentiles. |
| | `cots_hmsc_joint_source_120m.tif` | 120 m | Component used: 1 enhanced, 2 clean multiscale, 3 clean native. |
| **HMSC + Cull history** | `cots_hmsc_history_combined_120m.tif` | 120 m | Forecast for 2027: the HMSC habitat probability updated by prior Cull records, (κ·p + W·H) / (κ + W). Equals HMSC wherever there is no nearby history. |
| | `cots_cull_history_kernel_120m.tif` | 120 m | Band 1: kernel-weighted prior high-density rate H (500 m Gaussian, all active reef-years ≤ 2026); band 2: evidence weight W. |

### Reading them together

The layers answer different questions and are **not on a common probability
scale**:

- The two legacy layers are *annual* probabilities. Survey year is their largest
  single driver, so their absolute values mostly reflect 2026 prevalence, and
  their within-reef ordering carries the fine-scale signal.
- The HMSC layer is *conditional on an active reef-year* and excludes year,
  reef identity, outbreak probability, effort and method by design. Multiply it
  by a separate reef-year outbreak probability if an annual risk is needed; do
  not let the outbreak layer reorder cells within a reef.
- For operational targeting within a reef, compare **within-reef percentiles**
  (fig02), not raw values.

We do not blend the layers into a single ensemble. Averaging an annual
probability with a conditional habitat probability mixes two estimands. Use
fig07 (within-reef agreement) and the four-reef panels instead to show where the
layers concur.

## Key results (run of 2026-10-01)

### HMSC convergence

All three components pass the convergence gate (minimum Beta ESS ≥ 400,
maximum upper PSRF ≤ 1.1). Before this run, no v4 component had passed.

| Component | Training rows | Min ESS | Max upper PSRF | Apparent AUC (target) | Share of mapped cells |
|---|---:|---:|---:|---:|---:|
| enhanced_joint | 14,120 | 705 | 1.012 | 0.600 | 24.6% |
| clean_joint | 15,512 | 1,645 | 1.020 | 0.576 | 71.8% |
| clean_native_joint | 15,571 | 409 | 1.018 | 0.565 | 3.6% |

Source: `outputs/tables/hmsc_*_production_summary.csv`, `hmsc_joint_coverage.csv`.

### Validation

Each model's own held-out validation uses a different protocol, so these
numbers are not a leaderboard:

| Layer | Held-out validation (own protocol) | Source |
|---|---|---|
| Clean + year | Out-of-fold AUC 0.616 (5-fold, 50 km blocks); within-reef-year AUC 0.548 in the v2 matched held-site benchmark | `model_performance_metrics_clean_20260820_1041.rds`; `../RANK_SIGNAL_FINDINGS.md` |
| RG + year | Out-of-fold AUC 0.690 (training 0.968, so strongly overfit in-sample); within-reef-year AUC 0.584 in the same benchmark | as above |
| HMSC joint (enhanced) | Within-reef-year AUC 0.570 held 3 km cells, 0.545 leave-reef, 0.515 forward (train ≤ 2023) | `../iterations/hmsc_corrected_v4/outputs/enhanced_joint_*_metrics.csv` (preliminary-length CV refits) |

**Matched independent check (fig03).** Every layer scores the same 117,318
Manta tows in 1,003 active reef-years. Within each reef-year, tows with strong
COTS evidence are compared against strict negatives. Manta is not a response in
any of the models.

| Layer | Within-reef-year AUC (95% CI) | Top-20% lift (95% CI) |
|---|---|---|
| Clean + year | 0.516 (0.500–0.532) | 1.11 (1.06–1.15) |
| RG + year | 0.576 (0.548–0.602) | 1.23 (1.18–1.29) |
| HMSC joint | 0.567 (0.549–0.585) | 1.15 (1.11–1.19) |

### Agreement and coverage

- Within reefs, the layers rank cells only moderately alike. Median Spearman
  correlations: clean vs RG 0.41, clean vs HMSC 0.33, RG vs HMSC 0.24 (fig07).
- Across GBR reef polygons, clean + year and HMSC cover 98% of 120 m cells;
  RG + year covers 61%. Across the whole predictor domain, RG covers only 26%
  of 30 m cells and the HMSC enhanced component 25% of 120 m cells (fig08).

## Forecasting with Cull history (09-10)

The static habitat layers reach within-reef AUCs of only about 0.57. Cull site
records show that high-density areas persist, and that the persistent pattern
lives at about 0.5-1 km, so the Cull midpoint error is not the main limit (see
`outputs/tables/forward_*`). The history layer turns that persistence into a
forecast.

**Validation is strictly forward.** For each test year T from 2022 to 2026,
every score uses only data from earlier years. The habitat model is the
clean_joint HMSC refit on rows before T. History comes only from earlier Cull
site-years. Scoring uses first-visit Cull site-years in active reef-years that
contain both outcomes: 4,170 site-years in 315 reef-years.

| Score | All | With nearby history (3,212) | Little or none (958) | New sites (1,664) |
|---|---:|---:|---:|---:|
| Habitat (HMSC refit before T) | 0.569 | 0.584 | 0.554 | 0.524 |
| Habitat + reef context | 0.582 | 0.600 | 0.520 | 0.514 |
| Effort only (prior dives at the site) | 0.636 | 0.656 | 0.491 | 0.513 |
| History kernel (fixed) | 0.654 | 0.705 | 0.465 | 0.476 |
| **Habitat + history (fixed, deployed)** | **0.681** (0.640-0.718) | **0.716** | **0.540** | **0.553** |

All values are within-reef-year AUC, with a 95% bootstrap interval over
reef-years for the deployed score. The "fixed" setting (σ = 500 m, no recency
decay, κ = 2) was written into the script before any results were seen. A
per-year nested re-tuning on T − 1 performed no better overall (0.676). It was
unstable, mostly choosing a 1-year memory, and it degraded where history is
sparse (fig09, `forward_nested_selection.csv`).

How to report this:

- **The gain comes from history, and history partly encodes operator
  behaviour.** Effort alone (how often crews dived a site) already scores 0.636.
  The combined layer favours places crews have found COTS before. That is
  operationally useful, but it is not independent evidence of habitat.
- **Where history is thin, the blend falls back to habitat.** At new sites it
  scores 0.553 against 0.524 for habitat alone. History is within reach of only
  8.6% of 120 m cells (W > 2, where history dominates, at 2.8%), so the combined
  layer is the HMSC layer across most of the domain.
- **Reef context did not help robustly.** The 0.5-1 km predictors from `08`
  were windward and seaward position, distance to reef edge, upwind open water,
  distance to channels, reef area, and eReefs current speed and flow
  steadiness. They added 0.013 to habitat AUC overall, with overlapping
  intervals. The effect varied across years (−0.04 to +0.07) and was worse at
  new sites and where history is sparse. They are therefore not in the
  production HMSC. The forward context refits also mixed worse (minimum ESS
  45-148).
- **Larval connectivity was not added.** It is a reef-level quantity, constant
  within a reef, so it cannot change within-reef ranking. It belongs in the
  reef-year activation model. No connectivity data are in the repository.
- **The ReefScan check is inconclusive** (`reefscan_presence_check.csv`). The
  237 presence-only detections fall on 29 reefs, mostly the Swains, 2024-26.
  Without survey tracks that recorded no detections, the test mostly reflects
  where ReefScan was towed. Detections sat below the within-reef median for
  clean + year (0.44) and HMSC (0.41), and near it for RG (0.49) and the
  combined layer (0.53). It needs ReefScan effort tracks to become a real test.

eReefs source: GBR1 v2.0 hydrodynamic model, whole-period mean (`all-one`),
from the AIMS THREDDS server (`thredds.ereefs.aims.gov.au`). Mean current
speed and u/v are taken at −9 m, filled upward from −5.35, −2.35 and −0.5 m
where the model is shallower. The raw slices are cached in
`cache/ereefs/` (downloaded 2026-10-01).

## Figures

All in `outputs/figures/`, 7.5 in wide at 300 dpi, built by `07` from the tables in `outputs/tables/` (run it last).

| Figure | Shows |
|---|---|
| `fig01_gbr_overview.png` | The three layers across the GBR |
| `fig02_four_reef_comparison.png` | Eyrie, Batt, Broadhurst, Bait: within-reef percentiles + HMSC uncertainty |
| `fig03_manta_independent_check.png` | Matched independent Manta check |
| `fig04_hmsc_target_effects.png` | HMSC habitat effects on the mapped Cull target |
| `fig05_hmsc_joint_coefficients.png` | HMSC coefficients across all 12 jointly fitted responses |
| `fig06_legacy_variable_importance.png` | XGBoost gain, clean vs RG |
| `fig07_layer_agreement.png` | Within-reef agreement between layers |
| `fig08_coverage_by_source.png` | Which model fills each cell |
| `fig09_forward_validation.png` | Forward validation: habitat, reef context, history, combined |
| `fig10_history_combined_four_reefs.png` | Four reefs: habitat vs Cull history vs combined, with history weight |

## Limitations to state in the report

1. **Fine-scale signal is weak in every layer.** The best within-reef AUCs are
   around 0.57–0.58. The layers prioritise search effort; they do not locate
   aggregations.
2. **HMSC forward transfer is near chance** (0.515 when trained ≤ 2023 and
   tested on later years). The production chains fix convergence, not transfer.
3. **The HMSC map is mostly a fallback map.** 72% of cells come from the
   clean-multiscale component, including nearly all of the Far Northern and
   Torres Strait sector, which reads uniformly high in fig01. Use the source
   layer when interpreting that region.
4. **The RG fallback composite joins two weakly related rankings.** Where both
   legacy layers exist, their Spearman correlation is 0.34, so expect visible
   seams at the RG coverage boundary.
5. **The Manta check covers only reefs with RG coverage** (all three layers
   must score each tow). Clean + year's weak result there may not hold where
   it is the only layer.
6. **Expert zone review is still pending** for all layers (see
   `../reports/CURRENT_BEST_MODEL.md`).

## Run order

From the repository root:

```powershell
Rscript DataProcessing/cots_sdm_v2/final/scripts/00_check_inputs.R            # verify all frozen inputs (minutes)
Rscript DataProcessing/cots_sdm_v2/final/scripts/01_layer_legacy_clean_year.R # verify + copy frozen clean layer
Rscript DataProcessing/cots_sdm_v2/final/scripts/02_layer_legacy_rg_year.R    # predict RG + year 2026 (~1 h)
Rscript DataProcessing/cots_sdm_v2/final/scripts/03_fit_hmsc_production.R --component=enhanced_joint
Rscript DataProcessing/cots_sdm_v2/final/scripts/03_fit_hmsc_production.R --component=clean_joint
Rscript DataProcessing/cots_sdm_v2/final/scripts/03_fit_hmsc_production.R --component=clean_native_joint
Rscript DataProcessing/cots_sdm_v2/final/scripts/04_build_support_grid_120m.R --workers=18
Rscript DataProcessing/cots_sdm_v2/final/scripts/05_layer_hmsc_joint.R
Rscript DataProcessing/cots_sdm_v2/final/scripts/06_compare_layers.R
Rscript DataProcessing/cots_sdm_v2/final/scripts/08_build_reef_context_120m.R  # reef geometry + eReefs (downloads ~0.5 GB once)
# forward habitat refits, one per cutoff year and variant (independent; run in parallel)
Rscript DataProcessing/cots_sdm_v2/final/scripts/09a_fit_forward_habitat.R --cutoff=2021 --variant=base   # ... 2021-2026 x base/context
Rscript DataProcessing/cots_sdm_v2/final/scripts/09b_validate_forward.R
Rscript DataProcessing/cots_sdm_v2/final/scripts/10_layer_history_combined.R
Rscript DataProcessing/cots_sdm_v2/final/scripts/07_build_report_figures.R      # last: uses outputs of 06-10
```

The three `03` fits are independent and can run in parallel (each uses 4
cores; about 50 minutes each). `03 --smoke=1` runs a seconds-long pipeline
test whose outputs are tagged `_SMOKE` and never feed a layer. `04` caches each
stage under `cache/` and resumes if interrupted.

## Folder layout

```
final/
  config.R        frozen inputs (path + SHA-256), chain lengths, seeds
  R/              common.R, hmsc.R, support_grid.R, rasters.R, figures.R, context.R, forward.R
  scripts/        00-10 (see Run order; 07 runs last)
  outputs/
    layers/       delivered GeoTIFFs + layer_manifest.csv   (tifs git-ignored)
    tables/       every number quoted in the report            (tracked)
    figures/      fig01-fig08 PNGs                              (tracked)
  models/         Hmsc production posteriors                    (git-ignored)
  cache/          120 m support grid, intermediate rasters      (git-ignored)
```

## Provenance

| Input | Produced by |
|---|---|
| Clean + year model and 2026 raster | `DataProcessing/fit_predict_cots_sdm.R` via `refit_fy2526_clean.R`, run `20260820_1041` |
| RG + year model | `DataProcessing/fit_predict_cots_sdm.R` via `refit_fy2526_rg.R`, run `20260820_1217` (its raster was never written; `02` predicts it) |
| Hmsc training tables (`*_fit_data_120m.rds`) | `cots_sdm_v2/iterations/hmsc_corrected_v4/scripts/04`, `04b`, `06b`, built on the v3 support-matched cells and the v4 active-reef-year conditioning and Manta evidence policy |
| Manta tow evidence and footprints | `hmsc_corrected_v4/scripts/01c`, `hmsc_geomorphic_v3` footprint links |

The decisions behind these choices are recorded in
`../reports/FINAL_MODEL_SYNTHESIS_20260914.md` and
`../reports/ALL_MODEL_REVIEW_20260911.md`.

### What changed relative to the exploratory work

- **RG + year 2026** was predicted for the first time.
- **HMSC chains** were lengthened from 500 samples × thin 5 (preliminary) to
  1,000 × thin 30 after 5,000 transient iterations, 4 chains. The responses,
  traits, terrain PCA basis and training rows are unchanged.
- **HMSC was mapped across the whole domain** for the first time (previously
  only four reefs). Each 120 m cell's predictors are computed exactly as for
  training: a 4×4 block mean of the 30 m layers, the geomorphic class at the
  cell centre, and point-wise multiscale Gaussian terrain. `04` checks this
  against the 12,736 training cells before it maps anything.
- **HMSC uncertainty** comes from 400 posterior draws, with fallback calibration
  applied to every draw.
