# COTS fine-scale SDM: final layers

This folder delivers the COTS layers and figures for the technical report. It
is self-contained: every input is a frozen file pinned by SHA-256 in
[`config.R`](config.R), and nothing here refits the exploratory models except
for lengthening the Hmsc MCMC chains and the forward-validation refits.

## Layer roles

| Role | Layer | Use it for |
|---|---|---|
| **Primary search-planning layer** | `cots_search_planning_2026_30m.tif` (+ `_priority_`) | Planning manta and cull effort *within* a reef: where to tow and dive first. RG + year with the clean + year fallback, at 30 m, ranked within each reef polygon. |
| Known aggregation sites | `cots_cull_history_kernel_120m.tif` | Briefing crews, especially new contractors, on where COTS have been found before. A record of past finds, not a habitat map. |
| Explanatory model | `cots_hmsc_joint_120m.tif` (+ uncertainty) | Showing how Cull, SALAD and RHIS evidence respond to habitat (fig04-05), with uncertainty. Not the operational map. |

All layers are in GDA2020 / MGA zone 55 and live in `outputs/layers/`, with
checksums in [`layer_manifest.csv`](outputs/layers/layer_manifest.csv).

| Layer | File | Resolution | Answers |
|---|---|---|---|
| **Primary search-planning layer** | `cots_search_planning_2026_30m.tif` | 30 m | Band 1: P(Cull CPUE >= 0.02) in 2026 from RG + year, with clean + year (calibrated) where RG layers are missing. Band 2: percentile of the cell within its GBR reef polygon. |
| | `cots_search_planning_priority_2026_30m.tif` | 30 m | 1 = top 10 % of the reef, 2 = 10-20 %, 3 = 20-30 %, 0 = rest. |
| **Clean + year** | `cots_legacy_clean_year2026_30m.tif` | 30 m | P(Cull CPUE >= 0.02) at a site surveyed in 2026. XGBoost on the 12-layer clean stack + survey year. Full clean-stack coverage. The operational incumbent. |
| **RG + year** | `cots_legacy_rg_year2026_30m.tif` | 30 m | Same target, with five ReefGuide layers added (wave period/height, turbidity, bathymetry, slope). NA where the RG layers are missing. |
| | `cots_legacy_rg_year2026_clean_fallback_30m.tif` (+ `_source`) | 30 m | RG + year where available, otherwise clean + year calibrated onto the RG scale. Source: 1 = RG, 2 = clean fallback. Band 1 of the primary layer. |
| **HMSC joint** | `cots_hmsc_joint_120m.tif` | 120 m | Posterior mean P(problematic-density COTS given an active reef-year): a static habitat score from a joint model of Cull, SALAD and RHIS evidence. |
| | `cots_hmsc_joint_uncertainty_120m.tif` | 120 m | Bands: posterior SD, 5th and 95th percentiles. |
| | `cots_hmsc_joint_source_120m.tif` | 120 m | Component used: 1 enhanced, 2 clean multiscale, 3 clean native. |
| **HMSC + Cull history** | `cots_hmsc_history_combined_120m.tif` | 120 m | Forecast for 2027: the HMSC habitat probability updated by prior Cull records, (kappa p + W H) / (kappa + W). Equals HMSC wherever there is no nearby history. |
| | `cots_cull_history_kernel_120m.tif` | 120 m | Band 1: kernel-weighted prior high-density rate H (500 m Gaussian, all active reef-years <= 2026); band 2: evidence weight W. |

The primary layer mixes two models within most reefs. RG layers cover the
whole reef polygon on only 54 reefs; 2,586 reefs are mixed (median 65 % of
cells from RG) and 920 use the clean fallback alone
(`tables/search_planning_reef_sources.csv`). The fallback is calibrated onto
the RG probability scale before ranking. However, the two models rank cells
quite differently where both exist (Spearman 0.34), so ranks near a coverage
boundary within a reef are less reliable.

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

## Operational efficiency (11-13)

**How much manta and cull effort could the primary layer save, and at what
cost?** In each scenario, effort is cut *within* each survey (manta) or
reef-year (cull) by skipping the lowest-ranked tows or sites. It is always
compared with skipping the same number at random, because any cut "saves"
effort; only the difference from random is the layer's contribution. Results
are in `tables/operations_*.csv` and fig11.

Why these numbers replace the earlier estimates
(`DataProcessing/SDM_Validation_Report.qmd`):
- **No random baseline:** the earlier report credited the model with the
  whole saving.
- **One threshold for the whole GBR:** it applied the 88th/95th percentile
  of the raster everywhere, so whole low-scoring reefs dropped out.
- **In-sample scoring:** it scored cull dives the model had been trained on.

Here, cull sites are scored by held-out predictions (each site excluded from
the model that scored it). Manta tows are independent of the layer's
training.

**At 60 % of effort kept, 2022-2026** (`operations_headline.csv`):

| | Primary layer | Random |
|---|---:|---:|
| Manta tows kept | 61 % (7,400 km of tows skipped) | 61 % |
| Strong-evidence tows kept | **67.6 %** | 60.9 % |
| Outbreak reef-years still detected (at least 1 strong tow) | 88.6 % | 87.0 % |
| COTS at the kept cull sites | **69.9 %** | 62.5 % |
| Ghost dives (0 COTS) avoided | **40.9 %** (727 dives) | 36.5 % (650) |
| Catch rate at the kept sites vs all | **+4.2 %** | -0.7 % |

The cull comparison covers 278 reef-years, 11,878 dives, 1,779 ghost dives
and 188,000 COTS. Over all years since 2018 the layer's advantage over random
is smaller: +1.7 points of COTS reached and +3.6 points of ghost dives
avoided.

How to report this defensibly:
- **The large effort saving is the cut itself, not the model.** Skipping 40 %
  of sites saves about 18,400 diver-hours over 2022-2026 (about AUD 17.5 M at
  the earlier report's 950 AUD/hour) whichever sites are skipped. It also
  leaves roughly 30 % of COTS unculled. The layer's contribution is choosing
  better sites to skip: about 14,000 more COTS reached and 77 more ghost dives
  avoided than random, for the same number of sites.
- **Manta surveillance can be thinned with little loss of reef-level
  detection, layer or not.** Outbreak reef-years have many positive tows, so
  60 % of tows still flags 87-89 % of them. The layer raises the share of
  strong-evidence *tows* kept (+6.7 points). However, that does **not**
  translate into finding more *aggregations*; see the next section. Do not
  use the layer to choose which manta tows to skip.
- **On reefs that have never been culled, the layer is no better than random
  for manta.** Part of its advantage elsewhere may come from where crews have
  historically dived.
- **Cull history plus habitat (forward-only) does better again on ghost
  dives:** 46.5 % avoided vs 36.4 % at random, and 72.5 % of COTS reached. It
  is the stronger tool on reefs with a cull record.
- **Costs:** "bottom time" is summed per dive in diver-minutes; the 950
  AUD/hour rate is carried over from the earlier report
  (`cfg$operations$diver_hour_cost_aud`). Confirm both before quoting dollars.

**Search-plan example (fig12, `search_plan_examples.csv`).** On Batt (manta
survey 4 Aug 2025; cull 2025) and Eyrie (manta survey 26 Oct 2023; cull 2024),
keeping the layer's top 60 % of tows retains 32/42 and 9/11 strong-evidence
tows. The planned dive sites sit on the top-ranked diveable cells (16-0.5 m,
at least 400 m apart), ordered into a short route. These reefs are in the
layer's training data, so fig12 illustrates the planning logic; the held-out
numbers above are the evidence.

## Manta thinning and aggregation detection (14)

**Question:** how many manta tows could be cut without missing COTS
aggregations?

**Confound and fix.** The layer is trained on Cull dives, and Cull crews dive
where manta found COTS or scars, so scoring tows with the deployed layer
partly scores it against its own source. Script 14 refits the RG + year
XGBoost (legacy hyperparameters) five times, each time leaving out a fifth of
the reefs. It then scores every reef's tows with the model that never saw that
reef's Cull data.

**What counts as an aggregation.** Within a survey (reef x date),
strong-evidence tows lying within 500 m of each other are linked into one
aggregation. An aggregation is found if at least one of its tows is kept,
since that tow triggers a closer search.

**Strategies compared,** each keeping the same number of tows per survey:
- the layer's top-ranked tows;
- random selection (exact expectation);
- systematic thinning: every n-th tow around the reef, averaged over starting
  points.

Intervals are 500 bootstrap resamples of surveys (`manta_aggregation_cv_*.csv`,
fig13).

**Aggregations of 2+ strong-evidence tows, all years** (1,900 aggregations in
1,082 surveys):

| Tows kept | Layer (leave-reef-out) | Layer (deployed, confounded) | Random | Systematic |
|---|---:|---:|---:|---:|
| 80 % | 96.2 % (95.3-97.2) | 96.9 % | 98.8 % | **99.7 %** |
| 70 % | 93.0 % | 93.7 % | 97.0 % | **98.7 %** |
| 60 % | 88.1 % (86.6-89.5) | 90.4 % | 94.0 % | **96.6 %** |
| 50 % | 81.8 % | 85.2 % | 89.7 % | **93.0 %** |

**Largest cut that still finds at least 95 % of aggregations:**

| Strategy | All years | 2022-2026 |
|---|---:|---:|
| Systematic | **40 %** (30-40 % with bootstrap assurance) | **40 %** (30 %) |
| Random | 30 % | 30 % |
| Layer | 20 % | 20 % (10 %) |

What this means for the report:
- **Tow every second or third segment, spread evenly around the reef.**
  Systematic thinning can cut about a third of manta tows while still finding
  at least 95 % of aggregations.
- **Ranking tows by the layer finds fewer aggregations than cutting at
  random.** The layer concentrates kept tows in a few high-scoring stretches,
  so aggregations in whole low-scoring stretches are missed. It does slightly
  better on isolated single strong tows (64 % vs 61 % found at 60 % kept),
  but those are not aggregations.
- **The confound is real but modest.** The deployed layer looks 1-3 points
  better than its cross-validated version (about 3 points at 50 % kept).
- **Where the layer belongs:** in deciding where to cull and where to search
  more closely *after* a detection, not in deciding which surveillance tows to
  drop.
- **Caveat:** systematic thinning assumes tows can be ordered around the reef.
  Tow order is not recorded, so tows are ordered by bearing from the survey
  centre, which suits perimeter surveys.

## Figures

All in `outputs/figures/`, 7.5 in wide at 300 dpi, built by `07` from the tables in `outputs/tables/` (run it last).

| Figure | Shows |
|---|---|
| `fig01_gbr_overview.png` | The three layers across the GBR |
| `fig02_four_reef_comparison.png` | Eyrie, Batt, Broadhurst, Lady Musgrave: within-reef percentiles + HMSC uncertainty |
| `fig02b_four_reef_comparison_imagery.png` | Same as fig02 over Esri World Imagery; ranks below the 70th percentile drawn near-transparent (tiles cached in `cache/tiles`) |
| `fig03_manta_independent_check.png` | Matched independent Manta check |
| `fig04_hmsc_target_effects.png` | HMSC habitat effects on the mapped Cull target |
| `fig05_hmsc_joint_coefficients.png` | HMSC coefficients across all 12 jointly fitted responses |
| `fig06_legacy_variable_importance.png` | XGBoost gain, clean vs RG |
| `fig07_layer_agreement.png` | Within-reef agreement between layers |
| `fig08_coverage_by_source.png` | Which model fills each cell |
| `fig09_forward_validation.png` | Forward validation: habitat, reef context, history, combined |
| `fig10_history_combined_four_reefs.png` | Four reefs: habitat vs Cull history vs combined, with history weight |
| `fig11_operational_efficiency.png` | Manta and cull effort cut by the primary layer vs at random (2022-2026) |
| `fig12_search_plans_batt_eyrie.png` | Example manta and dive plans on Batt and Eyrie over satellite imagery |
| `fig13_manta_aggregation_cv.png` | Manta thinning vs aggregation detection: layer (leave-reef-out and deployed), random, systematic |

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
Rscript DataProcessing/cots_sdm_v2/final/scripts/11_layer_search_planning.R     # primary layer: within-reef percentile + priority
Rscript DataProcessing/cots_sdm_v2/final/scripts/12_operational_efficiency.R    # manta / cull effort scenarios vs random
Rscript DataProcessing/cots_sdm_v2/final/scripts/13_figure_search_plans.R       # fig12 (Batt, Eyrie)
Rscript DataProcessing/cots_sdm_v2/final/scripts/14_manta_aggregation_cv.R      # leave-reef-out manta thinning vs aggregations (~30 min first run)
Rscript DataProcessing/cots_sdm_v2/final/scripts/07_build_report_figures.R      # last: uses outputs of 06-14
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
  scripts/        00-14 (see Run order; 07 runs last)
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
