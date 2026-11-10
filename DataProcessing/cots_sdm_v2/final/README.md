# COTS fine-scale SDM: final layers

## Goal

Predict **where on a reef COTS are most likely to be found**, so the COTS
Control Program can:

1. **Manta tow less** without missing COTS aggregations.
2. **Search fewer cull sites** without missing problematic densities
   (CPUE ≥ 0.02). This means fewer dives where COTS are absent (ghost dives)
   or below the action threshold (low-density dives).

Every layer and analysis below is judged against these two uses, on held-out
data, against a no-information baseline at the same effort. Guidance for
agents working in the repository is in [`AGENTS.md`](../../../AGENTS.md).

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
| **Manta model (within reef-year)** | `cots_manta_within_30m.tif` | 30 m | Band 1: relative log-odds that a manta tow over the cell records strong COTS evidence, against the reef-year average (comparable within a reef, not a probability). Band 2: percentile within the GBR reef polygon. XGBoost on the RG stack, no survey year, reef-year offset; trained on 198,852 tows (script 18). NA outside RG coverage; values away from tow paths are extrapolations. Best for ranking manta tows (section 17). |
| **Exploration 2027** | `cots_exploration_2027_30m.tif` | 30 m | Band 1: habitat percentile (primary layer + manta ensemble) on every eligible reef cell, i.e. 0.5-16 m deep, > 100 m from designated cull sites and > 150 m from manta tows of the last 5 seasons (NA elsewhere). Bands 2-3: distance to recent tows and to cull sites (m). Reefs with designated cull sites only (script 21). |
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

## Higher-density management reassessment (23; 16 October 2026)

The original analyses use CPUE >=0.02. Script 23 separately evaluates the
**unchanged** saved rankings against **CPUE >0.04**, preserving <0.02 positive
dives as low-density, zero dives as ghost, and 0.02–0.04 inclusive as intermediate.
It uses existing held-out rankings from 15–17 and the forward evidence scores
from 19; no model is refitted and no 2027 priorities change. Habitat fits in the
forward comparison remain spatially held out, not restricted to prior years.

- Forward evidence test, 2019–2026: at 95% >0.04 retention, candidate site-years
  searched are 90.2% (95% CI 85.4–93.4) for habitat plus recent evidence, 94.0%
  (89.7–95.5) for habitat alone, and 95.6% (95.3–95.9) for random selection.
  The combined ranking avoids 8.9% (6.2–12.9) of low-density and 11.2%
  (8.0–16.8) of ghost dives. At 99% high-density retention it searches 98.2%
  (97.2–98.7) of candidate site-years. These are retrospective trade-offs,
  not agreed tolerances or demonstrated field savings.
- Stopping at the top 80% rank cutoff retains 89.7% (86.4–92.7) of >0.04 dives
  in that forward test. Fixed classes are not validated stopping rules.
- Effort is the actual selected count divided by all candidate site-years,
  including never-dived sites. Random expected outcomes match the site count on
  each reef. This differs slightly from nominal rank cutoffs on small reefs or
  tied scores. Matched comparisons interpolate the first target crossing; 500
  paired reef-year bootstrap resamples repeat that interpolation.
- The full cleaned cull source spans 2002–2026, with 49,130 valid records across
  361 reefs; records matched to the designated register cover 4,299 sites on
  339 reefs (`cull_record_scope.csv`). Each evaluation uses its stated subset.

Outputs are `cull_high_density_fixed.csv`, `_matched.csv`, `_differences.csv`,
`_composition.csv`, `_priority_classes.csv`, `_input_manifest.csv`, and
`cull_record_scope.csv`, all in `outputs/tables/`. New inputs are pinned under
`assessment_*` in config.R; refreshing upstream caches requires explicitly
reviewing those pins. Assertions verify recovered forward sample counts,
all-site rank denominators, monotonic retention and complete retention at full
effort. Figure 20 is built by 07; use `--high_density_only=1` after 23.

The trial's proposed 10–15% search and low-density/empty-dive reductions still
require a numerical tolerance for missed >0.04 dives, observations from new
searches, and accounting for annual or pre-closure checks. No annual vessel-time
saving or validated fifth-wave suppression benefit is inferred.

## Coral cover and a full refit at CPUE > 0.04 (24; 19 October 2026)

Script 24 asks two questions of the cull-site ranking: does tow coral cover
help, and does refitting every part at >0.04 beat the deployed ≥0.02 model?
The 2027 product does not use coral and is unchanged.

**Method.**
- Each tow's hard coral code, and its recently dead coral code (recorded
  from Nov 2018), is set to its bin midpoint (README sheet of the manta
  workbook, pinned as `manta_workbook`). Last season's values are smoothed onto
  each site with the 300 m evidence kernel.
- 52 % of candidate site-years had a tow nearby (mean nearby hard coral
  17.6 %).
- Refit at each threshold: the leave-reef-out habitat XGBoost, the
  threshold-based evidence terms, and the evidence GLM.
- Forward test as in 19 (2019–2026, 861 reef-years). Every ranking is scored
  against both definitions; 200 paired reef-year bootstraps.
- At 0.02 the evidence terms reproduce 16 exactly. The base model reproduces
  19's `cull_all` (AUC 0.576 vs 0.577; 80.5 % vs 81.1 % of sites at 90 %
  retention).

**Coral.**
- More coral nearby, more problematic dives, with or without strong manta
  tows nearby.
  - Without strong tows, >0.04 share rises from 27 % (<10 % coral) to 39 %
    (≥30 %); with them, from 50 % to 66 % (`threshold_coral_bins.csv`).
  - Coefficient per 10 % cover: +0.11 at 0.02 (p = 0.03) and +0.15 at 0.04
    (p = 0.003), with reef-year-clustered SEs.
- The coral × strong-manta interaction is negative and not significant, and
  recently dead coral adds nothing. The effect is additive.
- Forward test, 0.02-trained, at 90 % retention: coral cuts sites searched by
  1.4 points (95 % CI 0.4–2.3) on the 0.02 target and by 1.5 (0.4–2.8) on the
  0.04 target. Low-density and ghost dives avoided do not change (≤1 point,
  CIs span 0).
- No gain at 95 % retention. AUC gains are under 0.01 (CIs span 0).
- Coral is habitat evidence, but it exists only where crews towed last
  season. COTS also remove coral, so high coral nearby may mark reefs where
  the outbreak has not yet peaked.

**Threshold.**
- Training at >0.04 is no better, even on the >0.04 target. At 90 %
  retention of >0.04 dives:
  - the 0.04-trained base model searches 3.0 points more sites (0.0–5.5)
    than the 0.02-trained one;
  - its habitat rank alone has a within reef-year AUC of 0.495, against
    0.563 for the 0.02-trained rank.
- Likely reason: >0.04 site-years are fewer (3,095 vs 4,591 of 10,028
  training site-years). The ≥0.02 response carries the same spatial signal
  with more positives.
- The deployed (0.02-trained) ranking keeps 90 % of >0.04 dives at 80.4 % of
  sites, against 80.5 % for 90 % of ≥0.02 dives.

Outputs: `threshold_coral_forward.csv`, `_differences.csv`,
`_coefficients.csv`, `_lr_tests.csv`, `_bins.csv`, `_checks.csv` in
`outputs/tables/`; per-dive forward ranks in `cache/threshold_coral_forward.rds`.

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

## Cull site selection: avoiding non-actionable dives (15)

**Question:** if crews dived only the higher-ranked *designated cull sites*
on each reef, what share of problematic dives (CPUE >= 0.02) would be missed,
and what share of low-density dives (COTS present but < 0.02) and ghost dives
(0 COTS) would be avoided?

**Method:**
- **Candidates:** every designated cull site on the reef (EoTR cull-site KMZ:
  15,884 sites on 576 reefs, of which 4,299 have ever been dived). Sites that
  were never dived because there was little evidence of COTS are included.
- **Matching:** 93.8 % of dives are matched to a site, by site name with a
  spatial fallback.
- **Held-out scores:** each reef is scored by the RG + year model refitted
  without that reef's Cull data (5-fold leave-reef-out). A site's score is the
  mean prediction over its polygon for that year, ranked among all designated
  sites on the reef.
- **Uncertainty:** 500 bootstrap resamples of reef-years
  (`cull_site_selection_*.csv`, fig14).

**Held out, 2018-2026** (29,813 dives in 918 reef-years; 61 % problematic,
21 % low-density, 18 % ghost):

| Designated sites dived | Problematic dives kept | COTS at problematic dives kept | Low-density dives avoided | Ghost dives avoided |
|---|---:|---:|---:|---:|
| top 80 % | 84.9 % (81.9-87.7) | 84.6 % | 16.1 % | 19.8 % |
| top 70 % | 77.2 % (72.9-80.8) | 72.0 % | 24.6 % | 29.5 % |
| top 60 % | 68.7 % (64.1-72.8) | 64.3 % | 35.2 % | 38.9 % |
| top 50 % | 59.8 % (55.5-64.5) | 55.4 % | 44.9 % | 48.2 % |

What this means for the report:
- **The layer cannot be used to rule out cull sites without losing
  problematic dives in roughly the same proportion.** With the top 70 % of
  sites, 75.6 % of dives still happen. A layer with no information would
  then keep 75.6 % of problematic dives and avoid 24.4 % of low-density
  dives; the layer keeps 77.2 % and avoids 24.6 %. The held-out within
  reef-year AUC for problematic vs low-density dives is only 0.57.
- **It is somewhat better at flagging ghost dives** (AUC 0.60): 29.5 %
  avoided vs 24.4 % with no information.
- **Why it fails here:** whether a site with COTS is above or below 0.02 CPUE
  in a given visit is driven mainly by outbreak stage, earlier culling and
  aggregation movement, not by static habitat. Static predictors cannot
  separate "some COTS" from "enough COTS".
- **The deployed layer looks better because it is in-sample** (AUC 0.61 and
  0.66). Report the held-out figures.
- **Never-dived candidates:** 43 % of site-years the layer ranks in a reef's
  top 30 % are at sites that have never been dived. The layer can point to
  untested sites worth a first look, but the evidence above says it should
  set the *order* of searching, not exclude sites.

## Habitat plus last season's evidence (16)

**Question:** same as 15 (search fewer designated sites on a reef while keeping
problematic dives), but the ranking may also use last season's Cull dives and
manta tows. Nothing from the season being searched is used.

**Method:**
- **Evidence terms** for every designated site, from season T-1 only:
  - Cull dives at the site: any problematic dive, or dived with none problematic.
  - Cull dives nearby: kernel (sigma 300 m) problematic rate and dive effort.
  - Manta tows nearby: kernel count of strong-evidence tows, and total tows.
- **Scores**, each ranked among all designated sites on the reef:
  - habitat only (held-out rank from 15);
  - evidence only (logistic model);
  - habitat + evidence (logistic model);
  - a simple rule: sites flagged last season first (problematic dive at the
    site, or a strong manta tow nearby), then by habitat.
- **Held out:** the logistic models are fitted on other reefs' dives (5-fold
  leave-reef-out).
- **Matched retention:** evidence moves both problematic and low-density dives
  up the ranking. So scores are compared when each keeps the *same* share of
  problematic dives, set in advance (80-95 %).
- **Uncertainty:** 500 bootstrap resamples of reef-years, paired across scores
  (`cull_site_evidence_*.csv`, fig15).
- **Sensitivity:** an evidence window of the last three seasons gives the same
  picture.

**Held out, last-season evidence, 2018-2026** (30,037 dives in 926
reef-years). Values are habitat only → habitat + evidence; brackets are the
95 % interval of the gain.

| Problematic dives kept | Designated sites searched | Low-density dives avoided | Ghost dives avoided | Non-problematic dive hours avoided |
|---|---:|---:|---:|---:|
| 90 % | 88 % → 82 % (−11 to −1) | 10 % → 16 % (+2 to +10) | 13 % → 21 % (+4 to +13) | 11 % → 19 % (+4 to +12) |
| 80 % | 73 % → 67 % (−11 to −3) | 22 % → 28 % (+2 to +10) | 27 % → 36 % (+5 to +14) | 23 % → 31 % (+4 to +13) |

**2022-2026** (16,971 dives in 628 reef-years). The gains are larger:

| Problematic dives kept | Designated sites searched | Low-density dives avoided | Ghost dives avoided |
|---|---:|---:|---:|
| 90 % | 88 % → 79 % | 10 % → 19 % (+6 to +13) | 12 % → 25 % (+8 to +17) |
| 80 % | 74 % → 61 % | 20 % → 32 % (+6 to +20) | 25 % → 41 % (+8 to +24) |

What this means for the report:
- **Last season's dives and manta tows are what let crews search less.** For
  the same problematic dives, the evidence-informed ranking searches about 6
  to 13 points fewer designated sites. It also avoids roughly 6 to 12 points
  more low-density dives and 8 to 16 points more ghost dives than habitat
  alone.
- **Most of that comes from the evidence, not the habitat layer.** Evidence
  alone performs almost the same as habitat + evidence. The habitat layer is
  the tie-breaker among sites with no recent evidence, and the only
  information on the 30 % of reef-years whose reef was not visited last
  season.
- **Separating "some COTS" from "enough COTS" is still hard.** Among dives
  actually made, the within reef-year AUC for problematic vs low-density dives
  rises only from 0.57 to 0.59. The gain comes mainly from:
  - ranking quiet and never-dived sites lower;
  - catching more ghost dives (the rule reaches AUC 0.65 vs ghost).
- **What the models learned:**
  - A strong manta tow nearby, and a high problematic rate among nearby dives,
    push a site up.
  - A site dived last season with nothing problematic is pushed down strongly.
  - A site that *was* problematic last season ranks lower than its
    neighbourhood suggests, consistent with culling having reduced it.
- **The simple rule** ("last season's flags first, then habitat") searches the
  fewest sites. It avoids fewer low-density dives than the fitted score, so
  use the fitted score for low-density avoidance.
- **Caveat:** "sites searched" credits a ranking for agreeing with where crews
  chose to dive. Crews already use last season's knowledge, so part of this
  gain is codifying current practice. Problematic densities at sites that were
  never dived cannot be verified.
- **Where nothing helps:** reefs not visited last season (11 % of dives,
  within reef-year AUC about 0.53 for every score).

## Training on manta tows instead of cull dives (17)

**Question:** the layers above are trained on Cull site-years (CPUE ≥ 0.02).
Would training on manta tows serve the goal better? Manta tows have about 20
times more rows and cover places that were searched but never dived. The
model would then be validated on the cull data.

**Models,** all with the same predictors (RG stack) and the legacy XGBoost
settings, refitted on the same 5 folds of reefs:
- **Cull-trained:** the deployed target (Cull site-year CPUE ≥ 0.02), plus
  survey year.
- **Manta-trained:** strong-evidence tow (direct COTS, ≥ 2 counted scars,
  legacy "common" scars) vs strict negative, plus survey year. Weak single
  scars are dropped, and predictors are averaged along each tow.
- **Manta-trained, within reef-year:** the same response, but each
  reef-year's strong-tow rate is a fixed offset and survey year is dropped.
  The trees can therefore only learn *where within a reef-year* COTS are,
  which is the ranking the programme uses.

**Held out, both goals** (`manta_response_*.csv`, fig16):

| | Cull-trained | Manta-trained | Manta, within reef-year |
|---|---:|---:|---:|
| Tows: within reef-year AUC, strong vs negative | 0.547 (0.516-0.574) | 0.575 (0.547-0.599) | **0.612** (0.578-0.636) |
| Tows: top-20 % lift | 1.15 | 1.27 | **1.33** |
| Aggregations found, 60 % of tows kept | 88.2 % | 88.4 % | **91.0 %** |
| Aggregations found, 80 % of tows kept | 96.7 % | 97.1 % | **97.8 %** |
| Cull: AUC problematic vs low-density (2018-2026) | **0.567** | 0.482 | 0.566 |
| Cull: AUC problematic vs ghost (2018-2026) | 0.583 | 0.504 | **0.590** |
| Cull: low-density dives avoided, 90 % of problematic kept | **11.4 %** | 7.0 % | 8.7 % |
| ... with last season's evidence (as in 16) | 15.6 % | 15.6 % | 15.4 % |

Tow counts: 198,852 tows (16,199 strong evidence) in 1,015 reef-years.
Aggregations: 1,904 groups of 2 or more strong tows. Cull: 30,037 dives at
designated sites. For comparison, at 60 % of tows kept, random thinning finds
94.0 % of aggregations and systematic thinning 96.5 %.

What this means for the report:
- **For the manta goal, the within reef-year manta model is the better
  layer.** It ranks held-out tows best (AUC 0.61 vs 0.55 for the cull-trained
  layer). It finds 2.8 points more aggregations at 60 % of tows kept
  (95 % CI +1.1 to +4.6).
- **It still does not beat thinning tows evenly.** Every model-ranked
  thinning finds fewer aggregations than random or systematic thinning,
  because ranking clusters the kept tows. The largest cut that still finds
  ≥ 95 % of aggregations is 20 % for every model, against 30 % (random) and
  40 % (systematic). The model's place in manta work is choosing *where to
  look harder* (for example, where to add or repeat tows), not which tows to
  drop.
- **For the cull goal, training on manta does not help.**
  - With survey year, the manta model ranks cull sites at chance (AUC 0.48
    and 0.50): 70 % of its gain goes on survey year, which predicts
    year-to-year prevalence, not location.
  - The within reef-year version matches the cull-trained layer on 2018-2026
    and is weaker in 2022-2026 (0.51 vs 0.55).
  - Once last season's dives and tows are added, all three habitat models
    give the same result. Recent evidence, not the habitat model, drives
    cull-site selection.
- **Why more data did not transfer:**
  - Manta tows see COTS and scars visible from the surface along the reef
    perimeter. A problematic dive measures how many COTS a diver finds in a
    10 ha site. The two are related but not the same quantity.
  - Manta coverage is the perimeter, so manta-model scores away from tow
    paths (for example, back reef) are extrapolations.
- **Interpretation:** the within reef-year manta model has the cleanest
  meaning: the relative chance of seeing strong COTS evidence at a location,
  given the reef is active that year. Its top predictors are wave period and
  height, bathymetry and slope. If a single habitat layer must serve both
  goals, it is the better candidate (better for manta, equal for cull). Its
  30 m layer is `cots_manta_within_30m.tif` (script 18; shown in fig02/02b).

**Ensembles (17).** Two ensembles of the cull-trained and within reef-year
manta ranks were tested: their average, and both ranks as separate terms in
the evidence model. Neither helped once last season's evidence was added.
At 80 % of problematic dives kept (2018-2026), ensemble + evidence searched
65 % of sites, against 67 % for cull + evidence (2022-2026: 59 % vs 62 %).
Low-density dives avoided were the same (27 %), and the intervals overlap
throughout.

## Year-by-year forward test of the site ranking (19)

**Method.** For each test season 2019-2026, the habitat + evidence model is
fitted only on dives from earlier seasons, then every designated site on
each reef dived that season is ranked. Three training windows (all earlier
seasons, last 5, last 3) and three habitat terms (cull, ensemble, both) are
compared with habitat alone (`site_ranking_forward*.csv`).

**Sites searched and low-density dives avoided at 90 % of problematic dives
kept.** Habitat alone → cull habitat + evidence, all earlier seasons:

| Test season | Dives | Sites searched | Low-density dives avoided | Within reef-year AUC, problematic vs low |
|---|---:|---|---|---|
| 2019 | 4,380 | 88 % → 85 % | 11 % → 14 % | 0.59 → 0.50 |
| 2020 | 3,311 | 75 % → 73 % | 23 % → 21 % | 0.61 → 0.59 |
| 2021 | 3,769 | 78 % → 80 % | 19 % → 17 % | 0.54 → 0.56 |
| 2022 | 4,109 | 80 % → 81 % | 13 % → 15 % | 0.55 → 0.61 |
| 2023 | 4,119 | 89 % → 74 % | 10 % → 21 % | 0.59 → 0.61 |
| 2024 | 3,910 | 88 % → 84 % | 7 % → 15 % | 0.55 → 0.61 |
| 2025 | 3,676 | 88 % → 84 % | 11 % → 18 % | 0.52 → 0.60 |
| 2026 (Jan-May) | 1,157 | 96 % → 76 % | 5 % → 25 % | 0.53 → 0.71 |
| **2019-2026** | 28,431 | **86 % → 81 %** | **11 % → 17 %** | **0.57 → 0.58** |

What this means:
- **The evidence model has worked every season since 2022, most clearly
  from 2023.** It did not help in 2019-2021, when fewer earlier seasons were
  available to learn from.
- **Training window:** all earlier seasons (0.81 of sites, 17 % low-density
  dives avoided) matched the last 5 (0.81, 17 %) and beat the last 3 (0.86,
  13 %).
- **Habitat term:** cull, ensemble and both performed the same, so the
  simpler cull-only model is used.

## Cull-site priority for 2027 (20)

Every EoTR designated cull site (15,884 sites on 576 reefs) is ranked within
its reef for the 2027 season. Outputs are in `outputs/site_ranking/`:

| File | Contents |
|---|---|
| `cots_cull_site_priority_2027.kmz` | The EoTR cull-site KMZ, unchanged in structure (one folder per reef; same site names, IDs and polygons). Each polygon is coloured by priority (A red, B orange, C yellow, D pale yellow, E grey), and its data and pop-up gain *Priority 2027*, *Rank on reef*, *P(problematic if dived)* and *Last season*. |
| `cots_cull_site_priority_2027.csv` | One row per site: Cull site ID, priority, rank on reef, score, the evidence behind it (dives and problematic dives last season, nearby problematic rate, dive effort, strong and total manta tows) and both habitat ranks. |
| `cots_cull_site_priority_2027_reefs.csv` | Per reef: sites, priority-A sites, sites flagged last season, and whether the reef was dived or towed last season. |
| `cots_cull_site_priority_classes_2027.csv` | What each class kept in the forward test. |

**Model** (`config.R` → `site_ranking`; coefficients in
`tables/site_ranking_2027_coefficients.csv`):
- Logistic model of "a dive at this site is problematic" on the cull habitat
  rank and last season's evidence.
- Fitted on 28,431 dives, 2019-2026.
- **Last season for 2027** is the most recent 12 months of data: 1 June 2025
  to 31 May 2026, 3,521 dives and 17,902 manta tows (1,978 strong). Calendar
  2026 holds only January-May, and the model was trained on a full season of
  previous-season evidence.
- **Habitat rank for 2027:** the deployed RG + year layer (clean fallback),
  averaged over each site polygon. Training used held-out ranks, so habitat
  carries a little more weight here than it earned. Its coefficient is small
  (+0.26 across the whole reef range), so the effect is minor.

**Priority classes:** share of the reef's designated sites, highest score
first (`priority_breaks`; five classes since 2 November 2026, when the old
50-80 % band C was split at 66 % so the protocol's dive rule can sit at A-C;
scores and ranks unchanged). Retention is from the forward test, 2019-2026:

| Searching through | Sites on the reef | Problematic dives kept | Low-density dives avoided | Ghost dives avoided |
|---|---:|---:|---:|---:|
| A | top 20 % | 33 % | 75 % | 79 % |
| A + B | top 50 % | 64 % | 42 % | 49 % |
| A + B + C | top 66 % | 79 % | 30 % | 36 % |
| A to D | top 80 % | 89 % | 18 % | 22 % |
| all (incl. E) | 100 % | 100 % | 0 % | 0 % |

How to use it:
- **Work down the classes; don't stop at A.** Class E (the bottom 20 % of a
  reef's sites) is the realistic candidate to skip or leave until last; the
  proposed protocol (35) also leaves D for dives without strong evidence.
  Skipping it gave up about 11 % of problematic dives and avoided 18 % of
  low-density and 22 % of ghost dives.
- **Most sites have no recent evidence.** 10,654 sites (67 %) sit on reefs
  with no dives or tows in the last 12 months, so they are ranked by habitat
  alone. On such reefs the ranking is close to chance (AUC about 0.53;
  section 16). Treat it as a first search order and update it once manta
  tows come in.
- **Evidence drives the top of the ranking:**
  - 678 sites had a problematic dive last season;
  - 809 have strong manta evidence nearby;
  - 293 were dived last season with nothing problematic, which strongly
    lowers their rank.
- **Reading *P(problematic if dived)*:** it is the model's chance that a dive
  *made* at the site is problematic. It is learned from dives crews chose to
  make, so use it to compare sites on a reef, not as an absolute risk for an
  unvisited reef.
- **Re-run each season.** Update `evidence_end` and `training_years` in
  `config.R`, then run 20.

## Density-aware priority classes (20b) — 6 November 2026

**Why.** The classes in 20 are shares of each reef's sites, so every reef
gets 20 % A sites whatever its COTS density: in 2027, 1,985 of the 2,951 A
sites were on reefs with no dives or tows in the last 12 months. 20b keeps
each reef's order from 20 but lets the reef's density set how many of its
sites reach each class.

**Method** (`scripts/20b_site_ranking_density_2027.R`):
- **Reef density:** logistic model of "a dive here is problematic" on the
  reef's previous-season share of problematic dives and of strong tows (each
  shrunk by half a dive or tow), log dive count, no-dive / no-tow indicators,
  the GBR-wide outbreak forecast for the previous year (Stewart et al.
  hybrid model; `data/gbrPreds2026Adj_20262408.csv`, hindcasts in
  `gbrPredsAdj_20262408.csv`, both now pinned in `config.R`), plus the site's
  within-reef percentile. Reef density = mean P over the reef's sites.
- **Unknown reefs** (no dives or tows in the last 12 months; 394 of 576
  reefs): the model's rows for such reefs are reef-years crews dived anyway,
  usually after tows had just found COTS, so it put them at a median 43 %
  problematic. Their P is scaled by the chance such a reef needed culling at
  all: of 872 cull-site reef-years towed with no evidence the season before
  (2012–2026), 19 % had a problematic dive (14 % at forecast < 0.02, about
  25 % above 0.05; logistic on the forecast). Median unknown reef: 7 %.
  A reef without a forecast would get the median forecast; none in 2027.
- **Scaled quotas:** cumulative class share = plogis(logit(break) + Δ),
  Δ = logit(reef density) − logit(typical culled reef, 0.35). A typical
  reef keeps A 20 / B 50 / C 66 / D 80 %; the order within a reef is
  unchanged.
- **Absolute cut-offs on P were tried first and rejected:** the within-reef
  slope is about 1 logit across a whole reef, so density swamped it. Every
  reef above about 30 % got all its sites in A–C (Batt and Eyrie 100 % A)
  and every unknown reef none, leaving crews no order on most reefs.

**Forward test** (`site_ranking_density_forward.csv`, `_by_season.csv`;
seasons 2020–2026, models fitted on earlier seasons, 200 reef-year
bootstraps; per-season thresholds; 802 reef-years, 24,051 dives):

| At the same total sites as within-reef A–C (65 %) | Within-reef (20) | Density-aware (scaled) | Difference |
|---|---:|---:|---:|
| > 0.04 dives kept | 81.1 % | 88.5 % | **+7.3 (2.7 to 11.0)** |
| Problematic dives kept | 80.6 % | 88.2 % | +7.6 (3.4 to 10.9) |
| Low-density dives avoided | 29.9 % | 23.3 % | −6.6 (−8.3 to −4.5) |
| Ghost dives avoided | 36.3 % | 36.1 % | −0.2 (−2.6 to 3.0) |
| Sites to keep 90 % of > 0.04 dives | 78.5 % | 70.3 % | **−8.2 (−16.2 to −1.4)** |
| Ghost avoided at 90 % of > 0.04 kept | 23.6 % | 31.6 % | +8.0 (−0.1 to 15.6) |

- Better in every season (> 0.04 kept at equal effort: 2020 96 vs 86 %,
  2024 79 vs 77 %, 2026 92 vs 87 %).
- **Absolute cut-offs** keep more (+13.8 points without the unknown-reef
  scaling) but 13–27 % of dived reef-years get no A–C site at all.
- **The forecast adds little where there are data:** evidence-only and
  evidence + forecast differ by about 1 point; it matters for unknown reefs,
  which this test cannot score fairly (every unknown reef-year in it was
  dived, usually because tows found COTS). Scaling them down is penalised
  here; in use, the first tows of a visit re-rank the reef (33).
- **Where the gain comes from:** mostly moving A–C sites between reefs, i.e.
  the which-reef decision, and dives counted are where crews chose to dive.
  Within a reef the order is the same as 20.

**2027 product** (`site_ranking/cots_cull_site_priority_density_2027.csv`,
`_reefs.csv`, `.kmz`; `tables/site_ranking_density_2027_classes.csv`,
`_coefficients.csv`):

| Reef density (expected share of problematic dives) | Reefs | A share | A–C share |
|---|---:|---:|---:|
| < 10 % (mostly unknown reefs) | 376 | 3 % | 21 % |
| 10–20 % | 21 | 5 % | 31 % |
| 20–35 % | 36 | 15 % | 58 % |
| 35–50 % | 61 | 27 % | 75 % |
| 50–70 % | 61 | 41 % | 84 % |
| > 70 % | 21 | 65 % | 93 % |

- A sites 1,956 (within-reef 2,951); A–C 6,206 (10,170). Unknown reefs drop
  from 6,831 A–C sites to 2,319; reefs with recent dives and tows rise from
  2,616 to 3,024. Batt (69 %): 47 A of 92; Eyrie (74 %): 19 of 33.
- The KMZ pop-up gives both classes, the rank, P and the reef density with
  its source. 21 adds the exploration zones to this KMZ too.
- **Status:** an alternative to the within-reef product, not yet adopted.
  Absolute P mixes "conditional on crews diving" (reefs with evidence) with
  "unconditional" (scaled unknown reefs); the classes for unknown reefs are
  a survey order, and should be refreshed after the first tows.

## Exploration zones: outside the cull sites and recent tow paths (21)

**Why.** Every ranking above learns from where crews have already dived and
towed. Used alone, it would keep consolidating effort on the same sites and
tow paths. The exploration layer points to high-ranked habitat that the
current search footprint does not cover.

**Check on held-out tows (17).** A tow counts as *new ground* if no tow lay
within 300 m in the previous 5 seasons. New ground found strong COTS
evidence more often than repeat ground (13.6 % vs 5.8 % of tows; partly
because new ground is often surveyed early in an outbreak). The habitat
models still ranked tows there (`exploration_new_ground_check.csv`):

| On new ground | Within reef-year AUC | Top-20 % lift |
|---|---|---|
| Cull-trained layer | 0.56 (0.51-0.60) | 1.11 |
| Manta model (within reef-year) | **0.62** (0.59-0.65) | **1.28** |
| Ensemble of the two | 0.61 (0.57-0.65) | 1.23 |

**Method:**
- **Habitat score:** the ensemble of the primary layer and the manta model
  (mean of the two within-reef percentiles, re-ranked within the reef).
  Where the manta layer is missing, the primary layer is used alone.
- **Eligible cells:** reef cells 0.5-16 m deep (clean DEM), more than 100 m
  from any designated cull site and more than 150 m from any manta tow in the
  last 5 seasons (June 2021 – May 2026).
- **Exploration zones (6 November 2026; replace the discrete areas):**
  - habitat percentile smoothed with a 5 × 5 cell (150 m) mean and
    re-ranked within the reef, so zones read as ground rather than
    scattered cells;
  - all eligible ground in the reef's top 40, 30, 20 and 10 % of habitat,
    as four nested bands (each band clipped inside the one below);
  - outlines closed (+30 m / −30 m buffer), simplified to 20 m to cut the
    cell stairs and rounded by Chaikin corner-cutting; holes under 0.5 ha
    filled, pieces under 0.5 ha dropped, clipped to the reef and kept
    outside cull sites, simplified to 8 m. (Kernel smoothing took ~35 s per
    reef on these large bands, about 5 hours in all.)
  - in the KMZ the bands are stacked semi-transparent purple fills, so the
    shading darkens continuously towards the best-ranked habitat.
- **Why zones, not areas.** The earlier product drew the top 10 discrete
  areas per reef (top 30 % habitat, split into ~9 ha chunks, ≥ 350 m apart).
  In review, nothing separated those areas from the surrounding unsearched
  habitat: the chunking and spacing rules decided where they fell. The zones
  show all rarely searched habitat, graded by rank, and leave the choice of
  where to look to the crew.

**Outputs:**
- `layers/cots_exploration_2027_30m.tif`, three bands:
  1. exploration score (150 m-smoothed habitat percentile on eligible cells);
  2. distance to the nearest manta tow in the last 5 seasons (m);
  3. distance to the nearest designated cull site (m).

  It covers reefs with designated cull sites only.
- `site_ranking/cots_exploration_zones_2027.csv` (one row per reef and band:
  area, pieces, largest piece and its centre, mean habitat percentile, depth
  range) and `.gpkg`. Zones on 483 reefs; median per reef 38 ha in the top
  10 % band (421 reefs), 68 ha top 20 %, 95 ha top 30 %, 128 ha top 40 %.
  The old `cots_exploration_areas_2027.*` files are removed.
- Both 2027 KMZs (20 and 20b) gain an **"Exploration zones 2027 (purple,
  darker = better habitat)"** folder: one sub-folder per reef, one
  placemark per band (pop-up: band, area, depth). Each band fills at about
  25 % opacity, so the nested bands stack towards the best habitat; the top
  10 % band has a thin white outline. KMZs grew from 12 to 18 MB.
- fig17 shows Batt and Eyrie.

What this means:
- **The zones are large on big, perimeter-towed reefs.** On Batt the top-40 %
  band is 2,600 ha, mostly the outer northern reef flat and lagoon edge
  beyond the tow paths; the top-10 % band is 440 ha in 45 pieces. Eyrie's
  best ground is the north-western reef flat and the northern back reef near
  the island, which no tow has crossed in 5 seasons.
- **Status:** these are places to look, not predictions. The manta model is
  extrapolated away from tow paths, which is exactly where these zones lie.
  It ranks new ground about as well as old (AUC 0.62), but no layer is
  strong.
- **Suggested use:** pick a scouting tow or dive in the darkest band that is
  reachable on the day, and log nil results so the next season's evidence
  includes them.
- **Settings:** clearances, depth band, band floors, smoothing and minimum
  piece size are set at the top of script 21. The script takes about 17
  minutes.

## Metacommunity calibration reporting snapshot (22; 15 October 2026)

This reporting step summarises archived ADRIA/COTSMod calibration outputs; it
does not refit or promote a model. Inputs from the sibling `ADRIA.jl` repository
are pinned in `config.R` (`ADRIA_DIR` can override its location). Script 22 checks
every new pin, writes reported numbers with `write_table()`, and keeps the plotted
LTMP observation series in the ignored cache. The earlier two-run comparison in
`report/external/metacommunity/` is preserved as superseded historical material.

- **Replicated candidates:** W, C3a and C6, five ecological–connectivity seed pairs,
  run `20261014T_48_candidates_1991_v2`. Each matches 7/7 observed historical peaks
  in every seed under score v2. Mean ± seed SD losses: W 6.79 ± 0.79, C3a 6.72 ±
  0.70, C6 7.13 ± 0.85 (`metacommunity_candidates.csv`). C3a has the larger Lizard
  first wave and better coral fit, but higher troughs. C6 trades a worse coral fit
  for lower troughs. These SDs are not confidence intervals.
- **Observation dependence:** G2′ passes at A½ = 30 adults/ha are W 5/5, C3a 0/5,
  C6 3/5; at A½ = 10 they are 0/5, 0/5 and 3/5. These conditional calibration gates
  do not establish independent predictive performance.
- **Historical culling:** run 42, primary allocation sigma 210 m. At Eyrie in
  2021–24, 88% of requested removals cannot be applied because simulated COTS are
  absent; the domain share is 70% (`metacommunity_culling_removal_summary.csv`).
  North Direction's second-wave maximum falls from 0.82 to 0.55 COTS/tow
  (`metacommunity_culling_diagnostic.csv`). This is a consistency diagnostic,
  not a control-effectiveness estimate.
- **Latest v0.4 diagnostics:** 15 October runs 08–13 extend through 2026, with one
  seed pair and culling off. Score v3 is distinct from v2. Original → corrected
  outbreak boundary: loss 15.61 → 12.63 with earlier coral thermal settings;
  11.34 → 11.34 with revised settings. The lowest exploratory chlorophyll-modified
  loss is 10.99 (`metacommunity_2026_diagnostics.csv`). Shared cull/SALAD evidence
  in the corrected boundary and assessment precludes an independent validation
  interpretation. These runs are not a newly adopted working configuration.

Held-out reefs/years, regional transfer and sensitivity gates remain incomplete.
No adaptive starting/stopping or primary-outbreak suppression benefits are claimed.
The report and manuscript describe this snapshot. Their Discussions outline the
planned historical assessment of culling from 2012 to the present and future
fifth-wave suppression scenarios, without claiming completed efficacy results. The pipeline and operational SDM results above are unchanged.

## Figures

All in `outputs/figures/`, 7.5 in wide at 300 dpi, built by `07` from the tables in `outputs/tables/` (run it last). Use `07 --metacommunity_only=1` for figures 18–19 alone after script 22.

| Figure | Shows |
|---|---|
| `fig01_gbr_overview.png` | The three layers across the GBR |
| `fig02_four_reef_comparison.png` | Eyrie, Batt, Broadhurst, Lady Musgrave: within-reef percentiles of clean + year, RG + year, the manta model, HMSC and HMSC + history, plus HMSC uncertainty |
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
| `fig14_cull_search_area.png` | Diving only the top-ranked designated cull sites: problematic vs low-density vs ghost dives kept (held out and deployed) |
| `fig15_cull_search_area_evidence.png` | Habitat vs habitat + last season's dives and manta tows, at matched problematic-dive retention: sites searched, low-density and ghost dives avoided |
| `fig16_manta_response_model.png` | Cull-trained vs manta-trained habitat models on both goals: tow ranking, aggregations found at 60 % of tows, low-density dives avoided |
| `fig17_search_2027_exploration.png` | Batt and Eyrie: 2027 cull-site priority, manta tows of the last 5 seasons and exploration zones (nested habitat bands) over imagery |
| `fig22_protocol_batt_eyrie.png` | Batt, Eyrie, Big Broadhurst and Lady Musgrave: the proposed protocol (35) under the density-aware v2 classes (37): classes, the follow-up tow zone (A-C + 100 m), exploration zones and the last 12 months of strong vs other tows |
| `fig33_protocol_density_tradeoff.png` | Tow zone (A-B vs A-C) and dive rule (A-B/A-C/A-D, three evidence requirements) trade-offs, within-reef vs density-aware v1 and v2 classes (36) |
| `fig34_ghost_per_miss.png` | Ghost dives avoided per dive above 0.04 missed: evidence-only rules, random, dive rules A-B/A-C/A-D under within-reef, v1 and v2 classes, and zone + rule protocols, with 95 % intervals (36) |
| `fig23_workflow_decision.png` | Decision diagram for a reef visit under the proposed workflow (35): open sites, tow zone, strong-evidence trigger, A–C dive rule |
| `fig18_metacommunity_candidates.png` | W, C3a and C6 historical trajectories: five-seed means and SD, with LTMP observations |
| `fig21_manta_surveillance_highlight.png` | Plain-language surveillance summary of the same aggregation retention results: habitat-priority, evenly spaced and random tows; detailed figure 13 retained |
| `fig20_high_density_retention.png` | Existing forward rankings re-evaluated at >0.04 CPUE: high-density retention versus actual candidate site-years searched, with same-effort random expectation and 95% bootstrap intervals |
| `fig19_metacommunity_2026_diagnostic.png` | Single-seed recent build-up: original/corrected boundaries and earlier/revised coral settings |

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
Rscript DataProcessing/cots_sdm_v2/final/scripts/15_cull_search_area_cv.R       # leave-reef-out ranking of designated cull sites (~40 min first run)
Rscript DataProcessing/cots_sdm_v2/final/scripts/16_cull_search_area_evidence.R # habitat + last season's dives and manta tows (~10 min; needs 15)
Rscript DataProcessing/cots_sdm_v2/final/scripts/17_manta_response_model.R      # cull- vs manta-trained habitat models on both goals (~15 min; needs 14-16)
Rscript DataProcessing/cots_sdm_v2/final/scripts/18_layer_manta_within.R        # manta within reef-year 30 m layer (~30 min; needs 11, 14)
Rscript DataProcessing/cots_sdm_v2/final/scripts/19_site_ranking_forward.R      # year-by-year forward test of the site ranking (~6 min; needs 16-17)
Rscript DataProcessing/cots_sdm_v2/final/scripts/20_site_ranking_2027.R         # 2027 cull-site priority: CSV + KMZ (~5 min; needs 18-19)
Rscript DataProcessing/cots_sdm_v2/final/scripts/20b_site_ranking_density_2027.R # density-aware classes: forward test, CSV + KMZ (~15 min; needs 19, 20)
Rscript DataProcessing/cots_sdm_v2/final/scripts/21_exploration_areas.R         # exploration layer + zones, added to the 2027 KMZs (needs 17, 20, 20b; rerun after 37 for the v2 KMZ)
Rscript DataProcessing/cots_sdm_v2/final/scripts/22_metacommunity_results.R      # frozen ADRIA calibration reporting snapshot
Rscript DataProcessing/cots_sdm_v2/final/scripts/23_cull_high_density_evaluation.R # fixed rankings at >0.04; no refits
Rscript DataProcessing/cots_sdm_v2/final/scripts/24_threshold_coral_refit.R       # coral terms + full refit at >0.04, forward test (~5 min; needs 15-17)
Rscript DataProcessing/cots_sdm_v2/final/scripts/25_calibration_report_snapshot.R # 21 October ADRIA snapshot, independent of script 22
Rscript DataProcessing/cots_sdm_v2/final/scripts/26_larval_concentration_report.R # review 06: concentration tables/map preparation
# 27-32: see their README sections (27 October 2026)
Rscript DataProcessing/cots_sdm_v2/final/scripts/33_deployment_validation.R      # deployment-realistic forward test (~25 min; needs 14-17)
Rscript DataProcessing/cots_sdm_v2/final/scripts/34_manta_zone_strategy.R        # option A: tow the zone, tows guide dives (~10 min; needs 33)
Rscript DataProcessing/cots_sdm_v2/final/scripts/35_protocol_statistics.R        # proposed protocol: dive rule, zone, trial planning (~3 min; needs 34)
Rscript DataProcessing/cots_sdm_v2/final/scripts/36_protocol_density.R           # protocol under density-aware classes; alternative densities and site-level layers, v2 (~30 min; needs 20b, 33, 34)
Rscript DataProcessing/cots_sdm_v2/final/scripts/37_site_ranking_density_v2_2027.R # density-aware v2 product: CSV + KMZ (needs 36); then rerun 21 to add exploration zones to its KMZ
Rscript DataProcessing/cots_sdm_v2/final/scripts/07_build_report_figures.R      # last: uses outputs of 06-21
```

The three `03` fits are independent and can run in parallel (each uses 4
cores; about 50 minutes each). `03 --smoke=1` runs a seconds-long pipeline
test whose outputs are tagged `_SMOKE` and never feed a layer. `04` caches each
stage under `cache/` and resumes if interrupted.

## Folder layout

```
final/
  config.R        frozen inputs (path + SHA-256), chain lengths, seeds
  R/              common.R, hmsc.R, support_grid.R, rasters.R, figures.R, context.R, forward.R, cull_sites.R, site_selection.R
  scripts/        00-23 (see Run order; 07 runs last)
  outputs/
    layers/       delivered GeoTIFFs + layer_manifest.csv   (tifs git-ignored)
    tables/       every number quoted in the report            (tracked)
    figures/      fig01-fig19 PNGs                              (tracked)
    site_ranking/ 2027 cull-site priority, within-reef (20), density-aware (20b) and density-aware v2 (37) (CSV, KMZ incl. exploration zones), exploration zones (CSV, GPKG)
  report/         results.qmd: Results section for the technical report
                  (quarto render report/results.qmd -> .docx and .html)
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

### Report review 04 structure (16 October 2026)

The technical report's main Methods and Results now focus on the operational
objectives with a model-purpose table and an evaluation table. Detailed model
specification is in Appendix G, tactical comparisons in Appendix H, and
metacommunity calibration/diagnostics in Appendix I. Original table and figure
data are retained; this is an editorial reorganisation, not an analysis update.
Figure 21 presents the existing manta-thinning result without model shorthand;
rebuild it with `07 --report_highlights_only=1`. The report title now reflects both
habitat suitability and metacommunity work, incorporating the author's review 03
title edits. The manuscript is unchanged.

### Metacommunity report update (21 October 2026)

Section 3.4 now uses the dated `2026-10-21_wwave_lf04` figure set from the
ADRIA sandbox. Figure 22 combines source Figures 1 (COTS) and 4 (coral), with
one reef per row. Other figures are replotted at portrait width with eight-point
labels, preserving the source summaries and bands; projection medians are
checked numerically against the frozen source table. Figure 23 shows recent diver/cull-density comparisons;
Figure 24 checks additional reefs. Figures 25–27 show exploratory effort-based
projections through 2050 and their outcomes. These supplementary figures are
in Appendix I; Appendix J preserves the superseded 15 October snapshot.
The overall title remains the wording restored from review 03; Section 3.4's
heading now reflects the COTS/coral comparison and control-strategy purpose.

The current calibration uses five base seed pairs and 65 runs per culling
state, spanning one-at-a-time parameter perturbations. Its bands are full
ranges, not confidence intervals. The working model includes updated coral,
chlorophyll and wave-settlement responses and a +0.04 initial fast-coral
adjustment at Lizard. Historical removals use reef-proportional allocation.
Projections use 450 continuous-history arms (five climate models × five seeds
× two pathways × nine culling levels). Increased effort reduces simulated
near-term peaks, but 2026 overshoot, outside supply and structural uncertainty
remain unresolved. This does not complete independent validation or establish
program benefits. Population/coral summaries in Figure 27 cover five reefs;
removals cover the full domain.

Run script **25**, then **07 --calibration_report_only=1** to rebuild without
the GBR rasters or rerunning script 22. New metadata pins in `config.R`
freeze member tables/figures and source data checksums; the builder checks
those checksums before combining the COTS/coral plots. Outputs are
`calibration_20261021_{summary,peaks}.csv`,
`projection_20261021_reef_outcomes.csv` and Figures 22–27 (PNG/PDF).
No biological fitting is performed and old pins/results are retained.
ADRIA base revision was `1d18e8bef3dca1ab33ad9b1add7548e2fffc27ea`, with local
uncommitted projection changes. SALAD/cull inputs remain in the sibling worktree;
derived figures are for local report review, with publication permissions to
confirm as recorded in the source metadata. The manuscript is unchanged.

### Review 06: hydrodynamics and larval connectivity (23 October 2026)

The technical report now has a separate Methods section with the requested
two paragraphs for Sev and Owen to complete: hydrodynamic model/particle
tracking and connectivity matrices/network analysis. A short third paragraph
describes the completed concentration-map comparison. New Results Section 3.4
reports the 2022/23 concentration patterns before the population-model results
(now 3.5). It covers 2,895 sites, 113 reefs and ten releases. Concentration does
not establish local retention, settlement or source-to-destination exchange;
the formal matrix/network analysis remains prospective.

Figure 28 is a portrait version of source lc2 (day-10 hotspot persistence).
Appendix K contains Figures 29–30, reformatted from lc4 (variation among releases)
and lc6 (cross-shelf patterns). Source figures in
`DataProcessing/connectivity_lizard/outputs/figures/` are unchanged. The main
findings are low persistent concentration around Lizard in this season,
release-to-release variation and differences among latitude bands. Domain-edge
effects remain an alternative explanation for the north–south gradient.

Rebuild with **26**, then **07 --larval_concentration_only=1**. The source
tables and all day-10 shapefile components, including the DBF containing the
concentration attributes, are pinned separately in `config.R`. Script 26
verifies that its day-10 reef concentration ranks match the source summary.
It writes `larval_concentration_*.csv` to the report tables and map data to
the ignored cache. No new hydrodynamic run or network analysis is performed.
Methods, figures, limitations and Work Plan output F distinguish this completed
descriptive work from the outstanding matrix, source-buffer and expert-method
details. Review 06 comments 26/45 are recorded in its response note; review 07
provides the next clean review copy and baselines. The manuscript is unchanged.

### Retention range and fallback surveys (26 October 2026)

Table 3, Section 3.1.1 and the detailed Appendix H comparison now include
80%, 90%, 95% and 99% retention of >0.04 dives. Script 23 adds the 80% target
to the same held-out rankings, exact site counts and 500 paired bootstrap
draws. All earlier rows in the fixed, matched and paired-difference tables
were verified unchanged. Models and the operational product were not refitted.

In the forward 2019–2026 test, the combined ranking searches 66.6% of sites
at 80% retention (95% interval 62.0–71.1), compared with 72.6% using habitat
alone and 80.4% at random. At 90%, it searches 80.1% (73.9–85.2), compared
with 85.8% and 90.5%. At 80% retention it avoids 28.7% of low-density dives
and 34.5% of ghost dives; at 90% it avoids 17.5% and 21.8%. These are initial
selection comparisons, with 20% or 10% of observed >0.04 dives omitted.
They are not predictions of annual effort or full-season detection.

Figure 20 now marks all four retention targets. Rebuild with script **23**,
then **07 --high_density_only=1**. Discussion and monitoring describe annual
full-coverage surveys, pre-closure checks and earlier evidence-triggered checks
as fallback options. They require assessment of detections, delays and all
follow-up effort. Initial-visit audits can estimate false negatives; detections
on later annual visits can include new populations. No recovery rate or
vessel-time saving has been inferred. Review 08 is the new clean review copy;
earlier annotated reviews and baselines remain unchanged.

Before the next review version, Table 3 was expanded to include the existing
`low_avoided` and `ghost_avoided` outcomes for every ranking and retention target,
alongside sites searched, with the same 95% intervals. The Results text now
summarises avoidance at all four targets. Low density remains positive CPUE
below 0.02, ghost means zero COTS and 0.02–0.04 remains a separate intermediate
class. This is a reporting amendment; model fits and result tables are unchanged.

### Review 08 integration (26 October 2026)

Applied the tracked Results edits and comments 52, 106, 110, 126, 163 and 164.
The main Results now begin with GBR habitat predictions and the four-reef map;
the map's embedded title, subtitle and footnote are removed, while its panel
labels/legends remain and imagery attribution is retained in Appendix H.
Rebuild the four-reef figures without the GBR overview with
`07 --four_reef_maps_only=1`. The numerical layers and result tables are unchanged.

The standalone main manta-thinning section and simplified figure were removed;
a short benchmark interpretation remains in the preceding section and detailed
comparisons remain in Appendix H. The practical emphasis is habitat plus recent
observations and crew judgement, without asserting that the cull-site model
beat the tow-thinning benchmark. The product and metacommunity prose follows the
author's shorter wording. Projections demonstrate model capability; the recent
population overshoot and possible culling bias remain to be tested.

Connectivity Results now state the proposed early-warning application of larval
concentration and network analysis in the initiation box. Management implications
retain five subsections but use eight short bullet paragraphs, about 450 source
words versus about 2,400 previously. Full management interpretation, scenario
design and field-trial requirements are preserved in new Appendix L. The existing
author-completion checks remain there and in their relevant Methods/appendices.
Review 08 and its baselines are preserved; review 09 is the next clean copy,
with a response note mapping all comments and tracked-edit themes. The manuscript
is unchanged.

## Cull-site and manta-effort clarification — 27 October 2026

“Sites searched” in the >0.04 assessment counts candidate designated cull
site-years selected for a dive visit, including never-dived sites. The report
now calls this **sites selected for cull dives**. It is not a manta-tow count,
nor a count of sites where COTS were actually removed.

New script **27_operational_cutoff_and_tow_effort.R** replays all ten forward
score columns with the original tie-breaking stream. Every archived dive rank
must match exactly before proceeding. The source coefficients are now explicitly
pinned in `config.R`; there is no model refit or change to the deployed product.

The product's A–C rule selects descending rank/site count ≤0.8, or
floor(0.8 × site count) sites per reef. It selects 78.9% of candidate site-years
and retains 89.1% of observed >0.04 dives over 2019–2026, avoiding 18.8% of
low-density dives and 23.5% of ghost dives. Thus deferring D initially misses
10.9% of observed high-density dives. These are new exact-product-cutoff
results, with uncertainty and same-effort baselines in
`outputs/tables/operational_priority_cutoff.csv`. They are distinct from the
older percentile-cutoff results in `cull_high_density_priority_classes.csv`.
The four existing matched-retention point estimates remain unchanged.

Script 27 also overlays valid recorded start-to-end tow segments on designated
sites. It calculates potentially avoidable **whole tows** and **straight-line
tow distance** at the same retention targets and at the A–C cutoff. Surveillance
outside site polygons, and in portions shared with any retained site, remains.
Source fallback or zero-length geometries are excluded and coverage is reported.
Exact random-site expectations account for overlapping polygons, without double
counting distance. Uncertainty uses 500 paired reef-year bootstrap draws.

At A–C, the conditional gross reductions are **11.6% of whole manta tows**
(95% interval 11.1–12.1) and **16.2% of straight-line distance** (15.6–16.9).
At 80%, 90%, 95% and 99% dive retention the combined ranking's whole-tow
reductions are 19.2%, 10.9%, 5.1% and 0.9%; distance reductions are 25.5%,
15.3%, 7.6% and 1.4%. The tow denominator is 125,375 valid segments on 774
of the 861 evaluated reef-years; 12 records lacked usable segment geometry.
All 240 comparable fixed-effort and 48 matched-retention point estimates
reproduce the existing dive/site results; reductions decrease to zero at full
site coverage. At the same A–C site count, random selection retains only
78.5% of >0.04 dives, versus 89.1% with the combined ranking. Raw reductions
in site or tow effort alone should not be attributed to ranking skill.

These are conditional spatial effort proxies: high-density dive retention
assumes selected cull visits still occur. They do not establish detection after
removing manta surveillance, feasible replacement routes, vessel-time savings,
annual savings after deferred-site checks, or exploration yield. A field trial
must assess those outcomes, including annual/pre-closure coverage.

New tables are `operational_priority_cutoff.csv`,
`operational_tow_effort_fixed.csv`, `operational_tow_effort_matched.csv` and
`operational_tow_effort_coverage.csv`. Run **27**, then
**07 --high_density_only=1**, then render the technical report. The geometry
cache is keyed to frozen geometry inputs, the evaluation reef-years and an
explicit geometry algorithm version. Bump that version when changing geometry
logic. Figure 20 marks A–C explicitly; Figure 17's report caption links its
classes to the Table 3 trade-off and the separate manta-effort table.

This amendment updates the main report before the next review iteration.
After the author confirmed it had no comments, review 09 and both clean
baselines were refreshed from the updated main report on 27 October. Previous
copies are archived in `cache/review_09_before_effort_refresh_20261027_093641/`.
Review 09 is now ready for feedback; no review 10 has been created.

## CCIP Bridging Project report template — 27 October 2026

The technical report and review 09 use the supplied
`report/CCIP Final Report Template - Bridging Projects.docx`. The original is
never overwritten. Formatting and section-content requirements are governed by
[`report/CCIP_TEMPLATE_GUIDE.md`](report/CCIP_TEMPLATE_GUIDE.md), linked from
`WRITING_BRIEF.md` for future Codex and Claude revisions. This is the scientific
technical-report template, distinct from the administrative final-project form.

The main report now follows the template's section order, one-page Project
Summary, dedicated rights-holder/stakeholder engagement subsection, explicit
Outcomes & Impact Pathway, ≤300-word conclusion, Data Accessibility statement
and Coral Reefs reference style. Main Results contains four key displays; the
full matched-retention tables, manta-effort comparisons and retention curve
are retained in Appendix H. No model was refitted or analytical result changed.
The new observation-footprint map `fig29_report_study_area.png` uses the frozen
cull archive and calibration-reef locations, and is rebuilt through script 07
with `--ccip_template_only=1`.

From the repository root, export with:

```powershell
& .\DataProcessing\cots_sdm_v2\final\report\scripts\render_ccip_report.ps1
```

This renders Quarto content, assembles a working copy of the source template,
refreshes native Word contents/caption/page fields, exports the matching PDF,
checks source-template identity and critical formatting, and updates the main
exports. Create the next numbered review and matching Word/Quarto baselines
only after verification; preserve annotated review inputs. The legacy
`-RefreshReview09` option must not overwrite the now-annotated review 09.
Existing exports are archived under `cache/ccip_template/build_*/`.
Inspect the page images before handing over a new review. A bare Quarto render
is only the content intermediate and does not assemble the final template.

The next feedback file is `report/reviews/technical_report_review_11.docx`,
with matching clean baselines. Annotated review 09 and its original baselines
are preserved. See `report/reviews/technical_report_review_09_responses.md`
for the integration record. The older `_CCIP` copies are historical exports,
not the current report. Close a target review before replacing its on-disk copy.

The template's cover, publication information, standard acknowledgements,
disclaimer, teal/cyan styles, Roman/Arabic numbering, footer, logos and final
gradient page are reused. Contact, affiliations, citation-author and copyright
slots are held in `report/ccip_report_metadata.json`; title/authors remain in the
Quarto YAML. The original template checksum is
`e88e42c58b01be83287e0848dc4ff68b3f5f4f0e3ccc27ba08793424f685a1f5`.

Completed engagement now includes the author-confirmed April visit and June
follow-up with Malu and multiple GBRMPA meetings, including development of the
cull-site-ranking concept. The extra observations support a separate future Far
Northern SDM. The 5 November 2026 delivery workshop and January–February 2027
trials aboard BPM and with Malu remain planned. Authors/contact/copyright,
public release links and the outstanding hydrodynamic specification still need
author completion before final publication.

## Reading the A–C trade-off: sites, dives and denominators — 27 October 2026

Script **28_explain_operational_dive_tradeoff.R** adds a descriptive check of
Table 1 using the frozen forward ranks. Run it after script 27. It writes
`outputs/tables/operational_deferred_dive_breakdown.csv` and reproduces the
published high-density, low-density and ghost-dive point estimates exactly.
It also writes `operational_deferred_dive_metrics.csv`, with 500 paired
reef-year bootstrap intervals for all-dive avoidance, intermediate retention
and deferred-dive composition. These additions are now integrated into Table 1
and its accompanying explanation; no existing ranking or result was refitted.

The A–C rule keeps approximately the top **80%**, not the top 20%, of candidate
sites. Its 78.9% candidate-site coverage retains 85.1% of all **recorded dives**:
4,248 of 28,431 observed dives are deferred (14.9%). Random selection at the
same candidate-site count defers an expected 6,139.6 dives (21.6%). Thus the
comparison matches candidate-site coverage, not observed dive count or time.
The ranking retains more of the sites that crews visited repeatedly; this
reflects operator behaviour as well as the evidence used to rank sites.

Of the model-deferred dives, 28.5% are ghosts, 26.7% have positive CPUE below
0.02, 11.9% have CPUE from 0.02 through 0.04 inclusive, and 32.9% are above
0.04. The corresponding random-site expected composition is 18.4%, 21.3%,
15.4% and 44.9%. Ghosts plus low-density dives therefore form 55.2% of the
deferred set, versus 39.7% under random site selection. These are shares of
**deferred dives**, not the shares of each CPUE class avoided in Table 1.

The missing intermediate class is also mostly retained: 11.6% is deferred,
versus 21.5% under random site selection. Adding it completes the accounting;
it does not create an additional efficiency gain. Table 1 demonstrates better
protection of observed high-density dives at the same candidate-site count,
but little extra avoidance of observed low-density/ghost dives. It does not
establish a large reduction in dive effort or search time. A comparison at
matched observed dive effort, and eventually actual operational time, is needed
to assess that claim separately. The revised report distinguishes candidate-site
coverage from observed-dive coverage and includes all four CPUE classes with a
clearly labelled deferred-dive composition. Retaining the intermediate group
preserves useful culling opportunities and is desirable to the program.

## Review 09 integration and ranking explanation — 27 October 2026

The main report incorporates all six tracked passages and both comments in
review 09. Summary and Results elevate low-density (18.8%) and ghost-dive
avoidance (23.5%); intermediate-dive retention is 88.4% [86.7–90.1].
The summary also reflects hydrodynamic/connectivity work, ADRIA integration,
the fifth-wave application, repeated GBRMPA engagement and exploration areas.
The possibility of additional tow reductions outside polygons remains untested.

Methods now gives the actual logistic score used for operational ranking,
with plain-language input definitions. Appendix G adds the Gaussian weight
(300 m scale, truncated at 900 m, at most 400 nearby records per stream),
no-evidence and missing-habitat handling, the archived 2027 coefficients,
the evidence window and within-reef A–D assignment. This is the ≥0.02-trained
site-evidence model, not the separate HMSC history blend; the >0.04 evaluation
still uses the existing scores. Two editable display equations reuse the
template's borderless equation-table component and live equation numbers.

Figure 1 adds the actual 2,895-site v0.4 domain from the population-model
GeoPackage, pinned in `config.R`. The archive now resides in the adjacent
`COTSMod-Calibration/store/data/` tree. The hydrodynamic mesh panel requested
in comment 65 is explicitly reserved for the figure to be supplied; it is not
fabricated or substituted with a different geometry. Rebuild Figure 1 with
script 07 `--ccip_template_only=1`.

Review 09 is preserved as the annotated source. Review 10 and clean baselines
are the next review round; the template source is unchanged.

## Three ranking alternatives and habitat percentiles — 27 October 2026

The incumbent remains the default. Scripts **29–32** evaluate three fixed
alternatives specified in
[`report/ranking_alternatives_analysis_plan.md`](report/ranking_alternatives_analysis_plan.md):

1. Replace the prior-site ≥0.02 indicator with any prior dive **>0.04**.
2. Add site CPUE as recency-weighted catch / recency-weighted bottom time,
   plus the effective recent-record count, using a fixed **90-day half-life**.
3. Use the **within-site 90th-percentile habitat prediction** instead of its
   mean in the recent-CPUE model.

All retain the ≥0.02 training response and the original quiet-site, nearby
dive and manta terms. Intermediate dives remain useful culling opportunities.
Evidence is strictly preceding-season; habitat retains the existing spatial
holdout and temporal limitation. Neither frozen inputs nor deployed products
were replaced.

The original habitat means used 16 sample points per site in validation.
The new comparison uses 2,201,095 polygon-cell/fallback records, median 139 per
site. Five refitted habitat folds reproduce archived sample-point means with
**zero error**. A current-model control using the same all-cell mean separates
denser sampling from changing the predictor. Twenty-one sites have no usable
predictor values and are flagged. Median, 75th, 90th and 95th percentile habitat
summaries are additional habitat-only diagnostics, not post-hoc combined models.

The forward cohort is unchanged: **28,431 dives, 861 reef-years and 32,444
candidate site-years**. Exact incumbent A–C point estimates replay unchanged.
New trade-off curves use 1% integer-site steps and first-crossing interpolation,
rather than the older 5% percentile-cutoff steps. Thus interpolated incumbent
values can differ slightly from Appendix H; its exact A–C point is identical.
New intervals use 500 common paired reef-year bootstrap draws.

| Ranking | >0.04 retained at A–C (%) | Dives avoided at 90% >0.04 retention (%) | Same measure, 2024–2026 (%) |
|---|---:|---:|---:|
| Current ranking | 89.1 | 13.8 | 14.8 |
| All-cell mean control | 89.4 | 13.9 | 15.1 |
| Prior >0.04 indicator | 89.5 | 13.9 | 14.9 |
| Indicator + recent CPUE | 90.5 | 14.7 | 16.6 |
| Recent CPUE + 90th-percentile habitat | 90.0 | 14.3 | 13.3 |

These are point estimates; full intervals and paired differences are in
`outputs/tables/ranking_candidate_assessment.csv` and
`ranking_candidate_differences.csv`. Recent CPUE has the best positive overall
and recent comparison. At 90% high-density retention its extra recorded-dive
avoidance is **0.89 points [−0.08, 1.97]** versus the incumbent and **0.77
[−0.11, 1.83]** versus the all-cell control. Recorded bottom-time avoidance is
14.8% versus 13.7%, a paired gain of **1.18 points [−0.01, 2.33]**. These
intervals include zero. At the same actual recorded-dive budget its high-density
retention gain is smaller: **0.54 points [−0.05, 1.31]**. At matched 90%
retention it also retains fewer intermediate dives: 88.1% versus 89.6%.

No alternative passes the pre-specified strong replacement gate. The recent-CPUE
option is exported **experimentally only**, not promoted. Gains vary among
seasons and are positive in five of eight years. The 90th-percentile combined
version is worse recently than the mean version. Habitat-only useful-dive AUC
is 0.576 for the mean, 0.589 for the median, 0.575 for the 75th percentile,
0.568 for the 90th and 0.560 for the 95th. This does not support favouring the
highest-valued habitat patches in this ranking test; it does not show that real
hotspots are absent. Median habitat is a potential separately specified next test.

### Experimental product

`outputs/site_ranking_candidates/` contains
`cots_cull_site_priority_2027_recent_cpue_experimental.csv`, `.gpkg` and `.kmz`,
plus a plain-language README. All **15,884** sites are included; KMZ matching
is 100%. The evidence window remains **1 June 2025–31 May 2026**. The model is
fitted on the existing deployment years 2019–2026. The probabilities concern
≥0.02 outcomes conditional on diving, not >0.04 outcomes or unvisited-site
occurrence. The coefficient export explicitly labels the old internal
`site_problematic_prev` alias as the new `site_high_prev` (>0.04) input.

The candidate changes which sites fall in class D at 122 locations, while
preserving each reef's class counts. Incumbent CSV/KMZ checksums are unchanged.
The candidate manifest records its checksums and experimental status. Use it
for comparison and field trials; annual checks, exploration, manta detection
and actual vessel-time savings still require measurement. Candidate selection
is exploratory and intervals are not adjusted for choosing among three models.

### Reproduction and report

Run **29 → 30 → 31 → 32 → 07 `--ranking_candidates_only=1`**, then the native
report export. Script 29 caches distinct cell reads and verifies a direct
extraction at five distributed polygons. Inputs and the deployed habitat layer
are pinned in `config.R`. Script 30 replays every archived forward rank before
fitting, checks unchanged A–C point estimates, and records annual results.
Script 31 exports only a separately labelled comparison product; script 32
checks geometry/candidate coverage, date/support metadata and product changes.

Figures **31** (site and recorded-dive trade-offs) and **32** (habitat-summary
AUC) are built by script 07, using `R/ranking_candidate_figures.R`. Full methods,
three comparison tables and both figures are in new Appendix M. Main Summary,
Methods, Results, Discussion and the product list identify the option while
preserving existing Table 1/Appendix H results. Review 11 is the new clean copy;
review 10 and all earlier feedback remain preserved. No commit or push was made.

## Deployment-realistic validation (33) — 30 October 2026

**Question.** Do the leave-reef-out statistics undersell the ranking? A
deployed model has seen every reef the program works on, and crews plan a
voyage knowing everything up to the day it leaves. Leave-reef-out answers "how
does this work on a reef the model has never seen"; the deployment question is
"how well does it rank the next visit, given all earlier data".

**Method.** Strictly forward in time, 2019–2026:
- For each test season T, every model is fitted only on data from before T, but
  from **all reefs**, so the test reef's own earlier seasons are in training.
- All rankings are scored on the same 28,431 dives, grouped by reef visit
  (reef × voyage; 2,976 visits in 861 reef-years).
- **Habitat:** leave-reef-out rank (15); forward rank (XGBoost refitted on
  site-years before T, all reefs); and the deployed layer, which has seen
  season T's own dives (a leaked ceiling, shown only to size the leak).
- **Evidence** (the site-evidence logistic model, fitted on dives 2013..T−1):
  previous calendar season (as deployed), or the **365 days before the voyage
  left port**. The second counts earlier voyages in the same season, but not
  the voyage's own manta survey or dives.
- 300 paired reef-year bootstraps (`deployment_validation*.csv`).

**Results, 2019–2026, at 90 % of > 0.04 dives kept:**

| Ranking | Sites selected | Low-density avoided | Ghost avoided | > 0.04 kept at A–C (top 80 %) | AUC, > 0.04 vs rest |
|---|---:|---:|---:|---:|---:|
| Habitat, leave-reef-out | 88.7 % | 9.6 % | 11.9 % | 84.2 % | 0.554 |
| Habitat, forward (sees own reef's past) | 82.5 % | 13.6 % | 17.0 % | 88.0 % | 0.556 |
| Habitat, in-sample (leaked ceiling) | 84.7 % | 13.7 % | 16.6 % | 86.6 % | 0.582 |
| **Current:** leave-reef-out habitat + previous season | 81.3 % | 16.5 % | 21.6 % | 89.0 % | 0.560 |
| Forward habitat + previous season | 82.2 % | 15.6 % | 20.0 % | 88.7 % | 0.563 |
| Leave-reef-out habitat + evidence to voyage start | 73.3 % | 20.4 % | 29.6 % | 94.0 % | 0.605 |
| **Forward habitat + evidence to voyage start** | **72.8 %** | **20.2 %** | **28.7 %** | **93.7 %** | **0.609** |

Paired against the current validation, forward habitat + evidence to voyage
start selects **8.4 points fewer sites** (95 % CI 2.0–15.5). It keeps **4.7
points more > 0.04 dives** at the A–C cutoff (2.5–7.3) and avoids 7.1 points
more ghost dives (0.5–13.6). The low-density gain (+3.7) has an interval that
includes zero. Season by season it has the higher AUC in all eight seasons,
and selects fewer sites in seven (2023 is the exception: 76 % vs 75 %; range
53–77 % vs 68–86 %).

**What this means:**
- **Seeing the reef's own history helps the habitat model on its own** (88.7 %
  → 82.5 % of sites). Once last season's evidence is in, it adds nothing
  (81.3 % vs 82.2 %), because the evidence terms already carry that history.
  The leave-reef-out statistics undersold the habitat layer alone, not the
  combined ranking.
- **The in-sample layer is no better than the forward one.** Letting the
  habitat model see the test season's own dives does not buy accuracy here.
  Same-season leakage is not inflating the deployed habitat layer much.
- **The real gain is fresher evidence.** Updating the evidence to each
  voyage's start, so earlier voyages in the same season count, is what moves
  the result. It is legitimate: crews know those outcomes when they plan.
- **Manta tows:** the forward within reef-year manta model ranks this season's
  tows the same as leave-reef-out (AUC 0.615 vs 0.611; lift 1.46 vs 1.42;
  difference CI spans 0).
- **Deployment implication:** re-rank before each voyage, with the evidence
  window ending at departure, rather than once a year. Script 20 currently
  fixes `evidence_end` for the season. The product is unchanged by this
  analysis.
- **Not tested:** using the voyage's own manta survey before its dives. That
  is current practice (manta triggers cull) and would score well, but it is
  the loop the program wants to reduce.
- **Note on numbers:** grouping by visit and a 2 % fraction grid mean the
  current-ranking values differ slightly from 19, 23 and 27 (for example
  81.3 % vs 80.1 % of sites at 90 % > 0.04 retention).
- **Fair baseline caveat.** These savings are measured against random site
  selection. Crews are not random: they dive open sites (last > 0.04) first
  and then follow the manta tows. The prior > 0.04 indicator alone scores
  about the same as the full ranking (29–32), and the gain from fresher
  evidence largely reproduces what crews already know when they leave port.
  Script 34 tests the ranking against current practice instead.

## Option A: tow only the priority zone, then let the tows guide (34) — 30 October 2026

**Question.** Under a realistic workflow, what does the ranking add beyond
what crews already do? Open sites are dived first under any system. The
remaining choice is where to tow. Option A tows only within a priority zone
(the top 80 % or 50 % of sites, plus open sites) and lets the tows guide the dives, as now. The alternative,
diving A sites without towing, is what 16–33 score; it cannot be checked
fairly, because recorded dives were themselves directed by tows.

**Rules (applied identically to every strategy):**
- **Open site:** the last dive in the 365 days before planning was > 0.04.
  Dived first under every strategy, so excluded from outcomes and always in
  the zone.
- **Guiding tows:** tows on the reef from 30 days before the visit's first
  dive to its last dive, in recorded order. The export row order is the
  towing order: consecutive tows are a median 421 m apart, against 1.9 km for
  a random order. A segment breaks on a new day or a > 1 km jump.
- **Own tows (corrected 2 November):** each tow belongs to the earliest visit
  whose window holds it. Follow-up visits come a median 15 days after the
  previous visit, so the 30-day window pulled in the earlier visit's sweep:
  75,425 of 180,414 window tows were counted twice. A follow-up visit now
  keeps those tows as evidence (they were towed then) but can skip, and
  counts as effort, only its own tows: 104,989 distinct tows.
- **Planning date:** the zone is set before the visit's own survey. Evidence
  runs to the earlier of departure and the first own tow, using 33's season-T
  coefficients (forward habitat + evidence). Before the fix, inherited tows
  pulled the planning date back before the earlier visit's dives, so sites
  that visit found above 0.04 were not treated as open.
- **Zone:** sites in the top `cut` of the reef plus open sites. A tow is in
  the zone if it passes within 100 m of a zone site.
- **Feasible skipping:** only runs of three or more consecutive out-of-zone
  tows are skipped. Shorter gaps are towed on the way past.
- **Trigger:** crews dive designated sites within 100 m of a kept
  strong-evidence tow. **Neighbour fill** also dives sites near the skipped
  tow next to a kept strong tow.
- **Recent-tow variant:** in-zone tows that repeat ground towed in the last
  365 days (midpoints within 100 m) may also be skipped. A skipped tow whose
  earlier tow was strong still triggers.
- **Baselines at the same tows saved:**
  - *Shifted:* the same skip pattern rotated to random positions in the
    visit's tow sequence (20 offsets).
  - *Block thinning:* k tows on and m off, averaged over phases.
  - Single-tow systematic thinning (14) is dropped: it is not feasible, and
    it gained from counting an aggregation as found from one tow.
- **Outcomes:** the visit's actual dives at non-open sites, relative to full
  towing; 300 paired reef-year bootstraps (`manta_zone_*.csv`).

**Coverage (`manta_zone_coverage.csv`):**
- **Sample:** 1,876 visits (726 reef-years), 104,989 distinct tows; 96 % pass
  within 100 m of a designated site.
- **> 0.04 dives:** of 7,920, **61 % are at open sites**. Of the rest, 59 %
  lie within 100 m of a strong guiding tow; the others came from crew
  judgement or other information.

**Results, 2019–2026, all visits, direct trigger, fresh tows (full towing = 100):**

| Strategy | Tows saved | > 0.04 dives found | Aggregations seen | Low-density incurred | Ghost incurred |
|---|---:|---:|---:|---:|---:|
| Tow near designated sites only (cut 100 %) | 1.7 % | 100.0 | 99.8 | 100.0 | 100.0 |
| Zone top 80 %, combined ranking | 7.8 % | 95.9 | 98.2 | 97.9 | 93.3 |
| Same skips, random place | 7.8 % | 95.8 | 97.2 | 96.2 | 93.5 |
| Zone top 80 %, habitat only | 8.1 % | 96.4 | 97.8 | 95.5 | 92.9 |
| **Zone top 50 %, combined** | **23.1 %** | **91.7** | **93.7** | **89.3** | **85.3** |
| Same skips, random place | 23.1 % | 89.1 | 91.5 | 86.3 | 83.0 |
| Blocks 9 on 3 off | 25.0 % | 92.5 | 91.3 | 86.6 | 85.4 |
| Zone, recent tows reused (cut 100 %) | 39.0 % | 98.0 | 81.1 (fresh) | 105.5 | 106.5 |

Zone minus the same skips at random (paired, `manta_zone_strategy_differences.csv`):
top 50 % combined, > 0.04 found +2.6 (−0.5 to 5.3), aggregations +2.1
(1.3–3.1); top 80 % combined, +0.1 (−2.3 to 2.0) and +1.0 (0.6–1.5).

**Follow-up vs first visit of the season (`manta_zone_absolute.csv`):**

| Visits | Strategy | Tows saved | > 0.04 dives not triggered | Ghost avoided | Aggregations seen |
|---|---|---:|---:|---:|---:|
| Follow-up (1,256; 55,452 own tows) | **Top-50 % zone** | **11,734 (21 %)** | **41 of 1,155** | **20 of 244** | **96.0 %** |
| | Same skips, random place | 11,734 | 82 | 24 | 94.0 % |
| | Blocks 9 on 3 off | 13,863 (25 %) | 68 | 26 | 93.1 % |
| First of season (620; 49,537 tows) | Top-50 % zone | 12,526 (25 %) | 108 of 646 | 48 | 84.4 % |
| | Blocks 9 on 3 off | 12,384 | 68 | 42 | 83.9 % |

**What this means:**
- **Use the zone only on follow-up visits.** There it halves the misses of
  random placement and beats feasible block thinning; on first visits blocks
  are better. Sweep first, then use the zone.
- **The zone saves towing, not ghost dives:** only 464 ghost dives over eight
  seasons were near a strong guiding tow at all.
- **Reusing recent tows** still saves about 40 % of tows for 98 % of > 0.04
  dives, but sees only 81 % of aggregations fresh and adds 5–7 % ghost dives.
  It needs a field trial before use.
- **The earlier (double-counted) version of this section** reported 22.9 %
  tows saved and 118 > 0.04 dives missed on follow-ups. Both came mostly from
  re-skipping the earlier visit's sweep. The direction of every comparison is
  unchanged.
- **Caveats:** tows saved counts tows, not boat time; a follow-up that
  inherits tows from an earlier follow-up keeps them even if that visit
  skipped them (slightly optimistic); visits without tows in the window
  (1,154) are excluded.

**Where ghost dives come from (`manta_zone_dive_pools.csv`;** every dive of the
1,876 visits by why it could have been chosen, scored with 33's ranking at
voyage start):

| Pool | Dives | > 0.04 | Ghost | Ghost rate | AUC > 0.04 vs ghost |
|---|---:|---:|---:|---:|---:|
| Open site | 7,273 | 4,817 | 313 | 4 % | 0.59 |
| Strong tow within 100 m | 3,752 | 1,769 | 456 | 12 % | 0.60 |
| **Weak-scar tow only** | **5,393** | **872** | **1,952** | **36 %** | **0.65** |
| Negative tows only | 1,037 | 163 | 460 | 44 % | 0.61 |
| No guiding tow within 100 m | 612 | 230 | 127 | 21 % | 0.65 |

59 % of ghost dives follow a weak-scar tow. The dive decision after weak or
no evidence, not the tow zone, is the lever for ghost dives; 35 scores it.

## Proposed tow and cull dive protocol (35) — 2 November 2026

Written up for GBRMPA as a protocol proposal (Claude Doc, 2 November 2026).

**Protocol scored** (recorded dives of 2019–2026 visits, 34's rules):
1. Dive open sites first.
2. Follow-up visits tow only the zone: A–C sites (top 66 %; A–B, top 50 %,
   until 6 November 2026), open sites and,
   with **season memory**, sites within 100 m of a strong tow on an earlier
   visit that season. First visits tow everything.
3. Dive sites within 100 m of a kept strong tow.
4. **Dive rule:** any other dive (weak scar, only negative tows, or no tow
   nearby, e.g. an adjacent site or long-unvisited ground) only at sites
   ranked A–C (top 66 %) at voyage start, at crew discretion.

A dive whose only strong tows were skipped counts as missed. Ranking classes
changed with this work (`config.R`, script 20): A 20 %, B 50 %, **C 66 %**,
**D 80 %**, E the rest. The old A–C (top 80 %) is now A–D.

**Results** (`protocol_statistics.csv`; 18,182 dives, 104,989 tows; 95 %
reef-year bootstrap):

| Option | Tows saved | Ghost avoided | Low-density avoided | > 0.04 missed | CPUE of missed > 0.04, mean (CI) | COTS not culled | Ghost per > 0.04 missed |
|---|---:|---:|---:|---:|---|---:|---:|
| Dive rule A–D | 0 | 18 % (611) | 12 % (475) | 1.4 % (111) | 0.097 (0.077–0.121) | 1.3 % | 5.5 |
| **Dive rule A–C** | 0 | **30 % (991)** | **19 % (776)** | **2.8 % (223)** | **0.099 (0.085–0.115)** | **2.6 %** | **4.4** |
| Same number at random | 0 | 24 % (806) | 18 % (738) | 5.1 % (401) | 0.128 (0.099–0.165) | — | 2.0 |
| Dive rule A–B | 0 | 41 % (1,370) | 27 % (1,108) | 4.5 % (359) | 0.097 (0.086–0.108) | 4.0 % | 3.8 |
| Evidence rule: skip weak-scar dives | 0 | 59 % (1,959) | 46 % | 11 % (874) | 0.102 | 9.9 % | 2.2 |
| Evidence rule: strong evidence only | 0 | 77 % (2,546) | 57 % | 16 % (1,267) | 0.128 | 14.5 % | 2.0 |
| Zone (A–C, top 66 %), follow-ups | 7.2 % | 0.3 % | 0.5 % | 0.4 % (29) | 0.75 (0.11–1.10) | 0.9 % | 0.3 |
| Zone with season memory | 6.7 % (13 % of follow-up tows) | 0.2 % | 0.3 % | 0.2 % (13) | 0.41 (0.08–1.06) | 0.4 % | 0.6 |
| **Zone with memory + dive rule A–C** | **6.7 %** | **30 % (999)** | **19 % (788)** | **3.0 % (236)** | **0.116 (0.089–0.150)** | **3.0 %** | **4.2** |

Zone rows use the A–C zone (6 November 2026). The earlier A–B zone (top
50 %) saved 11 % of all tows (20 % of follow-up tows with season memory) and
missed 41 > 0.04 dives (22 with season memory); with the dive rule it missed
245. Widening the zone to A–C gives up about 7 points of follow-up tow
savings to halve the zone's misses again.

- **The dive rule carries the protocol.** A–C is the middle setting: from
  A–D to A–C each extra > 0.04 dive missed avoids 3.4 ghost dives; from A–C
  to A–B, 2.8 (random 2.0). Missed > 0.04 dives are marginal (median CPUE
  0.069, IQR 0.052–0.113; all > 0.04 dives: mean 0.187, median 0.10).
- **No-tow dives are productive.** Skipping every dive with no tow nearby
  would lose 230 > 0.04 dives to avoid 127 ghost dives (38 % of them were
  > 0.04). The dive rule therefore leaves them to crew discretion within A–C.
- **Zone misses come from one site.** Of the 29 > 0.04 dives the plain A–C
  zone missed (`protocol_zone_missed_dives.csv`, `_reefs.csv`), 14 were at
  Horseshoe Reef (22-104a) site 11 in April 2019 (CPUE 0.38–2.7, 52 % of the
  missed COTS). A strong tow on the season's first visit 70 days earlier had
  flagged it, but the ranking placed it bottom. Season memory removes them
  (29 → 13) for 0.8 points of follow-up tows saved. Trunk Reef site 126
  (2019, CPUE 1.7–1.8, 24 % of missed COTS) had no earlier flag.
- **Field trial (`protocol_voyage_summary.csv`, `protocol_trial_power.csv`;
  2024–2026 voyages):** a median voyage has 3 reefs, 2 follow-ups, 152 own
  tows, 25 dives and 2 dives the A–C rule would skip (0–8). One voyage shows
  the rule's effect (lower > 0.04 share among skipped dives) only 29 % of the
  time; about 25 voyages give 95 %. Proposed: one live voyage (zone-first
  towing, flag-and-defer dives) for feasibility, then a season of shadow
  recording with the ranking frozen before each voyage.
- **Caveats:** scored on dives crews made; missed counted cautiously; COTS
  not culled remain until the next visit; tows, not boat time.

34 now also scores a 0.66 zone cut (A–C), so the report's random-position and
block baselines exist for it: on follow-up visits the A–C zone (no memory)
saves 13.5 % of tows and misses 29 > 0.04 dives; the same tows skipped at
random positions miss 49; 9-on-3-off blocks save 25 % and miss 68.

## Protocol under density-aware classes; A–B vs A–C zone (36) — 9 November 2026

`scripts/36_protocol_density.R` replays 35's protocol with the zone and dive
rule taken from either the within-reef classes (20) or the density-aware
classes (20b), for 2020–2026 (20b's forward test starts in 2020; 1,683 visits,
1,119 follow-up, 14,755 dives). Each visit keeps its voyage-start order
(33, `fwd_voyage`); density-aware class shares are
`plogis(qlogis(break) + delta)`, with delta from the reef's dives and tows in
the 365 days before the visit's planning date (earlier visits that season
count) and the forecast for the year before, using 20b's model and
unknown-reef discount refitted on seasons before T. Dive-rule evidence
variants: **discretion** (35's rule: weak, negative-only and no-tow dives only
at A–x), **weak-scar only** (only weak-scar dives filtered), **weak scar
needed** (negative-only and no-tow dives skipped everywhere, weak-scar dives
only at A–x). 300 reef-year bootstraps; paired density − within differences.
20b's helpers moved to `R/density_classes.R`; 20b caches its training inputs
for 36 (outputs unchanged, verified byte-identical).

Outputs: `protocol_density_statistics.csv` (every rule set, both scopes, CIs),
`_differences.csv` (paired), `_tradeoff.csv` (compact), `_zone_size.csv`,
`_totals.csv`; fig33.

| Classes | Setting | Follow-up tows saved | > 0.04 missed | Ghost avoided | Low avoided | Ghost per miss |
|---|---|---:|---:|---:|---:|---:|
| within-reef | Zone A–B (memory) | 19.7 % | 20 | 0.6 % | 0.6 % | 0.9 |
| within-reef | Zone A–C (memory) | 12.7 % | 11 | 0.3 % | 0.3 % | 0.7 |
| density-aware | Zone A–B (memory) | 14.3 % | 4 | 0.3 % | 0.2 % | 2.5 |
| density-aware | Zone A–C (memory) | 8.5 % | 0 | 0.1 % | 0.0 % | — |
| within-reef | Dive rule A–C, discretion | 0 | 206 (3.7 %) | 30.1 % | 19.5 % | 4.3 |
| density-aware | Dive rule A–C, discretion | 0 | 166 (3.0 %) | 33.5 % | 19.1 % | 6.0 |
| within-reef | Dive rule A–B, discretion | 0 | 321 (5.8 %) | 41.6 % | 27.8 % | 3.9 |
| density-aware | Dive rule A–B, discretion | 0 | 232 (4.2 %) | 44.4 % | 25.8 % | 5.7 |
| within-reef | Dive rule A–C, weak scar needed | 0 | 455 (8.2 %) | 40.6 % | 28.5 % | 2.7 |
| density-aware | Dive rule A–C, weak scar needed | 0 | 423 (7.6 %) | 44.6 % | 28.7 % | 3.1 |
| within-reef | Zone A–C + dive rule A–C | 12.7 % | 217 (3.9 %) | 30.3 % | 19.8 % | 4.2 |
| density-aware | Zone A–B + dive rule A–C | 14.3 % | 170 (3.1 %) | 33.8 % | 19.4 % | 5.9 |

Random skipping within the same pool gives about 2.1 ghost per miss.

- **Zone.** Reefs revisited within a season are usually dense (median delta
  0.69 on follow-ups), so density-aware zones are larger: the A–B zone holds a
  median 68 % of a follow-up reef's sites (within-reef A–C 66 %). The
  density-aware A–B zone saves about as many follow-up tows as within-reef A–C
  (14.3 % vs 12.7 %) and misses 4 > 0.04 dives instead of 11. Paired
  differences at the same cut: tows saved −4 to −5 points, misses −0.2 to
  −0.4 points (CIs exclude zero). With density-aware classes, A–B is the zone
  to use; A–C saves only 8.5 %.
- **Dive rule.** At A–C, density-aware classes skip the same share of dives
  (13.9 %) and avoid more ghost dives (+3.4 points, CI 0.9–5.5), with fewer
  misses (−0.7 points, CI −1.4 to +0.2, not significant alone; −0.9, CI −1.6
  to −0.0, with the zone). Low-density avoidance is unchanged.
- **Requiring at least a weak scar** (skipping no-evidence dives) adds misses
  almost at the random rate (about 3 ghost per miss): no-tow dives are
  productive (as in 35). Filtering weak-scar dives only is a milder version of
  the discretion rule with a similar ghost-per-miss ratio.
- **Candidate protocol:** density-aware A–B zone with season memory plus the
  A–C dive rule at discretion: 14.3 % of follow-up tows saved, 170 > 0.04
  dives missed (3.1 %), 33.8 % ghost and 19.4 % low-density dives avoided.
  The within-reef A–C version: 12.7 %, 217 (3.9 %), 30.3 %, 19.8 %.
- **Caveats:** 2020's typical density came from one training season; the
  density at planning uses earlier visits' dives, which encode where crews
  chose to dive; hindcast forecasts are probably in-sample.

## Density-aware v2: recent density and each site's own record (36, 37) — 9 November 2026

**Why.** Field checks of the 20b classes:
- Big Broadhurst scored 0.44. Its 12-month window held 42 problematic dives
  out of 108, nearly all from June–October 2025; since January 2026 there
  were 2 problematic out of 36.
- John Brewer had 4 ghost dives and 91 tows without a strong scar, yet scored
  0.26. 20b's model is fitted on dived reef-years only, so it predicts
  P(problematic | crews dive).
- The outbreak forecast was not the cause: both reefs' 2026 forecasts are
  below 0.03.
- Batt Reef's north-eastern flank went from 24 A–C sites (within-reef) to 43,
  because the scaled quota lifts a dense reef's whole order.

**Alternatives tested in 36.** All use the same replay (2020–2026, 1,683
visits), with each layer's class cut-offs calibrated on the training seasons
to the within-reef classes' pooled effort.
- *Densities* (fitted on cull visits at their planning dates):
  - 12 m reef-level model (control);
  - recent, each dive and tow weighted 0.5^(age / 90 d);
  - recent without the forecast;
  - **recent two-part**: P(a survey occasion needs culling | evidence) ×
    share problematic if dived. Survey occasions are a cull-site reef's dive
    and tow dates split at gaps of more than 10 days, so tow-only surveys
    count.
- *Ways to keep a dense reef's sparse parts low* (each on the recent two-part
  density):
  - lift halved;
  - no lift where the site was searched clean in the last year;
  - site model: P(dive problematic) on reef density, logit rank and the
    site's own record of dives within 300 m and linked tows, recent and from
    1–6 years before;
  - site visit model: the same terms, fitted on every designated site with the
    outcome "problematic dive here on this visit", so never-dived sites count;
  - site + 0.5 × visit model;
  - quota version: two-part quota, with sites ordered by the site + visit
    score;
  - **density-aware v2**: site + visit on reefs with dives or tows in the year
    before planning, the quota version elsewhere.
- *Outputs:*
  - `density_methods_accuracy.csv`: between-visit AUC, with a half-life
    sensitivity;
  - `density_layers_matched.csv`: matched effort with paired CIs;
  - `protocol_density_*.csv` (all layers), `density_layer_thresholds.csv`,
    `density_site_model_coefficients.csv`;
  - `density_layers_2027_sites.csv` and `_totals.csv`,
    `density_methods_2027_reefs.csv`;
  - `density_clean_tow_sites.csv`, `density_batt_flank_2027.csv`,
    `density_site_suppression.csv`.

**Between-visit ranking of dives (pooled AUC):**
- 20b density: 0.699; 12 m: 0.734; recent two-part: 0.740.
- Half-lives from 60 to 180 days: 0.737–0.740. The forecast adds nothing.

**Matched effort** (`density_layers_matched.csv`; bootstrap intervals in the
CSV):

| Layer | Ghost avoided at 3 % missed | Low avoided at 3 % missed | > 0.04 missed at 14 % skipped | Zone misses at 12 % tows saved |
|---|---:|---:|---:|---:|
| within-reef | 25.7 % | 16.6 % | 3.73 % | 10.4 |
| density-aware v1 (20b) | 33.7 % | 19.3 % | 3.02 % | 2.4 |
| recent two-part | 38.3 % | 19.2 % | 2.71 % | 2.9 |
| site visit model | 29.1 % | 18.0 % | 3.44 % | 6.8 |
| site + 0.5 × visit | 38.6 % | 22.4 % | 2.53 % | 1.5 |
| quota, site + visit order | 37.6 % | 18.5 % | 2.85 % | 3.5 |
| **density-aware v2** | **40.1 %** | **22.3 %** | **2.43 %** | **1.5** |

v2 minus v1 (paired):

| Measure | Difference | 95 % CI |
|---|---:|---|
| Ghost avoided | +6.4 points | +0.7 to +12.3 |
| Low-density avoided | +3.0 points | −0.3 to +6.9 |
| > 0.04 missed | −0.58 points | −1.13 to +0.01 |
| Zone misses | −0.9 | −7.0 to +3.9 |

**Operating points** (v2 classes are stricter, so the same cut skips more):

| Rule | > 0.04 missed | Ghost avoided | Low avoided | Ghost per miss |
|---|---:|---:|---:|---:|
| v2 dive rule A–D | 123 (2.2 %) | 33.9 % | 18.0 % | 8.2 |
| v2 dive rule A–C | 230 (4.1 %) | 49.0 % | 28.3 % | 6.3 |
| v1 dive rule A–C (reference) | 166 (3.0 %) | 33.5 % | 19.1 % | 6.0 |

- The v2 A–B zone saves 24.1 % of follow-up tows with 9 misses; v1 A–B saves
  14.3 % with 4.
- v2 zone A–B with rule A–D: 132 missed, 34.7 % ghost avoided.

**Findings that shaped v2:**
- **Weight 0.5 was chosen on the test seasons.** The exploratory grid ran from
  0.5 to 5. Larger weights pushed long-clean flanks down (Batt NE flank: 23
  A–C at w = 5) but lost misses and tow-zone performance. The margin is
  therefore optimistic.
- **Long-clean ground is not safe ground.** On dense reefs, sites with at
  least 5 tows and no strong scar or dive in the previous 1–6 years were dived
  405 times; 51 % of those dives were problematic and 132 were above 0.04
  (`density_clean_tow_sites.csv`). Every layer lifts Batt's NE flank; the
  flank had 0.4 % strong tows among its 1,337 tows since 2019, against 6.2 %
  on the SW flank (`density_batt_flank_2027.csv`).
- **The site models over-allocate where there is no evidence.** On reefs with
  no recent dives or tows (forecast only), they gave about 60 % A–C against 27
  % under the two-part quota, and about 1 in 5 such reefs needed culling.
  Hence the hybrid.

**Untowed sites neutral (10 November).**
- *Problem.* Tongue showed it clearly: 155 of its 173 A–B sites had no tow or
  dive within 100 m in the year (lagoon patches). The 74 towed sites, including
  11 near last season's strong tows, sat in C–E. With no tow data, a site
  entered the models with a strong-tow rate of 50 % (logit 0), far above a
  typical towed site. John Brewer's lagoon showed the same effect.
- *Fix.* An untowed site gets the mean recent strong-tow rate of the
  visit's / reef's towed sites (as if towed with the reef's average result);
  its 1–6 year record still counts.
- *Rejected first:* shifting untowed scores to the towed sites' median. On
  dense reefs it lifted ground crews skip because it is usually clean: Batt's
  NE flank, 38 of its 49 sites untowed, went to 39 A.
- *Replay cost* (`density-aware v2, untowed neutral` vs `density-aware v2`):

| Measure | v2 | Untowed neutral |
|---|---:|---:|
| Dive rule A–D: missed / ghost per miss | 123 / 8.2 | 123 / 8.2 |
| Ghost avoided at 3 % missed (matched) | 40.1 % | 40.3 % |
| Zone A–B: follow-up tows saved | 24.1 % | 25.0 % |
| Zone A–B: > 0.04 missed | 9 | 16 (+7, CI 0 to +17) |
| Zone A–B + rule A–D: missed / ghost per miss | 132 / 7.8 | 139 / 7.5 |

  Some untowed ground did hold aggregations. Adopted as the v2 product.

**Habitat entered directly (10 November; tested, not adopted).**
- *Method.* The untowed-neutral v2 with the site's habitat percentile on its
  reef (33's forward cull layer per season; 2027: 20's deployed rank) and
  habitat × untowed in both site models.
- *Fitted effect (2026 fit).* Habitat is +0.35 across the full range in the
  dive model and about 0 in the visit model. Habitat × untowed is −0.51:
  where there are no tows the slope is slightly negative.
- *Untowed sites* (`density_untowed_auc.csv`). Within-visit AUC among 1,969
  dives on 314 visits:

| Score | AUC (problematic) | AUC (> 0.04) |
|---|---:|---:|
| Habitat alone | 0.546 | 0.558 |
| Voyage-start rank | 0.554 | 0.532 |
| v2 untowed neutral | 0.609 | 0.599 |
| v2 untowed neutral + habitat | 0.607 | 0.605 |

- *Replay.* Matched-effort results are unchanged. Zone A–B misses 14 against
  16.
- *Batt.* Batt's NE flank is unchanged at 43 B, although its habitat
  percentile averages 0.31 against 0.72 on the SW flank: the reef's density
  dominates.
- *Conclusion.* Habitat is already reflected through the rank, and the sites'
  older records separate untowed sites better than habitat does.

**Where on the reef density comes from: v3 and v3b (10 November; tested, not
adopted).** The question was whether a site's density could come from nearby
evidence rather than the whole reef, so that a dense flank (Batt's SW) stops
lifting a quiet one (the NE).
- *Quick check.* Tow evidence pooled over designated sites with a Gaussian
  kernel ranks dives within a visit better than a site's own 100 m tows:
  within-visit AUC 0.62 vs 0.58 in 2020–23, 0.63 vs 0.58 in 2024–26 (500 m
  kernel, older tows included).
- *Aspect.* Weighting the kernel by aspect similarity, the direction the nearest
  reef edge faces from the GBR outlines, gave no consistent gain. It handled
  opposite flanks correctly: Batt NE 60°, SW −122°.
- *v3 (in 36).* The untowed-neutral v2 plus a local tow-density term, with and
  without aspect. The term was small: 0.09–0.12 in the dive model and 0.21 in
  the visit model. 2027 results:
  - Batt's NE flank: 47 A–C against 48;
  - Tongue: A–C rose from 148 to 172–177;
  - John Brewer: A–C rose from 20 to 29–30.
  - Not run through the full replay.
- *v3b (quick test, scratch).* Each site's density from the recent two-part
  model applied to kernel-weighted local dives and tows (500 m), shrunk to the
  reef with 2 pseudo-dives and 5 pseudo-tows, refitted strictly forward.
  - Within-visit AUC was almost the same as v2: 0.707 vs 0.704 (2020–23) and
    0.683 vs 0.678 (2024–26); no gain at untowed sites.
  - Batt: local density put the NE flank at 0.47 against 0.64 on the SW, yet
    the NE flank kept 44 of 49 sites in A–C (v2 48).
  - Tongue fell to 94 A–C, but John Brewer rose to 47 and the total to 6,425.
- *Conclusion.* Stopped. Local density separates Batt's flanks only modestly,
  and the dives do not reward pulling the NE flank down. The season's first
  full tow and re-ranking before each voyage remain the way to drop clean
  ground.

**2027 product (37, untowed-neutral v2).**
- Writes `cots_cull_site_priority_density_v2_2027.csv`, `.kmz` (21 adds
  exploration zones) and `_reefs.csv`.
- 676 A and 6,110 A–C sites, against 6,206 A–C for v1 and 10,170 within-reef.

| Reef | v2 density | v1 density | v2 A–C | v1 A–C |
|---|---:|---:|---:|---:|
| Batt | 0.78 | 0.69 | 90 of 92 | 82 |
| Tongue | 0.37 | 0.47 | 148 of 270 (A–B 28) | 206 |
| Big Broadhurst | 0.14 | 0.44 | 43 of 145 | 107 |
| John Brewer | 0.05 | 0.26 | 20 of 101 | 57 |

**Caveats:**
- The weight was chosen on the test seasons.
- The visit model encodes where crews chose to dive.
- Reefs with tows but no dives get 64 % A–C under v2, and the
  replay scores them only where crews later dived.
- Never-dived reefs cannot be tested.

## Why the old SDM validation ROC looked better (2 November 2026)

`DataProcessing/SDM_Validation_Report.qmd` showed AUC 0.865 for the cull
XGBoost Reef Guide model and the weighted ensemble on 1,535 2025 cull dives
(top 30 reefs). Checked with `scratchpad` scripts, not part of the pipeline:

| Same 1,535 dives (2025) | Pooled AUC | Within-reef AUC |
|---|---:|---:|
| RG model, at dive GPS | 0.865 | 0.767 |
| Ensemble | 0.865 | 0.769 |
| Each reef's mean RG score only | 0.792 | 0.5 |
| Clean (12-predictor) model | 0.660 | 0.623 |

- **Pooled across reefs:** reef identity alone scores 0.79; that is most of
  the pooled curve.
- **In-sample:** the 2025 raster was fitted on data including all 2025 dives
  (training AUC 0.968 in 02's metrics; out-of-fold 0.690).
- **Point memorisation:** scored at the cull-site location instead of each
  dive's GPS point, the same raster's within-reef AUC on 2025 dives falls to
  0.563.
- **Unseen 2026 dives:** at dive points (March–May 2026, 440 dives) pooled
  0.79 (reef mean alone 0.81), within reef 0.656; at site locations
  (1,082 dives) 0.606 vs 0.603 for the forward habitat layer and 0.74–0.76
  for the current combined ranking (33).
- **Conclusion:** the RG layers are already the habitat input to the product
  (RG + year, clean fallback). The old curve reflects between-reef
  differences and in-sample site memory, which the ranking now takes from
  dive and tow history explicitly. No better habitat map is being missed.
