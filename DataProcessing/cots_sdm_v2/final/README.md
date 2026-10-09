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
| `cots_cull_site_priority_2027.kmz` | The EoTR cull-site KMZ, unchanged in structure (one folder per reef; same site names, IDs and polygons). Each polygon is coloured by priority (A red, B orange, C yellow, D grey), and its data and pop-up gain *Priority 2027*, *Rank on reef*, *P(problematic if dived)* and *Last season*. |
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
first. Retention is from the forward test, 2019-2026:

| Searching through | Sites on the reef | Problematic dives kept | Low-density dives avoided | Ghost dives avoided |
|---|---:|---:|---:|---:|
| A | top 20 % | 33 % | 75 % | 79 % |
| A + B | top 50 % | 64 % | 42 % | 49 % |
| A + B + C | top 80 % | 89 % | 18 % | 22 % |
| all (incl. D) | 100 % | 100 % | 0 % | 0 % |

How to use it:
- **Work down the classes; don't stop at A.** Class D (the bottom 20 % of a
  reef's sites) is the realistic candidate to skip or leave until last.
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

## Exploration areas: outside the cull sites and recent tow paths (21)

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
- **Candidate areas:**
  - eligible cells in the reef's top 30 % of habitat;
  - one-cell gaps closed through eligible cells, so nearby pockets join;
  - patches over about 9 ha split into compact chunks (k-means on cell
    centres; about one cull site each);
  - outlines rounded (+45 m / -30 m buffer, kernel smoothing, holes under
    1 ha filled, simplified to 8 m), then clipped to the reef and kept
    outside cull sites;
  - areas of at least 0.5 ha, ranked by mean habitat percentile;
  - the top 10 per reef polygon, at least 350 m apart.

**Outputs:**
- `layers/cots_exploration_2027_30m.tif`, three bands:
  1. exploration score (habitat percentile on eligible cells);
  2. distance to the nearest manta tow in the last 5 seasons (m);
  3. distance to the nearest designated cull site (m).

  It covers reefs with designated cull sites only.
- `site_ranking/cots_exploration_areas_2027.csv` and `.gpkg`: 3,807 areas on
  488 reefs, most with the full 10. Median area 8 ha, median habitat
  percentile 0.84, median 170 m from a cull site. 898 areas are on reefs
  with no manta tow in 5 seasons. Outlines are simplified to 8 m (area
  change < 1 %) to keep the KMZ small (12 MB).
- The 2027 KMZ gains an **"Exploration areas 2027 (purple)"** folder (one
  sub-folder per reef). Each area's pop-up gives its rank, area, habitat
  percentile, depth and distance to the nearest cull site and recent tow.
- fig17 shows Batt and Eyrie.

What this means:
- **On well-surveyed reefs the gaps are small.** On Batt and Eyrie, the cull
  sites and 5 seasons of tows already cover most of the high-habitat
  perimeter. On Batt the areas are outer-edge pockets and the ends of the
  reef. Eyrie also gets a group on the northern back reef near the island,
  which no tow has crossed in 5 seasons.
- **The larger opportunities are reefs with little recent towing.**
- **Status:** these are places to look, not predictions. The manta model is
  extrapolated away from tow paths, which is exactly where these areas lie.
  It ranks new ground about as well as old (AUC 0.62), but no layer is
  strong.
- **Suggested use:** add one or two exploration areas to a reef's manta
  plan, or as a single scouting dive per visit. Log the result so the next
  season's evidence includes it.
- **Settings:** clearances, depth band, habitat floor, chunk size, spacing
  and areas per reef are set at the top of script 21. The script takes
  about 30 minutes.

## Figures

All in `outputs/figures/`, 7.5 in wide at 300 dpi, built by `07` from the tables in `outputs/tables/` (run it last).

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
| `fig17_search_2027_exploration.png` | Batt and Eyrie: 2027 cull-site priority, manta tows of the last 5 seasons and exploration areas over imagery |

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
Rscript DataProcessing/cots_sdm_v2/final/scripts/21_exploration_areas.R         # exploration layer + areas, added to the 2027 KMZ (~30 min; needs 17, 20)
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
  scripts/        00-21 (see Run order; 07 runs last)
  outputs/
    layers/       delivered GeoTIFFs + layer_manifest.csv   (tifs git-ignored)
    tables/       every number quoted in the report            (tracked)
    figures/      fig01-fig17 PNGs                              (tracked)
    site_ranking/ 2027 cull-site priority (CSV, KMZ incl. exploration areas), exploration areas (CSV, GPKG)
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
