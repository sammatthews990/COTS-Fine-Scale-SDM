# Hmsc visual assessment and joint-modelling roadmap

## Decision

**Do not promote the current Hmsc surface. Retain the frozen legacy
clean-year model as the current best operational map.**

The rejected successor was not Hmsc: it was the static multiscale XGBoost
candidate. The maps below are newly derived from the fitted production Hmsc
model and use exactly the same four reef comparisons. They are a diagnostic
of what Hmsc has learned about transferable fine-scale habitat.

## What was mapped

The mapped response is Cull high density (pooled CPUE at least 0.02). The
surface uses posterior-mean fixed effects for 2026, fixes Cull and manta
effort at the median matched-data protocol, and excludes reef and 200 m cell
random effects. It therefore asks: before a survey, and without being told
reef identity, how does the fitted Hmsc model rank habitat within a reef?

The fixed protocol is 600 Cull bottom-time units, 2 Cull dives, 2 manta tows, 393.1 m manta distance and a hard-coral observation fraction of 1.00.

The full held-out Hmsc score cannot be attached to this deployable map.
The full detection-conditioned fit achieved overall AUC 0.857 and within-reef AUC 0.756, but 0.6% of Cull high-density variance was assigned to multiscale habitat, 25.2% to Cull effort, and 69.7% to reef/cell identity effects.
As a check, the simpler held-out static-habitat ablation had AUC 0.527, whereas effort plus year alone had AUC 0.775.

## Four-reef visual comparison

All first and second panels use within-reef percentiles; this reveals rank
structure without treating the compressed raw Hmsc scale as a failure. The
third panel is Hmsc minus incumbent percentile. Crosses are historical Cull
locations and are not response symbols.

| Reef | Raw 5%-median-95% | Rank correlation | Mean absolute rank difference |
|---|---|---|---|
| Eyrie Reef (14-118) | 0.056-0.063-0.074 | -0.147 | 0.371 |
| Batt Reef (16-029) | 0.055-0.061-0.073 | -0.168 | 0.378 |
| (Big) Broadhurst Reef (No 1) (18-100a) | 0.053-0.066-0.085 | -0.089 | 0.334 |
| Bait Reef (19-137) | 0.051-0.062-0.080 | -0.053 | 0.339 |

The rank correlation is negative on every candidate reef. This does not by
itself prove the incumbent is ecologically correct, but it confirms that the
current Hmsc habitat component has learned a substantially different spatial
ordering and requires expert scrutiny rather than metric-led promotion.

### Eyrie Reef (14-118)

![Hmsc fixed-protocol versus clean 2026: Eyrie Reef (14-118)](figures/hmsc_vs_incumbent_01_eyrie.png)

### Batt Reef (16-029)

![Hmsc fixed-protocol versus clean 2026: Batt Reef (16-029)](figures/hmsc_vs_incumbent_02_batt.png)

### (Big) Broadhurst Reef (No 1) (18-100a)

![Hmsc fixed-protocol versus clean 2026: (Big) Broadhurst Reef (No 1) (18-100a)](figures/hmsc_vs_incumbent_03_broadhurst.png)

### Bait Reef (19-137)

![Hmsc fixed-protocol versus clean 2026: Bait Reef (19-137)](figures/hmsc_vs_incumbent_04_bait.png)

## Why the current joint model is not yet delivering a useful habitat map

1. The four responses are observation channels for one focal animal
   (Cull direct presence/high density and manta direct/scar presence), not a
   genuine biological community. Shared survey history can create residual
   association without identifying a habitat mechanism.
2. Reef identity, cell identity and realised effort explain almost all of
   the fitted Cull high-density variation. Those terms help reconstruct the
   historical sample but cannot create transferable 30-120 m habitat detail.
3. The observation methods have different footprints, effort and error.
   Aggregating them to 200 m cell-years does not explicitly model imperfect
   detection, scar false positives, depletion after culling, or preferential
   selection of sites.
4. The terrain representation is sparse: each terrain family is represented
   at only one selected support. It does not yet implement the paper's full
   DEM-first Gaussian multiscale derivation and scale-selection workflow.
5. Random-effect convergence remains inadequate (maximum upper PSRF 1.631; minimum ESS 48.5), so near-perfect residual correlations must not be interpreted as established ecological coupling.

### Current Cull high-density habitat coefficients

| Predictor | Posterior mean | Support positive |
|---|---|---|
| BPI_s30m | -0.036 | 0.029 |
| SLO_s60m | -1.087 | 0.078 |
| VRM_s60m | -17.702 | 0.312 |
| DEM_s30m | -0.006 | 0.236 |
| SVF_s120m | -6.292 | 0.332 |
| coral_suitability_mean | 0.646 | 1.000 |

Posterior support is one-sided for coral suitability, but weak for most
terrain terms. These coefficients should be treated as conditional
associations within the current detection-conditioned specification.

## Data already available locally

The combined 2026 Cull-Manta-RHIS workbook is much richer than the current
four-response Hmsc table. It already contains Cull size cohorts; manta
direct counts, scar counts, tow distance and coral-cover classes; and RHIS
adult/juvenile COTS counts, scars by coral growth form, live coral by growth
form, recent mortality, visibility, depth and survey habitat.

| Method | Rows | Years | Reefs | Named sites | Rows with COTS/scar evidence |
|---|---|---|---|---|---|
| Cull |  49,150 | 2002-02-22 to 2026-05-29 | 360 |  4,828 | 39,661 |
| Manta | 214,320 | 2012-11-29 to 2026-06-01 | 734 | 12,341 | 20,070 |
| RHIS |  28,273 | 2018-06-02 to 2026-05-20 | 604 |  2,907 |  2,671 |

| Comparison | Named units | Units repeated across days | Shared reef-days |
|---|---|---|---|
| Cull named site |  4,826 | 2,976 | - |
| Manta named site | 12,475 | 9,307 | - |
| RHIS named site |  1,646 |   777 | - |
| Cull-Manta same reef/day | - | - |   951 |
| Cull-RHIS same reef/day | - | - | 1,969 |
| Manta-RHIS same reef/day | - | - | 1,583 |

These are preliminary compatibility counts, not proof of matched sampling.
Same reef/day is deliberately weaker than same footprint and protocol. The
large number of repeated named sites and shared reef-days nevertheless shows
that a visit-level alignment audit is worthwhile before seeking new data.

The older ReefScan script expects a data/ReefScan CSV directory, but that
directory is not present in the current repository. Its previous ingestion
also discarded survey date, effort, segment/image identity, confidence, coral
composition and separate direct/scar outcomes. Those fields are required if
ReefScan is to strengthen detection modelling rather than add pseudo-replicates.

## Data that would make a joint approach materially stronger

| Priority | Compatible sampling unit and fields | Role in model |
|---|---|---|
| 1 | Individual Cull dive/track: date-time, observer, pre/post control, counts by size, bottom time, footprint, visibility, depth, scars, coral taxa/cover | Separates abundance from detection and removal |
| 1 | Manta tow segment/track: direct count, scar count, tow time/distance, observer, visibility, coral cover/taxa and georeferenced footprint | Method-specific detection and scar error |
| 1 | ReefScan image/segment effort with validated direct-COTS and scar labels, negatives, confidence and coral composition | Independent observation channel with quantified classification error |
| 1 | Repeat visits to the same cells before intervention, ideally paired across methods | Identifies non-detection versus true absence |
| 2 | Adult and juvenile COTS, fresh scars, susceptible coral taxa/cover, recent coral mortality and standardized predator-guild observations on compatible units | Creates a genuinely multivariate ecological response |
| 2 | Culling history, time since first detection/outbreak, removals and revisit interval | Models depletion and outbreak phase |
| 2 | DEM-first Gaussian terrain at 30/60/120/240/360 m plus geomorphic zone, wave/current exposure and larval connectivity | Multiscale habitat and dispersal mechanism |
| 3 | SST/heat stress, cyclones, chlorophyll/productivity, turbidity and time-varying coral condition aligned before each visit | Environmental state and disturbance context |

A response matrix becomes genuinely useful when the columns are compatible
ecological states observed on common units. Adult COTS, juvenile COTS, fresh
scars, susceptible coral groups and coral mortality are stronger candidates
than treating multiple detections of the same animal as if they were species.

## Recommended modelling sequence

1. **H1: transferable fixed-effect benchmark.** Refit without realised effort
   and reef/cell identity, using complete Gaussian multiscale terrain and
   coral composition. Validate held-out cells, reefs and future years.
2. **H2: observation-aware visit model.** Retain individual visits and add
   method-specific detection, effort, observer, scar error and Cull-removal
   processes. If Hmsc cannot express the required observation model, use an
   integrated multivariate model in Stan, NIMBLE or TMB and keep Hmsc as a
   benchmark rather than forcing the data into Hmsc.
3. **H3: random-structure sensitivity.** Compare reef only, cell only, reef
   plus spatial field, and regularised temporal structures. Require stable
   chains before interpreting residual associations.
4. **H4: genuine multivariate ecology.** Add compatible ReefScan, coral taxon
   and COTS life-stage responses; test whether shared latent gradients improve
   out-of-sample Cull high-density ordering, not merely fit to sampled reefs.
5. **H5: prospective operational test.** Pre-register fixed-effort predictions
   for stratified random sites within reefs, then score discrimination and
   calibration without using realised survey effort.

Every iteration keeps four gates: within-reef rank validation, leave-reef
transfer, forward-year validation, and two-to-four reef expert map review.
A higher global AUC cannot override a failed spatial-plausibility gate.

## Immediate next implementation

The next bounded experiment should be H1: refit the Hmsc response structure
with the expanded Gaussian multiscale predictors, no reef/cell random effects
in the deployable habitat score, and no realised effort. In parallel, build a
visit-level data audit showing which Cull, manta and ReefScan records can be
aligned on cell, date and pre-intervention status. H2 should start only after
that audit identifies enough repeated compatible visits to estimate detection.

## Reproducibility

- scripts/23_fit_hmsc_joint_production.R fits the production Hmsc model.
- scripts/26_summarise_hmsc_joint.R summarizes performance and parameters.
- scripts/27_diagnose_hmsc_effort_signal.R separates habitat and effort signal.
- scripts/32_build_hmsc_reef_maps.R creates these fixed-protocol maps.
- scripts/34_audit_joint_response_data.R audits available joint responses.
- scripts/33_build_hmsc_joint_roadmap.R regenerates this report.
- Hmsc method: https://doi.org/10.1111/2041-210X.13345
