# Geomorphology-led Hmsc v3: iteration 03

Last evaluated: 2026-09-08

## Decision

**NOT PROMOTED. Legacy clean + year remains the incumbent benchmark.**

The new 120 m Hmsc surface is a credible independent habitat-ranking
candidate, but it is not yet a production successor. It passes the broad MCMC
mixing check and held-cell within-reef validation, but narrowly misses the ESS
gate, has weak forward-year transfer, and still requires expert review of the
four reef maps. The legacy predictor, probability and rank were not used in
training, feature selection, shrinkage or validation.

## What changed

- Created this isolated v3 directory; existing legacy and v2 artifacts remain
  untouched.
- Intersected manta and SALAD line segments with every crossed 60 m and 120 m
  cell. Survey length/area is allocated by within-cell line fraction.
- Used saved SALAD time-interpolated positions for localized COTS/scar events.
- Tested 1-, 3- and 7-day visit-alignment windows.
- Averaged native predictors over four 30 m subcells at 60 m and sixteen at
  120 m. Gaussian terrain covers 30, 60, 120, 240 and 360 m.
- Reduced the multiscale basis to one within-family PC for each of DEM, slope,
  eastness, northness, horizontal curvature, vertical curvature, BPI and local
  SD in the final Hmsc candidate.
- Restricted Hmsc to six COTS responses: Cull problematic density and direct
  presence, SALAD density and direct evidence, and RHIS adult and juvenile
  evidence. Mapped Acropora layers are covariates. Coral mortality is excluded.
- Fitted scar-persistence and Cull-depletion models outside Hmsc.
- Kept adjusted reef-year outbreak probability outside every habitat formula.

## Geometry and alignment audit

Manta supplied 214,319 usable start locations and 184,282 usable lines. SALAD
supplied 722 usable tracks and 2,917 valid time-interpolated event positions.
Endpoint distance was more than three times recorded survey length for 76
manta records and 37 SALAD tracks; these were flagged and reduced to start
points rather than allowed to cross implausibly large parts of the GBR.

| Support | Manta footprint-cell links | SALAD links | Cull points | RHIS points |
|---:|---:|---:|---:|---:|
| 60 m | 890,022 | 6,236 | 49,130 | 28,273 |
| 120 m | 552,108 | 3,500 | 49,130 | 28,273 |

| Support | Window | Alignment units | Multi-method units |
|---:|---:|---:|---:|
| 60 m | 1 day | 790,948 | 12,544 |
| 60 m | 3 days | 788,199 | 14,274 |
| 60 m | 7 days | 785,710 | 15,350 |
| 120 m | 1 day | 417,067 | 13,922 |
| 120 m | 3 days | 413,960 | 15,879 |
| 120 m | 7 days | 411,460 | 17,028 |

The 3-day window is the current working compromise. Window choice changes the
number of cross-method units, but none creates enough repeated localized point
visits to identify a fine-cell occupancy/detection model reliably.

## Predictor coverage

| Support | Candidate cells | Complete base support | Geomorphology observed | Complete Gaussian terrain |
|---:|---:|---:|---:|---:|
| 60 m | 26,464 | 24,875 | 24,045 | 26,267 |
| 120 m | 20,920 | 20,210 | 18,396 | 20,777 |

Cells outside a source raster remain missing and are excluded. No synthetic
geomorphic class, bathymetry or exposure value is imputed.

## Terrain-basis screening

Nested spatial shrinkage screens compared geomorphology/base covariates,
one-PC-per-family, up-to-two-PC family bases and all-scale lasso. The naive
pooled labels transferred poorly: best held-cell AUC was approximately 0.524,
while several leave-reef/forward fits shrank to an intercept. These screens do
not support a pooled-label successor. They motivated the smaller joint model
and observation-process separation.

## Formal occupancy/detection diagnostic

The first cell-level occupancy fit was rejected because manta-dominated zero
selection left an almost all-positive localized sample. After rebuilding
controls from localized Cull/RHIS/SALAD observations, 60 m and 120 m fits still
failed convergence: only 323 and 337 sites, respectively, contained repeated
localized visits. This is an identifiability failure at fine-cell support, not
a successful detection correction.

Manta positives were not assigned to every crossed cell. Their exact position
along a tow is unknown. A formal manta-inclusive occupancy model must operate
at a footprint/coarser site unit with footprint-weighted habitat covariates.

## Separate observation processes

Both auxiliary models converged.

- Scar persistence used 5,041 localized RHIS/SALAD observations across 393
  reefs. Scar probability declined with log days since prior direct evidence
  (coefficient -0.096, p=0.013). Post-Cull status alone was not significant.
- Cull depletion used 49,130 dives across 361 reefs. Catch rate declined with
  log days after the first local Cull record (coefficient -0.199, p<0.001),
  after offsetting by bottom time and modelling reef/cell random intercepts.

These results support keeping scars and Cull intervention outside the static
habitat score.

## Reduced Hmsc result at 120 m

The production-length fit used four chains, 2,000 retained samples per chain,
thin 5 and a 2,000-iteration transient. Maximum upper PSRF was 1.053. Minimum
beta ESS was 354, narrowly below the conservative 400 gate.

| Evaluation | Overall AUC | Within-reef-year AUC | Probability SD |
|---|---:|---:|---:|
| Apparent fit | 0.597 | 0.602 | not used as promotion evidence |
| Held 3 km cell blocks, 5-fold | 0.565 | 0.572 | 0.087 |
| Leave reef out, 5-fold | 0.493 | 0.565 | 0.099 |
| Forward: train through 2023 | 0.512 | 0.522 | 0.096 |

The leave-reef contrast is informative: overall discrimination fails, while
within-reef ranking remains above chance. Static habitat is therefore behaving
more like a local search-priority layer than a reef-outbreak model. Reef-year
outbreak probability should activate reefs separately and must not reorder
cells within a reef.

## Four-reef expert review

All maps are rendered at the fitted 120 m support. Percentiles reveal spatial
ordering; raw probability is shown separately. Colour stretching is not an
assessment metric. The legacy panel is context only and rank agreement is not
a promotion gate.

### Eyrie Reef

![Eyrie legacy, Hmsc percentile, raw probability and geomorphology](figures/v3_hmsc_120m_eyrie.png)

### Batt Reef

![Batt legacy, Hmsc percentile, raw probability and geomorphology](figures/v3_hmsc_120m_batt.png)

### Broadhurst Reef

![Broadhurst legacy, Hmsc percentile, raw probability and geomorphology](figures/v3_hmsc_120m_broadhurst.png)

### Bait Reef

![Bait legacy, Hmsc percentile, raw probability and geomorphology](figures/v3_hmsc_120m_bait.png)

Expert acceptance is pending for all four reefs. Descriptive Spearman
correlation with legacy ranges from 0.080 to 0.525 and is not used to accept or
reject the geomorphic model.

## Harmonised next steps

1. Complete blinded/operational expert review of all four 120 m panels. Record
   concrete accepted and implausible geomorphic features for each reef.
2. Build the manta-inclusive footprint occupancy model at an observation-scale
   site unit. Use line-weighted habitat covariates, method-specific detection,
   and repeated visits; do not spread a positive tow across every crossed cell.
3. Fit the reduced Hmsc response block at 60 m using the identical partitions,
   then compare 60 m and 120 m maps at fitted support.
4. Extend the final 120 m chain only if expert review is plausible; target
   beta ESS >=400 and upper PSRF <=1.1.
5. Improve forward transfer with independently mapped ecological covariates,
   especially connectivity/exposure and compatible susceptible-coral data.
   Keep method, effort, reef identity, year and outbreak state out of habitat X.
6. Validate the separate outbreak activation layer for reef-year calibration
   and transfer. It may change absolute activation between reefs but must never
   change within-reef habitat order.
7. Repeat held-cell, leave-reef, forward-year and four-reef maps for every new
   candidate. Promote only if convergence, independent transfer and expert
   operational acceptance all pass.

