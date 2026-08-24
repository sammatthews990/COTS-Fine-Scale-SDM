# Implementation plan

## Objective

Build two complementary map products and keep their interpretations explicit:

- a stable fine-scale habitat prior; and
- a target-year operational risk surface informed by reef-scale outbreak
  dynamics.

A later hurdle model will estimate high density conditional on confirmed COTS
presence.

## Phase 0 — foundations and controlled benchmark

Status: complete.

- Create an isolated project directory and immutable input configuration.
- Build stable raw-record, site, site-year, and first-visit identifiers.
- Replace `CPUE_max` as a target with effort-aware pooled CPUE targets.
- Extract static predictors without losing row-to-geometry correspondence.
- Import reef-year outbreak probabilities and construct lag-only features.
- Assign whole reefs, rather than individual rows, to spatial folds.
- Reserve 2026 for genuine forward validation.
- Compare feature ablations with identical rows and model settings.

Decision gate: determine how much discrimination is due to static habitat,
calendar year, and lagged outbreak state when evaluated on the same data.

Outcome: calendar year alone adds little under controlled comparisons. Static
predictors do not currently transfer well to unseen reef clusters, while lagged
outbreak history materially improves spatial and rolling-origin discrimination.
See `BENCHMARK_FINDINGS.md`.

## Phase 1 — static habitat prior

Status: next implementation phase.

- Confirm the preferred observation unit: first visit, first campaign, or a
  hierarchical count model using every dive.
- Replace single 30 m site points with polygon summaries or actual dive tracks.
- Model counts with bottom-time effort, or fit an explicitly effort-standardised
  threshold model.
- Add site and reef structure so repeated observations do not act as independent
  samples.
- Test alternative history summaries such as the fraction of adequately
  surveyed years above threshold. Do not use the unadjusted historical maximum.
- Add repeated manta non-detections to the presence stage with tow effort and
  survey-method effects. Keep manta/cull conflicts explicit and do not duplicate
  correlated tow rows as unweighted pseudo-absences.
- Validate on completely unseen reefs and report extrapolation/uncertainty.

Decision gate: accept a habitat prior only if it transfers to unseen reefs and
is stable across forecast cut-off years.

## Phase 2 — dynamic outbreak risk

Status: adjusted hindcast/forecast bridge implemented provisionally; provenance
audit still required.

- Treat the reef-level outbreak model as a reef-year prior, not as fine-scale
  habitat.
- Audit whether each historical probability is an out-of-fold hindcast or an
  in-sample fitted value.
- Use lagged outbreak probability, recent trend, and multi-year maximum/mean.
- Add time since last positive survey, time since culling, recent culling
  effort, and connectivity/distance-to-outbreak features.
- Use rolling-origin evaluation: train through year `t-1`, predict year `t`.

Decision gate: dynamic features must improve forward-year calibration and COTS
capture at a fixed operational survey budget.

## Phase 3 — conditional high-density hurdle

- Stage A: probability of presence from standardised survey data.
- Stage B: probability of CPUE >= 0.02 conditional on presence.
- Combine the two stages only when an unconditional high-density probability is
  needed.
- Compare classification with a negative-binomial or zero-inflated count model
  using log(bottom time) as an offset.

## Phase 4 — spatial prediction and delivery

- Fit final models only after the above decision gates are met.
- Generate separate habitat-prior and target-year operational rasters.
- Publish uncertainty and environmental-extrapolation layers beside the mean.
- Evaluate rankings using top-k capture, COTS missed, ghost dives avoided,
  calibration, ROC-AUC, PR-AUC, Brier score, and log loss.

## Outbreak-layer provenance requirement

The configured files provide adjusted 1991-2025 hindcasts and a separate 2026
forecast. Before the current-year bridge is accepted for performance claims,
record whether each value is:

- a genuine forecast made before that year's observations;
- an out-of-fold hindcast; or
- an in-sample fitted/hindcast value.

Only the first two are valid for performance evaluation. Lagging reduces direct
target-year leakage but does not fix leakage inside the upstream outbreak model;
that provenance still needs to be established.
