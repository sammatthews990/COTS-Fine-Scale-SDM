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

Status: clean/RG, categorical geomorphic, site-aggregated manta-zero, and
annual-response ranking ablations implemented. The lifetime "ever high"
response is now diagnostic rather than the proposed final target. Annual
site-year observations are the primary response for the next refit.

- Confirm the preferred observation unit: first visit, first campaign, or a
  hierarchical count model using every dive.
- Replace single 30 m site points with polygon summaries or actual dive tracks.
- Model counts with bottom-time effort, or fit an explicitly effort-standardised
  threshold model.
- Add site and reef structure so repeated observations do not act as independent
  samples.
- Retain the full annual sequence instead of collapsing a site to "ever high".
  Weight site-years so a repeatedly visited site does not dominate the fit.
- Treat year as a prevalence or calibration adjustment, separate from the
  static spatial score. Do not let a gain in between-year discrimination count
  as evidence of better within-reef habitat ranking.
- Keep alternative history summaries, such as the fraction of adequately
  surveyed years above threshold, as sensitivity analyses. Do not use the
  unadjusted historical maximum.
- Add repeated manta non-detections to the presence stage with tow effort and
  survey-method effects. Keep manta/cull conflicts explicit and do not duplicate
  correlated tow rows as unweighted pseudo-absences.
- Validate on completely unseen reefs and report extrapolation/uncertainty.
- Report within-reef discrimination separately from between-reef
  discrimination so a reef-selection gain cannot be mistaken for a habitat
  gain.

Decision gate: accept a habitat prior only if it ranks held-out sites within
the same reef-year, captures more high-density observations in the top survey
fraction, transfers to unseen reefs, and is stable across forecast cut-off
years. Probability spread and map colour scaling are not acceptance metrics.

The first scale decomposition shows that the unconstrained outbreak bridge
improves overall 2026 AUC while reducing within-reef AUC. See
VISUAL_CHECK_FINDINGS.md. This motivates an explicit two-stage or constrained
cross-scale architecture.

The follow-up site-habitat ablation deliberately excludes every reef-level
predictor and validates held-out cull sites within reefs. RG improves
within-reef AUC from about 0.553 to 0.626. Geomorphic class adds no measurable
incremental value. Full-strength manta augmentation harms RG; a 0.10
effective-weight version reaches about 0.630 and remains a sensitivity model
pending repeated-seed and stricter-zero checks. See SITE_HABITAT_FINDINGS.md.

RG is unavailable at 566 otherwise clean-covered cull sites across 167 reefs.
A within-reef RG-with-clean-fallback test retains complete clean coverage and
achieves AUC 0.650 overall and 0.620 within reef, compared with 0.572 and 0.545
for clean-only. This is the provisional delivery layer while RG coverage is
audited.

The annual-response rank audit retains the full precision of every prediction.
Within-reef-year AUC is exactly invariant when raw probabilities are replaced
by within-group percentile ranks, demonstrating that compressed values are not
being treated as ties or discarded. On held-out sites, RG without year gives
the strongest annual within-reef-year result (AUC 0.592; top-20% lift 1.147).
Adding year improves overall AUC but not the within-reef-year ordering. The
legacy year raster has much stronger descriptive agreement with 2026 sites,
but must be reproduced with leakage-safe held-out predictions before that
comparison can be accepted as validation. See RANK_SIGNAL_FINDINGS.md.

A follow-up fixed legacy-spec check used the saved clean-year hyperparameters
and the legacy all-dive annual pooled response, while holding every site out
across all of its years. Year increased overall AUC from 0.570 to 0.662, but
within-reef-year AUC only from 0.555 to 0.568; top-20% lift decreased from 1.104
to 1.095. The descriptive legacy raster score of 0.730 therefore does not
transfer to held-out sites. A nested or pre-registered RG comparison is the next
model-fitting gate.

The rank-preserving calibration stage is now implemented. Training-fold year
intercepts improve RG overall AUC from 0.624 to 0.684 and Brier score from 0.239
to 0.225, while within-fold reef-year rankings are exactly unchanged. The
static RG score and time calibration are therefore separate model products.

Across three independent site-fold assignments, annual RG has mean
within-reef-year AUC 0.590 (SD 0.004), versus 0.541 (SD 0.011) for clean. RG is
better on all three splits and retains a higher mean top-20% lift (1.129 versus
1.085). A blind expert-review application has been generated using percentile
surfaces, with the model key held separately.

The legacy clean-plus-year model is now frozen and exactly reproducible on a
2,000-cell audit. The first 30/60/120 m support test shows that multiscale
topography plus year modestly improves held-out within-reef-year AUC over the
30 m baseline (0.574 versus 0.568). The static multiscale candidate reaches
0.586 and the best top-20% lift, but has much lower overall AUC because it does
not encode outbreak/year state. This reinforces the two-component architecture.
The current scale layers are block means of derived 30 m variables; a DEM-first
Gaussian pyramid and re-derived terrain metrics remain required.

A joint four-response Hmsc workflow is also implemented on matched 200 m
Cull/manta cell-years. It includes method effort covariates, repeated-cell and
reef random levels, four longer chains, and five-fold cell-blocked validation.
Hmsc is treated as a joint evidence model rather than an explicit
occupancy/detection model.

The completed Hmsc run has high detection-conditioned discrimination (Cull
high-density AUC 0.857 overall and 0.756 within reef-year), but it fails the
habitat interpretation gate. An effort-only blocked GLM reaches within-reef AUC
0.777, while static habitat reaches only 0.543. Within reef-year, Hmsc
probability correlates 0.889 with Cull dive count and 0.845 with bottom time.
The Hmsc result therefore diagnoses preferential sampling/observation effort;
it is not a candidate pre-survey habitat layer. A subsequent model must
separate ecological state from detection or predict at standardised effort.

## Two-part fixed-effort iteration 01

Status: completed; statistical gates passed, candidate not promoted pending
outbreak provenance and blinded expert map review.

The first explicitly separated model was evaluated on 9,321 site-years, 4,627
sites and 346 reefs using three repeated five-fold, within-reef spatial-site
partitions. Four static candidates excluded reef identity, year, outbreak
probability and effort. Multiscale topography fitted to annual high-density
responses was the stable winner:

- mean overall AUC 0.582;
- mean within-reef-year AUC 0.577 (SD 0.0007; minimum 0.577);
- mean top-20% capture 0.254 and lift 1.111.

The incumbent legacy clean-year validation values are within-reef-year AUC
0.568 and top-20% lift 1.095. The new static candidate therefore passes the
pre-registered statistical thresholds, but the gain is small and does not
override spatial expert evidence.

A nested rank-locked bridge adds same-year reef outbreak state. It raises
overall AUC to 0.743 while leaving within-reef ordering at AUC 0.577. The
fixed-effort pooled layer uses two dives and 600 bottom-time units and gives
the same overall/within AUC (0.743/0.577), with slightly worse Brier score than
the state bridge (0.211 versus 0.208). It is useful as a standardised
observation scenario, not as evidence of improved habitat.

Supplying realised effort produces AUC 0.824 overall and 0.782 within
reef-year. This deliberately retained diagnostic confirms the Hmsc finding:
realised survey effort contains strong preferential operational information
and must never enter a pre-survey map.

The deployment contract now has three named products:

1. static fine-scale habitat rank;
2. reef-year outbreak-state-adjusted absolute probability; and
3. fixed-protocol observation probability.

State and fixed-effort components have zero within-reef rank violations. The
machine-readable successor registry retains the legacy clean-year model as
current best until the outbreak provenance audit and blinded expert map gate
also pass. The detailed decision is in
`reports/TWO_PART_FIXED_EFFORT_ITERATION_01.md`; the rolling incumbent record
is `reports/CURRENT_BEST_MODEL.md`.

Every future iteration report and current-best update must include two to four
named reef panels comparing the exact saved candidate with the incumbent
raster. Panels use matched within-reef percentile scales and a candidate-minus-
incumbent rank-difference map. Report generation must fail when the visual-map
manifest is missing, refers to a different candidate, or contains fewer than
two or more than four reefs. Visual review remains a human promotion gate and
is never inferred from colour range or rank correlation alone.

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
- Evaluate fine-scale rankings using within-reef-year pairwise concordance,
  within-reef percentile contrast, top-k capture, COTS missed, and ghost dives
  avoided. Evaluate calibration, ROC-AUC, PR-AUC, Brier score, and log loss
  separately; do not use visual colour stretch or raw prediction range as a
  performance criterion.

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
