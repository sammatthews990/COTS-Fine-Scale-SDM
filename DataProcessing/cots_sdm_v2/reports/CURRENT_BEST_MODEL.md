# Current best COTS fine-scale model

Last evaluated: 2026-09-14

> Focused Hmsc v2 update: the reduced geomorphic 360 m posterior now passes
> convergence and retains the preferred coherent geomorphic map structure, but
> was not promoted because leave-reef and forward transfer remain weak. See
> [iteration 03](HMSC_V2_GEOMORPHIC_G360_ITERATION_03.md).

> Final synthesis: five publication-ready figures now consolidate blocked
> validation, density validation, variable importance, partial dependence and
> four-reef spatial review. See the
> [final model synthesis](FINAL_MODEL_SYNTHESIS_20260914.md).

> Iteration 07 update: continuous Cull/SALAD density and RosettaCOTS joint
> models were not promoted. See [the current dated decision](CURRENT_BEST_MODEL_20260911.md)
> and [all-model review](ALL_MODEL_REVIEW_20260911.md). Legacy clean + year
> remains the incumbent.

## Current decision

**Legacy clean + year remains the incumbent.**

Iteration 06 adds a matched 60 m geomorphology-led Hmsc candidate and the first
integrated shared-state, multi-method occupancy/detection model. Neither passes
all statistical, transfer, coverage and expert-map gates.

## Current development candidates

| Candidate | Held within | Leave-reef within | Forward within | Decision |
|---|---:|---:|---:|---|
| v4 enhanced Hmsc 120 m | 0.569-0.570 | 0.545 | 0.515 | development baseline; not promoted |
| v4 enhanced Hmsc 60 m | 0.560 | 0.546 | 0.543 | expert review; convergence fails |
| v5 shared-state 120 m | 0.604 | 0.579 | 0.504 | correct framework; forward/convergence fail |
| v2 reduced geomorphic 360 m | 0.558 | 0.540 | 0.531 | converged structural benchmark; not promoted |

The 60 m result suggests that observation support matters: it substantially
improves forward transfer while remaining competitive on unseen reefs, but its
posterior has minimum ESS 63.7 and upper PSRF 1.336. Longer chains are deferred
until experts accept its spatial surface.

The v5 model is a real observation-process model rather than further response
augmentation. Cull, Manta, SALAD and RHIS observe a shared latent problematic-
density state through separate detection functions; scars and Cull timing are
handled separately. Its strong held and leave-reef scores are outweighed by
chance forward-year performance and incomplete optimization/identifiability.

## Coverage and deployment rule

The development baseline remains the 120 m geomorphology-led enhanced Hmsc
surface with a fixed hierarchy:

1. enhanced geomorphology, exposure, Acropora and multiscale terrain;
2. clean Acropora plus multiscale terrain;
3. clean native stack;
4. explicit missing where the source stack itself is absent.

Fallbacks are calibrated for seam consistency and never overwrite an available
enhanced prediction. A prediction-source layer is mandatory. Reef identity,
year, outbreak probability, survey effort and method remain excluded from the
deployable within-reef habitat score.

## Reports and expert review

- [Iteration 06: shared-state and 60 m results](../iterations/occupancy_detection_v5/reports/OCCUPANCY_DETECTION_V5_ITERATION_06.md)
- [Expert-zone review instructions](../iterations/occupancy_detection_v5/reports/EXPERT_ZONE_REVIEW_20260910.md)
- [Editable expert-zone record](../iterations/occupancy_detection_v5/reports/expert_zone_review_20260910.csv)
- [Iteration 05: 200 m support and Manta augmentation](../iterations/hmsc_corrected_v4/reports/HMSC_CORRECTED_V4_ITERATION_05_SUPPORT_AND_MANTA.md)
- [Focused iteration 03: converged reduced 360 m geomorphic Hmsc](HMSC_V2_GEOMORPHIC_G360_ITERATION_03.md)
- [Final synthesis and five final figures](FINAL_MODEL_SYNTHESIS_20260914.md)

The iteration 06 report embeds four reef comparisons for both new candidates:
Eyrie, Batt, Broadhurst and Bait. Expert decisions remain pending and must cite
specific plausible or implausible zones; matching legacy rank is not required.

## Harmonised next steps

1. Complete the expert-zone record for the old and new panel families,
   including the converged reduced 360 m structural control. Judge named
   geomorphic zones on their operational plausibility, not incumbent rank.
2. Extend Hmsc chains only if the 60 m surface is accepted; require ESS >=400
   and upper PSRF <=1.1 plus held-cell, leave-reef and forward acceptance.
3. Refine the shared-state likelihood: actual scar effort, scar recency/error,
   stricter fold convergence, and a joint rather than composite Manta process.
4. Add a dynamic reef-year activation process that cannot alter within-reef
   rank; keep the supplied outbreak layer separate from static habitat.
5. Prepare leakage-safe pre-treatment RHIS branching/table coral condition and
   acquire independent current/connectivity data.
6. Apply the clean fallback hierarchy before any shared-state candidate can be
   considered deployable.
