# Current best COTS fine-scale model

Last evaluated: 2026-09-11

## Current decision

**Legacy clean + year remains the production incumbent. No iteration 07
successor was promoted.**

The final continuous-density experiment retained native Cull counts/effort and
exact SALAD track counts/areas, and used a verified per-hectare RosettaCOTS link.
It showed that replacing `CPUE >= 0.02` with abundance does not solve the main
transfer problem. Cull density has effectively zero within-reef rank under all
blocked tests. SALAD-native 60 m has promising forward-year rank (rho 0.218)
but fails held-cell (-0.019) and leave-reef (0.042) replication. Joint and
method-balanced Rosetta models conflict with SALAD and are rejected.

## Retained portfolio

| Role | Model | Status |
|---|---|---|
| Operational incumbent | legacy clean + year | frozen and accepted |
| Static-habitat development baseline | v4 enhanced Hmsc 120 m plus clean fallbacks | not promoted; expert/convergence/forward gates |
| Support sensitivity for expert review | v4 enhanced Hmsc 60 m | not promoted; convergence/expert gates |
| Correct future framework | v5 multi-method shared state | not promoted; forward/identification gates |
| Density research layer | v6 SALAD-native 60 m | not promoted; held/leave/coverage gates |

The density maps have enhanced coverage of only 16% at Eyrie, 50% at Batt, 20%
at Broadhurst and 11% at Bait. A clean fallback was deliberately not extended
to a surface that fails transfer. The established v4 hierarchy remains the
deployment rule for any future accepted geomorphic score: enhanced, then clean
multiscale, then clean native, with a prediction-source layer.

## Required decision before more modelling

Complete the four-reef zone-level review across legacy, v2/v4, v5 and density
panels. If no challenger is explicitly accepted, close this modelling cycle and
design a prospective fixed-effort, pre-cull validation survey rather than fit
more variants to the same observational data.

## Current reports

- [Iteration 07 continuous density](../iterations/density_rosetta_v6/reports/DENSITY_ROSETTA_V6_ITERATION_07.md)
- [All-model review](ALL_MODEL_REVIEW_20260911.md)
- [Iteration 06 shared state and 60 m Hmsc](../iterations/occupancy_detection_v5/reports/OCCUPANCY_DETECTION_V5_ITERATION_06.md)
- [Iteration 05 support and Manta](../iterations/hmsc_corrected_v4/reports/HMSC_CORRECTED_V4_ITERATION_05_SUPPORT_AND_MANTA.md)

## Harmonised next steps

1. Record accepted/rejected zones on all four reefs; colour stretch and incumbent
   rank agreement are not acceptance criteria.
2. Choose whether field operations should use legacy alone or legacy plus one
   explicitly experimental geomorphic layer.
3. Extend Hmsc chains only for an expert-accepted surface; require ESS >=400,
   upper PSRF <=1.1 and all blocked transfer gates.
4. Revisit shared-state abundance only with new repeat fixed-effort data,
   independent connectivity/current information, scar-age/error data and
   temporally prior coral condition.
5. Keep outbreak activation separate from within-reef habitat rank and retain
   the clean fallback hierarchy for every deployable candidate.
