# Current best COTS fine-scale model

Last evaluated: 2026-09-09

## Current decision

**Legacy clean + year remains the incumbent.**

The corrected conditional Hmsc v4 candidate is not promoted. It restores the
intended active-reef-year estimand, incorporates strong manta evidence at tow
support, restores the clean fallback, and recovers much of the operationally
preferred Hmsc v2 static-joint geography. It still fails forward transfer,
preliminary convergence, and expert-map acceptance gates.

## Latest candidate: corrected conditional Hmsc v4

| Evaluation | Overall AUC | Within-reef-year AUC |
|---|---:|---:|
| Apparent, enhanced joint | 0.600 | 0.598 |
| Held 3 km cell, five-fold | 0.571 | 0.570 |
| Leave reef out, five-fold | 0.514 | 0.545 |
| Forward, train through 2023 | 0.513 | 0.515 |

Strong manta positives are direct COTS, at least two counted scars, or legacy
common scars. Single/present scars remain a sensitivity definition. The manta
model operates at whole-tow support and is not copied into crossed cells.

The enhanced COTS-only and joint-resource target predictions are effectively
identical. The preliminary enhanced posterior fails convergence, and forward
ranking remains weak. The candidate therefore remains diagnostic.

The restored coverage hierarchy is enhanced, clean multiscale, then clean
native. It covers 100% of the clean DEM domain at Batt, Broadhurst and Bait and
99.2% at Eyrie; 19 Eyrie edge cells lack even native clean predictors.

Full methods, audits, maps and next steps are in
`../iterations/hmsc_corrected_v4/reports/HMSC_CORRECTED_V4_ITERATION_04.md`.

## Required expert review panels

### Eyrie

![Eyrie v4 final comparison](../iterations/hmsc_corrected_v4/reports/figures/v4_final_coverage_comparison_eyrie.png)

### Batt

![Batt v4 final comparison](../iterations/hmsc_corrected_v4/reports/figures/v4_final_coverage_comparison_batt.png)

### Broadhurst

![Broadhurst v4 final comparison](../iterations/hmsc_corrected_v4/reports/figures/v4_final_coverage_comparison_broadhurst.png)

### Bait

![Bait v4 final comparison](../iterations/hmsc_corrected_v4/reports/figures/v4_final_coverage_comparison_bait.png)

## Harmonised next steps

1. Record expert acceptance/rejection and concrete reef-zone observations for
   all four v4 panels.
2. Build and fit the corrected exact-footprint 200 m version, then compare it
   with 120 m at fitted support.
3. Default to the smaller COTS-only response block unless expert review shows a
   repeatable benefit from the nearly identical joint-resource surface.
4. Extend Hmsc chains only for an expert-accepted map and require beta ESS
   >=400 and upper PSRF <=1.1.
5. Improve forward transfer with independent connectivity/exposure and
   compatible coral-condition data while keeping outbreak state separate.
6. Maintain enhanced, clean multiscale and clean native outputs as explicit
   source layers; never silently impute ReefGuide gaps.
