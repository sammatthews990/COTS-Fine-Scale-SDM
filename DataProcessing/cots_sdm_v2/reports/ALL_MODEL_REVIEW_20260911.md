# COTS fine-scale SDM: all-model review

Last updated: 2026-09-11

> A later, figure-led synthesis using the same decision framework is available
> in [FINAL_MODEL_SYNTHESIS_20260914.md](FINAL_MODEL_SYNTHESIS_20260914.md).

## Executive decision

**Legacy clean + year remains the only accepted operational incumbent.**

The programme has now tested annual and year-agnostic classification,
multiscale Maxent, Hmsc joint responses, geomorphology-led reduced Hmsc at 60,
120 and 200 m, whole-track Manta evidence, a formal multi-method shared-state
model, and continuous Cull/SALAD density with RosettaCOTS. No challenger has
passed convergence, held-cell, leave-reef, forward-year, coverage and expert-map
gates together.

That does not mean the later work failed to add insight. It separates three
useful but currently incompatible components:

- the legacy clean + year model contains strong operationally recognised
  within-reef spatial structure and full clean-stack coverage;
- v2/v4 geomorphology-led Hmsc surfaces are easier to interpret as habitat and
  introduced defensible multiscale/support handling;
- v5 provides the correct direction for separating latent state from
  method-specific detection, scars and Cull intervention.

The continuous-density experiment confirms that the main limitation is not the
`0.02` threshold alone. Static habitat cannot currently extract stable
within-reef abundance rank from Cull CPUE, and Cull and SALAD do not combine
cleanly through a fixed cross-method calibration.

## Comparable candidate summary

Classification AUC and density Spearman correlations are different quantities
and should not be placed on a single leaderboard.

| Candidate | Held within | Leave-reef within | Forward within | Main strength | Blocking issue |
|---|---:|---:|---:|---|---|
| Legacy clean + year | AUC 0.568 | not historical gate | accepted historical output | operationally accepted pattern, full clean coverage | annual/state effects entangled with habitat |
| Two-part fixed effort v1 | AUC 0.577 | incomplete provenance | not decisive | first habitat/state separation | provenance and expert gates failed |
| Hmsc v2 static joint | AUC 0.545 | AUC 0.535 | AUC 0.570 | interpretable broad geomorphic pattern | ESS 17.7, PSRF 4.02; transfer weak |
| Hmsc v4 enhanced 120 m | AUC 0.569-0.570 | AUC 0.545 | AUC 0.515 | corrected support, reduced basis, fallback architecture | convergence and forward/expert gates |
| Hmsc v4 enhanced 60 m | AUC 0.560 | AUC 0.546 | AUC 0.543 | better forward transfer at finer fitted support | ESS 63.7, PSRF 1.336; expert pending |
| Hmsc v4 enhanced 200 m | AUC 0.541 | AUC 0.528 | AUC 0.528 | smoother broad support | loses held/leave signal and coverage |
| Hmsc + strong Manta response | AUC 0.569 | AUC 0.546 | AUC 0.515 | validates Manta as evidence | target map unchanged; convergence worse |
| Shared-state v5 120 m | AUC 0.604 | AUC 0.579 | AUC 0.504 | proper method/detection framework | chance forward result; incomplete identification |
| Cull density 60 m | rho -0.025 | rho -0.006 | rho -0.004 | preserves count and effort | no within-reef abundance rank |
| SALAD density 60 m | rho -0.019 | rho 0.042 | rho 0.218 | promising recent forward signal | held/leave not replicated; only 553 tracks |
| Balanced Rosetta joint 60 m | rho -0.035 | rho -0.051 | rho 0.026 | native method likelihoods, verified units | incompatible method signals |

## What the iterations established

### Why the year model can look better

The annual model can borrow strong reef-year prevalence and sampling-state
information. That helps answer “given COTS activity this year, where did survey
and removal encounter them?” and can produce an operationally convincing narrow
probability range. A year-agnostic habitat model asks a harder transfer question:
which cells are consistently suitable across outbreaks, methods and years?
Preferential survey placement, repeated culling, imperfect zeros, scar lag and
changing abundance state weaken that static signal.

The legacy output should therefore remain frozen and reproducible, but its
linear predictor must not be used as an offset, teacher or rank constraint for
new habitat models. Visual agreement is context, not an objective function.

### Multiscale terrain and fitted support matter

Full 30-360 m terrain bases were unstable and mixed poorly. Within-family PCA
and shrinkage preserve broad geomorphic form more defensibly. Matched 60 m and
120 m fits outperform the 200 m detour operationally; 60 m improved forward
classification while 120 m remains the better established development support.
No output should be presented at 30 m merely because the source raster is 30 m
when observations represent tracks or approximately 10 ha footprints.

### More responses are not automatically more information

Appending Manta or abundant coral responses to Hmsc barely altered the target
surface. The v5 shared-state likelihood was the important conceptual advance:
Cull, Manta, SALAD and RHIS need separate observation models around a common
ecological state. Scar persistence/error and Cull depletion cannot be treated
as ordinary habitat responses.

### Density was not the missing fix

SALAD-native density retains more useful abundance information than a binary
threshold and deserves continued use as independent evidence. Cull CPUE does
not show stable within-reef monotonic rank under blocked validation. The
RosettaCOTS calibration is valuable for translating standardised monitoring
estimates with uncertainty, but it cannot by itself reconcile historical Cull
targeting, track-scale SALAD counts, scars and depletion.

## Current operational and development positions

### Operational incumbent

Keep the frozen **legacy clean + year** product. Preserve its exact data,
feature construction, random seed/model object and 2026 prediction entry point.
It remains accepted because domain experts recognise actionable site-scale
patterns, not because later models were forced to reproduce its rank.

### Best static-habitat development line

Retain **Hmsc v4 120 m enhanced with clean multiscale and clean-native
fallbacks** as the best engineered development baseline. The **60 m Hmsc** is
the most relevant alternative for expert review because its forward score is
better and its fitted support is closer to the desired operational scale. Do
not extend either chain until a specific reef-zone review accepts its surface.

### Best modelling framework

Retain **v5 shared state** as the framework to revisit only when better data are
available. It correctly separates habitat, reef-year activation and observation
processes, but the current version is not a better map.

### Density branch

Archive **v6 SALAD-native 60 m** as a research candidate and independent
abundance validation layer. Reject Cull-native and Rosetta-joint density maps as
successors. Do not spend effort adding full-coverage fallback to a rejected
surface.

## Required final expert review

Review Eyrie, Batt, Broadhurst and Bait across three panel families:

1. legacy versus v2/v4 geomorphic maps;
2. 60 m Hmsc and v5 shared-state maps;
3. v6 Cull/SALAD/joint density maps.

For every reef, record named or drawn zones that are plausible, missed,
spuriously prioritised or operationally unusable. Assess broad geomorphic form,
along-reef continuity, inner/outer and windward/leeward contrasts, accessible
search bands and false hotspots. Do not assess colour range, raw calibration or
similarity to legacy rank.

## Harmonised next steps

1. Complete the zone-level expert review before any further fitting.
2. Make an explicit operational choice: retain legacy alone, or retain legacy
   plus one labelled experimental geomorphic layer for field evaluation.
3. If the 60 m or 120 m Hmsc surface is accepted, run production chains only
   for that surface and require ESS >=400, upper PSRF <=1.1, held-cell,
   leave-reef and forward acceptance.
4. If no challenger is accepted, close the current modelling cycle and design
   prospective validation: fixed-effort, pre-cull repeat visits across
   stratified geomorphic zones, with compatible direct COTS and scar recording.
5. Revisit shared-state abundance only with repeat detection data, known search
   area/time, scar-age information and enough new reefs/years for honest forward
   testing. Add independent connectivity/current and temporally prior coral
   condition when available.
6. Keep reef-year outbreak probability separate from within-reef rank. It may
   activate a reef-year habitat surface but must never reorder cells.
7. Any future deployable surface must use the enhanced -> clean multiscale ->
   clean native fallback hierarchy and provide an explicit prediction-source
   layer.

## Reports used in this review

- `reports/HMSC_V2_VISIT_MULTISCALE_ITERATION_02.md`
- `iterations/hmsc_corrected_v4/reports/HMSC_CORRECTED_V4_ITERATION_05_SUPPORT_AND_MANTA.md`
- `iterations/occupancy_detection_v5/reports/OCCUPANCY_DETECTION_V5_ITERATION_06.md`
- `iterations/density_rosetta_v6/reports/DENSITY_ROSETTA_V6_ITERATION_07.md`
