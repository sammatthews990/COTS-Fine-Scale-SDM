# Current best COTS fine-scale model

Last evaluated: 2026-09-09

## Current decision

**Legacy clean + year remains the incumbent.**

The corrected v4 work has now tested active-reef-year conditioning, strong
Manta evidence, multiscale terrain, 120 m and 200 m fitted support, method-
specific Hmsc response blocks, and explicit coverage fallbacks. No candidate
passes convergence, forward transfer, independent transfer and expert-map
acceptance together.

## Current development baseline

The best v4 development surface remains the 120 m geomorphology-led enhanced
Hmsc output with a fixed coverage hierarchy:

1. enhanced geomorphology, exposure, Acropora and multiscale terrain;
2. clean Acropora plus multiscale terrain;
3. clean native stack;
4. explicit missing where the source stack itself is absent.

Fallbacks are calibrated only for seam consistency and never overwrite an
available enhanced prediction. Reef identity, year, outbreak probability,
survey effort and method are excluded from the deployable habitat score.

| Candidate | Held within-reef AUC | Leave-reef within AUC | Forward within AUC | Decision |
|---|---:|---:|---:|---|
| v4 enhanced 120 m | 0.569-0.570 | 0.545 | 0.515 | not promoted |
| v4 enhanced 200 m PCA2 | 0.541 | 0.528 | 0.528 | reject support change |
| v4 120 m + strong Manta response | 0.569 | 0.546 | 0.515 | no material change |

The 200 m version is slightly better forward but substantially worse on held
cells and unseen reefs. Its enhanced coverage is only 35-72% on the four review
reefs, so it cannot replace the 120 m fallback architecture.

Strong Manta evidence is useful independently: whole-tow enhanced habitat
achieves within-reef AUC 0.560 held, 0.544 leave-reef and 0.576 forward. However,
adding it as a simple Hmsc response changes the Cull map almost not at all
(four-reef rank correlations 0.9986-0.9994) and worsens convergence. Retain it
as a separate evidence/activation layer until a formal shared-state,
method-specific detection model is fitted.

Detailed methods and results:

- [Iteration 04: corrected conditional 120 m model](../iterations/hmsc_corrected_v4/reports/HMSC_CORRECTED_V4_ITERATION_04.md)
- [Iteration 05: 200 m support and strong-Manta Hmsc test](../iterations/hmsc_corrected_v4/reports/HMSC_CORRECTED_V4_ITERATION_05_SUPPORT_AND_MANTA.md)

## Required expert review

### 120 m versus 200 m

- [Eyrie](../iterations/hmsc_corrected_v4/reports/figures/v4_support_comparison_eyrie_120m_200m.png)
- [Batt](../iterations/hmsc_corrected_v4/reports/figures/v4_support_comparison_batt_120m_200m.png)
- [Broadhurst](../iterations/hmsc_corrected_v4/reports/figures/v4_support_comparison_broadhurst_120m_200m.png)
- [Bait](../iterations/hmsc_corrected_v4/reports/figures/v4_support_comparison_bait_120m_200m.png)

### COTS-only versus strong-Manta Hmsc

- [Eyrie](../iterations/hmsc_corrected_v4/reports/figures/v4_manta_hmsc_comparison_eyrie.png)
- [Batt](../iterations/hmsc_corrected_v4/reports/figures/v4_manta_hmsc_comparison_batt.png)
- [Broadhurst](../iterations/hmsc_corrected_v4/reports/figures/v4_manta_hmsc_comparison_broadhurst.png)
- [Bait](../iterations/hmsc_corrected_v4/reports/figures/v4_manta_hmsc_comparison_bait.png)

## Harmonised next steps

1. Record specific expert-accepted and rejected reef zones from both panel sets.
2. Retain 120 m and test a matched 60 m model; do not pursue 200 m further
   unless expert review reveals an operational advantage.
3. Implement a shared latent COTS state with separate Cull, Manta, SALAD and
   RHIS detection models, scar persistence/error, and Cull depletion.
4. Keep strong Manta positives as direct COTS, >=2 counted scars, or legacy
   common scars; retain single/present scars only as sensitivity evidence.
5. Add independent connectivity/exposure and compatible susceptible-coral
   condition to improve forward transfer while keeping outbreak activation
   separate from within-reef rank.
6. Extend Hmsc chains only for an expert-accepted surface and require ESS >=400,
   upper PSRF <=1.1, held-cell, leave-reef and forward-year acceptance.
