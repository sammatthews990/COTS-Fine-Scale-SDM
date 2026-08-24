# Initial benchmark findings

## Scope

These results are an ablation benchmark, not a final SDM. The four models use
the same 9,936 complete site-year rows, first-visit pooled-CPUE target, folds,
XGBoost settings, and static predictors. Only the inclusion of calendar year
and lagged reef outbreak history changes.

The current-year outbreak probability is excluded. Outbreak inputs are limited
to lags 1-3, their mean, maximum, and trend.

## Controlled benchmark

| Features | Unseen reef-cluster ROC-AUC | Unseen reef-cluster PR-AUC | Unseen reef-cluster Brier | Forward-2026 ROC-AUC | Forward-2026 PR-AUC | Forward-2026 Brier |
|---|---:|---:|---:|---:|---:|---:|
| Static | 0.512 | 0.449 | 0.257 | 0.705 | 0.663 | 0.220 |
| Static + year | 0.531 | 0.497 | 0.263 | 0.709 | 0.681 | 0.219 |
| Static + outbreak history | 0.653 | 0.577 | 0.244 | 0.745 | 0.657 | 0.207 |
| Static + year + outbreak history | 0.648 | 0.578 | 0.247 | 0.734 | 0.673 | 0.209 |

Interpretation:

- The current static predictor/point-extraction combination does not transfer
  reliably to unseen reef clusters. It should not yet be labelled a stable
  habitat prior.
- Calendar year provides only a small spatial-transfer uplift when the target,
  rows, folds, and fitting controls are held constant. Its large importance in
  the legacy year model is therefore not evidence of habitat suitability.
- Lagged outbreak history supplies the strongest transferable signal in this
  benchmark. It represents dynamic reef state, not fine-scale habitat.
- Forward 2026 is a different task from unseen-reef validation. Many test reefs
  and sites have earlier observations, so static predictors can rank within a
  familiar monitoring domain even though they transfer poorly to new reefs.
- With the adjusted hindcasts, outbreak history also materially improves the
  aggregate 2026 result. Its benefit is still judged over several forecast
  origins rather than from one year.

The deterministic geographic folds are currently imbalanced (181-2,689
complete benchmark rows per held-out fold). That should be improved before a
final model comparison. It does not explain the main result: static fold AUCs
are consistently 0.47-0.55, whereas static plus outbreak-history fold AUCs are
0.64-0.72.

## Rolling-origin benchmark, 2022-2026

| Features | ROC-AUC | PR-AUC | Brier | Log loss | Capture in top 20% |
|---|---:|---:|---:|---:|---:|
| Static | 0.596 | 0.463 | 0.238 | 0.668 | 25.7% |
| Static + year | 0.590 | 0.454 | 0.234 | 0.660 | 24.9% |
| Static + outbreak history | 0.696 | 0.540 | 0.213 | 0.615 | 31.5% |
| Static + year + outbreak history | 0.685 | 0.541 | 0.214 | 0.617 | 30.5% |

The adjusted outbreak-history gain is real enough to continue developing, but
it varies by year. It is strongest in 2023 and 2026, useful in 2022 and 2024,
and modest in 2025. Calendar year again adds no pooled rolling benefit.

The rolling output also tags each prediction as a previously seen or new reef.
New-reef sample sizes are small (48-215 per test year), so those subgroup
metrics are diagnostic rather than definitive.

## Target construction diagnostics

- There are 10,028 site-years; 4,189 contain more than one survey date.
- Restricting the target to the first visit changes 442 site-year labels (4.4%).
- Using any single-dive maximum instead of pooled CPUE changes 654 labels
  (6.5%); every change is a max-positive but pooled-negative case.
- The probability of ever exceeding 0.02 on a single dive rises from 15.6% for
  sites with one lifetime dive to 80.6% at 6-10 dives, 99.0% at 21-50 dives,
  and 100% above 50 dives.

This confirms that an unadjusted historical maximum is mainly a cumulative
sampling opportunity measure. The static habitat stage needs either an
effort-aware repeated-survey model or a history statistic with explicit
denominator and uncertainty.

## Outbreak-layer status

The configured adjusted hindcast contains 135,135 reef-years for 3,861 reefs
from 1991-2025. The separate adjusted forecast contains all 3,861 reefs for
2026 and supplies a 0.45 outbreak cutoff. Reef names exactly match 359 of the
361 cull reef names. The unmatched labels are `Intereefal Sites` and
`Joe Baker Reef (20-374)`.

The provisional operational bridge trains on same-year adjusted hindcasts and
deploys the 2026 forecast. On the 615 common 2026 rows:

| Operational features | ROC-AUC | PR-AUC | Brier | Top-20% capture |
|---|---:|---:|---:|---:|
| Raw reef forecast | 0.658 | 0.570 | 0.415 | 27.2% |
| Static only | 0.705 | 0.663 | 0.220 | 30.7% |
| Static + adjusted outbreak | 0.740 | 0.691 | 0.209 | 34.6% |
| Static + year + adjusted outbreak | 0.740 | 0.712 | 0.209 | 34.3% |

The raw reef forecast ranks reasonably but is not calibrated to site-level
high-density CPUE, which is a different outcome. The fitted bridge is therefore
the appropriate architecture, subject to the provenance caveat below.

Before same-year outbreak probability is allowed into a forecast or validation
model, its upstream provenance must identify it as a genuine pre-survey forecast
or an out-of-fold hindcast. Lagging prevents direct same-year leakage but cannot
repair leakage inside the upstream outbreak model.

## Manta negative-evidence status

The long-history manta file contains 343,940 usable observations across 822
reefs. A conservative rule requires every informative tow at a reef to record
zero COTS and absent feeding scars, with at least 10 tows in at least three
survey years. It identifies:

- 23 repeatedly negative manta reefs supported by 4,294 tows;
- 3,649 fine-scale-eligible tow midpoints after coordinate, coral-presence, and
  tow-distance checks;
- 15 reefs with manta-only, uncontradicted evidence;
- 3 reefs whose cull observations are also all zero; and
- 5 reefs with positive cull dives elsewhere or at other times.

The five conflicts demonstrate why manta non-detection is not automatically a
reef-wide true absence. The prepared data retain the method conflict and are
not appended to the high-density response. Manta negatives belong in the
presence stage with survey method and effort represented explicitly. Only the
18 uncontradicted reefs are flagged as static-prior negative candidates, and
the three cross-method-supported reefs form the strongest subset.

## Current decisions

1. Do not revive the lifetime `CPUE_max >= 0.02` agnostic target.
2. Keep habitat prior and outbreak state as separate named products.
3. Continue with lagged outbreak history for leakage-safe evaluation and the
   adjusted same-year hindcast/forecast bridge for provisional operations.
4. Do not use calendar year as a substitute for a dynamic ecological layer.
5. Build the next static model around repeat-survey effort and hierarchical
   site/reef structure, then validate it on unseen reefs.
6. Add manta evidence to the presence stage only, with survey-method effects;
   do not turn repeated tow rows into unweighted high-density pseudo-absences.
