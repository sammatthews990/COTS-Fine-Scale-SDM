# Visual comparison and prediction-scale findings

## What was compared

The visual inspector uses the same Lizard, Cairns, Townsville, and Whitsunday
windows as the legacy master map. All probability layers use the same default
0.20-0.80 display scale and the difference layers use a fixed -0.30 to +0.30
scale.

The inspector contains:

- the legacy 2026 year model;
- the legacy agnostic model;
- the v2 static-only model;
- the adjusted 2026 reef outbreak forecast;
- the unconstrained v2 static-plus-outbreak XGBoost model;
- an experimental additive, rank-preserving bridge;
- operational-minus-static differences; and
- 2026 first-visit observations and repeated manta-negative reefs.

Predictors are averaged from 30 m to 150 m for visual QA only. Fitted models,
site-level validation, and the source rasters are unchanged.

## Visual findings

The legacy agnostic surface is effectively flat. Across the four windows its
regional 5th-95th percentile range is approximately 0.499-0.521. It therefore
cannot provide a useful spatial ranking even where its absolute probabilities
look plausible. This agrees with the earlier finding that a lifetime maximum
target is dominated by cumulative survey opportunity.

The v2 static surface contains substantially more within-reef structure. Its
regional 5th-95th percentile range is approximately 0.316-0.637, compared with
about 0.394-0.507 for the legacy year layer. This extra contrast is not, by
itself, proof of transferable habitat suitability: unseen-reef spatial
validation for the static model remains weak.

The unconstrained operational surface visibly retains fine-scale structure but
is strongly shifted by reef-scale forecast probability. Its cell-level
correlation with the outbreak layer is 0.68-0.73 across the four windows,
whereas correlation with the v2 static layer is only 0.25-0.40. Reef-shaped
steps are expected because one forecast value is assigned to each reef before
being combined with 30 m predictors.

At mapped 2026 observations, mean v2 static prediction is higher at
high-density sites than below-threshold sites in all four windows. The same is
true for both operational candidates. The legacy year layer separates the two
classes in Lizard, Cairns, and Whitsunday, but reverses slightly in the small
Townsville subset. The legacy agnostic difference is negligible everywhere.

## Why the overall score improvement needs decomposition

The 615 common 2026 site-years cover 82 reefs; 41 reefs contain both outcome
classes. A pair-weighted AUC within those mixed reefs gives a direct answer to
the question: given that the survey team is on this reef, can the model rank
the high-density sites above the low-density sites?

| Model | Overall ROC-AUC | Between-reef-only ROC-AUC | Within-reef ROC-AUC |
|---|---:|---:|---:|
| Raw reef forecast | 0.658 | 0.658 | 0.500 |
| Static only | 0.705 | 0.740 | 0.652 |
| Static + outbreak XGBoost | 0.740 | 0.747 | 0.610 |
| Static + year + outbreak XGBoost | 0.740 | 0.740 | 0.638 |
| Additive rank-preserving bridge | 0.726 | 0.720 | 0.652 |

The unconstrained static-plus-outbreak model has the best overall ranking, but
its within-reef AUC is lower than static-only. Its improvement is therefore a
reef-selection improvement, not a better fine-scale habitat ranking. Allowing
XGBoost to interact a reef-constant prior with static predictors can change
pixel ordering inside a reef.

The experimental additive bridge was fitted to out-of-time 2022-2025 static
predictions and same-year adjusted hindcasts. It combines the two scores on the
logit scale with positive coefficients, so it cannot reverse the static
ranking within a reef. It raises overall AUC from 0.705 to 0.726 and preserves
within-reef AUC at 0.652, but it has worse 2026 PR-AUC and Brier score than the
unconstrained model. It is a useful architectural comparison, not yet the
selected production model.

## Important forecast-scale check

The adjusted 2026 forecast is strongly concentrated near zero:

- median 0.031;
- 95th percentile 0.116;
- 99th percentile 0.203; and
- only 6 of 3,861 reefs at or above the supplied 0.45 cutoff.

The forecast still ranks observed 2026 sites above chance, but raw forecast
probabilities are not calibrated to the site-level high-density target. This
explains why the additive surface is strongly depressed on the shared map
scale. The hindcast/forecast provenance and the meaning of the 0.45 cutoff
must be confirmed before probability magnitudes are interpreted operationally.

## Decision and next implementation steps

No single combined raster should currently be called the habitat map.

1. Keep three separately named products: static habitat candidate, reef-year
   outbreak prior, and combined operational survey score.
2. Audit whether adjusted historical outbreak values are genuine forecasts,
   out-of-fold hindcasts, or in-sample fitted values.
3. Fit an explicit two-stage model:
   - reef/presence stage, using outbreak history, survey effort, and
     conservative manta negative evidence; and
   - within-reef high-density stage, using static fine-scale predictors and
     conditioning on evidence that COTS are present or the reef is active.
4. Replace arbitrary cross-scale XGBoost interactions with either an explicit
   probability product or a constrained monotone/additive bridge.
5. Improve the static stage with repeat-survey effort, a count model with
   bottom-time offset, site/reef structure, and polygon or tow-track predictor
   summaries.
6. Report overall, between-reef, and within-reef metrics in every rolling-year
   and unseen-reef validation. A combined score must improve reef selection
   without silently degrading the conditional within-reef ranking.
7. Continue to show uncertainty, environmental extrapolation, and the reef
   prior beside any operational surface.

The additive bridge and scale-decomposition outputs are generated by
scripts/08_fit_additive_bridge.R. The matched-scale maps and summaries are
generated by scripts/07_build_visual_comparison.R.
