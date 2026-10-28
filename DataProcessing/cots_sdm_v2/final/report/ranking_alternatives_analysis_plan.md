# Three alternatives to the current cull-site ranking

Specified before fitting or inspecting candidate performance, 27 October 2026.
The existing product and pinned inputs are preserved. Scripts 29–30 write new
results separately; no candidate is promoted solely because it fits better.

## Candidate definitions

All three retain the existing logistic response (a dive with CPUE ≥0.02),
preserving the value of intermediate 0.02–0.04 dives. The primary management
evaluation remains retention of dives strictly above 0.04. Existing quiet-site,
nearby dive, manta and effort predictors remain unless specified below.

1. **High-density indicator:** replace the direct-site prior ≥0.02 indicator
   with any prior-season dive strictly above 0.04. The quiet-site indicator
   still means visited with no ≥0.02 dive, so intermediate finds are not
   reclassified as quiet searches.
2. **Recent CPUE:** candidate 1 plus the temporally weighted site CPUE and
   effective recent dive count. CPUE is weighted catch divided by weighted
   bottom time; records receive weight `2^(-age_days/90)`. Transform CPUE as
   `log(1 + CPUE/0.02)` and evidence amount as `log(1 + sum(weights))`.
   The fixed 90-day half-life is not selected using test-year outcomes.
3. **Habitat hotspots:** candidate 2 with a site's within-polygon 90th-percentile
   habitat prediction replacing its mean, then ranked among sites on its reef.
   The 90th percentile is fixed before assessment. Habitat-only diagnostics
   also compare mean, median, 75th, 90th and 95th percentiles; they do not select
   a winner retrospectively for candidate 3.

## Controls and habitat support

Retain the archived current ranking as the incumbent. A second current-model
control uses the same complete 30 m cell extraction as the new candidates,
with the original ≥0.02 indicator and the mean habitat prediction. This
separates denser polygon sampling from changing the indicator or summary.
The original spatial evaluation used about 16 regular points per polygon.
All new summaries use the same cells intersecting each polygon, with a
representative-point fallback if needed. Coverage and cell counts are reported.

Refit only the cull habitat model, with the existing hyperparameters and the
same five reef folds used in script 17. Check its predictions at the archived
sample points against saved held-out values before computing full-cell
summaries. No test reef contributes to its own habitat fit. As in the current
analysis, habitat fits are spatially held out, not restricted to earlier years;
the evidence regression and its predictors are strictly forward in time.

## Validation and decision criteria

For each test year 2019–2026, train the evidence regression only on earlier
years. Build each year's evidence from its preceding season; no current-year
observations enter its predictors. Rank all designated sites, including
never-dived ones, on the same evaluated reef-years as the incumbent. Report
2024–2026 separately to check whether overall gains persist recently.

Report within reef-year discrimination and high-density top-20% lift, the exact
A–C cutoff, and matched 80%, 90%, 95% and 99% high-density retention. Separate
candidate sites, recorded dive counts and recorded bottom time. Also compare
high-density retention at the same expected recorded dive/bottom-time budget.
Track intermediate retention, useful-dive retention (≥0.02), low-density and
ghost avoidance, and the composition of deferred dives. Bottom time is an
observed effort proxy, not a measured saving in search or vessel time.
Trade-off curves use integer top-site counts on a fixed 1% cutoff grid and
first-crossing interpolation. All models use that same grid; the exact A–C
point must reproduce the published incumbent, independently of interpolation.

Uncertainty uses 500 paired resamples of whole reef-years and paired differences
from both controls. Do not attribute raw reductions to ranking skill. A strong
candidate should improve the primary matched-90% trade-off by at least one
percentage point, have a paired interval supporting improvement, and retain
the direction of gain in 2024–2026. Report any loss of intermediate/useful dives
or contrasting metrics explicitly. Smaller or uncertain gains remain
experimental options rather than a replacement recommendation.

If supported, export a separate candidate site layer and coefficient table,
preserving the incumbent and labelling its status. Document the forward-test
choice before fitting a final all-years model or using the latest evidence.
