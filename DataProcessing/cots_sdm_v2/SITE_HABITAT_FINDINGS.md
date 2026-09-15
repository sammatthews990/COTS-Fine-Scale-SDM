# Site-habitat ablation findings

## Decision

The fine-scale habitat stage should not contain reef-level outbreak
probabilities, reef identity, or calendar year. Those variables answer which
reef is likely to contain COTS, while this stage must answer where on a selected
reef problematic density is most likely.

The lifetime cull-site model is now a diagnostic ablation rather than the
proposed final habitat model. Collapsing a site's history to "ever high" loses
the repeated annual contrast that appears to be central to the useful legacy
year model and gives more opportunity for survey history to affect the label.
The RG family remains the most promising static predictor family, but it should
be refitted against the annual response before acceptance. Categorical
geomorphic classes are not retained. Site-aggregated manta zeros are not yet
retained in the final high-density model: a very light weight is marginally
promising, but the gain is small and stronger weighting is harmful.

## Controlled experiment

The response is defined once per cull site. A site is positive if any surveyed
year has first-visit, pooled CPUE at or above 0.02. This avoids the legacy
maximum-over-dives response, although an ever positive target can still be
affected by the number of surveyed years.

All eight ablations use the same 4,234 cull sites with complete clean and RG
coverage: 2,118 positive and 2,116 negative sites on 339 reefs. Validation
holds out cull sites within reefs, keeps observations in the same 200 m cell
together, and never uses manta sites as test observations.

| Features | Training evidence | Overall AUC | Within-reef AUC | Between-only AUC | Brier |
|---|---|---:|---:|---:|---:|
| Clean | Cull only | 0.569 | 0.553 | 0.589 | 0.249 |
| Clean + geomorphic | Cull only | 0.570 | 0.554 | 0.588 | 0.249 |
| RG | Cull only | **0.655** | **0.626** | 0.674 | **0.232** |
| RG + geomorphic | Cull only | 0.654 | 0.626 | 0.675 | 0.232 |
| Clean | Cull + manta | 0.577 | 0.552 | 0.592 | 0.270 |
| Clean + geomorphic | Cull + manta | 0.582 | 0.553 | 0.597 | 0.269 |
| RG | Cull + manta | 0.652 | 0.622 | 0.660 | 0.249 |
| RG + geomorphic | Cull + manta | 0.650 | 0.618 | 0.661 | 0.249 |

The default manta rows in the table have total effective weight equal to the
cull-negative weight. They improve the clean model's overall score without
improving within-reef ranking, and they reduce RG performance and calibration.
That is the same scale-mismatch warning seen with reef priors, although much
smaller.

## Manta-zero construction and weight sensitivity

Manta observations are aggregated to 200 m midpoint cells, approximately the
median tow length. A candidate zero requires at least five informative tows in
at least two years, hard coral recorded on every retained tow, and no COTS or
scar evidence. Candidates within 250 m of a positive cull site or positive
manta cell are excluded.

- 7,345 candidate repeated-negative cells met the effort rule.
- 658 conflicted locally with a positive cull site.
- 3,031 conflicted locally with positive manta evidence.
- 4,144 retained cells on 270 reefs also had complete predictor coverage.

The total manta contribution was then capped relative to the number of cull
negatives. With RG, the within-reef AUCs at total weight ratios 0.10, 0.25,
0.50, and 1.00 were 0.630, 0.628, 0.622, and 0.622. The 0.10 result is only
about 0.004 above cull-only RG. This is too small to claim a robust improvement
without repeated-seed or bootstrap uncertainty and stricter survey-method
sensitivity checks.

The practical interpretation is that manta non-detection is useful evidence of
local absence, but it is not automatically a true zero for the cull-derived
ever high-density response. Its most defensible role remains a separate
presence stage or a low-weight prior.

## Geomorphic result

The geomorphic raster was treated as categorical and one-hot encoded, including
an explicit Unknown class. Raw high-density prevalence differs among classes,
but adding the classes changes within-reef AUC by less than 0.002 after the
continuous clean or RG predictors are included. The classes are therefore
descriptive, not an accepted incremental predictor in this model.

## What drives the RG improvement

The leading RG cull-only importance values are wave period (gain 0.165), wave
height (0.135), turbidity (0.065), eastness (0.058), sky-view factor (0.057),
RG slope (0.056), BPI (0.052), VRM (0.050), and RG bathymetry (0.046).
Wave exposure is the clearest additional source of within-reef signal.

Feature importance is not a causal effect. Partial dependence or SHAP curves,
spatially stratified by reef and management region, are still needed to confirm
that the learned wave relationships are ecologically plausible and not edge or
sampling artefacts.

## Visual QA

The regional inspector now contains three additional layers:

- Site habitat: clean, cull only
- Site habitat: RG, cull only
- Site habitat: RG plus light manta

Across the four visual regions, RG and RG plus light manta are almost identical
(raster correlation 0.978-0.984), consistent with the small validation change.
RG is much less correlated with the legacy agnostic map (0.273-0.577) than the
clean refit is (0.581-0.760), so it adds appreciable within-reef spatial
structure rather than simply reproducing the old surface.

The map overlay currently shows 2026 first-visit observations. These are useful
for a quick spatial sanity check, but the small regional sampleespecially
Townsville and Whitsundaymust not replace the cross-validated result.

## RG coverage and fallback

RG has complete values at 4,234 of 4,800 clean-covered cull sites. The 566
clean-only sites occur on 167 reefs and have a slightly lower positive
prevalence (0.447 versus 0.500). This is a deployment concern: an RG-only map
would have gaps at sites and pixels where the clean stack is available.

A controlled five-fold fallback test evaluates every clean-covered cull site.
The model uses RG where all RG layers exist and a clean-only model elsewhere.
It achieves overall AUC 0.650, within-reef AUC 0.620, and Brier 0.233, compared
with 0.572, 0.545, and 0.249 for clean-only across the same 4,800 sites. The
hybrid is therefore the provisional complete-coverage map product. It is now
the default layer in the regional visual inspector.

## Recommendation and next gate

1. Retain RG with clean fallback as a provisional diagnostic map, not the
   accepted habitat product; retain RG cull-only as the common-domain ablation.
2. Keep reef outbreak probability outside this model and combine it only in a
   separate operational reef-selection product.
3. Do not retain geomorphic class in the provisional model.
4. Keep the 0.10 manta version as a sensitivity layer, not the default.
5. Refit on annual site-year outcomes with site-history weighting, then repeat
   the RG versus RG plus light-manta comparison across several seeds and
   stricter manta definitions.
6. Use within-reef-year AUC, percentile contrast, and fixed-budget capture as
   the primary fine-scale tests. Do not assess the model from colour scaling or
   the width of its probability range.
