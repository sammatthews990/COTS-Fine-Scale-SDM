# Continuous rank-signal findings

## What the 0.49-0.53 pattern means

A narrow probability interval can contain a fully resolved and operationally
useful ordering. The assessment therefore uses every prediction at full numeric
precision. It does not bin values, judge contrast from a colour ramp, or reward
a model merely for producing a wider range.

As an explicit check, each model's predictions were replaced by exact
within-reef-year percentile ranks. The within-reef-year AUC was unchanged for
all four annual models. The same invariance held for every 2026 comparison
layer. A monotonic display stretch can make an ordering easier to see, but it
cannot improve these rank metrics.

The fine-scale question is now stated as:

> Among sites surveyed on the same reef in the same year, does the model rank
> sites with CPUE at or above 0.02 above the other surveyed sites?

## Why the legacy year and agnostic models differ

They are not the same habitat model with and without a year column.

- The legacy year model keeps one pooled response per site and year and labels
  annual pooled CPUE at or above 0.02.
- The legacy agnostic model collapses the full record at a site and, by
  default, labels it from the maximum individual-survey CPUE.
- The first v2 site-habitat model instead labels a site if any first-visit year
  was above the threshold.

The two history-collapsed responses discard the sequence of positive and
negative years. They also make a site's label more likely to become positive as
survey history accumulates. This target difference is a more important
explanation than the presence of the Year predictor alone.

For a map predicted at one fixed year, Year has the same value at every pixel.
It cannot directly create within-reef spatial variation. It can still change
the fitted tree structure through interactions during training, and it can
absorb large temporal changes in prevalence. Those effects can let the static
predictors learn a cleaner spatial ordering, but they must be separated from a
simple between-year performance gain.

## Held-out annual-response benchmark

The benchmark contains 8,801 annual first-visit site records from 4,234 sites
with complete clean and RG predictors. All years from one site remain in one
validation fold. Training rows are weighted by the reciprocal of the number of
years observed at that site. The primary comparison is within 448 mixed
reef-years containing both positive and negative sites.

| Model | Overall AUC | Within reef-year AUC | Top-20% capture | Top-20% lift | Positive-negative percentile gap |
|---|---:|---:|---:|---:|---:|
| Annual clean, static | 0.558 | 0.546 | 0.254 | 1.108 | 0.037 |
| Annual clean + year | 0.643 | 0.548 | 0.248 | 1.080 | 0.036 |
| Annual RG, static | 0.624 | **0.592** | **0.263** | **1.147** | **0.063** |
| Annual RG + year | **0.691** | 0.584 | 0.261 | 1.136 | 0.056 |

Year produces a large overall AUC gain, but essentially no clean-model gain and
a modest RG loss within the same reef-year. That is evidence that most of the
year gain concerns temporal or between-group discrimination, not improved
fine-scale habitat ordering. RG contributes genuine fine-scale signal, although
the current score is not yet strong enough to accept the new map on metrics
alone.

## Fixed legacy-spec held-out check

The saved legacy clean-year fit used 1,500 trees, depth 5, learning rate
0.00137, minimum child weight 20, and effectively zero loss reduction. Those
fixed settings were reapplied to the legacy annual pooled response, with all
years from each site held in the same validation fold. A paired model omitted
only Year.

| Model | Overall AUC | Within reef-year AUC | Top-20% capture | Top-20% lift | Percentile gap |
|---|---:|---:|---:|---:|---:|
| Legacy settings, annual clean static | 0.570 | 0.555 | 0.253 | 1.104 | 0.042 |
| Legacy settings, annual clean + year | 0.662 | 0.568 | 0.251 | 1.095 | 0.037 |

This check uses 9,940 annual pooled records at 4,800 held-out sites. Year again
has a large effect on overall discrimination, but only a small effect on
within-reef-year AUC and a slight negative effect on fixed-budget capture and
percentile separation. The useful spatial pattern is therefore not explained
by the Year column itself.

The legacy hyperparameters were originally selected using the full legacy
dataset, rather than selected independently inside every training fold. This
fixed-parameter transfer is still a stronger test than evaluating the fitted
raster at its own training sites, but a fully nested retune would be required
for a final unbiased performance estimate.

## Descriptive 2026 comparison

At 615 sites observed in 2026, the extracted legacy year raster agrees much
more strongly with the observed within-reef ordering than the legacy agnostic
raster or the current v2 operational static model.

| Layer | Within-reef AUC | Top-20% capture | Top-20% lift | Percentile gap |
|---|---:|---:|---:|---:|
| Legacy agnostic raster | 0.577 | 0.249 | 1.073 | 0.042 |
| Legacy year raster | **0.730** | **0.318** | **1.373** | **0.163** |
| v2 operational static | 0.652 | 0.284 | 1.223 | 0.098 |
| v2 static + outbreak | 0.610 | 0.269 | 1.159 | 0.069 |
| v2 static + year + outbreak | 0.638 | 0.274 | 1.180 | 0.084 |

This supports the operational observation that the legacy year surface
contains meaningful fine-scale ordering. It does not yet establish an unbiased
performance estimate: the legacy raster was fitted from the current workbook,
which contains 2026 observations, so these extracted predictions may be partly
in sample. The comparison is deliberately labelled descriptive in the output.

The gap between the legacy year's descriptive within-reef AUC (0.730) and its
held-out-site legacy-spec result (0.568) is too large to ignore. It is consistent
with substantial in-sample optimism or learned local detail that does not
transfer to unsurveyed sites. It does not prove that the expert-recognised
pattern is ecologically wrong; it means that pattern needs independent evidence
before it is promoted as the habitat model.

The outbreak layers again improve the broad overall score but weaken the
within-reef ordering. They should remain outside the fine-scale habitat score.

## Modelling decision and next experiment

The next candidate should preserve the annual response while separating two
components:

1. a static spatial score used to rank sites within a reef; and
2. a year-specific prevalence/calibration adjustment that cannot reorder sites
   on the same reef.

The fixed legacy year formulation has now been reproduced under held-out-site
validation. Its large descriptive advantage did not survive. The next refit
should use a static RG spatial score and a separately estimated year prevalence
adjustment, followed by a nested or pre-registered tuning comparison. The
legacy raster should remain in the blind expert review as a benchmark, not as
the supervised target for the replacement model.

## Rank-preserving year calibration

A year-specific intercept calibrator has now been added on top of the static
clean and RG scores. Each held-out fold learns its year shifts from training
sites only. The transformation adds a constant to the log-odds for all sites in
the same year, so a deployed model cannot reorder sites within a reef-year.

| Spatial score | Calibration | Overall AUC | Brier | Log loss |
|---|---|---:|---:|---:|
| Clean static | None | 0.558 | 0.251 | 0.698 |
| Clean static | Year intercept | 0.639 | 0.238 | 0.671 |
| RG static | None | 0.624 | 0.239 | 0.671 |
| RG static | Year intercept | **0.684** | **0.225** | **0.642** |

Within each model fold and reef-year, the calibrated and uncalibrated AUCs are
exactly identical: 0.561 for clean and 0.609 for RG. When all out-of-fold
predictions are pooled, very small differences appear because sites on one reef
were predicted by different fold-specific calibrators. That is a cross-fitting
artefact, not a change that the final single model would make to within-year
rankings.

This architecture therefore gives the useful broad temporal adjustment without
allowing Year to create or erase fine-scale spatial structure. A forecast-year
intercept still needs an external, prospectively available prevalence estimate;
it must not be estimated from target-year COTS observations.

## Rank stability across site splits

The annual clean and RG static models were evaluated across three independent
within-reef site-fold assignments.

| Model | Mean within reef-year AUC | SD | Range | Mean top-20% lift |
|---|---:|---:|---:|---:|
| Annual clean static | 0.541 | 0.011 | 0.529-0.548 | 1.085 |
| Annual RG static | **0.590** | **0.004** | **0.586-0.593** | **1.129** |

RG is better on every split and is less sensitive to the fold assignment. This
makes RG the preferred annual spatial score for the next tuning gate, while the
clean model remains necessary as a coverage fallback.

## Blind spatial-plausibility review

`outputs/expert_blind_review.html` contains five anonymous surfaces across the
Lizard, Cairns, Townsville, and Whitsunday regions. All surfaces are converted
to regional percentile ranks and shown on one fixed display scale. Model names,
observations, and outbreak probabilities are absent.

The review saves four 1-5 operational scores, confidence, and notes for every
region-surface combination and exports a CSV. The identity key is a separate
file and should remain hidden until independent scoring is complete. The
protocol is documented in `EXPERT_REVIEW_PROTOCOL.md`.

Acceptance should require all of the following:

- better held-out within-reef-year pairwise concordance;
- better high-density capture in a fixed top survey fraction;
- stable rankings across folds, seeds, and forecast cut-off years;
- a prespecified blind review on reefs with strong operational knowledge; and
- separate calibration checks when the output is presented as a probability.

The expert review should compare anonymous within-reef rank surfaces, not ask
reviewers to score colour intensity or probability spread.
