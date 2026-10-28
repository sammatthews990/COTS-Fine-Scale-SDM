# Table 3 effort and operational cutoff clarification

This amendment is now included in review 09. After the author confirmed it had
no comments, the Word review and both clean baselines were refreshed from the
updated main report on 27 October. Previous copies are archived in
`../../cache/review_09_before_effort_refresh_20261027_093641/`.

“Sites searched” meant designated cull sites selected for a dive visit, including
never-dived candidates. It did not count manta tows or sites at which COTS were
actually removed. The revised label is **sites selected for cull dives**.

The main Results now distinguish two decisions:

- Retain a specified percentage of observed >0.04 CPUE dives: the existing
  80%, 90%, 95% and 99% matched-retention comparisons remain.
- Retain priority classes A–C and defer D: the exact product cutoff now has
  explicit rows in Table 3 and a marked point in the retention figure. Figure 3's
  caption defines the same classes and links to the measured trade-off.

The latter is a fixed site cutoff, rather than an 80% dive-retention target.
Integer rounding means A–C selects floor(0.8 × number of sites) on each reef.
Historical forward rankings retain their original tie order; the deployed
product uses habitat to break ties and its current evidence to rank sites.

Both **whole manta tows** and **tow distance** are reported, as requested, in a
companion table. Whole tows are potentially avoidable only when the complete
recorded segment lies within deferred-only sites. Distance can be deferred
within such sites while retaining the rest of the tow. Portions outside
designated sites or shared with any retained site remain surveyed.

These spatial estimates pair with cull-dive retention conditional on selected
dives still taking place. They do not demonstrate that removing the associated
manta surveillance would leave detection unchanged. Endpoint segments also
cannot establish actual vessel routes or time savings. Annual/pre-closure
checks, additional exploration and their detection yield remain part of the
field assessment.

Reproducible sources: script `27_operational_cutoff_and_tow_effort.R` and tables
`operational_priority_cutoff.csv`, `operational_tow_effort_matched.csv`,
`operational_tow_effort_fixed.csv` and `operational_tow_effort_coverage.csv`.
The script must reproduce every saved forward rank and the existing
matched-retention site point estimates before writing results. Original
threshold-evaluation tables and deployed rankings remain unchanged.

## Retaining A–C and deferring D in the historical test

| Metric | Estimate (%) | 95% interval (%) |
|---|---:|---:|
| Cull sites selected for dives | 78.9 | 78.8–79.0 |
| Dives above 0.04 CPUE retained | 89.1 | 85.4–92.5 |
| Low-density dives avoided | 18.8 | 17.4–20.3 |
| Ghost dives avoided | 23.5 | 21.8–25.3 |
| Whole manta tows potentially avoided | 11.6 | 11.1–12.1 |
| Straight-line manta distance potentially avoided | 16.2 | 15.6–16.9 |

The initially omitted share of observed >0.04 dives is 10.9%. The tow estimates
are gross spatial reductions before annual/pre-closure checks and added exploration.
Random selection of the same number of sites retains 78.5% of >0.04 dives.
These results cover the historical evaluation, rather than forecasting 2027 performance.

Validation: all 240 comparable fixed-effort and 48 matched-retention point estimates
match the original results; every saved forward rank is reproduced exactly.
Tow reductions are monotonic and become zero when all sites are retained.
