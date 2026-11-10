# Review 11: three alternative cull-site rankings

This version adds a new analysis to the main report. It preserves review 10
and all earlier feedback/baselines; no new review 10 comments were requested
or integrated. The supplied CCIP template remains unchanged.

## What was assessed

1. Replace the prior-site ≥0.02 indicator with any previous-season dive >0.04.
2. Add recency-weighted site CPUE and its effective record count, with a fixed
   90-day half-life.
3. Use the within-site 90th-percentile habitat value instead of its mean in
   the recent-CPUE version.

The ≥0.02 fitted response remains unchanged so intermediate dives retain value.
All three were fitted forward in time, with the same reef-held-out habitat
folds, candidate sites and observed dive cohort. Habitat summaries used complete
30 m cells and a mean control at that same support. Median, 75th, 90th and
95th percentile habitat-only diagnostics were also assessed without choosing
a new combined model retrospectively from their results.

## Findings and decision

At the same A–C site coverage, the current ranking retains 89.1% of >0.04 dives;
the three alternatives retain 89.5%, 90.5% and 90.0%. At 90% high-density
retention, their recorded-dive avoidance is 13.9%, 14.7% and 14.3%, versus
13.8% for the incumbent in this new comparison. The recent-CPUE gain is
0.89 percentage points [paired 95% interval −0.08 to 1.97]. At the same actual
recorded-dive budget, its high-density-retention gain is smaller: 0.54 points
[−0.05 to 1.31]. Intermediate-dive trade-offs and recent-period results are
reported explicitly in Appendix M.

No alternative passed the pre-specified strong replacement criterion.
Recent CPUE is supplied as a separately labelled experimental option for
field comparison. The current ranking and its CSV/KMZ remain unchanged.
The high-percentile hypothesis was not supported by the 90th-percentile
combined model; median habitat was more promising as a diagnostic for a later
test. The result does not rule out real habitat hotspots.

The new matched-retention calculations use 1% integer-site steps and new common
paired bootstrap draws. Earlier report tables use the original 5% percentile
steps. Small changes in interpolated incumbent values or Monte Carlo intervals
come from that comparison protocol; its exact A–C point estimates reproduce
the published values. Existing Table 1 and Appendix H results are preserved.

## Report changes and products

- Summary, Methods, Results and Discussion identify the experimental comparison
  without replacing the existing operational product.
- Appendix M gives definitions, the recency equation, fixed-site and matched
  retention/effort tables, two figures, uncertainty and the decision.
- Appendix A adds the experimental delivery. The output is in
  `../../outputs/site_ranking_candidates/` as CSV, GeoPackage and KMZ, with
  checksums in `ranking_candidate_product_manifest.csv`.
- The writing brief records the interpretation for future Codex/Claude passes.

The original review 09 hydrodynamic mesh request remains pending its supplied
image. Other author-completion and model-validation flags are preserved.

## Export checks

The main DOCX and PDF were rendered through the CCIP template workflow and
passed its checks: 108 pages, no blank pages or field errors, preserved critical
styles and branding, and three editable display equations. Page images were
inspected, including all new equations, tables and figures. The three short
comparison tables remain together with their captions. The source template and
annotated review 09 have unchanged checksums. Review 11 and its DOCX/QMD
baselines match the final main source/export; review 10 is preserved.
