# Review 09 responses

Annotated input: `technical_report_review_09.docx` (SHA-256
`71345b264249632e8f92ddd24bf14779e16a705e4423c174f241e41a3d0f94ac`).
It and its original baselines remain unchanged. The main Quarto report is the
authoring source; review 10 is the next clean review copy. Extraction and
pre-integration backups are in `../../cache/review_09_integration_20261027_132711/`
and `../../cache/review_09_feedback.json`.

## Comments

| ID | Feedback | Response |
|---|---|---|
| 17 | Fill low-density and ghost-dive percentages in the summary and elevate them in Results. | Incorporated: 18.8% of low-density and 23.5% of ghost dives avoided at the exact A–C cutoff, read from the existing results table. The Results opening highlights these outcomes. Table 1 now also reports all-dive avoidance, retention of 0.02–0.04 CPUE dives and all four deferred-dive classes with explicit denominators and paired uncertainty. |
| 65 | Add hydrodynamic mesh panel C and the population/connectivity domain as panel D in Figure 1. | Panel D added from the actual v0.4 population-model GeoPackage: 2,895 polygons on 113 reefs, with its checksum pinned in config. Panel C is reserved and explicitly labelled “to be supplied”, as the comment indicates the mesh figure is forthcoming. The caption distinguishes this domain from completed fine-scale transport coupling. The mesh image is the remaining author input for this comment. |

## Tracked passages

- **Summary opening:** added the hydrodynamic/connectivity strand and prototype
  population model alongside the habitat and operational work.
- **Summary habitat/efficiency:** retained the author's focus on subsequent
  action, recent evidence and the approximate top-80% rule. Filled the two XX
  percentages from source CSVs, repaired punctuation and distinguished the
  78.9% actual site count from 89.1% high-density retention. Useful intermediate
  culling opportunities and annual follow-up are explicit.
- **Summary surveillance/exploration:** foregrounded the distance estimate and
  exploration opportunities. The estimate keeps outside-polygon tows, including
  nearby ones; further savings are described as possible and untested rather
  than an established underestimate of achievable operational savings.
- **Summary population modelling:** incorporated ADRIA integration, historical
  irruptions, representation of the emerging fifth wave and exploratory 2050
  scenarios. Retained the density-mismatch and validation qualifications.
- **Summary engagement:** incorporated regular GBRMPA meetings, delivery and
  in-water planning, with completed Malu visits and planned workshop/trials
  distinguished.
- **Methods observation footprint:** clarified the habitat-observation and
  population-calibration context of Figure 1, alongside the domain panel. The
  full archive footprint is not described as exactly the training subset of
  every habitat fit.

## Additional instructions from chat

The main Methods now includes the actual fitted logistic ranking equation.
Appendix G defines all inputs and distance weighting, lists the archived 2027
coefficients, and explains training, absent evidence, logistic transformation,
ties and within-reef priority classes. The ≥0.02 training response is distinct
from the >0.04 evaluation target; no model was refitted.

Table 1 includes recorded-dive avoidance and the deferred-dive breakdown.
Retention of 0.02–0.04 CPUE dives is described as preserving useful culling
opportunities. Methods, Results and Discussion distinguish candidate-site
coverage, recorded dive counts, CPUE-class denominators and deferred composition.
Existing retention results and conditional tow estimates remain unchanged.

The writing brief records these preferences for future Codex and Claude passes.
The supplied CCIP template is unchanged; equations reuse its borderless table
layout with editable Word mathematics and live numbering.
