# AGENTS.md

Guidance for coding agents (and people) working in this repository.

## Goal

Predict **where on a reef crown-of-thorns starfish (COTS) are most likely to
be found**, so the COTS Control Program can:

1. **Manta tow less** without missing COTS aggregations.
2. **Search fewer cull sites** without missing problematic densities
   (CPUE ≥ 0.02 COTS per minute of bottom time). This means avoiding dives
   where COTS are absent (ghost dives) or present below the action threshold
   (low-density dives).

Every model and analysis is judged against these two uses. A model is useful
only if it beats a no-information baseline on held-out data at a fixed level
of effort.

## Current state (2026-10-07)

- **Habitat layers rank locations within a reef only weakly.** Held-out
  within reef-year AUCs are about 0.55–0.61.
- **Cull-site selection.** The gain comes from the previous season's cull
  dives and manta tows, combined with a habitat rank (scripts 16–17). For
  the same problematic dives, this searches 6–13 points fewer designated
  sites and avoids more low-density and ghost dives.
- **Manta tows.**
  - Thinning tows evenly (systematic) finds more aggregations than dropping
    the lowest-ranked tows by any model (script 14).
  - A manta-trained habitat model fitted within reef-years ranks tows best
    (AUC 0.61). It does not improve cull-site selection (script 17).
- **Operational product.** `final/outputs/site_ranking/` holds the 2027
  priority ranking of every designated cull site (CSV and KMZ, script 20).
  It uses the cull habitat rank plus the last 12 months of dives and tows.
  - Model choices are in `config.R` → `site_ranking`, chosen by the
    year-by-year forward test (script 19).
  - Ensembles with the manta model added nothing.
  - Re-run script 20 each season with an updated `evidence_end`, then 21.
  - Script 21 adds exploration areas to the same KMZ: high-ranked habitat
    outside cull sites and away from recent tows. Rankings learned from past
    effort otherwise keep consolidating search on the same places.
- Details, numbers and caveats are in `DataProcessing/cots_sdm_v2/final/README.md`.

## Where things are

- `DataProcessing/cots_sdm_v2/final/` holds the delivered pipeline: layers,
  validation and report figures. **Read `final/README.md` first.** It records
  every result, caveat and the run order.
- `DataProcessing/cots_sdm_v2/final/report/results.qmd` is the Results section
  of the technical report.
- Everything else under `DataProcessing/cots_sdm_v2/` (iterations, artifacts,
  reports) is exploratory history. It is due to be moved to
  `cots_sdm_v2/archive/`; the user will do this move. The pipeline's
  `config.R` finds inputs in either location.
- `DataProcessing/SDM_Validation_Report.qmd` and `SDMReefValidation.qmd` are
  superseded earlier estimates. They had no random baseline and used in-sample
  scoring.
- Raw data are in `data/`. This includes the EoTR designated cull-site KMZ and
  the ReefScan detections (presence-only, no effort tracks).

## Rules for analyses

- **Frozen inputs.** `final/config.R` pins every input by SHA-256. Do not
  regenerate or replace a pinned file silently. Add new inputs to `config.R`
  with their checksum.
- **Rank within reefs.** Operational decisions are made within a reef, so
  report within reef-year metrics: AUC, top-20 % lift, share kept at a fixed
  effort. Survey year and reef identity dominate raw probabilities.
- **Held out, never in-sample.**
  - The cull-trained layers were fitted on Cull dives, and Cull crews dive
    where manta tows found COTS. Score each reef with a model refitted without
    that reef's data (leave-reef-out), as in scripts 14–17.
  - Validation that uses history must be strictly forward: only data from
    years before the test year.
- **Always compare with a baseline at the same effort.** Use random selection
  (exact expectation) and, for manta tows, systematic thinning. Any cut saves
  effort; only the difference from the baseline is the model's contribution.
- **Matched retention.** When a score moves problematic and low-density sites
  up together, compare scores at the same share of problematic dives kept.
- **Say when a score encodes operator behaviour.** Cull history and dive
  effort reflect where crews chose to dive. Report them as such, not as
  habitat evidence.
- **Designated cull sites are the candidates.** Cull-site selection ranks all
  designated sites on a reef (the KMZ), including never-dived ones, not only
  places that were dived.

## Code conventions

- R scripts, run from the repository root:
  `Rscript DataProcessing/cots_sdm_v2/final/scripts/NN_name.R`.
  Each script starts with a header comment stating its question, method and
  validation.
- Scripts source `final/R/common.R` and call `initialise_final()`. Shared
  helpers live in `final/R/`; for example, `cull_sites.R` provides designated
  sites and matched dives.
- Write every reported number to `final/outputs/tables/` with `write_table()`.
  Put figures in `07_build_report_figures.R`, which runs last, using the
  palette in `R/figures.R`. Intermediate results go in `final/cache/`
  (git-ignored).
- When a result changes, update `final/README.md` (results, figure list, run
  order) and, where relevant, `report/results.qmd`.

## Git

- Work on branch `tidy/final-pipeline`. Commit or push only when asked.
