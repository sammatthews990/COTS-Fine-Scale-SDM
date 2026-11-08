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
  - Script 21 adds exploration zones to the same KMZ: all high-ranked
    habitat outside cull sites and away from recent tows, as nested shaded
    bands (top 40/30/20/10 %). They replaced the top-10 discrete areas on
    6 November: nothing set those areas apart from the habitat around them.
    Rankings learned from past effort otherwise keep consolidating search on
    the same places.
- **Density-aware classes (script 20b, 6 November).** Within-reef classes
  give every reef 20 % A sites whatever its density. 20b keeps each reef's
  order and scales its class shares by reef density (last-season dives and
  tows; the GBR-wide forecast for reefs with neither, discounted by the 19 %
  of such reefs that historically needed culling). Forward test at equal
  sites: +7.3 points of > 0.04 dives kept (CI 2.7–11.0), 8 points fewer
  sites for 90 % retention. The gain is mostly between reefs. Separate
  product (`*_density_2027.*`); not yet adopted.
- **Deployment-realistic validation (script 33, 30 October).** Forward in
  time, with every reef's own history in training. The combined ranking gains
  little from that, but refreshing the evidence to each voyage's departure
  (earlier voyages in the season count) selects about 8 points fewer sites at
  90 % > 0.04 retention. It suggests re-ranking before each voyage; the
  product is unchanged. Those savings are against random selection; crews
  already dive open sites (last > 0.04) first, which alone scores about the
  same.
- **Option A, tow only the zone (script 34).** Crews dive open sites first,
  tow only the A–C zone, then let the tows guide the dives. This is the fair
  test against current practice; diving A sites without towing cannot be
  checked, because recorded dives were directed by tows.
  - Each tow counts once, on the visit that made it (fixed 2 November: the
    30-day window had double-counted the earlier visit's sweep).
  - On follow-up visits (after an earlier visit that season), the top-50 %
    zone saves 21 % of tows and misses 41 > 0.04 dives, against 82 for random
    placement and 68 for 9-on-3-off blocks. On first visits blocks are
    better. The zone avoids few ghost dives.
  - Reusing tows from the last 12 months saves about 40 % with little loss of
    > 0.04 dives, but adds ghost dives and stale triggers. It needs a field
    trial.
- **Proposed protocol (script 35, 2 November).** Written up for GBRMPA.
  - Priority classes are now A 20 %, B 50 %, C 66 %, D 80 %, E the rest
    (`config.R`, script 20); the old A–C (top 80 %) is now A–D.
  - Dive rule: without strong manta evidence, dive only sites ranked A–C, at
    crew discretion. It avoids 30 % of ghost and 19 % of low-density dives
    for 2.8 % of > 0.04 dives missed (mean CPUE 0.099); 4.4 ghost dives per
    miss against 2.0 at random.
  - The follow-up tow zone is A–C (top 66 %) since 6 November (was A–B).
    With season memory (sites near strong tows earlier that season) it saves
    13 % of follow-up tows and misses 13 > 0.04 dives (A–B: 20 %, 22).
  - One voyage cannot test the dive rule (about two flagged dives); about 25
    voyages of shadow records can.
- **Old SDM_Validation_Report ROC (AUC 0.865)** came from pooling across reefs,
  in-sample dives and point memorisation. On unseen 2026 dives the old RG
  map matches the current forward habitat layer. See the README.
- Details, numbers and caveats are in `DataProcessing/cots_sdm_v2/final/README.md`.
- **27 October ranking experiment.** Scripts 29–32 test three alternatives,
  including recent CPUE and polygon habitat percentiles. The original ranking
  remains the default. `outputs/site_ranking_candidates/` holds a separately
  labelled experimental recent-CPUE option; no candidate passed the strong
  replacement gate. Read the README comparison before changing deployment.

## Where things are

- `DataProcessing/cots_sdm_v2/final/` holds the delivered pipeline: layers,
  validation and report figures. **Read `final/README.md` first.** It records
  every result, caveat and the run order.
- `DataProcessing/cots_sdm_v2/final/report/technical_report.qmd` is the current
  technical-report source; `report/results.qmd` is an earlier Results draft.
  **Read `report/WRITING_BRIEF.md` and `report/CCIP_TEMPLATE_GUIDE.md` before
  report edits.** The supplied CCIP Bridging Project DOCX is immutable and
  governs section focus, length and formatting. Use
  `report/scripts/render_ccip_report.ps1` for the template-based Word/PDF export;
  a bare Quarto render produces only its content intermediate.
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
