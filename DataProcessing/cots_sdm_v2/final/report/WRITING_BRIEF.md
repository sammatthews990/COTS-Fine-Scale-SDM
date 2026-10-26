# Writing brief: COTS fine-scale SDM

Read by the `marine-science-writing` skill before drafting. Update it when a
decision changes. Items marked `[TO CONFIRM]` are for the authors to supply.

## Deliverables

| # | Document | Skill type | File (in this folder) | Venue / recipient | Due |
|---|---|---|---|---|---|
| 1 | CCIP Final Project Report | `ccip-report` | `ccip_final_report.qmd` | GBRF, COTS Control Innovation Program | [TO CONFIRM] |
| 2 | Technical report (CCIP Output A, public-facing) | `technical-report` | `technical_report.qmd` | GBRF / GBRMPA COTS Control Program | [TO CONFIRM] |
| 3 | Journal article, Original Research | `article-full` | `manuscript_frontiers.qmd` | Frontiers in Marine Science, Research Topic 81777 "Operationalising Reef Interventions" | 28 Feb 2027 |
| 4 | (optional) Short-form or Brief Research Report | `article-short` / `article-brief` | — | [TO CONFIRM] | — |

Scaffold each from this folder, e.g.:

```bash
Rscript ~/.claude/skills/marine-science-writing/scripts/new_document.R \
  --type=article-full --out=manuscript_frontiers.qmd --tables=../outputs/tables
```

`results.qmd` (existing) is a draft Results section for the technical report;
reuse its text and inline-number pattern.

## Sources of truth

- `../README.md`: goal, methods, every result with the table behind it,
  limitations to state, figure list, provenance.
- `../../../../AGENTS.md`: current state of findings (dated) and the goal.
- `../outputs/tables/*.csv`: every number in every document comes from here
  via `tab()` / `pick()`; never type a result.
- `../outputs/figures/fig01–fig17*.png`: built by `../scripts/07` (and 13);
  7.5 in (190 mm) wide at 300 dpi. **Frontiers needs ≤ 180 mm and TIFF/JPEG**:
  re-export at 180 mm (or 85 mm for single-column) for submission.
- `../outputs/site_ranking/`: the operational 2027 cull-site priority product
  (scripts 19–21).

## Goal and framing (keep consistent across all documents)

The project predicts where on a reef COTS are most likely to be found, so the
COTS Control Program can (1) manta tow less without missing aggregations and
(2) search fewer cull sites without missing problematic densities
(CPUE ≥ 0.02 COTS per minute of bottom time). Every result is judged on
held-out data against a no-information baseline at the same effort.

Key messages (detail and numbers in README/AGENTS.md):

1. Habitat layers rank locations within a reef only weakly (held-out
   within reef-year AUC ~0.55–0.61). They set search order; they do not
   locate aggregations.
2. For manta surveillance, evenly spaced (systematic) thinning keeps more
   aggregations than dropping low-ranked tows by any model (script 14).
3. For cull-site selection, the gain comes from last season's cull and manta
   evidence plus a habitat rank (scripts 16, 19): fewer designated sites
   searched and more low-density/ghost dives avoided for the same problematic
   dives; it has worked every season since 2022.
4. Cull history forecasts well but partly encodes where crews choose to dive
   (scripts 09–10); report it as "known aggregation sites", not habitat.
5. Operational product: the 2027 priority ranking of all designated cull sites
   plus exploration areas outside the current search footprint (scripts 20–21).

Decisions already made — do not re-open in the text:
- Primary within-reef layer: RG + year with clean + year fallback; HMSC is
  explanatory only; no blended ensemble.
- History is a separate "known aggregation sites" layer.
- Layers set the *order* of search, not exclusion of sites.
- Always report held-out numbers and the difference from random/no-information.

Do not overclaim:
- Reef-level outbreak detection by the layers is barely above random.
- HMSC forward transfer is near chance (0.515).
- Connectivity is not in these models. A reef-level hybrid forecast (Stewart
  et al. preprint) added little beyond last season's cull record when
  checked; a fine-scale hydrodynamic connectivity pilot (Port Douglas–Lizard)
  is future work.
- Expert zone review is still pending.

## Prior work to cite (all keys in references.bib, DOIs verified)

- Reef-level COTS SDMs: `matthews2020connectivity` (larval connectivity and
  water quality explain outbreak distribution); `stewart2026hybrid` (hybrid
  machine-learning reef-year outbreak forecast, EcoEvoRxiv preprint; also the
  metacommunity model's external larval supply).
- Fine-scale coral SDMs: `guillaume2025coral` (multiscale terrain; the coral
  habitat predictors are layers of this kind).
- Metacommunity lineage: `matthews2020cotsmod` (COTSMod), implemented in
  COTSMod.jl (`cotsmod2026`) within ADRIA (`adriacalibration2026`) with
  CoralBlox coral dynamics (`dealmeida2026coralblox`, bioRxiv preprint).
- Cues from other COTS models: MICE (`morello2014model`, `plaganyi2014mice`;
  prey switching, diet-linked starvation), CoCoNet (`condie2021interventions`;
  juvenile predation, Allee form), ReefMod-GBR (`bozec2022cumulative`;
  post-outbreak collapse, connectivity matrices).

## Research Topic fit (Frontiers 81777)

Frame the article as operational evaluation of a decision-support tool for
reef pest control at program scale: performance against current practice and
a null baseline, what that means for deployment, and how it will be
monitored. Topic keywords that apply: decision support; monitoring and
evaluation; reef interventions; scaling strategies.

**Conflict of interest:** Samuel Matthews is a topic editor of this Research
Topic. Declare it in the cover letter and the Conflict of interest statement
(standard Frontiers wording in the skill's `frontiers_marine_science.md`);
Frontiers will assign an independent handling editor.

## Authors and admin

- Authors, order, affiliations, CRediT roles: [TO CONFIRM]
- CCIP project code and title, contract dates, Work Plan deliverables, metrics
  targets, budget, WHS, Traditional Owner engagement: [TO CONFIRM — supply the
  CCIP Work Plan]
- Funding statement: COTS Control Innovation Program (Reef Trust, AIMS, GBRF)
  — [TO CONFIRM project code and wording]
- Data availability: COTS Control Program cull and manta data are held by
  GBRMPA (restricted); AIMS LTMP; ReefGuide layers — [TO CONFIRM access
  statements and any repository DOI for layers and code]
- Language: Australian English for CCIP and technical report; Frontiers
  article [TO CONFIRM: American, or state Australian English].

## Writing preference: explain model differences in words

Use shorthand model names sparingly, especially in the metacommunity narrative.
Describe the relevant difference when it matters: for example, "the version
started with more preferred coral at Lizard Island" or "the version with lower
background coral mortality". Readers should not need to remember W, C3a, C6,
score-version numbers, domain-version numbers or run identifiers to follow the
argument. Keep those identifiers in Methods, tables and source notes where they
are needed for traceability; introduce them once with their meaning. In the
Discussion, explain processes, evidence and management consequences directly.
This is the author's preference, recorded on 15 October 2026.

## Discussion outline confirmed by the author

1. Combine recent observations with habitat for daily search planning. Sources:
   `cull_site_evidence_matched.csv`, `site_ranking_forward.csv` and the current
   cull-site Results. Describe habitat's role in ordering places lacking evidence
   without reducing the project to a comparison of model labels.
2. Preserve spatial coverage in manta surveillance. Sources:
   `manta_aggregation_cv_summary.csv`, `manta_aggregation_cv_max_cut.csv` and the
   manta Results. Separate retrospective thinning from new route design.
3. Explain what the metacommunity prototype currently contributes: diagnosing
   population, coral and observation mismatches, with historical culling as a
   consistency check. Sources: `metacommunity_candidates.csv`,
   `metacommunity_culling_removal_summary.csv` and
   `metacommunity_2026_diagnostics.csv`. Peak matching and conditional qualitative
   checks do not establish benefits of adaptive control policies.
4. Include the author’s added emphasis: assess the culling program from 2012 to
   the present using paired historical simulations, and project suppression of
   the fifth primary outbreak wave and its local/downstream coral consequences.
   These are proposed assessments, not completed policy results.
5. Define the next evidence needed: field trials of the ranking and exploration
   areas, independent observations and held-out validation of the metacommunity
   model. Retain the distinction between demonstrated search gains and proposed
   strategic uses. Sources: current Methods, Results and pipeline limitations.

## Glossary (use these terms consistently)

- **Problematic dive:** cull dive with CPUE ≥ 0.02 COTS per minute of bottom time.
- **Low-density dive:** COTS found, CPUE < 0.02. **Ghost dive:** no COTS found.
- **Strong-evidence tow:** manta tow with direct COTS, ≥ 2 counted scars or
  legacy "common" scars. **Aggregation:** strong-evidence tows within 500 m in
  one survey.
- **Within reef-year AUC:** discrimination among locations on the same reef in
  the same year (0.5 = none).
- **Held-out / leave-reef-out / forward:** scored by a model that never saw
  that reef (or that season).
- **Designated cull site:** a site in the EoTR cull-site KMZ.

## Methods pass (15 October 2026)

Methods in the technical report and Frontiers manuscript now describe the SDM
pipeline and the metacommunity model, including observed historical culling.
Sources were the final scripts and helpers, and the local COTSMod.jl / ADRIA
development worktrees. Separate bibliography entries record their commits
and the uncommitted ADRIA changes; these are not final release citations.

Validation wording distinguishes grouped reef holdout from temporal holdout:
script 19 trains the evidence regression on preceding seasons but uses spatially
held-out habitat ranks from script 17, whose fits are not restricted in time.
The matched check of deployed layers is described separately. Forward ranking
intervals use 200 resamples; the main spatial comparisons use 500.

The Results reconciliation below supersedes the earlier 12 October comparison. Do not describe prescribed historical catches
as an implemented adaptive control policy. Final model release, package versions,
connectivity attribution and source-data access details remain to be confirmed.

On Windows, set `$env:LC_ALL = 'English_Australia.utf8'` before Quarto renders.
The inherited `C.UTF-8` locale is unavailable to this R installation and causes
some inline Unicode punctuation to print as character codes.

## Metacommunity Results reconciliation (15 October 2026)

The report and manuscript now use script 22 and `outputs/tables/metacommunity_*.csv`.
The former `report/external/metacommunity/` comparison is preserved as superseded.
All 15 source files are pinned in config.R (ADRIA_DIR defaults to the sibling
ADRIA.jl repository). Figures 18–19 are built by script 07; use
`--metacommunity_only=1` to rebuild them without the SDM figures.

Evidence is separated into the replicated W/C3a/C6 comparison (14 October run 48),
historical removals (run 42), and single-seed v0.4 diagnostics through 2026
(15 October runs 08–13). All three candidates match 7/7 historical peaks in each
of five seed pairs under score v2. This is calibration, and G2′ passes depend on
the observation curve. At A½ 30: W 5/5, C3a 0/5, C6 3/5; at A½ 10: W 0/5,
C3a 0/5, C6 3/5. Do not say that no candidate passes any qualitative gate, or that
passing such a gate completes validation.

Historical culling reveals missing recent populations (Eyrie 2021–24: 88% of
requested removals unmet at sigma 210 m). Latest diagnostics have culling off.
Score v3 includes the diver-anchored recent-density term and is not comparable
numerically with v2. Corrected outbreak inputs share cull/SALAD evidence with the
assessment, so these comparisons are not independent validation. There are no
validated results for adaptive starting/stopping rules or primary suppression.

Only the metacommunity Results, dependent Methods text and contradictory summary
claims were reconciled. The Discussion pass below addresses interpretation and intended strategic uses.
Final source attribution, licensing/access statements and software-release details
remain outstanding.

## Discussion pass (15 October 2026)

The author confirmed the outline and requested explicit treatment of the model's
intended historical assessment of culling from 2012 to the present and projections
of fifth-primary-wave suppression and future coral consequences. These uses are
described as planned paired comparisons. The Discussion explains controls,
outcomes, matched effort for future active strategies, external-boundary limits
on attribution, and the need for an effort/detection component for policy rules.
No new simulation results or efficacy estimates were introduced.

Both documents explain differences between metacommunity configurations in words;
candidate identifiers, score-version labels and run codes remain in technical
Methods/Results material. The report's management implications, limitations and
next steps were revised together; the manuscript Discussion is connected prose.
Existing model-fit results remain the frozen 15 October reporting snapshot.

Open author decisions: historical assessment endpoint and complete record coverage,
projection horizon, operational budgets, intervention rules, climate scenarios,
trial design and monitoring ownership. These do not imply an approved simulation
experiment. The next writing pass should align the Abstract/Executive summary and
Introduction/Background with the completed Discussion, within their word limits.

## Word review 01 integrated (16 October 2026)

Scope: technical report section 4 only. Review input is preserved in `reviews/`;
its original text matched the saved clean baseline. The manuscript and the frozen
model results were not revised in this pass. The author confirmed January 2027
and refined the proposed field-trial criteria: 10–15% fewer sites searched by
manta tow and 10–15% fewer low-density/empty dives, with negligible missed dives
above 0.04 CPUE. Deferred sites retain annual surveillance or a check before the
reef closes. The existing validation uses >=0.02 and must not be relabelled as a
test at >0.04. The numerical meaning of "negligible" remains to be agreed.

Recurring preferences learned from this review:
- Explain comparison groups and consequences concretely. Replace phrases such
  as "operational comparator" with what crews would normally do.
- Describe larval boundary assumptions in terms of outside reefs and the actual
  counterfactual, with a named source where relevant.
- Focus strategic recommendations on decisions from now, including both
  reallocation of existing resources and additional vessel capacity.
- Discuss primary-outbreak suppression separately from local protection during
  secondary outbreaks, and distinguish local from downstream coral benefits.
- Use "strategies robust to deep uncertainty" for the stated decision objective;
  do not claim robustness has already been demonstrated.
- Keep transient development diagnoses out of management implications when
  requested; retain the underlying scientific requirements and Results evidence.
- The author may draft the intended final assessment in past tense ahead of its
  completion. Preserve that editorial intent with an explicit completion check;
  do not add unverified numerical benefits or call an unfinished result validated.

Local scenario facts are author-supplied, not general style rules: control periods
2014–15 and 2020–2026, increased resources in 2026–2028, one or two additional
vessels, and southern-domain checks at Irene Reef and Startle Reefs. The full
historical assessment and record coverage through 2026 still require confirmation.
Both Codex and Claude skills now route Word feedback integration through the same
section-scoped extractor and instruct future agents to learn correction patterns
without generalising individual edits into universal rules.

## Word review 02: Background and objective alignment (16 October 2026)

Integrated sections 1.1–1.5 and dependent management discussion in the technical
report. The former 1.5 is now woven into the Background or relocated to Methods;
the metacommunity cross-reference follows the Methods subsection. The manuscript
is unchanged. Methods/Results analysis work is queued in
`reviews/technical_report_review_02_next_iteration.md`.

Author preferences and decisions:
- Use **irruption** throughout report narration. Preserve published titles,
  source-variable names and lookup strings for traceability.
- Frame efficiency alongside efficacy: effort released from repeated searches
  should allow searches of overlooked reef habitat, not reinforce the footprint.
- Distinguish <0.02 low-density dives, the 0.02–0.04 intermediate group, and >0.04
  dives whose omission carries the management penalty. Retain >=0.02 only as the
  existing analysis definition until the separate evaluation is run. The author
  accepts some missed intermediate dives; do not claim those omissions are
  already demonstrated to be safe. Ecological thresholds depend on coral cover.
- Keep model purpose in Background and implementation/calibration detail in
  Methods. Refine management questions within the contractual objectives without
  silently presenting new wording as an amended Work Plan.
- Explain baselines with the actual search choices and consequences.
- Strategic objectives include stopping where local coral cannot be saved,
  continued control of larval-source reefs despite limited local benefit, and
  resources needed for fifth-primary-wave suppression and downstream protection.

Annual activity of approximately 1,300 sites across 200 reefs is author-supplied
from the dashboard; the precise period/site figure needs a final source snapshot.
No new performance numbers or model runs were introduced.

## Word review 03 and threshold reassessment (16 October 2026)

Integrated all review 03 comments and tracked prose edits. Methods now expand
ODMAP before use, separate the two CPUE thresholds, and describe the author-
confirmed 2002–2026 source period. Site-matching detail remains in code; the
report identifies 4,299 matched sites on 339 reefs (361 reefs in the wider cleaned
record set). Objective 3 now describes site-level strategies for individual reefs,
following the author's wording, without treating the revised list as a quotation
of the contractual Work Plan. Word's untracked minus-sign normalisation in units
was treated as formatting, not a scientific revision.

Script 23 re-evaluates saved rankings from 15–17 and 19, without refitting or
altering the 2027 product. New input pins preserve the raw records, site register
and ranking snapshots. Actual selected candidate site-years are counted, including
never-dived sites; random expectations use the same site count per reef. Targets
90%, 95% and 99% are sensitivity levels, not agreed definitions of negligible.
Paired uncertainty uses 500 reef-year resamples and every interpolation is repeated
inside each draw. Original >=0.02 tables/figures are retained with explicit labels.

New findings: in the forward 2019–2026 comparison, 95% >0.04 retention requires
90.2% of candidate site-years with the combined ranking; 99% requires 98.2%.
These results do not establish the proposed 10–15% adoption gains at a negligible
miss rate. Numerical claims in the report read the new tables inline. Figure 20
and the high-density class summary show this trade-off. The frozen 15 October
metacommunity results remain unchanged; scenario design is now described explicitly
as planned work. The manuscript has not been revised in this pass.

Additional preferences: explain reporting frameworks and spell out acronyms on
first use; keep matching implementation details in code when they distract from
the study design; distinguish all-record coverage from matched evaluation samples.

Open work requiring new information remains in the updated review 02 revision
queue: missed-dive tolerance, field outcomes and annual effort accounting,
dashboard site-count period, connectivity attribution and delivery, final model
release and agreed operational scenario settings. Do not turn those items into
completed results.

## Review 04: write for marine ecology readers (16 October 2026)

This is a major editorial preference, applying to Methods and Results as well as
Discussion: introduce each model's purpose and inputs before naming it. Use plain
descriptions in the main text. Main Methods/Results should present the minimum
evidence needed to answer the objectives; put equations, fitting settings, model
codes, alternative comparisons and diagnostic histories in appendices. A table
of model purposes and a table of evaluation decisions help readers follow the
study. The annual models are described as full-extent (13 predictors including
year) and limited-extent (18 predictors including year); original software codes
are mapped in Appendix G. This preference supplements the earlier guidance on
avoiding unexplained metacommunity shorthand.

The missed title edits from review 03 are restored: “Where and when to search for
COTS on a reef”, with a subtitle covering fine-scale habitat suitability and
metacommunity modelling for tactical and strategic decisions on the GBR.
For future Word reviews, inspect title/subtitle and other front matter separately
as well as extracting the named body sections; a section-only extraction does
not include title-page edits. Align the cover metadata and suggested citation.

Methods comments and tracked edits supplied source-program dates, the operational
0.04 CPUE criterion for reefs below 40% coral cover, coverage limits of ReefGuide
layers in deeper reef areas, and source citations for SALAD, RHIS, ReefScan and
habitat maps. Code confirms ReefScan was used as a separate presence-only check,
not a fitted response in the delivered joint model; retain that factual distinction
despite the tracked grouping of these sources. Programme date ranges are not the
sample period of every individual model fit. Predictor product versions and
resampling details still need confirmation.

Main Methods and Results were condensed by approximately 58% and 73%, respectively,
with detailed material retained in Appendices G–I. The original analyses, thresholds,
numbers and operational products are unchanged. A simpler surveillance figure
uses the same existing aggregation table; the original figure remains in Appendix H.
The main report retains three figures; detailed maps and model comparisons remain
available in its appendices. The manuscript was not rewritten in this pass.

### Calibration figure update — 21 October 2026

Use the dated 21 October W-wave-lf04 figure set, rather than the older figures
at the top of the ADRIA calibration-report directory. Section 3.4 pairs COTS and
coral by reef; its heading is “Historical COTS and coral dynamics provide a basis
for testing control strategies”. Main text should explain progress and remaining
population/coral mismatches in plain language. Technical configuration codes and
sensitivity details belong in Appendix I. Preserve review 03's overall title.

Appendix I contains current recent-density, additional-reef and projection
figures. Earlier October 15 comparisons are archived in Appendix J. Effort-based
projections are now implemented, so avoid saying they remain entirely prospective;
operational validation, source-reef assumptions and domain-wide coral consequences
remain unfinished. Distinguish full calibration ranges from projection quantiles
and confidence intervals. Do not interpret paired historical overlays as a
completed independent program-impact assessment.

### Review 06 — hydrodynamics and connectivity (23 October 2026)

Comments 26 and 45 request a separate hydrodynamic/connectivity account before
the metacommunity Methods and Results. Leave two explicit author-completion
paragraphs for Sev and Owen: hydrodynamic development and particle tracking;
then matrix construction and network analysis. Do not infer their technical
settings from the concentration figures. Treat the concentration-map analysis
as preliminary descriptive evidence, distinct from source-to-destination
connectivity, local retention and settlement. The local network analysis plan
awaits the matrix; proposed network results must remain prospective. The figure
location was confirmed as `DataProcessing/connectivity_lizard/outputs/figures`.
The main report uses a portrait version of lc2; Appendix K uses lc4 and lc6.
The population-model Results move to Section 3.5 after the new dispersal section.

### Retention options and fallback coverage — 26 October 2026

Table 3 and Section 3.1.1 show 80%, 90%, 95% and 99% retention of >0.04 dives,
with matched habitat and random baselines and paired reef-year intervals.
Script 23 was extended to calculate 80% retention with the same rankings and
500 bootstrap draws; existing targets remain numerically unchanged. Do not
confound a stricter retention target with a change in CPUE threshold, or frame
small near-complete-retention savings as an absence of ranking benefit.

Discuss annual full-coverage surveys including all deferred sites, pre-closure
checks and earlier checks triggered by new evidence as options for field testing.
These are proposed ways to find false negatives, not evaluated recovery guarantees.
Account for detection limits, delayed control and all follow-up effort; initial
site reductions are not established annual vessel-time savings. The choice of
retention, acceptable omissions and maximum delay requires agreement with the
program. Preserve existing review copies and provide review 08 with clean baselines.

Before the next version, the author requested low-density and ghost cull-dive
avoidance in Table 3 and the text. Show these existing matched-results metrics
for all three rankings and all four retention targets, with uncertainty. Define
the denominators: designated sites searched, low-density dives avoided and ghost
dives avoided are different quantities. Keep the positive-below-0.02 definition
of low density and the separate 0.02–0.04 intermediate class explicit.

### Review 08 — direct Results and brief management implications

Use the author's revisions as a tone guide: lead with the practical finding,
keep model names descriptive, and avoid repeating the same caveat or baseline
values in adjacent paragraphs. Retain uncertainty and scope where they are
needed, with fuller detail in the appendices. Do not interpret a preference for
confident prose as permission to assert unmeasured culling efficacy or forecasts.

For the five management subsections specifically, the author requests short
bullet paragraphs with appendix pointers; 4.1 can be longer. This is an explicit
project preference overriding the generic prose-first guidance for that section.
Appendix L retains detailed management interpretation and assessment design.

The four-reef habitat map belongs in the main habitat Results, with no embedded
overall title, subtitle or footnote. Keep necessary panel labels, legends and
imagery attribution (the latter in the appendix). Standalone thinning results
and the simplified thinning figure are removed from main Results; the benchmark
remains in Appendix H. Explain practical searches through crew judgement and
recent evidence without claiming the cull-site ranking beat the tow benchmark.

Bring out the prospective early-warning use of larval concentration plus network
analysis in the initiation box. Describe the fifth-wave simulations and 2050
scenarios as demonstrations of capability; recent population overshoot and its
possible culling bias still require testing. Review 08 remains the annotated
input; review 09 is the next clean review copy.
