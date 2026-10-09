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
