# Frontiers submission pack

For *Frontiers in Marine Science*, Original Research, Research Topic 81777
"Operationalising Reef Interventions". Deadline 28 February 2027.

## Scope statement (for the submission form, 200 words)

Reef interventions are moving from trials to program-scale operations, yet few
studies test whether decision support improves those operations. This manuscript
evaluates fine-scale decision support for crown-of-thorns starfish control on the
Great Barrier Reef, one of the few reef interventions already running at program
scale. It fits the Coral Reef Research section and the Research Topic's themes of
decision-support frameworks and monitoring and evaluation. We test a site ranking
built from 30 m habitat models and recent surveillance on held-out reefs and
forward in time, against random selection at equal effort and against the triggers
crews already use. Replaying eight seasons of voyages shows where a ranking adds
value: applied to every site it mainly protects productive dives, whereas applied
only to weak-evidence dives, with density-aware classes, it avoids more empty dives
and misses fewer high-density ones. We report negative results plainly, including weak within-reef habitat
discrimination and the poor performance of model-ranked surveillance. We also
present a site-resolved population model as a calibrated prototype for strategic
scenarios, with its limits stated. The work offers a transferable approach to
evaluating decision support within the operational decision it is meant to change,
and sets out the field trial needed to measure realized benefits.

## Cover letter: required declaration

Samuel Matthews is a Topic Editor of this Research Topic and an author of this
manuscript, and requests assignment to an independent handling editor. [TO CONFIRM:
add Takuya Iwanaga if an author; Takuya Iwanaga is also a Topic Editor.]

## Must be done before submission

- [ ] Author list, order, full affiliations with city, and CRediT roles.
- [ ] Data availability: GBRMPA access statement and contact; repository and DOI
      for layers, ranking and tables.
- [ ] Frozen, citable releases of the R pipeline, COTSMod.jl and ADRIA. The
      reference list currently cites development revisions, one with uncommitted
      changes, which reviewers cannot reproduce.
- [ ] Ethics or permit statement, exact funding wording, and acknowledgments,
      including Traditional Owner wording in the authors' own words.
- [ ] Generative AI statement: tools, versions and purposes.
- [ ] GBRMPA agreement to publish results from restricted program records.

## Supplementary material

`manuscript_frontiers_supplementary.qmd` renders Supplementary Tables S1–S8 to DOCX
and PDF, numbered in the order the main text first cites them: habitat predictors
(S1), habitat-model training and settings (S2), the within-reef ranking coefficients
(S3), the density-aware v2 site models (S4), delivered products (S5), layer agreement
(S6), the season-by-season forward test (S7) and the class-set comparison (S8).
Elevation and turbidity sources were confirmed on 12 November. Table S5 still needs
the repository and DOI.

## Figures

`scripts/export_frontiers_figures.py` writes the nine figures, in order of
appearance, to `submission/figures/` as 180 mm wide, 300 dpi TIFF files.

- [ ] Remove the titles and subtitles burned into Figures 2–6, and the footnote
      in Figure 1, at source; the captions carry this text. Figure 7 has panel
      letters burned in, which Frontiers adds at typesetting.
- [ ] Check text is at least 8 pt at final size after the reduction to 180 mm.
- [ ] Add alt text for each figure on the submission form.
- [ ] Figures with Esri imagery have been left out, because the article is CC-BY.
      Figure 1 (`fig35_manuscript_gbr_and_workflow`, 12 November) shows the GBR-wide
      operational layer (RG + year with the clean + year fallback) with a 30 m
      Eyrie inset, and the 2027 workflow on four reefs over Sentinel-2 scenes
      (open Copernicus data, credited in the caption). It replaced the
      three-layer GBR overview (report fig01) and serves as the study-area map.

## Scientific risks to weigh

- The strategic half is retained (author decision, 11 November). It now reports
  checks on five reefs outside the calibration score, but has no validated
  control result. Martin and Linnet are the more independent check; Yonge, Carter
  and No Name informed the wave-exposure term, so expect reviewers to ask for that
  term to be tested on reefs not used in its design.
- All search results are retrospective. The January–February 2027 trial could
  strengthen the paper if early results arrive before the deadline.
- The recommended cut-offs and the visit-model weight were chosen on the test seasons;
  the Discussion says so.
