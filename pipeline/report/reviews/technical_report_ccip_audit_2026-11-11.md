# Technical report audit against the CCIP template (11 November 2026)

An outside read of `technical_report.qmd` against the CCIP Bridging Project
template and `CCIP_TEMPLATE_GUIDE.md`. The appendix consolidation below has been
applied. Main-body findings are recommendations and have not been applied,
except where the consolidation required a changed cross-reference.

**Update, later on 11 November.** The author resolved findings 1, 2, 3, 7 and 8:
the Summary now states its baseline and intended benefits, the product is density-aware
v2 with a dive rule at A–D and a tow zone at A–C, the workshop outcome is recorded, and
Methods and Results were shortened. Findings 5 and 6 were addressed in the text (matrix
results stated as pending; a Surveillance section). Finding 4 and the administrative
gaps in finding 10 remain open. The table below is the audit as first written.

## Template compliance

| Section | Template guidance | Current | Verdict |
|---|---|---|---|
| Project Summary | One page; each strand, findings, benefits for people and environment | 349 words, one page | Passes on length. Headline lacks a baseline and benefits are not stated (findings 1, 8) |
| Context & Objectives | 2–3 pages; ends with three objectives as bullets | 3 pages | Passes |
| Methods | 2–3 pages of text plus 1–3 figures; domain map; permits; engagement subsection | 7 pages, 2 figures | Over length. Figure 1 has an empty panel; permits unconfirmed (findings 4, 7) |
| Results | 3–4 pages of text plus 2–4 displays | 9 pages, 4 displays | Display count passes; text over length (finding 7) |
| Discussion & Outputs | 3–4 pages; ends with outputs and their status | 5 pages | Slightly over. Outputs list present |
| Outcomes & Impact Pathway | 1–2 pages; Work Plan outcomes 2 and 4 | 2 pages | Passes, apart from the stale workshop date (finding 3) |
| Conclusions & Next Steps | ≤300 words | 249 words | Passes |
| Data Accessibility | Actual links, FAIR registration | Placeholders only | Not ready for release |
| Formatting | Template styles, fields, furniture | Export checks pass, template hash unchanged | Passes |
| Appendices | Detail that supports the main text | Was 15 appendices, 89 pages; now 5, about 41 pages | Addressed below |

`check_ccip_report.py` warns on Context, Discussion and Outcomes page counts but
does not test Methods or Results, so those two overruns have never been flagged.

## Findings, most important first

1. **The Summary reports the site rule without its baseline.** It says the A–D
   cutoff avoided 18.8% of low-density and 23.5% of ghost dives. Random selection
   of the same number of sites avoids 21.6% and 21.9%. On low-density dives the
   ranking is worse than random, and the Discussion says so, but the Summary reads
   as a gain. The real gain at A–D is retention (89.1% against 78.5%). Suggested
   replacement for the second half of the sentence: "…retained 89.1% of observed
   dives above 0.04 COTS per minute of bottom time, against 78.5% for random sites,
   but avoided no more low-density or empty dives than random selection." The
   workflow sentence that follows should also give its random comparison (29.8%
   ghost avoided and 2.8% missed, against 24.2% and 5.1%).
2. **A reader cannot tell which product to use.** The report presents within-reef
   classes (called the default), density-aware v1, v2 and a recent-CPUE candidate.
   Results say v2 "did better again", Figure 2 shows within-reef classes and the
   workflow map shows v2. The January–February trial needs one named product and
   one operating point (v2 at A–D or at A–C).
3. **Dated statements have gone stale.** The Summary, Methods and Outcomes describe
   a workshop "planned for 5 November 2026". Confirm whether it happened and what
   it agreed. `check_ccip_report.py` requires the phrase "5 November 2026", so
   keep the date when rewording.
4. **Figure 1 has an empty panel.** Panel C is "reserved for the hydrodynamic mesh
   figure, to be supplied". The hydrodynamic Methods are now written, so the panel
   should be supplied or removed.
5. **Connectivity Methods have no matching Results.** Methods describe building
   connectivity matrices, but Results stop at larval concentration and the
   population model still uses reef-level transport. Work Plan output F promises
   four reef clusters; the report covers the Lizard region. Either report the
   matrices or say plainly that they are not yet complete and why.
6. **Objective 2 is not answered directly.** It promises optimised, standardised
   surveillance paths with a comparison against perimeter surveillance. The report
   offers a thinning benchmark and a follow-up tow zone, and no new routes. Results
   are also organised by model strand, with surveillance inside the cull-site
   section. A short statement per objective, and a surveillance sub-heading, would
   let GBRF see what was delivered against each one.
7. **Methods and Results are about twice the guided length.** The main candidates
   to move or cut are the two paragraphs narrating density-aware v1 and v2 in
   Methods, the two number-dense density paragraphs in Results, and the
   hydrodynamic detail, which could sit in Appendix E.
8. **The Summary does not state benefits for people or the environment**, which
   the template asks for. Fewer empty dives and released effort for unsearched
   habitat could be framed as intended benefits still to be measured.
9. **Comparison with crew practice is under-signalled.** Savings against random
   are not savings against what crews do. Open sites hold 61% of dives above 0.04.
   The workflow replay is the fair test, and the main text could say so once,
   early in Results.
10. **Administrative gaps remain.** 22 `[TO CONFIRM]` items remain, including
    authors, organisations, permits, repository, DOI, licence and FAIR registration,
    plus the cover and contact fields in `ccip_report_metadata.json`.

## Appendix consolidation

| New appendix | Purpose | Drawn from |
|---|---|---|
| A: Outputs, Work Plan status and remaining actions | What was delivered and what is left | Old A, E, F, L (checklist) |
| B: Model specification and evaluation design | How to reproduce the models and tests | Old B, C (one sentence), G |
| C: Detailed search results | Full tables behind the main Results | Old D (part), H, M |
| D: On-water workflow, density-aware classes and field trial | The operational proposal | Old L (monitoring), N |
| E: Larval dispersal and population modelling | The strategic strand | Old G (metacommunity), I, K, L (limits and design), O |

Removed, with the underlying tables retained in `pipeline/outputs/tables/` and the
full text in git history:

- Appendix J (the superseded 15 October metacommunity comparison) and the earlier
  primary-layer effort analysis, example search plans and four-reef history map in
  Appendix D.
- Figures that repeated a table beside them, and the HMSC effect plots.
- Tables that duplicated another or kept archived definitions, such as the
  percentile-cutoff classes.
- Three of the four ranking-alternative displays, and the two SSP projection time
  series (the effects figure covers both pathways).

The former Appendix L paragraph drafted in the past tense ("We compared paired
simulations…") was removed. That assessment is incomplete, and Appendix E now
describes it as planned. Restore it if the assessment finishes before release.

One figure error was corrected: the retention curve (now in Appendix C) still
labelled its cutoff "A-C only (defer D)" from before the five-class relabel.
