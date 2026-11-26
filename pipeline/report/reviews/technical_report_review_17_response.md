# Review 17 response

Integrated on 26 November 2026 into `technical_report.qmd`.

- The fresh export passes the CCIP checks. The Summary is 359 rendered words; the only warning is the earlier one that the Discussion runs to 5 pages.
- It includes the metapopulation concept figure (Figure 3) and its references, which review 17 did not have.
- Next copy for feedback: `technical_report_review_18.docx`. The annotated review 17 is unchanged.

## Tracked changes

All accepted, with these repairs:

- **Cover.**
  - Isobel Ryan⁴ is added, with affiliation 4 as the Great Barrier Reef Marine Park Authority, Townsville. This is in the YAML author line, `ccip_report_metadata.json` and the suggested citation ("… Iwanaga T, Ryan I").
  - The subtitle now ends "AIMS, JCU, QUT, GBRMPA".
- **"replayed" → "back-tested"** throughout the report (c126), in Methods, Results, Discussion, Appendix D and the captions.
- **Methods, density-aware classes.** "(Figure 2)" became a live reference to the 2027 product figure. The deleted sentence about each site's own record is gone.
- **Methods, products.** The Eye-on-the-Reef sentence now reads: "To operationalise them most effectively, the pipeline would need access to the Eye on the Reef API, or would need to run on GBRMPA's Databricks system."
- **Methods, population model.** Your text was rewritten in grammatical form:
  - "delineated with a newly developed habitat-partitioning algorithm (Haller-Bull et al. 2026)";
  - "COTSMod.jl, developed for this project";
  - "Because computational constraints delayed the fine-scale connectivity networks within the project, larvae still move between reefs through reef-level connections …";
  - "Initial projections, developed to demonstrate capability, compared levels of control effort through 2050 to test their effect on suppressing primary outbreaks in a changing climate".
- **Engagement.**
  - Malu Ventures feedback on Far Northern sea country is added.
  - The habitat models were extended into the Torres Strait ("Straight" corrected).
  - Trials are now on 6–10 January 2027 everywhere: Summary, Methods, Discussion and Outcomes. The checker phrase in `scripts/check_ccip_report.py` was updated to match.
  - To keep the Summary within 360 words, its clause reads "vessel trials run 6–10 January 2027".
- **Results, habitat.** Now reads: "separate high-density COTS habitat from low-density habitat on the same reef consistently better than chance, but not well enough to direct culling or surveillance effort by themselves."
- **Results, site selection (para 217).** Your edits are kept, but the paragraph moved to Appendix C (see c122).
  - One fidelity repair: the ranking avoided 18.8% of low-density and 23.5% of ghost dives against 21.6% and 21.9% at random. That is about the same as random, not more.
  - The text now reads "protecting productive dives while avoiding about as many empty ones as random selection".
- **Results, density-aware classes.** The first-version comparison is deleted as you asked.
  - The tuning caveat is kept as one clause ("one weight … was tuned on these seasons, so the margin is somewhat optimistic").
  - "v2" is dropped from the main text, Table 2 and the Discussion; Appendix D still separates v1 and v2.
- **Results, recommended settings.** "(CPUE below 0.02)" and "only" are added, and the tow-zone sentence is deleted.
- **Results, surveillance.**
  - Your sentences now read: "In the proposed workflow, the first visit to a reef each season keeps a full sweep of the reef perimeter. Some of the tow effort saved on follow-up visits is also expected to go to exploring the suggested habitat outside the traditional tow paths."
  - "It is important to note" was dropped (voice rules).
  - The "Seasonal products and exploration" heading is removed, as you deleted it.
- **Results, population model.**
  - "for example" and "overstate the resources needed to control irruptions" are added.
  - The wave-exposure sentence keeps your "added to reduce recruitment on exposed reefs" and restores "designed from their records". AGENTS.md requires that distinction.
- **Figure 6 caption (calibration reefs).** "vermillion" became "orange". The same change was made in the other population-model captions and the new concept figure, so the colour is named the same way everywhere.

## Comments

| ID | Comment | Response |
|---|---|---|
| c62 | A Results placeholder to fill | The Methods sentence now points to a new Results section, "Larval connectivity patterns across the Initiation Box". |
| c66 | ReefPartitionUniversal citation | Added as given: Haller-Bull V, Grier B, Iwanaga T (2026), R package 1.0.2, github.com/open-aims/ReefPartitionUniversal. Initials only, as supplied. |
| c70 | Might be addressed by final hand-in | No change. If the fine-scale matrices enter the population model before submission, edit that sentence and Appendix E. |
| c77 | More Pratchett papers on feeding preferences | Added Pratchett (2007, *Pacific Science* 61: 113–120), on feeding preferences under controlled food availability, and Pratchett (2010, *Coral Reefs* 29: 717–725), on Lizard Island coral change during the 1995–99 outbreak. Both are in the text and Table 1; both were already in Mendeley. |
| c78 | What does "remembered" cover mean? | Now: "Adult survival depends on a running average of preferred-coral cover over recent years, rather than on the current year's cover alone". Table 1 matches. This follows COTSMod's exponential memory of preferred cover. |
| c79 | Cowan 2017 for fish; Wolfe or Desbiens for decorator crabs | "predators, from fish (Cowan et al. 2017) to rubble-dwelling decorator crabs (Desbiens et al. 2023), remove most small juveniles …". Cowan et al. (2017, *Diversity* 9: 7) is added to Table 1 too. Desbiens et al. 2023 (Wolfe is a co-author) covers the crab. |
| c87 | Keep the table borders | You applied Word's "Table Grid2". `build_ccip_report.py` now gives every main-body table that style, so Tables 1 and 2 keep full borders on every export. Appendix tables are unchanged. **Say if you want them gridded too.** |
| c122 | Present only the workflow back-test | The site-selection subsection is removed from Results; its numbers, with your edits, are now the first paragraph of Appendix C's "Recent evidence and cull-site selection". The coverage point (67% of sites on reefs with no recent data) moved into the habitat paragraph, and the open-sites point (61% of dives above 0.04) opens the workflow results. In Methods, the site-selection comparison is cut to one pointer to Appendix C. The workflow back-test now has three paragraphs: what was re-run, what the workflow keeps and changes, and how each setting was scored against the random baseline. |
| c123 | Remove the recent-CPUE sentence | Removed from Results; the comparison stays in Appendix C. |
| c126 | "replayed" → "back-tested" | Done throughout. |
| c131 | Not necessary detail | The Big Broadhurst and John Brewer example is removed; it remains in Appendix D. |
| c140 | Add a column of manta tows avoided | Table 2 has a new "Follow-up tows saved (%)" column. Dive rules alone save 0; the recommended row saves 16.7% [13.5–20.3]. The build script now fits six-column tables to the portrait page in the main body. |
| c163, c164 | Placeholder figure for pathways and clustering; 1–2 short placeholder paragraphs | Added the new Results section with two short paragraphs, `[TO WRITE]`/`[TO CONFIRM]` markers for the findings, and a placeholder figure (`assets/placeholder_connectivity_network.png`, panels a pathways and b clusters). The paragraphs describe the planned comparison against a distance-only control but state no results. |
| c177 | The evaluation metrics aren't fair; the fits looked good | Rewritten to lead with what worked: peaks within two years, and survey years ordered better than the no-model series. It then says the scores for individual irruption years were "mixed, partly because a peak shifted by a year or two counts as both a missed and a false irruption year". That explains why the true skill statistic is harsh without dropping the result. Appendix E's table still reports the scores. |

## Other changes made to keep the template checks passing

- **Results displays.** The new placeholder figure made five Results displays; the template allows 2–4. The larval-concentration hotspot figure moved to Appendix E, beside the concentration methods. The Results text still cites it.
- **Concept figure.** It is now cited from the population-model Methods, the Results calibration paragraph and Appendix E.

## Added after review 18 was issued (26 November)

- The workflow diagram (Figure 7, `fig23_workflow_decision.png`) was redrawn as a staged decision tree in the IPM operational-strategy format:
  - three bands: before the voyage, on arrival, and for each remaining site;
  - decision diamonds;
  - new steps in blue, plus a dashed box for exploring with the tows saved.
- It stays in the Discussion, and is now also cited in the Results paragraph on the recommended settings ("(Figure 7; Table 2)").
- The figure is built by `R/workflow_figure.R`; to rebuild it alone, run `07_build_report_figures.R --workflow_diagram_only=1`.
- Review 18 was already open, so it shows the older diagram.

## Still open

- Copyright owner on the cover.
- Back-cover contact position title.
- Campus cities for JCU and QUT.
- **The Summary still reports the site-selection result** ("Deferring the lowest-ranked 20% of sites retained 89.1% … against 78.5% at random"), which now sits in Appendix C. Say if you want it replaced with a workflow number.
- The connectivity network results and figure.
