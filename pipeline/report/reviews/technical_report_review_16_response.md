# Review 16 response: Summary, Context and Methods, plus a voice pass on Methods and Results

Integrated 25 November 2026 into `technical_report.qmd`. The fresh export passes
the CCIP checks: the Summary is 360 rendered words (the limit), and the only
warning is the earlier one that the Discussion runs to 5 pages. Next copy for
feedback: `technical_report_review_17.docx`. The annotated review 16 is unchanged.

## Tracked changes

All accepted, with these proofreading repairs:

- **Summary.**
  - "the generated 30 m GBR habitat maps" became "the 30 m GBR habitat maps",
    because the Summary needed every word.
  - "engagement with traditional owner-led contractors (Malu Ventures)
    supported a future Far Northern specific model" became "engagement with
    Traditional Owner-led contractors (Malu Ventures) supported a future Far
    Northern-specific model".
  - Your rewrite of the population-model sentences is kept, with one clause
    restored: the three almost COTS-free reefs stayed free of false irruptions
    only "with a wave-exposure term designed from their records". AGENTS.md
    requires that distinction, because agreement on those reefs checks the
    term's design and is not validation.
  - "support important tactical and strategic decision-making" became "provide
    a basis for important tactical and strategic decisions". The population
    model is not yet validated for strategic use.
  - **Your call:** whether "provide a basis for" is acceptable.
  - "keeps diving … but ranks the remaining dives" became "still dives … but
    ranks the rest", to fit 360 words.
- **Effort shares, Summary and Context.**
  - "20–40%" with "David Williamson pers comms" is now cited as
    "(Great Barrier Reef Marine Park Authority 2026, D. Williamson pers. comm.)"
    in both places, so the Summary and Context agree.
- **Context.**
  - The IPM sentence now reads "links decisions at individual sites, reefs and
    voyages with regional strategic priorities set by the annual prioritisation
    process", citing Matthews et al. (2026, *Biological Conservation* 321:
    111910; comment c41). The entry was checked against Crossref.
  - "Manta tow can therefore" became "Manta tows can therefore".
  - "to more efficiently direct" became "to direct … more efficiently".
  - "at the same reduction of effort" became "at the same reduction in effort".
- **Methods, first section.**
  - "srategies", "rnaking" and "dens ities" were corrected.
  - "Model probabilities were converted to site ranking which were assessed …
    with dives … formed" became two grammatical sentences.
  - "(>10)" for common scars is kept as "(more than 10)".
- **Methods headings.**
  - "COTS-Coral meta-population model" became "COTS–coral metapopulation
    model", matching the title's spelling.
  - "Engagement with traditional owners and Management Stakeholders" was
    **not** applied. The CCIP template requires the heading "Engagement with
    Rights-holders and Management Stakeholders", and the export check fails
    without it. Malu Ventures is now described as "Traditional Owner-led
    contractors" in that section instead.
- **Trials.**
  - Deleting the Malu trial in the Far Northern GBR and Torres Strait was also
    applied to Outcomes & Impact Pathway, so trials are planned only aboard
    BPM in Cairns throughout.
  - **Check** that this is right.
- **Connectivity.**
  - The modellers' paragraph you moved into Methods is now in past tense, with
    the grammar fixed and "Appendix X" changed to Appendix E.
  - The duplicate copy in Appendix E was replaced with a pointer to Methods.
- **Accepted as written.**
  - Deletion of the permits placeholder (c87).
  - Deletion of the workshop-participants placeholder.
  - Deletion of the progress-report sentence.
  - The two sentences cut from the Figure 1 caption.
  - "Testing a proposed on-water workflow" and "Cull and surveillance search
    products".

## Comments

| ID | Comment | Response |
|---|---|---|
| c41 | Biol. Conserv. DOI | Cited as Matthews et al. (2026), doi 10.1016/j.biocon.2026.111910, for the annual prioritisation process. |
| c57 | Methods barely mention how the habitat models were built and tested | New paragraph at the end of the section. It covers the XGBoost models (full-extent and ReefGuide versions combined into one operational layer), the manta and HMSC models, why reefs are scored by models fitted without them (grouped five-fold validation), and how within reef-year AUC and top-20% lift are read. Details remain in Appendix B. |
| c86 | "open-AIMS" is not a full citation; name the repository | The text now names ReefGuide.jl, COTSMod.jl, ADRIA.jl and CoralBlox.jl, and the ReefGuide entry has a year, so it no longer cites as "open-AIMS n.d.". |
| c87 | Permits sentence unnecessary | Deleted. |
| c118 | Describe the density-aware layer behind the rankings | New paragraph after the class definitions: within-reef classes give every reef the same share of A sites, so the density-aware classes scale those shares with recent reef density, giving low-density reefs fewer A sites. |
| c119 | Hard to parse | Rewritten as "Two safeguards kept each test out of sample", followed by the two safeguards. |
| c139 | What are "density quotas"? No mention of the Stewart models for reefs without evidence | Now explained: reefs with no dives or tows in the past year take their density from the GBR-wide reef forecast (Stewart et al. 2026a), turned into the chance that a survey needs culling and the share of problematic dives if it does. That density sets how many of their sites fall in each class, and sites are ordered by habitat rank. Checked against `scripts/37_site_ranking_density_v2_2027.R` and `R/density_classes.R`. |
| c172 | Two paragraphs on the key mechanisms, with literature, bringing detail back from the appendix | Two paragraphs. The first covers model structure: the three packages, reef-level transport and calibration reefs. The second walks through grazing, starvation, fecundity, the Allee effect, the predator pit, larval supply and wave-exposed settlement, each with its source. |
| c173 | A table of the mechanisms and their literature support | New Table 1 in Methods: mechanism, how the model represents it, and support. Settlement on wave-exposed reefs is marked "designed from LTMP records; not independently tested", because no field study supports it. |

New references, each checked against Crossref and added to Mendeley ("Claude
added" for the five not already there):

- De'ath & Moran 1998
- Keesing & Lucas 1992
- Rogers & Plagányi 2022
- Caballes et al. 2016
- Babcock & Mundy 1992
- Rogers et al. 2017
- Keesing et al. 2018
- Desbiens et al. 2023
- McCallum 1987
- Fabricius et al. 2010
- Matthews et al. 2026 (*Biological Conservation*)

## Voice pass, Methods and Results

These sections were rewritten with the marine-science-writing voice rules
(`voice.md`):

- At most two intervals per sentence.
- Each Results paragraph opens with its claim and ends with what it means.
- Named reefs carry the density example: Big Broadhurst and John Brewer, and
  the five reefs outside the population model's calibration.
- Plain comparisons are used, such as "about 4 times the random rate".
- A nested-bracket sentence was removed.
- Terms such as "site rule" and "season memory" are now defined where they are
  first used, because your Methods edits removed their earlier definitions.

Every number is still drawn inline from the output tables.

Length: Methods now runs from page 5 to page 10, including two figures and the
new table. That is over the template's guide of 2–3 pages of text plus figures,
because of the requested additions. Results stays about the same length.

## Added after review 17 was issued (26 November)

The new Figure 3 is a conceptual diagram of the metapopulation model, placed
after Table 1 in Methods:

- **a.** The 2,895-site domain.
- **b.** The COTS life cycle and its coral feedbacks.
- **c.** The annual step order, and which of ADRIA.jl, CoralBlox.jl and COTSMod.jl
  runs each step.

The figure is built by `R/metapop_concept_figure.R`; to rebuild it alone, run
`07_build_report_figures.R --metapop_concept_only=1`. Review 17 was already being
annotated, so it does not contain this figure. Figures after it are renumbered by
one.

## Still open

- Copyright owner on the cover.
- Back-cover contact position title.
- Campus cities for JCU and QUT.
- Whether the Malu trial removal also applies to Outcomes (applied here).
- The other report `[TO CONFIRM]` items.
