# Review 18 response

Integrated on 26 November 2026 into `technical_report.qmd`.

- The fresh export passes the CCIP checks. The Summary is 359 rendered words.
- The Discussion warning remains, at 5 pages. The Discussion text is now 1,359 words, down from 1,900, but the workflow diagram takes most of a page.
- The new Figure 7 workflow diagram is in this export; review 18 showed the old one.
- Next copy for feedback: `technical_report_review_19.docx`. The annotated review 18 is unchanged.

## Tracked changes

All accepted, with these repairs:

- **Surveillance.**
  - Now reads "following feedback from GBRMPA".
  - Now reads "between full surveillance sweeps".
  - Your follow-up definition became its own sentence: "Follow-up visits go to reefs under active culling, which the IPM requires to be re-towed every six weeks."
- **Connectivity placeholder.** Now reads "the strongest links between reefs and sites".
- **Discussion opening.**
  - "can vastly improve" became "substantially improved". Keeping 90% of dives above 0.04 needed 80% of sites with the ranking against 91% at random, which is a clear gain but not a vast one.
  - "predicted habitat" is kept.
- **Discussion, exploration zones.** "repeated detections" and "surveillance and culling footprint" are kept, and the zones now cite Figure 2.
- **Discussion, workflow.** Your decision-support-system sentences were rewritten in grammatical form and cited to Matthews et al. 2026a (Conservation Intelligence, c105): "Recent developments to that system provide a web map, updated daily, that colour-codes the importance of sites, a ready platform for the site ranking and exploration zones developed here."

## Comments

| ID | Comment | Response |
|---|---|---|
| c27 | Figure 2 (2027 product) should be in Results | Moved to Results, after the paragraph on the delivered ranking and exploration zones. The Methods and Discussion still cite it. At your request the habitat-map figure stays in Results, so Results now has five displays; `check_ccip_report.py` now allows five. |
| c87 | Cite Pratchett on outbreak dynamics following available coral | Pratchett (2005, *Coral Reefs* 24: 453–462), on outbreak dynamics at Lizard Island in 1995–99, now supports "aggregations moving on as they deplete their preferred coral". It was already in Mendeley. |
| c88 | Hard to parse | Now: "The ranking is therefore best read as crew knowledge organised by habitat, not as independent proof that a highly ranked site holds COTS." |
| c95 | Focus on field trials and recording exploration | Now: "Field trials will be needed to establish whether these zones find more COTS than expanding the search area at random, and that test depends on recording exploration effort and nil results as carefully as finds." |
| c98 | Section much shorter | "Retention, surveillance and the cost of deferral" is now two paragraphs instead of five. The first covers the trade-off for the recommended settings: ghost dives avoided per miss against random, follow-up tows saved, and habitat for exploration rather than tow thinning. The second covers what happens to deferred sites and the January trial. The site-rule detail, the tow-geometry paragraph and the long thinning discussion were removed; they remain in Appendices C and D. |
| c99 | Focus on the recommended product; trade-offs of density-aware classes, the dive rule and the tow zone; tone pass | The workflow section now has five parts: (1) what the workflow is and where it plugs in; (2) the scar-trigger question; (3) a new paragraph on why density-aware classes sharpen the choice, the A–D dive-rule result, the tow zone's different benefit, and the A–C alternative as a management choice; (4) one limits paragraph; (5) the trial design. |
| c100 | Unnecessary detail | The recent-CPUE and habitat-percentile paragraph is removed from the Discussion. The comparison stays in Appendix C. |
| c103 | Too dense; apply the tone skill | Rewritten as above. The A–B versus A–C zone comparison and the missed-dive counting caveats now sit in Appendix D. |
| c105 | Matthews et al. 2026 (dashboard paper) | Cited for GBRMPA's decision support system. |

## Voice pass on the Discussion

Every Discussion section was rewritten with `voice.md`:

- each paragraph opens with its claim, gives the evidence, and ends with what it means for crews or managers;
- caveats are gathered into one paragraph per section;
- numbers are cut to the headline ones.

The population-model section is now four paragraphs. It keeps the distinction about the wave-exposure term and the warning about excess COTS in 2026. The skill's checker flagged one Discussion sentence over 40 words, now split.

## Still open

- **The Summary still reports the site-selection result** (89.1% against 78.5% at random), which now sits in Appendix C. Do you want a workflow number there instead?
- Copyright owner, back-cover position title, and the JCU and QUT campus cities.
- The connectivity network results and figure.

## Follow-up (26 November)

- **Summary.** It now leads with the workflow back-test, as you asked.
  - The habitat-alone result is one sentence.
  - The site-selection numbers (89.1% against 78.5%) are dropped from the Summary; they remain in Appendix C.
  - The workflow paragraph adds low-density dives avoided (18.2%) and gives the random baseline as ghost dives per miss (8.2 against 2.1).
  - The Summary is 347 rendered words.
- **Conclusions.** The second paragraph now leads with the workflow result instead of the site-selection result.
