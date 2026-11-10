# Review 12 responses

Annotated input: `technical_report_review_12.docx` (SHA-256
`d6d3e61f7f6f03fa9b77a034f3074cc1502f57f798774b80d023852d2484cb26`).
It and `technical_report_review_12_baseline.*` remain unchanged. The tracked
text matched the baseline apart from the tracked edits. The main Quarto report is
the authoring source.

Review 12 was generated before the Appendix N edit of 2 November ("Dive rule
A–C + follow-up tow zone" row and the sentence on additive zone effects). That
edit is kept.

## Comments

| ID | Feedback | Response |
|---|---|---|
| 36 | Check the two sections added by Sev/Owen; add references to the reference list if missing. | Incorporated. Five references were verified against Crossref/DataCite and added to `references.bib`: Stewart et al. (2026, *Coral Reefs* 45:1581–1599, doi:10.1007/s00338-026-02818-x), Pratchett et al. (2017, *Diversity* 9:2), Moneghetti et al. (2019, *Ecology* 100:e02730), Engwirda (2018, *J. Comput. Phys.* 375:155–176) and Beaman (2010, MTSRF report). Text checks and open items are below. |

## Tracked passages

- **Title:** now "CCIP-BR-13: Integrative fine-scale habitat suitability and
  metapopulation modelling for targeted Crown-of-Thorns Starfish management".
  The old title moved to the subtitle. The suggested citation and the export
  checker's required phrase were updated to match.
- **Water movement and larval connections (Methods):** the old
  larval-concentration paragraph and the `[TO WRITE]` placeholder were deleted,
  as tracked. Sev and Owen's two subsections were inserted with light edits:
  - "save in netcdf" became "saved in NetCDF format";
  - "This forms a network … can be formed" was repaired;
  - $C_{i,j}(t)$ is now an inline equation;
  - "Connectivity Modelling" is now sentence case;
  - the observation datasets are a list;
  - "crown-of-thorns starfish (COTS)" became "COTS", since it is defined in
    the Summary;
  - citations now use bib keys.
- **Appendix K placeholder:** "alongside the Methods placeholders" was removed,
  since the Methods text now exists.
- **Reference casing (not requested):** the CSL sentence-cases titles. It had
  rendered "acanthaster cf. Solaris … great barrier reef" in four entries,
  including SALAD. Species and place names are now brace-protected in the bib.
  Technical-report entries still render without their institution; that is
  existing CSL behaviour and was left unchanged.

Export: `cache/ccip_template/build_20261105_143514`, 121 pages. The only
warning is the existing one: Discussion is 5 pages against 3–4 guidance.
The next clean review copy is `technical_report_review_13.docx`.

## Follow-up answers (5 November)

- **Bathymetry:** now cites Beaman (2017), the 30 m Geoscience Australia
  dataset. Beaman (2010) is removed from the bib.
- **"Appendix X":** Sev and Owen's text is now Appendix O, "Larval connectivity
  methods", with these edits:
  - the tense is now past;
  - "settlement radii is" became "no settlement radius was applied", and
    "larva are" became "larvae were";
  - "Figure 2 in the main text" became "Figure 2 of that paper", meaning
    Stewart et al. 2026, because the report's own Figure 2 is the Batt/Eyrie
    map.
- **Engwirda (2018):** confirmed.
- **Tense:** the hydrodynamic forcing sentences and the connectivity subsection
  now use the past tense. The domain extent and the description of the matrix
  stay in the present tense.
- **Larval-concentration method:** a three-sentence version of the deleted
  paragraph is back at the end of the particle-transport subsection, pointing
  to Appendix K.
- **Accepted as is:** no model skill, no tool references and unnamed rivers.

Export: `cache/ccip_template/build_20261105_144418`, 122 pages. The review 13
copies were refreshed from this export.

## Original open items (answered above)

- **Model skill is not reported.** Performance against IMOS, AIMS and EcoRAP
  observations is described but no skill is given. Add a sentence or
  appendix table if available.
- **References for the remaining tools.** COMPAS–EMS, eReefs GBR4 H4p0, BARRA
  and the Source catchment model have no references. The two rivers are not
  named.
- **Tense.** The connectivity subsection is in the present tense ("are then
  used", "We begin"); the rest of Methods is in the past tense. Left as written.
- **Larval-concentration method.** With the old paragraph deleted, Methods no
  longer describes the larval-concentration ranking. The dispersal Results and
  Appendix K still report it (ages 5, 9 and 10 days, 10 releases). Appendix K
  summarises the method briefly. Re-add one sentence in Methods if wanted.
