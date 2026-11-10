# CCIP Bridging Project technical-report template

This is the project-specific authority for `technical_report.qmd`, its Word/PDF
exports and the current review copy. Read it with `WRITING_BRIEF.md`. It applies
to both Codex and Claude work on this report. The separate administrative CCIP
Final Project Report form is a different deliverable.

## Retained source and evidence

- Source: `CCIP Final Report Template - Bridging Projects.docx` in this folder.
- SHA-256: `e88e42c58b01be83287e0848dc4ff68b3f5f4f0e3ccc27ba08793424f685a1f5`.
- Native Word rendering: 15 pages, eight section definitions. Inspection and
  package inventories are in `../cache/ccip_template/`.
- Never overwrite or resave the supplied template. Build from a working copy.
  Verify its checksum before and after every export.
- The template's instructional text and four comments govern content but are
  removed from the filled report. Its standard funding acknowledgement,
  Traditional Owner acknowledgement, disclaimer and branding are retained.
- Correct the template's spelling of “Discussion”; retain its intended heading.

## Content and length contract

| Section | Required focus | Template guidance |
|---|---|---|
| Project Summary | Summarise each strand, headline findings and benefits for people/environment; plain language | Maximum one page; draft target ≤360 words, confirm actual Word pagination |
| Context & Objectives | Challenge, prior work, literature and Work Plan research priorities; conclude with the three specific objectives as bullets | 2–3 pages |
| Methods | Data, modelling and validation supporting the objectives; applicable ethics/permits; dedicated **Engagement with Rights-holders and Management Stakeholders** subsection | 2–3 pages of text plus 1–3 key figures; modelling observation/domain map included |
| Results | Most important findings answering the objectives, with uncertainty and comparable baselines | 3–4 pages of text plus 2–4 key figures/tables/images **in total**; detailed comparisons in appendices |
| Discussion & Outputs | Interpret results against literature; limitations, assumptions and unexpected findings; conclude with a bullet list of final outputs and their actual status | 3–4 pages |
| Outcomes & Impact Pathway | Work Plan CCIP Outcomes, users, route from outputs to use, engagement and remaining actions | 1–2 pages of text; figures optional |
| Conclusions & Next Steps | Significant change enabled, principal findings and remaining research/actions | ≤300 words |
| Acknowledgements | Non-author contributions and any funding beyond the standard front-matter statement | Concise |
| Data Accessibility | Actual repositories/access links, FAIR registration, custodians and justified restrictions | Public links/DOIs must be supplied before final release; never invent them |
| References | Coral Reefs author–date style, alphabetic bibliography and full DOI links where available | `csl/coral-reefs-ccip.csl` |
| Appendices | Detailed methods, full comparison tables, validation, archived snapshots and scenario/field designs | Preserve scientifically necessary detail and distinguish old definitions |

Page ranges guide focus and density; do not pad sections or shrink template type
to hit a page count. Validate text length separately from maps, tables and
captions. Insert each visual directly after its first substantive citation.

The current Work Plan marks **CCIP Outcomes 2 and 4**. Do not claim Outcomes 1
or 3 were selected because related activities also support them. Research
priorities are **Decisions for Coral Resilience** and **Detection & Monitoring**.

## Formatting contract

The source DOCX controls appearance, including its teal/cyan heading colours;
generic document-style defaults must not replace them.

- Reuse the cover table, its 22 pt Arial title, 12 pt author treatment,
  affiliations and the original program/partner logos. Keep the existing report
  title; the long administrative subtitle need not be forced into the cover slot.
- Reuse the publication-information page and final gradient/contact/branding page.
  Fill author, contact and copyright slots only from supplied facts; otherwise
  retain explicit author-completion placeholders.
- Main A4 portrait: 11,900 ×16,840 twips; margins left/right 1,418 and
  top/bottom 1,701 twips. Cover section top margin is 2,268 twips. Preserve
  section-specific header/footer distances and the final-page variant.
- Main text uses template `BodyText1`/`BodyText`: Arial 11 pt. `BodyText1`
  has exact 14 pt lines and 12 pt after spacing; do not silently substitute
  generic paragraph spacing. Mathematical displays use expanding line height.
- Use actual Heading 1/2/3 styles. Heading 1 is 15 pt bold capitals in source
  teal (`31849B` plus its theme attributes); Heading 2 is 14 pt and Heading 3
  is 12 pt. Preserve their spacing and keep-with-next behaviour.
- Project Summary and appendix titles use `ReferencesandAppendix`; contents
  headings use `Contents`. The template's main headings are unnumbered.
- Use source `ListBullet`, `ListBullet2`, `NumberedList`, `TableHeading`,
  `TableText` and `Caption` styles. Table structure follows the source sample
  table; adapt column widths to content without shrinking body type. Wide
  appendix tables may reuse the template's landscape section pattern.
- Captions use the source 9 pt style and native Word `SEQ Figure` / `SEQ Table`
  fields. Use live REF fields for numerical cross-references and native contents,
  list-of-figures and list-of-tables fields. Refresh them in Word before PDF export.
- Project Summary begins at Arabic page 1; front matter uses the source Roman
  numbering. Footers identify `CCIP-BR-13`. Preserve cover first-page treatment
  and final-page furniture.
- Preserve source styles, theme, numbering definitions, original media, custom
  XML and opaque package parts. Add only relationships/styles needed by report
  content. Source guide/example pages are not report content.

Native Word may rename image parts and remove unused screenshots/maps belonging
to the discarded style-guide examples. This is expected cleanup. Verify the
retained cover/back-page branding by image-content checksum, not the ZIP member
name; the branding itself must remain byte-identical. The source DOCX is always
checked byte-for-byte. Word's field-cache/metadata normalisation is permitted;
critical style properties, page geometry and recurring furniture must not change.

## Editable slots and persistent metadata

The immutable template hash makes the following body-element locators stable.

| Source location in `word/document.xml` | Treatment in the filled copy |
|---|---|
| Body element 0, cover table rows 0–2 | Title and author text from Quarto metadata; affiliations from `ccip_report_metadata.json`; source row/run formatting reused |
| Paragraphs 4–6, 8 and 17 | Contact fields, citation authors/year/title and copyright; citation page count is a native NUMPAGES field |
| Paragraphs 10, 13 and 15 | Fixed funding, Traditional Owner and disclaimer text preserved |
| Contents and caption-list fields, elements 19–59 | Remove instruction shapes; populate live TOC/SEQ/REF fields and source list styles |
| Acronym heading/table, elements 60–64 | Fill the source table component with abbreviations used in this report |
| Main placeholders, elements 65–136 | Replace with report content; preserve the source portrait section pattern and Arabic page-1 restart |
| Guide/examples, elements 137–183 | Omit instructional content; reuse documented table, heading and landscape/portrait patterns where needed |
| Back page, elements 184–194 | Preserve gradient, logos and final-section geometry; fill the contact slot from report metadata |

Title/author changes belong in `technical_report.qmd`; affiliations, contacts,
copyright and citation-author spelling belong in `ccip_report_metadata.json`.
Never re-enter template-guidance pages as report content. Keep mathematical
displays editable, using expanding line height. Scale unusually tall figures
proportionally so the complete figure and caption fit within the source margins.
Move a wide table's caption/introduction with it when using a landscape section.
Express a requested section break as `pageBreakBefore` on its heading in the
assembled DOCX. Do not leave a separate break-only paragraph immediately before
a heading that already starts a page; it can generate an empty page.
After Word updates the contents and caption lists, compact only their empty
closing-field paragraph if it overflows. Preserve all visible entry formatting.
Keep short operational, coefficient and ranking-comparison tables together
with their captions rather than splitting a few rows across pages.

## Evidence and review safeguards

Keep results read from the existing CSV tables. Distinguish selected cull sites,
whole manta tows, tow distance, dive outcomes and annual effort. The exact A–C
rule differs from archived percentile cutoffs. Main Results uses one concise
operational table; the full retention tables and curve remain in Appendix H.
No model is refitted merely to reformat this report.

Engagement facts confirmed by the author on 27 October 2026: an April visit to
Malu and June follow-up; observations for a separate future Far Northern SDM;
multiple GBRMPA meetings, including the cull-site-ranking concept; a delivery
workshop planned for 5 November 2026; and trials planned aboard BPM in Cairns
and with Malu in the Far Northern GBR/Torres Strait in January–February 2027.
Distinguish completed engagement from planned trials and uncompleted model work.

Display equations reuse the source's borderless equation-table component,
with editable native Word mathematics and live equation numbers. Merge the
two example-expression cells when a single wider expression is needed; retain
the source alignment and allow expanding row height.

Use `scripts/render_ccip_report.ps1` for the complete export. A bare Quarto
render supplies the content intermediate, not the finished template assembly.
Native Word export controls the PDF so its layout matches the DOCX. Save and
inspect page images, check section lengths/display counts and the immutable
template hash. Create the next numbered review and both clean baselines only
after the final main document passes checks. Preserve annotated review inputs:
review 09 is the preserved feedback source; review 11 adds the subsequent ranking
experiment while review 10 is retained. Do not
use `-RefreshReview09` to overwrite the annotated input.

The packaged `render_docx.py` cannot run in this Windows environment because
its `pdf2image` dependency is absent. The verified local workflow uses a separate
native Word automation instance, followed by Poppler page images. It does not
attach to, save or close the user's existing Word documents.
