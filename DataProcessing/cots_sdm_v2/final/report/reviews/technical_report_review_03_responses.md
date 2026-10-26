# Review 03 responses and Methods/Results integration

All annotated comments and prose edits were integrated into the Quarto source.
The annotated review 03 is preserved. Its Methods/Results baseline differences
outside tracked prose were Word's minus-sign normalisation in unit exponents;
the source retains correctly rendered mathematical units.

| Comment or edit | Response |
|---|---|
| Executive-summary tracked edit | Added “habitat suitability” to the delivered-layer description. Summary now also reports the new >0.04 retention trade-off. |
| Objective 3 tracked edit | Repaired the incomplete wording to “Refine site-level management strategies for individual reefs through predicted suitable habitat and metapopulation modelling.” The list is presented as project objectives, rather than a verbatim contractual quotation. |
| Background tracked deletion | Removed the sentence explaining that the management questions were refinements of the Work Plan. Scientific status remains explicit in the surrounding text. |
| 32 | Spelled out Overview, Data, Model, Assessment and Prediction (ODMAP), explaining what each part contributes and where the protocol applies. |
| 35 | Removed the detailed register-matching algorithm from the Methods narrative. It remains in the analysis code. |
| 39 | Counted the relevant scope from the source records: 4,299 matched designated sites on 339 reefs; 361 reefs in the wider cleaned source. Report numbers read `cull_record_scope.csv` inline. |
| 41 | Incorporated the confirmed 2002–2026 full source period, consistent with the cleaned records. Data-access wording remains an author check. |
| 42 | Distinguished the original ≥0.02 modelling response from the >0.04 omission criterion; defined ghost, positive low-density, intermediate and high-density outcomes and the exact boundaries. |

The previous revision queue is now implemented wherever the existing evidence
permits it. Script 23 evaluates saved spatial and forward rankings without model
refits or changes to the 2027 product. New tables, Figure 20 and a high-density
class summary report fixed effort, matched retention, paired uncertainty,
intermediate outcomes and per-season results. Original ≥0.02 results are retained
and labelled. Methods explain candidate-site effort, random expectations,
bootstrap/interpolation, field and annual coverage accounting, and the design of
the proposed historical and future metacommunity assessments.

The new forward result is materially more restrictive: retaining 95% of >0.04
dives requires about 90% of candidate site-years, while 99% requires about 98%.
Neither target defines “negligible”, and no field adoption target is declared met.
The report, README, brief, Discussion and recommendations reflect that result.

The appended implementation status in `technical_report_review_02_next_iteration.md`
identifies remaining work requiring new information: field observations, a
numerical omission tolerance, annual effort accounting, connectivity attribution
and delivery, final model release/effort modelling, agreed scenario settings,
dashboard source period and access statements. The metacommunity Results remain
the frozen 15 October snapshot. The manuscript is unchanged in this pass.

Continue review in `technical_report_review_04.docx` with Track Changes enabled.
