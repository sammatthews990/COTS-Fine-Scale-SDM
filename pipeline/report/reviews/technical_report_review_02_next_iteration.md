# Methods and Results revision queue — review 02

The Background and dependent management discussion have been revised. This note
records work for the next iteration; no models, thresholds in existing result
tables, or reported performance estimates were changed in this pass.

## Methods

1. Separate the response used for existing model fitting and validation (CPUE
   ≥0.02) from the revised operational loss criterion (>0.04 COTS per minute of
   bottom time). Specify three positive-density groups: <0.02, 0.02–0.04 inclusive,
   and >0.04; retain ghost dives as a separate zero group. Confirm the treatment
   of exactly 0.02 in legacy versus new bins rather than silently changing labels.
   Some intermediate dives can be missed; the penalty for missed high-density
   dives applies above 0.04. Agree its numerical tolerance and whether that is a
   share of dives, sites, aggregations or COTS. Do not infer aggregation retention
   from dive retention. Coral-dependent ecological targets need separate wording.
2. Re-evaluate held-out rankings at the revised >0.04 outcome before considering
   any refit. Compare retention and effort at matched effort and matched
   high-density retention, with random selection and systematic manta thinning
   where appropriate. Keep spatial and temporal validation limitations explicit.
3. Define site searches, tow effort, dive effort, and annual coverage separately.
   Record added exploration searches, deferred repeat searches, and annual or
   pre-closure checks. Distinguish redistributed effort from total effort saved.
   Compare field performance with crews' normal planning, not only random choice.
4. Describe why perimeter-defined sites can leave lagoon bommie fields, shallow
   crests and isolated patches unsearched. Specify how routes and exploration
   areas will be tested; retrospective thinning is not validation of a new route.
5. The former Background 1.5 has been removed. Its purpose and broad spatial scope
   are integrated into 1.2–1.4; implementation detail belongs in the metacommunity
   Methods subsection, which now retains the `sec-metacommunity` anchor. Check
   relocated habitable-area, disturbance, site-allocation and coupling descriptions
   against the final model release. Confirm connectivity attribution (annual
   2018–2022 matrices), output F delivery, and habitat-coupling plans. The current
   model still uses reef-level transport, not implemented sub-reef connectivity.
6. Specify paired 2012–present culling/no-culling histories and external-boundary
   assumptions once finalised. Future scenarios must include stopping when local
   coral cannot be saved versus continuing to reduce larval export, when to start
   or intensify control on source reefs, reallocating existing effort, and extra
   vessel capacity. Specify catchability/effort, feasible removal constraints,
   climate sequences, seeds and uncertainty before claiming operational efficacy.

## Results

1. Preserve all current ≥0.02 results as legacy analyses. Add >0.04 retention and
   missed-dive outcomes, low-density (<0.02) and ghost dives avoided, and the
   intermediate group retained/deferred, with uncertainty and baselines. Revisit
   scripts 15–17 and 19, relevant tables/figures and class summaries together.
   Do not relabel current gains as evidence for the revised criterion. Only
   refit or change deployed priorities after an explicit analysis decision.
2. Report the 10–15% site-search and low-density/empty-dive targets as proposed
   trial adoption criteria until tested. Keep the January 2027 trial distinct
   from retrospective validation. Quantify how later coverage checks and added
   exploration affect seasonal effort and missed high-density populations.
3. Distinguish evidence for ranking new ground from evidence that proposed
   exploration areas discover more COTS. The latter requires observations; the
   existing comparison is partly confounded by irruption stage.
4. Keep calibration peak matching, historical-removal diagnostics and independent
   validation separate. Update the frozen 15 October metacommunity snapshot only
   through pinned source inputs and script 22. Future historical impact and
   fifth-wave scenarios need results for populations, larval export, local coral,
   downstream coral and required resources, including infeasible removals.
5. After new analyses, update README, the executive summary, recommendations,
   Appendix E and dependent Discussion claims. Follow the separate review scope
   before propagating to the manuscript.

## Source and editorial checks

- The approximate 1,300 sites/year and 200 actioned reefs/year are author-supplied
  dashboard summaries. Save the reporting period and dashboard extract before
  finalisation. Register totals are a different quantity. The public dashboard
  landing page reports 150–250 reefs annually but does not expose the site total.
- Use irruption in narrative. Keep original published reference titles, named
  dashboard categories, source-variable names and R lookup strings unchanged.
- The official Work Plan objective wording is preserved apart from terminology;
  operational refinements are explained below it, not represented as contractual
  amendments. Confirm exact quoted Work Plan wording for the final submission.
- The numerical definition of negligible missed >0.04 dives remains open.

## Implementation status after review 03 — 16 October 2026

The queue above is retained as the original request. The following work is now
incorporated in the main Quarto source and regenerated report:

- **Methods 1–4:** thresholds and exact boundary treatment are explicit; script
  23 evaluates unchanged rankings at >0.04; candidate-site, tow, dive and annual
  effort are distinguished; the trial records exploration and deferred-site checks.
- **Methods 5:** the relocated model description retains the frozen source
  snapshot. No final model release or new fine-scale transport implementation
  has been assumed. Source attribution and delivery checks remain open below.
- **Methods 6:** a dedicated scenario-design subsection specifies paired
  histories, external-boundary alternatives, future starting/stopping and
  source-reef controls, resource comparisons and required effort/catchability.
  These remain planned experiments, not completed simulations.
- **Results 1:** new fixed-effort, matched-retention, paired-difference,
  intermediate-outcome and season-specific tables are in `outputs/tables/`.
  Script 23 uses the existing scores from 15–17 and 19 without changing their
  original >=0.02 tables or deployed rankings. Figure 20 and a separate class
  table report the new criterion. Same-effort random expectations and 500 paired
  reef-year bootstrap intervals are included. Manta systematic-thinning results
  remain a separate aggregation outcome, not a dive-CPUE comparison.
- **Results 2–3:** proposed adoption targets and the distinction between new
  ground evidence and untested exploration areas are explicit. Annual effort
  and new-area yields cannot be calculated without field observations.
- **Results 4:** calibration, diagnostics and validation remain separate. The
  metacommunity snapshot is unchanged; no program-impact or suppression benefits
  were invented or inferred from it.
- **Results 5:** README, summary, Discussion, recommendations and Appendix E have
  been aligned with the new assessment. The manuscript remains outside this pass.

Still open: a numerical negligible-omission tolerance and its operational unit;
field outcomes and seasonal effort accounting; output F/connectivity attribution,
delivery and habitat-coupling plans; final model release and effort/catchability
specification; scenario horizons, budgets and forcing agreed with partners; the
dashboard's approximate annual site-count period and supporting extract; data
access and remaining source attribution. These require new information or work,
rather than editorial integration of the existing evidence.
