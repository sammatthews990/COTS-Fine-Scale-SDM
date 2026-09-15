# Four-reef expert-zone review

Status: **awaiting expert annotations**

The comparison panels are complete, but accepted and rejected reef zones cannot
be inferred from model scores or map similarity. They require ecological and
operational judgement. Record zones in the accompanying CSV using recognizable
reef-sector descriptions and, where possible, approximate coordinates.

For each panel:

1. list zones whose relative priority is operationally plausible;
2. list zones whose priority is implausible or would misdirect a search team;
3. state whether geomorphic boundaries and within-zone texture are useful;
4. assign `accept`, `reject`, or `needs_discussion` only after the zone notes;
5. do not accept or reject solely because a candidate matches the incumbent.

The four panel families are:

- `support_120m_vs_200m`: legacy, Hmsc v2, v4 120 m hybrid and v4 200 m;
- `manta_augmentation`: legacy, Hmsc v2, v4 hybrid, COTS-only Hmsc and
  strong-Manta-response Hmsc;
- `support_60m_vs_120m`: legacy, Hmsc v2, v4 120 m hybrid and v4 60 m;
- `shared_state_vs_previous`: legacy, Hmsc v2, v4 120 m hybrid and the new
  shared-state 120 m candidate.

The 60 m panels are in
`../../hmsc_corrected_v4/reports/figures/v4_support_comparison_*_60m_120m.png`.
The shared-state panels are in `figures/shared_state_comparison_*.png`.

Until annotations are supplied, every expert gate remains `pending` and no
candidate can be promoted.
