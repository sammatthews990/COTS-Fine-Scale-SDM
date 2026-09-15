# Blinded expert review protocol

## Purpose

This review tests whether candidate fine-scale surfaces contain spatial
patterns that are credible to people with operational knowledge of the reefs.
It is a secondary acceptance gate alongside held-out observations. It does not
replace cross-validation and it must not be based on probability range or map
colour contrast.

## Display standardisation

Every candidate is converted to a percentile rank across valid pixels in the
display region. The same fixed colour scale is then used for every surface.
This exposes spatial ordering while removing any advantage from a wider or more
visually convenient raw probability range.

No COTS observations, model names, or reef outbreak probabilities are shown in
the blind review. A separate key is written beside the review file and should
be withheld until scores are exported.

## Candidate set

The current five-surface set contains:

- the legacy year surface;
- the legacy agnostic surface;
- the lifetime RG-with-clean-fallback surface;
- the annual-response RG-with-clean-fallback surface; and
- the existing v2 static operational surface.

The order is randomised once and held constant among regions and reviewers.

## Review questions

For every region and anonymous surface, score from 1 (poor) to 5 (strong):

- known high-density areas receive relatively high rank;
- fine-scale reef pattern is ecologically plausible;
- implausible or operationally misleading hotspots are uncommon; and
- the ordering would be useful for prioritising dives within a selected reef.

Reviewers also record confidence and free-text notes identifying specific
flanks, habitats, or false hotspots.

## Independence rules

- Complete the review before opening `expert_blind_review_key.csv`.
- Use at least two reviewers where possible.
- Reviewers should score independently before discussing disagreements.
- Do not alter model settings in response to individual review regions.
- If expert feedback guides a model change, evaluate that change on a separate
  preregistered set of reefs or a later observation period.

## Decision summary

Summarise scores by anonymous surface first. Reveal identities only after the
summary is frozen. A candidate should not be accepted when its metric gain is
small and its blind operational score is materially worse than the legacy-year
benchmark.
