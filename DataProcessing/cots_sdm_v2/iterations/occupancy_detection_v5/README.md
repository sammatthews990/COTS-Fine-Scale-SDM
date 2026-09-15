# COTS occupancy/detection v5

This isolated iteration develops a shared ecological-state model with separate
observation processes for Cull, Manta, SALAD and RHIS. It does not modify the
legacy clean + year incumbent or the corrected Hmsc v4 artifacts.

Core principles:

- static habitat estimates the shared within-reef ecological state;
- reef-year outbreak probability is an external activation layer only;
- method, effort, visibility, scars and Cull history belong to observation or
  intervention processes and cannot enter the deployable habitat score;
- positive Manta tows remain whole-track observations and are never copied to
  each crossed cell;
- strong Manta evidence is direct COTS, at least two counted scars, or legacy
  common scars; single/present scars are sensitivity evidence only;
- promotion requires held-cell, leave-reef, forward-year and expert-map gates.

The v3 `unmarked` model is retained as a predecessor only. It excluded positive
Manta tracks, had very few repeated sites, estimated occupancy close to one and
did not converge.
