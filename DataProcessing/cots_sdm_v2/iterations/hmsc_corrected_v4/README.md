# Hmsc corrected v4

This iteration is an isolated correction of the geomorphology-led Hmsc work.
It does not modify the frozen legacy clean-year model, Hmsc v2, or Hmsc v3.

The primary changes are:

- restore conditioning on reef-years with credible COTS/scar evidence;
- classify manta evidence into strong positive, weak scar, and strict negative;
- use manta evidence at whole-footprint support rather than assigning a tow
  positive to every crossed cell;
- restore scars, continuous SALAD density, and susceptible-coral responses;
- compare actual selected terrain scales with one/two-PC and full shrinkage
  representations;
- fit enhanced and clean Hmsc surfaces separately, with the clean surface used
  as a calibrated fallback wherever ReefGuide/geomorphic layers are absent;
- compare candidates at their fitted support against Hmsc v2 and the frozen
  clean-year incumbent.

Run scripts in numeric order. The first implementation uses the already
audited exact 120 m v3 footprint links, then adds a matched 200 m support once
the corrected response/conditioning specification passes the initial gates.
