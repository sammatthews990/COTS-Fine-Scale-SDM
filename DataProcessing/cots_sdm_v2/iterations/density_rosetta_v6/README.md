# COTS density and RosettaCOTS v6

This isolated final candidate family predicts continuous COTS density rather
than optimizing a binary CPUE threshold.

Primary candidates:

- native Cull counts with dive-minute offsets, producing expected Cull CPUE;
- native SALAD `Total` counts with swept-area offsets, producing expected COTS
  plus scars per hectare;
- an integrated SALAD-equivalent density model. SALAD uses its native count and
  area likelihood. Cull uses the published RosettaCOTS Cull-from-SALAD
  calibration inside a separate negative-binomial observation process.

All habitat predictors remain static and deployable. Reef, year, outbreak
probability, method, realised effort and Cull history are excluded from the
habitat score. The fitted estimand is conditional on an active reef-year;
outbreak activation remains separate.

RosettaCOTS 1.0.2 is vendored at commit
`ef4e08e7566f3c849a9da905733e8b4cb06ca5de`. Its calibration uncertainty is
reported and tested; converted point estimates are not treated as truth.
