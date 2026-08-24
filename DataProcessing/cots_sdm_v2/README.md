# COTS fine-scale SDM v2

This directory is an isolated refactor of the COTS modelling workflow. It does
not source, overwrite, or modify the legacy model scripts or their outputs.

The refactor separates three questions that were previously mixed together:

1. **Static habitat prior** — where can high-density COTS be supported?
2. **Dynamic outbreak risk** — which reefs are likely to be occupied or
   outbreaking in a target year?
3. **Conditional high density** — given COTS presence, where is CPUE likely to
   exceed 0.02?

The first implemented milestone is a controlled benchmark. All candidate
models use the same site-year target, rows, spatial folds, and XGBoost settings.
Only the feature set changes:

- static fine-scale habitat predictors;
- static predictors plus calendar year;
- static predictors plus lagged reef outbreak probability;
- static predictors plus both year and lagged outbreak probability.

Lagged adjusted-hindcast features are the leakage-safe default. A separate
provisional operational experiment trains with same-year adjusted hindcasts
through 2025 and deploys the adjusted 2026 forecast. It remains explicitly
provisional until upstream hindcast and forecast provenance is confirmed.

## Run from the repository root

```powershell
Rscript DataProcessing/cots_sdm_v2/scripts/00_check_inputs.R
Rscript DataProcessing/cots_sdm_v2/scripts/01_prepare_data.R
Rscript DataProcessing/cots_sdm_v2/scripts/02_fit_benchmark.R
Rscript DataProcessing/cots_sdm_v2/scripts/03_fit_rolling_origin.R
Rscript DataProcessing/cots_sdm_v2/scripts/04_diagnose_targets.R
Rscript DataProcessing/cots_sdm_v2/scripts/05_prepare_manta_negatives.R
Rscript DataProcessing/cots_sdm_v2/scripts/06_fit_operational_2026.R
Rscript DataProcessing/cots_sdm_v2/tests/smoke_test.R
```

`01_prepare_data.R` writes derived tables under `artifacts/`. The benchmark
writes fold predictions and metrics under `outputs/`, and fitted fold models
under `models/`. These generated files are ignored by Git.

## Data contracts

- Cull observations are read from the `Cull` sheet of the configured workbook.
- `site_year_all` pools all dives at a site within a calendar year.
- `site_year_first_visit` pools only dives on the first surveyed date at a site
  within a calendar year, reducing within-year depletion bias.
- `high_density` is based on pooled CPUE (`total COTS / bottom time >= 0.02`).
  `CPUE_max` is retained only as a diagnostic and is never the primary target.
- Predictor extraction preserves an immutable `site_id` and extraction ID.
- All observations from one reef are assigned to the same spatial fold.
- Forward-year validation trains only on years before the test year.
- Adjusted hindcasts (1991-2025) and the adjusted 2026 forecast are ingested as
  separate products; the forecast never enters lag-only benchmark training.
- Repeated manta non-detections are prepared as presence-stage evidence with
  method conflicts and effort retained. They are not appended as unweighted
  pseudo-absences to the high-density target.

See `IMPLEMENTATION_PLAN.md` for the staged roadmap and decision gates, and
`BENCHMARK_FINDINGS.md` for the first controlled and rolling-origin results.
