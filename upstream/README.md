# Upstream: how the frozen inputs were made

Provenance only. These scripts built the files that `pipeline/config.R` pins
in the data store's `frozen/` and `predictors/` folders. They are kept as they
were run, with their original paths (`DataProcessing/...`, as at tag
`v1.0-2027-products`), and the pipeline never runs them. Rebuilding an input
means restoring its raw sources from the local data backup, rerunning the
script, and re-pinning the result in `config.R` explicitly.

| Folder | Built | Store path |
|---|---|---|
| `predictors/` | 30 m predictor stacks: `spatialData_PreProcessing.R` (terrain variables from the ACA/DeepReef DEMs, coral MaxEnt rasters, ReefGuidance layers → combined stack), `build_clean_predictor_stack.R` (12-layer clean stack), `rebuild_full_predictor_stack.R` (extended RG stack) | `predictors/` |
| `legacy_sdm/` | Clean + year and RG + year xgboost SDMs: `fit_predict_cots_sdm.R` via `refit_fy2526_clean.R` (run `20260820_1041`) and `refit_fy2526_rg.R` (run `20260820_1217`); `run_all_refits.ps1` | `frozen/legacy_sdm/` |
| `cots_sdm_v2/scripts/` | `01_prepare_data.R`: cleaned cull dives, site-years and first-visit site-years from the CSIRO workbook, manta history and SALAD; `13_assess_rank_signal.R`: out-of-fold annual ranks | `frozen/prepared/` |
| `cots_sdm_v2/iterations/hmsc_corrected_v4/` | Hmsc training tables (`04`, `04b`, `06`), corrected 120 m cells (`01a`) and manta footprints (`01c`) | `frozen/hmsc_v4/` |
| `cots_sdm_v2/iterations/hmsc_geomorphic_v3/` | Tow footprint links (`footprint_links_120m.rds`) | `frozen/hmsc_v3/` |
| `cots_sdm_v2/R/`, `config/` | Shared helpers and configuration sourced by the scripts above | |
| `cots_sdm_v2/reports/` | Exploratory decision record (model reviews, final synthesis) cited by `pipeline/README.md` | |
| `connectivity_lizard/` | Lizard larval-concentration tables (scripts `01`, `02`) | `frozen/connectivity_lizard/` |

Inputs not built here: the five `frozen/final_pins/` files were written by
pipeline scripts 15–20 and 02; the ADRIA/COTSMod snapshot came from the
ADRIA.jl and COTSMod-Calibration repositories; `raw/` holds source data as
received.
