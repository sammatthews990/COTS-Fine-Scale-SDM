# COTS Fine Scale SDM

Fine-scale habitat and search-planning layers for crown-of-thorns starfish
(COTS) control on the Great Barrier Reef (GBRF CCIP Bridging Project
CCIP-BR-13). The models predict where on a reef COTS are most likely to be
found, so the COTS Control Program can manta tow less and search fewer cull
sites without missing problematic densities.

## Layout

| Folder | Contents |
|---|---|
| [`pipeline/`](pipeline/) | The delivered pipeline: configuration, scripts, tables, figures, the 2027 cull-site rankings and the technical report. **Start with [`pipeline/README.md`](pipeline/README.md)**, which records every result, caveat and the run order. |
| [`upstream/`](upstream/) | Provenance only: the scripts that built the frozen inputs (predictor stacks, legacy SDM refits, data preparation, Hmsc training tables, Lizard connectivity) and the exploratory decision reports. Not run by the pipeline. |
| [`AGENTS.md`](AGENTS.md) | Goals, current state and rules for anyone (person or coding agent) working in the repository. |

## Data

Inputs are not stored in the repository. They are in the project data store,
the CCIP-BR-13 SharePoint `Data` folder synced through OneDrive (set
`COTS_DATA_DIR` to use a copy elsewhere): `raw/`, `predictors/`, `frozen/` and
`products/`, about 12 GB, every file listed with its SHA-256 in
`DATA_MANIFEST.csv`. See "Data store" in [`pipeline/README.md`](pipeline/README.md).

## Running

From the repository root, with R 4.5:

```powershell
Rscript pipeline/scripts/00_check_inputs.R   # verify every pinned input in the data store
```

then the scripts in the order given in `pipeline/README.md` ("Run order").

## History

Tag `v1.0-2027-products` holds the full repository, including all exploratory
iterations, as it produced the technical report and the 2027 products. The
reorganisation into `pipeline/`, `upstream/` and the external data store
followed on 10 November 2026.
