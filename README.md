# Capstone: scRNA-seq UMAP + clustering (GSE288434)

An R/Seurat workflow and interactive Shiny dashboard for exploring CAR-T single-cell RNA-seq data from NCBI GEO accession [GSE288434](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE288434).

The command-line pipeline performs QC, normalization, dimensionality reduction, clustering, UMAP, and marker detection. The dashboard then compares inferred CD4 and CD8 populations using adjustable clustering parameters.

## Pipeline quick start

1. Put 10x folders under `data/raw/S01` … `S06` (see course notes / GEO download).
2. `cp metadata/sample_metadata.csv.example metadata/sample_metadata.csv` and edit labels.
3. `Rscript run_scrna_pipeline.R`
4. Outputs: `results/<timestamp>/` (plots, CSVs, `seurat_processed.rds`).

Edit paths and thresholds in [config.yaml](config.yaml). Full pipeline design and pseudocode are in [docs/PIPELINE.md](docs/PIPELINE.md).

## Interactive parameter playground

Launch from the project root:

```bash
R -e "shiny::runApp('app', launch.browser = TRUE)"
```

The dashboard:

- Infers CD4 versus CD8 lineage from `CD4`, `CD8A`, and `CD8B` expression
- Computes separate lineage-specific PCA, clusters, and UMAPs
- Compares Run A (inferred CD4) with Run B (inferred CD8)
- Adjusts cells per lineage, clustering resolution, dimensions, and random seed
- Colors UMAPs by cluster, donor, treatment, or sample
- Calculates cluster marker genes and labels clusters with top markers
- Compares a fixed T-cell expression panel with side-by-side violin plots
- Exports reproducible UMAP bundles with PNGs and `parameters.txt`

Saved dashboard bundles are written to:

```text
parameter-playground/Run-A-YYYYMMDD-HHMMSS/
parameter-playground/Run-B-YYYYMMDD-HHMMSS/
```

Each bundle contains UMAPs colored by cluster, donor, treatment, and sample, plus the settings used to generate them. See [docs/APP.md](docs/APP.md) for complete instructions.

## Documentation

- [Pipeline design](docs/PIPELINE.md)
- [Shiny app guide](docs/APP.md)
- [Assignment methods paragraph](docs/METHODS.md)

## Data and generated outputs

Large files are intentionally excluded from Git:

- `data/`
- `results/`
- `parameter-playground/`

Download the original matrices from GEO rather than committing them to this repository.
