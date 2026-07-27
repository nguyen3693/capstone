# Capstone: scRNA-seq UMAP + clustering (GSE288434)

Single-command Seurat pipeline: QC → clusters → UMAP → marker tables.

## Quick start

1. Put 10x folders under `data/raw/S01` … `S06` (see course notes / GEO download).
2. `cp metadata/sample_metadata.csv.example metadata/sample_metadata.csv` and edit labels.
3. `Rscript run_scrna_pipeline.R`
4. Outputs: `results/<timestamp>/` (plots, CSVs, `seurat_processed.rds`).

Full design and pseudocode: [docs/PIPELINE.md](docs/PIPELINE.md).  
**Methods paragraph for your write-up:** [docs/METHODS.md](docs/METHODS.md).  
Edit paths and thresholds in [config.yaml](config.yaml).
