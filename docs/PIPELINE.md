# scRNA-seq pipeline (MVP)

## Goal

One command runs QC → clustering → UMAP → marker export for 10x-style count matrices (GSE288434).

## Folder layout

```text
capstone/
├── config.yaml                 # thresholds and paths (edit this)
├── run_scrna_pipeline.R        # single entry point
├── metadata/
│   ├── sample_metadata.csv.example
│   └── sample_metadata.csv     # you create (not committed if you prefer)
├── data/raw/                   # gitignored; 10x folders S01 … S06
│   ├── S01/
│   │   ├── barcodes.tsv.gz
│   │   ├── features.tsv.gz
│   │   └── matrix.mtx.gz
│   └── ...
├── R/
│   ├── pipeline_utils.R
│   ├── 01_load_merge.R
│   ├── 02_qc_filter.R
│   ├── 03_normalize_cluster.R
│   └── 04_markers_export.R
├── results/                    # gitignored or commit only small summaries
│   └── <run_id>/
│       ├── plots/
│       ├── tables/
│       ├── seurat_processed.rds
│       ├── config_used.yaml
│       └── run_summary.txt
└── docs/PIPELINE.md
```

## Before first run

1. Organize GEO matrices into `data/raw/S01` … `S06` (standard 10x names).
2. Copy metadata template:
   ```bash
   cp metadata/sample_metadata.csv.example metadata/sample_metadata.csv
   ```
3. Confirm `donor` and `treatment` labels match the paper (DMSO vs venetoclax / VX800).
4. Activate your R environment and install packages: `Seurat`, `ggplot2`, `dplyr`, `yaml`, `readr`.

## Run

From project root:

```bash
Rscript run_scrna_pipeline.R
```

Optional custom config path:

```bash
Rscript run_scrna_pipeline.R path/to/config.yaml
```

## Pipeline stages (biology + code)

| Stage | What it does | Module |
|-------|----------------|--------|
| Load / merge | Read 10x counts per sample; merge; add sample metadata | `01_load_merge.R` |
| QC | Mitochondrial %; violin plots; filter low-quality cells | `02_qc_filter.R` |
| Normalize & embed | Log normalize, HVGs, scale, PCA, graph clustering, UMAP | `03_normalize_cluster.R` |
| Markers | Wilcoxon markers per cluster; CSV export | `04_markers_export.R` |

## Pseudocode (full flow)

```text
READ config.yaml
CREATE results/<run_id>/

meta <- READ sample_metadata.csv
FOR each row in meta:
    counts <- Read10X(data/raw/<folder_name>)
    obj_sample <- CreateSeuratObject(counts)
    TAG obj_sample with donor, treatment, sample_id
MERGE all obj_sample -> obj

obj <- ADD percent.mt
PLOT qc violin -> results/.../qc_violin_pre_filter.png
obj <- FILTER by min/max genes, counts, %MT from config

obj <- NormalizeData
obj <- FindVariableFeatures
obj <- ScaleData (regress percent.mt)
obj <- RunPCA
obj <- FindNeighbors
obj <- FindClusters(resolution from config)
obj <- RunUMAP

PLOT UMAP by cluster, treatment, sample
markers <- FindAllMarkers
WRITE cluster_markers_all.csv, cluster_markers_top.csv
SAVE seurat_processed.rds
WRITE run_summary.txt
```

## MVP vs later (v2)

**In scope now**

- One config file
- Multi-sample merge with treatment labels
- Standard Seurat workflow
- Reproducible output folder per run

**Later**

- Integration / batch correction (e.g. SCT, Harmony)
- Automatic cell-type labels from marker panels
- Condition-level DE (pseudobulk)
- CLI flags (`--samples S01,S02`)
- HTML QC report (R Markdown)

## Assignment-friendly aims (simplified)

1. Implement a reproducible pipeline that processes GSE288434 10x data from a single config.
2. Produce UMAP visualizations colored by cluster and treatment.
3. Export cluster marker genes and document QC filtering in `run_summary.txt`.

## Troubleshooting

| Issue | Check |
|-------|--------|
| `Sample folder not found` | `folder_name` in metadata matches directory under `data/raw` |
| `Metadata file not found` | Copy `.example` to `sample_metadata.csv` |
| Very few cells after QC | Loosen `max_percent_mt` or `min_counts` in `config.yaml` |
| One giant cluster | Increase `cluster_resolution` (e.g. 0.6–1.0) |
