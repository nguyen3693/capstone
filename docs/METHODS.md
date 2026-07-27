# Methods (capstone)

## Paragraph for your report

Single-cell RNA sequencing data were obtained from NCBI GEO accession GSE288434 (10x Genomics Chromium 5′ v2). Raw gene–cell count matrices were downloaded from GEO supplementary files and organized into per-sample directories in 10x format (`barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz`). Analysis was performed in R (v4.6) using Seurat. For each of six samples (three healthy donors × vehicle [DMSO] or venetoclax [VX800] during CAR-T manufacturing), cells were filtered to retain those with 500–30,000 UMIs, 200–6,000 detected genes, and ≤15% mitochondrial reads. To reduce memory use on a local workstation, up to 3,000 cells per sample were randomly subsampled after filtering (random seed 42), then matrices were merged into a single object with sample, donor, and treatment metadata. Counts were log-normalized, and 1,500 highly variable genes were selected. Data were scaled, reduced with PCA (30 components), clustered with Louvain modularity optimization (resolution 0.5), and visualized with UMAP (30 PCA dimensions). Cluster marker genes were identified with Wilcoxon rank-sum tests (Seurat `FindAllMarkers`; positive markers only, min.pct = 0.25, log₂FC threshold 0.25). All steps were run from a single configuration file via `Rscript run_scrna_pipeline.R`; outputs included QC plots, UMAP figures, marker tables, and a processed Seurat object.

## Key parameters (from `config.yaml`)

| Step | Setting |
|------|---------|
| QC | `nFeature_RNA` 200–6000; `nCount_RNA` 500–30000; `%MT` ≤ 15 |
| Subsample | 3,000 cells/sample (enabled) |
| HVGs | 1,500 |
| PCA / UMAP dims | 30 |
| Clustering | resolution 0.5 |
| Markers | Wilcoxon; top 25 per cluster exported |

Adjust wording if you verify sample-level treatment labels on the GEO page or change subsampling/QC in config.
