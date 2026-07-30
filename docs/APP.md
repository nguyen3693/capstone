# Parameter playground (Shiny)

Compare inferred CD4 and CD8 populations side by side without re-running QC from raw 10x data.

## What it does

1. Loads a saved `results/*/seurat_processed.rds` from the pipeline.
2. Infers lineage from `CD4` versus mean `CD8A`/`CD8B` expression, then splits the cells before dimensionality reduction.
3. Uses **Run A for inferred CD4** and **Run B for inferred CD8**, with the same requested cell count for balance.
4. Recomputes lineage-specific variable genes → scaling → PCA → neighbors → clusters → UMAP for each run (progress bar + estimated time left).
5. Shows side-by-side UMAPs colored by cluster, treatment, sample, or donor.
6. Computes positive marker genes on demand, labels each cluster with its top genes, and exports the marker results as CSV.
7. Generates side-by-side CD4 and CD8 expression panels, with each run's fixed-gene violin plots stacked in one column.
8. Saves reproducible Run A or Run B UMAP bundles directly under `parameter-playground/` and optionally downloads them as ZIP files.

The inferred lineage display compares `CD4` expression with mean `CD8A`/`CD8B` expression. It is transcriptome-derived because this GEO series does not include FACS CD4/CD8 labels.

Marker genes describe how a cluster differs from the other cells; they are not automatically validated cell-type annotations.

## Requirements

R packages: `shiny`, `Seurat`, `ggplot2`, `patchwork` (plus what the pipeline already uses).

```r
install.packages(c("shiny", "patchwork"))
```

## Launch

From the **project root** (`capstone/`):

```bash
R -e "shiny::runApp('app', launch.browser = TRUE)"
```

Or in RStudio: open `app/app.R` and click **Run App**.

## Suggested first comparison

1. Click **Load object** (defaults to newest `seurat_processed.rds`).
2. Choose the number of cells per lineage.
3. Run A (CD4): resolution `0.5`, dims `30`.
4. Run B (CD8): resolution `0.5`, dims `30`.
5. Color by **Clusters**, then switch to **Treatment**.
6. Ask: which transcriptional states are present within CD4 versus CD8 cells, and how are DMSO/VX800 cells distributed?
7. Select a run, choose how many genes to show, and click **Compute cluster markers**.
8. Review the marker-labeled UMAP and download the complete CSV if needed.
9. After both runs finish, click **Generate CD4 vs CD8 expression panels** to compare the fixed genes side by side.
10. Use **Save Run A/B outputs** to write a timestamped project folder, or use the adjacent ZIP button for a browser download.

## Exported UMAP bundle

Each export creates `parameter-playground/Run-A-YYYYMMDD-HHMMSS/` or the corresponding Run B folder containing:

- `umap_clusters.png`
- `umap_donor.png`
- `umap_treatment.png`
- `umap_sample.png`
- `parameters.txt`

`parameters.txt` records the lineage, source RDS, cell counts, subsampling status, resolution, PCA/UMAP dimensions, random seed, number of clusters, runtime, and workflow. The entire `parameter-playground/` directory is gitignored.

## Notes

- This app does **not** re-filter QC; it infers CD4/CD8 lineage and can randomly subsample each lineage from an already-processed object.
- CD4/CD8 calls are transcriptome-derived, not FACS-validated labels.
- Reclustering 18k cells may take ~30–90 seconds per run on a laptop.
- Marker calculation can take several minutes and only runs when requested.
- The fixed violin panel is `CD4`, `CD8A`, `GZMA`, `GZMB`, `IL7R`, `IL2RA`, `TCF7`, `FOXP3`, `TIGIT`, `CCR7`, `BCL2`, and `MCL1`.
- The app reports panel genes that are missing from the loaded dataset.
- Keep `results/` gitignored; the app only needs a local RDS.
