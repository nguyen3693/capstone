# Parameter playground (Shiny)

Compare clustering / UMAP settings side by side without re-running QC from raw 10x data.

## What it does

1. Loads a saved `results/*/seurat_processed.rds` from the pipeline.
3. Choose **total cells** (subsample for speed) and set **resolution**, **PCA/UMAP dims**, and **seed** for **Run A** and **Run B**.
4. Recomputes neighbors → clusters → UMAP for each run (progress bar + estimated time left).
5. Optionally split into **inferred CD4+ vs CD8+** UMAPs (from `CD4` vs `CD8A`/`CD8B` expression; this GEO series has no FACS lineage labels).

## Suggested first comparison

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
2. Run A: resolution `0.5`, dims `30`.
3. Run B: resolution `0.8`, dims `30`.
4. Color by **Clusters**, then switch to **Treatment**.
5. Ask: does the DMSO-heavy island stay intact? do clusters split/merge?

## Notes

- This app does **not** re-filter QC or re-subsample; it explores clustering geometry on an already-processed object.
- Reclustering 18k cells may take ~30–90 seconds per run on a laptop.
- Keep `results/` gitignored; the app only needs a local RDS.
