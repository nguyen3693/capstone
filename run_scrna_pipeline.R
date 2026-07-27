#!/usr/bin/env Rscript
# Single-entry scRNA-seq pipeline (Seurat).
# Usage (from project root):
#   Rscript run_scrna_pipeline.R
#   Rscript run_scrna_pipeline.R config.yaml

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config.yaml"

source("R/pipeline_utils.R")
source("R/01_load_merge.R")
source("R/02_qc_filter.R")
source("R/03_normalize_cluster.R")
source("R/04_markers_export.R")

`%||%` <- function(x, y) if (is.null(x)) y else x

require_config_packages()
cfg <- load_config(config_path)

run_info <- make_run_dir(cfg$paths$output_dir, cfg$project$run_id)
run_dir <- run_info$run_dir
message("Run directory: ", run_dir)

meta <- read_sample_metadata(cfg$paths$metadata_file, cfg$samples$include %||% character())

# --- Pipeline stages (pseudocode map) ---
# 1) Load each 10x folder -> Seurat object; merge; attach donor/treatment metadata
# 2) QC metrics + violin plots; filter cells by config thresholds
# 3) Normalize -> HVG -> scale -> PCA -> neighbors -> clusters -> UMAP
# 4) FindAllMarkers; export CSV + UMAP figures + optional RDS

loaded <- load_and_merge_samples(
  cfg$paths$data_dir,
  meta,
  qc = cfg$qc,
  subsample = cfg$subsample,
  seed = cfg$analysis$random_seed
)
obj <- loaded$obj
cells_before <- loaded$cells_before
cells_after <- ncol(obj)
message("Cells in merged object: ", cells_after, " (raw barcodes loaded: ", cells_before, ")")

obj <- add_qc_metrics(obj)
save_plot(
  plot_qc_violin(obj),
  file.path(run_dir, "plots", "qc_violin_post_filter.png"),
  width = cfg$outputs$plot_width,
  height = cfg$outputs$plot_height,
  dpi = cfg$outputs$dpi
)

obj <- run_standard_workflow(obj, cfg$analysis)

umap_plots <- plot_umap_panels(obj)
for (nm in names(umap_plots)) {
  save_plot(
    umap_plots[[nm]],
    file.path(run_dir, "plots", paste0("umap_", nm, ".", cfg$outputs$plot_format)),
    width = cfg$outputs$plot_width,
    height = cfg$outputs$plot_height,
    dpi = cfg$outputs$dpi
  )
}

markers <- find_cluster_markers(obj, cfg$markers)
export_marker_tables(markers, run_dir, cfg$markers$top_n_per_cluster)

if (isTRUE(cfg$outputs$save_rds)) {
  saveRDS(obj, file.path(run_dir, "seurat_processed.rds"))
}

yaml::write_yaml(cfg, file.path(run_dir, "config_used.yaml"))
write_run_summary(
  run_dir,
  cfg,
  meta,
  extra = list(
    cells_raw_barcodes = cells_before,
    cells_merged = cells_after,
    subsample_enabled = ifelse(!is.null(cfg$subsample) && isTRUE(cfg$subsample$enabled), "true", "false")
  )
)

message("Done. Outputs in: ", run_dir)
