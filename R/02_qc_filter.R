# Stage 2: QC metrics, plots, and filtering

add_qc_metrics <- function(obj) {
  obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = "^MT-")
  obj
}

plot_qc_violin <- function(obj) {
  Seurat::VlnPlot(
    obj,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    group.by = "sample_id",
    ncol = 3,
    pt.size = 0
  )
}

downsample_by_sample <- function(obj, max_cells_per_sample, seed = 42) {
  set.seed(seed)
  keep_cells <- unlist(
    tapply(colnames(obj), obj$sample_id, function(cells) {
      if (length(cells) <= max_cells_per_sample) {
        cells
      } else {
        sample(cells, max_cells_per_sample)
      }
    }),
    use.names = FALSE
  )
  subset(obj, cells = keep_cells)
}

filter_cells <- function(obj, qc) {
  subset(
    obj,
    subset =
      nFeature_RNA >= qc$min_features &
      nFeature_RNA <= qc$max_features &
      nCount_RNA >= qc$min_counts &
      nCount_RNA <= qc$max_counts &
      percent.mt <= qc$max_percent_mt
  )
}
