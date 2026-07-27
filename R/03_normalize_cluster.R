# Stage 3: Normalize, HVGs, PCA, neighbors, clusters, UMAP

run_standard_workflow <- function(obj, analysis) {
  set.seed(analysis$random_seed)
  regress_vars <- analysis$scale_regress_vars
  if (is.null(regress_vars)) {
    regress_vars <- character()
  }

  obj <- Seurat::NormalizeData(obj, normalization.method = analysis$normalization_method)
  obj <- Seurat::FindVariableFeatures(obj, selection.method = "vst", nfeatures = analysis$n_variable_features)
  if (length(regress_vars) == 0) {
    obj <- Seurat::ScaleData(obj)
  } else {
    obj <- Seurat::ScaleData(obj, vars.to.regress = regress_vars)
  }
  obj <- Seurat::RunPCA(obj, npcs = analysis$pca_dims, verbose = FALSE)
  obj <- Seurat::FindNeighbors(obj, dims = 1:analysis$pca_dims)
  obj <- Seurat::FindClusters(obj, resolution = analysis$cluster_resolution)
  obj <- Seurat::RunUMAP(obj, dims = 1:analysis$umap_dims, verbose = FALSE)
  obj
}

plot_umap_panels <- function(obj) {
  p1 <- Seurat::DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
    ggplot2::ggtitle("Clusters")
  p2 <- Seurat::DimPlot(obj, reduction = "umap", group.by = "treatment") +
    ggplot2::ggtitle("Treatment")
  p3 <- Seurat::DimPlot(obj, reduction = "umap", group.by = "sample_id") +
    ggplot2::ggtitle("Sample")
  list(clusters = p1, treatment = p2, sample = p3)
}
