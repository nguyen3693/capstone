# Stage 1: Per-sample QC/downsample, then merge to one Seurat object (single counts layer)

load_and_merge_samples <- function(data_dir, meta, qc = NULL, subsample = NULL, seed = 42) {
  matrices <- list()
  cells_before <- 0L

  for (i in seq_len(nrow(meta))) {
    sample_id <- meta$sample_id[i]
    folder <- file.path(data_dir, meta$folder_name[i])
    if (!dir.exists(folder)) {
      stop("Sample folder not found: ", folder, call. = FALSE)
    }

    counts <- Seurat::Read10X(data.dir = folder)
    cells_before <- cells_before + ncol(counts)

    if (!is.null(qc)) {
      obj <- Seurat::CreateSeuratObject(
        counts = counts,
        project = sample_id,
        min.cells = 3,
        min.features = 200
      )
      obj <- add_qc_metrics(obj)
      obj$sample_id <- sample_id
      obj <- filter_cells(obj, qc)
      if (!is.null(subsample) && isTRUE(subsample$enabled)) {
        obj <- downsample_by_sample(obj, subsample$max_cells_per_sample, seed)
      }
      counts <- Seurat::GetAssayData(obj, layer = "counts")
    }

    colnames(counts) <- paste0(sample_id, "_", colnames(counts))
    matrices[[sample_id]] <- counts
  }

  common_genes <- Reduce(intersect, lapply(matrices, rownames))
  if (length(common_genes) == 0) {
    stop("No shared genes across samples.", call. = FALSE)
  }

  matrices <- lapply(matrices, function(m) m[common_genes, , drop = FALSE])
  merged_counts <- do.call(cbind, matrices)

  obj <- Seurat::CreateSeuratObject(
    counts = merged_counts,
    project = "GSE288434",
    min.cells = 3,
    min.features = 200
  )

  sample_map <- sub("_.*", "", colnames(obj))
  obj$sample_id <- sample_map
  obj$donor <- meta$donor[match(sample_map, meta$sample_id)]
  obj$treatment <- meta$treatment[match(sample_map, meta$sample_id)]

  list(obj = obj, cells_before = cells_before)
}
