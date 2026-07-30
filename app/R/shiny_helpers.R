# Helpers for the Shiny parameter playground

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || (is.character(x) && !nzchar(x))) y else x

default_rds_candidates <- function(project_root = "..") {
  results_dir <- file.path(project_root, "results")
  if (!dir.exists(results_dir)) {
    return(character())
  }
  rds <- list.files(
    results_dir,
    pattern = "^seurat_processed\\.rds$",
    recursive = TRUE,
    full.names = TRUE
  )
  # Prefer newest full-run style folders first
  rds[order(file.info(rds)$mtime, decreasing = TRUE)]
}

load_seurat_rds <- function(path) {
  if (!file.exists(path)) {
    stop("RDS not found: ", path, call. = FALSE)
  }
  obj <- readRDS(path)
  if (!inherits(obj, "Seurat")) {
    stop("File is not a Seurat object: ", path, call. = FALSE)
  }
  if (!"pca" %in% names(obj@reductions)) {
    stop("Object has no PCA reduction. Re-run the pipeline first.", call. = FALSE)
  }
  obj
}

# Randomly keep up to n_cells. If sample_id exists, sample roughly evenly across samples.
subsample_cells <- function(obj, n_cells, seed = 42) {
  n_cells <- as.integer(n_cells)
  n_total <- ncol(obj)
  if (is.na(n_cells) || n_cells < 100) {
    stop("n_cells must be at least 100", call. = FALSE)
  }
  if (n_cells >= n_total) {
    return(list(obj = obj, n_cells = n_total, subsampled = FALSE))
  }

  set.seed(as.integer(seed))
  if ("sample_id" %in% colnames(obj[[]])) {
    keep <- unlist(
      tapply(colnames(obj), obj$sample_id, function(cells) {
        n_per <- max(1L, floor(n_cells * length(cells) / n_total))
        if (length(cells) <= n_per) cells else sample(cells, n_per)
      }),
      use.names = FALSE
    )
    # Top up or trim to exact n_cells
    if (length(keep) > n_cells) {
      keep <- sample(keep, n_cells)
    } else if (length(keep) < n_cells) {
      leftover <- setdiff(colnames(obj), keep)
      keep <- c(keep, sample(leftover, min(length(leftover), n_cells - length(keep))))
    }
  } else {
    keep <- sample(colnames(obj), n_cells)
  }

  list(
    obj = subset(obj, cells = keep),
    n_cells = length(keep),
    subsampled = TRUE
  )
}

fmt_secs <- function(secs) {
  secs <- max(0, as.numeric(secs))
  if (secs < 60) {
    sprintf("~%.0fs left", secs)
  } else {
    sprintf("~%.1f min left", secs / 60)
  }
}

# Recompute neighbors, clusters, and UMAP from an existing processed object.
# progress_cb(detail, value 0-1, eta_seconds) is optional (used by Shiny).
recluster_umap <- function(obj,
                           resolution = 0.5,
                           dims = 30,
                           seed = 42,
                           n_cells = NULL,
                           lineage = "All",
                           progress_cb = NULL) {
  report <- function(detail, value, eta = NA_real_) {
    if (is.function(progress_cb)) {
      progress_cb(detail = detail, value = value, eta = eta)
    }
  }

  t0 <- proc.time()[["elapsed"]]
  stage_times <- numeric()

  dims <- as.integer(dims)
  if (dims < 2) {
    stop("dims must be >= 2", call. = FALSE)
  }

  # Stage weights for ETA (rough relative cost)
  weights <- c(
    lineage = 0.05,
    subsample = 0.05,
    preprocess = 0.25,
    pca = 0.15,
    neighbors = 0.15,
    clusters = 0.10,
    umap = 0.25
  )
  done_w <- 0

  update_eta <- function(just_finished_stage) {
    elapsed <- proc.time()[["elapsed"]] - t0
    stage_times[[just_finished_stage]] <<- elapsed - sum(stage_times)
    # Estimate remaining from remaining weight × observed pace
    remaining_w <- sum(weights[names(weights) != just_finished_stage &
                                 !(names(weights) %in% names(stage_times))])
    # Pace: elapsed so far / weight completed
    completed_w <- sum(weights[names(stage_times)])
    if (completed_w <= 0) {
      return(NA_real_)
    }
    pace <- elapsed / completed_w
    remaining_w * pace
  }

  lineage <- match.arg(lineage, c("All", "CD4", "CD8"))
  if (lineage != "All") {
    report(paste("Inferring and selecting", lineage, "cells..."), done_w, NA_real_)
    obj <- assign_cd4_cd8(obj)
    lineage_cells <- colnames(obj)[as.character(obj$t_lineage) == lineage]
    if (length(lineage_cells) < 100) {
      stop(
        sprintf("Only %d inferred %s cells were found; at least 100 are required.", length(lineage_cells), lineage),
        call. = FALSE
      )
    }
    obj <- subset(obj, cells = lineage_cells)
    done_w <- done_w + weights[["lineage"]]
    report(
      sprintf("Selected %s inferred %s cells", format(ncol(obj), big.mark = ","), lineage),
      done_w,
      update_eta("lineage")
    )
  } else {
    weights[["lineage"]] <- 0
  }

  n_available <- ncol(obj)
  n_used <- n_available
  subsampled <- FALSE
  if (!is.null(n_cells) && as.integer(n_cells) < ncol(obj)) {
    report("Subsampling cells...", done_w, NA_real_)
    sub <- subsample_cells(obj, n_cells, seed)
    obj <- sub$obj
    n_used <- sub$n_cells
    subsampled <- sub$subsampled
    done_w <- done_w + weights[["subsample"]]
    report(
      sprintf("Subsampled to %s cells", format(n_used, big.mark = ",")),
      done_w,
      update_eta("subsample")
    )
  } else {
    # Skip subsample stage weight
    weights[["subsample"]] <- 0
  }

  set.seed(as.integer(seed))

  report("Finding lineage-specific variable genes and scaling...", done_w, NA_real_)
  obj <- Seurat::FindVariableFeatures(
    obj,
    selection.method = "vst",
    nfeatures = min(1500L, nrow(obj)),
    verbose = FALSE
  )
  obj <- Seurat::ScaleData(
    obj,
    features = Seurat::VariableFeatures(obj),
    verbose = FALSE
  )
  done_w <- done_w + weights[["preprocess"]]
  report("Preprocessing done. Recomputing PCA...", done_w, update_eta("preprocess"))

  max_dims <- min(dims, ncol(obj) - 1L, length(Seurat::VariableFeatures(obj)))
  if (max_dims < 2) {
    stop("Not enough cells or variable genes to compute PCA.", call. = FALSE)
  }
  dims <- max_dims
  obj <- Seurat::RunPCA(
    obj,
    features = Seurat::VariableFeatures(obj),
    npcs = dims,
    verbose = FALSE
  )
  done_w <- done_w + weights[["pca"]]
  report("PCA done. Finding neighbors...", done_w, update_eta("pca"))

  obj <- Seurat::FindNeighbors(obj, dims = 1:dims, verbose = FALSE)
  done_w <- done_w + weights[["neighbors"]]
  report("Neighbors done. Clustering...", done_w, update_eta("neighbors"))

  obj <- Seurat::FindClusters(obj, resolution = as.numeric(resolution), verbose = FALSE)
  done_w <- done_w + weights[["clusters"]]
  report("Clusters done. Running UMAP...", done_w, update_eta("clusters"))

  obj <- Seurat::RunUMAP(obj, dims = 1:dims, verbose = FALSE)
  done_w <- 1
  report("UMAP done.", done_w, 0)

  list(
    obj = obj,
    n_clusters = length(unique(Seurat::Idents(obj))),
    n_cells = n_used,
    n_lineage_available = n_available,
    lineage = lineage,
    subsampled = subsampled,
    resolution = as.numeric(resolution),
    dims = dims,
    seed = as.integer(seed),
    elapsed_sec = proc.time()[["elapsed"]] - t0
  )
}

plot_umap <- function(obj, group_by = "seurat_clusters", title = NULL, label = TRUE) {
  meta_cols <- colnames(obj[[]])
  if (!group_by %in% meta_cols) {
    stop("Metadata column not found: ", group_by, call. = FALSE)
  }

  p <- Seurat::DimPlot(
    obj,
    reduction = "umap",
    group.by = group_by,
    label = isTRUE(label) && group_by == "seurat_clusters",
    repel = TRUE
  ) + ggplot2::theme_minimal()

  if (!is.null(title) && nzchar(title)) {
    p <- p + ggplot2::ggtitle(title)
  }
  p
}

cluster_treatment_table <- function(obj) {
  if (!all(c("seurat_clusters", "treatment") %in% colnames(obj[[]]))) {
    return(NULL)
  }
  tab <- table(Cluster = obj$seurat_clusters, Treatment = obj$treatment)
  as.data.frame.matrix(tab)
}

param_label <- function(resolution, dims, seed, n_cells = NULL) {
  base <- sprintf("res=%.2f | dims=%d | seed=%d", resolution, dims, seed)
  if (!is.null(n_cells)) {
    base <- paste0(base, " | cells=", format(as.integer(n_cells), big.mark = ","))
  }
  base
}

# Infer CD4 vs CD8 from transcript levels (not FACS). Uses normalized data layer.
assign_cd4_cd8 <- function(obj, cd4_gene = "CD4", cd8_genes = c("CD8A", "CD8B"), min_score = 0.25) {
  feats <- rownames(obj)
  if (!cd4_gene %in% feats) {
    stop("CD4 gene not found in object.", call. = FALSE)
  }
  cd8_genes <- cd8_genes[cd8_genes %in% feats]
  if (length(cd8_genes) == 0) {
    stop("No CD8 genes (CD8A/CD8B) found in object.", call. = FALSE)
  }

  mat <- tryCatch(
    Seurat::GetAssayData(obj, layer = "data"),
    error = function(e) Seurat::GetAssayData(obj, slot = "data")
  )

  cd4 <- as.numeric(mat[cd4_gene, ])
  cd8 <- as.numeric(Matrix::colMeans(mat[cd8_genes, , drop = FALSE]))

  lineage <- rep("Ambiguous", length(cd4))
  lineage[cd4 > cd8 & cd4 >= min_score] <- "CD4"
  lineage[cd8 > cd4 & cd8 >= min_score] <- "CD8"

  obj$cd4_score <- cd4
  obj$cd8_score <- cd8
  obj$t_lineage <- factor(lineage, levels = c("CD4", "CD8", "Ambiguous"))
  obj
}

# Find positive marker genes for every generated cluster and create readable labels.
compute_cluster_markers <- function(obj,
                                    top_n = 5,
                                    min_pct = 0.25,
                                    logfc_threshold = 0.25) {
  if (!"seurat_clusters" %in% colnames(obj[[]])) {
    stop("No generated clusters found. Compute Run A or Run B first.", call. = FALSE)
  }

  top_n <- max(1L, as.integer(top_n))
  obj <- Seurat::SetIdent(obj, value = "seurat_clusters")
  markers <- Seurat::FindAllMarkers(
    obj,
    only.pos = TRUE,
    min.pct = as.numeric(min_pct),
    logfc.threshold = as.numeric(logfc_threshold),
    test.use = "wilcox",
    verbose = FALSE
  )

  if (nrow(markers) == 0) {
    stop("No marker genes passed the current thresholds.", call. = FALSE)
  }

  fc_col <- if ("avg_log2FC" %in% colnames(markers)) "avg_log2FC" else "avg_logFC"
  markers <- markers[order(markers$cluster, -markers[[fc_col]]), , drop = FALSE]
  by_cluster <- split(markers, markers$cluster)
  top_markers <- do.call(
    rbind,
    lapply(by_cluster, function(x) utils::head(x, top_n))
  )
  rownames(top_markers) <- NULL

  summaries <- lapply(by_cluster, function(x) {
    genes <- utils::head(x$gene, top_n)
    data.frame(
      cluster = as.character(x$cluster[[1]]),
      top_genes = paste(genes, collapse = ", "),
      stringsAsFactors = FALSE
    )
  })
  cluster_summary <- do.call(rbind, summaries)
  rownames(cluster_summary) <- NULL

  label_map <- stats::setNames(
    paste0("Cluster ", cluster_summary$cluster, ": ", cluster_summary$top_genes),
    cluster_summary$cluster
  )
  annotations <- unname(label_map[as.character(obj$seurat_clusters)])
  obj$cluster_annotation <- factor(annotations, levels = unname(label_map))

  list(
    obj = obj,
    all_markers = markers,
    top_markers = top_markers,
    cluster_summary = cluster_summary,
    fc_column = fc_col
  )
}

plot_annotated_umap <- function(marker_result, title = "Clusters labeled by top marker genes") {
  Seurat::DimPlot(
    marker_result$obj,
    reduction = "umap",
    group.by = "cluster_annotation",
    label = TRUE,
    repel = TRUE
  ) +
    ggplot2::ggtitle(title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")
}

# Plot a fixed T-cell gene panel across all generated clusters.
create_expression_panel_violin <- function(
    obj,
    genes = c(
      "CD4", "CD8A", "GZMA", "GZMB", "IL7R", "IL2RA",
      "TCF7", "FOXP3", "TIGIT", "CCR7", "BCL2", "MCL1"
    )) {
  present_genes <- genes[genes %in% rownames(obj)]
  missing_genes <- setdiff(genes, present_genes)
  if (length(present_genes) == 0) {
    stop("None of the requested expression-panel genes are present.", call. = FALSE)
  }

  p <- Seurat::VlnPlot(
    obj,
    features = present_genes,
    group.by = "seurat_clusters",
    pt.size = 0,
    ncol = 1,
    combine = TRUE
  ) &
    ggplot2::theme_minimal() &
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  p <- p + patchwork::plot_annotation(
    title = "T-cell expression panel across generated clusters"
  )

  list(
    plot = p,
    genes = present_genes,
    missing_genes = missing_genes
  )
}
