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
                           progress_cb = NULL) {
  report <- function(detail, value, eta = NA_real_) {
    if (is.function(progress_cb)) {
      progress_cb(detail = detail, value = value, eta = eta)
    }
  }

  t0 <- proc.time()[["elapsed"]]
  stage_times <- numeric()

  dims <- as.integer(dims)
  max_pc <- ncol(obj[["pca"]])
  if (dims < 2) {
    stop("dims must be >= 2", call. = FALSE)
  }
  if (dims > max_pc) {
    dims <- max_pc
  }

  # Stage weights for ETA (rough relative cost)
  weights <- c(subsample = 0.05, neighbors = 0.25, clusters = 0.15, umap = 0.55)
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

  n_used <- ncol(obj)
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
    weights <- weights / sum(weights)
  }

  set.seed(as.integer(seed))

  report("Finding neighbors...", done_w, NA_real_)
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

lineage_counts <- function(obj) {
  if (!"t_lineage" %in% colnames(obj[[]])) {
    return(NULL)
  }
  as.data.frame(table(Lineage = obj$t_lineage), stringsAsFactors = FALSE)
}

plot_umap_lineage_split <- function(obj, group_by = "seurat_clusters", label = TRUE) {
  if (!"t_lineage" %in% colnames(obj[[]])) {
    obj <- assign_cd4_cd8(obj)
  }

  cd4_cells <- colnames(obj)[obj$t_lineage == "CD4"]
  cd8_cells <- colnames(obj)[obj$t_lineage == "CD8"]

  if (length(cd4_cells) < 10 || length(cd8_cells) < 10) {
    stop(
      sprintf(
        "Not enough cells after split (CD4=%d, CD8=%d). Lower min score or include more cells.",
        length(cd4_cells), length(cd8_cells)
      ),
      call. = FALSE
    )
  }

  p_cd4 <- plot_umap(
    subset(obj, cells = cd4_cells),
    group_by = group_by,
    title = sprintf("CD4+ inferred (n=%s)", format(length(cd4_cells), big.mark = ",")),
    label = label
  )
  p_cd8 <- plot_umap(
    subset(obj, cells = cd8_cells),
    group_by = group_by,
    title = sprintf("CD8+ inferred (n=%s)", format(length(cd8_cells), big.mark = ",")),
    label = label
  )
  list(cd4 = p_cd4, cd8 = p_cd8, n_cd4 = length(cd4_cells), n_cd8 = length(cd8_cells))
}
