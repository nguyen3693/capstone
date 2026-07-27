# Stage 4: Cluster markers and export tables

find_cluster_markers <- function(obj, markers_cfg) {
  Seurat::FindAllMarkers(
    obj,
    only.pos = markers_cfg$only_pos,
    min.pct = markers_cfg$min_pct,
    logfc.threshold = markers_cfg$logfc_threshold,
    test.use = markers_cfg$test_use
  )
}

top_markers_per_cluster <- function(marker_df, top_n = 25) {
  marker_df |>
    dplyr::group_by(cluster) |>
    dplyr::slice_max(order_by = avg_log2FC, n = top_n, with_ties = FALSE) |>
    dplyr::ungroup()
}

export_marker_tables <- function(marker_df, run_dir, top_n = 25) {
  all_path <- file.path(run_dir, "tables", "cluster_markers_all.csv")
  top_path <- file.path(run_dir, "tables", "cluster_markers_top.csv")
  readr::write_csv(marker_df, all_path)
  readr::write_csv(top_markers_per_cluster(marker_df, top_n), top_path)
  invisible(list(all = all_path, top = top_path))
}
