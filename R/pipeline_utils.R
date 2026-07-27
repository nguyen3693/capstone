# Shared helpers for run_scrna_pipeline.R

require_config_packages <- function() {
  pkgs <- c("Seurat", "ggplot2", "dplyr", "yaml", "readr")
  missing <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing) > 0) {
    stop(
      "Install missing packages: ", paste(missing, collapse = ", "),
      "\nExample: install.packages(c('yaml', 'readr')); install.packages('Seurat')",
      call. = FALSE
    )
  }
  invisible(pkgs)
}

load_config <- function(path = "config.yaml") {
  if (!file.exists(path)) {
    stop("Config not found: ", path, call. = FALSE)
  }
  yaml::read_yaml(path)
}

make_run_dir <- function(output_dir, run_id = NULL) {
  if (is.null(run_id) || !nzchar(run_id)) {
    run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
  }
  run_dir <- file.path(output_dir, run_id)
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(run_dir, "plots"), showWarnings = FALSE)
  dir.create(file.path(run_dir, "tables"), showWarnings = FALSE)
  list(run_id = run_id, run_dir = run_dir)
}

read_sample_metadata <- function(metadata_file, include_samples = character()) {
  if (!file.exists(metadata_file)) {
    stop(
      "Metadata file not found: ", metadata_file,
      "\nCopy metadata/sample_metadata.csv.example to metadata/sample_metadata.csv and edit.",
      call. = FALSE
    )
  }
  meta <- readr::read_csv(metadata_file, show_col_types = FALSE)
  required <- c("sample_id", "folder_name", "donor", "treatment")
  if (!all(required %in% names(meta))) {
    stop("Metadata must contain columns: ", paste(required, collapse = ", "), call. = FALSE)
  }
  if (length(include_samples) > 0) {
    meta <- meta[meta$sample_id %in% include_samples, , drop = FALSE]
  }
  if (nrow(meta) == 0) {
    stop("No samples left after filtering metadata.", call. = FALSE)
  }
  meta
}

write_run_summary <- function(run_dir, cfg, meta, extra = list()) {
  summary_path <- file.path(run_dir, "run_summary.txt")
  lines <- c(
    paste("project:", cfg$project$name),
    paste("n_samples:", nrow(meta)),
    paste("samples:", paste(meta$sample_id, collapse = ", ")),
    paste("qc:", paste(names(cfg$qc), cfg$qc, sep = "=", collapse = "; ")),
    paste("cluster_resolution:", cfg$analysis$cluster_resolution),
    paste("seed:", cfg$analysis$random_seed)
  )
  if (length(extra) > 0) {
    lines <- c(lines, paste(names(extra), extra, sep = ": "))
  }
  writeLines(lines, summary_path)
  invisible(summary_path)
}

save_plot <- function(p, path, width = 8, height = 6, dpi = 150) {
  ggplot2::ggsave(path, plot = p, width = width, height = height, dpi = dpi)
  invisible(path)
}
