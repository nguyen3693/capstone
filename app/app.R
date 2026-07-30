# Capstone Shiny dashboard: compare Seurat clustering/UMAP parameter settings
# Launch from project root:
#   R -e "shiny::runApp('app')"
# Or from the app folder:
#   shiny::runApp()

library(shiny)
library(Seurat)
library(ggplot2)
library(patchwork)

source("R/shiny_helpers.R", local = TRUE)

project_root <- normalizePath("..")
rds_choices <- default_rds_candidates(project_root)
default_rds <- if (length(rds_choices)) rds_choices[[1]] else ""

ui <- fluidPage(
  titlePanel("scRNA-seq parameter playground (GSE288434)"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Data"),
      selectInput(
        "rds_path",
        "Processed Seurat RDS",
        choices = rds_choices,
        selected = default_rds
      ),
      actionButton("load_btn", "Load object", class = "btn-primary"),
      hr(),
      h4("Cells to include"),
      sliderInput(
        "n_cells",
        "Cells per lineage (random subsample; balanced by sample when possible)",
        min = 500,
        max = 3000,
        value = 2000,
        step = 100
      ),
      helpText("Run A uses inferred CD4 cells; Run B uses inferred CD8 cells. The same requested count keeps comparisons balanced."),
      hr(),
      h4("Shared display"),
      selectInput(
        "color_by",
        "Color UMAP by",
        choices = c(
          "Clusters" = "seurat_clusters",
          "Treatment" = "treatment",
          "Sample" = "sample_id",
          "Donor" = "donor"
        ),
        selected = "seurat_clusters"
      ),
      checkboxInput("label_clusters", "Label clusters (when coloring by cluster)", TRUE),
      hr(),
      h4("Run A — inferred CD4"),
      sliderInput("res_a", "Cluster resolution", min = 0.1, max = 1.5, value = 0.5, step = 0.1),
      sliderInput("dims_a", "PCA / UMAP dims", min = 5, max = 40, value = 30, step = 1),
      numericInput("seed_a", "Random seed", value = 42, min = 1, step = 1),
      actionButton("run_a", "Compute Run A", class = "btn-success"),
      hr(),
      h4("Run B — inferred CD8"),
      sliderInput("res_b", "Cluster resolution", min = 0.1, max = 1.5, value = 0.8, step = 0.1),
      sliderInput("dims_b", "PCA / UMAP dims", min = 5, max = 40, value = 30, step = 1),
      numericInput("seed_b", "Random seed", value = 42, min = 1, step = 1),
      actionButton("run_b", "Compute Run B", class = "btn-success"),
      hr(),
      h4("Export UMAP bundles"),
      actionButton("save_run_a", "Save Run A outputs", class = "btn-primary"),
      downloadButton("download_run_a", "Download Run A ZIP"),
      actionButton("save_run_b", "Save Run B outputs", class = "btn-primary"),
      downloadButton("download_run_b", "Download Run B ZIP"),
      helpText(
        "Save writes directly to parameter-playground/ in this project. ",
        "ZIP downloads the same folder through your browser."
      ),
      hr(),
      h4("Cluster annotations"),
      selectInput(
        "marker_run",
        "Annotate which run?",
        choices = c("Run A (CD4)" = "A", "Run B (CD8)" = "B"),
        selected = "A"
      ),
      sliderInput(
        "top_marker_n",
        "Top genes shown per cluster",
        min = 1,
        max = 10,
        value = 5,
        step = 1
      ),
      actionButton("compute_markers", "Compute cluster markers", class = "btn-info"),
      downloadButton("download_markers", "Download marker CSV"),
      actionButton("compute_violin", "Generate CD4 vs CD8 expression panels", class = "btn-warning"),
      downloadButton("download_violin", "Download violin PNG"),
      helpText(
        "Marker calculation is separate because it can take several minutes. ",
        "The expression comparison uses a fixed T-cell gene panel after both runs are computed. ",
        "Marker genes describe a cluster but do not automatically prove its cell type."
      ),
      hr(),
      helpText(
        "Loads a saved pipeline object, then re-runs neighbors → clusters → UMAP ",
        "with your settings. Progress bar shows stage + estimated time left."
      )
    ),
    mainPanel(
      width = 9,
      textOutput("status"),
      textOutput("export_status"),
      tags$hr(),
      fluidRow(
        column(6, h3("Run A — inferred CD4"), plotOutput("plot_a", height = "480px"), verbatimTextOutput("summary_a")),
        column(6, h3("Run B — inferred CD8"), plotOutput("plot_b", height = "480px"), verbatimTextOutput("summary_b"))
      ),
      tags$hr(),
      h4("Cluster × treatment counts"),
      fluidRow(
        column(6, tableOutput("table_a")),
        column(6, tableOutput("table_b"))
      ),
      tags$hr(),
      h3("Cluster marker annotations"),
      textOutput("marker_status"),
      plotOutput("marker_umap", height = "620px"),
      h4("Top associated genes"),
      tableOutput("marker_summary"),
      h4("Marker details"),
      tableOutput("marker_table"),
      tags$hr(),
      h3("CD4 vs CD8 T-cell gene-expression panels"),
      textOutput("violin_status"),
      fluidRow(
        column(6, h4("Run A — inferred CD4"), plotOutput("violin_plot_a", height = "2200px")),
        column(6, h4("Run B — inferred CD8"), plotOutput("violin_plot_b", height = "2200px"))
      )
    )
  )
)

server <- function(input, output, session) {
  base_obj <- reactiveVal(NULL)
  run_a <- reactiveVal(NULL)
  run_b <- reactiveVal(NULL)
  marker_result <- reactiveVal(NULL)
  violin_result <- reactiveVal(NULL)
  status_msg <- reactiveVal("Load a Seurat RDS to begin.")
  export_msg <- reactiveVal("No output bundle saved yet.")

  output$status <- renderText(status_msg())
  output$export_status <- renderText(export_msg())

  get_run_result <- function(run_label) {
    if (identical(run_label, "A")) run_a() else run_b()
  }

  save_bundle <- function(run_label) {
    result <- get_run_result(run_label)
    if (is.null(result)) {
      stop("Compute Run ", run_label, " before exporting.", call. = FALSE)
    }
    save_run_output_bundle(
      run_result = result,
      run_label = run_label,
      output_root = file.path(project_root, "parameter-playground"),
      source_rds = input$rds_path
    )
  }

  run_with_progress <- function(run_name, obj, resolution, dims, seed, n_cells, lineage) {
    withProgress(
      message = paste(run_name, "in progress"),
      value = 0,
      {
        setProgress(0, detail = "Starting...")
        out <- recluster_umap(
          obj,
          resolution = resolution,
          dims = dims,
          seed = seed,
          n_cells = n_cells,
          lineage = lineage,
          progress_cb = function(detail, value, eta) {
            detail_txt <- detail
            if (!is.na(eta)) {
              detail_txt <- paste0(detail, "  |  ", fmt_secs(eta))
            }
            setProgress(value = max(0, min(1, value)), detail = detail_txt)
          }
        )
        setProgress(1, detail = "Done")
        out
      }
    )
  }

  observeEvent(input$load_btn, {
    path <- input$rds_path
    if (is.null(path) || !nzchar(path)) {
      status_msg("No RDS path selected. Run the pipeline first to create results/*/seurat_processed.rds.")
      return()
    }
    status_msg(paste("Loading", basename(dirname(path)), "..."))
    tryCatch({
      withProgress(message = "Loading Seurat object", value = 0.3, {
        setProgress(0.3, detail = "Reading RDS (may take a bit)...")
        obj <- load_seurat_rds(path)
        obj <- assign_cd4_cd8(obj)
        setProgress(1, detail = "Loaded")
      })
      base_obj(obj)
      run_a(NULL)
      run_b(NULL)
      marker_result(NULL)
      violin_result(NULL)
      export_msg("No output bundle saved for the newly loaded object.")

      lineage_counts <- table(obj$t_lineage)
      n_cd4 <- as.integer(lineage_counts[["CD4"]])
      n_cd8 <- as.integer(lineage_counts[["CD8"]])
      max_lineage_cells <- min(n_cd4, n_cd8)
      updateSliderInput(
        session,
        "n_cells",
        min = min(500, max_lineage_cells),
        max = max_lineage_cells,
        value = min(2000, max_lineage_cells),
        step = max(50, round(max_lineage_cells / 25))
      )

      status_msg(sprintf(
        "Loaded %s cells (%s inferred CD4; %s inferred CD8). Choose cells per lineage, then compute both runs.",
        format(ncol(obj), big.mark = ","),
        format(n_cd4, big.mark = ","),
        format(n_cd8, big.mark = ",")
      ))
    }, error = function(e) {
      status_msg(paste("Load failed:", conditionMessage(e)))
    })
  })

  observeEvent(input$run_a, {
    obj <- base_obj()
    if (is.null(obj)) {
      status_msg("Load an object first.")
      return()
    }
    status_msg("Computing Run A...")
    tryCatch({
      out <- run_with_progress(
        "Run A",
        obj,
        resolution = input$res_a,
        dims = input$dims_a,
        seed = input$seed_a,
        n_cells = input$n_cells,
        lineage = "CD4"
      )
      run_a(out)
      marker_result(NULL)
      violin_result(NULL)
      export_msg("Run A changed; save again to export its current parameters.")
      status_msg(sprintf(
        "Run A done in %.1fs | %s | clusters = %d",
        out$elapsed_sec,
        param_label(out$resolution, out$dims, out$seed, out$n_cells),
        out$n_clusters
      ))
    }, error = function(e) {
      status_msg(paste("Run A failed:", conditionMessage(e)))
    })
  })

  observeEvent(input$run_b, {
    obj <- base_obj()
    if (is.null(obj)) {
      status_msg("Load an object first.")
      return()
    }
    status_msg("Computing Run B...")
    tryCatch({
      out <- run_with_progress(
        "Run B",
        obj,
        resolution = input$res_b,
        dims = input$dims_b,
        seed = input$seed_b,
        n_cells = input$n_cells,
        lineage = "CD8"
      )
      run_b(out)
      marker_result(NULL)
      violin_result(NULL)
      export_msg("Run B changed; save again to export its current parameters.")
      status_msg(sprintf(
        "Run B done in %.1fs | %s | clusters = %d",
        out$elapsed_sec,
        param_label(out$resolution, out$dims, out$seed, out$n_cells),
        out$n_clusters
      ))
    }, error = function(e) {
      status_msg(paste("Run B failed:", conditionMessage(e)))
    })
  })

  observeEvent(input$save_run_a, {
    tryCatch({
      bundle <- withProgress(message = "Saving Run A output bundle", value = 0.2, {
        setProgress(0.2, detail = "Rendering four UMAP PNGs...")
        out <- save_bundle("A")
        setProgress(1, detail = "Saved")
        out
      })
      export_msg(paste("Saved Run A outputs to:", bundle$directory))
    }, error = function(e) {
      export_msg(paste("Run A export failed:", conditionMessage(e)))
    })
  })

  observeEvent(input$save_run_b, {
    tryCatch({
      bundle <- withProgress(message = "Saving Run B output bundle", value = 0.2, {
        setProgress(0.2, detail = "Rendering four UMAP PNGs...")
        out <- save_bundle("B")
        setProgress(1, detail = "Saved")
        out
      })
      export_msg(paste("Saved Run B outputs to:", bundle$directory))
    }, error = function(e) {
      export_msg(paste("Run B export failed:", conditionMessage(e)))
    })
  })

  output$download_run_a <- downloadHandler(
    filename = function() {
      paste0("Run-A-", format(Sys.time(), "%Y%m%d-%H%M%S"), ".zip")
    },
    content = function(file) {
      bundle <- save_bundle("A")
      zip_run_output_bundle(bundle, file)
      export_msg(paste("Saved and downloaded Run A bundle:", bundle$directory))
    },
    contentType = "application/zip"
  )

  output$download_run_b <- downloadHandler(
    filename = function() {
      paste0("Run-B-", format(Sys.time(), "%Y%m%d-%H%M%S"), ".zip")
    },
    content = function(file) {
      bundle <- save_bundle("B")
      zip_run_output_bundle(bundle, file)
      export_msg(paste("Saved and downloaded Run B bundle:", bundle$directory))
    },
    contentType = "application/zip"
  )

  observeEvent(input$compute_markers, {
    out <- if (identical(input$marker_run, "B")) run_b() else run_a()
    if (is.null(out)) {
      status_msg(paste("Compute Run", input$marker_run, "before calculating markers."))
      return()
    }

    status_msg(paste("Computing cluster markers for Run", input$marker_run, "..."))
    tryCatch({
      result <- withProgress(
        message = paste("Finding markers for Run", input$marker_run),
        value = 0,
        {
          setProgress(0.1, detail = "Comparing each cluster with all other cells...")
          markers <- compute_cluster_markers(
            out$obj,
            top_n = input$top_marker_n
          )
          setProgress(0.95, detail = "Building labels and tables...")
          markers$source_run <- input$marker_run
          markers$source_lineage <- out$lineage
          markers$top_n <- as.integer(input$top_marker_n)
          setProgress(1, detail = "Done")
          markers
        }
      )
      marker_result(result)
      violin_result(NULL)
      status_msg(sprintf(
        "Marker annotations ready for Run %s (%s): %d clusters, top %d genes shown.",
        result$source_run,
        result$source_lineage,
        nrow(result$cluster_summary),
        result$top_n
      ))
    }, error = function(e) {
      marker_result(NULL)
      violin_result(NULL)
      status_msg(paste("Marker calculation failed:", conditionMessage(e)))
    })
  })

  observeEvent(input$compute_violin, {
    out_a <- run_a()
    out_b <- run_b()
    if (is.null(out_a) || is.null(out_b)) {
      status_msg("Compute both Run A (CD4) and Run B (CD8) before generating the expression comparison.")
      return()
    }

    tryCatch({
      panels <- withProgress(
        message = "Building CD4 vs CD8 expression panels",
        value = 0,
        {
          setProgress(0.1, detail = "Building inferred CD4 panel...")
          panel_a <- create_expression_panel_violin(out_a$obj)
          setProgress(0.55, detail = "Building inferred CD8 panel...")
          panel_b <- create_expression_panel_violin(out_b$obj)
          setProgress(1, detail = "Done")
          list(a = panel_a, b = panel_b)
        }
      )
      violin_result(panels)
      missing <- union(panels$a$missing_genes, panels$b$missing_genes)
      message <- paste0(
        "CD4 vs CD8 expression panels ready. Genes shown: ",
        paste(panels$a$genes, collapse = ", ")
      )
      if (length(missing) > 0) {
        message <- paste0(message, ". Missing: ", paste(missing, collapse = ", "))
      }
      status_msg(message)
    }, error = function(e) {
      violin_result(NULL)
      status_msg(paste("Violin plot failed:", conditionMessage(e)))
    })
  })

  output$plot_a <- renderPlot({
    out <- run_a()
    req(out)
    obj <- out$obj
    plot_umap(
      obj,
      group_by = input$color_by,
      title = paste("A — inferred CD4:", param_label(out$resolution, out$dims, out$seed, out$n_cells)),
      label = input$label_clusters
    )
  })

  output$plot_b <- renderPlot({
    out <- run_b()
    req(out)
    obj <- out$obj
    plot_umap(
      obj,
      group_by = input$color_by,
      title = paste("B — inferred CD8:", param_label(out$resolution, out$dims, out$seed, out$n_cells)),
      label = input$label_clusters
    )
  })

  output$summary_a <- renderText({
    out <- run_a()
    if (is.null(out)) {
      return("Not computed yet.")
    }
    paste0(
      "Population: inferred ", out$lineage, "\n",
      "Clusters: ", out$n_clusters, "\n",
      "Cells used: ", format(out$n_cells, big.mark = ","),
      if (isTRUE(out$subsampled)) " (subsampled)" else " (all loaded)", "\n",
      param_label(out$resolution, out$dims, out$seed), "\n",
      sprintf("Elapsed: %.1fs", out$elapsed_sec)
    )
  })

  output$summary_b <- renderText({
    out <- run_b()
    if (is.null(out)) {
      return("Not computed yet.")
    }
    paste0(
      "Population: inferred ", out$lineage, "\n",
      "Clusters: ", out$n_clusters, "\n",
      "Cells used: ", format(out$n_cells, big.mark = ","),
      if (isTRUE(out$subsampled)) " (subsampled)" else " (all loaded)", "\n",
      param_label(out$resolution, out$dims, out$seed), "\n",
      sprintf("Elapsed: %.1fs", out$elapsed_sec)
    )
  })

  output$table_a <- renderTable({
    out <- run_a()
    req(out)
    cluster_treatment_table(out$obj)
  }, rownames = TRUE)

  output$table_b <- renderTable({
    out <- run_b()
    req(out)
    cluster_treatment_table(out$obj)
  }, rownames = TRUE)

  output$marker_status <- renderText({
    result <- marker_result()
    if (is.null(result)) {
      return("Choose Run A or B and click “Compute cluster markers.”")
    }
    sprintf(
      "Run %s annotations use the top %d positive marker genes per cluster. These are descriptive gene signatures, not validated cell-type names.",
      paste(result$source_run, paste0("(", result$source_lineage, ")")),
      result$top_n
    )
  })

  output$marker_umap <- renderPlot({
    result <- marker_result()
    req(result)
    plot_annotated_umap(
      result,
      title = paste("Run", result$source_run, "— inferred", result$source_lineage, "clusters labeled by top marker genes")
    )
  })

  output$marker_summary <- renderTable({
    result <- marker_result()
    req(result)
    result$cluster_summary
  }, striped = TRUE, hover = TRUE, spacing = "s")

  output$marker_table <- renderTable({
    result <- marker_result()
    req(result)
    markers <- result$top_markers
    preferred <- c(
      "cluster", "gene", result$fc_column,
      "pct.1", "pct.2", "p_val_adj"
    )
    markers[, preferred[preferred %in% colnames(markers)], drop = FALSE]
  }, striped = TRUE, hover = TRUE, spacing = "xs", digits = 3)

  output$download_markers <- downloadHandler(
    filename = function() {
      result <- marker_result()
      run_name <- if (is.null(result)) input$marker_run else result$source_run
      lineage <- if (is.null(result)) "" else paste0("_", tolower(result$source_lineage))
      paste0("run_", run_name, lineage, "_cluster_markers.csv")
    },
    content = function(file) {
      result <- marker_result()
      req(result)
      utils::write.csv(result$all_markers, file, row.names = FALSE)
    }
  )

  output$violin_status <- renderText({
    panels <- violin_result()
    if (is.null(panels)) {
      return("Compute both runs, then click “Generate CD4 vs CD8 expression panels.”")
    }
    text <- paste0("Genes shown in both panels: ", paste(panels$a$genes, collapse = ", "))
    missing <- union(panels$a$missing_genes, panels$b$missing_genes)
    if (length(missing) > 0) {
      text <- paste0(text, ". Missing from dataset: ", paste(missing, collapse = ", "))
    }
    text
  })

  output$violin_plot_a <- renderPlot({
    panels <- violin_result()
    req(panels)
    panels$a$plot
  })

  output$violin_plot_b <- renderPlot({
    panels <- violin_result()
    req(panels)
    panels$b$plot
  })

  output$download_violin <- downloadHandler(
    filename = function() {
      "cd4_vs_cd8_t_cell_expression_panels.png"
    },
    content = function(file) {
      panels <- violin_result()
      req(panels)
      combined <- panels$a$plot | panels$b$plot
      ggplot2::ggsave(
        filename = file,
        plot = combined,
        width = 16,
        height = 24,
        dpi = 150
      )
    },
    contentType = "image/png"
  )
}

shinyApp(ui = ui, server = server)
