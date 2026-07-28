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
        "Total cells (random subsample; balanced by sample when possible)",
        min = 500,
        max = 18000,
        value = 6000,
        step = 500
      ),
      helpText("Lower = faster. Uses the same cell count for Run A and Run B so comparisons stay fair."),
      hr(),
      h4("Shared display"),
      selectInput(
        "color_by",
        "Color UMAP by",
        choices = c(
          "Clusters" = "seurat_clusters",
          "Treatment" = "treatment",
          "Sample" = "sample_id",
          "Donor" = "donor",
          "Inferred CD4/CD8" = "t_lineage"
        ),
        selected = "seurat_clusters"
      ),
      checkboxInput("label_clusters", "Label clusters (when coloring by cluster)", TRUE),
      hr(),
      h4("CD4+ / CD8+ split"),
      checkboxInput(
        "show_lineage_split",
        "Show separate UMAPs for inferred CD4+ vs CD8+",
        value = TRUE
      ),
      selectInput(
        "lineage_from_run",
        "Split which run?",
        choices = c("Run A" = "A", "Run B" = "B"),
        selected = "A"
      ),
      sliderInput(
        "lineage_min_score",
        "Min gene score to call CD4/CD8 (Ambiguous if below)",
        min = 0.05,
        max = 1.0,
        value = 0.25,
        step = 0.05
      ),
      helpText(
        "This GEO scRNA-seq series does not include FACS CD4/CD8 labels. ",
        "Lineage is inferred from CD4 vs mean(CD8A, CD8B) expression."
      ),
      hr(),
      h4("Run A"),
      sliderInput("res_a", "Cluster resolution", min = 0.1, max = 1.5, value = 0.5, step = 0.1),
      sliderInput("dims_a", "PCA / UMAP dims", min = 5, max = 40, value = 30, step = 1),
      numericInput("seed_a", "Random seed", value = 42, min = 1, step = 1),
      actionButton("run_a", "Compute Run A", class = "btn-success"),
      hr(),
      h4("Run B"),
      sliderInput("res_b", "Cluster resolution", min = 0.1, max = 1.5, value = 0.8, step = 0.1),
      sliderInput("dims_b", "PCA / UMAP dims", min = 5, max = 40, value = 30, step = 1),
      numericInput("seed_b", "Random seed", value = 42, min = 1, step = 1),
      actionButton("run_b", "Compute Run B", class = "btn-success"),
      hr(),
      helpText(
        "Loads a saved pipeline object, then re-runs neighbors → clusters → UMAP ",
        "with your settings. Progress bar shows stage + estimated time left."
      )
    ),
    mainPanel(
      width = 9,
      textOutput("status"),
      tags$hr(),
      fluidRow(
        column(6, h3("Run A"), plotOutput("plot_a", height = "480px"), verbatimTextOutput("summary_a")),
        column(6, h3("Run B"), plotOutput("plot_b", height = "480px"), verbatimTextOutput("summary_b"))
      ),
      tags$hr(),
      h4("Cluster × treatment counts"),
      fluidRow(
        column(6, tableOutput("table_a")),
        column(6, tableOutput("table_b"))
      ),
      tags$hr(),
      h3("CD4+ vs CD8+ UMAPs (transcriptome-inferred)"),
      textOutput("lineage_status"),
      tableOutput("lineage_table"),
      fluidRow(
        column(6, plotOutput("plot_cd4", height = "420px")),
        column(6, plotOutput("plot_cd8", height = "420px"))
      )
    )
  )
)

server <- function(input, output, session) {
  base_obj <- reactiveVal(NULL)
  run_a <- reactiveVal(NULL)
  run_b <- reactiveVal(NULL)
  status_msg <- reactiveVal("Load a Seurat RDS to begin.")

  output$status <- renderText(status_msg())

  run_with_progress <- function(run_name, obj, resolution, dims, seed, n_cells) {
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
        setProgress(1, detail = "Loaded")
      })
      base_obj(obj)
      run_a(NULL)
      run_b(NULL)

      n <- ncol(obj)
      updateSliderInput(
        session,
        "n_cells",
        min = min(500, n),
        max = n,
        value = min(6000, n),
        step = max(100, round(n / 40))
      )

      status_msg(sprintf(
        "Loaded %s cells from %s. Choose how many cells to include, then Compute Run A / Run B.",
        format(n, big.mark = ","),
        path
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
        n_cells = input$n_cells
      )
      run_a(out)
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
        n_cells = input$n_cells
      )
      run_b(out)
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

  output$plot_a <- renderPlot({
    out <- run_a()
    req(out)
    obj <- out$obj
    if (identical(input$color_by, "t_lineage")) {
      obj <- assign_cd4_cd8(obj, min_score = input$lineage_min_score)
    }
    plot_umap(
      obj,
      group_by = input$color_by,
      title = paste("A:", param_label(out$resolution, out$dims, out$seed, out$n_cells)),
      label = input$label_clusters
    )
  })

  output$plot_b <- renderPlot({
    out <- run_b()
    req(out)
    obj <- out$obj
    if (identical(input$color_by, "t_lineage")) {
      obj <- assign_cd4_cd8(obj, min_score = input$lineage_min_score)
    }
    plot_umap(
      obj,
      group_by = input$color_by,
      title = paste("B:", param_label(out$resolution, out$dims, out$seed, out$n_cells)),
      label = input$label_clusters
    )
  })

  output$summary_a <- renderText({
    out <- run_a()
    if (is.null(out)) {
      return("Not computed yet.")
    }
    paste0(
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

  lineage_obj <- reactive({
    req(isTRUE(input$show_lineage_split))
    out <- if (identical(input$lineage_from_run, "B")) run_b() else run_a()
    req(out)
    assign_cd4_cd8(out$obj, min_score = input$lineage_min_score)
  })

  lineage_plots <- reactive({
    obj <- lineage_obj()
    plot_umap_lineage_split(obj, group_by = input$color_by, label = input$label_clusters)
  })

  output$lineage_status <- renderText({
    if (!isTRUE(input$show_lineage_split)) {
      return("Enable “Show separate UMAPs…” to compute CD4/CD8 split.")
    }
    out <- if (identical(input$lineage_from_run, "B")) run_b() else run_a()
    if (is.null(out)) {
      return(paste("Compute", input$lineage_from_run, "first, then the CD4/CD8 UMAPs will appear."))
    }
    obj <- lineage_obj()
    counts <- table(obj$t_lineage)
    n_cd4 <- if ("CD4" %in% names(counts)) as.integer(counts[["CD4"]]) else 0L
    n_cd8 <- if ("CD8" %in% names(counts)) as.integer(counts[["CD8"]]) else 0L
    n_amb <- if ("Ambiguous" %in% names(counts)) as.integer(counts[["Ambiguous"]]) else 0L
    sprintf(
      "From %s | CD4=%s | CD8=%s | Ambiguous=%s (Ambiguous cells hidden from the two plots)",
      input$lineage_from_run,
      format(n_cd4, big.mark = ","),
      format(n_cd8, big.mark = ","),
      format(n_amb, big.mark = ",")
    )
  })

  output$lineage_table <- renderTable({
    req(isTRUE(input$show_lineage_split))
    obj <- lineage_obj()
    tab <- table(Lineage = obj$t_lineage, Treatment = obj$treatment)
    as.data.frame.matrix(tab)
  }, rownames = TRUE)

  output$plot_cd4 <- renderPlot({
    req(isTRUE(input$show_lineage_split))
    lineage_plots()$cd4
  })

  output$plot_cd8 <- renderPlot({
    req(isTRUE(input$show_lineage_split))
    lineage_plots()$cd8
  })
}

shinyApp(ui = ui, server = server)
