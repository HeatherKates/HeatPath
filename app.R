library(shiny)
library(heatmaply)
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(purrr)
library(shinyWidgets)
library(tibble)
library(digest)

# Load gene sets
GO <- readRDS("data/gene_sets_human.rds")

ui <- fluidPage(
  titlePanel("Pathway Heatmaps (GO / KEGG)"),
  sidebarLayout(
    sidebarPanel(
      fileInput("logcpm_file", "Upload logCPM Matrix (.csv)", accept = ".csv"),
      fileInput("sampleinfo_file", "Upload Sample Info (.csv)", accept = ".csv"),
      selectInput("db_select", "Select Database", choices = names(GO), selected = "GO"),
      uiOutput("pathway_selector"),
      textInput("sample_desc", "Plot Title (optional)", value = ""),
      actionButton("generate_plot", "Generate Plot"),
      width = 3
    ),
    mainPanel(
      uiOutput("heatmap_outputs")
    )
  )
)

server <- function(input, output, session) {
  
  logcpm_data <- reactive({
    req(input$logcpm_file)
    df <- read_csv(input$logcpm_file$datapath)
    mat <- as.data.frame(df)
    gene_ids <- mat[[1]]
    keep <- !is.na(gene_ids) & gene_ids != ""
    mat <- mat[keep, , drop = FALSE]
    rownames(mat) <- make.unique(mat[[1]])
    mat[[1]] <- NULL
    as.matrix(mat)
  })
  
  sample_info <- reactive({
    req(input$sampleinfo_file)
    df <- read_csv(input$sampleinfo_file$datapath)
    colnames(df)[1:2] <- c("sample_name", "group")
    df <- df %>% filter(sample_name %in% colnames(logcpm_data()))
    df
  })
  
  output$pathway_selector <- renderUI({
    req(input$db_select)
    pathways <- names(GO[[input$db_select]])
    pickerInput("selected_pathways", "Select Pathway(s)", 
                choices = pathways, multiple = TRUE, 
                options = list(`live-search` = TRUE, size = 10))
  })
  
  plots <- reactiveVal(list())
  plot_titles <- reactiveVal(character())
  plot_ids <- reactiveVal(character())
  
  observeEvent(input$generate_plot, {
    req(input$selected_pathways)
    mat <- isolate(logcpm_data())
    sample_ann <- isolate(sample_info())
    current_plots <- plots()
    current_titles <- plot_titles()
    current_ids <- plot_ids()
    
    new_outputs <- lapply(seq_along(input$selected_pathways), function(i) {
      pathway <- input$selected_pathways[[i]]
      uid <- paste0("plot_", digest::digest(paste(pathway, Sys.time(), runif(1))))
      plot_title <- if (nzchar(input$sample_desc)) input$sample_desc else paste("Heatmap of Gene Expression (row-scaled logCPM) in pathway", pathway)
      
      local({
        this_uid <- uid
        this_title <- plot_title
        this_pathway <- pathway
        
        output[[this_uid]] <- renderPlotly({
          genes <- GO[[input$db_select]][[this_pathway]]
          matched <- intersect(rownames(mat), genes)
          sub_mat <- mat[matched, , drop = FALSE]
          
          if (nrow(sub_mat) < 2) {
            validate(paste("Too few genes in pathway:", this_pathway))
          }
          
          sub_mat <- t(scale(t(sub_mat)))
          ann_df <- sample_ann %>% column_to_rownames("sample_name")
          ann_df <- ann_df[colnames(sub_mat), , drop = FALSE]
          
          heatmaply(sub_mat, 
                    col_side_colors = ann_df,
                    colors = colorRampPalette(c("red", "white", "blue"))(256),
                    main = this_title,
                    plot_method = "plotly",
                    scale_fill_gradient_fun = ggplot2::scale_fill_gradient2,
                    margins = c(60, 120, 40, 20),
                    key.title = "Z-score",
                    fontsize_row = 10,
                    dendrogram = "both",
                    showticklabels = c(TRUE, nrow(sub_mat) <= 20),
                    grid_color = "transparent",
                    cellnote = NULL,
                    row_text_angle = 0,
                    column_text_angle = 45,
                    labCol = colnames(sub_mat),
                    labRow = if (nrow(sub_mat) <= 20) rownames(sub_mat) else NULL,
                    fontsize_col = 10,
                    fontsize = 12)
        })
      })
      
      list(uid = uid, title = plot_title)
    })
    
    plots(c(current_plots, lapply(new_outputs, function(x) output[[x$uid]])))
    plot_titles(c(current_titles, sapply(new_outputs, `[[`, "title")))
    plot_ids(c(current_ids, sapply(new_outputs, `[[`, "uid")))
  })
  
  output$heatmap_outputs <- renderUI({
    ids <- plot_ids()
    titles <- plot_titles()
    if (length(ids) == 0) return(NULL)
    
    tagList(lapply(seq_along(ids), function(i) {
      plotlyOutput(ids[i])
      actionButton(paste0("remove_plot_", i), paste("Remove Plot", i))
    }))
  })
  
  observe({
    req(plot_ids())
    lapply(seq_along(plot_ids()), function(i) {
      observeEvent(input[[paste0("remove_plot_", i)]], {
        plots_list <- plots()
        titles_list <- plot_titles()
        ids_list <- plot_ids()
        
        keep <- setdiff(seq_along(ids_list), i)
        plots(plots_list[keep])
        plot_titles(titles_list[keep])
        plot_ids(ids_list[keep])
      }, ignoreInit = TRUE)
    })
  })
}

shinyApp(ui, server)
