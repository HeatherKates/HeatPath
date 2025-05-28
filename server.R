library(shiny)
library(heatmaply)
library(dplyr)
library(tidyr)
library(shinyWidgets)
library(tibble)
library(digest)

GO <- readRDS("data/gene_sets_human.rds")

server <- function(input, output, session) {
  updateSelectInput(session, "db_select", choices = names(GO), selected = names(GO)[1])
  
  observeEvent(input$db_select, {
    pathways <- names(GO[[input$db_select]])
    updatePickerInput(session, "selected_pathways", choices = pathways, selected = character(0))
  })
  
  logcpm_data <- reactive({
    req(input$logcpm_file)
    df <- read.csv(input$logcpm_file$datapath, check.names = FALSE)
    gene_ids <- df[[1]]
    keep <- !is.na(gene_ids) & gene_ids != ""
    df <- df[keep, , drop = FALSE]
    rownames(df) <- make.unique(as.character(df[[1]]))
    df[[1]] <- NULL
    as.matrix(df)
  })
  
  sample_info <- reactive({
    req(input$sampleinfo_file)
    df <- read.csv(input$sampleinfo_file$datapath, check.names = FALSE)
    colnames(df)[1] <- c("sample_name")
    df <- df %>% filter(sample_name %in% colnames(logcpm_data()))
    df
  })
  
  plot_data <- reactiveVal(list())
  plotted_pathways <- reactiveVal(character())
  
  observeEvent(input$generate_plot, {
    req(input$selected_pathways)
    mat <- logcpm_data()
    ann_df <- sample_info()
    
    current_plots <- plot_data()
    previous <- plotted_pathways()
    new_pathways <- setdiff(input$selected_pathways, previous)
    if (length(new_pathways) == 0) return()
    
    for (pathway in new_pathways) {
      uid <- paste0("plot_", digest::digest(paste0(pathway, Sys.time(), runif(1))))
      plot_title <- if (nzchar(input$sample_desc)) input$sample_desc else paste("Heatmap:", pathway)
      
      genes <- GO[[input$db_select]][[pathway]]
      matched <- intersect(rownames(mat), genes)
      sub_mat_raw <- mat[matched, , drop = FALSE]
      sub_mat <- t(scale(t(sub_mat_raw)))
      
      if (nrow(sub_mat) < 2) next
      
      ann_df2 <- sample_info() %>% column_to_rownames("sample_name")
      ann_df2 <- ann_df2[colnames(sub_mat), , drop = FALSE]
      
      # Store matrix BEFORE rendering
      current_plots[[uid]] <- list(title = plot_title, pathway = pathway, matrix = sub_mat_raw)
      
      local({
        local_uid <- uid
        local_title <- plot_title
        
        output[[local_uid]] <- renderPlotly({
          heatmaply(sub_mat,
                    col_side_colors = ann_df2,
                    colors = colorRampPalette(c("red", "white", "blue"))(256),
                    main = local_title,
                    plot_method = "plotly",
                    scale_fill_gradient_fun = ggplot2::scale_fill_gradient2,
                    margins = c(60, 120, 60, 60),
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
        
        output[[paste0("download_matrix_", local_uid)]] <- downloadHandler(
          filename = function() {
            paste0(gsub("[^a-zA-Z0-9]", "_", local_title), ".csv")
          },
          content = function(file) {
            plot_list <- plot_data()
            write.csv(plot_list[[local_uid]]$matrix, file)
          }
        )
      })
    }
    
    plot_data(current_plots)
    plotted_pathways(union(previous, new_pathways))
    updatePickerInput(session, "selected_pathways", selected = character(0))
    updateTextInput(session, "sample_desc", value = "")
  })
  
  output$heatmap_outputs <- renderUI({
    all_plots <- plot_data()
    if (length(all_plots) == 0) return(NULL)
    
    tagList(
      lapply(seq_along(all_plots), function(i) {
        uid <- names(all_plots)[i]
        pathway <- all_plots[[i]]$pathway
        tagList(
          tags$h4(paste("Heatmap:", pathway)),
          plotlyOutput(uid, height = "500px"),
          downloadButton(paste0("download_matrix_", uid), "Download matrix as CSV"),
          actionButton(paste0("remove_plot_", uid), "Remove Plot"),
          tags$hr()
        )
      })
    )
  })
  
  observe({
    lapply(names(plot_data()), function(uid) {
      observeEvent(input[[paste0("remove_plot_", uid)]], {
        current <- plot_data()
        removed_pathway <- current[[uid]]$pathway
        current[[uid]] <- NULL
        plot_data(current)
        current_paths <- plotted_pathways()
        plotted_pathways(setdiff(current_paths, removed_pathway))
      }, ignoreInit = TRUE, once = TRUE)
    })
  })
}
