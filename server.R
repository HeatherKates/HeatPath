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

GO <- readRDS("data/gene_sets_human.rds")

server <- function(input, output, session) {
  # Update the database select choices
  updateSelectInput(session, "db_select", choices = names(GO), selected = names(GO)[1])
  
  # Dynamically update the pathway choices based on db_select
  observeEvent(input$db_select, {
    pathways <- names(GO[[input$db_select]])
    updatePickerInput(session, "selected_pathways",
                      choices = pathways,
                      selected = character(0))  # Reset selection

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
    
    # Only plot new selections
    new_pathways <- setdiff(input$selected_pathways, previous)
    if (length(new_pathways) == 0) return()
    
    for (pathway in new_pathways) {
      uid <- paste0("plot_", digest::digest(paste0(pathway, Sys.time(), runif(1))))
      plot_title <- if (nzchar(input$sample_desc)) input$sample_desc else paste("Heatmap:", pathway)
      
      genes <- GO[[input$db_select]][[pathway]]
      matched <- intersect(rownames(mat), genes)
      sub_mat <- mat[matched, , drop = FALSE]
      
      if (nrow(sub_mat) < 2) next  # skip invalid plots
      
      sub_mat <- t(scale(t(sub_mat)))
      ann_df2 <- ann_df %>% column_to_rownames("sample_name")
      ann_df2 <- ann_df2[colnames(sub_mat), , drop = FALSE]
      
      local({
        local_uid <- uid
        local_title <- plot_title
        local_mat <- sub_mat
        local_ann <- ann_df2
        
        output[[local_uid]] <- renderPlotly({
          heatmaply(local_mat,
                    col_side_colors = local_ann,
                    colors = colorRampPalette(c("red", "white", "blue"))(256),
                    main = local_title,
                    plot_method = "plotly",
                    scale_fill_gradient_fun = ggplot2::scale_fill_gradient2,
                    margins = c(60, 120, 60, 60),
                    key.title = "Z-score",
                    fontsize_row = 10,
                    dendrogram = "both",
                    showticklabels = c(TRUE, nrow(local_mat) <= 20),
                    grid_color = "transparent",
                    cellnote = NULL,
                    row_text_angle = 0,
                    column_text_angle = 45,
                    labCol = colnames(local_mat),
                    labRow = if (nrow(local_mat) <= 20) rownames(local_mat) else NULL,
                    fontsize_col = 10,
                    fontsize = 12) # Re-centers vertically
        })
      })
      
      current_plots[[uid]] <- list(title = plot_title, pathway = pathway)
    }
    
    # Update stored plots and plotted pathways
    plot_data(current_plots)
    plotted_pathways(union(previous, new_pathways))
    
    # Reset picker selection (this is key to fixing the UI behavior)
    updatePickerInput(session, "selected_pathways", selected = character(0))
    # Reset plot title input
    updateTextInput(session, "sample_desc", value = "")
  })
  
  output$heatmap_outputs <- renderUI({
    all_plots <- plot_data()
    if (length(all_plots) == 0) return(NULL)
    
    tagList(
      lapply(seq_along(all_plots), function(i) {
        uid <- names(all_plots)[i]
        pathway <- all_plots[[i]]$pathway
        custom_title <- all_plots[[i]]$title
        tagList(
          tags$h4(paste("Heatmap:", pathway)),  # Top-left title
          plotlyOutput(uid, height = "500px"),
          actionButton(paste0("remove_plot_", uid), "Remove Plot"),
          tags$hr()
        )
      })
    )
    
  })
  
  # Remove logic
  observe({
    lapply(names(plot_data()), function(uid) {
      observeEvent(input[[paste0("remove_plot_", uid)]], {
        current <- plot_data()
        removed_pathway <- current[[uid]]$pathway
        current[[uid]] <- NULL
        plot_data(current)
        
        # Update plotted_pathways too
        current_paths <- plotted_pathways()
        plotted_pathways(setdiff(current_paths, removed_pathway))
      }, ignoreInit = TRUE, once = TRUE)
    })
  })
}
