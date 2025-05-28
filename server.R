library(shiny)
library(heatmaply)
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(purrr)
library(tibble)
library(digest)

GO <- readRDS("data/gene_sets_human.rds")

function(input, output, session) {
  
  updateSelectInput(session, "db_select", choices = names(GO), selected = "GO")
  
  logcpm_data <- reactive({
    req(input$logcpm_file)
    df <- read_csv(input$logcpm_file$datapath)
    gene_ids <- df[[1]]
    keep <- !is.na(gene_ids) & gene_ids != ""
    df <- df[keep, , drop = FALSE]
    rownames(df) <- make.unique(df[[1]])
    df[[1]] <- NULL
    as.matrix(df)
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
    pickerInput("selected_pathways", "Select Pathway(s)", 
                choices = names(GO[[input$db_select]]), 
                multiple = TRUE,
                options = list(`live-search` = TRUE, size = 10))
  })
  
  plot_titles <- reactiveVal(character())
  plot_ids <- reactiveVal(character())
  
  observeEvent(input$generate_plot, {
    req(input$selected_pathways)
    mat <- isolate(logcpm_data())
    sample_ann <- isolate(sample_info())
    current_titles <- plot_titles()
    current_ids <- plot_ids()
    
    new_outputs <- lapply(seq_along(input$selected_pathways), function(i) {
      pathway <- input$selected_pathways[[i]]
      uid <- paste0("plot_", digest::digest(paste(pathway, Sys.time(), runif(1))))
      plot_title <- if (nzchar(input$sample_desc)) input$sample_desc else paste("Heatmap:", pathway)
      
      local({
        this_uid <- uid
        this_title <- plot_title
        this_pathway <- pathway
        
        output[[this_uid]] <- renderPlotly({
          genes <- GO[[input$db_select]][[this_pathway]]
          matched <- intersect(rownames(mat), genes)
          sub_mat <- mat[matched, , drop = FALSE]
          
          validate(need(nrow(sub_mat) >= 2, paste("Too few genes in pathway:", this_pathway)))
          
          sub_mat <- t(scale(t(sub_mat)))
          ann_df <- sample_ann %>% column_to_rownames("sample_name")
          ann_df <- ann_df[colnames(sub_mat), , drop = FALSE]
          
          heatmaply(sub_mat, 
                    col_side_colors = ann_df,
                    colors = colorRampPalette(c("red", "white", "blue"))(256),
                    main = this_title,
                    plot_method = "plotly"
          )
        })
      })
      
      list(uid = uid, title = plot_title)
    })
    
    plot_titles(c(current_titles, sapply(new_outputs, `[[`, "title")))
    plot_ids(c(current_ids, sapply(new_outputs, `[[`, "uid")))
  })
  
  output$heatmap_outputs <- renderUI({
    ids <- plot_ids()
    titles <- plot_titles()
    if (length(ids) == 0) return(NULL)
    
    tagList(lapply(seq_along(ids), function(i) {
      tagList(
        h4(titles[i]),
        plotlyOutput(ids[i]),
        actionButton(paste0("remove_plot_", i), "Remove")
      )
    }))
  })
  
  observe({
    req(plot_ids())
    lapply(seq_along(plot_ids()), function(i) {
      observeEvent(input[[paste0("remove_plot_", i)]], {
        titles_list <- plot_titles()
        ids_list <- plot_ids()
        
        keep <- seq_along(ids_list) != i
        plot_titles(titles_list[keep])
        plot_ids(ids_list[keep])
      }, ignoreInit = TRUE)
    })
  })
}
