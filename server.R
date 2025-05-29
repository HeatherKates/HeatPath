library(shiny)
library(heatmaply)
library(dplyr)
library(tidyr)
library(shinyWidgets)
library(tibble)
library(digest)
id_converter <- readr::read_tsv("data/msigdbr_id_converter.txt", show_col_types = FALSE)

GO <- readRDS("data/gene_sets_human.rds")

server <- function(input, output, session) {
  output$download_example_logcpm <- downloadHandler(
    filename = function() {
      "example_logCPM.csv"
    },
    content = function(file) {
      file.copy("test_data/test_lcpm.csv", file)
    }
  )
  
  output$download_example_sampleinfo <- downloadHandler(
    filename = function() {
      "example_sample_info.csv"
    },
    content = function(file) {
      file.copy("test_data/test_sample_info.csv", file)
    }
  )
  
  updateSelectInput(session, "db_select", choices = names(GO), selected = names(GO)[1])
  
  observeEvent(input$db_select, {
    updateTextInput(session, "pathway_search", value = "")
    output$search_results <- renderUI(NULL)
  })
  
  observeEvent(input$search_button, {
    req(input$pathway_search)
    db <- input$db_select
    search_term <- tolower(input$pathway_search)
    metadata <- GO[[db]]$metadata
    
    matches <- metadata %>%
      filter(grepl(search_term, tolower(searchable_text), fixed = TRUE))
    
    if (nrow(matches) == 0) {
      output$search_results <- renderUI(tags$p("No pathways found."))
    } else {
      output$search_results <- renderUI({
        checkboxGroupInput("selected_pathways", "Select Pathways:", 
                           choices = setNames(matches$set_id, matches$name_clean),
                           selected = NULL)
      })
    }
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
  output$download_all_html <- downloadHandler(
    filename = function() {
      paste0("heatmaps_", Sys.Date(), ".html")
    },
    content = function(file) {
      temp_rmd <- tempfile(fileext = ".Rmd")
      temp_dir <- tempdir()
      plots <- plot_data()
      
      if (length(plots) == 0) {
        writeLines("<h3>No plots to export.</h3>", file)
        return()
      }
      
      # Save each scaled matrix to a separate .rds file
      matrix_paths <- lapply(names(plots), function(uid) {
        mat_path <- file.path(temp_dir, paste0(uid, ".rds"))
        saveRDS(t(scale(t(plots[[uid]]$matrix))), mat_path)
        list(title = plots[[uid]]$title, rds_path = mat_path)
      })
      
      # Create the Rmd file content
      rmd_lines <- c(
        "---",
        "title: 'All Heatmaps'",
        "output: html_document",
        "---",
        "",
        "```{r setup, include=FALSE}",
        "library(heatmaply)",
        "library(grDevices)",
        "```",
        ""
      )
      
      for (entry in matrix_paths) {
        title <- entry$title
        abs_path <- normalizePath(entry$rds_path, winslash = "/")  # works on all OS
        rmd_lines <- c(
          rmd_lines,
          paste0("## ", title),
          "",
          "```{r, echo=FALSE, message=FALSE, warning=FALSE}",
          paste0("mat <- readRDS(\"", abs_path, "\")"),
          "heatmaply(mat,",
          "  plot_method = 'plotly',",
          "  colors = colorRampPalette(c('red', 'white', 'blue'))(256)",
          ")",
          "```",
          ""
        )
      }
      
      writeLines(rmd_lines, temp_rmd)
      
      rmarkdown::render(
        input = temp_rmd,
        output_file = file,
        envir = new.env(),
        quiet = TRUE
      )
    }
  )
  
  
  
  plotted_pathways <- reactiveVal(character())
  observeEvent(input$pathway_file, {
    req(input$pathway_file)
    
    db <- input$db_select
    uploaded_ids <- readLines(input$pathway_file$datapath, warn = FALSE) %>%
      trimws() %>%
      discard(~ .x == "")
    
    # Filter the converter to just this database
    db_converter <- id_converter %>% filter(db == db)
    
    # Match uploaded external IDs to msigdbr IDs
    matched <- db_converter %>% filter(external_id %in% uploaded_ids)
    unmatched <- setdiff(uploaded_ids, matched$external_id)
    
    if (nrow(matched) == 0) {
      output$search_results <- renderUI(tags$p("No matching pathways found in MSigDB for the uploaded IDs."))
    } else {
      # Make named vector: msigdbr_id = label
      id_choices <- setNames(matched$msigdbr_id, matched$gs_name)
      
      output$search_results <- renderUI({
        checkboxGroupInput("selected_pathways", "Select Pathways from Uploaded File:",
                           choices = id_choices, selected = matched$msigdbr_id)
      })
    }
    
    if (length(unmatched) > 0) {
      showNotification(
        paste0("The following ", length(unmatched), " IDs were not matched and excluded: ", paste(unmatched, collapse = ", ")),
        type = "warning",
        duration = 10
      )
    }
    
  })
  
  
  
  observeEvent(input$generate_plot, {
    output$search_results <- renderUI(NULL)
    
    req(input$selected_pathways)
    mat <- logcpm_data()
    ann_df <- sample_info()
    
    current_plots <- plot_data()
    previous <- plotted_pathways()
    new_pathways <- setdiff(input$selected_pathways, previous)
    if (length(new_pathways) == 0) return()
    
    for (pathway in new_pathways) {
      uid <- paste0("plot_", digest::digest(paste0(pathway, Sys.time(), runif(1))))
      # Look up human-readable name from metadata if available
      db_meta <- GO[[input$db_select]]$metadata
      pretty_name <- db_meta$name_clean[match(pathway, db_meta$set_id)]
      
      plot_title <- if (nzchar(input$sample_desc)) input$sample_desc else paste("Heatmap:", pretty_name %||% pathway)
      
      
      genes <- GO[[input$db_select]]$genes[[pathway]]
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
    
    updateCheckboxGroupInput(session, "selected_pathways", selected = character(0))
    updateTextInput(session, "sample_desc", value = "")
    updateTextInput(session, "pathway_search", value = "")  # Optional
    
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
