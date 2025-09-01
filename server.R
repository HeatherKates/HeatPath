library(shiny)
library(heatmaply)
library(dplyr)
library(tidyr)
library(shinyWidgets)
library(tibble)
library(digest)
library(DT)

GO <- readRDS("data/gene_sets_human.rds")

library(biomaRt)  # Add this to your library imports

# Add this function at the top of your file
detect_organism_and_add_symbols <- function(gene_ids) {
  # Detect organism from Ensembl ID pattern
  if (any(grepl("^ENSG", gene_ids))) {
    organism <- "human"
    dataset <- "hsapiens_gene_ensembl"
  } else if (any(grepl("^ENSMUSG", gene_ids))) {
    organism <- "mouse"
    dataset <- "mmusculus_gene_ensembl"
  } else {
    return(NULL)  # Unsupported organism
  }
  tryCatch({
    # Use biomaRt to get gene symbols
    mart <- useMart("ensembl", dataset = dataset)
    # Get gene symbols
    gene_info <- getBM(
      attributes = c("ensembl_gene_id", "external_gene_name"),
      filters = "ensembl_gene_id",
      values = gene_ids,
      mart = mart
    )
    return(list(
      organism = organism,
      gene_symbols = setNames(gene_info$external_gene_name, gene_info$ensembl_gene_id)
    ))
  }, error = function(e) {
    cat("Error fetching gene symbols:", e$message, "\n")
    return(NULL)
  })
}

# Add reactive values for demo data
demo_count_data <- reactiveVal(NULL)
demo_meta_data <- reactiveVal(NULL)
demo_loaded <- reactiveVal(FALSE)

server <- function(input, output, session) {
  # Update the database select choices
  updateSelectInput(session, "db_select", choices = names(GO), selected = names(GO)[1])
  
  # Dynamically update the pathway choices based on db_select
  # Populate pathways when database is selected
  observeEvent(input$db_select, {
    req(input$db_select)
    
    if ("metadata" %in% names(GO[[input$db_select]])) {
      # New structure with metadata
      db_meta <- GO[[input$db_select]]$metadata
      # Use gs_id as value, name_clean as display text for searching
      pathway_choices <- setNames(db_meta$gs_id, db_meta$name_clean)
    } else {
      # Old flat structure  
      pathway_choices <- names(GO[[input$db_select]])
    }
    
    updatePickerInput(session, "selected_pathways",
                      choices = pathway_choices,
                      selected = character(0))
  })
  
  logcpm_data <- reactive({
    if (demo_loaded()) {
      return(demo_count_data())
    }
    
    req(input$logcpm_file)
    df <- read.csv(input$logcpm_file$datapath, check.names = FALSE)
    gene_ids <- df[[1]]
    keep <- !is.na(gene_ids) & gene_ids != ""
    df <- df[keep, , drop = FALSE]
    gene_ids <- gene_ids[keep]
    
    # Check if we have ENSEMBL IDs and convert to symbols
    if (any(grepl("^ENS", gene_ids))) {
      showNotification("ENSEMBL IDs detected. Converting to gene symbols...", 
                       type = "message", duration = 3)
      
      conversion_result <- detect_organism_and_add_symbols(gene_ids)
      
      if (!is.null(conversion_result)) {
        showNotification(paste("Detected", conversion_result$organism, "data. Conversion complete!"), 
                         type = "message", duration = 3)
        
        # Map ENSEMBL IDs to symbols
        gene_symbols <- conversion_result$gene_symbols[gene_ids]
        
        # Remove genes without symbol matches and duplicates
        valid_symbols <- !is.na(gene_symbols) & gene_symbols != ""
        df <- df[valid_symbols, , drop = FALSE]
        gene_symbols <- gene_symbols[valid_symbols]
        
        # Handle duplicate symbols by making them unique
        rownames(df) <- make.unique(gene_symbols)
      } else {
        showNotification("Could not convert ENSEMBL IDs. Using original IDs.", 
                         type = "warning", duration = 5)
        rownames(df) <- make.unique(as.character(gene_ids))
      }
    } else {
      # Assume gene symbols, use as-is
      rownames(df) <- make.unique(as.character(gene_ids))
    }
    
    df[[1]] <- NULL
    as.matrix(df)
  })
  
  sample_info <- reactive({
    if (demo_loaded()) {
      return(demo_meta_data())
    }
    
    req(input$sampleinfo_file)
    df <- read.csv(input$sampleinfo_file$datapath, check.names = FALSE)
    colnames(df)[1] <- c("sample_name")
    df <- df %>% filter(sample_name %in% colnames(logcpm_data()))
    df
  })
  
  plot_data <- reactiveVal(list())
  plotted_pathways <- reactiveVal(character())
  
  # Demo data loader
  observeEvent(input$load_demo, {
    tryCatch({
      # Define demo data paths
      demo_dir <- "demo_data"
      count_file <- file.path(demo_dir, "GSE100688_logCPM_data_2025-09-01.csv")
      meta_file_path <- file.path(demo_dir, "GSE100688_filtered_metadata.csv")
      
      # Check if demo files exist
      if (!file.exists(count_file) || !file.exists(meta_file_path)) {
        showNotification("Demo data files not found. Please ensure demo_data folder contains the required files.",
                         type = "error", duration = 5)
        return()
      }
      
      # Load count matrix (ENSEMBL format from reference app)
      count_df <- read.csv(count_file, check.names = FALSE)
      gene_ids <- count_df[[1]]  # First column should be ENSEMBL IDs
      
      # Remove first column and set up matrix
      count_df[[1]] <- NULL
      
      # Convert ENSEMBL to symbols
      showNotification("Loading demo data and converting ENSEMBL IDs to gene symbols...", 
                       type = "message", duration = 3)
      
      conversion_result <- detect_organism_and_add_symbols(gene_ids)
      
      if (!is.null(conversion_result)) {
        gene_symbols <- conversion_result$gene_symbols[gene_ids]
        
        # Keep only genes with valid symbols
        valid_symbols <- !is.na(gene_symbols) & gene_symbols != ""
        count_df <- count_df[valid_symbols, , drop = FALSE]
        gene_symbols <- gene_symbols[valid_symbols]
        
        # Set rownames to gene symbols
        rownames(count_df) <- make.unique(gene_symbols)
        count_matrix_final <- as.matrix(count_df)
      } else {
        showNotification("Could not convert ENSEMBL IDs in demo data.", type = "error", duration = 5)
        return()
      }
      
      # Load metadata and format for current app
      meta_df <- read.csv(meta_file_path, check.names = FALSE)
      colnames(meta_df)[1] <- "sample_name"
      if (ncol(meta_df) > 1) {
        colnames(meta_df)[2] <- "group"
      }
      
      # Filter metadata to match samples in count matrix
      meta_df <- meta_df %>% filter(sample_name %in% colnames(count_matrix_final))
      
      # Set demo data
      demo_count_data(count_matrix_final)
      demo_meta_data(meta_df)
      demo_loaded(TRUE)
      
      showNotification("Demo data loaded successfully! 🎉", type = "message", duration = 3)
      
    }, error = function(e) {
      showNotification(paste("Error loading demo data:", e$message), type = "error", duration = 5)
    })
  })
  
  # Demo info modal
  observeEvent(input$info_demo, {
    showModal(modalDialog(
      title = "Demo Data Information",
      div(
        p("Demo data obtained for study ", strong("GSE100688"), " from ",
          strong("GREIN: GEO RNA-seq Experiments Interactive Navigator")),
        br(),
        p("The data contains ENSEMBL gene IDs which will be automatically converted to gene symbols."),
        br(),
        p("Visit GREIN:",
          a("https://dev.ilincs.org/apps/grein/?gse=GSE100688",
            href = "https://dev.ilincs.org/apps/grein/?gse=GSE100688",
            target = "_blank",
            style = "color: #007bff;"))
      ),
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })
  
  observeEvent(input$generate_plot, {
    req(input$selected_pathways)
    mat <- logcpm_data()
    ann_df <- sample_info()
    current_plots <- plot_data()
    previous <- plotted_pathways()
    
    new_pathways <- setdiff(input$selected_pathways, previous)
    
    for (pathway in new_pathways) {
      # Get genes and display name based on structure
      if ("metadata" %in% names(GO[[input$db_select]])) {
        # New structure with metadata
        genes <- GO[[input$db_select]]$genes[[pathway]]
        db_meta <- GO[[input$db_select]]$metadata
        display_name <- db_meta$name_clean[db_meta$gs_id == pathway]
        if (length(display_name) == 0) display_name <- pathway
      } else {
        # Old flat structure
        genes <- GO[[input$db_select]][[pathway]]
        display_name <- pathway
      }
      
      genes_in_data <- intersect(genes, rownames(mat))
      
      if (length(genes_in_data) >= 3) {
    
        sub_mat <- mat[genes_in_data, , drop = FALSE]
        sub_mat <- sub_mat[, ann_df$sample_name, drop = FALSE]
        
        # Scale the data
        scaled_mat <- t(scale(t(sub_mat)))
        
        ann_colors <- rainbow(length(unique(ann_df$group)))
        names(ann_colors) <- unique(ann_df$group)
        ann_list <- list(group = ann_colors)
        
        ann <- ann_df %>%
          column_to_rownames("sample_name") %>%
          .[colnames(scaled_mat), , drop = FALSE]
        
        plot_title <- if (input$sample_desc != "") {
          paste(input$sample_desc, "-", display_name)
        } else {
          display_name  # Use clean display name instead of pathway ID
        }
        
        
        uid <- digest(paste(pathway, input$sample_desc, Sys.time()), algo = "md5")
        
        local({
          local_mat <- scaled_mat
          local_ann <- ann
          local_title <- plot_title
          local_uid <- uid
          local_sub_mat <- sub_mat  # Store original matrix for download
          
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
                      fontsize = 12)
          })
          
          # Add download handler for this specific plot
          output[[paste0("download_matrix_", local_uid)]] <- downloadHandler(
            filename = function() {
              paste0(gsub("[^a-zA-Z0-9]", "_", local_title), ".csv")
            },
            content = function(file) {
              write.csv(local_sub_mat, file, row.names = TRUE)
            }
          )
        })
        
        current_plots[[uid]] <- list(title = plot_title, pathway = pathway, matrix = sub_mat)
      }
    }
    
    # Update stored plots and plotted pathways
    plot_data(current_plots)
    plotted_pathways(union(previous, new_pathways))
    
    # Reset picker selection
    updatePickerInput(session, "selected_pathways", selected = character(0))
    # Reset plot title input
    updateTextInput(session, "sample_desc", value = "")
  })

  
  # Add this flag to show/hide the preview panel
  output$files_uploaded <- reactive({
    demo_loaded() || (!is.null(input$logcpm_file) && !is.null(input$sampleinfo_file))
  })
  outputOptions(output, "files_uploaded", suspendWhenHidden = FALSE)
  
  # Add preview tables

    output$lcpm_preview <- DT::renderDataTable({
      # This will work for both uploaded files and demo data
      mat <- logcpm_data()
      req(mat)
      
      # Convert matrix to data frame with gene names as first column
      preview_df <- data.frame(
        Gene = rownames(mat),
        mat,
        check.names = FALSE
      )
      
      # Show first 1000 rows to avoid performance issues
      if (nrow(preview_df) > 1000) {
        preview_df <- preview_df[1:1000, ]
      }
      
      DT::datatable(
        preview_df,
        options = list(
          scrollX = TRUE,
          scrollY = "400px",
          pageLength = 10,
          lengthMenu = c(10, 25, 50, 100),
          dom = 'frtip',
          autoWidth = TRUE,
          columnDefs = list(
            list(width = "120px", targets = 0),  # Fixed width for gene names
            list(width = "80px", targets = 1:(ncol(preview_df)-1))  # Fixed width for expression values
          )
        ),
        rownames = FALSE,
        caption = paste("Showing", nrow(preview_df), "genes")
      ) %>%
        DT::formatRound(columns = 2:ncol(preview_df), digits = 3)
    }, server = FALSE)
    
    output$metadata_preview <- DT::renderDataTable({
      # This will work for both uploaded files and demo data
      meta_df <- sample_info()
      req(meta_df)
      
      DT::datatable(
        meta_df,
        options = list(
          scrollX = TRUE,
          scrollY = "400px",
          pageLength = 10,
          lengthMenu = c(10, 25, 50, 100),
          dom = 'frtip',
          autoWidth = TRUE
        ),
        rownames = FALSE,
        caption = paste("Showing", nrow(meta_df), "samples")
      )
    }, server = FALSE)
  
  output$heatmap_outputs <- renderUI({
    all_plots <- plot_data()
    if (length(all_plots) == 0) return(NULL)
    
    tagList(
      lapply(seq_along(all_plots), function(i) {
        uid <- names(all_plots)[i]
        pathway <- all_plots[[i]]$pathway
        custom_title <- all_plots[[i]]$title
        tagList(
          tags$h4(paste("Heatmap:", pathway)),
          plotlyOutput(uid, height = "500px"),
          div(style = "margin: 10px 0;",
              downloadButton(paste0("download_matrix_", uid), "Download matrix as CSV", 
                             class = "btn btn-success btn-sm"),
              actionButton(paste0("remove_plot_", uid), "Remove Plot", 
                           class = "btn btn-danger btn-sm", 
                           style = "margin-left: 10px;")
          ),
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
