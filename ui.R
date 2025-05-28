library(shiny)
library(shinyjs)
library(shinyWidgets)

ui <- fluidPage(
  useShinyjs(),
  
  # Load fonts & custom styles
  tags$head(
    tags$link(rel = "stylesheet", href = "https://fonts.googleapis.com/css2?family=Roboto:wght@300;400;700&display=swap"),
    tags$style(HTML("
      body {
        font-family: 'Roboto', sans-serif;
        background-color: white;
        color: #323b3f;
      }
      .app-header {
        background-color: #005496;
        color: white;
        padding: 20px;
        font-size: 22px;
        font-weight: bold;
      }
      .section-header {
        background-color: #005496;
        color: white;
        padding: 10px;
        margin-top: 20px;
        font-size: 18px;
        font-weight: bold;
      }
      .btn-primary {
        background-color: #005496 !important;
        border-color: #005496 !important;
        color: white !important;
      }
      .nav-tabs > li > a {
        color: #005496;
        font-weight: bold;
      }
      .nav-tabs > li.active > a {
        background-color: #e5f0fa;
        border-color: #005496;
      }
    "))
  ),
  
  # Top UFHCC-branded header
  div(class = "app-header", "UFHCC BCB-SR: Pathway Heatmaps"),
  
  # Main app content with About tab
  tabsetPanel(
    tabPanel("Pathway Heatmaps",
             div(class = "section-header", "Pathway Heatmap Controls"),
             
             fileInput("logcpm_file", "Upload logCPM Matrix (.csv)", accept = ".csv"),
             fileInput("sampleinfo_file", "Upload Sample Info (.csv)", accept = ".csv"),
             selectInput("db_select", "Select Database", choices = NULL),
             
             pickerInput(
               inputId = "selected_pathways",
               label = "Select Pathway(s)",
               choices = NULL,
               multiple = TRUE,
               options = list(
                 `live-search` = TRUE,
                 `size` = 10,
                 `selected-text-format` = "count > 3",
                 `actions-box` = TRUE,
                 `style` = "btn-primary",
                 `iconBase` = "fas",
                 `tickIcon` = "fas fa-check text-primary"
               )
             ),
             
             textInput("sample_desc", "Plot Title (optional)", value = ""),
             actionButton("generate_plot", "Generate Plot", class = "btn btn-primary btn-block"),
             
             div(class = "section-header", "Heatmap Output"),
             downloadButton("download_all_html", "Download all plots as HTML file", class = "btn btn-primary", style = "margin-top: 10px;"),
             
             uiOutput("heatmap_outputs")
    ),
    
    tabPanel("About",
             h3("What is This Tool?"),
             p("This app allows users to generate heatmaps of gene expression for selected pathways using RNA-seq data. It supports GO, KEGG, and Reactome gene sets, and provides interactive plotting and annotation options."),
             h3("Usage Instructions"),
             tags$ol(
               tags$li("Upload logCPM expression matrix with gene symbols or Ensembl IDs."),
               tags$li("Upload sample information CSV with 'sample_name' and an annotation column. If >1 annotation columns are included, only the first will be used."),
               tags$li("Select a pathway database and choose one or more pathways."),
               tags$li("Optionally add a title for each plot and generate plots interactively."),
               tags$li("You may download the lcpm counts matrix for the genes and samples in the heatmap."),
               tags$li("Use 'Remove Plot' buttons to remove unwanted plots from the display.")
             ),
             h3("What data are plotted?"),
             tags$ol(
               tags$li("For each selected pathway, the app extracts the corresponding genes from your expression matrix and generates a heatmap."),
               
               tags$li("To make the heatmaps easier to interpret, the expression values are standardized within each gene: this means converting the logCPM values into Z-scores by subtracting the mean and dividing by the standard deviation across samples. This allows you to see relative changes in expression for each gene across the samples, regardless of their absolute expression levels."),
               
               tags$li("Red indicates higher-than-average expression for that gene in a sample, while blue indicates lower-than-average expression. White indicates expression close to the gene’s average. This makes it easier to visually identify patterns of up- or down-regulation across samples.")
    ))
  ),
  
  # Footer
  tags$footer(
    tags$hr(),
    div(style = "text-align: center; color: #999;",
        "Contact: hkates@ufl.edu | UF Health Cancer Center Bioinformatics Shared Resource"
    ),
    style = "margin-top: 40px;"
  )
)
