library(shiny)
library(shinyjs)
library(shinyWidgets)
library(DT)

ui <- fluidPage(
  useShinyjs(),
  
  # Load fonts & custom styles
  tags$head(
    tags$link(rel = "stylesheet", href = "https://fonts.googleapis.com/css2?family=Roboto:wght@300;400;700&display=swap"),
    tags$style(HTML("
      body {
        font-family: 'Roboto', sans-serif;
        background-color: #f8f9fa;
        color: #323b3f;
        margin: 0;
        padding: 20px;
      }
      .app-header {
        background: linear-gradient(135deg, #162c55, #BA5216);
        color: white;
        padding: 30px 20px;
        position: relative;
        margin-bottom: 20px;
        border-radius: 8px;
      }
      .header-content {
  display: flex;
  justify-content: center;
  align-items: center;
  position: relative;
}
.header-text {
  flex: 1;
  text-align: center;
}

      .header-title {
        font-size: 32px;
        font-weight: bold;
        margin: 0 0 10px 0;
      }
      .header-subtitle {
        font-size: 18px;
        opacity: 0.9;
        margin: 0;
        line-height: 1.4;
      }
      .header-org {
        font-size: 16px;
        font-weight: 500;
        margin: 15px 0 0 0;
      }
.header-logo {
  position: absolute;
  right: 0;
  margin-left: 0;
}
      .header-logo img {
        height: 100px;
        width: auto;
        padding: 15px;
        background: rgba(255, 255, 255, 0.7);
        border-radius: 50%;
        backdrop-filter: blur(10px);
        border: 2px solid rgba(255, 255, 255, 0.3);
        box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2);
      }
      .section-header {
        background: linear-gradient(135deg, #162c55, #BA5216);
        color: white;
        padding: 15px 20px;
        margin: 0 0 20px 0;
        border-radius: 8px;
        font-size: 20px;
        font-weight: bold;
        text-align: center;
      }
      .wellPanel {
        background-color: white;
        border: 1px solid #dee2e6;
        border-radius: 8px;
        padding: 20px;
        margin-bottom: 20px;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
      }
      .btn-primary {
        background-color: #162c55 !important;
        border-color: #162c55 !important;
        color: white !important;
      }
      .btn-primary:hover {
        background-color: #0f1f3d !important;
        border-color: #0f1f3d !important;
      }
      .nav-tabs > li > a {
        color: #162c55;
        font-weight: bold;
      }
      .nav-tabs > li.active > a {
        background-color: #e8f4fd;
        border-color: #162c55;
        color: #162c55;
      }
      .demo-section {
        background-color: #162c55;
        border-radius: 8px;
        padding: 20px;
        text-align: center;
        margin-bottom: 25px;
      }
      .demo-header {
        display: flex;
        align-items: center;
        justify-content: center;
        margin-bottom: 10px;
      }
      .demo-title {
        color: white;
        margin: 0;
        font-weight: bold;
        font-size: 18px;
      }
      .demo-description {
        color: white;
        margin: 5px 0 15px 0;
        font-size: 14px;
      }
      .section-title {
        color: #162c55;
        margin-top: 0;
        margin-bottom: 20px;
        font-size: 24px;
        font-weight: bold;
      }
      .control-section {
        margin-bottom: 25px;
      }
      .control-label {
        font-weight: 600;
        color: #333;
        margin-bottom: 8px;
      }
      .main-content {
        /* Removed max-width constraint to allow full screen usage */
        margin: 0;
        padding: 0;
      }
    "))
  ),
  
  # Top header with gradient background and logo
  div(class = "app-header",
      div(class = "header-content",
          div(class = "header-text",
              h1("HeatPath", class = "header-title"),
              p("Generate pathway gene expression heatmaps from counts by searching pathway databases", 
                class = "header-subtitle"),
              p("UF Health Cancer Center Biostatistics and Computational Biology Shared Resource Bioinformatics Unit", 
                class = "header-org")
          ),
          div(class = "header-logo",
              img(src = "logo.png", alt = "UF Health Logo")
          )
      )
  ),
  
  # Main content with tabs - removed width constraints
  div(class = "main-content",
      
      # Tab panel structure
      tabsetPanel(
        tabPanel("Pathway Heatmaps",
                 
                 # Main section header
                 wellPanel(
                   style = "text-align: center; background: linear-gradient(135deg, #162c55, #BA5216); color: white; border: none; margin-bottom: 30px;",
                   h2("Create Pathway Heatmaps", style = "margin: 15px 0; font-weight: bold;"),
                   p("Upload your logCPM data and generate interactive heatmaps", style = "font-size: 16px; margin: 10px 0;")
                 ),
                 
                 # Two-column layout
                 fluidRow(
                   # Left column - Upload controls
                   column(4,
                          wellPanel(
                            style = "background-color: white; border: 1px solid #dee2e6;",
                            h4("Upload Files", class = "section-title"),
                            
                            # Demo data section
                            div(class = "demo-section",
                                div(class = "demo-header",
                                    h5("Try HeatPath with Demo Data", class = "demo-title"),
                                    actionButton("info_demo", "", icon = icon("info-circle"),
                                                 class = "btn-sm",
                                                 style = "margin-left: 10px; padding: 2px 6px; background: rgba(255,255,255,0.2); border: 1px solid rgba(255,255,255,0.3); color: white;")
                                ),
                                p("Load example RNA-seq data instantly", class = "demo-description"),
                                actionButton("load_demo", "Load Demo Data",
                                             class = "btn-light",
                                             style = "font-weight: bold; color: #162c55; border: 2px solid white;",
                                             icon = icon("upload"))
                            ),
                            
                            # File upload sections
                            div(class = "control-section",
                                div(class = "control-label", "Normalized Count Matrix (logCPM)"),
                                fileInput("logcpm_file", NULL, accept = ".csv", width = "100%")
                            ),
                            
                            div(class = "control-section",
                                div(class = "control-label", "Sample Metadata"),
                                fileInput("sampleinfo_file", NULL, accept = ".csv", width = "100%")
                            ),
                            
                            # Pathway selection
                            div(class = "control-section",
                                div(class = "control-label", "Select Database"),
                                selectInput("db_select", NULL, choices = NULL, width = "100%")
                            ),
                            
                            div(class = "control-section",
                                div(class = "control-label", "Select Pathway(s)"),
                                pickerInput(
                                  inputId = "selected_pathways",
                                  label = NULL,
                                  choices = NULL,
                                  multiple = TRUE,
                                  options = list(
                                    `live-search` = TRUE,
                                    `size` = 10,
                                    `selected-text-format` = "count > 3",
                                    `actions-box` = TRUE,
                                    `style` = "btn-outline-primary",
                                    `iconBase` = "fas",
                                    `tickIcon` = "fas fa-check text-primary"
                                  ),
                                  width = "100%"
                                )
                            ),
                            
                            div(class = "control-section",
                                div(class = "control-label", "Plot Title (optional)"),
                                textInput("sample_desc", NULL, value = "", width = "100%")
                            ),
                            
                            actionButton("generate_plot", "Generate Plot", 
                                         class = "btn btn-primary btn-lg",
                                         style = "width: 100%; margin-top: 10px;")
                          )
                   ),
                   
                   # Right column - Data preview
                   column(8,
                          conditionalPanel(
                            condition = "output.files_uploaded",
                            wellPanel(
                              style = "background-color: white; border: 1px solid #dee2e6;",
                              h4("Preview Your Data", class = "section-title"),
                              tabsetPanel(id = "data_preview",
                                          tabPanel("Count Matrix",
                                                   br(),
                                                   DT::dataTableOutput("lcpm_preview")
                                          ),
                                          tabPanel("Metadata",
                                                   br(),
                                                   DT::dataTableOutput("metadata_preview")
                                          )
                              )
                            )
                          )
                   )
                 ),
                 
                 # Heatmap output section
                 wellPanel(
                   style = "background-color: white; border: 1px solid #dee2e6; margin-top: 30px;",
                   h4("Heatmap Output", class = "section-title"),
                   uiOutput("heatmap_outputs")
                 )
        ),
        
        # About tab
        tabPanel("About",
                 wellPanel(
                   h3("What is HeatPath?", style = "color: #162c55;"),
                   p("HeatPath allows users to generate heatmaps of gene expression for selected pathways using RNA-seq data. It supports GO, KEGG, and Reactome gene sets, and provides interactive plotting and annotation options."),
                   
                   h3("Usage Instructions", style = "color: #162c55;"),
                   tags$ol(
                     tags$li("Upload logCPM expression matrix with gene symbols or Ensembl IDs."),
                     tags$li("Upload sample information CSV with 'sample_name' and 'group' columns."),
                     tags$li("Select a pathway database and choose one or more pathways."),
                     tags$li("Optionally add a title for each plot and generate plots interactively."),
                     tags$li("Use download buttons to save matrix data and remove buttons to clean up plots.")
                   ),
                   
                   h3("What data are plotted?", style = "color: #162c55;"),
                   tags$ol(
                     tags$li("For each selected pathway, the app extracts the corresponding genes from your expression matrix and generates a heatmap."),
                     tags$li("To make the heatmaps easier to interpret, the expression values are standardized within each gene: this means converting the logCPM values into Z-scores by subtracting the mean and dividing by the standard deviation across samples."),
                     tags$li("Red indicates higher-than-average expression for that gene in a sample, while blue indicates lower-than-average expression. White indicates expression close to the gene's average.")
                   )
                 )
        )
      )
  ),
  
  # Footer
  tags$footer(
    tags$hr(),
    div(style = "text-align: center; color: #999; padding: 20px;",
        "Contact: hkates@ufl.edu | UF Health Cancer Center Bioinformatics Shared Resource"
    )
  )
)