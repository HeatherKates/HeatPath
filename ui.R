library(shiny)
library(shinyWidgets)

fluidPage(
  titlePanel("Pathway Heatmaps (GO / KEGG)"),
  sidebarLayout(
    sidebarPanel(
      fileInput("logcpm_file", "Upload logCPM Matrix (.csv)", accept = ".csv"),
      fileInput("sampleinfo_file", "Upload Sample Info (.csv)", accept = ".csv"),
      selectInput("db_select", "Select Database", choices = NULL),  # Populated in server
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
