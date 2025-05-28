library(shiny)
library(shinyWidgets)


fluidPage(
  
  titlePanel("Pathway Heatmaps (GO / KEGG)"),
  sidebarLayout(
    sidebarPanel(
      fileInput("logcpm_file", "Upload logCPM Matrix (.csv)", accept = ".csv"),
      fileInput("sampleinfo_file", "Upload Sample Info (.csv)", accept = ".csv"),
      selectInput("db_select", "Select Database", choices = NULL),  # dynamically populated in server
      pickerInput(
        inputId = "selected_pathways",
        label = "Select Pathway(s)",
        choices = NULL,  # Choices will be updated dynamically via updatePickerInput()
        multiple = TRUE,
        options = list(
          `live-search` = TRUE,
          `size` = 10,
          `selected-text-format` = "count > 3",
          `actions-box` = TRUE,
          `style` = "btn-primary"
        )
      )
      ,
      
      textInput("sample_desc", "Plot Title (optional)", value = ""),
      actionButton("generate_plot", "Generate Plot"),
      width = 3
    ),
    mainPanel(
      uiOutput("heatmap_outputs")
    )
  )
)
