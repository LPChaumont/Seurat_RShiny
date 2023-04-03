# Load required packages and custom functions
source("global.R")

ui <- dashboardPage(
  dashboardHeader(title = "scRNAseq Analysis"), # header
  dashboardSidebar( # sidebar menu
    sidebarMenu(
      useShinyjs(), # enable shinyjs
      id = "tab",
      menuItem("Home Page", tabName = "home", icon = icon("list")), # home page sidebar tab
      menuItem("scRNAseq Analyzer", tabName = "input", icon = icon("play")), #input sidebar tab
      conditionalPanel(condition = "input.tab == 'input'", # conditional panel for input sidebar tab
                       div(
                         fileInput("file", "Upload File", multiple = FALSE, accept = c(".rds")), # file input
                         actionButton("reset", "Reset", icon = icon("undo"), # reset button
                                      style = "color: #fff; background-color: #dc3545; width: 87.25%"), 
                         actionButton("run", "Run", icon = icon("play"), # run button
                                      style = "color: #fff; background-color: #28a745; width: 87.25%")
                       )
      )
    )
  ),
  dashboardBody( # main content
    tabItems(
      tabItem(tabName = "home", # refers to tab in sidebar
              mainPanel(includeMarkdown("./markdown/home_page.md"), width = 12)
      ),
      tabItem(tabName = "input",
              tabsetPanel(id = "main_tabs",
                          tabPanel("Instructions", # instructions panel
                                   includeMarkdown("./markdown/instructions.md")
                          )
              )
      )
    )
  )
)

server <- function(input, output, session) {
  options(shiny.maxRequestSize = 300*1024^2) # set max upload size (default: 5 MB)
  
  shinyjs::disable("run") # Disable "Run" by default
  
  # enable/disable run button based on file upload
  observe({
    if (is.null(input$file) != TRUE) {
      shinyjs::enable("run")
    } else {
      shinyjs::disable("run")
    }
  })
  
  # reset tabs and disable "Run" button when reset is clicked
  observeEvent( 
    input$reset,{
      shinyjs::reset("file")
      shinyjs::disable("run")
      removeTab("main_tabs", "UMAP") # remove "UMAP" tab
      removeTab("main_tabs", "Gene Expression") # remove "Gene Expression" tab
    }
  )
  
  # create UMAP and feature plot when "Run" button is clicked
  observeEvent( 
    input$run,{
      shinyjs::disable("run") # disable run button
      
      show_modal_spinner(text = "Preparing plots ...") # show loading spinner while plots are generated
      
      # load Seurat object from uploaded file
      obj <- load_seurat_obj(input$file$datapath) 
      if (is.vector(obj)) { 
        # display error message if file is not a valid Seurat object
        showModal(modalDialog( 
          title = "Error with file",
          HTML("<h5>There is an error with the file you uploaded. See below for more details.</h5><br>",
               paste(unlist(obj)), collapse = "<br><br>")
        ))
        shinyjs::enable("run") # re-enable "Run" button
        
      } else { 
        # if the file was loaded successfully create plots and download buttons
        
        # render UMAP plot
        output$umap <- renderPlot({ 
          create_metadata_umap(obj, input$metadata_col)
        })
        
        # render feature plot
        output$featurePlot <- renderPlot({ 
          create_feature_plot(obj, input$gene)
        })
        
        # download button for UMAP
        output$download_umap <- downloadHandler( 
          filename = function(){
            paste0(input$metadata_col, "_UMAP.png")
          },
          content = function(file){
            plot <- create_metadata_umap(obj, input$metadata_col)
            ggsave(filename = file, width = 10, height = 5, type = "cairo")
          }
        )
        
        # download button for Gene Expression
        output$downloadFeaturePlot <- downloadHandler( 
          filename = function(){
            paste0(input$gene, "_feature_plot.png")
          },
          content = function(file){
            plot <- create_feature_plot(obj, input$gene)
            ggsave(filename = file, width = 10, height = 5, type = "cairo")
          }
        )
        
        # insert UMAP tab
        insertTab(
          inputId = "main_tabs",
          tabPanel(
            "UMAP",
            fluidRow(
              column(
                width = 8,
                plotOutput(outputId = "umap"),
                downloadButton("download_umap", "Download UMAP")
              ),
              column(
                width = 4,
                selectizeInput("metadata_col", "Metadata Column", choices = NULL)
              )
            )
          )
        )
        # improves performance and efficiency, see ?selectizeInput
        updateSelectizeInput(session, "metadata_col",choices = colnames(obj@meta.data), server = TRUE)
        
        
        # insert Gene Expression tab
        insertTab(
          inputId = "main_tabs",
          tabPanel(
            "Gene Expression",
            fluidRow(
              column(
                width = 8,
                plotOutput(outputId = "featurePlot"),
                downloadButton("downloadFeaturePlot", "Download Feature Plot")
              ),
              column(
                width = 4,
                selectizeInput("gene", "Genes",
                               choices = NULL)
              )
            )
          )
        )
        # # improves performance and efficiency, see ?selectizeInput
        updateSelectizeInput(session, "gene",choices = rownames(obj), server = TRUE)
        # remove loading spinner
        remove_modal_spinner()
        shinyjs::enable("run")
      }
    }
  )
}

shinyApp(ui, server)
