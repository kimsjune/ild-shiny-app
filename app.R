library(shiny)
library(shinyjs)
library(bslib)

library(standR)
library(ggplot2)
library(ggh4x)

spe_ruv <- readRDS("data/spe_ruv.rds")

annotation_condition <- c("fibroblast_IPF","fibrosis_IPF","neutral_IPF",
                          "inflammatory_NSIP","central_NSIP","peripheral_NSIP")

# A named list organizes annotations by name (condition)
annotation_condition_list <- list(
  IPF = c("fibroblast_IPF", "fibrosis_IPF", "neutral_IPF", "lymphoid_IPF"),
  NSIP = c("central_NSIP", "inflammatory_NSIP", "peripheral_NSIP", "airway_NSIP"),
  CHP = c("granuloma_CHP", "inflammatory_CHP", "fibrosis_CHP", "neutral_CHP"),
  UNC = c("fibroblast_UNC", "lymphoid_UNC", "fibrosis_UNC", "neutral_UNC"),
  Normal = c("peripheral_NOR", "central_NOR", "airway_NOR", "pleura_NOR")
)







# Define UI ----
ui <- page_navbar(
  useShinyjs(),
  title = "ILD spatial transcriptomics visualization",
  #titlePanel("ILD spatial transcriptomics visualization"),
  sidebar = sidebar(
    
    selectInput(inputId = "anno_type_select",
                label = "Choose annotation-condition:",
                choice = annotation_condition_list,
                multiple = T # Allow >=1 ROIs to be chosen
                ),
    actionButton("run", "Run",
                 style = 'dispaly: inline-block; padding: 4px'),

    #uiOutput("run", style = 'display: inline-block; padding: 4px')
    # actionButton("run", "Run",
    #              style = 'dispaly: inline-block; padding: 4px'),
    # actionButton("toggle_PCAcustom", "Show/hide options",
    #              style = 'display: inline-block; padding: 4px'),
    # br(),
    # 
    # div(id = "PCAcustom", style = 'padding: 4px',
    #   uiOutput("customization")
    # )
  ),
   
  nav_panel(
    "Introduction",
    textOutput("introduction"),
    

  ),
  

  nav_panel(
    "PCA",
    # uiOutput("text_out",
    #          style="width:100px;"), 
    uiOutput("pca",
             style = 'padding: 4px'),
    
    actionButton("toggle_PCAcustom", "Show/hide options",
                 style = 'display: inline-block; padding: 4px'),
    br(),
    
    div(id = "PCAcustom", style = "display: inline-block;",
        uiOutput("customization")
    )
   # uiOutput("debug")
  ),
  
  nav_panel(
    "Volcano plot"
  )
  
      
    
    
    
    
  )
  
  


# Define server logic ----
server <- function(input, output, session) {
  # some introduction
  output$introduction <- renderText(
    "This app allows users to interactively visualize spatial transcripomics data from Kim et al. 2024"
  )
  
  # This reacts when when it's blank
  ROIs <- reactive({input$anno_type_select}
  )

  
  # output$run <- renderUI({
  #   actionButton("run",
  #                label="Run")
  # })
  
  observeEvent(input$toggle_PCAcustom,{
    shinyjs::toggle(id="PCAcustom")
  })
  
  # Reacts when "run" button is clicked, just to return the chosen ROIs from input$anno_type_select
  # This intermediate is required to wait for the run to be clicked before any processing
  reactiveRun <- eventReactive(input$run, {
    return(ROIs())
  })
  

  
  # Some text output for debugging
  output$text_out <- renderUI({
    renderText({
      ROIs_text <- c()
      grepl_text <- c()
      paste0(ROIs(), collapse = "|")
     
     
     
    })
  })
  
  # Opens a grouped widgets to pick shape and colour
  output$customization <- renderUI({
    # Must be initialized first
    shapes_colours_pca <- list()
    
    # seq_along() counts the total # of elements
    for (i in seq_along(reactiveRun())) {
      shapes_colours_pca[[i]] <- 
        #wellPanel( # wellPanel() groups radio and text nicely
        # but it also breaks inline-block style
        div(
          radioButtons(inputId = paste0("shape_",reactiveRun()[i]),
                       # For each ROI, I need to address its shape and colour separately
                       # ex/ "shape_fibroblast_IPF"
                       
                       label = paste0("Pick a shape for ",reactiveRun()[i]),
                       
                       # Each shape is matched to ggplot2 shapes
                       choices = list(
                         "Circle" = 21,
                         "Square" = 22,
                         "Diamond" = 23,
                         "Triangle_up" = 24,
                         "Triangle_down" = 25)
                       ),
          textInput(inputId = paste0("colour_",reactiveRun()[i]),
                    label = paste0("Pick a colour for ",reactiveRun()[i]),
                    value=sample(c("black","blue","pink3","purple","orange",
                                   "darkgreen","maroon","turquoise3"),1)
                    ),
          style = "display: inline-block;"
        )
   # )
    
    }

    shapes_colours_pca# Show the actual widget,

    })
  

  # debugging use
  output$debug <- renderUI({
    ROIshapes <- list()
    ROIcolours <- list()
    for (i in seq_along(ROIs())) {
      ROIshapes[i] <- input[[paste0("shape_",reactiveRun()[i])]]
      ROIcolours[i] <- input[[paste0("colour_",reactiveRun()[i])]]
    }
    
    renderText({
      as.integer(unlist(shapes))
      print(unlist(ROIcolours))
    })
  })
  
  # PCA plot pt 1
  # Calculate PCA once for any selection of ROIs 
  # <<- required for global scope used in pt 2
  # observeEvent(reactiveRun(),{
  #   # paste() with collapse to generate "reactiveRun()[1] | ..." required to grepl 
  #   # ... annotation-condition from spe_ruv
  #   spe_ruv_subset <<- spe_ruv[,grepl(paste(reactiveRun(), collapse = "|"), spe_ruv$anno_type)]
  #   spe_ruv_subset <<- scater::runPCA(spe_ruv_subset)
  #   pca_ruv_results_subset <<- reducedDim(spe_ruv_subset, "PCA")
  # })
  
  
  spe_ruv_subset <- eventReactive(reactiveRun(),{
    spe_ruv_subset <- spe_ruv[,grepl(paste(reactiveRun(), collapse = "|"), spe_ruv$anno_type)]
    spe_ruv_subset <- scater::runPCA(spe_ruv_subset)
    return(spe_ruv_subset)
  })
  
  pca_ruv_results_subset <- eventReactive(spe_ruv_subset(),{
    return(reducedDim(spe_ruv_subset(), "PCA"))
  })
  
  # PCA plot pt 2
  output$pca <- renderUI({
    # Initialize
    ROIshapes <- list()
    ROIcolours <- list()
    
    for (i in seq_along(reactiveRun())) {
      # Using shape and colour from the customization UI as inputs
      ROIshapes[i] <- input[[paste0("shape_",reactiveRun()[i])]]
      ROIcolours[i] <- input[[paste0("colour_",reactiveRun()[i])]]
    }
    
    # paste() with collapse to generate "reactiveRun()[1] | ..." required to grepl
    # ... annotation-condition from spe_ruv
    # spe_ruv_subset <- spe_ruv[,grepl(paste(reactiveRun(), collapse = "|"), spe_ruv$anno_type)]
    # spe_ruv_subset<- scater::runPCA(spe_ruv_subset)
    # pca_ruv_results_subset<- reducedDim(spe_ruv_subset, "PCA")

      
    renderPlot({
      # withProgress(message="Making plot", value=0, {
      #   n <- 10
      #   
      #   for (i in 1:n) {
      #     incProgress(1/n)
      #   }
      #  
      # })
      
    
      drawPCA(spe_ruv_subset(), precomputed=pca_ruv_results_subset())+
      #geom_point(colour="black", pch=21, size=2.5)+
      geom_point(aes(shape=anno_type, fill=anno_type), size=3, colour="black", stroke=0.5)+
      scale_shape_manual("Annotation-condition",
        values = as.integer(unlist(ROIshapes)) # as.integer() crucial
      )+
      scale_fill_manual("Annotation-condition",
        values = unlist(ROIcolours) # colour() is a base function, so must be avoided
      )+
      scale_y_continuous(
        labels=scales::number_format(accuracy=0.1))+
      scale_x_continuous(
        labels=scales::number_format(accuracy=0.1))+
      theme_bw()+
      theme(panel.grid.minor=element_blank(),
            panel.grid.major=element_blank(),
            axis.text = element_text(color="black",size=14),
            axis.line = element_blank(),
            axis.ticks = element_line(colour="black"),
            axis.title = element_text(size=14),
            #legend.position="none",
            legend.title=element_text(size=14, vjust=0.5, hjust=0.5, face='bold', family='sans'),
            #legend.title=element_blank(),
            legend.text=element_text(size=14, vjust=0.5, hjust=0, face='bold', family='sans'),
            plot.margin=unit(c(1,1,1,1),"mm"),
            plot.background=element_rect(fill="transparent", colour=NA),
            panel.border = element_rect(colour="black", linewidth=0.4),
            panel.background=element_rect(fill="transparent", colour=NA),
            legend.background=element_rect(fill="transparent", colour=NA),
            legend.box.background=element_rect(fill="transparent", colour=NA),
            legend.key=element_rect(fill="transparent", colour=NA),
            legend.position = 'bottom',
            aspect.ratio=1)+
      guides(fill=guide_legend(nrow=2, byrow=T, title.position="top"), shape=guide_legend(nrow=2,byrow=T, title.position="left"))+
        # similar to set_panel_size() from egg(), but this way there is no background grid visible
        force_panelsizes(rows = unit(2.5, "in"),
                         cols = unit(2.5, "in"))
      

      
    })
  })
  
  
  output$volcano <- renderUI({
    
  })
  
}
  
  


# Run the app ----
shinyApp(ui = ui, server = server)