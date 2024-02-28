library(shiny)
library(shinyjs)
library(bslib)
library(shinylive)
library(shinycssloaders)
library(DT)

library(standR)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(ggrepel)
library(ggpp)
library(ggh4x)
library(limma)

library(tidyverse)
library(tibble)
set.seed(0)
spe_ruv <- readRDS("data/spe_ruv.rds")
design <- readRDS("data/design.rds")
fit <- readRDS("data/fit.rds")

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


options(pillar.sigfig = 4)



# Define UI ----
ui <- bslib::page_navbar(
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
   
  bslib::nav_panel(
    "Introduction",
    textOutput("introduction"),
    

  ),
  

  bslib::nav_panel(
    "PCA",
    # uiOutput("text_out",
    #          style="width:100px;"), 
    uiOutput("pca",
             style = 'padding: 4px') %>% withSpinner() ,
    
    actionButton("toggle_PCAcustom", "Show/hide options",
                 style = 'display: inline-block; padding: 4px'),
    br(),
    
    div(id = "PCAcustom", style = "display: inline-block;",
        uiOutput("customization") 
    )
   # uiOutput("debug")
  ),
  
  bslib::nav_panel(
    "Table",
    uiOutput("table"),
    downloadButton('downloadTable', "Save table")
  ),
  
  bslib::nav_panel(
    "Volcano",
    uiOutput("volcano") %>% withSpinner()
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
  
# 
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
  
  # Subset data based on user input
  spe_ruv_subset <- eventReactive(reactiveRun(),{
    spe_ruv_subset <- spe_ruv[,grepl(paste(reactiveRun(), collapse = "|"), spe_ruv$anno_type)]
    spe_ruv_subset <- scater::runPCA(spe_ruv_subset)
    return(spe_ruv_subset)
  })
  
  # PCA plot pt 1
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
  
  
  # Contrasts
  contrast <- eventReactive(reactiveRun(),{
    comparisons <- list()
    
    # Find the total number of combinations and iterate through
    for (i in 1:(choose(length(reactiveRun()), 2))){
      comparisons[i] <- 
        # Textual description
        # noquote(
        #   paste0(
        #     # First vs second elements in i-th pair
        #     combn(reactiveRun(), 2, simplify=F)[[i]][1],
        #     "vs",
        #     combn(reactiveRun(), 2, simplify=F)[[i]][2]
        #   )
        # ),
        # "=",
        # Mathematical description
        # this will suffice
        noquote(
          paste0(
            combn(reactiveRun(), 2, simplify=F)[[i]][1],
            "-",
            combn(reactiveRun(), 2, simplify=F)[[i]][2]
          )
        )
      
    }
    
    con <- makeContrasts(
      # Must use as.character()
      contrasts=as.character(unlist(comparisons)),
      levels = colnames(design)
    )
    
    colnames(con) <- sub("-", "_vs_", colnames(con))
    
    return(con)
  })
  
  
  efit <- reactive({
    fit_contrast <- contrasts.fit(fit, contrasts = contrast())
    efit <- eBayes(fit_contrast, robust = TRUE)
    return(efit)
  })
  
  topTabDT <- reactive({
    # fit_contrast <- contrasts.fit(fit, contrasts = contrast())
    # efit <- eBayes(fit_contrast, robust = TRUE)
    

    dt <- topTable(efit(), coef=c(1:ncol(contrast())), n=Inf, p.value=0.05, adjust.method="BH", lfc=1) %>%
      tibble::rownames_to_column(., var="Gene") %>%
      select(!c('ProbeName','GeneID',  'HUGOSymbol', 'ProbeDisplayName', 'Accessions', 'GenomeBuild', 'AnalyteType', 'CodeClass', 'ProbePool',
                'TargetGroup', 'genes_lowCount_overNsamples')) 
    # mutate(across(which(is.numeric))) is not compatible with renderDT which is a 'datatable' object.
    
    # this part cannot be piped together with the previous section. ncol(dt) does not evaluate
    # columns = -c(1:2) does not work
    # keep 4 sigfigs for all columns up to n except for the first two
    dt_sigfigs <- dt %>% datatable() %>%  formatSignif(columns=c(3:ncol(dt)), digits=4)
    return(dt_sigfigs)  
      
  })
  
  
  topTabDF <- reactive({
    # fit_contrast <- contrasts.fit(fit, contrasts = contrast())
    # efit <- eBayes(fit_contrast, robust = TRUE)
    
    
    dt <- topTable(efit(), coef=c(1:ncol(contrast())), n=Inf, p.value=0.05, adjust.method="BH", lfc=1) %>%
      tibble::rownames_to_column(., var="Gene") %>%
      select(!c('ProbeName','GeneID',  'HUGOSymbol', 'ProbeDisplayName', 'Accessions', 'GenomeBuild', 'AnalyteType', 'CodeClass', 'ProbePool',
                'TargetGroup', 'genes_lowCount_overNsamples')) 
    # mutate(across(which(is.numeric))) is not compatible with renderDT which is a 'datatable' object.
    
    # this part cannot be piped together with the previous section. ncol(dt) does not evaluate
    # columns = -c(1:2) does not work
    # keep 4 sigfigs for all columns up to n except for the first two
    # dt_sigfigs <- dt %>% datatable() %>%  formatSignif(columns=c(3:ncol(dt)), digits=4)
    return(dt)  
    
  })
  
  output$table <- renderUI({
    renderDataTable(
      topTabDT()
    )
  })
  
  output$downloadTable <- downloadHandler(
    filename = "output.csv",
    content = function(file) {write.table(topTabDF(), file, sep=",", row.names = F)}
  )
  
  
  plotHeight <- reactive(350* ncol(contrast()))
  
  
  
  output$volcano <- renderUI({
    volcanoTable <- list()
    volcanoDF <- list()
    plots <- list()
    
    # for (i in 1:ncol(contrast())){
    #   volcanoTable[i] <- paste0("hello",i)
    # }
    # 
    # renderText({unlist(volcanoTable)})
    
    for (i in 1:ncol(contrast())) {
      volcanoTable[[i]] <- topTable(efit(), coef= i , n=Inf)
      volcanoDF[[i]] <- data.frame(Target.name = rownames(volcanoTable[[i]]),
                                 cbind(logFC = volcanoTable[[i]]$logFC,
                                       PValue = volcanoTable[[i]]$adj.P.Val))
      volcanoDF[[i]]$de <- "NO"
      volcanoDF[[i]]$de[volcanoDF[[i]]$logFC >=1 & volcanoDF[[i]]$PValue < 0.05] <- "UP"
      volcanoDF[[i]]$de[volcanoDF[[i]]$logFC <= -1 & volcanoDF[[i]]$PValue < 0.05] <- "DN"
      volcanoDF[[i]]$deLab <- NA
      volcanoDF[[i]]$deLab[volcanoDF[[i]]$de!="NO"] <- volcanoDF[[i]]$Target.name[volcanoDF[[i]]$de != "NO"]
      
     
      
      # Because i is used inside ggplot function, and for loops have no separate variable scope
      # local() is used
      plots[[i]] <- local({
        # makes i a local variable for ggplot
        i <- i
        ggplot(data=volcanoDF[[i]],
                               aes(x=logFC,
                                 y=-log10(PValue),
                                 # x=sample(seq),
                                 # y=logFC,
                                 col=de,
                                 label=deLab))+
        geom_point()+
        theme_bw()+
        theme(
          axis.ticks = element_line(colour="black"),
          panel.border = element_rect(colour="black"),
          text=element_text(size=14, color="black"),
          axis.text=element_text(size=14, color="black"),
          plot.margin=unit(c(1,1,1,1),"cm"),
          plot.background=element_rect(fill="transparent", colour=NA),
          panel.background=element_rect(fill="transparent", colour=NA),
          panel.grid = element_blank(),
          legend.background=element_rect(fill="transparent", colour=NA),
          legend.box.background=element_rect(fill="transparent", colour=NA),
          legend.key=element_rect(fill="transparent", colour=NA),
          legend.position="none")+
        geom_vline(xintercept=c(-1,1), col="black",
                   linetype=3)+
        geom_hline(yintercept=-log10(0.05), col="black",
                   linetype=3)+
        # geom_text_repel(
        #   segment.colour="black",
        #   segment.linetype=3,
        #   data=volcanoTable[[i]] %>%
        #     filter(logFC>=1|logFC<=-1),
        #   aes(label=deLab),
        #   #position=position_nudge_center(direction="radial",x=3,y=3, center_x=0, center_y=5),
        #   min.segment.length=0,
        #   max.overlaps=4)+
        scale_color_manual(values=c("darkblue", "grey75", "red"),
                           breaks=c("DN","NA","UP"))+
        xlab(expression(log[2]~fold~change))+
        ylab(expression(-log[10]~italic(P)~value))+
        # scale_y_continuous(
        #   labels = scales::label_number(accuracy = 0.1),
        #   breaks = function(x) pretty(floor(seq(0, max(x+1)*1.5))))+
        theme(aspect.ratio=1, rect=element_rect(fill="transparent"))+
        ggtitle(colnames(contrast())[i])
      })
      
    }
    
    # renderPlot({
    #   for (i in 1:ncol(contrast())){
    #    volcanoPlot[[i]]
    #   }
    # })
    # renderPlot({
    #   
    #   volcanoPlot
    #   
    #   
    # })
    
    # Plots as a list can be grouped into a grid and output
    renderPlot(
      #grid.arrange(grobs = plots)
      plot_grid(plotlist = plots, ncol=1, align="v"),
      height = plotHeight()
    )
    

  
  })
  
}
  
  


# Run the app ----
shinyApp(ui = ui, server = server)