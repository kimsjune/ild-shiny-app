library(renv)
library(shiny)
library(shinymanager)

library(dplyr)
library(shinyjs)
library(bslib)
library(shinylive)
library(shinycssloaders)
library(DT)
library(mathjaxr)
library(colorspace)

library(standR)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(ggrepel)
library(ggpp)
library(ggh4x)
library(limma)
library(statmod)


library(ComplexHeatmap)
library(circlize)

library(tidyverse)
library(tibble)
set.seed(0)
options(warn=-1)
withMathJax()
spe_ruv <- readRDS("data/spe_ruv.rds")
design <- readRDS("data/design.rds")
fit <- readRDS("data/fit.rds")

theme_set(theme_cowplot(font_size=16))

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


credentials <- data.frame(
  user = c("shiny"), # mandatory
  password = c("ild"), # mandatory
#  start = c("2019-04-15"), # optinal (all others)
#  expire = c(NA, "2019-12-31"),
  admin = FALSE,
  comment = "Simple and secure authentification mechanism 
  for single ‘Shiny’ applications.",
  stringsAsFactors = FALSE
)


# Define UI ----
ui <- bslib::page_navbar(
  tags$head(tags$link(rel="shortcut icon", href="favicon.ico/lung.png")),
  
  tags$style(
    "
    .introText {
    padding-left: 100px;
    padding-right: 100px;
    }
    "
  ),

  useShinyjs(),
  title = "ILD spatial transcriptomics visualization",
  id = "main",
  fillable = T,
  #titlePanel("ILD spatial transcriptomics visualization"),
  sidebar = sidebar(
    
    selectInput(inputId = "anno_type_select",
                label = "Choose annotation-condition:",
                choice = annotation_condition_list,
                multiple = T # Allow >=1 ROIs to be chosen
    ),
    numericInput(inputId = "lfc",
                 label = "log2 fold change cutoff",
                 value = 1
                 ),
    
    actionButton("run", "Run",
                 style = 'dispaly: inline-block; padding: 4px'),
    helpText(
      tags$p("IPF: idiopathic pulmonary fibrosis"),
      tags$p("NSIP: non-specific interstitial pneumonia"),
      tags$p("CHP: chronic hypersensitivity pneumonitis"),
      tags$p("UNC: unclassified"),
      tags$p("NOR:", em("bona fide"), "normal"),
      hr(),
      tags$p("neutral: uninvolved"),
      tags$p("fibroblast: fibroblastic foci")
    ),
    hr(),
    tags$div(
      tags$p(
        a(shiny::icon("github"), " ",
          style = "padding: 10px; text-decoration: none;",
          href = "https://github.com/rstudio/shiny")
      )
    )

    


    
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
    id = "introPage",
    value = "introPage",
    #htmlOutput("introduction"),
    # tags$div(class = "row",
            tags$div(class = "introText",
                     tags$h3("Background"),
                     p("Fibrosing interstitial lung diseases (ILDs) constitute a diverse group of scarring disorders of the lungs that cause progressive respiratory failure. The most common fibrosing ILDs are idiopathic pulmonary fibrosis (IPF), chronic hypersensitivity pneumonitis (CHP) and non-specific interstitial pneumonia (NSIP). Despite some therapeutic advances, fibrosing ILDs are the leading indication for lung transplantation worldwide. Furthermore, the inherent diversity of fibrosing ILDs poses a challenge to diagnostic agreement and clinical management. Thus, there is an urgent need to dissect their molecular underpinnings to develop new diagnostic tools and better match patients with specific treatments."),
                     p("Recently, gene expression profiles have been leveraged to delineate these diseases from each other, and to gain insight into the mechanism of disease progression. However, the current literature falls short in several important aspects. First, traditional bulk tissue RNA-seq data cannot depict the complexity of IPF and CHP that are characterized by regional, temporal and cellular heterogeneity. Indeed, recent single cell RNA-seq (scRNA-seq) studies have shown that diverse cell types are found in the lung. Second, the use of explant lung samples from advanced disease is often overlooked. Deregulated gene expression profiles in end-stage disease are unlikely to be actionable targets to reverse disease progression, especially for fibrosis. Lastly, there is a relative paucity of gene expression data for CHP and NSIP compared to IPF.")
            ),
            tags$div(class = "introText",
                     tags$h3("Scope"),
                     p("This app allows users to interactively visualize spatial transcriptomics data. Select any number of annotation-ILD subtype from the left bar, and hit run. Four output types are available with a few customization options: PCA, table, Volcano and Heatmap."),
                     tags$h3("Rationale"),
                     p("This app was created to facilitate open access to our data and biological interpretation. This is more challenging for spatial transcriptomics data compared to traditional bulk RNA-seq. In bulk RNA-seq, there are generally fewer groups and logical pairwise comparisons. For example, there might be 2-3 treatment groups across one or two categorical variables such as genotype. The number of comparisons is limited, and the full breadth of analysis can be conveyed in a manuscript. On the other hand, there are far more biological groups, and subsequently more comparisons, that can be made in spatial transcriptomic data. For example, there are 20 biological groups in our data. It is often reasonable to choose more than two groups per analysis, such as comparing all 4 annotations within one ILD subtype. Theoretically, there are 1 048 555 combinations of k choices (where k = 2~20) in total. Although most of these may not be biologically meaningful, even their subset cannot be represented easily in a manuscript. With this app, clinicians and scientists can easily explore the data and draw conclusions."),
                    # withMathJax(p(style ="inline-block",  "\\(\\displaystyle\\sum_{k=2}^{20}\\binom{20}{k}\\)"))
                    ),
    #       ),
    # tags$div(class = "row",
             tags$div(class = "introText",
                      # tags$div(
                        tags$h3("Conditions and annotations"),
                      tags$div(style = "text-align: center;",
                        tags$img(src = "images/alluvial_plot.svg", alt = "alluvial")
                      )
                              
                      ),
             tags$div(class = "introText",
                      # tags$div(
                        tags$h3("Citation"),
                        p("If you use this resource, please cite ", 
                          a(href="TBD.com", "Kim et. al. 2024")),
                              # ),
                      # tags$div(
                      #  tags$h3("Links"),
                      
                        # p(
                        #   a(shiny::icon("github"), " ", 
                        #     style = "padding: 10px; text-decoration: none;",
                        #     href = "https://github.com/rstudio/shiny")
                          # a(shiny::icon("linkedin-in")," ",  
                          #   style = "padding: 10px; text-decoration: none;",
                          #   href= "https://www.linkedin.com/in/joon-kim-7a140b90/")
                       #   )
                      )
                      
               
             
      
    

    
    
  ),
  
  
  # bslib::nav_panel(
  #   "PCA",
  #   # uiOutput("text_out",
  #   #          style="width:100px;"),
  #   uiOutput("pca",
  #            style = 'padding: 4px') %>% withSpinner() ,
  #   
  #   actionButton("toggle_PCAcustom", "Show/hide options",
  #                style = 'display: inline-block; padding: 4px'),
  #   br(),
  #   
  #   div(id = "PCAcustom", style = "display: inline-block;",
  #       uiOutput("customization")
  #   )
  #   # uiOutput("debug")
  # ),
  
  # bslib::nav_panel(
  #   "Table",
  #   uiOutput("table"),
  #   downloadButton('downloadTable', "Save table")
  # ),
  
  # bslib::nav_panel(
  #   "Volcano", 
  #   column(3,
  #          div(
  #     sliderInput(inputId = "maxOverlap",
  #                 "Density of gene names to show",
  #                 min = 0, max= 20, value =8),
  #     numericInput(inputId = "delabSize",
  #                  "Gene label size", value = 6),
  #     checkboxInput(inputId = "toggle_customRange",
  #                   "Customize x & y range"
  #                   ),
  #     style = "display: inline-block;"
  #   ))
  #   ,
  #   div(id = "show_customRange",
  #       uiOutput("customRange"), style = "display: inline-block;"
  #     
  #   ),
  # 
  #   uiOutput("volcano") %>% withSpinner()
  # ),
    bslib::nav_panel(
    # "Viz.",
    # navset_card_tab(
      #nav_panel(
        "PCA",
        value = "PCA",
        uiOutput("text_out",
                 style="width:100px;"),
        uiOutput("pca",
                 style = "padding: 4px;") %>% withSpinner(type=4) ,
        
        actionButton("toggle_PCAcustom", "Show/hide options",
                     style = "display: inline-block; padding: 4px;"),
        br(),
        
        div(id = "PCAcustom", style = "display: inline-block;",
            uiOutput("customization")
        )
        
      ),
  
      
      nav_panel(
        "Table",
        uiOutput("table") %>% withSpinner(type=4),
        downloadButton('downloadTable', "Save table")
      ),
      
      nav_panel(
        "Volcano", 
        layout_sidebar(

        #  div(
        #    sliderInput(inputId = "maxOverlap",
        #                "Density of gene names to show",
        #                min = 0, max= 20, value =8),
        #    numericInput(inputId = "delabSize",
        #                 "Gene label size", value = 6),
        #    checkboxInput(inputId = "toggle_customRange",
        #                  "Customize x & y range"
        #    ),
        #    style = "display: inline-block;"
        #  ),
        #  
        # 
        # 
        # div(id = "show_customRange", style = "display: inline-block;",
        #     uiOutput("customRange")
        #     
        # ),
        sidebar =     accordion(
          
          accordion_panel(
            "Volcano options",
            div(
              sliderInput(inputId = "maxOverlap",
                          "Density of gene names to show",
                          min = 0, max= 20, value =8),
              numericInput(inputId = "delabSize",
                           "Gene label size", value = 6),
              checkboxInput(inputId = "toggle_customRange",
                            "Customize x & y range"
              ),
              style = "display: inline-block;"
            ),
            
            
            
            div(id = "show_customRange", style = "display: inline-block;",
                uiOutput("customRange")
                
            )
          ),
          
          # accordion_panel(
          #   "Download",
          #   # radioButtons(inputId = "downloadVolcanoType",
          #   #              label = "format",
          #   #              inline = T,
          #   #              choices = c("png" = ".png",
          #   #                          "svg" = ".svg",
          #   #                          "tiff" = ".tiff",
          #   #                          "pdf" = ".pdf")),
          #   downloadButton("downloadVolcano",
          #                  "Download plots")
          # )
        ),
        
        uiOutput("volcanoUI") %>% withSpinner(type=4),
        downloadButton("downloadVolcano")
      )
      ),
      
      nav_panel(
        "Heatmap",
        layout_sidebar(
          sidebar = 
            accordion(
              accordion_panel(
                "Heatmap options",
                div(

                  numericInput(inputId = "top_n_genes",
                               "Number of genes to show", value = 50),
                  # numericInput(inputId = "heatmap_fc",
                  #              "log2 FC to test", value = 1),
                  textInput("heatmap_col", "Heatmap colour scheme\n (see appendix for options)",
                            value = "Inferno"),

                  sliderInput("heatmap_range",
                              "Z score range", value = c(-2,2), min=-5, max=5),
                  radioButtons("heatmap_size",
                                     "Plot display size",
                                     choices = c(
                                       "x-small" = 9,
                                       "small" = 12, 
                                       "medium" = 16,
                                       "large" = 20),
                               selected = 12),
                  numericInput(inputId = "heatmap_fontsize",
                               "Gene name size", value = 12)


                  
                )
              )
            ),
          uiOutput("heatmapUI") %>% withSpinner(type=4),
          downloadButton("downloadHeatmap")
        )
      ),
  
  
  bslib::nav_panel(
    "Appendix",
    div(
    tags$p("These are all possible colours schmes for the heatmap. Enter the names on the top exactly without quotes."),
    tags$img(src = "images/hcl.svg", alt = "hcl_palette")
    )
  )
  
  # bslib::nav_panel(
  #   "Credit",
  #   div(
  #     tags$p("Amin Manji and Meggie Vo helped with beta testing, and grammar, respectively."),
  #     tags$p("Funded by the AMOSO foundation and PSI),
  #     tags$img(src = "images/amoso-logo.png", alt = "amoso logo")
  #   )
  # )
  
  
  # nav_spacer(),
  # 
  # nav_menu(
  #   title = "Links",
  #   nav_item(
  #     tags$a(
  #       shiny::icon("github"), "", href = "https://github.com/rstudio/shiny"))
  # )
  
  
  
  
  
  
)





ui <- secure_app(ui)


# Define server logic ----
server <- function(input, output, session) {
  observeEvent(input$run,{
    print(paste(input$main))

    
  })
  observeEvent(input$run,{
    if (input$main == "introPage"){
      nav_select(id = "main", selected = "PCA")
    }
  })
  
  
  
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
      # ROIs_text <- c()
      # grepl_text <- c()
      # paste0(ROIs(), collapse = "|")
      # colnames(lcpm_subset_scale())
      # head(colData(spe_ruv[,grepl(paste(reactiveRun(), collapse = "|"), spe_ruv$anno_type)]))
      
      # 
     #  head(colData(spe_ruv[,grepl(paste(reactiveRun(), collapse = "|"), spe_ruv$anno_type)])[colnames(lcpm_subset_scale()), "anno_type"])
     # colnames(lcpm_subset_scale_2())
      
    })
    # renderTable(
    #   colData(spe_ruv[,grepl(paste(reactiveRun(), collapse = "|"), spe_ruv$anno_type)])[colnames()]
    # )
  })
  
  # observeEvent(reactiveRun(),{
  #   #print(c(unlist(reactiveRun())))
  #   # print(str(unlist(reactiveRun())))
  #   # print(unlist(input$anno_type_select))
  #  # print(topTabDF() %>% slice_head(n=50))
  #   #print(unlist(reactiveRun()))
  #  # print(length(reactiveRun()))
  # })
  
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
        # use factor() to disable alphabetical reordering in the legend
        geom_point(aes(shape=factor(anno_type, levels = unlist(reactiveRun())), 
                       fill=factor(anno_type, levels= unlist(reactiveRun()))), size=3, colour="black", stroke=0.5)+
        
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
              axis.text = element_text(color="black",size=16),
              axis.line = element_blank(),
              axis.ticks = element_line(colour="black"),
              axis.title = element_text(size=16),
              #legend.position="none",
              legend.title=element_text(size=16, vjust=0.5, hjust=0.5, face='bold', family='sans'),
              #legend.title=element_blank(),
              legend.text=element_text(size=16, vjust=0.5, hjust=0, face='bold', family='sans'),
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
  
  lfc <- reactive({
    input$lfc
  })
  
  ## Redundant, so deprecated
  # topTabDT <- reactive({
  #   # fit_contrast <- contrasts.fit(fit, contrasts = contrast())
  #   # efit <- eBayes(fit_contrast, robust = TRUE)
  #   dt <- topTable(efit(), coef=c(1:ncol(contrast())), n=Inf, p.value=0.05, sort.by = "P", adjust.method="BH", lfc=lfc()) %>%
  #     tibble::rownames_to_column(., var="Gene") %>%
  #     select(!c('ProbeName','GeneID', 'HUGOSymbol', 'ProbeDisplayName', 'Accessions', 'GenomeBuild', 'AnalyteType', 'CodeClass', 'ProbePool','TargetGroup', 'genes_lowCount_overNsamples'))
  #   # mutate(across(which(is.numeric))) is not compatible with renderDT which is a 'datatable' object.
  #   
  #   # this part cannot be piped together with the previous section. ncol(dt) does not evaluate
  #   # columns = -c(1:2) does not work
  #   # keep 4 sigfigs for all columns up to n except for the first two
  #   dt_sigfigs <- dt %>% datatable() %>%  formatSignif(columns=c(3:ncol(dt)), digits=4)
  #   return(dt_sigfigs)
  #   
  # })
  
  
  topTabDF <- reactive({
    # If there are more than two groups, must sort by F, not P...
    # but it might sort by default? could be redundant

    if (length(reactiveRun()) > 2 ) {
      dt <- topTable(efit(), coef=c(1:ncol(contrast())), n=Inf, p.value=0.05, sort.by = "F", adjust.method="BH", lfc=lfc()) %>%
        tibble::rownames_to_column(., var="Gene") %>%
        select(!c('ProbeName','GeneID',  'HUGOSymbol', 'ProbeDisplayName', 'Accessions', 'GenomeBuild', 'AnalyteType', 'CodeClass', 'ProbePool', 'TargetGroup', 'genes_lowCount_overNsamples'))
    # mutate(across(which(is.numeric))) is not compatible with renderDT which is a 'datatable' object.
    
    # this part cannot be piped together with the previous section. ncol(dt) does not evaluate
    # columns = -c(1:2) does not work
    # keep 4 sigfigs for all columns up to n except for the first two
    # dt_sigfigs <- dt %>% datatable() %>%  formatSignif(columns=c(3:ncol(dt)), digits=4)
    }
    else {
      dt <- topTable(efit(), coef=c(1:ncol(contrast())), n=Inf, p.value=0.05, sort.by = "F", adjust.method="BH", lfc=lfc()) %>%
        tibble::rownames_to_column(., var="Gene") %>%
        select(!c('ProbeName','GeneID',  'HUGOSymbol', 'ProbeDisplayName', 'Accessions', 'GenomeBuild', 'AnalyteType', 'CodeClass', 'ProbePool', 'TargetGroup', 'genes_lowCount_overNsamples'))
    }
    return(dt)
    
  })
  
  output$table <- renderUI({
    renderDataTable(
      # No need to create another reactive expression for the dataTABLE version of dataframe
      topTabDF() %>% datatable() %>%  formatSignif(columns=c(3:ncol(topTabDF())), digits=4)
    )
  })
  
  output$downloadTable <- downloadHandler(
    filename = "output.csv",
    content = function(file) {write.table(topTabDF(), file, sep=",", row.names = F)}
  )
  
  
  plotHeight <- reactive(350* ncol(contrast()))
  
  maxOverlap <- reactive(input$maxOverlap)
  
  customX <- reactive(input$customX)
  
  customY <- reactive(input$customY)
  
  delabSize <- reactive(input$delabSize)
  
  output$customRange <- renderUI({
    div(
      sliderInput(inputId = "customX",
                  "logFC", min = -10, max = 10, dragRange = F, value = c(-4,4)
      ),
      numericInput(inputId = "customY",
                "-log P value", value = 10)
    )
  })
  
  
  observeEvent(input$toggle_customRange, {
    shinyjs::toggle("show_customRange")
  })
  

  # Create volcano plots as a reactive expression
  # Used to renderUI/Plot AND downloadHandler
  # renderUI output cannot be used to download because it becomes HTML?
  volcano <- reactive({
    volcanoTable <- list()
    volcanoDF <- list()
    plots <- list()
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
    

    # If not using custom range, determine maximum absolute FC for each plot
        if (input$toggle_customRange == F) {
  
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
                aspect.ratio=1, rect=element_rect(fill="transparent"),
                axis.ticks = element_line(colour="black"),
                panel.border = element_rect(colour="black"),
                text=element_text(size=20, color="black"),
                axis.text=element_text(size=20, color="black"),
                axis.title = element_text(size=20, color="black"),
                plot.margin=unit(c(1,1,1,1),"mm"),
                plot.background=element_rect(fill="transparent", colour=NA),
                plot.title = element_text(size =20),
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
              geom_text_repel(
                size = delabSize(),
                segment.colour="black",
                segment.linetype=3,
                data=volcanoDF[[i]] %>%
                  filter(logFC>=1|logFC<=-1),
                aes(label=deLab),
                #position=position_nudge_center(direction="radial",x=3,y=3, center_x=0, center_y=5),
                min.segment.length=0,
                max.overlaps = maxOverlap())+
              scale_color_manual(values=c("darkblue", "grey75", "red"),
                                 breaks=c("DN","NA","UP"))+
              xlab(expression(log[2]~fold~change))+
              ylab(expression(-log[10]~italic(P)~value))+
              # scale_y_continuous(
              #   labels = scales::label_number(accuracy = 0.1),
              #   breaks = function(x) pretty(floor(seq(0, max(x+1)*1.5))))+
              
              labs(title = colnames(contrast())[i])+
              xlim(
                0- max(abs(volcanoTable[[i]] %>% select(logFC))),
                max(abs(volcanoTable[[i]] %>% select(logFC)))
              )
          })
          
        
      } else {
          
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
                text=element_text(size=16, color="black"),
                title = element_text(size=16, hjust=0.5),
                axis.text=element_text(size=16, color="black"),
                plot.margin=unit(c(1,1,1,1),"mm"),
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
              geom_text_repel(
                size = delabSize(),
                segment.colour="black",
                segment.linetype=3,
                data=volcanoDF[[i]] %>%
                  filter(logFC>=1|logFC<=-1),
                aes(label=deLab),
                #position=position_nudge_center(direction="radial",x=3,y=3, center_x=0, center_y=5),
                min.segment.length=0,
                max.overlaps = maxOverlap())+
              scale_color_manual(values=c("darkblue", "grey75", "red"),
                                 breaks=c("DN","NA","UP"))+
              xlab(expression(log[2]~fold~change))+
              ylab(expression(-log[10]~italic(P)~value))+
              # scale_y_continuous(
              #   labels = scales::label_number(accuracy = 0.1),
              #   breaks = function(x) pretty(floor(seq(0, max(x+1)*1.5))))+
              theme(aspect.ratio=1, rect=element_rect(fill="transparent"))+
              labs(title = colnames(contrast())[i])+
              xlim(
                customX()[1],
                customX()[2]
              )+
              ylim(
                0, customY()
              )
          })
      }
    }
    
    return(plots)
  })
  
  volcanoPlots <- reactive({
    return(plot_grid(plotlist = volcano(), ncol=1, align="v", axis="lr"))
    
  })
  
  output$volcanoUI <- renderUI({
      # Plots as a list can be grouped into a grid and output
      renderPlot(
        #grid.arrange(grobs = plots)
        plot_grid(plotlist = volcano(), ncol=1, align="v",  axis ='lr'),
        height = plotHeight()
      )
  })
  
  output$downloadVolcano <- downloadHandler(
    filename = function() {
      paste0("volcano_",Sys.Date(),".png")
    },
    content = function(file) { 
      save_plot(volcanoPlots(),filename = file)
    },
    contentType = "image/png"

  )
  
  

  
  top_n_genes <- reactive({
    input$top_n_genes
  })
  
  heatmap_col <- reactive({
    input$heatmap_col
  })
  
  heatmap_range <- reactive({
    input$heatmap_range
  })
  
  heatmap_size <- reactive({
    input$heatmap_size
  })
  
  heatmap_fontsize <-reactive({
    input$heatmap_fontsize
  })


  
  
  lcpm_subset_scale <- reactive({
    mydata <- list()
    for (i in seq_along(reactiveRun())) {
      mydata[[i]] <- assay(spe_ruv_subset(),2)[, colData(spe_ruv_subset())$anno_type==reactiveRun()[i]]
    }


    lcpm_subset_scale <- t(scale(t(data.frame(mydata))))

    return(lcpm_subset_scale)
  })
  
  colnames4heatmap <- reactive({
    mydata <- list()
    for (i in seq_along(reactiveRun())) {
      mydata[[i]] <- assay(spe_ruv_subset(),2)[, colData(spe_ruv_subset())$anno_type==reactiveRun()[i]]
    }

    
    return(colnames(do.call(cbind,mydata)))
  })
  

  
  
  lcpm_subset_scale_topGenes <- reactive({
    ## BEWARE! top_n() reorders rows by some column value. Must use slice_head() to pick first n rows
    lcpm_subset_scale_topGenes <- lcpm_subset_scale()[topTabDF()  %>% slice_head(n=top_n_genes()) %>% select(Gene) %>% unlist %>% unname,]
    
    return(lcpm_subset_scale_topGenes)
  })
  
  
  heatmap <- reactive({
    col_fun <- colorRamp2(c(heatmap_range()[1],0,heatmap_range()[2]), hcl_palette= heatmap_col())
    
    chm <- Heatmap(lcpm_subset_scale_topGenes(),
                   cluster_columns = F,
                   col = col_fun,
                   # width = unit(dim(lcpm_subset_scale_topGenes())[2]*15, "mm"),
                   # height = unit(dim(lcpm_subset_scale_topGenes())[1]*15, "mm"),
                   width = unit(as.numeric(heatmap_size())/2 * dim(lcpm_subset_scale_topGenes())[2],"mm"),
                   height = unit(as.numeric(heatmap_size())/2 * dim(lcpm_subset_scale_topGenes())[1],"mm"),
                   heatmap_legend_param = list(border = "black",
                                               title = "Z score",
                                               title_gp = gpar(fontsize= heatmap_fontsize(),fontface='plain', fontfamily='sans'),
                                               labels_gp = gpar(fontsize= heatmap_fontsize(), fontface='plain', fontfamily='sans'),
                                               legend_height = unit(3 * as.numeric(heatmap_size()), "mm")),
                   top_annotation = HeatmapAnnotation(
                     foo= anno_block(gp = gpar(lty=0, fill="transparent"), 
                                     labels = unlist(reactiveRun()),
                                     labels_gp = gpar(col="black", fontsize=14, fontfamily='sans',fontface='bold'),
                                     labels_rot=0, labels_just = "center", labels_offset = unit(4.5,"mm"))
                   ),
                   
                   border_gp =  gpar(col="black", lwd=0.2),
                   row_names_gp = gpar(fontfamily = 'sans', fontface = 'italic', fontsize = heatmap_fontsize()),
                   show_column_names = F,
         
                   # top_annotation= HeatmapAnnotation(
                   #   foo = anno_block(
                   #     gp = gpar(lty=0, fill="transparent"),
                   #     labels = unlist(reactiveRun()),
                   #     labels_gp = gpar(col="black", fontsize=7, fontfamily = "sans", fontface = "bold"),
                   #     labels_rot = 20, labels_just = "center", labels_offset = unit(4,"mm")
                   #   )
                   # ),
                   column_split = rep(LETTERS[seq_along(reactiveRun())],
                                      # times = as.numeric(unname(table(colData(spe_ruv_subset())[colnames(lcpm_subset_scale()), "anno_type"])))
                                      times = as.numeric(unname(table(colData(spe_ruv_subset())[colnames4heatmap(),"anno_type"])))
                   ),
                   column_title = NULL
    )
    
    return(chm)
  })
  
  output$heatmapUI <- renderUI({
    
    renderPlot(
      heatmap(),
      height = 1.5 * as.numeric(heatmap_size())  * dim(lcpm_subset_scale_topGenes())[1]
      # height = as.numeric(heatmap_size()) * dim(lcpm_subset_scale_topGenes())[1],
      # width = as.numeric(heatmap_size()) * dim(lcpm_subset_scale_topGenes())[2]
      )
      
  })
  
  output$downloadHeatmap <- downloadHandler(
    filename = function() {
      paste0("heatmap_",Sys.Date(),".png")
    },
    content = function(file) {
      png(file, width = dim(lcpm_subset_scale_topGenes())[2]*50, height = dim(lcpm_subset_scale_topGenes())[1]*50, bg="transparent")
      draw(heatmap())
      dev.off()

  
    },
    contentType = "image/png"
  )
  
  # check_credentials returns a function to authenticate users
  res_auth <- secure_server(
    check_credentials = check_credentials(credentials)
  )
  
  output$auth_output <- renderPrint({
    reactiveValuesToList(res_auth)
  })  
  

  
}




# Run the app ----
shinyApp(ui = ui, server = server)