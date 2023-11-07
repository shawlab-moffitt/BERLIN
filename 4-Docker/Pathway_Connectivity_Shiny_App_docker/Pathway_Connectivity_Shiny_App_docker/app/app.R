####----Install and load packages----####

packages <- c("shiny","shinyjqui","pheatmap","RColorBrewer",
              "ggdendro","factoextra","dplyr","DT","viridis","readr","ggrepel",
              "shinycssloaders","stringr","tools","plotly","reshape2","ggpubr","gridExtra")

installed_packages <- packages %in% rownames(installed.packages())
#if (any(installed_packages == FALSE)) {
#  install.packages(packages[!installed_packages])
#}
invisible(lapply(packages, library, character.only = TRUE))
#bioconductor packages
bioCpacks <- c("clusterProfiler")
installed_packages_BIOC <- bioCpacks %in% rownames(installed.packages())
#if (any(installed_packages_BIOC == FALSE)) {
#  BiocManager::install(bioCpacks[!installed_packages_BIOC], ask = F)
#}
invisible(lapply(bioCpacks, library, character.only = TRUE))


####---- Read in Files----####  

GeneSet_File <- "GeneSet_List_HS_v5.RData"

PreSet_CoxH_Ranked_Feature_file <- "1-Prep_Input_Files/Output/GSE116256_AMLscRNA_SingleCellCluster_withTCGA_GeneSet.txt"





##--Example Data--##

ExamplePathway_File <- "Pathways_Of_Interest.txt"

ExampleAnnotation_File <- "Gene_Annoation.txt"

# R Data list load function for naming
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
gs <- loadRData(GeneSet_File)
names(gs) <- gsub("[[:punct:]]",".",names(gs))

if (is.character(PreSet_CoxH_Ranked_Feature_file)) {
  if (file.exists(PreSet_CoxH_Ranked_Feature_file)) {
    PreSet_CoxH_Ranked_Features <- as.data.frame(readr::read_delim(PreSet_CoxH_Ranked_Feature_file, delim = '\t', col_names = T, comment = "#"))
  } else { PreSet_CoxH_Ranked_Features <- NULL }
} else { PreSet_CoxH_Ranked_Features <- NULL }


#learn Jaccard index function
jaccard <- function(a, b) {  
  a = unique(a)
  b = unique(b)
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (1-(intersection/union))
}

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}


#increase file upload size
options(shiny.maxRequestSize=5000*1024^2)



####----UI----####

ui <- 
  navbarPage("{ Jaccard Pathway Connectivity Index }",
             
             ####----Pathway Connectivity----####
             
             tabPanel("Pathway Connectivity",
                      fluidPage(
                        title = "Pathway Connectivity",
                        sidebarLayout(
                          sidebarPanel(
                            width = 3,
                            tabsetPanel(
                              
                              ####----Pathway Parameters----####
                              
                              tabPanel("Pathway Parameters",
                                       p(),
                                       h4("Upload Pathways of Interest"),
                                       fluidRow(
                                         column(9,
                                                fileInput("UserPathwayFile","Upload File")
                                         ),
                                         column(3,
                                                checkboxInput("HeaderCheckIntra","Header",value = T)
                                         )
                                       ),
                                       fluidRow(
                                         column(6, style = 'margin-top:-15px;',
                                                uiOutput("rendUseExpData")
                                         ),
                                         column(6, style = 'margin-top:-15px;',
                                                uiOutput("rendDownloadLink")
                                                )
                                       ),
                                       p(),
                                       uiOutput("rendTopFeatureSelect"),
                                       uiOutput("rendHRLine1"),
                                       uiOutput("rendClusterHeader"),
                                       fluidRow(
                                         column(6,
                                                uiOutput("rendClustMethodIntra"),
                                                uiOutput("rendViewClustTabIntra")
                                         ),
                                         column(6,
                                                uiOutput("rendNumClusters"),
                                                uiOutput("renddnldClustTabIntra")
                                         )
                                       ),
                                       uiOutput("rendClustTabIntra"),
                                       uiOutput("rendHRLine2"),
                                       uiOutput("rendSIFheader"),
                                       fluidRow(
                                         column(6,
                                                uiOutput("rendJaccDistCutoff")
                                                ),
                                         column(6,
                                                uiOutput("rendPrevSIF")
                                                )
                                       ),
                                       uiOutput("rendSIFPreview"),
                                       uiOutput("renddnldSIFTabIntra")
                              ),
                              
                              ####----Pathway Figure Parameters----####
                              
                              tabPanel("Figure Parameters",
                                       h4("Heatmap Parameters"),
                                       selectInput("ColorPaletteIntra", "Select Color Palette:",
                                                   choices = c("Red/Blue" = "original",
                                                               "OmniBlueRed" = "OmniBlueRed",
                                                               "LightBlue/BlackRed" = "LightBlueBlackRed",
                                                               "Green/Black/Red" = "GreenBlackRed",
                                                               "Yellow/Green/Blue" = "YlGnBu","Inferno" = "Inferno",
                                                               "Viridis" = "Viridis","Plasma" = "Plasma",
                                                               "Reds" = "OrRd","Blues" = "PuBu","Greens" = "Greens")),
                                       fluidRow(
                                         column(6,
                                                checkboxInput("HeatColNamesIntra","Show Heatmap Column Names", value = F),
                                                numericInput("HeatColFontIntra", "Heatmap Column Font Size:",
                                                             min = 5, max = 75,
                                                             value = 12, step = 1),
                                                numericInput("HeatColDendHeight","Column Dendrogram Height:",
                                                             value = 50, step = 1)
                                         ),
                                         column(6,
                                                checkboxInput("HeatRowNamesIntra","Show Heatmap Row Names", value = F),
                                                numericInput("HeatRowFontIntra", "Heatmap Row Font Size:",
                                                             min = 5, max = 75,
                                                             value = 9, step = 1),
                                                numericInput("HeatRowDendHeight","Row Dendrogram Height:",
                                                             value = 50, step = 1)
                                         )
                                       ),
                                       shiny::hr(),
                                       h4("Connectivity Visualization Parameters"),
                                       selectInput("ConnView","View Connectivity as:",
                                                   choices = c("Phylogeny" = "phylogenic","Dendrogram" = "rectangle","Circular" = "circular")),
                                       uiOutput("rendPhyloLayout"),
                                       fluidRow(
                                         column(6,
                                                checkboxInput("ShowConnLabels","Show Labels:",value = T)
                                         ),
                                         column(6,
                                                numericInput("ConnFontSize","Font Size:",
                                                             value = 0.6, step = 0.1)
                                         )
                                       )
                              )
                            )
                          ),
                          
                          ####----Pathway Main Panel----####
                          
                          mainPanel(
                            tabsetPanel(
                              
                              tabPanel("Jaccard Pathway Connectivity Table",
                                       p(),
                                       div(DT::dataTableOutput("JaccTableIntra"), style = "font-size:12px"),
                                       uiOutput("renddnldJaccTabIntra")
                              ),
                              
                              tabPanel("Heatmap",
                                       shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput('JaccHeatmapIntra', width = "100%", height = "800px")), type = 6),
                                       downloadButton("dnldJaccHeatmapIntraPDF","Download as PDF"),
                                       downloadButton("dnldJaccHeatmapIntraSVG","Download as SVG")
                              ),
                              
                              tabPanel("Clustering",
                                       uiOutput("rendJaccDendo"),
                                       downloadButton("dnldJaccDendoPDF","Download as PDF"),
                                       downloadButton("dnldJaccDendoSVG","Download as SVG")
                              ),
                              
                              tabPanel("Gene Clusters and Annoation",
                                       fluidRow(
                                         column(3,
                                                fileInput("UserAnnotationFile","Upload Annotation File",
                                                          accept = c(".tsv",".txt",".csv"))
                                                ),
                                         column(8,
                                                h2(),
                                                actionButton("useExpAnnoData","Load Example Data"),
                                                tags$a(href="http://shawlab.science/shiny/PATH_SURVEYOR_ExampleData/Pathway_Connectivity_App/GeneAnnoationFile_Upload/", "Download example data", target='_blank'),
                                                )
                                         ),
                                       div(DT::dataTableOutput("ClusterTabAnno"), style = "font-size:12px"),
                                       uiOutput("renddnldClusterTabAnno")
                              )
                            )
                          )
                        )
                      )
                      
                      
             ),
             
             ####----Matrix Connectivity----####
             
             tabPanel("Matrix Clustering",
                      fluidPage(
                        title = "Matrix Clustering",
                        sidebarLayout(
                          sidebarPanel(
                            width = 3,
                            tabsetPanel(
                              
                              ####----Matrix Cluster Parameters----####
                              
                              tabPanel("Cluster Parameters",
                                       p(),
                                       fileInput("UserMatrixFile","Upload Matrix", accept = c(".txt",".tsv",".csv")),
                                       uiOutput("rendHRLine3"),
                                       #shiny::hr(),
                                       uiOutput("rendClusterParamsHeader2"),
                                       #h4("Clustering Parameters"),
                                       fluidRow(
                                         column(6,
                                                uiOutput("rendClustMethodMat")
                                                #selectInput("ClustMethodMat","Method:",
                                                #            choices = c("ward.D", "complete", "ward.D2", "single", "average", "mcquitty", "median", "centroid"))
                                         ),
                                         column(6,
                                                uiOutput("rendNumClustersMat")
                                                #numericInput("NumClustersMat", step = 1, label = "Number of Clusters", value = 10)
                                         )
                                       ),
                                       uiOutput("rendViewClustTabMat"),
                                       #checkboxInput("ViewClustTabMat","View Cluster Results Table"),
                                       uiOutput("rendClustTabMat"),
                                       uiOutput("renddnldClustTabMat"),
                                       #downloadButton("dnldClustTabMat","Download Cluster Result"),
                                       uiOutput("rendHRLine4"),
                                       #shiny::hr(),
                                       uiOutput("rendSIFdnldHeader2"),
                                       #h4("SIF Download"),
                                       uiOutput("rendPrevSIFMat"),
                                       #checkboxInput("PrevSIFMat","Preview SIF File",value = F),
                                       uiOutput("rendSIFPreviewMat"),
                                       uiOutput("renddnldSIFTabMat"),
                                       #downloadButton("dnldSIFTabMat","Download SIF File")
                              ),
                              
                              ####----Matrix Figure Parameters----####
                              
                              tabPanel("Figure Parameters",
                                       p(),
                                       h4("Heatmap Parameters"),
                                       selectInput("ColorPaletteMat", "Select Color Palette:",
                                                   choices = c("Red/Blue" = "original",
                                                               "OmniBlueRed" = "OmniBlueRed",
                                                               "LightBlue/BlackRed" = "LightBlueBlackRed",
                                                               "Green/Black/Red" = "GreenBlackRed",
                                                               "Yellow/Green/Blue" = "YlGnBu","Inferno" = "Inferno",
                                                               "Viridis" = "Viridis","Plasma" = "Plasma",
                                                               "Reds" = "OrRd","Blues" = "PuBu","Greens" = "Greens")),
                                       fluidRow(
                                         column(6,
                                                checkboxInput("HeatColNamesMat","Show Heatmap Column Names", value = T),
                                                numericInput("HeatColFontMat", "Heatmap Column Font Size:",
                                                             min = 5, max = 75,
                                                             value = 12, step = 1),
                                                numericInput("HeatColDendHeightMat","Column Dendrogram Height:",
                                                             value = 50, step = 1)
                                         ),
                                         column(6,
                                                checkboxInput("HeatRowNamesMat","Show Heatmap Row Names", value = T),
                                                numericInput("HeatRowFontMat", "Heatmap Row Font Size:",
                                                             min = 5, max = 75,
                                                             value = 9, step = 1),
                                                numericInput("HeatRowDendHeightMat","Row Dendrogram Height:",
                                                             value = 50, step = 1)
                                         )
                                       ),
                                       shiny::hr(),
                                       h4("Connectivity Visualization Parameters"),
                                       selectInput("ConnViewMat","View Connectivity as:",
                                                   choices = c("Phylogeny" = "phylogenic","Dendrogram" = "rectangle","Circular" = "circular")),
                                       uiOutput("rendPhyloLayoutMat"),
                                       fluidRow(
                                         column(6,
                                                checkboxInput("ShowConnLabelsMat","Show Labels:",value = T)
                                         ),
                                         column(6,
                                                numericInput("ConnFontSizeMat","Font Size:",
                                                             value = 0.4, step = 0.1)
                                         )
                                       )
                              )
                            )
                          ),
                          
                          ####----Matrix Main Panel----####
                          
                          mainPanel(
                            tabsetPanel(
                              tabPanel("Table Preview",
                                       p(),
                                       div(DT::dataTableOutput("UserMatrixPreview"), style = "font-size:12px")
                              ),
                              tabPanel("Heatmap",
                                       p(),
                                       shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput('UserMatrixHeatmap', width = "100%", height = "800px")), type = 6),
                                       downloadButton("dnldUserMatrixHeatmapPDF","Download as PDF"),
                                       downloadButton("dnldUserMatrixHeatmapSVG","Download as SVG")
                              ),
                              tabPanel("Clustering",
                                       uiOutput("rendUserMatixCLustering"),
                                       downloadButton("dnldUserMatixCLusteringPDF","Download as PDF"),
                                       downloadButton("dnldUserMatixCLusteringSVG","Download as SVG")
                              )
                            )
                          )
                        )
                      )
             )
  )






server <- function(input, output, session) {
  
  ####----Jaccard Connectivity----####
  
  ####----Render UI----####
  
  output$rendUseExpData <- renderUI({
    
    if (is.null(PreSet_CoxH_Ranked_Features)) {
      actionButton("UseExpData","Load Example Data")
    }
    
  })
  
  output$rendDownloadLink <- renderUI({
    
    if (is.null(PreSet_CoxH_Ranked_Features)) {
      tags$a(href="http://shawlab.science/shiny/PATH_SURVEYOR_ExampleData/Pathway_Connectivity_App/PathwayFile_Upload/", "Download example data", target='_blank')
    }
    
  })
  
  output$rendHRLine1 <- renderUI({
    
    
    req(isTruthy(input$UserPathwayFile) | isTruthy(input$UseExpData) | isTruthy(PreSet_CoxH_Ranked_Features))
    shiny::hr()
    
  })
  
  output$rendHRLine2 <- renderUI({
    
    
    req(isTruthy(input$UserPathwayFile) | isTruthy(input$UseExpData) | isTruthy(PreSet_CoxH_Ranked_Features))
    shiny::hr()
    
  })
  
  output$renddnldClustTabIntra <- renderUI({
    
    
    req(isTruthy(input$UserPathwayFile) | isTruthy(input$UseExpData) | isTruthy(PreSet_CoxH_Ranked_Features))
    downloadButton("dnldClustTabIntra","Cluster Results")
    
  })
  
  output$renddnldSIFTabIntra <- renderUI({
    
    
    req(isTruthy(input$UserPathwayFile) | isTruthy(input$UseExpData) | isTruthy(PreSet_CoxH_Ranked_Features))
    downloadButton("dnldSIFTabIntra","SIF File")
    
  })
  
  output$rendClusterHeader <- renderUI({
    
    
    req(isTruthy(input$UserPathwayFile) | isTruthy(input$UseExpData) | isTruthy(PreSet_CoxH_Ranked_Features))
    h4("Clustering Parameters")
    
  })
  
  output$rendClustMethodIntra <- renderUI({
    
    
    req(isTruthy(input$UserPathwayFile) | isTruthy(input$UseExpData) | isTruthy(PreSet_CoxH_Ranked_Features))
    selectInput("ClustMethodIntra","Method:",
                choices = c("ward.D", "complete", "ward.D2", "single", "average", "mcquitty", "median", "centroid"))
    
  })
  
  output$rendNumClusters <- renderUI({
    
    
    req(isTruthy(input$UserPathwayFile) | isTruthy(input$UseExpData) | isTruthy(PreSet_CoxH_Ranked_Features))
    numericInput("NumClusters", step = 1, label = "Number of Clusters", value = 10)
    
  })
  
  output$rendViewClustTabIntra <- renderUI({
    
    
    req(isTruthy(input$UserPathwayFile) | isTruthy(input$UseExpData) | isTruthy(PreSet_CoxH_Ranked_Features))
    checkboxInput("ViewClustTabIntra","View Cluster Results")
    
  })
  
  output$rendSIFheader <- renderUI({
    
    
    req(isTruthy(input$UserPathwayFile) | isTruthy(input$UseExpData) | isTruthy(PreSet_CoxH_Ranked_Features))
    h4("SIF Download")
    
  })
  
  output$rendJaccDistCutoff <- renderUI({
    
    
    req(isTruthy(input$UserPathwayFile) | isTruthy(input$UseExpData) | isTruthy(PreSet_CoxH_Ranked_Features))
    numericInput("JaccDistCutoff","Jaccard Distance Cutoff",
                 min = 0,max = 1, step = 0.1, value = 0.9, width = "200px")
    
  })
  
  output$rendPrevSIF <- renderUI({
    
    
    req(isTruthy(input$UserPathwayFile) | isTruthy(input$UseExpData) | isTruthy(PreSet_CoxH_Ranked_Features))
    checkboxInput("PrevSIF","Preview SIF File",value = F)
    
  })
  
  output$rendTopFeatureSelect <- renderUI({
    
    
    req(isTruthy(input$UserPathwayFile) | isTruthy(input$UseExpData) | isTruthy(PreSet_CoxH_Ranked_Features))
    ranked_file <- user_file_loaded()
    FileName <- fileInName()
    ext <- tools::file_ext(FileName)
    if (length(ext) > 0) {
      if (ext == "txt" | ext == "tsv") {
        if (nrow(ranked_file) > 2) {
          numericInput("TopFeatureSelect","Number of Top Features for Clustering:",
                       value = 50, min = 3, width = "300px")
        }
      }
      else if (ext == "RData") {
        if (length(ranked_file) > 2) {
          numericInput("TopFeatureSelect","Number of Top Features for Clustering:",
                       value = 50, min = 3, width = "300px")
        }
      }
    }
    
  })
  
  output$rendClustTabIntra <- renderUI({
    
    
    req(isTruthy(input$UserPathwayFile) | isTruthy(input$UseExpData) | isTruthy(PreSet_CoxH_Ranked_Features))
    if (!is.null(input$ViewClustTabIntra)) {
      if (input$ViewClustTabIntra == TRUE) {
        
        div(DT::dataTableOutput("ClustTabIntra"), style = "font-size:10px; overflow-X: scroll")
        
      }
    }
    
  })
  
  output$rendSIFPreview <- renderUI({
    
    
    req(isTruthy(input$UserPathwayFile) | isTruthy(input$UseExpData) | isTruthy(PreSet_CoxH_Ranked_Features))
    if (!is.null(input$PrevSIF)) {
      if (input$PrevSIF == TRUE) {
        
        #DT::dataTableOutput("SIFPreview")
        div(DT::dataTableOutput("SIFPreview"), style = "font-size:10px; overflow-X: scroll")
        
      }
    }
    
  })
  
  output$renddnldJaccTabIntra <- renderUI({
    
    
    req(isTruthy(input$UserPathwayFile) | isTruthy(input$UseExpData) | isTruthy(PreSet_CoxH_Ranked_Features))
    downloadButton("dnldJaccTabIntra","Download Jaccard Connectivity Table")
    
  })
  
  output$renddnldClusterTabAnno <- renderUI({
    
    req(input$UserAnnotationFile)
    downloadButton("dnldClusterTabAnno","Download Cluster Annotation Table")
    
  })
  
  output$rendJaccDendo <- renderUI({
    
    VisType <- input$ConnView
    if (VisType == "rectangle" | VisType == "phylogenic") {
      
      shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotly::plotlyOutput('JaccDendoIntra', width = "100%", height = "900px")), type = 6)
      
    }
    else if (VisType == "circular") {
      
      shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput('JaccDendoIntraCirc', width = "100%", height = "900px")), type = 6)
      
    }
    
  })
  
  output$rendPhyloLayout <- renderUI({
    
    if (input$ConnView == "phylogenic") {
      
      selectInput("PhyloLayout","Choose Phylogenic Tree Layout:",
                  choices = c("Auto Layout" = "layout.auto",
                              "Force Directed (DrL) Layout" = "layout_with_drl",
                              "Tree Layout" = "layout_as_tree",
                              "GEM Layout" = "layout.gem",
                              "Mulitdimensional Scaling Layout" = "layout.mds",
                              "Large Graph Layout" = "layout_with_lgl"))
      
    }
    
  })
  
  ####----Reactives----####
  
  user_file_loaded <- reactiveVal()
  user_anno <- reactiveVal()
  
  fileInName <- reactiveVal()
  fileInName_base <- reactiveVal()
  
  observeEvent(input$UseExpData, {
    input <- as.data.frame(readr::read_delim(ExamplePathway_File, delim = '\t', col_names = T))
    input[,1] <- gsub("[[:punct:]]",".",input[,1])
    user_file_loaded(input)
    
    ext <- tools::file_ext(ExamplePathway_File)
    fileInName(ExamplePathway_File)
    Example_File_Name <- gsub(paste0(".",ext), "", basename(ExamplePathway_File))
    fileInName_base(Example_File_Name)
  })
  
  observe({
    
    if (!is.null(PreSet_CoxH_Ranked_Features)) {
      ranked_file <- PreSet_CoxH_Ranked_Features
      user_file_loaded(ranked_file)
      fileInName(PreSet_CoxH_Ranked_Feature_file)
      ext <- tools::file_ext(PreSet_CoxH_Ranked_Feature_file)
      Example_File_Name <- gsub(paste0(".",ext), "", basename(PreSet_CoxH_Ranked_Feature_file))
      fileInName_base(Example_File_Name)
    } else {
      header_check <- input$HeaderCheckIntra
      gs.u <- input$UserPathwayFile
      fileInName(gs.u$datapath)
      ext <- tools::file_ext(gs.u$datapath)
      req(gs.u)
      #validate(need(ext == c("tsv","txt","gmt"), "Please upload .tsv, .txt, or .gmt file"))
      
      if (ext == "gmt") {
        gmt <- clusterProfiler::read.gmt(gs.u$datapath)
        ranked_file <- data.frame(GeneSet = unique(gmt[,1]))
        ranked_file[,1] <- gsub("[[:punct:]]",".",ranked_file[,1])
      }
      if (ext == "RData") {
        loadRData <- function(fileName){
          #loads an RData file, and returns it
          load(fileName)
          get(ls()[ls() != "fileName"])
        }
        ranked_file <- loadRData(gs.u$datapath)
      }
      if (ext == "tsv" | ext == "txt") {
        ranked_file <- as.data.frame(readr::read_delim(gs.u$datapath, delim = '\t', col_names = header_check, comment = "#"))
        ranked_file[,1] <- gsub("[[:punct:]]",".",ranked_file[,1])
      }
      
      user_file_loaded(ranked_file)
    }
    
  })
  
  ## User Upload gene set list reactive
  user_gs <- reactive({
    
    
    req(isTruthy(input$UserPathwayFile) | isTruthy(input$UseExpData) | isTruthy(PreSet_CoxH_Ranked_Features))
    
    if (!is.null(PreSet_CoxH_Ranked_Features)) {
      ranked_file <- PreSet_CoxH_Ranked_Features
      if (ncol(ranked_file) == 2) {
        ranked_file <- ranked_file[which(ranked_file[,2] != "NA"),]
        gs_u2 <- list()
        for (i in unique(ranked_file[,1])){
          gs_u2[[i]] <- ranked_file[ranked_file[,1] == i,][,2]
        }
        gs_u2
      }
      else {
        TopFeat <- input$TopFeatureSelect
        if (isTruthy(TopFeat)) {
          if (TopFeat > 0) {
            paths <- as.vector(ranked_file[,1])
            paths <- gsub("[[:punct:]]",".",paths)
            gs_u <- gs[unique(paths)]
            if (TopFeat > length(gs_u)) {
              TopFeat <- length(gs_u)
            }
            gs_u2 <- gs_u[c(1:TopFeat)]
            gs_u2 <- Filter(Negate(is.null), gs_u2)
            gs_u2
          }
        }
      }
    } else {
      
      ranked_file <- user_file_loaded()
      FileName <- fileInName()
      ext <- tools::file_ext(FileName)
      if (ext == "gmt") {
        gmt <- clusterProfiler::read.gmt(gs.u$datapath)
        gmt <- gmt[which(gmt[,2] != "NA"),]
        gs_u2 <- list()
        for (i in unique(gmt[,1])){
          gs_u2[[i]] <- gmt[gmt[,1] == i,]$gene
        }
        gs_u2
      } else if (ext == "RData") {
        gs_u2 <- ranked_file
        if (is.null(input$TopFeatureSelect)) {
          TopFeat <- 50
          gs_u2 <- gs_u2[1:TopFeat]
          gs_u2
        } else {
          TopFeat <- input$TopFeatureSelect
          if (TopFeat > length(gs_u2)) {
            TopFeat <- length(gs_u2)
          }
          gs_u2 <- gs_u2[1:TopFeat]
          gs_u2
        }
      } else if (ext == "txt" | ext == "tsv") {
        if (ncol(ranked_file) == 2) {
          ranked_file <- ranked_file[which(ranked_file[,2] != "NA"),]
          genesets <- unique(ranked_file[,1])
          if (is.null(input$TopFeatureSelect)) {
            TopFeat <- 50
            showGS <- genesets[1:TopFeat]
            ranked_file <- ranked_file[which(ranked_file[,1] %in% showGS),]
          } else {
            TopFeat <- input$TopFeatureSelect
            if (TopFeat > nrow(ranked_file)) {
              TopFeat <- nrow(ranked_file)
            }
            showGS <- genesets[1:TopFeat]
            ranked_file <- ranked_file[which(ranked_file[,1] %in% showGS),]
          }
          gs_u2 <- list()
          for (i in unique(ranked_file[,1])){
            gs_u2[[i]] <- ranked_file[ranked_file[,1] == i,][,2]
          }
          gs_u2
        } else {
          if (is.null(input$TopFeatureSelect)) {
            TopFeat <- 50
            paths <- as.vector(ranked_file[,1])
            paths <- gsub("[[:punct:]]",".",paths)
            gs_u <- gs[unique(paths)]
            gs_u2 <- gs_u[c(1:TopFeat)]
            gs_u2 <- Filter(Negate(is.null), gs_u2)
            gs_u2
          } else {
            TopFeat <- input$TopFeatureSelect
            paths <- as.vector(ranked_file[,1])
            paths <- gsub("[[:punct:]]",".",paths)
            gs_u <- gs[unique(paths)]
            if (TopFeat > length(gs_u)) {
              TopFeat <- length(gs_u)
            }
            gs_u2 <- gs_u[c(1:TopFeat)]
            gs_u2 <- Filter(Negate(is.null), gs_u2)
            gs_u2
          }
        }
      }
      
      #gs_u2
    }
    
  })
  
  observeEvent(input$useExpAnnoData, {
    data <- as.data.frame(readr::read_delim(ExampleAnnotation_File, delim = '\t', col_names = T, comment = "#"))
    user_anno(data)
  })
  
  observe({
    
    gs.u <- input$UserAnnotationFile
    ext <- tools::file_ext(gs.u$datapath)
    req(gs.u)
    validate(need(ext == c("tsv","txt","csv"), "Please upload .tsv, .txt, or .csv file"))
    if (ext == "txt" | ext == "tsv") {
      
      df <- as.data.frame(readr::read_delim(gs.u$datapath, delim = '\t', col_names = T, comment = "#"))
      
    }
    else if (ext == "csv") {
      
      df <- as.data.frame(readr::read_delim(gs.u$datapath, delim = ',', col_names = T))
      
    }
    
    user_anno(df)
    
  })
  
  ## Jaccard Matrix Generation
  jacc_react_Intra <- reactive({
    
    if (length(user_gs()) > 0) {
      
      gs1 <- user_gs() #User Gene Sets Input
      gs2 <- user_gs() #User Gene Sets Input
      
      jac_list <- list()
      for (i in 1:length(gs1)) {
        
        jac_vect <- c()
        
        for (j in 1:length(gs2)) {
          
          gs1_genes <- gs1[[i]]
          gs2_genes <- gs2[[j]]
          jac_index <- jaccard(gs1_genes,gs2_genes)
          jac_vect <- c(jac_vect,jac_index)
          
        }
        
        jac_list[[paste(names(gs1)[i],sep = "")]] <- jac_vect
        
      }
      
      jac_df_Intra <- do.call(rbind,jac_list)
      colnames(jac_df_Intra) <- names(gs2)
      
      jac_df_Intra
      
    }
    
  })
  
  ClusterTabAnno_react <- reactive({
    
    if (!is.null(user_anno())) {
      
      ## user annotation table
      anno_df <- user_anno()
      
    }
    
    
    ## user gs list
    gs_u <- user_gs()
    
    ## cluster table
    ## Data matrix
    jacc_df <- jacc_react_Intra()
    jacc_mat <- as.matrix(jacc_df)
    clust_me <- input$ClustMethodIntra
    hm_res <- pheatmap::pheatmap(jacc_mat,
                                 clustering_method = clust_me,
                                 silent = T)
    jacc_df <- jacc_react_Intra()
    cut_k <- input$NumClusters
    clustTab <- data.frame(Pathways = rownames(jacc_df),
                           Cluster = cutree(hm_res$tree_row,k = cut_k))
    
    ## table of gene set with genes
    nameVector <- unlist(mapply(function(x,y){ rep(y, length(x)) }, gs_u, names(gs_u)))
    resultDF <- cbind.data.frame(Pathways = nameVector,Genes = unlist(gs_u))
    rownames(resultDF) <- 1:nrow(resultDF)
    
    ## double check punct matches
    clustTab$Pathways <- gsub("[[:punct:]]",".",clustTab$Pathways)
    resultDF$Pathways <- gsub("[[:punct:]]",".",resultDF$Pathways)
    
    ## Merge geneset and cluster table
    df2 <- merge(resultDF,clustTab, by = "Pathways")
    
    if (!is.null(user_anno())) {
      
      colnames(anno_df)[1] <- "Genes"
      df3 <- merge(df2,anno_df,by = "Genes", all.x = T)
      
    }
    else if (is.null(user_anno())) {
      
      df3 <- df2
      
    }
    
    
    df4 <- df3 %>%
      relocate(Pathways,Cluster,Genes)
    
    df5 <- df4[order(df4$Pathways),]
    df5
    
  })
  
  SIF_tab_PAthConn <- reactive({
    
    jacc_df <- jacc_react_Intra()
    dist_cutoff <- input$JaccDistCutoff
    jacc_df <- as.data.frame(jacc_df)
    jacc_df_melt <- melt(jacc_df)
    jacc_cols <- colnames(jacc_df)
    jacc_df_melt$Pathway_B <- rep_len(jacc_cols,length.out = nrow(jacc_df_melt))
    colnames(jacc_df_melt)[c(1,2)] <- c("Pathway_A","Jaccard_Distance")
    jacc_df_melt <- jacc_df_melt[which(jacc_df_melt$Jaccard_Distance <= dist_cutoff),]
    jacc_df_melt <- jacc_df_melt[which(jacc_df_melt$Pathway_A != jacc_df_melt$Pathway_B),]
    jacc_df_melt <- jacc_df_melt[!duplicated(t(apply(jacc_df_melt,1,sort))),]
    jacc_df_melt
    
  })
  
  ####----Data Tables----####
  
  ## Jaccard Matrix Table - INTRA pathway
  output$JaccTableIntra <- DT::renderDataTable({
    
    req(isTruthy(input$UserPathwayFile) | isTruthy(input$UseExpData) | isTruthy(PreSet_CoxH_Ranked_Features))
    ## Jaccard Table
    jacc_df <- as.data.frame(jacc_react_Intra())
    if (sum(dim(jacc_df)) > 0) {
      colnames(jacc_df) <- gsub("[[:punct:]]", " ",colnames(jacc_df))
      
      DT::datatable(jacc_df,
                    extensions = "FixedColumns",
                    options = list(lengthMenu = c(10,20,50,100,1000,5000,10000),
                                   pageLength = 20,
                                   scrollX = TRUE,
                                   autoWidth = TRUE,
                                   fixedColumns = list(leftColumns = 1)),
                    selection=list(mode = "multiple")) %>%
        formatRound(columns = c(1:ncol(jacc_df)), digits = 4)
    }
    
    
  })
  
  ## Cluster Table
  output$ClustTabIntra <- DT::renderDataTable({
    
    ## Data matrix
    jacc_df <- jacc_react_Intra()
    jacc_mat <- as.matrix(jacc_df)
    clust_me <- input$ClustMethodIntra
    
    hm_res <- pheatmap::pheatmap(jacc_mat,
                                 clustering_method = clust_me,
                                 silent = T)
    
    jacc_df <- jacc_react_Intra()
    cut_k <- input$NumClusters
    clustTab <- data.frame(Pathways = rownames(jacc_df),
                           Cluster = cutree(hm_res$tree_row,k = cut_k))
    
    DT::datatable(clustTab,
                  options = list(keys = TRUE,
                                 searchHighlight = TRUE,
                                 scrollX = TRUE,
                                 pageLength = 10,
                                 lengthMenu = c("10", "25", "50", "100")),
                  rownames = F)
    
  })
  
  output$ClusterTabAnno <- DT::renderDataTable({
    
    df5 <- ClusterTabAnno_react()
    
    DT::datatable(df5,
                  options = list(lengthMenu = c(10,20,50,100,1000,5000,10000),
                                 pageLength = 20,
                                 scrollX = TRUE),
                  rownames = F,
                  selection=list(mode = "multiple"))
    
    
  })
  
  output$SIFPreview <- DT::renderDataTable({
    
    jacc_df_melt <- SIF_tab_PAthConn()
    colnames(jacc_df_melt) <- gsub("_"," ",colnames(jacc_df_melt))
    DT::datatable(jacc_df_melt,
                  options = list(keys = TRUE,
                                 searchHighlight = TRUE,
                                 pageLength = 10,
                                 lengthMenu = c("10", "25", "50", "100")),
                  rownames = F) %>%
      formatRound(columns = 2, digits = 4)
    
  })
  
  ####----Plots----####
  
  ## Jaccard Heatmap - INTRA
  JaccHeatmapIntra_react <- reactive({
    
    ## User inputs
    row_font <- input$HeatRowFontIntra
    col_font <- input$HeatColFontIntra
    row_name <- input$HeatRowNamesIntra
    col_name <- input$HeatColNamesIntra
    col_pall <- input$ColorPaletteIntra
    clust_me <- input$ClustMethodIntra
    row_dend <- input$HeatRowDendHeight
    col_dend <- input$HeatColDendHeight
    
    ## Data matrix
    jacc_df <- jacc_react_Intra()
    jacc_mat <- as.matrix(jacc_df)
    
    ## Color choice
    minimum = min(jacc_mat)
    maximum = max(jacc_mat)
    bk = c(seq(minimum,minimum/2, length=100), seq(minimum/2,maximum/2,length=100),seq(maximum/2,maximum,length=100))
    #Heatmap color
    col_sets <- c("OrRd","PuBu","Greens","YlGnBu")
    if (col_pall == "original") {
      HeatMap_Colors <- c("dark blue","blue","white","red", "dark red")
      hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
    }
    else if (col_pall %in% col_sets) {
      HeatMap_Colors <- brewer.pal(n = 5, col_pall)
      hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
    }
    else if (col_pall == "Inferno") {
      hmcols <- inferno(500)
    }
    else if (col_pall == "Viridis") {
      hmcols <- viridis::viridis(500)
    }
    else if (col_pall == "Plasma") {
      hmcols <- plasma(500)
    }
    else if (col_pall == "OmniBlueRed") {
      hmcols<- colorRampPalette(c("#1984C5", "#22A7F0", "#63BFF0", "#A7D5ED", "#E2E2E2", "#E1A692", "#DE6E56", "#E14B31", "#C23728"))(length(bk)-1)
    }
    else if (col_pall == "LightBlueBlackRed") {
      hmcols<- colorRampPalette(c("#34C5FD","black","red"))(length(bk)-1)
    }
    else if (col_pall == "GreenBlackRed") {
      hmcols<- colorRampPalette(c("green","black","red"))(length(bk)-1)
    }
    heat <- pheatmap::pheatmap(jacc_mat,
                               fontsize_row = row_font,
                               fontsize_col = col_font,
                               show_rownames = row_name,
                               show_colnames = col_name,
                               treeheight_col = col_dend,
                               treeheight_row = row_dend,
                               angle_col = 90,
                               border_color = NA,
                               clustering_method = clust_me,
                               color = hmcols)
    heat
    
  })
  
  output$JaccHeatmapIntra <- renderPlot({
    
    heat <- JaccHeatmapIntra_react()
    heat
    
  })
  
  JaccDendoIntra_react <- reactive({
    
    VisType <- input$ConnView
    
    if (VisType == "rectangle" | VisType == "phylogenic") {
      
      jacc_mat <- jacc_react_Intra()
      clust_me <- input$ClustMethodIntra
      nclust <- input$NumClusters
      LableChoice <- input$ShowConnLabels
      FontSize <- input$ConnFontSize
      
      
      hc <- hclust(dist(jacc_mat), method = clust_me)
      hc_dend <- as.dendrogram(hc)
      dend_data <- dendro_data(hc_dend)
      lab_order <- label(dend_data)[3]
      lab_vect <- lab_order$label
      clust_df <- data.frame(Pathways = rownames(jacc_mat),cluster = cutree(hc,k = nclust))
      clust_df2 <- clust_df[match(lab_vect,clust_df$Pathways),]
      colnames(clust_df2)[1] <- "label"
      clust_df2$cluster <- as.factor(clust_df2$cluster)
      
      dend_data[["labels"]] <- merge(dend_data[["labels"]],clust_df2, by="label")
      
      if (VisType == "rectangle") {
        
        FontSize <- FontSize * 6
        g <- ggdendro::ggdendrogram(dend_data,
                                    rotate = TRUE,
                                    theme_dendro = TRUE,
                                    leaf_labels = F,
                                    segments = T,
                                    labels = F) + 
          labs(title = "") +
          scale_y_continuous(expand = c(0.4, 0))
        if (LableChoice == TRUE) {
          g <- g + geom_text(data=label(dend_data), aes(x, y, label=label, hjust=1, color=cluster),size=FontSize)
        }
        
        gp <- plotly::ggplotly(g)
        
      }
      else if (VisType == "phylogenic") {
        
        if (is.null(input$PhyloLayout) == TRUE) {
          PhyloLayout <- "layout.auto"
        }
        else if (is.null(input$PhyloLayout) == FALSE) {
          PhyloLayout <- input$PhyloLayout
        }
        
        
        if (LableChoice == TRUE) {
          k = nclust
          gp <- plotly::ggplotly(
            factoextra::fviz_dend(
              hc_dend,
              cex = FontSize,
              k = nclust,
              show_labels = TRUE,
              type = "phylogenic",
              phylo_layout = PhyloLayout)
          )
          for(i in seq(k, 1)){
            gp$x$data[[i+1]]$text <- gp$x$data[[i+1+k]]$text
          }
          
        }
        else if (LableChoice == FALSE) {
          k = nclust
          gp <- plotly::ggplotly(
            factoextra::fviz_dend(
              hc,
              cex = FontSize,
              k = k,
              show_labels = TRUE,
              type = "phylogenic",
              phylo_layout = PhyloLayout)
          )
          gp$x$data[[1]]$text <- NA 
          for(i in seq(k, 1)){
            gp$x$data[[i+1]]$text <- gp$x$data[[i+1+k]]$text
            gp$x$data[[i + 1 + k]] <- NULL
          }
          
        }
        
      }
      gp
    }
    
  })
  
  output$JaccDendoIntra <- renderPlotly({
    
    gp <- JaccDendoIntra_react()
    gp
    
  })
  
  JaccDendoIntraCirc_react <- reactive({
    
    VisType <- input$ConnView
    
    if (VisType == "circular") {
      
      jacc_mat <- jacc_react_Intra()
      clust_me <- input$ClustMethodIntra
      nclust <- input$NumClusters
      LableChoice <- input$ShowConnLabels
      FontSize <- input$ConnFontSize
      
      hc <- hclust(dist(jacc_mat), method = clust_me)
      hc_dend <- as.dendrogram(hc)
      f <- factoextra::fviz_dend(hc_dend, k = nclust, cex = FontSize, color_labels_by_k = T, type = "circular", show_labels = LableChoice) +
        theme(legend.position="none") +
        labs(title = "") +
        theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.title.x = element_blank(),
              axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank())
      
      
      f
      
    }
    
  })
  
  output$JaccDendoIntraCirc <- renderPlot({
    
    f <- JaccDendoIntraCirc_react()
    f
    
  })
  
  
  ####----Matrix Clustering----####
  
  ####----Render UI----####
  
  output$rendHRLine3 <- renderUI({
    
    req(input$UserMatrixFile)
    shiny::hr()
    
  })
  
  output$rendClusterParamsHeader2 <- renderUI({
    
    req(input$UserMatrixFile)
    h4("Clustering Parameters")
    
  })
  
  output$rendClustMethodMat <- renderUI({
    
    req(input$UserMatrixFile)
    selectInput("ClustMethodMat","Method:",
                choices = c("ward.D", "complete", "ward.D2", "single", "average", "mcquitty", "median", "centroid"))
    
  })
  
  output$rendNumClustersMat <- renderUI({
    
    req(input$UserMatrixFile)
    numericInput("NumClustersMat", step = 1, label = "Number of Clusters", value = 10)
    
  })
  
  output$rendViewClustTabMat <- renderUI({
    
    req(input$UserMatrixFile)
    checkboxInput("ViewClustTabMat","View Cluster Results Table")
    
  })
  
  output$renddnldClustTabMat <- renderUI({
    
    req(input$UserMatrixFile)
    downloadButton("dnldClustTabMat","Download Cluster Result")
    
  })
  
  output$rendHRLine4 <- renderUI({
    
    req(input$UserMatrixFile)
    shiny::hr()
    
  })
  
  output$rendSIFdnldHeader2 <- renderUI({
    
    req(input$UserMatrixFile)
    h4("SIF Download")
    
  })
  
  output$rendPrevSIFMat <- renderUI({
    
    req(input$UserMatrixFile)
    checkboxInput("PrevSIFMat","Preview SIF File",value = F)
    
  })
  
  output$renddnldSIFTabMat <- renderUI({
    
    req(input$UserMatrixFile)
    downloadButton("dnldSIFTabMat","Download SIF File")
    
  })
  
  output$rendClustTabMat <- renderUI({
    
    req(input$UserMatrixFile)
    if (input$ViewClustTabMat == TRUE) {
      
      div(DT::dataTableOutput("ClustTabMat"), style = "font-size:10px; height:400px; overflow-X: scroll")
      
    }
    
  })
  
  output$rendSIFPreviewMat <- renderUI({
    
    req(input$UserMatrixFile)
    if (input$PrevSIFMat == TRUE) {
      
      div(DT::dataTableOutput("SIFPreviewMat"), style = "font-size:10px; height:400px; overflow-X: scroll")
      
    }
    
  })
  
  output$rendUserMatixCLustering <- renderUI({
    
    VisType <- input$ConnViewMat
    if (VisType == "rectangle" | VisType == "phylogenic") {
      
      shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotly::plotlyOutput('JaccDendoMat', width = "100%", height = "900px")), type = 6)
      
    }
    else if (VisType == "circular") {
      
      shinycssloaders::withSpinner(shinyjqui::jqui_resizable(plotOutput('JaccDendoMatCirc', width = "100%", height = "900px")), type = 6)
      
    }
    
  })
  
  output$rendPhyloLayoutMat <- renderUI({
    
    if (input$ConnViewMat == "phylogenic") {
      
      selectInput("PhyloLayoutMat","Choose Phylogenic Tree Layout:",
                  choices = c("Auto Layout" = "layout.auto",
                              "Force Directed (DrL) Layout" = "layout_with_drl",
                              "Tree Layout" = "layout_as_tree",
                              "GEM Layout" = "layout.gem",
                              "Mulitdimensional Scaling Layout" = "layout.mds",
                              "Large Graph Layout" = "layout_with_lgl"))
      
    }
    
  })
  
  ####----Reactives----####
  
  user_matrix_loaded <- reactive({
    
    gs.u <- input$UserMatrixFile
    ext <- tools::file_ext(gs.u$datapath)
    req(gs.u)
    validate(need(ext == c("tsv","txt","csv"), "Please upload .tsv, .txt, or .csv file"))
    
    if (ext == "csv") {
      matrix.u <- as.data.frame(readr::read_delim(gs.u$datapath, delim = ',', col_names = T))
    }
    else {
      matrix.u <- as.data.frame(readr::read_delim(gs.u$datapath, delim = '\t', col_names = T))
    }
    
    
    matrix.u
    
  })
  
  
  ####----Data Tables----####
  
  output$UserMatrixPreview <- DT::renderDataTable({
    
    df <- user_matrix_loaded()
    
    DT::datatable(df,
                  extensions = "FixedColumns",
                  options = list(lengthMenu = c(10,20,50,100,1000,5000,10000),
                                 pageLength = 20,
                                 scrollX = TRUE,
                                 autoWidth = TRUE,
                                 fixedColumns = list(leftColumns = 1)),
                  rownames = F,
                  selection=list(mode = "multiple")) %>%
      formatRound(columns = c(2:ncol(df)), digits = 4)
    
    
  })
  
  ## Cluster Table
  output$ClustTabMat <- DT::renderDataTable({
    
    ## Data matrix
    df <- user_matrix_loaded()
    rownames(df) <- df[,1]
    df <- df[,-1]
    df_mat <- as.matrix(df)
    clust_me <- input$ClustMethodMat
    
    hm_res <- pheatmap::pheatmap(df_mat,
                                 clustering_method = clust_me,
                                 silent = T)
    
    cut_k <- input$NumClustersMat
    clustTab <- data.frame(Pathways = rownames(df),
                           Cluster = cutree(hm_res$tree_row,k = cut_k))
    
    DT::datatable(clustTab,
                  options = list(keys = TRUE,
                                 searchHighlight = TRUE,
                                 pageLength = 10,
                                 lengthMenu = c("10", "25", "50", "100")),
                  rownames = F)
    
  })
  
  output$SIFPreviewMat <- DT::renderDataTable({
    
    df <- user_matrix_loaded()
    rownames(df) <- df[,1]
    df <- df[,-1]
    df <- as.data.frame(df)
    df_melt <- melt(df)
    df_cols <- colnames(df)
    df_melt$Sample_B <- rep_len(df_cols,length.out = nrow(df_melt))
    colnames(df_melt)[c(1,2)] <- c("Sample_A","Distance")
    df_melt <- df_melt[which(df_melt$Sample_A != df_melt$Sample_B),]
    df_melt <- df_melt[!duplicated(t(apply(df_melt,1,sort))),]
    DT::datatable(df_melt,
                  options = list(keys = TRUE,
                                 searchHighlight = TRUE,
                                 pageLength = 10,
                                 lengthMenu = c("10", "25", "50", "100")),
                  rownames = F)
    
  })
  
  ####----Plots----####
  
  UserMatrixHeatmap_react <- reactive({
    
    ## User inputs
    row_font <- input$HeatRowFontMat
    col_font <- input$HeatColFontMat
    row_name <- input$HeatRowNamesMat
    col_name <- input$HeatColNamesMat
    col_pall <- input$ColorPaletteMat
    clust_me <- input$ClustMethodMat
    row_dend <- input$HeatRowDendHeightMat
    col_dend <- input$HeatColDendHeightMat
    
    ## Data matrix
    jacc_df <- user_matrix_loaded()
    rownames(jacc_df) <- jacc_df[,1]
    jacc_df <- jacc_df[,-1]
    jacc_mat <- as.matrix(jacc_df)
    
    ## Color choice
    minimum = min(jacc_mat)
    maximum = max(jacc_mat)
    bk = c(seq(minimum,minimum/2, length=100), seq(minimum/2,maximum/2,length=100),seq(maximum/2,maximum,length=100))
    #Heatmap color
    col_sets <- c("OrRd","PuBu","Greens","YlGnBu")
    if (col_pall == "original") {
      HeatMap_Colors <- c("dark blue","blue","white","red", "dark red")
      hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
    }
    else if (col_pall %in% col_sets) {
      HeatMap_Colors <- brewer.pal(n = 5, col_pall)
      hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
    }
    else if (col_pall == "Inferno") {
      hmcols <- inferno(500)
    }
    else if (col_pall == "Viridis") {
      hmcols <- viridis::viridis(500)
    }
    else if (col_pall == "Plasma") {
      hmcols <- plasma(500)
    }
    else if (col_pall == "OmniBlueRed") {
      hmcols<- colorRampPalette(c("#1984C5", "#22A7F0", "#63BFF0", "#A7D5ED", "#E2E2E2", "#E1A692", "#DE6E56", "#E14B31", "#C23728"))(length(bk)-1)
    }
    else if (col_pall == "LightBlueBlackRed") {
      hmcols<- colorRampPalette(c("#34C5FD","black","red"))(length(bk)-1)
    }
    else if (col_pall == "GreenBlackRed") {
      hmcols<- colorRampPalette(c("green","black","red"))(length(bk)-1)
    }
    heat <- pheatmap::pheatmap(jacc_mat,
                               fontsize_row = row_font,
                               fontsize_col = col_font,
                               show_rownames = row_name,
                               show_colnames = col_name,
                               treeheight_col = col_dend,
                               treeheight_row = row_dend,
                               angle_col = 90,
                               border_color = NA,
                               clustering_method = clust_me,
                               color = hmcols)
    heat
    
  })
  
  ## Jaccard Heatmap - INTRA
  output$UserMatrixHeatmap <- renderPlot({
    
    heat <- UserMatrixHeatmap_react()
    heat
    
  })
  
  JaccDendoMat_react <- reactive({
    
    VisType <- input$ConnViewMat
    
    if (VisType == "rectangle" | VisType == "phylogenic") {
      
      jacc_mat <- user_matrix_loaded()
      rownames(jacc_mat) <- jacc_mat[,1]
      jacc_mat <- jacc_mat[,-1]
      clust_me <- input$ClustMethodMat
      nclust <- input$NumClustersMat
      LableChoice <- input$ShowConnLabelsMat
      FontSize <- input$ConnFontSizeMat
      
      
      hc <- hclust(dist(jacc_mat), method = clust_me)
      hc_dend <- as.dendrogram(hc)
      dend_data <- dendro_data(hc_dend)
      lab_order <- label(dend_data)[3]
      lab_vect <- lab_order$label
      clust_df <- data.frame(Pathways = rownames(jacc_mat),cluster = cutree(hc,k = nclust))
      clust_df2 <- clust_df[match(lab_vect,clust_df$Pathways),]
      colnames(clust_df2)[1] <- "label"
      clust_df2$cluster <- as.factor(clust_df2$cluster)
      
      dend_data[["labels"]] <- merge(dend_data[["labels"]],clust_df2, by="label")
      
      if (VisType == "rectangle") {
        
        FontSize <- FontSize * 6
        g <- ggdendro::ggdendrogram(dend_data,
                                    rotate = TRUE,
                                    theme_dendro = TRUE,
                                    leaf_labels = F,
                                    segments = T,
                                    labels = F) + 
          labs(title = "") +
          scale_y_continuous(expand = c(0.4, 0))
        if (LableChoice == TRUE) {
          g <- g + geom_text(data=label(dend_data), aes(x, y, label=label, hjust=1, color=cluster),size=FontSize)
        }
        
        gp <- plotly::ggplotly(g)
        
      }
      else if (VisType == "phylogenic") {
        
        if (is.null(input$PhyloLayout) == TRUE) {
          PhyloLayout <- "layout.auto"
        }
        else if (is.null(input$PhyloLayout) == FALSE) {
          PhyloLayout <- input$PhyloLayout
        }
        
        
        if (LableChoice == TRUE) {
          k = nclust
          gp <- plotly::ggplotly(
            factoextra::fviz_dend(
              hc_dend,
              cex = FontSize,
              k = nclust,
              show_labels = TRUE,
              type = "phylogenic",
              phylo_layout = PhyloLayout)
          )
          for(i in seq(k, 1)){
            gp$x$data[[i+1]]$text <- gp$x$data[[i+1+k]]$text
          }
          
        }
        else if (LableChoice == FALSE) {
          k = nclust
          gp <- plotly::ggplotly(
            factoextra::fviz_dend(
              hc,
              cex = FontSize,
              k = k,
              show_labels = TRUE,
              type = "phylogenic",
              phylo_layout = PhyloLayout)
          )
          gp$x$data[[1]]$text <- NA 
          for(i in seq(k, 1)){
            gp$x$data[[i+1]]$text <- gp$x$data[[i+1+k]]$text
            gp$x$data[[i + 1 + k]] <- NULL
          }
          
        }
        
      }
      gp
    }
    
  })
  
  output$JaccDendoMat <- renderPlotly({
    
    gp <- JaccDendoMat_react()
    gp
    
  })
  
  JaccDendoMatCirc_react <- reactive({
    
    VisType <- input$ConnViewMat
    
    if (VisType == "circular") {
      
      jacc_mat <- user_matrix_loaded()
      clust_me <- input$ClustMethodMat
      nclust <- input$NumClustersMat
      LableChoice <- input$ShowConnLabelsMat
      FontSize <- input$ConnFontSizeMat
      
      hc <- hclust(dist(jacc_mat), method = clust_me)
      hc_dend <- as.dendrogram(hc)
      f <- factoextra::fviz_dend(hc_dend, k = nclust, cex = FontSize, color_labels_by_k = T, type = "circular", show_labels = LableChoice) +
        theme(legend.position="none") +
        labs(title = "") +
        theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.title.x = element_blank(),
              axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank())
      
      
      f
      
    }
    
  })
  
  output$JaccDendoMatCirc <- renderPlot({
    
    f <- JaccDendoMatCirc_react()
    f
    
  })
  
  
  
  
  ####----Downloads----####
  
  ####----Jaccard----####
  
  output$dnldJaccHeatmapIntraPDF <- downloadHandler(
    filename = function() {
      gs.u <- input$UserPathwayFile
      TopFeat <- input$TopFeatureSelect
      file_name <- gs.u$name
      file_name <- tools::file_path_sans_ext(gs.u)
      paste(file_name,"_Top",TopFeat,"_JaccardHeatmap.pdf", sep = '')
    },
    content = function(file) {
      heat <- JaccHeatmapIntra_react()
      ggsave(file,heat, width = 10, height = 12)
    }
  )
  output$dnldJaccHeatmapIntraSVG <- downloadHandler(
    filename = function() {
      gs.u <- input$UserPathwayFile
      TopFeat <- input$TopFeatureSelect
      file_name <- gs.u$name
      file_name <- tools::file_path_sans_ext(gs.u)
      paste(file_name,"_Top",TopFeat,"_JaccardHeatmap.svg", sep = '')
    },
    content = function(file) {
      heat <- JaccHeatmapIntra_react()
      ggsave(file,heat, width = 10, height = 12)
    }
  )
  output$dnldJaccDendoPDF <- downloadHandler(
    filename = function() {
      gs.u <- input$UserPathwayFile
      TopFeat <- input$TopFeatureSelect
      file_name <- gs.u$name
      file_name <- tools::file_path_sans_ext(gs.u)
      paste(file_name,"_Top",TopFeat,"_JaccardCluster.pdf", sep = '')
    },
    content = function(file) {
      dendo <- JaccDendoIntra_react()
      ggsave(file, dendo, width = 10, height = 8)
    }
  )
  output$dnldJaccDendoSVG <- downloadHandler(
    filename = function() {
      gs.u <- input$UserPathwayFile
      TopFeat <- input$TopFeatureSelect
      file_name <- gs.u$name
      file_name <- tools::file_path_sans_ext(gs.u)
      paste(file_name,"_Top",TopFeat,"_JaccardCluster.svg", sep = '')
    },
    content = function(file) {
      dendo <- JaccDendoIntra_react()
      ggsave(file, dendo, width = 10, height = 8)
    }
  )
  ## Download Cluster Table - INTRA
  output$dnldClustTabIntra <- downloadHandler(
    filename = function() {
      
      gs.u <- input$UserPathwayFile
      TopFeat <- input$TopFeatureSelect
      file_name <- gs.u$name
      file_name <- tools::file_path_sans_ext(gs.u)
      paste(file_name,"_Top",TopFeat,"_Cluster_Result.txt", sep = '')
      
    },
    content = function(file) {
      
      ## Data matrix
      jacc_df <- jacc_react_Intra()
      jacc_mat <- as.matrix(jacc_df)
      clust_me <- input$ClustMethodIntra
      
      hm_res <- pheatmap::pheatmap(jacc_mat,
                                   clustering_method = clust_me,
                                   silent = T)
      
      jacc_df <- jacc_react_Intra()
      cut_k <- input$NumClusters
      clustTab <- data.frame(Pathways = rownames(jacc_df),
                             Cluster = cutree(hm_res$tree_row,k = cut_k))
      readr::write_delim(clustTab,file,delim = '\t')
      
    }
  )
  output$dnldSIFTabIntra <- downloadHandler(
    filename = function() {
      
      gs.u <- input$UserPathwayFile
      TopFeat <- input$TopFeatureSelect
      file_name <- gs.u$name
      file_name <- tools::file_path_sans_ext(gs.u)
      paste(file_name,"_Top",TopFeat,"_sif_dend.sif", sep = '')
      
    },
    content = function(file) {
      
      jacc_df <- jacc_react_Intra()
      dist_cutoff <- input$JaccDistCutoff
      jacc_df <- as.data.frame(jacc_df)
      jacc_df_melt <- melt(jacc_df)
      jacc_cols <- colnames(jacc_df)
      jacc_df_melt$Pathway_B <- rep_len(jacc_cols,length.out = nrow(jacc_df_melt))
      colnames(jacc_df_melt)[c(1,2)] <- c("Pathway_A","Jaccard_Distance")
      jacc_df_melt <- jacc_df_melt[which(jacc_df_melt$Jaccard_Distance <= dist_cutoff),]
      jacc_df_melt <- jacc_df_melt[which(jacc_df_melt$Pathway_A != jacc_df_melt$Pathway_B),]
      jacc_df_melt <- jacc_df_melt[!duplicated(t(apply(jacc_df_melt,1,sort))),]
      jacc_df_melt <- SIF_tab_PAthConn()
      readr::write_delim(jacc_df_melt,file,delim = '\t')
      
      
    }
  )
  
  ## Download Jaccard output - Intra
  output$dnldJaccTabIntra <- downloadHandler(
    filename = function() {
      
      gs.u <- input$UserPathwayFile
      TopFeat <- input$TopFeatureSelect
      file_name <- gs.u$name
      file_name <- tools::file_path_sans_ext(gs.u)
      paste(file_name,"_Top",TopFeat,"_Jaccard_Connectivity.txt", sep = '')
      
    },
    content = function(file) {
      
      jacc_df <- jacc_react_Intra()
      jacc_df <- as.data.frame(jacc_df)
      jacc_df$Pathways <- rownames(jacc_df)
      jacc_df <- jacc_df %>%
        relocate(Pathways)
      readr::write_delim(jacc_df,file,delim = '\t')
      
    }
  )
  
  ## Download Jaccard output - Intra
  output$dnldClusterTabAnno <- downloadHandler(
    filename = function() {
      
      gs.u <- input$UserPathwayFile
      TopFeat <- input$TopFeatureSelect
      file_name <- gs.u$name
      file_name <- tools::file_path_sans_ext(gs.u)
      paste(file_name,"_Top",TopFeat,"_Cluster_Annotation.txt", sep = '')
      
    },
    content = function(file) {
      
      df <- ClusterTabAnno_react()
      readr::write_delim(df,file,delim = '\t')
      
    }
  )
  
  ####----Matrix Clustering----####
  
  output$dnldUserMatrixHeatmapPDF <- downloadHandler(
    filename = function() {
      gs.u <- input$UserMatrixFile
      file_name <- gs.u$name
      file_name <- tools::file_path_sans_ext(gs.u)
      paste(file_name,"_Heatmap.pdf", sep = '')
    },
    content = function(file) {
      heat <- UserMatrixHeatmap_react()
      ggsave(file, heat, width = 10, height = 12)
    }
  )
  output$dnldUserMatrixHeatmapSVG <- downloadHandler(
    filename = function() {
      gs.u <- input$UserMatrixFile
      file_name <- gs.u$name
      file_name <- tools::file_path_sans_ext(gs.u)
      paste(file_name,"_Heatmap.svg", sep = '')
    },
    content = function(file) {
      heat <- UserMatrixHeatmap_react()
      ggsave(file, heat, width = 10, height = 12)
    }
  )
  output$dnldUserMatixCLusteringPDF <- downloadHandler(
    filename = function() {
      gs.u <- input$UserMatrixFile
      file_name <- gs.u$name
      file_name <- tools::file_path_sans_ext(gs.u)
      paste(file_name,"_ClusterFigure.pdf", sep = '')
    },
    content = function(file) {
      dendo <- JaccDendoIntra_react()
      ggsave(file, dendo, width = 10, height = 8)
    }
  )
  output$dnldUserMatixCLusteringSVG <- downloadHandler(
    filename = function() {
      gs.u <- input$UserMatrixFile
      file_name <- gs.u$name
      file_name <- tools::file_path_sans_ext(gs.u)
      paste(file_name,"_ClusterFigure.svg", sep = '')
    },
    content = function(file) {
      dendo <- JaccDendoIntra_react()
      ggsave(file, dendo, width = 10, height = 8)
    }
  )
  ## Download Cluster Table - INTRA
  output$dnldClustTabMat <- downloadHandler(
    filename = function() {
      
      gs.u <- input$UserMatrixFile
      file_name <- gs.u$name
      file_name <- tools::file_path_sans_ext(gs.u)
      paste(file_name,"_Cluster_Result.txt", sep = '')
      
    },
    content = function(file) {
      
      ## Data matrix
      df <- user_matrix_loaded()
      rownames(df) <- df[,1]
      df <- df[,-1]
      df_mat <- as.matrix(df)
      clust_me <- input$ClustMethodMat
      
      hm_res <- pheatmap::pheatmap(df_mat,
                                   clustering_method = clust_me,
                                   silent = T)
      
      cut_k <- input$NumClustersMat
      clustTab <- data.frame(Pathways = rownames(df),
                             Cluster = cutree(hm_res$tree_row,k = cut_k))
      readr::write_delim(clustTab,file,delim = '\t')
      
    }
  )
  
  output$dnldSIFTabMat <- downloadHandler(
    filename = function() {
      
      gs.u <- input$UserMatrixFile
      file_name <- gs.u$name
      file_name <- tools::file_path_sans_ext(gs.u)
      paste(file_name,"_sif_dend.sif", sep = '')
      
    },
    content = function(file) {
      
      df <- user_matrix_loaded()
      rownames(df) <- df[,1]
      df <- df[,-1]
      df <- as.data.frame(df)
      df_melt <- melt(df)
      df_cols <- colnames(df)
      df_melt$Sample_B <- rep_len(df_cols,length.out = nrow(df_melt))
      colnames(df_melt)[c(1,2)] <- c("Sample_A","Distance")
      df_melt <- df_melt[which(df_melt$Sample_A != df_melt$Sample_B),]
      df_melt <- df_melt[!duplicated(t(apply(df_melt,1,sort))),]
      readr::write_delim(df_melt,file,delim = '\t')
      
      
    }
  )
  
}


shinyApp(ui,server)

















