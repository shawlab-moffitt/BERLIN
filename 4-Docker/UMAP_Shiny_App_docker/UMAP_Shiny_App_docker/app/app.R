

##################################
##                              ##
## Pre Calculate UMAP Shiny App ##
##                              ##
##################################



####----User File Input----####

Project_Name <- 'GSE116256 AMLscRNA'

Expression_File <- 'Example_Input_Data/GSE116256_AMLscRNA_1000_RNA_normalized_counts.txt'

Meta_File <- 'Example_Input_Data/GSE116256_AMLscRNA_1000_RNA_metafile.txt'

PreSelect_UMAP1 <- 'UMAP_1'
PreSelect_UMAP2 <- 'UMAP_2'
PreSelect_Annotation1 <- 'seurat_clusters'
PreSelect_Annotation2 <- 'CellType'


GeneSet_File <- 'GeneSet_Data/GeneSet_List_HS_v5.RData'
GeneSet_CatTab_File <- 'GeneSet_Data/GeneSet_CatTable_v5.zip'





####----Install and load packages----####

packages <- c("shiny","shinythemes","shinyjqui","pheatmap","RColorBrewer","umap","shinyjs","cowplot","patchwork",
              "ggdendro","factoextra","dplyr","DT","viridis","readr","tidyverse","ggrepel","ggVennDiagram","ggtree",
              "shinycssloaders","stringr","tools","plotly","reshape2","ggpubr","gridExtra","scales")

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
invisible(lapply(packages, library, character.only = TRUE))
#bioconductor packages
bioCpacks <- c("clusterProfiler","GSVA","slingshot","SingleCellExperiment")
installed_packages_BIOC <- bioCpacks %in% rownames(installed.packages())
if (any(installed_packages_BIOC == FALSE)) {
  BiocManager::install(bioCpacks[!installed_packages_BIOC], ask = F)
}
invisible(lapply(bioCpacks, library, character.only = TRUE))


#increase file upload size
options(shiny.maxRequestSize=5000*1024^2)

####----Read In Files----####

##--Meta Data--##

meta_ext <- tools::file_ext(Meta_File)
if (meta_ext == "csv") {
  metaIn <- as.data.frame(read_delim(Meta_File, delim = ',', col_names = T))
}
if (meta_ext != "csv") {
  metaIn <- as.data.frame(read_delim(Meta_File, delim = '\t', col_names = T))
}
colnames(metaIn)[1] <- "SampleName"
metaIn[,1] <- gsub("[[:punct:]]","_",metaIn[,1])
# Get variables for selection feature
SamplesToSelect <- c("Select All Samples",metaIn[,1])
anno_options_og <- colnames(metaIn)[2:ncol(metaIn)]
anno_options <- anno_options_og

##--Expression Data--##

# Check extension and read file
expr_ext <- tools::file_ext(Expression_File)
if (expr_ext == "csv") {
  expr <- as.data.frame(read_delim(Expression_File, delim = ',', col_names = T))
}
if (expr_ext != "csv") {
  expr <- as.data.frame(read_delim(Expression_File, delim = '\t', col_names = T))
}
colnames(expr)[-1] <- gsub("[[:punct:]]","_",colnames(expr)[-1])
# Check that expression columns are numeric
isChar <- unname(which(sapply(expr, function(x) is.character(x))))
isChar <-  isChar[-1]
if (length(isChar) > 0) {
  expr[isChar] <- sapply(expr[isChar],as.numeric)
}
# Check for duplicated symbols in first column
colnames(expr)[1] <- "Symbol"
if (TRUE %in% duplicated(expr[,1])) {
  expr <- expr %>%
    group_by(Symbol) %>%
    summarise_all(max)
}
# Make symbols rownames
expr <- as.data.frame(expr)
rownames(expr) <- expr[,1]
expr <- expr[,-1]

# Get variables for selection feature
GeneSymbols <- rownames(expr)

## Make sure Samples are in both expr and meta
sampSame <- intersect(colnames(expr),metaIn[,1])
expr <- expr[,sampSame]
metaIn <- metaIn[which(metaIn[,1] %in% sampSame),]

expr <- expr[,metaIn$SampleName]

expr_mat <- as.matrix(expr)
sce <- SingleCellExperiment(assays = List(counts = expr_mat))
assays(sce)$norm <- assays(sce)$counts

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
geneset <- loadRData(GeneSet_File)
GeneSet_CatTab <- as.data.frame(read_delim(GeneSet_CatTab_File,delim = '\t', col_names = T))
GeneSet_Cats <- unique(GeneSet_CatTab[,1])


Project_Name2 <- gsub("[[:punct:]]",".",Project_Name)
Project_Name2 <- gsub(" ",".",Project_Name)


if (is.null(PreSelect_UMAP1) || PreSelect_UMAP1 == "") {
  PreSelect_UMAP1 <- colnames(metaIn)[2]
}
if (is.null(PreSelect_UMAP2) || PreSelect_UMAP2 == "") {
  PreSelect_UMAP2 <- colnames(metaIn)[3]
}
if (is.null(PreSelect_Annotation1) || PreSelect_Annotation1 == "") {
  PreSelect_Annotation1 <- NULL
}
if (is.null(PreSelect_Annotation2) || PreSelect_Annotation2 == "") {
  PreSelect_Annotation2 <- NULL
}



####----Functions----####


g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

## High-Low
highlow = function(mat) {
  new_mat = mat;
  new_mat[mat > quantile(as.numeric(mat), na.rm = T)[3]] = "High";
  new_mat[mat <= quantile(as.numeric(mat), na.rm = T)[3]] = "Low";
  return (new_mat)
}

#slingMST_me <- function(x, as.df = FALSE) {
#  
#  if(!as.df){
#    return(metadata(x)$mst)
#  }else{
#    dfs <- lapply(seq_along(metadata(x)$lineages), function(l){
#      lin <- metadata(x)$lineages[[l]]
#      mst <- metadata(x)$mst
#      centers <- do.call(rbind, V(mst)$coordinates)
#      rownames(centers) <- V(mst)$name
#      return(data.frame(centers[lin,], Order = seq_along(lin), 
#                        Lineage = l, Cluster = lin))
#    })
#    return(do.call(rbind, dfs))
#  }
#  
#}



####----UI----####

ui <- 
  navbarPage(paste0("{ ",Project_Name," UMAP }"),
             
             tabPanel("UMAP",
                      fluidPage(
                        title = "UMAP",
                        sidebarLayout(
                          sidebarPanel(width = 3,
                                       tabsetPanel(
                                         id = "mainSidebar",
                                         
                                         ####----Sidebar Data Param Panel----####
                                         
                                         tabPanel("Data Parameters",
                                                  p(),
                                                  h4("Plot Data"),
                                                  fluidRow(
                                                    column(6,
                                                           selectInput("SelectPreCalc1","Select X-Axis Coordinate:",
                                                                       choices = colnames(metaIn)[2:ncol(metaIn)], selected = PreSelect_UMAP1)
                                                    ),
                                                    column(6,
                                                           selectInput("SelectPreCalc2","Select Y-Axis Coordinate:",
                                                                       choices = colnames(metaIn)[2:ncol(metaIn)], selected = PreSelect_UMAP2)
                                                    )
                                                  ),
                                                  fluidRow(
                                                    column(5, style = "margin-top:-15px;padding-right:2px",
                                                           checkboxInput("AddTrajLines","Add Trajectory", value = F)
                                                    ),
                                                    column(7, style = "margin-top:-24px;padding-left:2px;",
                                                           uiOutput("rendCurveOrLin")
                                                    )
                                                  ),
                                                  fluidRow(
                                                    column(6,
                                                           uiOutput("rendTrajClustSelect")
                                                    ),
                                                    column(6,
                                                           uiOutput("rendTrajClustStart")
                                                    )
                                                  ),
                                                  checkboxInput("AnnoPathway","Annotate Pathway Score", value = F),
                                                  uiOutput("rendPathTabsetTab"),
                                                  uiOutput("rendViewGeneSetGenes"),
                                                  div(DT::dataTableOutput("GenesInGeneSetTab"), style = "font-size:10px"),
                                                  hr(),
                                                  h4("Annotate Samples"),
                                                  uiOutput("rendUMAPsampSelect"),
                                                  value = 1),
                                         
                                         ####----Sidebar Figure Param Panel----####
                                         
                                         tabPanel("Figure Parameters",
                                                  p(),
                                                  selectInput("UMAPcolors","Continuous Variable Color Palette:",
                                                              choices = c("Red" = "Reds","Blue" = "Blues","Purple" = "Purples",
                                                                          "Green" = "Greens","Orange" = "Oranges","Grey" = "Greys",
                                                                          "Blue/Green" = "BuGn","Yellow/Green" = "YlGn","Red/Purple" = "RdPu")),
                                                  selectInput("UMAPorientation","UMAP Plot Trio Orientation",choices = c("Side-by-Side","Stacked")),
                                                  fluidRow(
                                                    column(6,
                                                           textInput("UMAPplotHeight","UMAP Plot Height:",value = "500px",
                                                                     placeholder = "e.g. '500px','auto','100%'")
                                                    ),
                                                    column(6,
                                                           textInput("UMAPplotWidth","UMAP Plot Width:",value = "100%",
                                                                     placeholder = "e.g. '450px','auto','100%'")
                                                    )
                                                  ),
                                                  hr(),
                                                  h4("Font Sizes"),
                                                  fluidRow(
                                                    column(3,
                                                           numericInput("UMAPtitleTextSize","Title:", value = 16)
                                                    ),
                                                    column(3,
                                                           numericInput("UMAPaxisTextSize","Axis:", value = 12)
                                                    ),
                                                    column(3,
                                                           numericInput("UMAPlegendTextSize","Legend:", value = 11)
                                                    ),
                                                    column(3,
                                                           numericInput("UMAPdotSize","Dot size:", value = 1)
                                                    )
                                                  ),
                                                  hr(),
                                                  fluidRow(
                                                    column(6, style = 'padding-right:2px;',
                                                           textInput("UMAPxAxisLim","X-Axis Limits: min,max",value = "")
                                                    ),
                                                    column(6, style = 'padding-right:2px;padding-left:2px;',
                                                           textInput("UMAPyAxisLim","Y-Axis Limits: min,max",value = "")
                                                    )
                                                  ),
                                                  value = 2)
                                       )
                          ),
                          mainPanel(
                            p(),
                            fluidRow(
                              column(4, style = 'padding-left:90px;',
                                     uiOutput("rendUMAPannotateSamps1"),
                                     style = 'border-right: 0.5px solid lightgray',
                              ),
                              column(4, style = 'padding-left:60px;',
                                     fluidRow(
                                       column(5, #style = 'padding-right:2px;',
                                              uiOutput("rendGeneSelection")
                                       ),
                                       column(2, #style = 'padding-right:0px;padding-left:2px;',
                                              checkboxInput("LogGeneSelection","Log2",value = F)
                                       ),
                                       column(5,
                                              uiOutput("rendGeneExprRange")
                                       )
                                     ),
                                     style = 'border-right: 0.5px solid lightgray',
                              ),
                              column(4, style = 'padding-left:70px;',
                                     uiOutput("rendUMAPannotateSamps2")
                              )
                            ),
                            p(),
                            uiOutput("rendUMAPplot_PreC_ALL"),
                            jqui_draggable(jqui_resizable(plotOutput('UMAPprecLegend', width = "100%", height = "110px"))),
                            p(),
                            downloadButton("dnldUMAP_SVG_PreC_clin1","Download Clinical Annotation UMAP 1"),
                            downloadButton("dnldUMAP_SVG_PreC_expr","Download Gene Expression UMAP"),
                            downloadButton("dnldUMAP_SVG_PreC_clin2","Download Clinical Annotation UMAP 2"),
                            p(),
                            div(DT::dataTableOutput("UMAP_PreC_CoordTable"), style = "font-size:12px")
                          )
                        )
                        #)
                        #)
                      ) 
             ),
             
             ####----Bar Plot Tab----####
             
             tabPanel("Feature Comparison",
                      fluidPage(
                        title = "Feature Comparison",
                        sidebarLayout(
                          sidebarPanel(
                            width = 3,
                            tabsetPanel(
                              tabPanel("Data Parameters",
                                       p(),
                                       uiOutput("rendBPsampSubset"),
                                       uiOutput("rendBPsampCriteria"),
                                       uiOutput("rendBPgroupCriteria"),
                                       fluidRow(
                                         column(8,
                                                uiOutput("rendBPgeneSelection")
                                         ),
                                         column(4,
                                                checkboxInput("log2barplot","Log2",value = F)
                                         )
                                       ),
                                       fluidRow(
                                         uiOutput("rendBPstatComp"),
                                         uiOutput("rendVplotsampledots"),
                                         uiOutput("rendErrorBP"),
                                         uiOutput("rendVplotDotSize")
                                       ),
                                       radioButtons("ViolinOrBoxP","View As:",choices = c("Violin Plot","Box Plot","Barplot","Stacked Barplot"), inline = T),
                                       uiOutput("rendBPfeatFill")
                              ),
                              tabPanel("Figure Parameters",
                                       p(),
                                       textInput("barplotColoCodes","Color Code(s):",value = "", placeholder = "HEX or R Color Code(s) (Space Delim)"),
                                       hr(),
                                       fluidRow(
                                         column(6,
                                                checkboxInput("barplotsampledots","Include Dot Annotation", value = F)
                                         ),
                                         column(6,
                                                numericInput("barplotDotSize","Dot Size:", value = 0.25, step = 0.25)
                                         )
                                       ),
                                       hr(),
                                       fluidRow(
                                         column(6,
                                                textInput("barPlotYlim","Y-Axis Limits", value = "", placeholder = "min,max"),
                                                selectInput("barxAxisOrient","X-Axis Label Orientation",
                                                            choices = c(45,90,0))
                                         ),
                                         column(6,
                                                numericInput("barplotYbreaks","Y-Axis Breaks",value = "", step = 0.5),
                                                selectInput("barplotXaxOrder","X-Axis Group Order",
                                                            choices = c("Not Specificed","Ascending","Descending"))
                                         )
                                       ),
                                       hr(),
                                       fluidRow(
                                         column(6,
                                                numericInput("barplot1TitleSize","Title Font Size:",
                                                             value = 20, step = 1)
                                         ),
                                         column(6,
                                                numericInput("barplot1AxisSize","Axis Font Size:",
                                                             value = 16, step = 1)
                                         )
                                       )
                              )
                            )
                          ),
                          mainPanel(
                            p(),
                            jqui_resizable(plotOutput('barplot', width = "100%", height = "600px")),
                            fluidRow(
                              downloadButton("dnldPlotSVG_Barplot","Download as SVG"),
                              downloadButton("dnldPlotPDF_Barplot","Download as PDF"),
                              p(),
                              div(DT::dataTableOutput("barplot_table"), style = "font-size:12px")
                            )
                          )
                        )
                      )
             ),
             
             ####----Enrichment test----####
             
             tabPanel("Enrichment Test",
                      fluidPage(
                        title = "Enrichment Test",
                        sidebarLayout(
                          sidebarPanel(
                            width = 3,
                            ## Allows all select inputs to be wide enough to read the contents
                            tags$head(
                              tags$style(HTML('.selectize-input {
                                                white-space: nowrap;
                                              }
                                              .selectize-dropdown {
                                                width: 400px !important;
                                              }'
                              )
                              )
                            ),
                            fluidRow(
                              column(6,
                                     selectInput("EnrichFeat1","Feature One:", choices = anno_options, selected = PreSelect_Annotation1),
                                     uiOutput("rendEnrichVar1")
                              ),
                              column(6,
                                     selectInput("EnrichFeat2","Feature Two:", choices = anno_options, selected = PreSelect_Annotation2),
                                     uiOutput("rendEnrichVar2")
                              )
                            ),
                            radioButtons("FisherTailChoice","",choices = c("Two-Tailed" = "two.sided","One-Tailed - Greater" = "greater","One-Tailed - Less" = "less"), inline = T)
                          ),
                          mainPanel(
                            withSpinner(jqui_resizable(plotOutput("EnrichVennPlot", height = "450px", width = "100%"))),
                            hr(),
                            htmlOutput("EnrichFisherText", style = "font-size:14px;"),
                            hr(),
                            fluidRow(
                              column(6,
                                     tableOutput("EnrichFisherTab")
                              ),
                              column(6,
                                     verbatimTextOutput("EnrichFisherOut")
                              )
                            )
                          )
                        )
                      )
             ),
             
             ####----Heatmap----####
             
             tabPanel("Heatmap",
                      sidebarLayout(
                        sidebarPanel(
                          width = 3,
                          p(),
                          uiOutput("rendHeatSampSubset"),
                          uiOutput("rendHeatSampCriteria"),
                          uiOutput("rendHeatGroupCriteria"),
                          uiOutput("rendHeatGroupSelect"),
                          uiOutput("rendHeatGeneSelect"),
                          fluidRow(
                            column(4,
                                   selectInput("ClusterMethodChoice","Cluster Method:",
                                               choices = c("ward.D", "ward.D2", "complete", "single", "average", "mcquitty", "median", "centroid"))
                                   ),
                            column(4,
                                   checkboxInput("ClusterCol","Cluster Cols", value = TRUE),
                                   checkboxInput("ClusterRow","Cluster Rows", value = TRUE)
                                   ),
                            column(4,
                                   checkboxInput("DotHeatChoice","As Dot Heatmap", value = F)
                                   )
                            )
                          #fluidRow(
                          #  column(4,
                          #         ),
                          #  column(8,
                          #         uiOutput("rendDotHeatSizeFactor")
                          #         )
                          #)
                          ),
                        mainPanel(
                          p(),
                          withSpinner(jqui_resizable(plotOutput("HeatmapAvg", height = "800px",width = "100%")), type = 6)
                        )
                      ))
  )


####----Server----####


server <- function(input, output, session) {
  
  ####----Render UI----####
  
  
  output$rendUMAPsampSelect <- renderUI({
    
    selectizeInput(
      "UMAPsampSelect", 
      label = "Select Samples:",
      choices = SamplesToSelect, 
      multiple = T,
      selected = "",
      options = list(delimiter = " ", create = T)
    )
    
  })
  
  output$rendTrajClustSelect <- renderUI({
    
    if (input$AddTrajLines) {
      selectInput("TrajClustSelect","Trajectory Clusters:", choices = anno_options, multiple = F, selected = PreSelect_Annotation1)
    }
    
  })
  
  output$rendCurveOrLin <- renderUI({
    
    if (input$AddTrajLines) {
      radioButtons("CurveOrLin","", choices = c("Curves","Min Spanning Tree"), inline = T)
    }
    
  })
  
  output$rendTrajClustStart <- renderUI({
    
    #meta <- meta_react()
    #print(meta[1:5,1:5])
    if (input$AddTrajLines) {
      if (!is.null(input$TrajClustSelect)) {
        StartChoices <- c(" ",unique(metaIn[,input$TrajClustSelect]))
        selectInput("TrajClustStart","Starting Cluster:", choices = StartChoices, multiple = F)
      }
    }
    
  })
  
  
  meta_react <- reactive({
    
    meta_update <- metaIn
    if (input$AddTrajLines) {
      sce_umap <- SlingShot_react()
      pseudo_df <- as.data.frame(slingPseudotime(sce_umap))
      pseudo_df$Pseudotime_Lineage_Avg <- rowMeans(pseudo_df, na.rm = T)
      pseudo_df$SampleName <- rownames(pseudo_df)
      meta_update <- merge(pseudo_df[,c(ncol(pseudo_df)-1,ncol(pseudo_df))],metaIn, by = "SampleName", all.x = T)
    }
    if (input$AnnoPathway) {
      gsSelected <- GeneSetSelect_react()
      if (length(names(gsSelected)) > 0) {
          if (input$UMAPannotateSamps1 == names(gsSelected) | input$UMAPannotateSamps2 == names(gsSelected)) {
            #if (input$UMAPannotateSamps1 == "Selected Geneset" | input$UMAPannotateSamps2 == "Selected Geneset") {
            if (length(unname(unlist(gsSelected)) %in% rownames(expr_mat)) > 0) {
              ssGSEA <- gsva(expr_mat,gsSelected,method = "ssgsea", verbose = FALSE)
              ssGSEA <- as.data.frame(t(ssGSEA))
              ssGSEA$SampleName <- rownames(ssGSEA)
              #colnames(ssGSEA)[1] <- "Selected Geneset"
              meta_update <- merge(ssGSEA,metaIn, by = "SampleName", all = T)
            }
          }
      }
    }
    meta_update
    
  })
  
  Anno1Select <- reactiveValues(Annotation1Selected = PreSelect_Annotation1)
  observeEvent(input$UMAPannotateSamps1, {
    Anno1Select$Annotation1Selected <- input$UMAPannotateSamps1
  })
  output$rendUMAPannotateSamps1 <- renderUI({
    
    if (input$AddTrajLines) {
      anno_options <- c("Pseudotime_Lineage_Avg",anno_options)
    }
    if (input$AnnoPathway) {
      gsSelected <- GeneSetSelect_react()
      geneset_name <- names(gsSelected)
      anno_options <- c(geneset_name,anno_options)
      #anno_options <- c("Selected Geneset",anno_options)
    }
    selectInput("UMAPannotateSamps1","Annotate Samples By:", choices = anno_options, multiple = F, selected = Anno1Select$Annotation1Selected)
    
  })
  Anno2Select <- reactiveValues(Annotation2Selected = PreSelect_Annotation2)
  observeEvent(input$UMAPannotateSamps2, {
    Anno2Select$Annotation2Selected <- input$UMAPannotateSamps2
  })
  output$rendUMAPannotateSamps2 <- renderUI({
    
    if (input$AddTrajLines) {
      anno_options <- c("Pseudotime_Lineage_Avg",anno_options)
    }
    if (input$AnnoPathway) {
      gsSelected <- GeneSetSelect_react()
      geneset_name <- names(gsSelected)
      anno_options <- c(geneset_name,anno_options)
    }
    selectInput("UMAPannotateSamps2","Annotate Samples By:", choices = anno_options, multiple = F, selected = Anno2Select$Annotation2Selected)
    
  })
  
  DiscContList_react <- reactive({
    
    ## Check which columns of covs are discrete or continuous
    meta <- meta_react()
    covs <- meta
    covsClass <- sapply(as_tibble(covs), typeof)
    covsClass_num <- names(covsClass[which(covsClass %in% c("integer","double"))])
    IntCovsaCols <- apply(covs[,covsClass_num],2,function(x) any(round(as.numeric(x)) != x))
    MaybeCont <- names(IntCovsaCols[which(IntCovsaCols == FALSE)])
    MaybeContData <- apply(covs[,MaybeCont],2,function(x) any(length(levels(as.factor(x)))<(nrow(covs)*0.60)))
    ContinuousData <- names(MaybeContData[which(MaybeContData == FALSE)])
    ContinuousData <- c(ContinuousData,names(IntCovsaCols[which(IntCovsaCols == TRUE)]))
    
    covs.continuous <- covs[,which(colnames(covs) %in% ContinuousData)]
    covs.discrete <- covs[,which(!colnames(covs) %in% ContinuousData)]
    
    DiscContList <- list(discreteCovs=colnames(covs.discrete), continuousCovs=colnames(covs.continuous))
    DiscContList
    
  })
  
  output$rendGeneSelection <- renderUI({
    
    selectizeInput("GeneSelection","Feature:", choices = GeneSymbols)
    
  })
  
  output$rendGeneExprRange <- renderUI({
    
    geneSelec <- input$GeneSelection
    LogChoice <- input$LogGeneSelection
    if (LogChoice == TRUE) {
      expr <- log2(as.matrix(expr) + 1)
    }
    exprVal <- expr[geneSelec,][!is.na(expr[geneSelec,])]
    #exprRange <- round(range(expr[geneSelec,], na.rm = T),2)
    exprRange <- round(range(exprVal, na.rm = T),2)
    exprRangeText <- paste(exprRange[1],",",exprRange[2],sep = "")
    textInput("GeneExprRange","Range:",value = exprRangeText, placeholder = "min,max")
    
  })
  
  output$rendUMAPplot_PreC_ALL <- renderUI({
    
    ph <- input$UMAPplotHeight
    pw <- input$UMAPplotWidth
    withSpinner(jqui_resizable(plotlyOutput('UMAPplot_PreC_ALL', width = pw, height = ph)), type = 6)
    
  })
  
  output$rendBPsampSubset <- renderUI({
    
    anno_options2 <- c("Select All",anno_options_og)
    selectInput("BPsampSubset","Subset Samples By:", choices = anno_options2)
    
  })
  
  output$rendBPsampCriteria <- renderUI({
    
    meta <- meta_react()
    SampSubset <- input$BPsampSubset
    if (!is.null(SampSubset)) {
      if (SampSubset != "Select All") {
        SubsetOptions <- unique(meta[,SampSubset])
        selectInput("BPsampCriteria",paste(SampSubset,"Criteria:"), choices = SubsetOptions)
      }
    }
    
  })
  
  output$rendBPgroupCriteria <- renderUI({
    
    DiscContList <- DiscContList_react()
    selectInput("BPgroupCriteria","Bar Plot Grouping Criteria:", choices = DiscContList$discreteCovs, selected = "seurat_clusters")
    
  })
  
  output$rendBPgeneSelection <- renderUI({
    
    meta <- meta_react()
    FeatureChoice <- c(GeneSymbols,colnames(meta)[-1])
    selectizeInput("BPgeneSelection","Select Feature:", choices = FeatureChoice)
    
  })
  
  #output$rendbarplotsampledots <- renderUI({
  #  
  #  SampSubset <- input$BPsampSubset
  #  if (!is.null(SampSubset)) {
  #    if (SampSubset == "Select All") {
  #      checkboxInput("barplotsampledots","Include Dot Annotation", value = F)
  #    }
  #    else if (SampSubset != "Select All") {
  #      checkboxInput("barplotsampledots","Include Dot Annotation", value = T)
  #    }
  #  }
  #  
  #})
  
  output$rendBPstatComp <- renderUI({
    
    if (input$ViolinOrBoxP %in% c("Violin Plot","Box Plot","Barplot")) {
      column(5,
             selectInput("BPstatComp","Stat Test Method:",
                         choices = c("none","wilcox.test","t.test","kruskal.test","anova"))
      )
    }
    
  })
  
  output$rendErrorBP <- renderUI({
    
    if (input$ViolinOrBoxP == "Barplot") {
      column(6,
             selectInput("ErrorBP","Error Bar Type",
                         choices = c("Standard Deviation","Standard Error","None"), width = "150%")
      )
    }
    
  })
  
  output$rendVplotsampledots <- renderUI({
    
    if (!input$ViolinOrBoxP %in% c("Barplot","Stacked Barplot")) {
      column(4, 
             checkboxInput("Vplotsampledots","Include Dot Annotation", value = F)
      )
    }
    
  })
  
  output$rendVplotDotSize <- renderUI({
    
    if (!input$ViolinOrBoxP %in% c("Barplot","Stacked Barplot")) {
      column(3,
             numericInput("VplotDotSize","Dot Size:", value =1, step = 0.25)
      )
      
    }
    
  })
  
  output$rendBPfeatFill <- renderUI({
    
    if (input$ViolinOrBoxP == "Stacked Barplot") {
      
      checkboxInput("BPfeatFill","View Barplot Fill as Percentage", value = F)
      
    }
    
  })
  
  output$rendEnrichVar1 <- renderUI({
    
    meta <- meta_react()
    Feat1 <- input$EnrichFeat1
    VarChoices <- unique(meta[,Feat1])
    selectInput("EnrichVar1", "Variable from Feature One:", choices = VarChoices, selected = VarChoices[1])
    
  })
  
  output$rendEnrichVar2 <- renderUI({
    
    meta <- meta_react()
    Feat2 <- input$EnrichFeat2
    VarChoices <- unique(meta[,Feat2])
    selectInput("EnrichVar2", "Variable from Feature Two:", choices = VarChoices, selected = VarChoices[1])
    
  })
  
  output$rendPathTabsetTab <- renderUI({
    
    if (input$AnnoPathway) {
      
      tabsetPanel(
        id = "GeneSetPanel",
        tabPanel("Provided Genesets",
                 p(),
                 selectInput("gsDatabase","Select Geneset Database:", choices = GeneSet_Cats),
                 div(DT::dataTableOutput("GeneSetTable"), style = "font-size:10px"),
                 value = 1
                 ),
        tabPanel("User Genesets",
                 p(),
                 fileInput("userGeneSet","Upload Geneset"),
                 div(DT::dataTableOutput("userGeneSetTable"), style = "font-size:10px"),
                 value = 2
                 )
      )
      
    }
    
  })
  
  output$rendViewGeneSetGenes <- renderUI({
    
    if (input$AnnoPathway) {
      checkboxInput("ViewGeneSetGenes","View genes in selected geneset", value = FALSE)
    }
    
  })
  
  GeneSetTable_React <- reactive({
    
    GS_database <- input$gsDatabase
    sub_tab <- GeneSet_CatTab[which(GeneSet_CatTab[,1] == GS_database),]
    new_tab <- sub_tab[,-c(1,2)]
    new_tab
    
  })
  
  ## Render Gene Set Selection Table
  output$GeneSetTable <- renderDataTable({
    
    GeneSetTable_sub <- GeneSetTable_React()
    DT::datatable(GeneSetTable_sub,
                  options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                 pageLength = 10,
                                 scrollX = T),
                  selection = list(mode = 'single', selected = 1),
                  rownames = F)
    
  })
  
  userGeneSet_react <- reactive({
    
    gs.u <- input$userGeneSet
    ext <- tools::file_ext(gs.u$datapath)
    req(gs.u)
    
    # If user provides GMT file
    if (ext == "gmt") {
      gs_u <- clusterProfiler::read.gmt(gs.u$datapath)
      colnames(gs_u) <- c("term","gene")
      gs_u
    } else if (ext == "RData") {
      gs_u <- loadRData(gs.u$datapath)
      gs_u
    } else {
      gs_u <- as.data.frame(readr::read_delim(gs.u$datapath, delim = '\t'))
      colnames(gs_u) <- c("term","gene")
      gs_u
    }
    
  })
  
  userGeneSet_select <- reactive({
    
    gs_u <- userGeneSet_react()
    gs_tab <- userGeneSetTable_react()
    geneset_select <- gs_tab[input$userGeneSetTable_rows_selected,1]
    gsU_name <- input$userGeneSet
    ext <- tools::file_ext(gsU_name$datapath)
    
    if (ext == "RData") {
      gs_u_select <- gs_u[geneset_select]
      print(head(gs_u_select))
      gs_u_select
    } else {
      tab_select <- gs_u[which(gs_u[,1] == geneset_select),]
      gs_u_select <- list()
      for (i in unique(tab_select[,1])){
        gs_u_select[[i]] <- tab_select[tab_select[,1] == i,]$gene
      }
      gs_u_select
    }
    
  })
  
  ## Render User Gene Set Table - Backend
  userGeneSetTable_react <- reactive({
    
    fileIn <- input$userGeneSet
    ext <- tools::file_ext(fileIn$datapath)
    gs_u <- userGeneSet_react()
    
    # If user provides GMT file
    if (ext == "gmt") {
      uGS_table <- as.data.frame(unique(gs_u[,1]))
      colnames(uGS_table)[1] <- "GeneSet"
    } else if (ext == "RData") {
      uGS_table <- as.data.frame(names(gs_u))
      colnames(uGS_table)[1] <- "GeneSet"
    } else {
      uGS_table <- as.data.frame(unique(gs_u[,1]))
      colnames(uGS_table)[1] <- "GeneSet"
    }
    uGS_table
    
  })
  
  ## Render User Gene Set Table
  output$userGeneSetTable <- renderDataTable({
    
    req(input$userGeneSet)
    uGS_table <- userGeneSetTable_react()
    DT::datatable(uGS_table,
                  options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                 pageLength = 10,
                                 scrollX = T),
                  selection = list(mode = 'single', selected = 1),
                  rownames = F)
    
    
  })
  
  GeneSetSelect_react <- reactive({
    if (input$AnnoPathway) {
      if (!is.null(input$GeneSetPanel)) {
        if (input$GeneSetPanel == 1) {
          gsTable <- GeneSetTable_React()
          GeneSet_Name <- gsTable[input$GeneSetTable_rows_selected,2]
          GeneSetSelect <- geneset[GeneSet_Name]
          GeneSetSelect
        } else if (input$GeneSetPanel == 2) {
          GeneSetSelect <- userGeneSet_select()
          GeneSetSelect
        }
      }
    }
    
  })
  
  output$rendHeatSampSubset <- renderUI({
    
    anno_options2 <- c("Select All",anno_options_og)
    selectInput("HeatSampSubset","Subset Samples By:", choices = anno_options2)
    
  })
  
  output$rendHeatSampCriteria <- renderUI({
    
    meta <- meta_react()
    SampSubset <- input$HeatSampSubset
    if (!is.null(SampSubset)) {
      if (SampSubset != "Select All") {
        SubsetOptions <- unique(meta[,SampSubset])
        selectInput("HeatSampCriteria",paste(SampSubset,"Criteria:"), choices = SubsetOptions)
      }
    }
    
  })
  
  output$rendHeatGroupCriteria <- renderUI({
    
    DiscContList <- DiscContList_react()
    selectInput("HeatGroupCriteria","Heatmap Grouping Criteria:", choices = DiscContList$discreteCovs, selected = "seurat_clusters")
    
  })
  
  output$rendHeatGroupSelect <- renderUI({
    
    GroupingCol <- input$HeatGroupCriteria
    meta <- meta_react()
    GroupCond <- unique(meta[,GroupingCol])
    selectInput("HeatGroupSelect","Group Coniditon:",
                choices = GroupCond, selected = GroupCond[c(1:length(GroupCond))],
                multiple = T)
    
  })
  
  output$rendHeatGeneSelect <- renderUI({
    
      selectizeInput(
        "HeatGeneSelect", 
        label = "Select or Paste (space delim) Features:",
        choices = GeneSymbols, 
        multiple = T,
        selected = GeneSymbols[1:4],
        options = list(delimiter = " ", create = T)
    )
    
  })
  
  
  ####----Reactives----####
  
  ##--Transform expression data if chosen--##
  exprTrans <- reactive({
    
    expr2 <- expr
    expr3 <- expr2 # save copy for naming columns
    if (input$LogExprMatrix == T) {
      expr2 <- log2(expr2 + 1)
    }
    if (input$NormExprMatrix == T) {
      expr2 = apply(expr2, 1, scale)
      expr2 = apply(expr2, 1, rev)
      colnames(expr2) <- colnames(expr3)
    }
    expr2 <- as.matrix(expr2[sort(rownames(expr2)),])
    
    expr2
    
  })
  
  ## UMAP Coordinate Table  - Pre-Calculated
  umap_plot_table_PreC_react <- reactive({
    
    meta <- meta_react()
    if (!is.null(input$SelectPreCalc1) && !is.null(input$SelectPreCalc2)) {
      umap1 <- input$SelectPreCalc1
      umap2 <- input$SelectPreCalc2
      
      tdata_fit_df <- meta[,c("SampleName",umap1,umap2)]
      if (ncol(tdata_fit_df) >= 3) {
        colnames(tdata_fit_df)[c(2,3)] <- c("UMAP1","UMAP2")
      }
      
      tdata_fit_df
    }
    
  })
  
  UMAP_PreC_CoordTable_react <- reactive({
    
    meta <- meta_react()
    tdata_fit_df <- umap_plot_table_PreC_react()
    
    ## Add Annotation column
    if (!is.null(input$UMAPannotateSamps1)) {
      if (input$UMAPannotateSamps1 != " ") {
        if (input$UMAPannotateSamps1 %in% colnames(meta)) {
          metaColanno1 <- input$UMAPannotateSamps1
          tdata_fit_df <- merge(tdata_fit_df,meta[,c("SampleName",metaColanno1)], by = "SampleName")
          if (metaColanno1 %in% DiscContList_react()$discreteCovs) {
            tdata_fit_df[,metaColanno1] <- as.factor(tdata_fit_df[,metaColanno1])
          }
        }
      }
    }
    
    ## Add Gene Expression Column
    LogChoice <- input$LogGeneSelection
    if (LogChoice == TRUE) {
      expr <- log2(as.matrix(expr) + 1)
    }
    if (!is.null(input$GeneSelection)) {
      metaColgene <- input$GeneSelection
      exprRange <- round(range(expr[metaColgene,], na.rm = T),2)
      expr2 <- as.data.frame(expr)
      expr_g <- expr2[metaColgene,]
      expr_g_t <- as.data.frame(t(expr_g))
      expr_g_t$SampleName <- rownames(expr_g_t)
      tdata_fit_df <- merge(tdata_fit_df,expr_g_t)
    }
    else if (is.null(input$GeneSelection)) {
      metaColgene <- rownames(expr)[1]
      exprRange <- round(range(expr[metaColgene,], na.rm = T),2)
      expr2 <- as.data.frame(expr)
      expr_g <- expr2[metaColgene,]
      expr_g_t <- as.data.frame(t(expr_g))
      expr_g_t$SampleName <- rownames(expr_g_t)
      tdata_fit_df <- merge(tdata_fit_df,expr_g_t)
    }
    
    if (!is.null(input$UMAPannotateSamps2)) {
      if (input$UMAPannotateSamps1 %in% colnames(meta)) {
        if (input$UMAPannotateSamps2 != " ") {
          metaColanno2 <- input$UMAPannotateSamps2
          tdata_fit_df <- merge(tdata_fit_df,meta[,c("SampleName",metaColanno2)], by = "SampleName")
          if (metaColanno2 %in% DiscContList_react()$discreteCovs) {
            tdata_fit_df[,metaColanno2] <- as.factor(tdata_fit_df[,metaColanno2])
          }
        }
      }
    }
    
    tdata_fit_df
    
    
  })
  
  BP_meta_subset <- reactive({
    
    meta <- meta_react()
    if (!is.null(input$BPsampSubset) && !is.null(input$BPgroupCriteria)) {
      if (input$BPgroupCriteria != " ") {
        sampSubset <- input$BPsampSubset
        sampGroupCriteria <- input$BPgroupCriteria
        if (sampSubset == "Select All") {
          metaSub <- meta
          metaSub <- metaSub %>%
            relocate(SampleName,sampGroupCriteria)
          metaSub
        }
        else if (sampSubset != "Select All") {
          sampSubsetCriteria <- input$BPsampCriteria
          metaSub <- meta[which(meta[,sampSubset] == sampSubsetCriteria),]
          metaSub <- metaSub %>%
            relocate(SampleName,sampGroupCriteria,sampSubset)
          metaSub
        }
        
      }
    }
    
  })
  
  ####----Data Tables----####
  
  output$UMAP_PreC_CoordTable <- DT::renderDataTable({
    
    tdata_fit_df <- UMAP_PreC_CoordTable_react()
    if (ncol(tdata_fit_df) >= 3) {
      umap1 <- input$SelectPreCalc1
      umap2 <- input$SelectPreCalc2
      colnames(tdata_fit_df)[c(2,3)] <- c(umap1,umap2)
      #table output,
      DT::datatable(tdata_fit_df,
                    options = list(keys = TRUE,
                                   searchHighlight = TRUE,
                                   pageLength = 20,
                                   scrollX = T,
                                   scrollY = T,
                                   lengthMenu = c("10", "20", "50", "100")
                    ),
                    rownames = F,
                    selection=list(mode = "multiple"))
    }
    
  })
  
  output$barplot_table <- DT::renderDataTable({
    
    #meta <- meta_react()
    if (!is.null(input$BPsampSubset) && !is.null(input$BPgroupCriteria)) {
      if (input$BPgroupCriteria != " ") {
        metaSub <- BP_meta_subset()
        geneSelected <- input$BPgeneSelection
        logchoice <- input$log2barplot
        exprSub <- expr[,metaSub[,1]]
        feature <- geneSelected
        DiscContList <- DiscContList_react()
        
        if (feature %in% rownames(exprSub)) {
          expr_gene <- as.data.frame(t(exprSub[feature,]))
          expr_gene$SampleName <- rownames(expr_gene)
        }
        if (!feature %in% rownames(exprSub)) {
          expr_gene <- as.data.frame(metaSub[,c("SampleName",feature)])
        }
        if (logchoice == T) {
          expr_gene[which(expr_gene[,1] < 0),1] <- 0
          expr_gene[,1] <- log2(expr_gene[,1] + 1)
          colnames(expr_gene)[1] <- paste(colnames(expr_gene)[1],"_Log2")
        }
        metaSub <- merge(expr_gene,metaSub,by = "SampleName", all = T)
        col_num <- ncol(metaSub)
        if (input$BPsampSubset == "Select All") {
          metaSub <- metaSub[,c(1,3,2,4:col_num)]
        }
        else if (input$BPsampSubset != "Select All") {
          metaSub <- metaSub[,c(1,3,4,2,5:col_num)]
        }
        #if (input$ViolinOrBoxP == "Stacked Barplot") {
        #  if (feature %in% DiscContList$continuousCovs | feature %in% rownames(expr)) {
        #    metaSub[,feature] <- highlow(as.numeric(metaSub[, which(colnames(metaSub) == feature)]))
        #  }
        #}
        #table output,
        DT::datatable(metaSub,
                      options = list(keys = TRUE,
                                     searchHighlight = TRUE,
                                     pageLength = 10,
                                     scrollX = T,
                                     lengthMenu = c("10", "20", "50", "100")
                      ),
                      rownames = F,
                      selection=list(mode = "multiple"))
      }
    }
    
    
  })
  
  output$GenesInGeneSetTab <- DT::renderDataTable({
    
    if (input$AnnoPathway) {
      if (!is.null(input$ViewGeneSetGenes)) {
        if (input$ViewGeneSetGenes) {
          geneset <- GeneSetSelect_react()
          gs_df <- as.data.frame(geneset)
          DT::datatable(gs_df, options = list(scrollY = T), rownames = F)
        }
      }
    }
    
    
  })
  
  ####----Plots----####
  
  #SlingShot_react <- reactive({
  #  
  #  umap1 <- input$SelectPreCalc1
  #  umap2 <- input$SelectPreCalc2
  #  TrajCluster <- input$TrajClustSelect
  #  TrajClusterStart <- input$TrajClustStart
  #  if (TrajClusterStart == " ") {
  #    TrajClusterStart <- NULL
  #  }
  #  #sce1 <- sce
  #  rd_umap <- metaIn[,c("SampleName",umap1,umap2)]
  #  rownames(rd_umap) <- rd_umap[,1]
  #  rd_umap <- rd_umap[,-1]
  #  rd_umap <- as.matrix(rd_umap)
  #  reducedDims(sce) <- SimpleList(UMAP = rd_umap)
  #  
  #  cluster_info <- metaIn[,TrajCluster]
  #  names(cluster_info) <- metaIn[,1]
  #  colData(sce)$cluster_info <- cluster_info
  #  
  #  sce_umap <- slingshot(sce, clusterLabels = 'cluster_info', reducedDim = 'UMAP', start.clus = TrajClusterStart)
  #  sce_umap
  #  
  #})
  
  
  SlingShot_react <- reactive({
    
    umap1 <- input$SelectPreCalc1
    umap2 <- input$SelectPreCalc2
    TrajCluster <- input$TrajClustSelect
    TrajClusterStart <- input$TrajClustStart
    if (TrajClusterStart == " ") {
      TrajClusterStart <- NULL
    }
    
    clusterGT1 <- table(metaIn[,TrajCluster])[which(table(metaIn[,TrajCluster]) > 1)]
    meta_sub <- metaIn[which(metaIn[,TrajCluster] %in% names(clusterGT1)),]
    cluster_info <- meta_sub[,TrajCluster]
    names(cluster_info) <- meta_sub[,1]
    
    sce_sub <- sce[,meta_sub$SampleName]
    
    rd_umap <- meta_sub[,c("SampleName",umap1,umap2)]
    rownames(rd_umap) <- rd_umap[,1]
    rd_umap <- rd_umap[,-1]
    rd_umap <- as.matrix(rd_umap)
    reducedDims(sce_sub) <- SimpleList(UMAP = rd_umap)
    
    colData(sce_sub)$cluster_info <- cluster_info
    
    sce_umap <- slingshot(sce_sub, clusterLabels = 'cluster_info', reducedDim = 'UMAP', start.clus = TrajClusterStart)
    sce_umap
    
  })
  
  ####----UMAP Clin1----####
  
  umap_plot_PreC_clin_react_base1 <- reactive({
    
    #meta <- meta_react()
    ## Variables
    sampSelected <- input$UMAPsampSelect         # Sample Names to annotate
    UMAPtitleText <- input$UMAPtitleTextSize     # Title text size
    UMAPaxisText <- input$UMAPaxisTextSize       # Axis text size
    UMAPlegendText <- input$UMAPlegendTextSize   # Legend text size
    UMAPdotSize <- input$UMAPdotSize             # Dot size
    metaColanno1 <- NULL
    metaColanno2 <- NULL
    
    umap1 <- input$SelectPreCalc1
    umap2 <- input$SelectPreCalc2
    
    plot_df <- UMAP_PreC_CoordTable_react()
    rownames(plot_df) <- plot_df[,1]
    if (input$UMAPannotateSamps1 != " ") {
      metaColanno1 <- input$UMAPannotateSamps1
    }
    if (input$UMAPannotateSamps2 != " ") {
      metaColanno2 <- input$UMAPannotateSamps2
    }
    
    if (!is.null(input$GeneSelection)) {
      metaColgene <- input$GeneSelection
    }
    if (is.null(input$GeneSelection)) {
      metaColgene <- rownames(expr)[1]
    }
    
    ## Meta Annotation Plot
    #if (input$UMAPannotateSamps1 != " " && input$UMAPannotateSamps2 != " ") {
      if (length(colnames(plot_df)[c(4,5,6)]) == 3) {
        colnames(plot_df)[c(4,5,6)] <- c("AnnoName1","GeneName","AnnoName2")
        k <- plot_df %>%
          ggplot(aes(UMAP1, UMAP2, colour=AnnoName1,
                     text = paste("</br> <b>Sample Name:</b> ", SampleName,
                                  "</br> <b>",metaColanno1,":</b> ", AnnoName1,
                                  "</br> <b>",metaColgene," Gene Expression:</b> ", round(GeneName,4),
                                  "</br> <b>",metaColanno2,":</b> ", AnnoName2,
                                  sep = "")))
      }
    #  colnames(plot_df)[c(4,5,6)] <- c("AnnoName1","GeneName","AnnoName2")
    #  k <- plot_df %>%
    #    ggplot(aes(UMAP1, UMAP2, colour=AnnoName1,
    #               text = paste("</br> <b>Sample Name:</b> ", SampleName,
    #                            "</br> <b>",metaColanno1,":</b> ", AnnoName1,
    #                            "</br> <b>",metaColgene," Gene Expression:</b> ", round(GeneName,4),
    #                            "</br> <b>",metaColanno2,":</b> ", AnnoName2,
    #                            sep = "")))
    #}
    #if (input$UMAPannotateSamps1 != " " && input$UMAPannotateSamps2 == " ") {
    #  colnames(plot_df)[c(4,5)] <- c("AnnoName1","GeneName")
    #  k <- plot_df %>%
    #    ggplot(aes(UMAP1, UMAP2, colour=AnnoName1,
    #               text = paste("</br> <b>Sample Name:</b> ", SampleName,
    #                            "</br> <b>",metaColanno1,":</b> ", AnnoName1,
    #                            "</br> <b>",metaColgene," Gene Expression:</b> ", round(GeneName,4),
    #                            sep = "")))
    #}
    #if (input$UMAPannotateSamps1 == " " && input$UMAPannotateSamps2 != " ") {
    #  colnames(plot_df)[c(4,5)] <- c("GeneName","AnnoName2")
    #  k <- plot_df %>%
    #    ggplot(aes(UMAP1, UMAP2,
    #               text = paste("</br> <b>Sample Name:</b> ", SampleName,
    #                            "</br> <b>",metaColgene," Gene Expression:</b> ", round(GeneName,4),
    #                            "</br> <b>",metaColanno2,":</b> ", AnnoName2,
    #                            sep = "")))
    #}
    #if (input$UMAPannotateSamps1 == " " && input$UMAPannotateSamps2 == " ") {
    #  colnames(plot_df)[4] <- "GeneName"
    #  k <- plot_df %>%
    #    ggplot(aes(UMAP1, UMAP2,
    #               text = paste("</br> <b>Sample Name:</b> ", SampleName,
    #                            "</br> <b>",metaColgene," Gene Expression:</b> ", round(GeneName,4),
    #                            sep = "")))
    #}
    
    k <- k + geom_point(shape = 19,
                        size = UMAPdotSize) +
      
      theme_minimal()
    
    if (input$UMAPannotateSamps1 != " ") {
      metaColanno1 <- input$UMAPannotateSamps1
      k <- k + labs(x = umap1,
                    y = umap2,
                    color = metaColanno1)
    }
    if (input$UMAPannotateSamps1 == " ") {
      k <- k + labs(x = umap1,
                    y = umap2)
    }
    
    if (input$UMAPannotateSamps1 != " ") {
      if (!metaColanno1 %in% DiscContList_react()$discreteCovs){
        myPalette <- colorRampPalette(rev(brewer.pal(9, input$UMAPcolors)))
        k <- k + scale_colour_gradientn(colours = rev(myPalette(100)))
      }
    }
    
    k <- k + theme(axis.text = element_text(size = UMAPaxisText),
                   axis.title = element_text(size = UMAPaxisText),
                   plot.title = element_text(size = UMAPtitleText),
                   legend.text=element_text(size=UMAPlegendText))
    
    if (length(sampSelected) > 0) {
      plotdata <- plot_df[sampSelected,]
      plotlabel <- rownames(plot_df[sampSelected,])
      k <- k + geom_point(data = plotdata,
                          aes(x = UMAP1, y = UMAP2),
                          pch = 1,
                          color = "black",
                          size = UMAPdotSize,
                          stroke = .3)
    }
    k
    
  })
  
  umap_plot_PreC_clin_react1 <- reactive({
    
    ## Variables
    sampSelected <- input$UMAPsampSelect         # Sample Names to annotate
    metaColanno1 <- NULL
    metaColanno2 <- NULL
    TrajLines <- input$AddTrajLines
    TrajCluster <- input$TrajClustSelect
    
    plot_df <- UMAP_PreC_CoordTable_react()
    rownames(plot_df) <- plot_df[,1]
    
    umap1 <- input$SelectPreCalc1
    umap2 <- input$SelectPreCalc2
    
    if (input$UMAPannotateSamps1 != " ") {
      metaColanno1 <- input$UMAPannotateSamps1
    }
    if (input$UMAPannotateSamps2 != " ") {
      metaColanno2 <- input$UMAPannotateSamps2
    }
    
    if (!is.null(input$GeneSelection)) {
      metaColgene <- input$GeneSelection
    }
    if (is.null(input$GeneSelection)) {
      metaColgene <- rownames(expr)[1]
    }
    
    
    k <- umap_plot_PreC_clin_react_base1()
    
    if (length(sampSelected) > 0) {
      plotdata <- plot_df[sampSelected,]
      plotlabel <- rownames(plot_df[sampSelected,])
      if (!is.null(metaColanno1) && !is.null(metaColanno2)) {
        colnames(plotdata)[c(4,5,6)] <- c("AnnoName1","GeneName","AnnoName2")
      }
      else if (!is.null(metaColanno1) && is.null(metaColanno2)) {
        colnames(plotdata)[c(4,5)] <- c("AnnoName1","GeneName")
      }
      else if (is.null(metaColanno1) && !is.null(metaColanno2)) {
        colnames(plotdata)[c(4,5)] <- c("GeneName","AnnoName2")
      }
      else if (is.null(metaColanno1) && is.null(metaColanno2)) {
        colnames(plotdata)[4] <- "GeneName"
      }
    }
    
    ## Make it plotly
    k2 <- ggplotly(k,
                   tooltip = "text") %>% 
      
      config(displayModeBar = F)  %>% 
      
      layout(font=list(color="#black"),
             xaxis=list(title=umap1,zeroline=F),
             yaxis=list(title=umap2,zeroline=F)) 
    
    k2 <- k2 %>%
      hide_legend() %>%
      hide_colorbar()
    
    if (length(sampSelected) > 0) {
      k2 <- k2 %>%
        add_annotations(x = plotdata$UMAP1,
                        y = plotdata$UMAP2,
                        text = rownames(plotdata),
                        showarrow = TRUE,
                        arrowhead = 4,
                        arrowsize = .5)
    }
    
    if (TrajLines) {
      if (!is.null(TrajCluster)) {
        
        sce_umap <- SlingShot_react()
        
        if (input$CurveOrLin == "Curves") {
          for (i in seq_along(slingCurves(sce_umap))) {
            curve_i <- slingCurves(sce_umap)[[i]]
            curve_i <- curve_i$s[curve_i$ord, ]
            colnames(curve_i) <- c("UMAP_1", "UMAP_2")
            curve_i <- as.data.frame(curve_i)
            k2 <- k2 %>% add_trace(data = curve_i, x = curve_i$UMAP_1, y = curve_i$UMAP_2, mode = "lines",
                                   hoverinfo='skip',
                                   line = list(color = 'black'))
          }
        } else {
          mst <- slingMST(sce_umap, as.df = TRUE)
          mst_arr <- mst %>% arrange(Order)
          k2 <- k2 %>% add_trace(data = mst %>% arrange(Order), x = mst_arr$UMAP_1, y = mst_arr$UMAP_2, mode = "lines",
                                 hoverinfo='skip',
                                 line = list(color = 'black'))
          k2 <- k2 %>% add_trace(data = mst %>% arrange(Order), x = mst_arr$UMAP_1, y = mst_arr$UMAP_2, mode = "markers",
                                 hoverinfo='skip',
                                 marker = list(size = 10, color = 'black'))
        }
      }
    }
    
    k2
    
  })
  
  
  #slingMST <- function(x, as.df = FALSE) {
  #  
  #  if(!as.df){
  #    return(metadata(x)$mst)
  #  }else{
  #    dfs <- lapply(seq_along(metadata(x)$lineages), function(l){
  #      lin <- metadata(x)$lineages[[l]]
  #      mst <- metadata(x)$mst
  #      centers <- do.call(rbind, V(mst)$coordinates)
  #      rownames(centers) <- V(mst)$name
  #      return(data.frame(centers[lin,], Order = seq_along(lin), 
  #                        Lineage = l, Cluster = lin))
  #    })
  #    return(do.call(rbind, dfs))
  #  }
  #  
  #}
  
  ####----UMAP Expr----####
  
  ## UMAP Plot  - PATH - Clin
  umap_plot_PreC_expr_react_base <- reactive({
    
    ## Variables
    sampSelected <- input$UMAPsampSelect         # Sample Names to annotate
    UMAPtitleText <- input$UMAPtitleTextSize     # Title text size
    UMAPaxisText <- input$UMAPaxisTextSize       # Axis text size
    UMAPlegendText <- input$UMAPlegendTextSize   # Legend text size
    UMAPdotSize <- input$UMAPdotSize             # Dot size
    metaColanno1 <- NULL
    metaColanno2 <- NULL
    
    umap1 <- input$SelectPreCalc1
    umap2 <- input$SelectPreCalc2
    
    plot_df <- UMAP_PreC_CoordTable_react()
    rownames(plot_df) <- plot_df[,1]
    if (input$UMAPannotateSamps1 != " ") {
      metaColanno1 <- input$UMAPannotateSamps1
    }
    if (input$UMAPannotateSamps2 != " ") {
      metaColanno2 <- input$UMAPannotateSamps2
    }
    if (!is.null(input$GeneSelection)) {
      metaColgene <- input$GeneSelection
    }
    if (is.null(input$GeneSelection)) {
      metaColgene <- rownames(expr)[1]
    }
    
    if (!is.null(input$GeneSelection)) {
      metaColgene <- input$GeneSelection
    }
    if (is.null(input$GeneSelection)) {
      metaColgene <- rownames(expr)[1]
    }
    LogChoice <- input$LogGeneSelection
    expr2 <- expr
    if (LogChoice == TRUE) {
      expr2 <- log2(as.matrix(expr) + 1)
    }
    exprRange <- round(range(expr2[metaColgene,], na.rm = T),2)
    
    if (!is.null(input$GeneExprRange)) {
      exprMin <- as.numeric(gsub(" ","",strsplit(as.character(input$GeneExprRange),",")[[1]][1]))
      exprMax <- as.numeric(gsub(" ","",strsplit(as.character(input$GeneExprRange),",")[[1]][2]))
    }
    else if (is.null(input$GeneExprRange)) {
      exprMin <- as.numeric(gsub(" ","",strsplit(as.character(exprRange),",")[[1]][1]))
      exprMax <- as.numeric(gsub(" ","",strsplit(as.character(exprRange),",")[[1]][2]))
    }
    myPalette <- colorRampPalette(rev(brewer.pal(9, input$UMAPcolors)))
    
    ## Expr Annotation Plot
    #if (input$UMAPannotateSamps1 != " " && input$UMAPannotateSamps2 != " ") {
      if (length(colnames(plot_df)[c(4,5,6)]) == 3) {
      colnames(plot_df)[c(4,5,6)] <- c("AnnoName1","GeneName","AnnoName2")
      k <- plot_df %>%
        ggplot(aes(UMAP1, UMAP2, colour=GeneName,
                   text = paste("</br> <b>Sample Name:</b> ", SampleName,
                                "</br> <b>",metaColanno1,":</b> ", AnnoName1,
                                "</br> <b>",metaColgene," Gene Expression:</b> ", round(GeneName,4),
                                "</br> <b>",metaColanno2,":</b> ", AnnoName2,
                                sep = "")))
    }
    #if (input$UMAPannotateSamps1 != " " && input$UMAPannotateSamps2 == " ") {
    #  colnames(plot_df)[c(4,5)] <- c("AnnoName1","GeneName")
    #  k <- plot_df %>%
    #    ggplot(aes(UMAP1, UMAP2, colour=GeneName,
    #               text = paste("</br> <b>Sample Name:</b> ", SampleName,
    #                            "</br> <b>",metaColanno1,":</b> ", AnnoName1,
    #                            "</br> <b>",metaColgene," Gene Expression:</b> ", round(GeneName,4),
    #                            sep = "")))
    #}
    #if (input$UMAPannotateSamps1 == " " && input$UMAPannotateSamps2 != " ") {
    #  colnames(plot_df)[c(4,5)] <- c("GeneName","AnnoName2")
    #  k <- plot_df %>%
    #    ggplot(aes(UMAP1, UMAP2, colour=GeneName,
    #               text = paste("</br> <b>Sample Name:</b> ", SampleName,
    #                            "</br> <b>",metaColgene," Gene Expression:</b> ", round(GeneName,4),
    #                            "</br> <b>",metaColanno2,":</b> ", AnnoName2,
    #                            sep = "")))
    #}
    #if (input$UMAPannotateSamps1 == " " && input$UMAPannotateSamps2 == " ") {
    #  colnames(plot_df)[4] <- "GeneName"
    #  k <- plot_df %>%
    #    ggplot(aes(UMAP1, UMAP2, colour=GeneName,
    #               text = paste("</br> <b>Sample Name:</b> ", SampleName,
    #                            "</br> <b>",metaColgene," Gene Expression:</b> ", round(GeneName,4),
    #                            sep = "")))
    #}
    
    k <- k + geom_point(shape = 19,
                        size = UMAPdotSize) +
      
      theme_minimal()
    
    k <- k + labs(x = umap1,
                  y = umap2,
                  color = metaColgene)
    k <- k + scale_colour_gradientn(colours = rev(myPalette(100)), limits=c(exprMin, exprMax), oob = scales::squish)
    
    
    k <- k + theme(axis.text = element_text(size = UMAPaxisText),
                   axis.title = element_text(size = UMAPaxisText),
                   plot.title = element_text(size = UMAPtitleText),
                   legend.text=element_text(size=UMAPlegendText))
    
    if (length(sampSelected) > 0) {
      plotdata <- plot_df[sampSelected,]
      plotlabel <- rownames(plot_df[sampSelected,])
      k <- k + geom_point(data = plotdata,
                          aes(x = UMAP1, y = UMAP2),
                          pch = 1,
                          color = "black",
                          size = UMAPdotSize,
                          stroke = .3)
    }
    k
    
  })
  
  umap_plot_PreC_expr_react <- reactive({
    
    ## Variables
    sampSelected <- input$UMAPsampSelect         # Sample Names to annotate
    metaColanno1 <- NULL
    metaColanno2 <- NULL
    TrajLines <- input$AddTrajLines
    TrajCluster <- input$TrajClustSelect
    
    plot_df <- UMAP_PreC_CoordTable_react()
    rownames(plot_df) <- plot_df[,1]
    
    umap1 <- input$SelectPreCalc1
    umap2 <- input$SelectPreCalc2
    
    if (input$UMAPannotateSamps1 != " ") {
      metaColanno1 <- input$UMAPannotateSamps1
    }
    if (input$UMAPannotateSamps2 != " ") {
      metaColanno2 <- input$UMAPannotateSamps2
    }
    if (!is.null(input$GeneSelection)) {
      metaColgene <- input$GeneSelection
    }
    if (is.null(input$GeneSelection)) {
      metaColgene <- rownames(expr)[1]
    }
    
    
    k <- umap_plot_PreC_expr_react_base()
    
    if (length(sampSelected) > 0) {
      plotdata <- plot_df[sampSelected,]
      plotlabel <- rownames(plot_df[sampSelected,])
      if (!is.null(metaColanno1) && !is.null(metaColanno2)) {
        colnames(plotdata)[c(4,5,6)] <- c("AnnoName1","GeneName","AnnoName2")
      }
      else if (!is.null(metaColanno1) && is.null(metaColanno2)) {
        colnames(plotdata)[c(4,5)] <- c("AnnoName1","GeneName")
      }
      else if (is.null(metaColanno1) && !is.null(metaColanno2)) {
        colnames(plotdata)[c(4,5)] <- c("GeneName","AnnoName2")
      }
      else if (is.null(metaColanno1) && is.null(metaColanno2)) {
        colnames(plotdata)[4] <- "GeneName"
      }
    }
    
    
    ## Make it plotly
    k2 <- ggplotly(k,
                   tooltip = "text") %>% 
      
      config(displayModeBar = F)  %>% 
      
      layout(font=list(color="#black"),
             xaxis=list(title=umap1,zeroline=F),
             yaxis=list(title=umap2,zeroline=F)) 
    
    k2 <- k2 %>%
      hide_legend() %>%
      hide_colorbar()
    
    if (length(sampSelected) > 0) {
      k2 <- k2 %>%
        add_annotations(x = plotdata$UMAP1,
                        y = plotdata$UMAP2,
                        text = rownames(plotdata),
                        showarrow = TRUE,
                        arrowhead = 4,
                        arrowsize = .5)
    }
    
    if (TrajLines) {
      if (!is.null(TrajCluster)) {
        
        sce_umap <- SlingShot_react()
        
        if (input$CurveOrLin == "Curves") {
          for (i in seq_along(slingCurves(sce_umap))) {
            curve_i <- slingCurves(sce_umap)[[i]]
            curve_i <- curve_i$s[curve_i$ord, ]
            colnames(curve_i) <- c("UMAP_1", "UMAP_2")
            curve_i <- as.data.frame(curve_i)
            k2 <- k2 %>% add_trace(data = curve_i, x = curve_i$UMAP_1, y = curve_i$UMAP_2, mode = "lines",
                                   hoverinfo='skip',
                                   line = list(color = 'black'))
          }
        } else {
          mst <- slingMST(sce_umap, as.df = TRUE)
          mst_arr <- mst %>% arrange(Order)
          k2 <- k2 %>% add_trace(data = mst %>% arrange(Order), x = mst_arr$UMAP_1, y = mst_arr$UMAP_2, mode = "lines",
                                 hoverinfo='skip',
                                 line = list(color = 'black'))
          k2 <- k2 %>% add_trace(data = mst %>% arrange(Order), x = mst_arr$UMAP_1, y = mst_arr$UMAP_2, mode = "markers",
                                 hoverinfo='skip',
                                 marker = list(size = 10, color = 'black'))
        }
      }
    }
    
    k2
    
  })
  
  ####----UMAP Clin2----####
  
  ## UMAP Plot  - PATH - Clin
  umap_plot_PreC_clin_react_base2 <- reactive({
    
    ## Variables
    sampSelected <- input$UMAPsampSelect         # Sample Names to annotate
    UMAPtitleText <- input$UMAPtitleTextSize     # Title text size
    UMAPaxisText <- input$UMAPaxisTextSize       # Axis text size
    UMAPlegendText <- input$UMAPlegendTextSize   # Legend text size
    UMAPdotSize <- input$UMAPdotSize             # Dot size
    metaColanno1 <- NULL
    metaColanno2 <- NULL
    
    umap1 <- input$SelectPreCalc1
    umap2 <- input$SelectPreCalc2
    
    plot_df <- UMAP_PreC_CoordTable_react()
    rownames(plot_df) <- plot_df[,1]
    if (input$UMAPannotateSamps1 != " ") {
      metaColanno1 <- input$UMAPannotateSamps1
    }
    if (input$UMAPannotateSamps2 != " ") {
      metaColanno2 <- input$UMAPannotateSamps2
    }
    
    if (!is.null(input$GeneSelection)) {
      metaColgene <- input$GeneSelection
    }
    if (is.null(input$GeneSelection)) {
      metaColgene <- rownames(expr)[1]
    }
    
    ## Meta Annotation Plot
    #if (input$UMAPannotateSamps1 != " " && input$UMAPannotateSamps2 != " ") {
      if (length(colnames(plot_df)[c(4,5,6)]) == 3) {
      colnames(plot_df)[c(4,5,6)] <- c("AnnoName1","GeneName","AnnoName2")
      k <- plot_df %>%
        ggplot(aes(UMAP1, UMAP2, colour=AnnoName2,
                   text = paste("</br> <b>Sample Name:</b> ", SampleName,
                                "</br> <b>",metaColanno1,":</b> ", AnnoName1,
                                "</br> <b>",metaColgene," Gene Expression:</b> ", round(GeneName,4),
                                "</br> <b>",metaColanno2,":</b> ", AnnoName2,
                                sep = "")))
    }
    #if (input$UMAPannotateSamps1 != " " && input$UMAPannotateSamps2 == " ") {
    #  colnames(plot_df)[c(4,5)] <- c("AnnoName1","GeneName")
    #  k <- plot_df %>%
    #    ggplot(aes(UMAP1, UMAP2,
    #               text = paste("</br> <b>Sample Name:</b> ", SampleName,
    #                            "</br> <b>",metaColanno1,":</b> ", AnnoName1,
    #                            "</br> <b>",metaColgene," Gene Expression:</b> ", round(GeneName,4),
    #                            sep = "")))
    #}
    #if (input$UMAPannotateSamps1 == " " && input$UMAPannotateSamps2 != " ") {
    #  colnames(plot_df)[c(4,5)] <- c("GeneName","AnnoName2")
    #  k <- plot_df %>%
    #    ggplot(aes(UMAP1, UMAP2, colour=AnnoName2,
    #               text = paste("</br> <b>Sample Name:</b> ", SampleName,
    #                            "</br> <b>",metaColgene," Gene Expression:</b> ", round(GeneName,4),
    #                            "</br> <b>",metaColanno2,":</b> ", AnnoName2,
    #                            sep = "")))
    #}
    #if (input$UMAPannotateSamps1 == " " && input$UMAPannotateSamps2 == " ") {
    #  colnames(plot_df)[4] <- "GeneName"
    #  k <- plot_df %>%
    #    ggplot(aes(UMAP1, UMAP2,
    #               text = paste("</br> <b>Sample Name:</b> ", SampleName,
    #                            "</br> <b>",metaColgene," Gene Expression:</b> ", round(GeneName,4),
    #                            sep = "")))
    #}
    
    k <- k + geom_point(shape = 19,
                        size = UMAPdotSize) +
      
      theme_minimal()
    
    if (input$UMAPannotateSamps2 != " ") {
      metaColanno2 <- input$UMAPannotateSamps2
      k <- k + labs(x = umap1,
                    y = umap2,
                    color = metaColanno2)
    }
    if (input$UMAPannotateSamps2 == " ") {
      k <- k + labs(x = umap1,
                    y = umap2)
    }
    
    if (input$UMAPannotateSamps2 != " ") {
      if (!metaColanno2 %in% DiscContList_react()$discreteCovs) {
        myPalette <- colorRampPalette(rev(brewer.pal(9, input$UMAPcolors)))
        k <- k + scale_colour_gradientn(colours = rev(myPalette(100)))
      }
    }
    
    k <- k + theme(axis.text = element_text(size = UMAPaxisText),
                   axis.title = element_text(size = UMAPaxisText),
                   plot.title = element_text(size = UMAPtitleText),
                   legend.text=element_text(size=UMAPlegendText))
    
    if (length(sampSelected) > 0) {
      plotdata <- plot_df[sampSelected,]
      plotlabel <- rownames(plot_df[sampSelected,])
      k <- k + geom_point(data = plotdata,
                          aes(x = UMAP1, y = UMAP2),
                          pch = 1,
                          color = "black",
                          size = UMAPdotSize,
                          stroke = .3)
    }
    k
    
  })
  
  umap_plot_PreC_clin_react2 <- reactive({
    
    ## Variables
    sampSelected <- input$UMAPsampSelect         # Sample Names to annotate
    metaColanno1 <- NULL
    metaColanno2 <- NULL
    TrajLines <- input$AddTrajLines
    TrajCluster <- input$TrajClustSelect
    
    plot_df <- UMAP_PreC_CoordTable_react()
    rownames(plot_df) <- plot_df[,1]
    
    umap1 <- input$SelectPreCalc1
    umap2 <- input$SelectPreCalc2
    
    if (input$UMAPannotateSamps1 != " ") {
      metaColanno1 <- input$UMAPannotateSamps1
    }
    if (input$UMAPannotateSamps2 != " ") {
      metaColanno2 <- input$UMAPannotateSamps2
    }
    
    if (!is.null(input$GeneSelection)) {
      metaColgene <- input$GeneSelection
    }
    if (is.null(input$GeneSelection)) {
      metaColgene <- rownames(expr)[1]
    }
    
    
    k <- umap_plot_PreC_clin_react_base2()
    
    if (length(sampSelected) > 0) {
      plotdata <- plot_df[sampSelected,]
      plotlabel <- rownames(plot_df[sampSelected,])
      if (!is.null(metaColanno1) && !is.null(metaColanno2)) {
        colnames(plotdata)[c(4,5,6)] <- c("AnnoName1","GeneName","AnnoName2")
      }
      else if (!is.null(metaColanno1) && is.null(metaColanno2)) {
        colnames(plotdata)[c(4,5)] <- c("AnnoName1","GeneName")
      }
      else if (is.null(metaColanno1) && !is.null(metaColanno2)) {
        colnames(plotdata)[c(4,5)] <- c("GeneName","AnnoName2")
      }
      else if (is.null(metaColanno1) && is.null(metaColanno2)) {
        colnames(plotdata)[4] <- "GeneName"
      }
    }
    
    ## Make it plotly
    k2 <- ggplotly(k,
                   tooltip = "text") %>% 
      
      config(displayModeBar = F)  %>% 
      
      layout(font=list(color="#black"),
             xaxis=list(title=umap1,zeroline=F),
             yaxis=list(title=umap2,zeroline=F)) 
    
    k2 <- k2 %>%
      hide_legend() %>%
      hide_colorbar()
    
    if (length(sampSelected) > 0) {
      k2 <- k2 %>%
        add_annotations(x = plotdata$UMAP1,
                        y = plotdata$UMAP2,
                        text = rownames(plotdata),
                        showarrow = TRUE,
                        arrowhead = 4,
                        arrowsize = .5)
    }
    
    if (TrajLines) {
      if (!is.null(TrajCluster)) {
        
        sce_umap <- SlingShot_react()
        
        if (input$CurveOrLin == "Curves") {
          for (i in seq_along(slingCurves(sce_umap))) {
            curve_i <- slingCurves(sce_umap)[[i]]
            curve_i <- curve_i$s[curve_i$ord, ]
            colnames(curve_i) <- c("UMAP_1", "UMAP_2")
            curve_i <- as.data.frame(curve_i)
            k2 <- k2 %>% add_trace(data = curve_i, x = curve_i$UMAP_1, y = curve_i$UMAP_2, mode = "lines",
                                   hoverinfo='skip',
                                   line = list(color = 'black'))
          }
        } else {
          mst <- slingMST(sce_umap, as.df = TRUE)
          mst_arr <- mst %>% arrange(Order)
          k2 <- k2 %>% add_trace(data = mst %>% arrange(Order), x = mst_arr$UMAP_1, y = mst_arr$UMAP_2, mode = "lines",
                                 hoverinfo='skip',
                                 line = list(color = 'black'))
          k2 <- k2 %>% add_trace(data = mst %>% arrange(Order), x = mst_arr$UMAP_1, y = mst_arr$UMAP_2, mode = "markers",
                                 hoverinfo='skip',
                                 marker = list(size = 10, color = 'black'))
        }
      }
    }
    
    k2
    
    
  })
  
  output$UMAPprecLegend <- renderPlot({
    
    if (!is.null(input$UMAPannotateSamps1) && !is.null(input$UMAPannotateSamps2)) {
      if (input$UMAPannotateSamps1 != " " && input$UMAPannotateSamps2 != " ") {
        k <- umap_plot_PreC_clin_react_base1()
        m <- umap_plot_PreC_expr_react_base()
        n <- umap_plot_PreC_clin_react_base2()
        legendk <- g_legend(k)
        legendm <- g_legend(m)
        legendn <- g_legend(n)
        legend_grid <- grid.arrange(legendk,legendm,legendn,nrow = 1)
      }
      else if (input$UMAPannotateSamps1 != " " && input$UMAPannotateSamps2 == " ") {
        k <- umap_plot_PreC_clin_react_base1()
        m <- umap_plot_PreC_expr_react_base()
        legendk <- g_legend(k)
        legendm <- g_legend(m)
        legend_grid <- grid.arrange(legendk,legendm,nrow = 1)
      }
      else if (input$UMAPannotateSamps1 == " " && input$UMAPannotateSamps2 != " ") {
        m <- umap_plot_PreC_expr_react_base()
        n <- umap_plot_PreC_clin_react_base2()
        legendm <- g_legend(m)
        legendn <- g_legend(n)
        legend_grid <- grid.arrange(legendm,legendn,nrow = 1)
      }
      if (input$UMAPannotateSamps1 == " " && input$UMAPannotateSamps2 == " ") {
        if (!is.null(input$GeneSelection)) {
          GeneSelec <- input$GeneSelection
        }
        if (is.null(input$GeneSelection)) {
          GeneSelec <- rownames(expr)[1]
        }
        if (length(GeneSelec) > 0) {
          m <- umap_plot_PreC_expr_react_base()
          legendm <- g_legend(m)
          legend_grid <- grid.arrange(legendm,nrow = 1)
        }
      }
      
      legend_grid
    }
    
  })
  
  UMAPplot_PreC_ALL_react <- reactive({
    
    k2 <- umap_plot_PreC_clin_react1()
    m2 <- umap_plot_PreC_expr_react()
    n2 <- umap_plot_PreC_clin_react2()
    
    if (input$UMAPannotateSamps1 != " ") {
      UMAPtitlek <- paste("Annotated by",input$UMAPannotateSamps1)
    }
    else if (input$UMAPannotateSamps1 == " ") {
      UMAPtitlek <- ""
    }
    
    if (!is.null(input$GeneSelection)) {
      GeneSelec <- input$GeneSelection
    }
    if (is.null(input$GeneSelection)) {
      GeneSelec <- rownames(expr)[1]
    }
    if (input$LogGeneSelection == T) {
      GeneSelecExpr <- paste(GeneSelec,"Expression (Log2)")
    }
    if (input$LogGeneSelection == F) {
      GeneSelecExpr <- paste(GeneSelec,"Expression")
    }
    UMAPtitlem <- paste("Annotated by",GeneSelecExpr)
    
    if (input$UMAPannotateSamps2 != " ") {
      UMAPtitlen <- paste("Annotated by",input$UMAPannotateSamps2)
    }
    else if (input$UMAPannotateSamps2 == " ") {
      UMAPtitlen <- ""
    }
    
    
    if (input$UMAPorientation == "Side-by-Side") {
      subplot_all_sbs <- subplot(k2, m2, n2, nrows = 1, shareY = TRUE, shareX = TRUE)
      plot_titles = list( 
        list( 
          x = 0,
          y = 1,
          text = UMAPtitlek,  
          xref = "x",  
          yref = "paper",  
          xanchor = "center",  
          yanchor = "bottom",  
          showarrow = FALSE 
        ),  
        list( 
          x = 0,
          y = 1,
          text = UMAPtitlem,  
          xref = "x2",
          yref = "paper",
          xanchor = "center",  
          yanchor = "bottom",    
          showarrow = FALSE 
        ),  
        list( 
          x = 0,
          y = 1,
          text = UMAPtitlen,  
          xref = "x3",
          yref = "paper",
          xanchor = "center",  
          yanchor = "bottom",  
          showarrow = FALSE 
        ))
    }
    else if (input$UMAPorientation == "Stacked") {
      subplot_all_sbs <- subplot(k2, m2, n2, nrows = 3, shareY = TRUE, shareX = TRUE, margin = 0.04)
      plot_titles = list( 
        list( 
          x = 0,
          y = 1,
          text = UMAPtitlek,  
          xref = "x",  
          yref = "paper",  
          xanchor = "center",  
          yanchor = "bottom",  
          showarrow = FALSE 
        ),  
        list( 
          x = 0,
          y = 0.65,
          text = UMAPtitlem,  
          xref = "x",
          yref = "paper",
          xanchor = "center",  
          yanchor = "bottom",    
          showarrow = FALSE 
        ),  
        list( 
          x = 0,
          y = 0.3,
          text = UMAPtitlen,  
          xref = "x",
          yref = "paper",
          xanchor = "center",  
          yanchor = "bottom",  
          showarrow = FALSE 
        ))
    }
    
    subplot_all_sbs <- subplot_all_sbs %>% layout(annotations = plot_titles)
    
    subplot_all_sbs
    
  })
  
  output$UMAPplot_PreC_ALL <- renderPlotly({
    
    p_all <- UMAPplot_PreC_ALL_react()
    p_all
    
  })
  
  ####----Bar Plot----####
  
  barplot_react <- reactive({
    
    if (!is.null(input$BPsampSubset)) {
      if (!is.null(input$BPgroupCriteria)) {
        
        title_font <- input$barplot1TitleSize            # Title font size
        Xaxis_font <- input$barplot1AxisSize              # Axis font size
        Yaxis_font <- input$barplot1AxisSize              # Axis font size
        anno_size <- input$Vplot1aAnnoSize
        hjust_orient <- 1                                # Initial hjust
        axis_orient <- as.numeric(input$barxAxisOrient)  # X-axis label orientation
        if (axis_orient == 0) {                          # Adjust hjust if orientation is 0
          hjust_orient <- 0.5
        }
        logchoice <- input$log2barplot                   # Log expression data option
        colorin <- input$barplotColoCodes                # Bar plot Color codes
        Vplotylim <- input$barPlotYlim                      # Y-limit
        Vplotybreaks <- input$barplotYbreaks                # Y-axis breaks
        dotsizein <- input$VplotDotSize
        StatMethod <- input$BPstatComp
        bar_type <- input$ErrorBP
        dotChoice <- input$Vplotsampledots
        VilOrBP <- input$ViolinOrBoxP
        DiscContList <- DiscContList_react()
        
        FeatMat <- expr
        metaSub <- BP_meta_subset()
        FeatMat <- FeatMat[,metaSub$SampleName]
        sampSubset <- input$BPsampSubset
        sampCriteria <- input$BPsampCriteria
        groupCriteria <- input$BPgroupCriteria
        featSelected <- input$BPgeneSelection
        sampSelected <- input$UMAPsampSelect
        metaSub2 <- metaSub[,c("SampleName",groupCriteria)]
        metaSub2[,groupCriteria] <- as.factor(metaSub2[,groupCriteria])
        
        if (length(featSelected) > 0){
          
          feature <- featSelected
          if (feature %in% rownames(FeatMat)) {
            feat_gene <- as.data.frame(t(FeatMat[feature,]))
            feat_gene$SampleName <- rownames(feat_gene)
          }
          if (!feature %in% rownames(FeatMat)) {
            feat_gene <- as.data.frame(metaSub[,c("SampleName",feature)])
          }
          feat_gene2 <- merge(feat_gene,metaSub2)
          feature_lab <- feature
          feattitle <- feature_lab
          if (logchoice == T) {
            feature_lab <- paste(feature_lab,"(Log2)")
            feattitle <- feature_lab
          }
          
          if (sampSubset == "Select All") {
            plottitle <- paste(feature_lab, "Across ",groupCriteria)
          }
          if (sampSubset != "Select All") {
            plottitle <- paste(feature_lab,"in",sampSubset,"-",sampCriteria,"Across",groupCriteria)
          }
          
          colnames(feat_gene2) <- c("SampleName","FeatureName","Type")
          
          if (VilOrBP == "Stacked Barplot") {
            #if (featSelected %in% DiscContList$continuousCovs | feature %in% rownames(expr)) {
            #  feat_gene2[,"FeatureName"] <- highlow(as.numeric(feat_gene2[,"FeatureName"]))
            #}
            feat_gene2[,"FeatureName"] <- as.factor(feat_gene2[,"FeatureName"])
            barp <- ggplot(feat_gene2, aes(fill = FeatureName, x = Type))
            FillChoice <- input$BPfeatFill
            if (FillChoice == FALSE) {
              barp <- barp + geom_bar() +
                theme_minimal()
            }
            if (FillChoice == TRUE) {
              barp <- barp + geom_bar(position = "fill") +
                theme_minimal()
            }
            if (colorin != "") {
              colorin <- strsplit(colorin," ")[[1]]
              if (length(colorin) == 1) {
                colorin <- rep(colorin,length(unique(metaSub[,groupCriteria])))
              }
              barp <- barp + scale_fill_manual(values=colorin)
            }
            barp <- barp + labs(x = groupCriteria,
                                y = feattitle,
                                title = plottitle)
            barp <- barp + theme(axis.text.x = element_text(angle = axis_orient, hjust = hjust_orient,size = Xaxis_font),
                                 axis.title.x = element_text(size = Xaxis_font),
                                 axis.text.y = element_text(size = Yaxis_font),
                                 axis.title.y = element_text(size = Yaxis_font),
                                 plot.title = element_text(size = title_font, margin=margin(0,0,30,0)))
            barp
            
          } 
          
          else if (VilOrBP == "Barplot") {
            if (logchoice == T) {
              #feat_gene2[which(feat_gene2[,2] < 0),2] <- 0
              feat_gene2[,2] <- log2(feat_gene2[,2] + 1)
            }
            se <- function(x) sd(x, na.rm = T)/sqrt(length(x))
            feat_gene_stats <- feat_gene2 %>%
              group_by(Type) %>%
              summarise_at("FeatureName",list(mean = mean, sd = sd, se = se))
            y_min <- 0
            y_max <- round(max((feat_gene2[,2]) + 1, na.rm = T))
            if (Vplotylim != "") {
              y_min <- as.numeric(gsub(" ","",strsplit(Vplotylim,",")[[1]][1]))
              y_max <- as.numeric(gsub(" ","",strsplit(Vplotylim,",")[[1]][2]))
            }
            
            if (input$barplotXaxOrder == "Descending"){
              barp <- ggplot(feat_gene_stats, aes(reorder(Type,desc(mean)), mean, fill=Type))
            }
            if (input$barplotXaxOrder == "Ascending"){
              barp <- ggplot(feat_gene_stats, aes(reorder(Type, mean), mean, fill=Type))
            }
            if (input$barplotXaxOrder == "Not Specificed"){
              barp <- ggplot(feat_gene_stats, aes(Type, mean, fill=Type))
            }
            if (bar_type != "None") {
              if (bar_type == "Standard Deviation") {
                barp <- barp + geom_errorbar(aes(ymin=0,ymax = mean+sd), size = 0.75, width = 0.5)
              }
              else if (bar_type == "Standard Error") {
                barp <- barp + geom_errorbar(aes(ymin=0,ymax = mean+se), size = 0.75, width = 0.5)
              }
            }
            barp <- barp + geom_bar(stat = "identity",
                                    width=0.75,
                                    size = 0.75,
                                    color="black",
                                    show.legend = FALSE) +
              theme_minimal()
            
            barp <- barp + labs(x = groupCriteria,
                                y = feattitle,
                                title = plottitle)
            barp <- barp + theme(axis.text.x = element_text(angle = axis_orient, hjust = hjust_orient),
                                 axis.text = element_text(size = Xaxis_font),
                                 axis.title = element_text(size = Xaxis_font),
                                 plot.title = element_text(size = title_font))
            if (Vplotylim == "") {
              if (is.na(Vplotybreaks)) {
                barp <- barp + scale_y_continuous(expand = c(0, 0))
              }
              else if (!is.na(Vplotybreaks)) {
                barp <- barp + scale_y_continuous(limits=c(0,y_max),expand = c(0, 0),breaks=seq(0,y_max,Vplotybreaks),oob = rescale_none)
              }
            }
            else if (Vplotylim != "") {
              if (is.na(Vplotybreaks)) {
                barp <- barp + scale_y_continuous(limits=c(y_min,y_max),expand = c(0, 0),oob = rescale_none)
              }
              else if (!is.na(Vplotybreaks)) {
                barp <- barp + scale_y_continuous(limits=c(y_min,y_max),expand = c(0, 0),breaks=seq(y_min,y_max,Vplotybreaks),oob = rescale_none)
              }
            }
            
            if (dotChoice == T) {
              barp <- barp + geom_dotplot(data = feat_gene2, aes(x=Type,y=GeneName),
                                          binaxis='y', stackdir='center',
                                          stackratio=1, dotsize=dotsizein, fill = "black")
            }
            
            if (colorin != "") {
              colorin <- strsplit(colorin," ")[[1]]
              if (length(colorin) == 1) {
                colorin <- rep(colorin,length(unique(metaSub2[,groupCriteria])))
              }
              barp <- barp + scale_fill_manual(values=colorin)
            }
            
            if (StatMethod != "None") {
              barp <- barp + stat_compare_means(data = feat_gene2, aes(x=Type,y=FeatureName), method = StatMethod)
            }
            
            barp <- barp + coord_cartesian(clip = "off")
            
          }
          
          else if (VilOrBP %in% c("Violin Plot","Box Plot")) {
            if (logchoice == T) {
              #feat_gene2[which(feat_gene2[,2] < 0),2] <- 0
              feat_gene2[,2] <- log2(feat_gene2[,2] + 1)
            }
            
            y_min <- 0
            y_max <- max(feat_gene2$FeatureName, na.rm = T) * 1.03
            if (Vplotylim != "") {
              y_min <- as.numeric(gsub(" ","",strsplit(Vplotylim,",")[[1]][1]))
              y_max <- as.numeric(gsub(" ","",strsplit(Vplotylim,",")[[1]][2]))
            }
            
            if (input$barplotXaxOrder == "Descending"){
              feat_gene2$Type <- reorder(feat_gene2$Type,-feat_gene2$FeatureName)
            }
            if (input$barplotXaxOrder == "Ascending"){
              feat_gene2$Type <- reorder(feat_gene2$Type,feat_gene2$FeatureName)
            }
            
            barp <- ggplot(data = feat_gene2, aes(x=Type,y=FeatureName, fill=Type))
            #if (VilOrBP == "Stacked Barplot") {
            #  barp <- ggplot(feat_gene2, aes(fill = FeatureName, x = Type))
            #}
            #if (VilOrBP != "Stacked Barplot") {
            #  barp <- ggplot(data = feat_gene2, aes(x=Type,y=FeatureName, fill=Type))
            #}
            if (StatMethod != "None") {
              barp <- barp + stat_compare_means(data = feat_gene2,
                                                aes(x=Type,y=FeatureName),
                                                method = StatMethod,
                                                label.x = 1,
                                                label.y = max(feat_gene2$FeatureName, na.rm = T) * 1.02)
            }
            if (VilOrBP == "Box Plot") {
              barp <- barp + geom_boxplot() +
                theme_minimal()
            }
            if (VilOrBP == "Violin Plot") {
              barp <- barp + geom_violin() +
                theme_minimal()
              barp <- barp + stat_summary(fun=median, geom="point", shape=23, size=dotsizein, color="black", bg = "black")
            }
            #if (VilOrBP == "Stacked Barplot") {
            #  FillChoice <- input$BPfeatFill
            #  if (FillChoice == FALSE) {
            #    barp <- barp + geom_bar() +
            #      theme_minimal()
            #  }
            #  else if (FillChoice == TRUE) {
            #    barp <- barp + geom_bar() +
            #      theme_minimal(position = "fill")
            #  }
            #}
            
            if (Vplotylim == "") {
              if (is.na(Vplotybreaks)) {
                #barp <- barp + scale_y_continuous(expansion(mult = c(0, 0.1))) ## expansion function causing issue with y-axis title
                barp <- barp + scale_y_continuous()
              }
              else if (!is.na(Vplotybreaks)) {
                #barp <- barp + scale_y_continuous(limits=c(0,y_max),breaks=seq(0,y_max,Vplotybreaks))
                barp <- barp + scale_y_continuous(limits=c(0,y_max),expansion(mult = c(0, 0.1)),breaks=seq(0,y_max,Vplotybreaks))
              }
            }
            else if (Vplotylim != "") {
              if (is.na(Vplotybreaks)) {
                #barp <- barp + scale_y_continuous(limits=c(y_min,y_max),expansion(mult = c(0, 0.1)))
                barp <- barp + scale_y_continuous(limits=c(y_min,y_max))
              }
              else if (!is.na(Vplotybreaks)) {
                #barp <- barp + scale_y_continuous(limits=c(y_min,y_max),expansion(mult = c(0, 0.1)),breaks=seq(y_min,y_max,Vplotybreaks))
                barp <- barp + scale_y_continuous(limits=c(y_min,y_max),breaks=seq(y_min,y_max,Vplotybreaks))
              }
            }
            
            if (dotChoice == T) {
              feat_gene3 <- feat_gene2
              feat_gene3$xj <- jitter(as.numeric(factor(feat_gene3$Type)))
              barp <- barp + geom_point(data = feat_gene3, aes(x=xj), col="gray44", size=dotsizein)
              if (length(sampSelected) > 0) {
                feat_gene3$SampLabel <- ifelse(feat_gene3[,1] %in% sampSelected, feat_gene3[,1], NA)
                barp <- barp + geom_text_repel(data = feat_gene3,
                                               aes(x = xj, y = FeatureName, label = SampLabel),
                                               size = anno_size,
                                               max.overlaps = Inf,
                                               color = "black",
                                               fontface = "bold",
                                               min.segment.length = unit(0, 'lines'),
                                               box.padding = unit(0.35, "lines"),
                                               point.padding = unit(0.3, "lines"),
                                               arrow = arrow(length = unit(0.015, "npc"))
                )
                feat_gene4 <- feat_gene3[!is.na(feat_gene3$SampLabel),]
                barp <- barp + geom_point(data = feat_gene4, aes(x=xj), col="darkred", size=dotsizein)
              }
            }
            
            if (colorin != "") {
              colorin <- strsplit(colorin," ")[[1]]
              if (length(colorin) == 1) {
                colorin <- rep(colorin,length(unique(metaSub[,groupCriteria])))
              }
              barp <- barp + scale_fill_manual(values=colorin)
            }
            barp <- barp + labs(x = groupCriteria,
                                y = feattitle,
                                title = plottitle)
            barp <- barp + theme(axis.text.x = element_text(angle = axis_orient, hjust = hjust_orient,size = Xaxis_font),
                                 axis.title.x = element_text(size = Xaxis_font),
                                 axis.text.y = element_text(size = Yaxis_font),
                                 axis.title.y = element_text(size = Yaxis_font),
                                 plot.title = element_text(size = title_font, margin=margin(0,0,30,0)),
                                 legend.position = 'none')
            barp
            
          }
          
        }
        
      }
    }
    
  })
  
  ##render boxplot
  #barplot_react <- reactive({
  #  
  #  if (!is.null(input$BPsampSubset) && !is.null(input$BPsampSubset)) {
  #    #if (input$BPsampSubset != " " && input$BPgroupCriteria != " ") {
  #      
  #      title_font <- input$barplot1TitleSize            # Title font size
  #      axis_font <- input$barplot1AxisSize              # Axis font size
  #      hjust_orient <- 1                                # Initial hjust
  #      axis_orient <- as.numeric(input$barxAxisOrient)  # X-axis label orientation
  #      if (axis_orient == 0) {                          # Adjust hjust if orientation is 0
  #        hjust_orient <- 0.5
  #      }
  #      bar_type <- input$errorbarplot                   # Error bar type
  #      logchoice <- input$log2barplot                   # Log expression data option
  #      colorin <- input$barplotColoCodes                # Bar plot Color codes
  #      bpylim <- input$barPlotYlim                      # Y-limit
  #      bpybreaks <- input$barplotYbreaks                # Y-axis breaks
  #      dotsizein <- input$barplotDotSize
  #      StatMethod <- input$BPstatComp
  #      dotChoice <- input$barplotsampledots
  #      #if (is.null(dotChoice)) {
  #      #  dotChoice <- FALSE
  #      #}
  #      metaSub <- BP_meta_subset()
  #      sampSubset <- input$BPsampSubset
  #      sampCriteria <- input$BPsampCriteria
  #      sampLable <- paste(sampSubset,sampCriteria,"Average Gene Expression")
  #      groupCriteria <- input$BPgroupCriteria
  #      geneSelected <- input$BPgeneSelection
  #      metaSub2 <- metaSub[,c("SampleName",groupCriteria)]
  #      metaSub2[,groupCriteria] <- as.factor(metaSub2[,groupCriteria])
  #      
  #      exprSub <- expr[,metaSub[,1]]
  #      
  #      if (input$BPsampSubset == "Select All") {
  #        sampLable <- "Average Gene Expression"
  #      }
  #      
  #      if (length(geneSelected) > 0){
  #        
  #        gene <- geneSelected
  #        
  #        expr_gene <- as.data.frame(t(exprSub[gene,]))
  #        expr_gene$SampleName <- rownames(expr_gene)
  #        expr_gene2 <- merge(expr_gene,metaSub2,by = "SampleName", all = T)
  #        plottitle <- paste(sampLable,"of",gene,"\nAcross",groupCriteria)
  #        genetitle <- paste(gene,"Average Expression")
  #        
  #        if (logchoice == T) {
  #          expr_gene2[which(expr_gene2[,2] < 0),2] <- 0
  #          expr_gene2[,2] <- log2(expr_gene2[,2] + 1)
  #          plottitle <- paste(sampLable,"(Log2) of",gene,"\nAcross",groupCriteria)
  #          genetitle <- paste(gene,"Average Expression (Log2)")
  #        }
  #        
  #        colnames(expr_gene2) <- c("SampleName","GeneName","Type")
  #        
  #        se <- function(x) sd(x)/sqrt(length(x))
  #        expr_gene_stats <- expr_gene2 %>%
  #          group_by(Type) %>%
  #          summarise_at("GeneName",list(mean = mean, sd = sd, se = se))
  #        
  #        y_min <- 0
  #        y_max <- round(max((expr_gene2[,2]) + 1))
  #        if (bpylim != "") {
  #          y_min <- as.numeric(gsub(" ","",strsplit(bpylim,",")[[1]][1]))
  #          y_max <- as.numeric(gsub(" ","",strsplit(bpylim,",")[[1]][2]))
  #        }
  #        
  #        if (input$barplotXaxOrder == "Descending"){
  #          barp <- ggplot(expr_gene_stats, aes(reorder(Type,desc(mean)), mean, fill=Type))
  #        }
  #        if (input$barplotXaxOrder == "Ascending"){
  #          barp <- ggplot(expr_gene_stats, aes(reorder(Type, mean), mean, fill=Type))
  #        }
  #        if (input$barplotXaxOrder == "Not Specificed"){
  #          barp <- ggplot(expr_gene_stats, aes(Type, mean, fill=Type))
  #        }
  #        if (bar_type != "None") {
  #          if (bar_type == "Standard Deviation") {
  #            barp <- barp + geom_errorbar(aes(ymin=0,ymax = mean+sd), size = 0.75, width = 0.5)
  #          }
  #          else if (bar_type == "Standard Error") {
  #            barp <- barp + geom_errorbar(aes(ymin=0,ymax = mean+se), size = 0.75, width = 0.5)
  #          }
  #        }
  #        barp <- barp + geom_bar(stat = "identity",
  #                                width=0.75,
  #                                size = 0.75,
  #                                color="black",
  #                                show.legend = FALSE) +
  #          theme_minimal()
  #        
  #        barp <- barp + labs(x = groupCriteria,
  #                            y = genetitle,
  #                            title = plottitle)
  #        barp <- barp + theme(axis.text.x = element_text(angle = axis_orient, hjust = hjust_orient),
  #                             axis.text = element_text(size = axis_font),
  #                             axis.title = element_text(size = axis_font),
  #                             plot.title = element_text(size = title_font))
  #        if (bpylim == "") {
  #          if (is.na(bpybreaks)) {
  #            barp <- barp + scale_y_continuous(expand = c(0, 0))
  #          }
  #          else if (!is.na(bpybreaks)) {
  #            barp <- barp + scale_y_continuous(limits=c(0,y_max),expand = c(0, 0),breaks=seq(0,y_max,bpybreaks),oob = rescale_none)
  #          }
  #        }
  #        else if (bpylim != "") {
  #          if (is.na(bpybreaks)) {
  #            barp <- barp + scale_y_continuous(limits=c(y_min,y_max),expand = c(0, 0),oob = rescale_none)
  #          }
  #          else if (!is.na(bpybreaks)) {
  #            barp <- barp + scale_y_continuous(limits=c(y_min,y_max),expand = c(0, 0),breaks=seq(y_min,y_max,bpybreaks),oob = rescale_none)
  #          }
  #        }
  #        
  #        if (dotChoice == T) {
  #          barp <- barp + geom_dotplot(data = expr_gene2, aes(x=Type,y=GeneName),
  #                                      binaxis='y', stackdir='center',
  #                                      stackratio=1, dotsize=dotsizein, fill = "black")
  #          #barp <- barp + geom_jitter(data = expr_gene2, aes(x=Type,y=GeneName),
  #          #                            fill = "black")
  #          #barp <- barp + geom_point(data=expr_gene2, aes(x=Type, y=GeneName),fill="black")
  #        }
  #        
  #        if (colorin != "") {
  #          colorin <- strsplit(colorin," ")[[1]]
  #          if (length(colorin) == 1) {
  #            colorin <- rep(colorin,length(unique(metaSub[,groupCriteria])))
  #          }
  #          barp <- barp + scale_fill_manual(values=colorin)
  #        }
  #        
  #        if (StatMethod != "None") {
  #          barp <- barp + stat_compare_means(data = expr_gene2, aes(x=Type,y=GeneName), method = StatMethod)
  #        }
  #        
  #        barp <- barp + coord_cartesian(clip = "off")
  #        
  #        
  #        
  #        barp
  #      }
  #    #}
  #  }
  #  
  #  
  #  
  #  
  #})
  
  output$barplot <- renderPlot({
    
    bpplot <- barplot_react()
    bpplot
    
  })
  
  ####----Enrichment Test----####
  
  output$EnrichFisherTab <- renderTable({
    
    meta <- meta_react()
    Feat1 <- input$EnrichFeat1
    Feat2 <- input$EnrichFeat2
    Var1 <- input$EnrichVar1
    Var2 <- input$EnrichVar2
    
    inBoth <- nrow(meta[which(meta[,Feat1] == Var1 & meta[,Feat2] == Var2),])
    inNone <- nrow(meta[which(meta[,Feat1] != Var1 & meta[,Feat2] != Var2),])
    inVar1 <- nrow(meta[which(meta[,Feat1] == Var1 & meta[,Feat2] != Var2),])
    inVar2 <- nrow(meta[which(meta[,Feat1] != Var1 & meta[,Feat2] == Var2),])
    
    FishMat <- as.data.frame(matrix(c(inBoth,inVar1,inVar2,inNone), nrow = 2))
    colnames(FishMat) <- c(paste0(Feat1," - ",Var1),paste0(Feat1," - Not ",Var1))
    rownames(FishMat) <- c(paste0(Feat2," - ",Var2),paste0(Feat2," - Not ",Var2))
    
    FishMat
    
  }, rownames = TRUE)
  
  Fisher_react <- reactive({
    
    meta <- meta_react()
    Feat1 <- input$EnrichFeat1
    Feat2 <- input$EnrichFeat2
    Var1 <- input$EnrichVar1
    Var2 <- input$EnrichVar2
    TailChoice <- input$FisherTailChoice
    
    inBoth <- nrow(meta[which(meta[,Feat1] == Var1 & meta[,Feat2] == Var2),])
    inNone <- nrow(meta[which(meta[,Feat1] != Var1 & meta[,Feat2] != Var2),])
    inVar1 <- nrow(meta[which(meta[,Feat1] == Var1 & meta[,Feat2] != Var2),])
    inVar2 <- nrow(meta[which(meta[,Feat1] != Var1 & meta[,Feat2] == Var2),])
    
    FishMat <- matrix(c(inBoth,inVar1,inVar2,inNone), nrow = 2)
    Fisher_react <- fisher.test(FishMat, alternative = TailChoice)
    Fisher_react
    
  })
  
  output$EnrichFisherOut <- renderText({
    
    Feat1 <- input$EnrichFeat1
    Feat2 <- input$EnrichFeat2
    Var1 <- input$EnrichVar1
    Var2 <- input$EnrichVar2
    
    text <- capture.output(Fisher_react())
    text[grep("data:  FishMat",text)] <- paste0("data:  ",Feat1," - ",Var1," vs. ",Feat2," - ",Var2)
    textOut <- paste(text,collapse = '\n')
    textOut
    
  })
  
  FisherVenn_react <- reactive({
    
    meta <- meta_react()
    Feat1 <- input$EnrichFeat1
    Feat2 <- input$EnrichFeat2
    Var1 <- input$EnrichVar1
    Var2 <- input$EnrichVar2
    
    inVar1 <- meta[which(meta[,Feat1] == Var1),1]
    inVar2 <- meta[which(meta[,Feat2] == Var2),1]
    
    VennList <- list(inVar1 = inVar1,
                     inVar2 = inVar2)
    names(VennList)[1] <- paste0(Feat1,":\n",Var1)
    names(VennList)[2] <- paste0(Feat2,":\n",Var2)
    VennObj <- Venn(VennList)
    VennData <- process_data(VennObj)
    
    ggplot() +
      # 1. region count layer
      geom_sf(aes(fill = count), data = venn_region(VennData), show.legend = F) +
      # 2. set edge layer
      geom_sf(aes(color = name), data = venn_setedge(VennData), show.legend = F, size = 1) +
      # 3. set label layer
      geom_sf_text(aes(label = name), data = venn_setlabel(VennData), size = 6) +
      # 4. region label layer
      geom_sf_label(aes(label = paste0(count)), 
                    data = venn_region(VennData),
                    size = 6,
                    alpha = 0.5) +
      scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")+
      scale_color_manual(values = c("#F4FAFE","#F4FAFE"))+
      theme_void()
    
  })
  
  output$EnrichVennPlot <- renderPlot({
    
    p <- FisherVenn_react()
    p
    
  })
  
  
  
  output$EnrichFisherText <- renderUI({
    
    meta <- meta_react()
    Feat1 <- input$EnrichFeat1
    Feat2 <- input$EnrichFeat2
    Var1 <- input$EnrichVar1
    Var2 <- input$EnrichVar2
    TailChoice <- input$FisherTailChoice
    FisherOut <- capture.output(Fisher_react())
    pval <- gsub("p-value [[:punct:]] ","",grep("^p-value ",FisherOut, value = T))
    if (as.numeric(pval) < 0.05) {
      SigWord <- "significantly"
    } else {SigWord <- "not significantly"}
    OR <- trimws(FisherOut[grep("^odds ratio",FisherOut)+1])
    if (OR < 1) {
      ORWord <- "depleted"
    } else {ORWord <- "enriched"}
    
    inBoth <- nrow(meta[which(meta[,Feat1] == Var1 & meta[,Feat2] == Var2),])
    inNone <- nrow(meta[which(meta[,Feat1] != Var1 & meta[,Feat2] != Var2),])
    inVar1 <- nrow(meta[which(meta[,Feat1] == Var1 & meta[,Feat2] != Var2),])
    inVar2 <- nrow(meta[which(meta[,Feat1] != Var1 & meta[,Feat2] == Var2),])
    
    HTML(paste("Fishers Exact test found that <b>",Feat1,"</b> - <b>",Var1,"</b> is <b>",SigWord,"</b> <b>",ORWord,
               "</b> in <b>",Feat2,"</b> - <b>",Var2,"</b>, based on a P.value of <b>",pval,"</b> and an odds ratio of <b>",OR,"</b>.", sep = ""))
    
  })
  
  
  ####----Heatmap----####
  
  AvgExprdf_react <- reactive({
    
    sampSubset <- input$HeatSampSubset
    sampCrit <- input$HeatSampCriteria
    sampGrouping <- input$HeatGroupCriteria
    sampGroupCond <- input$HeatGroupSelect
    GenesSelect <- input$HeatGeneSelect
    meta <- meta_react()

    meta <- meta[which(meta[,sampGrouping] %in% sampGroupCond),]
    
    if (sampSubset != "Select All") {
      meta <- meta[which(meta[,sampSubset] == sampCrit),]
      expr <- expr[,meta[,1]]
    }
    
    AvgExprDF <- data.frame(rownames(expr))
    for (i in unique(meta[,sampGrouping])) {
      samples <- meta[which(meta[,sampGrouping] == i),1]
      if (length(samples) <= 1) {
        AvgExprDF[,paste(sampGrouping,i, sep = "_")] <- expr[,samples]
      }
      else if (length(samples) > 1) {
        AvgExprDF[,paste(sampGrouping,i, sep = "_")] <- rowMeans(expr[,samples])
        
      }
    }
    rownames(AvgExprDF) <- AvgExprDF[,1]
    AvgExprDF <- AvgExprDF[,-1]
    AvgExprDF
    
  })
  
  AvgExprDotdf_react <- reactive({
    
    AvgExprDF <- AvgExprdf_react()
    GenesSelect <- input$HeatGeneSelect
    sampSubset <- input$HeatSampSubset
    sampCrit <- input$HeatSampCriteria
    sampGrouping <- input$HeatGroupCriteria
    sampGroupCond <- input$HeatGroupSelect
    sizeFactor <- input$DotHeatSizeFactor
    meta <- meta_react()
    meta <- meta[which(meta[,sampGrouping] %in% sampGroupCond),]
    if (sampSubset != "Select All") {
      meta <- meta[which(meta[,sampSubset] == sampCrit),]
    }
    
    AvgExprDF <- AvgExprDF[which(rownames(AvgExprDF) %in% GenesSelect),]
    AvgExprDF <- AvgExprDF[which(rowSums(AvgExprDF) > 0),]
    GenesSelect <- rownames(AvgExprDF)
    
    AvgExprDF2 <- AvgExprDF
    AvgExprDF2 <- log2(AvgExprDF2 + 1)
    zAvgExprDF2 <- apply(AvgExprDF2, 1, scale)
    zAvgExprDF2 <- apply(zAvgExprDF2, 1, rev)
    colnames(zAvgExprDF2) <- names(AvgExprDF2)
    AvgExprDF2 <- as.matrix(zAvgExprDF2)
    AvgExprDF2[is.na(AvgExprDF2)] <- 0
    
    AvgExprDF2 <- as.data.frame(AvgExprDF2)
    AvgExprDF2$gene <- rownames(AvgExprDF2)
    AvgExprDF_melt <- melt(AvgExprDF2[GenesSelect,])
    colnames(AvgExprDF_melt) <- c("gene",sampGrouping,"Expression")
    
    TotalCells_PerCluster <- as.data.frame(table(meta[,sampGrouping]))
    TotalCells_PerCluster[,1] <- paste(sampGrouping,TotalCells_PerCluster[,1], sep = "_")
    colnames(TotalCells_PerCluster) <- c(sampGrouping,"TotalCells_PerCluster")
    AvgExprDF_melt2 <- merge(AvgExprDF_melt,TotalCells_PerCluster)
    
    cellExpr <- function(vec) {
      prop <- sum(vec>0)
      return(prop)
    }
    
    AvgExprDF_melt3 <- AvgExprDF_melt2[1,]
    AvgExprDF_melt3$TotalClusterCells_Expressed <- NA
    for (cluster in unique(meta[,sampGrouping])) {
      AvgExprDF_melt2_sub <- AvgExprDF_melt2[which(AvgExprDF_melt2[,sampGrouping] == paste0(sampGrouping,"_",cluster)),]
      exprSub <- expr[GenesSelect,meta[which(meta[,sampGrouping] == cluster),1]]
      TotalClusterCells_Expressed <- as.data.frame(apply(exprSub,1,cellExpr))
      AvgExprDF_melt2_sub2 <- merge(AvgExprDF_melt2_sub,TotalClusterCells_Expressed, by.x = "gene", by.y = 0)
      colnames(AvgExprDF_melt2_sub2)[ncol(AvgExprDF_melt2_sub2)] <- "TotalClusterCells_Expressed"
      AvgExprDF_melt3 <- rbind(AvgExprDF_melt3,AvgExprDF_melt2_sub2)
      
    }
    AvgExprDF_melt3 <- AvgExprDF_melt3[-1,]
    
    AvgExprDF_melt3
    
  })
  
  HeatmapDotAvg_react <- reactive({
    
    AvgExprDF_melt3 <- AvgExprDotdf_react()
    sampGrouping <- input$HeatGroupCriteria
    sampGroupCond <- input$HeatGroupSelect
    ClusterMethod <- input$ClusterMethodChoice
    ClusterColChoice <- input$ClusterCol
    ClusterRowChoice <- input$ClusterRow
    GenesSelect <- input$HeatGeneSelect
    
    mat <- AvgExprDF_melt3 %>% 
      select(-TotalCells_PerCluster, -TotalClusterCells_Expressed) %>%  # drop unused columns to faciliate widening
      pivot_wider(names_from = !!sym(sampGrouping), values_from = Expression) %>% 
      data.frame() # make df as tibbles -> matrix annoying
    row.names(mat) <- mat$gene  # put gene in `row`
    mat <- mat[,-1] #drop gene column as now in rows
    clust <- hclust(dist(mat %>% as.matrix()), method = ClusterMethod) # hclust with distance matrix
    v_clust <- hclust(dist(mat %>% as.matrix() %>% t()), method = ClusterMethod) # hclust with distance matrix
    
    ddgram_row <- as.dendrogram(clust) # create dendrogram
    ggtree_plot <- ggtree::ggtree(ddgram_row) +
      geom_tiplab(align = T, linetype = "solid")
    ggtree_plot[["data"]][["label"]] <- rep("",length(ggtree_plot[["data"]][["label"]]))
    
    ddgram_col <- as.dendrogram(v_clust)
    ggtree_plot_col <- ggtree(ddgram_col) +
      layout_dendrogram() +
      geom_tiplab(align = T, linetype = "solid")
    ggtree_plot_col[["data"]][["label"]] <- rep("",length(ggtree_plot_col[["data"]][["label"]]))
    
    if (ClusterColChoice & ClusterRowChoice) {
      AvgExprDF_melt4 <- AvgExprDF_melt3 %>%
        mutate(`% Expressing` = (TotalClusterCells_Expressed/TotalCells_PerCluster) * 100,
               gene = factor(gene, levels = clust$labels[clust$order]),
               Grouping = factor(!!sym(sampGrouping), levels = v_clust$labels[v_clust$order]))
    } else if (ClusterColChoice & !ClusterRowChoice) {
      AvgExprDF_melt4 <- AvgExprDF_melt3 %>%
        mutate(`% Expressing` = (TotalClusterCells_Expressed/TotalCells_PerCluster) * 100,
               Grouping = factor(!!sym(sampGrouping), levels = v_clust$labels[v_clust$order]))
      AvgExprDF_melt4[,2] <- factor(AvgExprDF_melt4[,2], levels=rev(unique(GenesSelect)))
    } else if (!ClusterColChoice & ClusterRowChoice) {
      AvgExprDF_melt4 <- AvgExprDF_melt3 %>%
        mutate(`% Expressing` = (TotalClusterCells_Expressed/TotalCells_PerCluster) * 100,
               gene = factor(gene, levels = clust$labels[clust$order]))
      ClusOrder <- paste0(sampGrouping,"_",sampGroupCond)
      AvgExprDF_melt4[,1] <- factor(AvgExprDF_melt4[,1], levels=unique(ClusOrder))
    } else if (!ClusterColChoice & !ClusterRowChoice) {
      AvgExprDF_melt4 <- AvgExprDF_melt3 %>%
        mutate(`% Expressing` = (TotalClusterCells_Expressed/TotalCells_PerCluster) * 100)
      AvgExprDF_melt4[,2] <- factor(AvgExprDF_melt4[,2], levels=rev(unique(GenesSelect)))
      ClusOrder <- paste0(sampGrouping,"_",sampGroupCond)
      AvgExprDF_melt4[,1] <- factor(AvgExprDF_melt4[,1], levels=unique(ClusOrder))
    }
    dotplot <- AvgExprDF_melt4 %>%
      ggplot(aes(x=AvgExprDF_melt4[,1], y = gene, color = Expression, size = `% Expressing`)) + 
      geom_point() +
      theme_minimal() +
      scale_colour_gradient2(low = "dark blue", high = "dark red", mid = "white") +
      theme(axis.text.x = element_text(angle = 90, size = 16),
            axis.text.y = element_text(size = 14),
            axis.title.x = element_blank(),
            axis.title.y = element_blank()) +
      scale_y_discrete(position = "right")
    
    if (ClusterColChoice & ClusterRowChoice) {
      plot_spacer() + plot_spacer() + ggtree_plot_col +
        ggtree_plot + plot_spacer() + dotplot + 
        plot_layout(ncol = 3, widths = c(0.7, -0.1, 4), heights = c(0.9, 4))
    } else if (ClusterColChoice & !ClusterRowChoice) {
      plot_spacer() + ggtree_plot_col +
        plot_spacer() + dotplot + 
        plot_layout(ncol = 2, widths = c(-0.1, 4), heights = c(0.9, 4))
    } else if (!ClusterColChoice & ClusterRowChoice) {
        ggtree_plot + plot_spacer() + dotplot + 
        plot_layout(ncol = 3, widths = c(0.7, -0.1, 4), heights = c(4))
    } else if (!ClusterColChoice & !ClusterRowChoice) {
      dotplot
    }
    
  })
  
  HeatmapAvg_react <- reactive({
    
    if (length(input$HeatGeneSelect) >= 2) {
      
      sampGrouping <- input$HeatGroupCriteria
      sampGroupCond <- input$HeatGroupSelect
      ClusterMethod <- input$ClusterMethodChoice
      ClusterColChoice <- input$ClusterCol
      ClusterRowChoice <- input$ClusterRow
      color_choice <- input$ColorPaletteHeat
      GenesSelect <- input$HeatGeneSelect
      meta <- meta_react()
      AvgExprDF <- AvgExprdf_react()
      
      AvgExprDF <- AvgExprDF[which(rownames(AvgExprDF) %in% GenesSelect),]
      AvgExprDF <- AvgExprDF[which(rowSums(AvgExprDF) > 0),]
      GenesSelect <- rownames(AvgExprDF)
      
      dataset <- AvgExprDF
      dataset <- log2(dataset + 1)
      zdataset <- apply(dataset, 1, scale)
      zdataset <- apply(zdataset, 1, rev)
      colnames(zdataset) <- names(dataset)
      dataset <- as.matrix(zdataset)
      dataset[is.na(dataset)] <- 0
      
      anno_meta <- data.frame(colnames(AvgExprDF))
      anno_meta$Type <- gsub("AvgExpr_","",anno_meta[,1])
      rownames(anno_meta) <- anno_meta[,1]
      anno_meta <- anno_meta[,-1,drop = F]
      
      dataset = dataset[apply(dataset, 1, function(x) !all(x==0)),]
      minimum = -5;
      maximum = 5;
      if (abs(min(dataset)) > abs(max(dataset))) {
        dataset[dataset < -abs(max(dataset))] = -abs(max(dataset))
      } else {
        dataset[dataset > abs(min(dataset))] = abs(min(dataset))
      }
      
      bk = c(seq(minimum,minimum/2, length=100), seq(minimum/2,maximum/2,length=100),seq(maximum/2,maximum,length=100))
      
      HeatMap_Colors <- c("dark blue","blue","white","red", "dark red")
      hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
      dataset <- dataset[match(GenesSelect, rownames(dataset)),]
      ClusOrder <- paste0(sampGrouping,"_",sampGroupCond)
      dataset <- dataset[,match(ClusOrder, colnames(dataset))]
      
      hm <- pheatmap::pheatmap(as.matrix(dataset),
                               cluster_col = ClusterColChoice,
                               cluster_row = ClusterRowChoice,
                               fontsize_row = 12,
                               fontsize_col = 12,
                               clustering_method = ClusterMethod,
                               color=hmcols,
                               annotation_col = anno_meta,
                               angle_col = 90,
                               border_color = NA)
      
    }
    
    
  })
  
  output$HeatmapAvg <- renderPlot({
    
    if (input$DotHeatChoice) {
      plot <- HeatmapDotAvg_react()
    } else {
      plot <- HeatmapAvg_react()
    }
    plot
    
  })
  
  
  
  ####----Downloaders----####
  
  ## Pre-calculated downloads
  
  output$dnldUMAP_SVG_PreC_clin1 <- downloadHandler(
    filename = function() {
      umap1 <- input$SelectPreCalc1
      umap2 <- input$SelectPreCalc2
      metaCol1 <- input$UMAPannotateSamps1
      metaCol1 <- gsub("[[:punct:]]",".",metaCol1)
      metaCol1 <- gsub(" ",".",metaCol1)
      paste(Project_Name2,"_",umap1,"_",umap2,"_",metaCol1,"Anno_UMAP.svg", sep = '')
    },
    content = function(file) {
      plot <- umap_plot_PreC_clin_react_base1()
      #TrajCluster <- input$TrajClustSelect
      #umap1 <- input$SelectPreCalc1
      #umap2 <- input$SelectPreCalc2
      #if (input$AddTrajLines) {
      #  if (!is.null(TrajCluster)) {
      #    expr_mat <- as.matrix(expr)
      #    sce <- SingleCellExperiment(assays = List(counts = expr_mat))
      #    assays(sce)$norm <- assays(sce)$counts
      #    
      #    rd_umap <- meta[,c("SampleName",umap1,umap2)]
      #    rownames(rd_umap) <- rd_umap[,1]
      #    rd_umap <- rd_umap[,-1]
      #    rd_umap <- as.matrix(rd_umap)
      #    reducedDims(sce) <- SimpleList(UMAP = rd_umap)
      #    
      #    cluster_info <- meta[,TrajCluster]
      #    names(cluster_info) <- meta[,1]
      #    colData(sce)$cluster_info <- cluster_info
      #    
      #    sce_umap <- slingshot(sce, clusterLabels = 'cluster_info', reducedDim = 'UMAP')
      #    
      #    for (i in seq_along(slingCurves(sce_umap))) {
      #      curve_i <- slingCurves(sce_umap)[[i]]
      #      curve_i <- curve_i$s[curve_i$ord, ]
      #      colnames(curve_i) <- c("UMAP_1", "UMAP_2")
      #      curve_i <- as.data.frame(curve_i)
      #      plot <- plot + geom_path(data = as.data.frame(curve_i), col = "black", linewidth = 1)
      #    }
      #  }
      #}
      ggsave(file,plot, width = 8, height = 8)
    }
  )
  
  output$dnldUMAP_SVG_PreC_expr <- downloadHandler(
    filename = function() {
      umap1 <- input$SelectPreCalc1
      umap2 <- input$SelectPreCalc2
      GeneSelected <- input$GeneSelection
      paste(Project_Name2,"_",umap1,"_",umap2,"_",GeneSelected,"Expr_UMAP.svg", sep = '')
    },
    content = function(file) {
      plot <- umap_plot_PreC_expr_react_base()
      ggsave(file,plot, width = 8, height = 8)
    }
  )
  
  output$dnldUMAP_SVG_PreC_clin2 <- downloadHandler(
    filename = function() {
      umap1 <- input$SelectPreCalc1
      umap2 <- input$SelectPreCalc2
      metaCol2 <- input$UMAPannotateSamps2
      metaCol2 <- gsub("[[:punct:]]",".",metaCol2)
      metaCol2 <- gsub(" ",".",metaCol2)
      paste(Project_Name2,"_",umap1,"_",umap2,"_",metaCol2,"Anno_UMAP.svg", sep = '')
    },
    content = function(file) {
      plot <- umap_plot_PreC_clin_react_base2()
      ggsave(file,plot, width = 8, height = 8)
    }
  )
  
  output$dnldPlotSVG_Barplot <- downloadHandler(
    filename = function() {
      umap1 <- input$SelectPreCalc1
      umap2 <- input$SelectPreCalc2
      sampSubset <- input$BPsampSubset
      sampCriteria <- input$BPsampCriteria
      groupCriteria <- input$BPgroupCriteria
      geneSelected <- input$BPgeneSelection
      
      sampSubset <- gsub("[[:punct:]]",".",sampSubset)
      sampSubset <- gsub(" ",".",sampSubset)
      sampCriteria <- gsub("[[:punct:]]",".",sampCriteria)
      sampCriteria <- gsub(" ",".",sampCriteria)
      groupCriteria <- gsub("[[:punct:]]",".",groupCriteria)
      groupCriteria <- gsub(" ",".",groupCriteria)
      
      paste(Project_Name2,"_",sampSubset,"_",sampCriteria,"_",groupCriteria,"_",geneSelected,"_Expression_BarPlot.svg", sep = '')
    },
    content = function(file) {
      plot <- barplot_react()
      ggsave(file,plot, width = 8, height = 8)
    }
  )
  
  output$dnldPlotPDF_Barplot <- downloadHandler(
    filename = function() {
      umap1 <- input$SelectPreCalc1
      umap2 <- input$SelectPreCalc2
      sampSubset <- input$BPsampSubset
      sampCriteria <- input$BPsampCriteria
      groupCriteria <- input$BPgroupCriteria
      geneSelected <- input$BPgeneSelection
      
      sampSubset <- gsub("[[:punct:]]",".",sampSubset)
      sampSubset <- gsub(" ",".",sampSubset)
      sampCriteria <- gsub("[[:punct:]]",".",sampCriteria)
      sampCriteria <- gsub(" ",".",sampCriteria)
      groupCriteria <- gsub("[[:punct:]]",".",groupCriteria)
      groupCriteria <- gsub(" ",".",groupCriteria)
      
      paste(Project_Name2,"_",sampSubset,"_",sampCriteria,"_",groupCriteria,"_",geneSelected,"_Expression_BarPlot.pdf", sep = '')
    },
    content = function(file) {
      plot <- barplot_react()
      ggsave(file,plot, width = 8, height = 8)
    }
  )
  
  
  
}




# Run the application 
shinyApp(ui = ui, server = server)


