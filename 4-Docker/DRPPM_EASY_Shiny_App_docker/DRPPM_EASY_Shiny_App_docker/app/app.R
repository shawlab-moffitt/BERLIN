


####----Install and load packages----####

packages <- c("shiny","shinythemes","shinyjqui","pheatmap","BiocManager","RColorBrewer",
              "dplyr","DT","ggplot2","ggpubr","reshape2","tibble","viridis","scales",
              "plotly","readr","enrichR","ggrepel","tidyr","tools","shinycssloaders")

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
invisible(lapply(packages, library, character.only = TRUE))
#bioconductor packages
bioCpacks <- c("clusterProfiler","GSVA","limma","enrichplot")
installed_packages_BIOC <- bioCpacks %in% rownames(installed.packages())
if (any(installed_packages_BIOC == FALSE)) {
  BiocManager::install(bioCpacks[!installed_packages_BIOC], ask = F)
}
invisible(lapply(bioCpacks, library, character.only = TRUE))




####----User Data Input----####

#Input desired project name for webpage - will be followed by 'Expression Analysis'
ProjectName <- 'GSE116256 AMLscRNA'


##--User Input File Names--##

#expression data
expr_file <- 'Example_Input_Data/GSE116256_AMLscRNA_1000_RNA_normalized_counts.txt'

#meta data
meta_file <- 'Example_Input_Data/GSE116256_AMLscRNA_1000_RNA_metafile.txt'
#Is there a header?
header <- TRUE

#If human: set TRUE
#If mouse: set FALSE
human <- TRUE

# Starting Feature
Feature_Selected <- 'seurat_clusters'

set.seed(10242022)










##--User Gene Set Input (Optional)--##

#write in the name of your gene set list for shiny UI
userGSlist_name <- 'LINCS L1000'

#path to your gene set file .gmt or .txt/.tsv
userGS_file <- 'Gene_Set_Data/LINCS_L1000_gsNsym_HS_v2.zip'
#Does gene set file have header?
header.gs <- TRUE

#path to your R data list object for ssGSEA
userRData_file <- 'Gene_Set_Data/LINCS_L1000_gs_HS_v2.RData'





####----Backend Data Input----####

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

if (human == TRUE) {
  #MSigDB gene set
  msigdb <- 'Gene_Set_Data/msigdb_gsNsym_HS_v2.zip'
  #MSigDB gene set FOR UI
  msigdb2 <- 'Gene_Set_Data/msigdb_gsNcat_HS_v2.txt'
  #gene set list for ssGSEA
  gs <- loadRData('Gene_Set_Data/msigdb_gs_HS_v2.RData')
  #Cytokine genes for human
  CTKgenes <- c("IL2","IL12A","IL12B","IL17A","IFNA1","IFNB1","IFNG","IFNGR","CD11b",
                "ITGAM","CD33","ENTPD1","ICOSLG","CD275","CD278","TNFSF9","TNFRSF9",
                "CD40","CD40LG","CD70","CD27","TNFSF18","TNFRSF18","TNFSF14","TNFRSF14",
                "TNFSF4","TNFRSF4","HLA-A","CD3","CEACAM1","CD80","CD86","CTLA4","CD276",
                "VTCN1","PVR","CD226","TIGIT","CD96","LGALS3","LGALS3BP","LGALS9","LGALS9C",
                "HAVCR2","HHLA2","TMIGD2","CD274","PDCD1LG2","PDCD1","VSIR")
}
if (human == FALSE) {
  #MSigDB gene set
  msigdb <- 'Gene_Set_Data/msigdb_gsNsym_MM.zip'
  #MSigDB gene set FOR UI
  msigdb2 <- 'Gene_Set_Data/msigdb_gsNcat_MM.tsv'
  #gene set list for ssGSEA 
  gs <- loadRData('Gene_Set_Data/msigdb_gs_MM.RData')
  #Cytokine genes for mouse
  CTKgenes <- c("Il2","Il12a","Il12b","Il17a","Ifna13","Ifnb1","Ifng","Ifngr1","Cd11b","Itgam",
                "Cd33","Entpd1","Icosl","Icos","Tnfsf9","Tnfrsf9","Cd40","Cd40lg","Cd70","Cd27",
                "Tnfsf18","Tnfrsf18","Tnfsf14","Tnfrsf14","Tnfsf4","Tnfrsf4","H2-K1","CD3G",
                "Ceacam1","Cd80","Cd86","Ctla4","Cd276","Vtcn1","Pvr","Cd226","Tigit","Cd96","Lgals3",
                "Lgals3bp","Lgals9","Lgals9c","Havcr2","Hhla2","Cd274","Pdcd1lg2","Pdcd1","Vsir")
}



####----Read and Manipulate Files----####


##read files

##reading expression data
expr <- as.data.frame(read_delim(expr_file, delim = '\t', col_names = T))
colnames(expr)[1] <- "Gene"
# Remove Expression with NA
expr <- expr %>%
  drop_na()
# Check that expression data is numeric
isChar <- unname(which(sapply(expr, function(x) is.character(x))))
isChar <-  isChar[-1]
if (length(isChar) > 0) {
  expr[isChar] <- sapply(expr[isChar],as.numeric)
}
# Remove Duplicate genes
if (TRUE %in% duplicated(expr[,1])) {
  expr <- expr %>%
    group_by(Gene) %>%
    summarise_all(max)
}
expr <- as.data.frame(expr)
# Make rownames <- genenames
row.names(expr) <- expr[,1]
expr <- expr[,-1]
expr = expr[order(row.names(expr)), ]
# Harmonize special characters in sample names
colnames(expr) <- gsub("[[:punct:]]", "_", colnames(expr))
A <- as.matrix(expr)

#gene list file from expression data
Gene <- rownames(expr)
geneList <- as.data.frame(Gene)


#meta
meta <- read.delim(meta_file, sep = '\t', header = header, strip.white = T)
meta[,1] <- gsub("[_.-]", "_", meta[,1])
colnames(meta)[1] <- "SampleName"

metagroups <- as.vector(levels(factor(meta[,2])))



#for heatmap sample selection
sampsames <- intersect(colnames(expr),meta[,1])
#ensure expression samples and meta are exact
expr <- expr[,sampsames]
meta <- meta[which(meta[,1] %in% sampsames),]


#boxplot choices based on meta groups
if (length(metagroups) == 2) {
  boxopt <- c("wilcox.test", "t.test", "none")
}
if (length(metagroups) >= 3) {
  boxopt <- c("kruskal.test", "anova", "none")
}


#MSigDB gene sets
msigdb.gsea <- as.data.frame(read_delim(msigdb, delim = '\t'))
gmt <- msigdb.gsea
#MSigDB gene sets FOR UI
msigdb.gsea2 <- as.data.frame(read_delim(msigdb2, delim = '\t'))


#tab2 User gene set
if (file_ext(userGS_file) == "gmt") {
  tab2 <- read.gmt(userGS_file)
}
if (file_ext(userGS_file) != "gmt") {
  tab2 <- as.data.frame(read_delim(userGS_file, col_names = T, delim = '\t'))
}
#tab2 back end
GeneSet2 <- as.data.frame(unique(tab2[,1]))
rownames(GeneSet2) <- 1:nrow(GeneSet2)
colnames(GeneSet2)[1] <- "Gene_Set"
#tab2 R Data list
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
gs2 <- loadRData(userRData_file)

## starting ssgsea heatmap genes
#gs_names_start <- c(grep("^HALLMARK",names(gs),value = T),grep("^HALLMARK",names(gs),value = T, invert = T))
#gs_names_start <- gs_names_start[1:3]

#CV function for variance
cv <- function(x){
  (sd(x)/mean(x))*100
}

##enrichR pathway choices
#enrichRtab <- listEnrichrDbs()
#enrichRchoice <- enrichRtab[,"libraryName"]

####----Shiny UI----####


shinytheme("sandstone")

ui <-
  
  navbarPage(paste("{",ProjectName,"Expression Analysis }", sep=" "),
             
             ####----Intro Tab----####
             
             tabPanel("Intro/Methodology",
                      fluidPage(
                        mainPanel(
                          h3("Introduction"),
                          p("With the rapid generation of large datasets as a result of advancing next-generation sequecing (NGS) technologies, the need for quick and reproducible data analysis tools has become paramount. The ability to analyze both genomic and proteomic data sets can allow for the detection of compelling genes and features that may be of interest for further analysis. The General Expression Analysis App hosted in R Shiny does not require an extensive computation background to run these analyses. With user input of expression and meta data, along with gene set we have sourced and provided, the user may produce useful anaylsis and visualizations on a web-based interface within minutes. This method of reproducible data analysis generates a number of visualization to view Gene Set Enrichment (GSEA) and differential gene expression analysis, as well as the ability to download various tables and gene sets that are produced by the App for further use. Additional apps are being developed for the EASY family. DRPPM-EASY-Integraction allows users to further analyze data obtained from the main DRPPM-EASY app and compare expression data between two matrices. DRPPM-EASY-CCLE integrates a sample selection tab which allows users to select expression and meta data from the Cancer Cell Line Encyclopedia (CCLE) based on cancer type or lineage and sample type for analysis within the DRPPM-EASY app. The flow chart below gives a layout of the EASY app's family infrastructure which we will describe in further detail."),
                          h3("Unsupervised Clustering Methods"),
                          p("Unsupervised clustering is performed by calculating the top variable genes through MAD, CV, or VAR as the variance measure. The choice of variance measure and number to top genes/probes can be determined by the user, as well as the number of clusters with k-means. The identified most variable genes can be downloaded as a table for the user. The “complete” method was chosen as the default algorithm for clustering, but other methods such as Ward, average, and centroid could be chosen based on the ‘hclust()’ function in R."),
                          h3("Differential Gene Expression Methods"),
                          p("Differential gene expression analysis is performed on the expression data between two groups defined by the provided metadata file. The samples chosen are log-transformed (log2 + 1). A model matrix is then designed for LIMMA linear modeling followed by empirical Bayes statistics to moderate gene variance and modeling of the global characteristic of the data. Of note, the current pipeline expects the matrix input quantification are pre-normalized, such as in CPM, TMM, or FPKM."),
                          h3("Gene Set Enrichment Analysis (GSEA) Methods"),
                          p("Gene Set Enrichment Analysis (GSEA) is performed through signal-to-noise calculations on the expression matrix. This calculation requires at least two types of sample groups and at least three samples for each grouping. The signal-to-noise calculation scales the difference of means by standard deviation and generates a ranked list of genes, effectively capturing the differences between phenotypes. The GSEA identifies the enrichment score (ES) for a gene set, which indicates the extent that the gene set is overrepresented at the beginning or end of the list of ranked genes. The enrichment score is calculated by walking down the ranked gene list and increasing the running-sum statistic when a gene is in the gene set and decreasing when it is not. The maximum deviation from zero that is encountered through walking the list is the ES, where a positive ES indicates the gene set is enriched at the top of the ranked list and a negative ES indicates the gene set is enriched at the bottom of the ranked list. A leading-edge gene list is provided displaying a subset of the genes in the gene set that contribute the most to the ES."),
                          h4("Gene Set Sources"),
                          p("The Molecular Signatures Database (MSigDB) is a collection of over 32,000 annotated gene sets divided into 9 major collections as well as various sub-collections [PMID: 16199517, PMID: 12808457]. The MSigDB  gene sets were downloaded with the msigdbr package and processed through R using Tidyverse packages to combine the gene sets into a data frame."),
                          p("The Cell Marker gene sets were derived from a comprehensive list of cell markers from cell types in various tissues obtained through manually curating over 100,000 published papers [PMID: 30289549]. This data was obtained as an extensive data frame and parsed to make gene sets with the tidyverse packages."),
                          p("The Library of Integrated Network-based Cellular Signatures (LINCS) L1000 gene sets were derived from the Connectivity Map (CMAP) projects data which treated 4 human cell lines with ~1300 drugs followed by profiling genome-wide mRNA expression following small molecule perturbation [PMID: 24906883].")
                        )
                      )
             ),
             
             ####----Data Exploration Tab----####
             
             tabPanel("Data Exploration",
                      fluidPage(
                        title = "Data Exploration",
                        sidebarLayout(
                          sidebarPanel(
                            width = 3,
                            
                            ####----MVG Heatmap----####
                            
                            conditionalPanel(condition = "input.dataset == '1'",
                                             conditionalPanel(condition = "input.datasetheat == '1'",
                                                              p(),
                                                              h4("Feature Parameters"),
                                                              uiOutput("rendMVGheatAnnoCol"),
                                                              #selectInput("MVGheatAnnoCol","Annotation Data",
                                                              #            choices = colnames(meta)[2:ncol(meta)]),
                                                              numericInput("NumFeatures", step = 1, label = "Number of Genes", value = 100),
                                                              hr(),
                                                              h4("Clustering Parameters"),
                                                              selectInput("VarianceMeasure", "Select Variance Measure",
                                                                          choices = c("MAD","CV","VAR")),
                                                              fluidRow(
                                                                column(6,
                                                                       checkboxInput("clustrowsMVG","Cluster Heatmap Rows", value = T)
                                                                ),
                                                                column(6,
                                                                       checkboxInput("clustcolsMVG","Cluster Heatmap Columns", value = T)
                                                                )
                                                              ),
                                                              uiOutput("ClusterMethodMVG"),
                                                              numericInput("NumClusters", step = 1, label = "Number of Clusters (Cut Tree with ~k)", value = 2),
                                                              h4("Download Cluster Result:"),
                                                              downloadButton("downloadClusters", "Download .tsv"),
                                                              h4("Download Most Variable Genes List:"),
                                                              downloadButton("MVGdownload", "Download MVG .tsv"),
                                                              downloadButton("MVGdownloadgmt", "Download MVG .gmt"),
                                                              hr(),
                                                              h4("Figure Parameters"),
                                                              selectInput("ColorPalette1", "Select Color Palette:",
                                                                          choices = c("Red/Blue" = "original",
                                                                                      "OmniBlueRed" = "OmniBlueRed",
                                                                                      "LightBlue/BlackRed" = "LightBlueBlackRed",
                                                                                      "Green/Black/Red" = "GreenBlackRed",
                                                                                      "Yellow/Green/Blue" = "YlGnBu","Inferno" = "Inferno",
                                                                                      "Viridis" = "Viridis","Plasma" = "Plasma",
                                                                                      "Reds" = "OrRd","Blues" = "PuBu","Greens" = "Greens")),
                                                              fluidRow(
                                                                column(6,
                                                                       checkboxInput("ShowColNames1","Show Heatmap Column Names", value = T),
                                                                       numericInput("heatmapFont2.c", "Heatmap Column Font Size:",
                                                                                    min = 5, max = 75,
                                                                                    value = 12, step = 1)
                                                                ),
                                                                column(6,
                                                                       checkboxInput("ShowRowNames1","Show Heatmap Row Names", value = T),
                                                                       numericInput("heatmapFont2.r", "Heatmap Row Font Size:",
                                                                                    min = 5, max = 75,
                                                                                    value = 9, step = 1)
                                                                )
                                                              )
                                             )
                            ),
                            
                            ####----Custom Heatmap----####
                            
                            conditionalPanel(condition = "input.dataset == '1'",
                                             conditionalPanel(condition = "input.datasetheat == '2'",
                                                              p(),
                                                              h4("Selection Parameters"),
                                                              uiOutput("rendCustomheatAnnoCol"),
                                                              #selectInput("CustomheatAnnoCol","Annotation Data",
                                                              #            choices = colnames(meta)[2:ncol(meta)]),
                                                              selectInput("heatmapGeneSelec","Gene Selection:",
                                                                          choices = sort(as.vector(geneList[,1])),
                                                                          multiple = T, selected = CTKgenes),
                                                              textInput("userheatgenes", "Text Input of Gene List (space delimited):", value = ""),
                                                              selectInput("userheatsamp2", "Samples Selection:",
                                                                          choices = sampsames, multiple = T, selected = sampsames),
                                                              h4("Clustering Parameters"),
                                                              fluidRow(
                                                                column(6,
                                                                       checkboxInput("ClusRowOptCust","Cluster Heatmap Rows", value = FALSE)
                                                                ),
                                                                column(6,
                                                                       checkboxInput("ClusColOptCust","Cluster Heatmap Columns", value = FALSE))
                                                              ),
                                                              uiOutput("rendClustMethodsCust"),
                                                              h4("Figure Parameters"),
                                                              selectInput("ColorPalette2", "Select Color Palette:",
                                                                          choices = c("Red/Blue" = "original",
                                                                                      "OmniBlueRed" = "OmniBlueRed",
                                                                                      "LightBlue/BlackRed" = "LightBlueBlackRed",
                                                                                      "Green/Black/Red" = "GreenBlackRed",
                                                                                      "Yellow/Green/Blue" = "YlGnBu","Inferno" = "Inferno",
                                                                                      "Viridis" = "Viridis","Plasma" = "Plasma",
                                                                                      "Reds" = "OrRd","Blues" = "PuBu","Greens" = "Greens")),
                                                              fluidRow(
                                                                column(6,
                                                                       checkboxInput("ShowRowNames2","Show Heatmap Row Names", value = T),
                                                                       numericInput("heatmapFont3.r", "Heatmap Row Font Size:",
                                                                                    min = 5, max = 75,
                                                                                    value = 12, step = 1)
                                                                ),
                                                                column(6,
                                                                       checkboxInput("ShowColNames2","Show Heatmap Column Names", value = T),
                                                                       numericInput("heatmapFont3.c", "Heatmap Column Font Size:",
                                                                                    min = 5, max = 75,
                                                                                    value = 12, step = 1)
                                                                )
                                                              )
                                             )
                            ),
                            
                            ####----DEG Heatmap (per sample)----####
                            
                            conditionalPanel(condition = "input.dataset == '1'",
                                             conditionalPanel(condition = "input.datasetheat == '3'",
                                                              h4("Selection Parameters"),
                                                              uiOutput("rendMetaColDEGHeat"),
                                                              uiOutput("rendcomparisonA2_h"),
                                                              uiOutput("rendcomparisonB2_h"),
                                                              #selectInput("MetaColDEGHeat","Meta Column:",
                                                              #            choices = colnames(meta)[2:ncol(meta)]),
                                                              #selectInput("comparisonA2_h", "Comparison: GroupA",
                                                              #            choices = metagroups, selected = metagroups[1]),
                                                              #selectInput("comparisonB2_h", "Comparison: GroupB",
                                                              #            choices = metagroups, selected = metagroups[2]),
                                                              hr(),
                                                              numericInput("fc_cutoff_h", "LogFC Threshold (Absolute Value)",
                                                                           min = 0, max = 5, step = 0.1, value = 1),
                                                              numericInput("p_cutoff_h", "Significance Threshold P.Value:",
                                                                           min = 0, max = 10, step = 0.01, value = 0.05),
                                                              verbatimTextOutput("GenesAboveCutoff1"),
                                                              numericInput("top_x_h", "Number of Top Hits to Show on Heatmap:",
                                                                           value = 100, min = 0),
                                                              hr(),
                                                              h4("Clustering Parameters"),
                                                              fluidRow(
                                                                column(6,
                                                                       checkboxInput("ClusRowOptDEG","Cluster Heatmap Rows", value = TRUE)
                                                                ),
                                                                column(6,
                                                                       checkboxInput("ClusColOptDEG","Cluster Heatmap Columns", value = TRUE))
                                                              ),
                                                              uiOutput("rendClustMethodsDEG"),
                                                              hr(),
                                                              h4("Figure Parameters"),
                                                              selectInput("ColorPalette3", "Select Color Palette:",
                                                                          choices = c("Red/Blue" = "original",
                                                                                      "OmniBlueRed" = "OmniBlueRed",
                                                                                      "LightBlue/BlackRed" = "LightBlueBlackRed",
                                                                                      "Green/Black/Red" = "GreenBlackRed",
                                                                                      "Yellow/Green/Blue" = "YlGnBu","Inferno" = "Inferno",
                                                                                      "Viridis" = "Viridis","Plasma" = "Plasma",
                                                                                      "Reds" = "OrRd","Blues" = "PuBu","Greens" = "Greens")),
                                                              fluidRow(
                                                                column(6,
                                                                       checkboxInput("ShowRowNames3","Show Heatmap Row Names", value = T),
                                                                       numericInput("heatmapFont3.r.deg", "Heatmap Row Font Size:",
                                                                                    min = 5, max = 75,
                                                                                    value = 12, step = 1)
                                                                ),
                                                                column(6,
                                                                       checkboxInput("ShowColNames3","Show Heatmap Column Names", value = T),
                                                                       numericInput("heatmapFont3.c.deg", "Heatmap Column Font Size:",
                                                                                    min = 5, max = 75,
                                                                                    value = 12, step = 1)
                                                                )
                                                              )
                                             )
                            ),
                            
                            ####----Average Gene Expression Heatmap - Custom----####
                            
                            conditionalPanel(condition = "input.dataset == '1'",
                                             conditionalPanel(condition = "input.datasetheat == '5'",
                                                              h4("Selection Criteria"),
                                                              uiOutput("rendAvgHeatMetaCol"),
                                                              uiOutput("rendSampCondSelectionCust"),
                                                              #selectInput("SampCondSelectionCust","Sample Coniditon Selection",
                                                              #            choices = metagroups, selected = metagroups[c(1:length(metagroups))],
                                                              #            multiple = T),
                                                              selectizeInput("avgheatmapGeneSelec","Gene Selection:",
                                                                             choices = sort(as.vector(geneList[,1])),
                                                                             multiple = T, selected = CTKgenes),
                                                              textInput("avguserheatgenes", "Text Input of Gene List (space delimited):", value = ""),
                                                              h4("Clustering Parameters"),
                                                              fluidRow(
                                                                column(6,
                                                                       checkboxInput("AvgClusRowOptCust","Cluster Heatmap Rows", value = FALSE)
                                                                ),
                                                                column(6,
                                                                       checkboxInput("AvgClusColOptCust","Cluster Heatmap Columns", value = FALSE))
                                                              ),
                                                              uiOutput("rendAvgClustMethodsCust"),
                                                              h4("Figure Parameters"),
                                                              selectInput("ColorPalette5", "Select Color Palette:",
                                                                          choices = c("Red/Blue" = "original",
                                                                                      "OmniBlueRed" = "OmniBlueRed",
                                                                                      "LightBlue/BlackRed" = "LightBlueBlackRed",
                                                                                      "Green/Black/Red" = "GreenBlackRed",
                                                                                      "Yellow/Green/Blue" = "YlGnBu","Inferno" = "Inferno",
                                                                                      "Viridis" = "Viridis","Plasma" = "Plasma",
                                                                                      "Reds" = "OrRd","Blues" = "PuBu","Greens" = "Greens")),
                                                              fluidRow(
                                                                column(6,
                                                                       checkboxInput("ShowRowNames5","Show Heatmap Row Names", value = T),
                                                                       numericInput("heatmapFont3.r.deg.avg.cust", "Heatmap Row Font Size:",
                                                                                    min = 5, max = 75,
                                                                                    value = 12, step = 1)
                                                                ),
                                                                column(6,
                                                                       checkboxInput("ShowColNames5","Show Heatmap Row Names", value = T),
                                                                       numericInput("heatmapFont3.c.deg.avg.cust", "Heatmap Column Font Size:",
                                                                                    min = 5, max = 75,
                                                                                    value = 12, step = 1)
                                                                )
                                                              )
                                             )
                            ),
                            
                            ####----Gene scatter plot side panel----####
                            
                            conditionalPanel(condition = "input.dataset == '2'",
                                             p(),
                                             uiOutput("rendScatterPlotMetaCol"),
                                             selectInput("scatterG1","Select Gene 1",
                                                         choices = Gene),
                                             selectInput("scatterG2","Select Gene 2",
                                                         choices = Gene, selected = Gene[2]),
                                             checkboxInput("logask","Log2 Transform Expression Data"),
                                             h4("Figure Parameters"),
                                             fluidRow(
                                               column(6,
                                                      numericInput("GeneScatterTitleSize","Title Font Size",
                                                                   value = 14, step = 1)
                                               ),
                                               column(6,
                                                      numericInput("GeneScatterAxisSize","Axis Lable Font Size",
                                                                   value = 12, step = 1))
                                             )
                            ),
                            
                            ####----Boxplot side panel----####
                            
                            conditionalPanel(condition = "input.dataset == '3'",
                                             p(),
                                             uiOutput("rendBoxPlot1MetaCol"),
                                             #selectInput("BoxPlot1MetaCol","Select Sample Condition:",
                                             #            choices = colnames(meta)[2:ncol(meta)]),
                                             h4("Gene Selection:"),
                                             div(DT::dataTableOutput("GeneListTable2"), style = "font-size:10px"),
                                             h4("Figure Parameters"),
                                             fluidRow(
                                               column(4,
                                                      numericInput("boxplot1TitleSize","Title Font Size:",
                                                                   value = 20, step = 1)
                                               ),
                                               column(4,
                                                      numericInput("boxplot1AxisSize","Axis Font Size:",
                                                                   value = 16, step = 1)
                                               ),
                                               column(4,
                                                      numericInput("boxplotDot2", "Dot Size:",
                                                                   min = 0, max = 5,
                                                                   value = 0.5, step = .25)
                                               )
                                             )
                            ),
                            conditionalPanel(condition = "input.dataset == '4'",
                                             tabsetPanel(
                                               tabPanel("Data Parameters",
                                                        p(),
                                                        h4("Transform Data"),
                                                        uiOutput("rendBarPlot1MetaCol"),
                                                        fluidRow(
                                                          column(6, style = 'padding-right:4px;',
                                                                 checkboxInput("log2barplot","Log2(expr+1) Expression",value = T)
                                                          ),
                                                          column(6, style = 'padding-left:4px;',
                                                                 selectInput("errorbarplot","Error Bar Type",
                                                                             choices = c("Standard Deviation","Standard Error","None"))
                                                          )
                                                        ),
                                                        hr(),
                                                        h4("Gene Selection:"),
                                                        div(DT::dataTableOutput("GeneListTableBarPlot"), style = "font-size:10px"),
                                               ),
                                               tabPanel("Figure Paramters",
                                                        p(),
                                                        textInput("barplotColoCodes","Color Code(s):",value = "", placeholder = "HEX or R Color Code(s) (Space Delim)"),
                                                        hr(),
                                                        fluidRow(
                                                          column(6,
                                                                 checkboxInput("barplotsampledots","Include Dot Annotation", value = T)
                                                          ),
                                                          column(6,
                                                                 numericInput("barplotDotSize","Dot Size:", value = 1, step = 0.5)
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
                                                                 radioButtons("barplotXaxOrder","X-Axis Group Order",
                                                                              choices = c("Ascending","Descending"))
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
                            )
                          ),
                          mainPanel(
                            tabsetPanel(
                              id = "dataset",
                              tabPanel("Heatmaps",
                                       tabsetPanel(
                                         id = "datasetheat",
                                         tabPanel("Unsupervised Clustering Heatmap",
                                                  withSpinner(jqui_resizable(plotOutput("heatmap1", width = "100%", height = "1000px")), type = 6),
                                                  fluidRow(
                                                    downloadButton("dnldPlotSVG_heat1","Download as SVG"),
                                                    downloadButton("dnldPlotPDF_heat1","Download as PDF")
                                                  ),
                                                  value = 1),
                                         tabPanel("Custom Heatmap",
                                                  withSpinner(jqui_resizable(plotOutput("heatmap2", width = "100%", height = "1000px")), type = 6),
                                                  fluidRow(
                                                    downloadButton("dnldPlotSVG_heat2","Download as SVG"),
                                                    downloadButton("dnldPlotPDF_heat2","Download as PDF")
                                                  ),
                                                  value = 2),
                                         tabPanel("DEG Heatmap",
                                                  withSpinner(jqui_resizable(plotOutput("heatmap3", width = "100%", height = "1000px")), type = 6),
                                                  fluidRow(
                                                    downloadButton("dnldPlotSVG_heat3","Download as SVG"),
                                                    downloadButton("dnldPlotPDF_heat3","Download as PDF")
                                                  ),
                                                  value = 3),
                                         tabPanel("Custom Genes Average Expression Heatmap",
                                                  withSpinner(jqui_resizable(plotOutput("avgheatmap1Cust", width = "100%", height = "1000px")), type = 6),
                                                  fluidRow(
                                                    downloadButton("dnldPlotSVG_heat5","Download as SVG"),
                                                    downloadButton("dnldPlotPDF_heat5","Download as PDF")
                                                  ),
                                                  value = 5)
                                       ),
                                       value = 1
                              ),
                              tabPanel("Gene Scatter Plot",
                                       p(),
                                       withSpinner(jqui_resizable(plotlyOutput("geneScatter0", height = "500px")), type = 6),
                                       fluidRow(
                                         downloadButton("dnldPlotSVG_scatter","Download as SVG"),
                                         downloadButton("dnldPlotPDF_scatter","Download as PDF")
                                       ),
                                       DT::dataTableOutput("geneScatterTable"),
                                       downloadButton("geneScatterDownload", "Download Non-log2 Transformed .tsv"),
                                       value = 2),
                              tabPanel("Box Plot",
                                       p(),
                                       #withSpinner(jqui_resizable(plotOutput('boxplot3', width = "100%", height = "600px")), type = 6),
                                       jqui_resizable(plotOutput('boxplot3', width = "100%", height = "600px")),
                                       fluidRow(
                                         downloadButton("dnldPlotSVG_exprBox","Download as SVG"),
                                         downloadButton("dnldPlotPDF_exprBox","Download as PDF")
                                       ),
                                       value = 3),
                              tabPanel("Bar Plot",
                                       p(),
                                       jqui_resizable(plotOutput('barplot', width = "100%", height = "600px")),
                                       fluidRow(
                                         downloadButton("dnldPlotSVG_exprBar","Download as SVG"),
                                         downloadButton("dnldPlotPDF_exprBar","Download as PDF")
                                       ),
                                       value = 4)
                              
                            )
                          )
                        )
                      )
             ),
             
             ####----Differential Expression Analysis----####
             
             tabPanel("Differential Expression Analysis",
                      fluidPage(
                        title = "Differential Expression Analysis",
                        sidebarLayout(
                          sidebarPanel(
                            
                            ####----Volcano and MA plot side panel----####
                            
                            conditionalPanel(condition = "input.datasettwo == '1' || input.datasettwo == '2'",
                                             h4("Condition Selection:"),
                                             uiOutput("rendDEGMetaCol"),
                                             uiOutput("rendcomparisonA2"),
                                             uiOutput("rendcomparisonB2"),
                                             #selectInput("comparisonA2", "Comparison: GroupA",
                                             #            choices = metagroups, selected = metagroups[1]),
                                             #selectInput("comparisonB2", "Comparison: GroupB",
                                             #            choices = metagroups, selected = metagroups[2]),
                                             h4("Threshold Parameters"),
                                             numericInput("fc_cutoff", "LogFC Threshold (Absolute Value)",
                                                          min = 0, max = 5, step = 0.1, value = 1),
                                             numericInput("p_cutoff", "Significance Threshold P.Value):",
                                                          min = 0, max = 10, step = 0.1, value = 0.05),
                                             h4("Gene Selection Parameters:"),
                                             numericInput("top_x", "Number of Top Hits:", value = 10,
                                                          min = 0),
                                             selectizeInput("userGeneSelec", "User Selected Hits:",
                                                            choices = sort(as.vector(geneList[,1])), multiple = T, selected = "-"),
                                             textInput("userGeneSelec2", "Text Input of Gene List (space or tab delimited):", value = ""),
                                             #textInput("userGeneSelec2", "Text Input of Gene List (space or tab delimited):",
                                             #          value = "Bcl11b Tcf7 Bcl6 Id3 Id2 Ahr Ahrr Prdm1 Kit Klrg1 Klrb1c Klrb1b Pdcd1 Ctla4 Btla Itga1 Itgam Itgax S1pr1 S1pr5 Ccr7 Klf2 Gzmb Gzma Bcl2 Prf1 Zeb2 Cx3Cr1 Bhlhe40 Bach2 Slamf6 Eomes Icosl Nr4a1 Sell Klrk1 Klrc2 Klrc3 Klrd1 Klre1 Fcer1g Ncr1 Lef1 Havcr2 CD244 Cxcr5 Il2 Zfp683 Ly6c2 Cd27 Cd5 Cd160 Lamp1 Il7r NT5E Entpd1 Tigit Cd38 P2rx7 Klf4 Elovl6 Cpt1a Slc27a Slc2a1 Slc2a3 Pgk1 Lipa Pnpla2 Tcf7 Bach2 Satb1 Ccr7 Il7r Jund Gzmc Gzmf Gzmm Gzma Prf1 Klrk1 Klrb1a Klrc3 Klrb1c Klrc2 Klra5 Klre1 Klrg1 Ncr1 Cd244a Itgax Cd7 Itga1 Tox Zeb2 Mafb Tcf4 Irf5 Hdac9 Ctla4 Il10 Arg1 Tgfbi Cd300a Vav3 Themis2"),
                                             h4("Font Sizes"),
                                             fluidRow(
                                               column(4,
                                                      numericInput("VolMAAxisSize","Axis Lable",
                                                                   value = 18, step = 1)
                                               ),
                                               column(4,
                                                      numericInput("VolMATickSize","Axis Tick",
                                                                   value = 16, step = 1)
                                               ),
                                               column(4,
                                                      numericInput("VolMAAnnoSize","Annotation Text",
                                                                   value = 6, step = 1)
                                               )
                                               #column(3,
                                               #       selectInput("VolMAAnnoFont","Annotation Font",
                                               #                   choices = c("plain","bold","italic","bold.italic"))
                                               #)
                                             ),
                                             h4("Download Parameters"),
                                             fluidRow(
                                               column(4,
                                                      numericInput("VolMAdnldHeight","Plot Download Height",
                                                                   value = 8, step = 1, min = 1)
                                               ),
                                               column(4,
                                                      numericInput("VolMAdnldWidth","Plot Download Width",
                                                                   value = 10, step = 1, min = 1)
                                               ),
                                               column(4,
                                                      selectInput("VolMAdnldSizeUnits","Height/Width Units",
                                                                  choices = c("in","cm","mm","px"))
                                               )
                                             )
                            ),
                            
                            ####----Boxplot side panel----####
                            
                            conditionalPanel(condition = "input.datasettwo == '3'",
                                             p(),
                                             uiOutput("rendBoxPlotMetaColSelec"),
                                             h4("Gene Selection:"),
                                             div(DT::dataTableOutput("GeneListTable"), style = "font-size:10px"),
                                             p(),
                                             fluidRow(
                                               column(4,
                                                      numericInput("boxplot2TitleSize","Title Font Size",
                                                                   value = 20, step = 1)
                                               ),
                                               column(4,
                                                      numericInput("boxplot2AxisSize","Axis Lable Font Size",
                                                                   value = 16, step = 1)
                                               ),
                                               column(4,
                                                      numericInput("boxplotDot", "Dot Size:",
                                                                   min = 0, max = 5,
                                                                   value = 0.5, step = .25)
                                               )
                                             )
                                             
                            ),
                            
                            ####----Average Gene Expression Scatter Plot----####
                            
                            conditionalPanel(condition = "input.datasettwo == '4'",
                                             h4("Condition Selection"),
                                             uiOutput("rendAvgExprScatterMetaCol"),
                                             uiOutput("rendcomparisonA2.avg"),
                                             uiOutput("rendcomparisonB2.avg"),
                                             #selectInput("comparisonA2.avg", "Comparison: GroupA",
                                             #            choices = metagroups, selected = metagroups[1]),
                                             #selectInput("comparisonB2.avg", "Comparison: GroupB",
                                             #            choices = metagroups, selected = metagroups[2]),
                                             checkboxInput("AvgExpLogFC","Log2 Transform Expression Data", value = TRUE),
                                             h4("Gene Annotation Selection"),
                                             selectInput("scatterGeneSelec", "Select Genes to Annotate in Plot",
                                                         choices = rownames(expr),
                                                         multiple = T,
                                                         selected = "-"),
                                             textInput("gsSelection2", "Text Input of Genes (space or tab delimited):", value = ""),
                                             h4("Figure Parameters"),
                                             fluidRow(
                                               column(6,
                                                      numericInput("AvgExprScatterTitleSize","Title Font Size",
                                                                   value = 18, step = 1)
                                               ),
                                               column(6,
                                                      numericInput("AvgExprScatterAxisSize","Axis Lable Font Size",
                                                                   value = 16, step = 1))
                                             ),
                                             uiOutput("hover_info3")
                            ),
                            
                            ###--EnrichR pathway side panel--##
                            #
                            #conditionalPanel(condition = "input.datasettwo == '5'",
                            #                 h4("Condition Selection:"),
                            #                 uiOutput("rendEnrichRPathMetaCol"),
                            #                 uiOutput("rendcomparisonA2.path"),
                            #                 uiOutput("rendcomparisonB2.path"),
                            #                 #selectInput("comparisonA2.path", "Comparison: GroupA",
                            #                 #            choices = metagroups, selected = metagroups[1]),
                            #                 #selectInput("comparisonB2.path", "Comparison: GroupB",
                            #                 #            choices = metagroups, selected = metagroups[2]),
                            #                 selectInput("SelectedPathway", "Select Pathway",
                            #                             choices = enrichRchoice),
                            #                 h4("Threshold Parameters:"),
                            #                 numericInput("pathpval", "P-Value Cutoff:",
                            #                              min = 0, value = 0.01),
                            #                 numericInput("pathFC", "Log2 Fold Change Cutoff:",
                            #                              min = 0, value = 0),
                            #                 h4("Figure Paramters"),
                            #                 fluidRow(
                            #                   column(6,
                            #                          numericInput("PathwayTitleSize","Title Font Size",
                            #                                       value = 18, step = 1)
                            #                   ),
                            #                   column(6,
                            #                          numericInput("PathwayAxisSize","Axis Lable Font Size",
                            #                                       value = 16, step = 1))
                            #                 )
                            #                 
                            #),
                            
                            ##--DEG table side panel--##
                            
                            conditionalPanel(condition = "input.datasettwo == '6'",
                                             h4("Condition Selection:"),
                                             uiOutput("rendDEGtableMetaCol"),
                                             uiOutput("rendcomparisonA2.DEG"),
                                             uiOutput("rendcomparisonB2.DEG"),
                                             #selectInput("comparisonA2.DEG", "Comparison: GroupA",
                                             #            choices = metagroups, selected = metagroups[1]),
                                             #selectInput("comparisonB2.DEG", "Comparison: GroupB",
                                             #            choices = metagroups, selected = metagroups[2]),
                                             h4("Download DEG Table as GMT File"),
                                             textInput("DEGfileName", "File Name for Download:",value = "DEGgeneSet"),
                                             numericInput("fc_cutoff2", "LogFC Threshold",
                                                          min = 0, max = 5, step = 0.1, value = 1),
                                             selectInput("UpDnChoice","Up-regulated or Down-regulated:",
                                                         choices = c("UpAndDown_Regulated","Up_Regulated","Down_Regulated")),
                                             numericInput("p_cutoff2", "Adj.P.Value Cutoff:",
                                                          min = 0, max = 10, step = 0.1, value = 0.05),
                                             numericInput("top_x2", "Number of Top Hits:", value = 100),
                                             downloadButton("DEGgmtDownload", "Download DEG .gmt")
                            )
                          ),
                          mainPanel(
                            tabsetPanel(
                              id = "datasettwo",
                              tabPanel("Volcano Plot",
                                       p(),
                                       verbatimTextOutput("VolGroupsText"),
                                       jqui_resizable(plotOutput('Volcano3', width = "800px", height = "550px",
                                                                 hover = hoverOpts("plot_hover", delay = 10, delayType = "debounce"))),
                                       fluidRow(
                                         downloadButton("dnldPlotSVG_vol","Download as SVG"),
                                         downloadButton("dnldPlotPDF_vol","Download as PDF")
                                       ),
                                       uiOutput("hover_info"),
                                       value = 1),
                              tabPanel("MA Plot",
                                       p(),
                                       verbatimTextOutput("MAGroupsText"),
                                       jqui_resizable(plotOutput('MAPlot1', width = "800px", height = "550px",
                                                                 hover = hoverOpts("plot_hover2", delay = 10, delayType = "debounce"))),
                                       fluidRow(
                                         downloadButton("dnldPlotSVG_MA","Download as SVG"),
                                         downloadButton("dnldPlotPDF_MA","Download as PDF")
                                       ),
                                       uiOutput("hover_info2"),
                                       value = 2),
                              tabPanel("Box Plots",
                                       p(),
                                       #withSpinner(jqui_resizable(plotOutput('boxplot1', width = "100%", height = "600px")), type = 6),
                                       jqui_resizable(plotOutput('boxplot1', width = "100%", height = "600px")),
                                       fluidRow(
                                         downloadButton("dnldPlotSVG_exprBox2","Download as SVG"),
                                         downloadButton("dnldPlotPDF_exprBox2","Download as PDF")
                                       ),
                                       value = 3),
                              tabPanel("Average Expression Scatter Plot",
                                       p(),
                                       jqui_resizable(plotOutput("AvggeneScatter2", height = "500px",
                                                                 hover = hoverOpts("plot_hover3", delay = 10, delayType = "debounce"))),
                                       #withSpinner(jqui_resizable(plotOutput("AvggeneScatter2", height = "500px")), type = 6),
                                       #withSpinner(jqui_resizable(plotOutput("AvggeneScatter2", height = "500px",
                                       #                                      hover = hoverOpts("plot_hover", delay = 10, delayType = "debounce"))), type = 6),
                                       fluidRow(
                                         downloadButton("dnldPlotSVG_AvgScatter","Download as SVG"),
                                         downloadButton("dnldPlotPDF_AvgScatter","Download as PDF")
                                       ),
                                       DT::dataTableOutput("AvggeneScatterTable"),
                                       downloadButton("AvggeneScatterDownload", "Download Non-log2 Transformed .tsv"),
                                       value = 4),
                              #tabPanel("Pathway Analysis",
                              #         p("If Pathway Enrichment Analysis shows errors please adjust the P.Value and Log FC parameters."),
                              #         uiOutput("UpRegPathLabel"),
                              #         verbatimTextOutput("upregpath_text"),
                              #         withSpinner(plotOutput('UpRegPathway1'), type = 6),
                              #         fluidRow(
                              #           downloadButton("dnldPlotSVG_upPath","Download as SVG"),
                              #           downloadButton("dnldPlotPDF_upPath","Download as PDF")
                              #         ),
                              #         div(DT::dataTableOutput("UpRegPathwayTable1"), style = "font-size:10px"),
                              #         downloadButton("UpRegPathDownload", "Download Table .tsv"),
                              #         downloadButton("UpRegPathDownloadgmt", "Download .gmt"),
                              #         uiOutput("DnRegPathLabel"),
                              #         verbatimTextOutput("downregpath_text"),
                              #         withSpinner(plotOutput('DnRegPathway1'), type = 6),
                              #         fluidRow(
                              #           downloadButton("dnldPlotSVG_dnPath","Download as SVG"),
                              #           downloadButton("dnldPlotPDF_dnPath","Download as PDF")
                              #         ),
                              #         div(DT::dataTableOutput("DnRegPathwayTable1"), style = "font-size:10px"),
                              #         downloadButton("DnRegPathDownload", "Download Table .tsv"),
                              #         downloadButton("DnRegPathDownloadgmt", "Download .gmt"),
                              #         p(),
                              #         value = 5),
                              tabPanel("DEG Table",
                                       p(),
                                       verbatimTextOutput("degtext"),
                                       div(DT::dataTableOutput("DEGtable1"), style = "font-size:12px"),
                                       downloadButton("DEGtableDownload", "Download DEG table .tsv"),
                                       value = 6)
                            )
                          )
                        )
                      )
             ),
             
             ####----GSEA Tab----####
             
             tabPanel("GSEA Analysis",
                      fluidPage(
                        title = "GSEA Analysis",
                        sidebarLayout(
                          sidebarPanel(
                            tabsetPanel(id = "GSEA",
                                        tabPanel("GSEA Parameters",
                                                 h4("Condition Selection:"),
                                                 uiOutput("rendGSEAmetaCol"),
                                                 uiOutput("rendcomparisonA"),
                                                 uiOutput("rendcomparisonB"),
                                                 #selectInput("comparisonA", "Comparison: GroupA",
                                                 #            choices = metagroups, selected = metagroups[1]),
                                                 #selectInput("comparisonB", "Comparison: GroupB",
                                                 #            choices = metagroups, selected = metagroups[2]),
                                                 h4("GSEA Threshold Parameter:"),
                                                 numericInput("userPval", "Pvalue Cutoff", value = 1.0, width = '100%'),
                                                 h4("ssGSEA Boxplot Parameters"),
                                                 selectInput("ssGSEAtype","Choose ssGSEA Method",
                                                             choices = c("ssgsea","gsva","zscore","plage")),
                                                 tabsetPanel(
                                                   id = "tables",
                                                   tabPanel("MSigDB Gene Sets",
                                                            div(DT::dataTableOutput("msigdbTable"), style = "font-size:10px; height:500px; overflow-X: scroll"),
                                                            value = 1),
                                                   tabPanel(paste(userGSlist_name),
                                                            div(DT::dataTableOutput("tab2table"), style = "font-size:10px; height:500px; overflow-X: scroll"),
                                                            value = 3),
                                                   tabPanel("Use your own gene set",
                                                            p(),
                                                            fluidRow(
                                                              column(9,
                                                                     uiOutput("user.gmt")
                                                              ),
                                                              column(3,
                                                                     checkboxInput("UserGSheaderCheck","Header",value = T)
                                                              )
                                                            ),
                                                            #uiOutput("user.gmt"),
                                                            uiOutput("user.RDataButton"),
                                                            uiOutput("RDataMessage"),
                                                            uiOutput("user.GStable"),
                                                            value = 5)
                                                 )
                                        ),
                                        tabPanel("ssGSEA Heatmap Parameters",
                                                 p(),
                                                 selectInput("ssGSEAtypeHeat","Choose ssGSEA Method",
                                                             choices = c("ssgsea","gsva","zscore","plage")),
                                                 h4("Gene Set Selection"),
                                                 uiOutput("rendssgseaHeatGS"),
                                                 h4("Sample Selection"),
                                                 selectInput("userheatsampSS", "",
                                                             choices = sampsames, multiple = T, selected = sampsames)),
                                        tabPanel("Figure Parameters",
                                                 p(),
                                                 h4("Heatmap Parameters"),
                                                 selectInput("ColorPalette_gseaHeat", "Select Color Palette:",
                                                             choices = c("Red/Blue" = "original",
                                                                         "OmniBlueRed" = "OmniBlueRed",
                                                                         "LightBlue/BlackRed" = "LightBlueBlackRed",
                                                                         "Green/Black/Red" = "GreenBlackRed",
                                                                         "Yellow/Green/Blue" = "YlGnBu","Inferno" = "Inferno",
                                                                         "Viridis" = "Viridis","Plasma" = "Plasma",
                                                                         "Reds" = "OrRd","Blues" = "PuBu","Greens" = "Greens")),
                                                 fluidRow(
                                                   column(6,
                                                          checkboxInput("ShowColNamesSSheat","Show Heatmap Column Names", value = T),
                                                          numericInput("heatmapFont1.c", "Heatmap Column Font Size:",
                                                                       min = 5, max = 75,
                                                                       value = 12, step = 1),
                                                          checkboxInput("clustcolsSSheat","Cluster Heatmap Columns", value = F)
                                                   ),
                                                   column(6,
                                                          checkboxInput("ShowRowNames1SSheat","Show Heatmap Row Names", value = T),
                                                          numericInput("heatmapFont1.r", "Heatmap Row Font Size:",
                                                                       min = 5, max = 75,
                                                                       value = 10, step = 1),
                                                          checkboxInput("clustrowsSSheat","Cluster Heatmap Rows", value = F)
                                                   )
                                                 ),
                                                 uiOutput("rendClusterMethodSSheat"),
                                                 h4("Box Plot Parameters"),
                                                 selectInput("boxplotcompare", "Boxplot Stat Compare Method:",
                                                             choices = c("none","wilcox.test","t.test","kruskal.test","anova")),
                                                 fluidRow(
                                                   column(4,
                                                          numericInput("gseaBoxTitleSize","Title Font Size",
                                                                       value = 20, step = 1)
                                                   ),
                                                   column(4,
                                                          numericInput("gseaBoxAxisSize","Axis Lable Font Size",
                                                                       value = 16, step = 1)
                                                   ),
                                                   column(4,
                                                          numericInput("boxplotDotss", "Dot Size:",
                                                                       min = 0, max = 5,
                                                                       value = 0.75, step = .25)
                                                   )
                                                 )
                                        )
                            )
                          ),
                          mainPanel(
                            tabsetPanel(
                              id = "datasetthree",
                              tabPanel("Enrichment Plot",
                                       h3("GSEA Enrichment Plot"),
                                       verbatimTextOutput("NESandPval"),
                                       withSpinner(plotOutput("enrichplot0", width = "500px", height = "450px"), type = 6),
                                       fluidRow(
                                         downloadButton("dnldPlotSVG_gsea","Download as SVG"),
                                         downloadButton("dnldPlotPDF_gsea","Download as PDF")
                                       ),
                                       h3("Leading Edge Genes (~Signal2Noise Ranking)"),
                                       downloadButton("LEGdownload", "Download .tsv"),
                                       div(DT::dataTableOutput("LeadingEdgeGenes"), style = "font-size:12px; height:500px; overflow-y: scroll"),
                                       value = 1),
                              tabPanel("GSEA Heatmap",
                                       withSpinner(jqui_resizable(plotOutput("heatmap0", width = "100%", height = "1000px")), type = 6),
                                       fluidRow(
                                         downloadButton("dnldPlotSVG_gseaHeat","Download as SVG"),
                                         downloadButton("dnldPlotPDF_gseaHeat","Download as PDF")
                                       ),
                                       value = 2),
                              tabPanel("Generate Enriched Signatures Table",
                                       p("Please note this may take several minutes depending on size and quantity of gene sets in GMT file."),
                                       uiOutput("genESTbutton"),
                                       p(),
                                       withSpinner(DT::dataTableOutput("enrich_sig_table_gen"), type = 6),
                                       downloadButton("enrich_sig_download.u","Download .tsv"),
                                       value = 4),
                              tabPanel("ssGSEA Boxplots",
                                       p(),
                                       withSpinner(jqui_resizable(plotOutput('boxplot2', width = "100%", height = "500px")), type = 6),
                                       fluidRow(
                                         downloadButton("dnldPlotSVG_gseaBox","Download as SVG"),
                                         downloadButton("dnldPlotPDF_gseaBox","Download as PDF")
                                       ),
                                       DT::dataTableOutput("ssGSEAtable"),
                                       downloadButton("ssGSEAdownload", "Download .tsv"),
                                       value = 5),
                              tabPanel("ssGSEA Heatmap",
                                       tabsetPanel(
                                         tabPanel("ssGSEA zScore Heatmap",
                                                  p(),
                                                  withSpinner(jqui_resizable(plotOutput('ssgseaheatmap', width = "100%", height = "800px")), type = 6),
                                                  fluidRow(
                                                    downloadButton("dnldPlotSVG_ssgseaHeat","Download as SVG"),
                                                    downloadButton("dnldPlotPDF_ssgseaHeat","Download as PDF")
                                                  )
                                         ),
                                         tabPanel("ssGSEA Raw Difference Heatmap",
                                                  p(),
                                                  withSpinner(jqui_resizable(plotOutput('ssgseaheatmap2', width = "100%", height = "800px")), type = 6),
                                                  fluidRow(
                                                    downloadButton("dnldPlotSVG_ssgseaHeat2","Download as SVG"),
                                                    downloadButton("dnldPlotPDF_ssgseaHeat2","Download as PDF")
                                                  )
                                         ),
                                         tabPanel("Average Raw Difference ssGSEA Heatmap",
                                                  p(),
                                                  withSpinner(jqui_resizable(plotOutput('ssgseaheatmap3', width = "100%", height = "800px")), type = 6),
                                                  fluidRow(
                                                    downloadButton("dnldPlotSVG_ssgseaHeat3","Download as SVG"),
                                                    downloadButton("dnldPlotPDF_ssgseaHeat3","Download as PDF")
                                                  )
                                         )
                                       ),
                                       value = 6)
                            )
                          )
                        )
                      )
             )
  )


####----Server----####


server <- function(input, output, session) {
  
  
  ####----Render UI----####
  
  output$rendssgseaHeatGS <- renderUI({
    
    if (input$tables == 1) {
      gs_names <- c(grep("^HALLMARK",names(gs),value = T),grep("^HALLMARK",names(gs),value = T, invert = T))
    }
    if (input$tables == 3) {
      gs_names <- names(gs2)
    }
    if (input$tables == 5) {
      gs_names <- names(RDataListGen())
    }
    # take NA out if number of gs is less than 50
    gs_selected <- gs_names[1:3]
    gs_selected <- gs_selected[!is.na(gs_selected)]
    selectInput("ssgseaHeatGS", "Select Gene Sets:",
                choices = gs_names, selected = gs_selected, multiple = T)
    
  })
  
  # Render cluster method selection for average heatmap
  output$ClusterMethodMVG <- renderUI({
    
    if (input$clustrowsMVG == TRUE | input$clustcolsMVG == TRUE) {
      
      selectInput("ClusteringMethod",
                  "Select Clustering Method",
                  choices = c("complete", "ward.D", "ward.D2", "single", "average", "mcquitty", "median", "centroid"))
      
    }
    
  })
  
  output$rendMetaColDEGHeat <- renderUI({
    
    if (ncol(meta) > 2){
      selectInput("MetaColDEGHeat","Meta Column:",
                  choices = colnames(meta)[2:ncol(meta)],
                  selected = Feature_Selected)
    }
    
  })
  
  output$rendMVGheatAnnoCol <- renderUI({
    
    if (ncol(meta) > 2) {
      selectInput("MVGheatAnnoCol","Select Annotation Column:",
                  choices = colnames(meta[,2:ncol(meta)]),
                  selected = Feature_Selected)
    }
    
  })
  
  output$rendCustomheatAnnoCol <- renderUI({
    
    if (ncol(meta) > 2) {
      selectInput("CustomheatAnnoCol","Select Annotation Column:",
                  choices = colnames(meta[,2:ncol(meta)]),
                  selected = Feature_Selected)
    }
    
  })
  
  output$rendBoxPlot1MetaCol <- renderUI({
    
    if (ncol(meta) > 2) {
      selectInput("BoxPlot1MetaCol","Select Annotation Column:",
                  choices = colnames(meta[,2:ncol(meta)]),
                  selected = Feature_Selected)
    }
    
  })
  
  output$rendBarPlot1MetaCol <- renderUI({
    
    if (ncol(meta) > 2) {
      selectInput("BarPlot1MetaCol","Select Annotation Column:",
                  choices = colnames(meta[,2:ncol(meta)]),
                  selected = Feature_Selected)
    }
    
  })
  
  output$rendAvgHeatMetaCol <- renderUI({
    
    if (ncol(meta) > 2) {
      selectInput("AvgHeatMetaCol","Select Meta Column:",
                  choices = colnames(meta[,2:ncol(meta)]),
                  selected = Feature_Selected)
    }
    
  })
  
  output$rendSampCondSelectionCust <- renderUI({
    
    if (ncol(meta) > 2) {
      metagroups_new <- as.vector(levels(factor(meta[,input$AvgHeatMetaCol])))
      selectInput("SampCondSelectionCust","Sample Coniditon Selection",
                  choices = metagroups_new, selected = metagroups_new[c(1:length(metagroups_new))],
                  multiple = T)
    }
    else if (ncol(meta) == 2) {
      selectInput("SampCondSelectionCust","Sample Coniditon Selection",
                  choices = metagroups, selected = metagroups[c(1:length(metagroups))],
                  multiple = T)
    }
    
  })
  
  output$rendcomparisonA2_h <- renderUI({
    
    if (ncol(meta) == 2){
      selectInput("comparisonA2_h", "Comparison: GroupA",
                  choices = metagroups, selected = metagroups[1])
    }
    else if (ncol(meta) > 2){
      metagroups_new <- as.vector(levels(factor(meta[,input$MetaColDEGHeat])))
      selectInput("comparisonA2_h", "Comparison: GroupA",
                  choices = metagroups_new, selected = metagroups_new[1])
    }
    
  })
  
  output$rendGSEAmetaCol <- renderUI({
    
    if (ncol(meta) > 2) {
      selectInput("GSEAmetaCol","Select Meta Column:",
                  choices = colnames(meta[,2:ncol(meta)]),
                  selected = Feature_Selected)
    }
    
  })
  
  output$rendcomparisonA <- renderUI({
    
    if (ncol(meta) == 2){
      selectInput("comparisonA", "Comparison: GroupA",
                  choices = metagroups, selected = metagroups[1])
    }
    else if (ncol(meta) > 2){
      metagroups_new <- as.vector(levels(factor(meta[,input$GSEAmetaCol])))
      selectInput("comparisonA", "Comparison: GroupA",
                  choices = metagroups_new, selected = metagroups_new[1])
    }
    
  })
  
  output$rendcomparisonB <- renderUI({
    
    if (ncol(meta) == 2){
      selectInput("comparisonB", "Comparison: GroupB",
                  choices = metagroups, selected = metagroups[2])
    }
    else if (ncol(meta) > 2){
      metagroups_new <- as.vector(levels(factor(meta[,input$GSEAmetaCol])))
      selectInput("comparisonB", "Comparison: GroupB",
                  choices = metagroups_new, selected = metagroups_new[2])
    }
    
  })
  
  output$rendcomparisonB2_h <- renderUI({
    
    if (ncol(meta) == 2){
      selectInput("comparisonB2_h", "Comparison: GroupB",
                  choices = metagroups, selected = metagroups[2])
    }
    else if (ncol(meta) > 2){
      metagroups_new <- as.vector(levels(factor(meta[,input$MetaColDEGHeat])))
      selectInput("comparisonB2_h", "Comparison: GroupB",
                  choices = metagroups_new, selected = metagroups_new[2])
    }
    
  })
  
  output$rendDEGMetaCol <- renderUI({
    
    if (ncol(meta) > 2) {
      selectInput("DEGMetaCol","Select Meta Column:",
                  choices = colnames(meta[,2:ncol(meta)]),
                  selected = Feature_Selected)
    }
    
  })
  
  output$rendcomparisonA2 <- renderUI({
    
    if (ncol(meta) == 2){
      selectInput("comparisonA2", "Comparison: GroupA",
                  choices = metagroups, selected = metagroups[1])
    }
    else if (ncol(meta) > 2){
      metagroups_new <- as.vector(levels(factor(meta[,input$DEGMetaCol])))
      selectInput("comparisonA2", "Comparison: GroupA",
                  choices = metagroups_new, selected = metagroups_new[1])
    }
    
  })
  
  output$rendcomparisonB2 <- renderUI({
    
    if (ncol(meta) == 2){
      selectInput("comparisonB2", "Comparison: GroupB",
                  choices = metagroups, selected = metagroups[2])
    }
    else if (ncol(meta) > 2){
      metagroups_new <- as.vector(levels(factor(meta[,input$DEGMetaCol])))
      selectInput("comparisonB2", "Comparison: GroupB",
                  choices = metagroups_new, selected = metagroups_new[2])
    }
    
  })
  
  output$rendAvgExprScatterMetaCol <- renderUI({
    
    if (ncol(meta) > 2) {
      selectInput("AvgExprScatterMetaCol","Select Meta Column:",
                  choices = colnames(meta[,2:ncol(meta)]),
                  selected = Feature_Selected)
    }
    
  })
  
  output$rendcomparisonA2.avg <- renderUI({
    
    if (ncol(meta) == 2){
      selectInput("comparisonA2.avg", "Comparison: GroupA",
                  choices = metagroups, selected = metagroups[1])
    }
    else if (ncol(meta) > 2){
      metagroups_new <- as.vector(levels(factor(meta[,input$AvgExprScatterMetaCol])))
      selectInput("comparisonA2.avg", "Comparison: GroupA",
                  choices = metagroups_new, selected = metagroups_new[1])
    }
    
  })
  
  output$rendcomparisonB2.avg <- renderUI({
    
    if (ncol(meta) == 2){
      selectInput("comparisonB2.avg", "Comparison: GroupB",
                  choices = metagroups, selected = metagroups[2])
    }
    else if (ncol(meta) > 2){
      metagroups_new <- as.vector(levels(factor(meta[,input$AvgExprScatterMetaCol])))
      selectInput("comparisonB2.avg", "Comparison: GroupB",
                  choices = metagroups_new, selected = metagroups_new[2])
    }
    
  })
  
  output$rendScatterPlotMetaCol <- renderUI({
    
    if (ncol(meta) > 2){
      selectInput("ScatterPlotMetaCol","Select Sample Condition",
                  choices = colnames(meta)[2:ncol(meta)],
                  selected = Feature_Selected)
    }
    
  })
  
  output$rendEnrichRPathMetaCol <- renderUI({
    
    if (ncol(meta) > 2){
      selectInput("EnrichRPathMetaCol","Select Sample Condition",
                  choices = colnames(meta)[2:ncol(meta)],
                  selected = Feature_Selected)
    }
    
  })
  
  output$rendcomparisonA2.path <- renderUI({
    
    if (ncol(meta) == 2){
      selectInput("comparisonA2.path", "Comparison: GroupA",
                  choices = metagroups, selected = metagroups[1])
    }
    else if (ncol(meta) > 2){
      metagroups_new <- as.vector(levels(factor(meta[,input$EnrichRPathMetaCol])))
      selectInput("comparisonA2.path", "Comparison: GroupA",
                  choices = metagroups_new, selected = metagroups_new[1])
    }
    
  })
  
  output$rendcomparisonB2.path <- renderUI({
    
    if (ncol(meta) == 2){
      selectInput("comparisonB2.path", "Comparison: GroupB",
                  choices = metagroups, selected = metagroups[2])
    }
    else if (ncol(meta) > 2){
      metagroups_new <- as.vector(levels(factor(meta[,input$EnrichRPathMetaCol])))
      selectInput("comparisonB2.path", "Comparison: GroupB",
                  choices = metagroups_new, selected = metagroups_new[2])
    }
    
  })
  
  output$rendDEGtableMetaCol <- renderUI({
    
    if (ncol(meta) > 2){
      selectInput("DEGtableMetaCol","Select Sample Condition",
                  choices = colnames(meta)[2:ncol(meta)],
                  selected = Feature_Selected)
    }
    
  })
  
  output$rendcomparisonA2.DEG <- renderUI({
    
    if (ncol(meta) == 2){
      selectInput("comparisonA2.DEG", "Comparison: GroupA",
                  choices = metagroups, selected = metagroups[1])
    }
    else if (ncol(meta) > 2){
      metagroups_new <- as.vector(levels(factor(meta[,input$DEGtableMetaCol])))
      selectInput("comparisonA2.DEG", "Comparison: GroupA",
                  choices = metagroups_new, selected = metagroups_new[1])
    }
    
  })
  
  output$rendcomparisonB2.DEG <- renderUI({
    
    if (ncol(meta) == 2){
      selectInput("comparisonB2.DEG", "Comparison: GroupB",
                  choices = metagroups, selected = metagroups[2])
    }
    else if (ncol(meta) > 2){
      metagroups_new <- as.vector(levels(factor(meta[,input$DEGtableMetaCol])))
      selectInput("comparisonB2.DEG", "Comparison: GroupB",
                  choices = metagroups_new, selected = metagroups_new[2])
    }
    
  })
  
  output$rendBoxPlotMetaColSelec <- renderUI({
    
    if (ncol(meta) > 2){
      selectInput("BoxPlotMetaColSelec","Select Sample Condition",
                  choices = colnames(meta)[2:ncol(meta)],
                  selected = Feature_Selected)
    }
    
  })
  
  # Render cluster method selection for average heatmap
  output$rendClustMethodsCust <- renderUI({
    
    if (input$ClusRowOptCust == TRUE | input$ClusColOptCust == TRUE) {
      
      selectInput("ClusteringMethodCust",
                  "Select Clustering Method",
                  choices = c("complete", "ward.D", "ward.D2", "single", "average", "mcquitty", "median", "centroid"))
      
    }
    
  })
  
  # Render cluster method selection for average heatmap
  output$rendClustMethodsDEG <- renderUI({
    
    if (input$ClusRowOptDEG == TRUE | input$ClusColOptDEG == TRUE) {
      
      selectInput("ClusteringMethod_degh",
                  "Select Clustering Method",
                  choices = c("complete", "ward.D", "ward.D2", "single", "average", "mcquitty", "median", "centroid"))
      
    }
    
  })
  
  output$rendAvgClustMethodsCust <- renderUI({
    
    if (input$AvgClusRowOptCust == TRUE | input$AvgClusColOptCust == TRUE) {
      
      selectInput("AvgClusteringMethodCust",
                  "Select Clustering Method",
                  choices = c("complete", "ward.D", "ward.D2", "single", "average", "mcquitty", "median", "centroid"))
      
    }
    
  })
  
  # Render cluster method selection for average heatmap
  output$rendClusterMethodSSheat <- renderUI({
    
    if (input$clustcolsSSheat == TRUE | input$clustrowsSSheat == TRUE) {
      
      selectInput("ClusterMethodSSheat",
                  "Select Clustering Method",
                  choices = c("complete", "ward.D", "ward.D2", "single", "average", "mcquitty", "median", "centroid"))
      
    }
    
  })
  
  
  
  #render enrich signature table generation button for msigdb data
  output$genESTbutton <- renderUI({
    if (input$tables == 1) {
      actionButton("GenerateEST", "Generate Enriched Signature Table with MSigDB Gene Sets")
    }
  })
  
  #render user gmt data upload if indicated
  output$user.gmt <- renderUI({
    fileInput("user.gmt.file", "Gene Set File (.gmt, .tsv, or .txt)", accept = c(".gmt",".tsv",".txt"))
  })
  
  #render gene set table based off gmt file given
  output$user.GStable <- renderUI({
    req(input$user.gmt.file)
    div(DT::dataTableOutput("GStable.u"), style = "font-size:10px; height:500px; overflow-X: scroll")
  })
  
  #Warning message for RData list generation
  output$RDataMessage <- renderUI({
    req(input$user.gmt.file)
    p("Generating RData list may take up to several minutes depending on size of GMT file.")
  })
  
  #render action button to create RData list for ssGSEA boxplots
  output$user.RDataButton <- renderUI({
    req(input$user.gmt.file)
    actionButton("user.RData.Gen", "Generate RData list for ssGSEA Boxplot and Heatmap")
  })
  
  #render UI for hover text in volcano plot
  output$hover_info <- renderUI({
    top2 <- topgenereact()
    df <- top2 %>%
      select(GeneName,logFC,P.Value,adj.P.Val)
    colnames(df)[3] <- "-log10(P.Value)"
    df$P.Value <- df$'-log10(P.Value)'
    df$`-log10(P.Value)` <- -log10(df$`-log10(P.Value)`)
    hover <- input$plot_hover
    point <- nearPoints(df, hover, threshold = 10, maxpoints = 1, addDist = FALSE)
    if (nrow(point) == 0) return(NULL)
    wellPanel(
      p(HTML(paste0("<b> Name: </b>", point[1], "<br/>",
                    "<b> Fold change: </b>", round(point[2], digits = 4), "<br/>",
                    "<b> P Value: </b>", point[5], "<br/>",
                    "<b> adj P Value: </b>", point[4], "<br/>",
                    NULL
      ))))
  })
  
  #render UI for hover text in MA plot
  output$hover_info2 <- renderUI({
    top2 <- topgenereact()
    df <- top2 %>%
      select(GeneName,AveExpr,logFC,P.Value,adj.P.Val)
    colnames(df)[4] <- "-log10(P.Value)"
    df$P.Value <- df$'-log10(P.Value)'
    df$`-log10(P.Value)` <- -log10(df$`-log10(P.Value)`)
    hover <- input$plot_hover2
    point <- nearPoints(df, hover, threshold = 10, maxpoints = 1, addDist = FALSE)
    if (nrow(point) == 0) return(NULL)
    wellPanel(
      p(HTML(paste0("<b> Name: </b>", point[1], "<br/>",
                    "<b> Fold change: </b>", round(point[3], digits = 4), "<br/>",
                    "<b> P Value: </b>", point[6], "<br/>",
                    "<b> adj P Value: </b>", point[5], "<br/>",
                    "<b> Avg. Expression: </b>", round(point[2], digits = 4), "<br/>",
                    NULL
      ))))
  })
  
  AvgExprReact <- reactive({
    
    A_choice <- input$comparisonA2.avg
    B_choice <- input$comparisonB2.avg
    log_choice <- input$AvgExpLogFC
    if (ncol(meta) > 2) {
      metacol <- input$AvgExprScatterMetaCol
    }
    else if (ncol(meta) == 2) {
      metacol <- colnames(meta)[2]
    }
    # Get A and B sample names
    A <- meta[which(meta[,metacol] == A_choice),1]
    B <- meta[which(meta[,metacol] == B_choice),1]
    # Make A and B expression Matrices
    mat_A <- expr[,A]
    mat_B <- expr[,B]
    logFC_text <- ""
    # Log if user designates
    if (log_choice == TRUE) {
      mat_A <- log2(mat_A + 1)
      mat_B <- log2(mat_B + 1)
      logFC_text <- " (log2 +1)"
    }
    # Get avg expression of each gene
    mat_A$AvgExpression_GroupA <- rowMeans(mat_A)
    mat_B$AvgExpression_GroupB <- rowMeans(mat_B)
    # Add Gene Names as column
    mat_A$GeneSymbol <- rownames(mat_A)
    mat_B$GeneSymbol <- rownames(mat_B)
    # Merge average columns
    AvgExpr_Table <- merge(mat_A[,c("AvgExpression_GroupA","GeneSymbol")],
                           mat_B[,c("AvgExpression_GroupB","GeneSymbol")],
                           by = "GeneSymbol",
                           all = T)
    AvgExpr_Table
    
  })
  
  #render UI for hover text in volcano plot
  output$hover_info3 <- renderUI({
    
    #top2 <- topgenereact()
    #df2 <- top2 %>%
    #  select(GeneName,AveExpr,logFC,P.Value,adj.P.Val)
    #colnames(df2)[4] <- "-log10(P.Value)"
    #df2$P.Value <- df2$'-log10(P.Value)'
    #df2$`-log10(P.Value)` <- -log10(df2$`-log10(P.Value)`)
    
    A_choice <- input$comparisonA2.avg
    B_choice <- input$comparisonB2.avg
    df <- AvgExprReact()
    hover <- input$plot_hover3
    point <- nearPoints(df, hover, threshold = 10, maxpoints = 1, addDist = FALSE)
    if (nrow(point) == 0) return(NULL)
    wellPanel(
      p(HTML(paste0("<b> Gene: </b>", point[1], "<br/>",
                    "<b> ",A_choice," Average Expression: </b>", round(point[2], digits = 4), "<br/>",
                    "<b> ",B_choice," Average Expression: </b>", round(point[3], digits = 4), "<br/>",
                    NULL
      ))))
  })
  
  ####----Reactives----####
  
  
  ssGSEA_Heat_GS <- reactive({
    
    if (input$tables == 1) {
      ## starting ssgsea heatmap genes
      gs_names_start <- c(grep("^HALLMARK",names(gs),value = T),grep("^HALLMARK",names(gs),value = T, invert = T))
      gs_names_start <- gs_names_start[1:3]
    }
    else if (input$tables == 3) {
      gs_names_start <- names(gs2)[1:3]
    }
    else if (input$tables == 5) {
      gs_names_start <- names(RDataListGen())[1:3]
    }
    gs_names_start
    
  })
  
  GeneratedMSigDBEST <- reactive({
    if (input$GenerateEST == TRUE) {
      if (ncol(meta) > 2) {
        metacol <- input$GSEAmetaCol
      }
      else if (ncol(meta) == 2) {
        metacol <- colnames(meta)[2]
      }
      groupA <- meta[which(meta[,metacol] == input$comparisonA),1]
      groupB <- meta[which(meta[,metacol] == input$comparisonB),1]
      ##----Signal-to-Noise Calculation----##
      A <- A + 0.00000001
      P = as.matrix(as.numeric(colnames(A) %in% groupA))
      n1 <- sum(P[,1])
      M1 <- A %*% P
      M1 <- M1/n1
      A2 <- A*A
      S1 <- A2 %*% P
      S1 <- S1/n1 - M1*M1 
      S1 <- sqrt(abs((n1/(n1-1)) * S1))
      P = as.matrix(as.numeric(colnames(A) %in% groupB))
      n2 <- sum(P[,1])
      M2 <- A %*% P
      M2 <- M2/n2
      A2 <- A*A
      S2 <- A2 %*% P
      S2 <- S2/n2 - M2*M2
      S2 <- sqrt(abs((n2/(n2-1)) * S2))
      rm(A2)
      # small sigma "fix" as used in GeneCluster
      S2 <- ifelse(0.2*abs(M2) < S2, S2, 0.2*abs(M2))
      S2 <- ifelse(S2 == 0, 0.2, S2)
      S1 <- ifelse(0.2*abs(M1) < S1, S1, 0.2*abs(M1))
      S1 <- ifelse(S1 == 0, 0.2, S1)
      M1 <- M1 - M2
      rm(M2)
      S1 <- S1 + S2
      rm(S2)
      s2n.matrix <- M1/S1
      ##----Reformatting----##
      s2n.df <- as.data.frame(s2n.matrix)
      s2n.df$GeneID <- rownames(s2n.df)
      rownames(s2n.df) <- NULL
      data <- dplyr::select(s2n.df, GeneID, V1)
      data.gsea <- data$V1
      names(data.gsea) <- as.character(data$GeneID)
      s2n.matrix.s <- sort(data.gsea, decreasing = T)
      ####----GSEA----####
      #perform GSEA
      gsea.res <- GSEA(s2n.matrix.s, TERM2GENE = gmt, verbose = F, pvalueCutoff = 1)
      #extract results and convert to tibble
      gsea.df <- gsea.res@result
      gsea.df
    }
  })
  
  #reactive to generate RData list
  RDataListGen <- eventReactive(input$user.RData.Gen, {
    gmt <- GStable.ubg()
    colnames(gmt) <- c("term","gene")
    RData.u <- list()
    for (i in unique(gmt[,1])){
      RData.u[[i]] <- gmt[gmt[,1] == i,]$gene
    }
    RData.u
  })
  
  #reactive for ssGSEA function
  ssGSEAfunc <- reactive({
    if (input$tables == 1) {
      GS <- gs[(msigdb.gsea2[input$msigdbTable_rows_selected,3])]
    }
    if (input$tables == 3) {
      GS <- gs2[(GeneSet2[input$tab2table_rows_selected,1])]
    }
    if (input$tables == 5) {
      GS <- RDataListGen()[(user_gs_mirror()[input$GStable.u_rows_selected,1])]
    }
    GSVA::gsva(A, GS, method = input$ssGSEAtype, verbose = F)
  })
  
  #create background GMT from user input gene set table
  GStable.ubg <- reactive({
    gmt.u <- input$user.gmt.file
    ext <- tools::file_ext(gmt.u$datapath)
    req(gmt.u)
    headerCheck <- input$UserGSheaderCheck
    validate(need(ext == c("gmt","tsv","txt"), "Please upload .gmt, .tsv, or .txt file"))
    if (ext == "gmt") {
      read.gmt(gmt.u$datapath)
    }
    else {
      as.data.frame(read_delim(gmt.u$datapath, delim = '\t', col_names = headerCheck))
    }
  })
  
  #gs mirror from user input for selection help on back end
  user_gs_mirror <- reactive({
    GeneSet <- as.data.frame(unique(GStable.ubg()[1]))
    rownames(GeneSet) <- 1:nrow(GeneSet)
    colnames(GeneSet)[1] <- "Gene_Set"
    GeneSet
  })
  
  #perform sig2noise calculation and create GSEA result from user chosen gene set
  datasetInput <- reactive({
    if (ncol(meta) > 2) {
      req(input$GSEAmetaCol)
      metacol <- input$GSEAmetaCol
    }
    else if (ncol(meta) == 2) {
      metacol <- colnames(meta)[2]
    }
    groupA <- meta[which(meta[,metacol] == input$comparisonA),1]
    groupB <- meta[which(meta[,metacol] == input$comparisonB),1]
    ##----Signal-to-Noise Calculation----##
    A <- A + 0.00000001
    P = as.matrix(as.numeric(colnames(A) %in% groupA))
    n1 <- sum(P[,1])
    M1 <- A %*% P
    M1 <- M1/n1
    A2 <- A*A
    S1 <- A2 %*% P
    S1 <- S1/n1 - M1*M1 #
    S1 <- sqrt(abs((n1/(n1-1)) * S1))
    P = as.matrix(as.numeric(colnames(A) %in% groupB))
    n2 <- sum(P[,1])
    M2 <- A %*% P
    M2 <- M2/n2
    A2 <- A*A
    S2 <- A2 %*% P
    S2 <- S2/n2 - M2*M2
    S2 <- sqrt(abs((n2/(n2-1)) * S2))
    rm(A2)
    # small sigma "fix" as used in GeneCluster
    S2 <- ifelse(0.2*abs(M2) < S2, S2, 0.2*abs(M2))
    S2 <- ifelse(S2 == 0, 0.2, S2)
    S1 <- ifelse(0.2*abs(M1) < S1, S1, 0.2*abs(M1))
    S1 <- ifelse(S1 == 0, 0.2, S1)
    M1 <- M1 - M2
    rm(M2)
    S1 <- S1 + S2
    rm(S2)
    s2n.matrix <- M1/S1
    ##----Reformatting----##
    s2n.df <- as.data.frame(s2n.matrix)
    s2n.df$GeneID <- rownames(s2n.df)
    rownames(s2n.df) <- NULL
    data <- dplyr::select(s2n.df, GeneID, V1)
    data.gsea <- data$V1
    names(data.gsea) <- as.character(data$GeneID)
    s2n.matrix.s <- sort(data.gsea, decreasing = T)
    ##----GSEA----##
    if (input$tables == 1){
      GSEA(s2n.matrix.s, TERM2GENE = gmt[which(gmt[,1] == as.character(msigdb.gsea2[input$msigdbTable_rows_selected,3])),],
           verbose = F, pvalueCutoff = input$userPval)
    }
    else if (input$tables == 3){
      GSEA(s2n.matrix.s, TERM2GENE = tab2[which(tab2[,1] == as.character(GeneSet2[input$tab2table_rows_selected,1])),],
           verbose = F, pvalueCutoff = input$userPval)
    }
    else if (input$tables == 5){
      GSEA(s2n.matrix.s, TERM2GENE = GStable.ubg()[which(GStable.ubg()[,1] == as.character(user_gs_mirror()[input$GStable.u_rows_selected,1])),],
           verbose = F, pvalueCutoff = input$userPval)
    }
  })
  
  
  #gmt_sub <- GStable.ubg[which(GStable.ubg[,1] == as.character(GeneSet[1,1])),]
  
  #top genes data frame reactive
  topgenereact <- reactive({
    #make group based on user input
    if (ncol(meta) > 2) {
      metacol <- input$DEGMetaCol
    }
    else if (ncol(meta) == 2) {
      metacol <- colnames(meta)[2]
    }
    A <- meta[which(meta[,metacol] == input$comparisonA2),1]
    B <- meta[which(meta[,metacol] == input$comparisonB2),1]
    #make top table
    #samples <- c(A,B)
    #mat <- expr[,which(colnames(expr) %in% samples)]
    mat <- expr[,c(A,B)]
    mat <- log2(mat + 1.0)
    groupAOther <- factor(c(rep("A", length(A)), rep("B", length(B))))
    designA <- model.matrix(~0 + groupAOther)
    fit <- lmFit(mat, design = designA)
    contrast.matrix <- makeContrasts(groupAOtherA - groupAOtherB, levels = designA)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    options(digits = 4)
    top1 <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
    top2 <- top1
    
    #top2 = read.delim("DEseq_BCL11BKO_CD8_Peritoneal_Tomas_analysis_20230406_max_v2.txt",
    #                  sep = '\t', header = T, strip.white = T)
    #colnames(top2)[7] <- "adj.P.Val"
    #rownames(top2) = top2$GeneName
    top2["GeneName"] <- rownames(top2)
    Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
    top2['group'] <- "NotSignificant"
    p_cut <- input$p_cutoff
    f_cut <- input$fc_cutoff
    top2[which(top2$P.Value < p_cut & abs(top2$logFC) < abs(f_cut)), "group"] <- "Significant"
    top2[which(top2$P.Value > p_cut & abs(top2$logFC) > abs(f_cut)), "group"] <- "FoldChange"
    top2[which(top2$P.Value < p_cut & abs(top2$logFC) > abs(f_cut)), "group"] <- "Significant&FoldChange"
    top2['FCgroup'] <- "NotSignificant"
    top2[which(abs(top2$logFC) > abs(f_cut)), "group2"] <- "FoldChange"
    top2
  })
  
  topgenereact2 <- reactive({
    
    # UI Inputs
    A_choice <- input$comparisonA2_h            #Comparison group A
    B_choice <- input$comparisonB2_h            #Comparison group B
    
    # Make Top table - Compares only 2 groups
    #make group based on user input
    if (ncol(meta) > 2) {
      metacol <- input$MetaColDEGHeat
    }
    else if (ncol(meta) == 2) {
      metacol <- colnames(meta)[2]
    }
    A <- meta[which(meta[,metacol] == A_choice),1]
    B <- meta[which(meta[,metacol] == B_choice),1]
    mat <- expr[,c(A,B)]
    mat <- log2(mat + 1.0)
    groupAOther <- factor(c(rep("A", length(A)), rep("B", length(B))))
    designA <- model.matrix(~0 + groupAOther)
    fit <- lmFit(mat, design = designA)
    contrast.matrix <- makeContrasts(groupAOtherA - groupAOtherB, levels = designA)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    options(digits = 4)
    top1 <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
    top1
    
  })
  
  ####----Data Exploration----####
  
  ####----MVG Heatmap----####
  
  ## MVG heatmap reactive
  MVGheatmap_react <- reactive ({
    
    top_probes <- input$NumFeatures
    row_names_choice <- input$ShowRowNames1
    col_names_choice <- input$ShowColNames1
    row_font <- input$heatmapFont2.r
    col_font <- input$heatmapFont2.c
    clust_method <- input$ClusteringMethod
    color_choice <- input$ColorPalette1
    var_type <- input$VarianceMeasure
    clust_cols_opt <- input$clustcolsMVG
    clust_rows_opt <- input$clustrowsMVG
    if (ncol(meta) > 2) {
      metacol <- input$MVGheatAnnoCol
    }
    else if (ncol(meta) == 2) {
      metacol <- colnames(meta)[2]
    }
    
    
    exp <- expr
    mad <- NULL
    var <- NULL
    cv <- NULL
    if (var_type == "MAD"){
      mad <- apply(log2(exp + 1), 1, mad)
      mad <- sort(mad, decreasing = T)
      mad <- head(mad, n = (top_probes +1))
      out <- cbind(names(mad), mad[names(mad)], exp[names(mad),])
      colnames(out) <- c("Gene", "MAD", colnames(exp))
      dataset <- exp[names(mad),]
      variable_gene_list <- names(mad)
    }
    else if (var_type == "VAR"){
      var <- apply(log2(exp + 1), 1, var)
      var <- sort(var, decreasing = T)
      var <- head(var, n = (top_probes +1))
      out <- cbind(names(var), var[names(var)], exp[names(var),])
      colnames(out) <- c("Gene", "VAR", colnames(exp))
      dataset <- exp[names(var),]
      variable_gene_list <- names(var)
    }
    else if (var_type == "CV"){
      cv <- apply(log2(exp + 1), 1, cv)
      cv <- sort(cv, decreasing = T)
      cv <- head(cv, n = (top_probes +1))
      out <- cbind(names(cv), cv[names(cv)], exp[names(cv),])
      colnames(out) <- c("Gene", "CV", colnames(exp))
      dataset <- exp[names(cv),]
      variable_gene_list <- names(cv)
    }
    
    dataset <- log2(dataset + 1)
    zdataset <- apply(dataset, 1, scale)
    zdataset <- apply(zdataset, 1, rev)
    colnames(zdataset) <- names(dataset)
    dataset <- as.matrix(zdataset)
    dataset[is.na(dataset)] <- 0
    
    dataset = dataset[apply(dataset, 1, function(x) !all(x==0)),]
    minimum = -5;
    maximum = 5;
    if (abs(min(dataset)) > abs(max(dataset))) {
      dataset[dataset < -abs(max(dataset))] = -abs(max(dataset))
    } else {
      dataset[dataset > abs(min(dataset))] = abs(min(dataset))
    }
    #type <- meta[,input$MVGheatAnnoCol]
    #meta2 <- as.data.frame(as.factor(type))
    meta2 <- meta[,c(colnames(meta)[1],metacol)]
    meta2[,2] <- as.factor(meta2[,2])
    meta2 <- meta2[order(meta2[,2]),]
    rownames(meta2) <- meta2[,1]
    meta2 <- meta2[,-1,drop = F]
    dataset <- dataset[,rownames(meta2)]
    bk = c(seq(minimum,minimum/2, length=100), seq(minimum/2,maximum/2,length=100),seq(maximum/2,maximum,length=100))
    #Heatmap color
    col_sets <- c("OrRd","PuBu","Greens","YlGnBu")
    if (color_choice == "original") {
      HeatMap_Colors <- c("dark blue","blue","white","red", "dark red")
      hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
    }
    else if (color_choice %in% col_sets) {
      HeatMap_Colors <- brewer.pal(n = 5, color_choice)
      hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
    }
    else if (color_choice == "Inferno") {
      hmcols <- inferno(500)
    }
    else if (color_choice == "Viridis") {
      hmcols <- viridis(500)
    }
    else if (color_choice == "Plasma") {
      hmcols <- plasma(500)
    }
    else if (color_choice == "OmniBlueRed") {
      hmcols<- colorRampPalette(c("#1984C5", "#22A7F0", "#63BFF0", "#A7D5ED", "#E2E2E2", "#E1A692", "#DE6E56", "#E14B31", "#C23728"))(length(bk)-1)
    }
    else if (color_choice == "LightBlueBlackRed") {
      hmcols<- colorRampPalette(c("#34C5FD","black","red"))(length(bk)-1)
    }
    else if (color_choice == "GreenBlackRed") {
      hmcols<- colorRampPalette(c("green","black","red"))(length(bk)-1)
    }
    hm <- pheatmap::pheatmap(dataset,
                             cluster_col = clust_cols_opt,
                             cluster_row = clust_rows_opt,
                             fontsize_row = row_font,
                             fontsize_col = col_font,
                             show_rownames = row_names_choice,
                             show_colnames = col_names_choice,
                             annotation_col = meta2,
                             clustering_method = clust_method,
                             color = hmcols,
                             border_color = NA)
    
    hm
    
  })
  
  ####----Custom Heatmap----####
  
  #Custom heat map reactive
  CustomHeatmap_react <- reactive({
    
    if (length(input$heatmapGeneSelec) >= 2 || length(input$userheatgenes) >= 1) {
      row_names_choice <- input$ShowRowNames2
      col_names_choice <- input$ShowColNames2
      clust_col_choice <- input$ClusColOptCust
      genelist.uih <- NULL
      genelist.ush <- NULL
      genelist.uih2 <- NULL
      genelist.ush <- input$heatmapGeneSelec
      genelist.uih <- unlist(strsplit(input$userheatgenes, " "))
      #genelist.uih2 <- unlist(strsplit(input$userheatgenes, "\t"))
      heatgenes <- c(genelist.ush,genelist.uih,genelist.uih2)
      heatgenes <- heatgenes[!is.na(heatgenes)]
      heatgenes <- heatgenes[!is.null(heatgenes)]
      usersamps <- input$userheatsamp2
      exp <- expr[heatgenes,usersamps]
      meta <- meta[which(meta[,1] %in% usersamps),]
      
      if (ncol(meta) > 2) {
        metacol <- input$CustomheatAnnoCol
      }
      else if (ncol(meta) == 2) {
        metacol <- colnames(meta)[2]
      }
      
      dataset <- exp
      dataset <- log2(dataset + 1)
      zdataset <- apply(dataset, 1, scale)
      zdataset <- apply(zdataset, 1, rev)
      colnames(zdataset) <- names(dataset)
      dataset <- as.matrix(zdataset)
      dataset[is.na(dataset)] <- 0
      
      
      dataset = dataset[apply(dataset, 1, function(x) !all(x==0)),]
      minimum = -5;
      maximum = 5;
      if (abs(min(dataset)) > abs(max(dataset))) {
        dataset[dataset < -abs(max(dataset))] = -abs(max(dataset))
      } else {
        dataset[dataset > abs(min(dataset))] = abs(min(dataset))
      }
      meta2 <- meta[,c(colnames(meta)[1],metacol)]
      meta2[,2] <- as.factor(meta2[,2])
      meta2 <- meta2[order(meta2[,2]),]
      rownames(meta2) <- meta2[,1]
      meta2 <- meta2[,-1,drop = F]
      dataset <- dataset[,rownames(meta2)]
      bk = c(seq(minimum,minimum/2, length=100), seq(minimum/2,maximum/2,length=100),seq(maximum/2,maximum,length=100))
      #Heatmap color
      color_choice <- input$ColorPalette2
      col_sets <- c("OrRd","PuBu","Greens","YlGnBu")
      if (color_choice == "original") {
        HeatMap_Colors <- c("dark blue","blue","white","red", "dark red")
        hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
      }
      else if (color_choice %in% col_sets) {
        HeatMap_Colors <- brewer.pal(n = 5, color_choice)
        hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
      }
      else if (color_choice == "Inferno") {
        hmcols <- inferno(500)
      }
      else if (color_choice == "Viridis") {
        hmcols <- viridis(500)
      }
      else if (color_choice == "Plasma") {
        hmcols <- plasma(500)
      }
      else if (color_choice == "OmniBlueRed") {
        hmcols<- colorRampPalette(c("#1984C5", "#22A7F0", "#63BFF0", "#A7D5ED", "#E2E2E2", "#E1A692", "#DE6E56", "#E14B31", "#C23728"))(length(bk)-1)
      }
      else if (color_choice == "LightBlueBlackRed") {
        hmcols<- colorRampPalette(c("#34C5FD","black","red"))(length(bk)-1)
      }
      else if (color_choice == "GreenBlackRed") {
        hmcols<- colorRampPalette(c("green","black","red"))(length(bk)-1)
      }
      dataset <- dataset[match(heatgenes, rownames(dataset)),]
      
      hm <- pheatmap::pheatmap(dataset,
                               cluster_col = clust_col_choice,
                               cluster_row = input$ClusRowOptCust,
                               fontsize_row = input$heatmapFont3.r,
                               fontsize_col = input$heatmapFont3.c,
                               show_rownames = row_names_choice ,
                               show_colnames = col_names_choice,
                               annotation_col = meta2,
                               clustering_method = input$ClusteringMethodCust,
                               color=hmcols,
                               border_color = NA)
      
    }
    
    hm
    
  })
  
  ####----DEG Heatmap----####
  
  # DEG heatmap reactive
  DEGHeatmap_react <- reactive({
    
    # UI Inputs
    row_names_choice <- input$ShowRowNames3     #choose to show row names or not
    col_names_choice <- input$ShowColNames3
    #A_choice <- input$comparisonA2_h            #Comparison group A
    #B_choice <- input$comparisonB2_h            #Comparison group B
    FC_cutoff <- input$fc_cutoff_h              #FC cutoff for top gene selection 
    P_cutoff <- input$p_cutoff_h                #P-value cutoff for top gene selections
    top_probes <- input$top_x_h                 #Number of top genes to show on heatmap
    clust_method <- input$ClusteringMethod_degh #Cluster method for heatmap
    col_font <- input$heatmapFont3.c.deg        #Heatmap column font size
    row_font <- input$heatmapFont3.r.deg        #Heatmap row font size
    color_choice <- input$ColorPalette3
    clust_row_choice <- input$ClusRowOptDEG
    clust_col_choice <- input$ClusColOptDEG
    
    top1 <- topgenereact2()
    
    if (ncol(meta) > 2) {
      metacol <- input$MetaColDEGHeat
    }
    else if (ncol(meta) == 2) {
      metacol <- colnames(meta)[2]
    }
    
    top_above_cutoff <- top1[which(top1$logFC > abs(FC_cutoff) & top1$P.Value < P_cutoff),]
    
    # Get gene list from Top table
    genelist <- rownames(top_above_cutoff)[c(1:top_probes)]
    
    # Heatmap Calculations
    dataset <- expr[which(rownames(expr) %in% genelist),]
    dataset <- log2(dataset + 1)
    zdataset <- apply(dataset, 1, scale)
    zdataset <- apply(zdataset, 1, rev)
    colnames(zdataset) <- names(dataset)
    dataset <- as.matrix(zdataset)
    dataset[is.na(dataset)] <- 0
    
    dataset = dataset[apply(dataset, 1, function(x) !all(x==0)),]
    minimum = -5;
    maximum = 5;
    if (abs(min(dataset)) > abs(max(dataset))) {
      dataset[dataset < -abs(max(dataset))] = -abs(max(dataset))
    } else {
      dataset[dataset > abs(min(dataset))] = abs(min(dataset))
    }
    meta2 <- meta[,c(colnames(meta)[1],metacol)]
    meta2[,2] <- as.factor(meta2[,2])
    meta2 <- meta2[order(meta2[,2]),]
    rownames(meta2) <- meta2[,1]
    meta2 <- meta2[,-1,drop = F]
    dataset <- dataset[,rownames(meta2)]
    bk = c(seq(minimum,minimum/2, length=100), seq(minimum/2,maximum/2,length=100),seq(maximum/2,maximum,length=100))
    #Heatmap color
    col_sets <- c("OrRd","PuBu","Greens","YlGnBu")
    if (color_choice == "original") {
      HeatMap_Colors <- c("dark blue","blue","white","red", "dark red")
      hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
    }
    else if (color_choice %in% col_sets) {
      HeatMap_Colors <- brewer.pal(n = 5, color_choice)
      hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
    }
    else if (color_choice == "Inferno") {
      hmcols <- inferno(500)
    }
    else if (color_choice == "Viridis") {
      hmcols <- viridis(500)
    }
    else if (color_choice == "Plasma") {
      hmcols <- plasma(500)
    }
    else if (color_choice == "OmniBlueRed") {
      hmcols<- colorRampPalette(c("#1984C5", "#22A7F0", "#63BFF0", "#A7D5ED", "#E2E2E2", "#E1A692", "#DE6E56", "#E14B31", "#C23728"))(length(bk)-1)
    }
    else if (color_choice == "LightBlueBlackRed") {
      hmcols<- colorRampPalette(c("#34C5FD","black","red"))(length(bk)-1)
    }
    else if (color_choice == "GreenBlackRed") {
      hmcols<- colorRampPalette(c("green","black","red"))(length(bk)-1)
    }
    # Heatmap
    hm <- pheatmap::pheatmap(dataset,
                             cluster_col = clust_col_choice,
                             cluster_row = clust_row_choice,
                             fontsize_row = row_font,
                             fontsize_col = col_font,
                             show_rownames = row_names_choice,
                             show_colnames = col_names_choice,
                             clustering_method = clust_method,
                             color=hmcols,
                             annotation_col = meta2,
                             angle_col = 90,
                             border_color = NA)
    hm
    
  })
  
  ####----Avg Expr Heatmap----####
  
  #Custom average heatmap reactive
  AvgExprHeatmap_react <- reactive({
    
    if (length(input$avgheatmapGeneSelec) >= 2 || length(input$avguserheatgenes) >= 1) {
      row_names_choice <- input$ShowRowNames5      #Chose to show row names or not
      col_names_choice <- input$ShowColNames5
      clust_row_opt <- input$AvgClusRowOptCust
      clust_col_opt <- input$AvgClusColOptCust
      group_choices <- input$SampCondSelectionCust
      col_font <- input$heatmapFont3.c.deg.avg.cust     #Heatmap column font size
      row_font <- input$heatmapFont3.r.deg.avg.cust     #Heatmap row font size
      clust_method <- input$AvgClusteringMethodCust
      genelist.uih <- NULL
      genelist.ush <- NULL
      genelist.uih2 <- NULL
      genelist.ush <- input$avgheatmapGeneSelec
      genelist.uih <- unlist(strsplit(input$avguserheatgenes, " "))
      #genelist.uih2 <- unlist(strsplit(input$userheatgenes, "\t"))
      heatgenes <- c(genelist.ush,genelist.uih,genelist.uih2)
      heatgenes <- heatgenes[!is.na(heatgenes)]
      heatgenes <- heatgenes[!is.null(heatgenes)]
      
      if (ncol(meta) > 2) {
        metacol <- input$AvgHeatMetaCol
      }
      else if (ncol(meta) == 2) {
        metacol <- colnames(meta)[2]
      }
      
      AvgExprDF <- data.frame(rownames(expr))
      
      for (i in group_choices) {
        samples <- meta[which(meta[,metacol] == i),1]
        if (length(samples) <= 1) {
          AvgExprDF[,paste("AvgExpr_",i, sep = "")] <- expr[,samples]
        }
        else if (length(samples) > 1) {
          AvgExprDF[,paste("AvgExpr_",i, sep = "")] <- rowMeans(expr[,samples])
        }
      }
      
      rownames(AvgExprDF) <- AvgExprDF[,1]
      AvgExprDF <- AvgExprDF[,-1]
      
      AvgExprDF <- AvgExprDF[which(rownames(AvgExprDF) %in% heatgenes),]
      
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
      #Heatmap color
      color_choice <- input$ColorPalette5
      col_sets <- c("OrRd","PuBu","Greens","YlGnBu")
      if (color_choice == "original") {
        HeatMap_Colors <- c("dark blue","blue","white","red", "dark red")
        hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
      }
      else if (color_choice %in% col_sets) {
        HeatMap_Colors <- brewer.pal(n = 5, color_choice)
        hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
      }
      else if (color_choice == "Inferno") {
        hmcols <- inferno(500)
      }
      else if (color_choice == "Viridis") {
        hmcols <- viridis(500)
      }
      else if (color_choice == "Plasma") {
        hmcols <- plasma(500)
      }
      else if (color_choice == "OmniBlueRed") {
        hmcols<- colorRampPalette(c("#1984C5", "#22A7F0", "#63BFF0", "#A7D5ED", "#E2E2E2", "#E1A692", "#DE6E56", "#E14B31", "#C23728"))(length(bk)-1)
      }
      else if (color_choice == "LightBlueBlackRed") {
        hmcols<- colorRampPalette(c("#34C5FD","black","red"))(length(bk)-1)
      }
      else if (color_choice == "GreenBlackRed") {
        hmcols<- colorRampPalette(c("green","black","red"))(length(bk)-1)
      }
      dataset <- dataset[match(heatgenes, rownames(dataset)),]
      
      hm <- pheatmap::pheatmap(as.matrix(dataset),
                               cluster_col = clust_col_opt,
                               cluster_row = clust_row_opt,
                               fontsize_row = row_font,
                               fontsize_col = col_font,
                               show_rownames = row_names_choice,
                               show_colnames = col_names_choice,
                               clustering_method = clust_method,
                               color=hmcols,
                               annotation_col = anno_meta,
                               angle_col = 90,
                               border_color = NA)
      
    }
    
    hm
    
  })
  
  ssgseaheatmap_react <- reactive({
    
    row_names_choice <- input$ShowRowNames1SSheat
    col_names_choice <- input$ShowColNamesSSheat
    row_font <- input$heatmapFont1.r
    col_font <- input$heatmapFont1.c
    if (ncol(meta) > 2) {
      metacol <- input$GSEAmetaCol
    }
    else if (ncol(meta) == 2) {
      metacol <- colnames(meta)[2]
    }
    if (is.null(input$ClusterMethodSSheat) == TRUE) {
      clust_method <- 'complete'
    }
    else if (is.null(input$ClusterMethodSSheat) == FALSE) {
      clust_method <- input$ClusterMethodSSheat
    }
    color_choice <- input$ColorPalette_gseaHeat
    clust_cols_opt <- input$clustcolsSSheat
    clust_rows_opt <- input$clustrowsSSheat
    scoreMethod <- input$ssGSEAtypeHeat
    
    if (is.null(input$ssgseaHeatGS) == TRUE) {
      #geneset_names <- gs_names_start
      geneset_names <- ssGSEA_Heat_GS()
    }
    else if (is.null(input$ssgseaHeatGS) == FALSE) {
      geneset_names <- input$ssgseaHeatGS
    }
    samples_chosen <- input$userheatsampSS
    
    A <- A[,samples_chosen]
    meta <- meta[which(meta[,1] %in% samples_chosen),]
    
    if (input$tables == 1) {
      GS <- gs[geneset_names]
    }
    if (input$tables == 3) {
      GS <- gs2[geneset_names]
    }
    if (input$tables == 5) {
      GS <- RDataListGen()[geneset_names]
    }
    ssgsea <- gsva(A, GS, method = scoreMethod, verbose = F)
    
    ssgsea2 = t(ssgsea)
    ssgsea3 = apply(ssgsea2, 2, scale);
    ssgsea4 = apply(ssgsea3, 1, rev)
    colnames(ssgsea4) = rownames(ssgsea2)
    
    neworder_gs <- rownames(ssgsea4)
    final_gs <- intersect(geneset_names,neworder_gs)
    
    ssgsea4 <- ssgsea4[final_gs,]
    
    minimum = min(ssgsea4)
    maximum = max(ssgsea4)
    if (abs(min(ssgsea4)) > abs(max(ssgsea4))) {
      ssgsea4[ssgsea4 < -abs(max(ssgsea4))] = -abs(max(ssgsea4))
    } else {
      ssgsea4[ssgsea4 > abs(min(ssgsea4))] = abs(min(ssgsea4))
    }
    bk = c(seq(minimum,minimum/2, length=100), seq(minimum/2,maximum/2,length=100),seq(maximum/2,maximum,length=100))
    #Heatmap color
    col_sets <- c("OrRd","PuBu","Greens","YlGnBu")
    if (color_choice == "original") {
      HeatMap_Colors <- c("dark blue","blue","white","red", "dark red")
      hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
    }
    else if (color_choice %in% col_sets) {
      HeatMap_Colors <- brewer.pal(n = 5, color_choice)
      hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
    }
    else if (color_choice == "Inferno") {
      hmcols <- inferno(500)
    }
    else if (color_choice == "Viridis") {
      hmcols <- viridis(500)
    }
    else if (color_choice == "Plasma") {
      hmcols <- plasma(500)
    }
    else if (color_choice == "OmniBlueRed") {
      hmcols<- colorRampPalette(c("#1984C5", "#22A7F0", "#63BFF0", "#A7D5ED", "#E2E2E2", "#E1A692", "#DE6E56", "#E14B31", "#C23728"))(length(bk)-1)
    }
    else if (color_choice == "LightBlueBlackRed") {
      hmcols<- colorRampPalette(c("#34C5FD","black","red"))(length(bk)-1)
    }
    else if (color_choice == "GreenBlackRed") {
      hmcols<- colorRampPalette(c("green","black","red"))(length(bk)-1)
    }
    
    meta2 <- meta[,c(colnames(meta)[1],metacol)]
    meta2[,2] <- as.factor(meta2[,2])
    meta2 <- meta2[order(meta2[,2]),]
    rownames(meta2) <- meta2[,1]
    meta2 <- meta2[,-1,drop = F]
    ssgsea4 <- ssgsea4[,rownames(meta2)]
    
    #meta <- meta[order(meta[,2]),]
    #type <- meta[,2]
    #meta2 <- as.data.frame(type)
    #rownames(meta2) <- meta[,1]
    
    hm <- pheatmap::pheatmap(ssgsea4,
                             cluster_col = clust_cols_opt,
                             cluster_row = clust_rows_opt,
                             fontsize_row = row_font,
                             fontsize_col = col_font,
                             show_rownames = row_names_choice,
                             show_colnames = col_names_choice,
                             annotation_col = meta2,
                             clustering_method = clust_method,
                             color = hmcols)
    
    hm
    
    
  })
  
  ssgseaheatmap2_react <- reactive({
    
    row_names_choice <- input$ShowRowNames1SSheat
    col_names_choice <- input$ShowColNamesSSheat
    row_font <- input$heatmapFont1.r
    col_font <- input$heatmapFont1.c
    if (ncol(meta) > 2) {
      metacol <- input$GSEAmetaCol
    }
    else if (ncol(meta) == 2) {
      metacol <- colnames(meta)[2]
    }
    if (is.null(input$ClusterMethodSSheat) == TRUE) {
      clust_method <- 'complete'
    }
    else if (is.null(input$ClusterMethodSSheat) == FALSE) {
      clust_method <- input$ClusterMethodSSheat
    }
    color_choice <- input$ColorPalette_gseaHeat
    clust_cols_opt <- input$clustcolsSSheat
    clust_rows_opt <- input$clustrowsSSheat
    scoreMethod <- input$ssGSEAtypeHeat
    
    if (is.null(input$ssgseaHeatGS) == TRUE) {
      #geneset_names <- gs_names_start
      geneset_names <- ssGSEA_Heat_GS()
    }
    else if (is.null(input$ssgseaHeatGS) == FALSE) {
      geneset_names <- input$ssgseaHeatGS
    }
    samples_chosen <- input$userheatsampSS
    
    A <- A[,samples_chosen]
    meta <- meta[which(meta[,1] %in% samples_chosen),]
    
    if (input$tables == 1) {
      GS <- gs[geneset_names]
    }
    if (input$tables == 3) {
      GS <- gs2[geneset_names]
    }
    if (input$tables == 5) {
      GS <- RDataListGen()[geneset_names]
    }
    ssgsea <- gsva(A, GS, method = scoreMethod, verbose = F, ssgsea.norm = F)
    
    SD=apply(ssgsea,1, sd, na.rm = TRUE) #get SD
    
    ssgsea2 = t(ssgsea)
    ssgsea3 = apply(ssgsea2, 2, scale);
    ssgsea4 = apply(ssgsea3, 1, rev)
    colnames(ssgsea4) = rownames(ssgsea2)
    
    ssgsea5 = ssgsea4 * SD #multiply zscore matrix by SD
    
    neworder_gs <- rownames(ssgsea5)
    final_gs <- intersect(geneset_names,neworder_gs)
    
    ssgsea5 <- ssgsea5[final_gs,]
    
    minimum = min(ssgsea5)
    maximum = max(ssgsea5)
    if (abs(min(ssgsea5)) > abs(max(ssgsea5))) {
      ssgsea5[ssgsea5 < -abs(max(ssgsea5))] = -abs(max(ssgsea5))
    } else {
      ssgsea5[ssgsea5 > abs(min(ssgsea5))] = abs(min(ssgsea5))
    }
    bk = c(seq(minimum,minimum/2, length=100), seq(minimum/2,maximum/2,length=100),seq(maximum/2,maximum,length=100))
    #Heatmap color
    col_sets <- c("OrRd","PuBu","Greens","YlGnBu")
    if (color_choice == "original") {
      HeatMap_Colors <- c("dark blue","blue","white","red", "dark red")
      hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
    }
    else if (color_choice %in% col_sets) {
      HeatMap_Colors <- brewer.pal(n = 5, color_choice)
      hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
    }
    else if (color_choice == "Inferno") {
      hmcols <- inferno(500)
    }
    else if (color_choice == "Viridis") {
      hmcols <- viridis(500)
    }
    else if (color_choice == "Plasma") {
      hmcols <- plasma(500)
    }
    else if (color_choice == "OmniBlueRed") {
      hmcols<- colorRampPalette(c("#1984C5", "#22A7F0", "#63BFF0", "#A7D5ED", "#E2E2E2", "#E1A692", "#DE6E56", "#E14B31", "#C23728"))(length(bk)-1)
    }
    else if (color_choice == "LightBlueBlackRed") {
      hmcols<- colorRampPalette(c("#34C5FD","black","red"))(length(bk)-1)
    }
    else if (color_choice == "GreenBlackRed") {
      hmcols<- colorRampPalette(c("green","black","red"))(length(bk)-1)
    }
    
    meta2 <- meta[,c(colnames(meta)[1],metacol)]
    meta2[,2] <- as.factor(meta2[,2])
    meta2 <- meta2[order(meta2[,2]),]
    rownames(meta2) <- meta2[,1]
    meta2 <- meta2[,-1,drop = F]
    ssgsea5 <- ssgsea5[,rownames(meta2)]
    
    #meta <- meta[order(meta[,2]),]
    #type <- meta[,2]
    #meta2 <- as.data.frame(type)
    #rownames(meta2) <- meta[,1]
    
    hm <- pheatmap::pheatmap(ssgsea5,
                             cluster_col = clust_cols_opt,
                             cluster_row = clust_rows_opt,
                             fontsize_row = row_font,
                             fontsize_col = col_font,
                             show_rownames = row_names_choice,
                             show_colnames = col_names_choice,
                             annotation_col = meta2,
                             clustering_method = clust_method,
                             color = hmcols)
    
    hm
    
    
  })
  
  ssgseaheatmap3_react <- reactive({
    
    row_names_choice <- input$ShowRowNames1SSheat
    col_names_choice <- input$ShowColNamesSSheat
    row_font <- input$heatmapFont1.r
    col_font <- input$heatmapFont1.c
    if (ncol(meta) > 2) {
      metacol <- input$GSEAmetaCol
    }
    else if (ncol(meta) == 2) {
      metacol <- colnames(meta)[2]
    }
    if (is.null(input$ClusterMethodSSheat) == TRUE) {
      clust_method <- 'complete'
    }
    else if (is.null(input$ClusterMethodSSheat) == FALSE) {
      clust_method <- input$ClusterMethodSSheat
    }
    color_choice <- input$ColorPalette_gseaHeat
    clust_cols_opt <- input$clustcolsSSheat
    clust_rows_opt <- input$clustrowsSSheat
    scoreMethod <- input$ssGSEAtypeHeat
    
    if (is.null(input$ssgseaHeatGS) == TRUE) {
      #geneset_names <- gs_names_start
      geneset_names <- ssGSEA_Heat_GS()
    }
    else if (is.null(input$ssgseaHeatGS) == FALSE) {
      geneset_names <- input$ssgseaHeatGS
    }
    samples_chosen <- input$userheatsampSS
    
    A <- A[,samples_chosen]
    meta <- meta[which(meta[,1] %in% samples_chosen),]
    
    if (input$tables == 1) {
      GS <- gs[geneset_names]
    }
    if (input$tables == 3) {
      GS <- gs2[geneset_names]
    }
    if (input$tables == 5) {
      GS <- RDataListGen()[geneset_names]
    }
    ssgsea <- gsva(A, GS, method = scoreMethod, verbose = F, ssgsea.norm = F)
    
    SD=apply(ssgsea,1, sd, na.rm = TRUE) #get SD
    
    ssgsea2 = t(ssgsea)
    ssgsea3 = apply(ssgsea2, 2, scale);
    ssgsea4 = apply(ssgsea3, 1, rev)
    colnames(ssgsea4) = rownames(ssgsea2)
    
    ssgsea5 = ssgsea4 * SD #multiply zscore matrix by SD
    
    group_choices <- unique(meta[,2])
    
    AvgssGSEADF <- data.frame(rownames(ssgsea5))
    
    for (i in group_choices) {
      samples <- meta[which(meta[,metacol] == i),1]
      if (length(samples) <= 1) {
        AvgssGSEADF[,paste("Avg_ssGSEA_",i, sep = "")] <- ssgsea5[,samples]
      }
      else if (length(samples) > 1) {
        AvgssGSEADF[,paste("Avg_ssGSEA_",i, sep = "")] <- rowMeans(ssgsea5[,samples])
      }
    }
    
    rownames(AvgssGSEADF) <- AvgssGSEADF[,1]
    AvgssGSEADF <- AvgssGSEADF[,-1]
    
    neworder_gs <- rownames(AvgssGSEADF)
    final_gs <- intersect(geneset_names,neworder_gs)
    
    AvgssGSEADF <- AvgssGSEADF[final_gs,]
    
    
    minimum = min(AvgssGSEADF)
    maximum = max(AvgssGSEADF)
    if (abs(min(AvgssGSEADF)) > abs(max(AvgssGSEADF))) {
      AvgssGSEADF[AvgssGSEADF < -abs(max(AvgssGSEADF))] = -abs(max(AvgssGSEADF))
    } else {
      AvgssGSEADF[AvgssGSEADF > abs(min(AvgssGSEADF))] = abs(min(AvgssGSEADF))
    }
    bk = c(seq(minimum,minimum/2, length=100), seq(minimum/2,maximum/2,length=100),seq(maximum/2,maximum,length=100))
    #Heatmap color
    col_sets <- c("OrRd","PuBu","Greens","YlGnBu")
    if (color_choice == "original") {
      HeatMap_Colors <- c("dark blue","blue","white","red", "dark red")
      hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
    }
    else if (color_choice %in% col_sets) {
      HeatMap_Colors <- brewer.pal(n = 5, color_choice)
      hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
    }
    else if (color_choice == "Inferno") {
      hmcols <- inferno(500)
    }
    else if (color_choice == "Viridis") {
      hmcols <- viridis(500)
    }
    else if (color_choice == "Plasma") {
      hmcols <- plasma(500)
    }
    else if (color_choice == "OmniBlueRed") {
      hmcols<- colorRampPalette(c("#1984C5", "#22A7F0", "#63BFF0", "#A7D5ED", "#E2E2E2", "#E1A692", "#DE6E56", "#E14B31", "#C23728"))(length(bk)-1)
    }
    else if (color_choice == "LightBlueBlackRed") {
      hmcols<- colorRampPalette(c("#34C5FD","black","red"))(length(bk)-1)
    }
    else if (color_choice == "GreenBlackRed") {
      hmcols<- colorRampPalette(c("green","black","red"))(length(bk)-1)
    }
    
    anno_meta <- data.frame(colnames(AvgssGSEADF))
    anno_meta[,metacol] <- gsub("Avg_ssGSEA_","",anno_meta[,1])
    rownames(anno_meta) <- anno_meta[,1]
    anno_meta <- anno_meta[,-1,drop = F]
    
    hm <- pheatmap::pheatmap(as.matrix(AvgssGSEADF),
                             cluster_col = clust_cols_opt,
                             cluster_row = clust_rows_opt,
                             fontsize_row = row_font,
                             fontsize_col = col_font,
                             show_rownames = row_names_choice,
                             show_colnames = col_names_choice,
                             angle_col = 90,
                             annotation_col = anno_meta,
                             clustering_method = clust_method,
                             color = hmcols)
    
    hm
    
    
  })
  
  ####----Gene Scatter----####
  
  #render gene expression comparison scatter plot
  output$geneScatter0 <- renderPlotly({
    
    title_font <- input$GeneScatterTitleSize
    axis_font <- input$GeneScatterAxisSize
    if (ncol(meta) > 2){
      req(input$ScatterPlotMetaCol)
      metacol <- input$ScatterPlotMetaCol
    }
    else if (ncol(meta) == 2){
      metacol <- colnames(meta)[2]
    }
    
    
    #log if user designates
    if (input$logask == TRUE) {
      expr <- log2(expr + 1)
    }
    #transpose
    expr_t <- as.data.frame(t(expr))
    #reorder rowname to match meta for merging
    samporder <- meta[,1]
    expr_t2 <- as.data.frame(expr_t[samporder,])
    #add type
    expr_t3 <- expr_t2 %>% 
      mutate(type = case_when(
        rownames(expr_t2) == meta[,1] ~ meta[,metacol],
      ))
    expr_t3 <- expr_t3 %>%
      relocate(type)
    expr_t3$type <- as.factor(expr_t3$type)
    colnames(expr_t3)[which(colnames(expr_t3) == "type")] <- metacol
    #user gene input
    gene1.u <- input$scatterG1
    gene2.u <- input$scatterG2
    #get columns and info based off user input
    gene1 <- expr_t3[,gene1.u]
    gene2 <- expr_t3[,gene2.u]
    if (input$logask == TRUE) {
      gene1.n <- paste(colnames(expr_t3)[which(colnames(expr_t3) == gene1.u)], " Expression (log2 +1)", sep = "")
      gene2.n <- paste(colnames(expr_t3)[which(colnames(expr_t3) == gene2.u)], " Expression (log2 +1)", sep = "")
    }
    else if (input$logask == FALSE) {
      gene2.n <- paste(colnames(expr_t3)[which(colnames(expr_t3) == gene2.u)], " Expression", sep = "")
      gene1.n <- paste(colnames(expr_t3)[which(colnames(expr_t3) == gene1.u)], " Expression", sep = "")
    }
    #plot
    p <- ggplot(expr_t3, aes(x = gene1, y = gene2,
                             color = expr_t3[,metacol],
                             text = paste("</br> Sample: ", rownames(expr_t3),
                                          "</br> ",metacol,": ", expr_t3[,metacol],
                                          sep = ""))) +
      geom_point() +
      theme_minimal() +
      labs(x = gene1.n, y = gene2.n,
           title = paste(gene1.n, " vs. ", gene2.n, sep = ''),
           color = metacol) +
      theme(axis.title = element_text(size = axis_font),
            plot.title = element_text(size = title_font))
    
    ggplotly(p, tooltip = 'text')
    
  })
  
  ####----Data Tables----####
  
  
  #render MSigDB gene set table
  output$msigdbTable <- DT::renderDataTable({
    DT::datatable(msigdb.gsea2,
                  selection = list(mode = 'single', selected = 1),
                  options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100")))
    
  })
  
  #create user input gene set table
  output$GStable.u <- DT::renderDataTable({
    gmt.u <- input$user.gmt.file
    ext <- tools::file_ext(gmt.u$datapath)
    req(gmt.u)
    validate(need(ext == c("gmt","tsv","txt"), "Please upload .gmt, .tsv, or .txt file"))
    if (ext == "gmt") {
      gmt.us <- read.gmt(gmt.u$datapath)
    }
    else {
      gmt.us <- as.data.frame(read_delim(gmt.u$datapath, delim = '\t'))
      #colnames(gmt.us) <- c("term","gene")
    }
    GeneSet <- as.data.frame(unique(gmt.us[1]))
    rownames(GeneSet) <- 1:nrow(GeneSet)
    colnames(GeneSet)[1] <- "Gene_Set"
    DT::datatable(GeneSet,
                  selection = list(mode = 'single', selected = 1),
                  options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100")))
  })
  
  #render tab2 gene set table
  output$tab2table <- DT::renderDataTable({
    DT::datatable(GeneSet2,
                  selection = list(mode = 'single', selected = 1),
                  options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100")))
    
  })
  
  #render leading edge genes list
  output$LeadingEdgeGenes <- DT::renderDataTable({
    if (input$tables == 1){
      if (length(input$msigdbTable_rows_selected) > 0){
        res <- datasetInput()
        gsea.df <- as.data.frame(res@result)
        GS <- as.character(msigdb.gsea2[input$msigdbTable_rows_selected,3])
        ## Subset core enriched genes
        genes1 <- as.matrix(gsea.df[which(gsea.df$Description==GS),"core_enrichment"])
        genes2 <- strsplit(genes1,"/")
        GeneSymbol <- as.data.frame(genes2, col.names = "GeneSymbol")
        GeneSymbol$Rank <- rownames(GeneSymbol)
        GeneSymbol <- GeneSymbol[,c("Rank","GeneSymbol")]
        DT::datatable(GeneSymbol, options = list(paging = F), rownames = F)
      }
    }
    else if (input$tables == 3){
      if (length(input$tab2table_rows_selected) > 0){
        res <- datasetInput()
        gsea.df <- as.data.frame(res@result)
        GS <- as.character(GeneSet2[input$tab2table_rows_selected,1])
        ## Subset core enriched genes
        genes1 <- as.matrix(gsea.df[which(gsea.df$Description==GS),"core_enrichment"])
        genes2 <- strsplit(genes1,"/")
        GeneSymbol <- as.data.frame(genes2, col.names = "GeneSymbol")
        GeneSymbol$Rank <- rownames(GeneSymbol)
        GeneSymbol <- GeneSymbol[,c("Rank","GeneSymbol")]
        DT::datatable(GeneSymbol, options = list(paging = F), rownames = F)
      }
    }
    else if (input$tables == 5){
      if (length(input$GStable.u_rows_selected) > 0){
        res <- datasetInput()
        gsea.df <- as.data.frame(res@result)
        GS <- as.character(user_gs_mirror()[input$GStable.u_rows_selected,1])
        ## Subset core enriched genes
        genes1 <- as.matrix(gsea.df[which(gsea.df$Description==GS),"core_enrichment"])
        genes2 <- strsplit(genes1,"/")
        GeneSymbol <- as.data.frame(genes2, col.names = "GeneSymbol")
        GeneSymbol$Rank <- rownames(GeneSymbol)
        GeneSymbol <- GeneSymbol[,c("Rank","GeneSymbol")]
        DT::datatable(GeneSymbol, options = list(paging = F), rownames = F)
      }
    }
  })
  
  #render pre-loaded enriched signatures table
  output$enrich_sig_table <- DT::renderDataTable({
    gsea.df <- as_tibble(get(ES_Tab_List[which(ES_Tab_List == paste("ES_table",match(input$SigTableChoice, SigNames),sep = ""))]))
    DT::datatable(gsea.df,
                  extensions = c("KeyTable", "FixedHeader"),
                  caption = "Enriched Signatures",
                  options = list(keys = T, searchHighlight = T, pageLength = 10, lengthMenu = c("10", "25", "50", "100"),scrollX = T)) %>%
      formatRound(columns = c(2:10), digits = 2)
  })
  
  #render user generated enriched signature table based off other gmt
  output$enrich_sig_table_gen <- DT::renderDataTable({
    if (input$tables == 3) {
      if (ncol(meta) > 2) {
        req(input$GSEAmetaCol)
        metacol <- input$GSEAmetaCol
      }
      else if (ncol(meta) == 2) {
        metacol <- colnames(meta)[2]
      }
      groupA <- meta[which(meta[,metacol] == input$comparisonA),1]
      groupB <- meta[which(meta[,metacol] == input$comparisonB),1]
      ##----Signal-to-Noise Calculation----##
      A <- A + 0.00000001
      P = as.matrix(as.numeric(colnames(A) %in% groupA))
      n1 <- sum(P[,1])
      M1 <- A %*% P
      M1 <- M1/n1
      A2 <- A*A
      S1 <- A2 %*% P
      S1 <- S1/n1 - M1*M1 
      S1 <- sqrt(abs((n1/(n1-1)) * S1))
      P = as.matrix(as.numeric(colnames(A) %in% groupB))
      n2 <- sum(P[,1])
      M2 <- A %*% P
      M2 <- M2/n2
      A2 <- A*A
      S2 <- A2 %*% P
      S2 <- S2/n2 - M2*M2
      S2 <- sqrt(abs((n2/(n2-1)) * S2))
      rm(A2)
      # small sigma "fix" as used in GeneCluster
      S2 <- ifelse(0.2*abs(M2) < S2, S2, 0.2*abs(M2))
      S2 <- ifelse(S2 == 0, 0.2, S2)
      S1 <- ifelse(0.2*abs(M1) < S1, S1, 0.2*abs(M1))
      S1 <- ifelse(S1 == 0, 0.2, S1)
      M1 <- M1 - M2
      rm(M2)
      S1 <- S1 + S2
      rm(S2)
      s2n.matrix <- M1/S1
      ##----Reformatting----##
      s2n.df <- as.data.frame(s2n.matrix)
      s2n.df$GeneID <- rownames(s2n.df)
      rownames(s2n.df) <- NULL
      data <- dplyr::select(s2n.df, GeneID, V1)
      data.gsea <- data$V1
      names(data.gsea) <- as.character(data$GeneID)
      s2n.matrix.s <- sort(data.gsea, decreasing = T)
      ##----GSEA----##
      gmt.i <- tab2
      gsea.res <- GSEA(s2n.matrix.s, TERM2GENE = gmt.i, verbose = F, pvalueCutoff = input$userPval)
      gsea.df <- as_tibble(gsea.res@result)
      ## displaying the GSEA results as interactive data table
      DT::datatable(gsea.df,
                    extensions = c("KeyTable", "FixedHeader"),
                    caption = "Enriched Signatures",
                    options = list(keys = T, searchHighlight = T, pageLength = 10, lengthMenu = c("10", "25", "50", "100"), scrollX = T)) %>%
        formatRound(columns = c(2:10), digits = 2)
    }
    else if (input$tables == 5) {
      if (ncol(meta) > 2) {
        metacol <- input$GSEAmetaCol
      }
      else if (ncol(meta) == 2) {
        metacol <- colnames(meta)[2]
      }
      groupA <- meta[which(meta[,metacol] == input$comparisonA),1]
      groupB <- meta[which(meta[,metacol] == input$comparisonB),1]
      ##----Signal-to-Noise Calculation----##
      A <- A + 0.00000001
      P = as.matrix(as.numeric(colnames(A) %in% groupA))
      n1 <- sum(P[,1])
      M1 <- A %*% P
      M1 <- M1/n1
      A2 <- A*A
      S1 <- A2 %*% P
      S1 <- S1/n1 - M1*M1 
      S1 <- sqrt(abs((n1/(n1-1)) * S1))
      P = as.matrix(as.numeric(colnames(A) %in% groupB))
      n2 <- sum(P[,1])
      M2 <- A %*% P
      M2 <- M2/n2
      A2 <- A*A
      S2 <- A2 %*% P
      S2 <- S2/n2 - M2*M2
      S2 <- sqrt(abs((n2/(n2-1)) * S2))
      rm(A2)
      # small sigma "fix" as used in GeneCluster
      S2 <- ifelse(0.2*abs(M2) < S2, S2, 0.2*abs(M2))
      S2 <- ifelse(S2 == 0, 0.2, S2)
      S1 <- ifelse(0.2*abs(M1) < S1, S1, 0.2*abs(M1))
      S1 <- ifelse(S1 == 0, 0.2, S1)
      M1 <- M1 - M2
      rm(M2)
      S1 <- S1 + S2
      rm(S2)
      s2n.matrix <- M1/S1
      ##----Reformatting----##
      s2n.df <- as.data.frame(s2n.matrix)
      s2n.df$GeneID <- rownames(s2n.df)
      rownames(s2n.df) <- NULL
      data <- dplyr::select(s2n.df, GeneID, V1)
      data.gsea <- data$V1
      names(data.gsea) <- as.character(data$GeneID)
      s2n.matrix.s <- sort(data.gsea, decreasing = T)
      ##----GSEA----##
      gmt.i <- GStable.ubg()
      gsea.res <- GSEA(s2n.matrix.s, TERM2GENE = gmt.i, verbose = F, pvalueCutoff = input$userPval)
      gsea.df <- as_tibble(gsea.res@result)
      ## displaying the GSEA results as interactive data table
      DT::datatable(gsea.df,
                    extensions = c("KeyTable", "FixedHeader"),
                    caption = "Enriched Signatures",
                    options = list(keys = T, searchHighlight = T, pageLength = 10, lengthMenu = c("10", "25", "50", "100"), scrollX = T)) %>%
        formatRound(columns = c(2:10), digits = 2)
    }
    else if (input$tables == 1) {
      if (input$GenerateEST == TRUE) {
        #extract results and convert to tibble
        gsea.df <- GeneratedMSigDBEST()
        gsea.df <- as_tibble(gsea.df)
        ## displaying the GSEA results as interactive data table
        DT::datatable(gsea.df,
                      extensions = c("KeyTable", "FixedHeader"),
                      caption = "Enriched Signatures",
                      options = list(keys = T, searchHighlight = T, pageLength = 10, lengthMenu = c("10", "25", "50", "100"), scrollX = T)) %>%
          formatRound(columns = c(2:10), digits = 2)
      }
    }
  })
  
  #render DEG table
  output$DEGtable1 <- DT::renderDataTable({
    if (ncol(meta) > 2) {
      metacol <- input$DEGtableMetaCol
    }
    else if (ncol(meta) == 2) {
      metacol <- colnames(meta)[2]
    }
    A <- meta[which(meta[,metacol] == input$comparisonA2.DEG),1]
    B <- meta[which(meta[,metacol] == input$comparisonB2.DEG),1]
    mat <- expr[,c(A,B)]
    mat <- log2(mat + 1.0)
    groupAOther <- factor(c(rep("A", length(A)), rep("B", length(B))))
    designA <- model.matrix(~0 + groupAOther)
    fit <- lmFit(mat, design = designA)
    contrast.matrix <- makeContrasts(groupAOtherA - groupAOtherB, levels = designA)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    options(digits = 4)
    top1 <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
    DT::datatable(top1, options = list(lengthMenu = c(50,100,1000, 5000, 10000), pageLength = 100, scrollX = TRUE),
                  selection=list(mode = "multiple"))
  })
  
  #render up regulated pathway enrichment data table
  output$UpRegPathwayTable1 <- DT::renderDataTable({
    adjp <- input$pathpval
    FC <- input$pathFC
    if (ncol(meta) > 2) {
      metacol <- input$EnrichRPathMetaCol
    }
    else if (ncol(meta) == 2) {
      metacol <- colnames(meta)[2]
    }
    A <- meta[which(meta[,metacol] == input$comparisonA2.path),1]
    B <- meta[which(meta[,metacol] == input$comparisonB2.path),1]
    mat <- expr[,c(A,B)]
    mat <- log2(mat + 1.0)
    groupAOther <- factor(c(rep("A", length(A)), rep("B", length(B))))
    designA <- model.matrix(~0 + groupAOther)
    fit <- lmFit(mat, design = designA)
    contrast.matrix <- makeContrasts(groupAOtherA - groupAOtherB, levels = designA)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    options(digits = 4)
    top1 <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
    genes <- rownames(top1)[which(top1$P.Value < adjp & top1$logFC > FC)]
    dbs <- listEnrichrDbs() 
    enrichRLive <- TRUE 
    if (is.null(dbs)) { 
      enrichRLive <- FALSE 
    }
    dbs <- input$SelectedPathway
    enriched <- enrichr(genes, dbs) # Plot top 20 GO-BP results ordered by P-value
    DT::datatable(enriched[[1]], options = list(lengthMenu = c(10,20,50,100,1000), pageLength = 10, scrollX = TRUE),
                  selection=list(mode = "multiple"))
  })
  
  #render down regulated pathway enrichment data table
  output$DnRegPathwayTable1 <- DT::renderDataTable({
    adjp <- input$pathpval
    FC <- input$pathFC
    if (ncol(meta) > 2) {
      metacol <- input$EnrichRPathMetaCol
    }
    else if (ncol(meta) == 2) {
      metacol <- colnames(meta)[2]
    }
    A <- meta[which(meta[,metacol] == input$comparisonA2.path),1]
    B <- meta[which(meta[,metacol] == input$comparisonB2.path),1]
    mat <- expr[,c(A,B)]
    mat <- log2(mat + 1.0)
    groupAOther <- factor(c(rep("A", length(A)), rep("B", length(B))))
    designA <- model.matrix(~0 + groupAOther)
    fit <- lmFit(mat, design = designA)
    contrast.matrix <- makeContrasts(groupAOtherA - groupAOtherB, levels = designA)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    options(digits = 4)
    top1 <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
    genes <- rownames(top1)[which(top1$P.Value < adjp & top1$logFC > FC)]
    dbs <- listEnrichrDbs() 
    enrichRLive <- TRUE 
    if (is.null(dbs)) { 
      enrichRLive <- FALSE 
    }
    dbs <- input$SelectedPathway
    enriched <- enrichr(genes, dbs) # Plot top 20 GO-BP results ordered by P-value
    DT::datatable(enriched[[1]], options = list(lengthMenu = c(10,20,50,100,1000), pageLength = 10, scrollX = TRUE),
                  selection=list(mode = "multiple"))
  })
  
  #render gene list table for boxplot selection
  output$GeneListTable <- DT::renderDataTable({
    DT::datatable(geneList,
                  selection = list(mode = 'single', selected = 1),
                  options = list(keys = TRUE,
                                 searchHighlight = TRUE,
                                 pageLength = 10,
                                 lengthMenu = c("10", "25", "50", "100")))
  })
  
  #render gene list table for boxplot selection
  output$GeneListTable2 <- DT::renderDataTable({
    DT::datatable(geneList,
                  selection = list(mode = 'single', selected = 1),
                  options = list(keys = TRUE,
                                 searchHighlight = TRUE,
                                 pageLength = 10,
                                 lengthMenu = c("10", "25", "50", "100")),
                  rownames = F)
  })
  
  #render gene list table for boxplot selection
  output$GeneListTableBarPlot <- DT::renderDataTable({
    DT::datatable(geneList,
                  selection = list(mode = 'single', selected = 1),
                  options = list(keys = TRUE,
                                 searchHighlight = TRUE,
                                 pageLength = 10,
                                 lengthMenu = c("10", "25", "50", "100")),
                  rownames = F)
  })
  
  #render gene scatter plot data table
  output$geneScatterTable <- DT::renderDataTable({
    #log if user designates
    if (input$logask == TRUE) {
      expr <- log2(expr + 1)
    }
    #transpose
    expr_t <- as.data.frame(t(expr))
    if (ncol(meta) > 2) {
      metacol <- input$ScatterPlotMetaCol
    }
    else if (ncol(meta) == 2) {
      metacol <- colnames(meta)[2]
    }
    #reorder rowname to match meta for merging
    samporder <- meta[,1]
    expr_t2 <- as.data.frame(expr_t[samporder,])
    #add type
    expr_t3 <- expr_t2 %>% 
      mutate(type = case_when(
        rownames(expr_t2) == meta[,1] ~ meta[,metacol],
      ))
    expr_t3 <- expr_t3 %>%
      relocate(type)
    colnames(expr_t3)[which(colnames(expr_t3) == "type")] <- metacol
    #user gene input
    gene1.u <- input$scatterG1
    gene2.u <- input$scatterG2
    #get columns and info based off user input
    gene1 <- expr_t3[,gene1.u]
    gene2 <- expr_t3[,gene2.u]
    if (input$logask == TRUE) {
      gene1.n <- paste(colnames(expr_t3)[which(colnames(expr_t3) == gene1.u)], " Expression (log2 +1)", sep = "")
      gene2.n <- paste(colnames(expr_t3)[which(colnames(expr_t3) == gene2.u)], " Expression (log2 +1)", sep = "")
    }
    else if (input$logask == FALSE) {
      gene2.n <- paste(colnames(expr_t3)[which(colnames(expr_t3) == gene2.u)], " Expression", sep = "")
      gene1.n <- paste(colnames(expr_t3)[which(colnames(expr_t3) == gene1.u)], " Expression", sep = "")
    }
    #make table
    Sample <- rownames(expr_t3)
    Type <- expr_t3[,metacol]
    gene1col <- expr_t3[,gene1.u]
    gene2col <- expr_t3[,gene2.u]
    scatterTab <- data.frame(Sample, Type, gene1col, gene2col)
    colnames(scatterTab)[c(2,3,4)] <- c(metacol,gene1.n, gene2.n)
    #table output,
    DT::datatable(scatterTab,
                  options = list(keys = TRUE,
                                 searchHighlight = TRUE,
                                 pageLength = 10,
                                 lengthMenu = c("10", "25", "50", "100")))
    
  })
  
  output$AvggeneScatterTable <- renderDataTable({
    
    A_choice <- input$comparisonA2.avg
    B_choice <- input$comparisonB2.avg
    log_choice <- input$AvgExpLogFC
    if (ncol(meta) > 2){
      metacol <- input$AvgExprScatterMetaCol
    }
    else if (ncol(meta) == 2){
      metacol <- colnames(meta)[2]
    }
    
    # Get A and B sample names
    A <- meta[which(meta[,metacol] == A_choice),1]
    B <- meta[which(meta[,metacol] == B_choice),1]
    
    # Make A and B expression Matrices
    mat_A <- expr[,A]
    mat_B <- expr[,B]
    
    logFC_text <- ""
    
    # Log if user designates
    if (log_choice == TRUE) {
      mat_A <- log2(mat_A + 1)
      mat_B <- log2(mat_B + 1)
      logFC_text <- " (log2 +1)"
    }
    
    # Get avg expression of each gene
    mat_A$AvgExpression_GroupA <- rowMeans(mat_A)
    mat_B$AvgExpression_GroupB <- rowMeans(mat_B)
    
    # Add Gene Names as column
    mat_A$GeneSymbol <- rownames(mat_A)
    mat_B$GeneSymbol <- rownames(mat_B)
    
    # Merge average columns
    AvgExpr_Table <- merge(mat_A[,c("AvgExpression_GroupA","GeneSymbol")],
                           mat_B[,c("AvgExpression_GroupB","GeneSymbol")],
                           by = "GeneSymbol",
                           all = T)
    
    colnames(AvgExpr_Table)[2] <- paste("Average Expression ", A_choice, sep = "")
    colnames(AvgExpr_Table)[3] <- paste("Average Expression ", B_choice, sep = "")
    
    #table output
    DT::datatable(AvgExpr_Table,
                  options = list(keys = TRUE,
                                 searchHighlight = TRUE,
                                 pageLength = 10,
                                 lengthMenu = c("10", "25", "50", "100")),
                  rownames = F)
    
  })
  
  #render ssGSEA table
  output$ssGSEAtable <- DT::renderDataTable({
    if (ncol(meta) > 2) {
      metacol <- input$GSEAmetaCol
    }
    else if (ncol(meta) == 2) {
      metacol <- colnames(meta)[2]
    }
    if (input$tables == 3) {
      if (length(input$tab2table_rows_selected) > 0){
        ssgsea <- ssGSEAfunc()
        ssgsea2 <- as.data.frame(t(ssgsea))
        samporder <- meta[,1]
        ssgsea3 <- as.data.frame(ssgsea2[samporder,])
        colnames(ssgsea3)[1] <- colnames(ssgsea2[1])
        rownames(ssgsea3) <- samporder
        ssgsea4 <- ssgsea3 %>% 
          mutate(type = case_when(
            rownames(ssgsea3) == meta[,1] ~ meta[,metacol],
          ))
        ssgsea4 <- ssgsea4 %>%
          relocate(type)
        ssgsea4$type <- as.factor(ssgsea4$type)
        colnames(ssgsea4)[which(colnames(ssgsea4) == "type")] <- metacol
        #table output
        DT::datatable(ssgsea4,
                      options = list(keys = TRUE,
                                     searchHighlight = TRUE,
                                     pageLength = 10,
                                     lengthMenu = c("10", "25", "50", "100")))
      }
    }
    else if (input$tables == 1) {
      if (length(input$msigdbTable_rows_selected) > 0){
        ssgsea <- ssGSEAfunc()
        ssgsea2 <- as.data.frame(t(ssgsea))
        samporder <- meta[,1]
        ssgsea3 <- as.data.frame(ssgsea2[samporder,])
        colnames(ssgsea3)[1] <- colnames(ssgsea2[1])
        rownames(ssgsea3) <- samporder
        ssgsea4 <- ssgsea3 %>% 
          mutate(type = case_when(
            rownames(ssgsea3) == meta[,1] ~ meta[,metacol],
          ))
        ssgsea4 <- ssgsea4 %>%
          relocate(type)
        ssgsea4$type <- as.factor(ssgsea4$type)
        colnames(ssgsea4)[which(colnames(ssgsea4) == "type")] <- metacol
        #table output
        DT::datatable(ssgsea4,
                      options = list(keys = TRUE,
                                     searchHighlight = TRUE,
                                     pageLength = 10,
                                     lengthMenu = c("10", "25", "50", "100")))
      }
    }
    else if (input$tables == 5) {
      if (length(input$GStable.u_rows_selected) > 0){
        ssgsea <- ssGSEAfunc()
        ssgsea2 <- as.data.frame(t(ssgsea))
        samporder <- meta[,1]
        ssgsea3 <- as.data.frame(ssgsea2[samporder,])
        colnames(ssgsea3)[1] <- colnames(ssgsea2[1])
        rownames(ssgsea3) <- samporder
        ssgsea4 <- ssgsea3 %>% 
          mutate(type = case_when(
            rownames(ssgsea3) == meta[,1] ~ meta[,metacol],
          ))
        ssgsea4 <- ssgsea4 %>%
          relocate(type)
        ssgsea4$type <- as.factor(ssgsea4$type)
        colnames(ssgsea4)[which(colnames(ssgsea4) == "type")] <- metacol
        #table output
        DT::datatable(ssgsea4,
                      options = list(keys = TRUE,
                                     searchHighlight = TRUE,
                                     pageLength = 10,
                                     lengthMenu = c("10", "25", "50", "100")))
      }
    }
  })
  
  
  ####----Plots----####
  
  
  #render GSEA plot
  output$enrichplot0 <- renderPlot({
    if (input$tables == 1){
      if (length(input$msigdbTable_rows_selected) > 0){
        res <- datasetInput()
        gsea.df <- as.data.frame(res@result)
        geneset <- as.character(msigdb.gsea2[input$msigdbTable_rows_selected,3])
        title <- geneset
        gseaplot2(res,
                  geneset,
                  title,
                  pvalue_table = F)
      }
    }
    else if (input$tables == 3){
      if (length(input$tab2table_rows_selected) > 0){
        res <- datasetInput()
        gsea.df <- as.data.frame(res@result)
        geneset <- as.character(GeneSet2[input$tab2table_rows_selected,1])
        title <- geneset
        gseaplot2(res,
                  geneset,
                  title,
                  pvalue_table = F)
      }
    }
    else if (input$tables == 5){
      if (length(input$GStable.u_rows_selected) > 0){
        res <- datasetInput()
        gsea.df <- as.data.frame(res@result)
        geneset <- as.character(user_gs_mirror()[input$GStable.u_rows_selected,1])
        title <- geneset
        gseaplot2(res,
                  geneset,
                  title,
                  pvalue_table = F)
      }
    }
  })
  
  #render heatmap
  output$heatmap0 <- renderPlot({
    
    if (input$tables == 1) {
      if (length(input$msigdbTable_rows_selected) > 0){
        if (ncol(meta) > 2) {
          req(input$GSEAmetaCol)
          metacol <- input$GSEAmetaCol
        }
        else if (ncol(meta) == 2) {
          metacol <- colnames(meta)[2]
        }
        groupA <- meta[which(meta[,metacol] == input$comparisonA),1]
        groupB <- meta[which(meta[,metacol] == input$comparisonB),1]
        res <- datasetInput()
        gsea.df <- as.data.frame(res@result)
        GS <- as.character(msigdb.gsea2[input$msigdbTable_rows_selected,3])
        genes1 <- as.matrix(gsea.df[which(gsea.df$Description==GS),"core_enrichment"])
        genes2 <- strsplit(genes1,"/")
        genes3 <- as.data.frame(genes2, col.names = "genes")
        gene_symbol <- genes3$genes
        ## convert expression matrix to numeric
        class(A) <- "numeric"
        ## Transforming data
        A <- A[,c(groupA,groupB)]
        exp.mat1 = log2(A + 1) # log
        exp.mat2 = apply(exp.mat1, 1, scale); # z score
        exp.mat3 = apply(exp.mat2, 1, rev); # transpose
        colnames(exp.mat3) = colnames(A) # set the column name
        exp.mat4 = exp.mat3[rowSums(is.na(exp.mat3[,])) == 0,] # remove NaN rows
        exp.mat5 = exp.mat4[rownames(exp.mat4) %in% gene_symbol,] # grab gene symbols
        # reassign data
        dataset <- exp.mat5
        ## generate color for pheatmap
        if (abs(min(dataset)) > abs(max(dataset))) {
          dataset[dataset < -abs(max(dataset))] = -abs(max(dataset))
        } else {
          dataset[dataset > abs(min(dataset))] = abs(min(dataset))
        }
        meta2 <- meta[which(meta[,1] %in% c(groupA,groupB)),]
        
        meta2 <- meta2[,c(colnames(meta2)[1],metacol)]
        meta2[,2] <- as.factor(meta2[,2])
        meta2 <- meta2[order(meta2[,2]),]
        rownames(meta2) <- meta2[,1]
        meta3 <- meta2[,-1,drop = F]
        dataset <- dataset[,rownames(meta2)]
        
        zscore_range = 10;
        minimum = -zscore_range;
        maximum = zscore_range;
        bk = c(seq(minimum,minimum/2, length=100), seq(minimum/2,maximum/2,length=100),seq(maximum/2,maximum,length=100))
        #Heatmap color
        color_choice <- input$ColorPalette_gseaHeat
        col_sets <- c("OrRd","PuBu","Greens","YlGnBu")
        if (color_choice == "original") {
          HeatMap_Colors <- c("dark blue","blue","white","red", "dark red")
          hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
        }
        else if (color_choice %in% col_sets) {
          HeatMap_Colors <- brewer.pal(n = 5, color_choice)
          hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
        }
        else if (color_choice == "Inferno") {
          hmcols <- inferno(500)
        }
        else if (color_choice == "Viridis") {
          hmcols <- viridis(500)
        }
        else if (color_choice == "Plasma") {
          hmcols <- plasma(500)
        }
        else if (color_choice == "OmniBlueRed") {
          hmcols<- colorRampPalette(c("#1984C5", "#22A7F0", "#63BFF0", "#A7D5ED", "#E2E2E2", "#E1A692", "#DE6E56", "#E14B31", "#C23728"))(length(bk)-1)
        }
        else if (color_choice == "LightBlueBlackRed") {
          hmcols<- colorRampPalette(c("#34C5FD","black","red"))(length(bk)-1)
        }
        else if (color_choice == "GreenBlackRed") {
          hmcols<- colorRampPalette(c("green","black","red"))(length(bk)-1)
        }
        pheatmap::pheatmap(dataset,            #data
                           cluster_cols = F,    #cluster columns - NO
                           cluster_row = F,     #cluster rows - YES
                           fontsize_col = input$heatmapFont1.c,   #column fontsize
                           fontsize_row = input$heatmapFont1.r,
                           show_rownames = T,  
                           show_colnames = T,
                           color=hmcols,
                           annotation_col = meta3)
      }
    }
    else if (input$tables == 3) {
      if (length(input$tab2table_rows_selected) > 0){
        if (ncol(meta) > 2) {
          req(input$GSEAmetaCol)
          metacol <- input$GSEAmetaCol
        }
        else if (ncol(meta) == 2) {
          metacol <- colnames(meta)[2]
        }
        groupA <- meta[which(meta[,metacol] == input$comparisonA),1]
        groupB <- meta[which(meta[,metacol] == input$comparisonB),1]
        res <- datasetInput()
        gsea.df <- as.data.frame(res@result)
        GS <- as.character(GeneSet2[input$tab2table_rows_selected,1])
        genes1 <- as.matrix(gsea.df[which(gsea.df$Description==GS),"core_enrichment"])
        genes2 <- strsplit(genes1,"/")
        genes3 <- as.data.frame(genes2, col.names = "genes")
        gene_symbol <- genes3$genes
        ## convert expression matrix to numeric
        class(A) <- "numeric"
        ## Transforming data
        A <- A[,c(groupA,groupB)]
        exp.mat1 = log2(A + 1) # log
        exp.mat2 = apply(exp.mat1, 1, scale); # z score
        exp.mat3 = apply(exp.mat2, 1, rev); # transpose
        colnames(exp.mat3) = colnames(A) # set the column name
        exp.mat4 = exp.mat3[rowSums(is.na(exp.mat3[,])) == 0,] # remove NaN rows
        exp.mat5 = exp.mat4[rownames(exp.mat4) %in% gene_symbol,] # grab gene symbols
        # reassign data
        dataset <- exp.mat5
        ## generate color for pheatmap
        if (abs(min(dataset)) > abs(max(dataset))) {
          dataset[dataset < -abs(max(dataset))] = -abs(max(dataset))
        } else {
          dataset[dataset > abs(min(dataset))] = abs(min(dataset))
        }
        meta2 <- meta[which(meta[,1] %in% c(groupA,groupB)),]
        meta2 <- meta2[,c(colnames(meta2)[1],metacol)]
        meta2[,2] <- as.factor(meta2[,2])
        meta2 <- meta2[order(meta2[,2]),]
        rownames(meta2) <- meta2[,1]
        meta3 <- meta2[,-1,drop = F]
        dataset <- dataset[,rownames(meta2)]
        zscore_range = 10;
        minimum = -zscore_range;
        maximum = zscore_range;
        bk = c(seq(minimum,minimum/2, length=100), seq(minimum/2,maximum/2,length=100),seq(maximum/2,maximum,length=100))
        #Heatmap color
        color_choice <- input$ColorPalette_gseaHeat
        col_sets <- c("OrRd","PuBu","Greens","YlGnBu")
        if (color_choice == "original") {
          HeatMap_Colors <- c("dark blue","blue","white","red", "dark red")
          hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
        }
        else if (color_choice %in% col_sets) {
          HeatMap_Colors <- brewer.pal(n = 5, color_choice)
          hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
        }
        else if (color_choice == "Inferno") {
          hmcols <- inferno(500)
        }
        else if (color_choice == "Viridis") {
          hmcols <- viridis(500)
        }
        else if (color_choice == "Plasma") {
          hmcols <- plasma(500)
        }
        else if (color_choice == "OmniBlueRed") {
          hmcols<- colorRampPalette(c("#1984C5", "#22A7F0", "#63BFF0", "#A7D5ED", "#E2E2E2", "#E1A692", "#DE6E56", "#E14B31", "#C23728"))(length(bk)-1)
        }
        else if (color_choice == "LightBlueBlackRed") {
          hmcols<- colorRampPalette(c("#34C5FD","black","red"))(length(bk)-1)
        }
        else if (color_choice == "GreenBlackRed") {
          hmcols<- colorRampPalette(c("green","black","red"))(length(bk)-1)
        }
        pheatmap::pheatmap(dataset,            #data
                           cluster_cols = F,    #cluster columns - NO
                           cluster_row = F,     #cluster rows - YES
                           fontsize_col = input$heatmapFont1.c,   #column fontsize
                           fontsize_row = input$heatmapFont1.r,
                           show_rownames = T,  
                           show_colnames = T,
                           color=hmcols,
                           annotation_col = meta3)
      }
    }
    else if (input$tables == 5) {
      if (length(input$GStable.u_rows_selected) > 0){
        if (ncol(meta) > 2) {
          req(input$GSEAmetaCol)
          metacol <- input$GSEAmetaCol
        }
        else if (ncol(meta) == 2) {
          metacol <- colnames(meta)[2]
        }
        groupA <- meta[which(meta[,metacol] == input$comparisonA),1]
        groupB <- meta[which(meta[,metacol] == input$comparisonB),1]
        res <- datasetInput()
        gsea.df <- as.data.frame(res@result)
        GS <- as.character(user_gs_mirror()[input$GStable.u_rows_selected,1])
        genes1 <- as.matrix(gsea.df[which(gsea.df$Description==GS),"core_enrichment"])
        genes2 <- strsplit(genes1,"/")
        genes3 <- as.data.frame(genes2, col.names = "genes")
        gene_symbol <- genes3$genes
        ## convert expression matrix to numeric
        class(A) <- "numeric"
        ## Transforming data
        A <- A[,c(groupA,groupB)]
        exp.mat1 = log2(A + 1) # log
        exp.mat2 = apply(exp.mat1, 1, scale); # z score
        exp.mat3 = apply(exp.mat2, 1, rev); # transpose
        colnames(exp.mat3) = colnames(A) # set the column name
        exp.mat4 = exp.mat3[rowSums(is.na(exp.mat3[,])) == 0,] # remove NaN rows
        exp.mat5 = exp.mat4[rownames(exp.mat4) %in% gene_symbol,] # grab gene symbols
        # reassign data
        dataset <- exp.mat5
        ## generate color for pheatmap
        if (abs(min(dataset)) > abs(max(dataset))) {
          dataset[dataset < -abs(max(dataset))] = -abs(max(dataset))
        } else {
          dataset[dataset > abs(min(dataset))] = abs(min(dataset))
        }
        meta2 <- meta[which(meta[,1] %in% c(groupA,groupB)),]
        meta2 <- meta2[,c(colnames(meta2)[1],metacol)]
        meta2[,2] <- as.factor(meta2[,2])
        meta2 <- meta2[order(meta2[,2]),]
        rownames(meta2) <- meta2[,1]
        meta3 <- meta2[,-1,drop = F]
        dataset <- dataset[,rownames(meta2)]
        zscore_range = 10;
        minimum = -zscore_range;
        maximum = zscore_range;
        bk = c(seq(minimum,minimum/2, length=100), seq(minimum/2,maximum/2,length=100),seq(maximum/2,maximum,length=100))
        #Heatmap color
        color_choice <- input$ColorPalette_gseaHeat
        col_sets <- c("OrRd","PuBu","Greens","YlGnBu")
        if (color_choice == "original") {
          HeatMap_Colors <- c("dark blue","blue","white","red", "dark red")
          hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
        }
        else if (color_choice %in% col_sets) {
          HeatMap_Colors <- brewer.pal(n = 5, color_choice)
          hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
        }
        else if (color_choice == "Inferno") {
          hmcols <- inferno(500)
        }
        else if (color_choice == "Viridis") {
          hmcols <- viridis(500)
        }
        else if (color_choice == "Plasma") {
          hmcols <- plasma(500)
        }
        else if (color_choice == "OmniBlueRed") {
          hmcols<- colorRampPalette(c("#1984C5", "#22A7F0", "#63BFF0", "#A7D5ED", "#E2E2E2", "#E1A692", "#DE6E56", "#E14B31", "#C23728"))(length(bk)-1)
        }
        else if (color_choice == "LightBlueBlackRed") {
          hmcols<- colorRampPalette(c("#34C5FD","black","red"))(length(bk)-1)
        }
        else if (color_choice == "GreenBlackRed") {
          hmcols<- colorRampPalette(c("green","black","red"))(length(bk)-1)
        }
        pheatmap::pheatmap(dataset,            #data
                           cluster_cols = F,    #cluster columns - NO
                           cluster_row = F,     #cluster rows - YES
                           fontsize_col = input$heatmapFont1.c,   #column fontsize
                           fontsize_row = input$heatmapFont1.r,
                           show_rownames = T,  
                           show_colnames = T,
                           color=hmcols,
                           annotation_col = meta3)
      }
    }
  })
  
  #render MVG heatmap - per sample - 1
  output$heatmap1 <- renderPlot({
    
    heat <- MVGheatmap_react()
    heat
    
  })
  
  #render custom heatmap - per sample - 2
  output$heatmap2 <- renderPlot({
    
    heat <- CustomHeatmap_react()
    heat
    
  })
  
  # Render DEG Expression heatmap - #3
  output$heatmap3 <- renderPlot({
    
    heat <- DEGHeatmap_react()
    heat
    
    
  })
  
  # Render Average Expression heatmap - cutom
  output$avgheatmap1Cust <- renderPlot({
    
    heat <- AvgExprHeatmap_react()
    heat
    
    
  })
  
  output$ssgseaheatmap <- renderPlot({
    
    heat <- ssgseaheatmap_react()
    heat
    
  })
  
  output$ssgseaheatmap2 <- renderPlot({
    
    heat <- ssgseaheatmap2_react()
    heat
    
  })
  
  output$ssgseaheatmap3 <- renderPlot({
    
    heat <- ssgseaheatmap3_react()
    heat
    
  })
  
  #render MA plot
  output$MAPlot1 <- renderPlot({
    #title_font <- input$VolMATitleSize
    axis_font <- input$VolMAAxisSize
    top2 <- topgenereact()
    #add color categories based on FC and pval
    top2['threshold'] <- "none"
    top2[which(top2$logFC > abs(input$fc_cutoff)), "threshold"] <- "up"
    top2[which(top2$logFC < -abs(input$fc_cutoff)), "threshold"] <- "down"
    upRed <- "lightcoral"
    dnBlue <- "cadetblue3"
    mdGray <- "gray70"
    top2u <- top2[order(top2[,1], decreasing = TRUE),]
    top2d <- top2[order(top2[,1]),]
    top_hits_up <- top2u[head(which(top2u$logFC > abs(input$fc_cutoff)), n = input$top_x),]
    top_hits_dn <- top2d[head(which(top2d$logFC < -abs(input$fc_cutoff)), n = input$top_x),]
    #select genes to label based on user selection
    genesel.s <- NULL
    genesel.t <- NULL
    genesel.u <- NULL
    genesel.s <- unlist(strsplit(input$userGeneSelec2, " "))
    genesel.t <- unlist(strsplit(input$userGeneSelec2, "\t"))
    genesel.u <- input$userGeneSelec
    genesel.text <- c(genesel.s,genesel.t,genesel.u)
    top2_selec <- top2 %>%
      filter(GeneName %in% genesel.text)
    x <- ggplot(data = top2, aes(x = AveExpr, y = logFC)) +
      geom_point(size = 2, shape = 16) +
      theme_light(base_size = 16)
    #colors
    x <- x + aes(color = threshold) +
      scale_color_manual(values = c("up" = upRed,"down" = dnBlue, "none" = mdGray))
    x <- x + geom_hline(yintercept = c(-abs(input$fc_cutoff),abs(input$fc_cutoff)), linetype="dashed", color="gray20")
    x <- x + geom_text_repel(
      data =  top_hits_up,
      aes(label = rownames(top_hits_up)),
      color="gray20",
      size = 6,
      nudge_x = 0.2,
      nudge_y=0.2,
      box.padding = unit(0.9, "lines"),
      point.padding = unit(.3+4*0.1, "lines"),
      max.overlaps = 50)
    x <- x + geom_text_repel(
      data =  top_hits_dn,
      aes(label = rownames(top_hits_dn)),
      color="gray20",
      size = 6,
      nudge_x = 0.2,
      nudge_y=0.2,
      box.padding = unit(0.9, "lines"),
      point.padding = unit(.3+4*0.1, "lines"),
      max.overlaps = 50)
    x <- x + geom_text_repel(
      data =  top2_selec,
      aes(label = rownames(top2_selec)),
      color="gray20",
      size = 6,
      nudge_x = 0.2,
      nudge_y=0.2,
      box.padding = unit(0.9, "lines"),
      point.padding = unit(.3+4*0.1, "lines"),
      max.overlaps = 50)
    #coloring selected points
    x <- x + geom_point(data = top2_selec,
                        aes(x = AveExpr, y = logFC),
                        pch = 21,
                        color = "black",
                        size = 2)
    x <- x + geom_point(data = top_hits_dn,
                        aes(x = AveExpr, y = logFC),
                        pch = 21,
                        color = "black",
                        size = 2)
    x <- x + geom_point(data = top_hits_up,
                        aes(x = AveExpr, y = logFC),
                        pch = 21,
                        color = "black",
                        size = 2)
    x <- x + theme(legend.position="none")
    x <- x + labs(x = "Average Expression", y = "log2FC")
    x <- x + theme(axis.text = element_text(size=18))
    x <- x + theme(axis.title = element_text(size = axis_font))
    x
  })
  
  #render boxplot
  output$boxplot1 <- renderPlot({
    
    title_font <- input$boxplot2TitleSize
    axis_font <- input$boxplot2AxisSize
    if (ncol(meta) > 2) {
      metacol <- input$BoxPlotMetaColSelec
    }
    else if (ncol(meta) == 2) {
      metacol <- colnames(meta)[2]
    }
    
    if (length(input$GeneListTable_rows_selected) > 0){
      gene <- geneList[input$GeneListTable_rows_selected, 1]
      min <- min(log2(expr[gene,] + 1.0))
      max <- max(log2(expr[gene,] + 1.0))
      meta_temp <- meta
      rownames(meta_temp) <- meta[,1]
      meta_temp <- meta_temp[,metacol,drop = F]
      meta_temp[,metacol] <- as.factor(meta_temp[,metacol])
      #meta_temp <- meta_temp %>%
      #  select(Group)
      data = merge(t(expr[gene,]), meta_temp, by=0)
      colnames(data) = c("SampleName", "GeneExpr", "Cluster")
      ggplot(data, aes(x=Cluster, y=log2(GeneExpr + 1.0))) +
        geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=1) +
        geom_dotplot(binaxis='y', stackdir='center', dotsize=input$boxplotDot) +
        stat_compare_means(label = "p.signif") + ylim(min * 0.9, max * 1.3) +
        theme_bw() +
        labs(title= paste(gene, "Expression (log2)")) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1),
              axis.title = element_text(size = axis_font),
              plot.title = element_text(size = title_font))
    }
  })
  
  #render boxplot
  output$boxplot3 <- renderPlot({
    
    title_font <- input$boxplot1TitleSize
    axis_font <- input$boxplot1AxisSize
    if (ncol(meta) > 2) {
      metacol <- input$BoxPlot1MetaCol
    }
    else if (ncol(meta) == 2) {
      metacol <- colnames(meta)[2]
    }
    
    if (length(input$GeneListTable2_rows_selected) > 0){
      gene <- geneList[input$GeneListTable2_rows_selected, 1]
      min <- min(log2(expr[gene,] + 1.0))
      max <- max(log2(expr[gene,] + 1.0))
      meta_temp <- meta
      rownames(meta_temp) <- meta[,1]
      meta_temp <- meta_temp[,metacol,drop = F]
      meta_temp[,metacol] <- as.factor(meta_temp[,metacol])
      #meta_temp <- meta_temp %>%
      #  select(Group)
      data = merge(t(expr[gene,]), meta_temp, by=0)
      colnames(data) = c("SampleName", "GeneExpr", "Cluster")
      ggplot(data, aes(x=Cluster, y=log2(GeneExpr + 1.0))) +
        geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=1) +
        geom_dotplot(binaxis='y', stackdir='center', dotsize=input$boxplotDot2) +
        stat_compare_means(label = "p.signif") + ylim(min * 0.9, max * 1.3) +
        theme_bw() +
        labs(title= paste(gene, "Expression (log2)")) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1),
              axis.title = element_text(size = axis_font),
              plot.title = element_text(size = title_font))
    }
  })
  
  ####----Bar Plot----####
  barplot_react <- reactive({
    
    title_font <- input$barplot1TitleSize            # Title font size
    axis_font <- input$barplot1AxisSize              # Axis font size
    hjust_orient <- 1                                # Initial hjust
    axis_orient <- as.numeric(input$barxAxisOrient)  # X-axis label orientation
    if (axis_orient == 0) {                          # Adjust hjust if orientation is 0
      hjust_orient <- 0.5
    }
    bar_type <- input$errorbarplot                   # Error bar type
    logchoice <- input$log2barplot                   # Log expression data option
    colorin <- input$barplotColoCodes
    colorout <- input$barplotColoCodesOut
    bpylim <- input$barPlotYlim                      # Y-limit
    bpybreaks <- input$barplotYbreaks                # Y-axis breaks
    dotsizein <- input$barplotDotSize
    
    if (ncol(meta) > 2) {
      metacol <- input$BarPlot1MetaCol
    }
    else if (ncol(meta) == 2) {
      metacol <- colnames(meta)[2]
    }
    meta_temp <- meta
    rownames(meta_temp) <- meta[,1]
    meta_temp <- meta_temp[,metacol,drop = F]
    meta_temp[,metacol] <- as.factor(meta_temp[,metacol])
    
    if (length(input$GeneListTableBarPlot_rows_selected) > 0){
      
      gene <- geneList[input$GeneListTableBarPlot_rows_selected, 1]
      
      expr_gene <- merge(as.data.frame(t(expr[gene,])), meta_temp, by=0)
      colnames(expr_gene)[1] <- "SampleName"
      
      expr_gene2 <- merge(expr_gene,meta, all = T)
      plottitle <- paste(gene,"Average Gene Expression Across",metacol)
      genetitle <- paste(gene,"Average Expression")
      
      if (logchoice == T) {
        expr_gene2[,gene] <- log2(expr_gene2[,gene] + 1)
        plottitle <- paste(gene,"Average Gene Expression (Log2) Across",metacol)
        genetitle <- paste(gene,"Average Expression (Log2)")
      }
      
      colnames(expr_gene2)[c(1,2,3)] <- c("SampleName","Type","GeneName")
      expr_gene2 <- expr_gene2 %>%
        relocate(SampleName,GeneName,Type)
      
      se <- function(x) sd(x)/sqrt(length(x))
      expr_gene_stats <- expr_gene2 %>%
        group_by(Type) %>%
        summarise_at("GeneName",funs(mean,sd,se))
      
      y_min <- 0
      y_max <- round(max((expr_gene2[,2]) + 1))
      if (bpylim != "") {
        y_min <- as.numeric(gsub(" ","",strsplit(bpylim,",")[[1]][1]))
        y_max <- as.numeric(gsub(" ","",strsplit(bpylim,",")[[1]][2]))
      }
      
      if (input$barplotXaxOrder == "Descending"){
        barp <- ggplot(expr_gene_stats, aes(reorder(Type,-mean, na.rm = TRUE), mean, fill=Type))
      }
      if (input$barplotXaxOrder == "Ascending"){
        barp <- ggplot(expr_gene_stats, aes(reorder(Type, mean, na.rm = TRUE), mean, fill=Type))
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
      
      barp <- barp + labs(x = metacol,
                          y = genetitle,
                          title = plottitle)
      barp <- barp + theme(axis.text.x = element_text(angle = axis_orient, hjust = hjust_orient),
                           axis.text = element_text(size = axis_font),
                           axis.title = element_text(size = axis_font),
                           plot.title = element_text(size = title_font))
      if (bpylim == "") {
        if (is.na(bpybreaks)) {
          barp <- barp + scale_y_continuous(expand = c(0, 0))
        }
        else if (!is.na(bpybreaks)) {
          barp <- barp + scale_y_continuous(limits=c(0,y_max),expand = c(0, 0),breaks=seq(0,y_max,bpybreaks),oob = rescale_none)
        }
      }
      else if (bpylim != "") {
        if (is.na(bpybreaks)) {
          barp <- barp + scale_y_continuous(limits=c(y_min,y_max),expand = c(0, 0),oob = rescale_none)
        }
        else if (!is.na(bpybreaks)) {
          barp <- barp + scale_y_continuous(limits=c(y_min,y_max),expand = c(0, 0),breaks=seq(y_min,y_max,bpybreaks),oob = rescale_none)
        }
      }
      
      if (input$barplotsampledots == T) {
        barp <- barp + geom_dotplot(data = expr_gene2, aes(x=Type,y=GeneName),
                                    binaxis='y', stackdir='center',
                                    stackratio=1, dotsize=dotsizein, fill = "black")
      }
      
      if (colorin != "") {
        colorin <- strsplit(colorin," ")[[1]]
        if (length(colorin) == 1) {
          colorin <- rep(colorin,length(unique(meta[,metacol])))
        }
        barp <- barp + scale_fill_manual(values=colorin)
      }
      barp <- barp + coord_cartesian(clip = "off")
      
      
      
      barp
    }
  })
  
  output$barplot <- renderPlot({
    
    bpplot <- barplot_react()
    bpplot
    
  })
  
  #render up regulated pathway enrichment plot
  output$UpRegPathway1 <- renderPlot({
    title_font <- input$PathwayTitleSize
    axis_font <- input$PathwayAxisSize
    adjp <- input$pathpval
    FC <- input$pathFC
    if (ncol(meta) > 2) {
      metacol <- input$EnrichRPathMetaCol
    }
    else if (ncol(meta) == 2) {
      metacol <- colnames(meta)[2]
    }
    A <- meta[which(meta[,metacol] == input$comparisonA2.path),1]
    B <- meta[which(meta[,metacol] == input$comparisonB2.path),1]
    mat <- expr[,c(A,B)]
    mat <- log2(mat + 1.0)
    groupAOther <- factor(c(rep("A", length(A)), rep("B", length(B))))
    designA <- model.matrix(~0 + groupAOther)
    fit <- lmFit(mat, design = designA)
    contrast.matrix <- makeContrasts(groupAOtherA - groupAOtherB, levels = designA)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    options(digits = 4)
    top1 <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
    genes <- rownames(top1)[which(top1$P.Value < adjp & top1$logFC > FC)]
    dbs <- listEnrichrDbs() 
    enrichRLive <- TRUE 
    if (is.null(dbs)) { 
      enrichRLive <- FALSE 
    }
    dbs <- input$SelectedPathway
    enriched <- enrichr(genes, dbs) # Plot top 20 GO-BP results ordered by P-value 
    plotEnrich(enriched[[1]], showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value") +
      theme(axis.title = element_text(size = axis_font),
            plot.title = element_text(size = title_font))
  })
  
  #render up regulated pathway enrichment plot
  output$DnRegPathway1 <- renderPlot({
    title_font <- input$PathwayTitleSize
    axis_font <- input$PathwayAxisSize
    adjp <- input$pathpval
    FC <- input$pathFC
    if (ncol(meta) > 2) {
      metacol <- input$EnrichRPathMetaCol
    }
    else if (ncol(meta) == 2) {
      metacol <- colnames(meta)[2]
    }
    A <- meta[which(meta[,metacol] == input$comparisonA2.path),1]
    B <- meta[which(meta[,metacol] == input$comparisonB2.path),1]
    mat <- expr[,c(A,B)]
    mat <- log2(mat + 1.0)
    groupAOther <- factor(c(rep("A", length(A)), rep("B", length(B))))
    designA <- model.matrix(~0 + groupAOther)
    fit <- lmFit(mat, design = designA)
    contrast.matrix <- makeContrasts(groupAOtherA - groupAOtherB, levels = designA)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    options(digits = 4)
    top1 <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
    genes <- rownames(top1)[which(top1$P.Value < adjp & top1$logFC > FC)]
    dbs <- listEnrichrDbs() 
    enrichRLive <- TRUE 
    if (is.null(dbs)) { 
      enrichRLive <- FALSE 
    }
    dbs <- input$SelectedPathway
    enriched <- enrichr(genes, dbs) # Plot top 20 GO-BP results ordered by P-value 
    plotEnrich(enriched[[1]], showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value") +
      theme(axis.title = element_text(size = axis_font),
            plot.title = element_text(size = title_font))
  })
  
  Volcano3_react <- reactive({
    #title_font <- input$VolMATitleSize
    axis_font <- input$VolMAAxisSize
    tick_font <- input$VolMATickSize
    anno_font <- input$VolMAAnnoSize
    anno_face <- input$VolMAAnnoFont
    top2 <- topgenereact()
    #add color categories based on FC and pval
    top2['threshold'] <- "none"
    top2[which(top2$logFC > abs(input$fc_cutoff) & top2$P.Value < input$p_cutoff), "threshold"] <- "up"
    top2[which(top2$logFC < -abs(input$fc_cutoff) & top2$P.Value < input$p_cutoff), "threshold"] <- "down"
    upRed <- "lightcoral"
    dnBlue <- "cadetblue3"
    mdGray <- "gray70"
    #select number of top hits based on input
    top_hits_up <- top2[head(which(top2$logFC > abs(input$fc_cutoff) & top2$P.Value < input$p_cutoff), n = input$top_x),]
    top_hits_dn <- top2[head(which(top2$logFC < -abs(input$fc_cutoff) & top2$P.Value < input$p_cutoff), n = input$top_x),]
    #select genes to label based on user selection
    genesel.s <- NULL
    genesel.t <- NULL
    genesel.u <- NULL
    genesel.s <- unlist(strsplit(input$userGeneSelec2, " "))
    genesel.t <- unlist(strsplit(input$userGeneSelec2, "\t"))
    genesel.u <- input$userGeneSelec
    genesel.text <- c(genesel.s,genesel.t,genesel.u)
    top2_selec <- top2 %>%
      filter(GeneName %in% genesel.text)
    #create plot
    x <- ggplot(data = top2, aes(x = logFC, y = -log10(P.Value))) +
      geom_point(size = 2, shape = 16) +
      theme_light(base_size = 16)
    #colors
    x <- x + aes(color = threshold) +
      scale_color_manual(values = c("up" = upRed,"down" = dnBlue, "none" = mdGray))
    #FC and pval lines
    x <- x + geom_vline(xintercept = c(-abs(input$fc_cutoff),abs(input$fc_cutoff)), linetype="dashed", color="gray20")
    x <- x + geom_hline(yintercept = -log10(input$p_cutoff), linetype="dashed", color="gray20")
    #label top hits if needed
    if (input$top_x > 0) { 
      x <- x + geom_text_repel(
        data =  top_hits_up,
        #aes(label = rownames(top_hits_up), fontface = anno_face),
        aes(label = rownames(top_hits_up)),
        size = anno_font,
        color="gray20",
        min.segment.length = 0,
        nudge_x = 0.2,
        nudge_y=0.2,
        box.padding = unit(0.9, "lines"),
        point.padding = unit(.3+4*0.1, "lines"),
        max.overlaps = Inf)
      x <- x + geom_text_repel(
        data =  top_hits_dn,
        #aes(label = rownames(top_hits_dn), fontface = anno_face),
        aes(label = rownames(top_hits_dn)),
        size = anno_font,
        color="gray20",
        min.segment.length = 0,
        nudge_x = 0.2,
        nudge_y=0.2,
        box.padding = unit(0.9, "lines"),
        point.padding = unit(.3+4*0.1, "lines"),
        max.overlaps = Inf)
    }
    x <- x + geom_text_repel(
      data =  top2_selec,
      #aes(label = rownames(top2_selec), fontface = anno_face),
      aes(label = rownames(top2_selec)),
      size = anno_font,
      color="gray20",
      min.segment.length = 0,
      nudge_x = 0.2,
      nudge_y=0.2,
      box.padding = unit(0.9, "lines"),
      point.padding = unit(.3+4*0.1, "lines"),
      max.overlaps = Inf)
    #coloring selected points
    x <- x + geom_point(data = top2_selec,
                        aes(x = logFC, y = -log10(P.Value)),
                        pch = 21,
                        color = "black",
                        size = 2)
    x <- x + geom_point(data = top_hits_dn,
                        aes(x = logFC, y = -log10(P.Value)),
                        pch = 21,
                        color = "black",
                        size = 2)
    x <- x + geom_point(data = top_hits_up,
                        aes(x = logFC, y = -log10(P.Value)),
                        pch = 21,
                        color = "black",
                        size = 2)
    #axis parameters
    x <- x + theme(legend.position="none")
    x <- x + theme(axis.text = element_text(size=tick_font),
                   axis.title = element_text(size = axis_font))
    #x <- x + theme(axis.title = element_text(size = axis_font))
    x <- x + labs(x = "log2FC", y = "-log10(P.Value)")
    x
    
    
  })
  
  #render volcano plot
  output$Volcano3 <- renderPlot({
    plot <- Volcano3_react()
    plot
  })
  
  #render ssGSEA boxplot
  output$boxplot2 <- renderPlot({
    
    title_font <- input$gseaBoxTitleSize
    axis_font <- input$gseaBoxAxisSize
    if (ncol(meta) > 2) {
      metacol <- input$GSEAmetaCol
    }
    else if (ncol(meta) == 2) {
      metacol <- colnames(meta)[2]
    }
    
    #if tab 2 selected
    if (input$tables == 3) {
      if (length(input$tab2table_rows_selected) > 0){
        ssgsea <- ssGSEAfunc()
        ssgsea2 <- as.data.frame(t(ssgsea))
        samporder <- meta[,1]
        ssgsea3 <- as.data.frame(ssgsea2[samporder,])
        colnames(ssgsea3)[1] <- colnames(ssgsea2[1])
        rownames(ssgsea3) <- samporder
        ssgsea4 <- ssgsea3 %>% 
          mutate(type = case_when(
            rownames(ssgsea3) == meta[,1] ~ meta[,metacol],
          ))
        ssgsea4 <- ssgsea4 %>%
          relocate(type)
        ssgsea4$type <- as.factor(ssgsea4$type)
        colnames(ssgsea4)[which(colnames(ssgsea4) == "type")] <- metacol
        if (input$boxplotcompare == "none") {
          ggplot(ssgsea4, aes(factor(ssgsea4[,metacol]), ssgsea4[,2], fill = ssgsea4[,metacol])) +
            geom_boxplot(width = 0.5, lwd = 1, fill = "white") +
            geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = input$boxplotDotss) +
            labs(x = "Group", y = paste(colnames(ssgsea4)[2], " Score", sep = ""),
                 title = paste("ssGSEA Expression Signature: ", colnames(ssgsea4)[2],sep = ""),
                 fill = metacol) +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1),
                  axis.text = element_text(size = 12),
                  axis.title = element_text(size = axis_font),
                  plot.title = element_text(size = title_font))
        }
        else {
          ggplot(ssgsea4, aes(factor(ssgsea4[,metacol]), ssgsea4[,2], fill = ssgsea4[,metacol])) +
            geom_boxplot(width = 0.5, lwd = 1, fill = "white") +
            geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = input$boxplotDotss) +
            labs(x = "Group", y = paste(colnames(ssgsea4)[2], " Score", sep = ""),
                 title = paste("ssGSEA Expression Signature: ", colnames(ssgsea4)[2],sep = ""),
                 fill = metacol) +
            theme_bw() +
            stat_compare_means(method = input$boxplotcompare) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1),
                  axis.text = element_text(size = 12),
                  axis.title = element_text(size = axis_font),
                  plot.title = element_text(size = title_font))
        }
      }
    }
    #if MSigDB tab selected
    else if (input$tables == 1) {
      if (length(input$msigdbTable_rows_selected) > 0){
        ssgsea <- ssGSEAfunc()
        ssgsea2 <- as.data.frame(t(ssgsea))
        samporder <- meta[,1]
        ssgsea3 <- as.data.frame(ssgsea2[samporder,])
        colnames(ssgsea3)[1] <- colnames(ssgsea2[1])
        rownames(ssgsea3) <- samporder
        ssgsea4 <- ssgsea3 %>% 
          mutate(type = case_when(
            rownames(ssgsea3) == meta[,1] ~ meta[,metacol],
          ))
        ssgsea4 <- ssgsea4 %>%
          relocate(type)
        ssgsea4$type <- as.factor(ssgsea4$type)
        colnames(ssgsea4)[which(colnames(ssgsea4) == "type")] <- metacol
        if (input$boxplotcompare == "none") {
          ggplot(ssgsea4, aes(factor(ssgsea4[,metacol]), ssgsea4[,2], fill = ssgsea4[,metacol])) +
            geom_boxplot(width = 0.5, lwd = 1, fill = "white") +
            geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = input$boxplotDotss) +
            labs(x = "Group", y = paste(colnames(ssgsea4)[2], " Score", sep = ""),
                 title = paste("ssGSEA Expression Signature: ",colnames(ssgsea4)[2],sep = ""),
                 fill = metacol) +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1),
                  axis.text = element_text(size = 12),
                  axis.title = element_text(size = axis_font),
                  plot.title = element_text(size = title_font))
        }
        else {
          ggplot(ssgsea4, aes(factor(ssgsea4[,metacol]), ssgsea4[,2], fill = ssgsea4[,metacol])) +
            geom_boxplot(width = 0.5, lwd = 1, fill = "white") +
            geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = input$boxplotDotss) +
            labs(x = "Group", y = paste(colnames(ssgsea4)[2], " Score", sep = ""),
                 title = paste("ssGSEA Expression Signature: ",colnames(ssgsea4)[2],sep = ""),
                 fill = metacol) +
            theme_bw() +
            stat_compare_means(method = input$boxplotcompare) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1),
                  axis.text = element_text(size = 12),
                  axis.title = element_text(size = axis_font),
                  plot.title = element_text(size = title_font))
        }
      }
    }
    #if user generated GS table selected
    else if (input$tables == 5) {
      if (length(input$GStable.u_rows_selected) > 0){
        ssgsea <- ssGSEAfunc()
        ssgsea2 <- as.data.frame(t(ssgsea))
        samporder <- meta[,1]
        ssgsea3 <- as.data.frame(ssgsea2[samporder,])
        colnames(ssgsea3)[1] <- colnames(ssgsea2[1])
        rownames(ssgsea3) <- samporder
        ssgsea4 <- ssgsea3 %>% 
          mutate(type = case_when(
            rownames(ssgsea3) == meta[,1] ~ meta[,metacol],
          ))
        ssgsea4 <- ssgsea4 %>%
          relocate(type)
        ssgsea4$type <- as.factor(ssgsea4$type)
        colnames(ssgsea4)[which(colnames(ssgsea4) == "type")] <- metacol
        if (input$boxplotcompare == "none") {
          ggplot(ssgsea4, aes(factor(ssgsea4[,metacol]), ssgsea4[,2], fill = ssgsea4[,metacol])) +
            geom_boxplot(width = 0.5, lwd = 1, fill = "white") +
            geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = input$boxplotDotss) +
            labs(x = "Group", y = paste(colnames(ssgsea4)[2], " Score", sep = ""),
                 title = paste("ssGSEA Expression Signature: ",colnames(ssgsea4)[2],sep = ""),
                 fill = metacol) +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1),
                  axis.text = element_text(size = 12),
                  axis.title = element_text(size = axis_font),
                  plot.title = element_text(size = title_font))
        }
        else {
          ggplot(ssgsea4, aes(factor(ssgsea4[,metacol]), ssgsea4[,2], fill = ssgsea4[,metacol])) +
            geom_boxplot(width = 0.5, lwd = 1, fill = "white") +
            geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = input$boxplotDotss) +
            labs(x = "Group", y = paste(colnames(ssgsea4)[2], " Score", sep = ""),
                 title = paste("ssGSEA Expression Signature: ",colnames(ssgsea4)[2],sep = ""),
                 fill = metacol) +
            theme_bw() +
            stat_compare_means(method = input$boxplotcompare) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1),
                  axis.text = element_text(size = 12),
                  axis.title = element_text(size = axis_font),
                  plot.title = element_text(size = title_font))
        }
      }
    }
  })
  
  
  
  # Render Average Gene Expression comparison plot
  output$AvggeneScatter2 <- renderPlot({
    
    A_choice <- input$comparisonA2.avg
    B_choice <- input$comparisonB2.avg
    log_choice <- input$AvgExpLogFC
    title_font <- input$AvgExprScatterTitleSize
    axis_font <- input$AvgExprScatterAxisSize
    
    #select genes to label based on user selection
    genesel.s <- NULL
    genesel.t <- NULL
    genesel.u <- NULL
    genesel.s <- unlist(strsplit(input$gsSelection2, " "))
    genesel.t <- unlist(strsplit(input$gsSelection2, "\t"))
    genesel.u <- input$scatterGeneSelec
    genesel.text <- c(genesel.s,genesel.t,genesel.u)
    
    AvgExpr_Table <- AvgExprReact()
    
    logFC_text <- ""
    
    # Log if user designates
    if (log_choice == TRUE) {
      logFC_text <- " (log2 +1)"
    }
    
    # Add group for above or below abline
    AvgExpr_Table <- AvgExpr_Table %>%
      mutate(ColorCol = case_when(
        AvgExpr_Table$AvgExpression_GroupB > AvgExpr_Table$AvgExpression_GroupA ~ paste(B_choice," Expr > ",A_choice,' Expr', sep = ''),
        AvgExpr_Table$AvgExpression_GroupA > AvgExpr_Table$AvgExpression_GroupB ~ paste(A_choice," Expr > ",B_choice,' Expr', sep = '')
      ))
    
    AvgExpr_Table_selec <- AvgExpr_Table %>%
      filter(GeneSymbol %in% genesel.text)
    
    # plot
    p <- ggplot(AvgExpr_Table, aes(x = AvgExpression_GroupA, y = AvgExpression_GroupB,
                                   color = ColorCol)) +
      geom_point() +
      theme_minimal() +
      theme(legend.position="none",
            axis.title = element_text(size = axis_font),
            plot.title = element_text(size = title_font)) +
      xlab(paste("Average Expression: ",A_choice,logFC_text,sep = ""))+
      ylab(paste("Average Expression: ",B_choice,logFC_text,sep = "")) +
      labs(title = paste("Average Expression",logFC_text,": ",A_choice," vs. ",B_choice,sep = "")) +
      scale_color_manual(values = c("cadetblue3","lightcoral"))
    
    p <- p + geom_text_repel(
      data =  AvgExpr_Table_selec,
      aes(label = GeneSymbol),
      size = 6,
      color="gray20",
      nudge_x = 0.2,
      nudge_y=0.2,
      box.padding = unit(0.9, "lines"),
      point.padding = unit(.3+4*0.1, "lines"),
      max.overlaps = 50)
    
    #coloring selected points
    p <- p + geom_point(data = AvgExpr_Table_selec,
                        aes(x = AvgExpression_GroupA, y = AvgExpression_GroupB),
                        pch = 21,
                        color = "black",
                        size = 2)
    
    p
    
  })
  
  
  ####----Download Handlers----####
  
  output$dnldPlotSVG_exprBar <- downloadHandler(
    filename = function() {
      paste(gsub(" ","",ProjectName),"_AvgExpression_Barplot.svg", sep = '')
    },
    content = function(file) {
      bpplot <- barplot_react()
      ggsave(file,bpplot, width = 10, height = 8)
    }
  )
  
  output$dnldPlotPDF_exprBar <- downloadHandler(
    filename = function() {
      paste(gsub(" ","",ProjectName),"_AvgExpression_Barplot.pdf", sep = '')
    },
    content = function(file) {
      bpplot <- barplot_react()
      ggsave(file,bpplot, width = 10, height = 8)
    }
  )
  
  output$dnldPlotSVG_ssgseaHeat <- downloadHandler(
    filename = function() {
      paste(gsub(" ","",ProjectName),"_ssGSEA_Heatmap.svg", sep = '')
    },
    content = function(file) {
      
      heat <- ssgseaheatmap_react()
      
      ggsave(file,heat, width = 15, height = 20)
    }
  )
  
  output$dnldPlotPDF_ssgseaHeat <- downloadHandler(
    filename = function() {
      paste(gsub(" ","",ProjectName),"_ssGSEA_Heatmap.pdf", sep = '')
    },
    content = function(file) {
      
      heat <- ssgseaheatmap_react()
      ggsave(file,heat, width = 15, height = 20)
    }
  )
  
  output$dnldPlotSVG_ssgseaHeat2 <- downloadHandler(
    filename = function() {
      paste(gsub(" ","",ProjectName),"_ssGSEA_Heatmap.svg", sep = '')
    },
    content = function(file) {
      
      heat <- ssgseaheatmap2_react()
      
      ggsave(file,heat, width = 15, height = 20)
    }
  )
  
  output$dnldPlotPDF_ssgseaHeat2 <- downloadHandler(
    filename = function() {
      paste(gsub(" ","",ProjectName),"_ssGSEA_Heatmap.pdf", sep = '')
    },
    content = function(file) {
      
      heat <- ssgseaheatmap2_react()
      ggsave(file,heat, width = 15, height = 20)
    }
  )
  
  output$dnldPlotSVG_ssgseaHeat3 <- downloadHandler(
    filename = function() {
      paste(gsub(" ","",ProjectName),"_ssGSEA_Heatmap.svg", sep = '')
    },
    content = function(file) {
      
      heat <- ssgseaheatmap3_react()
      
      ggsave(file,heat, width = 15, height = 20)
    }
  )
  
  output$dnldPlotPDF_ssgseaHeat3 <- downloadHandler(
    filename = function() {
      paste(gsub(" ","",ProjectName),"_ssGSEA_Heatmap.pdf", sep = '')
    },
    content = function(file) {
      
      heat <- ssgseaheatmap3_react()
      ggsave(file,heat, width = 15, height = 20)
    }
  )
  
  
  output$dnldPlotSVG_heat1 <- downloadHandler(
    filename = function() {
      paste(gsub(" ","",ProjectName),"_MVG_Heatmap.svg", sep = '')
    },
    content = function(file) {
      
      heat <- MVGheatmap_react()
      
      ggsave(file,heat, width = 10, height = 20)
    }
  )
  
  output$dnldPlotSVG_heat2 <- downloadHandler(
    filename = function() {
      paste(gsub(" ","",ProjectName),"_Custom_Heatmap.svg", sep = '')
    },
    content = function(file) {
      
      heat <- CustomHeatmap_react()
      ggsave(file,heat, width = 10, height = 20)
    }
  )
  
  output$dnldPlotSVG_heat3 <- downloadHandler(
    filename = function() {
      paste(gsub(" ","",ProjectName),"_DEGExpr_Heatmap.svg", sep = '')
    },
    content = function(file) {
      
      heat <- DEGHeatmap_react()
      
      ggsave(file,heat, width = 10, height = 20)
    }
  )
  
  output$dnldPlotSVG_heat5 <- downloadHandler(
    filename = function() {
      paste(gsub(" ","",ProjectName),"_CustomAvgExpr_Heatmap.svg", sep = '')
    },
    content = function(file) {
      
      heat <- AvgExprHeatmap_react()
      
      ggsave(file,heat, width = 10, height = 20)
      
    }
  )
  
  output$dnldPlotSVG_scatter <- downloadHandler(
    filename = function() {
      paste(gsub(" ","",ProjectName),"_2GeneExpression_ScatterPlot.svg", sep = '')
    },
    content = function(file) {
      #log if user designates
      if (input$logask == TRUE) {
        expr <- log2(expr + 1)
      }
      #transpose
      expr_t <- as.data.frame(t(expr))
      #reorder rowname to match meta for merging
      samporder <- meta[,1]
      expr_t2 <- as.data.frame(expr_t[samporder,])
      #add type
      expr_t3 <- expr_t2 %>% 
        mutate(type = case_when(
          rownames(expr_t2) == meta[,1] ~ meta[,2],
        ))
      expr_t3 <- expr_t3 %>%
        relocate(type)
      #user gene input
      gene1.u <- input$scatterG1
      gene2.u <- input$scatterG2
      #get columns and info based off user input
      gene1 <- expr_t3[,gene1.u]
      gene2 <- expr_t3[,gene2.u]
      if (input$logask == TRUE) {
        gene1.n <- paste(colnames(expr_t3)[which(colnames(expr_t3) == gene1.u)], " Expression (log2 +1)", sep = "")
        gene2.n <- paste(colnames(expr_t3)[which(colnames(expr_t3) == gene2.u)], " Expression (log2 +1)", sep = "")
      }
      else if (input$logask == FALSE) {
        gene2.n <- paste(colnames(expr_t3)[which(colnames(expr_t3) == gene2.u)], " Expression", sep = "")
        gene1.n <- paste(colnames(expr_t3)[which(colnames(expr_t3) == gene1.u)], " Expression", sep = "")
      }
      #plot
      p <- ggplot(expr_t3, aes(x = gene1, y = gene2,
                               color = type)) +
        geom_point() +
        theme_minimal() +
        labs(x = gene1.n, y = gene2.n,
             title = paste(gene1.n, " vs. ", gene2.n, sep = ''))
      ggsave(file,p, width = 10, height = 8)
    }
  )
  
  output$dnldPlotSVG_exprBox <- downloadHandler(
    filename = function() {
      paste(gsub(" ","",ProjectName),"_Expression_Boxplot.svg", sep = '')
    },
    content = function(file) {
      if (length(input$GeneListTable2_rows_selected) > 0){
        gene <- geneList[input$GeneListTable2_rows_selected, 1]
        min <- min(log2(expr[gene,] + 1.0))
        max <- max(log2(expr[gene,] + 1.0))
        meta_temp <- meta
        rownames(meta_temp) <- meta[,1]
        meta_temp <- meta_temp %>%
          select(Group)
        data = merge(t(expr[gene,]), meta_temp, by=0)
        colnames(data) = c("SampleName", "GeneExpr", "Cluster")
        x <- ggplot(data, aes(x=Cluster, y=log2(GeneExpr + 1.0))) +
          geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=1) +
          geom_dotplot(binaxis='y', stackdir='center', dotsize=input$boxplotDot2) +
          stat_compare_means(label = "p.signif") + ylim(min * 0.9, max * 1.3) +
          theme_bw() +
          labs(title= paste(gene, "Expression (log2)")) +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1),
                text = element_text(size = input$boxplot.font2))
        ggsave(file,x, width = 10, height = 8)
      }
    }
  )
  
  output$dnldPlotSVG_vol <- downloadHandler(
    filename = function() {
      paste(gsub(" ","",ProjectName),"_VolcanoPlot.svg", sep = '')
    },
    content = function(file) {
      
      x <- Volcano3_react()
      PlotH <- input$VolMAdnldHeight
      PlotW <- input$VolMAdnldWidth
      PlotU <- input$VolMAdnldSizeUnits
      ggsave(file,x, width = PlotW, height = PlotH, units = PlotU)
    }
  )
  
  output$dnldPlotSVG_MA <- downloadHandler(
    filename = function() {
      paste(gsub(" ","",ProjectName),"_MA_Plot.svg", sep = '')
    },
    content = function(file) {
      top2 <- topgenereact()
      #add color categories based on FC and pval
      top2['threshold'] <- "none"
      top2[which(top2$logFC > abs(input$fc_cutoff)), "threshold"] <- "up"
      top2[which(top2$logFC < -abs(input$fc_cutoff)), "threshold"] <- "down"
      upRed <- "lightcoral"
      dnBlue <- "cadetblue3"
      mdGray <- "gray70"
      top2u <- top2[order(top2[,1], decreasing = TRUE),]
      top2d <- top2[order(top2[,1]),]
      top_hits_up <- top2u[head(which(top2u$logFC > abs(input$fc_cutoff)), n = input$top_x),]
      top_hits_dn <- top2d[head(which(top2d$logFC < -abs(input$fc_cutoff)), n = input$top_x),]
      #select genes to label based on user selection
      genesel.s <- NULL
      genesel.t <- NULL
      genesel.u <- NULL
      genesel.s <- unlist(strsplit(input$userGeneSelec2, " "))
      genesel.t <- unlist(strsplit(input$userGeneSelec2, "\t"))
      genesel.u <- input$userGeneSelec
      genesel.text <- c(genesel.s,genesel.t,genesel.u)
      top2_selec <- top2 %>%
        filter(GeneName %in% genesel.text)
      x <- ggplot(data = top2, aes(x = AveExpr, y = logFC)) +
        geom_point(size = 2, shape = 16) +
        theme_light(base_size = 16)
      #colors
      x <- x + aes(color = threshold) +
        scale_color_manual(values = c("up" = upRed,"down" = dnBlue, "none" = mdGray))
      x <- x + geom_hline(yintercept = c(-abs(input$fc_cutoff),abs(input$fc_cutoff)), linetype="dashed", color="gray20")
      x <- x + geom_text_repel(
        data =  top_hits_up,
        aes(label = rownames(top_hits_up)),
        color="gray20",
        size = 6,
        nudge_x = 0.2,
        nudge_y=0.2,
        box.padding = unit(0.9, "lines"),
        point.padding = unit(.3+4*0.1, "lines"),
        max.overlaps = 50)
      x <- x + geom_text_repel(
        data =  top_hits_dn,
        aes(label = rownames(top_hits_dn)),
        color="gray20",
        size = 6,
        nudge_x = 0.2,
        nudge_y=0.2,
        box.padding = unit(0.9, "lines"),
        point.padding = unit(.3+4*0.1, "lines"),
        max.overlaps = 50)
      x <- x + geom_text_repel(
        data =  top2_selec,
        aes(label = rownames(top2_selec)),
        color="gray20",
        size = 6,
        nudge_x = 0.2,
        nudge_y=0.2,
        box.padding = unit(0.9, "lines"),
        point.padding = unit(.3+4*0.1, "lines"),
        max.overlaps = 50)
      #coloring selected points
      x <- x + geom_point(data = top2_selec,
                          aes(x = AveExpr, y = logFC),
                          pch = 21,
                          color = "black",
                          size = 2)
      x <- x + geom_point(data = top_hits_dn,
                          aes(x = AveExpr, y = logFC),
                          pch = 21,
                          color = "black",
                          size = 2)
      x <- x + geom_point(data = top_hits_up,
                          aes(x = AveExpr, y = logFC),
                          pch = 21,
                          color = "black",
                          size = 2)
      x <- x + theme(legend.position="none")
      x <- x + labs(x = "Average Expression", y = "log2FC")
      x <- x + theme(axis.text = element_text(size=18))
      x <- x + theme(axis.title = element_text(size=24))
      ggsave(file,x, width = 10, height = 8)
    }
  )
  
  output$dnldPlotSVG_exprBox2 <- downloadHandler(
    filename = function() {
      paste(gsub(" ","",ProjectName),"_ExpressionBoxplot.svg", sep = '')
    },
    content = function(file) {
      if (length(input$GeneListTable_rows_selected) > 0){
        gene <- geneList[input$GeneListTable_rows_selected, 1]
        min <- min(log2(expr[gene,] + 1.0))
        max <- max(log2(expr[gene,] + 1.0))
        meta_temp <- meta
        rownames(meta_temp) <- meta[,1]
        meta_temp <- meta_temp %>%
          select(Group)
        data = merge(t(expr[gene,]), meta_temp, by=0)
        colnames(data) = c("SampleName", "GeneExpr", "Cluster")
        x <- ggplot(data, aes(x=Cluster, y=log2(GeneExpr + 1.0))) +
          geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=1) +
          geom_dotplot(binaxis='y', stackdir='center', dotsize=input$boxplotDot) +
          stat_compare_means(label = "p.signif") + ylim(min * 0.9, max * 1.3) +
          theme_bw() +
          labs(title= paste(gene, "Expression (log2)")) +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1),
                text = element_text(size = input$boxplot.font))
        ggsave(file,x, width = 10, height = 8)
      }
    }
  )
  
  output$dnldPlotSVG_AvgScatter <- downloadHandler(
    filename = function() {
      paste(gsub(" ","",ProjectName),"_AvgExpression_ScatterPlot.svg", sep = '')
    },
    content = function(file) {
      A_choice <- input$comparisonA2.avg
      B_choice <- input$comparisonB2.avg
      log_choice <- input$AvgExpLogFC
      
      #select genes to label based on user selection
      genesel.s <- NULL
      genesel.t <- NULL
      genesel.u <- NULL
      genesel.s <- unlist(strsplit(input$gsSelection2, " "))
      genesel.t <- unlist(strsplit(input$gsSelection2, "\t"))
      genesel.u <- input$scatterGeneSelec
      genesel.text <- c(genesel.s,genesel.t,genesel.u)
      
      if (ncol(meta) > 2){
        metacol <- input$AvgExprScatterMetaCol
      }
      else if (ncol(meta) == 2){
        metacol <- colnames(meta)[2]
      }
      
      # Get A and B sample names
      A <- meta[which(meta[,metacol] == A_choice),1]
      B <- meta[which(meta[,metacol] == B_choice),1]
      
      # Make A and B expression Matrices
      mat_A <- expr[,A]
      mat_B <- expr[,B]
      
      logFC_text <- ""
      
      # Log if user designates
      if (log_choice == TRUE) {
        mat_A <- log2(mat_A + 1)
        mat_B <- log2(mat_B + 1)
        logFC_text <- " (log2 +1)"
      }
      
      # Get avg expression of each gene
      mat_A$AvgExpression_GroupA <- rowMeans(mat_A)
      mat_B$AvgExpression_GroupB <- rowMeans(mat_B)
      
      # Add Gene Names as column
      mat_A$GeneSymbol <- rownames(mat_A)
      mat_B$GeneSymbol <- rownames(mat_B)
      
      # Merge average columns
      AvgExpr_Table <- merge(mat_A[,c("AvgExpression_GroupA","GeneSymbol")],
                             mat_B[,c("AvgExpression_GroupB","GeneSymbol")],
                             by = "GeneSymbol",
                             all = T)
      
      # Add group for above or below abline
      AvgExpr_Table <- AvgExpr_Table %>%
        mutate(ColorCol = case_when(
          AvgExpr_Table$AvgExpression_GroupB > AvgExpr_Table$AvgExpression_GroupA ~ paste(B_choice," Expr > ",A_choice,' Expr', sep = ''),
          AvgExpr_Table$AvgExpression_GroupA > AvgExpr_Table$AvgExpression_GroupB ~ paste(A_choice," Expr > ",B_choice,' Expr', sep = '')
        ))
      
      AvgExpr_Table_selec <- AvgExpr_Table %>%
        filter(GeneSymbol %in% genesel.text)
      
      # plot
      p <- ggplot(AvgExpr_Table, aes(x = AvgExpression_GroupA, y = AvgExpression_GroupB,
                                     color = ColorCol)) +
        geom_point() +
        theme_minimal() +
        theme(legend.position="none") +
        xlab(paste("Average Expression: ",A_choice,logFC_text,sep = ""))+
        ylab(paste("Average Expression: ",B_choice,logFC_text,sep = "")) +
        labs(title = paste("Average Expression",logFC_text,": ",A_choice," vs. ",B_choice,sep = "")) +
        scale_color_manual(values = c("cadetblue3","lightcoral"))
      
      p <- p + geom_text_repel(
        data =  AvgExpr_Table_selec,
        aes(label = GeneSymbol),
        size = 6,
        color="gray20",
        nudge_x = 0.2,
        nudge_y=0.2,
        box.padding = unit(0.9, "lines"),
        point.padding = unit(.3+4*0.1, "lines"),
        max.overlaps = 50)
      
      #coloring selected points
      p <- p + geom_point(data = AvgExpr_Table_selec,
                          aes(x = AvgExpression_GroupA, y = AvgExpression_GroupB),
                          pch = 21,
                          color = "black",
                          size = 2)
      
      ggsave(file,p, width = 10, height = 8)
    }
  )
  
  output$dnldPlotSVG_upPath <- downloadHandler(
    filename = function() {
      paste(gsub(" ","",ProjectName),"_UpRegulated_Pathways.svg", sep = '')
    },
    content = function(file) {
      adjp <- input$pathpval
      FC <- input$pathFC
      if (ncol(meta) > 2) {
        metacol <- input$EnrichRPathMetaCol
      }
      else if (ncol(meta) == 2) {
        metacol <- colnames(meta)[2]
      }
      A <- meta[which(meta[,metacol] == input$comparisonA2.path),1]
      B <- meta[which(meta[,metacol] == input$comparisonB2.path),1]
      mat <- expr[,c(A,B)]
      mat <- log2(mat + 1.0)
      groupAOther <- factor(c(rep("A", length(A)), rep("B", length(B))))
      designA <- model.matrix(~0 + groupAOther)
      fit <- lmFit(mat, design = designA)
      contrast.matrix <- makeContrasts(groupAOtherA - groupAOtherB, levels = designA)
      fit2 <- contrasts.fit(fit, contrast.matrix)
      fit2 <- eBayes(fit2)
      options(digits = 4)
      top1 <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
      genes <- rownames(top1)[which(top1$P.Value < adjp & top1$logFC > FC)]
      dbs <- listEnrichrDbs() 
      enrichRLive <- TRUE 
      if (is.null(dbs)) { 
        enrichRLive <- FALSE 
      }
      dbs <- input$SelectedPathway
      enriched <- enrichr(genes, dbs) # Plot top 20 GO-BP results ordered by P-value 
      x <- plotEnrich(enriched[[1]], showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value")
      ggsave(file,x, width = 10, height = 8)
    }
  )
  
  output$dnldPlotSVG_dnPath <- downloadHandler(
    filename = function() {
      paste(gsub(" ","",ProjectName),"_DownRegulated_Pathways.svg", sep = '')
    },
    content = function(file) {
      adjp <- input$pathpval
      FC <- input$pathFC
      if (ncol(meta) > 2) {
        metacol <- input$EnrichRPathMetaCol
      }
      else if (ncol(meta) == 2) {
        metacol <- colnames(meta)[2]
      }
      A <- meta[which(meta[,metacol] == input$comparisonA2.path),1]
      B <- meta[which(meta[,metacol] == input$comparisonB2.path),1]
      mat <- expr[,c(A,B)]
      mat <- log2(mat + 1.0)
      groupAOther <- factor(c(rep("A", length(A)), rep("B", length(B))))
      designA <- model.matrix(~0 + groupAOther)
      fit <- lmFit(mat, design = designA)
      contrast.matrix <- makeContrasts(groupAOtherA - groupAOtherB, levels = designA)
      fit2 <- contrasts.fit(fit, contrast.matrix)
      fit2 <- eBayes(fit2)
      options(digits = 4)
      top1 <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
      genes <- rownames(top1)[which(top1$P.Value < adjp & top1$logFC > FC)]
      dbs <- listEnrichrDbs() 
      enrichRLive <- TRUE 
      if (is.null(dbs)) { 
        enrichRLive <- FALSE 
      }
      dbs <- input$SelectedPathway
      enriched <- enrichr(genes, dbs) # Plot top 20 GO-BP results ordered by P-value 
      x <- plotEnrich(enriched[[1]], showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value")
      ggsave(file,x, width = 10, height = 8)
    }
  )
  
  output$dnldPlotSVG_gsea <- downloadHandler(
    filename = function() {
      paste(gsub(" ","",ProjectName),"_GSEA_EnrichmentPlot.svg", sep = '')
    },
    content = function(file) {
      if (input$tables == 1){
        if (length(input$msigdbTable_rows_selected) > 0){
          res <- datasetInput()
          gsea.df <- as.data.frame(res@result)
          geneset <- as.character(msigdb.gsea2[input$msigdbTable_rows_selected,3])
          title <- geneset
          x <- gseaplot2(res,
                         geneset,
                         title,
                         pvalue_table = F)
        }
      }
      else if (input$tables == 3){
        if (length(input$tab2table_rows_selected) > 0){
          res <- datasetInput()
          gsea.df <- as.data.frame(res@result)
          geneset <- as.character(GeneSet2[input$tab2table_rows_selected,1])
          title <- geneset
          x <- gseaplot2(res,
                         geneset,
                         title,
                         pvalue_table = F)
        }
      }
      else if (input$tables == 5){
        if (length(input$GStable.u_rows_selected) > 0){
          res <- datasetInput()
          gsea.df <- as.data.frame(res@result)
          geneset <- as.character(user_gs_mirror()[input$GStable.u_rows_selected,1])
          title <- geneset
          x <- gseaplot2(res,
                         geneset,
                         title,
                         pvalue_table = F)
        }
      }
      ggsave(file,x, width = 10, height = 8)
    }
  )
  
  output$dnldPlotSVG_gseaHeat <- downloadHandler(
    filename = function() {
      paste(gsub(" ","",ProjectName),"_GSEA_Heatmap.svg", sep = '')
    },
    content = function(file) {
      if (input$tables == 1) {
        if (length(input$msigdbTable_rows_selected) > 0){
          groupA <- meta[which(meta[,2] == input$comparisonA),1]
          groupB <- meta[which(meta[,2] == input$comparisonB),1]
          res <- datasetInput()
          gsea.df <- as.data.frame(res@result)
          GS <- as.character(msigdb.gsea2[input$msigdbTable_rows_selected,3])
          genes1 <- as.matrix(gsea.df[which(gsea.df$Description==GS),"core_enrichment"])
          genes2 <- strsplit(genes1,"/")
          genes3 <- as.data.frame(genes2, col.names = "genes")
          gene_symbol <- genes3$genes
          ## convert expression matrix to numeric
          class(A) <- "numeric"
          ## Transforming data
          A <- A[,c(groupA,groupB)]
          exp.mat1 = log2(A + 1) # log
          exp.mat2 = apply(exp.mat1, 1, scale); # z score
          exp.mat3 = apply(exp.mat2, 1, rev); # transpose
          colnames(exp.mat3) = colnames(A) # set the column name
          exp.mat4 = exp.mat3[rowSums(is.na(exp.mat3[,])) == 0,] # remove NaN rows
          exp.mat5 = exp.mat4[rownames(exp.mat4) %in% gene_symbol,] # grab gene symbols
          # reassign data
          dataset <- exp.mat5
          ## generate color for pheatmap
          if (abs(min(dataset)) > abs(max(dataset))) {
            dataset[dataset < -abs(max(dataset))] = -abs(max(dataset))
          } else {
            dataset[dataset > abs(min(dataset))] = abs(min(dataset))
          }
          meta2 <- meta[which(meta[,1] %in% c(groupA,groupB)),]
          meta2 <- meta2[order(meta2[,2]),]
          type <- meta2[,2]
          meta3 <- as.data.frame(type)
          rownames(meta3) <- meta2[,1]
          zscore_range = 10;
          minimum = -zscore_range;
          maximum = zscore_range;
          bk = c(seq(minimum,minimum/2, length=100), seq(minimum/2,maximum/2,length=100),seq(maximum/2,maximum,length=100))
          #hmcols<- colorRampPalette(c("dark blue","blue","white","red", "dark red"))(length(bk)-1)
          #Heatmap color
          color_choice <- input$ColorPalette_gseaHeat
          col_sets <- c("OrRd","PuBu","Greens","YlGnBu","Spectral")
          if (color_choice == "original") {
            HeatMap_Colors <- c("dark blue","blue","white","red", "dark red")
            hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
          }
          else if (color_choice %in% col_sets) {
            HeatMap_Colors <- brewer.pal(n = 5, color_choice)
            hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
          }
          else if (color_choice == "Inferno") {
            hmcols <- inferno(500)
          }
          else if (color_choice == "Viridis") {
            hmcols <- viridis(500)
          }
          else if (color_choice == "Plasma") {
            hmcols <- plasma(500)
          }
          pdf(NULL)
          x <- pheatmap::pheatmap(dataset,            #data
                                  cluster_cols = F,    #cluster columns - NO
                                  cluster_row = F,     #cluster rows - YES
                                  fontsize_col = input$heatmapFont1.c,   #column fontsize
                                  fontsize_row = input$heatmapFont1.r,
                                  show_rownames = T,  
                                  show_colnames = T,
                                  color=hmcols,
                                  annotation_col = meta3)
        }
      }
      else if (input$tables == 3) {
        if (length(input$tab2table_rows_selected) > 0){
          groupA <- meta[which(meta[,2] == input$comparisonA),1]
          groupB <- meta[which(meta[,2] == input$comparisonB),1]
          res <- datasetInput()
          gsea.df <- as.data.frame(res@result)
          GS <- as.character(GeneSet2[input$tab2table_rows_selected,1])
          genes1 <- as.matrix(gsea.df[which(gsea.df$Description==GS),"core_enrichment"])
          genes2 <- strsplit(genes1,"/")
          genes3 <- as.data.frame(genes2, col.names = "genes")
          gene_symbol <- genes3$genes
          ## convert expression matrix to numeric
          class(A) <- "numeric"
          ## Transforming data
          A <- A[,c(groupA,groupB)]
          exp.mat1 = log2(A + 1) # log
          exp.mat2 = apply(exp.mat1, 1, scale); # z score
          exp.mat3 = apply(exp.mat2, 1, rev); # transpose
          colnames(exp.mat3) = colnames(A) # set the column name
          exp.mat4 = exp.mat3[rowSums(is.na(exp.mat3[,])) == 0,] # remove NaN rows
          exp.mat5 = exp.mat4[rownames(exp.mat4) %in% gene_symbol,] # grab gene symbols
          # reassign data
          dataset <- exp.mat5
          ## generate color for pheatmap
          if (abs(min(dataset)) > abs(max(dataset))) {
            dataset[dataset < -abs(max(dataset))] = -abs(max(dataset))
          } else {
            dataset[dataset > abs(min(dataset))] = abs(min(dataset))
          }
          meta2 <- meta[which(meta[,1] %in% c(groupA,groupB)),]
          meta2 <- meta2[order(meta2[,2]),]
          type <- meta2[,2]
          meta3 <- as.data.frame(type)
          rownames(meta3) <- meta2[,1]
          zscore_range = 10;
          minimum = -zscore_range;
          maximum = zscore_range;
          bk = c(seq(minimum,minimum/2, length=100), seq(minimum/2,maximum/2,length=100),seq(maximum/2,maximum,length=100))
          #hmcols<- colorRampPalette(c("dark blue","blue","white","red", "dark red"))(length(bk)-1)
          #Heatmap color
          color_choice <- input$ColorPalette_gseaHeat
          col_sets <- c("OrRd","PuBu","Greens","YlGnBu","Spectral")
          if (color_choice == "original") {
            HeatMap_Colors <- c("dark blue","blue","white","red", "dark red")
            hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
          }
          else if (color_choice %in% col_sets) {
            HeatMap_Colors <- brewer.pal(n = 5, color_choice)
            hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
          }
          else if (color_choice == "Inferno") {
            hmcols <- inferno(500)
          }
          else if (color_choice == "Viridis") {
            hmcols <- viridis(500)
          }
          else if (color_choice == "Plasma") {
            hmcols <- plasma(500)
          }
          pdf(NULL)
          x <- pheatmap::pheatmap(dataset,            #data
                                  cluster_cols = F,    #cluster columns - NO
                                  cluster_row = F,     #cluster rows - YES
                                  fontsize_col = input$heatmapFont1.c,   #column fontsize
                                  fontsize_row = input$heatmapFont1.r,
                                  show_rownames = T,  
                                  show_colnames = T,
                                  color=hmcols,
                                  annotation_col = meta3)
        }
      }
      else if (input$tables == 5) {
        if (length(input$GStable.u_rows_selected) > 0){
          groupA <- meta[which(meta[,2] == input$comparisonA),1]
          groupB <- meta[which(meta[,2] == input$comparisonB),1]
          res <- datasetInput()
          gsea.df <- as.data.frame(res@result)
          GS <- as.character(user_gs_mirror()[input$GStable.u_rows_selected,1])
          genes1 <- as.matrix(gsea.df[which(gsea.df$Description==GS),"core_enrichment"])
          genes2 <- strsplit(genes1,"/")
          genes3 <- as.data.frame(genes2, col.names = "genes")
          gene_symbol <- genes3$genes
          ## convert expression matrix to numeric
          class(A) <- "numeric"
          ## Transforming data
          A <- A[,c(groupA,groupB)]
          exp.mat1 = log2(A + 1) # log
          exp.mat2 = apply(exp.mat1, 1, scale); # z score
          exp.mat3 = apply(exp.mat2, 1, rev); # transpose
          colnames(exp.mat3) = colnames(A) # set the column name
          exp.mat4 = exp.mat3[rowSums(is.na(exp.mat3[,])) == 0,] # remove NaN rows
          exp.mat5 = exp.mat4[rownames(exp.mat4) %in% gene_symbol,] # grab gene symbols
          # reassign data
          dataset <- exp.mat5
          ## generate color for pheatmap
          if (abs(min(dataset)) > abs(max(dataset))) {
            dataset[dataset < -abs(max(dataset))] = -abs(max(dataset))
          } else {
            dataset[dataset > abs(min(dataset))] = abs(min(dataset))
          }
          meta2 <- meta[which(meta[,1] %in% c(groupA,groupB)),]
          meta2 <- meta2[order(meta2[,2]),]
          type <- meta2[,2]
          meta3 <- as.data.frame(type)
          rownames(meta3) <- meta2[,1]
          zscore_range = 10;
          minimum = -zscore_range;
          maximum = zscore_range;
          bk = c(seq(minimum,minimum/2, length=100), seq(minimum/2,maximum/2,length=100),seq(maximum/2,maximum,length=100))
          #hmcols<- colorRampPalette(c("dark blue","blue","white","red", "dark red"))(length(bk)-1)
          #Heatmap color
          color_choice <- input$ColorPalette_gseaHeat
          col_sets <- c("OrRd","PuBu","Greens","YlGnBu","Spectral")
          if (color_choice == "original") {
            HeatMap_Colors <- c("dark blue","blue","white","red", "dark red")
            hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
          }
          else if (color_choice %in% col_sets) {
            HeatMap_Colors <- brewer.pal(n = 5, color_choice)
            hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
          }
          else if (color_choice == "Inferno") {
            hmcols <- inferno(500)
          }
          else if (color_choice == "Viridis") {
            hmcols <- viridis(500)
          }
          else if (color_choice == "Plasma") {
            hmcols <- plasma(500)
          }
          pdf(NULL)
          x <- pheatmap::pheatmap(dataset,            #data
                                  cluster_cols = F,    #cluster columns - NO
                                  cluster_row = F,     #cluster rows - YES
                                  fontsize_col = input$heatmapFont1.c,   #column fontsize
                                  fontsize_row = input$heatmapFont1.r,
                                  show_rownames = T,  
                                  show_colnames = T,
                                  color=hmcols,
                                  annotation_col = meta3)
        }
      }
      ggsave(file,x, width = 10, height = 20)
    }
  )
  
  output$dnldPlotSVG_gseaBox <- downloadHandler(
    filename = function() {
      paste(gsub(" ","",ProjectName),"_ssGSEA_Boxplot.svg", sep = '')
    },
    content = function(file) {
      #if tab 2 selected
      if (input$tables == 3) {
        if (length(input$tab2table_rows_selected) > 0){
          ssgsea <- ssGSEAfunc()
          ssgsea2 <- as.data.frame(t(ssgsea))
          samporder <- meta[,1]
          ssgsea3 <- as.data.frame(ssgsea2[samporder,])
          colnames(ssgsea3)[1] <- colnames(ssgsea2[1])
          rownames(ssgsea3) <- samporder
          ssgsea4 <- ssgsea3 %>% 
            mutate(type = case_when(
              rownames(ssgsea3) == meta[,1] ~ meta[,2],
            ))
          ssgsea4 <- ssgsea4 %>%
            relocate(type)
          if (input$boxplotcompare == "none") {
            x <- ggplot(ssgsea4, aes(factor(type), ssgsea4[,2], fill = type)) +
              geom_boxplot(width = 0.5, lwd = 1, fill = "white") +
              geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = input$boxplotDotss) +
              labs(x = "Group", y = paste(colnames(ssgsea4)[2], " Score", sep = ""),
                   title = paste("ssGSEA Expression Signature: ", colnames(ssgsea4)[2],sep = "")) +
              theme_bw() +
              theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1),
                    text = element_text(size = input$boxplotFontss))
          }
          else {
            x <- ggplot(ssgsea4, aes(factor(type), ssgsea4[,2], fill = type)) +
              geom_boxplot(width = 0.5, lwd = 1, fill = "white") +
              geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = input$boxplotDotss) +
              labs(x = "Group", y = paste(colnames(ssgsea4)[2], " Score", sep = ""),
                   title = paste("ssGSEA Expression Signature: ", colnames(ssgsea4)[2],sep = "")) +
              theme_bw() +
              stat_compare_means(method = input$boxplotcompare) +
              theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1),
                    text = element_text(size = input$boxplotFontss))
          }
        }
      }
      #if MSigDB tab selected
      else if (input$tables == 1) {
        if (length(input$msigdbTable_rows_selected) > 0){
          ssgsea <- ssGSEAfunc()
          ssgsea2 <- as.data.frame(t(ssgsea))
          samporder <- meta[,1]
          ssgsea3 <- as.data.frame(ssgsea2[samporder,])
          colnames(ssgsea3)[1] <- colnames(ssgsea2[1])
          rownames(ssgsea3) <- samporder
          ssgsea4 <- ssgsea3 %>% 
            mutate(type = case_when(
              rownames(ssgsea3) == meta[,1] ~ meta[,2],
            ))
          ssgsea4 <- ssgsea4 %>%
            relocate(type)
          if (input$boxplotcompare == "none") {
            x <- ggplot(ssgsea4, aes(factor(type), ssgsea4[,2], fill = type)) +
              geom_boxplot(width = 0.5, lwd = 1, fill = "white") +
              geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = input$boxplotDotss) +
              labs(x = "Group", y = paste(colnames(ssgsea4)[2], " Score", sep = ""),
                   title = paste("ssGSEA Expression Signature: ",colnames(ssgsea4)[2],sep = "")) +
              theme_bw() +
              theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1),
                    text = element_text(size = input$boxplotFontss))
          }
          else {
            x <- ggplot(ssgsea4, aes(factor(type), ssgsea4[,2], fill = type)) +
              geom_boxplot(width = 0.5, lwd = 1, fill = "white") +
              geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = input$boxplotDotss) +
              labs(x = "Group", y = paste(colnames(ssgsea4)[2], " Score", sep = ""),
                   title = paste("ssGSEA Expression Signature: ",colnames(ssgsea4)[2],sep = "")) +
              theme_bw() +
              stat_compare_means(method = input$boxplotcompare) +
              theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1),
                    text = element_text(size = input$boxplotFontss))
          }
        }
      }
      #if user generated GS table selected
      else if (input$tables == 5) {
        if (length(input$GStable.u_rows_selected) > 0){
          ssgsea <- ssGSEAfunc()
          ssgsea2 <- as.data.frame(t(ssgsea))
          samporder <- meta[,1]
          ssgsea3 <- as.data.frame(ssgsea2[samporder,])
          colnames(ssgsea3)[1] <- colnames(ssgsea2[1])
          rownames(ssgsea3) <- samporder
          ssgsea4 <- ssgsea3 %>% 
            mutate(type = case_when(
              rownames(ssgsea3) == meta[,1] ~ meta[,2],
            ))
          ssgsea4 <- ssgsea4 %>%
            relocate(type)
          if (input$boxplotcompare == "none") {
            x <- ggplot(ssgsea4, aes(factor(type), ssgsea4[,2], fill = type)) +
              geom_boxplot(width = 0.5, lwd = 1, fill = "white") +
              geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = input$boxplotDotss) +
              labs(x = "Group", y = paste(colnames(ssgsea4)[2], " Score", sep = ""),
                   title = paste("ssGSEA Expression Signature: ",colnames(ssgsea4)[2],sep = "")) +
              theme_bw() +
              theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1),
                    text = element_text(size = input$boxplotFontss))
          }
          else {
            x <- ggplot(ssgsea4, aes(factor(type), ssgsea4[,2], fill = type)) +
              geom_boxplot(width = 0.5, lwd = 1, fill = "white") +
              geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = input$boxplotDotss) +
              labs(x = "Group", y = paste(colnames(ssgsea4)[2], " Score", sep = ""),
                   title = paste("ssGSEA Expression Signature: ",colnames(ssgsea4)[2],sep = "")) +
              theme_bw() +
              stat_compare_means(method = input$boxplotcompare) +
              theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1),
                    text = element_text(size = input$boxplotFontss))
          }
        }
      }
      ggsave(file,x, width = 10, height = 8)
    }
  )
  
  output$dnldPlotPDF_heat1 <- downloadHandler(
    filename = function() {
      paste(gsub(" ","",ProjectName),"_MVG_Heatmap.pdf", sep = '')
    },
    content = function(file) {
      
      heat <- MVGheatmap_react()
      
      ggsave(file,heat, width = 10, height = 20)
    }
  )
  
  output$dnldPlotPDF_heat2 <- downloadHandler(
    filename = function() {
      paste(gsub(" ","",ProjectName),"_Custom_Heatmap.pdf", sep = '')
    },
    content = function(file) {
      
      heat <- CustomHeatmap_react()
      ggsave(file,heat, width = 10, height = 20)
      
    }
  )
  
  output$dnldPlotPDF_heat3 <- downloadHandler(
    filename = function() {
      paste(gsub(" ","",ProjectName),"_DEGExpr_Heatmap.pdf", sep = '')
    },
    content = function(file) {
      
      heat <- DEGHeatmap_react()
      
      ggsave(file,heat, width = 10, height = 20)
    }
  )
  
  output$dnldPlotPDF_heat5 <- downloadHandler(
    filename = function() {
      paste(gsub(" ","",ProjectName),"_CustomAvgExpr_Heatmap.pdf", sep = '')
    },
    content = function(file) {
      
      heat <- AvgExprHeatmap_react()
      
      ggsave(file,heat, width = 10, height = 20)
      
    }
  )
  
  output$dnldPlotPDF_scatter <- downloadHandler(
    filename = function() {
      paste(gsub(" ","",ProjectName),"_2GeneExpression_ScatterPlot.pdf", sep = '')
    },
    content = function(file) {
      #log if user designates
      if (input$logask == TRUE) {
        expr <- log2(expr + 1)
      }
      #transpose
      expr_t <- as.data.frame(t(expr))
      #reorder rowname to match meta for merging
      samporder <- meta[,1]
      expr_t2 <- as.data.frame(expr_t[samporder,])
      #add type
      expr_t3 <- expr_t2 %>% 
        mutate(type = case_when(
          rownames(expr_t2) == meta[,1] ~ meta[,2],
        ))
      expr_t3 <- expr_t3 %>%
        relocate(type)
      #user gene input
      gene1.u <- input$scatterG1
      gene2.u <- input$scatterG2
      #get columns and info based off user input
      gene1 <- expr_t3[,gene1.u]
      gene2 <- expr_t3[,gene2.u]
      if (input$logask == TRUE) {
        gene1.n <- paste(colnames(expr_t3)[which(colnames(expr_t3) == gene1.u)], " Expression (log2 +1)", sep = "")
        gene2.n <- paste(colnames(expr_t3)[which(colnames(expr_t3) == gene2.u)], " Expression (log2 +1)", sep = "")
      }
      else if (input$logask == FALSE) {
        gene2.n <- paste(colnames(expr_t3)[which(colnames(expr_t3) == gene2.u)], " Expression", sep = "")
        gene1.n <- paste(colnames(expr_t3)[which(colnames(expr_t3) == gene1.u)], " Expression", sep = "")
      }
      #plot
      p <- ggplot(expr_t3, aes(x = gene1, y = gene2,
                               color = type)) +
        geom_point() +
        theme_minimal() +
        labs(x = gene1.n, y = gene2.n,
             title = paste(gene1.n, " vs. ", gene2.n, sep = ''))
      ggsave(file,p, width = 10, height = 8)
    }
  )
  
  output$dnldPlotPDF_exprBox <- downloadHandler(
    filename = function() {
      paste(gsub(" ","",ProjectName),"_Expression_Boxplot.pdf", sep = '')
    },
    content = function(file) {
      if (length(input$GeneListTable2_rows_selected) > 0){
        gene <- geneList[input$GeneListTable2_rows_selected, 1]
        min <- min(log2(expr[gene,] + 1.0))
        max <- max(log2(expr[gene,] + 1.0))
        meta_temp <- meta
        rownames(meta_temp) <- meta[,1]
        meta_temp <- meta_temp %>%
          select(Group)
        data = merge(t(expr[gene,]), meta_temp, by=0)
        colnames(data) = c("SampleName", "GeneExpr", "Cluster")
        x <- ggplot(data, aes(x=Cluster, y=log2(GeneExpr + 1.0))) +
          geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=1) +
          geom_dotplot(binaxis='y', stackdir='center', dotsize=input$boxplotDot2) +
          stat_compare_means(label = "p.signif") + ylim(min * 0.9, max * 1.3) +
          theme_bw() +
          labs(title= paste(gene, "Expression (log2)")) +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1),
                text = element_text(size = input$boxplot.font2))
        ggsave(file,x, width = 10, height = 8)
      }
    }
  )
  
  output$dnldPlotPDF_vol <- downloadHandler(
    filename = function() {
      paste(gsub(" ","",ProjectName),"_VolcanoPlot.pdf", sep = '')
    },
    content = function(file) {
      
      x <- Volcano3_react()
      PlotH <- input$VolMAdnldHeight
      PlotW <- input$VolMAdnldWidth
      PlotU <- input$VolMAdnldSizeUnits
      ggsave(file,x, width = PlotW, height = PlotH, units = PlotU)
    }
  )
  
  output$dnldPlotPDF_MA <- downloadHandler(
    filename = function() {
      paste(gsub(" ","",ProjectName),"_MA_Plot.pdf", sep = '')
    },
    content = function(file) {
      top2 <- topgenereact()
      #add color categories based on FC and pval
      top2['threshold'] <- "none"
      top2[which(top2$logFC > abs(input$fc_cutoff)), "threshold"] <- "up"
      top2[which(top2$logFC < -abs(input$fc_cutoff)), "threshold"] <- "down"
      upRed <- "lightcoral"
      dnBlue <- "cadetblue3"
      mdGray <- "gray70"
      top2u <- top2[order(top2[,1], decreasing = TRUE),]
      top2d <- top2[order(top2[,1]),]
      top_hits_up <- top2u[head(which(top2u$logFC > abs(input$fc_cutoff)), n = input$top_x),]
      top_hits_dn <- top2d[head(which(top2d$logFC < -abs(input$fc_cutoff)), n = input$top_x),]
      #select genes to label based on user selection
      genesel.s <- NULL
      genesel.t <- NULL
      genesel.u <- NULL
      genesel.s <- unlist(strsplit(input$userGeneSelec2, " "))
      genesel.t <- unlist(strsplit(input$userGeneSelec2, "\t"))
      genesel.u <- input$userGeneSelec
      genesel.text <- c(genesel.s,genesel.t,genesel.u)
      top2_selec <- top2 %>%
        filter(GeneName %in% genesel.text)
      x <- ggplot(data = top2, aes(x = AveExpr, y = logFC)) +
        geom_point(size = 2, shape = 16) +
        theme_light(base_size = 16)
      #colors
      x <- x + aes(color = threshold) +
        scale_color_manual(values = c("up" = upRed,"down" = dnBlue, "none" = mdGray))
      x <- x + geom_hline(yintercept = c(-abs(input$fc_cutoff),abs(input$fc_cutoff)), linetype="dashed", color="gray20")
      x <- x + geom_text_repel(
        data =  top_hits_up,
        aes(label = rownames(top_hits_up)),
        color="gray20",
        size = 6,
        nudge_x = 0.2,
        nudge_y=0.2,
        box.padding = unit(0.9, "lines"),
        point.padding = unit(.3+4*0.1, "lines"),
        max.overlaps = 50)
      x <- x + geom_text_repel(
        data =  top_hits_dn,
        aes(label = rownames(top_hits_dn)),
        color="gray20",
        size = 6,
        nudge_x = 0.2,
        nudge_y=0.2,
        box.padding = unit(0.9, "lines"),
        point.padding = unit(.3+4*0.1, "lines"),
        max.overlaps = 50)
      x <- x + geom_text_repel(
        data =  top2_selec,
        aes(label = rownames(top2_selec)),
        color="gray20",
        size = 6,
        nudge_x = 0.2,
        nudge_y=0.2,
        box.padding = unit(0.9, "lines"),
        point.padding = unit(.3+4*0.1, "lines"),
        max.overlaps = 50)
      #coloring selected points
      x <- x + geom_point(data = top2_selec,
                          aes(x = AveExpr, y = logFC),
                          pch = 21,
                          color = "black",
                          size = 2)
      x <- x + geom_point(data = top_hits_dn,
                          aes(x = AveExpr, y = logFC),
                          pch = 21,
                          color = "black",
                          size = 2)
      x <- x + geom_point(data = top_hits_up,
                          aes(x = AveExpr, y = logFC),
                          pch = 21,
                          color = "black",
                          size = 2)
      x <- x + theme(legend.position="none")
      x <- x + labs(x = "Average Expression", y = "log2FC")
      x <- x + theme(axis.text = element_text(size=18))
      x <- x + theme(axis.title = element_text(size=24))
      ggsave(file,x, width = 10, height = 8)
    }
  )
  
  output$dnldPlotPDF_exprBox2 <- downloadHandler(
    filename = function() {
      paste(gsub(" ","",ProjectName),"_ExpressionBoxplot.pdf", sep = '')
    },
    content = function(file) {
      if (length(input$GeneListTable_rows_selected) > 0){
        gene <- geneList[input$GeneListTable_rows_selected, 1]
        min <- min(log2(expr[gene,] + 1.0))
        max <- max(log2(expr[gene,] + 1.0))
        meta_temp <- meta
        rownames(meta_temp) <- meta[,1]
        meta_temp <- meta_temp %>%
          select(Group)
        data = merge(t(expr[gene,]), meta_temp, by=0)
        colnames(data) = c("SampleName", "GeneExpr", "Cluster")
        x <- ggplot(data, aes(x=Cluster, y=log2(GeneExpr + 1.0))) +
          geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=1) +
          geom_dotplot(binaxis='y', stackdir='center', dotsize=input$boxplotDot) +
          stat_compare_means(label = "p.signif") + ylim(min * 0.9, max * 1.3) +
          theme_bw() +
          labs(title= paste(gene, "Expression (log2)")) +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1),
                text = element_text(size = input$boxplot.font))
        ggsave(file,x, width = 10, height = 8)
      }
    }
  )
  
  output$dnldPlotPDF_AvgScatter <- downloadHandler(
    filename = function() {
      paste(gsub(" ","",ProjectName),"_AvgExpression_ScatterPlot.pdf", sep = '')
    },
    content = function(file) {
      A_choice <- input$comparisonA2.avg
      B_choice <- input$comparisonB2.avg
      log_choice <- input$AvgExpLogFC
      
      #select genes to label based on user selection
      genesel.s <- NULL
      genesel.t <- NULL
      genesel.u <- NULL
      genesel.s <- unlist(strsplit(input$gsSelection2, " "))
      genesel.t <- unlist(strsplit(input$gsSelection2, "\t"))
      genesel.u <- input$scatterGeneSelec
      genesel.text <- c(genesel.s,genesel.t,genesel.u)
      
      if (ncol(meta) > 2){
        metacol <- input$AvgExprScatterMetaCol
      }
      else if (ncol(meta) == 2){
        metacol <- colnames(meta)[2]
      }
      
      # Get A and B sample names
      A <- meta[which(meta[,metacol] == A_choice),1]
      B <- meta[which(meta[,metacol] == B_choice),1]
      
      # Make A and B expression Matrices
      mat_A <- expr[,A]
      mat_B <- expr[,B]
      
      logFC_text <- ""
      
      # Log if user designates
      if (log_choice == TRUE) {
        mat_A <- log2(mat_A + 1)
        mat_B <- log2(mat_B + 1)
        logFC_text <- " (log2 +1)"
      }
      
      # Get avg expression of each gene
      mat_A$AvgExpression_GroupA <- rowMeans(mat_A)
      mat_B$AvgExpression_GroupB <- rowMeans(mat_B)
      
      # Add Gene Names as column
      mat_A$GeneSymbol <- rownames(mat_A)
      mat_B$GeneSymbol <- rownames(mat_B)
      
      # Merge average columns
      AvgExpr_Table <- merge(mat_A[,c("AvgExpression_GroupA","GeneSymbol")],
                             mat_B[,c("AvgExpression_GroupB","GeneSymbol")],
                             by = "GeneSymbol",
                             all = T)
      
      # Add group for above or below abline
      AvgExpr_Table <- AvgExpr_Table %>%
        mutate(ColorCol = case_when(
          AvgExpr_Table$AvgExpression_GroupB > AvgExpr_Table$AvgExpression_GroupA ~ paste(B_choice," Expr > ",A_choice,' Expr', sep = ''),
          AvgExpr_Table$AvgExpression_GroupA > AvgExpr_Table$AvgExpression_GroupB ~ paste(A_choice," Expr > ",B_choice,' Expr', sep = '')
        ))
      
      AvgExpr_Table_selec <- AvgExpr_Table %>%
        filter(GeneSymbol %in% genesel.text)
      
      # plot
      p <- ggplot(AvgExpr_Table, aes(x = AvgExpression_GroupA, y = AvgExpression_GroupB,
                                     color = ColorCol)) +
        geom_point() +
        theme_minimal() +
        theme(legend.position="none") +
        xlab(paste("Average Expression: ",A_choice,logFC_text,sep = ""))+
        ylab(paste("Average Expression: ",B_choice,logFC_text,sep = "")) +
        labs(title = paste("Average Expression",logFC_text,": ",A_choice," vs. ",B_choice,sep = "")) +
        scale_color_manual(values = c("cadetblue3","lightcoral"))
      
      p <- p + geom_text_repel(
        data =  AvgExpr_Table_selec,
        aes(label = GeneSymbol),
        size = 6,
        color="gray20",
        nudge_x = 0.2,
        nudge_y=0.2,
        box.padding = unit(0.9, "lines"),
        point.padding = unit(.3+4*0.1, "lines"),
        max.overlaps = 50)
      
      #coloring selected points
      p <- p + geom_point(data = AvgExpr_Table_selec,
                          aes(x = AvgExpression_GroupA, y = AvgExpression_GroupB),
                          pch = 21,
                          color = "black",
                          size = 2)
      
      ggsave(file,p, width = 10, height = 8)
    }
  )
  
  output$dnldPlotPDF_upPath <- downloadHandler(
    filename = function() {
      paste(gsub(" ","",ProjectName),"_UpRegulated_Pathways.pdf", sep = '')
    },
    content = function(file) {
      adjp <- input$pathpval
      FC <- input$pathFC
      if (ncol(meta) > 2) {
        metacol <- input$EnrichRPathMetaCol
      }
      else if (ncol(meta) == 2) {
        metacol <- colnames(meta)[2]
      }
      A <- meta[which(meta[,metacol] == input$comparisonA2.path),1]
      B <- meta[which(meta[,metacol] == input$comparisonB2.path),1]
      mat <- expr[,c(A,B)]
      mat <- log2(mat + 1.0)
      groupAOther <- factor(c(rep("A", length(A)), rep("B", length(B))))
      designA <- model.matrix(~0 + groupAOther)
      fit <- lmFit(mat, design = designA)
      contrast.matrix <- makeContrasts(groupAOtherA - groupAOtherB, levels = designA)
      fit2 <- contrasts.fit(fit, contrast.matrix)
      fit2 <- eBayes(fit2)
      options(digits = 4)
      top1 <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
      genes <- rownames(top1)[which(top1$P.Value < adjp & top1$logFC > FC)]
      dbs <- listEnrichrDbs() 
      enrichRLive <- TRUE 
      if (is.null(dbs)) { 
        enrichRLive <- FALSE 
      }
      dbs <- input$SelectedPathway
      enriched <- enrichr(genes, dbs) # Plot top 20 GO-BP results ordered by P-value 
      x <- plotEnrich(enriched[[1]], showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value")
      ggsave(file,x, width = 15, height = 8)
    }
  )
  
  output$dnldPlotPDF_dnPath <- downloadHandler(
    filename = function() {
      paste(gsub(" ","",ProjectName),"_DownRegulated_Pathways.pdf", sep = '')
    },
    content = function(file) {
      adjp <- input$pathpval
      FC <- input$pathFC
      if (ncol(meta) > 2) {
        metacol <- input$EnrichRPathMetaCol
      }
      else if (ncol(meta) == 2) {
        metacol <- colnames(meta)[2]
      }
      A <- meta[which(meta[,metacol] == input$comparisonA2.path),1]
      B <- meta[which(meta[,metacol] == input$comparisonB2.path),1]
      mat <- expr[,c(A,B)]
      mat <- log2(mat + 1.0)
      groupAOther <- factor(c(rep("A", length(A)), rep("B", length(B))))
      designA <- model.matrix(~0 + groupAOther)
      fit <- lmFit(mat, design = designA)
      contrast.matrix <- makeContrasts(groupAOtherA - groupAOtherB, levels = designA)
      fit2 <- contrasts.fit(fit, contrast.matrix)
      fit2 <- eBayes(fit2)
      options(digits = 4)
      top1 <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
      genes <- rownames(top1)[which(top1$P.Value < adjp & top1$logFC > FC)]
      dbs <- listEnrichrDbs() 
      enrichRLive <- TRUE 
      if (is.null(dbs)) { 
        enrichRLive <- FALSE 
      }
      dbs <- input$SelectedPathway
      enriched <- enrichr(genes, dbs) # Plot top 20 GO-BP results ordered by P-value 
      x <- plotEnrich(enriched[[1]], showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value")
      ggsave(file,x, width = 15, height = 8)
    }
  )
  
  output$dnldPlotPDF_gsea <- downloadHandler(
    filename = function() {
      paste(gsub(" ","",ProjectName),"_GSEA_EnrichmentPlot.pdf", sep = '')
    },
    content = function(file) {
      if (input$tables == 1){
        if (length(input$msigdbTable_rows_selected) > 0){
          res <- datasetInput()
          gsea.df <- as.data.frame(res@result)
          geneset <- as.character(msigdb.gsea2[input$msigdbTable_rows_selected,3])
          title <- geneset
          x <- gseaplot2(res,
                         geneset,
                         title,
                         pvalue_table = F)
        }
      }
      else if (input$tables == 3){
        if (length(input$tab2table_rows_selected) > 0){
          res <- datasetInput()
          gsea.df <- as.data.frame(res@result)
          geneset <- as.character(GeneSet2[input$tab2table_rows_selected,1])
          title <- geneset
          x <- gseaplot2(res,
                         geneset,
                         title,
                         pvalue_table = F)
        }
      }
      else if (input$tables == 5){
        if (length(input$GStable.u_rows_selected) > 0){
          res <- datasetInput()
          gsea.df <- as.data.frame(res@result)
          geneset <- as.character(user_gs_mirror()[input$GStable.u_rows_selected,1])
          title <- geneset
          x <- gseaplot2(res,
                         geneset,
                         title,
                         pvalue_table = F)
        }
      }
      ggsave(file,x, width = 10, height = 8)
    }
  )
  
  output$dnldPlotPDF_gseaHeat <- downloadHandler(
    filename = function() {
      paste(gsub(" ","",ProjectName),"_GSEA_Heatmap.pdf", sep = '')
    },
    content = function(file) {
      if (input$tables == 1) {
        if (length(input$msigdbTable_rows_selected) > 0){
          groupA <- meta[which(meta[,2] == input$comparisonA),1]
          groupB <- meta[which(meta[,2] == input$comparisonB),1]
          res <- datasetInput()
          gsea.df <- as.data.frame(res@result)
          GS <- as.character(msigdb.gsea2[input$msigdbTable_rows_selected,3])
          genes1 <- as.matrix(gsea.df[which(gsea.df$Description==GS),"core_enrichment"])
          genes2 <- strsplit(genes1,"/")
          genes3 <- as.data.frame(genes2, col.names = "genes")
          gene_symbol <- genes3$genes
          ## convert expression matrix to numeric
          class(A) <- "numeric"
          ## Transforming data
          A <- A[,c(groupA,groupB)]
          exp.mat1 = log2(A + 1) # log
          exp.mat2 = apply(exp.mat1, 1, scale); # z score
          exp.mat3 = apply(exp.mat2, 1, rev); # transpose
          colnames(exp.mat3) = colnames(A) # set the column name
          exp.mat4 = exp.mat3[rowSums(is.na(exp.mat3[,])) == 0,] # remove NaN rows
          exp.mat5 = exp.mat4[rownames(exp.mat4) %in% gene_symbol,] # grab gene symbols
          # reassign data
          dataset <- exp.mat5
          ## generate color for pheatmap
          if (abs(min(dataset)) > abs(max(dataset))) {
            dataset[dataset < -abs(max(dataset))] = -abs(max(dataset))
          } else {
            dataset[dataset > abs(min(dataset))] = abs(min(dataset))
          }
          meta2 <- meta[which(meta[,1] %in% c(groupA,groupB)),]
          meta2 <- meta2[order(meta2[,2]),]
          type <- meta2[,2]
          meta3 <- as.data.frame(type)
          rownames(meta3) <- meta2[,1]
          zscore_range = 10;
          minimum = -zscore_range;
          maximum = zscore_range;
          bk = c(seq(minimum,minimum/2, length=100), seq(minimum/2,maximum/2,length=100),seq(maximum/2,maximum,length=100))
          #hmcols<- colorRampPalette(c("dark blue","blue","white","red", "dark red"))(length(bk)-1)
          #Heatmap color
          color_choice <- input$ColorPalette_gseaHeat
          col_sets <- c("OrRd","PuBu","Greens","YlGnBu","Spectral")
          if (color_choice == "original") {
            HeatMap_Colors <- c("dark blue","blue","white","red", "dark red")
            hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
          }
          else if (color_choice %in% col_sets) {
            HeatMap_Colors <- brewer.pal(n = 5, color_choice)
            hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
          }
          else if (color_choice == "Inferno") {
            hmcols <- inferno(500)
          }
          else if (color_choice == "Viridis") {
            hmcols <- viridis(500)
          }
          else if (color_choice == "Plasma") {
            hmcols <- plasma(500)
          }
          pdf(NULL)
          x <- pheatmap(dataset,            #data
                        cluster_cols = F,    #cluster columns - NO
                        cluster_row = F,     #cluster rows - YES
                        fontsize_col = input$heatmapFont1.c,   #column fontsize
                        fontsize_row = input$heatmapFont1.r,
                        show_rownames = T,  
                        show_colnames = T,
                        color=hmcols,
                        annotation_col = meta3)
        }
      }
      else if (input$tables == 3) {
        if (length(input$tab2table_rows_selected) > 0){
          groupA <- meta[which(meta[,2] == input$comparisonA),1]
          groupB <- meta[which(meta[,2] == input$comparisonB),1]
          res <- datasetInput()
          gsea.df <- as.data.frame(res@result)
          GS <- as.character(GeneSet2[input$tab2table_rows_selected,1])
          genes1 <- as.matrix(gsea.df[which(gsea.df$Description==GS),"core_enrichment"])
          genes2 <- strsplit(genes1,"/")
          genes3 <- as.data.frame(genes2, col.names = "genes")
          gene_symbol <- genes3$genes
          ## convert expression matrix to numeric
          class(A) <- "numeric"
          ## Transforming data
          A <- A[,c(groupA,groupB)]
          exp.mat1 = log2(A + 1) # log
          exp.mat2 = apply(exp.mat1, 1, scale); # z score
          exp.mat3 = apply(exp.mat2, 1, rev); # transpose
          colnames(exp.mat3) = colnames(A) # set the column name
          exp.mat4 = exp.mat3[rowSums(is.na(exp.mat3[,])) == 0,] # remove NaN rows
          exp.mat5 = exp.mat4[rownames(exp.mat4) %in% gene_symbol,] # grab gene symbols
          # reassign data
          dataset <- exp.mat5
          ## generate color for pheatmap
          if (abs(min(dataset)) > abs(max(dataset))) {
            dataset[dataset < -abs(max(dataset))] = -abs(max(dataset))
          } else {
            dataset[dataset > abs(min(dataset))] = abs(min(dataset))
          }
          meta2 <- meta[which(meta[,1] %in% c(groupA,groupB)),]
          meta2 <- meta2[order(meta2[,2]),]
          type <- meta2[,2]
          meta3 <- as.data.frame(type)
          rownames(meta3) <- meta2[,1]
          zscore_range = 10;
          minimum = -zscore_range;
          maximum = zscore_range;
          bk = c(seq(minimum,minimum/2, length=100), seq(minimum/2,maximum/2,length=100),seq(maximum/2,maximum,length=100))
          #hmcols<- colorRampPalette(c("dark blue","blue","white","red", "dark red"))(length(bk)-1)
          #Heatmap color
          color_choice <- input$ColorPalette_gseaHeat
          col_sets <- c("OrRd","PuBu","Greens","YlGnBu","Spectral")
          if (color_choice == "original") {
            HeatMap_Colors <- c("dark blue","blue","white","red", "dark red")
            hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
          }
          else if (color_choice %in% col_sets) {
            HeatMap_Colors <- brewer.pal(n = 5, color_choice)
            hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
          }
          else if (color_choice == "Inferno") {
            hmcols <- inferno(500)
          }
          else if (color_choice == "Viridis") {
            hmcols <- viridis(500)
          }
          else if (color_choice == "Plasma") {
            hmcols <- plasma(500)
          }
          pdf(NULL)
          x <- pheatmap(dataset,            #data
                        cluster_cols = F,    #cluster columns - NO
                        cluster_row = F,     #cluster rows - YES
                        fontsize_col = input$heatmapFont1.c,   #column fontsize
                        fontsize_row = input$heatmapFont1.r,
                        show_rownames = T,  
                        show_colnames = T,
                        color=hmcols,
                        annotation_col = meta3)
        }
      }
      else if (input$tables == 5) {
        if (length(input$GStable.u_rows_selected) > 0){
          groupA <- meta[which(meta[,2] == input$comparisonA),1]
          groupB <- meta[which(meta[,2] == input$comparisonB),1]
          res <- datasetInput()
          gsea.df <- as.data.frame(res@result)
          GS <- as.character(user_gs_mirror()[input$GStable.u_rows_selected,1])
          genes1 <- as.matrix(gsea.df[which(gsea.df$Description==GS),"core_enrichment"])
          genes2 <- strsplit(genes1,"/")
          genes3 <- as.data.frame(genes2, col.names = "genes")
          gene_symbol <- genes3$genes
          ## convert expression matrix to numeric
          class(A) <- "numeric"
          ## Transforming data
          A <- A[,c(groupA,groupB)]
          exp.mat1 = log2(A + 1) # log
          exp.mat2 = apply(exp.mat1, 1, scale); # z score
          exp.mat3 = apply(exp.mat2, 1, rev); # transpose
          colnames(exp.mat3) = colnames(A) # set the column name
          exp.mat4 = exp.mat3[rowSums(is.na(exp.mat3[,])) == 0,] # remove NaN rows
          exp.mat5 = exp.mat4[rownames(exp.mat4) %in% gene_symbol,] # grab gene symbols
          # reassign data
          dataset <- exp.mat5
          ## generate color for pheatmap
          if (abs(min(dataset)) > abs(max(dataset))) {
            dataset[dataset < -abs(max(dataset))] = -abs(max(dataset))
          } else {
            dataset[dataset > abs(min(dataset))] = abs(min(dataset))
          }
          meta2 <- meta[which(meta[,1] %in% c(groupA,groupB)),]
          meta2 <- meta2[order(meta2[,2]),]
          type <- meta2[,2]
          meta3 <- as.data.frame(type)
          rownames(meta3) <- meta2[,1]
          zscore_range = 10;
          minimum = -zscore_range;
          maximum = zscore_range;
          bk = c(seq(minimum,minimum/2, length=100), seq(minimum/2,maximum/2,length=100),seq(maximum/2,maximum,length=100))
          #hmcols<- colorRampPalette(c("dark blue","blue","white","red", "dark red"))(length(bk)-1)
          #Heatmap color
          color_choice <- input$ColorPalette_gseaHeat
          col_sets <- c("OrRd","PuBu","Greens","YlGnBu","Spectral")
          if (color_choice == "original") {
            HeatMap_Colors <- c("dark blue","blue","white","red", "dark red")
            hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
          }
          else if (color_choice %in% col_sets) {
            HeatMap_Colors <- brewer.pal(n = 5, color_choice)
            hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
          }
          else if (color_choice == "Inferno") {
            hmcols <- inferno(500)
          }
          else if (color_choice == "Viridis") {
            hmcols <- viridis(500)
          }
          else if (color_choice == "Plasma") {
            hmcols <- plasma(500)
          }
          pdf(NULL)
          x <- pheatmap(dataset,            #data
                        cluster_cols = F,    #cluster columns - NO
                        cluster_row = F,     #cluster rows - YES
                        fontsize_col = input$heatmapFont1.c,   #column fontsize
                        fontsize_row = input$heatmapFont1.r,
                        show_rownames = T,  
                        show_colnames = T,
                        color=hmcols,
                        annotation_col = meta3)
        }
      }
      ggsave(file,x, width = 10, height = 20)
    }
  )
  
  output$dnldPlotPDF_gseaBox <- downloadHandler(
    filename = function() {
      paste(gsub(" ","",ProjectName),"_ssGSEA_Boxplot.pdf", sep = '')
    },
    content = function(file) {
      #if tab 2 selected
      if (input$tables == 3) {
        if (length(input$tab2table_rows_selected) > 0){
          ssgsea <- ssGSEAfunc()
          ssgsea2 <- as.data.frame(t(ssgsea))
          samporder <- meta[,1]
          ssgsea3 <- as.data.frame(ssgsea2[samporder,])
          colnames(ssgsea3)[1] <- colnames(ssgsea2[1])
          rownames(ssgsea3) <- samporder
          ssgsea4 <- ssgsea3 %>% 
            mutate(type = case_when(
              rownames(ssgsea3) == meta[,1] ~ meta[,2],
            ))
          ssgsea4 <- ssgsea4 %>%
            relocate(type)
          if (input$boxplotcompare == "none") {
            x <- ggplot(ssgsea4, aes(factor(type), ssgsea4[,2], fill = type)) +
              geom_boxplot(width = 0.5, lwd = 1, fill = "white") +
              geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = input$boxplotDotss) +
              labs(x = "Group", y = paste(colnames(ssgsea4)[2], " Score", sep = ""),
                   title = paste("ssGSEA Expression Signature: ", colnames(ssgsea4)[2],sep = "")) +
              theme_bw() +
              theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1),
                    text = element_text(size = input$boxplotFontss))
          }
          else {
            x <- ggplot(ssgsea4, aes(factor(type), ssgsea4[,2], fill = type)) +
              geom_boxplot(width = 0.5, lwd = 1, fill = "white") +
              geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = input$boxplotDotss) +
              labs(x = "Group", y = paste(colnames(ssgsea4)[2], " Score", sep = ""),
                   title = paste("ssGSEA Expression Signature: ", colnames(ssgsea4)[2],sep = "")) +
              theme_bw() +
              stat_compare_means(method = input$boxplotcompare) +
              theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1),
                    text = element_text(size = input$boxplotFontss))
          }
        }
      }
      #if MSigDB tab selected
      else if (input$tables == 1) {
        if (length(input$msigdbTable_rows_selected) > 0){
          ssgsea <- ssGSEAfunc()
          ssgsea2 <- as.data.frame(t(ssgsea))
          samporder <- meta[,1]
          ssgsea3 <- as.data.frame(ssgsea2[samporder,])
          colnames(ssgsea3)[1] <- colnames(ssgsea2[1])
          rownames(ssgsea3) <- samporder
          ssgsea4 <- ssgsea3 %>% 
            mutate(type = case_when(
              rownames(ssgsea3) == meta[,1] ~ meta[,2],
            ))
          ssgsea4 <- ssgsea4 %>%
            relocate(type)
          if (input$boxplotcompare == "none") {
            x <- ggplot(ssgsea4, aes(factor(type), ssgsea4[,2], fill = type)) +
              geom_boxplot(width = 0.5, lwd = 1, fill = "white") +
              geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = input$boxplotDotss) +
              labs(x = "Group", y = paste(colnames(ssgsea4)[2], " Score", sep = ""),
                   title = paste("ssGSEA Expression Signature: ",colnames(ssgsea4)[2],sep = "")) +
              theme_bw() +
              theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1),
                    text = element_text(size = input$boxplotFontss))
          }
          else {
            x <- ggplot(ssgsea4, aes(factor(type), ssgsea4[,2], fill = type)) +
              geom_boxplot(width = 0.5, lwd = 1, fill = "white") +
              geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = input$boxplotDotss) +
              labs(x = "Group", y = paste(colnames(ssgsea4)[2], " Score", sep = ""),
                   title = paste("ssGSEA Expression Signature: ",colnames(ssgsea4)[2],sep = "")) +
              theme_bw() +
              stat_compare_means(method = input$boxplotcompare) +
              theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1),
                    text = element_text(size = input$boxplotFontss))
          }
        }
      }
      #if user generated GS table selected
      else if (input$tables == 5) {
        if (length(input$GStable.u_rows_selected) > 0){
          ssgsea <- ssGSEAfunc()
          ssgsea2 <- as.data.frame(t(ssgsea))
          samporder <- meta[,1]
          ssgsea3 <- as.data.frame(ssgsea2[samporder,])
          colnames(ssgsea3)[1] <- colnames(ssgsea2[1])
          rownames(ssgsea3) <- samporder
          ssgsea4 <- ssgsea3 %>% 
            mutate(type = case_when(
              rownames(ssgsea3) == meta[,1] ~ meta[,2],
            ))
          ssgsea4 <- ssgsea4 %>%
            relocate(type)
          if (input$boxplotcompare == "none") {
            x <- ggplot(ssgsea4, aes(factor(type), ssgsea4[,2], fill = type)) +
              geom_boxplot(width = 0.5, lwd = 1, fill = "white") +
              geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = input$boxplotDotss) +
              labs(x = "Group", y = paste(colnames(ssgsea4)[2], " Score", sep = ""),
                   title = paste("ssGSEA Expression Signature: ",colnames(ssgsea4)[2],sep = "")) +
              theme_bw() +
              theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1),
                    text = element_text(size = input$boxplotFontss))
          }
          else {
            x <- ggplot(ssgsea4, aes(factor(type), ssgsea4[,2], fill = type)) +
              geom_boxplot(width = 0.5, lwd = 1, fill = "white") +
              geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = input$boxplotDotss) +
              labs(x = "Group", y = paste(colnames(ssgsea4)[2], " Score", sep = ""),
                   title = paste("ssGSEA Expression Signature: ",colnames(ssgsea4)[2],sep = "")) +
              theme_bw() +
              stat_compare_means(method = input$boxplotcompare) +
              theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1),
                    text = element_text(size = input$boxplotFontss))
          }
        }
      }
      ggsave(file,x, width = 10, height = 8)
    }
  )
  
  output$downloadClusters <- downloadHandler(
    filename = function() {
      paste(ProjectName,"_ClusterResults.tsv", sep = '')
    },
    content = function(file) {
      top_probes <- input$NumFeatures
      col_labels <- colnames(expr)
      isexpr <- rowSums(expr[,col_labels] > 1) >= length(col_labels)
      exp <- expr[isexpr,]
      mad <- NULL
      var <- NULL
      cv <- NULL
      var_type <- input$VarianceMeasure
      if (var_type == "MAD"){
        mad <- apply(log2(exp + 1), 1, mad)
        mad <- sort(mad, decreasing = T)
        mad <- head(mad, n = (top_probes +1))
        out <- cbind(names(mad), mad[names(mad)], exp[names(mad),])
        colnames(out) <- c("Gene", "MAD", colnames(exp))
        dataset <- exp[names(mad),]
        variable_gene_list <- names(mad)
      }
      if (var_type == "VAR"){
        var <- apply(log2(exp + 1), 1, var)
        var <- sort(var, decreasing = T)
        var <- head(var, n = (top_probes +1))
        out <- cbind(names(var), var[names(var)], exp[names(var),])
        colnames(out) <- c("Gene", "VAR", colnames(exp))
        dataset <- exp[names(var),]
        variable_gene_list <- names(var)
      }
      if (var_type == "CV"){
        cv <- apply(log2(exp + 1), 1, cv)
        cv <- sort(cv, decreasing = T)
        cv <- head(cv, n = (top_probes +1))
        out <- cbind(names(cv), cv[names(cv)], exp[names(cv),])
        colnames(out) <- c("Gene", "VAR", colnames(exp))
        dataset <- exp[names(cv),]
        variable_gene_list <- names(cv)
      }
      dataset <- log2(dataset + 1)
      zdataset <- apply(dataset, 1, scale)
      zdataset <- apply(zdataset, 1, rev)
      colnames(zdataset) <- names(dataset)
      dataset <- as.matrix(zdataset)
      dataset[is.na(dataset)] <- 0
      dataset = dataset[apply(dataset[,-1], 1, function(x) !all(x==0)),]
      minimum = -5;
      maximum = 5;
      if (abs(min(dataset)) > abs(max(dataset))) {
        dataset[dataset < -abs(max(dataset))] = -abs(max(dataset))
      } else {
        dataset[dataset > abs(min(dataset))] = abs(min(dataset))
      }
      results2 = hclust(dist(t(dataset)), method = input$ClusteringMethod)
      m = sort(cutree(results2, k=input$NumClusters))
      output = cbind(colnames(m), as.matrix(m))
      colnames(output) = c("Cluster")
      output <- as.data.frame(output)
      output$SampleName <- rownames(output)
      output <- output %>%
        relocate(SampleName)
      write_tsv(as.data.frame(output),file)
    }
  )
  
  #download gene scatter plot expression data
  output$geneScatterDownload <- downloadHandler(
    filename = function() {
      #user gene input
      gene1.u <- input$scatterG1
      gene2.u <- input$scatterG2
      paste(gene1.u, "_vs_", gene2.u, "_Expression.tsv", sep = "")
    },
    content = function(file) {
      #transpose
      expr_t <- as.data.frame(t(expr))
      #reorder rowname to match meta for merging
      samporder <- meta[,1]
      expr_t2 <- as.data.frame(expr_t[samporder,])
      #add type
      expr_t3 <- expr_t2 %>% 
        mutate(type = case_when(
          rownames(expr_t2) == meta[,1] ~ meta[,2],
        ))
      expr_t3 <- expr_t3 %>%
        relocate(type)
      #user gene input
      gene1.u <- input$scatterG1
      gene2.u <- input$scatterG2
      #get columns and info based off user input
      gene1 <- expr_t3[,gene1.u]
      gene1.n <- paste(colnames(expr_t3)[which(colnames(expr_t3) == gene1.u)], " Expression", sep = "")
      gene2 <- expr_t3[,gene2.u]
      gene2.n <- paste(colnames(expr_t3)[which(colnames(expr_t3) == gene2.u)], " Expression", sep = "")
      #make table
      Sample <- rownames(expr_t3)
      Type <- expr_t3[,'type']
      gene1col <- expr_t3[,gene1.u]
      gene2col <- expr_t3[,gene2.u]
      scatterTab <- data.frame(Sample, Type, gene1col, gene2col)
      colnames(scatterTab)[c(3,4)] <- c(gene1.n, gene2.n)
      write_tsv(scatterTab, file)
    }
  )
  
  output$AvggeneScatterDownload <- downloadHandler(
    filename = function() {
      A_choice <- gsub(" ","",input$comparisonA2.avg)
      B_choice <- gsub(" ","",input$comparisonB2.avg)
      log_choice <- input$AvgExpLogFC
      paste("AvgExpr_",A_choice,"_vs_",B_choice,".tsv",sep = "")
    },
    content = function(file) {
      A_choice <- input$comparisonA2.avg
      B_choice <- input$comparisonB2.avg
      log_choice <- input$AvgExpLogFC
      
      if (ncol(meta) > 2){
        metacol <- input$AvgExprScatterMetaCol
      }
      else if (ncol(meta) == 2){
        metacol <- colnames(meta)[2]
      }
      
      # Get A and B sample names
      A <- meta[which(meta[,metacol] == A_choice),1]
      B <- meta[which(meta[,metacol] == B_choice),1]
      
      # Make A and B expression Matrices
      mat_A <- expr[,A]
      mat_B <- expr[,B]
      
      # Get avg expression of each gene
      mat_A$AvgExpression_GroupA <- rowMeans(mat_A)
      mat_B$AvgExpression_GroupB <- rowMeans(mat_B)
      
      # Add Gene Names as column
      mat_A$GeneSymbol <- rownames(mat_A)
      mat_B$GeneSymbol <- rownames(mat_B)
      
      # Merge average columns
      AvgExpr_Table <- merge(mat_A[,c("AvgExpression_GroupA","GeneSymbol")],
                             mat_B[,c("AvgExpression_GroupB","GeneSymbol")],
                             by = "GeneSymbol",
                             all = T)
      
      A_choice_lab <- gsub(" ","",A_choice)
      B_choice_lab <- gsub(" ","",B_choice)
      
      colnames(AvgExpr_Table)[2] <- paste("Average_Expression_", A_choice_lab, sep = "")
      colnames(AvgExpr_Table)[3] <- paste("Average_Expression_", B_choice_lab, sep = "")
      
      write_tsv(AvgExpr_Table, file)
    }
  )
  
  #download ssGSEA score table
  output$ssGSEAdownload <- downloadHandler(
    filename = function() {
      if (input$tables == 1) {
        paste('ssGSEAscore_',names(gs[(msigdb.gsea2[input$msigdbTable_rows_selected,3])]),'.tsv',sep = '')
      }
      else if (input$tables == 3) {
        paste('ssGSEAscore_',names(gs2[(GeneSet2[input$tab2table_rows_selected,1])]),'.tsv',sep = '')
      }
      else if (input$tables == 5) {
        paste('ssGSEAscore_',names(RDataListGen()[(GStable.ubg()[input$GStable.u_rows_selected,1])]),'.tsv',sep = '')
      }
    },
    content = function(file) {
      if (input$tables == 1) {
        GS <- gs[(msigdb.gsea2[input$msigdbTable_rows_selected,3])]
      }
      else if (input$tables == 3) {
        GS <- gs2[(GeneSet2[input$tab2table_rows_selected,1])]
      }
      else if (input$tables == 5) {
        GS <- RDataListGen()[(GStable.ubg()[input$GStable.u_rows_selected,1])]
      }
      ssgsea <- ssGSEAfunc()
      ssgsea2 <- as.data.frame(t(ssgsea))
      samporder <- meta[,1]
      ssgsea3 <- as.data.frame(ssgsea2[samporder,])
      colnames(ssgsea3)[1] <- colnames(ssgsea2[1])
      rownames(ssgsea3) <- samporder
      ssgsea4 <- ssgsea3 %>% 
        mutate(type = case_when(
          rownames(ssgsea3) == meta[,1] ~ meta[,2],
        ))
      ssgsea4 <- ssgsea4 %>%
        relocate(type)
      ssgsea4$sample <- rownames(ssgsea4)
      ssgsea5 <- ssgsea4[,c(3,1,2)]
      rownames(ssgsea5) <- 1:nrow(ssgsea5)
      write_tsv(ssgsea5, file)
    }
  )
  
  #render download button for DEG GMT
  output$DEGgmtDownload <- downloadHandler(
    filename = function() {
      paste(input$DEGfileName,".gmt",sep = "")
    },
    content = function(file) {
      if (ncol(meta) > 2) {
        metacol <- input$DEGtableMetaCol
      }
      else if (ncol(meta) == 2) {
        metacol <- colnames(meta)[2]
      }
      A <- meta[which(meta[,metacol] == input$comparisonA2.DEG),1]
      B <- meta[which(meta[,metacol] == input$comparisonB2.DEG),1]
      mat <- expr[,c(A,B)]
      mat <- log2(mat + 1.0)
      groupAOther <- factor(c(rep("A", length(A)), rep("B", length(B))))
      designA <- model.matrix(~0 + groupAOther)
      fit <- lmFit(mat, design = designA)
      contrast.matrix <- makeContrasts(groupAOtherA - groupAOtherB, levels = designA)
      fit2 <- contrasts.fit(fit, contrast.matrix)
      fit2 <- eBayes(fit2)
      options(digits = 4)
      top1 <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
      if (input$UpDnChoice == "UpAndDown_Regulated"){
        genes <- rownames(top1)[which(abs(top1$logFC) > abs(input$fc_cutoff2) & top1$adj.P.Val < input$p_cutoff2)]
        genes.h <- t(as.data.frame(head(genes, n=input$top_x2)))
        genes.h.gmt <- data.frame(paste(input$top_x2,input$UpDnChoice,"DEgenes", sep = ""),
                                  paste(input$UpDnChoice,"DEgenes",sep = ""),
                                  genes.h)
      }
      else if (input$UpDnChoice == "Up_Regulated"){
        genes <- rownames(top1)[which(top1$logFC > input$fc_cutoff2 & top1$adj.P.Val < input$p_cutoff2)]
        genes.h <- t(as.data.frame(head(genes, n=input$top_x2)))
        genes.h.gmt <- data.frame(paste(input$top_x2,input$UpDnChoice,"DEgenes", sep = ""),
                                  paste(input$UpDnChoice,"DEgenes",sep = ""),
                                  genes.h)
      }
      else if (input$UpDnChoice == "Down_Regulated"){
        genes <- rownames(top1)[which(top1$logFC < -abs(input$fc_cutoff2) & top1$adj.P.Val < input$p_cutoff2)]
        genes.h <- t(as.data.frame(head(genes, n=input$top_x2)))
        genes.h.gmt <- data.frame(paste(input$top_x2,input$UpDnChoice,"DEgenes", sep = ""),
                                  paste(input$UpDnChoice,"DEgenes",sep = ""),
                                  genes.h)
      }
      write_delim(genes.h.gmt, file, delim = '\t', col_names = F)
    }
  )
  
  #render download button for upreg gene pathways GMT
  output$UpRegPathDownloadgmt <- downloadHandler(
    filename = function() {
      paste("UpReg_",input$SelectedPathway,"_pathway.gmt",sep = "")
    },
    content = function(file) {
      top1 <- topgenereact()
      genes <- rownames(top1)[which(top1$P.Value < adjp & top1$logFC > FC)]
      dbs <- listEnrichrDbs() 
      enrichRLive <- TRUE 
      if (is.null(dbs)) { 
        enrichRLive <- FALSE 
      }
      dbs <- input$SelectedPathway
      enriched <- enrichr(genes, dbs) # Plot top 20 GO-BP results ordered by P-value
      df <- enriched[[1]]
      path.genes <- as.data.frame(df[,'Genes'])
      genes.n <- c()
      for (i in path.genes[,1]) {
        g <- strsplit(i,";")
        for (j in g){
          genes.n <- c(genes.n, j)
        }
      }
      rm(g)
      genes.n <- t(as.data.frame(unique(genes.n)))
      gmt.df <- data.frame(paste("UpReg_",input$SelectedPathway,"_PathwayGenes",sep = ""),
                           input$SelectedPathway,genes.n)
      rownames(gmt.df) <- NULL
      write_delim(gmt.df, file, delim = '\t', col_names = F)
    }
  )
  
  #render download button for dnreg gene pathways GMT
  output$DnRegPathDownloadgmt <- downloadHandler(
    filename = function() {
      paste("DnReg_",input$SelectedPathway,"_pathway.gmt",sep = "")
    },
    content = function(file) {
      top1 <- topgenereact()
      genes <- rownames(top1)[which(top1$P.Value < adjp & top1$logFC > FC)]
      dbs <- listEnrichrDbs() 
      enrichRLive <- TRUE 
      if (is.null(dbs)) { 
        enrichRLive <- FALSE 
      }
      dbs <- input$SelectedPathway
      enriched <- enrichr(genes, dbs) # Plot top 20 GO-BP results ordered by P-value
      df <- enriched[[1]]
      path.genes <- as.data.frame(df[,'Genes'])
      genes.n <- c()
      for (i in path.genes[,1]) {
        g <- strsplit(i,";")
        for (j in g){
          genes.n <- c(genes.n, j)
        }
      }
      rm(g)
      genes <- t(as.data.frame(unique(genes.n)))
      gmt.df <- data.frame(paste("DnReg_",input$SelectedPathway,"_PathwayGenes",sep = ""),
                           input$SelectedPathway,genes.n)
      rownames(gmt.df) <- NULL
      write_delim(gmt.df, file, delim = '\t', col_names = F)
    }
  )
  
  #render download button for upreg gene pathways
  output$UpRegPathDownload <- downloadHandler(
    filename = function() {
      paste("UpReg_",input$SelectedPathway,"_pathway.tsv",sep = "")
    },
    content = function(file) {
      top1 <- topgenereact()
      genes <- rownames(top1)[which(top1$P.Value < adjp & top1$logFC > FC)]
      dbs <- listEnrichrDbs() 
      enrichRLive <- TRUE 
      if (is.null(dbs)) { 
        enrichRLive <- FALSE 
      }
      dbs <- input$SelectedPathway
      enriched <- enrichr(genes, dbs) # Plot top 20 GO-BP results ordered by P-value
      df <- enriched[[1]]
      write_delim(df, file, delim = '\t')
    }
  )
  
  #render download button for dnreg gene pathways
  output$DnRegPathDownload <- downloadHandler(
    filename = function() {
      paste("DnReg_",input$SelectedPathway,"_pathway.tsv",sep = "")
    },
    content = function(file) {
      top1 <- topgenereact()
      genes <- rownames(top1)[which(top1$P.Value < adjp & top1$logFC > FC)]
      dbs <- listEnrichrDbs() 
      enrichRLive <- TRUE 
      if (is.null(dbs)) { 
        enrichRLive <- FALSE 
      }
      dbs <- input$SelectedPathway
      enriched <- enrichr(genes, dbs) # Plot top 20 GO-BP results ordered by P-value
      df <- enriched[[1]]
      write_delim(df, file, delim = '\t')
    }
  )
  
  #render Most Variable Gene download button
  output$MVGdownload <- downloadHandler(
    filename = function() {
      top_probes <- input$NumFeatures
      paste(top_probes,"_Most_Variable_Genes", ".tsv", sep = "")
    },
    content = function(file){
      top_probes <- input$NumFeatures
      col_labels <- colnames(expr)
      isexpr <- rowSums(expr[,col_labels] > 1) >= length(col_labels)
      exp <- expr[isexpr,]
      mad <- NULL
      var <- NULL
      cv <- NULL
      var_type <- input$VarianceMeasure
      if (var_type == "MAD"){
        mad <- apply(log2(exp + 1), 1, mad)
        mad <- sort(mad, decreasing = T)
        mad <- head(mad, n = (top_probes +1))
        out <- cbind(names(mad), mad[names(mad)], exp[names(mad),])
        colnames(out) <- c("Gene", "MAD", colnames(exp))
        dataset <- exp[names(mad),]
        variable_gene_list <- names(mad)
      }
      if (var_type == "VAR"){
        var <- apply(log2(exp + 1), 1, var)
        var <- sort(var, decreasing = T)
        var <- head(var, n = (top_probes +1))
        out <- cbind(names(var), var[names(var)], exp[names(var),])
        colnames(out) <- c("Gene", "VAR", colnames(exp))
        dataset <- exp[names(var),]
        variable_gene_list <- names(var)
      }
      if (var_type == "CV"){
        cv <- apply(log2(exp + 1), 1, cv)
        cv <- sort(cv, decreasing = T)
        cv <- head(cv, n = (top_probes +1))
        out <- cbind(names(cv), cv[names(cv)], exp[names(cv),])
        colnames(out) <- c("Gene", "VAR", colnames(exp))
        dataset <- exp[names(cv),]
        variable_gene_list <- names(cv)
      }
      variable_gene_list <- as.data.frame(variable_gene_list)
      colnames(variable_gene_list)[1] <- "Genes"
      variable_gene_list$Rank <- rownames(variable_gene_list)
      variable_gene_list <- variable_gene_list %>%
        select(Rank, Genes)
      write_delim(variable_gene_list, file, delim = '\t')
    }
  )
  
  #render Most Variable Gene GMT download button
  output$MVGdownloadgmt <- downloadHandler(
    filename = function() {
      top_probes <- input$NumFeatures
      paste(top_probes,"_Most_Variable_Genes", ".gmt", sep = "")
    },
    content = function(file){
      top_probes <- input$NumFeatures
      col_labels <- colnames(expr)
      isexpr <- rowSums(expr[,col_labels] > 1) >= length(col_labels)
      exp <- expr[isexpr,]
      mad <- NULL
      var <- NULL
      cv <- NULL
      var_type <- input$VarianceMeasure
      if (var_type == "MAD"){
        mad <- apply(log2(exp + 1), 1, mad)
        mad <- sort(mad, decreasing = T)
        mad <- head(mad, n = (top_probes +1))
        out <- cbind(names(mad), mad[names(mad)], exp[names(mad),])
        colnames(out) <- c("Gene", "MAD", colnames(exp))
        dataset <- exp[names(mad),]
        variable_gene_list <- names(mad)
      }
      if (var_type == "VAR"){
        var <- apply(log2(exp + 1), 1, var)
        var <- sort(var, decreasing = T)
        var <- head(var, n = (top_probes +1))
        out <- cbind(names(var), var[names(var)], exp[names(var),])
        colnames(out) <- c("Gene", "VAR", colnames(exp))
        dataset <- exp[names(var),]
        variable_gene_list <- names(var)
      }
      if (var_type == "CV"){
        cv <- apply(log2(exp + 1), 1, cv)
        cv <- sort(cv, decreasing = T)
        cv <- head(cv, n = (top_probes +1))
        out <- cbind(names(cv), cv[names(cv)], exp[names(cv),])
        colnames(out) <- c("Gene", "VAR", colnames(exp))
        dataset <- exp[names(cv),]
        variable_gene_list <- names(cv)
      }
      variable_gene_list <- as.data.frame(variable_gene_list)
      colnames(variable_gene_list)[1] <- "Genes"
      genes.n <- c()
      for (i in variable_gene_list[,1]) {
        g <- strsplit(i,";")
        for (j in g){
          genes.n <- c(genes.n, j)
        }
      }
      rm(g)
      genes <- t(as.data.frame(unique(genes.n)))
      gmt.df <- data.frame(paste(top_probes,"_Most_Variable_Genes", sep = ""),
                           "Most_Variable_Genes",genes)
      rownames(gmt.df) <- NULL
      write_delim(gmt.df, file, delim = '\t', col_names = F)
    }
  )
  
  #download button for leading edge genes
  output$LEGdownload <- downloadHandler(
    filename = function() {
      if (input$tables == 1){
        if (length(input$msigdbTable_rows_selected) > 0){
          GS <- as.character(msigdb.gsea2[input$msigdbTable_rows_selected,3])
        }
      }
      else if (input$tables == 3){
        if (length(input$tab2table_rows_selected) > 0){
          GS <- as.character(GeneSet2[input$tab2table_rows_selected,1])
        }
      }
      else if (input$tables == 5){
        if (length(input$GStable.u_rows_selected) > 0){
          GS <- as.character(user_gs_mirror()[input$GStable.u_rows_selected,1])
        }
      }
      paste(GS,"_leading_edge_genes", ".tsv", sep = "")
    },
    content = function(file){
      if (input$tables == 1){
        if (length(input$msigdbTable_rows_selected) > 0){
          GS <- as.character(msigdb.gsea2[input$msigdbTable_rows_selected,3])
        }
      }
      else if (input$tables == 3){
        if (length(input$tab2table_rows_selected) > 0){
          GS <- as.character(GeneSet2[input$tab2table_rows_selected,1])
        }
      }
      else if (input$tables == 5){
        if (length(input$GStable.u_rows_selected) > 0){
          GS <- as.character(user_gs_mirror()[input$GStable.u_rows_selected,1])
        }
      }
      res <- datasetInput()
      gsea.df <- as.data.frame(res@result)
      ## Subset core enriched genes
      genes1 <- as.matrix(gsea.df[which(gsea.df$Description==GS),"core_enrichment"])
      genes2 <- strsplit(genes1,"/")
      GeneSymbol <- as.data.frame(genes2, col.names = "GeneSymbol")
      GeneSymbol$Rank <- rownames(GeneSymbol)
      GeneSymbol <- GeneSymbol[,c("Rank","GeneSymbol")]
      write_delim(GeneSymbol, file, delim = '\t')
    })
  
  #render DEG table download button
  output$DEGtableDownload <- downloadHandler(
    filename = function() {
      paste("DEG_Table_",Sys.Date(), ".tsv", sep = "")
    },
    content = function(file){
      if (ncol(meta) > 2) {
        metacol <- input$DEGtableMetaCol
      }
      else if (ncol(meta) == 2) {
        metacol <- colnames(meta)[2]
      }
      A <- meta[which(meta[,metacol] == input$comparisonA2.DEG),1]
      B <- meta[which(meta[,metacol] == input$comparisonB2.DEG),1]
      mat <- expr[,c(A,B)]
      mat <- log2(mat + 1.0)
      groupAOther <- factor(c(rep("A", length(A)), rep("B", length(B))))
      designA <- model.matrix(~0 + groupAOther)
      fit <- lmFit(mat, design = designA)
      contrast.matrix <- makeContrasts(groupAOtherA - groupAOtherB, levels = designA)
      fit2 <- contrasts.fit(fit, contrast.matrix)
      fit2 <- eBayes(fit2)
      options(digits = 4)
      top1 <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
      top1$Gene <- rownames(top1)
      top1 <- top1 %>%
        relocate(Gene)
      write_delim(top1, file, delim = '\t')
    }
  )
  
  #download button for Enrich Sig Table
  output$enrich_sig_download <- downloadHandler(
    filename = function() {
      groupA <- input$comparisonA
      groupB <- input$comparisonB
      paste(input$SigTableChoice,".tsv", sep = "")
    },
    content = function(file){
      gsea.df <- as_tibble(get(ES_Tab_List[which(ES_Tab_List == paste("ES_table",match(input$SigTableChoice, SigNames),sep = ""))]))
      write_delim(gsea.df, file, delim = '\t')
    })
  
  #download button for user enriched signature table
  output$enrich_sig_download.u <- downloadHandler(
    filename = function() {
      groupA <- input$comparisonA
      groupB <- input$comparisonB
      paste("Enrich_Sig_Table_",groupA,"vs",groupB,".tsv", sep = "")
    },
    content = function(file) {
      if (input$tables == 3) {
        groupA <- meta[which(meta[,2] == input$comparisonA),1]
        groupB <- meta[which(meta[,2] == input$comparisonB),1]
        ##----Signal-to-Noise Calculation----##
        A <- A + 0.00000001
        P = as.matrix(as.numeric(colnames(A) %in% groupA))
        n1 <- sum(P[,1])
        M1 <- A %*% P
        M1 <- M1/n1
        A2 <- A*A
        S1 <- A2 %*% P
        S1 <- S1/n1 - M1*M1 
        S1 <- sqrt(abs((n1/(n1-1)) * S1))
        P = as.matrix(as.numeric(colnames(A) %in% groupB))
        n2 <- sum(P[,1])
        M2 <- A %*% P
        M2 <- M2/n2
        A2 <- A*A
        S2 <- A2 %*% P
        S2 <- S2/n2 - M2*M2
        S2 <- sqrt(abs((n2/(n2-1)) * S2))
        rm(A2)
        # small sigma "fix" as used in GeneCluster
        S2 <- ifelse(0.2*abs(M2) < S2, S2, 0.2*abs(M2))
        S2 <- ifelse(S2 == 0, 0.2, S2)
        S1 <- ifelse(0.2*abs(M1) < S1, S1, 0.2*abs(M1))
        S1 <- ifelse(S1 == 0, 0.2, S1)
        M1 <- M1 - M2
        rm(M2)
        S1 <- S1 + S2
        rm(S2)
        s2n.matrix <- M1/S1
        ##----Reformatting----##
        s2n.df <- as.data.frame(s2n.matrix)
        s2n.df$GeneID <- rownames(s2n.df)
        rownames(s2n.df) <- NULL
        data <- dplyr::select(s2n.df, GeneID, V1)
        data.gsea <- data$V1
        names(data.gsea) <- as.character(data$GeneID)
        s2n.matrix.s <- sort(data.gsea, decreasing = T)
        ##----GSEA----##
        gmt.i <- tab2
        gsea.res <- GSEA(s2n.matrix.s, TERM2GENE = gmt.i, verbose = F, pvalueCutoff = input$userPval)
        gsea.df <- as.data.frame(gsea.res@result)
        write_delim(gsea.df, file, delim = '\t')
      }
      else if (input$tables == 5) {
        groupA <- meta[which(meta[,2] == input$comparisonA),1]
        groupB <- meta[which(meta[,2] == input$comparisonB),1]
        ##----Signal-to-Noise Calculation----##
        A <- A + 0.00000001
        P = as.matrix(as.numeric(colnames(A) %in% groupA))
        n1 <- sum(P[,1])
        M1 <- A %*% P
        M1 <- M1/n1
        A2 <- A*A
        S1 <- A2 %*% P
        S1 <- S1/n1 - M1*M1 
        S1 <- sqrt(abs((n1/(n1-1)) * S1))
        P = as.matrix(as.numeric(colnames(A) %in% groupB))
        n2 <- sum(P[,1])
        M2 <- A %*% P
        M2 <- M2/n2
        A2 <- A*A
        S2 <- A2 %*% P
        S2 <- S2/n2 - M2*M2
        S2 <- sqrt(abs((n2/(n2-1)) * S2))
        rm(A2)
        # small sigma "fix" as used in GeneCluster
        S2 <- ifelse(0.2*abs(M2) < S2, S2, 0.2*abs(M2))
        S2 <- ifelse(S2 == 0, 0.2, S2)
        S1 <- ifelse(0.2*abs(M1) < S1, S1, 0.2*abs(M1))
        S1 <- ifelse(S1 == 0, 0.2, S1)
        M1 <- M1 - M2
        rm(M2)
        S1 <- S1 + S2
        rm(S2)
        s2n.matrix <- M1/S1
        ##----Reformatting----##
        s2n.df <- as.data.frame(s2n.matrix)
        s2n.df$GeneID <- rownames(s2n.df)
        rownames(s2n.df) <- NULL
        data <- dplyr::select(s2n.df, GeneID, V1)
        data.gsea <- data$V1
        names(data.gsea) <- as.character(data$GeneID)
        s2n.matrix.s <- sort(data.gsea, decreasing = T)
        ##----GSEA----##
        gmt.i <- GStable.ubg()
        gsea.res <- GSEA(s2n.matrix.s, TERM2GENE = gmt.i, verbose = F, pvalueCutoff = input$userPval)
        gsea.df <- as.data.frame(gsea.res@result)
        write_delim(gsea.df, file, delim = '\t')
      }
      else if (input$tables == 1){
        if (input$GenerateEST == TRUE) {
          gsea.df <- GeneratedMSigDBEST()
          gsea.df <- as.data.frame(gsea.df)
          write_delim(gsea.df, file, delim = '\t')
        }
      }
    }
  )
  
  
  ####----Text----####
  
  
  #NES and Pval output
  output$NESandPval <- renderText({
    if (input$tables == 1){
      if (length(input$msigdbTable_rows_selected) > 0){
        res <- datasetInput()
        gsea.df <- as.data.frame(res@result)
        GS = as.character(msigdb.gsea2[input$msigdbTable_rows_selected,3])
        NES = gsea.df$NES[which(gsea.df[,'ID']==GS)]
        Pval = gsea.df$pvalue[which(gsea.df[,'ID']==GS)]
        NES.o <- paste0("NES: ", NES)
        Pval.o <- paste0("Pvalue: ", Pval)
        if (NES > 0){
          UpOrDown <- paste("Based on the normalized enrichment score above, the", GS, "gene set is upregulated in", input$comparisonA , "group.")
        }
        else {
          UpOrDown <- paste("Based on the normalized enrichment score above, the", GS, "gene set is downregulated in", input$comparisonA , "group.")
        }
        paste(NES.o, Pval.o, UpOrDown, sep = '\n')
      }
      else if (length(input$msigdbTable_rows_selected) == 0){
        paste("Please select gene set from side panel table to begin.", sep = '')
      }
    }
    else if (input$tables == 3){
      if (length(input$tab2table_rows_selected) > 0){
        res <- datasetInput()
        gsea.df <- as.data.frame(res@result)
        GS = as.character(GeneSet2[input$tab2table_rows_selected,1])
        NES = gsea.df$NES[which(gsea.df[,'ID']==GS)]
        Pval = gsea.df$pvalue[which(gsea.df[,'ID']==GS)]
        NES.o <- paste0("NES: ", NES)
        Pval.o <- paste0("Pvalue: ", Pval)
        if (NES > 0){
          UpOrDown <- paste("Based on the normalized enrichment score above, the", GS, "gene set is upregulated in", input$comparisonA , "group.")
        }
        else {
          UpOrDown <- paste("Based on the normalized enrichment score above, the", GS, "gene set is downregulated in", input$comparisonA , "group.")
        }
        paste(NES.o, Pval.o, UpOrDown, sep = '\n')
      }
      else if (length(input$tab2table_rows_selected) == 0){
        paste("Please select gene set from side panel table to begin.", sep = '')
      }
    }
    else if (input$tables == 5){
      if (length(input$GStable.u_rows_selected) > 0){
        res <- datasetInput()
        gsea.df <- as.data.frame(res@result)
        GS <- as.character(user_gs_mirror()[input$GStable.u_rows_selected,1])
        NES = gsea.df$NES[which(gsea.df[,'ID']==GS)]
        Pval = gsea.df$pvalue[which(gsea.df[,'ID']==GS)]
        NES.o <- paste0("NES: ", NES)
        Pval.o <- paste0("Pvalue: ", Pval)
        if (NES > 0){
          UpOrDown <- paste("Based on the normalized enrichment score above, the", GS, "gene set is upregulated in", input$comparisonA , "group.")
        }
        else {
          UpOrDown <- paste("Based on the normalized enrichment score above, the", GS, "gene set is downregulated in", input$comparisonA , "group.")
        }
        paste(NES.o, Pval.o, UpOrDown, sep = '\n')
      }
      else if (length(input$GStable.u_rows_selected) == 0){
        paste("Please select gene set from side panel table to begin.", sep = '')
      }
    }
  })
  
  output$GenesAboveCutoff1 <- renderText({
    
    # UI Inputs
    A_choice <- input$comparisonA2_h            #Comparison group A
    B_choice <- input$comparisonB2_h            #Comparison group B
    FC_cutoff <- input$fc_cutoff_h              #FC cutoff for top gene selection 
    P_cutoff <- input$p_cutoff_h                #P-value cutoff for top gene selections
    
    # Make Top table - Compares only 2 groups
    #make group based on user input
    A <- meta[which(meta[,2] == A_choice),1]
    B <- meta[which(meta[,2] == B_choice),1]
    mat <- expr[,c(A,B)]
    mat <- log2(mat + 1.0)
    groupAOther <- factor(c(rep("A", length(A)), rep("B", length(B))))
    designA <- model.matrix(~0 + groupAOther)
    fit <- lmFit(mat, design = designA)
    contrast.matrix <- makeContrasts(groupAOtherA - groupAOtherB, levels = designA)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    options(digits = 4)
    top1 <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
    
    top_above_cutoff <- top1[which(abs(top1$logFC) >= abs(FC_cutoff) & top1$P.Value <= P_cutoff),]
    
    Num_Genes_Above_Cutoff1 <- length(rownames(top_above_cutoff))
    
    paste("Number of Genes Above FC and P.Value Cutoffs: ",Num_Genes_Above_Cutoff1,sep = "")
    
  })
  
  output$GenesAboveCutoff2 <- renderText({
    
    # UI Inputs
    A_choice <- input$comparisonA2_ha            #Comparison group A
    B_choice <- input$comparisonB2_ha            #Comparison group B
    FC_cutoff <- input$fc_cutoff_ha              #FC cutoff for top gene selection 
    P_cutoff <- input$p_cutoff_ha                #P-value cutoff for top gene selections
    
    # Make Top table - Compares only 2 groups
    #make group based on user input
    A <- meta[which(meta[,2] == A_choice),1]
    B <- meta[which(meta[,2] == B_choice),1]
    mat <- expr[,c(A,B)]
    mat <- log2(mat + 1.0)
    groupAOther <- factor(c(rep("A", length(A)), rep("B", length(B))))
    designA <- model.matrix(~0 + groupAOther)
    fit <- lmFit(mat, design = designA)
    contrast.matrix <- makeContrasts(groupAOtherA - groupAOtherB, levels = designA)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    options(digits = 4)
    top1 <- topTable(fit2, coef = 1, n = 300000, sort = "p", p.value = 1.0, adjust.method = "BH")
    
    top_above_cutoff <- top1[which(abs(top1$logFC) >= abs(FC_cutoff) & top1$P.Value <= P_cutoff),]
    
    Num_Genes_Above_Cutoff2 <- length(rownames(top_above_cutoff))
    
    paste("Number of Genes Above FC and P.Value Cutoffs: ",Num_Genes_Above_Cutoff2,sep = "")
    
  })
  
  output$VolGroupsText <- renderText({
    paste("This volcano plot is comparing group A: ",input$comparisonA2, " and group B: ",input$comparisonB2, ".\nGenes with a positive log fold change are upregulated in the ",
          input$comparisonA2, " group.\nGenes with a negative log fold change are upregulated in the ",input$comparisonB2, " group.", sep = "")
  })
  
  output$MAGroupsText <- renderText({
    paste("This MA plot is comparing group A: ",input$comparisonA2, " and group B: ",input$comparisonB2, ".\nGenes with a positive log fold change are upregulated in the ",
          input$comparisonA2, " group.\nGenes with a negative log fold change are upregulated in the ",input$comparisonB2, " group.", sep = "")
  })
  
  output$upregpath_text <- renderText({
    paste("Genes in these enriched terms are upregulated in group A: ", input$comparisonA2.path, " group.", sep = "")
  })
  
  output$downregpath_text <- renderText({
    paste("Genes in these enriched terms are upregulated in group B: ", input$comparisonB2.path," group", sep = "")
  })
  
  output$degtext <- renderText({
    paste("This table represents differentially expressed genes when comparing group A: ",input$comparisonA2.DEG," and group B: ",input$comparisonB2.DEG,
          ".\nGenes with a positive logFC, indicate an upregulation in group A: ", input$comparisonA2.DEG,
          ".\nGenes with a negative logFC, indicate an upregulation in group B: ", input$comparisonB2.DEG,
          ".\nThe 'AveExpr' column represents the log transformed average expression between group A: ",input$comparisonA2.DEG," and group B: ",input$comparisonB2.DEG,".", sep = "")
  })
  
  output$UpRegPathLabel <- renderUI({
    
    FC <- input$pathFC
    h3(paste("(Up-regulated pathway (> ",FC," logFC)",sep = ""))
    
  })
  
  output$DnRegPathLabel <- renderUI({
    
    FC <- input$pathFC
    h3(paste("Down-regulated pathway (> -",FC," logFC)",sep = ""))
    
  })
}


# Run the application 
shinyApp(ui = ui, server = server)