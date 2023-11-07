


## Set local Github repository as working directory
setwd("~/R/BERLIN")

## Desired project name
Project_Name <- "GSE116256_AMLscRNA"

## Specify the output folder path
output_folder <- "3-R_Shiny_Viz_Applications/Pathway_Connectivity_Shiny_App/1-Prep_Input_Files/Output/"

## Path to Single Cell Cluster Gene set derived from "2-Single_Cell_Post_Processing.R"
sc_Pipeline_GeneSet_File <- "3-R_Shiny_Viz_Applications/Pathway_Connectivity_Shiny_App/1-Prep_Input_Files/Input/GSE116256_AMLscRNA_SingleCell_Cluster_DEG_GeneSet.txt"

## Path to provided TCGA Gene Set
Geneset_2_file <- "3-R_Shiny_Viz_Applications/Pathway_Connectivity_Shiny_App/1-Prep_Input_Files/Input/PMID30827681_AML_Lineages.tsv"


####----Load Package----####
library(readr)


####----Read in Files----####
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

if (tools::file_ext(sc_Pipeline_GeneSet_File) %in% c("txt","tsv","zip","gz")) {
  sc_Pipeline_GeneSet <- as.data.frame(read_delim(sc_Pipeline_GeneSet_File,delim = '\t', col_names = T))
  text1 <- TRUE
} else if (tools::file_ext(sc_Pipeline_GeneSet_File) %in% "RData") {
  sc_Pipeline_GeneSet <- loadRData(sc_Pipeline_GeneSet_File)
  text1 <- FALSE
} else {
  print("Please check file name/format of single cell pipeline geneset")
}

if (tools::file_ext(Geneset_2_file) %in% c("txt","tsv","zip","gz")) {
  Geneset_2 <- as.data.frame(read_delim(Geneset_2_file,delim = '\t', col_names = T))
  colnames(Geneset_2) <- c("term","gene")
  text2 <- TRUE
} else if (tools::file_ext(Geneset_2_file) %in% "RData") {
  Geneset_2 <- loadRData(Geneset_2_file)
  text2 <- FALSE
} else {
  print("Please check file name/format of geneset 2")
}


####----Combine Gene Sets----####

if (!text1 & !text2) {
  New_GeneSet <- c(sc_Pipeline_GeneSet,Geneset_2)
  save(New_GeneSet, file = paste0(output_folder,"/",Project_Name,"_SingleCellCluster_withTCGA_GeneSet.RData"))
}
if (text1 & text2) {
  New_GeneSet <- rbind(sc_Pipeline_GeneSet,Geneset_2)
  write_delim(New_GeneSet, paste0(output_folder,"/",Project_Name,"_SingleCellCluster_withTCGA_GeneSet.txt"), delim = '\t')
}
if (text1 != text2) {
  print("Please provide geneset files of the same format")
}





