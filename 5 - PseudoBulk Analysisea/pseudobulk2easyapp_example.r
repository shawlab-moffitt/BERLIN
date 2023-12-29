#### Installation ####
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("celldex")

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SingleR")


devtools::install_github("immunogenomics/presto")

options(future.globals.maxSize = 8000 * 1024^2)

####---- User Input ----####

## Set local Github repository as working directory
setwd("/Users/4472414/Projects/ADIOS_Immunology/GSE123813_PD1_TumorSpecificTcells")

## Desired project name
Project_Name <- "GSE123813_bcc"

## Define organism

organism <- "Human"
# Define the path to the count file
Count_file <- "/Users/4472414/Projects/ADIOS_Immunology/GSE123813_PD1_TumorSpecificTcells/GSE123813_bcc_scRNA_counts.txt"

# Define path to clinical/meta data (OPTIONAL)
meta_file <- "/Users/4472414/Projects/ADIOS_Immunology/GSE123813_PD1_TumorSpecificTcells/GSE123813_bcc_tcell_metadata.txt"

# Specify the output folder path
output_folder <- "/Users/4472414/Projects/ADIOS_Immunology/GSE123813_PD1_TumorSpecificTcells/GSE123813_bcc_output"

seed <- 42

resolution <- c(0.5,1,1.5,2)


####---- Load Packages ----####

packages <- c("dplyr","Seurat","patchwork","readr","DoubletFinder","celldex","Matrix",
              "fields","SeuratDisk","SeuratData","SingleR","SingleCellExperiment")
invisible(lapply(packages, library, character.only = TRUE))




####---- Run Script ----####

set.seed(seed)

# Create an output directory if it doesn't already exist
if (!file.exists(output_folder)) {
  dir.create(output_folder)
}

if (file.exists(meta_file)) {
  meta <- as.data.frame(read_delim(meta_file, delim = '\t', col_names = T))
  meta[,1] <- gsub("[[:punct:]]",".",meta[,1])
  colnames(meta)[1] <- "Cell"
} else {meta <- NULL}

# seurat_h5_unprocessed
seurat_unprocessed_h5_path = file.path(output_folder, paste0(Project_Name,"_unprocessed_h5friendly"));

if (file.exists(file)) {
  # Read the count file
  counts <- read.csv2(Count_file, sep = "", header = TRUE, row.names = 1)
  colnames(counts) <- gsub("[[:punct:]]",".",colnames(counts))

  # Define a function for quality control (QC) using Seurat
  qc.seurat <- function(seurat, species, nFeature) {
    mt.pattern <- case_when(
      species == "Human" ~ "^MT-",
      species == "Mouse" ~ "^mt-",
      TRUE ~ "^MT-"
    )
    ribo.pattern <- case_when(
      species == "Human" ~ "^RP[LS]",
      species == "Mouse" ~ "^Rp[ls]",
      TRUE ~ "^RP[LS]"
    )

    # Calculate percentage of mitochondrial and ribosomal genes
    seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = mt.pattern, assay = "RNA")
    seurat[["percent.rp"]] <- PercentageFeatureSet(seurat, pattern = ribo.pattern, assay = "RNA")

    # Filter cells based on QC criteria
    seurat[, seurat[["percent.mt"]] <= 20 & seurat[["nFeature_RNA"]] >= nFeature]
  }

  # Create a Seurat object and perform QC
  seurat_obj <- CreateSeuratObject(counts = counts)
  orig_seurat_obj = seurat_obj #
  #seurat_obj = orig_seurat_obj
  #head(seurat_obj@meta.data)

  meta_out <- seurat_obj@meta.data
  meta_out$Cell <- rownames(meta_out)
  meta_out <- meta_out %>% relocate(Cell)
  if (!is.null(meta)) {
    meta_out <- merge(meta,meta_out, by = "Cell", all.y = T)
    meta_obj <- meta_out
    rownames(meta_obj) <- meta_obj$Cell
    meta_obj <- meta_obj[,-1]
    seurat_obj@meta.data <- meta_obj
  }


  head(seurat_obj@meta.data)



  # Write seurat H5 file



} else {
  seurat_obj <- LoadH5Seurat(seurat_unprocessed_h5_file)

}

seurat_obj@meta.data = seurat_obj@meta.data[which(seurat_obj@meta.data$patient != "NA"),]
seurat_obj@meta.data = seurat_obj@meta.data[which(seurat_obj@meta.data$cluster != "NA"),]
seurat_obj@meta.data = seurat_obj@meta.data[which(seurat_obj@meta.data$treatment != "NA"),]
Idents(object=seurat_obj) <- "orig.ident";

SaveH5Seurat(seurat_obj, filename = seurat_unprocessed_h5_path,
             overwrite = TRUE, verbose = TRUE)
# seurat_obj <- qc.seurat(seurat_obj, "Human", 500) # not sure if this is necesary

# Normalize data using LogNormalize method
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

# Calculate cell cycle scores
seurat_obj <- CellCycleScoring(object = seurat_obj, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)

# Find highly variable features
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)


# Scale data by regressing out unwanted sources of variation
seurat_obj <- ScaleData(seurat_obj, vars.to.regress = c("nFeature_RNA", "percent.mt"), verbose = FALSE)

head(seurat_obj@meta.data)
agg = AggregateExpression(seurat_obj, return.seurat = T, group.by = c('patient', 'treatment', 'cluster'), normalization.method = "LogNormalize",scale.factor = 10000)
head(agg@assays$RNA$data)

cnames = colnames(agg@assays$RNA$data)
m1 = do.call(rbind, strsplit(cnames, "_"))

meta_text = cbind(cnames, m1)
colnames(meta_text) = c("ID", "Patient", "Treatment", "CellType")
head(meta_text)
treat_celltype = paste(meta_text[, "Treatment"], meta_text[,"CellType"], sep="_")
meta_text = cbind(meta_text, treat_celltype);
colnames(meta_text)
head(meta_text)

matrix_file = paste(output_folder, Project_Name, "_logNormalized_data.txt", sep="");
write.table(2^ agg@assays$RNA$data, file = matrix_file, sep="\t", row.names = TRUE,col.names=NA)
head(agg@assays$RNA$data)
meta_file = paste(output_folder, Project_Name, "_meta.txt", sep="");
write.table(meta_text, file = meta_file, sep="\t", row.names=F)

hist(agg@assays$RNA$data[,1], breaks=100)



### Write out the matrix ###
outputmatrix = matrix_file
outputmeta = meta_file

# generation of the EASY app
# It'll be a good idea to double check the files now before proceeding with generating the apps
UPPERCASE_organism = toupper(organism)
if (UPPERCASE_organism == "HUMAN") {

  output_dir <- file.path(output_folder, paste(Project_Name, "_Human_EASY_App", sep=""))
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
  } else {
    #print("Dir already exists!")
  }
  R_script_output_dir <- file.path(output_folder, paste(Project_Name, "_Human_EASY_App/R", sep=""))

  if (!dir.exists(R_script_output_dir)){
    dir.create(R_script_output_dir)
  } else {
    #print("Dir already exists!")
  }


  outputmatrix = paste(output_dir, "/", Project_Name, "_logNormalized_data.txt", sep="")
  outputmeta = paste(output_dir, "/", Project_Name, "_meta.txt", sep="")

  write.table(2^ agg@assays$RNA$data, file = outputmatrix, sep="\t", row.names = TRUE,col.names=NA)
  write.table(meta_text, file=outputmeta, row.names=FALSE, sep="\t")

  hs_LINC1000_url = "https://raw.githubusercontent.com/shawlab-moffitt/GEO_DataDownload_Example/e065189a65f861024ac10f683f1595d44f479ffa/EASYAppFiles/Template_Human_EASY_App/LINCS_L1000_gsNsym_HS_v2.zip";
  hs_linc1000_rdata_url = "https://github.com/shawlab-moffitt/GEO_DataDownload_Example/raw/main/EASYAppFiles/Template_Human_EASY_App/LINCS_L1000_gs_HS_v2.RData";
  hs_app_url = "https://raw.githubusercontent.com/shawlab-moffitt/GEO_DataDownload_Example/main/EASYAppFiles/Template_Human_EASY_App/app.R";
  hs_msigdb_url = "https://raw.githubusercontent.com/shawlab-moffitt/GEO_DataDownload_Example/main/EASYAppFiles/Template_Human_EASY_App/msigdb_gsNcat_HS_v2.txt"
  hs_msigdb_zip_url = "https://github.com/shawlab-moffitt/GEO_DataDownload_Example/raw/main/EASYAppFiles/Template_Human_EASY_App/msigdb_gsNsym_HS_v2.zip"
  hs_msigdb_rdata_url = "https://github.com/shawlab-moffitt/GEO_DataDownload_Example/raw/main/EASYAppFiles/Template_Human_EASY_App/msigdb_gs_HS_v2.RData"

  internalR_url = "https://raw.githubusercontent.com/shawlab-moffitt/GEO_DataDownload_Example/main/EASYAppFiles/Template_Human_EASY_App/R/internal.R";
  loginR_url = "https://raw.githubusercontent.com/shawlab-moffitt/GEO_DataDownload_Example/main/EASYAppFiles/Template_Human_EASY_App/R/login.R"
  logoutR_url = "https://raw.githubusercontent.com/shawlab-moffitt/GEO_DataDownload_Example/main/EASYAppFiles/Template_Human_EASY_App/R/logout.R"
  runExampleR_url = "https://raw.githubusercontent.com/shawlab-moffitt/GEO_DataDownload_Example/main/EASYAppFiles/Template_Human_EASY_App/R/runExample.R"


  output_file_LINC1000 = paste(output_dir, "/LINCS_L1000_gsNsym_HS_v2.zip", sep="")
  output_file_rdata = paste(output_dir, "/LINCS_L1000_gs_HS_v2.RData", sep="")
  output_file_app = paste(output_dir, "/app.R", sep="")
  output_file_msigdb = paste(output_dir, "/msigdb_gsNcat_HS_v2.txt", sep="")
  output_file_msigdbzip = paste(output_dir, "/msigdb_gsNsym_HS_v2.zip", sep="")
  output_file_msigdbrdata = paste(output_dir, "/msigdb_gs_HS_v2.RData", sep="")
  output_file_internalR = paste(R_script_output_dir, "/internal.R", sep="")
  output_file_loginR = paste(R_script_output_dir, "/login.R", sep="")
  output_file_logoutR = paste(R_script_output_dir, "/logout.R", sep="")
  output_file_runExampleR = paste(R_script_output_dir, "/runExample.R", sep="")


  download.file(hs_LINC1000_url, output_file_LINC1000)
  download.file(hs_linc1000_rdata_url, output_file_rdata)
  download.file(hs_app_url, output_file_app);
  download.file(hs_msigdb_url, output_file_msigdb);
  download.file(hs_msigdb_zip_url, output_file_msigdbzip);
  download.file(hs_msigdb_rdata_url, output_file_msigdbrdata);

  download.file(internalR_url, output_file_internalR);
  download.file(loginR_url, output_file_loginR);
  download.file(logoutR_url, output_file_logoutR);
  download.file(runExampleR_url, output_file_runExampleR);

  # replace the expr_and meta file
  tx  <- readLines(output_file_app)
  tx2 <- gsub(pattern = "expr_file <- ''", replace = paste("expr_file <- '", Project_Name, "_logNormalized_data.txt", "'", sep=""), x = tx)
  tx3 <- gsub(pattern = "meta_file <- ''", replace = paste("meta_file <- '", Project_Name, "_meta.txt", "'", sep=""), x = tx2)
  tx4 <- gsub(pattern = "ProjectName <- ''", replace = paste("ProjectName <- '", Project_Name, " Study", "'", sep=""), x = tx3)
  writeLines(tx4, con=output_file_app)


}
# write mouse samples
if (UPPERCASE_organism == "MOUSE") {

  output_dir <- file.path(outputpath, paste(Project_Name, "_Mouse_EASY_App", sep=""))
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
  } else {
    #print("Dir already exists!")
  }
  R_script_output_dir <- file.path(outputpath, paste(Project_Name, "_Mouse_EASY_App/R", sep=""))

  if (!dir.exists(R_script_output_dir)){
    dir.create(R_script_output_dir)
  } else {
    #print("Dir already exists!")
  }



  outputmatrix = paste(output_dir, "/", Project_Name, "_logNormalized_data.txt", sep="")
  outputmeta = paste(output_dir, "/", Project_Name, "_meta.txt", sep="")

  write.table(2^ agg@assays$RNA$data, file = outputmatrix, sep="\t", row.names = TRUE,col.names=NA)
  write.table(meta_text, file=outputmeta, row.names=FALSE, sep="\t")

  CellMarker_GS_MM_rdata_url = "https://github.com/shawlab-moffitt/GEO_DataDownload_Example/raw/main/EASYAppFiles/Template_Mouse_EASY_App/CellMarker_GS_MM.RData";
  CellMarker_gsNsym_MM_tsv_url = "https://raw.githubusercontent.com/shawlab-moffitt/GEO_DataDownload_Example/main/EASYAppFiles/Template_Mouse_EASY_App/CellMarker_gsNsym_MM.tsv";
  mm_app_url = "https://raw.githubusercontent.com/shawlab-moffitt/GEO_DataDownload_Example/main/EASYAppFiles/Template_Mouse_EASY_App/app.R";
  mm_msigdb_url = "https://raw.githubusercontent.com/shawlab-moffitt/GEO_DataDownload_Example/main/EASYAppFiles/Template_Mouse_EASY_App/msigdb_gsNcat_MM.tsv"
  mm_msigdb_zip_url = "https://github.com/shawlab-moffitt/GEO_DataDownload_Example/raw/main/EASYAppFiles/Template_Mouse_EASY_App/msigdb_gsNsym_MM.zip"
  mm_msigdb_rdata_url = "https://github.com/shawlab-moffitt/GEO_DataDownload_Example/raw/main/EASYAppFiles/Template_Mouse_EASY_App/msigdb_gs_MM.RData"

  internalR_url = "https://raw.githubusercontent.com/shawlab-moffitt/GEO_DataDownload_Example/main/EASYAppFiles/Template_Mouse_EASY_App/R/internal.R";
  loginR_url = "https://raw.githubusercontent.com/shawlab-moffitt/GEO_DataDownload_Example/main/EASYAppFiles/Template_Mouse_EASY_App/R/login.R"
  logoutR_url = "https://raw.githubusercontent.com/shawlab-moffitt/GEO_DataDownload_Example/main/EASYAppFiles/Template_Mouse_EASY_App/R/logout.R"
  runExampleR_url = "https://raw.githubusercontent.com/shawlab-moffitt/GEO_DataDownload_Example/main/EASYAppFiles/Template_Mouse_EASY_App/R/runExample.R"


  output_file_CellMarker_rdata = paste(output_dir, "/CellMarker_GS_MM.RData", sep="")
  output_file_CellMarker_tsv = paste(output_dir, "/CellMarker_gsNsym_MM.tsv", sep="")
  output_file_app = paste(output_dir, "/app.R", sep="")
  output_file_msigdb = paste(output_dir, "/msigdb_gsNcat_MM.tsv", sep="")
  output_file_msigdbzip = paste(output_dir, "/msigdb_gsNsym_MM.zip", sep="")
  output_file_msigdbrdata = paste(output_dir, "/msigdb_gs_MM.RData", sep="")
  output_file_internalR = paste(R_script_output_dir, "/internal.R", sep="")
  output_file_loginR = paste(R_script_output_dir, "/login.R", sep="")
  output_file_logoutR = paste(R_script_output_dir, "/logout.R", sep="")
  output_file_runExampleR = paste(R_script_output_dir, "/runExample.R", sep="")


  download.file(CellMarker_GS_MM_rdata_url, output_file_CellMarker_rdata)
  download.file(CellMarker_gsNsym_MM_tsv_url, output_file_CellMarker_tsv)
  download.file(mm_app_url, output_file_app);
  download.file(mm_msigdb_url, output_file_msigdb);
  download.file(mm_msigdb_zip_url, output_file_msigdbzip);
  download.file(mm_msigdb_rdata_url, output_file_msigdbrdata);

  download.file(internalR_url, output_file_internalR);
  download.file(loginR_url, output_file_loginR);
  download.file(logoutR_url, output_file_logoutR);
  download.file(runExampleR_url, output_file_runExampleR);

  # replace the expr_and meta file
  tx  <- readLines(output_file_app)
  tx2 <- gsub(pattern = "expr_file <- ''", replace = paste("expr_file <- '", Project_Name, "_logNormalized_data.txt", "'", sep=""), x = tx)
  tx3 <- gsub(pattern = "meta_file <- ''", replace = paste("meta_file <- '", Project_Name, "_meta.txt", "'", sep=""), x = tx2)
  tx4 <- gsub(pattern = "ProjectName <- ''", replace = paste("ProjectName <- '", Project_Name, " Study", "'", sep=""), x = tx3)
  writeLines(tx4, con=output_file_app)


}

