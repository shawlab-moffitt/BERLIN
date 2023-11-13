


####---- User Input ----####

## Set local Github repository as working directory
#setwd("C:/Users/Administrator/Desktop/BERLIN_02")

## Desired project name
Project_Name <- "GSE116256_AMLscRNA"

# Define the path to the count file
#Count_file <- "2-Single_Cell_RNAseq_Pipeline/Input/Count_File_Example_GSE116256_AML921A-D0.txt.gz"
Count_file <- "Input/Count_File_Example_GSE116256_AML921A-D0.txt.gz"

# Define path to clinical/meta data (OPTIONAL)
#meta_file <- "2-Single_Cell_RNAseq_Pipeline/Input/Meta_File_Example_GSE116256_AML921A-D0.txt.gz"
meta_file <- "Input/Meta_File_Example_GSE116256_AML921A-D0.txt.gz"


# Specify the output folder path
#output_folder <- "2-Single_Cell_RNAseq_Pipeline/Output/Single_Cell_RNAseq_Output/"
output_folder <- "Output/Single_Cell_RNAseq_Output"

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

seurat_obj <- qc.seurat(seurat_obj, "Human", 500)

# Normalize data using LogNormalize method
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

# Calculate cell cycle scores
seurat_obj <- CellCycleScoring(object = seurat_obj, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)

# Find highly variable features
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# Scale data by regressing out unwanted sources of variation
seurat_obj <- ScaleData(seurat_obj, vars.to.regress = c("nFeature_RNA", "percent.mt"), verbose = FALSE)

# Perform dimensionality reduction and clustering
seurat_obj <- RunPCA(seurat_obj, npcs = 30, verbose = FALSE, seed.use = seed)
seurat_obj <- RunTSNE(seurat_obj, reduction = "pca", dims = 1:30, seed.use = seed)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10, verbose = FALSE, seed.use = seed)
seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = resolution)

# Find doublets using DoubletFinder
suppressMessages(require(DoubletFinder))
nExp <- round(ncol(seurat_obj) * 0.04)  # expect 4% doublets
seurat_obj <- doubletFinder_v3(seurat_obj, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)

# Load reference databases from celldex
hpca.ref <- celldex::HumanPrimaryCellAtlasData()
dice.ref <- celldex::DatabaseImmuneCellExpressionData()
blueprint.ref <- celldex::BlueprintEncodeData()
monaco.ref <- celldex::MonacoImmuneData()
northern.ref <- celldex::NovershternHematopoieticData()

# Convert Seurat object to SingleCellExperiment
sce <- as.SingleCellExperiment(DietSeurat(seurat_obj))

# Auto-annotate cell types using reference databases from celldex
hpca.main <- SingleR(test = sce, assay.type.test = 1, ref = hpca.ref, labels = hpca.ref$label.main)
hpca.fine <- SingleR(test = sce, assay.type.test = 1, ref = hpca.ref, labels = hpca.ref$label.fine)
dice.main <- SingleR(test = sce, assay.type.test = 1, ref = dice.ref, labels = dice.ref$label.main)
dice.fine <- SingleR(test = sce, assay.type.test = 1, ref = dice.ref, labels = dice.ref$label.fine)
blue.main <- SingleR(test = sce, assay.type.test = 1, ref = blueprint.ref, labels = blueprint.ref$label.main)
blue.fine <- SingleR(test = sce, assay.type.test = 1, ref = blueprint.ref, labels = blueprint.ref$label.fine)
monaco.main <- SingleR(test = sce, assay.type.test = 1, ref = monaco.ref, labels = monaco.ref$label.main)
monaco.fine <- SingleR(test = sce, assay.type.test = 1, ref = monaco.ref, labels = monaco.ref$label.fine)
northern.main <- SingleR(test = sce, assay.type.test = 1, ref = northern.ref, labels = northern.ref$label.main)
northern.fine <- SingleR(test = sce, assay.type.test = 1, ref = northern.ref, labels = northern.ref$label.fine)

# Add the celldex annotations to the metadata of the Seurat object
seurat_obj@meta.data$hpca.main <- hpca.main$pruned.labels
seurat_obj@meta.data$hpca.fine <- hpca.fine$pruned.labels
seurat_obj@meta.data$dice.main <- dice.main$pruned.labels
seurat_obj@meta.data$dice.fine <- dice.fine$pruned.labels
seurat_obj@meta.data$monaco.main <- monaco.main$pruned.labels
seurat_obj@meta.data$monaco.fine <- monaco.fine$pruned.labels
seurat_obj@meta.data$northern.main <- northern.main$pruned.labels
seurat_obj@meta.data$northern.fine <- northern.fine$pruned.labels
seurat_obj@meta.data$blue.main <- blue.main$pruned.labels
seurat_obj@meta.data$blue.fine <- blue.fine$pruned.labels

# Add UMAP coordinates to the metadata
UMAP <- as.data.frame(Embeddings(object = seurat_obj[["umap"]]))
seurat_obj@meta.data$UMAP_1 <- UMAP$UMAP_1
seurat_obj@meta.data$UMAP_2 <- UMAP$UMAP_2

# Add tSNE coordinates to the metadata
tsne <- as.data.frame(Embeddings(object = seurat_obj[["tsne"]]))
seurat_obj@meta.data$tSNE_1 <- tsne$tSNE_1
seurat_obj@meta.data$tSNE_2 <- tsne$tSNE_2


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

# Write metadata to a file
write_delim(meta_out, file = file.path(output_folder, paste0(Project_Name,"_metafile_with_annotation.txt")),
            delim = "\t")

# Write seurat H5 file 
SaveH5Seurat(seurat_obj, filename = file.path(output_folder, paste0(Project_Name,"_h5friendly")),
             overwrite = TRUE, verbose = TRUE)

# write out count data
for (i in names(seurat_obj@assays)) {
  
  # Retrieve scaled_counts, raw counts, normalized data, scaled counts, and metadata
  scaled_counts <- as.data.frame(seurat_obj@assays[[i]]@scale.data)
  raw_counts <- as.data.frame(seurat_obj@assays[[i]]@counts)
  normalized_data <- as.data.frame(seurat_obj@assays[[i]]@data)
  
  # Add row names as a column for Alyssa Umap app
  scaled_counts <- tibble::rownames_to_column(scaled_counts, var = "Genes")
  raw_counts <- tibble::rownames_to_column(raw_counts, var = "Genes")
  normalized_data <- tibble::rownames_to_column(normalized_data, var = "Genes")
  
  # Write scaled_counts,raw_counts, normalized_data, and metadata to separate files
  write.table(scaled_counts, file = paste0(output_folder, "/", Project_Name,"_",i,"_ScaledCounts.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  write.table(raw_counts, file = paste0(output_folder, "/", Project_Name,"_",i,"_raw_counts.txt"),
              sep = "\t", row.names = FALSE, col.names = TRUE)
  write.table(normalized_data, file = paste0(output_folder, "/", Project_Name,"_",i,"_normalized_counts.txt"),
              sep = "\t", row.names = FALSE, col.names = TRUE)
}




####---- User Input ----####

## Set local Github repository as working directory
#setwd("C:/Users/Administrator/Desktop/BERLIN_02")

## Desired project name
Project_Name <- "GSE116256_AMLscRNA"

## Seurat H5 file
#h5_file <- "2-Single_Cell_RNAseq_Pipeline/Output/Single_Cell_RNAseq_Output/GSE116256_AMLscRNA_h5friendly.h5seurat"
h5_file <- "Output/Single_Cell_RNAseq_Output/GSE116256_AMLscRNA_h5friendly.h5seurat"
     
## Specify the output folder path
#output_folder <- "2-Single_Cell_RNAseq_Pipeline/Output/Single_Cell_PostProcessing_Output/"
output_folder <- "Output/Single_Cell_PostProcessing_Output/"

## Optional column to use as marker for differential expression analysis
Marker_Column <- "CellType"

## number for random cells to be subset for UMAP App
num_cells <- 1000

## Set seed
seed <- 28



####---- Load Packages ----####
# Define the list of packages to install
packages <- c("dplyr", "Seurat", "patchwork", "readr", "SingleR", "DoubletFinder",
              "Matrix", "fields", "SeuratDisk", "SeuratData", "openxlsx")
# Load in the necessary packages 
invisible(lapply(packages, library, character.only = TRUE))



####---- Run Script ----####

set.seed(seed)

# Load the file as a Seurat object using the Load H5Seurat function
seurat_obj <- SeuratDisk::LoadH5Seurat(h5_file,
                                       assays = c("integrated", "RNA"),
                                       reductions = c("pca", "tsne", "umap", "rpca"),
                                       graphs = NULL,
                                       images = NULL,
                                       meta.data = T,
                                       commands = F,
                                       verbose = F)


# Create an output directory if it doesn't already exist
if (!file.exists(output_folder)) {
  dir.create(output_folder)
}

####---- Find markers and Gene Expression ----#### 

# Generate the Markers for each clusters 
# Set the Idents 
Idents(seurat_obj) <- "seurat_clusters"

# get the number of cells in each cluster
cluster_sizes <- table(Idents(seurat_obj))

# only include clusters with more than 2 cells (FindMaker will only work on clusters that have at least 3 cells )
clusters <- names(cluster_sizes)[cluster_sizes > 2]

wb <- createWorkbook()

SinglCell_Cluster_gs <- data.frame(term = "term",gene = "gene")

SeuratOutput_Folder <- paste0(output_folder,"/Seurat_Cluster_DEG/")
dir.create(SeuratOutput_Folder,recursive = T)

# loop through each cluster and find marker genes
for (cluster_id in clusters) {
  
  print(paste("Finding marker genes for cluster", cluster_id))
  
  # use the FindMarker() function to get the top 10 marker genes
  unfiltered_markers <- FindMarkers(object = seurat_obj, ident.1 = cluster_id, min.pct = .25 )
  upreg_markers <- unfiltered_markers[unfiltered_markers$p_val_adj < 0.05 & unfiltered_markers$avg_log2FC > 0, ]
  dnreg_markers <- unfiltered_markers[unfiltered_markers$p_val_adj < 0.05 & unfiltered_markers$avg_log2FC < 0, ]
  
  unfiltered_markers <- tibble::rownames_to_column(unfiltered_markers, var = "Genes")
  upreg_markers <- tibble::rownames_to_column(upreg_markers, var = "Genes")
  dnreg_markers <- tibble::rownames_to_column(dnreg_markers, var = "Genes")
  
  if (nrow(upreg_markers) > 0) {
    SinglCell_Cluster_gs <- rbind(SinglCell_Cluster_gs,
                                  data.frame(term = paste0(Project_Name,"_cluster_", cluster_id, "_UP"),
                                             gene = upreg_markers$Genes))
  }
  
  
  # write the marker genes to a separate file 
  unfilter_output_file <- file.path(SeuratOutput_Folder, paste0(Project_Name,"_cluster_", cluster_id, "_unfilt_gene_expression.txt"))
  upreg_output_file <- file.path(SeuratOutput_Folder, paste0(Project_Name,"_cluster_", cluster_id, "_upreg_gene_expression.txt"))
  dnreg_output_file <- file.path(SeuratOutput_Folder, paste0(Project_Name,"_cluster_", cluster_id, "_dnreg_gene_expression.txt"))
  
  write.table(unfiltered_markers, file = unfilter_output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  write.table(upreg_markers, file = upreg_output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  write.table(dnreg_markers, file = dnreg_output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  # add a new sheet to the workbook for unfiltered markers and write the marker genes to the sheet
  unfiltered_sheet_name <- paste0("Cluster_", cluster_id, "_unfilt")
  addWorksheet(wb, unfiltered_sheet_name)
  writeData(wb, unfiltered_sheet_name, unfiltered_markers, row.names = FALSE)  # include row.names in the output
  setColWidths(wb, unfiltered_sheet_name, cols = 1:ncol(unfiltered_markers), widths = "auto") # Autofit columns
  
  # add a new sheet to the workbook for upregulated markers and write the marker genes to the sheet
  upregulated_sheet_name <- paste0("Cluster_", cluster_id, "_upreg")
  addWorksheet(wb, upregulated_sheet_name)
  writeData(wb, upregulated_sheet_name, upreg_markers, row.names = FALSE)  # include row.names in the output
  setColWidths(wb, upregulated_sheet_name, cols = 1:ncol(upreg_markers), widths = "auto") # Autofit columns
  
  # add a new sheet to the workbook for upregulated markers and write the marker genes to the sheet
  dnregulated_sheet_name <- paste0("Cluster_", cluster_id, "_dnreg")
  addWorksheet(wb, dnregulated_sheet_name)
  writeData(wb, dnregulated_sheet_name, dnreg_markers, row.names = FALSE)  # include row.names in the output
  setColWidths(wb, dnregulated_sheet_name, cols = 1:ncol(dnreg_markers), widths = "auto") # Autofit columns
  
}

SinglCell_Cluster_gs <- SinglCell_Cluster_gs[-1,]

# Save the workbook to an Excel file
xlsx_output_file <- file.path(SeuratOutput_Folder, paste0(Project_Name,"_SeuratCluster_DEG.xlsx"))
saveWorkbook(wb, xlsx_output_file)


if (all(c(!is.na(Marker_Column),!is.null(Marker_Column)),Marker_Column!="")) {
  MarkerOutput_Folder <- paste0(output_folder,"/",Marker_Column,"_Cluster_DEG/")
  dir.create(MarkerOutput_Folder,recursive = T)
  # Generate marker list  by cell types
  wb2 <-  createWorkbook()
  Idents(seurat_obj) <- Marker_Column
  cell_types <- unique(seurat_obj@meta.data[,Marker_Column])
  
  for(cells in cell_types){
    
    if (length(seurat_obj@meta.data[which(seurat_obj@meta.data[,Marker_Column] == cells),Marker_Column]) >= 3) {
      
      print(paste("Finding marker genes for", cells))
      
      # use the FindMarker() function to get the top 10 marker genes
      unfiltered_markers <- FindMarkers(object = seurat_obj, ident.1 = cells, min.pct = .25 )
      upreg_markers <- unfiltered_markers[unfiltered_markers$p_val_adj < 0.05 & unfiltered_markers$avg_log2FC > 0, ]
      dnreg_markers <- unfiltered_markers[unfiltered_markers$p_val_adj < 0.05 & unfiltered_markers$avg_log2FC < 0, ]
      
      unfiltered_markers <- tibble::rownames_to_column(unfiltered_markers, var = "Genes")
      upreg_markers <- tibble::rownames_to_column(upreg_markers, var = "Genes")
      dnreg_markers <- tibble::rownames_to_column(dnreg_markers, var = "Genes")
      
      if (nrow(upreg_markers) > 0) {
        SinglCell_Cluster_gs <- rbind(SinglCell_Cluster_gs,
                                      data.frame(term = paste0(Project_Name,"_",cells, "_UP"),
                                                 gene = upreg_markers$Genes))
      }
      
      # write the marker genes to a separate file 
      unfilter_output_file <- file.path(MarkerOutput_Folder, gsub("[[:punct:]]",".",gsub(" ","_",paste0(Project_Name,"_",cells, "_unfilt_gene_expression.txt"))))
      upreg_output_file <- file.path(MarkerOutput_Folder, gsub("[[:punct:]]",".",gsub(" ","_",paste0(Project_Name,"_",cells, "_upreg_gene_expression.txt"))))
      dnreg_output_file <- file.path(MarkerOutput_Folder, gsub("[[:punct:]]",".",gsub(" ","_",paste0(Project_Name,"_",cells, "_dnreg_gene_expression.txt"))))
      
      write.table(unfiltered_markers, file = unfilter_output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
      write.table(upreg_markers, file = upreg_output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
      write.table(dnreg_markers, file = dnreg_output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
      
      # add a new sheet to the workbook for unfiltered markers and write the marker genes to the sheet
      unfiltered_sheet_name <- paste0(cells, "_unfilt")
      addWorksheet(wb2, unfiltered_sheet_name)
      writeData(wb2, unfiltered_sheet_name, unfiltered_markers, row.names = FALSE)  # include row.names in the output
      setColWidths(wb2, unfiltered_sheet_name, cols = 1:ncol(unfiltered_markers), widths = "auto") # Autofit columns
      
      # add a new sheet to the workbook for upregulated markers and write the marker genes to the sheet
      upregulated_sheet_name <- paste0(cells, "_upreg")
      addWorksheet(wb2, upregulated_sheet_name)
      writeData(wb2, upregulated_sheet_name, upreg_markers, row.names = FALSE)  # include row.names in the output
      setColWidths(wb2, upregulated_sheet_name, cols = 1:ncol(upreg_markers), widths = "auto") # Autofit columns
      
      # add a new sheet to the workbook for upregulated markers and write the marker genes to the sheet
      dnregulated_sheet_name <- paste0(cells, "_dnreg")
      addWorksheet(wb2, dnregulated_sheet_name)
      writeData(wb2, dnregulated_sheet_name, dnreg_markers, row.names = FALSE)  # include row.names in the output
      setColWidths(wb2, dnregulated_sheet_name, cols = 1:ncol(dnreg_markers), widths = "auto") # Autofit columns
      
    }
    
  }
  
  
  # Save the workbook to an Excel file
  xlsx_output_cell_file <- file.path(MarkerOutput_Folder, paste0(Project_Name,"_",Marker_Column,"_DEG.xlsx"))
  saveWorkbook(wb2, xlsx_output_cell_file)
  
}




####---- Random Cell matrix for Umap app ----####

ShinyAppOutput_Folder <- paste0(output_folder,"/RShiny_App_Input_Files/")
dir.create(ShinyAppOutput_Folder,recursive = T)

# Randomly sample cells
random_cells <- sample(1:ncol(seurat_obj), size = num_cells, replace = FALSE)
subset_seurat_obj <- seurat_obj[,random_cells]

# Iterate over columns of the subset Seurat object
for (i in names(subset_seurat_obj@assays)) {
  
  # Retrieve scaled_counts, raw counts, normalized data, scaled counts, and metadata
  scaled_counts <- as.data.frame(subset_seurat_obj@assays[[i]]@scale.data)
  raw_counts <- as.data.frame(subset_seurat_obj@assays[[i]]@counts)
  normalized_data <- as.data.frame(subset_seurat_obj@assays[[i]]@data)
  metadata <- as.data.frame(subset_seurat_obj@meta.data)
  
  # Add row names as a column for Alyssa Umap app
  scaled_counts <- tibble::rownames_to_column(scaled_counts, var = "Genes")
  raw_counts <- tibble::rownames_to_column(raw_counts, var = "Genes")
  normalized_data <- tibble::rownames_to_column(normalized_data, var = "Genes")
  metadata <- tibble::rownames_to_column(metadata, var = "Cell")
  
  # Write scaled_counts,raw_counts, normalized_data, and metadata to separate files
  write.table(scaled_counts, file = paste0(ShinyAppOutput_Folder, "/", Project_Name,"_",num_cells,"_",i,"_ScaledCounts.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  write.table(raw_counts, file = paste0(ShinyAppOutput_Folder, "/", Project_Name,"_",num_cells,"_",i,"_raw_counts.txt"),
              sep = "\t", row.names = FALSE, col.names = TRUE)
  write.table(normalized_data, file = paste0(ShinyAppOutput_Folder, "/", Project_Name,"_",num_cells,"_",i,"_normalized_counts.txt"),
              sep = "\t", row.names = FALSE, col.names = TRUE)
  write.table(metadata, file = paste0(ShinyAppOutput_Folder, "/", Project_Name,"_",num_cells,"_",i,"_metafile.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

####----Write out gene set----####

SinglCell_Cluster_gs_list <- list()
for (i in unique(SinglCell_Cluster_gs[,1])){
  SinglCell_Cluster_gs_list[[i]] <- SinglCell_Cluster_gs[SinglCell_Cluster_gs[,1] == i,][,2]
}

save(SinglCell_Cluster_gs_list,file = paste0(ShinyAppOutput_Folder,"/",Project_Name,"_SingleCell_Cluster_DEG_GeneSet.RData"))
write_delim(SinglCell_Cluster_gs,paste0(ShinyAppOutput_Folder,"/",Project_Name,"_SingleCell_Cluster_DEG_GeneSet.txt"), delim = '\t')



####---- User Input ----####

## Set local Github repository as working directory
#setwd("BERLIN")

## Load Functions from Customized ISCVAM H5 File
#source("2-Single_Cell_RNAseq_Pipeline/R_Scripts/R/write.h5.R")
source("R/write.h5.R")

## Desired project name
Project_Name <- "GSE116256_AMLscRNA"

## Seurat H5 file
#h5_file <- "2-Single_Cell_RNAseq_Pipeline/Output/Single_Cell_RNAseq_Output/GSE116256_AMLscRNA_h5friendly.h5seurat"
h5_file <- "Output/Single_Cell_RNAseq_Output/GSE116256_AMLscRNA_h5friendly.h5seurat"
            
## Specify the output folder path
#output_folder <- "2-Single_Cell_RNAseq_Pipeline/Output/Generate_ISCVA_H5_Output"

output_folder <- "Output/Generate_ISCVA_H5_Output/"

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}



####---- Load Packages ----####
# Define the list of packages to install
packages <- c("dplyr", "Seurat", "rhdf5", "SeuratDisk")
# Load in the necessary packages 
invisible(lapply(packages, library, character.only = TRUE))



####---- Run Script ----####

# Load H5 Seurat Compliant File
seurat_obj <- LoadH5Seurat(h5_file,
                           assays = c("integrated", "RNA"),
                           reductions = c("pca", "tsne", "umap", "rpca"),
                           graphs = NULL,
                           images = NULL,
                           meta.data = TRUE,
                           commands = FALSE,
                           verbose = FALSE)

# Create Seurat.list
seurat.list <- list(all = list(seurat = seurat_obj, covs = seurat_obj@meta.data))

# Convert Seurat.list to H5 File
output_h5_file_name <- file.path(output_folder, paste0(Project_Name,"_ISCVA.h5"))
write.h5(seurat.list = seurat.list, fn = output_h5_file_name)

# Check Dimensions of the H5 File
h5 <- H5Fopen(output_h5_file_name)
h5ls(h5)

# Check Covs in ISCVA Format
colnames(h5$'/artifacts/all/covs')

# Check DiscreteCovs and ContinuousCovs
h5$'/artifacts/all/discreteCovs'
h5$'/artifacts/all/continuousCovs'

# Close the H5 File
h5closeAll()