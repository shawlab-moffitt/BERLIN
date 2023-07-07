



####---- User Input ----####

## Set local Github repository as working directory
setwd("BERLIN")

## Desired project name
Project_Name <- "GSE116256_AMLscRNA"

# Define the path to the count file
Count_file <- "2-Single_Cell_RNAseq_Pipeline/Input/Count_File_Example_GSE116256.txt"

# Define path to clinical/meta data (OPTIONAL)
meta_file <- "2-Single_Cell_RNAseq_Pipeline/Input/Meta_File_Example_GSE116256.txt"

# Specify the output folder path
output_folder <- "2-Single_Cell_RNAseq_Pipeline/Output/Single_Cell_RNAseq_Output/"






####---- Load Packages ----####

packages <- c("dplyr","Seurat","patchwork","readr","DoubletFinder","celldex","Matrix",
              "fields","SeuratDisk","SeuratData","SingleR","SingleCellExperiment")
invisible(lapply(packages, library, character.only = TRUE))






####---- Run Script ----####

# Create an output directory if it doesn't already exist
if (!file.exists(output_folder)) {
  dir.create(output_folder)
}

if (file.exists(meta_file)) {
  meta <- read.csv2(meta_file, sep = "", header = TRUE)
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
seurat_obj <- RunPCA(seurat_obj, npcs = 30, verbose = FALSE, seed.use = 42)
seurat_obj <- RunTSNE(seurat_obj, reduction = "pca", dims = 1:30, seed.use = 1)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10, verbose = FALSE, seed.use = 42)
seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = c(0.10, 0.15, 0.25,0.75))

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
seurat_obj@meta.data$tSNE_2 <- tsne$tSNE_1

meta_out <- seurat_obj@meta.data
meta_out$Cell <- rownames(meta_out)
meta_out <- meta_out %>% relocate(Cell)
if (!is.null(meta)) {
  meta_out <- merge(meta,meta_out, by = "Cell", all.y = T)
  seurat_obj@meta.data <- meta_out
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











