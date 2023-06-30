# The H5 script is taking the H5 file and loading it into a seurat object and outputing: 
## Gene expression txt files for each cluster and excel file: containing unfiltered and filtered by p_val_adj < 0.05
## Gene expression txt file for cell types and excel file: containing unfiltered and filtered by p_val_adj < 0.05
## Subsetting 1000 matrix and outputing, raw, normalized, scaled and metadata 



####---- User Input ----####

## Set local Github repository as working directory
setwd("BERLIN")

## Desired project name
Project_Name <- "GSE116256_AMLscRNA"

## Seurat H5 file
h5_file <- "2-Single_Cell_RNAseq_Pipeline/Input/GSE116256_AMLscRNA_h5friendly.h5seurat"

## Specify the output folder path
output_folder <- "2-Single_Cell_RNAseq_Pipeline/Output/"

## Set seed and number for random cells to be subset for UMAP App
set.seed(28)
num_cells <- 1000



####---- Load Packages ----####
# Define the list of packages to install
packages <- c("dplyr", "Seurat", "patchwork", "readr", "SingleR", "DoubletFinder",
              "Matrix", "fields", "SeuratDisk", "SeuratData", "openxlsx")
# Load in the necessary packages 
invisible(lapply(packages, library, character.only = TRUE))



####---- Run Script ----####

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

# Assign the full output directory path to a variable
full_output_dir <- paste0(getwd(),"/", output_folder, sep = "/")

####---- Find markers and Gene Expression ----#### 

# Generate the Markers for each clusters 
# Set the Idents 
Idents(seurat_obj) <- "seurat_clusters"


# get the number of cells in each cluster
cluster_sizes <- table(Idents(seurat_obj))

# only include clusters with more than 2 cells (FindMaker will only work on clusters that have at least 3 cells )
clusters <- names(cluster_sizes)[cluster_sizes > 2]

wb <- createWorkbook()

# loop through each cluster and find marker genes
for (cluster_id in clusters) {
  
  print(paste("Finding marker genes for cluster", cluster_id))
  
  # use the FindMarker() function to get the top 10 marker genes
  unfiltered_markers <- FindMarkers(object = seurat_obj, ident.1 = cluster_id, min.pct = .25 )
  upreg_markers <- unfiltered_markers[unfiltered_markers$p_val_adj < 0.05, ]
  
  unfiltered_markers <- tibble::rownames_to_column(unfiltered_markers, var = "Genes")
  upreg_markers <- tibble::rownames_to_column(upreg_markers, var = "Genes")
  
  
  ########### separate up and down reg
  
  # write the marker genes to a separate file 
  unfilter_output_file<- file.path(full_output_dir, paste0(Project_Name,"_cluster_", cluster_id, "_unfilt_gene_expression.txt"))
  upreg_output_file <- file.path(full_output_dir, paste0(Project_Name,"_cluster_", cluster_id, "_upreg_gene_expression.txt"))
  
  
  write.table(unfiltered_markers, file = unfilter_output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  write.table(upreg_markers, file = upreg_output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
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
  
}

# Save the workbook to an Excel file
xlsx_output_file <- file.path(full_output_dir, paste0(Project_Name,"_marker_genes_new.xlsx"))
saveWorkbook(wb, xlsx_output_file)




# Generate marker list  by cell types
wb2 <-  createWorkbook()


Idents(seurat_obj) <- "seurat_clusters_gabby_annotation"

cell_types <- unique(seurat_obj@meta.data$seurat_clusters_gabby_annotation)


for(cells in cell_types){
  print(paste("Finding marker genes for", cells))
  
  # use the FindMarker() function to get the top 10 marker genes
  unfiltered_markers <- FindMarkers(object = seurat_obj, ident.1 = cells, min.pct = .25 )
  upreg_markers <- unfiltered_markers[unfiltered_markers$p_val_adj < 0.05, ]
  
  unfiltered_markers <- tibble::rownames_to_column(unfiltered_markers, var = "Genes")
  upreg_markers <- tibble::rownames_to_column(upreg_markers, var = "Genes")
  
  # write the marker genes to a separate file 
  unfilter_output_file<- file.path(full_output_dir, paste0(Project_Name,"_",cells, "_unfilt_gene_expression.txt"))
  upreg_output_file <- file.path(full_output_dir, paste0(Project_Name,"_",cells, "_upreg_gene_expression.txt"))
  
  
  write.table(unfiltered_markers, file = unfilter_output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  write.table(upreg_markers, file = upreg_output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
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
  
  
}

# Save the workbook to an Excel file
xlsx_output_cell_file <- file.path(full_output_dir, paste0(Project_Name,"_cell_marker_genes.xlsx"))
saveWorkbook(wb2, xlsx_output_cell_file)



####---- Random 1000 Cell matrix for Umap app ----####

folder_path <- paste0(full_output_dir, "/", num_cells, "_matrix")
dir.create(folder_path, showWarnings = FALSE)


# Randomly sample cells
random_cells <- sample(1:nrow(seurat_obj), size = num_cells, replace = FALSE)
subset_seurat_obj <- seurat_obj[random_cells, ]

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
  metadata <- tibble::rownames_to_column(metadata, var = "Barcodes")
  
  # Write scaled_counts,raw_counts, normalized_data, and metadata to separate files
  write.table(scaled_counts, file = paste0(folder_path, "/", Project_Name,"_",num_cells,"_",i,"_ScaledCounts.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  write.table(raw_counts, file = paste0(folder_path, "/", Project_Name,"_",num_cells,"_",i,"_raw_counts.txt"),
              sep = "\t", row.names = FALSE, col.names = TRUE)
  write.table(normalized_data, file = paste0(folder_path, "/", Project_Name,"_",num_cells,"_",i,"_normalized_counts.txt"),
              sep = "\t", row.names = FALSE, col.names = TRUE)
  write.table(metadata, file = paste0(folder_path, "/", Project_Name,"_",num_cells,"_",i,"_metafile.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}


#####---- Find All Markers ----####
# Set the identity of the cells in the Seurat object
Idents(seurat_obj) <- "seurat_clusters_gabby_annotation"

# Find markers for each cluster with a minimum percentage threshold
all.markers <- FindAllMarkers(seurat_obj, min.pct = 0.1)

# Calculate the difference between pct.1 and pct.2 and add a new column
all.markers$diff_pct <- all.markers$pct.1 - all.markers$pct.2

# Specify the file path for saving the all markers data
all_marker_file_path <- file.path(full_output_dir, paste0(Project_Name,"_all_markers.txt"))

all.markers <- all.markers %>% relocate("gene")

# Write the all markers data to a tab-separated file
write.table(all.markers, file = all_marker_file_path, 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

