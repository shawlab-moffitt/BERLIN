



#' Perform dimension reduction and clustering analysis on Seurat object
#'
#' @param object Seurat object. If NULL, must supply count data.
#' @param assay String. Name of the initial/default assay.
#' @param seed Numeric. Random seed, default is 42. Setting to NULL will remove seed.
#' @param verbose Boolean. To show progress, TRUE, else FALSE.
#' @param resolution Numeric vector. Above one to obtain a larger number of communities, and below one to obtain a smaller number.
#' @param pca_npcs Numeric. Total number of principal components to compute for PCA plot. Default is 30.
#' @param tsne_dims Numeric sequence. Dimensions to use as input features for TSNE plot. Default 1:30.
#' @param umap_dims  Numeric sequence. Dimensions to use as input features for UMAP. Default 1:10.
#' @param neighbor_dims Numeric sequence. Dimension of reduction to use as input to find neighbors.
#' @param doublet_pN Numeric. The number of generated artificial doublets, expressed as a proportion of the merged real-artificial data. Default is set to 0.25.
#' @param doublet_pK Numeric. The principal component neighborhood size to compute pANN, expressed as a proportion of the merged real-artificial data. Default is 0.09.
#' @param doublet_prop Numeric. Estimated proportion of homotypic doublets. Default is 0.04.
#' @param doublet_PCs Numeric. Number of statistically-significant principal components.
#'
#' @return A Seurat object.
#' @export
#'
#'


berlin_cluster <- function(object = NULL, assay = "RNA", seed = 42, verbose = TRUE, resolution = c(0.5,1,1.5,2),
                                pca_npcs = 30, tsne_dims = 1:30, umap_dims = 1:10, neighbor_dims = 1:30) {


  if (is.null(object)) stop("Please supply Seurat object")

  if (!assay %in% names(object)) stop("Assay input is not found in object")
  DefaultAssay(object) <- assay

  object <- Seurat::RunPCA(object, npcs = pca_npcs, verbose = verbose, seed.use = seed)
  object <- Seurat::RunTSNE(object, reduction = "pca", dims = tsne_dims, seed.use = seed)
  object <- Seurat::RunUMAP(object, dims = umap_dims, verbose = verbose, seed.use = seed)
  object <- Seurat::FindNeighbors(object, reduction = "pca", dims = neighbor_dims, verbose = verbose)
  object <- Seurat::FindClusters(object, resolution = resolution, verbose = verbose)

  # Add UMAP coordinates to the metadata
  UMAP <- as.data.frame(Embeddings(object = object[["umap"]]))
  object <- Seurat::AddMetaData(object,metadata = UMAP)
  # Add tSNE coordinates to the metadata
  tsne <- as.data.frame(Embeddings(object = object[["tsne"]]))
  object <- Seurat::AddMetaData(object,metadata = tsne)

  return(object)

}
