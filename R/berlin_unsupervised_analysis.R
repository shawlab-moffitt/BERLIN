



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
#'
#' @return A Seurat object.
#' @export
#'
#'


berlin_unsupervised_analysis <- function(object = NULL, assay = "RNA", seed = 42, verbose = TRUE, resolution = c(0.5,1,1.5,2),
                                         run_pca = TRUE, run_umap = TRUE, run_TSNE = FALSE,
                                         pca_npcs = 30, tsne_dims = 1:30, umap_dims = 1:10, neighbor_dims = 1:30) {

  call.string <- deparse(expr = sys.calls()[[1]])
  func_name <- "berlin_unsupervised_analysis"
  time.stamp <- Sys.time()
  argg <- c(as.list(environment()))
  argg <- Filter(function(x) any(is.numeric(x) | is.character(x)), argg)

  if (is.null(object)) stop("Please supply Seurat object")

  if (!assay %in% names(object)) stop("Assay input is not found in object")
  SeuratObject::DefaultAssay(object) <- assay

  if (run_pca) {
    object <- Seurat::RunPCA(object, npcs = pca_npcs, verbose = verbose, seed.use = seed)
  }
  if (run_umap) {
    object <- Seurat::RunUMAP(object, dims = umap_dims, verbose = verbose, seed.use = seed)
    # Add UMAP coordinates to the metadata
    UMAP <- as.data.frame(SeuratObject::Embeddings(object = object[["umap"]]))
    object <- Seurat::AddMetaData(object,metadata = UMAP)
  }
  if (run_TSNE) {
    object <- Seurat::RunTSNE(object, reduction = "pca", dims = tsne_dims, seed.use = seed)
    # Add tSNE coordinates to the metadata
    tsne <- as.data.frame(SeuratObject::Embeddings(object = object[["tsne"]]))
    object <- Seurat::AddMetaData(object,metadata = tsne)
  }
  object <- Seurat::FindNeighbors(object, reduction = "pca", dims = neighbor_dims, verbose = verbose)
  object <- Seurat::FindClusters(object, resolution = resolution, verbose = verbose)


  slot(object = object, name = "commands")[[func_name]] <- list(name = func_name,
                                                                time.stamp = time.stamp,
                                                                call.string = call.string,
                                                                params = argg)

  return(object)

}
