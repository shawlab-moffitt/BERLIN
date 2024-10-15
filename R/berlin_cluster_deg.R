



#' Perform differential gene expression on clusters in a Seurat object
#'
#' @param object Seurat object. If NULL, must supply count data.
#' @param assay String. Name of the initial/default assay.
#' @param cluster_col String. Column name from meta data of primary identified clusters. Default to 'seurat_clusters'.
#' @param meta_cluster_col String. Column name from meta data of an additional column to use for clustering and differential gene expression. Optional.
#' @param min.pct Numeric. Only test genes that are detected in a minimum fraction of cells. Default is 0.25.
#' @param adjpval_cutoff Numeric. An adjusted p-value value cutoff for identifying significantly expressed genes. default is 0.05.
#' @param log2fc_cutoff Numeric. An absolute value of a log2 fold change cutoff to identify up and down regulated genes. Defualt is 0.
#' @param verbose Boolean. To show progress, TRUE, else FALSE.
#'
#' @return A list of data frames containing DGE of clusters and a geneset of the markers of each cluster.
#' @export
#'
#'



berlin_cluster_deg <- function(object = NULL, assay = "RNA", cluster_col = "seurat_clusters", meta_cluster_col = NULL, min.pct = 0.25,
                               adjpval_cutoff = 0.05, log2fc_cutoff = 0, verbose = TRUE) {

  if (is.null(object)) stop("Please supply Seurat object as input.")

  if (!assay %in% names(object)) stop("Assay input is not found in object")
  DefaultAssay(object) <- assay

  SinglCell_Cluster_gs <- data.frame(term = character(0),gene = character(0))
  cluster_marker_list <- list()

  if (!is.null(cluster_col)) {
    if (!cluster_col %in% names(object[[]])) {
      stop("cluster_col input not found in object meta data.")
    } else {
      Idents(object) <- cluster_col
      # get the number of cells in each cluster
      cluster_sizes <- table(Idents(object))
      # only include clusters with more than 2 cells (FindMaker will only work on clusters that have at least 3 cells )
      clusters <- names(cluster_sizes)[cluster_sizes > 2]

      for (id in clusters) {
        unfiltered_markers <- Seurat::FindMarkers(object = object, ident.1 = id, min.pct = min.pct )
        upreg_markers <- unfiltered_markers[unfiltered_markers$p_val_adj < adjpval_cutoff & unfiltered_markers$avg_log2FC > log2fc_cutoff,]
        dnreg_markers <- unfiltered_markers[unfiltered_markers$p_val_adj < adjpval_cutoff & unfiltered_markers$avg_log2FC < log2fc_cutoff,]
        unfiltered_markers <- tibble::rownames_to_column(unfiltered_markers, var = "Genes")
        upreg_markers <- tibble::rownames_to_column(upreg_markers, var = "Genes")
        dnreg_markers <- tibble::rownames_to_column(dnreg_markers, var = "Genes")
        cluster_marker_list[[paste0(cluster_col,"_Cluster_",id,"_Unfiltered_DEG")]] <- unfiltered_markers
        cluster_marker_list[[paste0(cluster_col,"_Cluster_",id,"_Upregulated_DEG")]] <- upreg_markers
        cluster_marker_list[[paste0(cluster_col,"_Cluster_",id,"_Downregulated_DEG")]] <- dnreg_markers
        if (nrow(upreg_markers) > 0) {
          SinglCell_Cluster_gs <- rbind(SinglCell_Cluster_gs,
                                        data.frame(term = rep(paste0(cluster_col,"_Cluster_", id, "_Upregulated"),length(upreg_markers$Genes)),
                                                   gene = upreg_markers$Genes))
        }
      }
    }

  }

  if (!is.null(meta_cluster_col)) {
    if (!meta_cluster_col %in% names(object[[]])) {
      stop("meta_cluster_col input not found in object meta data.")
    } else {
      Idents(object) <- meta_cluster_col
      # get the number of cells in each cluster
      meta_cluster_col_sizes <- table(Idents(object))
      # only include clusters with more than 2 cells (FindMaker will only work on clusters that have at least 3 cells )
      meta_cluster_cols <- names(meta_cluster_col_sizes)[meta_cluster_col_sizes > 2]

      for (id in meta_cluster_cols) {
        unfiltered_markers <- Seurat::FindMarkers(object = object, ident.1 = id, min.pct = min.pct )
        upreg_markers <- unfiltered_markers[unfiltered_markers$p_val_adj < adjpval_cutoff & unfiltered_markers$avg_log2FC > log2fc_cutoff,]
        dnreg_markers <- unfiltered_markers[unfiltered_markers$p_val_adj < adjpval_cutoff & unfiltered_markers$avg_log2FC < log2fc_cutoff,]
        cluster_marker_list[[paste0(meta_cluster_col,"_",id,"_Unfiltered_DEG")]] <- unfiltered_markers
        cluster_marker_list[[paste0(meta_cluster_col,"_",id,"_Upregulated_DEG")]] <- upreg_markers
        cluster_marker_list[[paste0(meta_cluster_col,"_",id,"_Downregulated_DEG")]] <- dnreg_markers
        unfiltered_markers <- tibble::rownames_to_column(unfiltered_markers, var = "Genes")
        upreg_markers <- tibble::rownames_to_column(upreg_markers, var = "Genes")
        dnreg_markers <- tibble::rownames_to_column(dnreg_markers, var = "Genes")
        if (nrow(upreg_markers) > 0) {
          SinglCell_Cluster_gs <- rbind(SinglCell_Cluster_gs,
                                        data.frame(term = paste0(meta_cluster_col,"_", id, "_Upregulated"),
                                                   gene = upreg_markers$Genes))
        }
      }
    }

  }

  output_list <- list(Cluster_Markers_DEG = cluster_marker_list,
                      Cluster_Markers_Geneset = SinglCell_Cluster_gs)
  return(output_list)


}
