



#' Perform differential gene expression on clusters in a Seurat object
#'
#' @param object Seurat object. If NULL, must supply count data.
#' @param assay String. Name of the initial/default assay. Default is RNA.
#' @param cluster_cols String or Vector. Column name from meta data of primary identified clusters. Default to 'seurat_clusters'. Multiple columns allowed.
#' @param min.pct Numeric. Only test genes that are detected in a minimum fraction of cells. Default is 0.25.
#' @param adjpval_cutoff Numeric. An adjusted p-value value cutoff for identifying significantly expressed genes. Default is 0.05.
#' @param log2fc_cutoff Numeric. An absolute value of a log2 fold change cutoff to identify up and down regulated genes. Default is 0.
#' @param verbose Boolean. To show progress, TRUE, else FALSE.
#'
#' @return A list of data frames containing DGE of clusters and a geneset of the markers of each cluster.
#' @export
#'
#'



berlin_cluster_deg <- function(object = NULL, assay = "RNA", cluster_cols = "seurat_clusters", min.pct = 0.25,
                               adjpval_cutoff = 0.05, log2fc_cutoff = 0, verbose = TRUE) {

  if (is.null(object)) stop("Please supply Seurat object as input.")

  if (!assay %in% names(object)) stop("Assay input is not found in object")
  SeuratObject::DefaultAssay(object) <- assay

  SinglCell_Cluster_gs <- data.frame(term = character(0),gene = character(0))
  cluster_marker_list <- list()

  for (cluster_col in cluster_cols) {
    if (verbose) {
      base::message(paste0("Finding markers on cluster column: ",cluster_col))
    }
    if (!is.null(cluster_col)) {
      if (!cluster_col %in% names(object[[]])) {
        stop(paste0(cluster_col," column not found in object meta data."))
      } else {
        SeuratObject::Idents(object) <- cluster_col
        # get the number of cells in each cluster
        cluster_sizes <- table(SeuratObject::Idents(object))
        # only include clusters with more than 2 cells (FindMaker will only work on clusters that have at least 3 cells )
        clusters <- names(cluster_sizes)[cluster_sizes > 2]
        i <- 1
        for (id in clusters) {
          if (verbose) {
            pb = utils::txtProgressBar(min = 0, max = length(clusters), initial = 0, style = 3)
            utils::setTxtProgressBar(pb,i)
          }
          unfiltered_markers <- Seurat::FindMarkers(object = object, ident.1 = id, min.pct = min.pct )
          upreg_markers <- unfiltered_markers[unfiltered_markers$p_val_adj < adjpval_cutoff & unfiltered_markers$avg_log2FC > log2fc_cutoff,]
          dnreg_markers <- unfiltered_markers[unfiltered_markers$p_val_adj < adjpval_cutoff & unfiltered_markers$avg_log2FC < log2fc_cutoff,]
          unfiltered_markers <- tibble::rownames_to_column(unfiltered_markers, var = "Genes")
          upreg_markers <- tibble::rownames_to_column(upreg_markers, var = "Genes")
          dnreg_markers <- tibble::rownames_to_column(dnreg_markers, var = "Genes")
          cluster_marker_list[[paste0(cluster_col,"_Cluster_",id,"_Unfiltered_Markers")]] <- unfiltered_markers
          #cluster_marker_list[[paste0(cluster_col,"_Cluster_",id,"_Upregulated_DEG")]] <- upreg_markers
          #cluster_marker_list[[paste0(cluster_col,"_Cluster_",id,"_Downregulated_DEG")]] <- dnreg_markers
          if (nrow(upreg_markers) > 0) {
            SinglCell_Cluster_gs <- rbind(SinglCell_Cluster_gs,
                                          data.frame(term = rep(paste0(cluster_col,"_Cluster_", id, "_Upregulated"),length(upreg_markers$Genes)),
                                                     gene = upreg_markers$Genes))
          }
          i <- i+1
        }
        close(pb)
      }

    }
  }

  output_list <- list(Cluster_Markers_Unfiltered = cluster_marker_list,
                      Cluster_Markers_Geneset = SinglCell_Cluster_gs)
  return(output_list)


}



#' Derive table of intersected genes between clusters and geneset
#'
#' @param object Seurat object. If NULL, must supply count data.
#' @param geneset Data.frame or list. Two column data frame with first column consisting of geneset names and second column of the genes withing that geneset. One gene per row, the geneset name can repeat in the first column. If list input, this would be a named list where the names are the geneset names and each element a vector of genes as strings.
#' @param cluster_col String. Column name from meta data of primary identified clusters. Default to 'seurat_clusters'.
#' @param assay String. Name of the initial/default assay. Default is RNA.
#' @param min.pct Numeric. Only test genes that are detected in a minimum fraction of cells. Default is 0.25.
#' @param adjpval_cutoff Numeric. An adjusted p-value value cutoff for identifying significantly expressed genes. Default is 0.05.
#' @param log2fc_cutoff Numeric. An absolute value of a log2 fold change cutoff to identify up and down regulated genes. Default is 0.
#' @param verbose Boolean. To show progress, TRUE, else FALSE.
#'
#' @return Data.frame
#' @export
#'


berlin_intersect <- function(object = NULL, geneset = NULL, cluster_col = "seurat_clusters", assay = "RNA", min.pct = 0.25,
                               adjpval_cutoff = 0.05, log2fc_cutoff = 0, verbose = TRUE) {


  if (is.null(object)) stop("Please supply Seurat object as input.")

  if (!assay %in% names(object)) stop("Assay input is not found in object")
  SeuratObject::DefaultAssay(object) <- assay

  if (!is.null(cluster_col)) {
    if (!cluster_col %in% names(object[[]])) {
      stop("cluster_col input not found in object meta data.")
    } else {
      SeuratObject::Idents(object) <- cluster_col
      # get the number of cells in each cluster
      cluster_sizes <- table(SeuratObject::Idents(object))
      # only include clusters with more than 2 cells (FindMaker will only work on clusters that have at least 3 cells )
      clusters <- names(cluster_sizes)[cluster_sizes > 2]
      if (verbose) {
        base::message(paste0("Finding markers for cluster column: ",cluster_col))
      }
      cluster_marker_list <- lapply(clusters, function(x) {
        markers <- Seurat::FindMarkers(object = object, ident.1 = x, min.pct = min.pct )
        upreg_markers <- rownames(markers[markers$p_val_adj < adjpval_cutoff & markers$avg_log2FC > log2fc_cutoff,])
        upreg_markers
      })
      names(cluster_marker_list) <- paste0(cluster_col,"_",clusters)

      if (is.data.frame(geneset) & ncol(geneset) == 2) {
        gs <- split(geneset[,2],geneset[,1])
      }
      if (verbose) {
        base::message("Intersecting cluster markers with geneset.")
      }
      cluster_marker_list_gsint <- lapply(cluster_marker_list, function(x) lapply(gs, function(y) intersect(x,y)))

      cluster_marker_list_gsint_df <- do.call(rbind,lapply(cluster_marker_list_gsint, function(x) {
        do.call(cbind,lapply(x, function(y) {
          out <- paste(y,collapse = ",")
          out[which(out=='')] <- NA
          out
        }))
      }))
      cluster_marker_list_gsint_df <- as.data.frame(cbind(Cluster = names(cluster_marker_list_gsint),cluster_marker_list_gsint_df))
      return(cluster_marker_list_gsint_df)

    }
  }
}


#' Generate dotplot heatmap of geneset gene expression within clusters
#'
#' @param objectSeurat object. If NULL, must supply count data.
#' @param geneset Data.frame or list. Two column data frame with first column consisting of geneset names and second column of the genes withing that geneset. One gene per row, the geneset name can repeat in the first column. If list input, this would be a named list where the names are the geneset names and each element a vector of genes as strings.
#' @param cluster_col String. Column name from meta data of primary identified clusters. Default to 'seurat_clusters'.
#' @param assay String. Name of the initial/default assay. Default is RNA.
#' @param min.pct Numeric. Only test genes that are detected in a minimum fraction of cells. Default is 0.25.
#' @param adjpval_cutoff Numeric. An adjusted p-value value cutoff for identifying significantly expressed genes. Default is 0.05.
#' @param log2fc_cutoff Numeric. An absolute value of a log2 fold change cutoff to identify up and down regulated genes. Default is 0.
#' @param verbose Boolean. To show progress, TRUE, else FALSE.
#' @param plot_title String. Desired plot title.
#' @param title_size Numeric. Title font size.
#' @param x_text_size Numeric. x-axis tick font size.
#' @param y_text_size Numeric. y-axis tick font size.
#' @param x_title_size Numeric. x-axis title font size.
#' @param strip_text_size Numeric. y-axis geneset annotation font size.
#' @param leg_title_size Numeric. Legend title font size.
#' @param leg_text_size Numeric. Legend text font size.
#'
#' @return ggplot object
#' @export
#'


berlin_dotplot <- function(object = NULL, geneset = NULL, cluster_col = "seurat_clusters", assay = "RNA", min.pct = 0.25,
                           adjpval_cutoff = 0.05, log2fc_cutoff = 0, verbose = TRUE, plot_title = NULL, title_size = 20,
                           x_text_size = 14, y_text_size = 14, x_title_size = 16, strip_text_size = 16,
                           leg_title_size = 16, leg_text_size = 14) {


  if (is.null(object)) stop("Please supply Seurat object as input.")

  if (!assay %in% names(object)) stop("Assay input is not found in object")
  SeuratObject::DefaultAssay(object) <- assay

  if (!is.null(cluster_col)) {
    if (!cluster_col %in% names(object[[]])) {
      stop("cluster_col input not found in object meta data.")
    } else {
      SeuratObject::Idents(object) <- cluster_col
      # get the number of cells in each cluster
      cluster_sizes <- table(SeuratObject::Idents(object))
      # only include clusters with more than 2 cells (FindMaker will only work on clusters that have at least 3 cells )
      clusters <- names(cluster_sizes)[cluster_sizes > 2]
      if (verbose) {
        base::message(paste0("Finding markers for cluster column: ",cluster_col))
      }
      cluster_marker_list <- lapply(clusters, function(x) {
        markers <- Seurat::FindMarkers(object = object, ident.1 = x, min.pct = min.pct )
        upreg_markers <- rownames(markers[markers$p_val_adj < adjpval_cutoff & markers$avg_log2FC > log2fc_cutoff,])
        upreg_markers
      })
      names(cluster_marker_list) <- paste0(cluster_col,"_",clusters)

      if (is.data.frame(geneset)) {
        if (ncol(geneset) == 2) {
          gs <- split(geneset[,2],geneset[,1])
        }
      } else {
        gs <- geneset
        geneset <- stack(gs)
        geneset <- geneset[,c(2,1)]
        colnames(geneset) <- c("term","gene")
      }
      if (verbose) {
        base::message("Intersecting cluster markers with geneset.")
      }
      cluster_marker_list_gsint <- lapply(cluster_marker_list, function(x) lapply(gs, function(y) intersect(x,y)))

      dp <- Seurat::DotPlot(object = object, features = unique(geneset[,2]))
      dp_data <- dp$data

      dp_data_gs <- merge(geneset,dp_data, by.x = "gene",by.y = "features.plot")

      dotplot <- dp_data_gs %>%
        ggplot2::ggplot(aes(x=id, y = gene, color = avg.exp.scaled, size = pct.exp)) +
        ggplot2::geom_point() +
        ggplot2::theme_bw() +
        ggplot2::scale_colour_gradient2(low = "dark blue", high = "dark red", mid = "white") +
        ggplot2::theme(axis.text.x = element_text(size = x_text_size),
                       axis.text.y = element_text(size = y_text_size),
                       axis.title.x = element_text(size = x_title_size),
                       axis.title.y = element_blank(),
                       plot.title = element_text(size = title_size),
                       legend.title = element_text(size = leg_title_size),
                       legend.text = element_text(size = leg_text_size),
                       strip.placement = "outside",
                       strip.background.y = element_blank(),
                       strip.text.y.left = element_text(face = "bold", hjust = 1, size = strip_text_size, angle = 0)) +
        ggplot2::scale_y_discrete(position = "left") +
        ggplot2::facet_grid(term ~ ., switch = "y", scales = "free", space = "free") +
        ggplot2::xlab(cluster_col) +
        ggplot2::ggtitle(plot_title)

      return(dotplot)
    }
  }
}
