



#' Perform differential gene expression on clusters in a Seurat object
#'
#' @param object Seurat object. If NULL, must supply count data.
#' @param assay String. Name of the initial/default assay. Default is RNA.
#' @param cluster_cols String or Vector. Column name from meta data of primary identified clusters. Default to 'seurat_clusters'. Multiple columns allowed.
#' @param min.pct Numeric. Only test genes that are detected in a minimum fraction of cells. Default is 0.25.
#' @param adjpval_cutoff Numeric. An adjusted p-value value cutoff for identifying significantly expressed genes. Default is 0.05.
#' @param log2fc_cutoff Numeric. An absolute value of a log2 fold change cutoff to identify up and down regulated genes. Default is 0.
#' @param verbose Boolean. To show progress, TRUE, else FALSE.
#' @param min.diff.pct Numeric. Only test genes that show a minimum difference in the fraction of detection between the two groups. Set to 0 by default.
#' @param return_seurat Boolean. If TRUE, a Seurat object will be the output with the results within the stored in the Misc slot.
#'
#' @return A list of data frames containing DGE of clusters and a geneset of the markers of each cluster.
#' @export
#'
#'


berlin_cluster_deg <- function(object = NULL, assay = "RNA", cluster_cols = "seurat_clusters", verbose = TRUE,
                               min.pct = 0.25, min.diff.pct = 0, adjpval_cutoff = 0.05, log2fc_cutoff = 1, return_seurat = FALSE) {


  call.string <- deparse(expr = sys.calls()[[1]])
  func_name <- "berlin_cluster_deg"
  time.stamp <- Sys.time()
  argg <- c(as.list(environment()))
  argg <- Filter(function(x) any(is.numeric(x) | is.character(x)), argg)

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

        if (verbose) {
          base::message(paste0("Finding all unfiltered markers on cluster column: ",cluster_col))
        }
        unfiltered_all_markers <- Seurat::FindAllMarkers(object = object, min.pct = 0, min.diff.pct = 0, logfc.threshold = 0, verbose = verbose)
        rownames(unfiltered_all_markers) <- NULL
        unfiltered_all_markers <- unfiltered_all_markers %>%
          relocate(cluster,gene)
        unfiltered_all_markers$cluster <- paste0(cluster_col,"_",unfiltered_all_markers$cluster)

        unfiltered_pct_markers <- unfiltered_all_markers[which(unfiltered_all_markers$pct.1 >= min.pct | unfiltered_all_markers$pct.2 >= min.pct),]

        upreg_markers <- unfiltered_pct_markers[unfiltered_pct_markers$p_val_adj < adjpval_cutoff & unfiltered_pct_markers$avg_log2FC > log2fc_cutoff,]
        dnreg_markers <- unfiltered_pct_markers[unfiltered_pct_markers$p_val_adj < adjpval_cutoff & unfiltered_pct_markers$avg_log2FC < log2fc_cutoff,]

        upreg_markers_gs <- data.frame(term = paste0(upreg_markers$cluster,"_UpReg"),gene = upreg_markers$gene)
        dnreg_markers_gs <- data.frame(term = paste0(dnreg_markers$cluster,"_DownReg"),gene = dnreg_markers$gene)

        markers_gs <- rbind(upreg_markers_gs,dnreg_markers_gs)
        cluster_marker_list[[cluster_col]] <- unfiltered_all_markers

        SeuratObject::Misc(object = object, slot = paste0(cluster_col,"_FindMarkers_Results")) <- as.data.frame(unfiltered_all_markers)
        SeuratObject::Misc(object = object, slot = paste0(cluster_col,"_FindMarkers_Geneset")) <- as.data.frame(markers_gs)
      }
    }
  }
  cluster_marker_list_out <- data.table::rbindlist(cluster_marker_list, idcol = "Cluster_Column")
  output_list <- list(FindMarkers_Results = cluster_marker_list_out,
                      FindMarkers_Geneset = markers_gs)
  if (!return_seurat) {
    return(output_list)
  } else {
    slot(object = object, name = "commands")[[func_name]] <- list(name = func_name,
                                                                  time.stamp = time.stamp,
                                                                  call.string = call.string,
                                                                  params = argg)
    return(object)
  }
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
#' @param markers_res Data.frame. Result data from Seurat's FindAllMarkers() function.
#' @param return_seurat Boolean. If TRUE, a Seurat object will be the output with the results within the stored in the Misc slot.
#'
#' @return Data.frame
#' @export
#'





berlin_intersect <- function(object = NULL, markers_res = NULL,geneset = NULL, cluster_col = "seurat_clusters", assay = "RNA",
                             min.pct = 0.25, adjpval_cutoff = 0.05, log2fc_cutoff = 1.0, verbose = TRUE, return_seurat = FALSE) {

  call.string <- deparse(expr = sys.calls()[[1]])
  func_name <- "berlin_intersect"
  time.stamp <- Sys.time()
  argg <- c(as.list(environment()))
  argg <- Filter(function(x) any(is.numeric(x) | is.character(x)), argg)

  if (is.null(object) & is.null(markers_res)) stop("Please supply Seurat object or FindAllMarkers result as input.")

  if (!assay %in% names(object)) stop("Assay input is not found in object")
  SeuratObject::DefaultAssay(object) <- assay

  if (is.null(markers_res)) {
    # Object with find markers res
    ## Check if find marker result is in object
    fm_res_name <- grep("FindMarkers_Results",names(Misc(object = object)), value = T)
    if (length(fm_res_name) > 0) {
      markers_df <- as.data.frame(Misc(object = object, slot = fm_res_name))
      if (any(grepl(paste0("^",cluster_col,"_"),markers_df[,1]))) {
        markers_df_cl <- markers_df[which(grepl(paste0("^",cluster_col,"_"),markers_df[,1])),]
        unfiltered_all_markers <- markers_df_cl
      } else {
        stop(paste0(cluster_col," markers not found in FindMarkers_Results within Seurat object"))
      }
    } else {
      # Object without find markers res
      if (verbose) {
        base::message(paste0("Finding all unfiltered markers on cluster column: ",cluster_col))
      }
      SeuratObject::Idents(object) <- cluster_col
      unfiltered_all_markers <- Seurat::FindAllMarkers(object = object, min.pct = 0, min.diff.pct = 0, logfc.threshold = 0, verbose = verbose)
      rownames(unfiltered_all_markers) <- NULL
      unfiltered_all_markers <- unfiltered_all_markers %>%
        relocate(cluster,gene)
      unfiltered_all_markers$cluster <- paste0(cluster_col,"_",unfiltered_all_markers$cluster)
    }
  } else {
    # find markers res provided
    if (!"cluster" %in% colnames(markers_res)) stop("'cluster' column not found in FindAllMarkers result table.")
    unfiltered_all_markers <- markers_res
  }

  unfiltered_pct_markers <- unfiltered_all_markers[which(unfiltered_all_markers$pct.1 >= min.pct | unfiltered_all_markers$pct.2 >= min.pct),]
  sig_upreg_markers <- unfiltered_pct_markers[unfiltered_pct_markers$p_val_adj < adjpval_cutoff & unfiltered_pct_markers$avg_log2FC > log2fc_cutoff,]
  clusters <- unique(sig_upreg_markers$cluster)

  if (is.data.frame(geneset)) {
    if (ncol(geneset) == 2) {
      colnames(geneset) <- c("term","gene")
    }
  } else {
    geneset <- stack(geneset)
    geneset <- geneset[,c(2,1)]
    colnames(geneset) <- c("term","gene")
  }

  geneset_sig_upreg <- merge(sig_upreg_markers,geneset)[,c(1,2,8)]
  clusters_miss <- setdiff(clusters,unique(geneset_sig_upreg$cluster))
  if (length(clusters_miss) > 0) {
    geneset_sig_upreg <- rbind(geneset_sig_upreg,
                               data.frame(gene = NA,
                                          cluster = clusters_miss,
                                          term = unique(geneset$term)))
  }

  geneset_sig_upreg$cluster <- factor(geneset_sig_upreg$cluster, levels = clusters)
  geneset_sig_upreg <- geneset_sig_upreg %>%
    group_by(cluster,term) %>%
    summarise(genes = paste0(gene[!is.na(gene)], collapse = ", "), .groups = "drop") %>%
    mutate(genes = replace(genes,genes=="",NA)) %>%
    unique() %>%
    pivot_wider(
      id_cols = cluster,
      names_from = term,
      values_from = genes
    ) %>%
    as.data.frame()

  if (!return_seurat) {
    return(geneset_sig_upreg)
  } else {
    SeuratObject::Misc(object = object, slot = paste0(cluster_col,"_Sig_UpReg_Intersect")) <- as.data.frame(geneset_sig_upreg)
    slot(object = object, name = "commands")[[func_name]] <- list(name = func_name,
                                                                  time.stamp = time.stamp,
                                                                  call.string = call.string,
                                                                  params = argg)
    return(object)
  }
}


#' Generate dotplot heatmap of geneset gene expression within clusters
#'
#' @param object Seurat object. If NULL, must supply count data.
#' @param geneset Data.frame or list. Two column data frame with first column consisting of geneset names and second column of the genes withing that geneset. One gene per row, the geneset name can repeat in the first column. If list input, this would be a named list where the names are the geneset names and each element a vector of genes as strings.
#' @param cluster_col String. Column name from meta data of primary identified clusters. Default to 'seurat_clusters'.
#' @param assay String. Name of the initial/default assay. Default is RNA.
#' @param min.pct Numeric. Only test genes that are detected in a minimum fraction of cells. Default is 0.25.
#' @param scaled Boolean. TRUE if plotting scaled expression data, FALSE if un-scaled.
#' @param exp_binary Boolean. If TRUE, size of dot will be binary, large if significantly upregulated, small if not. If FALSE, dot size will be based on percent of cells expressing gene.
#' @param log_exp Boolean. Only relevant if 'scaled' argument equals FALSE. If TRUE, average expression value will be logged with log2.
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
#' @param markers_res Data.frame. Result data from Seurat's FindAllMarkers() function.
#' @param return_seurat Boolean. If TRUE, a Seurat object will be the output with the results within the stored in the Misc slot.
#'
#' @return ggplot object
#' @export
#'


berlin_dotplot <- function(object = NULL, markers_res = NULL, geneset = NULL, cluster_col = "seurat_clusters", assay = "RNA", min.pct = 0.25,
                           scaled = FALSE, exp_binary = TRUE, log_exp = TRUE, adjpval_cutoff = 0.05, log2fc_cutoff = 1.0, verbose = TRUE,
                           return_seurat = FALSE, plot_title = NULL, title_size = 20, x_text_size = 14, y_text_size = 14,
                           x_title_size = 16, strip_text_size = 16, leg_title_size = 16, leg_text_size = 14) {

  call.string <- deparse(expr = sys.calls()[[1]])
  func_name <- "berlin_dotplot"
  time.stamp <- Sys.time()
  argg <- c(as.list(environment()))
  argg <- Filter(function(x) any(is.numeric(x) | is.character(x)), argg)

  if (is.null(object) & is.null(markers_res)) stop("Please supply Seurat object or FindAllMarkers result as input.")

  if (!assay %in% names(object)) stop("Assay input is not found in object")
  SeuratObject::DefaultAssay(object) <- assay

  #if (!is.null(cluster_col)) {
    #if (!cluster_col %in% names(object[[]])) {
    #  stop("cluster_col input not found in object meta data.")
    #} else {


      if (!assay %in% names(object)) stop("Assay input is not found in object")
      SeuratObject::DefaultAssay(object) <- assay

      if (is.null(markers_res)) {
        # Object with find markers res
        ## Check if find marker result is in object
        fm_res_name <- grep("FindMarkers_Results",names(Misc(object = object)), value = T)
        if (length(fm_res_name) > 0) {
          markers_df <- as.data.frame(Misc(object = object, slot = fm_res_name))
          if (any(grepl(paste0("^",cluster_col,"_"),markers_df[,1]))) {
            markers_df_cl <- markers_df[which(grepl(paste0("^",cluster_col,"_"),markers_df[,1])),]
            unfiltered_all_markers <- markers_df_cl
          } else {
            stop(paste0(cluster_col," markers not found in FindMarkers_Results within Seurat object"))
          }
        } else {
          # Object without find markers res
          if (verbose) {
            base::message(paste0("Finding all unfiltered markers on cluster column: ",cluster_col))
          }
          SeuratObject::Idents(object) <- cluster_col
          unfiltered_all_markers <- Seurat::FindAllMarkers(object = object, min.pct = 0, min.diff.pct = 0, logfc.threshold = 0, verbose = verbose)
          rownames(unfiltered_all_markers) <- NULL
          unfiltered_all_markers <- unfiltered_all_markers %>%
            relocate(cluster,gene)
          unfiltered_all_markers$cluster <- paste0(cluster_col,"_",unfiltered_all_markers$cluster)
        }
      } else {
        # find markers res provided
        if (!"cluster" %in% colnames(markers_res)) stop("'cluster' column not found in FindAllMarkers result table.")
        unfiltered_all_markers <- markers_res
      }


      unfiltered_pct_markers <- unfiltered_all_markers[which(unfiltered_all_markers$pct.1 >= min.pct | unfiltered_all_markers$pct.2 >= min.pct),]
      sig_upreg_markers <- unfiltered_pct_markers[unfiltered_pct_markers$p_val_adj < adjpval_cutoff & unfiltered_pct_markers$avg_log2FC > log2fc_cutoff,]
      clusters <- unique(sig_upreg_markers$cluster)

      if (is.data.frame(geneset)) {
        if (ncol(geneset) == 2) {
          colnames(geneset) <- c("term","gene")
        }
      } else {
        geneset <- stack(geneset)
        geneset <- geneset[,c(2,1)]
        colnames(geneset) <- c("term","gene")
      }

      geneset_sig_upreg <- merge(sig_upreg_markers,geneset)[,c(1,2,8)]
      clusters_miss <- setdiff(clusters,unique(geneset_sig_upreg$cluster))
      if (length(clusters_miss) > 0) {
        geneset_sig_upreg <- rbind(geneset_sig_upreg,
                                   data.frame(gene = NA,
                                              cluster = clusters_miss,
                                              term = unique(geneset$term)))
      }
      geneset_sig_upreg$cluster <- factor(geneset_sig_upreg$cluster, levels = clusters)
      geneset_sig_upreg[,"Significantly Upregulated"] <- TRUE

      ## Reformat intersected list to data frame
      #cluster_gs_overlap <- as.data.frame(data.table::rbindlist(lapply(cluster_marker_list_gsint,function(x) {
      #  stack(x)
      #}), fill = T, use.names = T, idcol = "Cluster"))
      #colnames(cluster_gs_overlap) <- c("Cluster","Gene.Symbol","Annotation")
      #cluster_gs_overlap[,"Significantly Upregulated"] <- TRUE

      # Get overall dp plot data
      #dp <- Seurat::DotPlot(object = object, features = unique(cluster_gs_overlap[,2]))
      dp <- Seurat::DotPlot(object = object, features = unique(geneset_sig_upreg$gene))
      dp_data <- dp$data
      dp_data$cluster <- paste0(cluster_col,"_",dp_data$id)
      dp_data$gene <- dp_data[,"features.plot"]

      # Merge overall data with geneset annotation
      dp_data_gs <- merge(geneset,dp_data, all.y = T)

      # Merge with sig upreg markers and annotation
      dp_data_gs2 <- merge(geneset_sig_upreg,dp_data_gs, all.y = T)
      dp_data_gs2[,"Significantly Upregulated"] <- ifelse(is.na(dp_data_gs2[,"Significantly Upregulated"]),FALSE,dp_data_gs2[,"Significantly Upregulated"])
      dp_data_gs2[,"Significantly Upregulated"] <- factor(dp_data_gs2[,"Significantly Upregulated"], levels = c("TRUE","FALSE"))

      # Select data to plot based on user input
      expr_col <- ifelse(scaled,"avg.exp.scaled","avg.exp")
      if (!scaled) {
        if (log_exp) {
          dp_data_gs2[,paste0(expr_col,"_log2")] <- log2(dp_data_gs2[,expr_col]+1)
          expr_col <- paste0(expr_col,"_log2")
        }
      }
      pct_col <- ifelse(exp_binary,"Significantly Upregulated","pct.exp")

      # Plot
      dotplot <- dp_data_gs2 %>%
        ggplot2::ggplot(ggplot2::aes(x=id, y = gene, color = !!sym(expr_col), size = !!sym(pct_col))) +
        ggplot2::geom_point() +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(size = x_text_size),
                       axis.text.y = ggplot2::element_text(size = y_text_size),
                       axis.title.x = ggplot2::element_text(size = x_title_size),
                       axis.title.y = ggplot2::element_blank(),
                       plot.title = ggplot2::element_text(size = title_size),
                       legend.title = ggplot2::element_text(size = leg_title_size),
                       legend.text = ggplot2::element_text(size = leg_text_size),
                       strip.placement = "outside",
                       strip.background.y = ggplot2::element_blank(),
                       strip.text.y.left = ggplot2::element_text(face = "bold", hjust = 1, size = strip_text_size, angle = 0)) +
        ggplot2::scale_y_discrete(position = "left") +
        ggplot2::facet_grid(term ~ ., switch = "y", scales = "free", space = "free") +
        ggplot2::xlab(cluster_col) +
        ggplot2::ggtitle(plot_title)
      if (exp_binary) {
        dotplot <- dotplot +
          ggplot2::scale_size_manual(values = c("TRUE"=5,"FALSE"=1))
      }
      if (scaled) {
        dotplot <- dotplot +
          ggplot2::scale_colour_gradient2(low = "dark blue", high = "dark red", mid = "white")
      } else {
        dotplot <- dotplot +
          ggplot2::scale_colour_gradient(low = "#a1d6ff", high = "#132B43")
      }

      if (!return_seurat) {
        return(dotplot)
      } else {
        SeuratObject::Misc(object = object, slot = paste0(cluster_col,"_Sig_UpReg_DotPlot")) <- dotplot
        slot(object = object, name = "commands")[[func_name]] <- list(name = func_name,
                                                                      time.stamp = time.stamp,
                                                                      call.string = call.string,
                                                                      params = argg)
        return(object)
      }

    #}
  #}
}
