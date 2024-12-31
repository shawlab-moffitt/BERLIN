


#' List available ScType tissue gene sets
#'
#' @return Printed text.
#' @export
#'
#' @examples
#' ScType_genesets()


ScType_genesets <- function() {
  ScType_Pos_Genesets <- c("adrenal","brain","eye","heart","immune","intestine","kidney",
                           "liver","lung","muscle","pancreas","spleen","stomach","thymus")
  ScType_Neg_Genesets <- c("immune")
  cat("Positive Indicators:\n")
  print(ScType_Pos_Genesets)
  cat("Negative Indicators:\n")
  print(ScType_Neg_Genesets)
}


#' Generate ScType encrichment score matrix from Seurat object
#'
#' @param object Seurat object.
#' @param scaled Boolean. TRUE if the data has already been scaled. If the data still needs to be scaled, FALSE and the function will scale the data.
#' @param geneset String. Input a tissue geneset name from the provided in the package (see details). Run ScType_genesets() function to view available genesets.
#' @param pos_geneset List object. A list of positive marker genesets. Each element in the list should be a vector of genes and the name is the geneset name.
#' @param neg_geneset List object. A list of negative marker genesets. Each element in the list should be a vector of genes and the name is the geneset name.
#' @param verbose Boolean. To show progress, TRUE, else FALSE.
#'
#' @return Matrix
#' @export
#'


berlin_sctype <- function(object = NULL, scaled = TRUE, geneset = "immune", pos_geneset = NULL, neg_geneset = NULL, verbose = TRUE) {

  # Check input data
  if (is.null(object)) stop("Please supply Seurat object as input.")
  if (scaled) {
    if ("scale.data" %in% SeuratObject::Layers(object)) {
      data <- object@assays$RNA$scale.data
    } else {
      stop("Layer 'scale.data' is empty. Add scaled data to set 'scaled' argument to FALSE.")
    }
  } else {
    data <- object@assays$RNA$data
  }
  # If scaled, check for scaled data
  if (scaled & (sum(dim(data)) == 0)) stop("Scaled argument is TRUE but the dimensions of the scaled data are zero.\nChange scale to FALSE or add scaled data.")
  # If not scaled, check for unscaled data
  if (!scaled & (sum(dim(data)) == 0)) stop("Matrix data not found in object.")
  if (!scaled) {
    if (verbose) {
      message("Scaling data")
    }
    scaleddata <- Seurat::ScaleData(object, vars.to.regress = c("nFeature_RNA", "percent.mt"), verbose = verbose)
    data <- scaleddata@assays$RNA$scale.data
  }

  # Get geneset for ScType
  ScType_Pos_Genesets <- c("adrenal","brain","eye","heart","immune","intestine","kidney",
                           "liver","lung","muscle","pancreas","placenta","spleen","stomach","thymus")
  ScType_Neg_Genesets <- c("immune")
  # If geneset argument not in set genesets
  if (!tolower(geneset) %in% c(ScType_Pos_Genesets,ScType_Neg_Genesets)) {
    # if user input genesets are also null
    if (is.null(pos_geneset) & is.null(neg_geneset) & !exists(geneset)) {
      stop("Geneset argument not found.\nPlease see available genesets and indicators with the ScType_genesets() funciton or\ninput your own positve and negative marker genesets.")
    # If user input genesets are NOT null
    } else if (!is.null(pos_geneset) | !is.null(neg_geneset)) {
      # If users positive marker geneset exists, load it in
      if (exists(pos_geneset)) {
        gs_pos <- pos_geneset
      } else {
        # If users positive marker geneset does NOT exists
        stop("Positive marker geneset indicated, but not found.")
      }
      # If users negative marker geneset exists, load it in
      if (exists(neg_geneset)) {
        gs_neg <- neg_geneset
      } else {
        # If users negative marker geneset does NOT exists
        stop("Negative marker geneset indicated, but not found.")
      }
    }
  # If geneset argument in set genesets
  } else {
    # If geneset selected is immune markers
    if (tolower(geneset) == "immune") {
      data(list = paste0("ScType_immune_system_PosIndicator_Geneset_DB"), package = "BERLIN")
      assign("gs_pos",get(paste0("ScType_immune_system_PosIndicator_Geneset_DB")))
      data(list = paste0("ScType_immune_system_NegIndicator_Geneset_DB"), package = "BERLIN")
      assign("gs_neg",get(paste0("ScType_immune_system_NegIndicator_Geneset_DB")))
    } else {
    # If geneset selected is NOT immune markers
      data(list = paste0("ScType_",tolower(geneset),"_PosIndicator_Geneset_DB"), package = "BERLIN")
      assign("gs_pos",get(paste0("ScType_",tolower(geneset),"_PosIndicator_Geneset_DB")))
      gs_neg <- NULL
    }
  }

  gs <- c(gs_pos,gs_neg)

  # Get marker sensitivity from genesets
  marker_stat = sort(table(unlist(gs_pos)), decreasing = T)
  marker_sensitivity = data.frame(score_marker_sensitivity = scales::rescale(as.numeric(marker_stat),to = c(0,1), from = c(length(gs_pos),1)),
                                  gene = names(marker_stat), stringsAsFactors = FALSE)

  # Keep only gs genes that are found in the input data
  intsercted_genes <- intersect(rownames(data),unique(unname(unlist(gs))))
  gs_pos <- lapply(gs_pos, function(x){
    return(x[x %in% intsercted_genes])
  })
  gs_neg <- lapply(gs_neg, function(x){
    return(x[x %in% intsercted_genes])
  })
  cell_markers_genes_score <- marker_sensitivity[which(marker_sensitivity[,2] %in% intsercted_genes),]
  cell_markers_genes_score <- cell_markers_genes_score[order(cell_markers_genes_score[,2]),]
  rownames(cell_markers_genes_score) <- cell_markers_genes_score[,2]

  data_gs <- data[intsercted_genes,]
  data_gs <- data_gs[order(rownames(data_gs)),]

  data_gs_weighted <- data_gs * cell_markers_genes_score[,1]

  # combine scores
  es = do.call("rbind",lapply(names(gs_pos), function(gss_){
    sapply(1:ncol(data_gs_weighted), function(j) {
      # For each geneset, take the genes of that geneset
      # Apply the genes to each barcode to get the scores
      gs_z = data_gs_weighted[gs_pos[[gss_]], j]            # Get the score for the positive markers
      gz_2 = data_gs_weighted[gs_neg[[gss_]], j] * -1       # Get the score for the negative markers and multiply by -1 ??
      sum_t1 = sum(gs_z, na.rm = T) / sqrt(length(gs_z))    # Get some sort of average score for the barcode of the positive and negative markers ??
      sum_t2 = sum(gz_2, na.rm = T) / sqrt(length(gz_2))
      if(is.na(sum_t2)){                                    # if no negative markers set to 0
        sum_t2 = 0
      }
      sum_t1 + sum_t2                                       # get the sum of the positive and negative scores
    })                                                      # This is done across all barcodes
  })                                                        # For each geneset
  )                                                         # The geneset scores for each barcode are rbinded to a line to make a data frame

  dimnames(es) = list(names(gs_pos), colnames(data_gs_weighted))
  es.max <- es[stats::complete.cases(es),]
  return(es.max)

}



#' Assign cluster to scType
#'
#' @param object Seurat object.
#' @param score Matrix. Results from berlin_scType function, of a matrix with rownames as scTypes and columns of each cell and their associated scType enrichment score.
#' @param meta data.frame. If no Seurat object available, please provide meta data with the cluster column and rownames matching score cell names.
#' @param cluster_col String. Column name of meta data column that annotated clusters. Defaults to "seurat_clusters".
#' @param scType_col String. Desired column name for column to store scType classification results. Defaults to "scType_classification".
#' @param ncell.pct Numeric. Checks whether the total scType score sum is greater than a fraction of the number of cells in that cluster. Number between 0 and 1, defaults to 0.25.
#' @param verbose Boolean. To show progress, TRUE, else FALSE.
#'
#' @return Seurat object with appended meta data.
#' @export
#'
#'


berlin_sctype_classify <- function(object = NULL, score = NULL, meta = NULL, cluster_col = "seurat_clusters",
                                   scType_cluster_col = NULL, scType_cell_col = NULL, ncell.pct = 0.25, verbose = TRUE) {

  call.string <- deparse(expr = sys.calls()[[1]])
  func_name <- "berlin_sctype_classify"
  time.stamp <- Sys.time()
  argg <- c(as.list(environment()))
  argg <- Filter(function(x) any(is.numeric(x) | is.character(x)), argg)

  if (is.null(object) & is.null(meta)) stop("Please provide Seurat object or meta data.")
  if (is.null(meta) & !is.null(object)) {
    meta <- object[[]]
  }
  if (is.null(cluster_col) & (!"seurat_clusters" %in% colnames(meta))) stop("Please provide cluster column name.")
  scType_cluster_col <- ifelse(is.null(scType_cluster_col),paste0("scType_",cluster_col,"_Classification"),scType_cluster_col)
  scType_cell_col <- ifelse(is.null(scType_cell_col),"scType_Cell_Classification",scType_cell_col)

  if (verbose) {
    message("Calculating scType cluster results.")
  }

  cL_results = do.call("rbind", lapply(unique(meta[,cluster_col]), function(cl){
    es.max.cl = sort(rowSums(score[ ,rownames(meta[meta[,cluster_col]==cl, ])]), decreasing = T)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(meta[,cluster_col]==cl)), 1)
  }))



  colnames(cL_results) <- c(cluster_col,scType_cluster_col,"scType Score Sum","ncells")
  rownames(cL_results) <- NULL

  cL_results[,scType_cluster_col][as.numeric(cL_results[,"scType Score Sum"]) < as.numeric(cL_results[,"ncells"])*ncell.pct] = "Unknown"

  meta <- cbind(barcode = rownames(meta),meta)
  col2move2 <- colnames(meta)[which(colnames(meta)==cluster_col)-1]
  meta2 <- merge(meta,cL_results, all.x = T, sort = F)
  meta2 <- meta2 %>%
    relocate(any_of(c(cluster_col)), .after = !!sym(col2move2)) %>%
    as.data.frame()


  if (verbose) {
    message("Calculating scType cell results.")
  }
  score_t <- as.data.frame(t(score))
  score_t[,scType_cell_col] <- apply(score_t,1,function(x) {
    paste(names(x[x==max(x)]), collapse = ", ")
  })
  scType_cell_class <- cbind(barcode = rownames(score_t),
                                  score_t[,scType_cell_col, drop = F])
  meta2 <- merge(meta2,scType_cell_class, all.x = T, sort = F)

  rownames(meta2) <- meta2[,1]
  meta2 <- meta2[,-1]
  meta3 <- meta2[rownames(meta),]

  object <- Seurat::AddMetaData(object,metadata = meta3)

  slot(object = object, name = "commands")[[func_name]] <- list(name = func_name,
                                                                time.stamp = time.stamp,
                                                                call.string = call.string,
                                                                params = argg)

  return(object)
}






#' Summarize scType cell annotation across clusters
#'
#' @param object Seurat object.
#' @param scType_cell_col String. Name of column containing predicted scType cell levele classification. If not supplied it will be predicted based on string matching.
#' @param cluster_col String. Name of column containing cluster information. Default to "seurat_clusters".
#' @param show_unassigned Boolean. If TRUE, unassigned or NA cell types will be included.
#'
#' @return A data.frame object.
#' @export
#'



berlin_sctype_summarize <- function(object = NULL, scType_cell_col = NULL, cluster_col = "seurat_clusters", show_unassigned = TRUE) {

  if (is.null(object)) stop("Please supply Seurat object")

  meta <- object[[]]
  if (is.null(scType_cell_col)) {
    scType_cell_col <- grep("^scType_Cell_Classification$", colnames(meta), value = T)
    if (length(celltype_col) > 0) {
      message(paste0("'scType_cell_col' argument is NULL\nWill use identified possible cell type column(s) '",paste0(scType_cell_col, collapse = ", "),"'"))
    } else {
      stop(paste0("'scType_cell_col' argument is NULL\nCannot indentify cell type column\nPlease provide"))
    }
  }

  if (is.null(cluster_col)) stop("Please supply cluster column name")

  cluster_celltype_df_reshape_list <- lapply(scType_cell_col,function(c) {
    cluster_celltype_df <- data.frame(
      cluster = meta[,cluster_col],
      cell_type = meta[,c]
    )
    # Convert NA values to a string so they are treated as a category
    cluster_celltype_df$cell_type <- ifelse(is.na(cluster_celltype_df$cell_type), "Unassigned/NA", cluster_celltype_df$cell_type)
    # Generate a table of cluster vs cell type counts
    cluster_celltype_count <- table(cluster_celltype_df$cluster, cluster_celltype_df$cell_type)
    # Convert counts to percentages
    cluster_celltype_percentage <- prop.table(cluster_celltype_count, margin = 1) * 100
    # Convert to df for easier viewing
    cluster_celltype_df <- as.data.frame(cluster_celltype_percentage)
    colnames(cluster_celltype_df) <- c('cluster', 'scType', 'percentage')
    # Reshape for easy viewing
    cluster_celltype_df_reshape <- reshape2::dcast(cluster_celltype_df, scType ~ cluster, value.var = "percentage")
    colnames(cluster_celltype_df_reshape)[-1] <- paste0("Percent_scTypeCells_",cluster_col,"_",colnames(cluster_celltype_df_reshape)[-1])
    if (!show_unassigned) {
      cluster_celltype_df_reshape <- cluster_celltype_df_reshape[which(cluster_celltype_df_reshape[,1] != "Unassigned/NA"),]
    }
    cluster_celltype_df_reshape <- data.frame(cluster_celltype_df_reshape)
  })
  cluster_celltype_df_reshape_out <- data.table::rbindlist(cluster_celltype_df_reshape_list,fill = T)

  return(cluster_celltype_df_reshape_out)

}
