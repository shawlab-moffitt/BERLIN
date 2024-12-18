


#' Generate scGate model from parameter data
#'
#' @param params A data.frame. The scGate model will be generated based on these parameters. More information about the format in details.
#' @param signatures A two column data frame. First column of signature names and second column of signature featurs/genes.
#'
#' @return A R List Object.
#' @export
#'
#' @details
#' The first column should indicate the cell type label. The following columns should be named to indicate the level and positive or negative indication of the signature.
#' The level should be a numeric number starting at 1 and the positive and negative indication will be determined by grepping the column name for 'pos' or 'neg'.
#' For example the column name 'Level2Pos' indicated a column of level 2 positive indication signatures. The user does not need to include the word 'level' and
#' can write out the word of 'positive' or 'negative' and put the level number after the indicator if preferred.
#'
#'



berlin_scgate_model <- function(params = NULL, signatures = NULL) {
  if (is.null(params)) stop("Please supply scGate model parameters")
  posneg_cols <- grep("pos|neg",colnames(params), ignore.case = T)
  if (length(posneg_cols) == 0) stop("Positive/Negative indicators not found in column names")
  leveln_cols <- grep("[0-9]",colnames(params))
  if (length(leveln_cols) == 0) stop("Level indicators not found in column names")

  if (is.null(signatures)) {
    data(list = paste0("scGate_Signatures"), package = "BERLIN")
    assign("scGate_Signatures",get(paste0("scGate_Models_Geneset")))
  }

  model_list <- lapply(seq_along(params[,1]), function(x) {

    # Get cell type
    cell_type <- params[x,1]

    # Apply function on each level/pos/neg
    temp_mod <- do.call("rbind",sapply(colnames(params)[-1], function(y) {
      leveln <- paste0("level",gsub("\\D+", "", y))
      posneg <- ifelse(grepl("pos",y,ignore.case = T),"positive","negative")
      sig_name <- strsplit(gsub(" ","",params[x,y]),",")[[1]]

      # Apply function on each signature name (sig_name)
      # Output is 4 values to make row of df
      rowz <- do.call("rbind",sapply(sig_name, function(z) {
        genes <- paste0(scGate_Signatures[which(scGate_Signatures[,1] == z),2], collapse = ";")
        rowmade <- c(leveln,posneg,z,genes)
        return(as.data.frame(rowmade))

      }))

    }))

    # Clean cell type model data frame
    temp_mod <- as.data.frame(temp_mod)
    colnames(temp_mod) <- c("levels","use_as","name","signature")
    rownames(temp_mod) <- NULL
    temp_mod[temp_mod==""] <- NA
    temp_mod <- temp_mod[complete.cases(temp_mod),]
    return(temp_mod)

  })
  names(model_list) <- params[,1]
  return(model_list)
}







#' Performs scGate function with a single model on a Seurat object
#'
#' @param object Seurat object. If NULL, must supply count data.
#' @param model A single scGate model, or a list of scGate models. Please see details for formatting.
#' @param model_name String. A name for the input model.
#' @param assay String. Name of the initial/default assay. Default to RNA.
#' @param ncores Numeric. Number of processors for parallel processing. Default to 1.
#' @param seed Numeric. Random seed, default is 42. Setting to NULL will remove seed.
#' @param verbose Boolean. To show progress, TRUE, else FALSE.
#'
#' @return A data.frame.
#' @export
#'


berlin_run_scgate <- function(object = NULL, model = NULL, model_name = NULL, assay = "RNA", pos.thr = 0.2, neg.thr = 0.2, ncores = 1, seed = 42, verbose = TRUE) {
  if (is.null(object)) stop("Please supply Seurat object")
  if (is.null(model)) stop("Please supply scGate model")
  if (is.data.frame(model)) {
    if (all(colnames(model) != c("levels","use_as","name","signature"))) stop("Please check model format.\nMust follow column names: levels, use_as, name, signature.")
  }
  if (is.data.frame(model) | is.data.frame(model[[1]])) {
    if (is.null(model_name)) stop("Please provide model name.")
  }

  original_meta_col_n <- ncol(object[[]])

  object_scgate <- scGate::scGate(object, model, assay = assay, pos.thr = pos.thr, neg.thr = neg.thr, ncores = ncores, seed = seed, verbose = verbose)

  meta <- object_scgate[[]]
  new_meta_cols_n <- ncol(meta)
  new_col_name <- paste0("scGate_", model_name, "_Celltype")
  colnames(meta)[c((original_meta_col_n+1):ncol(meta))] <- paste0(colnames(meta)[c((original_meta_col_n+1):ncol(meta))],new_col_name)
  pure_columns <- grep("^is.pure", colnames(meta), value = TRUE)
  if (length(pure_columns) == 1) {
    meta[[new_col_name]] <- ifelse(meta[,pure_columns] == "Pure",model_name,NA)
  } else {
    pure_temp <- sapply(pure_columns, function(c) {
      feature <- sub("^is.pure_", "", c)
      feature <- sub(new_col_name, "", feature)
      col <- meta[,c]
      new_col <- ifelse(meta[,c] == "Pure",feature,NA)
      return(new_col)
    })
    pure_comb <- apply(pure_temp,1, function(x) {
      paste(x[!is.na(x)],collapse = ", ")
    })
    meta[,new_col_name] <- pure_comb
    meta[which(meta[,new_col_name] == ""),new_col_name] <- NA
  }
  return(meta)

}


#' Performs scGate function on a Seurat object with one or multiple models
#'
#' @param object Seurat object. If NULL, must supply count data.
#' @param assay String. Name of the initial/default assay. Default to RNA.
#' @param model List object or data.frame. A single scGate model, or a list of scGate models. Please see \link[scGate]{scGate} for further formatting parameters.
#' @param model_name String. A desired name for the input model. Required if only running a single model.
#' @param species String. Denoting if the data is from "human" or "mouse". Required if using package provided models.
#' @param ncores Numeric. Number of processors for parallel processing. Default to 1.
#' @param seed Numeric. Random seed, default is 42. Setting to NULL will remove seed.
#' @param verbose Boolean. To show progress, TRUE, else FALSE.
#'
#' @return Seurat object.
#' @export
#'
#'
#'


berlin_scgate <- function(object = NULL, assay = "RNA", model = NULL, model_name = NULL, pos.thr = 0.2, neg.thr = 0.2,
                          species = "human", ncores = 1, seed = 42, verbose = TRUE) {

  if (is.null(object)) stop("Please supply Seurat object")

  if (is.null(model)) {
    # detect species
    species_detected <- detect_species(Features(object))
    # if different than input species, notify user
    if (species_detected != species) {
      message(paste0("Species detected does not equeal the species argument. Will be treating data as ",species_detected," data."))
    }
    scGate_models_DB <-  get_scGateDB()
    model <- scGate_models_DB[[species_detected]]
  } else {
    if (is.data.frame(model)) {
      if (all(colnames(model) != c("levels","use_as","name","signature"))) stop("Please check model format.\nMust follow column names: levels, use_as, name, signature.")
    }
  }

  # Create a list to store new meta data
  if (is.data.frame(model) | is.data.frame(model[[1]])) {
    new_meta <- berlin_run_scgate(object = object, model = model, model_name = model_name, assay = assay, ncores = ncores, seed = seed, verbose = verbose)
  } else {
    model_names <- names(model)
    new_metas <- sapply(model_names, function(m) {
      model_in <- model[[m]]
      new_meta <- berlin_run_scgate(object = object, model = model_in, model_name = m, assay = assay, ncores = ncores, seed = seed, verbose = verbose)
      new_meta <- cbind(barcode = rownames(new_meta),new_meta)
      new_meta
    })
    new_meta <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, all = TRUE),
                       new_metas)
  }

  if (verbose) {
    message("Adding new meta data")
  }
  object <- Seurat::AddMetaData(object,metadata = new_meta)

  return(object)

}


#' Summarize scGate annotation across clusters
#'
#' @param object Seurat object. If NULL, must supply count data.
#' @param celltype_col String. Name of column containing predicted cell type names. If not supplied it will be predicted based on string matching.
#' @param cluster_col String. Name of column containing cluster information.
#'
#' @return A data.frame object.
#' @export
#'



berlin_scgate_summarize <- function(object = NULL, celltype_col = NULL, cluster_col = "seurat_clusters", show_unassigned = FALSE) {

  if (is.null(object)) stop("Please supply Seurat object")

  meta <- object[[]]
  if (is.null(celltype_col)) {
    celltype_col <- grep("^scGate_.*\\_Celltype$", colnames(meta), value = T)
    celltype_col <- grep("multiscGate",celltype_col, invert = T, value = T)
    if (length(celltype_col) > 0) {
      message(paste0("'celltype_col' argument is NULL\nWill use identified possible cell type column(s) '",paste0(celltype_col, collapse = ", "),"'"))
    } else {
      stop(paste0("'celltype_col' argument is NULL\nCannot indentify cell type column\nPlease provide"))
    }
  }

  if (is.null(cluster_col)) stop("Please supply cluster column name")

  cluster_celltype_df_reshape_list <- lapply(celltype_col,function(c) {
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
    colnames(cluster_celltype_df) <- c('cluster', 'cell_type', 'percentage')
    # Reshape for easy viewing
    cluster_celltype_df_reshape <- reshape2::dcast(cluster_celltype_df, cell_type ~ cluster, value.var = "percentage")
    colnames(cluster_celltype_df_reshape)[-1] <- paste0("Percent_Celltype_",cluster_col,"_",colnames(cluster_celltype_df_reshape)[-1])
    if (show_unassigned) {
      cluster_celltype_df_reshape <- cluster_celltype_df_reshape[which(cluster_celltype_df_reshape[,1] != "Unassigned/NA"),]
    }
    cluster_celltype_df_reshape <- cbind(data.frame(Model = rep(c,nrow(cluster_celltype_df_reshape)),
                                                    cluster_celltype_df_reshape))
  })
  cluster_celltype_df_reshape_out <- rbindlist(cluster_celltype_df_reshape_list,fill = T)

  return(cluster_celltype_df_reshape_out)

}

