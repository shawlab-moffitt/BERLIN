

#' Performs ScGate function with a single model on a Seurat object
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


berlin_run_scgate <- function(object = NULL, model = NULL, model_name = NULL, assay = "RNA", ncores = 1, seed = 42, verbose = TRUE) {
  if (is.null(object)) stop("Please supply Seurat object")
  if (is.null(model)) stop("Please supply ScGate model")
  if (is.data.frame(model)) {
    if (colnames(model) != c("levels","use_as","name","signature")) stop("Please check model format.\nMust follow column names: levels, use_as, name, signature.")
  }
  if (is.data.frame(model) | is.data.frame(model[[1]])) {
    if (is.null(model_name)) stop("Please provide model name.")
  }

  original_meta_col_n <- ncol(object[[]])

  object_scgate <- scGate::scGate(object, model, assay = assay, ncores = ncores, seed = seed, verbose = verbose)

  meta <- object_scgate[[]]
  new_meta_cols_n <- ncol(meta)
  colnames(meta)[c((original_meta_col_n+1):ncol(meta))] <- paste0(colnames(meta)[c((original_meta_col_n+1):ncol(meta))],"_ScGate_", model_name, "_CellType")
  pure_columns <- grep("^is.pure", colnames(meta), value = TRUE)
  new_col_name <- paste0("scGate_", model_name, "_Celltype")
  if (length(pure_columns) == 1) {
    meta[[new_col_name]] <- ifelse(meta[,pure_columns] == "Pure",model_name,NA)
  } else {
    pure_temp <- sapply(pure_columns, function(c) {
      feature <- sub("^is.pure_", "", c)
      feature <- sub(paste0("_ScGate_", model_name, "_CellType"), "", feature)
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


#' Performs ScGate function on a Seurat object with one or multiple models
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


berlin_scgate <- function(object = NULL, assay = "RNA", model = NULL, model_name = NULL, species = "human", ncores = 1, seed = 42, verbose = TRUE) {

  if (is.null(object)) stop("Please supply Seurat object")

  if (is.null(model)) {
    # detect species
    species_detected <- detect_species(Features(object))
    # if different than input species, notify user
    if (species_detected != species) {
      message(paste0("Species detected does not equeal the species argument. Will be treating data as ",species_detected," data."))
    }
    scGate_models_DB <-  get_scGateDB()
    models <- scGate_models_DB[[species_detected]]
  } else {
    if (is.data.frame(model)) {
      if (colnames(model) != c("levels","use_as","name","signature")) stop("Please check model format.\nMust follow column names: levels, use_as, name, signature.")
    }
    if (is.data.frame(model) | is.data.frame(model[[1]])) {
      if (is.null(model_name)) stop("Please provide model name.")
    }
  }

  # Create a list to store new meta data
  if (is.data.frame(model) | is.data.frame(model[[1]])) {
    new_meta <- berlin_run_scgate(object = object, model = model, model_name = model_name, assay = assay, ncores = ncores, seed = seed, verbose = verbose)
  } else {
    model_names <- names(model)
    new_metas <- sapply(model_names[1:3], function(m) {
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
