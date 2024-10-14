



#' Detect human or mouse data
#'
#' @param genes Character vector of gene symbols.
#'
#' @return String stating the species detected.
#' @export
#'
#' @examples
#' x <- c("MYC","BRAF","TP53")
#' y <- c("Myc","Braf","Tp53")
#' detect_species(x)
#' detect_species(y)


detect_species <- function(genes) {
  # check the capitalization pattern
  is_human_gene <- function(gene) {
    return(grepl("^[A-Z]+$", gene))
  }
  is_mouse_gene <- function(gene) {
    return(grepl("^[A-Z][a-z]*$", gene))
  }
  # Count the number of human and mouse gene patterns
  human_count <- sum(sapply(genes, is_human_gene))
  mouse_count <- sum(sapply(genes, is_mouse_gene))
  # Determine the majority match
  if (human_count > mouse_count) {
    return("human")
  } else if (mouse_count > human_count) {
    return("mouse")
  } else {
    return("undetermined") # If counts are equal or if there are no matches
  }
}

#' Convert excel dates to gene symbols
#'
#' @param vec Character vector of gene symbols.
#' @param mouse Boolean. If data is from mouse, use TRUE, else FALSE.
#'
#' @return Character vector of gene symbols.
#' @export
#'
#' @examples
#' x <- c("1-Dec","MYC","TP53")
#' y <- c("1-Dec","Myc","Tp53")
#' date_to_gene(x)
#' date_to_gene(y, mouse = TRUE)


date_to_gene <- function(vec,mouse = FALSE) {
  if (is.null(vec)) exit("Please provide vector argument")
  if (!mouse) {
    vec <- gsub('(\\d+)-(Mar)','MARCH\\1',vec)
    vec <- gsub('(\\d+)-(Sep)','SEPT\\1',vec)
    vec <- gsub('(\\d+)-(Dec)','DEC\\1',vec)
  } else {
    vec <- gsub('(\\d+)-(Mar)','March\\1',vec)
    vec <- gsub('(\\d+)-(Sep)','Sept\\1',vec)
    vec <- gsub('(\\d+)-(Dec)','Dec\\1',vec)
  }
  return(vec)
}

#' Filter Seurat object by percent mitochondria or number of features
#'
#' @param seurat Seurat object.
#' @param species String denoting if the data is from "human" or "mouse".
#' @param percent_mt Numeric. Percent mitochondria to filter out. Default is NULL.
#' @param minFeatures Numeric. Minimum number of features a barcode/cell must have. Default to 500.
#'
#' @return Seurat object of filtered features.
#' @export
#'


berlin_filter <- function(seurat = NULL, species = "human", percent_mt = NULL, minFeatures = 500) {
  species <- tolower(species)
  mt.pattern <- case_when(
    species == "human" ~ "^MT-",
    species == "mouse" ~ "^mt-",
    TRUE ~ "^MT-"
  )
  ribo.pattern <- case_when(
    species == "human" ~ "^RP[LS]",
    species == "mouse" ~ "^Rp[ls]",
    TRUE ~ "^RP[LS]"
  )

  # Calculate percentage of mitochondrial and ribosomal genes
  seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = mt.pattern, assay = "RNA")
  seurat[["percent.rp"]] <- PercentageFeatureSet(seurat, pattern = ribo.pattern, assay = "RNA")

  # Filter cells based on QC criteria
  if (!is.null(percent_mt)) {
    seurat[,seurat[["percent.mt"]] <= percent_mt & seurat[["nFeature_RNA"]] >= minFeatures]
  } else {
    seurat[,seurat[["nFeature_RNA"]] >= minFeatures]
    }
}


#' Perform quality control and initial formatting of Seurat object
#'
#' @param object Seurat object. If NULL, must supply count data.
#' @param counts A matrix, data.frame, or dgCMatrix. The rows are unique features. The feature name can be the rowname or the first column. Each column is a uniquely labeled barcode/cell.
#' @param meta A data.frame. The rows are barcode/cell names, and the column are any additionally supplied metadata. Row names in the metadata need to match the column names of the counts matrix.
#' @param assay String. Name of the initial/default assay.
#' @param project_name String. Desired project name.
#' @param species String. Denoting if the data is from "human" or "mouse".
#' @param percent_mt Numeric. Percent mitochondria to filter out. Default is NULL.
#' @param minFeatures Numeric. Minimum number of features a barcode/cell must have. Default to 500.
#' @param varFeatures Numeric. Number of features to select as top variable features. Default to 2000.
#' @param verbose Boolean. To show progress, TRUE, else FALSE.
#' @param logScale Numeric. Set the scale to log normalize data. Defualt to 10000.
#' @param remove_duplicates Boolean. To remove duplicate features that may be found in the data, use TRUE, else FALSE and the function will notify you if duplicate features are found.
#'
#' @return Seurat object.
#' @export
#'
berlin_qc <- function(object = NULL, counts = NULL, meta = NULL, assay = "RNA", project_name = "BERLIN_Project", species = "human",
                      percent_mt = NULL, minFeatures = 500, varFeatures = 2000, verbose = TRUE, logScale = 10000, remove_duplicates = TRUE) {

  if (is.null(object) & is.null(counts)) stop("Please supply Seurat object or counts matrix")

  if (!is.null(counts) & is.null(object)) {
    if (class(counts)[1] != "dgCMatrix") {
      # If first column is gene symbols
      if (is.character(counts[,1])) {
        # Check for duplicate genes
        counts_dup <- counts[which(counts[,1] %in% counts[,1][duplicated(counts[,1])]),]
        if (nrow(counts_dup) > 0) {
          if (!remove_duplicates) {
            stop("Duplicate features found. Remove duplicates or set 'remove_duplicates' argument to TRUE.")
          }
        }
        counts_nondup <- counts[which(!counts[,1] %in% counts[,1][duplicated(counts[,1])]),]
        if (nrow(counts_dup) > 0) {
          counts_dup <- counts_dup %>%
            group_by(!!sym(colnames(counts)[1])) %>%
            summarise_all(max) %>%
            as.data.frame()
        }
        counts <- rbind(counts_dup,counts_nondup)
        rownames(counts) <- counts[,1]
        counts <- as.matrix(counts[,-1])
      }
    }
    # detect species
    species_detected <- detect_species(rownames(counts))
    # if different than input species, notify user
    if (species_detected != species) {
      print(paste0("Species detected does not equeale the species argument. Will be treating data as ",species_detected," data."))
    }
    # Check that gene symbols are not in excel date format
    rownames(counts) <- date_to_gene(rownames(counts), ifelse(species_detected == "mouse",TRUE,FALSE))
    object <- CreateSeuratObject(counts = counts, assay = assay, project = project_name, meta.data = meta)
  }
  if (!is.null(object)) {
    if (!(class(object) == "Seurat")) stop("Input object must be Seurat object.")
    # detect species
    species_detected <- detect_species(Features(object))
    # if different than input species, notify user
    if (species_detected != species) {
      print(paste0("Species detected does not equeale the species argument. Will be treating data as ",species_detected," data."))
    }
  }

  if (!assay %in% names(object)) stop("Assay input is not found in object")
  DefaultAssay(object) <- assay
  object <- berlin_filter(object, species_detected, percent_mt, minFeatures)
  object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = logScale, verbose = verbose)
  object <- CellCycleScoring(object = object, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes, verbose = verbose)
  object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = varFeatures, verbose = verbose)
  object <- ScaleData(object, vars.to.regress = c("nFeature_RNA", "percent.mt"), verbose = verbose)

  if (!is.null(meta)) {
    object <- AddMetaData(object,metadata = meta)
  }

  return(object)
}
