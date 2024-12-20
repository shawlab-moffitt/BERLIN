


#' Run UMAP Shiny App
#'
#' @param launch.browser Boolean logic to launch application in browser. Default to TRUE.
#' @return An R Shiny Application.
#' @export
#' @rawNamespace import(shiny, except=c(dataTableOutput, renderDataTable))
#' @import shinycssloaders
#' @import shinythemes
#' @importFrom shinyjs useShinyjs
#' @import shinyjqui
#' @import pheatmap
#' @import RColorBrewer
#' @import umap
#' @import cowplot
#' @import patchwork
#' @import ggdendro
#' @import factoextra
#' @import dplyr
#' @import viridis
#' @import readr
#' @import ggVennDiagram
#' @import ggtree
#' @import stringr
#' @import tools
#' @import plotly
#' @import reshape2
#' @import gridExtra
#' @import tidyverse
#' @import DT
#' @import ggplot2
#' @import ggpubr
#' @import ggrepel
#' @import scales
#' @import data.table
#'
runUMAPapp <- function(object = NULL, counts = NULL, meta = NULL, n_cells = 2000, assay = "RNA", save_data = TRUE, project_name = "BERLIN_Project",
                       umap1_col = "UMAP_1", umap2_col = "UMAP_2", anno1_col = "seurat_clusters", anno2_col = NULL, anno3_col = NULL,
                       seed = 42, species = "human", remove_duplicates = TRUE, launch.browser = TRUE, species_detected = "human", verbose = TRUE) {

  set.seed(42)
  if (is.null(object) & is.null(counts)) stop("Please supply Seurat object or counts matrix")

  if (!is.null(counts) & is.null(object)) {
    if (verbose) {
      message("Reading count data")
    }
      # If first column is gene symbols
      if (is.character(counts[,1])) {
        # Check for duplicate genes
        counts_dup <- counts[which(counts[,1] %in% counts[,1][duplicated(counts[,1])]),]
        if (nrow(counts_dup) > 0) {
          if (verbose) {
            message("Reducing duplicate features")
          }
          if (!remove_duplicates) {
            stop("Duplicate features found. Remove duplicates or set 'remove_duplicates' argument to TRUE.")
          }
        }
        counts_nondup <- counts[which(!counts[,1] %in% counts[,1][duplicated(counts[,1])]),]
        if (nrow(counts_dup) > 0) {
          counts_dup <- counts_dup %>%
            dplyr::group_by(!!sym(colnames(counts)[1])) %>%
            dplyr::summarise_all(max) %>%
            as.data.frame()
        }
        counts <- rbind(counts_dup,counts_nondup)
      }
    # detect species
    species_detected <- detect_species(counts[,1])
    # if different than input species, notify user
    if (species_detected != species) {
      if (verbose) {
        message(paste0("Species detected does not equeal the species argument. Will be treating data as ",species_detected," data."))
      }
    }
    # Check that gene symbols are not in excel date format
    counts[,1] <- date_to_gene(counts[,1], ifelse(species_detected == "mouse",TRUE,FALSE))
  }
  if (is.null(meta) & is.null(object)) {
    stop("No meta data provided")
  }
  if (!is.null(object)) {
    if (verbose) {
      message("Extracting count data from Seurat object")
    }
    counts <- as.data.frame(object[[assay]]$data)
    counts <- cbind(gene = rownames(counts),counts)
    meta <- object[[]]
    if (nrow(meta) == 0) {
      stop("No meta data found in Seurat object")
    }
    meta <- cbind(data.frame(SampleName = rownames(meta)),
                  meta)
  }
  sampsame <- intersect(meta[,1], colnames(counts)[-1])
  if (length(sampsame) == 0) {
    stop("Sample name column missing from meta data\nShould be first column")
  }
  meta <- meta[which(meta[,1] %in% sampsame),]
  counts <- counts[,which(colnames(counts) %in% c("gene",sampsame))]

  if ((ncol(counts)-1) > n_cells) {
    if (verbose) {
      message(paste0("Subsetting ",n_cells," random cells from data for app."))
    }
    random_cells <- sample(2:ncol(counts), size = n_cells, replace = FALSE)
    counts <- counts[,c("gene",random_cells)]

    meta <- meta[which(meta[,1] %in% random_cells),]
  }

  if (save_data) {
    if (verbose) {
      message(paste0("Saving app input data to: ",getwd()))
    }
    data.table::fwrite(counts,paste0(project_name,"_NormCounts_",Sys.Date(),".txt"), sep = '\t', na = NA)
    data.table::fwrite(meta,paste0(project_name,"_MetaData_",Sys.Date(),".txt"), sep = '\t', na = NA)
  }

  if (verbose) {
    message("Starting up application.")
  }
  shiny::shinyOptions(project_name = project_name, counts = counts, meta = meta, species_detected = species_detected,
               umap1_col = umap1_col, umap2_col = umap2_col,
               anno1_col = anno1_col, anno2_col = anno2_col, anno3_col = anno3_col)


  appDir <- system.file("UMAP_App", package = "BERLIN")
  if (appDir == "") {
    stop("Could not find UMAP_App Try re-installing `BERLIN`.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal", launch.browser = launch.browser)
}


