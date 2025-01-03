

#' Write out BERLIN Annotation Summary
#'
#' @param object Seurat object. If NULL, must supply meta data.
#' @param meta A data.frame. The rows are barcode/cell names, and the column are any additionally supplied metadata. Row names in the metadata need to match the column names of the counts matrix.
#' @param file String. Desired file name with directory for output.
#' @param excel Boolean. If TRUE, file will be save as excel workbook file.
#' @param anno_cols Character Vector. Vector of column names to be subset and written out to file. List of default columns in details.
#' @param return_seurat Boolean. If TRUE, a Seurat object will be the output with the results within the stored in the Misc slot.
#'
#' @return data frame
#' @export
#'
#' @details
#' Default annotation columns: c("scType_seurat_clusters_Classification","scType Score Sum","ncells","scType_Cell_Classification",
#' "scGate_Myeloid_Macrohphage_Celltype","seurat_clusters","RNA_snn_res.0.5","RNA_snn_res.1",
#' "RNA_snn_res.1.5","RNA_snn_res.2","dice.main","dice.fine", "monaco.main","monaco.fine",
#' "northern.main","northern.fine","blue.main","blue.fine")
#'
#'




berlin_anno_summary <- function(object = NULL, meta = NULL, file = NULL, excel = TRUE, return_seurat = FALSE,
                                anno_cols = c("^scType_.*\\_Classification$","scType Score Sum","ncells","scType_Cell_Classification",
                                              "^scGate_.*\\_Celltype$","seurat_clusters","RNA_snn_res.0.5","RNA_snn_res.1",
                                              "RNA_snn_res.1.5","RNA_snn_res.2","dice.main","dice.fine", "monaco.main","monaco.fine",
                                              "northern.main","northern.fine","blue.main","blue.fine")) {

  call.string <- deparse(expr = sys.calls()[[1]])
  func_name <- "berlin_anno_summary"
  time.stamp <- Sys.time()
  argg <- c(as.list(environment()))
  argg <- Filter(function(x) any(is.numeric(x) | is.character(x)), argg)

  if (is.null(object) & is.null(meta)) stop("Please supply Seurat object or meta data")

  if (is.null(meta)) {
    meta <- object[[]]
  }

  cols_found <- grep(paste(anno_cols,collapse = "|"), colnames(meta), value = T, ignore.case = T)

  meta_anno <- meta[,cols_found]
  scGate_col <- grep("^scGate_.*\\_Celltype$", colnames(meta_anno), value = T)
  scGate_col <- grep("multiscGate",scGate_col, invert = T, value = T)
  scType_col <- grep("^scType_.*\\_Classification$", colnames(meta_anno), value = T)
  scType_col <- grep("scType_Cell_Classification",scType_col,invert = T, value = T)

  meta_anno <- meta_anno %>%
    relocate(any_of(c(scType_col,"scType Score Sum","ncells","scType_Cell_Classification",scGate_col,"seurat_clusters")))
  meta_anno <- cbind(cell = rownames(meta_anno),
                     meta_anno)


  if (excel) {
    if (is.null(file)) {
      file <- paste0(getwd(),"/BERLIN_Annotation_Summary_",Sys.Date(),".xlsx")
    }
    writexl::write_xlsx(meta_anno, file)
  } else {
    if (is.null(file)) {
      file <- paste0(getwd(),"/BERLIN_Annotation_Summary_",Sys.Date(),".txt")
    }
    write.table(meta_anno,file, sep = '\t', row.names = F)
  }

  if (!return_seurat) {
    return(meta_anno)
  } else {
    SeuratObject::Misc(object = object, slot = "BERLIN_Annotation_Summary") <- meta_anno
    slot(object = object, name = "commands")[[func_name]] <- list(name = func_name,
                                                                  time.stamp = time.stamp,
                                                                  call.string = call.string,
                                                                  params = argg)
    return(object)
  }

}
