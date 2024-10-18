#' Perform SingleR annotation of Seurat object
#'
#' @param object Seurat object. If NULL, must supply count data.
#' @param verbose Boolean. To show progress, TRUE, else FALSE.
#'
#' @return A Seurat object.
#' @export
#'
berlin_annotate <- function(object = NULL, verbose = TRUE) {

  if (is.null(object)) stop("Please supply Seurat object as input.")

  # Load reference databases from celldex
  if (verbose) {
    message("Loading reference databases from celldex")
  }
  hpca.ref <- celldex::HumanPrimaryCellAtlasData()
  dice.ref <- celldex::DatabaseImmuneCellExpressionData()
  blueprint.ref <- celldex::BlueprintEncodeData()
  monaco.ref <- celldex::MonacoImmuneData()
  northern.ref <- celldex::NovershternHematopoieticData()

  # Convert Seurat object to SingleCellExperiment
  if (verbose) {
    message("Converting Seurat object to SingleCellExperiment")
  }
  sce <- Seurat::as.SingleCellExperiment(Seurat::DietSeurat(object))

  # Auto-annotating cell types using reference databases from celldex
  if (verbose) {
    message("Auto-annotating cell types using reference databases from celldex")
    pb = txtProgressBar(min = 0, max = 10, initial = 0, style = 3)
    hpca.main <- SingleR::SingleR(test = sce, assay.type.test = 1, ref = hpca.ref, labels = hpca.ref$label.main)
    setTxtProgressBar(pb,1)
    hpca.fine <- SingleR::SingleR(test = sce, assay.type.test = 1, ref = hpca.ref, labels = hpca.ref$label.fine)
    setTxtProgressBar(pb,2)
    dice.main <- SingleR::SingleR(test = sce, assay.type.test = 1, ref = dice.ref, labels = dice.ref$label.main)
    setTxtProgressBar(pb,3)
    dice.fine <- SingleR::SingleR(test = sce, assay.type.test = 1, ref = dice.ref, labels = dice.ref$label.fine)
    setTxtProgressBar(pb,4)
    blue.main <- SingleR::SingleR(test = sce, assay.type.test = 1, ref = blueprint.ref, labels = blueprint.ref$label.main)
    setTxtProgressBar(pb,5)
    blue.fine <- SingleR::SingleR(test = sce, assay.type.test = 1, ref = blueprint.ref, labels = blueprint.ref$label.fine)
    setTxtProgressBar(pb,6)
    monaco.main <- SingleR::SingleR(test = sce, assay.type.test = 1, ref = monaco.ref, labels = monaco.ref$label.main)
    setTxtProgressBar(pb,7)
    monaco.fine <- SingleR::SingleR(test = sce, assay.type.test = 1, ref = monaco.ref, labels = monaco.ref$label.fine)
    setTxtProgressBar(pb,8)
    northern.main <- SingleR::SingleR(test = sce, assay.type.test = 1, ref = northern.ref, labels = northern.ref$label.main)
    setTxtProgressBar(pb,9)
    northern.fine <- SingleR::SingleR(test = sce, assay.type.test = 1, ref = northern.ref, labels = northern.ref$label.fine)
    setTxtProgressBar(pb,10)
    close(pb)
  } else {
    hpca.main <- SingleR::SingleR(test = sce, assay.type.test = 1, ref = hpca.ref, labels = hpca.ref$label.main)
    hpca.fine <- SingleR::SingleR(test = sce, assay.type.test = 1, ref = hpca.ref, labels = hpca.ref$label.fine)
    dice.main <- SingleR::SingleR(test = sce, assay.type.test = 1, ref = dice.ref, labels = dice.ref$label.main)
    dice.fine <- SingleR::SingleR(test = sce, assay.type.test = 1, ref = dice.ref, labels = dice.ref$label.fine)
    blue.main <- SingleR::SingleR(test = sce, assay.type.test = 1, ref = blueprint.ref, labels = blueprint.ref$label.main)
    blue.fine <- SingleR::SingleR(test = sce, assay.type.test = 1, ref = blueprint.ref, labels = blueprint.ref$label.fine)
    monaco.main <- SingleR::SingleR(test = sce, assay.type.test = 1, ref = monaco.ref, labels = monaco.ref$label.main)
    monaco.fine <- SingleR::SingleR(test = sce, assay.type.test = 1, ref = monaco.ref, labels = monaco.ref$label.fine)
    northern.main <- SingleR::SingleR(test = sce, assay.type.test = 1, ref = northern.ref, labels = northern.ref$label.main)
    northern.fine <- SingleR::SingleR(test = sce, assay.type.test = 1, ref = northern.ref, labels = northern.ref$label.fine)
  }

  # Add the celldex annotations to the metadata of the Seurat object
  object@meta.data$hpca.main <- hpca.main$pruned.labels
  object@meta.data$hpca.fine <- hpca.fine$pruned.labels
  object@meta.data$dice.main <- dice.main$pruned.labels
  object@meta.data$dice.fine <- dice.fine$pruned.labels
  object@meta.data$monaco.main <- monaco.main$pruned.labels
  object@meta.data$monaco.fine <- monaco.fine$pruned.labels
  object@meta.data$northern.main <- northern.main$pruned.labels
  object@meta.data$northern.fine <- northern.fine$pruned.labels
  object@meta.data$blue.main <- blue.main$pruned.labels
  object@meta.data$blue.fine <- blue.fine$pruned.labels

  return(object)
}
