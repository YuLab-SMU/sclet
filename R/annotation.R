#' Run SingleR for cell type annotation
#' 
#' Wrapper for SingleR to perform automatic cell type annotation.
#' 
#' @title RunSingleR
#' @param object A SingleCellExperiment object
#' @param ref Reference dataset (e.g. from celldex). If NULL, will try to use HumanPrimaryCellAtlasData.
#' @param labels Labels in the reference dataset. If NULL, tries to use 'label.main'.
#' @param assay.type Assay to use for annotation. Default is "logcounts".
#' @param ... Additional arguments passed to SingleR::SingleR
#' @return A SingleCellExperiment object with annotation added to colData
#' @export
RunSingleR <- function(object, ref = NULL, labels = NULL, assay.type = "logcounts", ...) {
    if (!requireNamespace("SingleR", quietly = TRUE)) {
        stop("Package 'SingleR' is needed for this function to work. Please install it.")
    }
    
    if (is.null(ref)) {
        if (!requireNamespace("celldex", quietly = TRUE)) {
            stop("Package 'celldex' is needed to load default reference. Please install it.")
        }
        message("Loading HumanPrimaryCellAtlasData from celldex...")
        ref <- celldex::HumanPrimaryCellAtlasData()
    }
    
    if (is.null(labels)) {
        if ("label.main" %in% colnames(SummarizedExperiment::colData(ref))) {
            labels <- ref$label.main
        } else {
            stop("Please specify 'labels' for the reference dataset.")
        }
    }
    
    pred <- SingleR::SingleR(test = object, ref = ref, labels = labels, assay.type.test = assay.type, ...)
    
    # Add predictions to colData
    SummarizedExperiment::colData(object)$SingleR_labels <- pred$labels
    SummarizedExperiment::colData(object)$SingleR_pruned.labels <- pred$pruned.labels
    
    message("Annotation added to colData columns: 'SingleR_labels' and 'SingleR_pruned.labels'")
    
    return(object)
}
