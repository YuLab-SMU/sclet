#' Aggregate expression for pseudobulk analysis
#' 
#' Aggregates expression values by groups (e.g. cluster and sample) to create a pseudobulk object.
#' 
#' @title AggregateExpression
#' @param object A SingleCellExperiment object
#' @param group_by A character vector of column names in colData(object) to group by. 
#'        Typically c("cluster", "sample").
#' @param assay.type Assay to aggregate. Default is "counts".
#' @param fun Aggregation function. Default is "sum".
#' @return A SummarizedExperiment object with aggregated counts. 
#'         The colData contains the unique combinations of the grouping variables.
#' @importFrom SummarizedExperiment SummarizedExperiment assay
#' @importFrom Matrix rowSums
#' @importFrom methods as
#' @importFrom scuttle aggregateAcrossCells
#' @export
AggregateExpression <- function(object, group_by, assay.type = "counts", fun = "sum") {
    
    if (!all(group_by %in% colnames(colData(object)))) {
        stop("Not all grouping variables are present in colData.")
    }
    
    # Use scrapper::aggregateAcrossCells.se if available, otherwise fallback to scuttle
    # It returns a SummarizedExperiment
    ids <- colData(object)[, group_by, drop=FALSE]
    
    # aggregateAcrossCells handles "sum" (counts) naturally
    # If other stats needed, we might need other approach, but for pseudobulk DE, sum is standard.
    if (fun != "sum") {
        warning("Currently only 'sum' is supported optimally. Using 'sum'.")
    }
    
    if (requireNamespace("scrapper", quietly = TRUE) && identical(assay.type, "counts") && identical(fun, "sum")) {
        ids_list <- stats::setNames(lapply(group_by, function(x) colData(object)[[x]]), group_by)
        agg_sce <- tryCatch(
            scrapper::aggregateAcrossCells.se(object, ids_list),
            error = function(e) NULL
        )
    } else {
        agg_sce <- NULL
    }

    if (is.null(agg_sce)) {
        agg_sce <- sclet_muffle_known_warnings(
            scuttle::aggregateAcrossCells(object, ids = ids, use.assay.type = assay.type, statistics = "sum"),
            patterns = c("'aggregateAcrossCells' is deprecated")
        )
    }
    
    return(agg_sce)
}


#' Find markers using pseudobulk DE (DESeq2)
#' 
#' Performs differential expression analysis on pseudobulk data using DESeq2.
#' 
#' @title FindMarkers_pseudobulk
#' @param object A SingleCellExperiment or SummarizedExperiment object (output of AggregateExpression)
#' @param design A formula for the design matrix (e.g. ~ condition).
#' @param contrast A character vector specifying the contrast (e.g. c("condition", "treated", "control")).
#' @param ... Additional arguments passed to DESeq2::DESeq
#' @return A data.frame with DE results
#' @export
FindMarkers_pseudobulk <- function(object, design, contrast, ...) {
    if (!requireNamespace("DESeq2", quietly = TRUE)) {
        stop("Package 'DESeq2' is needed for pseudobulk DE. Please install it.")
    }
    
    dds <- tryCatch(
        DESeq2::DESeqDataSet(object, design = design),
        error = function(e) NULL
    )
    if (is.null(dds)) {
        counts_mat <- assay(object, "counts")
        if (!is.matrix(counts_mat)) counts_mat <- as.matrix(counts_mat)
        mode(counts_mat) <- "integer"
        dds <- DESeq2::DESeqDataSetFromMatrix(
            countData = counts_mat,
            colData = colData(object),
            design = design
        )
    }
    
    # Run DESeq
    dds <- DESeq2::DESeq(dds, ...)
    
    # Get results
    res <- DESeq2::results(dds, contrast = contrast)
    
    return(as.data.frame(res))
}
