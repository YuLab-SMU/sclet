#' @include state.R
NULL

#' Run DoubletFinder to identify and filter doublets
#'
#' @title RunDoubletFinder
#' @param sce a SingleCellExperiment object
#' @param ... additional parameters passed to \code{scDblFinder::scDblFinder}
#' @return an updated SingleCellExperiment object with doublet scores and classes
#' @export
#' @importFrom SummarizedExperiment colData colData<- assayNames
RunDoubletFinder <- function(sce, ...) {
    if (!requireNamespace("scDblFinder", quietly = TRUE)) {
        stop("Please install the 'scDblFinder' package to use RunDoubletFinder.")
    }
    
    # Check if 'counts' assay exists
    if (!"counts" %in% SummarizedExperiment::assayNames(sce)) {
        stop("The 'counts' assay is required for scDblFinder.")
    }
    
    prev_state <- sclet_get_state(sce)
    
    # Run scDblFinder
    # scDblFinder::scDblFinder returns a SCE with several colData columns:
    # scDblFinder.class, scDblFinder.score, scDblFinder.weighted, etc.
    res <- scDblFinder::scDblFinder(sce, ...)
    
    # Transfer the required columns to the original sce object
    cd <- SummarizedExperiment::colData(res)
    SummarizedExperiment::colData(sce)$scDblFinder.score <- cd$scDblFinder.score
    SummarizedExperiment::colData(sce)$scDblFinder.class <- cd$scDblFinder.class
    
    sce <- sclet_restore_state(sce, prev_state)
    
    # Log the command and update analysis state
    sce <- sclet_log_command(
        sce,
        "RunDoubletFinder",
        params = list(...)
    )
    
    doublet_count <- sum(cd$scDblFinder.class == "doublet", na.rm = TRUE)
    singlet_count <- sum(cd$scDblFinder.class == "singlet", na.rm = TRUE)
    
    sce <- sclet_set_analysis_state(
        object = sce,
        type = "preprocess",
        id = "doublet_finder",
        method = "scDblFinder",
        inputs = list(assay = "counts"),
        artifacts = list(
            colData = c("scDblFinder.score", "scDblFinder.class")
        ),
        params = list(...),
        summary = list(
            doublet_count = doublet_count,
            singlet_count = singlet_count
        )
    )
    
    return(sce)
}

#' Run decontX to remove ambient RNA contamination
#'
#' @title RunDecontX
#' @param sce a SingleCellExperiment object
#' @param assay the name of the assay to use as input. Defaults to "counts".
#' @param ... additional parameters passed to \code{celda::decontX}
#' @return an updated SingleCellExperiment object with 'decontXcounts' assay
#' @export
#' @importFrom SummarizedExperiment assay assay<- colData colData<- assayNames
RunDecontX <- function(sce, assay = "counts", ...) {
    if (!requireNamespace("celda", quietly = TRUE)) {
        stop("Please install the 'celda' package to use RunDecontX.")
    }
    
    if (!assay %in% SummarizedExperiment::assayNames(sce)) {
        stop("Assay '", assay, "' not found in assays(sce).")
    }
    
    prev_state <- sclet_get_state(sce)
    
    counts_mat <- SummarizedExperiment::assay(sce, assay)
    
    # Run decontX
    # celda::decontX returns a list containing 'decontXcounts', 'contamination', etc.
    res <- celda::decontX(counts_mat, ...)
    
    # Store the decontaminated matrix into a new assay
    SummarizedExperiment::assay(sce, "decontXcounts") <- res$decontXcounts
    
    # Optional: store contamination scores in colData
    SummarizedExperiment::colData(sce)$decontX_contamination <- res$contamination
    
    sce <- sclet_restore_state(sce, prev_state)
    
    # Update layers
    sce <- sclet_set_layer(
        sce,
        name = "decontXcounts",
        assay = "decontXcounts",
        role = "cleaned",
        params = list(source = assay, ...)
    )
    
    # Log the command and update analysis state
    sce <- sclet_log_command(
        sce,
        "RunDecontX",
        params = list(assay = assay, ...)
    )
    
    sce <- sclet_set_analysis_state(
        object = sce,
        type = "preprocess",
        id = "decontX",
        method = "celda::decontX",
        inputs = list(assay = assay),
        artifacts = list(
            assay = "decontXcounts",
            colData = "decontX_contamination"
        ),
        params = list(assay = assay, ...),
        summary = list(
            mean_contamination = mean(res$contamination, na.rm = TRUE)
        )
    )
    
    return(sce)
}
