#' @include state.R
NULL

#' Run imputation to denoise single-cell data
#'
#' @title RunImputation
#' @param sce a SingleCellExperiment object
#' @param method the imputation method to use. Currently only supports "ALRA" (zero-preserving, SVD-based imputation).
#' @param assay the assay to use as input. Defaults to "logcounts".
#' @param ... additional parameters passed to the underlying imputation method (e.g., \code{ALRA::alra}).
#' @return an updated SingleCellExperiment object with imputed data in a new assay (e.g., 'alra_imputed')
#' @export
#' @importFrom SummarizedExperiment assay assay<- assayNames
RunImputation <- function(sce, method = "ALRA", assay = "logcounts", ...) {
    method <- match.arg(toupper(method), choices = c("ALRA"))
    
    if (!assay %in% SummarizedExperiment::assayNames(sce)) {
        stop("Assay '", assay, "' not found in assays(sce).")
    }
    
    prev_state <- sclet_get_state(sce)
    
    if (method == "ALRA") {
        if (!requireNamespace("ALRA", quietly = TRUE)) {
            stop("Please install the 'ALRA' package via remotes::install_github('KlugerLab/ALRA')")
        }
        
        mat <- SummarizedExperiment::assay(sce, assay)
        
        # ALRA expects cells as rows and genes as columns.
        # Convert to a plain matrix before calling ALRA, because ALRA's
        # internal matrix operations can drop dimensions on sparse/delayed inputs.
        mat_t <- as.matrix(t(mat))
        
        # We wrap this in a message because ALRA can be verbose and take time
        message("Running ALRA imputation...")
        
        # ALRA can sometimes throw errors if the matrix is purely sparse, 
        # but rsvd (which ALRA uses) usually handles dgCMatrix well.
        # To be safe and compatible with standard ALRA pipeline:
        res <- ALRA::alra(mat_t, ...)
        
        # The completed matrix is A_norm_rank_k_cor_sc
        imputed_mat <- t(res$A_norm_rank_k_cor_sc)
        rownames(imputed_mat) <- rownames(sce)
        colnames(imputed_mat) <- colnames(sce)
        
        target_assay <- "alra_imputed"
        SummarizedExperiment::assay(sce, target_assay) <- imputed_mat
    }
    
    sce <- sclet_restore_state(sce, prev_state)
    
    # Update layers
    sce <- sclet_set_layer(
        sce,
        name = target_assay,
        assay = target_assay,
        role = "imputed",
        params = list(source = assay, method = method, ...)
    )
    
    # Log command and update analysis state
    sce <- sclet_log_command(
        sce,
        "RunImputation",
        params = list(method = method, assay = assay, ...)
    )
    
    sce <- sclet_set_analysis_state(
        object = sce,
        type = "preprocess",
        id = "imputation",
        method = method,
        inputs = list(assay = assay),
        artifacts = list(assay = target_assay),
        params = list(method = method, assay = assay, ...),
        summary = list()
    )
    
    return(sce)
}
