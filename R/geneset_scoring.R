#' Run Gene Set Scoring
#'
#' This function calculates gene set or pathway activity scores for individual cells using 
#' methods like AUCell, UCell, or GSVA. The resulting scores are added to the `colData` of the 
#' SingleCellExperiment object and registered in the analysis-state contract.
#'
#' @param sce A SingleCellExperiment object.
#' @param gene_sets A named list of character vectors containing gene symbols.
#' @param method Character. The scoring method to use. Options are "AUCell", "UCell", or "GSVA". Defaults to "UCell".
#' @param assay_use Character. The assay to use for scoring. Defaults to "counts" (often preferred for AUCell/UCell) or "logcounts".
#' @param ncores Integer. Number of cores for parallel processing. Defaults to 1.
#' @param name Character. Name of the scoring run, used as prefix for colData columns. Defaults to "Score".
#' @param ... Additional arguments passed to the underlying scoring function.
#' 
#' @return A SingleCellExperiment object with gene set scores added to `colData(sce)` and state registered.
#' 
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom SummarizedExperiment assay colData colData<-
#' @export
RunGeneSetScoring <- function(sce, gene_sets, method = c("UCell", "AUCell", "GSVA"), 
                              assay_use = "counts", ncores = 1, name = "Score", ...) {
    
    method <- match.arg(method)
    
    if (!is.list(gene_sets) || is.null(names(gene_sets))) {
        stop("gene_sets must be a named list of character vectors.")
    }
    
    if (!assay_use %in% assayNames(sce)) {
        stop(sprintf("Assay '%s' not found in the SingleCellExperiment object.", assay_use))
    }
    
    expr_mat <- assay(sce, assay_use)
    
    scores <- NULL
    
    if (method == "UCell") {
        if (!requireNamespace("UCell", quietly = TRUE)) {
            stop("Please install 'UCell' package to use method='UCell'.")
        }
        message("Running UCell scoring...")
        scores <- UCell::ScoreSignatures_UCell(expr_mat, features = gene_sets, ncores = ncores, name = "", ...)
        
    } else if (method == "AUCell") {
        if (!requireNamespace("AUCell", quietly = TRUE)) {
            stop("Please install 'AUCell' package to use method='AUCell'.")
        }
        message("Running AUCell scoring...")
        cells_rankings <- AUCell::AUCell_buildRankings(expr_mat, nCores = ncores, plotStats = FALSE, verbose = FALSE)
        cells_AUC <- AUCell::AUCell_calcAUC(gene_sets, cells_rankings, nCores = ncores, ...)
        scores <- t(as.matrix(AUCell::getAUC(cells_AUC)))
        
    } else if (method == "GSVA") {
        if (!requireNamespace("GSVA", quietly = TRUE)) {
            stop("Please install 'GSVA' package to use method='GSVA'.")
        }
        message("Running GSVA scoring...")
        
        # In GSVA 1.50+, parameter specification is via gsvaParam
        if (utils::packageVersion("GSVA") >= "1.50.0") {
            param <- GSVA::gsvaParam(expr_mat, gene_sets, ...)
            gsva_res <- GSVA::gsva(param)
        } else {
            gsva_res <- GSVA::gsva(expr_mat, gene_sets, parallel.sz = ncores, ...)
        }
        scores <- t(as.matrix(gsva_res))
    }
    
    # Prefix colnames with name
    colnames(scores) <- paste0(name, "_", colnames(scores))
    
    # Add to colData
    for (cn in colnames(scores)) {
        colData(sce)[[cn]] <- scores[, cn]
    }
    
    # Register state
    metadata(sce)$sclet$analyses$geneset_scoring <- list(
        method = method,
        assay_use = assay_use,
        gene_sets = names(gene_sets),
        score_columns = colnames(scores),
        timestamp = Sys.time()
    )
    
    metadata(sce)$sclet$active$geneset_scoring <- "geneset_scoring"
    
    message(sprintf("Gene set scoring completed using %s. Added %d columns to colData.", method, ncol(scores)))
    return(sce)
}
