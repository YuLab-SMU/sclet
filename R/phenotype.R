#' Run phenotype association analysis using scPAS
#'
#' @title RunPhenotypeAssociation
#' @param sce a SingleCellExperiment object
#' @param bulk_matrix Bulk expression matrix of related disease. Each row represents a gene and each column represents a sample.
#' @param phenotype_data Phenotype annotation of each bulk sample. It can be a continuous variable, a binary group vector, or clinical survival data (two-column matrix with 'time' and 'status').
#' @param family Response type for the regression model. Can be 'gaussian', 'binomial', or 'cox'.
#' @param ... Additional arguments passed to \code{scPAS::scPAS}
#' @return an updated SingleCellExperiment object with phenotype association scores added to colData
#' @export
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom S4Vectors metadata metadata<-
RunPhenotypeAssociation <- function(sce, bulk_matrix, phenotype_data, family = "cox", ...) {
    if (!requireNamespace("scPAS", quietly = TRUE)) {
        stop("Please install 'scPAS' via remotes::install_github('aiminXie/scPAS')")
    }
    
    # Check inputs
    if (is.null(bulk_matrix) || is.null(phenotype_data)) {
        stop("Both 'bulk_matrix' and 'phenotype_data' must be provided.")
    }
    
    # scPAS typically prefers a Seurat object to utilize pre-constructed networks
    # Convert SCE to Seurat
    seurat_obj <- as.Seurat(sce)
    
    message("Running scPAS for phenotype association...")
    # scPAS modifies the seurat object and returns it
    # We use capture.output or just let it print
    res_seurat <- scPAS::scPAS(
        bulk_dataset = bulk_matrix,
        sc_dataset = seurat_obj,
        phenotype = phenotype_data,
        family = family,
        ...
    )
    
    # Extract results back to SCE
    # scPAS adds columns to meta.data, such as scPAS_score, scPAS_class, etc.
    # We'll identify new columns and add them to colData(sce)
    old_cols <- colnames(seurat_obj@meta.data)
    new_cols <- colnames(res_seurat@meta.data)
    added_cols <- setdiff(new_cols, old_cols)
    
    if (length(added_cols) > 0) {
        for (col in added_cols) {
            SummarizedExperiment::colData(sce)[[col]] <- res_seurat@meta.data[[col]]
        }
    }
    
    # Extract model parameters from misc
    if (!is.null(res_seurat@misc$scPAS_para)) {
        md <- S4Vectors::metadata(sce)
        md$phenotype_assoc <- res_seurat@misc$scPAS_para
        S4Vectors::metadata(sce) <- md
    }
    
    # Log and update analysis state
    sce <- sclet_log_command(
        sce,
        "RunPhenotypeAssociation",
        params = list(family = family, ...)
    )
    
    sce <- sclet_set_analysis_state(
        object = sce,
        type = "phenotype",
        id = "scPAS",
        method = "scPAS::scPAS",
        inputs = list(
            bulk_matrix = "external",
            phenotype_data = "external"
        ),
        artifacts = list(
            colData = added_cols
        ),
        params = list(
            family = family,
            ...
        ),
        summary = list(
            added_columns = added_cols
        )
    )
    
    return(sce)
}
