#' add QC metrics
#' 
#' @title QCMetrics
#' @param object a 'SingleCellExperiment' object
#' @return update object with two QC metrics ('nFeature_RNA' and 'nCount_RNA')
#' @importFrom scuttle perCellQCMetrics
#' @export
QCMetrics <- function(object) {
    qc_metrics <- scuttle::perCellQCMetrics(object)
    SummarizedExperiment::colData(object)$nFeature_RNA <- qc_metrics$detected
    SummarizedExperiment::colData(object)$nCount_RNA <- qc_metrics$sum
    return(object)
}

#' calculate percentage of genes that matched a pattern
#' 
#' @title PercentageFeatureSet
#' @param object a SingleCellExperiment object
#' @param pattern the pattern to search for
#' @param feature search pattern from feature of the object. If NULL, use the `rownames(object)`, otherwise, use `rowData(object)[[feature]]`
#' @return the percentage of each cell
#' @export
PercentageFeatureSet <- function(object, pattern = NULL, feature=NULL) {
    if (is.null(feature)) {
        features <- rownames(object)
    } else {
        features <- rowData(object)[[feature]]
    }

    has_pattern <- grep(pattern, features)
    qc_metrics <- scuttle::perCellQCMetrics(object, subsets = list(pattern=has_pattern))
    return(qc_metrics$subsets_pattern_percent)
}

#' normalize data
#' 
#' @title NormalizeData
#' @param object a SingleCellExperiment object
#' @param scale.factor scale factor
#' @return a SingleCellExperiment object with 'logcounts' assay created/updated.
#' @importFrom SummarizedExperiment 'assay<-'
#' @importFrom yulab.utils get_fun_from_pkg
#' @importFrom SummarizedExperiment assay
#' @importFrom methods as
#' @export
NormalizeData <- function(object, scale.factor = 10000) {
    # No longer force conversion to dgCMatrix. scuttle::logNormCounts handles DelayedMatrix efficiently.
    # if (is(assay(object, "counts"), "DelayedMatrix")) {
    #    assay(object, "counts") <- as(assay(object, "counts"), "dgCMatrix")
    # }
    
    # Use scuttle::logNormCounts which is standard for SCE
    libsize <- MatrixGenerics::colSums(SummarizedExperiment::assay(object, "counts"))
    size_factors <- libsize / scale.factor
    object <- scuttle::logNormCounts(
        object,
        size.factors = size_factors,
        center.size.factors = FALSE,
        name = "logcounts"
    )
    
    return(object)
}

#' identify variable features
#' 
#' @title FindVariableFeatures
#' @param object a SingleCellExperiment object
#' @param nfeatures number of features to be selected as highly variable features
#' @param method one of 'seurat' or 'scran'
#' @param ... additional parameters for 'method = 'scran', see also `scran::modelGeneVar()`
#' @return an updated SingleCellExperiment object with identified highly variable features
#' @importFrom SingleCellExperiment counts
#' @importFrom SingleCellExperiment logcounts
#' @importFrom SummarizedExperiment 'rowData<-'
#' @importFrom stats loess
#' @export
FindVariableFeatures <- function(object, nfeatures = 2000, method = "seurat", ...) {
    method <- match.arg(method, c("seurat", "scran"))

    if (method == "seurat") {
        if (!requireNamespace("Seurat", quietly = TRUE)) {
             warning("Seurat is not installed. Falling back to scran method.")
             hvf.info <- scran::modelGeneVar(object, ...)
             method <- "scran"
        } else {
             hvf.info <- FindVariableFeatures_seurat(object)
        }
    } else {
        hvf.info <- scran::modelGeneVar(object, ...)
    }

    rd <- rowData(object)
    SummarizedExperiment::rowData(object) <- cbind(rd, hvf.info[rownames(rd),])
    object@metadata$nVariableFeatures <- nfeatures
    object@metadata$hvgmethod <- method
    object@metadata$hvgcols <- names(hvf.info)

    return(object)
}

FindVariableFeatures_seurat <- function(object) {
    # Check for Seurat namespace again to be safe
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Seurat package is required for method='seurat'")
    }

    sce <- object
    if (!is.null(sce@metadata$nVariableFeatures)) {
        rd <- rowData(sce)
        items <- c("mean", "variance", "variance.expected", "variance.standardized")
        if (all(items %in% names(rd))) {
            return(rd[, items])
        }
 
    }

    object <- counts(sce)
    SparseRowVar2 <- get_fun_from_pkg('Seurat', "SparseRowVar2")
    SparseRowVarStd <- get_fun_from_pkg('Seurat', "SparseRowVarStd")

    ## taken from Seurat
    clip.max <- sqrt(x = ncol(x = object))
    hvf.info <- data.frame(mean = Matrix::rowMeans(object))
    hvf.info$variance <- SparseRowVar2(
      mat = object,
      mu = hvf.info$mean,
      display_progress = TRUE
    )
    hvf.info$variance.expected <- 0
    hvf.info$variance.standardized <- 0
    not.const <- hvf.info$variance > 0
    fit <- loess(
      formula = log10(x = variance) ~ log10(x = mean),
      data = hvf.info[not.const, ],
      span = .3
    )
    hvf.info$variance.expected[not.const] <- 10 ^ fit$fitted
    hvf.info$variance.standardized <- SparseRowVarStd(
      mat = object,
      mu = hvf.info$mean,
      sd = sqrt(hvf.info$variance.expected),
      vmax = clip.max,
      display_progress = TRUE
    )
    ## end of seurat code

    return(hvf.info)
}


#' get variable features
#' 
#' @title VariableFeatures
#' @param object a SingleCellExperiment object
#' @param method one of 'seurat' or 'scran'
#' @param ... additional parameters for 'method = "scran"', see also `scran::getTopHVGs()`
#' @return highly variable features
#' @export
VariableFeatures <- function(object, method = "seurat", ...) {
    method <- match.arg(method, c("seurat", "scran"))

    nfeatures <- object@metadata$nVariableFeatures
    if (is.null(nfeatures)) {
        stop("You should run 'FindVariableFeatures' first.")
    }

    if (method == "scran") {
        res <- scran::getTopHVGs(object, n=nfeatures, ...)
    } else {
        res <- getTopHVGs_seurat(object, n=nfeatures)
    }

    return(res)
}

getTopHVGs_seurat <- function(object, n) {
    if (is(object, "SingleCellExperiment")) {
        d <- rowData(object)
        d <- d[rownames(d) %in% rownames(object), ]
    } else {
        d <- object
    }
    i <- order(d$variance.standardized, decreasing = TRUE)
    d <- d[i, ]
    res <- rownames(d)[1:n]
    return(res)
}

#' scale data
#' 
#' @title ScaleData
#' @param object a SingleCellExperiment object
#' @param features selected features to be scaled, default is all features
#' @param assay selected assay to be scaled
#' @return an updated SingleCellExperiment object with 'scaled' assay
#' @export
#' @importFrom stats sd
ScaleData <- function(object, features = NULL, assay = "logcounts") {
    if (is.null(features)) {
        # all genes
        features <- rownames(object)
    } 

    # For DelayedArray, we must avoid densifying the matrix.
    # scater::logNormCounts already provides logcounts.
    # Standard scaling (subtract mean, divide sd) on sparse/HDF5 matrix results in a DENSE matrix.
    # This is bad for memory.
    # However, 'scaled' slot in Seurat is typically dense. 
    # If we want to support big data, we should use DelayedMatrix to represent the scaled matrix lazily.
    
    mat <- assay(object, assay)[features, , drop=FALSE]
    
    # Calculate row means and sds efficiently
    # DelayedMatrixStats or MatrixGenerics should be used
    if (is(mat, "DelayedMatrix")) {
        # DelayedArray logic
        rm <- MatrixGenerics::rowMeans(mat)
        rs <- MatrixGenerics::rowSds(mat)
    } else {
        # In-memory logic (could be sparse or dense)
        if (is(mat, "dgCMatrix")) {
             rm <- Matrix::rowMeans(mat)
             # Sparse matrix rowSds calculation can be tricky if not careful, 
             # but MatrixGenerics handles it or we use sparseMatrixStats if available.
             # Fallback to apply for now if not available, but let's assume MatrixGenerics works.
             rs <- MatrixGenerics::rowSds(mat)
        } else {
             rm <- rowMeans(mat)
             rs <- apply(mat, 1, stats::sd)
        }
    }
    
    rs[rs == 0] <- 0.01
    
    # Create a DelayedMatrix that represents the scaled data
    # (mat - rm) / rs
    # In DelayedArray, we can do arithmetic directly.
    
    if (is(mat, "DelayedMatrix")) {
        # This creates a DelayedMatrix (lazy evaluation)
        # We need to broadcast the subtraction and division.
        # DelayedArray supports standard broadcasting for vector-matrix ops if dimensions match?
        # Usually it matches on columns? No, R matrices match on columns (recycling).
        # mat is genes x cells. rm is genes.
        # (mat - rm) needs rm to be repeated for each cell.
        
        # DelayedArray arithmetic:
        # We can use sweep-like operations or just rely on broadcasting if implemented correctly.
        # Alternatively, create a DelayedMatrix backend.
        
        scaled <- (mat - rm) / rs
        
    } else {
        # For memory matrices, we might still want to return a DelayedMatrix if it's too big?
        # But if it's in memory, maybe just compute it.
        # NOTE: Seurat ScaleData forces dense. 
        # Here we try to keep it compatible but lazy if possible.
        
        # If the user explicitly wants HDF5/BigData, they should have used Read10X(backend="HDF5").
        # So here, if input is DelayedMatrix, output is DelayedMatrix.
        # If input is memory, output is memory (dense, as scaling destroys sparsity).
        
        scaled <- (mat - rm) / rs
    }
    
    SummarizedExperiment::assay(object, "scaled") <- scaled
    return(object)
}
