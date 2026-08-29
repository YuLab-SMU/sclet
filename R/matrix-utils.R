sclet_matrix_colSums <- function(x, na.rm = FALSE, dims = 1L) {
    if (methods::is(x, "DelayedArray")) {
        return(MatrixGenerics::colSums(x, na.rm = na.rm))
    }
    if (methods::is(x, "Matrix")) {
        return(Matrix::colSums(x, na.rm = na.rm))
    }
    base::colSums(x, na.rm = na.rm, dims = dims)
}

sclet_matrix_rowSums <- function(x, na.rm = FALSE, dims = 1L) {
    if (methods::is(x, "DelayedArray")) {
        return(MatrixGenerics::rowSums(x, na.rm = na.rm))
    }
    if (methods::is(x, "Matrix")) {
        return(Matrix::rowSums(x, na.rm = na.rm))
    }
    base::rowSums(x, na.rm = na.rm, dims = dims)
}

sclet_matrix_colMeans <- function(x, na.rm = FALSE, dims = 1L) {
    if (methods::is(x, "DelayedArray")) {
        return(MatrixGenerics::colMeans(x, na.rm = na.rm))
    }
    if (methods::is(x, "Matrix")) {
        return(Matrix::colMeans(x, na.rm = na.rm))
    }
    base::colMeans(x, na.rm = na.rm, dims = dims)
}

sclet_matrix_rowMeans <- function(x, na.rm = FALSE, dims = 1L) {
    if (methods::is(x, "DelayedArray")) {
        return(MatrixGenerics::rowMeans(x, na.rm = na.rm))
    }
    if (methods::is(x, "Matrix")) {
        return(Matrix::rowMeans(x, na.rm = na.rm))
    }
    base::rowMeans(x, na.rm = na.rm, dims = dims)
}

#' @importFrom MatrixGenerics rowSds
sclet_matrix_rowSds <- function(x, na.rm = FALSE, dims = 1L) {
    if (methods::is(x, "DelayedArray")) {
        return(MatrixGenerics::rowSds(x, na.rm = na.rm))
    }
    if (methods::is(x, "Matrix")) {
        return(MatrixGenerics::rowSds(x, na.rm = na.rm))
    }
    apply(x, 1, stats::sd, na.rm = na.rm)
}

sclet_as_dgCMatrix <- function(x) {
    if (methods::is(x, "dgCMatrix")) {
        return(x)
    }
    if (methods::is(x, "DelayedArray")) {
        return(methods::as(x, "dgCMatrix"))
    }
    if (methods::is(x, "Matrix")) {
        return(methods::as(x, "dgCMatrix"))
    }
    methods::as(Matrix::Matrix(x, sparse = TRUE), "dgCMatrix")
}

# BiocNeighbors 2.0.x-2.4.x (Bioconductor 3.20-3.22, R 4.4/4.5) stopped
# coercing X to an ordinary matrix and only dispatch on `matrix`, so
# queryKNN()/findKNN() fail with 'unable to find an inherited method for
# function ... for signature X = "dgCMatrix"' (YuLab-SMU/sclet#25).
# Sparse-matrix support returned in BiocNeighbors 2.6.0 (Bioconductor 3.23),
# which needs R >= 4.6 and is out of reach for users on older R versions.
# Detect support through S4 dispatch instead of pinning version numbers, so
# Bioc devel (2.5.x) and future releases keep the sparse path automatically.
sclet_biocneighbors_accepts_sparse <- function() {
    if (!requireNamespace("BiocNeighbors", quietly = TRUE)) {
        return(FALSE)
    }
    tryCatch(
        {
            # findKNN() and queryKNN() share the same coercion behaviour, so
            # probing one generic covers both call sites.
            methods::selectMethod(
                BiocNeighbors::queryKNN,
                signature = c(X = "dgCMatrix", BNPARAM = "missing")
            )
            TRUE
        },
        error = function(e) FALSE
    )
}

# Prepare a cells x features matrix for BiocNeighbors::queryKNN()/findKNN().
# Keep sparse inputs sparse when the installed BiocNeighbors can dispatch on
# them; otherwise fall back to an ordinary matrix, materializing the dense
# block only on the affected BiocNeighbors versions.
sclet_as_knn_matrix <- function(x) {
    if (is.matrix(x) || sclet_biocneighbors_accepts_sparse()) {
        return(x)
    }
    cli::cli_alert_info(
        "BiocNeighbors < 2.6 does not accept sparse matrices; ",
        "converting to a dense matrix (may increase memory usage)"
    )
    as.matrix(x)
}

sclet_extract_cell_feature_matrix <- function(object, assay_name, features = NULL) {
    mat <- SummarizedExperiment::assay(object, assay_name)
    if (!is.null(features)) {
        features <- intersect(features, rownames(object))
        if (length(features) == 0) {
            stop("No overlapping features found in assay '", assay_name, "'.")
        }
        mat <- mat[features, , drop = FALSE]
    }
    sclet_as_dgCMatrix(Matrix::t(mat))
}