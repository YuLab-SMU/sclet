sclet_matrix_colSums <- function(x, na.rm = FALSE, dims = 1L) {
    if (methods::is(x, "DelayedArray")) {
        return(MatrixGenerics::colSums(x, na.rm = na.rm, dims = dims))
    }
    if (methods::is(x, "Matrix")) {
        return(Matrix::colSums(x, na.rm = na.rm, dims = dims))
    }
    base::colSums(x, na.rm = na.rm, dims = dims)
}

sclet_matrix_rowSums <- function(x, na.rm = FALSE, dims = 1L) {
    if (methods::is(x, "DelayedArray")) {
        return(MatrixGenerics::rowSums(x, na.rm = na.rm, dims = dims))
    }
    if (methods::is(x, "Matrix")) {
        return(Matrix::rowSums(x, na.rm = na.rm, dims = dims))
    }
    base::rowSums(x, na.rm = na.rm, dims = dims)
}

sclet_matrix_colMeans <- function(x, na.rm = FALSE, dims = 1L) {
    if (methods::is(x, "DelayedArray")) {
        return(MatrixGenerics::colMeans(x, na.rm = na.rm, dims = dims))
    }
    if (methods::is(x, "Matrix")) {
        return(Matrix::colMeans(x, na.rm = na.rm, dims = dims))
    }
    base::colMeans(x, na.rm = na.rm, dims = dims)
}

sclet_matrix_rowMeans <- function(x, na.rm = FALSE, dims = 1L) {
    if (methods::is(x, "DelayedArray")) {
        return(MatrixGenerics::rowMeans(x, na.rm = na.rm, dims = dims))
    }
    if (methods::is(x, "Matrix")) {
        return(Matrix::rowMeans(x, na.rm = na.rm, dims = dims))
    }
    base::rowMeans(x, na.rm = na.rm, dims = dims)
}

sclet_matrix_rowSds <- function(x, na.rm = FALSE, dims = 1L) {
    if (methods::is(x, "DelayedArray")) {
        return(MatrixGenerics::rowSds(x, na.rm = na.rm, dims = dims))
    }
    if (methods::is(x, "Matrix")) {
        return(MatrixGenerics::rowSds(x, na.rm = na.rm, dims = dims))
    }
    apply(x, 1, stats::sd, na.rm = na.rm)
}

sclet_as_dgCMatrix <- function(x) {
    if (methods::is(x, "dgCMatrix")) {
        return(x)
    }
    if (methods::is(x, "Matrix")) {
        return(methods::as(x, "dgCMatrix"))
    }
    methods::as(Matrix::Matrix(x, sparse = TRUE), "dgCMatrix")
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