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
