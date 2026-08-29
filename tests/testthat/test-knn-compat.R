library(testthat)
library(sclet)
library(SingleCellExperiment)

# Regression tests for YuLab-SMU/sclet#25:
# BiocNeighbors 2.0.x-2.4.x cannot dispatch queryKNN()/findKNN() on
# dgCMatrix inputs, so sclet must fall back to ordinary matrices there
# while keeping the sparse path on BiocNeighbors >= 2.6.

test_that("sparse probe returns a consistent flag", {
    supported <- sclet:::sclet_biocneighbors_accepts_sparse()
    expect_true(is.logical(supported) && length(supported) == 1 && !is.na(supported))
})

test_that("sclet_as_knn_matrix keeps sparse input when supported", {
    m <- methods::as(Matrix::rsparsematrix(20, 10, density = 0.3), "dgCMatrix")
    out <- sclet:::sclet_as_knn_matrix(m)
    if (sclet:::sclet_biocneighbors_accepts_sparse()) {
        expect_identical(out, m)
        expect_s4_class(out, "dgCMatrix")
    } else {
        expect_true(is.matrix(out))
    }
    expect_equal(as.matrix(out), as.matrix(m))
})

test_that("sclet_as_knn_matrix passes ordinary matrices through untouched", {
    m <- matrix(rnorm(100), nrow = 10)
    expect_identical(sclet:::sclet_as_knn_matrix(m), m)
})

test_that("sclet_as_knn_matrix falls back to dense on old BiocNeighbors", {
    testthat::local_mocked_bindings(
        sclet_biocneighbors_accepts_sparse = function() FALSE,
        .package = "sclet"
    )
    m <- methods::as(Matrix::rsparsematrix(20, 10, density = 0.3), "dgCMatrix")
    out <- sclet:::sclet_as_knn_matrix(m)
    expect_true(is.matrix(out))
    expect_false(methods::is(out, "dgCMatrix"))
    expect_equal(out, as.matrix(m))
})

test_that("RunKNNPredict survives old-BiocNeighbors dense fallback (#25)", {
    testthat::local_mocked_bindings(
        sclet_biocneighbors_accepts_sparse = function() FALSE,
        .package = "sclet"
    )
    ref_counts <- matrix(rpois(1000, lambda = 1), ncol = 50)
    rownames(ref_counts) <- paste0("Gene", 1:20)
    colnames(ref_counts) <- paste0("RefCell", 1:50)
    ref <- SingleCellExperiment(
        assays = list(logcounts = Matrix::Matrix(ref_counts, sparse = TRUE))
    )
    colLabels(ref) <- sample(c("RefTypeA", "RefTypeB"), 50, replace = TRUE)

    query_counts <- matrix(rpois(1000, lambda = 1), ncol = 50)
    rownames(query_counts) <- paste0("Gene", 1:20)
    colnames(query_counts) <- paste0("QueryCell", 1:50)
    query <- SingleCellExperiment(assays = list(counts = query_counts))
    query <- NormalizeData(query)
    query <- FindVariableFeatures(query, nfeatures = 10)
    # mirror the reported issue exactly: sparse genes x cells layers
    SummarizedExperiment::assay(query, "logcounts") <- Matrix::Matrix(
        as.matrix(SummarizedExperiment::assay(query, "logcounts")),
        sparse = TRUE
    )

    query <- RunKNNPredict(query, ref = ref, labels = "label", k = 3)
    expect_true("knn_map_labels" %in% colnames(colData(query)))
    expect_true("knn_map_score" %in% colnames(colData(query)))
})
