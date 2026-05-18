library(testthat)
library(sclet)
library(SingleCellExperiment)

test_that("RunKNNPredict and ProjectionPlot work correctly", {
    # Create reference
    ref_counts <- matrix(rpois(1000, lambda = 1), ncol=50)
    rownames(ref_counts) <- paste0("Gene", 1:20)
    colnames(ref_counts) <- paste0("RefCell", 1:50)
    ref <- SingleCellExperiment(assays = list(logcounts = ref_counts))
    colLabels(ref) <- sample(c("RefTypeA", "RefTypeB"), 50, replace = TRUE)
    reducedDim(ref, "UMAP") <- matrix(rnorm(100), ncol=2)

    # Create query
    query_counts <- matrix(rpois(1000, lambda = 1), ncol=50)
    rownames(query_counts) <- paste0("Gene", 1:20)
    colnames(query_counts) <- paste0("QueryCell", 1:50)
    query <- SingleCellExperiment(assays = list(counts = query_counts))
    query <- NormalizeData(query)
    query <- FindVariableFeatures(query, nfeatures=10)
    reducedDim(query, "UMAP") <- matrix(rnorm(100), ncol=2)

    # Test RunKNNPredict
    query <- RunKNNPredict(query, ref = ref, labels = "label", k = 3)
    
    expect_true("knn_map_labels" %in% colnames(colData(query)))
    expect_true("knn_map_score" %in% colnames(colData(query)))
    
    # Check mapping state
    map_state <- get_mapping(query)
    expect_false(is.null(map_state))
    expect_equal(map_state$method, "KNNPredict")
    expect_equal(map_state$artifacts$labels_col, "knn_map_labels")

    # Test ProjectionPlot
    p <- ProjectionPlot(query, ref, reduction="UMAP")
    expect_s3_class(p, "ggplot")
})
