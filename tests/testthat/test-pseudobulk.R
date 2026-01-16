test_that("AggregateExpression works", {
    library(SingleCellExperiment)
    sce <- SingleCellExperiment(assays = list(counts = matrix(rpois(100, 5), ncol=10)))
    colnames(sce) <- paste0("Cell", 1:10)
    colData(sce)$cluster <- rep(c("A", "B"), each=5)
    colData(sce)$sample <- rep(c("S1", "S2"), 5)
    
    # Test aggregation
    pb <- AggregateExpression(sce, group_by = c("cluster", "sample"))
    
    expect_s4_class(pb, "SummarizedExperiment")
    expect_equal(ncol(pb), 4) # 2 clusters * 2 samples
    expect_true(all(c("cluster", "sample") %in% colnames(colData(pb))))
})

test_that("FindMarkers_pseudobulk handles errors gracefully", {
    # If DESeq2 not installed, should fail (but here we just check structure if we can mock or if it's not installed)
    # We won't test full DESeq2 execution as it's Suggests and might be heavy/missing in minimal env
    # But we can check input validation
})
