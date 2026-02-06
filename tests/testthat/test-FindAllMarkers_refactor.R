test_that("FindAllMarkers works with refactored presto implementation", {
    # Check if presto is installed, otherwise skip or expect installation
    if (!requireNamespace("presto", quietly = TRUE)) {
        skip("presto not installed")
    }

    # Create a mock SCE
    set.seed(123)
    n_genes <- 100
    n_cells <- 200
    counts <- matrix(rpois(n_genes * n_cells, lambda = 5), nrow = n_genes)
    rownames(counts) <- paste0("Gene", 1:n_genes)
    colnames(counts) <- paste0("Cell", 1:n_cells)
    
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
    
    # We need logcounts for FindAllMarkers (it uses assay(object, "logcounts"))
    # Manual log normalization to avoid dependency loop issues in test if sclet not fully loaded
    libsize <- colSums(counts)
    size_factors <- libsize / 10000
    logcounts <- log1p(t(t(counts) / size_factors))
    SummarizedExperiment::assay(sce, "logcounts") <- logcounts
    
    # Add clusters
    clusters <- factor(rep(c("A", "B", "C"), length.out = n_cells))
    SingleCellExperiment::colLabels(sce) <- clusters
    
    # Run FindAllMarkers
    # Note: sclet package needs to be loaded or functions available. 
    # In devtools::test(), they are available.
    
    markers <- FindAllMarkers(sce, min.pct = 0, logfc.threshold = 0, return.thresh = 1)
    
    expect_s3_class(markers, "data.frame")
    
    # Check column names
    expected_cols <- c("gene", "cluster", "avg_log2FC", "pct.1", "pct.2", "pval", "padj")
    expect_true(all(expected_cols %in% colnames(markers)))
    
    expect_true(nrow(markers) > 0)
    
    # Check validity of values
    expect_true(all(markers$pval >= 0 & markers$pval <= 1))
    expect_true(all(markers$pct.1 >= 0 & markers$pct.1 <= 1))
    expect_true(all(markers$pct.2 >= 0 & markers$pct.2 <= 1))
    
    # Check if cluster column contains correct levels
    expect_true(all(unique(markers$cluster) %in% c("A", "B", "C")))
})
