library(testthat)
library(sclet)
library(SingleCellExperiment)
library(ggplot2)

test_that("CellDimPlot and FeatureDimPlot work correctly", {
    counts <- matrix(rpois(10000, lambda = 1), ncol=100)
    rownames(counts) <- paste0("Gene", 1:100)
    colnames(counts) <- paste0("Cell", 1:100)
    sce <- SingleCellExperiment(assays = list(counts = counts))
    colLabels(sce) <- sample(c("TypeA", "TypeB"), 100, replace = TRUE)
    sce <- RunStandardPipeline(sce, nfeatures=50, npcs=10, dims=1:5, resolution=0.8)

    # Test CellDimPlot with active ident
    p1 <- CellDimPlot(sce)
    expect_s3_class(p1, "ggplot")
    
    # Test CellDimPlot with specific column
    sce$batch <- sample(c("B1", "B2"), 100, replace = TRUE)
    p2 <- CellDimPlot(sce, group.by = "batch")
    expect_s3_class(p2, "ggplot")
    
    # Test FeatureDimPlot
    p3 <- FeatureDimPlot(sce, features = c("Gene1", "Gene2"))
    expect_s3_class(p3, "ggplot")
})

test_that("GroupHeatmap and CellStatPlot work correctly", {
    counts <- matrix(rpois(10000, lambda = 1), ncol=100)
    rownames(counts) <- paste0("Gene", 1:100)
    colnames(counts) <- paste0("Cell", 1:100)
    sce <- SingleCellExperiment(assays = list(counts = counts))
    colLabels(sce) <- sample(c("TypeA", "TypeB"), 100, replace = TRUE)
    sce$batch <- sample(c("B1", "B2"), 100, replace = TRUE)
    sce <- RunStandardPipeline(sce, nfeatures=50, npcs=10, dims=1:5, resolution=0.8)

    # Test GroupHeatmap
    p4 <- GroupHeatmap(sce, features = c("Gene1", "Gene2"))
    expect_true(inherits(p4, "Heatmap") || inherits(p4, "ggplot") || inherits(p4, "pheatmap"))
    
    # Test CellStatPlot
    p5 <- CellStatPlot(sce)
    expect_s3_class(p5, "ggplot")
    
    p6 <- CellStatPlot(sce, split.by = "batch")
    expect_s3_class(p6, "ggplot")
})
