library(testthat)
library(sclet)
library(SingleCellExperiment)

test_that("RunIntegration fastMNN works correctly", {
    counts <- matrix(rpois(2000, lambda = 1), ncol=100)
    rownames(counts) <- paste0("Gene", 1:20)
    colnames(counts) <- paste0("Cell", 1:100)
    sce <- SingleCellExperiment(assays = list(counts = counts))
    sce$batch <- sample(c("B1", "B2"), 100, replace = TRUE)
    
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, nfeatures=10)
    
    suppressWarnings({
        sce <- RunIntegration(sce, method = "fastMNN", batch = "batch", PARAM = batchelor::FastMnnParam(d=5))
    })
    
    # The corrected embedding must be preserved (fastMNN's intended artifact)
    expect_true("corrected" %in% reducedDimNames(sce))
    expect_equal(sclet_get_active_reduction(sce), "corrected")
    expect_equal(sclet:::sclet_get_active_state(sce, "integration"), "fastmnn")
    
    # Original expression data must be preserved (not replaced by HVG-only 'reconstructed')
    expect_true("logcounts" %in% assayNames(sce))
    expect_equal(nrow(sce), 20)
    
    expect_true("corrected" %in% Layers(sce))
    expect_equal(sclet_get_active_layer(sce), "corrected")
    
    state <- get_analysis_context(sce)
    expect_equal(state$active$layer, "corrected")
    
    # With correct.all = FALSE (default) the integration is embedding-based: the
    # corrected layer resolves to the source expression for quantitative analyses.
    expect_true(isFALSE(get_integration(sce, "artifacts")$corrected_expression))
    expect_true(isTRUE(get_integration(sce, "artifacts")$embedding_only))
    expect_equal(get_integration(sce, "artifacts")$reduction, "corrected")
    
    # Check downstream
    sce <- RunPCA(sce, ncomponents=5)
    expect_equal(get_analysis_context(sce)$active$reduction, "PCA")
    
    # PCA should depend on integration batchcorrect
    pca_state <- sclet:::sclet_get_state_record(sce, "reduction", "pca")
    expect_equal(pca_state$inputs$integration$id, "fastmnn")
})

test_that("RunIntegration Harmony works correctly", {
    skip_if_not_installed("harmony")
    
    counts <- matrix(rpois(2000, lambda = 1), ncol=100)
    rownames(counts) <- paste0("Gene", 1:20)
    colnames(counts) <- paste0("Cell", 1:100)
    sce <- SingleCellExperiment(assays = list(counts = counts))
    sce$batch <- sample(c("B1", "B2"), 100, replace = TRUE)
    
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, nfeatures=10)
    sce <- ScaleData(sce)
    sce <- RunPCA(sce, ncomponents=5)
    
    sce <- RunIntegration(sce, method = "Harmony", batch = "batch")
    
    expect_true("HARMONY" %in% reducedDimNames(sce))
    expect_equal(sclet_get_active_reduction(sce), "HARMONY")
    
    # Check downstream
    sce <- RunUMAP(sce, dims=1:2)
    expect_equal(get_analysis_context(sce)$active$reduction, "UMAP")
    
    # UMAP should depend on Harmony
    umap_state <- sclet:::sclet_get_state_record(sce, "reduction", "umap")
    expect_equal(umap_state$inputs$reduction, "HARMONY")
    expect_equal(umap_state$inputs$integration, "harmony")
})
