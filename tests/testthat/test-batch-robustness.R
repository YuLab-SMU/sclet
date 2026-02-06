
test_that("BatchRemover and subsequent clustering workflow works", {
    skip_if_not_installed("scuttle")
    skip_if_not_installed("batchelor")
    skip_if_not_installed("scran")
    
    # Mock data
    set.seed(123)
    sce1 <- scuttle::mockSCE(ncells = 50, ngenes = 200)
    sce1 <- scuttle::logNormCounts(sce1)
    sce1$batch <- "batch1"
    
    sce2 <- scuttle::mockSCE(ncells = 50, ngenes = 200)
    sce2 <- scuttle::logNormCounts(sce2)
    sce2$batch <- "batch2"
    
    # Pre-process (simulate user workflow)
    # We use scran method by default now
    sce1 <- FindVariableFeatures(sce1, nfeatures = 50)
    sce2 <- FindVariableFeatures(sce2, nfeatures = 50)
    
    # Merge
    pbmc_list <- list(pbmc1 = sce1, pbmc2 = sce2)
    combined <- sce_merge(pbmc_list)
    
    expect_s4_class(combined, "SingleCellExperiment")
    expect_true("batch" %in% names(colData(combined)))
    
    # Batch Correction
    rm_batch <- BatchRemover(combined)
    
    # Check if BatchRemover returns a valid SCE
    expect_s4_class(rm_batch, "SingleCellExperiment")
    
    # The critical failure point was running VariableFeatures on this result
    # because it might look for assays that don't exist or are named differently
    expect_no_error({
        vf <- VariableFeatures(rm_batch)
    })
    
    expect_true(length(vf) > 0)
    
    # Test Clustering workflow on this object
    expect_no_error({
        # We use "reconstructed" assay if it exists (standard for batchelor output)
        assay_to_use <- "logcounts"
        if ("reconstructed" %in% assayNames(rm_batch)) {
            assay_to_use <- "reconstructed"
        }
        
        rm_batch <- runPCA(rm_batch, subset_row = vf, exprs_values = assay_to_use)
        rm_batch <- FindNeighbors(rm_batch, dims = 1:5)
        rm_batch <- FindClusters(rm_batch, resolution = 0.5)
    })
    
    expect_true(!is.null(colLabels(rm_batch)))
})

test_that("VariableFeatures handles missing assays gracefully (scran method)", {
    skip_if_not_installed("scuttle")
    
    set.seed(42)
    sce <- scuttle::mockSCE(ncells = 20, ngenes = 100)
    sce <- scuttle::logNormCounts(sce)
    
    # Run modelGeneVar (via FindVariableFeatures)
    sce <- FindVariableFeatures(sce, method="scran", nfeatures=20)
    
    # Simulate a scenario where logcounts/counts are missing or not found
    # This mimics the state where `assay(sce, i="logcounts")` failed in the bug report
    assay(sce, "counts") <- NULL 
    assay(sce, "logcounts") <- NULL 
    
    # VariableFeatures should still work because it should use rowData
    expect_no_error({
        vf <- VariableFeatures(sce)
    })
    expect_equal(length(vf), 20)
})
