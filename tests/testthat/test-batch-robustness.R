
test_that("BatchRemover and subsequent clustering workflow works", {
    skip_if_not_installed("scuttle")
    skip_if_not_installed("batchelor")
    skip_if_not_installed("scran")
    
    # Mock data
    set.seed(123)
    sce1 <- scuttle::mockSCE(ncells = 50, ngenes = 200)
    sce1 <- NormalizeData(sce1)
    sce1$batch <- "batch1"
    
    sce2 <- scuttle::mockSCE(ncells = 50, ngenes = 200)
    sce2 <- NormalizeData(sce2)
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
    expect_true(has_integration(combined, id = "merged_inputs"))
    expect_equal(get_integration(combined, id = "merged_inputs")$method, "sce_merge")
    expect_equal(get_integration(combined, id = "merged_inputs")$summary$n_inputs, 2)
    
    # Batch Correction
    rm_batch <- BatchRemover(combined)
    
    # Check if BatchRemover returns a valid SCE
    expect_s4_class(rm_batch, "SingleCellExperiment")
    expect_true("corrected" %in% Layers(rm_batch))
    expect_equal(DefaultLayer(rm_batch), "corrected")
    expect_equal(sclet_get_layer(rm_batch, "corrected")$role, "corrected")
    expect_equal(sclet_get_active_state(rm_batch, "integration"), "batchcorrect")
    expect_equal(get_integration(rm_batch, "method"), "batchelor::batchCorrect")
    expect_equal(get_integration(rm_batch, "artifacts")$layer, "corrected")
    expect_equal(get_integration(rm_batch, "inputs")$merge$id, "merged_inputs")
    
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
    sce <- NormalizeData(sce)
    
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

test_that("BatchRemover preserves scrapper HVG workflow", {
    skip_if_not_installed("scuttle")
    skip_if_not_installed("batchelor")
    skip_if_not_installed("scrapper")

    set.seed(202)
    sce1 <- scuttle::mockSCE(ncells = 36, ngenes = 120)
    sce1 <- NormalizeData(sce1)
    sce1$batch <- "batch1"
    sce1 <- FindVariableFeatures(sce1, nfeatures = 30, method = "scrapper")

    sce2 <- scuttle::mockSCE(ncells = 34, ngenes = 120)
    sce2 <- NormalizeData(sce2)
    sce2$batch <- "batch2"
    sce2 <- FindVariableFeatures(sce2, nfeatures = 30, method = "scrapper")

    combined <- sce_merge(list(pbmc1 = sce1, pbmc2 = sce2))
    rm_batch <- BatchRemover(combined, nHVG = 30)

    expect_equal(get_hvg(rm_batch, "method"), "scrapper")
    expect_equal(get_integration(rm_batch, "method"), "batchelor::batchCorrect")
    expect_true("corrected" %in% Layers(rm_batch))

    expect_no_error({
        vf <- VariableFeatures(rm_batch)
    })
    expect_equal(length(vf), 30)
    expect_true(all(vf %in% rownames(rm_batch)))
})
