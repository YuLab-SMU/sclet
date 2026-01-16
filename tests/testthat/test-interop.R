test_that("Interop functions check dependencies", {
    # If zellkonverter is missing, ReadH5AD should fail
    if (!requireNamespace("zellkonverter", quietly = TRUE)) {
        expect_error(ReadH5AD("dummy.h5ad"), "zellkonverter")
    }
    
    # If Seurat is missing, as.Seurat should fail
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        expect_error(as.Seurat(NULL), "Seurat")
    }
})
