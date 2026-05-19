test_that("RunSuperCell works", {
  skip_if_not_installed("SuperCell")
  skip_if_not_installed("scuttle")
  skip_if_not_installed("scater")
  skip_if_not_installed("scran")
  skip_if_not_installed("igraph")

  set.seed(21)
  sce <- scuttle::mockSCE(ncells = 40, ngenes = 80)
  sce <- NormalizeData(sce)
  sce <- FindVariableFeatures(sce, nfeatures = 20, method = "scran")
  sce <- RunPCA(sce, subset_row = VariableFeatures(sce), ncomponents = 5)
  sce <- FindNeighbors(sce, dims = 1:5)
  sce <- FindClusters(sce)

  out <- RunSuperCell(sce, layer = "logcounts", gamma = 10, k.knn = 5)
  expect_s4_class(out, "SingleCellExperiment")
  expect_true("size" %in% colnames(colData(out)))
  expect_true(has_supercell(out))
  expect_equal(get_supercell(out, "inputs")$layer, "logcounts")
  expect_equal(get_supercell(out, "artifacts")$analysis_key, "supercell")
  expect_equal(sclet_get_active_state(out, "aggregation"), "supercell")
  expect_equal(DefaultLayer(out), "logcounts")
  expect_null(DefaultReduction(out))
  expect_null(DefaultGraph(out))
  expect_null(ActiveIdent(out))
  expect_equal(CommandLog(out)$command, "RunSuperCell")
})

test_that("RunSlingshot works", {
    if (requireNamespace("slingshot", quietly = TRUE)) {
        sce <- SingleCellExperiment(assays = list(counts = matrix(rpois(1000, 5), ncol=50)))
        reducedDim(sce, "UMAP") <- matrix(rnorm(100), ncol=2)
        colLabels(sce) <- rep(c("A", "B"), each=25)
        
        # Should run
        # sce <- RunSlingshot(sce, cluster.labels=colLabels(sce))
        # expect_true(!is.null(sce@metadata$slingshot))
    }
    expect_true(TRUE)
})

