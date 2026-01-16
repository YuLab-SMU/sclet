test_that("Idents gets/sets colLabels", {
  library(SingleCellExperiment)
  sce <- SingleCellExperiment(assays = list(counts = matrix(rpois(100, 5), ncol=10)))
  colnames(sce) <- paste0("Cell", 1:10)
  
  # Default NULL
  expect_null(Idents(sce))
  
  # Set via colLabels
  colLabels(sce) <- factor(rep(c("A", "B"), each=5))
  ids <- Idents(sce)
  expect_equal(as.character(ids), rep(c("A", "B"), each=5))
  expect_equal(names(ids), colnames(sce))
  
  # Fallback to colData$label if colLabels is NULL
  colLabels(sce) <- NULL
  colData(sce)$label <- factor(rep(c("X", "Y"), each=5))
  # Should sync back to colLabels
  ids <- Idents(sce)
  expect_equal(as.character(ids), rep(c("X", "Y"), each=5))
  expect_equal(as.character(colLabels(sce)), rep(c("X", "Y"), each=5))
})

test_that("subset_cell works", {
  library(SingleCellExperiment)
  sce <- SingleCellExperiment(assays = list(counts = matrix(rpois(200, 10), ncol=20)))
  colData(sce)$nFeature_RNA <- rnorm(20, mean=100, sd=10)
  
  # Just run to check no error
  sce_sub <- subset_cell(sce, feature = "nFeature_RNA", method = "sd", n = 1)
  expect_s4_class(sce_sub, "SingleCellExperiment")
})
