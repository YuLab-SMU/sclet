test_that("Idents requires an explicit active identity source", {
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
  
  # Non-default colData labels no longer become active identities implicitly
  colLabels(sce) <- NULL
  colData(sce)$manual_label <- factor(rep(c("X", "Y"), each=5))
  expect_null(Idents(sce))

  ActiveIdent(sce) <- "manual_label"
  ids <- Idents(sce)
  expect_equal(as.character(ids), rep(c("X", "Y"), each=5))
})

test_that("subset_cell works", {
  library(SingleCellExperiment)
  sce <- SingleCellExperiment(assays = list(counts = matrix(rpois(200, 10), ncol=20)))
  colData(sce)$nFeature_RNA <- rnorm(20, mean=100, sd=10)
  
  # Just run to check no error
  sce_sub <- subset_cell(sce, feature = "nFeature_RNA", method = "sd", n = 1)
  expect_s4_class(sce_sub, "SingleCellExperiment")
})

test_that("FindVariableFeatures legacy seurat method redirects to scran", {
  library(SingleCellExperiment)
  skip_if_not_installed("scuttle")
  skip_if_not_installed("scran")

  sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
  sce <- NormalizeData(sce)
  expect_warning(
    sce <- FindVariableFeatures(sce, nfeatures = 10, method = "seurat"),
    "deprecated"
  )
  expect_equal(sclet::get_hvg(sce, "method"), "scran")
})

test_that("legacy seurat HVG metadata no longer becomes the default path", {
  library(SingleCellExperiment)
  skip_if_not_installed("scuttle")

  sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
  rowData(sce)$mean <- seq_len(nrow(sce))
  rowData(sce)$variance.standardized <- rev(seq_len(nrow(sce)))
  sce <- sclet:::sclet_set_hvg_state(
    sce,
    nfeatures = 10,
    method = "seurat",
    hvgcols = c("mean", "variance.standardized")
  )

  expect_error(
    VariableFeatures(sce),
    "No Bioconductor-native HVG statistics found"
  )
  expect_length(VariableFeatures(sce, method = "seurat"), 10)
})
