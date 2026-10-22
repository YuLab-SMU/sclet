## The README quickstart, step 1-2.
##
## The quickstart is the first thing a user runs and none of it was covered:
## `Read10X()`, `PercentageFeatureSet()`, `subset_feature()` and `QCMetrics()`
## appeared in zero test files (measured 2026-10-22, 105 of 268 own exports
## untested). A break in any of them is a break on the user's first command,
## and nothing would have said so.
##
## These tests run the steps in the order the README runs them, on inputs small
## enough to be fast: a real 10x directory for the reader, a mock object for
## the rest.

# A minimal CellRanger-style `filtered_feature_bc_matrix/`: three files that
# `DropletUtils::read10xCounts()` reads. Small enough to write inline, real
# enough that the parser is actually exercised.
write_tiny_10x <- function(dir, n_genes = 20L, n_cells = 5L, with_mito = TRUE) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    set.seed(42)
    counts <- matrix(
        stats::rpois(n_genes * n_cells, lambda = 5),
        nrow = n_genes, ncol = n_cells
    )
    genes <- paste0("GENE", seq_len(n_genes))
    if (with_mito && n_genes >= 4L) {
        genes[1:4] <- paste0("MT-", c("ND1", "ND2", "CO1", "ATP6"))
    }
    barcodes <- paste0("CELL", seq_len(n_cells), "-1")
    Matrix::writeMM(Matrix::Matrix(counts, sparse = TRUE),
                    file.path(dir, "matrix.mtx"))
    writeLines(paste0(genes, "\t", genes), file.path(dir, "features.tsv"))
    writeLines(barcodes, file.path(dir, "barcodes.tsv"))
    dir
}

test_that("Read10X reads a 10x directory and computes the QC columns", {
    skip_if_not_installed("DropletUtils")
    dir <- write_tiny_10x(file.path(tempdir(), "tiny10x"))
    on.exit(unlink(dir, recursive = TRUE), add = TRUE)

    sce <- Read10X(dir)

    expect_s4_class(sce, "SingleCellExperiment")
    expect_equal(ncol(sce), 5L)
    expect_equal(nrow(sce), 20L)
    expect_true("counts" %in% SummarizedExperiment::assayNames(sce))
    # The README says Read10X also computes these; that claim is the thing
    # being tested, not just that a matrix came back.
    expect_true(all(c("nFeature_RNA", "nCount_RNA") %in%
                    colnames(SummarizedExperiment::colData(sce))))
    expect_equal(
        as.numeric(SummarizedExperiment::colData(sce)$nCount_RNA),
        as.numeric(Matrix::colSums(SummarizedExperiment::assay(sce, "counts")))
    )
})

test_that("Read10X refuses a directory that is not a 10x matrix", {
    skip_if_not_installed("DropletUtils")
    empty <- file.path(tempdir(), "not_a_10x_dir")
    dir.create(empty, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(empty, recursive = TRUE), add = TRUE)
    expect_error(suppressWarnings(Read10X(empty)))
})

test_that("PercentageFeatureSet reproduces the per-cell mitochondrial fraction", {
    skip_if_not_installed("scuttle")
    sce <- scuttle::mockSCE(ncells = 10, ngenes = 50)
    rownames(sce)[1:5] <- paste0("MT-", c("ND1", "ND2", "CO1", "CO2", "ATP6"))

    sce$percent.mt <- PercentageFeatureSet(sce, "^MT-")

    counts <- SummarizedExperiment::assay(sce, "counts")
    mito <- Matrix::colSums(counts[grepl("^MT-", rownames(sce)), , drop = FALSE])
    total <- Matrix::colSums(counts)
    expect_equal(as.numeric(sce$percent.mt), as.numeric(mito / total * 100))
    expect_true(all(sce$percent.mt >= 0 & sce$percent.mt <= 100))
})

test_that("PercentageFeatureSet on an object with no mitochondrial genes is zero", {
    skip_if_not_installed("scuttle")
    sce <- scuttle::mockSCE(ncells = 10, ngenes = 50)
    rownames(sce) <- paste0("GENE", seq_len(nrow(sce)))
    expect_true(all(PercentageFeatureSet(sce, "^MT-") == 0))
})

test_that("subset_feature drops genes seen in too few cells, as the README uses it", {
    skip_if_not_installed("scuttle")
    sce <- scuttle::mockSCE(ncells = 10, ngenes = 50)
    # Make five genes detectable in a single cell only, so mincell = 3 drops
    # them and keeps everything else.
    counts <- SummarizedExperiment::assay(sce, "counts")
    sparse_rows <- seq_len(5)
    counts[sparse_rows, ] <- 0
    counts[sparse_rows, 1] <- 10
    SummarizedExperiment::assay(sce, "counts") <- counts

    before <- nrow(sce)
    sce <- subset_feature(sce, mincell = 3, peek = FALSE)

    expect_lt(nrow(sce), before)
    detected <- Matrix::rowSums(counts > 0)
    expect_equal(nrow(sce), sum(detected >= 3))
    expect_false(any(rownames(counts)[sparse_rows] %in% rownames(sce)))
})

test_that("subset_feature refuses a non-SingleCellExperiment", {
    expect_error(subset_feature(matrix(1:4, 2), mincell = 1), "SingleCellExperiment")
})

test_that("QCMetrics writes the same two columns Read10X does", {
    skip_if_not_installed("scrapper")
    skip_if_not_installed("scuttle")
    sce <- scuttle::mockSCE(ncells = 10, ngenes = 50)
    rownames(sce) <- paste0("GENE", seq_len(nrow(sce)))

    sce <- QCMetrics(sce)

    expect_true(all(c("nFeature_RNA", "nCount_RNA") %in%
                    colnames(SummarizedExperiment::colData(sce))))
    counts <- SummarizedExperiment::assay(sce, "counts")
    expect_equal(as.numeric(sce$nFeature_RNA), as.numeric(Matrix::colSums(counts > 0)))
    expect_equal(as.numeric(sce$nCount_RNA), as.numeric(Matrix::colSums(counts)))
})

test_that("the quickstart steps run in order and leave a usable object", {
    # The whole point of the quickstart is that a reader can paste it and get
    # somewhere. This runs steps 1-3 on a mock object (a real Read10X is
    # covered above) and asserts the object is in the state the README claims
    # by step 3: a UMAP, cluster labels, and a Status() that reports them.
    skip_if_not_installed("scuttle")
    sce <- scuttle::mockSCE(ncells = 100, ngenes = 300)
    rownames(sce)[1:5] <- paste0("MT-", c("ND1", "ND2", "CO1", "CO2", "ATP6"))

    sce[["percent.mt"]] <- PercentageFeatureSet(sce, "^MT-")
    sce <- subset_feature(sce, mincell = 3, peek = FALSE)
    sce <- RunBasicWorkflow(sce, n_features = 50, n_pcs = 5,
                            cluster_resolution = 0.5, verbose = FALSE)

    expect_true(all(c("PCA", "UMAP") %in% SingleCellExperiment::reducedDimNames(sce)))
    expect_false(is.null(SingleCellExperiment::colLabels(sce)))
    status <- Status(sce)
    expect_true(is.list(status))
    expect_true("health" %in% names(status))
})
