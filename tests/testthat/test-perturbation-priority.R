## `RunPerturbationPriority()` with `group_col = "colLabels"`.
##
## The special value materialises `colLabels(sce)` into a temporary colData
## column and points `group_col` at it. The column is written to `sce` AFTER
## `colData(sce)` has already been read into `cd`, so the "does this column
## exist" check that follows was reading the pre-write snapshot and refused the
## column the branch had just created. Measured 2026-10-21: the call failed with
## "A valid `group_col` must be provided that exists in colData(sce)." on every
## object, and no test covered the branch.
##
## These tests drive the real function rather than a mocked one, because the bug
## is precisely in the interaction between the two. They are skipped when Augur
## is absent (GitHub-only), so the suite stays runnable without it.

skip_if_not_installed("Augur")

# Augur requires at least 20 cells per cell type PER CONDITION, so a fixture has
# to be a balanced 2x2 grid of (cell type x condition), not one cell type that
# happens to sit in two conditions.
fixture <- function(n_per_cell = 20L, n_genes = 100L) {
    grid <- expand.grid(ct = c("a", "b"), cond = c("x", "y"), stringsAsFactors = FALSE)
    grid <- grid[rep(seq_len(nrow(grid)), each = n_per_cell), , drop = FALSE]
    n <- nrow(grid)
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(rpois(n_genes * n, 5), nrow = n_genes, ncol = n))
    )
    SingleCellExperiment::colLabels(sce) <- grid$ct
    SummarizedExperiment::colData(sce)$condition <- grid$cond
    sce
}

test_that("group_col = 'colLabels' is accepted and used as the grouping", {
    sce <- fixture()

    out <- RunPerturbationPriority(sce, condition_col = "condition",
        group_col = "colLabels")

    expect_s4_class(out, "SingleCellExperiment")
    expect_true(has_perturbation_priority(out))

    # The record must carry the resolved column, not the magic value: a reader
    # cannot resolve "colLabels" again from a saved object.
    rec <- get_perturbation_priority(out)
    expect_identical(rec$group_col, "sclet_priority_group_temp")

    # Both cluster levels must have been scored, which is only true if the
    # grouping actually reached Augur.
    expect_true(is.data.frame(rec$AUC))
    expect_setequal(as.character(rec$AUC$cell_type), c("a", "b"))
})

test_that("a real colData column still works and is not renamed", {
    sce <- fixture()
    SummarizedExperiment::colData(sce)$grp <- SingleCellExperiment::colLabels(sce)

    out <- RunPerturbationPriority(sce, condition_col = "condition", group_col = "grp")
    rec <- get_perturbation_priority(out)
    expect_identical(rec$group_col, "grp")
    expect_setequal(as.character(rec$AUC$cell_type), c("a", "b"))
})

test_that("group_col = NULL falls back to the active identity", {
    sce <- fixture()
    out <- RunPerturbationPriority(sce, condition_col = "condition", group_col = NULL)
    rec <- get_perturbation_priority(out)
    # ActiveIdent() resolves the clustering, which the fixture wrote into
    # colLabels; the resolved name is whatever the object calls it.
    expect_true(nzchar(rec$group_col))
    expect_setequal(as.character(rec$AUC$cell_type), c("a", "b"))
})

test_that("an unknown group_col is still refused", {
    sce <- fixture()
    expect_error(
        RunPerturbationPriority(sce, condition_col = "condition", group_col = "nope"),
        "group_col"
    )
})

test_that("a single-level condition_col is still refused", {
    sce <- fixture()
    SummarizedExperiment::colData(sce)$condition <- "only_one"
    expect_error(
        RunPerturbationPriority(sce, condition_col = "condition",
            group_col = "colLabels"),
        "at least two distinct condition"
    )
})

## A `SingleCellExperiment` assay is often a `DelayedMatrix` (TENxPBMCData and
## every HDF5Array-backed object hand one back), and Augur accepts neither a
## DelayedMatrix nor a SummarizedExperiment -- only Seurat, monocle, a sparse
## matrix, a numeric matrix or a numeric data frame. The wrapper used to pass the
## assay straight through, so `RunPerturbationPriority()` failed on exactly the
## objects the book demonstrates it on. Measured 2026-10-21: "input must be
## Seurat, monocle, sparse matrix, numeric matrix, or numeric data frame".
test_that("a DelayedMatrix assay is densified before reaching Augur", {
    skip_if_not_installed("DelayedArray")
    sce <- fixture()
    # Wrap the assay in a DelayedArray: same values, an unacceptable class.
    counts <- SummarizedExperiment::assay(sce, "counts")
    SummarizedExperiment::assay(sce, "counts") <-
        DelayedArray::DelayedArray(counts)
    expect_s4_class(SummarizedExperiment::assay(sce, "counts"), "DelayedArray")

    out <- RunPerturbationPriority(sce, condition_col = "condition",
        group_col = "colLabels")

    rec <- get_perturbation_priority(out)
    expect_setequal(as.character(rec$AUC$cell_type), c("a", "b"))
})
