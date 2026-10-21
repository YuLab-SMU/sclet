## `AIAuditAnalysisChain()` must accept the LOGICAL layer name.
##
## The package's vocabulary for an integration result is the layer, not the
## assay: `RunIntegration(method = "fastMNN", correct.all = TRUE)` registers a
## layer called `"corrected"` and backs it with an assay whose name comes from
## batchelor (`"reconstructed"`). The audit required the value to be present in
## `assayNames()`, so the documented call
##
##     AIAuditAnalysisChain(sce, features, raw_layer = "counts",
##                          integrated_layer = "corrected")
##
## failed with "An integrated assay is required for AIAuditAnalysisChain()" on
## an object that had just been integrated successfully. Measured 2026-10-21
## while running the book chapter, which is the only thing that exercised it.
##
## `sclet_ai_call()` is mocked through the documented
## `options(sclet.ai.call = ...)` hook, so no model is contacted.

skip_if_not_installed("batchelor")

fake_call <- function(task, context, ...) {
    list(answer = "ok", proposed_actions = list())
}

make_integrated <- function() {
    set.seed(11)
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(rpois(60 * 40, 5), nrow = 60, ncol = 40))
    )
    SummarizedExperiment::colData(sce)$batch <- rep(c("b1", "b2"), each = 20)
    sce <- sclet::NormalizeData(sce)
    sce <- sclet::FindVariableFeatures(sce, method = "scran")
    sclet::RunIntegration(sce, method = "fastMNN", batch = "batch",
        correct.all = TRUE)
}

test_that("the logical 'corrected' layer resolves to its backing assay", {
    sce <- make_integrated()

    # The layer exists, and the assay behind it is NOT called "corrected" --
    # that mismatch is the whole point of the test.
    expect_true("corrected" %in% sclet::Layers(sce))
    expect_false("corrected" %in% SummarizedExperiment::assayNames(sce))

    old <- options(sclet.ai.call = fake_call)
    on.exit(options(old), add = TRUE)

    audit <- AIAuditAnalysisChain(
        sce,
        features = head(rownames(sce), 3),
        raw_layer = "counts",
        integrated_layer = "corrected"
    )
    expect_identical(audit$answer, "ok")
})

test_that("an explicit assay name still works", {
    sce <- make_integrated()
    # Any assay other than the raw one exercises the did-not-need-resolving
    # branch; the name itself is not the point.
    expect_true("logcounts" %in% SummarizedExperiment::assayNames(sce))

    old <- options(sclet.ai.call = fake_call)
    on.exit(options(old), add = TRUE)

    audit <- AIAuditAnalysisChain(
        sce,
        features = head(rownames(sce), 3),
        raw_layer = "counts",
        integrated_layer = "logcounts"
    )
    expect_identical(audit$answer, "ok")
})

test_that("an object with no integration is still refused", {
    mat <- matrix(1, nrow = 3, ncol = 4,
        dimnames = list(paste0("g", 1:3), NULL))
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = mat))
    SummarizedExperiment::colData(sce)$group <- rep(c("a", "b"), each = 2)

    expect_error(
        AIAuditAnalysisChain(sce, features = rownames(sce),
            raw_layer = "counts", integrated_layer = "corrected"),
        "integrated assay"
    )
})
