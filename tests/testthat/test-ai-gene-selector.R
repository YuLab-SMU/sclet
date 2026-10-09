# Gene selectors must resolve, or fail with a message that names the problem.
#
# `features` (scale_data, run_integration) and `subset_row` (run_pca) are
# declared `character_vector`, so a string passes validation. Model output
# expresses selections as R range notation in a string -- "1:20" -- which is
# unambiguous but not a gene name. Passed straight through, the action subsets
# by a name that does not exist and base R raises an opaque error, reproduced
# on the real pipeline:
#
#   ScaleData(object, features = "1:20", assay = "logcounts")
#   #> Error: subscript out of bounds
#
# That names neither the parameter nor the problem, and because it happens
# after validation has passed, the repair loop never sees it.

selector_test_sce <- function() {
    counts <- matrix(1:40, nrow = 10L, ncol = 4L,
        dimnames = list(paste0("gene_", 1:10), paste0("cell_", 1:4)))
    SingleCellExperiment::SingleCellExperiment(list(counts = counts))
}

test_that("range notation in a string resolves to row names", {
    sce <- selector_test_sce()
    sel <- sclet:::sclet_ai_gene_selector
    expect_identical(sel(sce, "1:3"), c("gene_1", "gene_2", "gene_3"))
    expect_identical(sel(sce, "1:2,5:6"), c("gene_1", "gene_2", "gene_5", "gene_6"))
    # The exact string that produced "subscript out of bounds".
    expect_length(sel(sce, "1:10"), 10L)
})

test_that("row names, positions and NULL all resolve", {
    sce <- selector_test_sce()
    sel <- sclet:::sclet_ai_gene_selector
    expect_identical(sel(sce, c("gene_2", "gene_4")), c("gene_2", "gene_4"))
    expect_identical(sel(sce, c(2L, 4L)), c("gene_2", "gene_4"))
    expect_identical(sel(sce, c(2, 4)), c("gene_2", "gene_4"))
    expect_null(sel(sce, NULL))
    expect_identical(sel(sce, list("gene_1", "gene_2")), c("gene_1", "gene_2"))
})

test_that("an unresolvable selector names the parameter and the accepted forms", {
    sce <- selector_test_sce()
    sel <- sclet:::sclet_ai_gene_selector
    err <- tryCatch(sel(sce, "banana"), error = function(e) conditionMessage(e))
    expect_true(grepl("features is not usable", err, fixed = TRUE))
    expect_true(grepl("banana", err, fixed = TRUE))
    expect_true(grepl("positions within 1..10", err, fixed = TRUE))
    expect_true(grepl("\"1:20\"", err, fixed = TRUE))

    # An out-of-range range is reported, not silently truncated.
    err2 <- tryCatch(sel(sce, "1:200"), error = function(e) conditionMessage(e))
    expect_true(grepl("1:200", err2, fixed = TRUE))

    # A partly-valid name list reports the offending entries.
    err3 <- tryCatch(sel(sce, c("gene_1", "nope")), error = function(e) conditionMessage(e))
    expect_true(grepl("nope", err3, fixed = TRUE))

    # The argument name is used, so the message points at the right parameter.
    err4 <- tryCatch(sel(sce, "nope", arg = "subset_row"), error = function(e) conditionMessage(e))
    expect_true(grepl("subset_row is not usable", err4, fixed = TRUE))
})

test_that("an object without row names is rejected clearly", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1:8, nrow = 4L))
    )
    err <- tryCatch(
        sclet:::sclet_ai_gene_selector(sce, "1:2"),
        error = function(e) conditionMessage(e)
    )
    expect_true(grepl("no row names", err, fixed = TRUE))
})

test_that("scale_data accepts the range notation that used to crash", {
    sce <- selector_test_sce()
    sce <- NormalizeData(sce)
    registry <- AIExecutionRegistry(list(
        AIDefaultExecutionRegistry(sce, include = "preprocess")[["scale_data"]]
    ))
    plan <- sclet:::new_sclet_ai_plan(
        task = "analysis_plan",
        actions = list(list(id = "s1", action = "scale_data", params = list(features = "1:3")))
    )
    validation <- ValidateAIPlan(plan, object = sce, registry = registry, strict = FALSE)
    expect_true(validation$valid)

    execution <- RunAIPlan(
        object = sce, plan = plan, registry = registry,
        dry_run = FALSE, confirm = TRUE
    )
    expect_equal(execution$status, "completed")
    expect_true("scaled" %in% SummarizedExperiment::assayNames(execution$object))
})

test_that("scale_data reports an unusable selector instead of crashing", {
    sce <- selector_test_sce()
    sce <- NormalizeData(sce)
    registry <- AIExecutionRegistry(list(
        AIDefaultExecutionRegistry(sce, include = "preprocess")[["scale_data"]]
    ))
    plan <- sclet:::new_sclet_ai_plan(
        task = "analysis_plan",
        actions = list(list(id = "s1", action = "scale_data", params = list(features = "banana")))
    )
    execution <- RunAIPlan(
        object = sce, plan = plan, registry = registry,
        dry_run = FALSE, confirm = TRUE
    )
    expect_equal(execution$status, "failed")
    # Not "subscript out of bounds".
    expect_true(grepl("features is not usable", execution$results[[1L]]$error, fixed = TRUE))
})
