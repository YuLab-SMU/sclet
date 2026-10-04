test_that("AIInvestigate uses bounded profile and is read-only", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrix(1, nrow = 3L, ncol = 4L)),
        colData = S4Vectors::DataFrame(batch = rep(c("batch-x", "batch-y"), each = 2L))
    )
    before <- sclet:::GetAIProfile(sce)$fingerprint
    captured <- NULL
    old <- options(sclet.ai.call = function(task, context, ...) {
        captured <<- list(task = task, context = context)
        list(answer = "Need confirm design semantics.", warnings = list(), findings = list(), recommendations = list())
    })
    on.exit(options(old), add = TRUE)

    result <- sclet:::AIInvestigate(sce, "Should I correct batch effects?")

    expect_s3_class(result, "sclet_ai_result")
    expect_true(isTRUE(result$metadata$read_only))
    expect_false(isTRUE(result$metadata$actions_executed))
    expect_true(isTRUE(result$metadata$clarification_required))
    expect_equal(captured$task, "Should I correct batch effects?")
    expect_false(captured$context$execution$allowed)
    expect_false(isTRUE(captured$context$privacy$complete_matrix_included))
    expect_false(isTRUE(captured$context$privacy$metadata_values_included))
    expect_true(any(vapply(result$warnings, function(x) identical(x$code, "clarification_required"), logical(1L))))
    expect_identical(before, sclet:::GetAIProfile(sce)$fingerprint)
})

test_that("AIInvestigate validates questions and never needs an execution registry", {
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = matrix(1, nrow = 2L, ncol = 2L)))
    expect_error(sclet:::AIInvestigate(sce, " "), "question")
    expect_error(sclet:::AIInvestigate(list(), "question"), "SingleCellExperiment")
})
