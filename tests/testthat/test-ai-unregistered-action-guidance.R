## "action is not registered" must say which group holds the action.
##
## The message was correct but a dead end: when a plan names a real action that
## the registry simply did not request, the caller was told only that it was
## "not registered" and reasonably concluded the package could not do it.
## Measured 2026-10-21 running the book: `RunAIAnalysis(goal = "find markers")`
## proposed `run_de_test` (a real action, in the `annotation` group) and the only
## feedback was the bare "not registered" line.

test_that("an action that exists in another group is named as such", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 4, ncol = 3))
    )
    expect_setequal(sclet_ai_action_groups(sce, "run_de_test"), "annotation")
    expect_setequal(sclet_ai_action_groups(sce, "normalize_data"), "preprocess")

    plan <- new_sclet_ai_plan(task = "t", actions = list(
        list(id = "s1", action = "run_de_test")
    ))
    registry <- AIDefaultExecutionRegistry(sce,
        include = c("read", "preprocess", "dimred", "graph", "cluster"))
    validation <- ValidateAIPlan(plan, object = sce, registry = registry)

    expect_false(validation$valid)
    msg <- paste(unlist(validation$errors), collapse = " ")
    expect_match(msg, "not registered", fixed = TRUE)
    expect_match(msg, "annotation", fixed = TRUE)
    expect_match(msg, "include", fixed = TRUE)
})

test_that("an action that exists nowhere keeps the original, actionable text", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 4, ncol = 3))
    )
    expect_length(sclet_ai_action_groups(sce, "no_such_action"), 0L)

    plan <- new_sclet_ai_plan(task = "t", actions = list(
        list(id = "s1", action = "no_such_action")
    ))
    registry <- AIDefaultExecutionRegistry(sce, include = "read")
    validation <- ValidateAIPlan(plan, object = sce, registry = registry)

    expect_false(validation$valid)
    msg <- paste(unlist(validation$errors), collapse = " ")
    expect_match(msg, "exact registry key", fixed = TRUE)
    # The "it exists in group X" branch must not fire for a genuinely unknown
    # name -- pointing at a group that does not contain it would be worse than
    # saying nothing.
    expect_false(grepl("it exists, but in the", msg, fixed = TRUE))
})

test_that("an action already in the registry validates without that error", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 4, ncol = 3))
    )
    SummarizedExperiment::colData(sce)$batch <- c("a", "a", "b")
    plan <- new_sclet_ai_plan(task = "t", actions = list(
        list(id = "s1", action = "normalize_data")
    ))
    registry <- AIDefaultExecutionRegistry(sce, include = "preprocess")
    validation <- ValidateAIPlan(plan, object = sce, registry = registry)
    msg <- paste(unlist(validation$errors), collapse = " ")
    expect_false(grepl("action is not registered", msg, fixed = TRUE))
})
