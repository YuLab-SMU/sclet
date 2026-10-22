## `RunAIAnalysis()` must plan against the registry it will execute against.
##
## It builds that registry from `include`, passes it to `ValidateAIPlan()` and
## to `ExecuteAIPlan()` -- and used to omit it from the `AIPlanAnalysis()` call
## in between. `AIPlanAnalysis()` then defaulted to `include = "all"` (39
## actions), the planner legitimately proposed an action from a group outside
## `include`, and validation against the narrower registry refused it.
##
## Measured 2026-10-21 with the book's own example: goal "find markers that
## distinguish the clusters", default `include`. The planner proposed
## `run_de_test` -- which IS a registry action, in the `annotation` group -- and
## the result was `invalid_plan`, "action is not registered: run_de_test". The
## planner had been told the action existed. The failure looked like a model
## problem and was a wiring problem.
##
## `sclet_ai_call()` is mocked through `options(sclet.ai.call = ...)`, so no
## model is contacted.

# `run_de_test` lives in the annotation group; this asserts the premise both
# tests below depend on, and fails loudly if the groups are ever reorganised.
test_that("run_de_test is in the annotation group, not in the beginner default", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 4, ncol = 3))
    )
    beginner <- c("read", "preprocess", "dimred", "graph", "cluster")
    expect_false("run_de_test" %in%
        names(AIDefaultExecutionRegistry(sce, include = beginner)))
    expect_true("run_de_test" %in%
        names(AIDefaultExecutionRegistry(sce, include = "annotation")))
})

test_that("a plan naming an included action validates and runs", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(rpois(30 * 20, 5), nrow = 30, ncol = 20,
            dimnames = list(paste0("g", 1:30), NULL)))
    )
    SummarizedExperiment::colData(sce)$batch <- rep(c("a", "b"), each = 10)
    sce <- NormalizeData(sce)

    # The planner proposes an action from the annotation group. With
    # `include` naming that group, planning and execution see the same
    # catalog, so the plan is accepted.
    old <- options(sclet.ai.call = function(task, context, ...) {
        list(
            answer = "Inspect the recorded quality checks.",
            proposed_actions = list(list(
                id = "qc",
                action = "check_qc",
                params = list()
            ))
        )
    })
    on.exit(options(old), add = TRUE)

    result <- RunAIAnalysis(
        sce,
        goal = "check QC",
        include = c("read", "preprocess", "dimred", "graph", "cluster", "annotation"),
        confirm = TRUE
    )
    expect_s3_class(result, "sclet_ai_analysis")
    expect_false(identical(result$status, "invalid_plan"))
    # The registry the planner saw and the registry that executed are one object.
    expect_identical(result$registry, {
        AIDefaultExecutionRegistry(
            sce,
            include = c("read", "preprocess", "dimred", "graph", "cluster", "annotation")
        )
    })
})

test_that("planning sees only the groups it will execute against", {
    # The inverse of the bug: if the planner is NOT given a group, its context
    # must not advertise that group's actions, so it cannot propose them and be
    # refused later. This asserts the catalog the planner is shown.
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 4, ncol = 3))
    )
    seen <- NULL
    old <- options(sclet.ai.call = function(task, context, ...) {
        seen <<- vapply(context$available_actions %||% list(),
            function(a) a$name %||% "", character(1))
        list(answer = "ok", proposed_actions = list())
    })
    on.exit(options(old), add = TRUE)

    RunAIAnalysis(
        sce,
        goal = "check QC",
        include = c("read", "preprocess", "dimred", "graph", "cluster"),
        confirm = FALSE
    )
    expect_true("check_qc" %in% seen)
    expect_false("run_de_test" %in% seen)
})
