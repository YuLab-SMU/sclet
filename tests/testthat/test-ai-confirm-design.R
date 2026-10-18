## Design roles as an action.
##
## `run_integration` refuses until ConfirmAIDesignSemantics() names the batch
## column, and no action could make that call: on full Kang (2026-10-18) every
## plan that asked to integrate was rejected with
## `design_semantics_not_confirmed`, repair rounds included.

test_that("confirm_design clears the integration design gate", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 10L, ncol = 8L,
            dimnames = list(paste0("g", 1:10), paste0("c", 1:8)))))
    SummarizedExperiment::colData(sce)$batch <- rep(c("b1", "b2"), each = 4L)
    SummarizedExperiment::colData(sce)$condition <- rep(c("ctrl", "stim"), 4L)
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "integration"))

    ## Before the design is confirmed the integration step is refused.
    integration <- list(id = "i1", action = "run_integration",
        params = list(method = "fastMNN", batch = "batch"))
    refusal <- registry$run_integration$prerequisites(sce, integration$params)
    expect_true(grepl("design_semantics_not_confirmed", refusal, fixed = TRUE))

    ## confirm_design is available and executes.
    expect_true("confirm_design" %in% names(registry))
    plan <- sclet:::new_sclet_ai_plan(task = "design",
        actions = list(list(id = "d1", action = "confirm_design",
            params = list(batch = "batch", condition = "condition"))))
    validation <- ValidateAIPlan(plan, object = sce, registry = registry, strict = TRUE)
    expect_true(validation$valid)
    executed <- ExecuteAIPlan(sce, plan, registry, validation = validation,
        dry_run = FALSE, confirmation = validation$confirmation_token)
    expect_equal(executed$status, "completed")

    ## The gate is now satisfied on the executed object.
    after <- registry$run_integration$prerequisites(executed$object,
        integration$params)
    expect_true(!grepl("design_semantics_not_confirmed", after, fixed = TRUE))
})

test_that("a plan that confirms the design can validate the integration step", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 10L, ncol = 8L,
            dimnames = list(paste0("g", 1:10), paste0("c", 1:8)))))
    SummarizedExperiment::colData(sce)$batch <- rep(c("b1", "b2"), each = 4L)
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "integration"))
    plan <- sclet:::new_sclet_ai_plan(
        task = "design then integrate",
        actions = list(
            list(id = "d1", action = "confirm_design", params = list(batch = "batch")),
            list(id = "i1", action = "run_integration",
                params = list(method = "fastMNN", batch = "batch"))
        ))
    validation <- ValidateAIPlan(plan, object = sce, registry = registry, strict = FALSE)
    expect_true(validation$valid)
    expect_false(any(grepl("design_semantics_not_confirmed", validation$errors)))
})

test_that("confirm_design refuses unknown columns and empty roles", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 5L, ncol = 4L)))
    SummarizedExperiment::colData(sce)$batch <- rep(c("b1", "b2"), each = 2L)
    registry <- AIDefaultExecutionRegistry(sce, include = "read")
    action <- registry$confirm_design

    none <- action$prerequisites(sce, list())
    expect_true(grepl("design_role_missing", none, fixed = TRUE))

    bad <- action$prerequisites(sce, list(batch = "sample_id"))
    expect_true(grepl("design_role_invalid", bad, fixed = TRUE))
})

test_that("a plan may reference step.output.reduction after integration", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 10L, ncol = 8L,
            dimnames = list(paste0("g", 1:10), paste0("c", 1:8)))))
    SummarizedExperiment::colData(sce)$batch <- rep(c("b1", "b2"), each = 4L)
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "integration"))
    plan <- sclet:::new_sclet_ai_plan(
        task = "integrate then neighbours",
        actions = list(
            list(id = "d1", action = "confirm_design", params = list(batch = "batch")),
            list(id = "i1", action = "run_integration",
                params = list(method = "fastMNN", batch = "batch")),
            list(id = "n1", action = "find_neighbors",
                params = list(reduction = "${i1.output.reduction}", dims = 1:5),
                depends_on = list("i1"))
        ))
    validation <- ValidateAIPlan(plan, object = sce, registry = registry, strict = FALSE)
    expect_false(any(grepl("value_reference_unresolvable", validation$errors)))
})
