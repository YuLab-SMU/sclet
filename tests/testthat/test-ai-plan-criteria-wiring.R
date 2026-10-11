## Tests for end-to-end wiring of plan-level success_criteria and stop_conditions
## through AIPlanAnalysis() -> new_sclet_ai_plan() -> ValidateAIPlan() -> RunAIPlan().
##
## Background:
##   Before these fixes, three links were broken so model-produced plans could
##   never carry criteria end to end:
##     1. The native structured-output schema did not declare the fields,
##        so the model had nowhere to put them.
##     2. new_sclet_ai_result() / validate_sclet_ai_result() dropped them.
##     3. sclet_ai_plan_once() never forwarded result$success_criteria or
##        result$stop_conditions into new_sclet_ai_plan().
##   Each test below targets one of those links so any future regression that
##   re-breaks a link will fail exactly one named test.
##
##   Deliberately NOT tested here: whether a criterion's `field` is meetable.
##   Validation does not judge that, by design -- an unfindable field is
##   reported at execution as `not_available`, and a criterion may legitimately
##   name any key of the action's declared output_schema (see
##   test-ai-velocity-execution.R, where a stop condition names
##   `required_states`). An earlier attempt to reject unmeetable fields at
##   validation broke those tests and was reverted.

test_that("mock-supplied success_criteria and stop_conditions reach the plan (link 3)", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(rpois(20, 5), nrow = 5, ncol = 4))
    )
    old <- options(sclet.ai.call = function(task, context, ...) {
        list(
            answer = "Plan with criteria.",
            proposed_actions = list(list(
                id = "step1",
                action = "test_sce_action_wired",
                params = list(),
                depends_on = character()
            )),
            success_criteria = list(list(
                id = "c1",
                description = "n_cells exists in step1 output",
                source = "action_output",
                step_id = "step1",
                check = "exists",
                field = "n_cells"
            )),
            stop_conditions = list(list(
                id = "s1",
                description = "stop if step1 output fails",
                source = "action_output",
                step_id = "step1",
                check = "exists",
                field = "n_cells"
            ))
        )
    })
    on.exit(options(old), add = TRUE)

    # Register an action whose output_summary actually carries referenceable
    # fields (n_cells, n_features, ...), which Status(object) / inspect_status
    # do not: their summaries collapse to list(class=..., length=...).
    test_action <- AIAction(
        name = "test_sce_action_wired",
        handler = function(object, params) object,
        returns = "sce",
        requires_confirmation = FALSE,
        mutates_object = FALSE,
        idempotent = TRUE
    )
    registry <- AIExecutionRegistry(list(test_action))
    plan <- AIPlanAnalysis(
        sce,
        model = NULL,
        structured_output = FALSE,
        registry = registry,
        max_repair_rounds = 0L
    )
    expect_s3_class(plan, "sclet_ai_plan")
    # The schema/result/forwarding links are all wired when these hold.
    expect_equal(length(plan$success_criteria), 1L)
    expect_equal(length(plan$stop_conditions), 1L)
    expect_equal(plan$success_criteria[[1L]]$id, "c1")
    expect_equal(plan$stop_conditions[[1L]]$id, "s1")
    expect_true(isTRUE(plan$validation$valid))
    # RunAIPlan exercises the assessment path that used to be unreachable.
    exec <- RunAIPlan(sce, plan, registry, dry_run = FALSE, record = FALSE)
    expect_equal(exec$status, "completed")
    expect_length(exec$success_assessment, 1L)
    expect_equal(exec$success_assessment[[1L]]$status, "met")
})

test_that("schema link declares success_criteria and stop_conditions (link 1)", {
    skip_if_not_installed("aisdk")
    schema <- sclet:::sclet_ai_native_result_schema()
    expect_s3_class(schema, "z_schema")
    expect_true("success_criteria" %in% names(schema$properties))
    expect_true("stop_conditions" %in% names(schema$properties))
})

test_that("result link carries both fields through new_sclet_ai_result (link 2)", {
    result <- sclet:::new_sclet_ai_result(
        task = "analysis_plan",
        answer = "ok",
        success_criteria = list(list(
            id = "c1", description = "x",
            source = "action_output", step_id = "step1", check = "exists"
        )),
        stop_conditions = list(list(
            id = "s1", description = "y",
            source = "action_output", step_id = "step1", check = "exists"
        ))
    )
    expect_s3_class(result, "sclet_ai_result")
    expect_equal(length(result$success_criteria), 1L)
    expect_equal(length(result$stop_conditions), 1L)
    # The same shape must hold when validate_sclet_ai_result() is the entry
    # point (it is used by RecordAIResult directly, bypassing the constructor).
    res2 <- list(
        task = "analysis_plan",
        findings = list(),
        evidence = list(),
        warnings = list(),
        recommendations = list(),
        proposed_actions = list(),
        success_criteria = list(list(id = "c1")),
        stop_conditions = list(list(id = "s1"))
    )
    class(res2) <- c("sclet_ai_result", "list")
    expect_true(isTRUE(sclet:::validate_sclet_ai_result(res2)))
})

test_that("direct-construction path still works (no regression on existing tests)", {
    # Regression guard for the file-level tests already covering
    # new_sclet_ai_plan() / ValidateAIPlan() / ExecuteAIPlan(): the same
    # construction shape must keep passing after the new fields and the new
    # unmeetable rule land.
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(rpois(20, 5), nrow = 5, ncol = 4))
    )
    test_action <- AIAction(
        name = "test_direct_action",
        handler = function(object, params) object,
        returns = "sce",
        requires_confirmation = FALSE,
        mutates_object = FALSE,
        idempotent = TRUE
    )
    registry <- AIExecutionRegistry(list(test_action))
    plan <- sclet:::new_sclet_ai_plan(
        task = "direct_construction",
        actions = list(list(id = "step1", action = "test_direct_action")),
        success_criteria = list(list(
            id = "c1",
            description = "n_cells exists",
            source = "action_output",
            step_id = "step1",
            check = "exists",
            field = "n_cells"
        )),
        stop_conditions = list(list(
            id = "s1",
            description = "step output exists",
            source = "action_output",
            step_id = "step1",
            check = "exists"
        ))
    )
    v <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_true(v$valid)
    exec <- ExecuteAIPlan(sce, plan, registry, validation = v, dry_run = FALSE, record = FALSE)
    expect_equal(exec$status, "completed")
    expect_equal(exec$success_assessment[[1L]]$status, "met")
})

test_that("mock supplying success_criteria alone still validates and assesses", {
    # Smaller wire: the schema/result/forwarding links hold even when the
    # mock only sets success_criteria and omits stop_conditions entirely.
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(rpois(20, 5), nrow = 5, ncol = 4))
    )
    old <- options(sclet.ai.call = function(task, context, ...) {
        list(
            answer = "Plan with one criterion.",
            proposed_actions = list(list(
                id = "step1",
                action = "test_sce_action_partial",
                params = list(),
                depends_on = character()
            )),
            success_criteria = list(list(
                id = "c1",
                description = "n_cells exists in step1 output",
                source = "action_output",
                step_id = "step1",
                check = "exists",
                field = "n_cells"
            ))
        )
    })
    on.exit(options(old), add = TRUE)
    test_action <- AIAction(
        name = "test_sce_action_partial",
        handler = function(object, params) object,
        returns = "sce",
        requires_confirmation = FALSE,
        mutates_object = FALSE,
        idempotent = TRUE
    )
    registry <- AIExecutionRegistry(list(test_action))
    plan <- AIPlanAnalysis(
        sce,
        model = NULL,
        structured_output = FALSE,
        registry = registry,
        max_repair_rounds = 0L
    )
    expect_equal(length(plan$success_criteria), 1L)
    expect_equal(length(plan$stop_conditions), 0L)
    expect_true(isTRUE(plan$validation$valid))
    exec <- RunAIPlan(sce, plan, registry, dry_run = FALSE, record = FALSE)
    expect_equal(exec$status, "completed")
    expect_length(exec$success_assessment, 1L)
    expect_equal(exec$success_assessment[[1L]]$status, "met")
})
