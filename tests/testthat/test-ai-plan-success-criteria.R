test_that("new_sclet_ai_plan defaults are backward compatible", {
    plan <- sclet:::new_sclet_ai_plan(task = "test_plan")
    expect_identical(plan$success_criteria, list())
    expect_identical(plan$candidate_routes, list())
    expect_null(plan$selected_route)
    expect_identical(plan$assumptions, list())
    expect_identical(plan$risks, list())
    expect_s3_class(plan, "sclet_ai_plan")
})

test_that("success criterion referencing unknown step_id fails validation", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 5, ncol = 4))
    )
    registry <- AIDefaultExecutionRegistry(sce, include = "read")
    plan <- sclet:::new_sclet_ai_plan(
        task = "bad_criterion",
        actions = list(list(id = "step1", action = "inspect_status")),
        success_criteria = list(list(
            id = "c1",
            description = "check nonexistent step",
            source = "action_output",
            step_id = "nonexistent_step",
            check = "exists"
        ))
    )
    v <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_false(v$valid)
    expect_true(any(grepl("success_criterion_unresolvable", v$errors)))
})

test_that("success criterion with check gte but no value fails validation", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 5, ncol = 4))
    )
    registry <- AIDefaultExecutionRegistry(sce, include = "read")
    plan <- sclet:::new_sclet_ai_plan(
        task = "incomplete_criterion",
        actions = list(list(id = "step1", action = "inspect_status")),
        success_criteria = list(list(
            id = "c1",
            description = "needs value",
            source = "action_output",
            step_id = "step1",
            check = "gte",
            field = "n_cells"
        ))
    )
    v <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_false(v$valid)
    expect_true(any(grepl("success_criterion_incomplete", v$errors)))
})

test_that("selected_route not in candidate_routes fails validation", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 5, ncol = 4))
    )
    registry <- AIDefaultExecutionRegistry(sce, include = "read")
    plan <- sclet:::new_sclet_ai_plan(
        task = "bad_route",
        actions = list(list(id = "step1", action = "inspect_status")),
        candidate_routes = list(list(id = "harmony", description = "Harmony integration")),
        selected_route = "x"
    )
    v <- ValidateAIPlan(plan, registry = registry)
    expect_false(v$valid)
    expect_true(any(grepl("selected_route_unknown", v$errors)))
})

test_that("valid plan with matching route and criteria passes validation", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 5, ncol = 4))
    )
    registry <- AIDefaultExecutionRegistry(sce, include = "read")
    plan <- sclet:::new_sclet_ai_plan(
        task = "good_plan",
        actions = list(list(id = "step1", action = "inspect_status")),
        candidate_routes = list(list(id = "harmony", description = "Harmony integration")),
        selected_route = "harmony",
        success_criteria = list(list(
            id = "c1",
            description = "step output exists",
            source = "action_output",
            step_id = "step1",
            check = "exists"
        )),
        assumptions = list("batch is technical"),
        risks = list("may over-correct")
    )
    v <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_true(v$valid)
})

test_that("dry run returns empty success_assessment", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 5, ncol = 4))
    )
    registry <- AIDefaultExecutionRegistry(sce, include = "read")
    plan <- sclet:::new_sclet_ai_plan(
        task = "dry_plan",
        actions = list(list(id = "step1", action = "inspect_status")),
        success_criteria = list(list(
            id = "c1",
            description = "check output",
            source = "action_output",
            step_id = "step1",
            check = "exists"
        ))
    )
    v <- ValidateAIPlan(plan, object = sce, registry = registry)
    exec <- ExecuteAIPlan(sce, plan, registry, validation = v, dry_run = TRUE)
    expect_identical(exec$success_assessment, list())
})

test_that("executed plan with existing output criterion returns met", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(rpois(20, 5), nrow = 5, ncol = 4))
    )
    test_action <- AIAction(
        name = "test_sce_action",
        handler = function(object, params) {
            object
        },
        returns = "sce",
        requires_confirmation = FALSE,
        mutates_object = FALSE,
        idempotent = TRUE
    )
    registry <- AIExecutionRegistry(list(test_action))
    plan <- sclet:::new_sclet_ai_plan(
        task = "met_criterion",
        actions = list(list(id = "step1", action = "test_sce_action")),
        success_criteria = list(list(
            id = "c1",
            description = "n_cells exists",
            source = "action_output",
            step_id = "step1",
            check = "exists",
            field = "n_cells"
        ))
    )
    v <- ValidateAIPlan(plan, registry = registry)
    expect_true(v$valid)
    exec <- ExecuteAIPlan(sce, plan, registry, validation = v, dry_run = FALSE, record = FALSE)
    expect_equal(exec$status, "completed")
    expect_length(exec$success_assessment, 1L)
    expect_equal(exec$success_assessment[[1L]]$status, "met")
})

test_that("failed step produces not_available criterion and status remains failed", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(rpois(20, 5), nrow = 5, ncol = 4))
    )
    failing_action <- AIAction(
        name = "test_failing_action",
        handler = function(object, params) {
            stop("deliberate test failure")
        },
        returns = "value",
        requires_confirmation = FALSE,
        mutates_object = FALSE,
        idempotent = FALSE
    )
    registry <- AIExecutionRegistry(list(failing_action))
    plan <- sclet:::new_sclet_ai_plan(
        task = "fail_criterion",
        actions = list(list(id = "step1", action = "test_failing_action")),
        success_criteria = list(list(
            id = "c1",
            description = "check output of failed step",
            source = "action_output",
            step_id = "step1",
            check = "exists"
        ))
    )
    v <- ValidateAIPlan(plan, registry = registry)
    expect_true(v$valid)
    exec <- ExecuteAIPlan(sce, plan, registry, validation = v, dry_run = FALSE, record = FALSE)
    expect_equal(exec$status, "failed")
    expect_length(exec$success_assessment, 1L)
    expect_equal(exec$success_assessment[[1L]]$status, "not_available")
})

test_that("evidence criterion with unwritten evidence_id returns not_available", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(rpois(20, 5), nrow = 5, ncol = 4))
    )
    test_action <- AIAction(
        name = "test_no_evidence_action",
        handler = function(object, params) {
            list(done = TRUE)
        },
        returns = "value",
        requires_confirmation = FALSE,
        mutates_object = FALSE,
        idempotent = TRUE
    )
    registry <- AIExecutionRegistry(list(test_action))
    plan <- sclet:::new_sclet_ai_plan(
        task = "evidence_criterion",
        actions = list(list(id = "step1", action = "test_no_evidence_action")),
        success_criteria = list(list(
            id = "c1",
            description = "evidence that was never written",
            source = "evidence",
            evidence_id = "nonexistent_evidence_id",
            check = "exists"
        ))
    )
    v <- ValidateAIPlan(plan, registry = registry)
    expect_true(v$valid)
    exec <- ExecuteAIPlan(sce, plan, registry, validation = v, dry_run = FALSE, record = FALSE)
    expect_equal(exec$status, "completed")
    expect_length(exec$success_assessment, 1L)
    expect_equal(exec$success_assessment[[1L]]$status, "not_available")
    expect_true(grepl("no evidence node found", exec$success_assessment[[1L]]$reason))
})

test_that("success criterion with check equals met when values match", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(rpois(20, 5), nrow = 5, ncol = 4))
    )
    test_action <- AIAction(
        name = "test_equals_action",
        handler = function(object, params) {
            object
        },
        returns = "sce",
        requires_confirmation = FALSE,
        mutates_object = FALSE,
        idempotent = TRUE
    )
    registry <- AIExecutionRegistry(list(test_action))
    plan <- sclet:::new_sclet_ai_plan(
        task = "equals_criterion",
        actions = list(list(id = "step1", action = "test_equals_action")),
        success_criteria = list(list(
            id = "c1",
            description = "n_cells equals 4",
            source = "action_output",
            step_id = "step1",
            check = "equals",
            field = "n_cells",
            value = 4L
        ))
    )
    v <- ValidateAIPlan(plan, registry = registry)
    expect_true(v$valid)
    exec <- ExecuteAIPlan(sce, plan, registry, validation = v, dry_run = FALSE, record = FALSE)
    expect_equal(exec$status, "completed")
    expect_equal(exec$success_assessment[[1L]]$status, "met")
})

test_that("success criterion with check gte not_met when value is below threshold", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(rpois(20, 5), nrow = 5, ncol = 4))
    )
    test_action <- AIAction(
        name = "test_gte_action",
        handler = function(object, params) {
            object
        },
        returns = "sce",
        requires_confirmation = FALSE,
        mutates_object = FALSE,
        idempotent = TRUE
    )
    registry <- AIExecutionRegistry(list(test_action))
    plan <- sclet:::new_sclet_ai_plan(
        task = "gte_criterion",
        actions = list(list(id = "step1", action = "test_gte_action")),
        success_criteria = list(list(
            id = "c1",
            description = "n_cells >= 100",
            source = "action_output",
            step_id = "step1",
            check = "gte",
            field = "n_cells",
            value = 100
        ))
    )
    v <- ValidateAIPlan(plan, registry = registry)
    expect_true(v$valid)
    exec <- ExecuteAIPlan(sce, plan, registry, validation = v, dry_run = FALSE, record = FALSE)
    expect_equal(exec$status, "completed")
    expect_equal(exec$success_assessment[[1L]]$status, "not_met")
})

test_that("duplicate success criterion ids are rejected", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 5, ncol = 4))
    )
    registry <- AIDefaultExecutionRegistry(sce, include = "read")
    plan <- sclet:::new_sclet_ai_plan(
        task = "dup_criteria",
        actions = list(list(id = "step1", action = "inspect_status")),
        success_criteria = list(
            list(id = "c1", description = "first", source = "action_output", step_id = "step1", check = "exists"),
            list(id = "c1", description = "duplicate", source = "action_output", step_id = "step1", check = "exists")
        )
    )
    v <- ValidateAIPlan(plan, registry = registry)
    expect_false(v$valid)
    expect_true(any(grepl("success_criterion ids must be unique", v$errors)))
})

test_that("failed step not_available reason reports the step did not complete", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(rpois(20, 5), nrow = 5, ncol = 4))
    )
    failing_action <- AIAction(
        name = "test_failing_action_reason",
        handler = function(object, params) {
            stop("deliberate test failure")
        },
        returns = "value",
        requires_confirmation = FALSE,
        mutates_object = FALSE,
        idempotent = FALSE
    )
    registry <- AIExecutionRegistry(list(failing_action))
    plan <- sclet:::new_sclet_ai_plan(
        task = "fail_criterion_reason",
        actions = list(list(id = "step1", action = "test_failing_action_reason")),
        success_criteria = list(list(
            id = "c1",
            description = "check output of failed step",
            source = "action_output",
            step_id = "step1",
            check = "exists"
        ))
    )
    v <- ValidateAIPlan(plan, registry = registry)
    exec <- ExecuteAIPlan(sce, plan, registry, validation = v, dry_run = FALSE, record = FALSE)
    expect_equal(exec$success_assessment[[1L]]$status, "not_available")
    expect_true(grepl("did not complete", exec$success_assessment[[1L]]$reason))
})

test_that("success criterion with check in is met and not_met correctly", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(rpois(20, 5), nrow = 5, ncol = 4))
    )
    test_action <- AIAction(
        name = "test_in_action",
        handler = function(object, params) object,
        returns = "sce",
        requires_confirmation = FALSE,
        mutates_object = FALSE,
        idempotent = TRUE
    )
    registry <- AIExecutionRegistry(list(test_action))
    plan_met <- sclet:::new_sclet_ai_plan(
        task = "in_criterion_met",
        actions = list(list(id = "step1", action = "test_in_action")),
        success_criteria = list(list(
            id = "c1",
            description = "n_cells is in expected set",
            source = "action_output",
            step_id = "step1",
            check = "in",
            field = "n_cells",
            value = c(4L, 5L, 6L)
        ))
    )
    v_met <- ValidateAIPlan(plan_met, registry = registry)
    expect_true(v_met$valid)
    exec_met <- ExecuteAIPlan(sce, plan_met, registry, validation = v_met, dry_run = FALSE, record = FALSE)
    expect_equal(exec_met$success_assessment[[1L]]$status, "met")

    plan_not_met <- sclet:::new_sclet_ai_plan(
        task = "in_criterion_not_met",
        actions = list(list(id = "step1", action = "test_in_action")),
        success_criteria = list(list(
            id = "c1",
            description = "n_cells is in an unrelated set",
            source = "action_output",
            step_id = "step1",
            check = "in",
            field = "n_cells",
            value = c(100L, 200L)
        ))
    )
    v_not_met <- ValidateAIPlan(plan_not_met, registry = registry)
    exec_not_met <- ExecuteAIPlan(sce, plan_not_met, registry, validation = v_not_met, dry_run = FALSE, record = FALSE)
    expect_equal(exec_not_met$success_assessment[[1L]]$status, "not_met")
})

test_that("NA-containing comparisons degrade to not_available instead of erroring", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(rpois(20, 5), nrow = 5, ncol = 4))
    )
    na_action <- AIAction(
        name = "test_na_action",
        handler = function(object, params) {
            list(metric = NA_real_)
        },
        returns = "value",
        requires_confirmation = FALSE,
        mutates_object = FALSE,
        idempotent = TRUE
    )
    registry <- AIExecutionRegistry(list(na_action))
    plan <- sclet:::new_sclet_ai_plan(
        task = "na_criterion",
        actions = list(list(id = "step1", action = "test_na_action")),
        success_criteria = list(list(
            id = "c1",
            description = "metric gte threshold, but metric is NA",
            source = "action_output",
            step_id = "step1",
            check = "gte",
            field = "metric",
            value = 1
        ))
    )
    v <- ValidateAIPlan(plan, registry = registry)
    expect_true(v$valid)
    exec <- expect_no_error(
        ExecuteAIPlan(sce, plan, registry, validation = v, dry_run = FALSE, record = FALSE)
    )
    expect_equal(exec$status, "completed")
    expect_equal(exec$success_assessment[[1L]]$status, "not_available")
})

test_that("evidence criterion with empty evidence_id is rejected", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 5, ncol = 4))
    )
    registry <- AIDefaultExecutionRegistry(sce, include = "read")
    plan <- sclet:::new_sclet_ai_plan(
        task = "empty_evidence",
        actions = list(list(id = "step1", action = "inspect_status")),
        success_criteria = list(list(
            id = "c1",
            description = "empty evidence ref",
            source = "evidence",
            evidence_id = "",
            check = "exists"
        ))
    )
    v <- ValidateAIPlan(plan, registry = registry)
    expect_false(v$valid)
    expect_true(any(grepl("success_criterion_unresolvable.*empty evidence_id", v$errors)))
})
