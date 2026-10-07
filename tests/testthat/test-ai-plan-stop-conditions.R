test_that("new_sclet_ai_plan defaults include empty stop_conditions and preserve P0.1 defaults", {
    plan <- sclet:::new_sclet_ai_plan(task = "test_plan")
    expect_identical(plan$stop_conditions, list())
    expect_identical(plan$assumptions, list())
    expect_identical(plan$candidate_routes, list())
    expect_null(plan$selected_route)
    expect_identical(plan$success_criteria, list())
    expect_identical(plan$risks, list())
    expect_s3_class(plan, "sclet_ai_plan")
})

test_that("stop condition with action_output referencing unknown step_id fails validation", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 5, ncol = 4))
    )
    registry <- AIDefaultExecutionRegistry(sce, include = "read")
    plan <- sclet:::new_sclet_ai_plan(
        task = "bad_stop",
        actions = list(list(id = "step1", action = "inspect_status")),
        stop_conditions = list(list(
            id = "s1",
            description = "stop if nonexistent step fails",
            source = "action_output",
            step_id = "nonexistent_step",
            check = "exists"
        ))
    )
    v <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_false(v$valid)
    expect_true(any(grepl("stop_condition_unresolvable", v$errors)))
})

test_that("stop condition with evidence source and empty evidence_id fails validation", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 5, ncol = 4))
    )
    registry <- AIDefaultExecutionRegistry(sce, include = "read")
    plan <- sclet:::new_sclet_ai_plan(
        task = "empty_evidence_stop",
        actions = list(list(id = "step1", action = "inspect_status")),
        stop_conditions = list(list(
            id = "s1",
            description = "stop on empty evidence ref",
            source = "evidence",
            evidence_id = "",
            check = "exists"
        ))
    )
    v <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_false(v$valid)
    expect_true(any(grepl("stop_condition_unresolvable", v$errors)))
})

test_that("stop condition with check gte but no field/value fails validation", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 5, ncol = 4))
    )
    registry <- AIDefaultExecutionRegistry(sce, include = "read")
    plan <- sclet:::new_sclet_ai_plan(
        task = "incomplete_stop",
        actions = list(list(id = "step1", action = "inspect_status")),
        stop_conditions = list(list(
            id = "s1",
            description = "needs field and value",
            source = "action_output",
            step_id = "step1",
            check = "gte"
        ))
    )
    v <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_false(v$valid)
    expect_true(any(grepl("stop_condition_incomplete", v$errors)))
})

test_that("well-formed stop_conditions on a valid plan passes validation", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 5, ncol = 4))
    )
    registry <- AIDefaultExecutionRegistry(sce, include = "read")
    plan <- sclet:::new_sclet_ai_plan(
        task = "good_stop",
        actions = list(list(id = "step1", action = "inspect_status")),
        stop_conditions = list(list(
            id = "s1",
            description = "stop if step fails",
            source = "action_output",
            step_id = "step1",
            check = "exists"
        ))
    )
    v <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_true(v$valid)
})

test_that("run_pca ncomponents larger than object dimensions fails data_scale check", {
    set.seed(1)
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(rpois(40, 5), nrow = 10, ncol = 4))
    )
    rownames(sce) <- paste0("g", seq_len(10))
    colnames(sce) <- paste0("c", seq_len(4))
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 5)
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "dimred"))
    plan <- sclet:::new_sclet_ai_plan(
        task = "too_many_pcs",
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
        actions = list(list(
            id = "pca", action = "run_pca",
            params = list(ncomponents = 50)
        ))
    )
    v <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_false(v$valid)
    expect_true(any(grepl("data_scale_incompatible", v$errors)))
})

test_that("run_pca ncomponents within range passes data_scale check", {
    set.seed(1)
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(rpois(40, 5), nrow = 10, ncol = 4))
    )
    rownames(sce) <- paste0("g", seq_len(10))
    colnames(sce) <- paste0("c", seq_len(4))
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 5)
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "dimred"))
    plan <- sclet:::new_sclet_ai_plan(
        task = "ok_pcs",
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
        actions = list(list(
            id = "pca", action = "run_pca",
            params = list(ncomponents = 2)
        ))
    )
    v <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_true(v$valid)
    expect_false(any(grepl("data_scale_incompatible", v$errors)))
})

test_that("data_scale check does not fire for non-run_pca actions", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 5, ncol = 4))
    )
    registry <- AIDefaultExecutionRegistry(sce, include = "read")
    plan <- sclet:::new_sclet_ai_plan(
        task = "no_pca",
        actions = list(list(id = "step1", action = "inspect_status")),
        stop_conditions = list()
    )
    v <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_true(v$valid)
    expect_false(any(grepl("data_scale_incompatible", v$errors)))
})

test_that("human_confirmations includes confirmed design from analysis_story", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrix(1, 6L, 8L)),
        colData = S4Vectors::DataFrame(
            batch_col = rep(c("a", "b"), each = 4L)
        )
    )
    sce <- ConfirmAIDesignSemantics(sce, design = list(batch = "batch_col"))
    registry <- AIDefaultExecutionRegistry(sce, include = "read")
    plan <- sclet:::new_sclet_ai_plan(
        task = "confirmed_design",
        actions = list(list(id = "step1", action = "inspect_status"))
    )
    v <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_true(length(v$human_confirmations) >= 1L)
    confirmed_entries <- Filter(function(x) identical(x$status, "confirmed"), v$human_confirmations)
    expect_true(length(confirmed_entries) >= 1L)
    expect_equal(confirmed_entries[[1L]]$role, "batch")
    expect_equal(confirmed_entries[[1L]]$column, "batch_col")
})

test_that("human_confirmations includes required_not_confirmed for unconfirmed design", {
    set.seed(1)
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(rpois(40, 5), nrow = 10, ncol = 4))
    )
    SummarizedExperiment::colData(sce)$batch <- c("a", "a", "b", "b")
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "integration"))
    plan <- sclet:::new_sclet_ai_plan(
        task = "unconfirmed_plan",
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
        actions = list(list(
            id = "integrate", action = "run_integration",
            params = list(batch = "batch", method = "fastMNN")
        ))
    )
    v <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_false(v$valid)
    req_entries <- Filter(function(x) identical(x$status, "required_not_confirmed"), v$human_confirmations)
    expect_true(length(req_entries) >= 1L)
})

test_that("dry run surfaces success_criteria_preview and stop_conditions_preview", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 5, ncol = 4))
    )
    registry <- AIDefaultExecutionRegistry(sce, include = "read")
    plan <- sclet:::new_sclet_ai_plan(
        task = "dry_preview",
        actions = list(list(id = "step1", action = "inspect_status")),
        success_criteria = list(list(
            id = "c1",
            description = "check output",
            source = "action_output",
            step_id = "step1",
            check = "exists"
        )),
        stop_conditions = list(list(
            id = "s1",
            description = "stop if step fails",
            source = "action_output",
            step_id = "step1",
            check = "exists"
        ))
    )
    v <- ValidateAIPlan(plan, object = sce, registry = registry)
    exec <- ExecuteAIPlan(sce, plan, registry, validation = v, dry_run = TRUE)
    expect_identical(exec$success_criteria_preview, plan$success_criteria)
    expect_identical(exec$stop_conditions_preview, plan$stop_conditions)
    expect_identical(exec$success_assessment, list())
})

test_that("dry run threads through human_confirmations from validation", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrix(1, 6L, 8L)),
        colData = S4Vectors::DataFrame(
            batch_col = rep(c("a", "b"), each = 4L)
        )
    )
    sce <- ConfirmAIDesignSemantics(sce, design = list(batch = "batch_col"))
    registry <- AIDefaultExecutionRegistry(sce, include = "read")
    plan <- sclet:::new_sclet_ai_plan(
        task = "dry_confirm",
        actions = list(list(id = "step1", action = "inspect_status"))
    )
    v <- ValidateAIPlan(plan, object = sce, registry = registry)
    exec <- ExecuteAIPlan(sce, plan, registry, validation = v, dry_run = TRUE)
    expect_identical(exec$human_confirmations, v$human_confirmations)
})

test_that("existing success_criteria and execution tests still pass unchanged", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 5, ncol = 4))
    )
    registry <- AIDefaultExecutionRegistry(sce, include = "read")
    plan <- sclet:::new_sclet_ai_plan(
        task = "regression_guard",
        actions = list(list(id = "step1", action = "inspect_status")),
        success_criteria = list(list(
            id = "c1",
            description = "step output exists",
            source = "action_output",
            step_id = "step1",
            check = "exists"
        ))
    )
    v <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_true(v$valid)
    exec <- ExecuteAIPlan(sce, plan, registry, validation = v, dry_run = TRUE)
    expect_equal(exec$status, "dry_run")
    expect_identical(exec$success_assessment, list())
})
