test_that("default execution registry exposes only safe read-only actions", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 3, ncol = 2))
    )
    registry <- AIDefaultExecutionRegistry(sce)
    expect_s3_class(registry, "sclet_ai_execution_registry")
    expect_equal(names(registry), c("inspect_status", "inspect_ledger", "check_qc"))
    expect_true(all(!vapply(registry, function(x) isTRUE(x$mutates_object), logical(1))))
    expect_true(all(!vapply(registry, function(x) isTRUE(x$requires_confirmation), logical(1))))

    plan <- new_sclet_ai_plan(
        task = "read_only_plan",
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
        actions = list(
            list(id = "status", action = "inspect_status"),
            list(id = "qc", action = "check_qc", depends_on = "status")
        )
    )
    validation <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_true(isTRUE(validation$valid))
    result <- ExecuteAIPlan(sce, plan, registry, validation = validation, dry_run = FALSE,
        confirmation = validation$confirmation_token)
    expect_equal(result$status, "completed")
    expect_true(isTRUE(result$recorded))
    expect_equal(nrow(result$object), nrow(sce))
    expect_equal(ncol(result$object), ncol(sce))
})

test_that("action input schemas reject unknown and invalid parameters", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 3, ncol = 2))
    )
    registry <- AIDefaultExecutionRegistry(sce)
    plan <- new_sclet_ai_plan(
        task = "bad_params",
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
        actions = list(list(
            id = "ledger",
            action = "inspect_ledger",
            params = list(detail = "not-a-detail", unexpected = TRUE)
        ))
    )
    validation <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_false(isTRUE(validation$valid))
    expect_true(any(grepl("unknown parameter", validation$errors, fixed = TRUE)))
    expect_true(any(grepl("must be one of", validation$errors, fixed = TRUE)))
})
test_that("AIPlanAnalysis creates a non-executing structured plan", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 4, ncol = 3))
    )
    old <- options(sclet.ai.call = function(task, context, ...) {
        list(
            answer = "Run the registered marker action after checking counts.",
            proposed_actions = list(list(
                id = "markers",
                action = "mark_cells",
                params = list(label = "review"),
                depends_on = character(),
                expected_outputs = list(column = "ai_plan_marker")
            ))
        )
    })
    on.exit(options(old), add = TRUE)

    plan <- AIPlanAnalysis(sce)
    expect_s3_class(plan, "sclet_ai_plan")
    expect_length(plan$actions, 1)
    expect_equal(plan$actions[[1]]$action, "mark_cells")
    expect_true(isTRUE(plan$requires_confirmation))
})

test_that("validated plans execute only registered actions after confirmation", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 4, ncol = 3))
    )
    registry <- AIExecutionRegistry(list(
        mark_cells = AIAction(
            "mark_cells",
            function(object, params) {
                cd <- SummarizedExperiment::colData(object)
                cd$ai_plan_marker <- rep(params$label, ncol(object))
                SummarizedExperiment::colData(object) <- cd
                object
            },
            prerequisites = function(object, params) {
                "counts" %in% SummarizedExperiment::assayNames(object)
            }
        )
    ))
    plan <- new_sclet_ai_plan(
        task = "test_plan",
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
        actions = list(list(
            id = "markers",
            action = "mark_cells",
            params = list(label = "review")
        ))
    )
    validation <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_true(isTRUE(validation$valid))
    expect_true(nzchar(validation$confirmation_token))

    dry <- ExecuteAIPlan(sce, plan, registry, validation = validation, dry_run = TRUE)
    expect_equal(dry$status, "dry_run")
    expect_false("ai_plan_marker" %in% colnames(SummarizedExperiment::colData(dry$object)))

    expect_error(
        ExecuteAIPlan(sce, plan, registry, validation = validation, dry_run = FALSE),
        class = "sclet_ai_confirmation_required"
    )

    executed <- ExecuteAIPlan(
        sce,
        plan,
        registry,
        validation = validation,
        dry_run = FALSE,
        confirmation = validation$confirmation_token
    )
    expect_equal(executed$status, "completed")
    expect_true("ai_plan_marker" %in% colnames(SummarizedExperiment::colData(executed$object)))
    expect_true(isTRUE(executed$recorded))
    expect_true(grepl("^ai_execution_", executed$execution_id))
    expect_true(!is.null(sclet_get_analysis(executed$object, executed$execution_id)))
    expect_true(length(sclet_get_commands(executed$object)) >= 1)
})

test_that("failed actions are recorded and stop subsequent execution", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 2, ncol = 2))
    )
    calls <- 0L
    registry <- AIExecutionRegistry(list(
        fail_action = AIAction(
            "fail_action",
            function(object, params) {
                calls <<- calls + 1L
                stop("intentional action failure")
            }
        ),
        never_action = AIAction(
            "never_action",
            function(object, params) {
                calls <<- calls + 1L
                object
            }
        )
    ))
    plan <- new_sclet_ai_plan(
        task = "failure_plan",
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
        actions = list(
            list(id = "fail", action = "fail_action"),
            list(id = "never", action = "never_action", depends_on = "fail")
        )
    )
    validation <- ValidateAIPlan(plan, object = sce, registry = registry)
    result <- ExecuteAIPlan(
        sce, plan, registry, validation = validation, dry_run = FALSE,
        confirmation = validation$confirmation_token
    )
    expect_equal(result$status, "failed")
    expect_equal(calls, 1L)
    expect_equal(result$results[[1]]$status, "failed")
    expect_true(grepl("intentional action failure", result$results[[1]]$error, fixed = TRUE))
    expect_true(isTRUE(result$recorded))
    execution <- sclet_get_analysis(result$object, result$execution_id)
    expect_equal(execution$summary$status, "failed")
})

test_that("validation rejects a stale plan fingerprint", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 2, ncol = 2))
    )
    registry <- AIExecutionRegistry(list(
        no_op = AIAction("no_op", function(object, params) object)
    ))
    plan <- new_sclet_ai_plan(
        task = "stale_plan",
        context_fingerprint = "not-the-current-fingerprint",
        actions = list(list(id = "step", action = "no_op"))
    )
    validation <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_false(isTRUE(validation$valid))
    expect_true(any(grepl("fingerprint", validation$errors, fixed = TRUE)))
    expect_error(
        ValidateAIPlan(plan, object = sce, registry = registry, strict = TRUE),
        class = "sclet_ai_invalid_plan"
    )
})
test_that("validation rejects unregistered actions and unmet prerequisites", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(logcounts = matrix(1, nrow = 2, ncol = 2))
    )
    registry <- AIExecutionRegistry(list(
        needs_counts = AIAction(
            "needs_counts",
            function(object, params) object,
            prerequisites = function(object, params) {
                if (!"counts" %in% SummarizedExperiment::assayNames(object)) {
                    return("counts assay is required")
                }
                TRUE
            }
        )
    ))
    plan <- new_sclet_ai_plan(
        task = "invalid_plan",
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
        actions = list(
            list(id = "missing", action = "not_registered"),
            list(id = "needs", action = "needs_counts")
        )
    )
    validation <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_false(isTRUE(validation$valid))
    expect_true(any(grepl("not registered", validation$errors, fixed = TRUE)))
    expect_true(any(grepl("counts assay is required", validation$errors, fixed = TRUE)))
    expect_error(
        ExecuteAIPlan(sce, plan, registry, dry_run = FALSE),
        class = "sclet_ai_invalid_plan"
    )
})
