test_that("RunBasicWorkflow runs the standard pipeline end-to-end with defaults", {
    set.seed(1)
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(rpois(50 * 40, lambda = 5), nrow = 50L, ncol = 40L,
            dimnames = list(paste0("gene_", seq_len(50L)), paste0("c", seq_len(40L)))))
    )
    result <- RunBasicWorkflow(sce, n_pcs = 5, verbose = FALSE)
    expect_true("PCA" %in% SingleCellExperiment::reducedDimNames(result))
    expect_true("UMAP" %in% SingleCellExperiment::reducedDimNames(result))
    expect_false(is.null(ActiveIdent(result)))
    expect_true("logcounts" %in% SummarizedExperiment::assayNames(result))
    expect_equal(ncol(result), ncol(sce))
    expect_equal(
        unname(SummarizedExperiment::assay(result, "counts")),
        unname(SummarizedExperiment::assay(sce, "counts"))
    )
    # the underlying functions register their own state; nothing is duplicated
    records <- sclet_get_state(result)$states$records %||% list()
    expect_true("pca" %in% names(records$reduction %||% list()))
    expect_true("louvain_clusters" %in% names(records$clustering %||% list()))
})

test_that("RunBasicWorkflow respects a restricted steps argument", {
    set.seed(1)
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(rpois(50 * 40, lambda = 5), nrow = 50L, ncol = 40L,
            dimnames = list(paste0("gene_", seq_len(50L)), paste0("c", seq_len(40L)))))
    )
    result <- RunBasicWorkflow(sce, steps = c("normalize", "variable_features"), verbose = FALSE)
    expect_true("logcounts" %in% SummarizedExperiment::assayNames(result))
    expect_false("PCA" %in% SingleCellExperiment::reducedDimNames(result))
    expect_false("UMAP" %in% SingleCellExperiment::reducedDimNames(result))
    expect_null(ActiveIdent(result))
})

test_that("RunBasicWorkflow rejects unknown steps and steps that need a missing reduction", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(rpois(50 * 20, lambda = 5), nrow = 50L, ncol = 20L))
    )
    expect_error(RunBasicWorkflow(sce, steps = "not_a_step"), "unknown step")
    expect_error(RunBasicWorkflow(sce, steps = "clusters", verbose = FALSE), "FindNeighbors")
    expect_error(RunBasicWorkflow(sce, n_pcs = 0), "n_pcs")
})

test_that("RunBasicWorkflow is not registered as an AI execution action", {
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1, 5L, 4L)))
    registry <- AIDefaultExecutionRegistry(sce, include = "all")
    expect_false("RunBasicWorkflow" %in% names(registry))
})

test_that("AskAI provides a bounded beginner question interface", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 4, ncol = 3))
    )
    old <- options(sclet.ai.call = function(task, context, ...) {
        list(answer = "The object contains counts and has not been normalized yet.")
    })
    on.exit(options(old), add = TRUE)

    result <- AskAI(sce, "What should I do first?")
    expect_s3_class(result, "sclet_ai_result")
    expect_equal(result$task, "copilot")
    expect_match(result$answer, "counts")
    expect_false(isTRUE(result$context$dataset$data$expression$counts$values_included))
})

test_that("RunAIAnalysis hides plan execution details behind a beginner workflow", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 4, ncol = 3))
    )
    old <- options(sclet.ai.call = function(task, context, ...) {
        list(
            answer = "Inspect the current status before proposing computation.",
            proposed_actions = list(list(
                id = "status",
                action = "inspect_status",
                params = list()
            ))
        )
    })
    on.exit(options(old), add = TRUE)

    dry <- RunAIAnalysis(
        sce,
        goal = "Help me understand the current object.",
        confirm = "no"
    )
    expect_s3_class(dry, "sclet_ai_analysis")
    expect_equal(dry$status, "dry_run")
    expect_equal(dry$preview$status, "dry_run")
    expect_equal(nrow(dry$object), nrow(sce))

    done <- RunAIAnalysis(
        sce,
        goal = "Help me understand the current object.",
        confirm = TRUE
    )
    expect_s3_class(done, "sclet_ai_analysis")
    expect_equal(done$status, "completed")
    expect_equal(done$execution$status, "completed")
    expect_true(any(grepl("ai_plan_", names(GetAnalysisLedger(done$object)$analyses))))
})

test_that("RunAIAnalysis returns validation errors instead of executing an invalid plan", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 4, ncol = 3))
    )
    old <- options(sclet.ai.call = function(task, context, ...) {
        list(
            answer = "Use an unavailable action.",
            proposed_actions = list(list(id = "bad", action = "not_registered"))
        )
    })
    on.exit(options(old), add = TRUE)

    result <- RunAIAnalysis(sce, goal = "Do something useful", confirm = TRUE)
    expect_equal(result$status, "invalid_plan")
    expect_null(result$execution)
    expect_true(any(grepl("not registered", result$validation$errors, fixed = TRUE)))
    # the report layer gained a structured clarification block without changing status
    expect_equal(result$report$status, "invalid_plan")
    expect_equal(result$report$clarification$status, "clarification_required")
    expect_true(result$report$clarification$n_questions >= 1L)
    expect_equal(result$report$errors, result$validation$errors)
})
