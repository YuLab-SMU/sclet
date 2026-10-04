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

test_that("sclet_ai_format_clarification maps a design-confirmation error to a structured question", {
    errors <- c(
        "prerequisites not met for action integrate: design_semantics_not_confirmed: call ConfirmAIDesignSemantics(object, design = list(batch = 'batch', ...)) before running integration"
    )
    result <- sclet:::sclet_ai_format_clarification(errors)
    expect_equal(result$status, "clarification_required")
    expect_true(length(result$questions) >= 1L)
    expect_equal(result$raw_errors, errors)
    question <- result$questions[[1L]]
    expect_equal(question$id, "design_batch")
    expect_true(question$recognized)
    expect_equal(question$blocked_action, "run_integration")
    expect_equal(question$related_function, "ConfirmAIDesignSemantics")
    expect_true(grepl("ConfirmAIDesignSemantics", question$text, fixed = TRUE))
    # the formatter must never answer the question or suggest a column
    expect_false(grepl("'batch'", question$text, fixed = TRUE))
})

test_that("sclet_ai_format_clarification keeps unknown errors and never guesses an answer", {
    errors <- c(
        "prerequisites not met for action t1: start_cluster_missing: root must be explicit.",
        "prerequisites not met for action a2: reference_missing: supply ref.",
        "prerequisites not met for action d1: labels_missing: supply labels.",
        "some completely unrecognized failure mode"
    )
    result <- sclet:::sclet_ai_format_clarification(errors)
    expect_equal(result$raw_errors, errors)
    ids <- vapply(result$questions, function(q) q$id, character(1L))
    expect_true(all(c("trajectory_root", "annotation_reference", "annotation_labels",
        "unclassified") %in% ids))
    unknown <- Filter(function(q) identical(q$id, "unclassified"), result$questions)
    expect_equal(length(unknown), 1L)
    expect_equal(unknown[[1L]]$source_error, "some completely unrecognized failure mode")
    expect_false(unknown[[1L]]$recognized)
    # every entry stays a question: none of them carries an answer value
    expect_true(all(vapply(result$questions, function(q) is.null(q$answer), logical(1L))))
})

test_that("sclet_ai_format_clarification returns a typed status when nothing is blocked", {
    result <- sclet:::sclet_ai_format_clarification(character())
    expect_equal(result$status, "no_clarification_needed")
    expect_length(result$questions, 0L)
    expect_equal(result$raw_errors, character())
})

test_that("RunAIAnalysis surfaces a structured clarification report on invalid_plan", {
    set.seed(1)
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(rpois(40, 5), nrow = 10L, ncol = 4L))
    )
    SummarizedExperiment::colData(sce)$batch <- c("a", "a", "b", "b")
    old <- options(sclet.ai.call = function(task, context, ...) {
        list(
            answer = "Integrate the batches.",
            proposed_actions = list(list(
                id = "integrate", action = "run_integration",
                params = list(batch = "batch", method = "fastMNN")
            ))
        )
    })
    on.exit(options(old), add = TRUE)

    result <- RunAIAnalysis(sce, goal = "Integrate batches", confirm = TRUE,
        include = c("read", "integration"))
    expect_equal(result$status, "invalid_plan")
    expect_true(any(grepl("design_semantics_not_confirmed", result$validation$errors)))
    expect_equal(result$report$clarification$status, "clarification_required")
    expect_true(length(result$report$clarification$questions) >= 1L)
    expect_equal(result$report$clarification$questions[[1L]]$id, "design_batch")
    expect_equal(result$report$clarification$raw_errors, result$validation$errors)
})

test_that("sclet_ai_record_clarification_response writes a user_decision evidence node", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, 2L, 4L)),
        colData = S4Vectors::DataFrame(batch = rep(c("a", "b"), 2L))
    )
    updated <- sclet:::sclet_ai_record_clarification_response(
        sce,
        question_id = "design_batch",
        answer = "batch",
        blocked_action = "run_integration"
    )
    ledger <- GetAnalysisLedger(updated, detail = "summary",
        include_artifacts = FALSE, include_data = FALSE)
    found <- Filter(function(x) identical(x$summary$kind, "user_decision"),
        ledger$state_records$ai_evidence)
    expect_true(length(found) >= 1L)
    node <- found[[1L]]$summary$values
    expect_true(node$question_is_known)
    expect_true(node$answer_column_present)
    expect_false(node$free_text_recorded)
    # the answer text itself is never written to the ledger
    expect_false(any(grepl("batch", node, fixed = TRUE)))
    # and the node is retrievable through the normal evidence query
    expect_true("ev:clarification_q1_a11" %in% names(sclet:::sclet_ai_evidence_get_all(updated)))
})

test_that("clarification recording never stores the free-text answer", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, 2L, 4L)),
        colData = S4Vectors::DataFrame(batch = rep(c("a", "b"), 2L))
    )
    updated <- sclet:::sclet_ai_record_clarification_response(
        sce, question_id = "design_batch", answer = "patient 42 has stage 3 disease"
    )
    node <- sclet:::sclet_ai_evidence_get_all(updated)[[1L]]
    expect_false(node$values$answer_column_present)
    dumped <- paste(capture.output(str(node)), collapse = " ")
    expect_false(grepl("patient", dumped))
    expect_false(grepl("stage 3", dumped))
})

test_that("ResolveAIClarifications refuses to prompt in a non-interactive session", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, 4L, 4L)),
        colData = S4Vectors::DataFrame(batch = rep(c("a", "b"), 2L))
    )
    clarification <- sclet:::sclet_ai_format_clarification(
        "prerequisites not met for action i: design_semantics_not_confirmed: x"
    )
    result <- ResolveAIClarifications(sce, clarification)
    expect_equal(result$status, "needs_interactive")
    expect_equal(result$resolved, character())
    expect_equal(result$skipped, "design_batch")
    # nothing was written: no confirmation, no evidence, object untouched
    expect_identical(result$object, sce)
    expect_null(sclet_get_state_record(result$object, "ai_design_confirmation", "design;batch=batch"))
    expect_length(sclet:::sclet_ai_evidence_get_all(result$object), 0L)
})

test_that("ResolveAIClarifications applies a design answer and records a user_decision", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, 4L, 4L)),
        colData = S4Vectors::DataFrame(batch = rep(c("a", "b"), 2L))
    )
    clarification <- sclet:::sclet_ai_format_clarification(
        "prerequisites not met for action i: design_semantics_not_confirmed: x"
    )
    result <- ResolveAIClarifications(sce, clarification, ask = function(q) "batch")
    expect_equal(result$status, "resolved")
    expect_equal(result$resolved, "design_batch")

    record <- sclet_get_state_record(result$object, "ai_design_confirmation", "design;batch=batch")
    expect_false(is.null(record))
    expect_equal(record$inputs$batch, "batch")

    decisions <- Filter(function(x) identical(x$summary$kind, "user_decision"),
        GetAnalysisLedger(result$object, detail = "summary",
            include_artifacts = FALSE, include_data = FALSE)$state_records$ai_evidence)
    expect_true(length(decisions) >= 1L)
    expect_equal(result$transcript[[1L]]$outcome, "resolved")
    expect_true(result$transcript[[1L]]$recorded)
})

test_that("ResolveAIClarifications never confirms an answer that is not a real column", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, 4L, 4L)),
        colData = S4Vectors::DataFrame(batch = rep(c("a", "b"), 2L))
    )
    clarification <- sclet:::sclet_ai_format_clarification(
        "prerequisites not met for action i: design_semantics_not_confirmed: x"
    )
    result <- ResolveAIClarifications(sce, clarification, ask = function(q) "not_a_column")
    expect_equal(result$status, "unresolved")
    expect_equal(result$invalid, "design_batch")
    expect_equal(result$resolved, character())
    # a rejected answer must not leave a confirmation or a decision behind
    expect_null(sclet_get_state_record(result$object, "ai_design_confirmation", "design;batch=not_a_column"))
    expect_length(sclet:::sclet_ai_evidence_get_all(result$object), 0L)
    expect_false(result$transcript[[1L]]$applied)
})

test_that("ResolveAIClarifications treats an empty answer as skipped, never auto-filled", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, 4L, 4L)),
        colData = S4Vectors::DataFrame(batch = rep(c("a", "b"), 2L))
    )
    clarification <- sclet:::sclet_ai_format_clarification(
        c(
            "prerequisites not met for action i: design_semantics_not_confirmed: x",
            "prerequisites not met for action a: reference_missing: y"
        )
    )
    result <- ResolveAIClarifications(sce, clarification, ask = function(q) "")
    expect_equal(result$status, "unresolved")
    expect_setequal(result$skipped, c("design_batch", "annotation_reference"))
    expect_equal(result$resolved, character())
    expect_identical(result$object, sce)
    expect_length(sclet:::sclet_ai_evidence_get_all(result$object), 0L)
})

test_that("ResolveAIClarifications uses apply_answer for questions with no built-in confirmation", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, 4L, 4L)),
        colData = S4Vectors::DataFrame(cluster = rep(c("A", "B"), 2L))
    )
    clarification <- sclet:::sclet_ai_format_clarification(
        "prerequisites not met for action t: start_cluster_missing: root required."
    )
    seen <- character()
    result <- ResolveAIClarifications(
        sce, clarification,
        ask = function(q) "A",
        apply_answer = function(object, question, answer) {
            seen <<- c(seen, answer)
            object
        }
    )
    expect_equal(result$status, "resolved")
    expect_equal(result$resolved, "trajectory_root")
    expect_equal(seen, "A")
    # recorded for the audit trail even though nothing was applied to the object
    expect_true(length(sclet:::sclet_ai_evidence_get_all(result$object)) >= 1L)
})

test_that("ResolveAIClarifications can resolve without recording evidence", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, 4L, 4L)),
        colData = S4Vectors::DataFrame(batch = rep(c("a", "b"), 2L))
    )
    clarification <- sclet:::sclet_ai_format_clarification(
        "prerequisites not met for action i: design_semantics_not_confirmed: x"
    )
    result <- ResolveAIClarifications(sce, clarification, ask = function(q) "batch", record = FALSE)
    expect_equal(result$status, "resolved")
    expect_false(result$transcript[[1L]]$recorded)
    expect_length(sclet:::sclet_ai_evidence_get_all(result$object), 0L)
    expect_false(is.null(sclet_get_state_record(result$object, "ai_design_confirmation", "design;batch=batch")))
})

test_that("resolving a clarification unblocks the plan that was previously rejected", {
    set.seed(1)
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(rpois(40, 5), nrow = 10L, ncol = 4L))
    )
    SummarizedExperiment::colData(sce)$batch <- c("a", "a", "b", "b")
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "integration"))
    plan <- sclet:::new_sclet_ai_plan(
        task = "integrate",
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
        actions = list(list(id = "i", action = "run_integration",
            params = list(batch = "batch", method = "fastMNN")))
    )
    before <- ValidateAIPlan(plan, object = sce, registry = registry, strict = FALSE)
    expect_true(any(grepl("design_semantics_not_confirmed", before$errors)))

    solved <- ResolveAIClarifications(
        sce,
        sclet:::sclet_ai_format_clarification(before$errors),
        ask = function(q) "batch"
    )
    expect_equal(solved$status, "resolved")

    # the design error is gone; only the expected fingerprint staleness remains
    after <- ValidateAIPlan(plan, object = solved$object, registry = registry, strict = FALSE)
    expect_false(any(grepl("design_semantics_not_confirmed", after$errors)))
})

test_that("ResolveAIClarifications handles an empty clarification and rejects bad callbacks", {
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1, 4L, 4L)))
    empty <- ResolveAIClarifications(sce, list(questions = list()))
    expect_equal(empty$status, "nothing_to_resolve")
    clarification <- sclet:::sclet_ai_format_clarification("start_cluster_missing: x")
    expect_error(ResolveAIClarifications(sce, clarification, ask = "not a function"), "ask must be")
    expect_error(
        ResolveAIClarifications(sce, clarification, apply_answer = "not a function"),
        "apply_answer must be"
    )
    expect_error(
        ResolveAIClarifications(list(counts = matrix(1)), clarification, ask = function(q) "a"),
        "SingleCellExperiment"
    )
})

test_that("clarification formatting does not alter ValidateAIPlan or prerequisites behavior", {
    set.seed(1)
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(rpois(40, 5), nrow = 10L, ncol = 4L))
    )
    SummarizedExperiment::colData(sce)$batch <- c("a", "a", "b", "b")
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "integration"))

    # (1) the AI self-declaring design confirmation is still rejected
    attack <- sclet:::new_sclet_ai_plan(
        task = "attacker_plan",
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
        actions = list(list(
            id = "integrate", action = "run_integration",
            params = list(batch = "batch", method = "fastMNN", .design_confirmed = TRUE)
        ))
    )
    attack_validation <- ValidateAIPlan(attack, object = sce, registry = registry, strict = FALSE)
    expect_false(attack_validation$valid)
    expect_true(any(grepl("design_semantics_not_confirmed", attack_validation$errors)))

    # (2) omitting the required parameter is still rejected
    missing_batch <- sclet:::new_sclet_ai_plan(
        task = "missing_param",
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
        actions = list(list(id = "integrate", action = "run_integration",
            params = list(method = "fastMNN")))
    )
    missing_validation <- ValidateAIPlan(missing_batch, object = sce, registry = registry, strict = FALSE)
    expect_false(missing_validation$valid)
    expect_true(any(grepl("missing required parameter", missing_validation$errors)))

    # (3) formatting is lossless and adds no decision of its own
    formatted <- sclet:::sclet_ai_format_clarification(attack_validation$errors)
    expect_equal(formatted$raw_errors, attack_validation$errors)
    ids <- vapply(formatted$questions, function(q) q$id, character(1L))
    # every rejection reason survives, recognized or not
    expect_true(all(c("design_batch", "unclassified") %in% ids))
    design_question <- Filter(function(q) identical(q$id, "design_batch"), formatted$questions)
    expect_equal(length(design_question), 1L)
    expect_equal(design_question[[1L]]$blocked_action, "run_integration")
    expect_equal(design_question[[1L]]$related_function, "ConfirmAIDesignSemantics")

    # (4) the formatter must not mutate the object or the validation result
    before <- GetAnalysisLedger(sce)$fingerprint
    invisible(sclet:::sclet_ai_format_clarification(attack_validation$errors))
    expect_equal(GetAnalysisLedger(sce)$fingerprint, before)
    repeat_plan <- sclet:::new_sclet_ai_plan(
        task = "attacker_plan",
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
        actions = attack$actions
    )
    expect_equal(
        ValidateAIPlan(repeat_plan, object = sce, registry = registry, strict = FALSE)$errors,
        attack_validation$errors
    )
})
