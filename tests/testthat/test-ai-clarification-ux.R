## Structured clarification: the whole surface.
##
## Owns every clarification test: mapping validation errors to questions,
## recording answers, resolving them interactively, the planner channel, and the
## regression guard that the judgment layer is untouched. The first-round blocks
## were moved here from test-ai-beginner.R on 2026-10-18 without changing an
## assertion; the beginner file keeps the facade (RunBasicWorkflow/AskAI/
## RunAIAnalysis) it is named for.
##
## The round-two additions pin:
##
##   * `cci_method` is a known question (method_missing / method_invalid);
##   * an accepted answer is written into the evidence node's summary when it is
##     a public vocabulary term (engine, celldex dataset, existing column,
##     anonymized cluster code), and that summary is what the planner reads back
##     as `analysis_story$user_decisions`;
##   * an answer the blocked action would reject is reported `invalid` and is
##     neither applied nor recorded.
##
## The judgment layer is untouched: one test asserts that validation still
## rejects exactly what it rejected before, and that formatting copies the
## errors rather than rewriting them.

clarification_test_object <- function() {
    SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, 2L, 4L)),
        colData = S4Vectors::DataFrame(
            batch = rep(c("a", "b"), 2L),
            cluster = factor(rep(c("c1", "c2"), 2L))
        )
    )
}

clarification_test_evidence <- function(object) {
    ledger <- GetAnalysisLedger(object, detail = "summary",
        include_artifacts = FALSE, include_data = FALSE)
    nodes <- ledger$state_records$ai_evidence %||% list()
    summaries <- lapply(nodes, function(record) record$summary)
    Filter(function(node) identical(node$kind, "user_decision"), summaries)
}


## -- first-round contract (moved from test-ai-beginner.R, 2026-10-18)
## The beginner file is about RunBasicWorkflow/AskAI/RunAIAnalysis; these
## blocks test the clarification surface itself, so they live with the
## rest of it. No assertion was changed in the move.

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
    actions <- sort(unique(names(AIDefaultExecutionRegistry(updated, include = "all"))))
    expected_id <- paste0(
        "ev:clarification_q",
        match("design_batch", sclet:::sclet_ai_clarification_question_ids()),
        "_a",
        match("run_integration", actions)
    )
    expect_true(expected_id %in% names(sclet:::sclet_ai_evidence_get_all(updated)))
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

## -- the engine question exists at all ---------------------------------------

test_that("method errors map to a run_cci engine question", {
    missing <- "prerequisites not met for action cc1: method_missing: run_cci requires an explicit engine choice."
    invalid <- "prerequisites not met for action cc2: method_invalid: run_cci accepts method = \"CellChat\", \"CellPhoneDB\" or \"NicheNet\"."
    ## One question per blocked choice, whichever of the two errors produced it.
    for (errors in list(missing, invalid, c(missing, invalid))) {
        result <- sclet:::sclet_ai_format_clarification(errors)
        expect_equal(result$status, "clarification_required")
        expect_length(result$questions, 1L)
        question <- result$questions[[1L]]
        expect_equal(question$id, "cci_method")
        expect_true(question$recognized)
        expect_equal(question$blocked_action, "run_cci")
        expect_equal(question$related_function, "run_cci")
        ## Names the choices, never chooses one.
        expect_true(grepl("CellChat", question$text, fixed = TRUE))
        expect_null(question$answer)
        ## Untouched raw errors: formatting is a copy, not a rewrite.
        expect_equal(result$raw_errors, errors)
    }
})

test_that("a run_cci plan without an engine is still rejected the same way", {
    sce <- clarification_test_object()
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "communication"))
    plan <- sclet:::new_sclet_ai_plan(
        task = "cell communication",
        actions = list(list(id = "cc1", action = "run_cci", params = list(group = "cluster")))
    )
    err <- tryCatch(ValidateAIPlan(plan, object = sce, registry = registry, strict = TRUE),
        error = function(e) e)
    expect_s3_class(err, "sclet_ai_invalid_plan")
    message <- conditionMessage(err)
    expect_true(grepl("method_missing", message, fixed = TRUE))
    ## ... and the presentation layer turns that exact message into the question
    ## without changing it.
    clarification <- sclet:::sclet_ai_format_clarification(message)
    expect_equal(clarification$questions[[1L]]$id, "cci_method")
    expect_equal(clarification$raw_errors, message)
})

## -- an accepted answer reaches the next planning round -----------------------

test_that("a public-vocabulary answer is recorded and visible to the planner", {
    sce <- clarification_test_object()
    updated <- sclet:::sclet_ai_record_clarification_response(
        sce, question_id = "cci_method", answer = "CellChat",
        blocked_action = "run_cci")

    nodes <- clarification_test_evidence(updated)
    expect_length(nodes, 1L)
    node <- nodes[[1L]]
    expect_true(node$values$answer_is_public_vocabulary)
    expect_equal(node$values$answer_code,
        match("CellChat", sclet:::sclet_ai_clarification_engines()))
    ## The numeric payload stays inside the shared evidence guard.
    expect_true(sclet:::sclet_ai_evidence_value_ok(node$values))

    ## The summary names the choice, and it is what the story carries.
    story <- GetAnalysisLedger(updated)$analysis_story
    expect_true(length(story$user_decisions) >= 1L)
    decision <- story$user_decisions[[1L]]$summary
    expect_true(grepl("CellChat", decision, fixed = TRUE))
    expect_true(grepl("cci_method", decision, fixed = TRUE))

    ## The planning context is the actual channel: this is what the model reads.
    context <- sclet:::sclet_ai_context(updated, target = "planning")
    summaries <- vapply(context$analysis_story$user_decisions,
        function(x) as.character(x$summary %||% ""), character(1L))
    expect_true(any(grepl("CellChat", summaries, fixed = TRUE)))
})

test_that("a free-text answer is recorded without absorbing the text", {
    sce <- clarification_test_object()
    secret <- "sample_42_patient_note"
    updated <- sclet:::sclet_ai_record_clarification_response(
        sce, question_id = "unclassified", answer = secret)

    node <- clarification_test_evidence(updated)[[1L]]
    expect_false(node$values$answer_is_public_vocabulary)
    expect_equal(node$values$answer_code, 0L)
    expect_false(node$values$free_text_recorded)

    expect_true(grepl("value withheld", node$summary, fixed = TRUE))
    ## The text is nowhere: not in the values, not in the story, not in context.
    story <- GetAnalysisLedger(updated)$analysis_story
    dumped <- paste(
        utils::capture.output(str(node)),
        utils::capture.output(str(story)),
        utils::capture.output(str(sclet:::sclet_ai_context(updated, target = "planning"))),
        collapse = "\n"
    )
    expect_false(grepl(secret, dumped, fixed = TRUE))
})

## -- unusable answers are reported, never recorded ---------------------------

test_that("ResolveAIClarifications rejects an engine the action does not register", {
    sce <- clarification_test_object()
    clarification <- sclet:::sclet_ai_format_clarification(
        "prerequisites not met for action cc1: method_missing: choose an engine.")

    bad <- ResolveAIClarifications(sce, clarification,
        ask = function(question) "LIANA", record = TRUE)
    expect_equal(bad$status, "unresolved")
    expect_equal(bad$invalid, "cci_method")
    expect_equal(bad$transcript[[1L]]$outcome, "invalid")
    expect_length(clarification_test_evidence(bad$object), 0L)

    good <- ResolveAIClarifications(sce, clarification,
        ask = function(question) "CellPhoneDB", record = TRUE)
    expect_equal(good$status, "resolved")
    expect_equal(good$resolved, "cci_method")
    nodes <- clarification_test_evidence(good$object)
    expect_length(nodes, 1L)
    expect_equal(nodes[[1L]]$values$answer_code,
        match("CellPhoneDB", sclet:::sclet_ai_clarification_engines()))
})

test_that("ResolveAIClarifications keeps the design-confirmation behaviour", {
    sce <- clarification_test_object()
    clarification <- sclet:::sclet_ai_format_clarification(
        "prerequisites not met for action integrate: design_semantics_not_confirmed: confirm the batch column.")

    bad <- ResolveAIClarifications(sce, clarification,
        ask = function(question) "not_a_column", record = TRUE)
    expect_equal(bad$status, "unresolved")
    expect_length(clarification_test_evidence(bad$object), 0L)

    good <- ResolveAIClarifications(sce, clarification,
        ask = function(question) "batch", record = TRUE)
    expect_equal(good$status, "resolved")
    ## The existing confirmation mechanism still runs, unchanged.
    confirmations <- sclet:::sclet_get_state_records(good$object, "ai_design_confirmation")
    expect_length(confirmations, 1L)
    expect_equal(confirmations[[1L]]$summary$design, list(batch = "batch"))
})

test_that("ResolveAIClarifications validates an annotation reference answer", {
    sce <- clarification_test_object()
    clarification <- sclet:::sclet_ai_format_clarification(
        "prerequisites not met for action a1: reference_missing: supply ref.")

    bad <- ResolveAIClarifications(sce, clarification,
        ask = function(question) "HELD - user-supplied celldex reference required")
    expect_equal(bad$status, "unresolved")
    expect_equal(bad$invalid, "annotation_reference")

    good <- ResolveAIClarifications(sce, clarification,
        ask = function(question) "celldex::MonacoImmuneData")
    expect_equal(good$status, "resolved")
    node <- clarification_test_evidence(good$object)[[1L]]
    expect_true(node$values$answer_is_public_vocabulary)
    expect_true(grepl("MonacoImmuneData", node$summary, fixed = TRUE))
})

## -- the prompt tells the planner that recorded choices are settled ----------

test_that("the planning prompt carries the user-decision rule", {
    line <- sclet:::sclet_ai_user_decision_prompt_line()
    expect_true(is.character(line) && length(line) == 1L)
    ## Names the field the planner must read...
    expect_true(grepl("analysis_story$user_decisions", line, fixed = TRUE))
    ## ... says the choice is settled, and forbids re-asking or substituting.
    expect_true(grepl("settled", line, fixed = TRUE))
    expect_true(grepl("do not ask", line, fixed = TRUE))
    expect_true(grepl("do not substitute", line, fixed = TRUE))
})

## -- the transcript distinguishes "applied" from "recorded" ------------------

test_that("a plan-parameter answer is recorded without being reported as applied", {
    sce <- clarification_test_object()
    clarification <- sclet:::sclet_ai_format_clarification(
        "prerequisites not met for action cc1: method_missing: choose an engine.")
    result <- ResolveAIClarifications(sce, clarification,
        ask = function(question) "NicheNet", record = TRUE)
    entry <- result$transcript[[1L]]
    expect_equal(entry$outcome, "resolved")
    ## Nothing on the object changed: the engine answer is a plan parameter.
    expect_false(entry$applied)
    ## But the decision IS recorded, which is what the planner reads next round.
    expect_true(entry$recorded)
    expect_true(identical(result$object, sce) ||
        length(clarification_test_evidence(result$object)) == 1L)

    ## A design confirmation IS applied, so both flags are true.
    design <- sclet:::sclet_ai_format_clarification(
        "prerequisites not met for action i: design_semantics_not_confirmed: x")
    confirmed <- ResolveAIClarifications(sce, design,
        ask = function(question) "batch", record = TRUE)
    expect_true(confirmed$transcript[[1L]]$applied)
    expect_true(confirmed$transcript[[1L]]$recorded)

    ## So is an answer handled by a caller-supplied applier.
    root <- sclet:::sclet_ai_format_clarification(
        "prerequisites not met for action t: start_cluster_missing: root required.")
    hooked <- ResolveAIClarifications(sce, root, ask = function(question) "cluster_1",
        apply_answer = function(object, question, answer) object, record = TRUE)
    expect_true(hooked$transcript[[1L]]$applied)
})

## -- one engine vocabulary, so the recorded answer is a real engine ----------

test_that("the clarification engine list is the one RunCCI actually accepts", {
    ## The public definition is RunCCI()'s formal default; match.arg() there is
    ## what rejects an unknown engine. The clarification helper must agree with
    ## it, or an answer could be validated as an engine the action refuses.
    accepted <- eval(formals(RunCCI)$method)
    expect_equal(sclet:::sclet_ai_clarification_engines(), accepted)
    ## An engine the action accepts is a valid answer; anything else is not.
    sce <- clarification_test_object()
    for (engine in accepted) {
        expect_null(sclet:::sclet_ai_clarification_answer_problem("cci_method", engine, sce))
    }
    expect_false(is.null(sclet:::sclet_ai_clarification_answer_problem(
        "cci_method", "LIANA", sce)))
})

## -- a second answer to the same question is not silently dropped -----------

test_that("re-answering a question records a second decision instead of colliding", {
    sce <- clarification_test_object()
    first <- sclet:::sclet_ai_record_clarification_response(
        sce, "cci_method", "CellChat", "run_cci")
    second <- sclet:::sclet_ai_record_clarification_response(
        first, "cci_method", "NicheNet", "run_cci")

    nodes <- clarification_test_evidence(second)
    expect_length(nodes, 2L)
    summaries <- vapply(nodes, function(node) as.character(node$summary), character(1L))
    expect_true(any(grepl("CellChat", summaries, fixed = TRUE)))
    expect_true(any(grepl("NicheNet", summaries, fixed = TRUE)))
    ids <- vapply(nodes, function(node) as.character(node$id), character(1L))
    expect_equal(length(unique(ids)), 2L)
    ## The first answer keeps the plain id; the second is suffixed. Compute the
    ## expected id instead of hardcoding it: both components depend on catalog
    ## order and on the registered action list, so a literal would break for a
    ## reason unrelated to this behaviour.
    question_code <- match("cci_method", sclet:::sclet_ai_clarification_question_ids())
    action_code <- match("run_cci",
        sort(unique(names(AIDefaultExecutionRegistry(sce, include = "all")))))
    base_id <- paste0("ev:clarification_q", question_code, "_a", action_code)
    expect_true(base_id %in% ids)
    expect_true(paste0(base_id, "_2") %in% ids)
    ## Both reach the planner, in order, so the latest choice is visible too.
    context <- sclet:::sclet_ai_context(second, target = "planning")
    ctx_summaries <- vapply(context$analysis_story$user_decisions,
        function(x) as.character(x$summary %||% ""), character(1L))
    expect_true(any(grepl("NicheNet", ctx_summaries, fixed = TRUE)))
})
