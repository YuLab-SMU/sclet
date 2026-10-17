## Structured clarification, round two (2026-10-17).
##
## The first implementation ("turn an error into a question, record the answer")
## left two things open that the roadmap named: the communication engine had no
## question at all, and an answer never reached the next planning round -- the
## caller had to thread it back into the goal text by hand. These tests pin the
## two additions:
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
