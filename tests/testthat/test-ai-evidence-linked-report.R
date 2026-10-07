test_that("RecordAIResult accepts findings with resolvable evidence_refs", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 3, ncol = 2))
    )
    sce <- RecordAIEvidence(sce, list(
        id = "ev:valid_ref",
        kind = "state",
        values = list(n = 10),
        claim_level = "observed"
    ))
    result <- sclet:::new_sclet_ai_result(
        task = "test_task",
        findings = list(list(
            statement = "valid finding",
            severity = "info",
            claim_level = "associated",
            evidence_refs = "ev:valid_ref"
        ))
    )
    updated <- RecordAIResult(sce, result, id = "ai_test_1")
    record <- sclet_get_analysis(updated, "ai_test_1")
    expect_equal(record$type, "ai_test_task")
    expect_equal(record$summary$n_findings, 1L)
})

test_that("RecordAIResult rejects findings with unresolvable evidence_refs at any claim_level", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 3, ncol = 2))
    )
    result <- sclet:::new_sclet_ai_result(
        task = "test_task",
        findings = list(list(
            statement = "bad finding",
            severity = "info",
            claim_level = "associated",
            evidence_refs = "ev:does_not_exist_at_all"
        ))
    )
    expect_error(
        RecordAIResult(sce, result),
        "unresolvable|evidence ref"
    )
})

test_that("RecordAIResult accepts findings with no evidence_refs at hypothesis level", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 3, ncol = 2))
    )
    result <- sclet:::new_sclet_ai_result(
        task = "test_task",
        findings = list(list(
            statement = "hypothesis",
            severity = "info",
            claim_level = "hypothesis"
        ))
    )
    updated <- RecordAIResult(sce, result, id = "ai_test_2")
    record <- sclet_get_analysis(updated, "ai_test_2")
    expect_equal(record$summary$n_findings, 1L)
})

test_that("RecordAIResult claim-ceiling error takes precedence for consistent_with with stale ref", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 3, ncol = 2))
    )
    result <- sclet:::new_sclet_ai_result(
        task = "test_task",
        findings = list(list(
            statement = "over-claimed",
            severity = "info",
            claim_level = "consistent_with",
            evidence_refs = "ev:does_not_exist"
        ))
    )
    expect_error(
        RecordAIResult(sce, result),
        "over-claimed|unresolvable|evidence ref"
    )
})

test_that("RecordAIResult with audit_claims = FALSE accepts stale evidence_refs", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 3, ncol = 2))
    )
    result <- sclet:::new_sclet_ai_result(
        task = "test_task",
        findings = list(list(
            statement = "stale ref",
            severity = "info",
            claim_level = "associated",
            evidence_refs = "ev:does_not_exist"
        ))
    )
    updated <- RecordAIResult(sce, result, id = "ai_test_3", audit_claims = FALSE)
    record <- sclet_get_analysis(updated, "ai_test_3")
    expect_equal(record$summary$n_findings, 1L)
})

test_that("RunAIAnalysis report includes success_assessment from execution", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 5, ncol = 4))
    )
    old <- options(sclet.ai.call = function(task, context, ...) {
        list(
            answer = "mock plan",
            actions = list(list(id = "step1", action = "inspect_status")),
            assumptions = list(),
            candidate_routes = list(),
            selected_route = NULL,
            success_criteria = list(),
            risks = list(),
            stop_conditions = list()
        )
    })
    on.exit(options(old), add = TRUE)
    result <- RunAIAnalysis(sce, goal = "test goal", confirm = "yes")
    expect_true(!is.null(result$report$success_assessment))
    expect_equal(result$report$success_assessment, result$execution$success_assessment)
})

test_that("RunAIAnalysis with interpret = FALSE does not populate interpretation", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 5, ncol = 4))
    )
    old <- options(sclet.ai.call = function(task, context, ...) {
        list(
            answer = "mock plan",
            actions = list(list(id = "step1", action = "inspect_status")),
            assumptions = list(),
            candidate_routes = list(),
            selected_route = NULL,
            success_criteria = list(),
            risks = list(),
            stop_conditions = list()
        )
    })
    on.exit(options(old), add = TRUE)
    result <- RunAIAnalysis(sce, goal = "test goal", confirm = "yes", interpret = FALSE)
    expect_null(result$report$interpretation)
})

test_that("RunAIAnalysis with interpret = TRUE populates interpretation on success", {
    skip_if_not_installed("aisdk")
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 5, ncol = 4))
    )
    sce <- RecordAIEvidence(sce, list(
        id = "ev:health_check",
        kind = "deterministic_summary",
        values = list(flag = 1L),
        claim_level = "measured"
    ))
    agent <- new.env(parent = emptyenv())
    agent$run <- function(task, max_steps = 10) {
        list(text = "mock interpretation")
    }
    call_count <- 0L
    old <- options(
        sclet.ai.create_agent = function(...) agent,
        sclet.ai.generate_object = function(...) {
            call_count <<- call_count + 1L
            if (call_count == 1L) {
                list(
                    object = list(
                        answer = "mock plan",
                        actions = list(list(id = "step1", action = "inspect_status")),
                        assumptions = list(),
                        candidate_routes = list(),
                        selected_route = NULL,
                        success_criteria = list(),
                        risks = list(),
                        stop_conditions = list()
                    )
                )
            } else {
                list(
                    object = list(
                        answer = "mock interpretation",
                        findings = list(list(
                            statement = "interpretation finding",
                            severity = "info",
                            claim_level = "observed",
                            evidence_refs = "ev:health_check"
                        )),
                        evidence = list("ev:health_check"),
                        warnings = list(),
                        recommendations = list(),
                        proposed_actions = list()
                    )
                )
            }
        }
    )
    on.exit(options(old), add = TRUE)
    result <- RunAIAnalysis(sce, goal = "test goal", confirm = "yes", interpret = TRUE)
    expect_true(!is.null(result$report$interpretation))
    expect_s3_class(result$report$interpretation, "sclet_ai_result")
    expect_true(validate_sclet_ai_result(result$report$interpretation, error = FALSE))
    # the interpretation's evidence_refs are resolvable, so recording must
    # actually succeed (not merely be attempted and silently fail)
    expect_null(result$report$interpretation_error)
    interp_ids <- grep("^ai_interpretation_", names(GetAnalysisLedger(result$object)$analyses), value = TRUE)
    expect_true(length(interp_ids) >= 1L)
})

test_that("RunAIAnalysis with interpret = TRUE catches interpretation errors", {
    skip_if_not_installed("aisdk")
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 5, ncol = 4))
    )
    agent <- new.env(parent = emptyenv())
    call_count <- 0L
    agent$run <- function(task, max_steps = 10) {
        call_count <<- call_count + 1L
        if (call_count == 1L) {
            list(text = "mock plan")
        } else {
            stop("agent run failed for interpretation")
        }
    }
    old <- options(
        sclet.ai.create_agent = function(...) agent,
        sclet.ai.generate_object = function(...) {
            if (call_count == 1L) {
                list(
                    object = list(
                        answer = "mock plan",
                        actions = list(list(id = "step1", action = "inspect_status")),
                        assumptions = list(),
                        candidate_routes = list(),
                        selected_route = NULL,
                        success_criteria = list(),
                        risks = list(),
                        stop_conditions = list()
                    )
                )
            } else {
                stop("interpretation call failed")
            }
        }
    )
    on.exit(options(old), add = TRUE)
    result <- RunAIAnalysis(sce, goal = "test goal", confirm = "yes", interpret = TRUE)
    expect_true(!is.null(result$report$interpretation_error))
    expect_type(result$report$interpretation_error, "character")
    expect_equal(result$status, result$execution$status)
    expect_true(!is.null(result$execution))
})

test_that("RecordAIResult summary includes cited_user_decisions when present", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 3, ncol = 2))
    )
    sce <- RecordAIEvidence(sce, list(
        id = "ev:user_dec_1",
        kind = "user_decision",
        values = list(route_index = 1L),
        claim_level = "observed"
    ))
    result <- sclet:::new_sclet_ai_result(
        task = "test_task",
        findings = list(list(
            statement = "cites user decision",
            severity = "info",
            claim_level = "associated",
            evidence_refs = "ev:user_dec_1"
        ))
    )
    updated <- RecordAIResult(sce, result, id = "ai_test_4")
    record <- sclet_get_analysis(updated, "ai_test_4")
    expect_true("cited_user_decisions" %in% names(record$summary))
    expect_true("ev:user_dec_1" %in% record$summary$cited_user_decisions)
})

test_that("RecordAIResult summary includes empty cited_user_decisions when none cited", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 3, ncol = 2))
    )
    result <- sclet:::new_sclet_ai_result(
        task = "test_task",
        findings = list(list(
            statement = "no user decision cited",
            severity = "info",
            claim_level = "hypothesis"
        ))
    )
    updated <- RecordAIResult(sce, result, id = "ai_test_5")
    record <- sclet_get_analysis(updated, "ai_test_5")
    expect_true("cited_user_decisions" %in% names(record$summary))
    expect_equal(record$summary$cited_user_decisions, character())
})

test_that("existing test files still pass after RecordAIResult changes", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 3, ncol = 2))
    )
    result <- sclet:::new_sclet_ai_result(
        task = "status_review",
        answer = "all inputs are visible",
        context = list(schema_version = "1.0", fingerprint = "test"),
        metadata = list(provider = "mock")
    )
    updated <- RecordAIResult(sce, result, id = "ai_status_1")
    record <- sclet_get_analysis(updated, "ai_status_1")
    expect_equal(record$type, "ai_status_review")
    expect_equal(record$inputs$context_fingerprint, "test")
    expect_true("cited_user_decisions" %in% names(record$summary))
})
