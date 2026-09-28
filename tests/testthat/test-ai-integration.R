test_that("GetAnalysisLedger returns a bounded, versioned AI view", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 5, ncol = 4))
    )
    sce <- sclet_set_analysis(
        sce, "augur",
        list(id = "augur", method = "Augur", AUC = data.frame(auc = 0.8))
    )
    sce <- sclet_set_analysis_state(
        sce, "priority", "augur", "Augur",
        inputs = list(group_col = "condition"),
        summary = list(best_cell_type = "T")
    )
    sce <- sclet_set_analysis(
        sce, "rare_cells",
        list(id = "rare_1", method = "density", labels = rep(FALSE, 4))
    )
    sce <- sclet_set_analysis_state(
        sce, "rare_cells", "rare_1", "density",
        artifacts = list(label_col = "rare_cluster")
    )

    ledger <- GetAnalysisLedger(sce, detail = "full", include_data = TRUE)
    expect_equal(ledger$schema_version, "1.0")
    expect_equal(ledger$dataset$n_cells, 4)
    expect_true(ledger$health$has_perturbation_priority)
    expect_true(ledger$health$has_rare_cells)
    expect_true(any(vapply(ledger$analyses, function(x) identical(x$type, "priority"), logical(1))))
    expect_true("priority" %in% names(ledger$state_records))
    expect_true("rare_cells" %in% names(ledger$state_records))
    expect_false(ledger$dataset$data$expression[["counts"]]$values_included)
    expect_true(nzchar(ledger$fingerprint))
})

test_that("legacy LLM context is generated from the AI ledger", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 3, ncol = 2))
    )
    sce <- sclet_set_analysis(
        sce, "state_priority",
        list(id = "state_priority", method = "state_priority", summary = list(ok = TRUE))
    )
    text <- SummarizeContextForLLM(sce)
    expect_type(text, "character")
    expect_true(grepl("SCLET AI LEDGER", text, fixed = TRUE))
    expect_true(grepl("state_priority", text, fixed = TRUE))
})

test_that("read-only AI tools expose deterministic handlers", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 3, ncol = 2))
    )
    tools <- sclet_ai_tool_registry(sce)
    expect_true(all(c("get_status", "get_ledger", "get_capabilities", "get_analysis_lineage", "get_quality_checks") %in% names(tools)))
    expect_true(all(vapply(tools, function(x) isTRUE(x$read_only), logical(1))))
    expect_equal(sclet_ai_call_tools(tools, "get_status")$n_commands, 0)
    expect_equal(length(sclet_ai_tool_registry(sce, mode = "execute")), 0)
    skip_if_not_installed("aisdk")
    native <- sclet_ai_native_tools(tools)
    expect_length(native, length(tools))
    expect_true(all(vapply(native, inherits, logical(1), what = "Tool")))
    expect_s3_class(sclet_ai_native_result_schema(), "z_schema")
})

test_that("mocked aisdk calls return structured AI results", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 3, ncol = 2))
    )
    old <- options(sclet.ai.call = function(task, context, ...) {
        list(
            answer = paste("mock answer for", task),
            findings = list(list(
                statement = "status was observed",
                severity = "info",
                claim_level = "observed",
                evidence_refs = "health"
            )),
            recommendations = list("inspect missing prerequisites")
        )
    })
    on.exit(options(old), add = TRUE)

    result <- AIStatus(sce)
    expect_s3_class(result, "sclet_ai_result")
    expect_equal(result$task, "status_review")
    expect_true(grepl("mock answer", result$answer, fixed = TRUE))
    expect_equal(validate_sclet_ai_result(result), TRUE)

    qc <- AIReviewQC(sce)
    expect_s3_class(qc, "sclet_ai_result")
    expect_equal(qc$task, "qc_review")
})

test_that("native tool loop can be followed by structured aisdk summary", {
    skip_if_not_installed("aisdk")
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 3, ncol = 2))
    )
    agent <- new.env(parent = emptyenv())
    agent$run <- function(task, max_steps = 10) {
        list(text = paste("tool loop completed", task))
    }
    old <- options(
        sclet.ai.create_agent = function(...) agent,
        sclet.ai.generate_object = function(...) {
            list(
                object = list(
                    answer = "structured summary",
                    findings = list(list(
                        statement = "the tool loop was completed",
                        severity = "info",
                        claim_level = "observed",
                        evidence_refs = "health"
                    )),
                    evidence = list("health"),
                    warnings = list(),
                    recommendations = list("continue review"),
                    proposed_actions = list()
                )
            )
        }
    )
    on.exit(options(old), add = TRUE)

    result <- AIStatus(sce)
    expect_s3_class(result, "sclet_ai_result")
    expect_equal(result$answer, "structured summary")
    expect_true(isTRUE(result$metadata$native_tools))
    expect_true(isTRUE(result$metadata$structured_output))
    expect_equal(result$findings[[1]]$claim_level, "observed")
})

test_that("structured summary failure falls back to the tool-loop response", {
    skip_if_not_installed("aisdk")
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 3, ncol = 2))
    )
    agent <- new.env(parent = emptyenv())
    agent$run <- function(task, max_steps = 10) {
        list(text = "tool response retained")
    }
    old <- options(
        sclet.ai.create_agent = function(...) agent,
        sclet.ai.generate_object = function(...) stop("schema endpoint unavailable")
    )
    on.exit(options(old), add = TRUE)

    result <- AIStatus(sce)
    expect_equal(result$answer, "tool response retained")
    expect_false(isTRUE(result$metadata$structured_output))
    expect_true(any(vapply(result$warnings, function(x) {
        identical(x$code, "structured_output_fallback")
    }, logical(1))))
})
test_that("AI result validation rejects unsupported causal claims", {
    expect_error(
        new_sclet_ai_result(
            "test",
            findings = list(list(claim_level = "causal"))
        ),
        class = "sclet_ai_invalid_output"
    )
})

test_that("AI results can be explicitly recorded without changing active view", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 3, ncol = 2))
    )
    result <- new_sclet_ai_result(
        "status_review",
        answer = "all inputs are visible",
        context = list(schema_version = "1.0", fingerprint = "test"),
        metadata = list(provider = "mock")
    )
    updated <- RecordAIResult(sce, result, id = "ai_status_1")
    record <- sclet_get_analysis(updated, "ai_status_1")
    expect_equal(record$type, "ai_status_review")
    expect_equal(record$inputs$context_fingerprint, "test")
    expect_equal(DefaultAssay(updated), "counts")
})
