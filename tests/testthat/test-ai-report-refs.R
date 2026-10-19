## Report-layer citation coverage.
##
## Real-data runs wrote most findings without an evidence id (3/13 on the Kang
## contrast run, 2026-10-18); an unreferenced finding cannot be distinguished
## from a fabricated one, so the explanation prompt now asks for a per-finding
## reference and the result carries a coverage count.

test_that("the explanation prompt requires a per-finding evidence reference", {
    prompt <- sclet:::sclet_ai_explanation_prompt()
    expect_true(grepl("evidence_refs", prompt, fixed = TRUE))
    expect_true(grepl("ev:", prompt, fixed = TRUE))
    expect_true(grepl("hypothesis", prompt, fixed = TRUE))
    expect_true(grepl("not", prompt, fixed = TRUE))
})

test_that("ref coverage counts findings with and without evidence ids", {
    result <- new_sclet_ai_result(
        "analysis_explanation",
        answer = "two findings",
        context = list(schema_version = "1.0", fingerprint = "test"),
        findings = list(
            list(statement = "cited", evidence_refs = "ev:findallmarkers"),
            list(statement = "uncited", evidence_refs = character())
        )
    )
    attached <- sclet:::sclet_ai_attach_ref_coverage(result)
    coverage <- attached$metadata$evidence_ref_coverage
    expect_equal(coverage$n_findings, 2L)
    expect_equal(coverage$n_with_refs, 1L)
    expect_equal(coverage$fraction, 0.5)
})

test_that("the explanation prompt lists the real evidence ids, bounded", {
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1, 5, 4)))
    sce <- RecordAIEvidence(sce, list(id = "ev:findallmarkers",
        kind = "deterministic_summary", claim_level = "associated",
        values = list(n_rows = 3L)))
    sce <- RecordAIEvidence(sce, list(id = "ev:annotation_monaco",
        kind = "deterministic_summary", claim_level = "consistent_with",
        values = list(n_rows = 2L)))
    lines <- sclet:::sclet_ai_evidence_prompt_lines(sce)
    expect_true(any(grepl("ev:findallmarkers", lines, fixed = TRUE)))
    expect_true(any(grepl("ev:annotation_monaco", lines, fixed = TRUE)))
    prompt <- sclet:::sclet_ai_explanation_prompt(lines)
    expect_true(grepl("Use these exact strings", prompt, fixed = TRUE))
    expect_true(grepl("ev:findallmarkers", prompt, fixed = TRUE))

    ## An object with no evidence adds no list; the rules still apply.
    empty <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1, 5, 4)))
    expect_length(sclet:::sclet_ai_evidence_prompt_lines(empty), 0L)
    expect_false(grepl("Use these exact strings",
        sclet:::sclet_ai_explanation_prompt(character()), fixed = TRUE))
})

test_that("the id list is capped and says how many were omitted", {
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1, 5, 4)))
    for (i in seq_len(25L)) {
        sce <- RecordAIEvidence(sce, list(id = paste0("ev:e", i),
            kind = "deterministic_summary", claim_level = "associated",
            values = list(n_rows = as.integer(i))))
    }
    lines <- sclet:::sclet_ai_evidence_prompt_lines(sce, max_ids = 10L)
    expect_equal(sum(grepl("^- ev:", lines)), 10L)
    expect_true(any(grepl("more ids exist", lines)))
})

test_that("low citation coverage produces an explicit warning", {
    threshold <- sclet:::sclet_ai_coverage_warning
    expect_match(threshold(list(n_findings = 9L, n_with_refs = 2L, fraction = 2/9)),
        "only 2 of 9 findings cite an evidence id")
    expect_null(threshold(list(n_findings = 9L, n_with_refs = 8L, fraction = 8/9)))
    expect_null(threshold(list(n_findings = 1L, n_with_refs = 0L, fraction = 0)))
    expect_null(threshold(NULL))
})

test_that("explanation asks for one revision when citations are missing", {
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1, 5, 4)))
    sce <- RecordAIEvidence(sce, list(id = "ev:findallmarkers",
        kind = "deterministic_summary", claim_level = "associated",
        values = list(n_rows = 3L)))
    calls <- 0L
    old <- options(sclet.ai.call = function(task, context, ...) {
        calls <<- calls + 1L
        if (calls == 1L) {
            return(new_sclet_ai_result("analysis_explanation",
                context = list(schema_version = "1.0", fingerprint = "test"),
                findings = list(
                    list(statement = "a", evidence_refs = character()),
                    list(statement = "b", evidence_refs = character()),
                    list(statement = "c", evidence_refs = character()))))
        }
        new_sclet_ai_result("analysis_explanation",
            context = list(schema_version = "1.0", fingerprint = "test"),
            findings = list(
                list(statement = "a", evidence_refs = "ev:findallmarkers"),
                list(statement = "b", evidence_refs = "ev:findallmarkers"),
                list(statement = "c", evidence_refs = "ev:findallmarkers")))
    })
    on.exit(options(old), add = TRUE)

    result <- AIExplainAnalysis(sce, model = "mock", structured_output = FALSE)
    expect_equal(calls, 2L)
    expect_equal(result$metadata$evidence_ref_coverage$fraction, 1)
    expect_equal(result$metadata$explanation_revision_rounds, 1L)
})

test_that("a second revision runs when the first still leaves most findings uncited", {
    ## 2026-10-19: four real-data records at 0.3, 0.56, 0.5 and 1.0. One round
    ## was not enough to put the report layer above the 0.8 target, so a second
    ## round runs only while coverage stays below half -- and stops as soon as
    ## a round stops helping.
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1, 5, 4)))
    sce <- RecordAIEvidence(sce, list(id = "ev:findallmarkers",
        kind = "deterministic_summary", claim_level = "associated",
        values = list(n_rows = 3L)))
    uncited <- function() list(statement = "asserted", evidence_refs = character())
    incited <- function() list(statement = "cited", evidence_refs = "ev:findallmarkers")
    calls <- 0L
    old <- options(sclet.ai.call = function(task, context, ...) {
        calls <<- calls + 1L
        findings <- switch(as.character(calls),
            "1" = list(uncited(), uncited(), uncited()),
            "2" = list(incited(), uncited(), uncited()),
            list(incited(), incited(), incited()))
        new_sclet_ai_result("analysis_explanation",
            context = list(schema_version = "1.0", fingerprint = "test"),
            findings = findings)
    })
    on.exit(options(old), add = TRUE)

    result <- AIExplainAnalysis(sce, model = "mock", structured_output = FALSE)
    expect_equal(calls, 3L)
    expect_equal(result$metadata$evidence_ref_coverage$fraction, 1)
    expect_equal(result$metadata$explanation_revision_rounds, 2L)
})

test_that("a revision that does not help ends the loop", {
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1, 5, 4)))
    sce <- RecordAIEvidence(sce, list(id = "ev:findallmarkers",
        kind = "deterministic_summary", claim_level = "associated",
        values = list(n_rows = 3L)))
    uncited <- function() list(statement = "asserted", evidence_refs = character())
    calls <- 0L
    old <- options(sclet.ai.call = function(task, context, ...) {
        calls <<- calls + 1L
        new_sclet_ai_result("analysis_explanation",
            context = list(schema_version = "1.0", fingerprint = "test"),
            findings = list(uncited(), uncited(), uncited()))
    })
    on.exit(options(old), add = TRUE)

    result <- AIExplainAnalysis(sce, model = "mock", structured_output = FALSE)
    ## First proposal plus exactly one round: a round with no improvement is
    ## not repeated.
    expect_equal(calls, 2L)
    expect_null(result$metadata$explanation_revision_rounds)
    expect_equal(result$metadata$evidence_ref_coverage$fraction, 0)
})
