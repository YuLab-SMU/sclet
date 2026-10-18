## Invalid evidence citations are dropped, not allowed to delete the report.
##
## Kang contrast run (2026-10-18) cited an execution id as evidence; the whole
## interpretation recording was refused and the session lost its findings.

test_that("unresolvable refs are stripped and the finding is downgraded", {
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1, 5, 4)))
    sce <- RecordAIEvidence(sce, list(id = "ev:real", kind = "deterministic_summary",
        claim_level = "associated", values = list(n = 1L)))
    result <- new_sclet_ai_result("analysis_explanation",
        context = list(schema_version = "1.0", fingerprint = "test"),
        findings = list(
            list(statement = "supported", evidence_refs = "ev:real"),
            list(statement = "bogus", evidence_refs = "ai_execution_plan_x",
                claim_level = "associated"),
            list(statement = "mixed", evidence_refs = c("ev:real", "sclet-ai-1"))
        ))
    cleaned <- sclet:::sclet_ai_sanitize_finding_refs(sce, result)
    expect_equal(cleaned$findings[[1L]]$evidence_refs, "ev:real")
    expect_equal(cleaned$findings[[1L]]$claim_level, NULL)  # untouched
    expect_equal(length(cleaned$findings[[2L]]$evidence_refs), 0L)
    expect_equal(cleaned$findings[[2L]]$claim_level, "hypothesis")
    expect_true(isTRUE(cleaned$findings[[2L]]$refs_dropped))
    expect_equal(cleaned$findings[[3L]]$evidence_refs, "ev:real")
    expect_true("ai_execution_plan_x" %in% cleaned$metadata$dropped_evidence_refs)
    expect_true("sclet-ai-1" %in% cleaned$metadata$dropped_evidence_refs)
})
