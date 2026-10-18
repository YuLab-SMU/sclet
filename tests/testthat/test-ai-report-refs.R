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
