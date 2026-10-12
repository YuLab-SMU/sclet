sclet_ai_test_evidence_object <- function() {
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrix(1L, nrow = 4L, ncol = 4L,
            dimnames = list(paste0("gene_", 1:4), paste0("cell_", 1:4))))
    )
    sce <- sclet:::RecordAIEvidence(sce, list(
        id = "ev:findallmarkers",
        kind = "deterministic_summary",
        values = list(n_rows = 100L, n_significant = 5L, top_up_ranks = c(1.5, 2.3)),
        claim_level = "associated"
    ))
    sce
}

test_that("qc context includes evidence ids", {
    sce <- sclet_ai_test_evidence_object()
    ctx <- sclet:::sclet_ai_context(sce, target = "qc", detail = "full")
    serialized <- paste(utils::capture.output(dput(ctx)), collapse = " ")
    expect_true(grepl("ev:findallmarkers", serialized, fixed = TRUE))
})

test_that("status context includes evidence ids", {
    sce <- sclet_ai_test_evidence_object()
    ctx <- sclet:::sclet_ai_context(sce, target = "status", detail = "full")
    serialized <- paste(utils::capture.output(dput(ctx)), collapse = " ")
    expect_true(grepl("ev:findallmarkers", serialized, fixed = TRUE))
})

test_that("planning context does not include evidence ids", {
    sce <- sclet_ai_test_evidence_object()
    ctx <- sclet:::sclet_ai_context(sce, target = "planning", detail = "full")
    serialized <- paste(utils::capture.output(dput(ctx)), collapse = " ")
    expect_false(grepl("ev:findallmarkers", serialized, fixed = TRUE))
})

test_that("evidence values are absent from qc context", {
    sce <- sclet_ai_test_evidence_object()
    ctx <- sclet:::sclet_ai_context(sce, target = "qc", detail = "full")
    serialized <- paste(utils::capture.output(dput(ctx)), collapse = " ")
    expect_false(grepl("n_rows", serialized, fixed = TRUE))
    expect_false(grepl("top_up_ranks", serialized, fixed = TRUE))
    expect_false(grepl("n_significant", serialized, fixed = TRUE))
})

test_that("evidence values are absent from status context", {
    sce <- sclet_ai_test_evidence_object()
    ctx <- sclet:::sclet_ai_context(sce, target = "status", detail = "full")
    serialized <- paste(utils::capture.output(dput(ctx)), collapse = " ")
    expect_false(grepl("n_rows", serialized, fixed = TRUE))
    expect_false(grepl("top_up_ranks", serialized, fixed = TRUE))
    expect_false(grepl("n_significant", serialized, fixed = TRUE))
})

test_that("report context still includes evidence values", {
    sce <- sclet_ai_test_evidence_object()
    ctx <- sclet:::sclet_ai_context(sce, target = "report", detail = "full")
    serialized <- paste(utils::capture.output(dput(ctx)), collapse = " ")
    expect_true(grepl("ev:findallmarkers", serialized, fixed = TRUE))
    expect_true(grepl("n_rows", serialized, fixed = TRUE))
})

test_that("evidence kind and claim_level survive in qc context", {
    sce <- sclet_ai_test_evidence_object()
    ctx <- sclet:::sclet_ai_context(sce, target = "qc", detail = "full")
    serialized <- paste(utils::capture.output(dput(ctx)), collapse = " ")
    expect_true(grepl("deterministic_summary", serialized, fixed = TRUE))
    expect_true(grepl("associated", serialized, fixed = TRUE))
})

test_that("non-completed evidence is excluded", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrix(1L, nrow = 2L, ncol = 2L))
    )
    sce <- sclet:::sclet_init_state(sce)
    state <- sclet:::sclet_get_state(sce)
    state$states$records$ai_evidence <- list(
        "ev:incomplete" = list(
            id = "ev:incomplete",
            type = "ai_evidence",
            status = "pending",
            method = "sclet_ai_evidence",
            summary = list(id = "ev:incomplete", kind = "test", status = "pending",
                values = list(n = 1L), claim_level = "measured")
        )
    )
    sce <- sclet:::sclet_set_state(sce, state)
    ctx <- sclet:::sclet_ai_context(sce, target = "qc", detail = "full")
    serialized <- paste(utils::capture.output(dput(ctx)), collapse = " ")
    expect_false(grepl("ev:incomplete", serialized, fixed = TRUE))
})

test_that("sanitizer produces no new warnings or redactions for evidence in qc/status", {
    sce <- sclet_ai_test_evidence_object()

    ctx_qc <- sclet:::sclet_ai_context(sce, target = "qc", detail = "full")
    san_qc <- sclet:::sclet_ai_sanitize_context(sce, ctx_qc, privacy = "standard")

    ctx_status <- sclet:::sclet_ai_context(sce, target = "status", detail = "full")
    san_status <- sclet:::sclet_ai_sanitize_context(sce, ctx_status, privacy = "standard")

    qc_ev_warnings <- san_qc$policy$warnings[grepl("ai_evidence|ev:", san_qc$policy$warnings)]
    status_ev_warnings <- san_status$policy$warnings[grepl("ai_evidence|ev:", san_status$policy$warnings)]
    expect_length(qc_ev_warnings, 0L)
    expect_length(status_ev_warnings, 0L)

    qc_ev_redacted <- san_qc$policy$redacted_fields[grepl("ai_evidence|ev:", san_qc$policy$redacted_fields)]
    status_ev_redacted <- san_status$policy$redacted_fields[grepl("ai_evidence|ev:", san_status$policy$redacted_fields)]
    expect_length(qc_ev_redacted, 0L)
    expect_length(status_ev_redacted, 0L)

    qc_serialized <- paste(utils::capture.output(dput(san_qc$payload)), collapse = " ")
    status_serialized <- paste(utils::capture.output(dput(san_status$payload)), collapse = " ")
    expect_false(grepl("n_rows", qc_serialized, fixed = TRUE))
    expect_false(grepl("n_rows", status_serialized, fixed = TRUE))
})

test_that("sclet_ai_target_matches includes ai_evidence for qc and status", {
    expect_true(sclet:::sclet_ai_target_matches("ai_evidence", "ev:x", "qc"))
    expect_true(sclet:::sclet_ai_target_matches("ai_evidence", "ev:x", "status"))
    expect_false(sclet:::sclet_ai_target_matches("ai_evidence", "ev:x", "planning"))
    expect_true(sclet:::sclet_ai_target_matches("ai_evidence", "ev:x", "report"))
})
