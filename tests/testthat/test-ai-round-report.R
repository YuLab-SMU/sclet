# M12-e: every execution round renders its report and carries the self-check.
#
# The claim: the report's self-check is what makes a failure visible, so a round
# that does not render has lost that check -- and the next round is then left
# reading a ledger that grows every round. So rendering moved from "an optional
# finishing move" to "the default end of a round".
#
# The other half of the claim is what the statuses must NOT collapse into:
#   - `failed` self-check is not `failed` execution; the run happened.
#   - `not_rendered` (nothing recorded) is not a pass.
#   - a broken render must never discard an execution that already ran.
#
# These tests drive sclet_ai_render_round_report() directly (it is the seam),
# so they need no provider and no real data.

m12e_recorded_sce <- function() {
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrix(1L, nrow = 4L, ncol = 4L))
    )
    result <- new_sclet_ai_result(
        answer = "test report",
        findings = list(list(
            statement = "the object has four cells",
            claim_level = "hypothesis",
            evidence_refs = character()
        )),
        task = "interpret"
    )
    RecordAIResult(sce, result, id = "ai_round_test")
}

test_that("a recorded round renders and its self-check passes", {
    sce <- m12e_recorded_sce()
    rendered <- sclet:::sclet_ai_render_round_report(sce, id = "ai_round_test")
    expect_identical(rendered$status, "passed")
    expect_true(!is.null(rendered$markdown))
    expect_true(nzchar(rendered$markdown))
    expect_true(isTRUE(rendered$self_check$ok))
    expect_null(rendered$message)
})

test_that("an unrecorded round reports not_rendered instead of passing", {
    # Nothing was recorded, so there is no report. That is not a pass, and the
    # reason has to be readable: the previous behaviour was silence.
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrix(1L, nrow = 4L, ncol = 4L))
    )
    rendered <- sclet:::sclet_ai_render_round_report(sce, id = "ai_missing")
    expect_identical(rendered$status, "not_rendered")
    expect_true(grepl("no recorded AI result", rendered$message))
    expect_null(rendered$markdown)
})

test_that("disabled rendering says so rather than pretending it passed", {
    sce <- m12e_recorded_sce()
    rendered <- sclet:::sclet_ai_render_round_report(sce, id = "ai_round_test",
        enabled = FALSE)
    expect_identical(rendered$status, "disabled")
    expect_null(rendered$markdown)
})

test_that("a render failure never raises, so an execution that ran is not discarded", {
    # The seam is called with something that cannot be rendered at all. It must
    # return a status, not throw: the analysis is the work, the report is a
    # view of it.
    rendered <- sclet:::sclet_ai_render_round_report(
        SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1L, 2L, 2L))),
        id = "definitely_not_here"
    )
    expect_true(rendered$status %in% c("not_rendered", "failed"))
    expect_true(is.character(rendered$message))
})

test_that("RunAIAnalysis exposes the two switches and documents their defaults", {
    args <- formals(RunAIAnalysis)
    expect_identical(args$render_report, TRUE)
    expect_identical(args$verify_before_continue, TRUE)
})
