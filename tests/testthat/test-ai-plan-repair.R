## Tests for AIPlanAnalysis() repair loop (spec: .dev/spec-ai-plan-repair.md)
##
## Mock-invocation note:
## ---------------------
## The `options(sclet.ai.call = ...)` mock branch in `sclet_ai_call()` is
## invoked as `mock(task, context, tools, schema, model, ...)` and returns
## before the real path assembles `prompt`. The mock therefore cannot
## observe the `system_prompt` (repair or otherwise). Scenarios 1-9
## distinguish first-call from repair-call by **invocation count** (the
## mock returns invalid on call 1, valid on call 2, etc.). Prompt
## *content* is covered by scenario 10, which tests the pure
## `sclet_ai_repair_prompt_block()` helper directly.

test_that("first proposal valid -> no repair (scenario 1)", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1L, nrow = 3L, ncol = 2L))
    )
    calls <- 0L
    old <- options(sclet.ai.call = function(task, context, ...) {
        calls <<- calls + 1L
        list(
            answer = "Inspect status only.",
            proposed_actions = list(list(
                id = "status",
                action = "inspect_status",
                params = list(),
                depends_on = character()
            ))
        )
    })
    on.exit(options(old), add = TRUE)

    plan <- AIPlanAnalysis(sce)
    expect_s3_class(plan, "sclet_ai_plan")
    expect_true(isTRUE(plan$validation$valid))
    expect_equal(plan$metadata$repair_rounds_used, 0L)
    expect_equal(plan$metadata$validation_error_counts, 0L)
    # A first attempt that was already valid has no first-attempt errors, but
    # the field must still exist so callers can rely on its shape.
    expect_identical(plan$metadata$first_attempt_errors, character())
    expect_equal(calls, 1L)
})

test_that("invalid then valid -> exactly one repair (scenario 2)", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1L, nrow = 3L, ncol = 2L))
    )
    calls <- 0L
    old <- options(sclet.ai.call = function(task, context, ...) {
        calls <<- calls + 1L
        if (calls == 1L) {
            list(
                answer = "Bad plan.",
                proposed_actions = list(list(
                    id = "bad",
                    action = "not_a_real_action",
                    params = list(),
                    depends_on = character()
                ))
            )
        } else {
            list(
                answer = "Fixed plan.",
                proposed_actions = list(list(
                    id = "status",
                    action = "inspect_status",
                    params = list(),
                    depends_on = character()
                ))
            )
        }
    })
    on.exit(options(old), add = TRUE)

    plan <- AIPlanAnalysis(sce)
    expect_true(isTRUE(plan$validation$valid))
    expect_equal(plan$metadata$repair_rounds_used, 1L)
    expect_equal(length(plan$metadata$validation_error_counts), 2L)
    expect_true(plan$metadata$validation_error_counts[1L] > 0L)
    expect_equal(plan$metadata$validation_error_counts[2L], 0L)
    # The first attempt's error text survives repair, and its length matches
    # the recorded first count -- the two views of the same fact must agree.
    expect_length(
        plan$metadata$first_attempt_errors,
        plan$metadata$validation_error_counts[1L]
    )
    expect_true(all(nzchar(plan$metadata$first_attempt_errors)))
    # ...and it is the *first* attempt's text, not the repaired plan's.
    expect_identical(plan$validation$errors, character())
    expect_equal(calls, 2L)
})

test_that("max_repair_rounds = 0L -> no repair attempted (scenario 3)", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1L, nrow = 3L, ncol = 2L))
    )
    calls <- 0L
    old <- options(sclet.ai.call = function(task, context, ...) {
        calls <<- calls + 1L
        list(
            answer = "Bad plan.",
            proposed_actions = list(list(
                id = "bad",
                action = "not_a_real_action",
                params = list(),
                depends_on = character()
            ))
        )
    })
    on.exit(options(old), add = TRUE)

    plan <- AIPlanAnalysis(sce, max_repair_rounds = 0L)
    expect_equal(plan$metadata$repair_rounds_used, 0L)
    expect_equal(calls, 1L)
    expect_false(isTRUE(plan$validation$valid))
})

test_that("never improves -> stops early after 2 calls, not 4 (scenario 4)", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1L, nrow = 3L, ncol = 2L))
    )
    calls <- 0L
    old <- options(sclet.ai.call = function(task, context, ...) {
        calls <<- calls + 1L
        list(
            answer = "Same bad plan every time.",
            proposed_actions = list(list(
                id = "bad",
                action = "not_a_real_action",
                params = list(),
                depends_on = character()
            ))
        )
    })
    on.exit(options(old), add = TRUE)

    plan <- AIPlanAnalysis(sce, max_repair_rounds = 3L)
    expect_equal(calls, 2L)
    expect_equal(plan$metadata$repair_rounds_used, 1L)
    expect_equal(plan$metadata$validation_error_counts, c(1L, 1L))
})

test_that("returns best attempt, not last (scenario 5)", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1L, nrow = 3L, ncol = 2L))
    )
    calls <- 0L
    bad3 <- function() list(
        answer = "3 errors.",
        proposed_actions = list(
            list(id = "a", action = "fake_a", params = list(), depends_on = character()),
            list(id = "b", action = "fake_b", params = list(), depends_on = character()),
            list(id = "c", action = "fake_c", params = list(), depends_on = character())
        )
    )
    bad1 <- function() list(
        answer = "1 error.",
        proposed_actions = list(list(
            id = "only",
            action = "fake_only",
            params = list(),
            depends_on = character()
        ))
    )
    bad5 <- function() list(
        answer = "5 errors.",
        proposed_actions = list(
            list(id = "a", action = "fake_a", params = list(), depends_on = character()),
            list(id = "b", action = "fake_b", params = list(), depends_on = character()),
            list(id = "c", action = "fake_c", params = list(), depends_on = character()),
            list(id = "d", action = "fake_d", params = list(), depends_on = character()),
            list(id = "e", action = "fake_e", params = list(), depends_on = character())
        )
    )
    old <- options(sclet.ai.call = function(task, context, ...) {
        calls <<- calls + 1L
        switch(
            as.character(calls),
            "1" = bad3(),
            "2" = bad1(),
            "3" = bad5(),
            bad5()
        )
    })
    on.exit(options(old), add = TRUE)

    plan <- AIPlanAnalysis(sce, max_repair_rounds = 3L)
    expect_equal(length(plan$actions), 1L)
    expect_equal(plan$actions[[1L]]$action, "fake_only")
    expect_equal(plan$metadata$validation_error_counts, c(3L, 1L, 5L))
    expect_equal(plan$metadata$repair_rounds_used, 2L)
})

test_that("repair executes nothing (scenario 6)", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1L, nrow = 3L, ncol = 2L))
    )
    calls <- 0L
    old <- options(sclet.ai.call = function(task, context, ...) {
        calls <<- calls + 1L
        if (calls == 1L) {
            list(
                answer = "Bad.",
                proposed_actions = list(list(
                    id = "bad", action = "not_a_real_action",
                    params = list(), depends_on = character()
                ))
            )
        } else {
            list(
                answer = "Good.",
                proposed_actions = list(list(
                    id = "status", action = "inspect_status",
                    params = list(), depends_on = character()
                ))
            )
        }
    })
    on.exit(options(old), add = TRUE)

    ledger_before <- GetAnalysisLedger(sce)
    plan <- AIPlanAnalysis(sce)
    ledger_after <- GetAnalysisLedger(sce)
    expect_true(isTRUE(plan$validation$valid))
    # Falsifiable "nothing was executed" check. The ledger fingerprint is
    # stable across repeated reads of an unmodified object (verified: two
    # reads 1.2s apart give the same value), but changes as soon as an
    # analysis or state record is written, so a real execution during
    # planning would fail this assertion.
    #
    # This replaces an earlier assertion that compared
    # `length(ledger_$execution$records)` before and after. GetAnalysisLedger()
    # has no `execution` field at all, so that comparison was NULL == NULL --
    # always true and incapable of failing.
    expect_identical(ledger_after$fingerprint, ledger_before$fingerprint)
    expect_equal(length(ledger_after$analyses), 0L)
    expect_equal(length(ledger_after$state_records), 0L)
    expect_null(plan$ai_result$execution)
})

test_that("out-of-range / non-integer max_repair_rounds errors (scenario 7)", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1L, nrow = 3L, ncol = 2L))
    )
    expect_error(AIPlanAnalysis(sce, max_repair_rounds = -1L),
                 class = "sclet_ai_invalid_arg")
    expect_error(AIPlanAnalysis(sce, max_repair_rounds = 99L),
                 class = "sclet_ai_invalid_arg")
    expect_error(AIPlanAnalysis(sce, max_repair_rounds = NA),
                 class = "sclet_ai_invalid_arg")
    expect_error(AIPlanAnalysis(sce, max_repair_rounds = "two"),
                 class = "sclet_ai_invalid_arg")
    expect_error(AIPlanAnalysis(sce, max_repair_rounds = 1.5),
                 class = "sclet_ai_invalid_arg")
    expect_error(AIPlanAnalysis(sce, max_repair_rounds = NULL),
                 class = "sclet_ai_invalid_arg")
    expect_error(AIPlanAnalysis(sce, max_repair_rounds = c(1L, 2L)),
                 class = "sclet_ai_invalid_arg")
    expect_error(AIPlanAnalysis(sce, max_repair_rounds = TRUE),
                 class = "sclet_ai_invalid_arg")
    # Non-finite values are the case that originally slipped through: the
    # range predicate compared against as.integer(Inf), which both warned
    # ("NAs introduced by coercion to integer range") and made the enclosing
    # `if` condition NA, producing a bare simpleError instead of the typed
    # sclet_ai_invalid_arg. Guard the warning explicitly, not just the class.
    expect_error(AIPlanAnalysis(sce, max_repair_rounds = Inf),
                 class = "sclet_ai_invalid_arg")
    expect_error(AIPlanAnalysis(sce, max_repair_rounds = -Inf),
                 class = "sclet_ai_invalid_arg")
    expect_error(AIPlanAnalysis(sce, max_repair_rounds = NaN),
                 class = "sclet_ai_invalid_arg")
    expect_warning(
        tryCatch(AIPlanAnalysis(sce, max_repair_rounds = Inf),
                 error = function(e) NULL),
        NA
    )
})

test_that("no-repair path keeps full plan shape (scenario 8)", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1L, nrow = 3L, ncol = 2L))
    )
    old <- options(sclet.ai.call = function(task, context, ...) {
        list(
            answer = "Inspect.",
            proposed_actions = list(list(
                id = "status", action = "inspect_status",
                params = list(), depends_on = character()
            ))
        )
    })
    on.exit(options(old), add = TRUE)

    plan <- AIPlanAnalysis(sce, max_repair_rounds = 0L)
    required <- c("plan_id", "task", "rationale", "actions", "validation",
                  "metadata", "assumptions", "candidate_routes", "selected_route",
                  "success_criteria", "risks", "stop_conditions")
    for (field in required) {
        expect_true(field %in% names(plan), info = paste("missing field:", field))
    }
    expect_s3_class(plan, "sclet_ai_plan")
    expect_true("repair_rounds_used" %in% names(plan$metadata))
    expect_true("validation_error_counts" %in% names(plan$metadata))
})

## Scenario 9 ("the full ai-* suite still passes") is deliberately NOT
## implemented as an in-file test.
##
## It was originally written as a test_that() that looped over every other
## test-ai-*.R file and ran it via testthat::test_file(). That is the wrong
## shape for this guard, for three reasons:
##   1. It re-ran the entire suite (including real scVelo / scDoubletFinder /
##      marker-test fixtures) inside a single mock-only unit-test file, making
##      the file take ~150s instead of ~1s, and that cost grows with the suite.
##   2. It downgraded a file-level *error* to `skip()`, so a genuine regression
##      that threw rather than failed an expectation would be reported as a
##      skip, not a failure.
##   3. Nested test_file() execution shares the same process, options and RNG
##      state as the outer run, so it can mask or create cross-test
##      interference that would not occur in a normal run.
##
## The full-suite regression guard is instead performed by running the suite
## itself (`devtools::test(filter = "ai-")`, or the equivalent per-file loop),
## which is what the spec's acceptance checklist calls for. See
## .dev/spec-ai-plan-repair.md section 4 scenario 9 and section 5.

test_that("sclet_ai_repair_prompt_block content (scenario 10)", {
    block <- sclet:::sclet_ai_repair_prompt_block(c(
        "action is not registered: fake_action",
        "parameter ncomponents has invalid type"
    ))
    expect_true(is.character(block))
    expect_true(length(block) > 0L)
    joined <- paste(block, collapse = "\n")

    expect_true(any(grepl("action is not registered: fake_action",
                          block, fixed = TRUE)))
    expect_true(any(grepl("parameter ncomponents has invalid type",
                          block, fixed = TRUE)))

    expect_true(grepl("rejected", joined, ignore.case = TRUE))
    expect_true(grepl("not executed|NOT executed|not run", joined, ignore.case = TRUE))

    expect_true(grepl("Do NOT invent|do not invent|Do not invent", joined))
    expect_true(grepl("DELETE|delete|omit|remove", joined))

    expect_true(grepl("input_schema", joined))
    expect_true(grepl("type|TYPE", joined))

    empty_null <- sclet:::sclet_ai_repair_prompt_block(NULL)
    expect_true(length(empty_null) == 0L)
    empty_char <- sclet:::sclet_ai_repair_prompt_block(character(0))
    expect_true(length(empty_char) == 0L)
})
