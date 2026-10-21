# M12-c: validation errors must say how to fix the plan, not only that it is wrong.
#
# The measured failure shapes this pins are the knowledge-type ones from the
# 2026-10-19 Kang reliability repeats (9 trials, 4 failures, all correctly
# refused by the validator):
#
#   1. a zero-action plan for a 43k-cell goal (`plan_empty`), which said only
#      "re-propose the plan" and so recurred;
#   2. a step whose `depends_on` still named a step the redundant-DE dedupe had
#      removed (`depends on unknown step(s)`), 2 of 9 trials -- the error named
#      the unknown id but neither of the two possible repairs.
#
# The rule these tests enforce: every refusal of a *shape* the model can fix
# itself carries a `fix:` clause naming what to change. The verdicts themselves
# are unchanged -- this is message quality, not a relaxed or tightened gate, so
# each case asserts both the verdict AND the guidance.

m12c_fixture <- function() {
    SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 3, ncol = 2))
    )
}

test_that("plan_empty states the floor a plan must clear, and still refuses", {
    sce <- m12c_fixture()
    registry <- AIDefaultExecutionRegistry(sce, include = "read")
    plan <- new_sclet_ai_plan(
        task = "m12c_empty",
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
        actions = list()
    )
    validation <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_false(validation$valid)

    error <- paste(validation$errors, collapse = " ")
    expect_true(grepl("plan_empty", error, fixed = TRUE))
    # the guidance must name the repair, not merely restate the failure
    expect_true(grepl("fix:", error, fixed = TRUE))
    expect_true(grepl("success_criteria", error, fixed = TRUE))
    expect_true(grepl("zero steps is not a way", error, fixed = TRUE))
})

test_that("a dependency on a removed step names both repairs it could take", {
    sce <- m12c_fixture()
    registry <- AIDefaultExecutionRegistry(sce, include = "all")
    plan <- new_sclet_ai_plan(
        task = "m12c_dangling_dep",
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
        actions = list(
            list(
                id = "step_1",
                action = "normalize_data",
                params = list(),
                depends_on = list()
            ),
            list(
                id = "step_2",
                action = "run_pca",
                params = list(ncomponents = 20L),
                # step_9 is not in the plan: the shape the dedupe produced
                depends_on = list("step_1", "step_9")
            )
        )
    )
    validation <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_false(validation$valid)

    dep_errors <- grep("depends on unknown step", validation$errors, value = TRUE)
    expect_true(length(dep_errors) >= 1L)
    dep_error <- paste(dep_errors, collapse = " ")
    # the unknown id must still be named, or the model cannot find it
    expect_true(grepl("step_9", dep_error, fixed = TRUE))
    # and both repairs must be stated
    expect_true(grepl("fix:", dep_error, fixed = TRUE))
    expect_true(grepl("remove those ids", dep_error, fixed = TRUE))
    expect_true(grepl("add the step back", dep_error, fixed = TRUE))
})

test_that("an unregistered action keeps its do-not-use-a-placeholder guidance", {
    sce <- m12c_fixture()
    registry <- AIDefaultExecutionRegistry(sce, include = "read")
    plan <- new_sclet_ai_plan(
        task = "m12c_not_run",
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
        actions = list(list(
            id = "step_1",
            action = "not_run",
            params = list(),
            depends_on = list()
        ))
    )
    validation <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_false(validation$valid)

    error <- paste(validation$errors, collapse = " ")
    expect_true(grepl("action is not registered: not_run", error, fixed = TRUE))
    # the 2026-10-19 failure was exactly this shape: a placeholder action used
    # to "drop" a step. The message must say to remove the step instead.
    expect_true(grepl("remove the step from the plan instead", error, fixed = TRUE))
})

test_that("an unknown step reference states how to repair it", {
    sce <- m12c_fixture()
    registry <- AIDefaultExecutionRegistry(sce, include = "read")
    plan <- new_sclet_ai_plan(
        task = "m12c_dangling_ref",
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
        actions = list(
            list(
                id = "step_1",
                action = "run_pca",
                params = list(ncomponents = 20L),
                depends_on = list()
            ),
            list(
                id = "step_2",
                action = "run_umap",
                params = list(reduction = "${step_7.output.reduction}"),
                depends_on = list("step_1")
            )
        )
    )
    validation <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_false(validation$valid)

    ref_errors <- grep("references unknown step", validation$errors, value = TRUE)
    expect_true(length(ref_errors) >= 1L)
    ref_error <- paste(ref_errors, collapse = " ")
    expect_true(grepl("step_7", ref_error, fixed = TRUE))
    expect_true(grepl("fix:", ref_error, fixed = TRUE))
    expect_true(grepl("pass a literal value", ref_error, fixed = TRUE))
})

# ---- M12-a follow-up: the contrast-level gap ---------------------------------
# The only failure in the 12 same-tree trials (commit 5947e4b) was
# `contrast_invalid: ident.1/ident.2 must be levels of 'condition': ctrl, stim`.
# Root cause was not the validator: the planning prompt sent group-valued
# parameters to `capabilities$identity_levels$codes`, which describes the
# CLUSTERING, while a `group` contrast takes levels of a colData column -- a set
# that had no representation in the context at all. These tests pin the fix.

test_that("contrast_levels exposes colData levels separately from identity levels", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrix(1L, nrow = 3L, ncol = 8L)),
        colData = S4Vectors::DataFrame(
            condition = rep(c("ctrl", "stim"), each = 4L),
            donor = rep(paste0("d", 1:4), times = 2L)
        )
    )
    caps <- sclet:::sclet_ai_capabilities(
        list(available = list(assays = "counts")),
        list(analyses = list(), workflows = list()),
        object = sce
    )
    expect_true(!is.null(caps$contrast_levels))
    columns <- vapply(caps$contrast_levels$columns, function(x) x$column, character(1))
    expect_true("condition" %in% columns)
    expect_true("donor" %in% columns)
    expect_true("identity_levels" %in% names(caps))
    expect_false(isTRUE(caps$contrast_levels$raw_values_included))

    condition <- Filter(function(x) identical(x$column, "condition"),
        caps$contrast_levels$columns)[[1L]]
    expect_identical(condition$n_levels, 2L)
    expect_setequal(unlist(condition$levels), c("ctrl", "stim"))
    expect_true(condition$levels_available)
})

test_that("a column with too many levels reports its count but withholds the values", {
    # A donor column with many levels is a design role, not a usable contrast.
    # The count stays visible so "withheld" is distinguishable from "empty".
    n <- 20L
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrix(1L, nrow = 2L, ncol = n)),
        colData = S4Vectors::DataFrame(subject = paste0("s", seq_len(n)))
    )
    levels_info <- sclet:::sclet_ai_build_contrast_levels(sce, max_levels = 12L)
    subject <- Filter(function(x) identical(x$column, "subject"),
        levels_info$columns)[[1L]]
    expect_identical(subject$n_levels, as.integer(n))
    expect_false(subject$levels_available)
    # `character(0)`, not NULL: the field is always a bounded vector so the
    # privacy layer passes it through, and `levels_available = FALSE` is what
    # says "withheld" rather than "none".
    expect_length(subject$levels, 0L)
})

test_that("the planning prompt names where group contrasts get their levels", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1L, nrow = 3L, ncol = 4L))
    )
    captured <- new.env(parent = emptyenv())
    fake_call <- function(task, context, ..., system_prompt = NULL) {
        captured$system_prompt <- system_prompt
        response <- list(
            answer = "mock",
            proposed_actions = list(list(id = "step1", action = "inspect_status")),
            success_criteria = list(),
            stop_conditions = list()
        )
        sclet:::sclet_ai_normalize_response(response, task = task,
            context = context, metadata = list(mock = TRUE))
    }
    testthat::with_mocked_bindings(
        {
            invisible(AIPlanAnalysis(sce, goal = "compare stimulated vs control", model = "mock"))
            prompt <- captured$system_prompt %||% ""
            # the old instruction pointed every group-valued parameter at the
            # clustering codes; the contrast case must be stated separately
            expect_true(grepl("contrast_levels", prompt, fixed = TRUE))
            expect_true(grepl("WITH `group = <column>`", prompt, fixed = TRUE))
            expect_true(grepl("WITHOUT a `group` parameter", prompt, fixed = TRUE))
        },
        sclet_ai_call = fake_call,
        .package = "sclet"
    )
})

test_that("contrast_invalid tells the model to copy a listed level", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1L, nrow = 3L, ncol = 8L))
    )
    SummarizedExperiment::colData(sce)$condition <- rep(c("ctrl", "stim"), each = 4L)
    registry <- AIDefaultExecutionRegistry(sce, include = "all")
    de <- registry[["run_de_test"]]
    expect_true(!is.null(de))
    result <- de$prerequisites(sce, list(
        group = "condition", ident.1 = "treat", ident.2 = "ctrl"
    ))
    expect_true(is.character(result))
    expect_true(grepl("contrast_invalid", result, fixed = TRUE))
    # the valid levels are named, the wrong one is named, and the repair says
    # "copy one of these" rather than leaving the model to guess
    expect_true(grepl("ctrl", result, fixed = TRUE))
    expect_true(grepl("stim", result, fixed = TRUE))
    expect_true(grepl("you wrote: treat", result, fixed = TRUE))
    expect_true(grepl("fix:", result, fixed = TRUE))
    expect_true(grepl("contrast_levels", result, fixed = TRUE))
})
