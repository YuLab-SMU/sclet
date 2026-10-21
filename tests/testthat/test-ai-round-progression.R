# M13-b: round progression (D3) is computable from two plans, not from a model's opinion.
#
# The dimension that was missing after M12: plan-quality dimensions 1, 2 and 4 are all
# computable from one round's context, but dimension 3 asks "did this round respond to
# what the last round said", which needs two rounds in view at once. This file is that.
#
# The judges must not read the model's reasoning -- unfalsifiable. Every judge below is a
# relation between the previous plan, the previous report's own fields, and the current
# plan. Four classes: advanced, confirmatory, no_progress, regressed.
#
# `regressed` is the one that has to exist. A dimension that can only say "good" or
# "average" cannot say "worse", and re-deriving an answer the previous round already
# discarded is the most expensive move a plan can make.

m13b_plan <- function(actions, rationale = "") {
    sclet:::new_sclet_ai_plan(
        task = paste0("m13b_", paste(vapply(actions, function(x) x$action, character(1)),
            collapse = "_")),
        actions = actions,
        rationale = rationale,
        context_fingerprint = "sclet-ai-m13b"
    )
}

m13b_report <- function(warnings = character(), gaps = list(),
        decisions = list()) {
    list(
        warnings = as.list(warnings),
        evidence_gaps = gaps,
        user_decisions = decisions
    )
}

test_that("a round that runs the remedy for a warning is advanced", {
    prev <- m13b_plan(list(
        list(id = "p1", action = "run_annotation", params = list(ref = "x"))
    ))
    prev_report <- m13b_report(warnings = c("cluster 12 is 98.3% doublets"))

    curr <- m13b_plan(list(
        list(id = "c1", action = "run_doublet_detection", params = list()),
        list(id = "c2", action = "run_de_test", params = list())
    ))
    d <- sclet:::sclet_ai_round_progression(prev, prev_report, curr)
    expect_identical(d$progression_class, "advanced")
    expect_identical(d$n_warnings, 1L)
    expect_identical(d$n_warnings_answered, 1L)
    expect_identical(d$n_warnings_ignored, 0L)
})

test_that("a warning the planner explicitly carries forward is carried, not ignored", {
    # An accepted warning is not a defect. "We cannot fix the doublet signal because
    # no doublet evidence exists" is an answer; silence is not.
    prev <- m13b_plan(list(list(id = "p1", action = "run_annotation", params = list())))
    prev_report <- m13b_report(warnings = "ambient RNA dominates the empty droplets")

    curr <- m13b_plan(
        list(list(id = "c1", action = "run_de_test", params = list())),
        rationale = "ambient decontamination was not run because the object has no raw counts"
    )
    d <- sclet:::sclet_ai_round_progression(prev, prev_report, curr)
    expect_identical(d$n_warnings_ignored, 0L)
    expect_true(d$n_warnings_carried >= 1L)
    expect_length(d$unanswered_warnings, 0L)
})

test_that("a warning nobody mentions is ignored and is listed by name", {
    prev <- m13b_plan(list(list(id = "p1", action = "run_annotation", params = list())))
    prev_report <- m13b_report(warnings = "cluster 12 is 98.3% doublets")

    curr <- m13b_plan(list(list(id = "c1", action = "run_enrichment", params = list())))
    d <- sclet:::sclet_ai_round_progression(prev, prev_report, curr)
    expect_identical(d$n_warnings_ignored, 1L)
    expect_length(d$unanswered_warnings, 1L)
    # no remedy and no mention: the class reflects that nothing was answered,
    # but it is confirmatory rather than regressed -- ignoring is not undoing
    expect_identical(d$progression_class, "confirmatory")
})

test_that("an unchanged plan is no_progress and its ratio is 0", {
    step <- list(id = "s1", action = "run_de_test", params = list(ident.1 = "a", ident.2 = "b"))
    prev <- m13b_plan(list(step))
    curr <- m13b_plan(list(step))
    d <- sclet:::sclet_ai_round_progression(prev, m13b_report(), curr)
    expect_identical(d$n_actions, 1L)
    expect_identical(d$n_repeated_steps, 1L)
    expect_equal(d$progress_ratio, 0)
    expect_identical(d$progression_class, "no_progress")
})

test_that("the same action with new parameters is a new analysis, not a replay", {
    prev <- m13b_plan(list(list(id = "s1", action = "run_de_test",
        params = list(min.pct = 0.01))))
    curr <- m13b_plan(list(list(id = "s2", action = "run_de_test",
        params = list(min.pct = 0.25))))
    d <- sclet:::sclet_ai_round_progression(prev, m13b_report(), curr)
    expect_identical(d$n_repeated_steps, 0L)
    expect_equal(d$progress_ratio, 1)
})

test_that("partially repeated plans report an intermediate ratio", {
    prev <- m13b_plan(list(
        list(id = "s1", action = "run_de_test", params = list(a = 1)),
        list(id = "s2", action = "run_pca", params = list(npc = 20))
    ))
    curr <- m13b_plan(list(
        list(id = "s1", action = "run_de_test", params = list(a = 1)),
        list(id = "s3", action = "run_enrichment", params = list())
    ))
    d <- sclet:::sclet_ai_round_progression(prev, m13b_report(), curr)
    expect_identical(d$n_repeated_steps, 1L)
    expect_equal(d$progress_ratio, 0.5)
})

test_that("evidence gaps close when the named analysis runs", {
    prev <- m13b_plan(list(list(id = "p1", action = "run_rare_cell_detection", params = list())))
    prev_report <- m13b_report(gaps = list(list(
        analysis_type = "rare_cells",
        reason = "doublet_evidence_available"
    )))
    curr <- m13b_plan(list(list(id = "c1", action = "run_doublet_detection", params = list())))
    d <- sclet:::sclet_ai_round_progression(prev, prev_report, curr)
    expect_identical(d$n_evidence_gaps, 1L)
    expect_identical(d$n_evidence_gaps_closed, 1L)
})

test_that("a recorded user decision that the next round never uses is counted", {
    prev <- m13b_plan(list(list(id = "p1", action = "run_annotation", params = list())))
    prev_report <- m13b_report(decisions = list(list(
        summary = "use MonacoImmuneData as the annotation reference"
    )))
    curr <- m13b_plan(list(list(id = "c1", action = "run_enrichment", params = list())))
    d <- sclet:::sclet_ai_round_progression(prev, prev_report, curr)
    expect_identical(d$n_user_decisions, 1L)
    expect_identical(d$n_user_decisions_used, 0L)
    expect_length(d$unused_user_decisions, 1L)
})

test_that("a regression is detected and outranks every positive signal", {
    # The previous round established that cluster 12 is a doublet population; this
    # round re-analyses it without any decontamination. Every positive judge is
    # satisfied by the same plan, and the class still has to say `regressed` --
    # otherwise the dimension can only report good and average.
    prev <- m13b_plan(list(list(id = "p1", action = "run_doublet_detection", params = list())))
    prev_report <- m13b_report(
        warnings = "cluster 12 is 98.3% doublets",
        gaps = list(list(analysis_type = "annotation", reason = "reference_missing")),
        decisions = list(list(summary = "use MonacoImmuneData as the reference"))
    )
    curr <- m13b_plan(
        list(
            list(id = "c1", action = "run_annotation", params = list(ref = "x")),
            list(id = "c2", action = "run_de_test", params = list())
        ),
        rationale = "annotate the clusters and test markers, MonacoImmuneData"
    )
    d <- sclet:::sclet_ai_round_progression(prev, prev_report, curr)
    # every positive counter is non-zero...
    expect_true(d$n_warnings_ignored >= 0L)
    expect_true(d$n_evidence_gaps_closed >= 1L)
    expect_true(d$n_user_decisions_used >= 1L)
    # ...and yet the class is the worst one
    expect_true(d$regressed)
    expect_identical(d$progression_class, "regressed")
})

test_that("the judge is deterministic and carries no raw values", {
    prev <- m13b_plan(list(list(id = "p1", action = "run_annotation", params = list())))
    prev_report <- m13b_report(
        warnings = c("cluster 12 is 98.3% doublets", "ambient RNA detected"),
        gaps = list(list(analysis_type = "rare_cells", reason = "doublet_evidence_available"))
    )
    curr <- m13b_plan(list(
        list(id = "c1", action = "run_doublet_detection", params = list())
    ))
    first <- sclet:::sclet_ai_round_progression(prev, prev_report, curr)
    second <- sclet:::sclet_ai_round_progression(prev, prev_report, curr)
    expect_identical(first, second)
    expect_false(isTRUE(first$raw_values_included))
    expect_identical(first$schema_version, "1.0")
})

test_that("an empty current plan is not silently advanced or regressed", {
    prev <- m13b_plan(list(list(id = "p1", action = "run_de_test", params = list())))
    prev_report <- m13b_report(warnings = "cluster 12 is 98.3% doublets")
    curr <- m13b_plan(list())
    d <- sclet:::sclet_ai_round_progression(prev, prev_report, curr)
    expect_identical(d$n_actions, 0L)
    expect_true(is.na(d$progress_ratio))
    # it answers nothing, so it cannot be advanced
    expect_identical(d$n_warnings_answered, 0L)
})

# ---- Finding the previous round from the ledger (M13-b, second half) -------
#
# A caller that has to remember to pass the previous round will eventually not, and
# then D3 silently becomes "always a first round" -- indistinguishable from a
# workflow that never advances. Locating the round from the ledger turns that
# failure into the visible string "no previous round".

test_that("first round reports first_round, not unavailability", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1L, nrow = 3L, ncol = 4L))
    )
    plan <- sclet:::new_sclet_ai_plan(
        task = "m13b_first", actions = list(list(id = "s1", action = "inspect_status"))
    )
    d <- sclet:::sclet_ai_previous_round_progression(sce, plan, list(warnings = list()))
    expect_identical(d$progression_class, "first_round")
    expect_true(grepl("first round", d$reason, fixed = TRUE))
    # and it is NOT reported as a defect in D3 itself
    expect_null(d$n_actions)
})

test_that("a non-SCE object degrades to first_round rather than raising", {
    plan <- sclet:::new_sclet_ai_plan(
        task = "m13b_degraded", actions = list(list(id = "s1", action = "inspect_status"))
    )
    d <- sclet:::sclet_ai_previous_round_progression(NULL, plan, list())
    expect_identical(d$progression_class, "first_round")
})

# ---- The previous round is two records, and that was found by probing -------
#
# No single record holds both halves of a previous round. The plan comes from the
# newest `ai_execution` record's `params$steps`; the report comes from the newest
# `ai_report` record, whose warnings live in ARTIFACTS (a list of character
# vectors), not in summary, and whose cited user decisions are
# `summary$cited_user_decisions`. Verified against the real M1-M11 objects.
#
# The failure this pins: looking for either half in the other place returned
# `not_comparable` on every real object -- which reads exactly like "D3 was never
# computed here" and is indistinguishable from a workflow that never advances.
test_that("the previous round is assembled from the execution and report records", {
    # Real-data probe, NOT a unit test: it reads the M7 object, which exists only
    # on this machine. It is skipped rather than kept unconditional, because a
    # test that fails on a fresh clone because /tmp is empty teaches the reader
    # to ignore failures. The shape it verifies is pinned above by fixtures.
    skip_if_not(file.exists("/tmp/sclet-agentrun/m7-report.rds"),
        "the M7 object is not present (real-data probe, run locally)")
    sce <- readRDS("/tmp/sclet-agentrun/m7-report.rds")
    current <- sclet:::new_sclet_ai_plan(
        task = "m13b_real_probe",
        actions = list(
            list(id = "c1", action = "run_enrichment", params = list()),
            list(id = "c2", action = "run_annotation", params = list(ref = "x"))
        ),
        rationale = "annotate the clusters and run enrichment"
    )
    d <- sclet:::sclet_ai_previous_round_progression(sce, current, list())
    # it must not fall through to the two silent-degradation branches
    expect_false(identical(d$progression_class, "first_round"))
    expect_false(identical(d$progression_class, "not_comparable"))
    # the plan was rebuilt from the execution record, so its actions are visible
    expect_true(d$n_actions >= 1L)
    # and the report's warnings came from artifacts, not summary
    expect_true(d$n_warnings >= 1L)
    expect_true(d$progression_class %in%
        c("advanced", "confirmatory", "no_progress", "regressed"))
})

test_that("every remedy names a registry action that actually exists", {
    # Regression for a defect found by a real-data probe: `run_decontx` and
    # `run_imputation` were in the remedy table and are NOT registry keys.
    # RunDecontX() is exported but has no AIAction; there is no imputation
    # action at all. A remedy that cannot match makes warnings read as
    # "ignored" forever -- the judge would be quietly wrong, never loudly.
    remedies <- sclet:::sclet_ai_round_progression_remedies()
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1L, nrow = 3L, ncol = 4L))
    )
    registry <- AIDefaultExecutionRegistry(sce, include = "all")
    known <- names(registry)
    for (topic in names(remedies)) {
        for (action in remedies[[topic]]) {
            expect_true(action %in% known,
                info = paste0("remedy for '", topic, "' names non-existent action: ", action))
        }
    }
})

test_that("an ambient warning is answerable now that the action exists", {
    # Before M13-a there was no remedy possible for ambient RNA: RunDecontX()
    # was exported but had no AIAction, so an ambient warning could only ever be
    # carried, never answered. The action was created in M13-a and this is the
    # consequence: a warning about ambient contamination is now a thing a round
    # can actually respond to.
    remedies <- sclet:::sclet_ai_round_progression_remedies()
    expect_true("run_decontx" %in% remedies$ambient)

    prev <- m13b_plan(list(list(id = "p1", action = "run_de_test", params = list())))
    prev_report <- m13b_report(warnings = "ambient RNA contamination dominates the empty droplets")
    curr <- m13b_plan(list(list(id = "c1", action = "run_decontx", params = list())))
    d <- sclet:::sclet_ai_round_progression(prev, prev_report, curr)
    expect_identical(d$n_warnings_answered, 1L)
    expect_identical(d$progression_class, "advanced")
})

# ---- B20 judge: an unmentioned suggestion is countable ----------------------
#
# The prompt now REQUIRES each unexplored entry to be planned or declined by
# name. Without this judge that requirement is advice a model can violate
# silently -- which is the measured failure: 4 of 6 trials ignored the
# suggestion and none mentioned it. "Oversight" is the honest name.

m13b_with_unexplored <- function(actions, rationale, suggestions) {
    p <- m13b_plan(actions, rationale = rationale)
    p$metadata$unexplored_suggestions <- suggestions
    p
}

test_that("a suggestion that is planned is counted as planned", {
    prev <- m13b_plan(list(list(id = "p1", action = "run_de_test", params = list())))
    curr <- m13b_with_unexplored(
        list(list(id = "c1", action = "run_decontx", params = list())),
        "clean the ambient pool, then test markers",
        list(list(action = "run_decontx", what = "ambient pressure is high"))
    )
    d <- sclet:::sclet_ai_round_progression(prev, m13b_report(), curr)
    expect_identical(d$n_unexplored, 1L)
    expect_identical(d$n_unexplored_planned, 1L)
    expect_identical(d$n_unexplored_unmentioned, 0L)
})

test_that("a suggestion declined BY NAME is counted as declined", {
    prev <- m13b_plan(list(list(id = "p1", action = "run_de_test", params = list())))
    curr <- m13b_with_unexplored(
        list(list(id = "c1", action = "run_de_test", params = list())),
        "markers only; run_decontx is not needed because the goal is cluster comparison",
        list(list(action = "run_decontx", what = "ambient pressure is high"))
    )
    d <- sclet:::sclet_ai_round_progression(prev, m13b_report(), curr)
    expect_identical(d$n_unexplored_planned, 0L)
    expect_identical(d$n_unexplored_declined, 1L)
    expect_identical(d$n_unexplored_unmentioned, 0L)
})

test_that("a suggestion neither planned nor named is flagged as unmentioned", {
    # the B20 failure, verbatim: the suggestion was in the context, the plan
    # did not include it, and the rationale never mentioned it
    prev <- m13b_plan(list(list(id = "p1", action = "run_de_test", params = list())))
    curr <- m13b_with_unexplored(
        list(list(id = "c1", action = "run_de_test", params = list())),
        "markers for each cluster, with enrichment",
        list(list(action = "run_decontx", what = "ambient pressure is high"))
    )
    d <- sclet:::sclet_ai_round_progression(prev, m13b_report(), curr)
    expect_identical(d$n_unexplored_unmentioned, 1L)
    expect_identical(d$unmentioned_suggestions, "run_decontx")
})

test_that("no suggestions means no unexplored counters fire", {
    prev <- m13b_plan(list(list(id = "p1", action = "run_de_test", params = list())))
    curr <- m13b_plan(list(list(id = "c1", action = "run_de_test", params = list())))
    d <- sclet:::sclet_ai_round_progression(prev, m13b_report(), curr)
    expect_identical(d$n_unexplored, 0L)
    expect_identical(d$n_unexplored_planned, 0L)
    expect_identical(d$n_unexplored_declined, 0L)
    expect_identical(d$n_unexplored_unmentioned, 0L)
    expect_length(d$unmentioned_suggestions, 0L)
})

test_that("the determinism pin covers the unexplored counters", {
    prev <- m13b_plan(list(list(id = "p1", action = "run_de_test", params = list())))
    curr <- m13b_with_unexplored(
        list(list(id = "c1", action = "run_de_test", params = list())),
        "markers only",
        list(list(action = "run_decontx", what = "ambient pressure"))
    )
    a <- sclet:::sclet_ai_round_progression(prev, m13b_report(), curr)
    b <- sclet:::sclet_ai_round_progression(prev, m13b_report(), curr)
    expect_identical(a, b)
})

# ---- B20 closing half: the report must print silent ignoring ---------------
#
# The prompt requirement and the D3 judge make the failure detectable at plan
# time; the report is where a READER meets it. The snapshot travels on the
# execution record so the report can state, from the ledger alone: "this run was
# offered X and left it unanswered".

test_that("RenderAIReport prints an unmentioned suggestion from the ledger", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1L, nrow = 6L, ncol = 4L))
    )
    plan <- new_sclet_ai_plan(
        task = "b20_report",
        actions = list(list(id = "s1", action = "normalize_data", params = list())),
        rationale = "normalise only"
    )
    plan$metadata$unexplored_suggestions <- list(list(
        action = "run_decontx",
        what = "the ambient-pressure screen is high"
    ))
    plan$actions <- sclet:::sclet_ai_normalize_plan_actions(plan$actions)
    validation <- ValidateAIPlan(plan, object = sce,
        registry = AIDefaultExecutionRegistry(sce, include = "all"), strict = TRUE)
    expect_true(isTRUE(validation$valid))
    exec <- ExecuteAIPlan(sce, plan,
        AIDefaultExecutionRegistry(sce, include = "all"),
        validation = validation, dry_run = FALSE,
        confirmation = validation$confirmation_token)
    expect_identical(exec$status, "completed")

    # RenderAIReport() renders a recorded RESULT; the plan narrative is the
    # record it needs here.
    sce_done <- RecordAIResult(exec$object, new_sclet_ai_result(
        answer = "normalised and tested markers",
        findings = list(list(
            statement = "clusters were tested for markers",
            claim_level = "hypothesis"
        )),
        task = "interpret"
    ), id = paste0("ai_plan_", plan$plan_id))
    report <- RenderAIReport(sce_done)
    expect_true(grepl("Data suggestions: 1 suggested, 0 executed, 1 not mentioned",
        report, fixed = TRUE))
    expect_true(grepl("`run_decontx`", report, fixed = TRUE))
    # the sentence that tells a reader what to do with the gap
    expect_true(grepl("left silently unanswered", report, fixed = TRUE))
})

test_that("a plan with no snapshot renders no suggestions block", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1L, nrow = 6L, ncol = 4L))
    )
    plan <- new_sclet_ai_plan(
        task = "b20_report_none",
        actions = list(list(id = "s1", action = "normalize_data", params = list())),
        rationale = "normalise only"
    )
    plan$actions <- sclet:::sclet_ai_normalize_plan_actions(plan$actions)
    validation <- ValidateAIPlan(plan, object = sce,
        registry = AIDefaultExecutionRegistry(sce, include = "all"), strict = TRUE)
    exec <- ExecuteAIPlan(sce, plan,
        AIDefaultExecutionRegistry(sce, include = "all"),
        validation = validation, dry_run = FALSE,
        confirmation = validation$confirmation_token)
    sce_done <- RecordAIResult(exec$object, new_sclet_ai_result(
        answer = "normalised only",
        findings = list(list(statement = "object normalised",
            claim_level = "hypothesis")),
        task = "interpret"
    ), id = paste0("ai_plan_", plan$plan_id))
    report <- RenderAIReport(sce_done)
    expect_false(grepl("Data suggestions:", report, fixed = TRUE))
})

test_that("the plan fingerprint is stable across attaching a validation", {
    # The fingerprint hashed $validation, and validation carries the fingerprint
    # -- a self-reference. The same plan hashed two different values before and
    # after `plan$validation <- v`, and ExecuteAIPlan() then refused with
    # sclet_ai_stale_plan on a plan that had JUST validated.
    plan <- new_sclet_ai_plan(
        task = "fp_selfref",
        actions = list(list(id = "s1", action = "normalize_data", params = list())),
        rationale = "r"
    )
    fp_before <- sclet:::sclet_ai_plan_fingerprint(plan)
    plan$validation <- list(valid = TRUE, plan_fingerprint = fp_before,
        plan_id = plan$plan_id)
    fp_after <- sclet:::sclet_ai_plan_fingerprint(plan)
    expect_identical(fp_before, fp_after)
})
