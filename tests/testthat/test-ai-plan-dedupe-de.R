## Redundant DE steps are removed, but never one something depends on.

dedupe_plan <- function(...) {
    sclet:::new_sclet_ai_plan(task = "dedupe", actions = list(...))
}

markers_step <- function(id = "m1") {
    list(id = id, action = "run_de_test", params = list(name = "findallmarkers"))
}
contrast_step <- function(id = "c1") {
    list(id = id, action = "run_de_test",
        params = list(group = "condition", ident.1 = "stim", ident.2 = "ctrl",
            covariates = "donor", name = "cond"))
}

test_that("all-cluster markers are dropped when a targeted contrast exists", {
    plan <- dedupe_plan(markers_step(), contrast_step())
    out <- sclet:::sclet_ai_deduplicate_de_steps(plan,
        "Compare stimulated versus control cells within donors.")
    expect_equal(length(out$actions), 1L)
    expect_equal(out$actions[[1L]]$action, "run_de_test")
    expect_equal(out$actions[[1L]]$id, "c1")
    expect_true(any(grepl("m1: run_de_test", out$metadata$removed_redundant_steps)))
})

test_that("markers stay when the goal asks for them", {
    plan <- dedupe_plan(markers_step(), contrast_step())
    out <- sclet:::sclet_ai_deduplicate_de_steps(plan,
        "Find marker genes for each cluster and compare stim vs control.")
    expect_equal(length(out$actions), 2L)
})

test_that("markers stay when a later step references them", {
    plan <- dedupe_plan(markers_step(),
        list(id = "e1", action = "run_enrichment",
            params = list(detest_id = "${m1.output.detest_id}")),
        contrast_step())
    out <- sclet:::sclet_ai_deduplicate_de_steps(plan,
        "Compare stimulated versus control cells within donors.")
    expect_equal(length(out$actions), 3L)
})

test_that("a plan with only markers is untouched", {
    plan <- dedupe_plan(markers_step())
    out <- sclet:::sclet_ai_deduplicate_de_steps(plan, "Find clusters.")
    expect_equal(length(out$actions), 1L)
    expect_null(out$metadata$removed_redundant_steps)
})

test_that("markers stay when a success criterion references them", {
    plan <- dedupe_plan(markers_step(), contrast_step())
    plan$success_criteria <- list(list(id = "c9",
        criterion = "marker table from ${m1.output.detest_id} is non-empty"))
    out <- sclet:::sclet_ai_deduplicate_de_steps(plan,
        "Compare stimulated versus control cells within donors.")
    expect_equal(length(out$actions), 2L)
    expect_null(out$metadata$removed_redundant_steps)
})

test_that("markers stay when a criterion names the step by id", {
    ## Kang reliability, 2026-10-19: the criterion carried `step_id = "m1"`,
    ## not a token, so the dedupe removed m1 and validation then failed with
    ## `success_criterion_unresolvable: c9 references unknown step m1`.
    plan <- dedupe_plan(markers_step(), contrast_step())
    plan$success_criteria <- list(list(id = "c9", source = "action_output",
        step_id = "m1", check = "exists", field = "n_genes"))
    out <- sclet:::sclet_ai_deduplicate_de_steps(plan,
        "Compare stimulated versus control cells within donors.")
    expect_equal(length(out$actions), 2L)
    expect_null(out$metadata$removed_redundant_steps)
})

test_that("markers stay when a later step depends on them", {
    ## Same run: the later step carried `depends_on = "m1"`, which lives on the
    ## action, not in its params, and the scan never looked at it.
    plan <- dedupe_plan(markers_step(),
        list(id = "e1", action = "run_enrichment",
            params = list(detest_id = "markers"), depends_on = "m1"),
        contrast_step())
    out <- sclet:::sclet_ai_deduplicate_de_steps(plan,
        "Compare stimulated versus control cells within donors.")
    expect_equal(length(out$actions), 3L)
    expect_null(out$metadata$removed_redundant_steps)
})
