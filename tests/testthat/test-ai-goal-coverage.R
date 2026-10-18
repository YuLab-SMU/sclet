## Goal coverage: a plan can be valid and still not answer the question.
##
## Measured on the first PBMC 10k real-data pilot (2026-10-18): the goal said
## "annotate the cell types", one run's plan omitted run_annotation entirely,
## validation passed, execution completed, and nothing in the pipeline or the
## report said so. These tests pin the deterministic detection, the single
## feedback round that tries to repair it, and the report field.

goal_coverage_plan <- function(goal, actions) {
    sclet:::new_sclet_ai_plan(
        task = "goal coverage test",
        actions = actions,
        metadata = list(goal = goal)
    )
}

test_that("coverage names the actions a goal asks for and what is missing", {
    coverage <- sclet:::sclet_ai_goal_coverage
    plan <- goal_coverage_plan(
        "Annotate the cell types and find marker genes for each cluster.",
        list(list(id = "s1", action = "run_annotation",
            params = list(ref = "x", labels = "y")))
    )
    result <- coverage(plan)
    expect_true("run_annotation" %in% result$required)
    expect_true("run_de_test" %in% result$required)
    expect_false("run_annotation" %in% result$missing)
    expect_true("run_de_test" %in% result$missing)

    ## A goal without analysis keywords reports nothing (no false positives).
    quiet <- coverage(goal_coverage_plan(
        "Normalize the counts, run PCA, build a graph and cluster the cells.",
        list(list(id = "s1", action = "normalize_data", params = list()))
    ))
    expect_length(quiet$required, 0L)
    expect_length(quiet$missing, 0L)

    ## No goal: nothing to check.
    none <- coverage(goal_coverage_plan(NULL,
        list(list(id = "s1", action = "normalize_data", params = list()))))
    expect_null(none$goal)
    expect_length(none$missing, 0L)
})

test_that("a goal gap triggers one feedback round and a complete valid plan replaces it", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(stats::rpois(40L * 20L, 5), nrow = 40L, ncol = 20L,
            dimnames = list(paste0("gene_", seq_len(40L)), paste0("c", seq_len(20L)))))
    )
    ref <- matrix(stats::rpois(20L * 10L, 10L), nrow = 20L, ncol = 10L,
        dimnames = list(paste0("gene_", seq_len(20L)), paste0("r", seq_len(10L))))
    labels <- rep(c("A", "B"), each = 5L)

    calls <- 0L
    old <- options(sclet.ai.call = function(task, context, ...) {
        calls <<- calls + 1L
        if (calls == 1L) {
            list(answer = "normalize only",
                proposed_actions = list(list(id = "s1", action = "normalize_data",
                    params = list())))
        } else {
            list(answer = "normalize and annotate",
                proposed_actions = list(
                    list(id = "s1", action = "normalize_data", params = list()),
                    list(id = "s2", action = "run_annotation",
                        params = list(ref = ref, labels = labels,
                            name = "annot_goal"))
                ))
        }
    })
    on.exit(options(old), add = TRUE)

    plan <- AIPlanAnalysis(sce, goal = "Annotate the cell types in this dataset.",
        model = "mock", registry = AIDefaultExecutionRegistry(sce, include = "all"))
    expect_equal(calls, 2L)
    expect_true("run_annotation" %in% vapply(plan$actions, function(a) a$action, character(1)))
    expect_equal(plan$metadata$goal_repair_rounds_used, 1L)
    expect_length(plan$metadata$goal_coverage$missing, 0L)
})

test_that("an invalid or no-better candidate is not accepted over the executable plan", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(stats::rpois(30L * 12L, 5), nrow = 30L, ncol = 12L,
            dimnames = list(paste0("gene_", seq_len(30L)), paste0("c", seq_len(12L)))))
    )
    calls <- 0L
    old <- options(sclet.ai.call = function(task, context, ...) {
        calls <<- calls + 1L
        if (calls == 1L) {
            list(answer = "normalize only",
                proposed_actions = list(list(id = "s1", action = "normalize_data",
                    params = list())))
        } else {
            ## Invalid candidate: an unregistered action name.
            list(answer = "bogus",
                proposed_actions = list(
                    list(id = "s1", action = "normalize_data", params = list()),
                    list(id = "s2", action = "not_registered", params = list())
                ))
        }
    })
    on.exit(options(old), add = TRUE)

    plan <- AIPlanAnalysis(sce, goal = "Annotate the cell types in this dataset.",
        model = "mock", registry = AIDefaultExecutionRegistry(sce, include = "all"))
    expect_equal(calls, 2L)
    expect_equal(length(plan$actions), 1L)
    expect_true(isTRUE(plan$validation$valid))
    expect_equal(plan$metadata$goal_repair_rounds_used, 0L)
    expect_true("run_annotation" %in% plan$metadata$goal_coverage$missing)
})

test_that("RunAIAnalysis surfaces the goal gap in its report", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(stats::rpois(30L * 12L, 5), nrow = 30L, ncol = 12L,
            dimnames = list(paste0("gene_", seq_len(30L)), paste0("c", seq_len(12L)))))
    )
    old <- options(sclet.ai.call = function(task, context, ...) {
        list(answer = "normalize only",
            proposed_actions = list(list(id = "s1", action = "normalize_data",
                params = list())))
    })
    on.exit(options(old), add = TRUE)

    result <- RunAIAnalysis(sce, goal = "Annotate the cell types in this dataset.",
        model = "mock", confirm = TRUE, include = "all",
        max_repair_rounds = 0L)
    expect_equal(result$status, "completed")
    expect_true("goal_coverage" %in% names(result$report))
    expect_true("run_annotation" %in% result$report$goal_coverage$missing)
})
