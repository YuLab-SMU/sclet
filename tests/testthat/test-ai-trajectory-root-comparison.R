## Trajectory multi-root comparison (Phase D / trajectory follow-up)
##
## Proves that compare_trajectory_roots() reads only completed trajectory state
## records plus bounded ev:trajectory_<name> evidence nodes, compares explicit
## root hypotheses as a read-only summary layer, and never chooses a preferred
## root or exposes per-cell values.

sclet_ai_test_root_object <- function() {
    counts <- matrix(1, nrow = 12L, ncol = 18L,
        dimnames = list(paste0("gene_", seq_len(12L)), paste0("cell_", seq_len(18L))))
    SingleCellExperiment::SingleCellExperiment(list(counts = counts))
}

sclet_ai_test_add_root_run <- function(
    sce,
    id,
    root_code = "cluster_1",
    start_group_code = 1L,
    n_lineages = 2L,
    n_groups = 3L,
    median = 0.5,
    iqr = 0.25,
    max_val = 1.0,
    reduction = "UMAP",
    group = "cluster"
) {
    sce <- sclet_set_analysis_state(
        object = sce,
        type = "trajectory",
        id = id,
        method = "slingshot",
        inputs = list(reduction = reduction, group = group),
        artifacts = list(analysis_key = "trajectory"),
        params = list(start_cluster = paste0("root_for_", id)),
        summary = list(n_lineages = n_lineages),
        active = FALSE
    )
    sce <- RecordAIEvidence(sce, list(
        id = paste0("ev:trajectory_", id),
        kind = "deterministic_summary",
        values = list(
            n_lineages = as.integer(n_lineages),
            n_groups = as.integer(n_groups),
            start_group_code = as.integer(start_group_code),
            group_codes = paste0("cluster_", seq_len(n_groups)),
            start_group = root_code,
            n_summarized_lineages = 1L,
            first_lineage_median = as.numeric(median),
            first_lineage_iqr = as.numeric(iqr),
            first_lineage_max = as.numeric(max_val),
            relative_ordering_only = TRUE,
            pseudotime_is_absolute_time = FALSE,
            raw_values_included = FALSE
        ),
        claim_level = "consistent_with"
    ), source = id)
    sce
}

test_that("compare_trajectory_roots returns not_available when no trajectory runs exist", {
    sce <- sclet_ai_test_root_object()
    result <- compare_trajectory_roots(sce)
    expect_equal(result$status, "not_available")
    expect_equal(result$reason, "no_trajectory_runs")
    expect_equal(result$n_runs, 0L)
    expect_null(result$recommendation)
    expect_false(result$raw_values_included)
})

test_that("compare_trajectory_roots returns not_available when only one usable run exists", {
    sce <- sclet_ai_test_root_object()
    sce <- sclet_ai_test_add_root_run(sce, "traj_one", root_code = "cluster_1")
    result <- compare_trajectory_roots(sce)
    expect_equal(result$status, "not_available")
    expect_equal(result$reason, "at_least_two_completed_trajectory_runs_required")
    expect_equal(result$n_runs, 1L)
    expect_false(result$raw_values_included)
})

test_that("compare_trajectory_roots compares two different root hypotheses", {
    sce <- sclet_ai_test_root_object()
    sce <- sclet_ai_test_add_root_run(sce, "traj_a",
        root_code = "cluster_1", start_group_code = 1L,
        n_lineages = 2L, iqr = 0.2, max_val = 1.0)
    sce <- sclet_ai_test_add_root_run(sce, "traj_b",
        root_code = "cluster_2", start_group_code = 2L,
        n_lineages = 3L, iqr = 0.3, max_val = 1.2)
    result <- compare_trajectory_roots(sce)

    expect_equal(result$status, "available")
    expect_equal(result$n_runs, 2L)
    expect_equal(result$n_roots, 2L)
    expect_equal(result$n_comparisons, 1L)
    expect_equal(result$claim_level, "consistent_with")
    expect_null(result$recommendation)
    expect_false(result$raw_values_included)
    expect_true(all(c("traj_a", "traj_b") %in% names(result$runs)))
    expect_equal(result$runs$traj_a$root_code, "cluster_1")
    expect_equal(result$runs$traj_b$root_code, "cluster_2")
    expect_equal(result$roots$cluster_1$run_count, 1L)
    expect_equal(result$roots$cluster_2$run_count, 1L)
    expect_false(result$comparisons$comparison_1$same_root)
    expect_equal(result$comparisons$comparison_1$root_difference,
        "different_recorded_root_code")
})

test_that("compare_trajectory_roots keeps same-root runs available without choosing a winner", {
    sce <- sclet_ai_test_root_object()
    sce <- sclet_ai_test_add_root_run(sce, "traj_same_1",
        root_code = "cluster_1", start_group_code = 1L,
        n_lineages = 2L, iqr = 0.25, max_val = 1.0)
    sce <- sclet_ai_test_add_root_run(sce, "traj_same_2",
        root_code = "cluster_1", start_group_code = 1L,
        n_lineages = 2L, iqr = 0.5, max_val = 1.0)
    result <- compare_trajectory_roots(sce)

    expect_equal(result$status, "available")
    expect_equal(result$n_roots, 1L)
    expect_true(result$comparisons$comparison_1$same_root)
    expect_equal(result$comparisons$comparison_1$root_difference,
        "same_recorded_root_code")
    expect_null(result$recommendation)
})

test_that("compare_trajectory_roots respects explicit ids filtering", {
    sce <- sclet_ai_test_root_object()
    sce <- sclet_ai_test_add_root_run(sce, "traj_keep_1", root_code = "cluster_1")
    sce <- sclet_ai_test_add_root_run(sce, "traj_keep_2", root_code = "cluster_2")
    sce <- sclet_ai_test_add_root_run(sce, "traj_drop", root_code = "cluster_3")
    result <- compare_trajectory_roots(sce, ids = c("traj_keep_1", "traj_keep_2"))

    expect_equal(result$status, "available")
    expect_equal(result$n_runs, 2L)
    expect_true(all(c("traj_keep_1", "traj_keep_2") %in% names(result$runs)))
    expect_false("traj_drop" %in% names(result$runs))
    expect_equal(result$n_roots, 2L)
})

test_that("compare_trajectory_roots skips trajectory runs without bounded evidence", {
    sce <- sclet_ai_test_root_object()
    sce <- sclet_ai_test_add_root_run(sce, "traj_with_ev", root_code = "cluster_1")
    sce <- sclet_set_analysis_state(
        object = sce,
        type = "trajectory",
        id = "traj_without_ev",
        method = "slingshot",
        inputs = list(reduction = "UMAP", group = "cluster"),
        artifacts = list(analysis_key = "trajectory"),
        summary = list(n_lineages = 2L),
        active = FALSE
    )
    result <- compare_trajectory_roots(sce)

    expect_equal(result$status, "not_available")
    expect_equal(result$reason, "at_least_two_completed_trajectory_runs_required")
    expect_equal(result$n_runs, 1L)
})

test_that("compare_trajectory_roots returns not_available for unknown requested ids", {
    sce <- sclet_ai_test_root_object()
    sce <- sclet_ai_test_add_root_run(sce, "traj_known", root_code = "cluster_1")
    result <- compare_trajectory_roots(sce, ids = c("traj_missing"))
    expect_equal(result$status, "not_available")
    expect_equal(result$reason, "no_requested_trajectory_runs")
    expect_equal(result$n_runs, 0L)
})

test_that("compare_trajectory_roots output is bounded and privacy-safe", {
    sce <- sclet_ai_test_root_object()
    sce <- sclet_ai_test_add_root_run(sce, "traj_priv_1", root_code = "cluster_1")
    sce <- sclet_ai_test_add_root_run(sce, "traj_priv_2", root_code = "cluster_2")
    result <- compare_trajectory_roots(sce)

    expect_equal(result$status, "available")
    all_names <- unique(c(
        names(result),
        unlist(lapply(result$runs, names), use.names = FALSE),
        unlist(lapply(result$roots, names), use.names = FALSE),
        unlist(lapply(result$comparisons, names), use.names = FALSE)
    ))
    bad <- vapply(all_names, sclet:::sclet_ai_evidence_bad_field, logical(1L))
    expect_false(any(bad),
        info = paste("Bad field names found:", paste(all_names[bad], collapse = ", ")))
    expect_true(is.list(result$caveats))
    expect_lte(length(result$caveats), 5L)
    expect_false(isTRUE(result$matching$cell_level_overlap_available))
    expect_true(isTRUE(result$matching$root_labels_anonymized))
    expect_null(result$recommendation)
})
