test_that("CompareAIAnalyses reports no routes without recorded integration analyses", {
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = matrix(1, 2L, 2L)))
    result <- sclet:::CompareAIAnalyses(sce)
    expect_equal(result$status, "not_available")
    expect_false(result$execution$performed)
})

test_that("CompareAIAnalyses compares recorded routes without executing", {
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = matrix(1, 2L, 4L)))
    sce <- sclet:::sclet_set_analysis_state(
        sce, type = "integration", id = "raw_baseline", method = "raw",
        summary = list(status = "completed", metrics = list(
            list(name = "biological_preservation", value = 1, direction = "higher_is_better",
                uncertainty = list(status = "reported")),
            list(name = "runtime_sec", value = 1, direction = "lower_is_better",
                uncertainty = list(status = "reported"))
        ))
    )
    sce <- sclet:::sclet_set_analysis_state(
        sce, type = "integration", id = "harmony_1", method = "harmony",
        summary = list(status = "completed", metrics = list(
            list(name = "batch_mixing", value = 0.8, direction = "higher_is_better"),
            list(name = "biological_preservation", value = 0.9, direction = "higher_is_better"),
            list(name = "runtime_sec", value = 2, direction = "lower_is_better")
        ))
    )
    result <- sclet:::CompareAIAnalyses(sce, criterion = "biological_preservation")
    expect_equal(result$status, "available")
    expect_equal(result$baseline, "raw_baseline")
    expect_length(result$routes, 2L)
    expect_true(length(result$metrics) >= 3L)
    expect_null(result$recommendation)
    expect_false(result$execution$performed)
    expect_false(grepl("counts|1 1 1", paste(utils::capture.output(str(result)), collapse = " ")))
})

test_that("CompareAIAnalyses orders routes by criterion", {
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = matrix(1, 2L, 4L)))
    sce <- sclet:::sclet_set_analysis_state(sce, type = "integration", id = "raw_baseline",
        method = "raw",
        summary = list(status = "completed", metrics = list(
            list(name = "biological_preservation", value = 0.5, direction = "higher_is_better"),
            list(name = "batch_mixing", value = 0.3, direction = "higher_is_better")
        )))
    sce <- sclet:::sclet_set_analysis_state(sce, type = "integration", id = "route_A",
        method = "mocked_A",
        summary = list(status = "completed", metrics = list(
            list(name = "biological_preservation", value = 0.95, direction = "higher_is_better"),
            list(name = "batch_mixing", value = 0.2, direction = "higher_is_better")
        )))
    sce <- sclet:::sclet_set_analysis_state(sce, type = "integration", id = "route_B",
        method = "mocked_B",
        summary = list(status = "completed", metrics = list(
            list(name = "biological_preservation", value = 0.7, direction = "higher_is_better"),
            list(name = "batch_mixing", value = 0.4, direction = "higher_is_better")
        )))
    r <- CompareAIAnalyses(sce, criterion = "biological_preservation")
    ord <- names(r$routes)
    # raw_baseline can be first or last; the others should have route_A before route_B
    non_baseline <- setdiff(ord, "raw_baseline")
    expect_equal(which(non_baseline == "route_A"), 1L)
})

test_that("CompareAIAnalyses detects batch_vs_bio tradeoff and refuses recommendation", {
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = matrix(1, 2L, 4L)))
    sce <- sclet:::sclet_set_analysis_state(sce, type = "integration", id = "raw_baseline",
        method = "raw",
        summary = list(status = "completed", metrics = list(
            list(name = "batch_mixing", value = 0.4, direction = "higher_is_better"),
            list(name = "biological_preservation", value = 0.8, direction = "higher_is_better")
        )))
    sce <- sclet:::sclet_set_analysis_state(sce, type = "integration", id = "harmony_conflict",
        method = "harmony",
        summary = list(status = "completed", metrics = list(
            list(name = "batch_mixing", value = 0.85, direction = "higher_is_better"),
            list(name = "biological_preservation", value = 0.6, direction = "higher_is_better")
        )))
    r <- CompareAIAnalyses(sce, criterion = "biological_preservation")
    tradeoffs <- r$tradeoffs
    expect_true(is.list(tradeoffs) && (length(tradeoffs) > 0L || !is.null(tradeoffs$status)))
    has_conflict <- any(vapply(tradeoffs, function(t) identical(t$type, "batch_vs_bio"), logical(1L)))
    expect_true(has_conflict)
    expect_null(r$recommendation)
    expect_false("best_route" %in% names(r))
})

test_that("CompareAIAnalyses normalizes missing metrics to not_available", {
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = matrix(1, 2L, 4L)))
    sce <- sclet:::sclet_set_analysis_state(sce, type = "integration", id = "raw_baseline",
        method = "raw",
        summary = list(status = "completed", metrics = list(
            list(name = "batch_mixing", value = 0.3, direction = "higher_is_better")
        )))
    r <- CompareAIAnalyses(sce, criterion = "biological_preservation")
    baseline_metrics_key <- grep("^raw_baseline::", names(r$metrics), value = TRUE)
    metric_names <- sub("^raw_baseline::", "", baseline_metrics_key)
    expect_true("cluster_stability" %in% metric_names)
    stab_key <- "raw_baseline::cluster_stability"
    m <- r$metrics[[stab_key]]
    expect_identical(m$uncertainty$status, "not_available")
    expect_true(is.na(m$value))
})

sclet_ai_test_rare_comparison_state <- function(object, id, summary_status = "completed") {
    sclet:::sclet_set_analysis_state(
        object, type = "rare_cells", id = id, method = "density",
        inputs = list(reduction = "PCA", dims = 1:3),
        summary = list(status = summary_status, n_rare_clusters = 2L)
    )
}

sclet_ai_test_rare_comparison_node <- function(object, source, label, size,
                                                signals, claim_level = "associated",
                                                low_confidence = TRUE) {
    RecordAIEvidence(
        object,
        list(
            id = paste0("ev:rare_", source, "_", label),
            kind = "deterministic_summary",
            values = list(
                population_label = label,
                population_size = as.integer(size),
                population_fraction = size / 100,
                n_independent_signals = as.integer(signals),
                low_confidence = isTRUE(low_confidence),
                raw_values_included = FALSE
            ),
            claim_level = claim_level
        ),
        source = source
    )
}

test_that("compare_rare_cell_evidence handles missing and single-run evidence", {
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1, 2L, 4L)))
    empty <- compare_rare_cell_evidence(sce)
    expect_equal(empty$status, "not_available")
    expect_equal(empty$reason, "no_rare_cell_evidence_records")

    sce <- sclet_ai_test_rare_comparison_state(sce, "rare_a")
    sce <- sclet_ai_test_rare_comparison_node(sce, "rare_a", "cluster_1", 3L, 1L)
    single <- compare_rare_cell_evidence(sce)
    expect_equal(single$status, "not_available")
    expect_equal(single$reason, "at_least_two_completed_rare_cell_runs_required")
    expect_equal(single$n_runs, 1L)
})

test_that("compare_rare_cell_evidence excludes non-completed runs and invalid labels", {
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1, 2L, 4L)))
    sce <- sclet_ai_test_rare_comparison_state(sce, "rare_failed", summary_status = "failed")
    sce <- sclet_ai_test_rare_comparison_node(sce, "rare_failed", "cluster_1", 3L, 1L)
    result <- compare_rare_cell_evidence(sce)
    expect_equal(result$status, "not_available")
    expect_equal(result$reason, "no_rare_cell_evidence_records")
    invalid_node <- list(
        kind = "deterministic_summary",
        source = "rare_failed",
        values = list(
            population_label = NA_character_,
            population_size = 3L,
            n_independent_signals = 1L
        )
    )
    expect_false(sclet:::sclet_ai_rare_cell_node_is_valid(invalid_node, "rare_failed"))
})

test_that("compare_rare_cell_evidence compares runs without inflating same-run support", {
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1, 2L, 4L)))
    sce <- sclet_ai_test_rare_comparison_state(sce, "rare_a")
    sce <- sclet_ai_test_rare_comparison_state(sce, "rare_b")
    sce <- sclet_ai_test_rare_comparison_node(sce, "rare_a", "cluster_1", 3L, 1L)
    sce <- sclet_ai_test_rare_comparison_node(sce, "rare_a", "cluster_2", 2L, 2L,
        claim_level = "consistent_with", low_confidence = FALSE)
    sce <- sclet_ai_test_rare_comparison_node(sce, "rare_b", "cluster_1", 4L, 2L,
        claim_level = "consistent_with", low_confidence = FALSE)
    sce <- sclet_ai_test_rare_comparison_node(sce, "rare_b", "cluster_3", 1L, 1L)

    result <- compare_rare_cell_evidence(sce)
    expect_equal(result$status, "available")
    expect_equal(result$n_runs, 2L)
    expect_equal(result$n_populations, 3L)
    expect_equal(result$n_recurring_populations, 1L)
    expect_equal(result$matching$method, "same_recorded_population_label")
    expect_false(result$matching$cell_level_overlap_available)
    expect_equal(result$independence$n_support_lines, 2L)
    expect_true(result$independence$populations_within_one_run_are_not_independent)
    expect_equal(result$populations$cluster_1$run_count, 2L)
    expect_true(result$populations$cluster_1$recurring)
    expect_setequal(result$populations$cluster_1$runs_present, c("run_1", "run_2"))
    expect_equal(result$runs$run_1$n_populations, 2L)
    expect_equal(result$runs$run_2$n_populations, 2L)
    expect_false(result$raw_values_included)
    expect_null(result$recommendation)
    expect_false(grepl("rare_a|rare_b", paste(utils::capture.output(str(result)), collapse = " ")))
})

test_that("compare_rare_cell_evidence can restrict completed runs", {
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1, 2L, 4L)))
    sce <- sclet_ai_test_rare_comparison_state(sce, "rare_a")
    sce <- sclet_ai_test_rare_comparison_state(sce, "rare_b")
    sce <- sclet_ai_test_rare_comparison_node(sce, "rare_a", "cluster_1", 3L, 1L)
    sce <- sclet_ai_test_rare_comparison_node(sce, "rare_b", "cluster_1", 4L, 2L,
        claim_level = "consistent_with", low_confidence = FALSE)
    result <- compare_rare_cell_evidence(sce, ids = c("rare_b", "missing"))
    expect_equal(result$status, "not_available")
    expect_equal(result$n_runs, 1L)
    expect_error(compare_rare_cell_evidence(sce, ids = 1), "ids must be")
})
