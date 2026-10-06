test_that("deterministic advanced diagnostics are bounded and handle missing inputs", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrix(1:24, nrow = 6L, ncol = 4L)),
        colData = S4Vectors::DataFrame(
            sample_id = c("secret-A", "secret-A", "secret-B", "secret-B"),
            condition = c("control", "control", "treated", "treated"),
            cluster = c("private-c1", "private-c1", "private-c2", "private-c2"),
            qc_metric = c(1, 2, 10, 12)
        )
    )
    SingleCellExperiment::reducedDim(sce, "PCA") <- cbind(PC1 = c(1, 2, 3, 4), PC2 = c(4, 3, 2, 1))

    qc <- sclet:::summarize_qc_by_group(sce, "sample_id")
    assoc <- sclet:::summarize_pca_metadata_association(sce, "condition")
    composition <- sclet:::summarize_cluster_sample_composition(sce, "sample_id")
    small <- sclet:::summarize_small_clusters(sce, threshold = 2L)
    ready <- sclet:::check_integration_readiness(sce)

    expect_equal(qc$status, "available")
    expect_equal(qc$metrics$qc_metric$by_group$group_1$mean, 1.5)
    expect_equal(assoc$status, "available")
    expect_equal(assoc$reduction, "PCA")
    expect_equal(composition$status, "available")
    expect_equal(dim(composition$counts), c(2L, 2L))
    expect_equal(small$n_small_clusters, 2L)
    expect_equal(ready$status, "clarification_required")
    expect_false(grepl("secret-A|private-c1", paste(capture.output(str(list(qc, assoc, composition, small))), collapse = " ")))
})

test_that("check_rare_cell_readiness requires cluster assignment and PCA", {
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1, 2L, 4L)))
    result <- sclet:::check_rare_cell_readiness(sce)
    expect_equal(result$status, "not_ready")
    expect_true("run_rare_cell_detection" %in% result$blocked_actions)
    expect_true(length(result$questions) >= 1L)
    expect_false(result$checks$has_pca)
})

test_that("check_rare_cell_readiness is ready with cluster and PCA and reports missing doublet evidence", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrix(1, 6L, 8L)),
        colData = S4Vectors::DataFrame(cluster = rep(c("c1", "c2"), each = 4L))
    )
    SingleCellExperiment::reducedDim(sce, "PCA") <- matrix(as.numeric(seq_len(16L)), nrow = 8L, ncol = 2L)
    result <- sclet:::check_rare_cell_readiness(sce)
    expect_equal(result$status, "ready_for_diagnostic")
    expect_true(result$checks$has_cluster_column)
    expect_false(result$checks$doublet_evidence_available)
    expect_true(grepl("doublet evidence missing", paste(result$notes, collapse = " ")))
})

test_that("summarize_small_cluster_evidence counts only available signals", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrix(1, 6L, 8L)),
        colData = S4Vectors::DataFrame(
            cluster = c("c1", "c1", rep("c2", 6L)),
            sample_id = rep(c("s1", "s2"), each = 4L),
            scDblFinder.class = c("singlet", "doublet", "singlet", "singlet",
                "doublet", "singlet", "singlet", "doublet")
        )
    )
    result <- sclet:::summarize_small_cluster_evidence(sce, "cluster", size_threshold = 2L)
    expect_equal(result$status, "available")
    expect_equal(result$n_small_clusters, 1L)
    expect_equal(result$size, 2L)
    expect_true(isTRUE(result$independent_signals$doublet$available))
    expect_true(isTRUE(result$independent_signals$sample_replication$available))
    expect_false(isTRUE(result$independent_signals$qc$available))
    expect_equal(result$independent_signals$qc$reason, "no_numeric_qc_columns")
    expect_false(isTRUE(result$independent_signals$marker$available))
    expect_equal(result$n_independent_signals_available, 2L)
    expect_false(result$raw_values_included)
    expect_false(grepl("c1|c2|s1|s2", paste(capture.output(str(result)), collapse = " ")))
})

test_that("summarize_small_cluster_evidence reports qc deviation without judging the population", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrix(1, 6L, 6L)),
        colData = S4Vectors::DataFrame(
            cluster = c("small", "small", "other", "other", "other", "other"),
            nCount_RNA = c(2, 3, 90, 95, 100, 105)
        )
    )
    result <- sclet:::summarize_small_cluster_evidence(sce, "cluster", size_threshold = 2L)
    expect_true(isTRUE(result$independent_signals$qc$available))
    expect_true(result$independent_signals$qc$deviates_from_other_cells)
    expect_equal(result$n_independent_signals_available, 1L)
    # the diagnostic must not contain any verdict about the population itself
    expect_false(any(c("real", "noise", "verdict", "is_rare") %in% names(result)))
})

sclet_ai_diag_rare_fixture <- function() {
    set.seed(1)
    counts <- matrix(rpois(50 * 40, lambda = 5), nrow = 50L, ncol = 40L,
        dimnames = list(paste0("gene_", seq_len(50L)), paste0("c", seq_len(40L))))
    counts[1:5, 1:6] <- counts[1:5, 1:6] + 60
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 30L)
    RunPCA(sce, ncomponents = 10L)
}

sclet_ai_diag_marker_credited <- function(summary) {
    vapply(summary$clusters, function(entry) isTRUE(entry$independent_signals$marker$available), logical(1L))
}

test_that("an all-cluster marker run credits populations individually, not object-wide", {
    skip_if_not_installed("BiocNeighbors")
    # regression: a FindAllMarkers run stores no per-cluster labels in its state
    # record, so it used to credit every population on the object at once
    sce <- sclet_ai_diag_rare_fixture()
    sce <- FindNeighbors(sce, dims = seq_len(5L), reduction = "PCA")
    sce <- FindClusters(sce, resolution = 1.2)
    sce <- RunDEtest(sce, name = "all_markers")

    rare <- RunRareCellDetection(sce, name = "rareX", k = 5, dims = 1:3, rare_threshold = 1000)
    summary <- sclet:::summarize_small_cluster_evidence(rare, cluster = "rare_cluster", size_threshold = 1000)
    credited <- sclet_ai_diag_marker_credited(summary)

    expect_equal(summary$status, "available")
    expect_gt(summary$n_small_clusters, 2L)
    expect_lt(sum(credited), length(credited))
    expect_gt(sum(credited), 0L)
    # the object really does hold one marker record; credit is still per population
    expect_equal(summary$clusters[[1L]]$independent_signals$marker$n_marker_records, 1L)
    for (index in which(credited)) {
        expect_true(summary$clusters[[index]]$independent_signals$marker$cluster_covered_by_marker_analysis)
        expect_gte(summary$clusters[[index]]$independent_signals$marker$n_significant_up_genes, 1L)
    }
    for (index in which(!credited)) {
        expect_equal(summary$clusters[[index]]$independent_signals$marker$reason,
            "no_marker_evidence_attributable_to_population")
    }
})

test_that("marker credit needs significant up genes attributable to the population", {
    skip_if_not_installed("BiocNeighbors")
    sce <- sclet_ai_diag_rare_fixture()
    sce <- RunRareCellDetection(sce, name = "rareB", k = 5, dims = 1:3, rare_threshold = 1000)
    ActiveIdent(sce) <- "rare_cluster"
    sce <- RunDEtest(sce, name = "rare_markers")

    # noise data with no marker passing the default adjusted-p-value cutoff
    strict <- sclet:::summarize_small_cluster_evidence(sce, cluster = "rare_cluster", size_threshold = 1000)
    expect_identical(sum(sclet_ai_diag_marker_credited(strict)), 0L)
    expect_identical(sum(vapply(strict$clusters, function(entry) entry$n_independent_signals_available,
        integer(1L))), 0L)

    # relaxing the cutoff to a level no real claim would use shows that exactly the
    # populations whose own group passes are the ones credited
    loose <- sclet:::summarize_small_cluster_evidence(sce, cluster = "rare_cluster",
        size_threshold = 1000, pvalue_cutoff = 0.5)
    credited <- sclet_ai_diag_marker_credited(loose)
    expect_identical(sum(credited), 1L)
    expect_gte(loose$clusters[[which(credited)[1L]]]$independent_signals$marker$n_significant_up_genes, 1L)

    expect_error(sclet:::summarize_small_cluster_evidence(sce, cluster = "rare_cluster", pvalue_cutoff = 5))
    expect_error(sclet:::summarize_small_cluster_evidence(sce, cluster = "rare_cluster", min_up_genes = 0L))
})

test_that("uninformative signals are not counted as independent evidence", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrix(1, 6L, 8L)),
        colData = S4Vectors::DataFrame(
            cluster = c("c1", "c1", rep("c2", 6L)),
            const_value = rep(7, 8L),
            all_missing = rep(NA_real_, 8L),
            scDblFinder.class = rep(NA_character_, 8L),
            sample_id = rep("sample_a", 8L)
        )
    )
    result <- sclet:::summarize_small_cluster_evidence(sce, "cluster", size_threshold = 2L)
    signals <- result$independent_signals
    expect_equal(result$n_independent_signals_available, 0L)
    expect_false(isTRUE(signals$qc$available))
    expect_equal(signals$qc$reason, "no_informative_qc_columns")
    expect_equal(signals$qc$n_qc_columns_screened, 2L)
    expect_false(isTRUE(signals$doublet$available))
    expect_equal(signals$doublet$reason, "scDblFinder_calls_missing_for_population_or_rest")
    expect_false(isTRUE(signals$sample_replication$available))
    expect_equal(signals$sample_replication$reason, "fewer_than_two_samples")
    expect_false(isTRUE(signals$marker$available))
    expect_equal(signals$marker$reason, "no_marker_evidence_attributable_to_population")
    expect_false(any(c("real", "noise", "verdict", "is_rare") %in% names(result)))
})

test_that("cells without a sample label are not counted as an extra sample group", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrix(1, 6L, 8L)),
        colData = S4Vectors::DataFrame(
            cluster = c("c1", "c1", rep("c2", 6L)),
            sample_id = c("s1", NA, "s2", "s2", "s1", "s2", "s1", "s2")
        )
    )
    result <- sclet:::summarize_small_cluster_evidence(sce, "cluster", size_threshold = 2L)
    replication <- result$independent_signals$sample_replication
    expect_true(replication$available)
    expect_equal(replication$n_samples_total, 2L)
    expect_equal(replication$n_samples_present, 1L)
    expect_equal(replication$n_cells_unlabeled, 1L)
    expect_false(isTRUE(replication$present_in_multiple_samples))
})

test_that("summarize_small_cluster_evidence returns not_available without cluster information", {
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1, 4L, 4L)))
    result <- sclet:::summarize_small_cluster_evidence(sce)
    expect_equal(result$status, "not_available")
    expect_equal(result$reason, "missing_cluster_column")
})

    test_that("check_trajectory_readiness requires cluster assignment and a reduction", {
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1, 2L, 4L)))
    result <- sclet:::check_trajectory_readiness(sce)
    expect_equal(result$status, "not_ready")
    expect_true("trajectory" %in% result$blocked_actions)
    expect_true(length(result$questions) >= 1L)
})

test_that("check_trajectory_readiness never recommends a specific start cluster", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrix(1, 6L, 8L)),
        colData = S4Vectors::DataFrame(cluster = rep(c("c1", "c2"), each = 4L))
    )
    SingleCellExperiment::reducedDim(sce, "UMAP") <- matrix(as.numeric(seq_len(16L)), nrow = 8L, ncol = 2L)
    result <- sclet:::check_trajectory_readiness(sce)
    expect_equal(result$status, "ready_for_diagnostic")
    expect_equal(result$checks$reduction_resolved, "UMAP")
    expect_equal(result$checks$n_clusters, 2L)
    flattened <- unlist(result, recursive = TRUE, use.names = TRUE)
    expect_false(any(grepl("recommend|suggest|candidate_root|start_cluster|root_|origin_cluster",
        names(flattened), ignore.case = TRUE)))
    expect_true(any(grepl("no root or start cluster is suggested", tolower(unlist(result)), fixed = FALSE)))
})

test_that("summarize_trajectory_cluster_order reports position only and never suggests a root", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrix(1, 6L, 6L)),
        colData = S4Vectors::DataFrame(cluster = c("c1", "c1", "c2", "c2", "c2", "c2"))
    )
    SingleCellExperiment::reducedDim(sce, "UMAP") <- cbind(
        UMAP_1 = c(1, 2, 8, 9, 10, 11),
        UMAP_2 = c(3, 4, 12, 13, 14, 15)
    )
    result <- sclet:::summarize_trajectory_cluster_order(sce, "cluster")
    expect_equal(result$status, "available")
    expect_equal(result$reduction, "UMAP")
    expect_equal(result$n_groups, 2L)
    expect_false(result$root_suggested)
    expect_equal(result$by_dimension$dim_1$group_1$mean, 1.5)
    expect_equal(result$by_dimension$dim_1$group_2$mean, 9.5)
    # no per-cell values and no ordering claim
    dumped <- paste(capture.output(str(result)), collapse = " ")
    expect_false(grepl("c1|c2", dumped))
})

test_that("summarize_trajectory_cluster_order returns not_available without a usable reduction", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrix(1, 6L, 4L)),
        colData = S4Vectors::DataFrame(cluster = rep(c("c1", "c2"), each = 2L))
    )
    result <- sclet:::summarize_trajectory_cluster_order(sce, "cluster")
    expect_equal(result$status, "not_available")
    expect_equal(result$reason, "reduction_not_available")
})

test_that("diagnostics return typed not_available when required inputs are absent", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrix(1, nrow = 2L, ncol = 2L))
    )
    expect_equal(sclet:::summarize_qc_by_group(sce, "missing")$status, "not_available")
    expect_equal(sclet:::summarize_pca_metadata_association(sce, "missing")$status, "not_available")
    expect_equal(sclet:::summarize_cluster_sample_composition(sce, "missing")$status, "not_available")
    expect_equal(sclet:::summarize_small_clusters(sce)$status, "not_available")
    expect_equal(sclet:::check_integration_readiness(sce)$status, "clarification_required")
})

test_that("integration readiness requires confirmed batch semantics", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrix(1, nrow = 2L, ncol = 4L)),
        colData = S4Vectors::DataFrame(batch = rep(c("x", "y"), each = 2), condition = rep(c("a", "b"), 2))
    )
    expect_equal(sclet:::check_integration_readiness(sce)$status, "clarification_required")
    expect_equal(sclet:::check_integration_readiness(sce, list(batch = "batch", condition = "condition"))$status, "ready_for_diagnostic")
    expect_equal(sclet:::check_integration_readiness(sce, list(batch = "condition", condition = "condition"))$status, "not_ready")
})

test_that("annotation readiness reports clarification_required or not_ready before clustering", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrix(rpois(20 * 6, lambda = 5), nrow = 20L, ncol = 6L))
    )
    SummarizedExperiment::assay(sce, "logcounts") <- log2(SummarizedExperiment::assay(sce, "counts") + 1)
    ready <- sclet:::check_annotation_readiness(sce)
    expect_true(ready$status %in% c("not_ready", "clarification_required"))
    expect_true("run_annotation" %in% ready$blocked_actions)
    expect_true("run_de_test" %in% ready$blocked_actions)
    expect_true(length(ready$questions) >= 1L)
})
