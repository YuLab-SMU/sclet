## End-to-end RunAIAnalysis() coverage for existing domain adapters (P4)
##
## Proves that the four existing domain adapters (marker/DE, annotation,
## rare-cell/doublet, trajectory) complete a real context -> plan -> execute ->
## report loop through the RunAIAnalysis() orchestrator, not just the
## ValidateAIPlan()/ExecuteAIPlan() layer.

sclet_ai_test_e2e_clustered_object <- function(n_cells = 60L) {
    set.seed(1L)
    counts <- matrix(rpois(50L * n_cells, lambda = 5L), nrow = 50L, ncol = n_cells,
        dimnames = list(paste0("gene_", seq_len(50L)), paste0("c", seq_len(n_cells))))
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 30L)
    sce <- RunPCA(sce, ncomponents = 10L)
    sce <- FindNeighbors(sce, dims = seq_len(5L), reduction = "PCA")
    sce <- FindClusters(sce, resolution = 1.2)
    sce
}

test_that("run_de_test survives RunAIAnalysis end to end with bounded evidence", {
    sce <- sclet_ai_test_e2e_clustered_object()
    old <- options(sclet.ai.call = function(task, context, ...) {
        list(
            answer = "Run an all-cluster marker test.",
            proposed_actions = list(list(
                id = "de1",
                action = "run_de_test",
                params = list(name = "e2e_markers")
            ))
        )
    })
    on.exit(options(old), add = TRUE)

    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "annotation"))
    result <- RunAIAnalysis(sce, goal = "Find markers for every cluster.",
        confirm = "yes", registry = registry)

    expect_equal(result$status, "completed")
    expect_equal(result$execution$status, "completed")

    evidence <- sclet:::sclet_ai_evidence_get_all(result$object)
    de_ev <- evidence[["ev:findallmarkers"]]
    expect_false(is.null(de_ev))
    expect_equal(de_ev$claim_level, "associated")

    n_cells <- ncol(result$object)
    per_cell <- any(vapply(de_ev$values, function(v) {
        is.atomic(v) && length(v) >= n_cells
    }, logical(1L)))
    expect_false(isTRUE(per_cell))
})

test_that("run_annotation survives RunAIAnalysis end to end and never overwrites Idents", {
    sce <- sclet_ai_test_e2e_clustered_object()
    ng <- 25L
    use_genes <- paste0("gene_", seq_len(ng))
    set.seed(2L)
    ref_mat <- matrix(rpois(ng * 10L, lambda = rep(c(5, 20), each = 5L)),
        nrow = ng, ncol = 10L,
        dimnames = list(use_genes, paste0("ref", seq_len(10L))))
    ref_labels <- rep(c("group_A", "group_B"), each = 5L)
    sce_sub <- sce[use_genes, , drop = FALSE]
    before_ident <- ActiveIdent(sce_sub)

    old <- options(sclet.ai.call = function(task, context, ...) {
        list(
            answer = "Annotate clusters against the supplied reference.",
            proposed_actions = list(list(
                id = "ann1",
                action = "run_annotation",
                params = list(ref = ref_mat, labels = ref_labels, name = "e2e_singler")
            ))
        )
    })
    on.exit(options(old), add = TRUE)

    registry <- AIDefaultExecutionRegistry(sce_sub, include = c("read", "annotation"))
    result <- RunAIAnalysis(sce_sub, goal = "Annotate clusters.",
        confirm = "yes", registry = registry)

    expect_equal(result$status, "completed")
    expect_equal(result$execution$status, "completed")
    expect_equal(ActiveIdent(result$object), before_ident)

    evidence <- sclet:::sclet_ai_evidence_get_all(result$object)
    ann_ev <- evidence[["ev:annotation_e2e_singler"]]
    expect_false(is.null(ann_ev))
    expect_true(ann_ev$claim_level %in% c("associated", "consistent_with"))
})

test_that("run_doublet_detection + run_rare_cell_detection survive RunAIAnalysis and never remove cells", {
    skip_if_not_installed("BiocNeighbors")
    skip_if_not_installed("scDblFinder")

    set.seed(11L)
    n_cells <- 60L
    counts <- matrix(rpois(50L * n_cells, lambda = 5L), nrow = 50L, ncol = n_cells,
        dimnames = list(paste0("gene_", seq_len(50L)), paste0("cell_", seq_len(n_cells))))
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 30L)
    sce <- RunPCA(sce, ncomponents = 10L)

    old <- options(sclet.ai.call = function(task, context, ...) {
        list(
            answer = "Score doublets then look for rare populations.",
            proposed_actions = list(
                list(id = "dbl1", action = "run_doublet_detection", params = list()),
                list(id = "rare1", action = "run_rare_cell_detection",
                    params = list(name = "e2e_rare", reduction = "PCA",
                        rare_threshold = 1000L, k = 5L, dims = c(1L, 2L, 3L)))
            )
        )
    })
    on.exit(options(old), add = TRUE)

    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "rare_cell"))
    # scDblFinder emits its own informational notice about small sample sizes on
    # a fixture this small (60 cells); this is benign upstream diagnostic output,
    # not a sclet behavior under test, so it is suppressed rather than left to
    # leak as an unasserted testthat warning.
    result <- suppressWarnings(RunAIAnalysis(sce, goal = "Find doublets and rare populations.",
        confirm = "yes", registry = registry))

    expect_equal(result$status, "completed")
    expect_equal(result$execution$status, "completed")
    expect_equal(ncol(result$object), ncol(sce))
    expect_false(is.null(SummarizedExperiment::colData(result$object)$scDblFinder.class))
})

test_that("run_trajectory survives RunAIAnalysis end to end with bounded evidence", {
    skip_if_not_installed("slingshot")

    set.seed(3L)
    n_cells <- 80L
    counts <- matrix(rpois(60L * n_cells, lambda = 8L), nrow = 60L, ncol = n_cells,
        dimnames = list(paste0("gene_", seq_len(60L)), paste0("cell_", seq_len(n_cells))))
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 30L)
    sce <- RunPCA(sce, ncomponents = 10L)
    set.seed(4L)
    embedding <- matrix(stats::rnorm(n_cells * 2L), ncol = 2L)
    embedding[1:30, 1] <- embedding[1:30, 1] - 4
    embedding[31:60, 2] <- embedding[31:60, 2] - 4
    embedding[61:80, 1] <- embedding[61:80, 1] + 3
    embedding[61:80, 2] <- embedding[61:80, 2] + 3
    SingleCellExperiment::reducedDim(sce, "UMAP") <- embedding
    SummarizedExperiment::colData(sce)$cluster <- rep(c("A", "B", "C"), times = c(30L, 30L, 20L))
    ActiveIdent(sce) <- "cluster"

    real_start <- as.character(Idents(sce))[1L]

    old <- options(sclet.ai.call = function(task, context, ...) {
        list(
            answer = "Infer a slingshot trajectory from the supplied start cluster.",
            proposed_actions = list(list(
                id = "traj1",
                action = "run_trajectory",
                params = list(group = "cluster", start_cluster = real_start,
                    reduction = "UMAP", name = "e2e_traj")
            ))
        )
    })
    on.exit(options(old), add = TRUE)

    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "trajectory"))
    result <- RunAIAnalysis(sce, goal = "Infer a trajectory.",
        confirm = "yes", registry = registry)

    expect_equal(result$status, "completed")
    expect_equal(result$execution$status, "completed")

    evidence <- sclet:::sclet_ai_evidence_get_all(result$object)
    traj_ev <- evidence[["ev:trajectory_e2e_traj"]]
    expect_false(is.null(traj_ev))
    expect_equal(traj_ev$claim_level, "consistent_with")
    expect_false(traj_ev$values$pseudotime_is_absolute_time)
    expect_true(traj_ev$values$relative_ordering_only)
})

test_that("a prerequisite failure surfaces as invalid_plan with structured clarification through RunAIAnalysis", {
    skip_if_not_installed("slingshot")

    set.seed(3L)
    n_cells <- 80L
    counts <- matrix(rpois(60L * n_cells, lambda = 8L), nrow = 60L, ncol = n_cells,
        dimnames = list(paste0("gene_", seq_len(60L)), paste0("cell_", seq_len(n_cells))))
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 30L)
    sce <- RunPCA(sce, ncomponents = 10L)
    set.seed(4L)
    embedding <- matrix(stats::rnorm(n_cells * 2L), ncol = 2L)
    SingleCellExperiment::reducedDim(sce, "UMAP") <- embedding
    SummarizedExperiment::colData(sce)$cluster <- rep(c("A", "B", "C"), times = c(30L, 30L, 20L))
    ActiveIdent(sce) <- "cluster"

    old <- options(sclet.ai.call = function(task, context, ...) {
        list(
            answer = "Infer a trajectory without specifying a root.",
            proposed_actions = list(list(
                id = "traj1",
                action = "run_trajectory",
                params = list(group = "cluster", reduction = "UMAP", name = "no_root")
            ))
        )
    })
    on.exit(options(old), add = TRUE)

    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "trajectory"))
    result <- RunAIAnalysis(sce, goal = "Infer a trajectory.",
        confirm = "yes", registry = registry)

    expect_equal(result$status, "invalid_plan")
    expect_null(result$execution)
    expect_equal(result$report$clarification$status, "clarification_required")
    question_ids <- vapply(result$report$clarification$questions,
        function(q) q$id, character(1L))
    expect_true("trajectory_root" %in% question_ids)
})

test_that("RunAIAnalysis with interpret = TRUE works for a domain other than inspect_status/run_integration", {
    sce <- sclet_ai_test_e2e_clustered_object()

    old <- options(sclet.ai.call = function(task, context, ...) {
        if (identical(task, "analysis_plan")) {
            list(
                answer = "Run an all-cluster marker test.",
                proposed_actions = list(list(
                    id = "de1",
                    action = "run_de_test",
                    params = list(name = "interp_markers")
                ))
            )
        } else if (identical(task, "analysis_explanation")) {
            list(
                answer = "Interpretation of the marker test.",
                findings = list(list(
                    statement = "Markers were computed for each cluster.",
                    severity = "info",
                    claim_level = "associated",
                    evidence_refs = "ev:findallmarkers"
                )),
                evidence = list("ev:findallmarkers"),
                warnings = list(),
                recommendations = list(),
                proposed_actions = list()
            )
        } else {
            list(answer = "Unhandled task.")
        }
    })
    on.exit(options(old), add = TRUE)

    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "annotation"))
    result <- RunAIAnalysis(sce, goal = "Find markers and interpret them.",
        confirm = "yes", registry = registry, interpret = TRUE)

    expect_equal(result$status, "completed")
    expect_false(is.null(result$report$interpretation))
    expect_s3_class(result$report$interpretation, "sclet_ai_result")
    expect_true(validate_sclet_ai_result(result$report$interpretation, error = FALSE))
    expect_null(result$report$interpretation_error)

    interp_ids <- grep("^ai_interpretation_",
        names(GetAnalysisLedger(result$object)$analyses), value = TRUE)
    expect_true(length(interp_ids) >= 1L)
})
