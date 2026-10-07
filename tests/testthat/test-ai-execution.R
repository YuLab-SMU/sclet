test_that("AIDefaultExecutionRegistry has integration group only when explicitly included", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 10, ncol = 8))
    )
    SummarizedExperiment::colData(sce)$batch <- rep(c("group_1", "group_2"), each = 4)

    read_only <- AIDefaultExecutionRegistry(sce, include = "read")
    expect_false("run_integration" %in% names(read_only))

    int_only <- AIDefaultExecutionRegistry(sce, include = "integration")
    expect_true("run_integration" %in% names(int_only))

    all_reg <- AIDefaultExecutionRegistry(sce, include = "all")
    expect_true("run_integration" %in% names(all_reg))
})

test_that("run_integration prerequisites rejects missing batch", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 5, ncol = 4))
    )
    SummarizedExperiment::colData(sce)$batch <- rep(letters[1:2], each = 2)
    sce <- ConfirmAIDesignSemantics(
        sce,
        design = list(batch = "batch")
    )
    reg <- AIDefaultExecutionRegistry(sce, include = "integration")
    action <- reg$run_integration

    result <- action$prerequisites(sce, list(method = "fastMNN"))
    expect_type(result, "character")
    expect_true(grepl("missing required parameter.*batch", result))
})

test_that("run_integration prerequisites rejects unavailable batch column", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 5, ncol = 4))
    )
    SummarizedExperiment::colData(sce)$sample <- rep(letters[1:2], each = 2)
    sce <- ConfirmAIDesignSemantics(
        sce,
        design = list(batch = "sample")
    )
    reg <- AIDefaultExecutionRegistry(sce, include = "integration")
    action <- reg$run_integration

    result <- action$prerequisites(sce, list(batch = "nonexistent_col", method = "fastMNN"))
    expect_type(result, "character")
    expect_true(grepl("batch column is not available", result))
})

test_that("run_integration prerequisites rejects unconfirmed design", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 5, ncol = 4))
    )
    SummarizedExperiment::colData(sce)$batch <- rep(letters[1:2], each = 2)
    reg <- AIDefaultExecutionRegistry(sce, include = "integration")
    action <- reg$run_integration

    result <- action$prerequisites(sce, list(batch = "batch", method = "fastMNN"))
    expect_type(result, "character")
    expect_true(grepl("design_semantics_not_confirmed", result))

    sce_confirmed <- ConfirmAIDesignSemantics(sce, design = list(batch = "batch"))
    result_ok <- action$prerequisites(sce_confirmed, list(batch = "batch", method = "fastMNN"))
    expect_true(isTRUE(result_ok))
})

test_that("run_integration prerequisites returns readable reason for missing harmony", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 5, ncol = 4))
    )
    SummarizedExperiment::colData(sce)$batch <- rep(letters[1:2], each = 2)
    sce <- ConfirmAIDesignSemantics(sce, design = list(batch = "batch"))
    reg <- AIDefaultExecutionRegistry(sce, include = "integration")
    action <- reg$run_integration

    real_req_ns <- base::requireNamespace
    result <- testthat::with_mocked_bindings(
        action$prerequisites(sce, list(batch = "batch", method = "Harmony")),
        .package = "base",
        requireNamespace = function(package, ...) {
            if (identical(package, "harmony")) return(FALSE)
            real_req_ns(package, ...)
        }
    )
    expect_type(result, "character")
    expect_true(grepl("optional_package_missing.*harmony", result))
})

test_that("run_integration action has standard schema fields", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 5, ncol = 4))
    )
    reg <- AIDefaultExecutionRegistry(sce, include = "integration")
    action <- reg$run_integration

    expect_s3_class(action, "sclet_ai_action")
    expect_equal(action$name, "run_integration")
    expect_true(action$requires_confirmation)
    expect_true(action$mutates_object)
    expect_false(action$idempotent)
    expect_true("integration" %in% action$allowed_state_types)
    expect_true("reduction" %in% action$allowed_state_types)
    expect_true(action$estimated_cost %in% c("low", "medium", "high"))

    schema <- action$input_schema
    expect_true(is.list(schema))
    expect_true("batch" %in% names(schema))
    expect_true(is.list(schema$batch) && isTRUE(schema$batch$required))
    expect_true("method" %in% names(schema))
    expect_true(is.character(schema$method) && length(schema$method) > 1L)
    expect_true(all(c("fastMNN", "Harmony", "scVI") %in% schema$method))
    expect_false(".design_confirmed" %in% names(schema))
})

test_that("run_integration rejects plans that only self-declare design confirmation", {
    set.seed(1)
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(rpois(40, 5), nrow = 10, ncol = 4))
    )
    SummarizedExperiment::colData(sce)$batch <- c("a", "a", "b", "b")
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "integration"))
    plan <- sclet:::new_sclet_ai_plan(
        task = "attacker_plan",
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
        actions = list(list(
            id = "integrate", action = "run_integration",
            params = list(batch = "batch", method = "fastMNN")
        ))
    )
    v <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_false(v$valid)
    expect_true(any(grepl("design_semantics_not_confirmed", v$errors)))
})

test_that("run_integration succeeds only after ConfirmAIDesignSemantics is called", {
    set.seed(1)
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(rpois(40, 5), nrow = 10, ncol = 4))
    )
    rownames(sce) <- paste0("g", seq_len(10))
    colnames(sce) <- paste0("c", seq_len(4))
    SummarizedExperiment::colData(sce)$batch <- c("a", "a", "b", "b")
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 5)
    sce <- RunPCA(sce, ncomponents = 2)

    sce <- ConfirmAIDesignSemantics(sce, design = list(batch = "batch"))

    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "integration"))
    plan <- sclet:::new_sclet_ai_plan(
        task = "confirmed_plan",
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
        actions = list(list(
            id = "integrate", action = "run_integration",
            params = list(batch = "batch", method = "fastMNN")
        ))
    )
    v <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_true(v$valid)

    executed <- testthat::with_mocked_bindings({
        ExecuteAIPlan(sce, plan, registry, validation = v, dry_run = FALSE,
            confirmation = v$confirmation_token)
    }, .package = "sclet", RunIntegration = function(object, ...) {
        args <- list(...)
        SingleCellExperiment::reducedDim(object, "correctedPCA") <- matrix(1:8, ncol = 2)
        sclet:::sclet_set_analysis_state(object, "integration", args$name %||% "fastmnn",
            method = "mocked_integration", inputs = list(batch = args$batch), active = FALSE)
    })
    expect_equal(executed$status, "completed")
    expect_false(is.null(SingleCellExperiment::reducedDim(executed$object, "correctedPCA")))
    expect_false(is.null(sclet:::sclet_get_state_record(executed$object, "integration", "fastmnn")))
})

test_that("stale design confirmation (fingerprint mismatch) is rejected", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(rpois(40, 5), nrow = 10, ncol = 4))
    )
    SummarizedExperiment::colData(sce)$batch <- c("a", "a", "b", "b")
    sce <- ConfirmAIDesignSemantics(sce, design = list(batch = "batch"))

    SummarizedExperiment::colData(sce)$plate <- factor(rep(1:2, each = 2))

    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "integration"))
    plan <- sclet:::new_sclet_ai_plan(
        task = "stale_plan",
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
        actions = list(list(
            id = "integrate", action = "run_integration",
            params = list(batch = "batch", method = "fastMNN")
        ))
    )
    v <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_false(v$valid)
    expect_true(any(grepl("design_semantics_not_confirmed", v$errors)))
})

test_that("stale design confirmation (same column name, values reassigned) is rejected", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(rpois(40, 5), nrow = 10, ncol = 4))
    )
    SummarizedExperiment::colData(sce)$batch <- c("a", "a", "b", "b")
    sce <- ConfirmAIDesignSemantics(sce, design = list(batch = "batch"))

    SummarizedExperiment::colData(sce)$batch <- c("T_cell", "T_cell", "B_cell", "B_cell")

    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "integration"))
    plan <- sclet:::new_sclet_ai_plan(
        task = "reassigned_values_plan",
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
        actions = list(list(
            id = "integrate", action = "run_integration",
            params = list(batch = "batch", method = "fastMNN")
        ))
    )
    v <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_false(v$valid)
    expect_true(any(grepl("design_semantics_not_confirmed", v$errors)))
})

test_that("AIDefaultExecutionRegistry isolates annotation group by default", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrix(1, nrow = 10L, ncol = 8L))
    )
    default_reg <- AIDefaultExecutionRegistry(sce)
    expect_false("run_annotation" %in% names(default_reg))
    expect_false("run_de_test" %in% names(default_reg))

    annot_reg <- AIDefaultExecutionRegistry(sce, include = "annotation")
    expect_true("run_annotation" %in% names(annot_reg))
    expect_true("run_de_test" %in% names(annot_reg))

    all_reg <- AIDefaultExecutionRegistry(sce, include = "all")
    expect_true("run_annotation" %in% names(all_reg))
    expect_true("run_de_test" %in% names(all_reg))
})

test_that("run_annotation prerequisites reject missing ref with an explicit warning about defaults", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrix(1, nrow = 10L, ncol = 4L))
    )
    reg <- AIDefaultExecutionRegistry(sce, include = "annotation")
    action <- reg$run_annotation
    reason <- action$prerequisites(sce, list())
    expect_type(reason, "character")
    expect_true(grepl("reference_missing", reason))
    expect_true(grepl("HumanPrimaryCellAtlasData", reason))
})

test_that("run_de_test end-to-end via ValidateAIPlan + ExecuteAIPlan writes active detest state", {
    set.seed(1L)
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(rpois(50 * 40, lambda = 5), nrow = 50L, ncol = 40L,
            dimnames = list(paste0("gene_", seq_len(50L)), paste0("c", seq_len(40L)))))
    )
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 30L)
    sce <- RunPCA(sce, ncomponents = 10L)
    sce <- FindNeighbors(sce, dims = seq_len(5L), reduction = "PCA")
    sce <- FindClusters(sce, resolution = 1.2)

    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "annotation"))
    plan <- sclet:::new_sclet_ai_plan(
        task = "marker_screen",
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
        actions = list(list(id = "de1", action = "run_de_test",
            params = list(name = "pipeline_detest")))
    )
    validation <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_true(validation$valid)
    executed <- ExecuteAIPlan(sce, plan, registry, validation = validation,
        dry_run = FALSE, confirmation = validation$confirmation_token)
    expect_equal(sclet_get_active_state(executed$object, "detest"), "pipeline_detest")
    rec <- sclet_get_state_record(executed$object, "detest", "pipeline_detest")
    expect_false(is.null(rec))
    expect_equal(rec$status, "completed")
})

test_that("AIDefaultExecutionRegistry has rare_cell group only when explicitly included", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 10L, ncol = 8L))
    )
    default_registry <- AIDefaultExecutionRegistry(sce)
    expect_false("run_rare_cell_detection" %in% names(default_registry))
    expect_false("run_doublet_detection" %in% names(default_registry))
    full_registry <- AIDefaultExecutionRegistry(sce, include = c("read", "rare_cell"))
    expect_true("run_rare_cell_detection" %in% names(full_registry))
    expect_true("run_doublet_detection" %in% names(full_registry))
    all_registry <- AIDefaultExecutionRegistry(sce, include = "all")
    expect_true("run_rare_cell_detection" %in% names(all_registry))
})

test_that("run_rare_cell_detection prerequisites reject a missing reduction", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 10L, ncol = 8L))
    )
    reg <- AIDefaultExecutionRegistry(sce, include = "rare_cell")
    reason <- reg$run_rare_cell_detection$prerequisites(sce, list())
    expect_type(reason, "character")
    expect_true(grepl("reduction_missing", reason))
    expect_true(grepl("PCA", reason))
})

test_that("run_doublet_detection prerequisites reject an object without counts", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(values = matrix(1, nrow = 10L, ncol = 8L))
    )
    reg <- AIDefaultExecutionRegistry(sce, include = "rare_cell")
    reason <- reg$run_doublet_detection$prerequisites(sce, list())
    expect_true(grepl("counts_assay_missing", reason))
})

sclet_ai_test_rare_object <- function(n_cells = 60L) {
    set.seed(11L)
    counts <- matrix(rpois(50L * n_cells, lambda = 5L), nrow = 50L, ncol = n_cells,
        dimnames = list(paste0("gene_", seq_len(50L)), paste0("cell_", seq_len(n_cells))))
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 30L)
    RunPCA(sce, ncomponents = 10L)
}

sclet_ai_test_rare_plan <- function(sce, registry, name, rare_threshold = 1000L) {
    plan <- sclet:::new_sclet_ai_plan(
        task = "rare_cell_review",
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
        actions = list(list(id = "rare1", action = "run_rare_cell_detection",
            params = list(name = name, rare_threshold = rare_threshold, k = 5L,
                dims = c(1L, 2L, 3L))))
    )
    validation <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_true(validation$valid)
    ExecuteAIPlan(sce, plan, registry, validation = validation,
        dry_run = FALSE, confirmation = validation$confirmation_token)
}

test_that("run_rare_cell_detection runs end to end and never removes cells", {
    skip_if_not_installed("BiocNeighbors")
    sce <- sclet_ai_test_rare_object()
    counts_before <- unname(SummarizedExperiment::assay(sce, "counts"))
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "rare_cell"))
    executed <- sclet_ai_test_rare_plan(sce, registry, "rareq_e2e")

    expect_equal(ncol(executed$object), ncol(sce))
    expect_equal(unname(SummarizedExperiment::assay(executed$object, "counts")), counts_before)
    expect_false(is.null(SummarizedExperiment::colData(executed$object)$rare_cluster))
    record <- sclet_get_state_record(executed$object, "rare_cells", "rareq_e2e")
    expect_false(is.null(record))
    expect_equal(record$status, "completed")
})

test_that("rare cluster with exactly 1 independent signal is recorded as associated with low_confidence flag", {
    skip_if_not_installed("BiocNeighbors")
    sce <- sclet_ai_test_rare_object()
    # only one signal class is available: doublet calls, with no QC metrics and no sample column
    SummarizedExperiment::colData(sce)$scDblFinder.class <- rep("singlet", ncol(sce))
    SummarizedExperiment::colData(sce)$nCount_RNA <- NULL
    SummarizedExperiment::colData(sce)$nFeature_RNA <- NULL
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "rare_cell"))
    executed <- sclet_ai_test_rare_plan(sce, registry, "rareq_one")

    evidence <- sclet:::sclet_ai_evidence_get_all(executed$object)
    expect_true(length(evidence) >= 1L)
    for (node in evidence) {
        expect_true(grepl("^ev:rare_rareq_one_", node$id))
        expect_equal(node$values$n_independent_signals, 1L)
        expect_equal(node$claim_level, "associated")
        expect_true(node$values$low_confidence)
    }
})

test_that("rare cluster with 2 independent signals is recorded as consistent_with", {
    skip_if_not_installed("BiocNeighbors")
    sce <- sclet_ai_test_rare_object()
    SummarizedExperiment::colData(sce)$scDblFinder.class <- rep(c("singlet", "doublet"), length.out = ncol(sce))
    SummarizedExperiment::colData(sce)$sample_id <- rep(c("sample_a", "sample_b"), length.out = ncol(sce))
    SummarizedExperiment::colData(sce)$nCount_RNA <- NULL
    SummarizedExperiment::colData(sce)$nFeature_RNA <- NULL
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "rare_cell"))
    executed <- sclet_ai_test_rare_plan(sce, registry, "rareq_two")

    evidence <- sclet:::sclet_ai_evidence_get_all(executed$object)
    expect_true(length(evidence) >= 1L)
    for (node in evidence) {
        expect_equal(node$values$n_independent_signals, 2L)
        expect_equal(node$claim_level, "consistent_with")
        expect_false(node$values$low_confidence)
    }
})

test_that("rare cluster with 0 independent signals does not record any evidence", {
    skip_if_not_installed("BiocNeighbors")
    sce <- sclet_ai_test_rare_object()
    SummarizedExperiment::colData(sce)$nCount_RNA <- NULL
    SummarizedExperiment::colData(sce)$nFeature_RNA <- NULL
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "rare_cell"))
    executed <- sclet_ai_test_rare_plan(sce, registry, "rareq_zero")

    ledger <- GetAnalysisLedger(executed$object, detail = "summary",
        include_artifacts = FALSE, include_data = FALSE)
    expect_length(ledger$state_records$ai_evidence, 0L)
    expect_length(sclet:::sclet_ai_evidence_get_all(executed$object), 0L)
    notes <- attr(executed$object, "sclet_ai_note")
    expect_true(length(notes) >= 1L)
    expect_true(all(grepl("no_evidence_recorded", notes)))
})

test_that("run_rare_cell_detection never records evidence from cluster size alone", {
    skip_if_not_installed("BiocNeighbors")
    sce <- sclet_ai_test_rare_object()
    SummarizedExperiment::colData(sce)$nCount_RNA <- NULL
    SummarizedExperiment::colData(sce)$nFeature_RNA <- NULL
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "rare_cell"))
    executed <- sclet_ai_test_rare_plan(sce, registry, "rareq_size_only")

    ledger <- GetAnalysisLedger(executed$object, detail = "summary",
        include_artifacts = FALSE, include_data = FALSE)
    state_records <- unlist(ledger$state_records %||% list(), recursive = FALSE, use.names = FALSE)
    for (record in c(ledger$analyses %||% list(), state_records)) {
        if (!is.list(record) || !identical(as.character(record$type %||% ""), "ai_evidence")) next
        expect_true(as.integer(record$summary$values$n_independent_signals %||% 0L) >= 1L)
        expect_true(record$summary$claim_level %in% c("associated", "consistent_with"))
    }
})

test_that("AIDefaultExecutionRegistry has trajectory group only when explicitly included", {
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1, 2L, 4L)))
    default_registry <- AIDefaultExecutionRegistry(sce)
    expect_false("run_trajectory" %in% names(default_registry))
    full_registry <- AIDefaultExecutionRegistry(sce, include = c("read", "trajectory"))
    expect_true("run_trajectory" %in% names(full_registry))
    all_registry <- AIDefaultExecutionRegistry(sce, include = "all")
    expect_true("run_trajectory" %in% names(all_registry))
})

test_that("the trajectory registry contains no cellrank, fate, spatial, multimodal or regvelo execution actions", {
    ## P5 velocity readiness slice intentionally adds run_velocity (see
    ## .dev/spec-p5-velocity-readiness.md); this guard is narrowed to the
    ## domains still explicitly paused by the roadmap (P5 recommended order:
    ## velocity -> trajectory/velocity interpretation -> CellRank/fate ->
    ## spatial -> multimodal).
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1, 2L, 4L)))
    registry <- AIDefaultExecutionRegistry(sce, include = "all")
    expect_length(names(registry), length(unique(names(registry))))
    expect_true("run_velocity" %in% names(registry))
    expect_false(any(grepl("cellrank|fate|spatial|multimodal|regvelo",
        names(registry), ignore.case = TRUE)))
})

sclet_ai_test_trajectory_object <- function(n_cells = 80L) {
    set.seed(3L)
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
    sce
}

sclet_ai_test_trajectory_plan <- function(sce, registry, params) {
    plan <- sclet:::new_sclet_ai_plan(
        task = "trajectory_inference",
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
        actions = list(list(id = "traj1", action = "run_trajectory", params = params))
    )
    ValidateAIPlan(plan, object = sce, registry = registry, strict = FALSE)
}

test_that("run_trajectory rejects plans without an explicit start_cluster through ValidateAIPlan", {
    skip_if_not_installed("slingshot")
    sce <- sclet_ai_test_trajectory_object()
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "trajectory"))

    omitted <- sclet_ai_test_trajectory_plan(sce, registry, list(group = "cluster"))
    expect_false(omitted$valid)
    expect_true(any(grepl("start_cluster_missing", omitted$errors)))

    nulled <- sclet_ai_test_trajectory_plan(sce, registry,
        list(group = "cluster", start_cluster = NULL))
    expect_false(nulled$valid)
    expect_true(any(grepl("start_cluster_missing", nulled$errors)))

    empty <- sclet_ai_test_trajectory_plan(sce, registry,
        list(group = "cluster", start_cluster = ""))
    expect_false(empty$valid)
    expect_true(any(grepl("start_cluster_missing", empty$errors)))
})

test_that("run_trajectory rejects an unknown start cluster and a missing reduction", {
    sce <- sclet_ai_test_trajectory_object()
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "trajectory"))

    unknown <- sclet_ai_test_trajectory_plan(sce, registry,
        list(group = "cluster", start_cluster = "NOT_A_CLUSTER"))
    expect_false(unknown$valid)
    expect_true(any(grepl("start_cluster_unknown", unknown$errors)))

    bad_group <- sclet_ai_test_trajectory_plan(sce, registry,
        list(group = "not_a_column", start_cluster = "A"))
    expect_false(bad_group$valid)
    expect_true(any(grepl("group_column_missing", bad_group$errors)))

    bad_reduction <- sclet_ai_test_trajectory_plan(sce, registry,
        list(group = "cluster", start_cluster = "A", reduction = "NOT_A_REDUCTION"))
    expect_false(bad_reduction$valid)
    expect_true(any(grepl("reduction_missing", bad_reduction$errors)))
})

test_that("run_trajectory succeeds end-to-end with an explicit start_cluster", {
    skip_if_not_installed("slingshot")
    sce <- sclet_ai_test_trajectory_object()
    idents_before <- as.character(Idents(sce))
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "trajectory"))
    validation <- sclet_ai_test_trajectory_plan(sce, registry,
        list(group = "cluster", start_cluster = "A", name = "traj_ok"))
    expect_true(validation$valid)
    executed <- ExecuteAIPlan(sce, validation$plan, registry, validation = validation,
        dry_run = FALSE, confirmation = validation$confirmation_token)
    expect_equal(executed$status, "completed")
    record <- sclet_get_state_record(executed$object, "trajectory", "traj_ok")
    expect_false(is.null(record))
    expect_equal(record$status, "completed")
    expect_equal(record$params$start_cluster, "A")
    expect_equal(as.character(Idents(executed$object)), idents_before)
})

test_that("run_trajectory records consistent_with evidence without per-cell pseudotime", {
    skip_if_not_installed("slingshot")
    sce <- sclet_ai_test_trajectory_object()
    n_cells <- ncol(sce)
    counts_before <- unname(SummarizedExperiment::assay(sce, "counts"))
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "trajectory"))
    plan <- sclet:::new_sclet_ai_plan(
        task = "trajectory_inference",
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
        actions = list(list(id = "traj1", action = "run_trajectory",
            params = list(group = "cluster", start_cluster = "A", name = "traj_ev")))
    )
    validation <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_true(validation$valid)
    executed <- ExecuteAIPlan(sce, plan, registry, validation = validation,
        dry_run = FALSE, confirmation = validation$confirmation_token)

    evidence <- sclet:::sclet_ai_evidence_get_all(executed$object)
    node <- evidence[["ev:trajectory_traj_ev"]]
    expect_false(is.null(node))
    expect_equal(node$claim_level, "consistent_with")
    expect_true(node$values$relative_ordering_only)
    expect_false(node$values$pseudotime_is_absolute_time)
    # the start cluster is recorded as an anonymized code, never a raw label
    expect_match(node$values$start_group, "^cluster_[0-9]+$")
    expect_false(any(grepl("cluster_label|group_label",
        names(node$values), ignore.case = TRUE)))
    # no per-cell pseudotime vector may appear anywhere in the values
    per_cell <- any(vapply(node$values, function(value) {
        is.atomic(value) && length(value) >= n_cells
    }, logical(1L)))
    expect_false(per_cell)
    # and no time-point style column is invented
    cd_names <- colnames(SummarizedExperiment::colData(executed$object))
    expect_false(any(grepl("^day[0-9]|timepoint|time_point", cd_names, ignore.case = TRUE)))
    expect_equal(ncol(executed$object), n_cells)
    expect_equal(unname(SummarizedExperiment::assay(executed$object, "counts")), counts_before)
})

test_that("run_trajectory is idempotent for identical inputs", {
    skip_if_not_installed("slingshot")
    sce <- sclet_ai_test_trajectory_object()
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "trajectory"))
    run_once <- function(name) {
        plan <- sclet:::new_sclet_ai_plan(
            task = "trajectory_inference",
            context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
            actions = list(list(id = "traj1", action = "run_trajectory",
                params = list(group = "cluster", start_cluster = "A", name = name)))
        )
        validation <- ValidateAIPlan(plan, object = sce, registry = registry)
        ExecuteAIPlan(sce, plan, registry, validation = validation,
            dry_run = FALSE, confirmation = validation$confirmation_token)$object
    }
    first <- run_once("traj_idem_1")
    second <- run_once("traj_idem_2")
    expect_equal(
        SummarizedExperiment::colData(first)$slingPseudotime_1,
        SummarizedExperiment::colData(second)$slingPseudotime_1
    )
})

test_that("run_annotation writes consistent_with evidence without per-cell vectors and never overwrites Idents", {
    set.seed(1L)
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(rpois(50 * 40, lambda = 5), nrow = 50L, ncol = 40L,
            dimnames = list(paste0("gene_", seq_len(50L)), paste0("c", seq_len(40L)))))
    )
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 30L)
    sce <- RunPCA(sce, ncomponents = 10L)
    sce <- FindNeighbors(sce, dims = seq_len(5L), reduction = "PCA")
    sce <- FindClusters(sce, resolution = 1.2)

    ng <- 25L
    use_genes <- paste0("gene_", seq_len(ng))
    ref_mat <- matrix(rpois(ng * 10L, lambda = rep(c(5, 20), each = 5L)),
        nrow = ng, ncol = 10L,
        dimnames = list(use_genes, paste0("ref", seq_len(10L))))
    ref_labels <- rep(c("group_A", "group_B"), each = 5L)
    sce_test <- sce[use_genes, , drop = FALSE]
    before_ident <- ActiveIdent(sce_test)

    registry <- AIDefaultExecutionRegistry(sce_test, include = c("read", "annotation"))
    plan <- sclet:::new_sclet_ai_plan(
        task = "annotate",
        context_fingerprint = GetAnalysisLedger(sce_test)$fingerprint,
        actions = list(list(id = "a1", action = "run_annotation",
            params = list(ref = ref_mat, labels = ref_labels, name = "pipeline_annot")))
    )
    validation <- ValidateAIPlan(plan, object = sce_test, registry = registry)
    expect_true(validation$valid)
    executed <- ExecuteAIPlan(sce_test, plan, registry, validation = validation,
        dry_run = FALSE, confirmation = validation$confirmation_token)
    expect_equal(ActiveIdent(executed$object), before_ident)

    evidence <- sclet:::sclet_ai_evidence_get_all(executed$object)
    expect_true(length(evidence) >= 1L)
    annotation_ev <- evidence[["ev:annotation_pipeline_annot"]]
    expect_false(is.null(annotation_ev))
    expect_true(annotation_ev$claim_level %in% c("associated", "consistent_with"))
    expect_false(annotation_ev$claim_level %in% c("measured", "observed"))

    n_cells <- ncol(executed$object)
    per_cell_found <- any(vapply(annotation_ev$values, function(v) {
        is.atomic(v) && length(v) >= n_cells
    }, logical(1L)))
    expect_false(isTRUE(per_cell_found))
})
