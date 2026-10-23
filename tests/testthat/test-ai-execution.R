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
    expect_false("dims" %in% names(schema))
    expect_false(".design_confirmed" %in% names(schema))
})

test_that("cached validation cannot authorize a modified plan", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 5L, ncol = 4L))
    )
    registry <- AIExecutionRegistry(list(
        mark_a = AIAction(
            name = "mark_a",
            handler = function(object, params) {
                SummarizedExperiment::colData(object)$marker <- rep("A", ncol(object))
                object
            },
            returns = "sce",
            mutates_object = TRUE,
            allowed_state_types = character()
        ),
        mark_b = AIAction(
            name = "mark_b",
            handler = function(object, params) {
                SummarizedExperiment::colData(object)$marker <- rep("B", ncol(object))
                object
            },
            returns = "sce",
            mutates_object = TRUE,
            allowed_state_types = character()
        )
    ))
    plan <- sclet:::new_sclet_ai_plan(
        task = "tamper",
        actions = list(list(id = "step1", action = "mark_a", params = list()))
    )
    validation <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_true(validation$valid)
    expect_match(validation$plan_fingerprint, "^sclet-ai-plan-")

    tampered <- plan
    tampered$actions[[1L]]$action <- "mark_b"
    expect_error(
        ExecuteAIPlan(sce, tampered, registry, validation = validation,
            dry_run = FALSE, confirmation = validation$confirmation_token,
            record = FALSE),
        class = "sclet_ai_stale_plan"
    )
    expect_null(SummarizedExperiment::colData(sce)$marker)
})

test_that("a valid raw (not pre-normalized) plan is not misrejected as stale", {
    # Regression: ValidateAIPlan fingerprinted the caller's raw plan while
    # ExecuteAIPlan fingerprints the plan AFTER normalizing plan$actions
    # (default-filling max_retries, requires_confirmation, params, ...).
    # The two digests differed for any not-yet-normalized plan, so every
    # valid plan failed with sclet_ai_stale_plan on a non-dry run.
    set.seed(1)
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(rpois(40, 5), nrow = 10, ncol = 4))
    )
    SummarizedExperiment::colData(sce)$celltype <- rep(c("a", "b"), each = 2)
    registry <- AIDefaultExecutionRegistry(sce, include = "read")
    action <- names(registry)[[1L]]
    raw_plan <- list(
        plan_id = "plan_raw_shape",
        task = "raw",
        # Deliberately NOT run through sclet_ai_normalize_plan_actions: no
        # requires_confirmation / max_retries / depends_on defaults.
        actions = list(list(id = "step1", action = action, params = list()))
    )
    validation <- ValidateAIPlan(raw_plan, object = sce, registry = registry)
    expect_true(validation$valid)
    executed <- ExecuteAIPlan(
        sce, raw_plan, registry, validation = validation,
        dry_run = FALSE, confirmation = validation$confirmation_token,
        record = FALSE
    )
    expect_s3_class(executed, "sclet_ai_execution")
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
    # A real subpopulation (genes 1-8 elevated in 6 cells): the corrected
    # seeding rule finds density PEAKS, so a pure-noise blob collapses into
    # one group with no rare populations and the claim-graded evidence tests
    # would have nothing to grade. The injection gives the algorithm a peak
    # to find, which is what the graded-evidence contract is about.
    counts[1:8, 1:6] <- counts[1:8, 1:6] + 40L
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

test_that("the full registry exposes no spatial, multimodal or regvelo computation action", {
    ## Each slice intentionally adds its own computation action as it lands
    ## (run_velocity, run_fate_analysis). The domains still paused by the
    ## roadmap are spatial, multimodal and regvelo -- but the guard must test
    ## that no *computation* action exists for them, not that the substring is
    ## absent. The spatial/multimodal batches deliberately shipped read-only
    ## catalog actions (inspect_*_readiness, summarize_*_evidence,
    ## compare_*), so a substring check now fails on exactly the actions that
    ## were designed to be safe. Mutation is what this guard is about.
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1, 2L, 4L)))
    registry <- AIDefaultExecutionRegistry(sce, include = "all")
    expect_length(names(registry), length(unique(names(registry))))
    expect_true("run_velocity" %in% names(registry))
    expect_true("run_fate_analysis" %in% names(registry))

    paused_domains <- grepl("spatial|multimodal|regvelo", names(registry), ignore.case = TRUE)
    mutation_actions <- names(registry)[paused_domains & grepl("^run_", names(registry))]
    # These domains were paused while their actions were unreviewed. Two have
    # since been reviewed and ship as driver-only mutation actions; anything else
    # appearing here is a regression of that rule.
    expect_setequal(
        mutation_actions,
        c("run_spatial_deconvolution", "run_regvelo")
    )

    ## The read-only catalog actions for those domains are expected, and none
    ## of them may mutate the object -- that is what makes them safe to ship
    ## ahead of the computation actions.
    catalog <- setdiff(names(registry)[paused_domains], mutation_actions)
    expect_true(length(catalog) > 0L)
    expect_false(any(vapply(registry[catalog], function(action) {
        isTRUE(action$mutates_object) || isTRUE(action$requires_confirmation)
    }, logical(1L))))
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

# ---- Regressions for the 2026-10-17 pbmc3k annotation benchmark -------------

# A tiny reference with a labels column, so the string-column path is exercised
# without any celldex download.
annotation_test_reference <- function(genes, n_ref = 10L, seed = 2L) {
    set.seed(seed)
    ref <- SingleCellExperiment::SingleCellExperiment(
        list(logcounts = matrix(
            stats::rpois(length(genes) * n_ref, lambda = rep(c(5, 20), each = n_ref / 2L)),
            nrow = length(genes), ncol = n_ref,
            dimnames = list(genes, paste0("ref", seq_len(n_ref)))
        ))
    )
    SummarizedExperiment::colData(ref)$label.main <- rep(c("group_A", "group_B"), each = n_ref / 2L)
    ref
}

annotation_test_object <- function(n_genes = 30L, n_cells = 40L, seed = 1L) {
    set.seed(seed)
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(stats::rpois(n_genes * n_cells, lambda = 5),
            nrow = n_genes, ncol = n_cells,
            dimnames = list(paste0("gene_", seq_len(n_genes)), paste0("c", seq_len(n_cells)))))
    )
    NormalizeData(sce)
}

test_that("RunSingleR resolves a labels string as a reference colData column", {
    sce <- annotation_test_object()
    ref <- annotation_test_reference(rownames(sce))

    ## The documented spelling. Before the fix this reached SingleR as a
    ## one-element label vector and died on
    ## `length(curlabels) == ncol(curref)`, which is how the 2026-10-17
    ## pbmc3k + MonacoImmuneData run failed after the whole chain had run.
    annotated <- RunSingleR(sce, ref = ref, labels = "label.main")
    labels <- SummarizedExperiment::colData(annotated)[["SingleR_labels"]]
    expect_equal(length(labels), ncol(sce))
    expect_true(all(labels %in% c("group_A", "group_B")))

    ## A string that is neither a column nor a one-sample label is refused here
    ## with the available names, instead of surfacing SingleR's assertion.
    err <- tryCatch(RunSingleR(sce, ref = ref, labels = "not_a_column"),
        error = function(e) e)
    expect_s3_class(err, "error")
    expect_true(grepl("label.main", conditionMessage(err), fixed = TRUE))
})

test_that("run_annotation executes from a labels string and records its evidence", {
    sce <- annotation_test_object()
    ref <- annotation_test_reference(rownames(sce))
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "annotation"))
    plan <- sclet:::new_sclet_ai_plan(
        task = "annotate",
        actions = list(list(id = "a1", action = "run_annotation",
            params = list(ref = ref, labels = "label.main", name = "annot_str")))
    )
    validation <- ValidateAIPlan(plan, object = sce, registry = registry, strict = TRUE)
    expect_true(validation$valid)
    executed <- ExecuteAIPlan(sce, validation$plan, registry, validation = validation,
        dry_run = FALSE, confirmation = validation$confirmation_token)
    expect_equal(executed$status, "completed")

    labels <- SummarizedExperiment::colData(executed$object)[["annot_str_labels"]]
    expect_equal(length(labels), ncol(sce))

    evidence <- sclet:::sclet_ai_evidence_get_all(executed$object)
    node <- evidence[["ev:annotation_annot_str"]]
    expect_false(is.null(node))
    expect_equal(node$claim_level, "consistent_with")
    ## Counts and fractions only: nothing with one entry per cell.
    expect_false(any(vapply(node$values, function(v) {
        is.atomic(v) && length(v) >= ncol(executed$object)
    }, logical(1L))))
})

test_that("run_annotation refuses to report success without its evidence node", {
    sce <- annotation_test_object()
    ref <- annotation_test_reference(rownames(sce))
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "annotation"))
    plan <- sclet:::new_sclet_ai_plan(
        task = "annotate",
        actions = list(list(id = "a1", action = "run_annotation",
            params = list(ref = ref, labels = "label.main", name = "annot_no_ev")))
    )
    validation <- ValidateAIPlan(plan, object = sce, registry = registry, strict = TRUE)
    expect_true(validation$valid)

    executed <- with_mocked_bindings(
        ExecuteAIPlan(sce, validation$plan, registry, validation = validation,
            dry_run = FALSE, confirmation = validation$confirmation_token),
        RecordAIEvidence = function(...) stop("simulated ledger refusal"),
        .package = "sclet"
    )
    expect_equal(executed$status, "failed")
    errors <- vapply(executed$results, function(step) step$error %||% "", character(1L))
    expect_true(any(grepl("simulated ledger refusal", errors, fixed = TRUE)))
})

test_that("the annotation evidence id is sanitized like every other record id", {
    ## A record id is caller-supplied; the evidence id must not echo it raw.
    expect_equal(sclet:::sclet_ai_annotation_evidence_id("annot.1 run"), "ev:annotation_annot_1_run")
    expect_equal(sclet:::sclet_ai_annotation_evidence_id("plain_id"), "ev:annotation_plain_id")
})

test_that("a celldex reference spelling is normalized but a placeholder is not", {
    ## Measured 2026-10-17: the model wrote the namespaced literal
    ## `<celldex::MonacoImmuneData>`; the handler's `grepl("Data$")` rejected it
    ## after the whole chain had executed, and execution never re-plans.
    normalize <- sclet:::sclet_celldex_reference_name
    expect_equal(normalize("MonacoImmuneData"), "MonacoImmuneData")
    expect_equal(normalize("<celldex::MonacoImmuneData>"), "MonacoImmuneData")
    expect_equal(normalize("celldex::HumanPrimaryCellAtlasData"), "HumanPrimaryCellAtlasData")
    expect_equal(normalize("\"celldex::MonacoImmuneData\""), "MonacoImmuneData")
    expect_equal(normalize("`celldex::MonacoImmuneData`"), "MonacoImmuneData")

    ## A placeholder, a typo and an NA are not references: the user's choice
    ## must be the only thing that can arrive.
    expect_null(normalize("HELD - user-supplied celldex reference required"))
    expect_null(normalize("MonacoImmune"))
    expect_null(normalize("MonacoImmuneData()"))
    expect_null(normalize(NA_character_))
    expect_null(normalize(c("MonacoImmuneData", "HumanPrimaryCellAtlasData")))
})

test_that("a namespaced celldex reference executes instead of failing after the chain", {
    sce <- annotation_test_object()
    ref <- annotation_test_reference(rownames(sce))
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "annotation"))
    plan <- sclet:::new_sclet_ai_plan(
        task = "annotate",
        actions = list(list(id = "a1", action = "run_annotation",
            params = list(ref = "<celldex::MonacoImmuneData>", labels = "label.main",
                name = "annot_namespaced")))
    )
    validation <- ValidateAIPlan(plan, object = sce, registry = registry, strict = TRUE)
    expect_true(validation$valid)

    ## The dataset lookup is the only thing mocked; the reference itself is a
    ## real SummarizedExperiment, so no download happens.
    executed <- with_mocked_bindings(
        ExecuteAIPlan(sce, validation$plan, registry, validation = validation,
            dry_run = FALSE, confirmation = validation$confirmation_token),
        sclet_celldex_dataset = function(name) function() ref,
        .package = "sclet"
    )
    expect_equal(executed$status, "completed")
    expect_true("ev:annotation_annot_namespaced" %in%
        names(sclet:::sclet_ai_evidence_get_all(executed$object)))
})

test_that("a placeholder reference is refused as MISSING at validation", {
    ## The choice belongs to the user, so the failure has to say "missing":
    ## `reference_invalid` has no clarification question, and on the PBMC 10k
    ## real-data run (2026-10-18) a placeholder plan was discarded without ever
    ## producing an `annotation_reference` question.
    sce <- annotation_test_object()
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "annotation"))
    plan <- sclet:::new_sclet_ai_plan(
        task = "annotate",
        actions = list(list(id = "a1", action = "run_annotation",
            params = list(ref = "HELD - user-supplied celldex reference required",
                labels = "label.main")))
    )
    err <- tryCatch(ValidateAIPlan(plan, object = sce, registry = registry, strict = TRUE),
        error = function(e) e)
    expect_s3_class(err, "sclet_ai_invalid_plan")
    expect_true(grepl("reference_missing", conditionMessage(err), fixed = TRUE))
    expect_true(grepl("Ask the user", conditionMessage(err), fixed = TRUE))
})

test_that("an angle-bracketed user-named reference is a placeholder, not a dataset", {
    sce <- annotation_test_object()
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "annotation"))
    plan <- sclet:::new_sclet_ai_plan(
        task = "annotate",
        actions = list(list(id = "a1", action = "run_annotation",
            params = list(ref = "<user-named celldex dataset>",
                labels = "<user-named labels column>")))
    )
    err <- tryCatch(ValidateAIPlan(plan, object = sce, registry = registry, strict = TRUE),
        error = function(e) e)
    expect_s3_class(err, "sclet_ai_invalid_plan")
    expect_true(grepl("reference_missing", conditionMessage(err), fixed = TRUE))
})

test_that("a placeholder labels column is refused as MISSING at validation", {
    sce <- annotation_test_object()
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "annotation"))
    plan <- sclet:::new_sclet_ai_plan(
        task = "annotate",
        actions = list(list(id = "a1", action = "run_annotation",
            params = list(ref = "celldex::MonacoImmuneData",
                labels = "<user-named labels column>")))
    )
    err <- tryCatch(ValidateAIPlan(plan, object = sce, registry = registry, strict = TRUE),
        error = function(e) e)
    expect_s3_class(err, "sclet_ai_invalid_plan")
    expect_true(grepl("labels_missing", conditionMessage(err), fixed = TRUE))
})

test_that("a near-miss celldex name is refused at validation, not in the handler", {
    ## Shape-based parsing cannot know which names celldex exports:
    ## `BlueprintENCODEData` normalizes fine (ends with "Data") and the real
    ## export is `BlueprintEncodeData`. Before the export check, that typo
    ## passed validation and failed inside the handler -- after the upstream
    ## steps had run.
    skip_if_not_installed("SingleR")
    skip_if_not_installed("celldex")
    expect_equal(
        sclet:::sclet_celldex_reference_name("BlueprintENCODEData"),
        "BlueprintENCODEData"
    )
    sce <- annotation_test_object()
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "annotation"))
    plan <- sclet:::new_sclet_ai_plan(
        task = "annotate",
        actions = list(list(id = "a1", action = "run_annotation",
            params = list(ref = "BlueprintENCODEData", labels = "label.main")))
    )
    err <- tryCatch(ValidateAIPlan(plan, object = sce, registry = registry, strict = TRUE),
        error = function(e) e)
    expect_s3_class(err, "sclet_ai_invalid_plan")
    expect_true(grepl("reference_invalid", conditionMessage(err), fixed = TRUE))
    expect_true(grepl("not a dataset exported by celldex",
        conditionMessage(err), fixed = TRUE))
})

# ---- Regressions for the 2026-10-16 benchmark findings (sclet-benchmark-ai) --

test_that("a zero-action plan is refused at validation and at execution", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 3, ncol = 2))
    )
    registry <- AIDefaultExecutionRegistry(sce, include = "read")
    plan <- new_sclet_ai_plan(
        task = "empty_plan_regression",
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
        actions = list()
    )
    validation <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_false(validation$valid)
    expect_true(any(grepl("plan_empty", validation$errors)))

    # The backstop: even a hand-built validation object that claims the empty
    # plan is valid must not execute to status="completed" with zero results,
    # which is what a repair attempt that returned nothing used to produce.
    forged <- validation
    forged$valid <- TRUE
    forged$confirmation_token <- "sclet-confirm-forged"
    expect_error(
        ExecuteAIPlan(sce, plan, registry, validation = forged,
            dry_run = FALSE, confirmation = forged$confirmation_token),
        class = "sclet_ai_empty_plan"
    )
})

test_that("run_milo accepts a PCA promised earlier in the same plan", {
    # The gate read reducedDimNames(object) only, so every chained milo plan
    # was refused with reduction_missing at validation although run_pca had
    # not executed yet -- the deferred projection 54248b3 wired into the
    # enrichment and fate gates but not this one.
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 3, ncol = 2))
    )
    registry <- AIDefaultExecutionRegistry(sce, include = "all")
    chained <- new_sclet_ai_plan(
        task = "milo_chain_regression",
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
        actions = list(
            list(id = "g1", action = "normalize_data"),
            list(id = "g2", action = "find_variable_features", depends_on = "g1"),
            list(id = "g3", action = "scale_data", depends_on = "g2"),
            list(id = "g4", action = "run_pca", params = list(ncomponents = 2L), depends_on = "g3"),
            list(id = "g5", action = "run_milo",
                params = list(sample_col = "s", condition_col = "c", contrasts = "a - b"),
                depends_on = "g4")
        )
    )
    v_chained <- ValidateAIPlan(chained, object = sce, registry = registry)
    expect_false(any(grepl("reduction_missing", v_chained$errors)))
    # the plan still fails on its unknown design columns -- that gate is the
    # point -- but the failure must name the columns, not the reduction
    expect_true(any(grepl("not_in_colData|_missing", v_chained$errors)))

    # Control: without a planned PCA the reduction gate still fires.
    bare <- new_sclet_ai_plan(
        task = "milo_bare_regression",
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
        actions = list(
            list(id = "g1", action = "run_milo",
                params = list(sample_col = "s", condition_col = "c", contrasts = "a - b"))
        )
    )
    v_bare <- ValidateAIPlan(bare, object = sce, registry = registry)
    expect_true(any(grepl("reduction_missing", v_bare$errors)))
})

test_that("run_de_test exposes detest_id as a referenceable output field", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(1, nrow = 3, ncol = 2))
    )
    registry <- AIDefaultExecutionRegistry(sce, include = "all")
    allowed <- sclet:::sclet_ai_referenceable_output_fields(registry[["run_de_test"]])
    expect_true("detest_id" %in% allowed)

    # The reference the benchmark measured as unreachable must now validate:
    # run_enrichment(detest_id = ${step.output.detest_id}) after run_de_test.
    chained <- new_sclet_ai_plan(
        task = "detest_id_regression",
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
        actions = list(
            list(id = "g1", action = "normalize_data"),
            list(id = "g2", action = "find_variable_features", depends_on = "g1"),
            list(id = "g3", action = "scale_data", depends_on = "g2"),
            list(id = "g4", action = "run_pca", params = list(ncomponents = 2L), depends_on = "g3"),
            list(id = "g5", action = "find_neighbors", params = list(dims = 1:2), depends_on = "g4"),
            list(id = "g6", action = "find_clusters", depends_on = "g5"),
            list(id = "g7", action = "run_de_test", params = list(name = "markers"), depends_on = "g6"),
            list(id = "g8", action = "run_enrichment",
                params = list(detest_id = "${g7.output.detest_id}", name = "enrich_go"),
                depends_on = "g7")
        )
    )
    v <- ValidateAIPlan(chained, object = sce, registry = registry)
    expect_false(any(grepl("value_reference_unresolvable", v$errors)))
    expect_false(any(grepl("detest_state_missing", v$errors)))
})

test_that("an unresolvable citation in the plan narrative does not abort the analysis", {
    ## Real-data run 7 (PBMC 10k, 2026-10-18): the plan narrative cited
    ## `sclet-ai-*` (AI result ids) rather than `ev:*`; RecordAIResult correctly
    ## refused, and because the plan-recording call was not wrapped the refusal
    ## escaped RunAIAnalysis, turning a completed 9/9 execution into a failed
    ## run. The analysis must survive; the unrecorded narrative is reported.
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(stats::rpois(30L * 12L, 5), nrow = 30L, ncol = 12L,
            dimnames = list(paste0("gene_", seq_len(30L)), paste0("c", seq_len(12L)))))
    )
    old <- options(sclet.ai.call = function(task, context, ...) {
        if (identical(task, "analysis_plan")) {
            return(new_sclet_ai_result(
                "analysis_plan",
                answer = "normalize the counts",
                context = list(schema_version = "1.0", fingerprint = "test"),
                findings = list(list(claim_level = "hypothesis",
                    evidence_refs = "sclet-ai-deadbeef")),
                proposed_actions = list(list(id = "s1", action = "normalize_data",
                    params = list()))
            ))
        }
        new_sclet_ai_result("analysis_explanation", answer = "ok",
            context = list(schema_version = "1.0", fingerprint = "test"))
    })
    on.exit(options(old), add = TRUE)

    result <- RunAIAnalysis(sce, goal = "Normalize the counts in this dataset.",
        model = "mock", confirm = TRUE, include = "all", interpret = FALSE,
        record = TRUE)
    expect_equal(result$status, "completed")
    ## The bad citation is stripped and the questionable finding downgraded, so
    ## the plan narrative is recorded instead of dropped. The dropped id is
    ## reported for audit rather than silently vanishing.
    expect_null(result$report$plan_record_error)
    expect_true(any(grepl("^ai_plan_", names(GetAnalysisLedger(result$object)$analyses))))
    recorded <- GetAnalysisLedger(result$object)$analyses[[
        grep("^ai_plan_", names(GetAnalysisLedger(result$object)$analyses))[[1L]]]]
    expect_true(grepl("sclet-ai-deadbeef", recorded$summary$dropped_evidence_refs,
        fixed = TRUE))
    expect_equal(recorded$summary$n_findings_with_refs, 0L)
})

test_that("the recorded round keeps the plan's topology", {
    # The failure this pins: the plan is a DAG, but the record used to be a
    # flat array (`depends_on` appeared nowhere in R/ai-execution.R), so a
    # finished round could not answer "which step consumed which earlier
    # step's output". `AnalysisGraph()` reads these fields; if a refactor
    # drops them the graph silently degrades to an unconnected trace.
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(stats::rpois(40, 5), nrow = 10L, ncol = 4L))
    )
    write_record <- function(name, id) {
        AIAction(
            name = name,
            handler = function(object, params) {
                sclet:::sclet_set_analysis_state(object, type = "preprocess",
                    id = id, method = "trace_pin", summary = list(n = 1))
            },
            returns = "sce",
            allowed_state_types = "preprocess"
        )
    }
    registry <- AIExecutionRegistry(list(
        step_a = write_record("step_a", "pin_a"),
        step_b = write_record("step_b", "pin_b")
    ))
    plan <- sclet:::new_sclet_ai_plan(
        task = "topology_pin",
        actions = list(
            list(id = "s1", action = "step_a"),
            list(id = "s2", action = "step_b", depends_on = "s1")
        ),
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint
    )
    validation <- ValidateAIPlan(plan, object = sce, registry = registry)
    executed <- ExecuteAIPlan(sce, plan, registry, validation = validation,
        dry_run = FALSE, confirmation = validation$confirmation_token)

    record <- sclet:::sclet_get_analysis(executed$object, executed$execution_id)
    steps <- record$params$steps
    expect_equal(length(steps), 2L)
    expect_identical(steps[[1L]]$depends_on, character())
    expect_identical(steps[[2L]]$depends_on, "s1")
    expect_identical(steps[[2L]]$record_ids, "preprocess::pin_b")
    expect_true(all(vapply(steps, function(s) is.numeric(s$attempts), logical(1L))))
})
