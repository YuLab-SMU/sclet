# The six mainlines that had no AI action: registry, gates, records, evidence.
#
# The contract these tests defend is the same for all six: an action either
# refuses with a named gate, or runs and leaves a record plus one evidence node
# whose values are numbers or logicals. "Registered" without a gate or an
# evidence node would be a name in a catalog, not a capability.

ml_sce <- function(ncells = 60, ngenes = 120) {
    sce <- scuttle::mockSCE(ncells = ncells, ngenes = ngenes)
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = SummarizedExperiment::assay(sce, "counts")),
        colData = S4Vectors::DataFrame(SummarizedExperiment::colData(sce))
    )
    sce
}

ml_sce_normalized <- function(...) {
    sce <- ml_sce(...)
    SummarizedExperiment::assay(sce, "logcounts") <-
        log1p(as.matrix(SummarizedExperiment::assay(sce, "counts")))
    sce
}

ml_action <- function(sce, action, include = "all") {
    AIDefaultExecutionRegistry(sce, include = include)[[action]]
}

ml_evidence <- function(object, id) {
    sclet:::sclet_ai_evidence_get_all(object)[[id]]
}

ml_gate <- function(action, object, params = list()) {
    action$prerequisites(object, params)
}

# A stub with the real signature, writing the record the real function writes.
ml_augur_stub <- function(auc = data.frame(cell_type = c("a", "b"), auc = c(0.57, 0.36))) {
    function(sce, condition_col, group_col = NULL, layer = NULL, assay = NULL,
             name = "augur", ...) {
        sclet:::sclet_set_analysis_state(
            sce,
            type = "priority", id = name, method = "Augur",
            artifacts = list(analysis_key = "augur", AUC = auc),
            summary = list()
        )
    }
}

ml_phenotype_stub <- function() {
    function(sce, bulk_matrix, phenotype_data, family = "cox", ...) {
        sclet:::sclet_set_analysis_state(
            sce,
            type = "phenotype", id = "scPAS", method = "scPAS::scPAS",
            summary = list(added_columns = c("col_a", "col_b"))
        )
    }
}

ml_scenic_stub <- function() {
    function(sce, tfs_path, motif_annotations_path, database_paths,
             assay_use = "counts", num_workers = 1L, seed = 123L,
             name = "scenic", output_dir = NULL, ...) {
        sclet:::sclet_set_state_record(
            sce, type = "scenic", id = name, active = TRUE,
            value = list(
                method = "pySCENIC",
                inputs = list(assay = assay_use),
                artifacts = list(altExp = "SCENIC_AUC"),
                params = list(num_workers = num_workers, seed = seed),
                summary = list(n_regulons = 7L, n_cells = ncol(sce)),
                created_at = Sys.time()
            )
        )
    }
}

ml_oracle_stub <- function() {
    function(sce, target_gene, perturbation_value = 0.0, base_grn_path,
             cluster_key = NULL, reduction = "PCA", ...) {
        sclet:::sclet_set_state_record(
            sce, type = "perturbation",
            id = paste0("celloracle_", target_gene), active = TRUE,
            value = list(
                method = "celloracle::Oracle",
                inputs = list(assay = "counts", reduction = reduction),
                artifacts = list(shift_reduction = "X_perturbation"),
                params = list(target_gene = target_gene,
                    perturbation_value = perturbation_value),
                summary = list(n_cells = ncol(sce)),
                created_at = Sys.time()
            )
        )
    }
}

test_that("the six mainline actions are registered under their own groups", {
    sce <- ml_sce()
    groups <- c(
        run_geneset_scoring = "geneset",
        run_supercell = "metacell",
        run_perturbation_priority = "priority",
        run_phenotype_association = "phenotype",
        run_scenic = "regulon",
        run_in_silico_perturbation = "perturbation"
    )
    all_actions <- AIDefaultExecutionRegistry(sce, include = "all")
    # 40 through the six mainline actions; the record-only milestone adds six
    # more (imputation, malignant classifier, two mappings, deconvolution, RegVelo).
    # 40 through the six mainline actions; the record-only milestone adds six
    # more (imputation, malignant classifier, two mappings, deconvolution, RegVelo).
    expect_equal(length(all_actions), 46L)
    expect_setequal(names(groups), intersect(names(groups), names(all_actions)))
    other_groups <- setdiff(unique(groups), character())
    for (action in names(groups)) {
        in_group <- AIDefaultExecutionRegistry(sce, include = groups[[action]])
        expect_true(action %in% names(in_group), info = action)
        # Exclusivity: an action registered in every group would pass the
        # membership check above and make every include request carry it.
        for (other in setdiff(other_groups, groups[[action]])) {
            expect_false(action %in% names(AIDefaultExecutionRegistry(sce, include = other)),
                info = paste(action, "leaked into", other))
        }
        descriptor <- all_actions[[action]]
        expect_identical(descriptor$returns, "sce")
        expect_true(descriptor$mutates_object)
        expect_true(descriptor$requires_confirmation)
        expect_false(descriptor$idempotent)
        expect_identical(isTRUE(descriptor$replaces_object),
            identical(action, "run_supercell"))
        expect_equal(descriptor$allowed_state_types[1L],
            c(run_geneset_scoring = "geneset_scoring",
              run_supercell = "aggregation",
              run_perturbation_priority = "priority",
              run_phenotype_association = "phenotype",
              run_scenic = "scenic",
              run_in_silico_perturbation = "perturbation")[[action]])
        expect_true("ai_evidence" %in% descriptor$allowed_state_types)
    }
})

test_that("only the Augur action is offered to the planner", {
    sce <- ml_sce()
    catalog <- sclet:::sclet_ai_describe_registry(
        AIDefaultExecutionRegistry(sce, include = "all")
    )
    expect_true("run_perturbation_priority" %in% names(catalog))
    hidden <- c("run_geneset_scoring", "run_supercell",
        "run_phenotype_association", "run_scenic", "run_in_silico_perturbation")
    expect_false(any(hidden %in% names(catalog)))
    # The registry still carries them for validation/execution.
    expect_true(all(hidden %in% names(AIDefaultExecutionRegistry(sce, include = "all"))))
})

test_that("every registered action is discoverable through the group hint", {
    sce <- ml_sce()
    actions <- names(AIDefaultExecutionRegistry(sce, include = "all"))
    for (action in actions) {
        groups <- sclet:::sclet_ai_action_groups(sce, action)
        expect_true(length(groups) >= 1L, info = action)
    }
})

test_that("each mainline refuses a fixture that lacks its input, naming the gate", {
    # gene sets
    geneset <- ml_action(ml_sce(), "run_geneset_scoring")
    expect_match(ml_gate(geneset, ml_sce(), list(gene_sets = list())),
        "gene_sets_empty")
    expect_match(ml_gate(geneset, ml_sce(),
        list(gene_sets = list(only = c("NOT_A_GENE_1", "NOT_A_GENE_2")))),
        "gene_sets_not_found")

    # metacells need a normalized assay, and enough cells for two metacells
    supercell <- ml_action(ml_sce(), "run_supercell")
    expect_match(ml_gate(supercell, ml_sce(), list()),
        "logcounts_assay_missing")
    small <- ml_sce_normalized(ncells = 20)
    expect_match(ml_gate(supercell, small, list(gamma = 20)),
        "too_few_metacells")

    # Augur: the condition column is the user's design
    augur <- ml_action(ml_sce(), "run_perturbation_priority")
    expect_match(ml_gate(augur, ml_sce(), list(condition_col = "nope")),
        "condition_col_missing")
    one_level <- ml_sce()
    SummarizedExperiment::colData(one_level)$only_one <- "a"
    expect_match(ml_gate(augur, one_level, list(condition_col = "only_one")),
        "condition_levels_insufficient")
    expect_match(ml_gate(augur, ml_sce(), list(condition_col = "Treatment")),
        "group_column_missing")

    # phenotype: scPAS is a required, absent backend here
    phenotype <- ml_action(ml_sce(), "run_phenotype_association")
    expect_match(ml_gate(phenotype, ml_sce(),
        list(bulk_matrix = matrix(1, 2, 2), phenotype_data = data.frame(a = 1))),
        "optional_package_missing: scPAS")

    # SCENIC: the three file inputs are the user's
    scenic <- ml_action(ml_sce(), "run_scenic")
    expect_match(ml_gate(scenic, ml_sce(),
        list(tfs_path = "nope.txt", motif_annotations_path = "nope2.txt",
            database_paths = "nope3.feather")),
        "scenic_files_missing")

    # CellOracle: target gene, base GRN, identity, reduction
    oracle <- ml_action(ml_sce(), "run_in_silico_perturbation")
    expect_match(ml_gate(oracle, ml_sce(),
        list(target_gene = "NOT_A_GENE", base_grn_path = "nope.txt")),
        "target_gene_not_found")
    grn <- tempfile(fileext = ".txt")
    writeLines("x", grn)
    on.exit(unlink(grn), add = TRUE)
    expect_match(ml_gate(oracle, ml_sce(),
        list(target_gene = rownames(ml_sce())[1L], base_grn_path = "nope.txt")),
        "base_grn_missing")
    expect_match(ml_gate(oracle, ml_sce(),
        list(target_gene = rownames(ml_sce())[1L], base_grn_path = grn)),
        "cluster_identity_missing")
})

test_that("the identity gate accepts what the callees accept and rejects what they reject", {
    sce <- ml_sce()
    SingleCellExperiment::colLabels(sce) <- rep(c("a", "b"), length.out = ncol(sce))
    sce <- RunPCA(sce, ncomponents = 5)
    grn <- tempfile(fileext = ".txt")
    on.exit(unlink(grn), add = TRUE)
    writeLines("x", grn)

    # Augur: ActiveIdent() is enough; the gate used to demand colLabels itself.
    augur <- ml_action(sce, "run_perturbation_priority")
    expect_true(isTRUE(ml_gate(augur, sce, list(condition_col = "Treatment"))))
    expect_match(ml_gate(augur, sce, list(condition_col = "Treatment", group_col = "nope")),
        "group_col_unknown")

    # CellOracle: the sentinel "colLabels" is not a colData column, so the gate
    # must resolve it to the real one ("label") and stay silent.
    oracle <- ml_action(sce, "run_in_silico_perturbation")
    target <- rownames(sce)[1L]
    expect_true(isTRUE(ml_gate(oracle, sce,
        list(target_gene = target, base_grn_path = grn))))
    expect_match(ml_gate(oracle, sce,
        list(target_gene = target, base_grn_path = grn, cluster_key = "nope")),
        "cluster_key_unknown")

    # ... and without any identity at all it still refuses.
    bare <- ml_sce()
    SingleCellExperiment::colLabels(bare) <- NULL
    expect_match(ml_gate(oracle, bare,
        list(target_gene = rownames(bare)[1L], base_grn_path = grn)),
        "cluster_identity_missing")

    # A reduction that does not exist is its own gate, not a Python failure.
    no_reduction <- ml_sce()
    SingleCellExperiment::colLabels(no_reduction) <- rep(c("a", "b"),
        length.out = ncol(no_reduction))
    expect_match(ml_gate(oracle, no_reduction,
        list(target_gene = rownames(no_reduction)[1L], base_grn_path = grn)),
        "reduction_missing")

    # An assay named by the user must exist.
    geneset <- ml_action(ml_sce(), "run_geneset_scoring")
    genes <- rownames(ml_sce())
    expect_match(ml_gate(geneset, ml_sce(),
        list(gene_sets = list(a = genes[1:3]), assay_use = "nope")),
        "assay_missing")
})

test_that("an object-replacing action cannot share a plan with other steps", {
    sce <- ml_sce_normalized(ncells = 60)
    sce <- RunPCA(sce, ncomponents = 5)
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "metacell"))
    plan <- sclet:::new_sclet_ai_plan(
        task = "supercell_mid_plan",
        actions = list(
            list(id = "s1", action = "inspect_status"),
            list(id = "s2", action = "run_supercell", params = list(gamma = 20))
        ),
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint
    )
    validation <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_false(isTRUE(validation$valid))
    expect_match(paste(unlist(validation$errors), collapse = " | "),
        "object_replacement_not_standalone")
})

test_that("gene-set scoring runs, records the analysis and a numeric evidence node", {
    skip_if_not_installed("UCell")
    sce <- ml_sce()
    genes <- rownames(sce)
    action <- ml_action(sce, "run_geneset_scoring")
    result <- action$handler(sce, list(
        gene_sets = list(set_a = genes[1:5], set_b = genes[6:10]),
        method = "UCell",
        name = "my_sets"
    ))
    record <- sclet:::sclet_get_state_records(result, "geneset_scoring")[["my_sets"]]
    expect_false(is.null(record))
    evidence <- ml_evidence(result, "ev:geneset_my_sets")
    expect_false(is.null(evidence))
    expect_identical(evidence$claim_level, "measured")
    expect_equal(evidence$values$n_sets, 2L)
    expect_equal(evidence$values$n_genes_used, 10L)
    expect_true(all(vapply(evidence$values, is.numeric, logical(1L))))
    expect_identical(evidence$source, "my_sets")
})

test_that("metacell aggregation runs and records its own, replaced, object", {
    skip_if_not_installed("SuperCell")
    sce <- ml_sce_normalized(ncells = 60)
    sce <- sclet:::sclet_set_analysis_state(sce, type = "clustering", id = "before_metacells",
        method = "fixture", summary = list(n = 1))
    before <- ncol(sce)
    action <- ml_action(sce, "run_supercell")
    result <- action$handler(sce, list(gamma = 20, k.knn = 5))
    expect_lt(ncol(result), before)
    evidence <- ml_evidence(result, "ev:supercell_supercell")
    expect_false(is.null(evidence))
    expect_equal(evidence$values$n_cells, before)
    expect_gte(evidence$values$n_metacells, 2)
    # The ledger really is replaced, which is why the action is driver-only.
    expect_null(sclet:::sclet_get_state_records(result, "clustering")[["before_metacells"]])
    # The evidence source resolves to the analysis record, not to the state
    # record that shares the bare id "supercell".
    source_of <- AnalysisGraph(result)$edges
    source_of <- source_of[source_of$kind == "source_of", , drop = FALSE]
    expect_true("record:analysis::supercell" %in% source_of$to)
})

test_that("Augur priority runs and records AUC aggregates at 'associated'", {
    sce <- ml_sce()
    action <- ml_action(sce, "run_perturbation_priority")
    result <- with_mocked_bindings(
        action$handler(sce, list(condition_col = "Treatment")),
        RunPerturbationPriority = ml_augur_stub(),
        .package = "sclet"
    )
    evidence <- ml_evidence(result, "ev:augur_augur")
    expect_false(is.null(evidence))
    expect_identical(evidence$claim_level, "associated")
    expect_equal(evidence$values$n_cell_types, 2L)
    expect_equal(evidence$values$auc_max, 0.57)
    expect_equal(evidence$values$auc_min, 0.36)
    expect_equal(evidence$values$auc_median, stats::median(c(0.57, 0.36)))
    expect_identical(evidence$source, "augur")
})

test_that("phenotype association records external inputs as 'observed'", {
    sce <- ml_sce()
    action <- ml_action(sce, "run_phenotype_association")
    result <- with_mocked_bindings(
        action$handler(sce, list(bulk_matrix = matrix(1, 2, 2),
            phenotype_data = data.frame(a = 1))),
        RunPhenotypeAssociation = ml_phenotype_stub(),
        .package = "sclet"
    )
    evidence <- ml_evidence(result, "ev:phenotype_association")
    expect_false(is.null(evidence))
    expect_identical(evidence$claim_level, "observed")
    expect_equal(evidence$values$n_columns_added, 2L)
    expect_true(isTRUE(evidence$values$external_inputs))
    expect_identical(evidence$source, "scPAS")
})

test_that("SCENIC records regulon counts and flags the Python backend", {
    sce <- ml_sce()
    action <- ml_action(sce, "run_scenic")
    tfs <- tempfile(fileext = ".txt")
    motif <- tempfile(fileext = ".tbl")
    database <- tempfile(fileext = ".feather")
    on.exit(unlink(c(tfs, motif, database)), add = TRUE)
    for (path in c(tfs, motif, database)) writeLines("x", path)
    result <- with_mocked_bindings(
        action$handler(sce, list(tfs_path = tfs, motif_annotations_path = motif,
            database_paths = database)),
        RunSCENIC = ml_scenic_stub(),
        .package = "sclet"
    )
    evidence <- ml_evidence(result, "ev:scenic_scenic")
    expect_false(is.null(evidence))
    expect_equal(evidence$values$n_regulons, 7L)
    expect_true(isTRUE(evidence$values$python_backend))
})

test_that("in-silico perturbation records a simulation at 'associated'", {
    sce <- ml_sce()
    target <- rownames(sce)[1L]
    action <- ml_action(sce, "run_in_silico_perturbation")
    grn <- tempfile(fileext = ".txt")
    on.exit(unlink(grn), add = TRUE)
    writeLines("x", grn)
    result <- with_mocked_bindings(
        action$handler(sce, list(target_gene = target, base_grn_path = grn,
            perturbation_value = 0.5)),
        RunInSilicoPerturbation = ml_oracle_stub(),
        .package = "sclet"
    )
    node_id <- paste0("ev:perturbation_",
        sclet:::sclet_ai_evidence_id_token(target))
    evidence <- ml_evidence(result, node_id)
    expect_false(is.null(evidence))
    expect_identical(evidence$claim_level, "associated")
    expect_true(isTRUE(evidence$values$simulated))
    expect_equal(evidence$values$perturbation_value, 0.5)
    expect_identical(evidence$source, paste0("celloracle_", target))
})

test_that("a mainline plan round-trips through validation, execution and the trace graph", {
    sce <- ml_sce()
    # Augur's grouping gate needs an identity; a plan cannot pass a cell-type
    # column that does not exist, so the fixture clusters first.
    SingleCellExperiment::colLabels(sce) <- rep(c("a", "b"), length.out = ncol(sce))
    registry <- AIDefaultExecutionRegistry(sce, include = "priority")
    plan <- sclet:::new_sclet_ai_plan(
        task = "mainline_priority",
        actions = list(list(id = "s1", action = "run_perturbation_priority",
            params = list(condition_col = "Treatment"))),
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint
    )
    validation <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_true(validation$valid)
    executed <- with_mocked_bindings(
        ExecuteAIPlan(sce, plan, registry, validation = validation,
            dry_run = FALSE, confirmation = validation$confirmation_token),
        RunPerturbationPriority = ml_augur_stub(),
        .package = "sclet"
    )
    expect_equal(executed$status, "completed")
    graph <- AnalysisGraph(executed$object)
    produced <- graph$edges[graph$edges$kind == "produced", , drop = FALSE]
    expect_true("record:priority::augur" %in% produced$to)
    expect_true("evidence:ev:augur_augur" %in% produced$to)
})

test_that("a user-owned input is refused by validation, not silently defaulted", {
    sce <- ml_sce()
    registry <- AIDefaultExecutionRegistry(sce, include = "geneset")
    plan <- sclet:::new_sclet_ai_plan(
        task = "missing_gene_sets",
        actions = list(list(id = "s1", action = "run_geneset_scoring",
            params = list(method = "UCell"))),
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint
    )
    validation <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_false(isTRUE(validation$valid))
    expect_true(any(grepl("gene_sets", unlist(validation$errors), fixed = TRUE)))
})

test_that("the default round still includes only the basic groups", {
    expect_identical(
        eval(formals(RunAIAnalysis)$include),
        c("read", "preprocess", "dimred", "graph", "cluster")
    )
})

test_that("the skill tells the planner which actions are not in its catalog", {
    block <- sclet:::sclet_ai_skill_context(NULL)
    expect_match(block, "Actions that are not in your catalog", fixed = TRUE)
    expect_match(block, "run_geneset_scoring", fixed = TRUE)
    expect_match(block, "run_supercell", fixed = TRUE)
})

test_that("the oracle action passes the real identity column, not the sentinel", {
    # Reverting `sclet_ai_mainline_identity_column()` to return "colLabels"
    # leaves every other test green while the callee receives a zero-length
    # label vector and fails inside Python. Capture what it is handed.
    captured <- new.env(parent = emptyenv())
    stub <- function(sce, target_gene, perturbation_value = 0.0, base_grn_path,
                     cluster_key = NULL, reduction = "PCA", ...) {
        captured$cluster_key <- cluster_key
        sclet:::sclet_set_state_record(
            sce, type = "perturbation", id = paste0("celloracle_", target_gene),
            active = TRUE,
            value = list(
                method = "celloracle::Oracle",
                inputs = list(assay = "counts", reduction = reduction),
                artifacts = list(shift_reduction = "X_perturbation"),
                params = list(target_gene = target_gene,
                    perturbation_value = perturbation_value),
                summary = list(n_cells = ncol(sce)),
                created_at = Sys.time()
            )
        )
    }
    grn <- tempfile(fileext = ".txt")
    on.exit(unlink(grn), add = TRUE)
    writeLines("x", grn)

    expect_identity <- function(sce, expected) {
        action <- ml_action(sce, "run_in_silico_perturbation")
        with_mocked_bindings(
            action$handler(sce, list(target_gene = rownames(sce)[1L],
                base_grn_path = grn)),
            RunInSilicoPerturbation = stub,
            .package = "sclet"
        )
        expect_identical(captured$cluster_key, expected)
    }

    by_labels <- ml_sce()
    SingleCellExperiment::colLabels(by_labels) <- rep(c("a", "b"),
        length.out = ncol(by_labels))
    by_labels <- RunPCA(by_labels, ncomponents = 5)
    expect_identity(by_labels, "label")

    by_column <- ml_sce()
    SummarizedExperiment::colData(by_column)$grp <- rep(c("a", "b"),
        length.out = ncol(by_column))
    ActiveIdent(by_column) <- "grp"
    by_column <- RunPCA(by_column, ncomponents = 5)
    expect_identity(by_column, "grp")
})

test_that("the priority gate sees an identity the same plan will create", {
    sce <- ml_sce()
    registry <- AIDefaultExecutionRegistry(sce, include = "all")
    plan <- sclet:::new_sclet_ai_plan(
        task = "cluster_then_priority",
        actions = list(
            list(id = "s1", action = "normalize_data"),
            list(id = "s2", action = "find_variable_features",
                params = list(nfeatures = 100L), depends_on = "s1"),
            list(id = "s3", action = "scale_data", depends_on = "s2"),
            list(id = "s4", action = "run_pca", params = list(ncomponents = 10L),
                depends_on = "s3"),
            list(id = "s5", action = "find_neighbors", params = list(dims = 1:10),
                depends_on = "s4"),
            list(id = "s6", action = "find_clusters", depends_on = "s5"),
            list(id = "s7", action = "run_perturbation_priority",
                params = list(condition_col = "Treatment"), depends_on = "s6")
        ),
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint
    )
    validation <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_true(validation$valid)
    expect_false(any(grepl("group_column_missing", unlist(validation$errors))))
})

test_that("an object-replacing round is declared in the graph, not silently truncated", {
    sce <- ml_sce_normalized(ncells = 60)
    registry <- AIDefaultExecutionRegistry(sce, include = "metacell")
    plan <- sclet:::new_sclet_ai_plan(
        task = "supercell_round",
        actions = list(list(id = "s1", action = "run_supercell", params = list(gamma = 20))),
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint
    )
    validation <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_true(validation$valid)
    executed <- ExecuteAIPlan(sce, plan, registry, validation = validation,
        dry_run = FALSE, confirmation = validation$confirmation_token)
    graph <- AnalysisGraph(executed$object)
    expect_true(isTRUE(graph$meta$object_replaced))
    # A round that replaced nothing says so.
    basic <- AnalysisGraph(ml_sce())
    expect_false(isTRUE(basic$meta$object_replaced))
})

test_that("a same-id re-run is a refreshed edge, not a stale producer", {
    # The Python seam is mocked, NOT RunInSilicoPerturbation: the real writer
    # must run, because the point of this test is that its record carries the
    # `created_at` the trace diff needs. Mocking the exported function instead
    # would make the test pass with that stamp removed.
    sce <- ml_sce()
    SingleCellExperiment::colLabels(sce) <- rep(c("a", "b"), length.out = ncol(sce))
    sce <- RunPCA(sce, ncomponents = 5)
    target <- rownames(sce)[1L]
    grn <- tempfile(fileext = ".txt")
    on.exit(unlink(grn), add = TRUE)
    writeLines("x", grn)
    registry <- AIDefaultExecutionRegistry(sce, include = "perturbation")
    plan <- sclet:::new_sclet_ai_plan(
        task = "perturb_twice",
        actions = list(
            list(id = "s1", action = "run_in_silico_perturbation",
                params = list(target_gene = target, base_grn_path = grn)),
            list(id = "s2", action = "run_in_silico_perturbation",
                params = list(target_gene = target, base_grn_path = grn),
                depends_on = "s1")
        ),
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint
    )
    validation <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_true(validation$valid)
    fake_basilisk <- function(env, fun, ...) {
        args <- list(...)
        # `emb` is cells x dims; the callee assigns rownames(co_res) <- cells.
        matrix(0, nrow = nrow(args$emb), ncol = 2L)
    }
    executed <- with_mocked_bindings(
        ExecuteAIPlan(sce, plan, registry, validation = validation,
            dry_run = FALSE, confirmation = validation$confirmation_token),
        sclet_basilisk_run = fake_basilisk,
        .package = "sclet"
    )
    expect_equal(executed$status, "completed")
    graph <- AnalysisGraph(executed$object)
    record_id <- paste0("record:perturbation::celloracle_", target)
    produced <- graph$edges$from[graph$edges$kind == "produced" & graph$edges$to == record_id]
    refreshed <- graph$edges$from[graph$edges$kind == "refreshed" & graph$edges$to == record_id]
    expect_length(produced, 1L)
    expect_length(refreshed, 1L)
    expect_identical(produced, paste0("step:", executed$execution_id, ":s1"))
    expect_identical(refreshed, paste0("step:", executed$execution_id, ":s2"))
    expect_setequal(graph_producers(graph, record_id)$ids,
        paste0("step:", executed$execution_id, ":s", 1:2))
})

test_that("a plan can satisfy the oracle reduction gate with its own earlier step", {
    sce <- ml_sce()
    target <- rownames(sce)[1L]
    grn <- tempfile(fileext = ".txt")
    on.exit(unlink(grn), add = TRUE)
    writeLines("x", grn)
    registry <- AIDefaultExecutionRegistry(sce, include = "all")
    plan <- sclet:::new_sclet_ai_plan(
        task = "pca_then_perturb",
        actions = list(
            list(id = "s1", action = "normalize_data"),
            list(id = "s2", action = "find_variable_features",
                params = list(nfeatures = 100L), depends_on = "s1"),
            list(id = "s3", action = "scale_data", depends_on = "s2"),
            list(id = "s4", action = "run_pca", params = list(ncomponents = 10L),
                depends_on = "s3"),
            list(id = "s5", action = "find_neighbors", params = list(dims = 1:10),
                depends_on = "s4"),
            list(id = "s6", action = "find_clusters", depends_on = "s5"),
            list(id = "s7", action = "run_in_silico_perturbation",
                params = list(target_gene = target, base_grn_path = grn),
                depends_on = "s6")
        ),
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint
    )
    validation <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_false(any(grepl("reduction_missing", unlist(validation$errors))))
    expect_false(any(grepl("cluster_identity_missing", unlist(validation$errors))))
})
