# The record-only mainlines: registered, gated, recorded, evidenced.
#
# Same contract as test-ai-registry-mainlines.R, one difference: most of these
# inputs live on the object, so several actions are planner-visible; the ones
# that need a SECOND object or a prior network are driver-only.

ro_sce <- function(ncells = 80, ngenes = 120) {
    set.seed(11)
    counts <- matrix(stats::rpois(ngenes * ncells, 6), nrow = ngenes, ncol = ncells,
        dimnames = list(paste0("gene_", seq_len(ngenes)), paste0("cell_", seq_len(ncells))))
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
    SummarizedExperiment::assay(sce, "logcounts") <- log1p(as.matrix(counts))
    SummarizedExperiment::colData(sce)$sample <- rep(c("s1", "s2", "s3", "s4"),
        length.out = ncells)
    SummarizedExperiment::colData(sce)$condition <- rep(c("ctrl", "stim"),
        length.out = ncells)
    sce
}

ro_ref <- function(ngenes = 120, ncells = 40) {
    set.seed(12)
    counts <- matrix(stats::rpois(ngenes * ncells, 6), nrow = ngenes, ncol = ncells,
        dimnames = list(paste0("gene_", seq_len(ngenes)), paste0("ref_", seq_len(ncells))))
    ref <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
    SummarizedExperiment::assay(ref, "logcounts") <- log1p(as.matrix(counts))
    SummarizedExperiment::colData(ref)$label <- rep(c("T", "B"), length.out = ncells)
    ref
}

ro_gate <- function(action, object, params = list()) {
    action$prerequisites(object, params)
}

ro_malignant_stub <- function(record_name = "malignant") {
    function(object, seed, method = c("iterative", "cnv"), cells = NULL,
             sample_col = "sample", assay = NULL, n_top = 50, padj = 0.01, ...,
             name = record_name) {
        stopifnot(is.list(seed), !is.null(seed$tumor), !is.null(seed$normal))
        sclet:::sclet_set_state_record(object, type = "malignancy", id = name,
            active = TRUE,
            value = list(
                method = "iterative", inputs = list(sample_col = sample_col),
                artifacts = list(columns = paste0(name, "_classification")),
                summary = list(n_malignant = 12, n_non_malignant = 8,
                    malignant_fraction = 0.6, iterations = 3, converged = TRUE),
                created_at = Sys.time()
            ))
    }
}

ro_mapping_stub <- function(method = "KNNPredict", record_name = "knn_map") {
    function(object, ref, labels = NULL, features = NULL, layer = NULL,
             k = 5, vars = NULL, ...) {
        # Mirror the real record's summary for the method: KNN writes
        # n_labels/n_reference_labels, Symphony writes n_query_cells/n_ref_cells.
        summary <- if (identical(method, "Symphony")) {
            list(n_query_cells = ncol(object), n_ref_cells = ncol(ref),
                n_common_genes = 100L)
        } else {
            list(n_labels = 2L, n_reference_labels = 2L)
        }
        sclet:::sclet_set_analysis_state(object, type = "mapping", id = record_name,
            method = method, summary = summary)
    }
}

ro_regvelo_stub <- function(record_name = "regvelo") {
    function(sce, grn, ...) {
        # The real callee builds the prior from a matrix/data.frame and rejects a
        # named list; the stub asserts the shape so a wrong gate cannot hide.
        stopifnot(is.matrix(grn) || methods::is(grn, "Matrix") || is.data.frame(grn))
        sclet:::sclet_set_state_record(sce, type = "velocity", id = record_name,
            active = TRUE,
            value = list(
                method = "regvelo::REGVELOVI",
                inputs = list(grn = list(n_targets = 40L, n_regulators = 6L, n_edges = 120L)),
                artifacts = list(velocity_assay = "velocity_regvelo"),
                summary = list(n_features = 40L, n_cells = ncol(sce), n_regulators = 6L),
                created_at = Sys.time()
            ))
    }
}

ro_deconv_stub <- function(record_name = "spatial_deconv") {
    function(sce_spatial, sce_ref, ref_group_key, ...) {
        sclet:::sclet_set_analysis_state(sce_spatial, type = "spatial",
            id = record_name, method = "cell2location",
            inputs = list(ref_group_key = ref_group_key, n_common_genes = 60L),
            artifacts = list(analysis_key = "spatial_deconv",
                columns = c("abund_T", "abund_B")))
    }
}

test_that("the record-only actions are registered under their own groups", {
    sce <- ro_sce()
    expected <- c(
        run_imputation = "preprocess",
        classify_malignant_epithelium = "malignant",
        map_reference_knn = "annotation",
        map_reference_symphony = "annotation",
        run_spatial_deconvolution = "spatial",
        run_regvelo = "velocity"
    )
    all_actions <- AIDefaultExecutionRegistry(sce, include = "all")
    expect_equal(length(all_actions), 46L)
    for (action in names(expected)) {
        expect_true(action %in% names(all_actions), info = action)
        expect_true(action %in% names(AIDefaultExecutionRegistry(sce,
            include = expected[[action]])), info = action)
        descriptor <- all_actions[[action]]
        expect_identical(descriptor$returns, "sce")
        expect_true(descriptor$mutates_object)
        expect_true(descriptor$requires_confirmation)
    }
    # `preprocess` gained one action, and nothing else in it changed.
    expect_setequal(names(AIDefaultExecutionRegistry(sce, include = "preprocess")),
        c("normalize_data", "find_variable_features", "scale_data",
          "run_decontx", "run_imputation"))
})

test_that("exactly two of the six are offered to the planner", {
    sce <- ro_sce()
    catalog <- sclet:::sclet_ai_describe_registry(
        AIDefaultExecutionRegistry(sce, include = "all"))
    expect_true("run_imputation" %in% names(catalog))
    expect_true("classify_malignant_epithelium" %in% names(catalog))
    hidden <- c("map_reference_knn", "map_reference_symphony", "run_regvelo",
        "run_spatial_deconvolution")
    expect_false(any(hidden %in% names(catalog)))
    expect_true(all(hidden %in% names(AIDefaultExecutionRegistry(sce, include = "all"))))
})

test_that("every registered action is discoverable through the group hint", {
    sce <- ro_sce()
    for (action in names(AIDefaultExecutionRegistry(sce, include = "all"))) {
        expect_true(length(sclet:::sclet_ai_action_groups(sce, action)) >= 1L,
            info = action)
    }
})

test_that("each action refuses a fixture that lacks its input, naming the gate", {
    # imputation: the assay it imputes must exist
    imp <- AIDefaultExecutionRegistry(ro_sce(), include = "preprocess")[["run_imputation"]]
    expect_match(ro_gate(imp, ro_sce(), list(assay = "nope")), "assay_missing")
    expect_true(isTRUE(ro_gate(imp, ro_sce(), list())))

    # malignant: two seed vectors, a sample column, and limma
    mal <- AIDefaultExecutionRegistry(ro_sce(), include = "malignant")[["classify_malignant_epithelium"]]
    genes <- rownames(ro_sce())
    expect_match(ro_gate(mal, ro_sce(), list(seed_tumor = genes[1:3])), "seed_missing")
    expect_match(ro_gate(mal, ro_sce(), list(seed_tumor = "x", seed_normal = "y")),
        "seed_not_found")
    expect_match(ro_gate(mal, ro_sce(), list(seed_tumor = genes[1:3],
        seed_normal = genes[4:6], sample_col = "nope")), "sample_column_missing")

    # mappings: a second SCE is required
    knn <- AIDefaultExecutionRegistry(ro_sce(), include = "annotation")[["map_reference_knn"]]
    expect_match(ro_gate(knn, ro_sce(), list()), "reference_missing")
    no_label <- ro_ref()
    SummarizedExperiment::colData(no_label)$label <- NULL
    expect_match(ro_gate(knn, ro_sce(), list(ref = no_label)), "annotation_levels_missing")

    # deconvolution: reference, group key, counts both sides
    dec <- AIDefaultExecutionRegistry(ro_sce(), include = "spatial")[["run_spatial_deconvolution"]]
    expect_match(ro_gate(dec, ro_sce(), list()), "reference_missing")
    expect_match(ro_gate(dec, ro_sce(), list(sce_ref = ro_ref(),
        ref_group_key = "nope")), "reference_group_missing")

    # RegVelo: spliced/unspliced and a GRN with regulators + W
    rv <- AIDefaultExecutionRegistry(ro_sce(), include = "velocity")[["run_regvelo"]]
    expect_match(ro_gate(rv, ro_sce(), list(grn = list())), "spliced_unspliced_missing")
    with_assays <- ro_sce()
    base_assay <- SummarizedExperiment::assay(with_assays, "counts")
    SummarizedExperiment::assay(with_assays, "spliced") <- base_assay
    SummarizedExperiment::assay(with_assays, "unspliced") <- base_assay
    expect_match(ro_gate(rv, with_assays, list(grn = list())), "grn_invalid")
    # A named list is what the callee REJECTS; the gate must reject it too.
    expect_match(ro_gate(rv, with_assays,
        list(grn = list(regulators = "gene_1", W = matrix(1)))), "grn_invalid")
    prior <- matrix(1, nrow = 2, ncol = 2,
        dimnames = list(c("gene_1", "gene_2"), c("gene_3", "gene_4")))
    expect_true(isTRUE(ro_gate(rv, with_assays, list(grn = prior))))
    edges <- data.frame(target = "gene_1", regulator = "gene_3")
    expect_true(isTRUE(ro_gate(rv, with_assays, list(grn = edges))))
})

test_that("imputation runs for real, records and evidences the real assay name", {
    # REAL RunImputation, no mock: ALRA is installed, and mocking it is exactly
    # what let a wrong assay name ("ALRA" vs the real "alra_imputed") ship.
    skip_if_not_installed("ALRA")
    action <- AIDefaultExecutionRegistry(ro_sce(), include = "preprocess")[["run_imputation"]]
    # ALRA's automatic k-selection (`ALRA::choose_k`) can pick a k that does not
    # fit the matrix, which fails on any shape without an explicit `k`. The
    # action forwards extra arguments, so a real run passes one.
    sce <- ro_sce(ncells = 160, ngenes = 90)
    result <- action$handler(sce, list(method = "ALRA", k = 5L))
    expect_true("alra_imputed" %in% SummarizedExperiment::assayNames(result))
    expect_false("ALRA" %in% SummarizedExperiment::assayNames(result))
    # The declared contract must name the same assay, or the output gate would
    # accept a run that wrote nothing the description promised.
    expect_identical(action$output_schema$required_assays, "alra_imputed")
    evidence <- sclet:::sclet_ai_evidence_get_all(result)[["ev:imputation"]]
    expect_false(is.null(evidence))
    expect_equal(evidence$values$n_cells, ncol(sce))
    expect_equal(evidence$values$n_genes, nrow(sce))
    expect_identical(evidence$claim_level, "measured")
    # A re-run on the same object refreshes the node instead of erroring.
    again <- action$handler(result, list(method = "ALRA", k = 5L))
    expect_false(is.null(
        sclet:::sclet_ai_evidence_get_all(again)[["ev:imputation"]]))
})

test_that("the malignant classifier takes two seed vectors and evidencing is 'associated'", {
    action <- AIDefaultExecutionRegistry(ro_sce(), include = "malignant")[["classify_malignant_epithelium"]]
    sce <- ro_sce()
    genes <- rownames(sce)
    result <- with_mocked_bindings(
        action$handler(sce, list(seed_tumor = genes[1:5], seed_normal = genes[6:10],
            sample_col = "sample", name = "malignant")),
        ClassifyMalignantEpithelium = ro_malignant_stub(),
        .package = "sclet"
    )
    evidence <- sclet:::sclet_ai_evidence_get_all(result)[["ev:malignant_iterative_malignant"]]
    expect_false(is.null(evidence))
    expect_identical(evidence$claim_level, "associated")
    expect_equal(evidence$values$n_malignant, 12)
    expect_true(isTRUE(evidence$values$converged))
    # The node must not collide with the call_malignant_cells node prefix.
    expect_false("ev:malignant_malignant" %in% names(
        sclet:::sclet_ai_evidence_get_all(result)))
})

test_that("both mappings transfer labels and evidence their label counts", {
    sce <- ro_sce()
    ref <- ro_ref()
    knn <- AIDefaultExecutionRegistry(sce, include = "annotation")[["map_reference_knn"]]
    out <- with_mocked_bindings(
        knn$handler(sce, list(ref = ref, name = "knn_map")),
        RunKNNPredict = ro_mapping_stub("KNNPredict", "knn_map"),
        .package = "sclet"
    )
    ev <- sclet:::sclet_ai_evidence_get_all(out)[["ev:mapping_knn_knn_map"]]
    expect_equal(ev$values$n_labels, 2L)
    expect_equal(ev$values$n_cells, ncol(sce))

    sym <- AIDefaultExecutionRegistry(sce, include = "annotation")[["map_reference_symphony"]]
    out2 <- with_mocked_bindings(
        sym$handler(sce, list(ref = ref, name = "sym_map")),
        RunSymphonyMapping = ro_mapping_stub("Symphony", "sym_map"),
        .package = "sclet"
    )
    ev2 <- sclet:::sclet_ai_evidence_get_all(out2)[["ev:mapping_symphony_sym_map"]]
    # The real Symphony record has n_query_cells / n_ref_cells, not n_labels.
    expect_equal(ev2$values$n_query_cells, ncol(sce))
    expect_equal(ev2$values$n_ref_cells, ncol(ref))
    expect_equal(ev2$values$n_cells, ncol(sce))
})

test_that("deconvolution records the spatial result with its common-gene count", {
    sce <- ro_sce()
    ref <- ro_ref()
    action <- AIDefaultExecutionRegistry(sce, include = "spatial")[["run_spatial_deconvolution"]]
    out <- with_mocked_bindings(
        action$handler(sce, list(sce_ref = ref, ref_group_key = "label",
            name = "spatial_deconv")),
        RunSpatialDeconvolution = ro_deconv_stub(),
        .package = "sclet"
    )
    ev <- sclet:::sclet_ai_evidence_get_all(out)[["ev:spatial_deconv_spatial_deconv"]]
    expect_false(is.null(ev))
    expect_equal(ev$values$n_common_genes, 60L)
    expect_equal(ev$values$n_columns, 2L)
    expect_identical(ev$claim_level, "measured")
})

test_that("RegVelo records its GRN counts", {
    action <- AIDefaultExecutionRegistry(ro_sce(), include = "velocity")[["run_regvelo"]]
    out <- with_mocked_bindings(
        action$handler(ro_sce(), list(grn = matrix(1, 2, 2,
            dimnames = list(c("gene_1", "gene_2"), c("gene_3", "gene_4"))))),
        RunRegVelo = ro_regvelo_stub(),
        .package = "sclet"
    )
    ev <- sclet:::sclet_ai_evidence_get_all(out)[["ev:velocity_regvelo_regvelo"]]
    expect_equal(ev$values$n_targets, 40L)
    expect_equal(ev$values$n_edges, 120L)
    expect_identical(ev$claim_level, "measured")
})

test_that("a driver plan for a mapping round-trips into the trace graph", {
    sce <- ro_sce()
    ref <- ro_ref()
    registry <- AIDefaultExecutionRegistry(sce, include = "annotation")
    plan <- sclet:::new_sclet_ai_plan(
        task = "knn_map",
        actions = list(list(id = "s1", action = "map_reference_knn",
            params = list(ref = ref, name = "knn_map"))),
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint
    )
    validation <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_true(validation$valid)
    executed <- with_mocked_bindings(
        ExecuteAIPlan(sce, plan, registry, validation = validation,
            dry_run = FALSE, confirmation = validation$confirmation_token),
        RunKNNPredict = ro_mapping_stub("KNNPredict", "knn_map"),
        .package = "sclet"
    )
    expect_equal(executed$status, "completed")
    graph <- AnalysisGraph(executed$object)
    produced <- graph$edges$to[graph$edges$kind == "produced"]
    expect_true("record:mapping::knn_map" %in% produced)
    expect_true("evidence:ev:mapping_knn_knn_map" %in% produced)
    expect_equal(graph$meta$unresolved_edges, 0L)
})

test_that("imputation can chain into the rest of a plan", {
    sce <- ro_sce()
    registry <- AIDefaultExecutionRegistry(sce, include = "all")
    plan <- sclet:::new_sclet_ai_plan(
        task = "impute_then_hvg",
        actions = list(
            list(id = "s1", action = "run_imputation"),
            list(id = "s2", action = "find_variable_features",
                params = list(nfeatures = 50L), depends_on = "s1")
        ),
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint
    )
    validation <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_true(validation$valid)
})

test_that("a seed supplied as one flat vector is refused", {
    sce <- ro_sce()
    registry <- AIDefaultExecutionRegistry(sce, include = "malignant")
    genes <- rownames(sce)
    plan <- sclet:::new_sclet_ai_plan(
        task = "flat_seed",
        actions = list(list(id = "s1", action = "classify_malignant_epithelium",
            params = list(seed = genes[1:6], sample_col = "sample"))),
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint
    )
    validation <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_false(isTRUE(validation$valid))
    expect_match(paste(unlist(validation$errors), collapse = " | "), "seed_tumor")
})

test_that("the skill describes the newly hidden actions", {
    block <- sclet:::sclet_ai_skill_context(NULL)
    for (id in c("map_reference_knn", "map_reference_symphony", "run_regvelo",
        "run_spatial_deconvolution")) {
        expect_match(block, id, fixed = TRUE)
    }
    expect_match(block, "seed_tumor", fixed = TRUE)
})
