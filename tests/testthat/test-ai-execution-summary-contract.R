# Tests for execution summary and referenceable output fields agreement
# These tests verify that sclet_ai_execution_summary() and
# sclet_ai_referenceable_output_fields() agree by construction.

test_that("the three reproduction criteria grade met after the fix", {
    skip_if_not_installed("scuttle")
    skip_if_not_installed("scran")
    skip_if_not_installed("igraph")

    set.seed(42)
    counts <- matrix(rpois(50 * 40, lambda = 5), nrow = 50, ncol = 40)
    rownames(counts) <- paste0("gene", seq_len(nrow(counts)))
    colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))

    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "preprocess", "dimred", "graph", "cluster", "annotation"))
    plan <- sclet:::new_sclet_ai_plan(
        task = "reproduction_plan",
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
        actions = list(
            list(id = "s1", action = "normalize_data", params = list(scale.factor = 10000)),
            list(id = "s2", action = "find_variable_features", params = list(nfeatures = 15), depends_on = "s1"),
            list(id = "s3", action = "scale_data", params = list(features = rownames(sce)), depends_on = "s2"),
            list(id = "s4", action = "run_pca", params = list(ncomponents = 5), depends_on = "s3"),
            list(id = "s5", action = "find_neighbors", params = list(dims = 1:5, k = 5), depends_on = "s4"),
            list(id = "s6", action = "find_clusters", params = list(resolution = 1.0), depends_on = "s5"),
            list(id = "s7", action = "run_de_test", params = list(name = "findallmarkers"), depends_on = "s6")
        ),
        success_criteria = list(
            list(
                id = "c_graph",
                description = "active graph is set after find_neighbors",
                source = "action_output",
                step_id = "s5",
                check = "exists",
                field = "active_graph"
            ),
            list(
                id = "c_reqstates",
                description = "required_states is declared in run_de_test output",
                source = "action_output",
                step_id = "s7",
                check = "exists",
                field = "required_states"
            ),
            list(
                id = "c_ident",
                description = "active_ident is set after find_clusters",
                source = "action_output",
                step_id = "s6",
                check = "exists",
                field = "active_ident"
            )
        )
    )

    v <- ValidateAIPlan(plan, object = sce, registry = registry)
    expect_true(v$valid, paste(v$errors, collapse = "; "))

    exec <- ExecuteAIPlan(
        sce, plan, registry,
        validation = v,
        dry_run = FALSE,
        confirmation = v$confirmation_token,
        record = FALSE
    )
    expect_equal(exec$status, "completed")
    expect_length(exec$success_assessment, 3L)

    for (crit in exec$success_assessment) {
        expect_equal(crit$status, "met", info = paste(crit$id, ":", crit$reason))
    }
})

test_that("find_neighbors summary carries all declared keys", {
    skip_if_not_installed("scuttle")
    skip_if_not_installed("scran")
    skip_if_not_installed("igraph")

    set.seed(1)
    sce <- scuttle::mockSCE(ncells = 40, ngenes = 120)
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 30, method = "scran")
    sce <- RunPCA(sce, subset_row = VariableFeatures(sce), ncomponents = 10)
    sce <- FindNeighbors(sce, dims = 1:5)

    registry <- AIDefaultExecutionRegistry(sce, include = "graph")
    descriptor <- registry$find_neighbors
    summary <- sclet:::sclet_ai_execution_summary(sce, descriptor)

    declared_keys <- names(descriptor$output_schema)
    declared_keys <- setdiff(declared_keys, "note")
    for (key in declared_keys) {
        expect_true(key %in% names(summary), info = paste("declared key", key, "missing from summary"))
    }
})

test_that("run_de_test summary carries declared required_states", {
    skip_if_not_installed("scuttle")
    skip_if_not_installed("scran")
    skip_if_not_installed("igraph")

    set.seed(2)
    sce <- scuttle::mockSCE(ncells = 40, ngenes = 120)
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 30, method = "scran")
    sce <- RunPCA(sce, subset_row = VariableFeatures(sce), ncomponents = 10)
    sce <- FindNeighbors(sce, dims = 1:5)
    sce <- FindClusters(sce)

    registry <- AIDefaultExecutionRegistry(sce, include = "all")
    descriptor <- registry$run_de_test
    summary <- sclet:::sclet_ai_execution_summary(sce, descriptor)

    expect_true("required_states" %in% names(summary))
})

test_that("note is absent from every registered action's summary", {
    skip_if_not_installed("scuttle")

    sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
    registry <- AIDefaultExecutionRegistry(sce, include = "all")

    for (action_name in names(registry)) {
        descriptor <- registry[[action_name]]
        summary <- sclet:::sclet_ai_execution_summary(sce, descriptor)
        expect_false("note" %in% names(summary), info = paste("note found in", action_name, "summary"))
    }
})

test_that("observed values beat declared values in summary", {
    skip_if_not_installed("scuttle")
    skip_if_not_installed("scran")
    skip_if_not_installed("igraph")

    set.seed(3)
    sce <- scuttle::mockSCE(ncells = 40, ngenes = 120)
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 30, method = "scran")
    sce <- RunPCA(sce, subset_row = VariableFeatures(sce), ncomponents = 10)
    sce <- FindNeighbors(sce, dims = 1:5)

    registry <- AIDefaultExecutionRegistry(sce, include = "graph")
    descriptor <- registry$find_neighbors

    descriptor_modified <- descriptor
    descriptor_modified$output_schema$active_graph <- "different_graph_name"

    summary <- sclet:::sclet_ai_execution_summary(sce, descriptor_modified)

    expect_equal(summary$active_graph, "knn_graph")
    expect_false(identical(summary$active_graph, "different_graph_name"))
})

test_that("privacy sanitizer handles broadened summary keys without new redactions", {
    skip_if_not_installed("scuttle")
    skip_if_not_installed("scran")
    skip_if_not_installed("igraph")

    set.seed(5)
    sce <- scuttle::mockSCE(ncells = 40, ngenes = 120)
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 30, method = "scran")
    sce <- RunPCA(sce, subset_row = VariableFeatures(sce), ncomponents = 10)
    sce <- FindNeighbors(sce, dims = 1:5)
    sce <- FindClusters(sce)

    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "preprocess", "dimred", "graph", "cluster"))
    descriptor <- registry$find_clusters
    summary <- sclet:::sclet_ai_execution_summary(sce, descriptor)

    context <- list(
        state = list(
            output_summary = summary
        )
    )

    sanitized <- sclet_ai_sanitize_context(sce, context, privacy = "standard")

    expect_false("active_ident" %in% sanitized$policy$redacted_fields)
    expect_false("active_graph" %in% sanitized$policy$redacted_fields)
    expect_false("required_states" %in% sanitized$policy$redacted_fields)

    expect_length(sanitized$policy$warnings, 0L)
})

test_that("referenceable fields and summary keys agree for all registered actions", {
    skip_if_not_installed("scuttle")
    skip_if_not_installed("scran")
    skip_if_not_installed("igraph")

    set.seed(4)
    sce <- scuttle::mockSCE(ncells = 40, ngenes = 120)
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 30, method = "scran")
    sce <- RunPCA(sce, subset_row = VariableFeatures(sce), ncomponents = 10)
    sce <- FindNeighbors(sce, dims = 1:5)
    sce <- FindClusters(sce)

    registry <- AIDefaultExecutionRegistry(sce, include = "all")

    for (action_name in names(registry)) {
        descriptor <- registry[[action_name]]
        referenceable <- sclet:::sclet_ai_referenceable_output_fields(descriptor)
        summary <- sclet:::sclet_ai_execution_summary(sce, descriptor)

        if (!identical(descriptor$returns, "sce")) {
            expect_identical(referenceable, character(), info = action_name)
            next
        }

        summary_keys <- names(summary)

        # `active_ident` and `active_graph` are the one deliberate exception,
        # and only in one direction: like active_assay/active_reduction they are
        # advertised unconditionally (a reference may resolve once the state
        # exists) but present only when the object actually has that state. That
        # exception is asserted explicitly after the loop rather than skipped
        # silently here, so a regression in those two fields cannot hide.
        for (field in setdiff(referenceable, c("active_ident", "active_graph"))) {
            expect_true(field %in% summary_keys, info = paste(action_name, ":", field, "in referenceable but not in summary"))
        }

        for (key in summary_keys) {
            expect_true(key %in% referenceable, info = paste(action_name, ":", key, "in summary but not in referenceable"))
        }
    }

    # This object has been clustered and has a graph, so both observed pointers
    # must genuinely be present -- not merely advertised.
    prepared <- sclet:::sclet_ai_execution_summary(sce, registry[["find_clusters"]])
    expect_true(all(c("active_ident", "active_graph") %in% names(prepared)))
    # ...and they must be observations, not the contract's declarations.
    expect_equal(prepared$active_ident, sclet:::sclet_get_active_ident(sce))
    expect_equal(prepared$active_graph, DefaultGraph(sce))
})

test_that("nested contract values are flattened to their names", {
    # A contract value such as list(clustering = "louvain_clusters") carries an
    # internal state-type label whose KEY matches the privacy classifier's group
    # pattern. Carried verbatim it made the sanitizer redact
    # `required_states.clustering` as "group labels do not meet the aggregation
    # threshold" -- a false positive on a state name, but a real warning on a
    # path that leads to the ledger. The key must survive (criteria name it) and
    # the nested key must not.
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(rpois(20, 5), nrow = 5, ncol = 4))
    )
    registry <- AIDefaultExecutionRegistry(sce, include = "all")
    summary <- sclet:::sclet_ai_execution_summary(sce, registry[["find_clusters"]])

    expect_true("required_states" %in% names(summary))
    expect_false(is.list(summary$required_states))
    expect_true(all(c("clustering") %in% summary$required_states))

    # The declared `note` prose must never be carried, for any action.
    for (action_name in names(registry)) {
        s <- sclet:::sclet_ai_execution_summary(sce, registry[[action_name]])
        if (!is.list(s)) next
        expect_false("note" %in% names(s), info = action_name)
    }
})

test_that("the exposable-contract-keys definition drives both functions", {
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = matrix(rpois(20, 5), nrow = 5, ncol = 4))
    )
    registry <- AIDefaultExecutionRegistry(sce, include = "all")
    for (action_name in names(registry)) {
        descriptor <- registry[[action_name]]
        keys <- sclet:::sclet_ai_exposable_contract_keys(descriptor)
        expect_false("note" %in% keys, info = action_name)
        if (!identical(descriptor$returns, "sce")) next
        # Every exposable key is advertised by the validator...
        expect_true(all(keys %in% sclet:::sclet_ai_referenceable_output_fields(descriptor)),
            info = action_name)
    }
})
