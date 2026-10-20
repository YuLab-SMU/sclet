## Four defects found by driving sclet with its own skills (inst/skills/*) as
## an external agent against real data -- the cached PBMC 10k fixture, 33,538
## genes x 11,769 cells -- on 2026-10-19. Each test is pinned to the artifact
## the defect corrupted, not to a printed message, because "the log said it
## worked" is exactly the claim these defects falsified.
##
##   1. The marker evidence node recorded hardcoded thresholds (0.05 / 0.1)
##      instead of the ones the run used, so a plan asking logfc.threshold =
##      0.25 produced a node describing a different analysis.
##   2. `FindNeighbors()` (via `set_dimred()`) and `RunUMAP()` left their
##      scratch `.dimred` reduction on the object forever: reducedDimNames()
##      after a plain pipeline was c("PCA", ".dimred"), the stray entry a
##      11,769 x 30 duplicate of the PCA subset.
##   3. `run_enrichment` accepted a plan with no `db` and silently applied the
##      engine default (GO), although the annotation database is the user's
##      modelling choice -- every sibling action with such a choice refuses.
##   4. The ledger did not expose the HVG selection (n / method) at all, and
##      the execution record kept no per-step parameters or timing, so a
##      reader could not reproduce or even describe the run from the ledger.

## A clustered object with real counts and a logcounts assay: enough for an
## all-cluster marker test, a graph rebuild and a UMAP.
regression_clustered_object <- function(n_cells = 60L) {
    set.seed(11L)
    counts <- matrix(stats::rpois(50L * n_cells, lambda = 5L),
        nrow = 50L, ncol = n_cells,
        dimnames = list(paste0("gene_", seq_len(50L)), paste0("c", seq_len(n_cells))))
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 30L)
    sce <- RunPCA(sce, ncomponents = 10L)
    sce <- FindNeighbors(sce, dims = seq_len(5L), reduction = "PCA")
    sce <- FindClusters(sce, resolution = 1.2)
    sce
}

regression_plan <- function(object, actions, task = "regression") {
    sclet:::new_sclet_ai_plan(
        task = task,
        context_fingerprint = GetAnalysisLedger(object)$fingerprint,
        actions = actions
    )
}

regression_execute <- function(object, plan, registry) {
    validation <- ValidateAIPlan(plan, object = object, registry = registry, strict = TRUE)
    expect_true(validation$valid)
    ExecuteAIPlan(object, plan, registry, validation = validation,
        dry_run = FALSE, confirmation = validation$confirmation_token)
}

## -- 1: marker evidence records the thresholds the run was given ------------
test_that("marker evidence records the run's own thresholds, not constants", {
    sce <- regression_clustered_object()
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "annotation"))
    plan <- regression_plan(sce, list(list(
        id = "m1", action = "run_de_test",
        params = list(name = "tuned_markers", min.pct = 0.2, logfc.threshold = 0.5)
    )))
    executed <- regression_execute(sce, plan, registry)
    expect_equal(executed$status, "completed")

    nodes <- sclet:::sclet_ai_evidence_get_all(executed$object)
    node <- nodes[["ev:findallmarkers"]]
    expect_false(is.null(node))
    expect_equal(node$values$logfc_threshold, 0.5)
    expect_equal(node$values$min_pct, 0.2)

    ## The default path records the defaults, so neither branch is a constant.
    default_sce <- regression_clustered_object()
    default_plan <- regression_plan(default_sce, list(list(
        id = "m2", action = "run_de_test", params = list(name = "default_markers")
    )))
    executed_default <- regression_execute(default_sce, default_plan, registry)
    default_nodes <- sclet:::sclet_ai_evidence_get_all(executed_default$object)
    default_node <- default_nodes[["ev:findallmarkers"]]
    expect_equal(default_node$values$logfc_threshold, 0.1)
    expect_equal(default_node$values$min_pct, 0.01)
})

## -- 2: no scratch reduction survives the graph or the embedding -----------
test_that("FindNeighbors and RunUMAP leave no .dimred on the object", {
    sce <- regression_clustered_object()
    expect_false(".dimred" %in% SingleCellExperiment::reducedDimNames(sce))

    ## The graph must still be built from the requested dims: compare it with
    ## the subset handed straight to bluster.
    manual <- bluster::makeSNNGraph(
        SingleCellExperiment::reducedDim(sce, "PCA")[, seq_len(5L), drop = FALSE],
        k = 10L, type = "rank"
    )
    graph <- sclet:::sclet_get_graph(sce, "knn_graph")
    expect_equal(igraph::ecount(graph), igraph::ecount(manual))
    expect_equal(igraph::vcount(graph), igraph::vcount(manual))

    umap <- RunUMAP(sce, dims = seq_len(5L), reduction = "PCA")
    expect_true("UMAP" %in% SingleCellExperiment::reducedDimNames(umap))
    expect_false(".dimred" %in% SingleCellExperiment::reducedDimNames(umap))
})

## -- 3: enrichment will not choose the database for the user ---------------
regression_enrichment_object <- function() {
    genes <- paste0("gene_", seq_len(20L))
    counts <- matrix(1L, nrow = length(genes), ncol = 6L,
        dimnames = list(genes, paste0("cell_", seq_len(6L))))
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
    markers <- data.frame(
        gene = genes,
        pval = c(rep(1e-9, 8L), rep(0.9, 12L)),
        avg_log2FC = c(rep(3.5, 8L), rep(-0.4, 12L)),
        cluster = rep(c("cluster_1", "cluster_2"), c(8L, 12L)),
        row.names = genes, stringsAsFactors = FALSE
    )
    sclet:::sclet_set_analysis_state(
        object = sce, type = "detest", id = "de_ai",
        method = "sclet::RunDEtest", inputs = list(),
        artifacts = list(result = markers),
        params = list(), summary = list(n_markers = nrow(markers)),
        active = TRUE
    )
}

test_that("run_enrichment refuses to pick db for an AI plan", {
    sce <- regression_enrichment_object()
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "enrichment"))

    missing_db <- regression_plan(sce, list(list(
        id = "en1", action = "run_enrichment", params = list(detest_id = "de_ai")
    )))
    validation <- ValidateAIPlan(missing_db, object = sce, registry = registry, strict = FALSE)
    expect_false(isTRUE(validation$valid))
    expect_true(any(grepl("db_missing", validation$errors, fixed = TRUE)))

    bad_db <- regression_plan(sce, list(list(
        id = "en2", action = "run_enrichment",
        params = list(detest_id = "de_ai", db = "MSigDB")
    )))
    bad_validation <- ValidateAIPlan(bad_db, object = sce, registry = registry, strict = FALSE)
    expect_false(isTRUE(bad_validation$valid))
    expect_true(any(grepl("db_invalid", bad_validation$errors, fixed = TRUE)))

    ## An explicit choice still validates; no clusterProfiler run is needed.
    good_db <- regression_plan(sce, list(list(
        id = "en3", action = "run_enrichment",
        params = list(detest_id = "de_ai", db = "GO")
    )))
    good_validation <- ValidateAIPlan(good_db, object = sce, registry = registry, strict = FALSE)
    expect_true(good_validation$valid)
})

## -- 4: the ledger carries the HVG choice and the per-step parameters ------
test_that("the ledger exposes the HVG selection the object carries", {
    sce <- regression_clustered_object()
    ledger <- GetAnalysisLedger(sce)
    expect_true("features" %in% names(ledger))
    expect_equal(ledger$features$hvg$n, 30L)
    expect_true(nzchar(as.character(ledger$features$hvg$method)))

    ## The sanitizer classifies it as an aggregate, not as an unknown field:
    ## an "unknown" verdict would silently drop it from an outbound payload.
    expect_equal(sclet:::sclet_ai_payload_field_kind("features", "features"),
        "aggregate")
})

test_that("the execution record keeps per-step parameters and timing", {
    sce <- regression_clustered_object()
    registry <- AIDefaultExecutionRegistry(sce, include = "all")
    plan <- regression_plan(sce, list(
        list(id = "n1", action = "find_neighbors",
            params = list(dims = 1:10, reduction = "PCA")),
        list(id = "c1", action = "find_clusters", params = list(resolution = 0.5),
            depends_on = "n1")
    ))
    executed <- regression_execute(sce, plan, registry)
    expect_equal(executed$status, "completed")

    ## The returned results carry the bounded parameters too, not only the id.
    expect_equal(executed$results[[1L]]$params$dims,
        list(type = "integer", length = 10L))
    expect_equal(executed$results[[1L]]$params$reduction, "PCA")

    ## The ledger record, in full detail, keeps the per-step rows.
    full <- GetAnalysisLedger(executed$object, detail = "full")
    execution_keys <- grep("^ai_execution_", names(full$analyses), value = TRUE)
    expect_gte(length(execution_keys), 1L)
    record <- full$analyses[[execution_keys[[length(execution_keys)]]]]
    steps <- record$params$steps
    expect_length(steps, 2L)
    expect_equal(steps[[1L]]$action, "find_neighbors")
    expect_equal(steps[[2L]]$params$resolution, 0.5)
    expect_true(is.numeric(steps[[1L]]$duration_sec))
    expect_gte(steps[[1L]]$duration_sec, 0)

    ## A long vector is bounded, so the record cannot grow with the object.
    expect_false(is.null(steps[[1L]]$params$dims$length))
    expect_null(steps[[1L]]$params$dims$preview)
})

## -- 5: Louvain clustering reproduces, and records the seed it used --------
## Two runs of the same PBMC 10k plan returned 13 and 12 clusters (2026-10-19):
## `igraph::cluster_louvain()` consumes random numbers, `find_clusters` had no
## seed parameter, and the record therefore described a result it could not
## reproduce. `sclet-provenance` promises the opposite ("a reader must be able
## to re-run this plan and get the same clusters").
test_that("find_clusters reproduces its labels and records the seed", {
    sce <- regression_clustered_object()
    first <- SingleCellExperiment::colLabels(FindClusters(sce, resolution = 0.8))
    second <- SingleCellExperiment::colLabels(FindClusters(sce, resolution = 0.8))
    expect_equal(as.character(first), as.character(second))

    tuned <- FindClusters(sce, resolution = 0.8, seed = 3L)
    record <- sclet:::sclet_get_state_record(tuned, "clustering", "louvain_clusters")
    expect_equal(record$params$seed, 3L)
    expect_equal(record$params$resolution, 0.8)

    ## And the seed travels through a plan's parameters into that record.
    registry <- AIDefaultExecutionRegistry(sce, include = "all")
    plan <- regression_plan(sce, list(list(id = "c1", action = "find_clusters",
        params = list(resolution = 0.4, seed = 9L))))
    executed <- regression_execute(sce, plan, registry)
    plan_record <- sclet:::sclet_get_state_record(executed$object, "clustering",
        "louvain_clusters")
    expect_equal(plan_record$params$seed, 9L)
    expect_equal(plan_record$params$resolution, 0.4)
})

test_that("clustering does not move the caller's random stream", {
    sce <- regression_clustered_object()
    set.seed(7L)
    before <- stats::runif(3L)
    set.seed(7L)
    invisible(FindClusters(sce, resolution = 0.8))
    after <- stats::runif(3L)
    expect_equal(before, after)
})
