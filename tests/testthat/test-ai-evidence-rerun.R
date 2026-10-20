## A re-run of the same analysis must not leave a stale evidence node behind.
##
## Measured 2026-10-21: a node id is derived from the analysis name
## (`ev:doublet_doublet_finder`), so running the same action twice reproduces
## it. RecordAIEvidence() refuses the duplicate and every recorder wrapped the
## call in tryCatch(error = function(e) object), so the second run returned an
## object whose colData said 0 doublets while the citable node still said 150.
## Seven recorders shared the shape; the fix is one shared function, pinned
## here, plus a routing assertion per recorder below.

rerun_test_object <- function(n = 300L, g = 50L) {
    set.seed(1)
    counts <- matrix(stats::rpois(n * g, 3), g, n,
        dimnames = list(paste0("g", seq_len(g)), paste0("c", seq_len(n))))
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
    sce <- sclet:::sclet_set_analysis_state(sce, type = "preprocess",
        id = "doublet_finder", method = "scDblFinder", inputs = list(),
        artifacts = list(), summary = list(), active = TRUE)
    sce$scDblFinder.class <- rep(c("singlet", "doublet"), length.out = n)
    sce$scDblFinder.score <- rep(c(0.1, 0.9), length.out = n)
    sce
}

test_that("an identical re-run is a no-op, not a duplicate error", {
    sce <- rerun_test_object()
    first <- sclet:::sclet_ai_record_doublet_evidence(sce, name = "doublet_finder")
    nodes <- sclet:::sclet_ai_evidence_get_all(first)
    expect_length(nodes, 1L)
    expect_equal(nodes[[1L]]$values$n_called_doublets, 150L)
    ## Same data, recorded again: no error, one node, same numbers.
    again <- sclet:::sclet_ai_record_doublet_evidence(first, name = "doublet_finder")
    expect_length(sclet:::sclet_ai_evidence_get_all(again), 1L)
    expect_equal(sclet:::sclet_ai_evidence_get_all(again)[[1L]]$values$n_called_doublets, 150L)
})

test_that("changed data without replace refuses loudly instead of keeping the old node", {
    sce <- rerun_test_object()
    first <- sclet:::sclet_ai_record_doublet_evidence(sce, name = "doublet_finder")
    ## The honest number is now 0; the node must not keep saying 150.
    SummarizedExperiment::colData(first)$scDblFinder.class <- "singlet"
    expect_error(
        sclet:::sclet_ai_record_doublet_evidence(first, name = "doublet_finder"),
        "evidence_id_conflict"
    )
})

test_that("replace = TRUE refreshes the node to the new numbers", {
    sce <- rerun_test_object()
    first <- sclet:::sclet_ai_record_doublet_evidence(sce, name = "doublet_finder")
    SummarizedExperiment::colData(first)$scDblFinder.class <- "singlet"
    refreshed <- sclet:::sclet_ai_record_doublet_evidence(first, name = "doublet_finder",
        replace = TRUE)
    nodes <- sclet:::sclet_ai_evidence_get_all(refreshed)
    expect_length(nodes, 1L)
    expect_equal(nodes[[1L]]$values$n_called_doublets, 0L)
    expect_equal(sum(refreshed$scDblFinder.class == "doublet"), 0L)
})

test_that("the shared recorder refuses a different node, and no recorder silently drops one", {
    ## The routing itself is the fix: a recorder that still calls
    ## RecordAIEvidence() inside tryCatch would fail these.
    recorders <- c(
        "sclet_ai_record_doublet_evidence", "sclet_ai_record_cnv_evidence",
        "sclet_ai_record_milo_evidence", "sclet_ai_record_velocity_evidence",
        "sclet_ai_record_trajectory_evidence", "sclet_ai_record_fate_evidence",
        "sclet_ai_record_enrichment_evidence", "sclet_ai_record_communication_evidence"
    )
    for (name in recorders) {
        fn <- getFromNamespace(name, "sclet")
        body_text <- paste(deparse(body(fn)), collapse = " ")
        expect_true(grepl("sclet_ai_record_evidence_for_run", body_text, fixed = TRUE),
            info = paste(name, "must route through the re-run-aware recorder"))
    }
})

## -- the code-to-label mapping is recoverable, and this pins how -----------
## Every group reaching an evidence node is anonymized to cluster_N by
## paste0("cluster_", seq_along(unique(<column values>))). Measured on the
## pancreas object (2026-10-20): the trajectory reported start_group
## "cluster_7" and inputs$group = "cluster", whose seventh level is
## "Ngn3 low EP" -- the root the plan actually asked for. Nothing documented
## that round-trip, so a reader could not resolve a code without guessing.
test_that("an anonymized group code resolves back to its label positionally", {
    sce <- rerun_test_object()
    sce$cluster <- factor(rep(c("Alpha", "Beta", "Delta", "Ductal", "Epsilon",
        "Ngn3 high EP", "Ngn3 low EP", "Pre-endocrine"), length.out = ncol(sce)))
    labels <- levels(factor(as.character(sce$cluster)))
    ## The derivation the evidence layer uses.
    expect_equal(paste0("cluster_", seq_along(labels))[[7L]], "cluster_7")
    expect_equal(labels[[7L]], "Ngn3 low EP")
    ## Re-levelling changes the mapping, which is why a code must never be
    ## carried across objects.
    subsetted <- sce[, as.character(sce$cluster) != "Alpha"]
    subsetted$cluster <- droplevels(subsetted$cluster)
    expect_false(identical(levels(factor(as.character(subsetted$cluster)))[[7L]],
        "Ngn3 low EP"))
})

## -- audit 2026-10-21: the marker node and the clarification transcript -----
## Two more instances of the same two shapes, found by scanning every
## tryCatch(error = function(e) <object>) site rather than by reading the
## backlog.

test_that("re-running marker detection with other thresholds does not leave a stale node", {
    ## ev:findallmarkers is a FIXED id, so a stricter re-run collides. Before
    ## the fix the first run's n_rows/min_pct stayed in the ledger while the
    ## record held the second run's table (measured: record 0 rows, node still
    ## n_rows = 2, min_pct = 0.01).
    set.seed(1)
    n <- 120L; g <- 40L
    counts <- matrix(stats::rpois(n * g, 6), g, n,
        dimnames = list(paste0("g", seq_len(g)), paste0("c", seq_len(n))))
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
    sce$cluster <- factor(rep(c("a", "b"), each = 60))
    ActiveIdent(sce) <- "cluster"
    sce <- NormalizeData(sce)
    sce <- RunPCA(sce, ncomponents = 10)

    run_round <- function(object, task, params) {
        registry <- AIDefaultExecutionRegistry(object, include = "all")
        plan <- sclet:::new_sclet_ai_plan(task = task,
            context_fingerprint = GetAnalysisLedger(object)$fingerprint,
            actions = list(list(id = "d1", action = "run_de_test", params = params)))
        validation <- ValidateAIPlan(plan, object = object, registry = registry,
            strict = FALSE)
        expect_true(isTRUE(validation$valid))
        ExecuteAIPlan(object, plan, registry, validation = validation, dry_run = FALSE,
            confirmation = validation$confirmation_token)
    }
    first <- run_round(sce, "markers1",
        list(min.pct = 0.01, logfc.threshold = 0.1, name = "markers"))
    node1 <- sclet:::sclet_ai_evidence_get_all(first$object)[["ev:findallmarkers"]]$values
    expect_equal(node1$min_pct, 0.01)

    ## Stricter thresholds: without replace the recorder refuses, so the stale
    ## node cannot be mistaken for the new run.
    strict <- list(min.pct = 0.9, logfc.threshold = 5, name = "markers")
    refused <- run_round(first$object, "markers2", strict)
    expect_equal(refused$status, "failed")
    failed_step <- Filter(function(r) identical(r$status, "failed"), refused$results)[[1L]]
    expect_true(grepl("evidence_id_conflict", failed_step$error, fixed = TRUE))

    ## With replace_evidence = TRUE the node follows the object.
    strict$replace_evidence <- TRUE
    refreshed <- run_round(first$object, "markers3", strict)
    expect_equal(refreshed$status, "completed")
    node2 <- sclet:::sclet_ai_evidence_get_all(refreshed$object)[["ev:findallmarkers"]]$values
    expect_equal(node2$min_pct, 0.9)
    detest <- sclet:::sclet_get_state_record(refreshed$object, "detest", "markers")
    expect_equal(as.integer(node2$n_rows), as.integer(nrow(detest$artifacts$result)))
})

test_that("the clarification transcript reports the recording outcome, not the intent", {
    ## Measured 2026-10-21 by mocking the recorder to fail: the transcript said
    ## recorded = TRUE while the ledger had 0 nodes. A transcript that claims a
    ## durable record the ledger does not have is how a later round reads a
    ## decision the ledger never stored.
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1L, 5L, 12L)))
    clarification <- list(questions = list(list(id = "q1", text = "pick",
        blocked_action = "run_cci", options = c("CellChat"))))

    broken <- with_mocked_bindings(
        ResolveAIClarifications(sce, clarification, ask = function(q) "CellChat",
            record = TRUE),
        sclet_ai_record_clarification_response = function(...) stop("recorder exploded"),
        .package = "sclet"
    )
    expect_false(isTRUE(broken$transcript[[1L]]$recorded))
    expect_true(grepl("recorder exploded", broken$transcript[[1L]]$record_error,
        fixed = TRUE))
    expect_length(sclet:::sclet_ai_evidence_get_all(broken$object), 0L)

    ## The happy path still reports TRUE, and the node really exists.
    ok <- ResolveAIClarifications(sce, clarification, ask = function(q) "CellChat",
        record = TRUE)
    expect_true(isTRUE(ok$transcript[[1L]]$recorded))
    expect_true(is.null(ok$transcript[[1L]]$record_error))
    expect_length(sclet:::sclet_ai_evidence_get_all(ok$object), 1L)
})

test_that("a success criterion blames the ledger read, not the analysis, when the read failed", {
    ## Same class as the two above: a swallowed error turned into a sentence
    ## about the analysis. "no evidence node found for ev:x" sends a caller
    ## looking for a step that never ran when the real problem was the read.
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1L, 5L, 12L)))
    criterion <- list(list(id = "c1", description = "x", source = "evidence",
        evidence_id = "ev:nope", check = "exists"))
    plan <- sclet:::new_sclet_ai_plan(task = "t",
        actions = list(list(id = "a1", action = "inspect_status")),
        success_criteria = criterion)

    failed_read <- with_mocked_bindings(
        sclet:::sclet_ai_assess_success_criteria(plan, list(), sce, character()),
        sclet_ai_evidence_get_all = function(...) stop("ledger unreadable"),
        .package = "sclet"
    )
    expect_equal(failed_read[[1L]]$status, "not_available")
    expect_true(grepl("could not read the evidence ledger", failed_read[[1L]]$reason,
        fixed = TRUE))
    expect_false(grepl("no evidence node found", failed_read[[1L]]$reason, fixed = TRUE))

    ## The genuinely-absent case keeps its own wording.
    absent <- sclet:::sclet_ai_assess_success_criteria(plan, list(), sce, character())
    expect_true(grepl("no evidence node found", absent[[1L]]$reason, fixed = TRUE))
})
