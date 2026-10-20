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
