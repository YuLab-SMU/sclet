## The report layer (M3): what a report may say, and what it may not.
##
## Three defects were fixed here and one layer was built:
##
##   1. `analysis_story$conflicts` fired from the second round of any analysis,
##      because the append-only AI log types (ai_execution / ai_evidence)
##      were treated like competing computational states.
##   2. The marker evidence node carried fields that were never measured --
##      `top_up_ranks` was 1:10 by construction, `top_up_fractions` was ten
##      NAs, two opaque integer codes had no consumer -- plus an `n_significant`
##      that restated the table's own pre-filter (12,559 of 12,572).
##   3. `RunDEtest()` did not record its own cutoffs, so `get_detest()` alone
##      could not say what produced the table.
##
## And `RenderAIReport()` renders a recorded result from the ledger and
## re-checks it: unresolved citations and claim-level drift are printed, not
## smoothed over.

render_report_object <- function(n_cells = 40L) {
    set.seed(3L)
    counts <- matrix(stats::rpois(20L * n_cells, lambda = 5L), nrow = 20L,
        ncol = n_cells, dimnames = list(paste0("gene_", seq_len(20L)),
            paste0("c", seq_len(n_cells))))
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, nfeatures = 10L)
    sce <- RunPCA(sce, ncomponents = 5L)
    sce <- FindNeighbors(sce, dims = seq_len(5L), reduction = "PCA")
    sce <- FindClusters(sce, resolution = 0.8, seed = 1L)
    ## Two identity groups are planted so the marker test has something to
    ## compare: a 20-gene random matrix can collapse to one cluster, and the
    ## point of these tests is the RECORD, not the clustering.
    SingleCellExperiment::colLabels(sce) <- factor(rep(c("g1", "g2"),
        length.out = ncol(sce)))
    ActiveIdent(sce) <- "colLabels"
    sce
}

render_report_evidence <- function(sce) {
    RecordAIEvidence(sce, list(
        id = "ev:render_test",
        kind = "deterministic_summary",
        values = list(n_rows = 12L, raw_values_included = FALSE),
        claim_level = "associated"
    ), source = NULL, parents = character(), scope = NULL)
}

## -- 1: append-only logs are not conflicts ---------------------------------
test_that("multiple executions are not reported as a conflict", {
    sce <- render_report_object()
    for (round in c("r1", "r2")) {
        sce <- sclet_set_analysis(sce, paste0("ai_execution_", round), list(
            id = paste0("ai_execution_", round), type = "ai_execution",
            method = "ExecuteAIPlan", summary = list(status = "completed")))
    }
    story <- GetAnalysisLedger(sce)$analysis_story
    exec_conflicts <- Filter(function(x) identical(x$type, "ai_execution"),
        story$conflicts)
    expect_length(exec_conflicts, 0L)

    ## The computational state types still conflict: two detest records mean
    ## the reader cannot tell which one is in force.
    sce <- sclet_set_analysis_state(sce, "detest", "one", "FindAllMarkers",
        inputs = list(), summary = list(n_markers = 1L), active = FALSE)
    sce <- sclet_set_analysis_state(sce, "detest", "two", "FindAllMarkers",
        inputs = list(), summary = list(n_markers = 2L), active = TRUE)
    story <- GetAnalysisLedger(sce)$analysis_story
    detest_conflicts <- Filter(function(x) identical(x$type, "detest"),
        story$conflicts)
    expect_gte(length(detest_conflicts), 1L)
})

## -- 2: the marker evidence says only what it measured ---------------------
test_that("marker evidence carries no unmeasured or restated fields", {
    sce <- render_report_object()
    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "annotation"))
    plan <- sclet:::new_sclet_ai_plan(
        task = "markers",
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
        actions = list(list(id = "m1", action = "run_de_test",
            params = list(name = "markers", min.pct = 0.2, logfc.threshold = 0.3)))
    )
    validation <- ValidateAIPlan(plan, object = sce, registry = registry, strict = TRUE)
    executed <- ExecuteAIPlan(sce, plan, registry, validation = validation,
        dry_run = FALSE, confirmation = validation$confirmation_token)
    node <- sclet:::sclet_ai_evidence_get_all(executed$object)[["ev:findallmarkers"]]
    values <- node$values

    ## Never measured -> gone.
    expect_null(values$top_up_ranks)
    expect_null(values$top_up_fractions)
    expect_null(values$annotation_method_code)
    expect_null(values$annotation_scope_code)
    ## The table is already filtered by the caller's cutoffs, so a second
    ## significance count restates the filter.
    expect_null(values$n_significant)

    ## What remains is measured, and the thresholds are the run's own.
    expect_true(values$n_rows >= 0L)
    expect_true(values$n_upregulated >= 0L)
    expect_equal(values$min_pct, 0.2)
    expect_equal(values$logfc_threshold, 0.3)
})

## -- 3: the detest record carries its own cutoffs --------------------------
test_that("a detest record records the cutoffs that produced it", {
    sce <- render_report_object()
    sce <- RunDEtest(sce, name = "with_cutoffs", min.pct = 0.25, logfc.threshold = 0.4)
    record <- get_detest(sce, "with_cutoffs")
    expect_equal(record$params$min.pct, 0.25)
    expect_equal(record$params$logfc.threshold, 0.4)

    ## The same for a pairwise test that relies on the defaults: the record
    ## must still say what those defaults were.
    ids <- levels(SingleCellExperiment::colLabels(sce))
    if (length(ids) >= 2L) {
        sce <- RunDEtest(sce, ident.1 = ids[[1L]], ident.2 = ids[[2L]],
            name = "pairwise")
        pair <- get_detest(sce, "pairwise")
        expect_equal(pair$params$min.pct, 0.01)
        expect_equal(pair$params$logfc.threshold, 0.1)
    }
})

## -- 4: RenderAIReport renders from the ledger and re-checks it ------------
test_that("RenderAIReport renders a recorded result and passes its own checks", {
    sce <- render_report_object()
    sce <- render_report_evidence(sce)
    result <- sclet:::new_sclet_ai_result(
        task = "render_report",
        answer = "One finding.",
        findings = list(list(
            statement = "markers were tested",
            claim_level = "associated",
            evidence_refs = list("ev:render_test")
        )),
        warnings = list("a caveat the reader needs"),
        recommendations = list("do the next thing")
    )
    sce <- RecordAIResult(sce, result, id = "ai_report_render_test")

    text <- RenderAIReport(sce)
    expect_type(text, "character")
    expect_match(text, "ai_report_render_test", fixed = TRUE)
    expect_match(text, "markers were tested", fixed = TRUE)
    expect_match(text, "ev:render_test", fixed = TRUE)
    expect_match(text, "associated", fixed = TRUE)
    expect_match(text, "## Verification", fixed = TRUE)
    expect_match(text, "all checks passed", fixed = TRUE)
    ## The warning the result carried survives into the ledger and the report.
    expect_match(text, "a caveat the reader needs", fixed = TRUE)
    expect_match(text, "do the next thing", fixed = TRUE)

    path <- tempfile(fileext = ".md")
    on.exit(unlink(path), add = TRUE)
    written <- RenderAIReport(sce, id = "ai_report_render_test", path = path)
    expect_true(file.exists(path))
    expect_equal(readLines(path, warn = FALSE), strsplit(written, "\n", fixed = TRUE)[[1L]])

    expect_error(RenderAIReport(sce, id = "does_not_exist"), "no recorded AI result")
})

test_that("RenderAIReport reports a claim level that no longer re-derives", {
    sce <- render_report_object()
    sce <- render_report_evidence(sce)
    result <- sclet:::new_sclet_ai_result(
        task = "render_drift",
        findings = list(list(
            statement = "a claim",
            claim_level = "associated",
            evidence_refs = list("ev:render_test")
        ))
    )
    recorded <- RecordAIResult(sce, result, id = "ai_report_drift")
    ## Tamper with the RECORDED finding, not the evidence: the report must
    ## notice that the stored level is no longer what the engine derives.
    ledger <- GetAnalysisLedger(recorded, detail = "full", include_artifacts = TRUE)
    expect_true("ai_report_drift" %in% names(ledger$analyses))
    record <- sclet:::sclet_get_analysis(recorded, "ai_report_drift")
    record$artifacts$findings[[1L]]$claim_level <- "consistent_with"
    tampered <- sclet_set_analysis(recorded, "ai_report_drift", record)

    text <- RenderAIReport(tampered, id = "ai_report_drift")
    expect_match(text, "Claim levels re-derived: FAIL", fixed = TRUE)
    expect_match(text, "problems found", fixed = TRUE)
})

test_that("RenderAIReport lists executed steps, so a read-only step is traceable", {
    sce <- render_report_object()
    registry <- AIDefaultExecutionRegistry(sce, include = "read")
    plan <- sclet:::new_sclet_ai_plan(
        task = "qc_step",
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
        actions = list(list(id = "q1", action = "check_qc"))
    )
    validation <- ValidateAIPlan(plan, object = sce, registry = registry, strict = TRUE)
    executed <- ExecuteAIPlan(sce, plan, registry, validation = validation,
        dry_run = FALSE, confirmation = validation$confirmation_token)
    result <- sclet:::new_sclet_ai_result(task = "after_qc", answer = "ok")
    executed$object <- RecordAIResult(executed$object, result, id = "ai_report_after_qc")

    text <- RenderAIReport(executed$object, id = "ai_report_after_qc")
    expect_match(text, "## Executed steps", fixed = TRUE)
    expect_match(text, "check_qc", fixed = TRUE)
    expect_match(text, "q1", fixed = TRUE)
})

## -- 5: an enrichment run off a DERIVED marker table keeps its parent ------
## Measured on the M1 report (2026-10-20): the KEGG enrichment ran on the
## id-converted copy of the marker table, whose record carries no evidence node
## of its own. The enrichment node therefore had no parents and the ceiling
## counted markers + GO + KEGG as TWO independent lines of support when they are
## one measurement annotated twice.
test_that("enrichment from a derived detest record keeps its parent chain", {
    ## Two groups with planted signal, so the marker table is non-empty at the
    ## default cutoffs (a random matrix with two groups can return zero rows,
    ## because FindAllMarkers drops pval > return.thresh before returning).
    set.seed(9L)
    n_cells <- 60L
    counts <- matrix(stats::rpois(20L * n_cells, lambda = 5L), nrow = 20L,
        ncol = n_cells, dimnames = list(paste0("gene_", seq_len(20L)),
            paste0("c", seq_len(n_cells))))
    groups <- rep(c("g1", "g2"), each = n_cells / 2L)
    counts[seq_len(5L), groups == "g1"] <- counts[seq_len(5L), groups == "g1"] * 30
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
    sce <- NormalizeData(sce)
    SingleCellExperiment::colLabels(sce) <- factor(groups)
    ActiveIdent(sce) <- "colLabels"

    registry <- AIDefaultExecutionRegistry(sce, include = c("read", "annotation"))
    plan <- sclet:::new_sclet_ai_plan(
        task = "markers",
        context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
        actions = list(list(id = "m1", action = "run_de_test",
            params = list(name = "markers_by_cluster")))
    )
    validation <- ValidateAIPlan(plan, object = sce, registry = registry, strict = TRUE)
    executed <- ExecuteAIPlan(sce, plan, registry, validation = validation,
        dry_run = FALSE, confirmation = validation$confirmation_token)
    sce <- executed$object

    markers <- get_detest(sce, "markers_by_cluster")$artifacts$result
    expect_gt(nrow(markers), 0L)
    derived <- markers
    derived$gene <- paste0("entrez_", seq_len(nrow(derived)))
    sce <- sclet_set_analysis_state(
        object = sce, type = "detest", id = "markers_entrez",
        method = "test derived copy", inputs = list(parent_detest = "markers_by_cluster"),
        artifacts = list(result = derived), params = list(),
        summary = list(n_markers = nrow(derived)), active = FALSE
    )

    parents <- sclet:::sclet_ai_enrichment_parent_nodes(sce, "markers_entrez")
    expect_equal(parents, "ev:findallmarkers")
})
