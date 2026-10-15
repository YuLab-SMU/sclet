## Enrichment contract (batch A of the single-cell re-scope): readiness, the
## run_enrichment action, bounded evidence and the DE lineage that keeps a
## claim ceiling from counting one measurement twice.
##
## Two properties are the point of this file, and both are checked in a way
## that can actually fail:
##
##   1. No gene symbol and no enrichment category (name or accession) may leave
##      these functions into evidence. The check is run against the terms and
##      genes a REAL clusterProfiler run actually returned, so it tracks
##      whatever came back rather than a hardcoded guess list.
##   2. Enrichment is derived FROM the DE marker table, so DE evidence plus
##      enrichment evidence must still count as ONE independent line of
##      support. A control with two genuinely unrelated sources asserts 2,
##      which is what stops the lineage assertion from passing vacuously.
##
## The compute is real: clusterProfiler and org.Hs.eg.db are available in this
## environment, so no Python/basilisk backend and no LLM is ever contacted.

## ---------------------------------------------------------------------------
## Fixtures
## ---------------------------------------------------------------------------

## Real human symbols drawn from two genuine biology buckets: a hemoglobin /
## erythrocyte set that scores strongly for oxygen-transport GO terms, and a
## proliferation set. The first cluster is the one that will enrich.
sclet_ai_test_enrichment_genes <- function() {
    c("HBB", "HBA1", "HBD", "HBE1", "HBG1", "HBG2", "AHSP", "ALAS2", "SLC4A1",
        "GYPA", "TP53", "ACTB", "GAPDH", "CDK1", "CCNB1", "MCM2", "PCNA",
        "MKI67", "BRCA1", "RAD51")
}

sclet_ai_test_enrichment_markers <- function(with_cluster = TRUE) {
    genes <- sclet_ai_test_enrichment_genes()
    markers <- data.frame(
        gene = genes,
        pval = c(rep(1e-12, 8L), rep(0.9, length(genes) - 8L)),
        avg_log2FC = c(rep(3.5, 8L), rep(-0.4, length(genes) - 8L)),
        row.names = genes,
        stringsAsFactors = FALSE
    )
    if (isTRUE(with_cluster)) {
        markers$cluster <- c(rep("cluster_1", 8L), rep("cluster_2", nrow(markers) - 8L))
    }
    markers
}

## An object carrying a completed detest record whose artifacts$result is the
## marker table, which is exactly what RunEnrichment() consumes.
sclet_ai_test_enrichment_object <- function(detest_id = "de_ai", with_cluster = TRUE,
                                          markers = sclet_ai_test_enrichment_markers(with_cluster)) {
    genes <- sclet_ai_test_enrichment_genes()
    counts <- matrix(1L, nrow = length(genes), ncol = 6L,
        dimnames = list(genes, paste0("cell_", seq_len(6L))))
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
    sce$cluster <- factor(rep(c("cluster_1", "cluster_2"), each = 3L))
    ActiveIdent(sce) <- "cluster"
    sclet_set_analysis_state(
        object = sce, type = "detest", id = detest_id,
        method = "sclet::RunDEtest", inputs = list(),
        artifacts = list(result = markers),
        params = list(), summary = list(n_markers = nrow(markers)),
        active = TRUE
    )
}

## A DE evidence node recorded against `detest_id`, standing in for what
## run_de_test registers, so the enrichment node has a real parent to name.
sclet_ai_test_enrichment_de_evidence <- function(object, detest_id = "de_ai") {
    RecordAIEvidence(
        object,
        list(
            id = paste0("ev:detest_", detest_id),
            kind = "deterministic_summary",
            values = list(n_rows = 20L, n_significant = 8L, raw_values_included = FALSE),
            claim_level = "associated"
        ),
        source = detest_id, parents = character(), scope = NULL
    )
}

sclet_ai_test_enrichment_registry <- function(object) {
    AIDefaultExecutionRegistry(object, include = c("read", "enrichment"))
}

sclet_ai_test_enrichment_plan <- function(object, registry, params) {
    plan <- sclet:::new_sclet_ai_plan(
        task = "enrichment",
        context_fingerprint = GetAnalysisLedger(object)$fingerprint,
        actions = list(list(id = "en1", action = "run_enrichment", params = as.list(params)))
    )
    ValidateAIPlan(plan, object = object, registry = registry, strict = FALSE)
}

## Every character value anywhere in a nested object, so a leak cannot hide in
## a list element nobody thought to verify.
sclet_ai_test_enrichment_characters <- function(x) {
    if (is.character(x)) return(x)
    if (is.list(x)) return(unlist(lapply(x, sclet_ai_test_enrichment_characters), use.names = FALSE))
    character()
}

## Identity tokens worth scanning for: real gene symbols only. Empty strings
## and "NA" are dropped because grepl("", x) is TRUE for every x and the
## literal "NA" would match the missing-value markers str() prints, so keeping
## them would make the leak check fail without meaning anything.
sclet_ai_test_enrichment_identity_tokens <- function(values) {
    values <- as.character(values[!is.na(values)])
    values <- unlist(strsplit(values, "/", fixed = TRUE), use.names = FALSE)
    values <- unique(trimws(values))
    values <- values[nzchar(values) & values != "NA"]
    values[grepl("^[A-Za-z][A-Za-z0-9.-]{1,}$", values)]
}

## -- Block 1: readiness is blocked with no detest, and names the fix -------
test_that("check_enrichment_readiness is blocked on a bare object and names RunDEtest()", {
    bare <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1L, 3L, 4L)))
    result <- check_enrichment_readiness(bare)
    expect_equal(result$status, "not_ready")
    expect_true("enrichment" %in% result$blocked_actions)
    expect_true(length(result$questions) > 0L)
    expect_true(any(grepl("RunDEtest", result$questions)))
    expect_false(isTRUE(result$checks$detest_state_available))
    expect_false(isTRUE(result$checks$raw_values_included))
})

## -- Block 2: readiness passes on a real marker table ---------------------
test_that("check_enrichment_readiness is ready once a detest table with clusters exists", {
    sce <- sclet_ai_test_enrichment_object()
    result <- check_enrichment_readiness(sce)
    expect_equal(result$status, "ready_for_diagnostic")
    expect_equal(length(result$blocked_actions), 0L)
    expect_true(result$checks$detest_state_available)
    expect_true(result$checks$detest_marker_column_has_cluster)
    expect_equal(result$checks$detest_id_resolved, "de_ai")
    expect_false(isTRUE(result$checks$raw_values_included))
    ## No modelling choice is offered.
    expect_null(result$recommendation)
})

## -- Block 3: the degeneracy gate -----------------------------------------
## RunEnrichment() sets every marker to "Cluster1" when the table has no
## cluster column: a warning, not an error, and a comparison of one pooled set
## against itself. Readiness must make that hazard visible rather than absorb
## it, because a planning agent cannot see the warning.
test_that("readiness blocks the single-pooled-group degeneracy instead of reporting ready", {
    sce <- sclet_ai_test_enrichment_object(with_cluster = FALSE)
    result <- check_enrichment_readiness(sce)
    expect_false(result$checks$detest_marker_column_has_cluster)
    expect_equal(result$status, "not_ready")
    expect_true("enrichment" %in% result$blocked_actions)
    expect_true(any(grepl("single pooled group|pooled", result$questions)))
    ## The question must NOT pick a cluster column for the user.
    expect_false(any(grepl("suggest|recommend|use column|建议", result$questions)))

    ## And the hazard is real, not hypothetical: prove the underlying function
    ## does run on this input and does collapse everything into one pooled
    ## group, which is precisely the meaningless comparison the gate blocks.
    degenerate <- suppressWarnings(
        RunEnrichment(sce, db = "GO", name = "degenerate_run")
    )
    table <- as.data.frame(get_enrichment(degenerate, id = "degenerate_run")$artifacts$result)
    expect_gt(nrow(table), 0L)
    expect_equal(length(unique(as.character(table$cluster))), 1L)
})

## -- Block 4: real enrichment end to end ---------------------------------
test_that("run_enrichment executes a real enrichment and records bounded evidence", {
    sce <- sclet_ai_test_enrichment_de_evidence(sclet_ai_test_enrichment_object())
    registry <- sclet_ai_test_enrichment_registry(sce)
    validation <- sclet_ai_test_enrichment_plan(sce, registry, list(db = "GO", name = "run1"))
    expect_true(validation$valid)

    executed <- ExecuteAIPlan(sce, validation$plan, registry, validation = validation,
        dry_run = FALSE, confirmation = validation$confirmation_token)
    expect_equal(executed$status, "completed")
    expect_true(has_enrichment(executed$object, "run1"))

    nodes <- sclet:::sclet_ai_evidence_get_all(executed$object)
    expect_true("ev:enrichment_run1" %in% names(nodes))
    expect_equal(nodes[["ev:enrichment_run1"]]$claim_level, "consistent_with")

    ## It really computed something: the result table is non-empty.
    table <- as.data.frame(get_enrichment(executed$object, id = "run1")$artifacts$result)
    expect_gt(nrow(table), 0L)
})

## -- Block 5: bounded evidence / privacy ---------------------------------
## The heart of the contract. The assertions are built from what the real run
## actually returned, so they track the data instead of a guessed list.
test_that("enrichment evidence leaks no gene symbol, category name or accession", {
    sce <- sclet_ai_test_enrichment_de_evidence(sclet_ai_test_enrichment_object())
    registry <- sclet_ai_test_enrichment_registry(sce)
    validation <- sclet_ai_test_enrichment_plan(sce, registry, list(db = "GO", name = "run2"))
    executed <- ExecuteAIPlan(sce, validation$plan, registry, validation = validation,
        dry_run = FALSE, confirmation = validation$confirmation_token)
    node <- sclet:::sclet_ai_evidence_get_all(executed$object)[["ev:enrichment_run2"]]
    dumped <- paste(capture.output(str(node)), collapse = "\n")

    result <- as.data.frame(get_enrichment(executed$object, id = "run2")$artifacts$result)
    real_genes <- unique(unlist(strsplit(as.character(result$geneID), "/"), use.names = FALSE))
    real_terms <- unique(as.character(result$Description))
    real_ids <- unique(as.character(result$ID))
    ## Sanity: the fixture really produced checkable identity content.
    expect_gt(length(real_genes), 0L)
    expect_gt(length(real_terms), 0L)
    expect_gt(length(real_ids), 0L)

    for (gene in real_genes) expect_false(grepl(gene, dumped, fixed = TRUE))
    for (term in real_terms) expect_false(grepl(term, dumped, fixed = TRUE))
    for (id in real_ids) expect_false(grepl(id, dumped, fixed = TRUE))
    ## Column names that would betray the payload shape are absent too.
    expect_false(grepl("geneID", dumped, fixed = TRUE))
    expect_false(grepl("Description", dumped, fixed = TRUE))

    ## Every character value anywhere in the node is an anonymized code.
    codes <- sclet_ai_test_enrichment_characters(node$values)
    expect_gt(length(codes), 0L)
    expect_true(all(grepl("^(group|cluster|condition|sample|route)_[0-9]+$", codes)))

    expect_false(isTRUE(node$values$raw_values_included))
    expect_true(isTRUE(node$values$enrichment_is_derived_from_markers))
    ## No per-cell / per-category vector smuggled through: nothing as long as
    ## the number of enriched rows it summarises.
    expect_lt(length(node$values$enriched_per_cluster), 1000L)

    ## The leak detector itself has teeth: a node that DID carry a gene is
    ## caught by the same predicate.
    tampered <- node$values
    tampered$decoy <- real_genes[[1L]]
    tampered_dump <- paste(capture.output(str(tampered)), collapse = "\n")
    expect_true(any(vapply(real_genes, function(g) grepl(g, tampered_dump, fixed = TRUE), logical(1L))))
})

## -- Block 6: the lineage assertion --------------------------------------
## DE evidence + enrichment evidence = ONE line, because enrichment
## re-describes the DE measurement. The control proves the assertion can fail.
test_that("DE plus enrichment evidence counts as one independent line of support", {
    sce <- sclet_ai_test_enrichment_de_evidence(sclet_ai_test_enrichment_object())
    registry <- sclet_ai_test_enrichment_registry(sce)
    validation <- sclet_ai_test_enrichment_plan(sce, registry, list(db = "GO", name = "run3"))
    executed <- ExecuteAIPlan(sce, validation$plan, registry, validation = validation,
        dry_run = FALSE, confirmation = validation$confirmation_token)
    object <- executed$object

    nodes <- sclet:::sclet_ai_evidence_get_all(object)
    enrichment <- nodes[["ev:enrichment_run3"]]
    ## The node names its DE parent explicitly.
    expect_true("ev:detest_de_ai" %in% enrichment$parents)

    ceiling <- check_claim_ceiling(object, evidence_ids = c("ev:detest_de_ai", "ev:enrichment_run3"))
    expect_equal(ceiling$n_independent_groups, 1L)

    ## CONTROL: the same DE node plus evidence from a genuinely unrelated
    ## source (different source, different dependency_group, no parentage)
    ## must count as two. Without this the assertion above could pass because
    ## the ceiling engine never separates anything.
    unrelated <- RecordAIEvidence(
        object,
        list(id = "ev:unrelated_control", kind = "deterministic_summary",
            values = list(n_observed = 5L, raw_values_included = FALSE),
            claim_level = "associated"),
        source = NULL, parents = character(), scope = list(control_marker = 1L)
    )
    control <- check_claim_ceiling(unrelated,
        evidence_ids = c("ev:detest_de_ai", "ev:unrelated_control"))
    expect_equal(control$n_independent_groups, 2L)
})

## -- Block 7: typed prerequisites ----------------------------------------
test_that("run_enrichment rejects a missing detest and a bad database", {
    bare <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1L, 3L, 4L)))
    registry <- sclet_ai_test_enrichment_registry(bare)
    ## ValidateAIPlan prefixes errors with "prerequisites not met for action
    ## en1: ", so the token is matched anywhere and never with a ^ anchor.
    missing <- sclet_ai_test_enrichment_plan(bare, registry, list())
    expect_false(missing$valid)
    expect_true(any(grepl("detest_state_missing", missing$errors)))

    sce <- sclet_ai_test_enrichment_object()
    bad_db <- sclet_ai_test_enrichment_plan(sce, sclet_ai_test_enrichment_registry(sce),
        list(db = "not_a_database"))
    expect_false(bad_db$valid)
    expect_true(any(grepl("db_invalid", bad_db$errors)))
    ## The refusal says it will not choose for the user.
    expect_true(any(grepl("never selected automatically|modelling choice", bad_db$errors)))

    bad_id <- sclet_ai_test_enrichment_plan(sce, sclet_ai_test_enrichment_registry(sce),
        list(detest_id = "nope"))
    expect_false(bad_id$valid)
    expect_true(any(grepl("detest_state_missing", bad_id$errors)))
})

## -- Block 8: the silent no-result path must not report success ----------
## RunEnrichment() warns and returns the object unchanged when nothing passes
## the cutoffs, writing no state. Reporting that as a completed step would
## tell a planning agent the analysis produced a result.
test_that("exhausted cutoffs surface a typed failure instead of a fake success", {
    sce <- sclet_ai_test_enrichment_de_evidence(sclet_ai_test_enrichment_object())
    registry <- sclet_ai_test_enrichment_registry(sce)
    validation <- sclet_ai_test_enrichment_plan(sce, registry,
        list(db = "GO", name = "empty_run", p_val_cutoff = 1e-300))
    expect_true(validation$valid)

    executed <- withCallingHandlers(
        ExecuteAIPlan(sce, validation$plan, registry, validation = validation,
            dry_run = FALSE, confirmation = validation$confirmation_token),
        warning = function(w) invokeRestart("muffleWarning")
    )
    expect_equal(executed$status, "failed")
    messages <- unlist(lapply(executed$results, function(step) step$error))
    expect_true(any(grepl("enrichment_no_result", messages)))
    ## Nothing was registered as evidence for a run that produced nothing.
    expect_false(has_enrichment(executed$object, "empty_run"))
    expect_false("ev:enrichment_empty_run" %in%
        names(sclet:::sclet_ai_evidence_get_all(executed$object)))
})

## -- Block 9: the read-only summary --------------------------------------
test_that("summarize_enrichment_evidence is bounded and offers no recommendation", {
    bare <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1L, 3L, 4L)))
    absent <- summarize_enrichment_evidence(bare)
    expect_equal(absent$status, "not_available")
    expect_equal(absent$reason, "no_enrichment_record")

    sce <- sclet_ai_test_enrichment_de_evidence(sclet_ai_test_enrichment_object())
    registry <- sclet_ai_test_enrichment_registry(sce)
    validation <- sclet_ai_test_enrichment_plan(sce, registry, list(db = "GO", name = "run4"))
    executed <- ExecuteAIPlan(sce, validation$plan, registry, validation = validation,
        dry_run = FALSE, confirmation = validation$confirmation_token)
    summary <- summarize_enrichment_evidence(executed$object, name = "run4")

    expect_equal(summary$status, "available")
    expect_null(summary$recommendation)
    expect_false(isTRUE(summary$raw_values_included))
    expect_true(summary$derived_from_evidence_node)
    expect_equal(summary$lineage_parents, "ev:detest_de_ai")
    expect_true(all(grepl("^(group|cluster|condition|sample|route)_[0-9]+$",
        sclet_ai_test_enrichment_characters(summary$cluster_codes))))

    ## Same identity ban as the evidence node.
    result <- as.data.frame(get_enrichment(executed$object, id = "run4")$artifacts$result)
    dumped <- paste(capture.output(str(summary)), collapse = "\n")
    for (gene in unique(unlist(strsplit(as.character(result$geneID), "/"), use.names = FALSE))) {
        expect_false(grepl(gene, dumped, fixed = TRUE))
    }
    for (term in unique(as.character(result$Description))) {
        expect_false(grepl(term, dumped, fixed = TRUE))
    }
})

## -- Block 10: strictly read-only ---------------------------------------
test_that("the enrichment diagnostics never mutate the object", {
    sce <- sclet_ai_test_enrichment_de_evidence(sclet_ai_test_enrichment_object())
    before_ledger <- GetAnalysisLedger(sce)$fingerprint
    before_evidence <- paste(sort(names(sclet:::sclet_ai_evidence_get_all(sce))), collapse = ",")

    check_enrichment_readiness(sce)
    check_enrichment_readiness(sce, detest_id = "de_ai")
    summarize_enrichment_evidence(sce, name = "nothing_here")
    check_claim_ceiling(sce)

    expect_equal(GetAnalysisLedger(sce)$fingerprint, before_ledger)
    expect_equal(paste(sort(names(sclet:::sclet_ai_evidence_get_all(sce))), collapse = ","),
        before_evidence)
})

## -- Block 11: the registry surface is opt-in and read-only where it should be
test_that("enrichment is an opt-in group and defaults hide it", {
    bare <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1L, 3L, 4L)))
    default_registry <- AIDefaultExecutionRegistry(bare)
    expect_false("run_enrichment" %in% names(default_registry))

    all_registry <- AIDefaultExecutionRegistry(bare, include = "all")
    expect_true("run_enrichment" %in% names(all_registry))

    descriptor <- all_registry$run_enrichment
    ## It mutates: that is the point of an analysis action, and it is why
    ## confirmation is required.
    expect_true(descriptor$mutates_object)
    expect_true(descriptor$requires_confirmation)
    expect_equal(descriptor$returns, "sce")
    expect_equal(descriptor$estimated_cost, "medium")
    expect_true(all(c("detest", "enrichment", "ai_evidence") %in% descriptor$allowed_state_types))
})
