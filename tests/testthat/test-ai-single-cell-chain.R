## Testing phase: does the ledger actually let an AI plan, execute and interpret
## an advanced single-cell chain?
##
## The batch suites (test-ai-enrichment.R, test-ai-communication.R,
## test-ai-interpretation-gaps.R) each unit-test one stage, and each builds the
## record it expects upstream to have produced. That leaves the seams untested.
## This file drives the chain for real, on the same object, and checks three
## things no single-stage test can:
##
##   1. the plan layer enforces the ordering by prerequisite rather than by us
##      guessing intent;
##   2. the enrichment evidence really links to the DE evidence the DE action
##      actually wrote -- `ev:findallmarkers` with source = the detest record,
##      which is NOT the shape the unit fixtures used;
##   3. the claim ceiling reports ONE line of support over the whole executed
##      chain, because enrichment re-describes the DE measurement instead of
##      measuring it again.
##
## DE and enrichment run for real (clusterProfiler + org.Hs.eg.db are installed
## and a real hemoglobin-gene marker table yields real GO hits). Only Milo and
## CellChat are mocked, for the reasons recorded in their own suites: their
## third-party fits fail on synthetic counts, so driving them would assert
## against those packages' data requirements rather than against this contract.

## ---------------------------------------------------------------------------
## Fixtures
## ---------------------------------------------------------------------------

## 40 real human symbols. The first eight are the erythrocyte/hemoglobin module
## and are strongly up in cluster_1, so a real GO run has something true to
## find rather than something we assert about noise.
sclet_ai_test_chain_genes <- function() {
    c("HBB", "HBA1", "HBD", "HBE1", "HBG1", "HBG2", "AHSP", "ALAS2",
        "TP53", "ACTB", "GAPDH", "IL6", "TNF", "CDK1", "CCNB1", "MCM2",
        "PCNA", "MKI67", "BRCA1", "ATM", "RAD51", "EGFR", "MET", "ALK",
        "NTRK1", "FGFR1", "KRAS", "CD74", "HLA-A", "CXCL12", "ACKR1",
        "CLDN5", "PECAM1", "COL1A1", "FN1", "VEGFA", "KRT8", "EPCAM",
        "C3", "CCR7")
}

sclet_ai_test_chain_object <- function(n_per_group = 25L) {
    genes <- sclet_ai_test_chain_genes()
    set.seed(9L)
    mu <- matrix(5, nrow = length(genes), ncol = 2L)
    mu[seq_len(8L), 1L] <- 120
    mu[seq_len(8L), 2L] <- 3
    counts <- matrix(NA_integer_, length(genes), 2L * n_per_group)
    for (j in seq_along(genes)) {
        counts[j, seq_len(n_per_group)] <- stats::rpois(n_per_group, mu[j, 1L])
        counts[j, seq.int(to = 2L * n_per_group, length.out = n_per_group)] <-
            stats::rpois(n_per_group, mu[j, 2L])
    }
    rownames(counts) <- genes
    colnames(counts) <- paste0("cell_", seq_len(ncol(counts)))
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
    sce$cluster <- factor(rep(c("cluster_1", "cluster_2"), each = n_per_group))
    ActiveIdent(sce) <- "cluster"
    sce
}

## Restores the provider stub, actively unsetting it when it was not set
## before: a no-op restore leaks the stub into every later suite.
sclet_ai_test_chain_stub <- function(handler) {
    previous <- getOption("sclet.ai.call", NULL)
    options(sclet.ai.call = handler)
    function() options(sclet.ai.call = previous)
}

## Values only, flattened. capture.output(str()) is not used for leak scans
## because it elides long character vectors, so a real leak could hide in an
## omitted section and the assertion would pass vacuously.
sclet_ai_test_chain_values_text <- function(node) {
    paste(vapply(node$values, function(v) paste(v, collapse = ","), character(1)),
        collapse = "|")
}

## Execute a multi-step plan through the real orchestrator layer.
sclet_ai_test_chain_run <- function(object, actions, groups) {
    registry <- AIDefaultExecutionRegistry(object, include = c("read", groups))
    plan <- sclet:::new_sclet_ai_plan(
        task = "chain",
        context_fingerprint = GetAnalysisLedger(object)$fingerprint,
        actions = actions
    )
    validation <- ValidateAIPlan(plan, object = object, registry = registry, strict = FALSE)
    list(registry = registry, validation = validation)
}

sclet_ai_test_chain_execute <- function(object, actions, groups) {
    prepared <- sclet_ai_test_chain_run(object, actions, groups)
    expect_true(prepared$validation$valid)
    ExecuteAIPlan(object, prepared$validation$plan, prepared$registry,
        validation = prepared$validation, dry_run = FALSE,
        confirmation = prepared$validation$confirmation_token)
}

## --------------------------------------------------------------------------
## 1. The plan layer enforces ordering through prerequisites
## --------------------------------------------------------------------------
test_that("enrichment cannot be planned before a marker test exists", {
    sce <- sclet_ai_test_chain_object()

    ## Standalone enrichment: nothing upstream, so the plan must be refused.
    alone <- sclet_ai_test_chain_run(sce, list(list(
        id = "enr1", action = "run_enrichment", params = list(db = "GO", name = "chain_enr"))),
        "enrichment")
    expect_false(alone$validation$valid)
    expect_true(any(grepl("detest_state_missing", alone$validation$errors)))
    ## The refusal is self-repairing: it names the action to run first.
    expect_true(any(grepl("RunDEtest", alone$validation$errors)))

    ## Readiness agrees with the gate rather than green-lighting it.
    ready <- check_enrichment_readiness(sce)
    expect_equal(ready$status, "not_ready")
    expect_true("enrichment" %in% ready$blocked_actions)
})

## --------------------------------------------------------------------------
## 2-5. The real chain, executed once
## --------------------------------------------------------------------------
test_that("a planned DE -> enrichment chain executes end to end", {
    skip_if_not_installed("clusterProfiler")
    skip_if_not_installed("org.Hs.eg.db")
    sce <- sclet_ai_test_chain_object()

    executed <- sclet_ai_test_chain_execute(sce, list(
        list(id = "de1", action = "run_de_test", params = list(name = "chain_de")),
        list(id = "enr1", action = "run_enrichment",
            params = list(detest_id = "chain_de", db = "GO", name = "chain_enr"))
    ), c("annotation", "enrichment"))

    expect_equal(executed$status, "completed")
    expect_true(has_detest(executed$object, "chain_de"))
    expect_true(has_enrichment(executed$object, "chain_enr"))

    ## Both stages really produced content, not empty records.
    markers <- get_detest(executed$object, "chain_de")$artifacts$result
    expect_gt(nrow(markers), 0L)
    enrichment <- as.data.frame(get_enrichment(executed$object, id = "chain_enr")$artifacts$result)
    expect_gt(nrow(enrichment), 0L)
})

test_that("the executed chain links enrichment evidence to the real DE evidence", {
    skip_if_not_installed("clusterProfiler")
    skip_if_not_installed("org.Hs.eg.db")
    sce <- sclet_ai_test_chain_object()
    executed <- sclet_ai_test_chain_execute(sce, list(
        list(id = "de1", action = "run_de_test", params = list(name = "chain_de")),
        list(id = "enr1", action = "run_enrichment",
            params = list(detest_id = "chain_de", db = "GO", name = "chain_enr"))
    ), c("annotation", "enrichment"))

    nodes <- sclet:::sclet_ai_evidence_get_all(executed$object)
    expect_true("ev:findallmarkers" %in% names(nodes))
    expect_true("ev:enrichment_chain_enr" %in% names(nodes))

    ## THE SEAM THIS TEST EXISTS FOR. The batch-A unit fixtures hand-made a DE
    ## evidence node; here the node was written by run_de_test itself, under
    ## its own id and source conventions. If the recorder's parent lookup had
    ## been tuned to the fixture instead of the real producer, this is where it
    ## breaks.
    parents <- nodes[["ev:enrichment_chain_enr"]]$parents
    expect_true("ev:findallmarkers" %in% as.character(parents))
})

test_that("one measurement described twice is still one line of support", {
    skip_if_not_installed("clusterProfiler")
    skip_if_not_installed("org.Hs.eg.db")
    sce <- sclet_ai_test_chain_object()
    executed <- sclet_ai_test_chain_execute(sce, list(
        list(id = "de1", action = "run_de_test", params = list(name = "chain_de")),
        list(id = "enr1", action = "run_enrichment",
            params = list(detest_id = "chain_de", db = "GO", name = "chain_enr"))
    ), c("annotation", "enrichment"))
    object <- executed$object

    ## Two evidence nodes, one underlying measurement.
    ceiling <- check_claim_ceiling(object,
        evidence_ids = c("ev:findallmarkers", "ev:enrichment_chain_enr"))
    expect_equal(ceiling$n_independent_groups, 1L)
    expect_equal(ceiling$ceiling_claim_level, "associated")

    ## And a claim that wants more than the evidence can carry is capped, with
    ## the request preserved so the movement is visible.
    verdict <- check_claim_ceiling(object, claims = list(list(
        claim_level = "consistent_with",
        evidence_refs = c("ev:findallmarkers", "ev:enrichment_chain_enr")
    )))$verdicts[[1L]]
    expect_equal(verdict$requested_claim_level, "consistent_with")
    expect_equal(verdict$claim_level, "associated")
    expect_equal(verdict$ceiling_status, "downgraded")

    ## CONTROL: the same two nodes plus one genuinely independent measurement
    ## does lift the ceiling, so the assertion above is not passing because
    ## nothing is ever counted as independent.
    independent <- RecordAIEvidence(object, list(
        id = "ev:independent_control", kind = "deterministic_summary",
        values = list(n_observed = 5L, raw_values_included = FALSE),
        claim_level = "associated"), source = NULL, parents = character(), scope = NULL)
    lifted <- check_claim_ceiling(independent, evidence_ids = c(
        "ev:findallmarkers", "ev:enrichment_chain_enr", "ev:independent_control"))
    expect_equal(lifted$n_independent_groups, 2L)
})

test_that("the executed chain leaks no gene or GO-term identity into evidence", {
    skip_if_not_installed("clusterProfiler")
    skip_if_not_installed("org.Hs.eg.db")
    sce <- sclet_ai_test_chain_object()
    executed <- sclet_ai_test_chain_execute(sce, list(
        list(id = "de1", action = "run_de_test", params = list(name = "chain_de")),
        list(id = "enr1", action = "run_enrichment",
            params = list(detest_id = "chain_de", db = "GO", name = "chain_enr"))
    ), c("annotation", "enrichment"))
    nodes <- sclet:::sclet_ai_evidence_get_all(executed$object)
    node <- nodes[["ev:enrichment_chain_enr"]]
    text <- sclet_ai_test_chain_values_text(node)

    ## The genes and terms below come from what the run actually produced, so
    ## this checks the real payload rather than a guessed one.
    markers <- get_detest(executed$object, "chain_de")$artifacts$result
    enrichment <- as.data.frame(get_enrichment(executed$object, id = "chain_enr")$artifacts$result)
    real_genes <- unique(c(rownames(markers),
        unlist(strsplit(as.character(enrichment$geneID), "/"), use.names = FALSE)))
    real_genes <- real_genes[nzchar(real_genes) & real_genes != "NA"]
    real_terms <- unique(as.character(enrichment$Description))
    real_ids <- unique(as.character(enrichment$ID))

    ## The identity really exists upstream: without this the scan below could
    ## pass simply because no gene or term was ever computed.
    expect_gt(length(real_genes), 0L)
    expect_gt(length(real_terms), 0L)

    for (gene in real_genes) expect_false(grepl(gene, text, fixed = TRUE))
    for (term in real_terms) expect_false(grepl(term, text, fixed = TRUE))
    for (id in real_ids) expect_false(grepl(id, text, fixed = TRUE))
    expect_false(isTRUE(node$values$raw_values_included))
    ## Nothing per-cell: the marker table is far wider than any summary vector.
    for (value in node$values) expect_true(length(value) < nrow(markers))
})

## --------------------------------------------------------------------------
## 6. The orchestrator, with a stubbed model proposing the chain
## --------------------------------------------------------------------------
test_that("RunAIAnalysis orchestrates a model-proposed DE -> enrichment chain", {
    skip_if_not_installed("clusterProfiler")
    skip_if_not_installed("org.Hs.eg.db")
    sce <- sclet_ai_test_chain_object()
    restore <- sclet_ai_test_chain_stub(function(task, context, ...) {
        list(
            answer = "Test every cluster for markers, then summarize them.",
            proposed_actions = list(
                list(id = "de1", action = "run_de_test", params = list(name = "chain_de")),
                list(id = "enr1", action = "run_enrichment",
                    params = list(detest_id = "chain_de", db = "GO", name = "chain_enr"))
            ),
            findings = list(list(
                text = "Erythrocyte-associated programs distinguish the clusters.",
                claim_level = "consistent_with",
                evidence_refs = c("ev:findallmarkers", "ev:enrichment_chain_enr")
            )),
            warnings = list(), recommendations = list()
        )
    })
    on.exit(restore(), add = TRUE)

    result <- RunAIAnalysis(sce, goal = "Describe what distinguishes the clusters.",
        confirm = "yes",
        registry = AIDefaultExecutionRegistry(sce, include = c("read", "annotation", "enrichment")))
    expect_equal(result$status, "completed")

    ## The model claimed "consistent_with" on two nodes that are one
    ## measurement. The ceiling it was audited against must cap that.
    audited <- sclet_ai_audit_result_claims(result$object, list(
        findings = list(list(
            text = "Erythrocyte-associated programs distinguish the clusters.",
            claim_level = "consistent_with",
            evidence_refs = c("ev:findallmarkers", "ev:enrichment_chain_enr")))))
    expect_equal(audited$findings[[1L]]$claim_level, "associated")
    expect_equal(audited$findings[[1L]]$requested_claim_level, "consistent_with")
})

## --------------------------------------------------------------------------
## 7. The plan-time pre-check, on a real chain
## --------------------------------------------------------------------------
test_that("the plan layer reports the ceiling the executed chain can carry", {
    skip_if_not_installed("clusterProfiler")
    skip_if_not_installed("org.Hs.eg.db")
    sce <- sclet_ai_test_chain_object()

    ## Before anything runs there is nothing to claim from.
    before <- sclet_ai_test_chain_run(sce, list(list(
        id = "de1", action = "run_de_test", params = list(name = "chain_de"))), "annotation")
    expect_equal(before$validation$claim_ceiling$status, "unsupported")

    executed <- sclet_ai_test_chain_execute(sce, list(
        list(id = "de1", action = "run_de_test", params = list(name = "chain_de")),
        list(id = "enr1", action = "run_enrichment",
            params = list(detest_id = "chain_de", db = "GO", name = "chain_enr"))
    ), c("annotation", "enrichment"))

    after <- sclet_ai_test_chain_run(executed$object, list(list(
        id = "de1", action = "run_de_test", params = list(name = "chain_de2"))), "annotation")
    expect_equal(after$validation$claim_ceiling$n_independent_groups, 1L)
    expect_equal(after$validation$claim_ceiling$ceiling_claim_level, "associated")
    ## Informational only: it must not have changed the verdict.
    expect_true(after$validation$valid)
})

## --------------------------------------------------------------------------
## 8. Read-only guarantee across the whole surface
## --------------------------------------------------------------------------
test_that("every interpretive helper leaves the executed chain untouched", {
    skip_if_not_installed("clusterProfiler")
    skip_if_not_installed("org.Hs.eg.db")
    sce <- sclet_ai_test_chain_object()
    executed <- sclet_ai_test_chain_execute(sce, list(
        list(id = "de1", action = "run_de_test", params = list(name = "chain_de")),
        list(id = "enr1", action = "run_enrichment",
            params = list(detest_id = "chain_de", db = "GO", name = "chain_enr"))
    ), c("annotation", "enrichment"))
    object <- executed$object

    fingerprint <- GetAnalysisLedger(object)$fingerprint
    evidence_before <- sort(names(sclet:::sclet_ai_evidence_get_all(object)))

    check_enrichment_readiness(object)
    check_de_readiness(object)
    check_communication_readiness(object)
    summarize_enrichment_evidence(object, name = "chain_enr")
    check_claim_ceiling(object)
    check_claim_ceiling(object, claims = list(list(
        claim_level = "measured", evidence_refs = "ev:findallmarkers")))

    expect_equal(GetAnalysisLedger(object)$fingerprint, fingerprint)
    expect_equal(sort(names(sclet:::sclet_ai_evidence_get_all(object))), evidence_before)
})

## --------------------------------------------------------------------------
## 9. Deferred validation is enforced at execution, not just waived
## --------------------------------------------------------------------------
## A chained plan must be plannable, so the marker-shape gate cannot run at
## validation time. That is only safe if the same gate then runs when the
## upstream step has really produced its table. Here the upstream step produces
## a table with no cluster column -- the exact input that makes RunEnrichment()
## quietly pool every marker into one group and still write a non-empty,
## successful-looking record -- so if the execution-time check were missing,
## this plan would complete and register evidence for a meaningless comparison.
test_that("a chained plan whose DE table lacks clusters fails at execution", {
    skip_if_not_installed("clusterProfiler")
    sce <- sclet_ai_test_chain_object()

    # Faithful to a real RunDEtest() marker table in every column EXCEPT
    # `cluster`, which is the condition under test. The stub must satisfy the
    # producer's own downstream code: run_de_test's evidence block reads the
    # table, and a thin stub dies inside that step with "subscript out of
    # bounds", so the chain never reaches run_enrichment and the assertion
    # below would fail for the wrong reason -- or, if the expectation were
    # loosened, pass for it.
    markers_without_cluster <- data.frame(
        pval = rep(1e-10, nrow(sce)),
        padj = rep(1e-9, nrow(sce)),
        avg_log2FC = rep(3, nrow(sce)),
        pct.1 = rep(1, nrow(sce)),
        pct.2 = rep(0, nrow(sce)),
        gene = rownames(sce),
        avgExpr = rep(2, nrow(sce)),
        statistic = rep(5, nrow(sce)),
        auc = rep(1, nrow(sce)),
        row.names = rownames(sce),
        stringsAsFactors = FALSE
    )

    actions <- list(
        list(id = "de1", action = "run_de_test", params = list(name = "bad_de")),
        list(id = "enr1", action = "run_enrichment",
            params = list(detest_id = "bad_de", db = "GO", name = "bad_enr"))
    )
    prepared <- sclet_ai_test_chain_run(sce, actions, c("annotation", "enrichment"))

    ## Planning must accept it: the upstream record is promised by an earlier
    ## step of the same plan, so refusing here would make chains unplannable.
    expect_true(prepared$validation$valid)

    executed <- with_mocked_bindings(
        ExecuteAIPlan(sce, prepared$validation$plan, prepared$registry,
            validation = prepared$validation, dry_run = FALSE,
            confirmation = prepared$validation$confirmation_token),
        RunDEtest = function(object, ...) {
            sclet_set_analysis_state(
                object = object, type = "detest", id = "bad_de",
                method = "FindAllMarkers", inputs = list(),
                artifacts = list(result = markers_without_cluster),
                params = list(), summary = list(n_markers = nrow(markers_without_cluster)),
                active = TRUE
            )
        },
        .package = "sclet"
    )

    ## Execution must refuse it, with the same typed reason the standalone path
    ## would have produced at validation time.
    expect_equal(executed$status, "failed")
    messages <- unlist(lapply(executed$results, function(step) step$error))
    expect_true(any(grepl("cluster_identity_missing", messages)))
    expect_false(has_enrichment(executed$object, "bad_enr"))
    expect_false("ev:enrichment_bad_enr" %in%
        names(sclet:::sclet_ai_evidence_get_all(executed$object)))
})

## --------------------------------------------------------------------------
## 10. The projection must not over-promise
## --------------------------------------------------------------------------
## Accepting a chained plan is only half the fix. If validation trusted any
## planned record it would wave through a step that consumes a record no
## earlier step produces, which is the same hole on the other side.
test_that("planned records are honoured only when an earlier step produces them", {
    sce <- sclet_ai_test_chain_object()
    # spliced/unspliced must be ASSAYS. Measured here: `sce$spliced <-
    # sce$counts` adds no assay at all on a SingleCellExperiment, so the
    # velocity step's own prerequisite failed first and this test proved
    # nothing about the projection. assays(sce)$x <- is the form that works.
    SummarizedExperiment::assays(sce)$spliced <- SummarizedExperiment::assays(sce)$counts
    SummarizedExperiment::assays(sce)$unspliced <- SummarizedExperiment::assays(sce)$counts
    expect_true(all(c("spliced", "unspliced") %in% SummarizedExperiment::assayNames(sce)))
    SingleCellExperiment::reducedDim(sce, "PCA") <- matrix(
        stats::rnorm(ncol(sce) * 5L), nrow = ncol(sce), ncol = 5L,
        dimnames = list(colnames(sce), NULL))

    registry <- AIDefaultExecutionRegistry(sce,
        include = c("read", "velocity", "fate"))
    validate <- function(actions) {
        plan <- sclet:::new_sclet_ai_plan(
            task = "chain",
            context_fingerprint = GetAnalysisLedger(sce)$fingerprint,
            actions = actions)
        ValidateAIPlan(plan, object = sce, registry = registry, strict = FALSE)
    }

    ## Chained: an earlier step really produces the velocity record.
    chained <- validate(list(
        list(id = "v1", action = "run_velocity", params = list(name = "vel_chain")),
        list(id = "f1", action = "run_fate_analysis", params = list(velocity_id = "vel_chain"))))
    expect_true(chained$valid)

    ## Same plan shape, but the id names something no step produces.
    typo <- validate(list(
        list(id = "v1", action = "run_velocity", params = list(name = "vel_chain")),
        list(id = "f1", action = "run_fate_analysis", params = list(velocity_id = "not_produced"))))
    expect_false(typo$valid)
    expect_true(any(grepl("velocity_id_unknown", typo$errors)))

    ## Nothing in the plan produces a velocity record at all.
    orphan <- validate(list(
        list(id = "q1", action = "check_qc", params = list()),
        list(id = "f1", action = "run_fate_analysis", params = list())))
    expect_false(orphan$valid)
    expect_true(any(grepl("velocity_record_missing", orphan$errors)))
})
