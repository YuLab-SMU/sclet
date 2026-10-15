## Batch C, parts 1 and 2: the doublet evidence node that was missing, and DE
## readiness for the stage that feeds enrichment.
##
## run_doublet_detection previously called scDblFinder and stopped: per-cell
## scores landed in colData and nothing reached the ledger, so a planning agent
## could execute the step and then find no evidence to cite for it. C1 closes
## that. scDblFinder is installed here, so these tests run the real classifier
## end to end rather than mocking it -- the only thing asserted about its output
## is the shape of what sclet registers, never a biological claim.
##
## C2 covers the silent hazard that degrades the whole DE -> enrichment chain:
## a one-population object makes RunDEtest() return an empty marker table
## without erroring, and that emptiness only shows up later as an enrichment run
## that recorded nothing.

## ---------------------------------------------------------------------------
## Fixtures
## ---------------------------------------------------------------------------

sclet_ai_test_gap_object <- function(n_cells = 90L, n_genes = 40L, groups = 2L) {
    counts <- matrix(rpois(n_genes * n_cells, lambda = 5L),
        nrow = n_genes, ncol = n_cells,
        dimnames = list(paste0("gene_", seq_len(n_genes)), paste0("cell_", seq_len(n_cells))))
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
    labels <- if (groups > 1L) {
        rep(paste0("cluster_", seq_len(groups)), length.out = n_cells)
    } else {
        rep("cluster_1", n_cells)
    }
    sce$cluster <- factor(labels)
    ActiveIdent(sce) <- "cluster"
    sce
}

## Every piece of text a node can carry. capture.output(str()) is deliberately
## NOT used for leak scans: it elides long character vectors, so a value that
## really leaked could sit in an omitted section and the assertion would pass
## vacuously. Same reasoning as test-ai-communication.R.
sclet_ai_test_gap_all_text <- function(x) {
    pieces <- character()
    collect <- function(value, name) {
        if (!is.null(name) && length(name)) pieces[[length(pieces) + 1L]] <<- as.character(name)
        if (is.character(value)) {
            pieces[[length(pieces) + 1L]] <<- paste(value, collapse = "\r")
        } else if (is.data.frame(value)) {
            for (column in colnames(value)) collect(value[[column]], column)
        } else if (is(value, "DataFrame") || is(value, "DFrame")) {
            for (column in colnames(value)) collect(value[[column]], column)
        } else if (is.matrix(value) || is.array(value)) {
            pieces[[length(pieces) + 1L]] <<- paste(as.character(value), collapse = "\r")
        } else if (is.list(value)) {
            if (!is.null(names(value))) {
                pieces[[length(pieces) + 1L]] <<- paste(names(value), collapse = "\r")
            }
            for (i in seq_along(value)) collect(value[[i]], names(value)[i])
        }
    }
    collect(x, NULL)
    paste(unlist(pieces, use.names = FALSE), collapse = "\n")
}

## -- C1: the doublet evidence node, real classifier ------------------------
test_that("run_doublet_detection now registers bounded evidence for its call", {
    object <- sclet_ai_test_gap_object()
    registry <- AIDefaultExecutionRegistry(object, include = c("read", "rare_cell"))
    plan <- sclet:::new_sclet_ai_plan(
        task = "doublets",
        context_fingerprint = GetAnalysisLedger(object)$fingerprint,
        actions = list(list(id = "dbl1", action = "run_doublet_detection", params = list()))
    )
    validation <- ValidateAIPlan(plan, object = object, registry = registry, strict = FALSE)
    expect_true(validation$valid)

    executed <- ExecuteAIPlan(object, validation$plan, registry, validation = validation,
        dry_run = FALSE, confirmation = validation$confirmation_token)
    expect_equal(executed$status, "completed")

    nodes <- sclet:::sclet_ai_evidence_get_all(executed$object)
    node <- nodes[["ev:doublet_doublet_finder"]]
    expect_false(is.null(node))
    ## One classifier's call is an association, never an observation of fact.
    expect_equal(node$claim_level, "associated")
    expect_equal(node$values$n_cells, as.integer(ncol(executed$object)))
    expect_true(node$values$n_called_doublets <= node$values$n_cells)
    expect_true(isTRUE(node$values$doublet_call_is_model_estimate))
    expect_true(isTRUE(node$values$no_cell_removed_or_filtered))
    expect_false(isTRUE(node$values$raw_values_included))
    ## The per-cell data really is still there for R-side use.
    expect_true("scDblFinder.class" %in% colnames(SummarizedExperiment::colData(executed$object)))
})

test_that("doublet evidence holds no per-cell vector, barcode or class label", {
    object <- sclet_ai_test_gap_object()
    scored <- suppressMessages(suppressWarnings(RunDoubletFinder(object)))
    recorded <- sclet:::sclet_ai_record_doublet_evidence(scored, name = "doublet_finder")
    node <- sclet:::sclet_ai_evidence_get_all(recorded)[["ev:doublet_doublet_finder"]]
    expect_false(is.null(node))

    values <- node$values
    ## Nothing may be per-cell: a per-cell vector would be exactly ncol() long.
    for (value in values) {
        expect_true(length(value) < ncol(recorded))
    }
    ## No character values at all: the class labels are the identity payload,
    ## and the code pattern would let anonymized values through, so banning the
    ## type entirely is the honest statement of this contract.
    expect_false(any(vapply(values, is.character, logical(1L))))
    ## The word "doublet" legitimately appears in field names and in the node
    ## id; what must not exist is a class label or barcode as a VALUE.
    dumped_values <- paste(vapply(values, function(v) paste(v, collapse = ","), character(1)),
        collapse = "|")
    expect_false(grepl("singlet", dumped_values, fixed = TRUE))
    expect_false(grepl("doublet", dumped_values, fixed = TRUE))
    for (barcode in colnames(object)[seq_len(5L)]) {
        expect_false(grepl(barcode, dumped_values, fixed = TRUE))
    }
    ## POSITIVE CONTROL: prove the scan finds these tokens when present. The
    ## class labels really do exist in the scored object (as a factor, which is
    ## why this is checked rather than assumed), and a tampered values list with
    ## a barcode in it is caught by the exact expression used above.
    class_values <- as.character(SummarizedExperiment::colData(scored)$scDblFinder.class)
    expect_true(any(grepl("singlet", class_values, fixed = TRUE)))
    tampered <- values
    tampered$decoy <- colnames(object)[[1L]]
    tampered_dump <- paste(vapply(tampered, function(v) paste(v, collapse = ","), character(1)),
        collapse = "|")
    expect_true(grepl(colnames(object)[[1L]], tampered_dump, fixed = TRUE))
    tampered$decoy <- "singlet"
    expect_true(grepl("singlet", paste(vapply(tampered, function(v) paste(v, collapse = ","),
        character(1)), collapse = "|"), fixed = TRUE))
})

test_that("re-running doublet detection neither errors nor duplicates evidence", {
    object <- sclet_ai_test_gap_object()
    once <- sclet:::sclet_ai_record_doublet_evidence(
        suppressMessages(suppressWarnings(RunDoubletFinder(object))), name = "doublet_finder")
    ids_after_one <- names(sclet:::sclet_ai_evidence_get_all(once))
    ## The action is idempotent = FALSE, so a second real run must degrade to a
    ## no-op in the recorder rather than failing on the stable evidence id.
    twice <- sclet:::sclet_ai_record_doublet_evidence(
        suppressMessages(suppressWarnings(RunDoubletFinder(once))), name = "doublet_finder")
    ids_after_two <- names(sclet:::sclet_ai_evidence_get_all(twice))
    expect_setequal(ids_after_one, ids_after_two)
    expect_equal(sum(grepl("^ev:doublet_", ids_after_two)), 1L)
})

test_that("doublet recorder is a no-op when no scores exist", {
    object <- sclet_ai_test_gap_object()
    fingerprint <- GetAnalysisLedger(object)$fingerprint
    unchanged <- sclet:::sclet_ai_record_doublet_evidence(object, name = "doublet_finder")
    expect_equal(GetAnalysisLedger(unchanged)$fingerprint, fingerprint)
    expect_length(sclet:::sclet_ai_evidence_get_all(unchanged), 0L)
})

## -- C2: DE readiness and the empty-marker-table hazard -------------------
test_that("check_de_readiness is ready with two populations and blocked with one", {
    two <- check_de_readiness(sclet_ai_test_gap_object(groups = 2L))
    expect_equal(two$status, "ready_for_diagnostic")
    expect_equal(two$checks$n_groups, 2L)
    expect_equal(length(two$blocked_actions), 0L)
    expect_false(isTRUE(two$checks$raw_values_included))

    one <- check_de_readiness(sclet_ai_test_gap_object(groups = 1L))
    expect_equal(one$status, "not_ready")
    expect_true("de" %in% one$blocked_actions)
    expect_equal(one$checks$n_groups, 1L)
    expect_true(length(one$questions) > 0L)
    expect_true(any(grepl("empty table|compare one group", one$questions)))
})

test_that("DE readiness reports the hazard without choosing a column or resolution", {
    ## A second, equally plausible group column exists; the diagnostic must
    ## report what it resolved to, and must not steer the user to the other one.
    object <- sclet_ai_test_gap_object(groups = 1L)
    object$second <- factor(rep("only", ncol(object)))
    result <- check_de_readiness(object)
    questions <- paste(result$questions, collapse = " ")
    steering <- c("instead", "prefer", "try ", "switch", "use the second",
        "better column", "resolution =", "set resolution")
    for (pattern in steering) {
        ## The disclaimer's own wording is allowed; an imperative is not.
        hits <- grep(pattern, result$questions, ignore.case = TRUE)
        for (hit in hits) {
            expect_true(grepl("no resolution or alternative column is suggested",
                result$questions[[hit]]),
                info = paste("steering language:", pattern))
        }
    }
    ## It names what it resolved (diagnostics), but never names a substitute.
    expect_equal(as.character(result$checks$group_column_resolved), "cluster")
    expect_false(any(grepl("second", result$questions, fixed = TRUE)))
})

test_that("DE readiness blocks an unresolvable grouping and tiny groups", {
    bare <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1L, 10L, 6L)))
    result <- check_de_readiness(bare)
    expect_equal(result$status, "not_ready")
    expect_false(result$checks$group_column_resolvable)
    expect_equal(result$checks$smallest_group_n, 0L)
    expect_true(any(grepl("cluster the data or set an active identity", result$questions)))

    ## Two groups but one has a single cell: no within-group variance to test.
    counts <- matrix(rpois(40L * 21L, 5L), nrow = 40L, ncol = 21L,
        dimnames = list(paste0("gene_", seq_len(40L)), paste0("cell_", seq_len(21L))))
    tiny <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
    tiny$cluster <- factor(c("cluster_1", rep("cluster_2", 20L)))
    ActiveIdent(tiny) <- "cluster"
    small <- check_de_readiness(tiny)
    expect_equal(small$status, "not_ready")
    expect_equal(small$checks$smallest_group_n, 1L)
    expect_equal(small$checks$n_groups, 2L)
    expect_true(any(grepl("within-group variance", small$questions)))
})

test_that("DE readiness reports an existing detest record without claiming for it", {
    object <- sclet_ai_test_gap_object()
    object <- sclet_set_analysis_state(
        object = object, type = "detest", id = "de_prior", method = "FindAllMarkers",
        inputs = list(), artifacts = list(result = data.frame(pval = 0.5)),
        params = list(), summary = list(n_markers = 1L), active = TRUE
    )
    result <- check_de_readiness(object)
    expect_true(result$checks$detest_record_available)
    expect_true("de_prior" %in% result$checks$available_detest_ids)
    ## Reporting that a record exists is not endorsing what it says, and the
    ## readiness summary never carries the marker table itself.
    expect_null(result$recommendation)
    expect_false(isTRUE(result$checks$raw_values_included))
    flattened <- unlist(result$checks, recursive = TRUE, use.names = FALSE)
    expect_false(any(vapply(result$checks, is.data.frame, logical(1L))))
    expect_false(any(grepl("pval", sclet_ai_test_gap_all_text(result$checks), fixed = TRUE)))
})

## -- C2/C1: read-only, as promised ---------------------------------------
test_that("the batch C diagnostics mutate nothing", {
    for (maker in list(function() sclet_ai_test_gap_object(groups = 1L),
            function() sclet_ai_test_gap_object(groups = 2L))) {
        object <- maker()
        fingerprint <- GetAnalysisLedger(object)$fingerprint
        evidence <- names(sclet:::sclet_ai_evidence_get_all(object))
        check_de_readiness(object)
        check_de_readiness(object, group = "cluster")
        invisible(capture.output(print(check_de_readiness(object)$status)))
        expect_equal(GetAnalysisLedger(object)$fingerprint, fingerprint)
        expect_equal(names(sclet:::sclet_ai_evidence_get_all(object)), evidence)
    }
})
