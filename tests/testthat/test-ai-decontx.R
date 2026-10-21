# run_decontx: the ambient-RNA action (M13-a, an explicit exception to "no new actions").
#
# Why this action exists at all: `RunDecontX()` has been exported in this package
# while having NO AIAction. That meant the single correct answer to ambient-RNA
# contamination -- clean before interpreting -- was not a step a planner could
# choose. M13 measured that gap and, deliberately, closes it rather than working
# around it.
#
# The rules this file pins are the ones that make the action honest:
#
#   1. it LABELS and CLEANS; it removes no cell and switches no active assay, so
#      downstream steps must name the cleaned assay explicitly;
#   2. it requires a real integer count matrix -- a log-normalised matrix under
#      `counts` would yield a contamination number about the wrong thing;
#   3. the evidence node records a DISTRIBUTION, never per-cell values, and says
#      the number is a model estimate.

m13a_sce <- function(n_genes = 200L, n_cells = 60L, seed = 42L) {
    set.seed(seed)
    # Two cell populations plus a shared ambient pool, which is the situation
    # decontX exists for. Integer counts, as the action requires.
    #
    # Gene index blocks are derived from n_genes rather than hard-coded: the
    # first version wrote `pop2[101:120, ]`, so any fixture under 120 genes died
    # with "subscript out of bounds" and the failure looked like a defect in the
    # action.
    per_pop <- max(10L, n_cells %/% 3L)
    block <- max(3L, n_genes %/% 6L)
    pop1_genes <- seq_len(block)
    pop2_genes <- seq.int(n_genes - block + 1L, n_genes)
    ambient_genes <- seq_len(min(5L, block))

    pop1 <- matrix(stats::rpois(n_genes * per_pop, 6), nrow = n_genes)
    pop2 <- matrix(stats::rpois(n_genes * per_pop, 2), nrow = n_genes)
    pop1[pop1_genes, ] <- pop1[pop1_genes, ] + 8L
    pop2[pop2_genes, ] <- pop2[pop2_genes, ] + 8L
    ambient <- matrix(stats::rpois(n_genes * per_pop, 1), nrow = n_genes)
    # ambient genes appear in every cell at low level
    ambient[ambient_genes, ] <- ambient[ambient_genes, ] + 3L

    counts <- cbind(pop1, pop2, ambient)
    storage.mode(counts) <- "integer"
    SingleCellExperiment::SingleCellExperiment(
        list(counts = Matrix::Matrix(counts, sparse = TRUE))
    )
}

test_that("run_decontx is registered and declares what it writes", {
    sce <- m13a_sce()
    registry <- AIDefaultExecutionRegistry(sce, include = "preprocess")
    expect_true("run_decontx" %in% names(registry))
    action <- registry[["run_decontx"]]
    expect_true(isTRUE(action$mutates_object))
    expect_true(isTRUE(action$requires_confirmation))
    expect_false(isTRUE(action$idempotent))
    expect_identical(action$output_schema$required_assays, "decontXcounts")
    expect_identical(action$output_schema$required_columns,
        "decontX_contamination")
})

test_that("a non-count assay is refused with the gate name, not silently used", {
    sce <- m13a_sce()
    # log-normalised values under the `counts` name: the shape that would make
    # decontX report a contamination number about the wrong matrix
    logcounts <- log1p(as.matrix(SummarizedExperiment::assay(sce, "counts")))
    SummarizedExperiment::assay(sce, "counts") <- logcounts
    registry <- AIDefaultExecutionRegistry(sce, include = "preprocess")
    verdict <- registry[["run_decontx"]]$prerequisites(sce, list())
    expect_true(is.character(verdict))
    expect_true(grepl("counts_assay_not_count", verdict, fixed = TRUE))
})

test_that("a missing assay is refused by name", {
    sce <- m13a_sce()
    registry <- AIDefaultExecutionRegistry(sce, include = "preprocess")
    verdict <- registry[["run_decontx"]]$prerequisites(sce, list(assay = "nope"))
    expect_true(is.character(verdict))
    expect_true(grepl("assay is not available: nope", verdict, fixed = TRUE))
})

test_that("run_decontx writes a cleaned assay, keeps the counts, and removes no cell", {
    skip_if_not(requireNamespace("celda", quietly = TRUE),
        "celda is not installed")
    sce <- m13a_sce()
    n_before <- ncol(sce)
    registry <- AIDefaultExecutionRegistry(sce, include = "preprocess")
    action <- registry[["run_decontx"]]
    out <- action$handler(sce, list(assay = "counts", name = "decontx_test"))

    # the two artefacts the action promises
    expect_true("decontXcounts" %in% SummarizedExperiment::assayNames(out))
    expect_true("decontX_contamination" %in%
        colnames(SummarizedExperiment::colData(out)))

    # the original matrix and the cell set are untouched: this action never
    # filters, which is the difference between cleaning and deleting
    expect_identical(ncol(out), n_before)
    expect_identical(dim(out), dim(sce))
    expect_true("counts" %in% SummarizedExperiment::assayNames(out))

    # contamination is a fraction per cell
    contamination <- SummarizedExperiment::colData(out)$decontX_contamination
    expect_true(all(contamination >= 0 & contamination <= 1, na.rm = TRUE))
})

test_that("the evidence node records a distribution, not per-cell values", {
    skip_if_not(requireNamespace("celda", quietly = TRUE),
        "celda is not installed")
    sce <- m13a_sce()
    registry <- AIDefaultExecutionRegistry(sce, include = "preprocess")
    out <- registry[["run_decontx"]]$handler(sce,
        list(assay = "counts", name = "decontx_ev"))

    nodes <- sclet:::sclet_ai_evidence_get_all(out)
    ids <- names(nodes)
    decontx_ids <- grep("^ev:decontx_", ids, value = TRUE)
    expect_true(length(decontx_ids) == 1L)

    node <- nodes[[decontx_ids[[1L]]]]
    values <- node$values %||% node$summary$values
    # bounded aggregates only
    for (field in c("n_cells", "contamination_mean", "contamination_median",
        "contamination_q95", "contamination_max")) {
        expect_true(field %in% names(values), info = field)
    }
    # and the two statements that keep the node honest
    expect_true(isTRUE(values$contamination_call_is_model_estimate))
    expect_true(isTRUE(values$no_cell_removed_or_filtered))
    expect_true(isTRUE(values$original_assay_preserved))
    expect_false(isTRUE(values$raw_values_included))
    # no per-cell vector leaked into the node
    for (value in values) {
        expect_true(length(value) == 1L,
            info = "an evidence value is not a single bounded aggregate")
    }
})

test_that("the node is not redacted by the privacy layer", {
    skip_if_not(requireNamespace("celda", quietly = TRUE),
        "celda is not installed")
    sce <- m13a_sce()
    registry <- AIDefaultExecutionRegistry(sce, include = "preprocess")
    out <- registry[["run_decontx"]]$handler(sce,
        list(assay = "counts", name = "decontx_priv"))

    context <- sclet_ai_context(out, target = "planning", detail = "full")
    sanitized <- sclet_ai_sanitize_context(out, context, privacy = "standard")
    redacted <- sanitized$policy$redacted_fields
    decontx_redacted <- redacted[grepl("decontx", redacted, ignore.case = TRUE)]
    expect_length(decontx_redacted, 0L)
})

test_that("decontX is deterministic, so an identical re-run is an idempotent no-op", {
    skip_if_not(requireNamespace("celda", quietly = TRUE),
        "celda is not installed")
    # Measured 2026-10-21 on the GSE158055 sample: two runs on the same input
    # gave identical contamination vectors (max abs diff 0). celda's decontX has
    # no random initialisation by default, so "same input, same numbers" is the
    # real behaviour -- and the recorder's no-op on identical values is therefore
    # the correct outcome, not a missed conflict.
    #
    # This matters because the opposite assumption is the tempting one: a
    # stochastic-looking estimator is assumed to vary, and a test asserting "a
    # re-run conflicts" would then be pinned to a false premise.
    sce <- m13a_sce(n_genes = 100L, n_cells = 30L)
    registry <- AIDefaultExecutionRegistry(sce, include = "preprocess")
    action <- registry[["run_decontx"]]

    first <- action$handler(sce, list(assay = "counts", name = "det"))
    second <- action$handler(sce, list(assay = "counts", name = "det"))
    expect_equal(
        SummarizedExperiment::colData(first)$decontX_contamination,
        SummarizedExperiment::colData(second)$decontX_contamination
    )
})

test_that("a changed input under the same node name is a conflict, not a silent overwrite", {
    skip_if_not(requireNamespace("celda", quietly = TRUE),
        "celda is not installed")
    # The complement of the test above: when the numbers DO differ, the recorder
    # must refuse rather than leave the earlier node silently describing a
    # different run. This is the B23 defect class.
    sce <- m13a_sce(n_genes = 120L, n_cells = 40L)
    registry <- AIDefaultExecutionRegistry(sce, include = "preprocess")
    action <- registry[["run_decontx"]]

    first <- action$handler(sce, list(assay = "counts", name = "det"))
    # same node name, different cells -> different aggregates
    subset <- first[, seq_len(20L)]
    conflict <- tryCatch(
        action$handler(subset, list(assay = "counts", name = "det")),
        error = function(e) conditionMessage(e)
    )
    expect_true(is.character(conflict))
    expect_true(grepl("evidence_id_conflict", conflict, fixed = TRUE))

    # and replace = TRUE is the documented way through
    refreshed <- tryCatch(
        action$handler(subset, list(assay = "counts", name = "det",
            replace_evidence = TRUE)),
        error = function(e) conditionMessage(e)
    )
    expect_false(is.character(refreshed))
})
