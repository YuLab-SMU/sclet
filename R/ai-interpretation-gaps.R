# Batch C: closing interpretation gaps for three single-cell stages that could
# previously run without leaving anything the AI layer can reason over.
#
#  - run_doublet_detection called scDblFinder and wrote per-cell colData, but
#    registered no evidence node at all, so a planning agent could execute it
#    and then find nothing in the ledger to cite.
#  - RunDEtest had no readiness check, so a single-group object (or one whose
#    smallest group has one cell) reached the DE engine and produced an empty
#    marker table, which then surfaced far away as an enrichment no-result.
#  - RunMilo takes sample and condition columns, i.e. design semantics, and the
#    existing confirmation gate is keyed on the batch role only.
#
# Shared discipline: aggregates only. No per-cell vector, no cell barcode, no
# class label, no sample or condition label enters evidence, because an evidence
# node is eligible to leave the process as provider context.

# ---------------------------------------------------------------------------
# C1: doublet evidence
# ---------------------------------------------------------------------------

# Bounded aggregates from scDblFinder's two colData columns. Returns NULL when
# the columns are absent, so a recorder never invents a verdict.
sclet_ai_doublet_facts <- function(object) {
    columns <- colnames(SummarizedExperiment::colData(object))
    score_column <- "scDblFinder.score"
    class_column <- "scDblFinder.class"
    if (!all(c(score_column, class_column) %in% columns)) return(NULL)
    coldata <- SummarizedExperiment::colData(object)
    calls <- as.character(coldata[[class_column]])
    calls <- calls[!is.na(calls)]
    if (!length(calls)) return(NULL)
    scores <- suppressWarnings(as.numeric(coldata[[score_column]]))
    scores <- scores[is.finite(scores)]
    quantiles <- if (length(scores)) {
        round(unname(stats::quantile(scores, c(0.25, 0.5, 0.75))), 6L)
    } else {
        rep(NA_real_, 3L)
    }
    n_cells <- length(calls)
    n_doublets <- sum(calls == "doublet")
    list(
        n_cells = as.integer(n_cells),
        n_called_doublets = as.integer(n_doublets),
        doublet_fraction = round(n_doublets / n_cells, 6L),
        score_quantile_25 = as.numeric(quantiles[[1L]]),
        score_median = as.numeric(quantiles[[2L]]),
        score_quantile_75 = as.numeric(quantiles[[3L]]),
        score_max = if (length(scores)) round(max(scores), 6L) else NA_real_
    )
}

#' Record a bounded evidence node for a doublet call
#'
#' Registers what scDblFinder concluded as aggregate counts and score quantiles
#' only. The per-cell class labels and per-cell scores stay in colData: a score
#' per cell is cell-level identity plus a model output, and an evidence node can
#' be sent to a provider as context. The claim level is deliberately
#' \code{"associated"} rather than \code{"observed"} or
#' \code{"consistent_with"}: one classifier's call, with no independent signal
#' corroborating it, is not a measurement of doubleness.
#'
#' @param object A \code{SingleCellExperiment} object.
#' @param name Preprocess record id carrying the doublet call.
#' @return The object with one \code{ev:doublet_<name>} node appended, or
#'   unchanged when no doublet columns exist or the node already exists.
#' @keywords internal
#' @noRd
sclet_ai_record_doublet_evidence <- function(object, name = "doublet_finder") {
    facts <- sclet_ai_doublet_facts(object)
    if (is.null(facts)) return(object)
    evidence <- list(
        id = paste0("ev:doublet_", name),
        kind = "deterministic_summary",
        values = c(facts, list(
            doublet_call_is_model_estimate = TRUE,
            no_cell_removed_or_filtered = TRUE,
            raw_values_included = FALSE
        )),
        claim_level = "associated"
    )
    tryCatch(
        RecordAIEvidence(object, evidence, source = name, parents = character(), scope = NULL),
        error = function(e) object
    )
}

# ---------------------------------------------------------------------------
# C2: differential-expression readiness
# ---------------------------------------------------------------------------

#' Check whether an object is ready for differential-expression testing
#'
#' Reports whether an active group identity resolves to at least two
#' populations with enough cells to test, and whether a detest record already
#' exists. It deliberately does not choose a group column, a clustering
#' resolution or a marker threshold: those are user decisions.
#'
#' @param object A \code{SingleCellExperiment} object.
#' @param group Optional group column name. When NULL it resolves from the
#'   active identity.
#' @return A read-only readiness summary list.
#' @export
check_de_readiness <- function(object, group = NULL) {
    if (!inherits(object, "SingleCellExperiment")) {
        stop("object must be a SingleCellExperiment", call. = FALSE)
    }
    assay_names <- SummarizedExperiment::assayNames(object)
    counts_available <- "counts" %in% assay_names
    columns <- colnames(SummarizedExperiment::colData(object))

    group_key <- group %||% tryCatch(ActiveIdent(object), error = function(e) NULL)
    if (!is.null(group_key) && !as.character(group_key) %in% columns &&
            !identical(as.character(group_key), "colLabels")) {
        group_key <- NULL
    }
    if (is.null(group_key)) {
        candidates <- grep("cluster|leiden|louvain|ident", columns,
            ignore.case = TRUE, value = TRUE)
        if (length(candidates)) group_key <- candidates[[1L]]
    }

    n_groups <- 0L
    smallest <- 0L
    if (!is.null(group_key)) {
        values <- tryCatch(
            if (identical(as.character(group_key), "colLabels"))
                as.character(SingleCellExperiment::colLabels(object))
            else as.character(SummarizedExperiment::colData(object)[[group_key]]),
            error = function(e) NULL
        )
        if (!is.null(values)) {
            counts <- as.integer(table(values[!is.na(values)]))
            n_groups <- as.integer(length(counts))
            smallest <- if (length(counts)) min(counts) else 0L
        }
    }
    group_resolvable <- !is.null(group_key) && n_groups > 0L
    enough_groups <- n_groups >= 2L
    enough_cells <- smallest >= 2L

    records <- tryCatch(sclet_get_state_records(object, "detest"), error = function(e) list())
    if (is.null(records)) records <- list()
    detest_ids <- names(records)
    de_package_available <- base::requireNamespace("limma", quietly = TRUE)

    questions <- character()
    if (!group_resolvable) {
        questions <- c(questions, paste0(
            "No group identity is resolvable, and a marker test compares groups against each other; ",
            "cluster the data or set an active identity first. Available colData columns: ",
            paste(columns, collapse = ", "), "."
        ))
    } else if (!enough_groups) {
        questions <- c(questions, paste0(
            "The resolved grouping '", group_key, "' has a single population, so every marker test would ",
            "compare one group against nothing and return an empty table. Re-cluster so at least two ",
            "populations exist; no resolution or alternative column is suggested here."
        ))
    }
    if (group_resolvable && enough_groups && !enough_cells) {
        questions <- c(questions, paste0(
            "The smallest population holds ", smallest, " cell(s); a two-group comparison with fewer than ",
            "two observations per group has no within-group variance to test, so a marker table can come ",
            "back empty or degenerate rather than erroring."
        ))
    }
    if (!counts_available) {
        questions <- c(questions, paste0(
            "No 'counts' assay is present; available assays: ",
            if (length(assay_names)) paste(assay_names, collapse = ", ") else "none", "."
        ))
    }
    if (!de_package_available) {
        questions <- c(questions,
            "The 'limma' package is not installed, which the marker test relies on.")
    }

    ready <- group_resolvable && enough_groups && enough_cells && de_package_available
    list(
        status = if (ready) "ready_for_diagnostic" else "not_ready",
        checks = list(
            counts_assay_available = counts_available,
            available_assays = assay_names,
            group_column_resolvable = group_resolvable,
            group_column_resolved = group_key,
            n_groups = n_groups,
            smallest_group_n = smallest,
            detest_record_available = length(detest_ids) > 0L,
            available_detest_ids = detest_ids,
            de_package_available = de_package_available,
            raw_values_included = FALSE
        ),
        blocked_actions = if (ready) character() else "de",
        questions = questions,
        notes = paste(
            "No group column, clustering resolution or significance threshold is chosen here. Marker",
            "results are statistical associations under the selected test and thresholds, not measured",
            "quantities, and they are the input to enrichment rather than a second independent line of",
            "evidence for it. No cell is removed or filtered by this diagnostic."
        )
    )
}
