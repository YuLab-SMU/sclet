sclet_ai_trajectory_root_run_records <- function(object, requested = NULL) {
    records <- sclet_get_state_records(object, "trajectory")
    if (!length(records)) return(list())
    keep <- vapply(records, function(record) {
        if (!is.list(record)) return(FALSE)
        record_status <- as.character(record$status %||% "")
        method <- as.character(record$method %||% "")
        analysis_key <- as.character(record$artifacts$analysis_key %||% "")
        identical(record_status, "completed") &&
            (identical(method, "slingshot") || identical(analysis_key, "trajectory"))
    }, logical(1L))
    records <- records[keep]
    if (!is.null(requested)) {
        requested <- unique(as.character(requested))
        records <- records[names(records) %in% requested]
    }
    records
}

sclet_ai_trajectory_root_node_is_valid <- function(node, run_ids) {
    if (!is.list(node)) return(FALSE)
    if (!identical(as.character(node$kind %||% ""), "deterministic_summary")) return(FALSE)
    source <- as.character(node$source %||% "")
    if (!source %in% run_ids) return(FALSE)
    values <- node$values %||% list()
    root <- as.character(values$start_group %||% "")
    is.list(values) && nzchar(root) &&
        length(values$n_lineages %||% numeric()) == 1L &&
        length(values$raw_values_included %||% logical()) == 1L
}

sclet_ai_trajectory_root_spread_ratio <- function(a, b) {
    if (!is.finite(a) || !is.finite(b)) return(NA_real_)
    lower <- min(a, b)
    upper <- max(a, b)
    if (!is.finite(lower) || lower <= 0) return(NA_real_)
    upper / lower
}

#' Compare completed trajectory runs across explicit root hypotheses
#'
#' This function reads only completed slingshot trajectory state records and
#' their bounded `ev:trajectory_<name>` evidence nodes. It compares multiple
#' trajectory runs generated from different explicit `start_cluster`
#' assumptions, summarizes how many runs support each anonymized root code, and
#' reports pairwise differences in lineage count and relative pseudotime spread.
#' It never reads per-cell pseudotime vectors, never chooses a preferred root,
#' and never writes a new evidence node.
#'
#' @param object A `SingleCellExperiment` object.
#' @param ids Optional completed trajectory analysis ids to compare. When
#'   `NULL`, all completed slingshot trajectory runs with bounded evidence are
#'   used.
#' @return A bounded list with `status`, per-run summaries, root-hypothesis
#'   summaries, pairwise comparisons, and explicit caveats. Fewer than two
#'   usable runs return `status = "not_available"` with a typed `reason`. The
#'   result contains no per-cell values, always sets `recommendation = NULL`,
#'   and always sets `raw_values_included = FALSE`.
#' @export
compare_trajectory_roots <- function(object, ids = NULL) {
    if (!inherits(object, "SingleCellExperiment")) {
        stop("object must be a SingleCellExperiment", call. = FALSE)
    }
    if (!is.null(ids) && (!is.character(ids) || any(!nzchar(ids)))) {
        stop("ids must be NULL or a non-empty character vector", call. = FALSE)
    }

    all_records <- sclet_ai_trajectory_root_run_records(object)
    selected_records <- sclet_ai_trajectory_root_run_records(object, requested = ids)
    if (!length(all_records)) {
        return(list(
            status = "not_available",
            reason = "no_trajectory_runs",
            n_runs = 0L,
            runs = list(),
            roots = list(),
            comparisons = list(),
            recommendation = NULL,
            raw_values_included = FALSE
        ))
    }
    if (!length(selected_records)) {
        return(list(
            status = "not_available",
            reason = "no_requested_trajectory_runs",
            n_runs = 0L,
            runs = list(),
            roots = list(),
            comparisons = list(),
            recommendation = NULL,
            raw_values_included = FALSE
        ))
    }

    evidence <- sclet_ai_evidence_get_all(object)
    evidence <- evidence[vapply(evidence, sclet_ai_trajectory_root_node_is_valid,
        logical(1L), run_ids = names(selected_records))]
    usable_runs <- names(selected_records)[vapply(names(selected_records), function(run_id) {
        any(vapply(evidence, function(node) identical(as.character(node$source %||% ""), run_id), logical(1L)))
    }, logical(1L))]

    if (length(usable_runs) < 2L) {
        return(list(
            status = "not_available",
            reason = if (!length(evidence)) {
                "no_trajectory_evidence_records"
            } else {
                "at_least_two_completed_trajectory_runs_required"
            },
            n_runs = as.integer(length(usable_runs)),
            runs = list(),
            roots = list(),
            comparisons = list(),
            recommendation = NULL,
            raw_values_included = FALSE
        ))
    }

    by_run <- lapply(usable_runs, function(run_id) {
        record <- selected_records[[run_id]]
        node <- evidence[[which(vapply(evidence, function(candidate) {
            identical(as.character(candidate$source %||% ""), run_id)
        }, logical(1L)))[[1L]]]]
        values <- node$values %||% list()
        list(
            run_id = run_id,
            root_code = as.character(values$start_group %||% "cluster_unknown"),
            start_group_code = as.integer(values$start_group_code %||% NA_integer_),
            n_lineages = as.integer(values$n_lineages %||% NA_integer_),
            n_groups = as.integer(values$n_groups %||% NA_integer_),
            n_summarized_lineages = as.integer(values$n_summarized_lineages %||% 0L),
            first_lineage_median = as.numeric(values$first_lineage_median %||% NA_real_),
            first_lineage_iqr = as.numeric(values$first_lineage_iqr %||% NA_real_),
            first_lineage_max = as.numeric(values$first_lineage_max %||% NA_real_),
            relative_spread = sclet_ai_tv_spread(
                values$first_lineage_iqr %||% NA_real_,
                values$first_lineage_max %||% NA_real_
            ),
            reduction = as.character(record$inputs$reduction %||% record$reduction %||% ""),
            group_column = as.character(record$inputs$group %||% record$cluster_key %||% "")
        )
    })
    names(by_run) <- usable_runs

    run_summaries <- lapply(by_run, function(run) {
        list(
            root_code = run$root_code,
            start_group_code = run$start_group_code,
            n_lineages = run$n_lineages,
            n_groups = run$n_groups,
            n_summarized_lineages = run$n_summarized_lineages,
            first_lineage_median = run$first_lineage_median,
            first_lineage_iqr = run$first_lineage_iqr,
            first_lineage_max = run$first_lineage_max,
            relative_spread = run$relative_spread,
            reduction = run$reduction,
            raw_values_included = FALSE
        )
    })

    root_codes <- sort(unique(vapply(by_run, function(run) run$root_code, character(1L))))
    root_summaries <- lapply(root_codes, function(root_code) {
        matching_runs <- by_run[vapply(by_run, function(run) identical(run$root_code, root_code), logical(1L))]
        lineage_counts <- vapply(matching_runs, function(run) run$n_lineages, integer(1L))
        spreads <- vapply(matching_runs, function(run) run$relative_spread, numeric(1L))
        finite_spreads <- spreads[is.finite(spreads)]
        list(
            root_code = root_code,
            run_count = as.integer(length(matching_runs)),
            run_fraction = round(length(matching_runs) / length(by_run), 4L),
            runs = names(matching_runs),
            n_lineages_min = as.integer(min(lineage_counts)),
            n_lineages_max = as.integer(max(lineage_counts)),
            n_lineages_mean = round(mean(lineage_counts), 4L),
            relative_spread_min = if (length(finite_spreads)) round(min(finite_spreads), 6L) else NA_real_,
            relative_spread_max = if (length(finite_spreads)) round(max(finite_spreads), 6L) else NA_real_,
            raw_values_included = FALSE
        )
    })
    names(root_summaries) <- root_codes

    pair_indices <- utils::combn(seq_along(by_run), 2L, simplify = FALSE)
    comparisons <- lapply(pair_indices, function(index_pair) {
        left <- by_run[[index_pair[[1L]]]]
        right <- by_run[[index_pair[[2L]]]]
        same_root <- identical(left$root_code, right$root_code)
        ratio <- sclet_ai_trajectory_root_spread_ratio(left$relative_spread, right$relative_spread)
        list(
            run_a = left$run_id,
            run_b = right$run_id,
            same_root = same_root,
            root_difference = if (same_root) "same_recorded_root_code" else "different_recorded_root_code",
            same_reduction = identical(left$reduction, right$reduction),
            same_group_column = identical(left$group_column, right$group_column),
            lineage_delta = as.integer(abs(left$n_lineages - right$n_lineages)),
            relative_spread_ratio = ratio,
            comparable_spread = is.finite(ratio),
            interpretation = if (same_root) {
                "The two runs use the same recorded root code; any differences summarize non-root parameters or bounded algorithmic variation, not a preferred root."
            } else {
                "The two runs use different recorded root codes; differences summarize alternative biological root assumptions and do not select a preferred root."
            },
            raw_values_included = FALSE
        )
    })
    names(comparisons) <- paste0("comparison_", seq_along(comparisons))

    list(
        status = "available",
        n_runs = as.integer(length(by_run)),
        runs = run_summaries,
        roots = root_summaries,
        n_roots = as.integer(length(root_summaries)),
        n_distinct_root_hypotheses = as.integer(length(root_summaries)),
        comparisons = comparisons,
        n_comparisons = as.integer(length(comparisons)),
        claim_level = "consistent_with",
        caveats = list(
            "Trajectory root is an explicit biological assumption supplied by the user; this function summarizes alternative completed runs but does not choose a preferred root.",
            "Comparisons use only bounded run-level summaries (lineage counts and relative pseudotime spread) and do not establish per-cell agreement, temporal truth, or biological correctness."
        ),
        independence = list(
            n_runs = as.integer(length(by_run)),
            n_root_hypotheses = as.integer(length(root_summaries)),
            note = "Each completed trajectory run is one explicit root-hypothesis line, not an independent biological replicate."
        ),
        matching = list(
            method = "same_object_completed_trajectory_runs",
            root_labels_anonymized = TRUE,
            cell_level_overlap_available = FALSE,
            interpretation = "Runs are compared only through recorded bounded summaries; no per-cell pseudotime vectors or cell-level overlap are read."
        ),
        recommendation = NULL,
        raw_values_included = FALSE
    )
}
