# Read-only readiness diagnostic and bounded lineage-driver summaries
#
# Everything in this file is deliberately privacy-withholding: the stored
# CellRank driver table carries gene symbols, and none of the functions below
# ever let a gene symbol, a per-cell value or a nominated "the" fate leave the
# package boundary. Only counts, the metric column NAME, anonymized route_N
# codes and aggregate statistics are emitted.

# Resolve the completed CellRank trajectory record used by every entry point in
# this file. Returns NULL when the id is unknown or the record is not CellRank.
sclet_ai_fate_drivers_record <- function(object, id = NULL) {
    records <- tryCatch(sclet_get_state_records(object, "trajectory"), error = function(e) list())
    if (is.null(records)) records <- list()
    if (is.null(id)) {
        id <- tryCatch(sclet_get_active_state(object, "trajectory"), error = function(e) NULL)
    }
    candidate <- NULL
    if (!is.null(id) && length(id) == 1L && nzchar(as.character(id))) {
        candidate <- records[[as.character(id)]]
    }
    if (is.null(candidate) && !is.null(id) && length(records) == 1L) {
        candidate <- records[[1L]]
    }
    if (is.null(candidate) || !is.list(candidate)) return(NULL)
    record_status <- as.character(candidate$status %||% "completed")
    if (!identical(record_status, "completed")) return(NULL)
    method <- as.character(candidate$method %||% "")
    analysis_key <- as.character(candidate$artifacts$analysis_key %||% "")
    is_cellrank <- identical(method, "CellRank") ||
        identical(analysis_key, "cellrank") ||
        identical(tolower(analysis_key), "cellrank")
    if (!is_cellrank) return(NULL)
    candidate
}

# The driver table lives in the state record artifacts. Fall back to the
# cellrank analysis record for objects written before artifacts carried it.
sclet_ai_fate_drivers_table <- function(object, record) {
    drivers <- record$artifacts$lineage_drivers %||% record$lineage_drivers
    if (is.null(drivers)) {
        analysis <- tryCatch(sclet_get_analysis(object, "cellrank"), error = function(e) NULL)
        drivers <- analysis$lineage_drivers
    }
    if (is.null(drivers)) return(NULL)
    if (!is.data.frame(drivers)) {
        drivers <- tryCatch(as.data.frame(drivers, stringsAsFactors = FALSE), error = function(e) NULL)
    }
    if (is.null(drivers) || !nrow(drivers)) return(NULL)
    drivers
}

# The user-visible metric column NAME (a plain identifier such as "corr"), or
# NULL when the table carries no numeric metric column at all.
sclet_ai_fate_drivers_metric <- function(drivers) {
    if (is.null(drivers)) return(NULL)
    tryCatch(sclet_pick_cellrank_driver_metric(drivers), error = function(e) NULL)
}

# Available lineage labels. Only used internally to build route_N codes; no
# lineage label is ever returned to the caller.
sclet_ai_fate_drivers_lineages <- function(record, drivers) {
    available <- NULL
    if (!is.null(drivers) && "lineage" %in% colnames(drivers)) {
        available <- unique(as.character(drivers$lineage))
    } else if (!is.null(record)) {
        available <- record$artifacts$fate_probability_names %||% record$fate_probability_names
        available <- as.character(available)
    }
    if (is.null(available)) available <- character()
    available[!is.na(available) & nzchar(available)]
}

# Deterministic anonymized code for a lineage label. Positions come from the
# available-lineage order, so the same label always yields the same route_N.
sclet_ai_fate_drivers_code <- function(lineage, available_lineages) {
    if (is.null(lineage) || is.null(available_lineages) || !length(available_lineages)) {
        return(NULL)
    }
    index <- match(as.character(lineage), as.character(available_lineages))
    if (is.na(index)) return(NULL)
    paste0("route_", index)
}

# The shared caveat block. Wordings are fixed here so no caller can drift into
# describing a driver as causal or a lineage as the fate of the dataset.
sclet_ai_fate_drivers_caveats <- function() {
    list(
        paste(
            "Lineage drivers are CellRank model estimates computed from a velocity kernel; they are derived from",
            "the fitted model, not measured from an experiment."
        ),
        paste(
            "The association reported here is with the model's absorption probabilities, so it is an association",
            "with a model output and NOT a causal relationship."
        ),
        paste(
            "The ranking is a modeling output, not a biological claim: a highly ranked driver is highly ranked",
            "under this model and this kernel only."
        ),
        paste(
            "Gene identity is deliberately withheld from this summary; only counts, the metric column name,",
            "anonymized route codes and aggregate statistics are returned."
        )
    )
}

# Bounded aggregate statistics of one lineage's ranked driver set. The bounded
# `driver_rank_codes` field name avoids the word "top", which would read as a
# nomination.
sclet_ai_fate_drivers_metric_summary <- function(absolute_values) {
    values <- as.numeric(absolute_values)
    values <- values[is.finite(values)]
    if (!length(values)) {
        return(list(
            n_observed = 0L,
            mean = NA_real_,
            quantile_25 = NA_real_,
            median = NA_real_,
            quantile_75 = NA_real_,
            max = NA_real_
        ))
    }
    quantiles <- round(unname(stats::quantile(values, c(0.25, 0.5, 0.75))), 4L)
    list(
        n_observed = as.integer(length(values)),
        mean = round(mean(values), 4L),
        quantile_25 = as.numeric(quantiles[[1L]]),
        median = as.numeric(quantiles[[2L]]),
        quantile_75 = as.numeric(quantiles[[3L]]),
        max = round(max(values), 4L)
    )
}

# Rank one lineage's rows by |metric| descending. A character `metric_value`
# (e.g. an unparsed "0.31" column) is coerced before ranking; gene identity is
# carried as a plain local vector that never enters a return value.
sclet_ai_fate_drivers_ranked <- function(subset_rows, metric) {
    values <- subset_rows[[metric]]
    if (!is.numeric(values)) {
        values <- suppressWarnings(as.numeric(as.character(values)))
    }
    absolute <- abs(values)
    order_key <- absolute
    order_key[!is.finite(order_key)] <- -Inf
    subset_rows <- subset_rows[order(-order_key, seq_len(nrow(subset_rows))), , drop = FALSE]
    list(
        rows = subset_rows,
        gene_ids = as.character(subset_rows$gene),
        absolute = abs(suppressWarnings(as.numeric(as.character(subset_rows[[metric]]))))
    )
}

# Jaccard overlap between the driver sets of two summarized lineages. Computed
# against gene identity internally and emitted as a single number, never as the
# shared identity itself.
sclet_ai_fate_drivers_jaccard <- function(left, right) {
    left <- unique(left[!is.na(left)])
    right <- unique(right[!is.na(right)])
    if (!length(left) || !length(right)) return(NA_real_)
    union_size <- length(union(left, right))
    if (!union_size) return(NA_real_)
    round(length(intersect(left, right)) / union_size, 4L)
}

#' Check whether an object is ready for lineage driver summary
#'
#' Reports whether a completed CellRank trajectory record with a lineage driver
#' table is present. This diagnostic deliberately never reports a gene name, a
#' driver ranking or a nominated fate: which lineages matter, and what a driver
#' score means biologically, is a modeling assumption only the user can make.
#' The result therefore contains only counts and column names.
#'
#' @param object A \code{SingleCellExperiment} object.
#' @param id Optional CellRank trajectory record id. When supplied, that exact
#'   completed record must be present.
#' @return A read-only readiness summary.
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' readiness <- check_fate_drivers_readiness(sce)
#' readiness$status
#' readiness$questions
#' @export
check_fate_drivers_readiness <- function(object, id = NULL) {
    if (!inherits(object, "SingleCellExperiment")) {
        stop("object must be a SingleCellExperiment", call. = FALSE)
    }
    records <- tryCatch(sclet_get_state_records(object, "trajectory"), error = function(e) list())
    if (is.null(records)) records <- list()
    requested_id <- if (is.null(id)) NULL else as.character(id)[[1L]]
    record <- sclet_ai_fate_drivers_record(object, id = id)
    record_available <- !is.null(record)
    resolved_id <- NULL
    if (record_available) {
        resolved_id <- as.character(record$id %||% requested_id %||% "")
        if (!nzchar(resolved_id)) resolved_id <- requested_id
    }

    drivers <- if (record_available) sclet_ai_fate_drivers_table(object, record) else NULL
    table_available <- !is.null(drivers) && is.data.frame(drivers) && nrow(drivers) > 0L
    has_lineage_column <- table_available && "lineage" %in% colnames(drivers)
    metric <- sclet_ai_fate_drivers_metric(drivers)
    available_lineages <- sclet_ai_fate_drivers_lineages(record, drivers)
    lineage_column <- if (has_lineage_column) as.character(drivers$lineage) else character()
    n_lineages <- as.integer(length(unique(lineage_column[!is.na(lineage_column) & nzchar(lineage_column)])))
    n_drivers <- if (table_available) as.integer(nrow(drivers)) else 0L

    # A declared-but-absent driver table is a real outcome: RunCellRank() only
    # computes drivers when a lineage driver plan is available, and it still
    # records the run as completed with has_lineage_drivers = FALSE.
    declared_available <- isTRUE(record$artifacts$has_lineage_drivers) || table_available
    driver_computation_ran <- record_available && declared_available && table_available

    questions <- character()
    if (!record_available) {
        available_trajectory_ids <- names(records)
        detail <- if (is.null(requested_id)) {
            if (length(available_trajectory_ids)) {
                paste0("Available trajectory records: ", paste(available_trajectory_ids, collapse = ", "), ".")
            } else {
                "No trajectory record is stored on this object."
            }
        } else {
            paste0(
                "Requested trajectory id '", requested_id, "' is not a completed CellRank record. ",
                if (length(available_trajectory_ids)) {
                    paste0("Available trajectory records: ", paste(available_trajectory_ids, collapse = ", "), ".")
                } else {
                    "No trajectory record is stored on this object."
                }
            )
        }
        questions <- c(questions, paste(
            "Lineage driver summary consumes a completed cell fate analysis; a completed",
            "RunCellFate() or RunCellRank() record is required, and none is available.", detail
        ))
    } else if (!driver_computation_ran) {
        questions <- c(questions, paste(
            "The CellRank record is present but its lineage driver computation did not run, so no driver table is",
            "available. Rerun the analysis with driver computation available for this record before summarizing drivers."
        ))
    }

    list(
        status = if (record_available) "ready_for_diagnostic" else "not_ready",
        checks = list(
            cellrank_record_available = record_available,
            cellrank_id_resolved = resolved_id,
            lineage_driver_table_available = table_available,
            driver_table_has_lineage_column = has_lineage_column,
            driver_metric_column = metric,
            n_lineages = n_lineages,
            n_drivers = n_drivers,
            raw_values_included = FALSE
        ),
        blocked_actions = if (record_available) character() else "fate_drivers",
        questions = questions,
        notes = paste(
            "No gene name, driver ranking or nominated fate is reported here: driver identity and the choice of",
            "lineage are modeling questions reserved for the user, and only counts and column names are returned.",
            "Lineage drivers are CellRank model estimates derived from the velocity kernel, not measured biology,",
            "and no cell is removed or filtered by this diagnostic."
        )
    )
}

#' Summarize the lineage drivers of one user-specified lineage
#'
#' Reads the stored CellRank lineage driver table for a single lineage named by
#' the user, ranks its drivers by the absolute value of the resolved driver
#' metric, and returns a bounded summary. Gene symbols exist in the underlying
#' table but never leave this function: only counts, the metric column name,
#' anonymized \code{route_N} codes and aggregate statistics of the absolute
#' metric are returned.
#'
#' @param object A \code{SingleCellExperiment} object.
#' @param id Optional CellRank trajectory record id.
#' @param lineage The lineage to summarize. When \code{NULL} and several
#'   lineages are available, no lineage is chosen for the user and the result
#'   asks for an explicit choice.
#' @param top_n Maximum number of drivers to summarize. Defaults to 10.
#' @return A read-only driver summary list with \code{status},
#'   \code{lineage_code}, \code{n_drivers_total}, \code{n_drivers_summarized},
#'   \code{driver_metric}, \code{driver_rank_codes},
#'   \code{driver_abs_metric_summary}, \code{claim_level}, \code{caveats} and
#'   \code{raw_values_included = FALSE}.
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' # No cell fate run is on record, so no lineage is chosen.
#' summary <- summarize_fate_lineage_drivers(sce)
#' summary$status
#' summary$reason
#' @export
summarize_fate_lineage_drivers <- function(object, id = NULL, lineage = NULL, top_n = 10L) {
    if (!inherits(object, "SingleCellExperiment")) {
        stop("object must be a SingleCellExperiment", call. = FALSE)
    }
    if (!is.null(lineage) && (!is.character(lineage) || length(lineage) != 1L || !nzchar(lineage))) {
        stop("lineage must be NULL or a single non-empty character value", call. = FALSE)
    }
    if (!is.numeric(top_n) || length(top_n) != 1L || !is.finite(top_n) || top_n < 1) {
        stop("top_n must be a single positive number", call. = FALSE)
    }
    top_n <- as.integer(top_n)

    record <- sclet_ai_fate_drivers_record(object, id = id)
    if (is.null(record)) {
        return(list(
            status = "not_available",
            reason = "no_cellrank_record",
            claim_level = "associated",
            caveats = sclet_ai_fate_drivers_caveats(),
            raw_values_included = FALSE
        ))
    }
    drivers <- sclet_ai_fate_drivers_table(object, record)
    metric <- sclet_ai_fate_drivers_metric(drivers)
    if (is.null(drivers) || is.null(metric)) {
        return(list(
            status = "not_available",
            reason = "no_lineage_driver_table",
            claim_level = "associated",
            caveats = sclet_ai_fate_drivers_caveats(),
            raw_values_included = FALSE
        ))
    }

    available_lineages <- sclet_ai_fate_drivers_lineages(record, drivers)
    resolved <- tryCatch(sclet_resolve_cellrank_lineage(
        list(lineage_drivers = drivers, fate_probability_names = available_lineages),
        lineage
    ), error = function(e) NULL)

    if (is.null(lineage) && !is.null(resolved) && length(available_lineages) > 1L) {
        # sclet_resolve_cellrank_lineage() falls back to the sole lineage only;
        # guard here so a multi-lineage object can never be summarized silently.
        resolved <- NULL
    }
    if (is.null(lineage) && is.null(resolved) && length(available_lineages) > 1L) {
        return(list(
            status = "clarification_required",
            n_available_lineages = as.integer(length(available_lineages)),
            question = paste0(
                "This object stores lineage drivers for ", length(available_lineages),
                " lineages; which one should be summarized? Name the lineage explicitly."
            ),
            claim_level = "associated",
            caveats = sclet_ai_fate_drivers_caveats(),
            raw_values_included = FALSE
        ))
    }
    if (is.null(resolved)) {
        return(list(
            status = "not_available",
            reason = "lineage_unknown",
            claim_level = "associated",
            caveats = sclet_ai_fate_drivers_caveats(),
            raw_values_included = FALSE
        ))
    }

    rows <- as.character(drivers$lineage)
    subset_rows <- drivers[!is.na(rows) & rows == as.character(resolved), , drop = FALSE]
    if (!nrow(subset_rows)) {
        return(list(
            status = "not_available",
            reason = "lineage_unknown",
            claim_level = "associated",
            caveats = sclet_ai_fate_drivers_caveats(),
            raw_values_included = FALSE
        ))
    }

    ranked <- sclet_ai_fate_drivers_ranked(subset_rows, metric)
    n_total <- as.integer(nrow(ranked$rows))
    n_summarized <- as.integer(min(top_n, n_total))
    keep <- seq_len(n_summarized)

    list(
        status = "available",
        lineage_code = sclet_ai_fate_drivers_code(resolved, available_lineages),
        n_drivers_total = n_total,
        n_drivers_summarized = n_summarized,
        driver_metric = as.character(metric),
        driver_rank_codes = paste0("route_", keep),
        driver_abs_metric_summary = sclet_ai_fate_drivers_metric_summary(ranked$absolute[keep]),
        claim_level = "associated",
        caveats = sclet_ai_fate_drivers_caveats(),
        raw_values_included = FALSE
    )
}

#' Compare lineage drivers across lineages
#'
#' Observational comparison of the stored CellRank driver table across at least
#' two lineages. The function summarizes how many drivers each anonymized
#' lineage carries, how large their absolute metric values are, and how much
#' the lineages' driver sets overlap. It never nominates, ranks or describes any
#' lineage as *the* fate, and no gene symbol, per-cell value or
#' *best lineage* field is emitted.
#'
#' @param object A \code{SingleCellExperiment} object.
#' @param id Optional CellRank trajectory record id.
#' @return A read-only comparison list with \code{status}, \code{n_lineages},
#'   \code{lineage_codes}, per-lineage summaries keyed by \code{route_N} codes,
#'   \code{distinctness}, \code{claim_level}, \code{caveats} and
#'   \code{raw_values_included = FALSE}.
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' comparison <- compare_fate_lineages(sce)
#' comparison$status
#' comparison$reason
#' @export
compare_fate_lineages <- function(object, id = NULL) {
    if (!inherits(object, "SingleCellExperiment")) {
        stop("object must be a SingleCellExperiment", call. = FALSE)
    }
    record <- sclet_ai_fate_drivers_record(object, id = id)
    drivers <- if (is.null(record)) NULL else sclet_ai_fate_drivers_table(object, record)
    metric <- sclet_ai_fate_drivers_metric(drivers)
    if (is.null(record) || is.null(drivers) || is.null(metric)) {
        return(list(
            status = "not_available",
            reason = "at_least_two_lineages_required",
            n_lineages = 0L,
            lineage_codes = character(),
            lineages = list(),
            distinctness = list(
                method = "jaccard_of_summarized_driver_sets",
                n_pairs = 0L,
                jaccard_min = NA_real_,
                jaccard_mean = NA_real_,
                jaccard_max = NA_real_
            ),
            claim_level = "associated",
            caveats = c(sclet_ai_fate_drivers_caveats(), paste(
                "At least two lineages with a driver table are required for an across-lineage comparison; none is",
                "selected or preferred when they are missing."
            )),
            recommendation = NULL,
            raw_values_included = FALSE
        ))
    }

    available_lineages <- sclet_ai_fate_drivers_lineages(record, drivers)
    rows <- as.character(drivers$lineage)
    lineage_labels <- unique(rows[!is.na(rows) & nzchar(rows)])
    lineage_labels <- available_lineages[available_lineages %in% lineage_labels]
    if (length(lineage_labels) < 2L) {
        return(list(
            status = "not_available",
            reason = "at_least_two_lineages_required",
            n_lineages = as.integer(length(lineage_labels)),
            lineage_codes = character(),
            lineages = list(),
            distinctness = list(
                method = "jaccard_of_summarized_driver_sets",
                n_pairs = 0L,
                jaccard_min = NA_real_,
                jaccard_mean = NA_real_,
                jaccard_max = NA_real_
            ),
            claim_level = "associated",
            caveats = c(sclet_ai_fate_drivers_caveats(), paste(
                "At least two lineages with a driver table are required for an across-lineage comparison; none is",
                "selected or preferred when they are missing."
            )),
            recommendation = NULL,
            raw_values_included = FALSE
        ))
    }

    top_n <- 10L
    per_lineage <- lapply(lineage_labels, function(label) {
        subset_rows <- drivers[rows == label, , drop = FALSE]
        ranked <- sclet_ai_fate_drivers_ranked(subset_rows, metric)
        keep <- seq_len(min(top_n, nrow(ranked$rows)))
        absolute <- ranked$absolute[keep]
        finite <- absolute[is.finite(absolute)]
        list(
            code = sclet_ai_fate_drivers_code(label, available_lineages),
            n_drivers = as.integer(nrow(ranked$rows)),
            n_drivers_summarized = as.integer(length(keep)),
            metric_mean = if (length(finite)) round(mean(finite), 4L) else NA_real_,
            metric_median = if (length(finite)) round(stats::median(finite), 4L) else NA_real_,
            metric_max = if (length(finite)) round(max(finite), 4L) else NA_real_,
            driver_ids = ranked$gene_ids[keep]
        )
    })
    codes <- vapply(per_lineage, function(entry) entry$code, character(1L))
    names(per_lineage) <- codes

    lineage_codes <- unname(codes)
    lineages <- lapply(per_lineage, function(entry) {
        list(
            n_drivers = entry$n_drivers,
            n_drivers_summarized = entry$n_drivers_summarized,
            metric_mean = entry$metric_mean,
            metric_median = entry$metric_median,
            metric_max = entry$metric_max,
            raw_values_included = FALSE
        )
    })
    names(lineages) <- lineage_codes

    pair_indices <- utils::combn(seq_along(per_lineage), 2L, simplify = FALSE)
    jaccard <- vapply(pair_indices, function(index_pair) {
        sclet_ai_fate_drivers_jaccard(
            per_lineage[[index_pair[[1L]]]]$driver_ids,
            per_lineage[[index_pair[[2L]]]]$driver_ids
        )
    }, numeric(1L))
    finite_jaccard <- jaccard[is.finite(jaccard)]

    list(
        status = "available",
        n_lineages = as.integer(length(per_lineage)),
        lineage_codes = lineage_codes,
        lineages = lineages,
        driver_metric = as.character(metric),
        distinctness = list(
            method = "jaccard_of_summarized_driver_sets",
            top_n = as.integer(top_n),
            n_pairs = as.integer(length(pair_indices)),
            jaccard_min = if (length(finite_jaccard)) round(min(finite_jaccard), 4L) else NA_real_,
            jaccard_mean = if (length(finite_jaccard)) round(mean(finite_jaccard), 4L) else NA_real_,
            jaccard_max = if (length(finite_jaccard)) round(max(finite_jaccard), 4L) else NA_real_,
            raw_values_included = FALSE
        ),
        claim_level = "associated",
        caveats = c(sclet_ai_fate_drivers_caveats(), paste(
            "This is an observational comparison of model outputs across lineages; no lineage is preferred, ranked",
            "or described as the fate of the dataset, and the comparison does not identify a target lineage."
        ), paste(
            "Overlapping driver sets are expected between related lineages, so a high overlap fraction is not",
            "evidence that the lineages are biologically identical."
        )),
        recommendation = NULL,
        raw_values_included = FALSE
    )
}
