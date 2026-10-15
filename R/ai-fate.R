# Readiness diagnostic and bounded evidence recorder for cell fate analysis

#' Check whether an object is ready for cell fate analysis
#'
#' Reports whether a completed velocity record, a usable embedding and a
#' cluster assignment are present for \code{RunCellFate()}. This diagnostic
#' deliberately never proposes, ranks or scores terminal states: which states
#' count as terminal, and what a fate probability means biologically, is a
#' modeling assumption only the user can make. The result therefore contains no
#' candidate terminal state, no ordering of clusters, and no fate ranking.
#'
#' @param object A \code{SingleCellExperiment} object.
#' @param reduction Optional preferred reduction name. Defaults to UMAP when
#'   available, then any usable embedding.
#' @param velocity_id Optional velocity record id. When supplied, that exact
#'   completed velocity record must be present.
#' @return A read-only readiness summary.
#' @export
check_fate_readiness <- function(object, reduction = NULL, velocity_id = NULL) {
    if (!inherits(object, "SingleCellExperiment")) {
        stop("object must be a SingleCellExperiment", call. = FALSE)
    }
    columns <- colnames(SummarizedExperiment::colData(object))
    velocity_records <- tryCatch(sclet_get_state_records(object, "velocity"), error = function(e) list())
    if (is.null(velocity_records)) velocity_records <- list()
    is_completed <- vapply(velocity_records, function(record) {
        identical(as.character(record$status %||% "completed"), "completed")
    }, logical(1L))
    completed_ids <- names(velocity_records)
    if (length(is_completed) && length(completed_ids) == length(is_completed)) {
        completed_ids <- completed_ids[is_completed]
    } else {
        completed_ids <- character()
    }
    velocity_record_available <- length(completed_ids) > 0L
    if (!is.null(velocity_id)) {
        velocity_record_available <- velocity_record_available &&
            as.character(velocity_id) %in% completed_ids
    }
    resolved_velocity_id <- NULL
    if (velocity_record_available) {
        if (is.null(velocity_id)) {
            active_id <- tryCatch(sclet_get_active_state(object, "velocity"), error = function(e) NULL)
            if (!is.null(active_id) && as.character(active_id) %in% completed_ids) {
                resolved_velocity_id <- as.character(active_id)
            } else {
                resolved_velocity_id <- completed_ids[[1L]]
            }
        } else {
            resolved_velocity_id <- as.character(velocity_id)
        }
    }
    assay_names <- SummarizedExperiment::assayNames(object)
    has_spliced_unspliced <- all(c("spliced", "unspliced") %in% assay_names)
    has_velocity_assay <- any(grepl("velocity", assay_names, ignore.case = TRUE))
    spliced_unspliced_or_velocity_assay_available <- has_spliced_unspliced || has_velocity_assay
    reductions <- tryCatch(SingleCellExperiment::reducedDimNames(object), error = function(e) character())
    resolved_reduction <- sclet_ai_diag_embedding_name(object, reduction)
    has_reduction <- !is.null(resolved_reduction)
    cluster_key <- tryCatch(ActiveIdent(object), error = function(e) NULL)
    if (!is.null(cluster_key) && identical(as.character(cluster_key), "colLabels")) {
        if (is.null(tryCatch(SingleCellExperiment::colLabels(object), error = function(e) NULL))) {
            cluster_key <- NULL
        }
    }
    if (!is.null(cluster_key) && !as.character(cluster_key) %in% columns &&
            !identical(as.character(cluster_key), "colLabels")) {
        cluster_key <- NULL
    }
    if (is.null(cluster_key)) {
        candidates <- grep("cluster|leiden|louvain|ident", columns, ignore.case = TRUE, value = TRUE)
        if (length(candidates)) cluster_key <- candidates[[1L]]
    }
    cluster_key_resolvable <- !is.null(cluster_key)
    questions <- character()
    if (!velocity_record_available) {
        detail <- if (is.null(velocity_id)) {
            paste0("Available velocity records: ", paste(completed_ids, collapse = ", "))
        } else {
            paste0("Requested velocity id '", as.character(velocity_id), "' is not a completed velocity record. ",
                "Available velocity records: ", paste(completed_ids, collapse = ", "))
        }
        questions <- c(questions, paste0(
            "Cell fate analysis consumes a completed RNA velocity estimate; no completed velocity record is available. ",
            detail
        ))
    }
    if (!spliced_unspliced_or_velocity_assay_available) {
        questions <- c(questions, paste0(
            "Cell fate analysis needs the spliced/unspliced counts or a stored velocity assay; the current assays are: ",
            paste(assay_names, collapse = ", ")
        ))
    }
    if (!has_reduction) {
        questions <- c(questions, paste0(
            "No usable embedding is available; run the dimensionality reduction first. ",
            "Available reductions: ", paste(reductions, collapse = ", ")
        ))
    }
    if (!cluster_key_resolvable) {
        questions <- c(questions, paste0(
            "No cluster assignment is available; cell fate analysis needs an explicit cluster column. ",
            "Available colData columns: ", paste(columns, collapse = ", ")
        ))
    }
    status <- if (!velocity_record_available || !has_reduction || !cluster_key_resolvable) {
        "not_ready"
    } else {
        "ready_for_diagnostic"
    }
    list(
        status = status,
        checks = list(
            velocity_record_available = velocity_record_available,
            velocity_id_resolved = resolved_velocity_id,
            spliced_unspliced_or_velocity_assay_available = spliced_unspliced_or_velocity_assay_available,
            reduction_resolved = resolved_reduction,
            cluster_key_resolvable = cluster_key_resolvable,
            cluster_key_resolved = cluster_key,
            available_reductions = reductions,
            available_velocity_ids = completed_ids,
            raw_values_included = FALSE
        ),
        blocked_actions = if (identical(status, "not_ready")) "fate" else character(),
        questions = questions,
        notes = paste(
            "No terminal state, fate or lineage is suggested, ranked or scored here: choosing which states count as",
            "terminal is a modeling assumption reserved for the user. Fate (absorption) probabilities are CellRank model",
            "estimates derived from the velocity field, not measured cell fate, and no cell is removed or filtered by",
            "this diagnostic."
        )
    )
}

# Record a bounded aggregate evidence node for a CellRank fate analysis.
# Only aggregate counts and quantiles enter the node; no per-cell vector is
# ever read into evidence, and no terminal state is ranked or recommended.
sclet_ai_record_fate_evidence <- function(object, name, cluster_key) {
    record <- tryCatch(sclet_get_state_record(object, "trajectory", name), error = function(e) NULL)
    if (is.null(record)) return(object)
    method <- as.character(record$method %||% "")
    analysis_key <- as.character(record$artifacts$analysis_key %||% "")
    is_cellrank <- identical(method, "CellRank") ||
        identical(analysis_key, "cellrank") ||
        identical(tolower(analysis_key), "cellrank")
    if (!is_cellrank) return(object)
    coldata <- SummarizedExperiment::colData(object)
    summary <- record$summary %||% list()
    n_terminal_states <- as.integer(summary$n_terminal_states %||% 0L)
    n_fate_dimensions <- as.integer(summary$n_fate_dimensions %||% 0L)
    n_lineage_drivers <- as.integer(summary$n_lineage_drivers %||% 0L)
    cluster_levels <- NULL
    if (!is.null(cluster_key) && is.character(cluster_key) && length(cluster_key) == 1L &&
            cluster_key %in% colnames(coldata)) {
        cluster_levels <- levels(as.factor(as.character(coldata[[cluster_key]])))
    }
    terminal_state_codes <- character()
    terminal_state_fractions <- numeric()
    terminal_state_column <- as.character(record$artifacts$terminal_state_col %||% "")
    if (length(terminal_state_column) == 1L && nzchar(terminal_state_column) &&
            terminal_state_column %in% colnames(coldata)) {
        labels <- as.character(coldata[[terminal_state_column]])
        labels <- labels[!is.na(labels)]
        labels <- labels[nzchar(labels)]
        if (length(labels)) {
            counts <- stats::na.omit(table(labels))
            if (length(counts)) {
                total <- sum(counts)
                labels_seen <- names(counts)
                index <- match(labels_seen, cluster_levels)
                if (is.null(cluster_levels) || any(is.na(index))) {
                    fallback <- seq_along(labels_seen)
                    index[is.na(index)] <- fallback[is.na(index)]
                }
                terminal_state_codes <- paste0("cluster_", index)
                terminal_state_fractions <- round(as.numeric(counts) / total, 4L)
            }
        }
    }
    fate_summary <- list()
    fate_columns <- as.character(record$artifacts$fate_probability_cols %||% character())
    fate_columns <- intersect(fate_columns, colnames(coldata))
    for (column in fate_columns) {
        values <- as.numeric(coldata[[column]])
        values <- values[is.finite(values)]
        if (!length(values)) next
        quantiles <- round(unname(stats::quantile(values, c(0.25, 0.5, 0.75))), 4L)
        fate_summary[[column]] <- list(
            n_observed = as.integer(length(values)),
            mean = round(mean(values), 4L),
            quantile_25 = as.numeric(quantiles[[1L]]),
            median = as.numeric(quantiles[[2L]]),
            quantile_75 = as.numeric(quantiles[[3L]]),
            max = round(max(values), 4L)
        )
    }
    values <- list(
        n_terminal_states = n_terminal_states,
        n_fate_dimensions = n_fate_dimensions,
        n_lineage_drivers = n_lineage_drivers,
        terminal_state_codes = terminal_state_codes,
        terminal_state_fractions = terminal_state_fractions,
        absorption_probability_is_model_estimate = TRUE,
        pseudotime_is_absolute_time = FALSE,
        raw_values_included = FALSE
    )
    if (length(fate_summary)) {
        first <- fate_summary[[1L]]
        values$n_summarized_fate_dimensions <- as.integer(length(fate_summary))
        values$first_fate_median <- first$median
        values$first_fate_iqr <- round(first$quantile_75 - first$quantile_25, 4L)
        values$first_fate_max <- first$max
    }
    evidence <- list(
        id = paste0("ev:fate_", name),
        kind = "deterministic_summary",
        values = values,
        claim_level = "consistent_with"
    )
    tryCatch(
        RecordAIEvidence(object, evidence, source = name, parents = character(), scope = NULL),
        error = function(e) object
    )
}
