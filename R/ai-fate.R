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
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' readiness <- check_fate_readiness(sce)
#' readiness$status
#' readiness$questions
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
sclet_ai_record_fate_evidence <- function(object, name, cluster_key, replace = FALSE) {
    # One predicate for "this is a usable CellRank result", shared with the
    # action handler: a duplicated copy drifted once already (it missed the
    # `status = "completed"` requirement and the lowercased analysis key, so a
    # lowercase-cased record could be recorded here and rejected there).
    resolved <- tryCatch(
        sclet_ai_resolve_cellrank_record(object, cellrank_id = name),
        error = function(e) list(id = NULL, record = NULL, message = conditionMessage(e))
    )
    record <- resolved$record
    if (is.null(resolved$id) || is.null(record)) return(object)
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
    # `terminal_state_fractions` is the share of the TERMINAL-STATE CELLS per
    # state, not the share of cells that will end there. CellRank labels a fixed
    # small set of cells as terminal, so on real data that share is close to
    # 1/k by construction (measured 2026-10-20 on the pancreas dataset: 30 / 30
    # / 30 of 3,696 cells, i.e. 0.333 each). A reader will read it as "one third
    # of the cells go to each fate", which is wrong. The honest destination
    # distribution is the argmax over the absorption probabilities across ALL
    # cells, recorded alongside it, plus the per-dimension means.
    n_terminal_state_cells <- if (length(terminal_state_column) == 1L &&
            nzchar(terminal_state_column) &&
            terminal_state_column %in% colnames(coldata)) {
        as.integer(sum(!is.na(coldata[[terminal_state_column]])))
    } else {
        0L
    }
    destination_fractions <- numeric()
    fate_dimension_means <- numeric()
    if (length(fate_summary)) {
        fate_dimension_means <- round(vapply(fate_summary, function(x) as.numeric(x$mean),
            numeric(1L)), 4L)
    }
    if (length(fate_columns) >= 2L) {
        mat <- do.call(cbind, lapply(fate_columns, function(cn) {
            as.numeric(coldata[[cn]])
        }))
        complete <- stats::complete.cases(mat)
        if (sum(complete) > 0L) {
            dest <- max.col(mat[complete, , drop = FALSE], ties.method = "first")
            destination_fractions <- round(as.numeric(tabulate(dest,
                nbins = length(fate_columns))) / sum(complete), 4L)
        }
    }
    values <- list(
        n_terminal_states = n_terminal_states,
        n_fate_dimensions = n_fate_dimensions,
        n_lineage_drivers = n_lineage_drivers,
        terminal_state_codes = terminal_state_codes,
        terminal_state_fractions = terminal_state_fractions,
        n_terminal_state_cells = as.integer(n_terminal_state_cells),
        # argmax over the absorption probabilities: where the cells actually go
        destination_fractions = destination_fractions,
        destination_fractions_are_argmax = TRUE,
        fate_dimension_mean_absorption = fate_dimension_means,
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
        id = sclet_ai_fate_evidence_id(name),
        kind = "deterministic_summary",
        values = values,
        claim_level = "consistent_with"
    )
    sclet_ai_record_evidence_for_run(object, evidence, source = name,
        parents = character(), scope = NULL, replace = replace)
}

# Sanitize a record id for use inside an evidence id.
#
# A record id can be a user-supplied string and therefore may carry a path, a
# file name or a private sample name. The evidence id must be stable and must
# never echo that text back, so anything outside [A-Za-z0-9_-] collapses to "_"
# and the result is truncated. Collisions are harmless here: an identical
# sanitized id means an identical record id, which is the idempotent case.
sclet_ai_evidence_id_token <- function(value) {
    token <- as.character(value %||% "")
    token <- token[[1L]]
    if (is.na(token) || !nzchar(token)) return("unknown")
    token <- gsub("[^A-Za-z0-9_-]+", "_", token)
    token <- gsub("_+", "_", token)
    token <- gsub("^_+|_+$", "", token)
    if (!nzchar(token)) return("unknown")
    substr(token, 1L, 48L)
}

# The one evidence id a fate record's bounded summary is stored under.
#
# Both the recorder and the action handler ask this function, so the id the
# handler verifies is by construction the id that was written. They disagreed
# once: the recorder wrote the raw record id while the handler looked up a
# sanitized one, so a record named e.g. "fate.1" recorded successfully and then
# failed the action with a spurious fate_evidence_missing.
sclet_ai_fate_evidence_id <- function(name) {
    paste0("ev:fate_", sclet_ai_evidence_id_token(name))
}

# Record the user's own explicit terminal-state confirmation for a CellRank run.
#
# This records a USER DECISION, not a finding about biology. The node's
# claim_level is therefore "observed" -- a statement of what the user chose --
# and never "consistent_with" or "associated", which would assert something
# about the data. The only thing entering `values` is the count of confirmed
# states, the anonymized state codes the user named, and boolean flags; the
# cellrank record id, the colData column text and any free-text rationale are
# never stored. Mirroring sclet_ai_record_clarification_response(), free text is
# deliberately absent so the payload stays inside sclet_ai_evidence_value_ok()
# without a special case.
#
# This function never picks, suggests, ranks or guesses a terminal state: it
# records exactly the codes it is handed, after dropping any code the stored
# result does not actually carry. An empty confirmation records nothing rather
# than inventing a state.
sclet_ai_record_fate_terminal_state_confirmation <- function(object, cellrank_id, terminal_state_codes) {
    if (!inherits(object, "SingleCellExperiment")) {
        stop("object must be a SingleCellExperiment", call. = FALSE)
    }
    # Resolve defensively: a confirmation that cannot be tied to a completed
    # CellRank record leaves the object untouched rather than writing an
    # orphaned decision into the ledger.
    record <- tryCatch(
        sclet_ai_resolve_cellrank_record(object, cellrank_id),
        error = function(e) list(id = NULL, record = NULL, message = "unresolved")
    )
    if (is.null(record$id) || is.null(record$record)) return(object)
    # The resolved id names a completed record; it is safe to use as `source`.
    source_id <- record$id
    codes <- sclet_ai_vector_param(terminal_state_codes, "character")
    codes <- codes[!is.na(codes)]
    codes <- trimws(codes)
    codes <- codes[nzchar(codes)]
    # Only anonymized codes may ever enter evidence, and only codes the stored
    # result really carries. Anything else is dropped -- never repaired, never
    # guessed, never replaced by a similar-looking state.
    codes <- codes[grepl("^(group|cluster|condition|sample|route)_[0-9]+$", codes)]
    valid <- sclet_ai_terminal_state_codes(object, record$record)
    codes <- unique(codes[codes %in% valid])
    if (!length(codes)) return(object)
    values <- list(
        n_confirmed_terminal_states = as.integer(length(codes)),
        confirmed_terminal_state_codes = codes,
        cellrank_record_confirmed = TRUE,
        raw_values_included = FALSE,
        free_text_recorded = FALSE
    )
    evidence <- list(
        id = paste0("ev:fate_terminal_", sclet_ai_evidence_id_token(source_id)),
        kind = "user_decision",
        values = values,
        claim_level = "observed"
    )
    sclet_ai_record_evidence_for_run(object, evidence, source = source_id,
        parents = character(), scope = NULL)
}
