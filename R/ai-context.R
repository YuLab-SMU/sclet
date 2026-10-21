#' Build a stable AI-facing analysis ledger
#'
#' `GetAnalysisLedger()` returns a bounded, versioned view of the current
#' `SingleCellExperiment` analysis state. It includes deterministic dataset
#' metadata, the active view, health flags from [Status()], analysis and state
#' records, and conservative read-only capability information. Expression
#' matrices and other large objects are represented by size/type references;
#' they are never copied into the ledger.
#'
#' @param object A `SingleCellExperiment` object.
#' @param detail Detail level. One of `"summary"` or `"full"`.
#' @param target Optional context target(s), such as `"status"`, `"qc"`,
#'   `"planning"`, `"markers"`, `"trajectory"`, `"velocity"`,
#'   `"perturbation"`, or `"report"`.
#' @param include_artifacts Logical. Include bounded artifact references in
#'   analysis records. Defaults to `FALSE`, which also hides the `findings` a
#'   recorded AI result keeps under `artifacts` -- read a result back with
#'   `include_artifacts = TRUE`, or its findings come back empty rather than
#'   absent (a silent-empty trap measured 2026-10-20).
#' @param include_data Logical. Include bounded metadata previews and assay
#'   dimensions. Complete expression matrices are never included. Defaults to
#'   `FALSE`.
#' @return A versioned, bounded AI context list.
#' @export
GetAnalysisLedger <- function(
    object,
    detail = c("summary", "full"),
    target = NULL,
    include_artifacts = FALSE,
    include_data = FALSE
) {
    detail <- match.arg(detail)
    sclet_ai_build_context(
        object = object,
        target = target,
        detail = detail,
        include_artifacts = isTRUE(include_artifacts),
        include_data = isTRUE(include_data)
    )
}

sclet_ai_context <- function(
    object,
    target = c(
        "status", "qc", "planning", "markers", "trajectory",
        "velocity", "perturbation", "report"
    ),
    detail = c("summary", "full"),
    registry = NULL
) {
    detail <- match.arg(detail)
    sclet_ai_build_context(
        object = object,
        target = target,
        detail = detail,
        registry = registry,
        include_artifacts = FALSE,
        include_data = FALSE
    )
}

sclet_ai_build_context <- function(
    object,
    target = NULL,
    detail = c("summary", "full"),
    registry = NULL,
    include_artifacts = FALSE,
    include_data = FALSE
) {
    detail <- match.arg(detail)
    status <- Status(object)
    state <- sclet_get_state(object)
    target <- sclet_ai_normalize_target(target)
    records <- sclet_ai_collect_records(
        object = object,
        state = state,
        detail = detail,
        target = target,
        include_artifacts = include_artifacts,
        include_data = include_data
    )

    dataset <- sclet_ai_dataset_summary(object, include_data = include_data)
    active_view <- status$active
    if (is.null(active_view)) {
        active_view <- list()
    }
    active_states <- tryCatch(
        sclet_get_active_state(object),
        error = function(e) list()
    )
    health <- status$health
    if (is.null(health)) {
        health <- list()
    }

    result <- list(
        schema_version = "1.0",
        # Which build produced this ledger. The package version alone is not
        # enough: during 2026-10-19/20 the installed tree sat two days behind
        # the source at the SAME version number, and only a schema difference
        # revealed it. `sclet_version` plus the install timestamp makes "which
        # code am I talking to" answerable from the ledger itself.
        sclet_version = tryCatch(
            as.character(utils::packageVersion("sclet")),
            error = function(e) NA_character_
        ),
        sclet_installed_at = tryCatch({
            path <- system.file("DESCRIPTION", package = "sclet")
            if (nzchar(path)) format(file.info(path)$mtime, "%Y-%m-%d %H:%M:%S") else NA_character_
        }, error = function(e) NA_character_),
        dataset = dataset,
        active_view = active_view,
        active_states = active_states,
        # The HVG selection is a result-shaping choice (`scran` and `seurat`
        # pick different genes) and it is stored only in the raw state, where
        # the ledger never looked: on PBMC 10k (2026-10-19) features$hvg held
        # n = 2000 / method = "scran" and nothing in the context said so, so a
        # next-round planner could not see, or reproduce, the choice in force.
        features = if (is.null(state$features)) {
            list()
        } else {
            sclet_ai_safe_value(state$features, include_data = include_data)
        },
        health = health,
        analyses = records$analyses,
        state_records = records$state_records,
        workflows = records$workflows,
        lineage = records$lineage,
        analysis_story = sclet_ai_build_analysis_story(object, records, status),
        capabilities = sclet_ai_capabilities(status, records, object = object),
        blocked_actions = sclet_ai_blocked_actions(health),
        quality_checks = list(
            health = health,
            mainline_missing = sclet_ai_or(health$mainline_missing, character())
        ),
        available_actions = sclet_ai_describe_registry(registry),
        warnings = records$warnings
    )
    result$fingerprint <- sclet_ai_fingerprint(result)
    result
}

sclet_ai_normalize_target <- function(target) {
    if (is.null(target)) {
        return(NULL)
    }
    target <- unique(tolower(trimws(as.character(target))))
    target[nzchar(target)]
}

sclet_ai_target_matches <- function(type, id, target) {
    if (is.null(target) || !length(target)) {
        return(TRUE)
    }
    known <- c(
        "status", "qc", "planning", "markers", "trajectory", "velocity",
        "perturbation", "report"
    )
    if (all(known %in% target) || "report" %in% target) {
        return(TRUE)
    }
    value <- tolower(paste(c(type, id), collapse = " "))
    aliases <- list(
        status = c("ai_evidence"),
        qc = c("qc", "preprocess", "clean", "doublet", "phenotype", "ai_evidence"),
        planning = c("workflow", "plan"),
        markers = c("marker", "enrichment", "detest", " da", "da ", "gene_set", "geneset"),
        trajectory = c("trajectory", "cellrank"),
        velocity = c("velocity", "regvelo"),
        perturbation = c("perturbation", "priority", "rare", "augur", "state_priority")
    )
    any(vapply(target, function(item) {
        if (item %in% names(aliases)) {
            if (!length(aliases[[item]])) {
                return(FALSE)
            }
            return(any(vapply(aliases[[item]], grepl, logical(1), x = value, fixed = TRUE)))
        }
        grepl(item, value, fixed = TRUE)
    }, logical(1)))
}

sclet_ai_evidence_record_is_completed <- function(record) {
    if (!is.list(record)) return(FALSE)
    record_status <- record$status
    if (is.character(record_status) && length(record_status) == 1L && record_status == "completed") {
        return(TRUE)
    }
    summary <- record$summary
    if (is.list(summary)) {
        summary_status <- summary$status
        if (is.character(summary_status) && length(summary_status) == 1L && summary_status == "completed") {
            return(TRUE)
        }
    }
    FALSE
}

sclet_ai_strip_evidence_for_interpretive <- function(normalized) {
    if (!is.list(normalized)) return(normalized)
    if (is.list(normalized$summary)) {
        keep <- c("id", "kind", "claim_level", "status")
        normalized$summary <- normalized$summary[intersect(names(normalized$summary), keep)]
    }
    normalized$artifacts <- NULL
    normalized$inputs <- NULL
    normalized
}

sclet_ai_target_exposes_full_evidence <- function(target) {
    # Internal ledger consumers call GetAnalysisLedger() without a target and
    # must see full evidence; only explicit interpretive targets are bounded.
    if (is.null(target) || !length(target)) return(TRUE)
    "report" %in% target
}

sclet_ai_collect_records <- function(
    object,
    state,
    detail,
    target,
    include_artifacts,
    include_data
) {
    entries <- list()
    keys <- character()
    warnings <- list()

    add_entry <- function(record, key, type, source) {
        if (is.null(record)) {
            return(invisible(NULL))
        }
        if (!is.list(record)) {
            record <- list(value = record)
        }
        id <- record$id
        if ((is.null(id) || !length(id)) && !identical(source, "legacy")) {
            id <- key
        }
        record_type <- record$type
        if (is.null(record_type) || !length(record_type)) {
            record_type <- type
        }
        if (identical(record_type, "ai_evidence") && !sclet_ai_evidence_record_is_completed(record)) {
            return(invisible(NULL))
        }
        entry_key <- as.character(key)
        if (!nzchar(entry_key)) {
            entry_key <- paste(record_type, id, sep = "::")
        }
        if (entry_key %in% keys) {
            entry_key <- paste(source, entry_key, sep = "::")
        }
        while (entry_key %in% keys) {
            entry_key <- paste0(entry_key, "_2")
        }
        match_id <- if (is.null(id)) key else id
        if (!sclet_ai_target_matches(record_type, match_id, target)) {
            return(invisible(NULL))
        }
        normalized <- sclet_ai_normalize_record(
            record = record,
            id = id,
            type = record_type,
            source = source,
            detail = detail,
            include_artifacts = include_artifacts,
            include_data = include_data
        )
        if (identical(record_type, "ai_evidence") && !sclet_ai_target_exposes_full_evidence(target)) {
            normalized <- sclet_ai_strip_evidence_for_interpretive(normalized)
        }
        entries[[entry_key]] <<- normalized
        keys <<- c(keys, entry_key)
        invisible(NULL)
    }

    registry <- state$analyses
    if (!is.list(registry)) {
        registry <- list()
    }
    if (length(registry)) {
        registry_names <- names(registry)
        if (is.null(registry_names)) {
            registry_names <- paste0("analysis_", seq_along(registry))
        }
        for (i in seq_along(registry)) {
            record <- registry[[i]]
            fallback_type <- if (is.list(record)) record$type else NULL
            if (is.null(fallback_type)) {
                fallback_type <- registry_names[[i]]
            }
            add_entry(record, registry_names[[i]], fallback_type, "analysis")
        }
    }

    accessor_context <- tryCatch(get_analysis_context(object), error = function(e) NULL)
    accessor_records <- if (is.list(accessor_context)) accessor_context$records else NULL
    if (is.list(accessor_records) && length(accessor_records)) {
        for (nm in names(accessor_records)) {
            add_entry(accessor_records[[nm]], nm, nm, "accessor")
        }
    }

    legacy_types <- c("trajectory", "milo", "aggregation", "batch")
    for (type in legacy_types) {
        legacy_record <- tryCatch(
            sclet_get_legacy_analysis_record(object, type),
            error = function(e) NULL
        )
        if (!is.null(legacy_record)) {
            add_entry(legacy_record, paste0("legacy_", type), type, "legacy")
        }
    }

    state_records <- state$states$records
    if (!is.list(state_records)) {
        state_records <- list()
    }
    state_view <- list()
    if (length(state_records)) {
        for (type in names(state_records)) {
            type_records <- state_records[[type]]
            if (!is.list(type_records)) {
                next
            }
            type_view <- list()
            record_names <- names(type_records)
            if (is.null(record_names)) {
                record_names <- paste0(type, "_", seq_along(type_records))
            }
            for (i in seq_along(type_records)) {
                id <- if (is.list(type_records[[i]])) type_records[[i]]$id else NULL
                if (is.null(id)) {
                    id <- record_names[[i]]
                }
                if (identical(type, "ai_evidence") && !sclet_ai_evidence_record_is_completed(type_records[[i]])) {
                    next
                }
                if (sclet_ai_target_matches(type, id, target)) {
                    normalized <- sclet_ai_normalize_record(
                        record = type_records[[i]],
                        id = id,
                        type = type,
                        source = "state",
                        detail = detail,
                        include_artifacts = include_artifacts,
                        include_data = include_data
                    )
                    if (identical(type, "ai_evidence") && !sclet_ai_target_exposes_full_evidence(target)) {
                        normalized <- sclet_ai_strip_evidence_for_interpretive(normalized)
                    }
                    type_view[[record_names[[i]]]] <- normalized
                }
                add_entry(
                    record = type_records[[i]],
                    key = paste(type, id, sep = "::"),
                    type = type,
                    source = "state"
                )
            }
            if (length(type_view)) {
                state_view[[type]] <- type_view
            }
        }
    }

    workflow_view <- list()
    if (length(entries)) {
        for (nm in names(entries)) {
            record <- entries[[nm]]
            type <- record$type
            id <- record$id
            is_workflow <- grepl("workflow", tolower(paste(nm, type, id))) ||
                identical(tolower(as.character(type)), "state_priority")
            if (isTRUE(is_workflow)) {
                workflow_view[[nm]] <- record
            }
        }
    }

    lineage <- list()
    if (length(entries)) {
        for (nm in names(entries)) {
            record <- entries[[nm]]
            parents <- record$parents
            if (is.null(parents) && is.list(record$inputs)) {
                parents <- record$inputs$parents
            }
            lineage[[nm]] <- list(
                id = record$id,
                type = record$type,
                parents = if (is.null(parents)) NULL else parents
            )
        }
    }

    list(
        analyses = entries,
        state_records = state_view,
        workflows = workflow_view,
        lineage = lineage,
        warnings = warnings
    )
}

sclet_ai_build_analysis_story <- function(object, records, status) {
    analyses <- records$analyses
    if (!is.list(analyses)) {
        analyses <- list()
    }
    lineage <- records$state_records
    if (!is.list(records$lineage)) {
        lineage_map <- list()
    } else {
        lineage_map <- records$lineage
    }

    timeline <- list()
    unordered_keys <- character()
    analysis_keys <- names(analyses)

    for (key in analysis_keys) {
        record <- analyses[[key]]
        created_at <- record$created_at
        if (is.null(created_at) || !length(created_at)) {
            unordered_keys <- c(unordered_keys, key)
            created_at_str <- NULL
        } else {
            created_at_str <- as.character(created_at)
        }

        parents <- lineage_map[[key]]$parents
        if (is.null(parents)) {
            parents <- character()
        } else {
            parents <- as.character(parents)
        }
        parents <- intersect(parents, analysis_keys)

        timeline[[length(timeline) + 1L]] <- list(
            key = key,
            type = as.character(record$type %||% ""),
            id = if (is.null(record$id)) NULL else as.character(record$id),
            created_at = created_at_str,
            status = if (is.null(record$status)) NULL else as.character(record$status),
            method = if (is.null(record$method)) NULL else as.character(record$method),
            depends_on = parents
        )
    }

    if (length(timeline) > 1L) {
        has_time <- vapply(timeline, function(x) !is.null(x$created_at), logical(1L))
        if (any(has_time)) {
            parse_time <- function(x) {
                if (is.null(x$created_at)) return(NA_real_)
                t <- tryCatch(as.POSIXct(x$created_at, tryFormats = c(
                    "%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M:%OS",
                    "%Y/%m/%d %H:%M:%S", "%Y-%m-%dT%H:%M:%S",
                    "%Y-%m-%dT%H:%M:%OS"
                )), error = function(e) NA_real_)
                if (is.na(t)) NA_real_ else as.numeric(t)
            }
            times <- vapply(timeline, parse_time, numeric(1L))
            with_time <- which(!is.na(times))
            without_time <- which(is.na(times))
            if (length(with_time) > 1L) {
                ord <- order(times[with_time])
                with_time <- with_time[ord]
            }
            timeline <- timeline[c(with_time, without_time)]
        }
    }

    state_records <- records$state_records
    if (!is.list(state_records)) {
        state_records <- list()
    }

    user_decisions <- list()
    # Read user decisions from the object's raw state, not from the
    # target-filtered record list.
    #
    # The target filter decides which records a given context includes, and
    # `ai_evidence` is classified under the `status` target -- so a planning
    # context used to receive an EMPTY user_decisions list even when the user
    # had answered a clarification question. That made the clarification loop
    # cosmetic: the answer was recorded, and then the next planning round could
    # not see it, which is precisely the gap the structured-clarification work
    # was meant to close.
    #
    # Raw state, not sclet_ai_evidence_get_all(): the latter builds a full
    # ledger, which builds this story, which would call it again. Only the id,
    # a timestamp and a bounded summary string are copied, so this adds no
    # cell-level or free-form data to the context.
    raw_state <- tryCatch(sclet_get_state(object), error = function(e) list())
    raw_evidence <- raw_state$states$records$ai_evidence %||% list()
    if (is.list(raw_evidence) && length(raw_evidence)) {
        for (nm in names(raw_evidence)) {
            rec <- raw_evidence[[nm]]
            summary <- if (is.list(rec)) rec$summary else NULL
            if (!is.list(summary)) summary <- list()
            if (!identical(summary$kind, "user_decision")) next
            decision_summary <- summary$summary
            if (is.null(decision_summary) || !is.character(decision_summary) || !length(decision_summary)) {
                decision_summary <- "user decision recorded"
            } else {
                decision_summary <- as.character(decision_summary)[1L]
            }
            user_decisions[[length(user_decisions) + 1L]] <- list(
                id = as.character(summary$id %||% nm),
                created_at = if (is.null(summary$created_at)) NULL else as.character(summary$created_at),
                summary = decision_summary
            )
        }
    }

    design_confirmations <- list()
    design_records <- state_records[["ai_design_confirmation"]]
    if (is.list(design_records) && length(design_records)) {
        for (nm in names(design_records)) {
            rec <- design_records[[nm]]
            summary <- rec$summary
            if (!is.list(summary)) summary <- list()
            design <- summary$design
            if (!is.list(design)) next
            confirmed_at <- summary$confirmed_at
            design_confirmations[[length(design_confirmations) + 1L]] <- list(
                id = as.character(rec$id %||% nm),
                created_at = if (is.null(confirmed_at)) NULL else as.character(confirmed_at),
                design = design
            )
        }
    }

    evidence_gaps <- list()
    has_rare_cells <- "rare_cells" %in% analysis_keys ||
        any(vapply(names(state_records), function(type) {
            type_records <- state_records[[type]]
            if (!is.list(type_records)) return(FALSE)
            any(vapply(type_records, function(r) {
                is.list(r) && identical(as.character(r$type %||% ""), "rare_cells")
            }, logical(1L)))
        }, logical(1L)))
    if (has_rare_cells) {
        readiness <- tryCatch(
            check_rare_cell_readiness(object),
            error = function(e) NULL
        )
        if (!is.null(readiness)) {
            doublet_available <- readiness$checks$doublet_evidence_available
            if (isFALSE(doublet_available)) {
                notes <- readiness$notes
                if (is.character(notes) && length(notes)) {
                    note <- notes[1L]
                } else {
                    note <- "doublet evidence missing"
                }
                evidence_gaps[[1L]] <- list(
                    analysis_type = "rare_cells",
                    gap = "doublet_evidence_available",
                    reason = note
                )
            }
        }
    }

    conflicts <- list()
    # A "conflict" means several records compete for one active slot, so the
    # reader cannot tell which is in force. That is true of the computational
    # state types (two detest records, two clusterings). It is NOT true of the
    # append-only AI log types: one ai_execution record is written per round and
    # one ai_evidence node per measurement, by design. Flagging those made the
    # conflicts field fire on every analysis with more than one round -- a
    # warning that always fires teaches the reader to ignore it (observed on
    # the first multi-round run, 2026-10-19).
    append_only_types <- c("ai_execution", "ai_evidence", "ai_design_confirmation")
    type_to_keys <- list()
    for (key in analysis_keys) {
        record <- analyses[[key]]
        type <- as.character(record$type %||% "")
        if (!nzchar(type)) next
        if (type %in% append_only_types) next
        if (is.null(type_to_keys[[type]])) {
            type_to_keys[[type]] <- character()
        }
        type_to_keys[[type]] <- c(type_to_keys[[type]], key)
    }
    for (type in names(type_to_keys)) {
        keys <- type_to_keys[[type]]
        if (length(keys) > 1L) {
            conflicts[[length(conflicts) + 1L]] <- list(
                type = type,
                keys = keys,
                reason = paste0("multiple ", type, " records exist; only the active one (if any) is in active_states")
            )
        }
    }

    list(
        schema_version = "1.1",
        timeline = timeline,
        n_steps = length(timeline),
        unordered_steps = unordered_keys,
        user_decisions = user_decisions,
        design_confirmations = design_confirmations,
        evidence_gaps = evidence_gaps,
        conflicts = conflicts,
        # M12-b: what the data can and cannot answer, and which analyses the
        # ledger has never run. Deterministic, no model call, no raw values.
        # `evidence_gaps` says "an analysis lacks evidence"; `data_support`
        # says what is available at all -- a planner needs both to choose a
        # next round instead of repeating the last one.
        data_support = sclet_ai_build_data_support(object, records),
        raw_values_included = FALSE
    )
}

sclet_ai_normalize_record <- function(
    record,
    id,
    type,
    source,
    detail,
    include_artifacts,
    include_data
) {
    if (!is.list(record)) {
        record <- list(value = record)
    }
    result <- list(
        id = id,
        type = type,
        source = source
    )
    fields <- c(
        "status", "method", "inputs", "summary", "parents", "created_at",
        "validation", "warnings", "software"
    )
    if (identical(detail, "full")) {
        fields <- c(fields, "params")
    }
    for (field in fields) {
        if (!is.null(record[[field]])) {
            result[[field]] <- sclet_ai_safe_value(
                record[[field]],
                include_data = include_data
            )
        }
    }
    if (isTRUE(include_artifacts) && !is.null(record$artifacts)) {
        result$artifacts <- sclet_ai_safe_value(
            record$artifacts,
            include_data = include_data
        )
    }
    if (identical(detail, "full")) {
        known <- c(
            "id", "type", "status", "method", "inputs", "summary", "parents",
            "created_at", "validation", "warnings", "software", "params",
            "artifacts"
        )
        extra <- setdiff(names(record), known)
        if (length(extra)) {
            for (field in sort(extra)) {
                result[[field]] <- sclet_ai_safe_value(
                    record[[field]],
                    include_data = include_data
                )
            }
        }
    }
    result
}

sclet_ai_dataset_summary <- function(object, include_data = FALSE) {
    col_data <- SummarizedExperiment::colData(object)
    row_data <- SummarizedExperiment::rowData(object)
    assay_names <- SummarizedExperiment::assayNames(object)
    modalities <- tryCatch(SingleCellExperiment::altExpNames(object), error = function(e) character())
    result <- list(
        n_cells = ncol(object),
        n_genes = nrow(object),
        assays = assay_names,
        layers = tryCatch(Layers(object), error = function(e) character()),
        reductions = SingleCellExperiment::reducedDimNames(object),
        modalities = modalities,
        columns = colnames(col_data),
        row_data_columns = colnames(row_data)
    )
    if (isTRUE(include_data)) {
        result$data <- list(
            included = TRUE,
            expression = lapply(assay_names, function(name) {
                value <- SummarizedExperiment::assay(object, name)
                list(
                    name = name,
                    class = class(value),
                    dimensions = dim(value),
                    values_included = FALSE
                )
            }),
            coldata_preview = sclet_ai_table_preview(col_data),
            rowdata_preview = sclet_ai_table_preview(row_data)
        )
        names(result$data$expression) <- assay_names
    }
    result
}

sclet_ai_table_preview <- function(value, max_rows = 10L, max_columns = 25L) {
    columns <- colnames(value)
    if (is.null(columns)) {
        columns <- character()
    }
    columns <- utils::head(columns, max_columns)
    preview <- list(
        class = class(value),
        n_rows = nrow(value),
        n_columns = ncol(value),
        columns = columns
    )
    if (length(columns) && nrow(value)) {
        row_ids <- seq_len(min(nrow(value), max_rows))
        preview$rows <- lapply(row_ids, function(i) {
            row <- lapply(columns, function(column) {
                sclet_ai_safe_value(value[[column]][[i]], include_data = TRUE)
            })
            names(row) <- columns
            row
        })
    }
    preview
}

sclet_ai_safe_value <- function(value, include_data = FALSE, depth = 0L) {
    if (is.null(value)) {
        return(NULL)
    }
    if (depth > 8L) {
        return(list(kind = "nested_value", class = class(value)))
    }
    if (is.matrix(value) || inherits(value, "Matrix") || inherits(value, "DelayedMatrix")) {
        return(list(
            kind = "matrix",
            class = class(value),
            dimensions = dim(value),
            values_included = FALSE
        ))
    }
    if (is.data.frame(value)) {
        return(sclet_ai_table_preview(value))
    }
    if (isS4(value)) {
        return(list(kind = "object", class = class(value)))
    }
    if (is.atomic(value)) {
        if (inherits(value, c("POSIXct", "POSIXlt", "Date"))) {
            return(as.character(value))
        }
        if (length(value) <= 100L) {
            return(value)
        }
        return(list(
            kind = "vector",
            class = class(value),
            length = length(value),
            values = if (isTRUE(include_data)) utils::head(value, 25L) else NULL
        ))
    }
    if (is.list(value)) {
        if (length(value) == 0L) {
            return(list())
        }
        result <- lapply(value, sclet_ai_safe_value, include_data = include_data, depth = depth + 1L)
        names(result) <- names(value)
        result
    } else {
        list(kind = "value", class = class(value))
    }
}

#' Describe an AIExecutionRegistry's action vocabulary for model-facing context
#'
#' Exposes only the fields a planning model needs to propose a valid,
#' correctly-parameterized action: `name`, `description`, `input_schema`,
#' `requires_confirmation`, `mutates_object`, `estimated_cost`. Deliberately
#' excludes `handler` and `prerequisites`, which are R closures and must
#' never be serialized into a prompt.
#'
#' @param registry An `sclet_ai_execution_registry`, or `NULL`.
#' @return A named list keyed by action name, or an empty list if `registry`
#'   is `NULL` or empty.
sclet_ai_describe_registry <- function(registry) {
    if (is.null(registry) || !length(registry)) {
        return(list())
    }
    actions <- lapply(registry, function(descriptor) {
        # `output_schema` is exposed because a planning model needs to know the
        # concrete names a step will *produce* -- e.g. that run_pca registers a
        # reduction called "PCA" -- in order to fill a downstream action's
        # `reduction`/`assay` parameter. Without it, a real model wrote the
        # descriptive phrase "PCA name from step_4" where the literal "PCA" was
        # required (measured: sclet-benchmark-ai results/ai_plan_repair.jsonl).
        #
        # `note` is stripped: every output_schema note is prose guidance, and
        # the privacy classifier hard-denies the field name `note`, so leaving
        # it in would redact it anyway and emit a spurious privacy warning on
        # every planning call (verified: all five actions that carry a note
        # produced 1 redaction + 1 warning each when it was included).
        output_schema <- descriptor$output_schema %||% list()
        output_schema[["note"]] <- NULL
        list(
            name = descriptor$name,
            description = descriptor$description,
            input_schema = descriptor$input_schema,
            output_schema = output_schema,
            requires_confirmation = isTRUE(descriptor$requires_confirmation),
            mutates_object = isTRUE(descriptor$mutates_object),
            estimated_cost = descriptor$estimated_cost %||% NA_character_
        )
    })
    names(actions) <- names(registry)
    actions
}

sclet_ai_capabilities <- function(status, records, object = NULL) {
    available <- status$available
    if (!is.list(available)) {
        available <- list()
    }
    result <- list(
        read_only_context = TRUE,
        assays = sclet_ai_or(available$assays, character()),
        layers = sclet_ai_or(available$layers, character()),
        reductions = sclet_ai_or(available$reductions, character()),
        graphs = sclet_ai_or(available$graphs, character()),
        recorded_analyses = names(records$analyses),
        recorded_workflows = names(records$workflows),
        execution = NA
    )
    if (!is.null(object)) {
        result$identity_levels <- sclet_ai_build_identity_levels(object)
        # Contrast levels for colData columns, separate from identity levels.
        #
        # Measured 2026-10-21 (M12-a batch, commit 5947e4b, the only failure in
        # 12 same-tree trials): a round-2 plan wrote `ident.1/ident.2` as levels
        # of `condition` that the object does not have, and died with
        # `contrast_invalid: ident.1/ident.2 must be levels of 'condition':
        # ctrl, stim`. The planning prompt told the model to take group-valued
        # parameters from `identity_levels$codes`, but that block describes the
        # ACTIVE IDENT (the clustering), not a design column -- and
        # `run_de_test(group = "condition")` contrasts levels of a colData
        # column, which had no representation in the context at all. So the
        # instruction pointed at the one place the answer could not be.
        result$contrast_levels <- sclet_ai_build_contrast_levels(object)
    }
    result
}

# Levels of the colData columns that a `group` contrast can name, for the
# planner to read rather than guess.
#
# Bounded on purpose: only low-cardinality columns are reported (a donor column
# with 16 levels is a design role, not a usable contrast on its own), and the
# count is reported even when the levels themselves are withheld, so "I did not
# show you these" and "there are none" stay distinguishable.
sclet_ai_build_contrast_levels <- function(object, max_levels = 12L) {
    coldata <- tryCatch(SummarizedExperiment::colData(object), error = function(e) NULL)
    if (is.null(coldata)) {
        return(list(columns = list(), note = "colData unavailable", raw_values_included = FALSE))
    }
    # The active identity column is the CLUSTERING, and its values are the
    # answer the analysis is supposed to discover -- `identity_levels` already
    # exposes them as anonymised codes precisely so they are not handed over as
    # text. Re-listing the same column here under its real labels would undo
    # that: measured 2026-10-21, a mock object with clusters named `treat1` /
    # `treat2` had those labels appear in `contrast_levels` after the block was
    # added. A clustering is not a design role, so it is excluded by name.
    ident_column <- tryCatch(ActiveIdent(object), error = function(e) NULL)
    if (is.null(ident_column) || !nzchar(ident_column)) ident_column <- NULL
    entries <- list()
    for (column in colnames(coldata)) {
        if (!is.null(ident_column) && identical(column, ident_column)) next
        values <- tryCatch(as.character(coldata[[column]]), error = function(e) NULL)
        if (is.null(values)) next
        values <- values[!is.na(values)]
        levels <- sort(unique(values))
        n_levels <- length(levels)
        if (n_levels < 2L) next
        usable <- n_levels <= max_levels
        entries[[length(entries) + 1L]] <- list(
            column = column,
            n_levels = as.integer(n_levels),
            # A character VECTOR, not as.list(): the privacy layer redacts an
            # unnamed list of scalars as "unnamed nested values are not
            # bounded", which silently removed every level and left the block
            # useless. A vector is bounded by n_levels and passes as-is.
            levels = if (usable) levels else character(),
            levels_available = usable
        )
    }
    # The field is `guidance`, not `note`: `note` is on the privacy layer's deny
    # list (`|note|comment`), so text explaining the block was redacted as a
    # sensitive field while the block stayed.
    list(
        columns = entries,
        max_levels = as.integer(max_levels),
        guidance = paste(
            "Use these for `group` contrasts (run_de_test with group = <column>",
            "and ident.1 / ident.2 as two of that column's levels).",
            "`identity_levels` is the clustering identity and is NOT what a",
            "condition contrast takes."
        ),
        raw_values_included = FALSE
    )
}

sclet_ai_build_identity_levels <- function(object) {
    mapping <- sclet_ai_ident_group_codes(object)
    if (is.null(mapping) || is.null(mapping$codes) || !length(mapping$codes)) {
        return(list(
            column = NULL,
            n_levels = 0L,
            aggregation_threshold = 10L,
            codes = list(),
            omitted_below_threshold = 0L,
            all_levels_available = TRUE,
            raw_values_included = FALSE
        ))
    }
    threshold <- 10L
    keep <- mapping$sizes >= threshold
    n_omitted <- sum(!keep)
    codes_list <- lapply(which(keep), function(i) {
        list(code = mapping$codes[[i]], n_cells = as.integer(mapping$sizes[[i]]))
    })
    list(
        column = mapping$col,
        n_levels = as.integer(sum(keep)),
        aggregation_threshold = threshold,
        codes = codes_list,
        omitted_below_threshold = as.integer(n_omitted),
        all_levels_available = n_omitted == 0L,
        raw_values_included = FALSE
    )
}

sclet_ai_ident_group_codes <- function(object) {
    col <- tryCatch(ActiveIdent(object), error = function(e) NULL)
    if (is.null(col) || !nzchar(col)) {
        return(list(col = NULL, codes = character(), labels = character(), sizes = integer(), n_groups = 0L))
    }
    # NOTE: `col` may legitimately be "colLabels" -- that is what FindClusters()
    # sets, and it is the most common way an object acquires an identity.
    # Idents() already resolves the colLabels case (R/clustering.R), so do NOT
    # bail out here: an earlier version returned an empty mapping for
    # "colLabels", which silently made group codes unusable on the standard
    # clustering path while still working on objects with an explicit colData
    # ident column. A `colLabels`-identified object has a perfectly usable
    # identity, so it must produce a mapping like any other.
    labels <- tryCatch(as.character(Idents(object)), error = function(e) NULL)
    if (is.null(labels) || !length(labels)) {
        return(list(col = col, codes = character(), labels = character(), sizes = integer(), n_groups = 0L))
    }
    sorted_labels <- sort(unique(labels))
    sizes <- as.integer(table(factor(labels, levels = sorted_labels)))
    codes <- paste0("group_", seq_along(sorted_labels))
    list(
        col = col,
        codes = codes,
        labels = sorted_labels,
        sizes = sizes,
        n_groups = as.integer(length(sorted_labels))
    )
}

sclet_ai_resolve_ident_value <- function(object, value) {
    if (is.null(value) || !length(value)) {
        return(value)
    }
    if (!is.character(value)) {
        return(value)
    }
    if (!all(grepl("^group_[0-9]+$", value))) {
        return(value)
    }
    mapping <- sclet_ai_ident_group_codes(object)
    if (is.null(mapping) || !length(mapping$codes)) {
        return(value)
    }
    code_indices <- as.integer(sub("^group_", "", value))
    resolved <- vapply(seq_along(value), function(i) {
        idx <- code_indices[[i]]
        if (is.na(idx) || idx < 1L || idx > length(mapping$labels)) {
            return(value[[i]])
        }
        mapping$labels[[idx]]
    }, character(1L))
    resolved
}


# Resolve a value against an arbitrary colData column (condition, treatment).
#
# Cluster identity codes look like group_1 and are resolved against the active
# identity; a condition contrast must match that column's own values instead.
# Unknown values are returned untouched so validation reports them rather than
# silently dropping the request.
sclet_ai_resolve_group_value <- function(object, group, value) {
    if (is.null(value) || !length(value)) return(value)
    if (!is.character(value)) return(value)
    if (is.null(group) || !group %in% colnames(SummarizedExperiment::colData(object))) {
        return(sclet_ai_resolve_ident_value(object, value))
    }
    levels <- unique(as.character(SummarizedExperiment::colData(object)[[group]]))
    value[vapply(value, function(v) as.character(v) %in% levels, logical(1L))]
}

sclet_ai_blocked_actions <- function(health) {
    list(
        execute_analysis = list(
            blocked = NA,
            reason = "execution capability depends on the AIExecutionRegistry supplied at call time, not on this read-only context"
        ),
        missing_mainline_steps = sclet_ai_or(health$mainline_missing, character()),
        unavailable_prerequisites = sclet_ai_or(health$mainline_missing, character())
    )
}

sclet_ai_or <- function(value, fallback) {
    if (is.null(value)) fallback else value
}

sclet_ai_fingerprint <- function(value) {
    raw <- serialize(value, connection = NULL, version = 2)
    if (!length(raw)) {
        return("sclet-ai-00000000")
    }
    weights <- seq_along(raw) %% 100003
    checksum <- sum((as.double(raw) + 1) * (weights + 1)) %% 2147483647
    sprintf("sclet-ai-%08x", as.integer(checksum))
}
