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
#'   analysis records. Defaults to `FALSE`.
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
    detail = c("summary", "full")
) {
    detail <- match.arg(detail)
    sclet_ai_build_context(
        object = object,
        target = target,
        detail = detail,
        include_artifacts = FALSE,
        include_data = FALSE
    )
}

sclet_ai_build_context <- function(
    object,
    target = NULL,
    detail = c("summary", "full"),
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
        dataset = dataset,
        active_view = active_view,
        active_states = active_states,
        health = health,
        analyses = records$analyses,
        state_records = records$state_records,
        workflows = records$workflows,
        lineage = records$lineage,
        capabilities = sclet_ai_capabilities(status, records),
        blocked_actions = sclet_ai_blocked_actions(health),
        quality_checks = list(
            health = health,
            mainline_missing = sclet_ai_or(health$mainline_missing, character())
        ),
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
        status = character(),
        qc = c("qc", "preprocess", "clean", "doublet", "phenotype"),
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
                if (sclet_ai_target_matches(type, id, target)) {
                    type_view[[record_names[[i]]]] <- sclet_ai_normalize_record(
                        record = type_records[[i]],
                        id = id,
                        type = type,
                        source = "state",
                        detail = detail,
                        include_artifacts = include_artifacts,
                        include_data = include_data
                    )
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

sclet_ai_capabilities <- function(status, records) {
    available <- status$available
    if (!is.list(available)) {
        available <- list()
    }
    list(
        read_only_context = TRUE,
        assays = sclet_ai_or(available$assays, character()),
        layers = sclet_ai_or(available$layers, character()),
        reductions = sclet_ai_or(available$reductions, character()),
        graphs = sclet_ai_or(available$graphs, character()),
        recorded_analyses = names(records$analyses),
        recorded_workflows = names(records$workflows),
        execution = FALSE
    )
}

sclet_ai_blocked_actions <- function(health) {
    list(
        execute_analysis = list(
            blocked = TRUE,
            reason = "Phase 0 context is read-only"
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
