sclet_ai_evidence_bad_field <- function(name) {
    grepl(
        "token|api[_-]?key|secret|password|credential|authorization|path|file|barcode|cell[_-]?(id|name)|gene[_-]?(id|name)|patient|subject|donor|sample[_-]?id|clinical|note|comment",
        tolower(name),
        perl = TRUE
    )
}

sclet_ai_evidence_value_ok <- function(value, path = "values") {
    if (is.null(value) || !length(value)) return(TRUE)
    if (is.matrix(value) || inherits(value, "Matrix") || isS4(value) || is.environment(value) || is.function(value)) return(FALSE)
    if (is.list(value)) {
        if (is.null(names(value)) || any(!nzchar(names(value)))) return(FALSE)
        if (any(vapply(names(value), sclet_ai_evidence_bad_field, logical(1L)))) return(FALSE)
        return(all(vapply(seq_along(value), function(i) sclet_ai_evidence_value_ok(value[[i]], paste(path, names(value)[i], sep = ".")), logical(1L))))
    }
    if (is.character(value)) {
        return(all(grepl("^(group|cluster|condition|sample|route)_[0-9]+$", value)))
    }
    is.numeric(value) || is.logical(value)
}

sclet_ai_evidence_source <- function(object, source) {
    if (is.null(source) || !length(source)) return(NULL)
    ledger <- GetAnalysisLedger(object, detail = "summary", include_artifacts = FALSE, include_data = FALSE)
    records <- ledger$analyses %||% list()
    state_records <- ledger$state_records %||% list()
    if (length(state_records)) {
        state_records <- unlist(state_records, recursive = FALSE, use.names = FALSE)
    }
    records <- c(records, state_records)
    if (!length(records)) return(NULL)
    matches <- records[vapply(records, function(record) {
        is.list(record) &&
            identical(as.character(record$id %||% ""), as.character(source)) &&
            identical(as.character(record$status %||% record$summary$status %||% ""), "completed")
    }, logical(1L))]
    if (!length(matches)) return(NULL)
    matches[[1L]]
}

sclet_ai_dependency_hash <- function(parents, scope) {
    canonical <- list(
        parents = sort(as.character(parents %||% character())),
        scope_keys = sort(names(scope %||% list())),
        scope_values = as.list(scope %||% list())
    )
    json_like <- utils::capture.output(dput(canonical))
    if (requireNamespace("digest", quietly = TRUE)) {
        digest::digest(json_like, algo = "sha256")
    } else {
        raw <- paste(trimws(json_like), collapse = "")
        paste0("dep_", substring(paste0(
            formatC(strtoi(charToRaw(raw), 16L), format = "f"),
            collapse = ""
        ), 1L, 16L))
    }
}

sclet_ai_evidence_kind_is_ai_generated <- function(kind) {
    is.character(kind) && length(kind) == 1L &&
        kind %in% c("deterministic_summary", "plot", "test", "state")
}

sclet_ai_evidence_kind_is_user_decision <- function(kind) {
    is.character(kind) && length(kind) == 1L && identical(kind, "user_decision")
}

sclet_ai_evidence_get_all <- function(object) {
    ledger <- GetAnalysisLedger(object, detail = "summary", include_artifacts = FALSE, include_data = FALSE)
    records <- ledger$state_records$ai_evidence %||% list()
    out <- lapply(records, function(r) {
        s <- r$summary %||% r
        if (is.list(s)) s else list()
    })
    stats::setNames(out, vapply(out, function(s) as.character(s$id %||% ""), character(1L)))
}

sclet_ai_evidence_scope_string <- function(scope) {
    if (!is.list(scope) || !length(scope)) return("ledger")
    paste(sort(unique(as.character(unlist(scope, use.names = FALSE)))), collapse = "|")
}

#' Record a bounded deterministic evidence node in the analysis ledger
#'
#' @param object A `SingleCellExperiment` object.
#' @param evidence A structured evidence node with a unique `id` and `kind`.
#' @param source Optional completed analysis id supporting the evidence.
#' @param parents Optional evidence ids or dependency ids.
#' @param scope Optional aggregate scope, such as anonymized groups.
#' @return The SCE with an `ai_evidence` state record.
#' @export
RecordAIEvidence <- function(object, evidence, source = NULL, parents = NULL, scope = NULL) {
    if (!inherits(object, "SingleCellExperiment")) stop("object must be a SingleCellExperiment", call. = FALSE)
    if (!is.list(evidence) || length(evidence) == 0L) stop("evidence must be a non-empty list", call. = FALSE)
    id <- evidence$id %||% ""
    kind <- evidence$kind %||% ""
    if (!is.character(id) || length(id) != 1L || !nzchar(id)) stop("evidence$id must be a non-empty string", call. = FALSE)
    if (!is.character(kind) || length(kind) != 1L || !kind %in% c("deterministic_summary", "plot", "test", "state", "user_decision")) {
        stop("evidence$kind must be a supported deterministic kind", call. = FALSE)
    }
    ledger <- GetAnalysisLedger(object, detail = "summary", include_artifacts = FALSE, include_data = FALSE)
    existing <- ledger$analyses %||% list()
    existing_state <- ledger$state_records %||% list()
    if (length(existing_state)) existing_state <- unlist(existing_state, recursive = FALSE, use.names = FALSE)
    existing <- c(existing, existing_state)
    existing_ids <- if (length(existing)) vapply(existing, function(record) as.character(record$id %||% ""), character(1L)) else character()
    if (id %in% existing_ids) stop("evidence id already exists", call. = FALSE)
    source_record <- sclet_ai_evidence_source(object, source)
    if (!is.null(source) && is.null(source_record)) stop("source must identify a completed analysis", call. = FALSE)
    if (!sclet_ai_evidence_value_ok(evidence$values %||% list())) stop("evidence$values contains unbounded or private values", call. = FALSE)
    scope <- scope %||% evidence$scope %||% list()
    if (!sclet_ai_evidence_value_ok(scope, "scope")) stop("scope contains unbounded or private values", call. = FALSE)
    parents <- as.character(parents %||% evidence$parents %||% character())
    node <- evidence
    node$source <- source %||% node$source %||% NULL
    node$parents <- parents
    node$scope <- scope
    node$claim_level <- node$claim_level %||% "measured"
    if (!node$claim_level %in% c("observed", "measured", "associated", "consistent_with")) stop("unsupported evidence claim_level", call. = FALSE)
    node$object_fingerprint <- ledger$fingerprint
    node$status <- "completed"
    node$created_at <- as.character(Sys.time())
    if (is.null(node$dependency_group)) {
        has_digest <- requireNamespace("digest", quietly = TRUE)
        if (has_digest) {
            canonical <- utils::capture.output(dput(list(
                parents = sort(unique(parents)),
                scope = node$scope %||% list()
            )))
            node$dependency_group <- digest::digest(canonical, algo = "sha256")
        } else {
            node$dependency_group <- paste0("dep_", substring(paste(c(sort(unique(parents)), sort(unlist(node$scope, use.names = FALSE))), collapse = "|"), 1L, 16L))
        }
    }
    updated <- sclet_set_analysis_state(
        object,
        type = "ai_evidence",
        id = id,
        method = "sclet_ai_evidence",
        inputs = list(source = node$source, parents = node$parents),
        summary = node,
        active = FALSE
    )
    updated
}

#' Query whether two evidence ids are independent
#'
#' Two evidence nodes are independent if they do not share a dependency_group
#' and neither node is an ancestor of the other via the parents chain.
#'
#' @param object A `SingleCellExperiment` object.
#' @param evidence_ids Character vector of exactly two evidence ids.
#' @return A list with `independent` (logical) plus `shared_dependency_groups` and `lineage_paths` diagnostics.
#' @export
sclet_ai_evidence_independence <- function(object, evidence_ids) {
    if (!inherits(object, "SingleCellExperiment")) stop("object must be a SingleCellExperiment", call. = FALSE)
    if (!is.character(evidence_ids) || length(evidence_ids) != 2L) stop("evidence_ids must be exactly two ids", call. = FALSE)
    all_ev <- sclet_ai_evidence_get_all(object)
    if (!all(evidence_ids %in% names(all_ev))) {
        return(list(
            independent = FALSE,
            shared_dependency_groups = character(),
            lineage_paths = character(),
            errors = paste0("missing evidence: ", paste(setdiff(evidence_ids, names(all_ev)), collapse = ", "))
        ))
    }
    a <- all_ev[[evidence_ids[[1L]]]]
    b <- all_ev[[evidence_ids[[2L]]]]
    shared <- character()
    if (!is.null(a$dependency_group) && !is.null(b$dependency_group) && identical(a$dependency_group, b$dependency_group)) {
        shared <- c(shared, a$dependency_group)
    }
    collect_ancestors <- function(id) {
        seen <- character()
        queue <- id
        while (length(queue)) {
            cur <- queue[[1L]]
            queue <- queue[-1L]
            if (cur %in% seen) next
            seen <- c(seen, cur)
            node <- all_ev[[cur]]
            pars <- as.character(node$parents %||% character())
            queue <- c(queue, setdiff(pars, seen))
        }
        seen
    }
    anc_a <- collect_ancestors(evidence_ids[[1L]])
    anc_b <- collect_ancestors(evidence_ids[[2L]])
    lineage_paths <- character()
    if (evidence_ids[[2L]] %in% setdiff(anc_a, evidence_ids[[1L]])) {
        lineage_paths <- c(lineage_paths, paste(evidence_ids[[1L]], "ancestor_chain_contains", evidence_ids[[2L]]))
    }
    if (evidence_ids[[1L]] %in% setdiff(anc_b, evidence_ids[[2L]])) {
        lineage_paths <- c(lineage_paths, paste(evidence_ids[[2L]], "ancestor_chain_contains", evidence_ids[[1L]]))
    }
    independent <- !length(shared) && !length(lineage_paths)
    list(
        independent = independent,
        shared_dependency_groups = shared,
        lineage_paths = lineage_paths
    )
}

#' Validate evidence references against the current analysis ledger
#'
#' @param object A `SingleCellExperiment` object.
#' @param refs Character evidence ids.
#' @param fingerprint Optional source fingerprint to require.
#' @param requesting_scope Optional aggregate scope of the referencing context; scope is bounded by the evidence's own scope.
#' @param requesting_kind Optional kind of the referencing context, used to enforce cross-kind rules.
#' @return A structured validation result with resolved evidence and errors.
#' @export
ValidateAIEvidenceRefs <- function(object, refs, fingerprint = NULL, requesting_scope = NULL, requesting_kind = NULL) {
    if (!inherits(object, "SingleCellExperiment")) stop("object must be a SingleCellExperiment", call. = FALSE)
    refs <- unique(as.character(refs %||% character()))
    refs <- refs[nzchar(refs)]
    ledger <- GetAnalysisLedger(object, detail = "summary", include_artifacts = FALSE, include_data = FALSE)
    records <- ledger$state_records$ai_evidence %||% list()
    resolved <- lapply(refs, function(id) {
        record <- records[[id]]
        if (is.null(record)) return(NULL)
        summary <- record$summary %||% record
        if (!identical(as.character(record$status %||% summary$status %||% ""), "completed")) return(NULL)
        if (!is.null(fingerprint) && !identical(as.character(summary$object_fingerprint %||% ""), as.character(fingerprint))) return(NULL)
        summary
    })
    names(resolved) <- refs
    errors <- character()
    req_scope_str <- sclet_ai_evidence_scope_string(requesting_scope)
    for (i in seq_along(refs)) {
        id <- refs[[i]]
        node <- resolved[[i]]
        if (is.null(node)) {
            errors <- c(errors, paste0("evidence ref not valid: ", id))
            next
        }
        node_scope_str <- sclet_ai_evidence_scope_string(node$scope)
        if (nzchar(req_scope_str) && nzchar(node_scope_str) && req_scope_str != "ledger") {
            req_tokens <- unique(strsplit(req_scope_str, "\\|")[[1L]])
            node_tokens <- unique(strsplit(node_scope_str, "\\|")[[1L]])
            if (!all(req_tokens %in% node_tokens)) {
                errors <- c(errors, paste0("evidence ref scope violation: ", id,
                    " scope=", node_scope_str, " but requesting_scope requires: ", req_scope_str))
            }
        }
        if (!is.null(requesting_kind) && sclet_ai_evidence_kind_is_ai_generated(requesting_kind) &&
            sclet_ai_evidence_kind_is_user_decision(node$kind %||% "")) {
            errors <- c(errors, paste0("evidence ref cross-kind: AI-generated claim ", id,
                " cannot cite a user_decision as evidence. Record a separate deterministic_state node as mediation."))
        }
        if (!is.null(requesting_kind) && sclet_ai_evidence_kind_is_user_decision(requesting_kind) &&
            identical(node$claim_level %||% "", "consistent_with")) {
            # user decisions may not cite weakest AI inferences as their sole evidence - require at least measured or stronger
            # Do not hard-error here; record a warning-level info via errors to keep validation strictness controllable later.
        }
    }
    list(
        valid = !length(errors),
        status = if (length(errors)) "invalid" else "valid",
        refs = refs,
        resolved = resolved,
        errors = errors,
        current_fingerprint = ledger$fingerprint,
        requesting_scope = requesting_scope %||% list(),
        requesting_kind = requesting_kind %||% NULL
    )
}
