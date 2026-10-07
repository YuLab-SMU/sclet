sclet_ai_payload_policy <- function(
    object,
    context,
    privacy = c("strict", "standard", "local")
) {
    privacy <- match.arg(privacy)
    scan <- sclet_ai_payload_scan(
        object = object,
        context = context,
        privacy = privacy,
        user_consent = FALSE
    )
    sclet_ai_payload_policy_from_scan(scan)
}

sclet_ai_sanitize_context <- function(
    object,
    context,
    privacy = c("strict", "standard", "local"),
    user_consent = FALSE
) {
    privacy <- match.arg(privacy)
    if (length(user_consent) != 1L || is.na(user_consent)) {
        stop("user_consent must be a single non-missing logical value")
    }
    scan <- sclet_ai_payload_scan(
        object = object,
        context = context,
        privacy = privacy,
        user_consent = isTRUE(user_consent)
    )
    list(
        payload = scan$payload,
        policy = sclet_ai_payload_policy_from_scan(scan)
    )
}

sclet_ai_payload_policy_from_scan <- function(scan) {
    payload <- sclet_ai_payload_canonical(scan$payload)
    serialized <- paste(utils::capture.output(dput(payload)), collapse = "")
    estimated_tokens <- max(
        1L,
        as.integer(ceiling(nchar(serialized, type = "bytes") / 4))
    )
    fingerprint <- sclet_ai_fingerprint(payload)
    list(
        allowed_fields = sort(unique(scan$allowed_fields)),
        redacted_fields = sort(unique(scan$redacted_fields)),
        aggregation_threshold = 10L,
        estimated_tokens = estimated_tokens,
        requires_user_consent = isTRUE(scan$requires_user_consent),
        warnings = unique(scan$warnings),
        payload_fingerprint = sub(
            "^sclet-ai-",
            "sclet-ai-payload-",
            fingerprint
        )
    )
}

sclet_ai_payload_scan <- function(
    object,
    context,
    privacy,
    user_consent = FALSE
) {
    invisible(object)
    state <- new.env(parent = emptyenv())
    state$allowed_fields <- character()
    state$redacted_fields <- character()
    state$warnings <- character()
    state$requires_user_consent <- FALSE

    add_allowed <- function(path) {
        if (nzchar(path)) {
            state$allowed_fields <- c(state$allowed_fields, path)
        }
    }
    add_redacted <- function(path, reason, consent = FALSE) {
        if (nzchar(path)) {
            state$redacted_fields <- c(state$redacted_fields, path)
        }
        if (isTRUE(consent)) {
            state$requires_user_consent <- TRUE
        }
        warning <- paste0("Redacted ", if (nzchar(path)) path else "context", ": ", reason)
        state$warnings <- c(state$warnings, warning)
    }
    unknown_field <- function(path) {
        if (identical(privacy, "strict")) {
            condition <- structure(
                list(message = paste0("Unknown outbound payload field: ", path)),
                class = c("sclet_ai_privacy_error", "error", "condition")
            )
            stop(condition)
        }
        add_redacted(path, "field is not allowlisted", consent = !identical(privacy, "local"))
        NULL
    }

    visit <- function(
        value,
        path = "",
        name = "",
        depth = 0L,
        parent_kind = "container"
    ) {
        if (depth > 8L) {
            add_redacted(path, "maximum nesting depth", consent = TRUE)
            return(NULL)
        }

        kind <- sclet_ai_payload_field_kind(name, path)
        if (identical(kind, "deny")) {
            add_redacted(path, sclet_ai_payload_deny_reason(name, path))
            return(NULL)
        }
        if (identical(kind, "metadata")) {
            add_redacted(path, "unconfirmed metadata values", consent = TRUE)
            return(NULL)
        }
        if (identical(kind, "unknown") &&
            depth >= 2L &&
            (identical(parent_kind, "aggregate") || identical(parent_kind, "container")) &&
            (is.numeric(value) || is.logical(value) || is.null(value) || is.list(value) ||
                (is.character(value) && length(value) <= 12L))) {
            kind <- "aggregate"
        }

        if (is.null(value)) {
            if (identical(kind, "unknown") && nzchar(path)) {
                return(unknown_field(path))
            }
            add_allowed(path)
            return(NULL)
        }
        if (is.matrix(value) || inherits(value, "Matrix") ||
            inherits(value, "DelayedMatrix") || isS4(value)) {
            add_redacted(path, "matrix or object values are never outbound")
            return(NULL)
        }
        if (is.environment(value) || is.function(value) || is.language(value)) {
            add_redacted(path, "unsupported value type")
            return(NULL)
        }
        if (is.data.frame(value)) {
            if (identical(kind, "unknown")) {
                return(unknown_field(path))
            }
            columns <- names(value)
            safe_columns <- columns[
                !vapply(
                    columns,
                    function(column) identical(
                        sclet_ai_payload_field_kind(column, path),
                        "deny"
                    ),
                    logical(1L)
                )
            ]
            result <- list(
                n_rows = as.integer(nrow(value)),
                n_columns = as.integer(ncol(value)),
                columns = safe_columns
            )
            add_allowed(path)
            return(result)
        }
        if (is.list(value)) {
            if (identical(kind, "unknown") && !isTRUE(user_consent)) {
                return(unknown_field(path))
            }
            if (is.null(names(value)) || any(!nzchar(names(value)))) {
                if (length(value) == 0L && !identical(kind, "unknown")) {
                    add_allowed(path)
                    return(value)
                }
                if (identical(kind, "unknown")) {
                    return(unknown_field(path))
                }
                add_redacted(path, "unnamed nested values are not bounded")
                return(NULL)
            }
            result <- lapply(seq_along(value), function(index) {
                child_name <- names(value)[[index]]
                child_path <- sclet_ai_payload_path(path, child_name)
                visit(
                    value[[index]],
                    child_path,
                    child_name,
                    depth + 1L,
                    parent_kind = kind
                )
            })
            names(result) <- names(value)
            keep <- !vapply(result, is.null, logical(1L))
            result <- result[keep]
            if (!length(result)) {
                return(NULL)
            }
            add_allowed(path)
            result
        } else if (is.atomic(value)) {
            if (identical(kind, "unknown") &&
                !(identical(privacy, "local") && isTRUE(user_consent))) {
                return(unknown_field(path))
            }
            result <- sclet_ai_payload_atomic(
                value = value,
                kind = kind,
                name = name,
                path = path,
                add_redacted = add_redacted,
                add_allowed = add_allowed
            )
            result
        } else {
            add_redacted(path, "unsupported value type")
            NULL
        }
    }

    if (!is.list(context)) {
        context <- list(context = context)
    }
    payload <- visit(context, path = "", name = "", depth = 0L)
    if (is.null(payload)) {
        payload <- list()
    }
    list(
        payload = payload,
        allowed_fields = state$allowed_fields,
        redacted_fields = state$redacted_fields,
        warnings = state$warnings,
        requires_user_consent = state$requires_user_consent
    )
}

sclet_ai_payload_atomic <- function(
    value,
    kind,
    name,
    path,
    add_redacted,
    add_allowed
) {
    if (identical(kind, "group")) {
        if (all(grepl("^group_[0-9]+$", as.character(value)))) {
            result <- unname(as.character(value))
            add_allowed(path)
            return(result)
        }
        counts <- table(as.character(value), useNA = "no")
        if (!length(counts) || any(counts < 10L)) {
            add_redacted(path, "group labels do not meet the aggregation threshold")
            return(NULL)
        }
        labels <- setNames(
            paste0("group_", seq_along(counts)),
            names(counts)
        )
        result <- unname(unname(labels[as.character(value)]))
        add_allowed(path)
        return(result)
    }

    if (sclet_ai_payload_is_count_field(name) &&
        is.numeric(value) && length(value)) {
        keep <- !is.na(value) & value >= 10
        if (!any(keep)) {
            add_redacted(path, "aggregate count is below the privacy threshold")
            return(NULL)
        }
        if (any(!keep)) {
            add_redacted(path, "some aggregate counts are below the privacy threshold")
        }
        result <- unname(value[keep])
        add_allowed(path)
        return(result)
    }

    if (is.character(value)) {
        if (sclet_ai_payload_contains_secret_or_path(value)) {
            add_redacted(path, "secret or private path detected")
            return(NULL)
        }
        if (identical(kind, "text")) {
            result <- sclet_ai_payload_bound_text(value)
            if (!length(result) || !nzchar(result[[1L]])) {
                add_redacted(path, "empty or unsafe text")
                return(NULL)
            }
            add_allowed(path)
            return(result)
        }
        if (identical(kind, "name") || identical(kind, "aggregate") ||
            identical(kind, "local")) {
            result <- sclet_ai_payload_bound_character(value)
            add_allowed(path)
            return(result)
        }
        add_redacted(path, "raw character values are not outbound")
        return(NULL)
    }

    if (is.factor(value)) {
        add_redacted(path, "factor levels may contain raw identifiers")
        return(NULL)
    }
    if (length(value) > 32L) {
        result <- list(
            kind = "bounded_vector",
            class = class(value),
            length = as.integer(length(value))
        )
    } else {
        result <- unname(value)
    }
    add_allowed(path)
    result
}

sclet_ai_payload_field_kind <- function(name, path) {
    field <- tolower(if (nzchar(name)) name else path)
    if (!nzchar(field)) {
        return("container")
    }
    if (grepl(
        "token|api[_-]?key|secret|password|credential|authorization|private[_-]?key|access[_-]?key|prompt|instruction|local[_-]?summary|(^|[_-])path($|[_-])|(^|[_-])file($|[_-])|barcode|cell[_-]?(id|name)|gene[_-]?(id|name)|patient|subject|donor|sample[_-]?id|clinical|note|comment",
        field,
        perl = TRUE
    )) {
        return("deny")
    }
    if (grepl("coldata|metadata[_-]?values|raw[_-]?metadata", field, perl = TRUE)) {
        return("metadata")
    }
    if (grepl("group|cluster", field, perl = TRUE) &&
        grepl("label|labels|groups?|clusters?|classes?", field, perl = TRUE) &&
        !grepl("count|size|number|composition", field, perl = TRUE)) {
        return("group")
    }
    if (grepl("goal|question|objective|user[_-]?request", field, perl = TRUE)) {
        return("text")
    }
    if (grepl(
        "assay|layer|reduction|graph|modality|feature[_-]?name|software|version|method|algorithm|parameter|(^|[_-])params?($|[_-])|status|state|schema|available|dimension|(^|[_-])class($|[_-])|summary|aggregate|count|number|composition|qc|diagnostic|profile|ledger|dataset|design|semantic|candidate|candidates|column|requested|scope|parents|claim[_-]?level|evidence|kind|source|uncertainty|privacy|execution|actions?[_-]?executed|capabilit|structure|cost|threshold|sparse|values?[_-]?included|raw[_-]?values?[_-]?included|n[_-]?(cells?|genes?|samples?|features?|clusters?|groups?)|(^|[_-])(active[_-]?view|active[_-]?states|analyses|health|state[_-]?records|workflows|lineage|blocked[_-]?actions|quality[_-]?checks|warnings|fingerprint|capabilities|analysis[_-]?story|analysis[_-]?type|timeline|n[_-]?steps|unordered[_-]?steps|user[_-]?decisions|conflicts|depends[_-]?on|created[_-]?at|confirmed[_-]?at|reason|gap|key|keys)($|[_-])",
        field,
        perl = TRUE
    )) {
        return("aggregate")
    }
    if (grepl("^names?$|^id$|^label$", field, perl = TRUE)) {
        return("name")
    }
    "unknown"
}

sclet_ai_payload_deny_reason <- function(name, path) {
    field <- tolower(if (nzchar(name)) name else path)
    if (grepl("token|api[_-]?key|secret|password|credential|authorization|private[_-]?key|access[_-]?key", field, perl = TRUE)) {
        return("credential or secret field")
    }
    if (grepl("prompt|instruction", field, perl = TRUE)) {
        return("prompt or instruction text")
    }
    if (grepl("path|file", field, perl = TRUE)) {
        return("local path or file reference")
    }
    if (grepl("barcode|cell|gene|patient|subject|donor|sample", field, perl = TRUE)) {
        return("raw identifier field")
    }
    "sensitive field"
}

sclet_ai_payload_is_count_field <- function(name) {
    grepl(
        "group|cluster|celltype|sample[_-]?(composition|counts?)|condition",
        tolower(name),
        perl = TRUE
    ) && grepl("count|size|number|n[_-]?|composition", tolower(name), perl = TRUE)
}

sclet_ai_payload_path <- function(parent, child) {
    if (!nzchar(parent)) child else paste(parent, child, sep = ".")
}

sclet_ai_payload_contains_secret_or_path <- function(value) {
    any(grepl(
        "(?i)(DEEPSEEK_API_KEY|[A-Z][A-Z0-9_]*(API[_-]?KEY|TOKEN|SECRET))\\s*[:=]\\s*[^[:space:]]+|bearer\\s+[A-Za-z0-9._-]+|(^|[[:space:]])(~[/\\\\]|/home/|/Users/|/tmp/|[A-Za-z]:[/\\\\])",
        as.character(value),
        perl = TRUE
    ))
}

sclet_ai_payload_bound_text <- function(value) {
    value <- as.character(value)[[1L]]
    value <- gsub("(?i)(DEEPSEEK_API_KEY|[A-Z][A-Z0-9_]*(API[_-]?KEY|TOKEN|SECRET))\\s*[:=]\\s*[^[:space:]]+", "[REDACTED]", value, perl = TRUE)
    value <- gsub("(^|[[:space:]])(~[/\\\\]|/home/|/Users/|/tmp/|[A-Za-z]:[/\\\\])[^[:space:]]+", "\\1[REDACTED_PATH]", value, perl = TRUE)
    substr(value, 1L, 500L)
}

sclet_ai_payload_bound_character <- function(value) {
    value <- as.character(value)
    value <- substr(value, 1L, 120L)
    unname(value)
}

sclet_ai_payload_canonical <- function(value) {
    if (is.list(value)) {
        if (!length(value)) {
            return(list())
        }
        if (!is.null(names(value))) {
            value <- value[order(names(value))]
            names(value) <- names(value)[order(names(value))]
        }
        return(lapply(value, sclet_ai_payload_canonical))
    }
    if (is.atomic(value)) {
        return(unname(value))
    }
    value
}
