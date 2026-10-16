#' Confirm the design-role semantics of colData columns
#'
#' Registers the user's own declaration of which colData column plays which
#' design role (batch, condition, subject, sample). Actions whose meaning
#' depends on the design -- run_integration, run_milo -- refuse to run until
#' this is called, and the confirmation is void when the declared column's
#' values are later reassigned. Roles are never inferred from column names:
#' a column named \code{sample} is not necessarily the sampling unit.
#'
#' This must be exported because it is part of the external-agent path: an
#' agent driving sclet without aisdk has to make this declaration itself
#' before any design-gated action can validate. See
#' \code{inst/skills/sclet-driver/SKILL.md}.
#'
#' @param object A SingleCellExperiment.
#' @param design A named list mapping roles (batch/condition/subject/sample)
#'   to colData column names, e.g. \code{list(batch = "batch_id")}.
#' @return The object with the confirmation recorded on its ledger.
#' @export
ConfirmAIDesignSemantics <- function(object, design) {
    if (!inherits(object, "SingleCellExperiment")) {
        stop("`object` must be a SingleCellExperiment.")
    }
    if (!is.list(design) || is.null(names(design))) {
        stop("`design` must be a named list mapping roles (batch/condition/subject/sample) to colData column names.")
    }
    known_roles <- c("batch", "condition", "subject", "sample")
    roles <- names(design)
    if (is.null(roles) || !all(nzchar(roles))) {
        stop("`design` must be a named list with non-empty names.")
    }
    unknown_roles <- setdiff(roles, known_roles)
    if (length(unknown_roles)) {
        stop(
            "Unsupported design role(s): ",
            paste(unknown_roles, collapse = ", "),
            ". Supported roles are: batch, condition, subject, sample."
        )
    }
    if (!length(roles)) {
        stop("`design` must map at least one role to a colData column name.")
    }
    cd <- SummarizedExperiment::colData(object)
    cd_names <- colnames(cd)
    for (role in roles) {
        col <- design[[role]]
        if (!is.character(col) || length(col) != 1L || is.na(col) || !nzchar(col)) {
            stop("Design role '", role, "' must be a single non-empty column name (character scalar).")
        }
        if (!col %in% cd_names) {
            stop("Design role '", role, "' points to column '", col, "' which is not present in colData(object). Available columns: ", paste(cd_names, collapse = ", "))
        }
    }
    design_sorted <- design[order(roles)]
    id <- paste(
        c("design", paste(sprintf("%s=%s", names(design_sorted), design_sorted), collapse = ",")),
        collapse = ";"
    )
    value_key <- sclet_ai_design_value_key(object, design_sorted)
    design_fp <- as.character(sclet_ai_design_fingerprint(object))
    ledger_fp_pre <- as.character(GetAnalysisLedger(object, detail = "summary", include_artifacts = FALSE, include_data = FALSE)$fingerprint)
    confirmed_at <- as.character(Sys.time())
    sclet_set_analysis_state(
        object = object,
        type = "ai_design_confirmation",
        id = id,
        method = "ConfirmAIDesignSemantics",
        inputs = design_sorted,
        summary = list(
            design = design_sorted,
            design_value_key = value_key,
            object_fingerprint = design_fp,
            ledger_fingerprint_before = ledger_fp_pre,
            confirmed_at = confirmed_at
        ),
        active = FALSE
    )
}

sclet_ai_design_value_key <- function(object, design_role_map) {
    cols <- unique(as.character(unlist(design_role_map)))
    cols <- cols[nzchar(cols)]
    if (!length(cols)) {
        return(NULL)
    }
    cd <- tryCatch(SummarizedExperiment::colData(object), error = function(e) NULL)
    if (is.null(cd) || !all(cols %in% colnames(cd))) {
        return(NULL)
    }
    hashes <- vapply(cols, function(col) {
        v <- cd[[col]]
        if (is.factor(v)) {
            paste(levels(v), as.integer(v), collapse = "|", sep = ":")
        } else if (is.atomic(v)) {
            if (length(v) <= 200L) {
                paste(as.character(v), collapse = "\u0001")
            } else {
                val <- as.character(v)
                n <- length(val)
                paste(
                    c("n=", as.character(n), "|head=", paste(val[seq_len(10L)], collapse = ","), "|tail=", paste(val[seq.int(to = n, length.out = 10L)], collapse = ","), "|uniq=", paste(unique(val), collapse = ",")),
                    collapse = ""
                )
            }
        } else {
            paste(class(v), collapse = ",")
        }
    }, character(1L), USE.NAMES = TRUE)
    out <- as.list(hashes)
    names(out) <- cols
    out
}

sclet_ai_design_fingerprint <- function(object) {
    tryCatch({
        cd <- SummarizedExperiment::colData(object)
        sclet_ai_fingerprint(list(
            n_cells = as.integer(ncol(object)),
            n_features = as.integer(nrow(object)),
            cd_names = as.character(colnames(cd)),
            has_rowdata = !is.null(tryCatch(
                SummarizedExperiment::rowData(object),
                error = function(e) NULL
            ))
        ))
    }, error = function(e) NULL)
}

# Resolve a design confirmation from the action's point of view.
#
# Passing `design = list(role = column, ...)` checks that this object carries a
# confirmation whose design map COVERS every role the action needs, on this
# fingerprint and with the same per-column value signature. `role %in%
# names(confirmed_map)` rather than equality, so a confirmation made for a
# richer design (batch + condition + sample) still satisfies an action needing
# only two of them: the user's declaration was not scoped to one action, and
# re-asking after a different step confirmed more roles would be theatre. The
# value signature is checked over the CONFIRMED map, so the check is that the
# columns the user attested to still hold the values they attested to.
#
# The named-argument form below is kept for run_integration's existing call.
sclet_ai_find_design_confirmation <- function(object, design = NULL, batch = NULL) {
    if (!inherits(object, "SingleCellExperiment")) {
        return(NULL)
    }
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
        return(NULL)
    }
    if (is.null(design)) {
        if (is.null(batch)) {
            return(NULL)
        }
        design <- list(batch = batch)
    }
    if (!is.list(design) || is.null(names(design)) || !length(design) ||
            any(!nzchar(as.character(names(design) %||% character())))) {
        return(NULL)
    }
    wanted_roles <- names(design)
    records <- tryCatch(
        sclet_get_state_records(object, "ai_design_confirmation"),
        error = function(e) list()
    )
    if (!length(records)) {
        return(NULL)
    }
    design_fp <- as.character(sclet_ai_design_fingerprint(object))
    for (id in names(records)) {
        rec <- records[[id]]
        summary <- rec$summary %||% list()
        if (!identical(as.character(summary$object_fingerprint %||% ""), design_fp)) {
            next
        }
        design_role_map <- summary$design %||% list()
        if (!all(wanted_roles %in% names(design_role_map))) {
            next
        }
        # A role the action names must point at the same column the user
        # confirmed it against, so confirming batch = "donor" cannot be
        # satisfied by a later claim that batch = "lane".
        mismatched <- wanted_roles[!vapply(wanted_roles, function(role) {
            identical(as.character(design_role_map[[role]] %||% ""),
                as.character(design[[role]] %||% ""))
        }, logical(1L))]
        if (length(mismatched)) {
            next
        }
        if (!is.null(summary$design_value_key)) {
            cur_key <- sclet_ai_design_value_key(object, design_role_map)
            if (is.null(cur_key) || !identical(cur_key, summary$design_value_key)) {
                next
            }
        }
        return(rec)
    }
    NULL
}

# The roles a design confirmation currently covers, for error messages that
# need to say what is still missing rather than restate the whole schema.
sclet_ai_design_confirmation_roles <- function(confirmation) {
    if (is.null(confirmation)) return(character())
    as.character(names(confirmation$summary$design %||% list()))
}
