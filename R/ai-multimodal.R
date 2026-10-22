# Read-only multimodal readiness diagnostic and bounded modality summaries
#
# Everything in this file is deliberately privacy-withholding and deliberately
# non-directive. A modality exists on the object in two different senses, and
# the distinction matters:
#
#   ATTACHED  -- the altExp exists, so the data is on the object. Discovered by
#                has_multimodal_adt() / has_multimodal_atac().
#   REGISTERED -- an analysis-state record of type "multimodal" exists, so the
#                mainline can reason about the modality without reaching into
#                object internals. Written by register_multimodal_adt() /
#                register_multimodal_atac(), which only RECORD dimensions.
#
# Registration computes nothing. It performs no integration, no batch
# correction, no normalization and no modality-specific analysis, so nothing
# below may describe a registered modality as integrated or corrected. Nothing
# below may either: recommend an integration method, rank modalities, name a
# feature or antibody, or emit a per-cell value. AltExp names and the assay
# platform labels "ADT" / "ATAC" are assay-platform names rather than data, so
# they may appear in a diagnostic; feature and antibody identity may not.

# Collect the registered multimodal state records. Falls back to an empty list
# rather than erroring on an object whose state layer is absent.
sclet_ai_multimodal_records <- function(object) {
    records <- tryCatch(sclet_get_state_records(object, "multimodal"), error = function(e) list())
    if (is.null(records)) records <- list()
    if (!is.list(records)) return(list())
    records[!vapply(records, is.null, logical(1L))]
}

# The analysis key a multimodal record was written under: "multimodal_adt",
# "multimodal_atac", or "multimodal_workflow" for the workflow record.
sclet_ai_multimodal_keys <- function(record) {
    if (!is.list(record)) return(character())
    keys <- as.character(record$artifacts$analysis_key %||% character())
    keys <- keys[!is.na(keys) & nzchar(keys)]
    # A record that carries no artifacts$analysis_key is still readable: its
    # method field was written by the same registration call.
    method <- as.character(record$method %||% "")
    if (!length(keys) && length(method) == 1L && !is.na(method) && nzchar(method)) {
        keys <- switch(method,
            altExp_ADT = "multimodal_adt",
            altExp_ATAC = "multimodal_atac",
            multimodal_mainline = "multimodal_workflow",
            character()
        )
    }
    keys
}

# The first record whose analysis key is `key`, or NULL. Search order is the
# record list order so the result is stable across calls.
sclet_ai_multimodal_record_for <- function(records, key) {
    for (record in records) {
        if (key %in% sclet_ai_multimodal_keys(record)) return(record)
    }
    NULL
}

# The workflow record written by RunMultimodalWorkflow(). Note carefully WHERE
# it lives: unlike the two registration functions, the workflow writes ONLY the
# analysis layer via sclet_set_analysis(), under the key "multimodal_workflow".
# It writes no state record of its own, so it must be read from the analyses
# rather than from sclet_get_state_records(). It is the only record that stores
# has_ADT / has_ATAC. The state-record scan is kept as a fallback so a future
# workflow that also registers a state record is picked up either way.
sclet_ai_multimodal_workflow_record <- function(object, records = NULL) {
    record <- tryCatch(sclet_get_analysis(object, "multimodal_workflow"), error = function(e) NULL)
    if (is.null(record) && !is.null(records)) {
        record <- sclet_ai_multimodal_record_for(records, "multimodal_workflow")
    }
    if (!is.list(record)) return(NULL)
    record
}

# TRUE when at least one multimodal record is registered, of any kind.
sclet_ai_multimodal_has_records <- function(records) {
    length(records) > 0L
}

# An integer count from a record field, or NA_integer_ when the record does not
# carry one. NA rather than 0 so "not recorded" never reads as "zero features".
sclet_ai_multimodal_count <- function(value) {
    if (is.null(value)) return(NA_integer_)
    value <- suppressWarnings(as.numeric(value))
    value <- value[is.finite(value)]
    if (!length(value)) return(NA_integer_)
    as.integer(value[[1L]])
}

# The registered ADT / ATAC record, preferring the one keyed by analysis_key.
sclet_ai_multimodal_registered <- function(records, modality) {
    key <- paste0("multimodal_", tolower(modality))
    !is.null(sclet_ai_multimodal_record_for(records, key))
}

# Platform labels, in the fixed ADT-then-ATAC order. "ADT" and "ATAC" name an
# assay platform, not a feature, so they are safe to return; ordering is fixed
# here rather than derived, so it cannot drift with record-insertion order.
sclet_ai_multimodal_kinds <- function(records) {
    kinds <- character()
    if (sclet_ai_multimodal_registered(records, "ADT")) kinds <- c(kinds, "ADT")
    if (sclet_ai_multimodal_registered(records, "ATAC")) kinds <- c(kinds, "ATAC")
    kinds
}

# Anonymized positional code for a modality label: group_1, group_2, ... based
# on position in the fixed kind order. The label itself is never returned.
sclet_ai_multimodal_code <- function(kind, kinds) {
    index <- match(as.character(kind), as.character(kinds))
    if (is.na(index)) return(NULL)
    paste0("group_", index)
}

# The shared caveat block for the bounded summaries. Wordings are fixed here so
# no caller can drift into describing a registered modality as integrated, or a
# bigger modality as a better one.
sclet_ai_multimodal_caveats <- function() {
    list(
        paste(
            "Modality registration records what was ATTACHED to the object and tracked by the mainline; it performs",
            "no integration, no batch correction and no modality-specific analysis."
        ),
        paste(
            "The counts reported here describe the stored altExp dimensions, not biological content and not data",
            "quality; they say how large the stored modality is, not how informative it is."
        ),
        paste(
            "No per-cell value, cell barcode, antibody name, feature or gene symbol, file path or sample identifier",
            "is included in this summary."
        ),
        paste(
            "No integration method is recommended here: how modalities should be combined is a modeling decision",
            "reserved for the user."
        )
    )
}

#' Check whether an object is ready for multimodal diagnostics
#'
#' Reports which multimodal modalities are attached to the object as altExps,
#' which of them carry a multimodal analysis-state record, and whether a
#' multimodal workflow record is present. This diagnostic deliberately never
#' says which modalities an object *should* have, never recommends an
#' integration method, and never ranks or scores a modality: which modalities
#' are expected, and how they should be combined, is a modeling assumption only
#' the user can make.
#'
#' The diagnostic distinguishes two separate states for each modality. A
#' modality is **attached** when an altExp such as \code{"ADT"} or
#' \code{"ATAC"} exists on the object. It is **registered** when a
#' \code{"multimodal"} analysis-state record exists for it. Registration makes
#' the modality trackable by the mainline; it computes nothing. An object can
#' therefore be reported as attached without being registered, and the returned
#' questions say so explicitly rather than asking for data to be re-attached.
#'
#' @param object A \code{SingleCellExperiment} object.
#' @param modalities Optional character vector of altExp names the caller is
#'   asking about. Used only to report which of the requested names are present;
#'   no requested name is ever reported as missing-in-error, and no name outside
#'   the caller's own request enters the result.
#' @return A read-only readiness summary list with elements \code{status},
#'   \code{checks}, \code{blocked_actions}, \code{questions} and \code{notes}.
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' readiness <- check_multimodal_readiness(sce)
#' readiness$status
#' readiness$questions
#' @export
check_multimodal_readiness <- function(object, modalities = NULL) {
    if (!inherits(object, "SingleCellExperiment")) {
        stop("object must be a SingleCellExperiment", call. = FALSE)
    }
    if (!is.null(modalities)) {
        if (!is.character(modalities) || any(is.na(modalities)) || any(!nzchar(modalities))) {
            stop("modalities must be NULL or a non-empty character vector", call. = FALSE)
        }
        modalities <- unique(modalities)
    }

    alt_exp_names <- tryCatch(
        as.character(SingleCellExperiment::altExpNames(object)),
        error = function(e) character()
    )
    if (is.null(alt_exp_names)) alt_exp_names <- character()

    adt_available <- tryCatch(has_multimodal_adt(object), error = function(e) FALSE)
    atac_available <- tryCatch(has_multimodal_atac(object), error = function(e) FALSE)
    adt_available <- isTRUE(adt_available)
    atac_available <- isTRUE(atac_available)
    n_modalities_available <- as.integer(sum(adt_available, atac_available))

    records <- sclet_ai_multimodal_records(object)
    adt_registered <- sclet_ai_multimodal_registered(records, "ADT")
    atac_registered <- sclet_ai_multimodal_registered(records, "ATAC")
    workflow_record_available <- !is.null(sclet_ai_multimodal_workflow_record(object, records))
    any_record <- sclet_ai_multimodal_has_records(records)

    # Requested modality names are narrowed to what actually exists, so the
    # diagnostic can report a requested-but-absent name without inventing one.
    requested_resolved <- NULL
    if (!is.null(modalities)) {
        requested_resolved <- alt_exp_names[alt_exp_names %in% modalities]
    }

    questions <- character()

    # Branch 1: the modality is NOT attached. The advice is to attach it.
    if (!adt_available) {
        questions <- c(questions, paste0(
            "No ADT modality is attached: no altExp named \"ADT\" exists on this object, so no ADT record can be ",
            "registered. Attach the modality first with altExp(sce, \"ADT\") <- adt_sce. ",
            "Available altExps: ", paste(alt_exp_names, collapse = ", ")
        ))
    }
    if (!atac_available) {
        questions <- c(questions, paste0(
            "No ATAC modality is attached: no altExp named \"ATAC\" exists on this object, so no ATAC record can be ",
            "registered. Attach the modality first with altExp(sce, \"ATAC\") <- atac_sce. ",
            "Available altExps: ", paste(alt_exp_names, collapse = ", ")
        ))
    }

    # Branch 2: the modality IS attached but NOT registered. The advice is
    # explicitly NOT to attach anything again -- it is to register, because
    # registration is what makes the already-present data trackable.
    if (adt_available && !adt_registered) {
        questions <- c(questions, paste0(
            "An ADT modality is attached (altExp \"ADT\" exists) but is not registered: no multimodal analysis-state ",
            "record describes it, so the mainline cannot track it. Registration is what makes the attached data ",
            "trackable, and it only records the stored dimensions -- it computes no integration or batch correction. ",
            "Register it with register_multimodal_adt(sce); do not re-attach data that is already present."
        ))
    }
    if (atac_available && !atac_registered) {
        questions <- c(questions, paste0(
            "An ATAC modality is attached (altExp \"ATAC\" exists) but is not registered: no multimodal analysis-state ",
            "record describes it, so the mainline cannot track it. Registration is what makes the attached data ",
            "trackable, and it only records the stored dimensions -- it computes no integration or batch correction. ",
            "Register it with register_multimodal_atac(sce); do not re-attach data that is already present."
        ))
    }

    # A caller-supplied request that matches no altExp is reported as such and
    # is never turned into a recommendation about what should exist.
    if (!is.null(requested_resolved) && length(requested_resolved) < length(modalities)) {
        missing_requested <- setdiff(modalities, requested_resolved)
        questions <- c(questions, paste0(
            "The requested modality name(s) ", paste(missing_requested, collapse = ", "),
            " do not match any altExp on this object. Available altExps: ",
            paste(alt_exp_names, collapse = ", "),
            ". No modality is recommended here; this reports only whether the name you asked about is present."
        ))
    }

    # Neither attached nor registered: genuinely nothing multimodal to work with.
    if (n_modalities_available == 0L && !any_record) {
        questions <- c(questions, paste(
            "No multimodal modality is attached and no multimodal state record exists, so there is nothing for a",
            "multimodal diagnostic to read. Which modalities an object should carry is a question reserved for the",
            "user; nothing is suggested or recommended here."
        ))
    }

    status <- if (n_modalities_available == 0L && !any_record) {
        "not_ready"
    } else {
        "ready_for_diagnostic"
    }

    list(
        status = status,
        checks = list(
            alt_exp_names = alt_exp_names,
            requested_modalities = modalities,
            requested_modalities_present = requested_resolved,
            adt_available = adt_available,
            atac_available = atac_available,
            n_modalities_available = n_modalities_available,
            adt_registered = adt_registered,
            atac_registered = atac_registered,
            workflow_record_available = workflow_record_available,
            raw_values_included = FALSE
        ),
        blocked_actions = if (identical(status, "not_ready")) "multimodal" else character(),
        questions = questions,
        notes = paste(
            "Attached and registered are different states. Attached means an altExp such as \"ADT\" or \"ATAC\" exists",
            "on the object; registered means a \"multimodal\" state record makes it trackable by the mainline.",
            "Registration only records what was attached and its stored dimensions -- it performs no integration, no",
            "batch correction and no modality-specific analysis. No modality is ranked or preferred, no downstream",
            "multimodal method is recommended, and no cell is removed or filtered by this diagnostic."
        )
    )
}

#' Summarize the registered multimodal modalities
#'
#' Returns a bounded, read-only summary of the registered multimodal state
#' records. Only counts, booleans, plain assay-platform labels and aggregate
#' statements are emitted. No cell barcode, antibody name, feature or gene
#' symbol, file path, sample identifier or per-cell vector is ever returned, and
#' no integration method is recommended. The function does not mutate the
#' object and writes no evidence node.
#'
#' @param object A \code{SingleCellExperiment} object.
#' @param name Optional multimodal record id. When supplied, only that record is
#'   summarized; when \code{NULL}, every multimodal record is summarized.
#' @return A read-only summary list with elements \code{status}, and when
#'   available \code{n_modalities}, \code{modality_kinds}, \code{modalities}
#'   (per-kind feature and cell counts), \code{has_ADT}, \code{has_ATAC},
#'   \code{claim_level}, \code{caveats} and \code{raw_values_included = FALSE}.
#'   When no modality is available and no record exists, \code{status =
#'   "not_available"} with \code{reason = "no_multimodal_record"}.
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' # No multimodal record is on record, so nothing is reported.
#' summary <- summarize_multimodal_evidence(sce)
#' summary$status
#' summary$reason
#' @export
summarize_multimodal_evidence <- function(object, name = NULL) {
    if (!inherits(object, "SingleCellExperiment")) {
        stop("object must be a SingleCellExperiment", call. = FALSE)
    }
    if (!is.null(name) && (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name))) {
        stop("name must be NULL or a single non-empty character value", call. = FALSE)
    }

    adt_available <- isTRUE(tryCatch(has_multimodal_adt(object), error = function(e) FALSE))
    atac_available <- isTRUE(tryCatch(has_multimodal_atac(object), error = function(e) FALSE))

    records <- sclet_ai_multimodal_records(object)
    if (!is.null(name)) {
        # An explicitly requested id narrows the record set; an unknown id
        # leaves it empty rather than silently falling back to every record.
        records <- records[names(records) %in% as.character(name)]
    }

    adt_record <- sclet_ai_multimodal_record_for(records, "multimodal_adt")
    atac_record <- sclet_ai_multimodal_record_for(records, "multimodal_atac")
    workflow_record <- sclet_ai_multimodal_workflow_record(object, records)

    has_record <- !is.null(adt_record) || !is.null(atac_record) || !is.null(workflow_record)
    if ((!adt_available && !atac_available && !has_record)) {
        return(list(
            status = "not_available",
            reason = "no_multimodal_record",
            n_modalities = 0L,
            modality_kinds = character(),
            modalities = list(),
            claim_level = "consistent_with",
            caveats = sclet_ai_multimodal_caveats(),
            recommendation = NULL,
            raw_values_included = FALSE
        ))
    }

    # Booleans come from the workflow record when it exists, since that is the
    # record which stores has_ADT / has_ATAC; otherwise they are derived from
    # the per-modality records that do exist.
    has_adt <- if (!is.null(workflow_record)) {
        isTRUE(workflow_record$summary$has_ADT)
    } else {
        !is.null(adt_record)
    }
    has_atac <- if (!is.null(workflow_record)) {
        isTRUE(workflow_record$summary$has_ATAC)
    } else {
        !is.null(atac_record)
    }

    modality_kinds <- character()
    if (has_adt || !is.null(adt_record)) modality_kinds <- c(modality_kinds, "ADT")
    if (has_atac || !is.null(atac_record)) modality_kinds <- c(modality_kinds, "ATAC")

    counts_for <- function(record) {
        if (is.null(record)) return(list(n_features = NA_integer_, n_cells = NA_integer_))
        list(
            n_features = sclet_ai_multimodal_count(record$artifacts$n_features),
            n_cells = sclet_ai_multimodal_count(record$artifacts$n_cells)
        )
    }
    per_kind <- list(
        ADT = counts_for(adt_record),
        ATAC = counts_for(atac_record)
    )
    per_kind <- per_kind[modality_kinds]
    modalities <- lapply(per_kind, function(entry) {
        list(
            n_features = entry$n_features,
            n_cells = entry$n_cells,
            raw_values_included = FALSE
        )
    })

    list(
        status = "available",
        n_modalities = as.integer(length(modality_kinds)),
        modality_kinds = modality_kinds,
        modalities = modalities,
        has_ADT = has_adt,
        has_ATAC = has_atac,
        workflow_record_available = !is.null(workflow_record),
        claim_level = "consistent_with",
        caveats = sclet_ai_multimodal_caveats(),
        recommendation = NULL,
        raw_values_included = FALSE
    )
}

#' Compare registered multimodal modalities
#'
#' Observationally compares the registered multimodal modalities. The function
#' reports how many modalities are registered, how many features and cells each
#' one stores, and nothing else. It never selects, ranks or prefers a modality,
#' never describes one as the informative or dominant modality, and never
#' recommends an integration method. Because different assay platforms have
#' different feature spaces, the counts compare stored dimensions rather than
#' biological signal.
#'
#' @param object A \code{SingleCellExperiment} object.
#' @return A read-only comparison list with elements \code{status}, and when
#'   available \code{n_modalities}, \code{modality_codes} (anonymized
#'   \code{group_N} codes only), a per-code count table, \code{claim_level},
#'   \code{caveats}, \code{recommendation = NULL} and
#'   \code{raw_values_included = FALSE}. Fewer than two registered modalities
#'   return \code{status = "not_available"} with
#'   \code{reason = "at_least_two_modalities_required"}.
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' comparison <- compare_multimodal_evidence(sce)
#' comparison$status
#' comparison$reason
#' @export
compare_multimodal_evidence <- function(object) {
    if (!inherits(object, "SingleCellExperiment")) {
        stop("object must be a SingleCellExperiment", call. = FALSE)
    }

    records <- sclet_ai_multimodal_records(object)
    adt_record <- sclet_ai_multimodal_record_for(records, "multimodal_adt")
    atac_record <- sclet_ai_multimodal_record_for(records, "multimodal_atac")

    # Only registered per-modality records count: the workflow record is a
    # container over them, not a modality of its own.
    kinds <- character()
    if (!is.null(adt_record)) kinds <- c(kinds, "ADT")
    if (!is.null(atac_record)) kinds <- c(kinds, "ATAC")

    insufficient <- function(n) {
        list(
            status = "not_available",
            reason = "at_least_two_modalities_required",
            n_modalities = as.integer(n),
            modality_codes = character(),
            modalities = list(),
            claim_level = "consistent_with",
            caveats = c(sclet_ai_multimodal_caveats(), paste(
                "At least two registered modalities are required for an across-modality comparison; none is selected,",
                "ranked or preferred when they are missing."
            )),
            recommendation = NULL,
            raw_values_included = FALSE
        )
    }
    if (length(kinds) < 2L) {
        return(insufficient(length(kinds)))
    }

    counts_for <- function(record) {
        list(
            n_features = sclet_ai_multimodal_count(record$artifacts$n_features),
            n_cells = sclet_ai_multimodal_count(record$artifacts$n_cells)
        )
    }
    per_kind <- list(ADT = counts_for(adt_record), ATAC = counts_for(atac_record))
    per_kind <- per_kind[kinds]
    codes <- vapply(kinds, function(kind) sclet_ai_multimodal_code(kind, kinds), character(1L))
    names(per_kind) <- codes

    modalities <- lapply(per_kind, function(entry) {
        list(
            n_features = entry$n_features,
            n_cells = entry$n_cells,
            raw_values_included = FALSE
        )
    })

    list(
        status = "available",
        n_modalities = as.integer(length(kinds)),
        modality_codes = unname(codes),
        modalities = modalities,
        claim_level = "consistent_with",
        caveats = c(sclet_ai_multimodal_caveats(), paste(
            "This compares the stored DIMENSIONS of different assay platforms rather than biological signal: a larger",
            "modality is not a better one, a modality with more features is not more informative, and no modality is",
            "preferred, ranked or recommended here."
        )),
        recommendation = NULL,
        raw_values_included = FALSE
    )
}
