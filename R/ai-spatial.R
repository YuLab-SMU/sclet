# Read-only spatial readiness, bounded evidence summaries and analysis comparison
#
# Everything in this file is deliberately privacy-withholding and purely
# observational. Spatial state records carry deconvolution-derived feature
# column names, per-spot model estimates and file-level inputs; none of the
# functions below ever let a cell/spot barcode, a sample identifier, a file
# path, a gene symbol, a raw matrix or any per-cell/per-spot vector leave the
# package boundary. Only counts, algorithm names, anonymized codes and fixed
# caveats are emitted. None of these functions records AI evidence -- they read
# existing state records and never call RecordAIEvidence().

# The analysis keys the spatial mainline writes into a state record's
# artifacts$analysis_key. This is the vocabulary every entry point in this file
# resolves against; it is never used to rank or prefer an analysis.
sclet_ai_spatial_keys <- function() {
    c("spatial_deconv", "colocalization", "spatial_niche")
}

# All spatial state records on the object, keyed by record id. A record that is
# not a list is dropped rather than coerced, so a malformed record can never be
# read as if it were a completed analysis.
sclet_ai_spatial_records <- function(object) {
    records <- tryCatch(sclet_get_state_records(object, "spatial"), error = function(e) list())
    if (is.null(records)) records <- list()
    records[vapply(records, is.list, logical(1L))]
}

# The record's analysis key, normalized to "" when the record does not declare
# one. Records written before artifacts carried the key fall back to the method
# only to classify the record -- never to identify a data column.
sclet_ai_spatial_key <- function(record) {
    if (is.null(record) || !is.list(record)) return("")
    key <- as.character(record$artifacts$analysis_key %||% "")
    key <- key[!is.na(key) & nzchar(key)]
    if (length(key)) return(key[[1L]])
    method <- tolower(as.character(record$method %||% ""))
    if (identical(method, "cell2location")) return("spatial_deconv")
    if (identical(method, "svp_globalbv")) return("colocalization")
    if (identical(method, "svp_lisa")) return("spatial_niche")
    ""
}

# Resolve the requested analysis against the records actually present. Returns
# a stable, key-ordered list keyed by analysis key: multiple records for the
# same analysis are kept as a list under that key, and the order of the keys is
# fixed by sclet_ai_spatial_keys() so downstream codes never shuffle.
sclet_ai_spatial_by_key <- function(object, analysis = NULL) {
    if (!is.null(analysis) && (!is.character(analysis) || length(analysis) != 1L || !nzchar(analysis))) {
        stop("analysis must be NULL or a single non-empty character value", call. = FALSE)
    }
    records <- sclet_ai_spatial_records(object)
    by_key <- list()
    for (record_id in names(records)) {
        key <- sclet_ai_spatial_key(records[[record_id]])
        if (!nzchar(key)) next
        by_key[[key]] <- c(by_key[[key]] %||% list(), stats::setNames(list(records[[record_id]]), record_id))
    }
    present <- sclet_ai_spatial_keys()
    present <- present[present %in% names(by_key)]
    by_key <- by_key[present]
    if (!is.null(analysis)) {
        by_key <- by_key[names(by_key) %in% as.character(analysis)]
    }
    by_key
}

# Does the modality expose spatial coordinates? SpatialExperiment stores them
# outside the SCE colData/assay slots and this package does not depend on that
# class, so the check is deliberately defensive: it looks for coordinates
# wherever the object actually puts them. Only the logical answer escapes --
# the coordinate values themselves are never read into a return value.
sclet_ai_spatial_coordinates <- function(object) {
    if (inherits(object, "SpatialExperiment")) return(TRUE)
    ## SpatialExperiment is an optional, undeclared dependency on purpose: the
    ## class and its accessor are probed only when actually installed, so the
    ## "::" call must sit behind requireNamespace() or R CMD check reports an
    ## undeclared import.
    coordinates <- NULL
    if (requireNamespace("SpatialExperiment", quietly = TRUE)) {
        coordinates <- tryCatch(SpatialExperiment::spatialCoords(object), error = function(e) NULL)
    }
    if (is.matrix(coordinates) || is.data.frame(coordinates)) return(TRUE)
    reductions <- tryCatch(SingleCellExperiment::reducedDimNames(object), error = function(e) character())
    if (any(grepl("spatial", reductions, ignore.case = TRUE))) return(TRUE)
    columns <- tryCatch(colnames(SummarizedExperiment::colData(object)), error = function(e) character())
    coordinate_names <- c("spatial_x", "spatial_y", "spatial_coords", "x", "y", "row", "col",
        "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres")
    any(tolower(columns) %in% coordinate_names)
}

# The counts assay gate shared by the readiness diagnostic. Spatial backends in
# this package consume a raw "counts" assay; the check reports availability only.
sclet_ai_spatial_has_counts <- function(object) {
    assays <- tryCatch(SummarizedExperiment::assayNames(object), error = function(e) character())
    "counts" %in% assays
}

# Bounded counts for one record. Every value is derived from stored sizes or
# from the DIMENSIONS of a stored artifact -- never from an artifact's contents.
# `n_cells` is the object's spot/cell count, not a per-spot vector.
sclet_ai_spatial_counts <- function(object, record) {
    artifacts <- record$artifacts %||% list()
    summary <- record$summary %||% list()
    features <- artifacts$features %||% record$features
    n_features <- length(features)
    if (!n_features) n_features <- summary$n_features %||% artifacts$n_features %||% 0L
    columns <- artifacts$columns %||% artifacts$features %||% record$columns
    n_columns <- length(columns)
    if (!n_columns) n_columns <- summary$n_columns %||% artifacts$n_columns %||% 0L
    if (!n_columns && !is.null(summary$n_features)) n_columns <- summary$n_features
    list(
        n_features = as.integer(n_features %||% 0L),
        n_columns = as.integer(n_columns %||% 0L),
        n_cells = as.integer(ncol(object))
    )
}

# The method identifier as a plain algorithm name ("cell2location", "SVP_LISA",
# ...). Algorithm names are not data, so they are the one string this file may
# echo verbatim. Anything longer than a bare identifier is dropped rather than
# risk carrying an embedded path or file name.
sclet_ai_spatial_method <- function(record) {
    method <- as.character(record$method %||% "")
    method <- method[!is.na(method) & nzchar(method)]
    if (!length(method)) return(NULL)
    method <- method[[1L]]
    if (grepl("[\\\\/]|\\.(h5|h5ad|rds|rda|csv|tsv|mtx|loom|zarr)", method, ignore.case = TRUE)) return(NULL)
    if (nchar(method) > 40L) return(NULL)
    method
}

# Deterministic anonymized code for one analysis key. The position comes from
# the fixed key order, so the same analysis always yields the same code and the
# code never encodes the record id, which is user-chosen and may be identifying.
# The `route_` prefix matches the code vocabulary used elsewhere in the package.
sclet_ai_spatial_code <- function(key, available_keys = NULL) {
    if (is.null(key) || is.null(available_keys) || !length(available_keys)) return(NULL)
    index <- match(as.character(key), as.character(available_keys))
    if (is.na(index)) return(NULL)
    paste0("route_", index)
}

# The shared caveat block. Wordings are fixed here so no caller can drift into
# describing a model estimate as a measurement or a comparison as a verdict.
sclet_ai_spatial_caveats <- function() {
    list(
        paste(
            "Spatial deconvolution proportions and niche assignments are model estimates, not measured",
            "cell-type composition and not ground-truth tissue architecture."
        ),
        paste(
            "Those estimates depend on the reference used, the chosen method and its parameters, so the same",
            "input analysed with a different reference or setting can yield different numbers."
        ),
        paste(
            "No per-cell or per-spot value is included here: only record-level counts, algorithm names and",
            "anonymized codes are returned."
        )
    )
}

#' Check whether an object is ready for spatial diagnostics
#'
#' Reports what spatial analyses are stored on the object and whether the
#' modality exposes the prerequisites a spatial record needs. This diagnostic
#' deliberately never recommends a spatial method: deconvolution,
#' colocalization and niche detection answer different questions, and which of
#' them a dataset needs is an analysis decision only the user can make. The
#' result therefore contains counts, availability flags and analysis keys only.
#'
#' @param object A \code{SingleCellExperiment} object.
#' @param analysis Optional analysis key to report on. One of
#'   \code{"spatial_deconv"}, \code{"colocalization"} or
#'   \code{"spatial_niche"}. When \code{NULL}, every analysis present is
#'   reported.
#' @return A read-only readiness summary list with \code{status},
#'   \code{checks}, \code{blocked_actions}, \code{questions} and \code{notes}.
#' @export
check_spatial_readiness <- function(object, analysis = NULL) {
    if (!inherits(object, "SingleCellExperiment")) {
        stop("object must be a SingleCellExperiment", call. = FALSE)
    }
    if (!is.null(analysis) && (!is.character(analysis) || length(analysis) != 1L || !nzchar(analysis))) {
        stop("analysis must be NULL or a single non-empty character value", call. = FALSE)
    }

    records <- sclet_ai_spatial_records(object)
    by_key <- sclet_ai_spatial_by_key(object)
    available_analyses <- names(by_key)
    n_spatial_records <- as.integer(length(records))
    record_available <- n_spatial_records > 0L
    has_coordinates <- sclet_ai_spatial_coordinates(object)
    has_counts <- sclet_ai_spatial_has_counts(object)

    requested <- if (is.null(analysis)) NULL else as.character(analysis)
    requested_available <- is.null(requested) || requested %in% available_analyses

    questions <- character()
    if (!record_available) {
        questions <- c(questions, paste(
            "No completed spatial analysis is stored on this object, so there is nothing to diagnose yet;",
            "a completed spatial analysis is required. The spatial entry points are",
            "RunSpatialDeconvolution(), RunSpatialColocalization() and RunSpatialNiche()."
        ))
    }
    if (!is.null(requested) && !requested_available) {
        questions <- c(questions, paste0(
            "The requested spatial analysis '", requested, "' is not available on this object. ",
            if (length(available_analyses)) {
                paste0("Available spatial analyses: ", paste(available_analyses, collapse = ", "), ".")
            } else {
                "No spatial analysis is stored on this object."
            }
        ))
    }
    if (!has_coordinates) {
        questions <- c(questions, paste(
            "No spatial coordinates were found in the modality: deconvolution proportions and niche",
            "assignments are only spatially interpretable when the spots carry their coordinates."
        ))
    }
    if (!has_counts) {
        questions <- c(questions, paste(
            "No 'counts' assay was found; the spatial backends in this package read a raw counts assay.",
            "Available assays: ", paste(tryCatch(SummarizedExperiment::assayNames(object),
                error = function(e) character()), collapse = ", ")
        ))
    }
    if ("spatial_deconv" %in% available_analyses) {
        questions <- c(questions, paste(
            "A spatial deconvolution record is present; downstream niche detection consumes its feature",
            "columns, so run niche detection on the same object or pass those columns explicitly."
        ))
    }

    list(
        status = if (record_available) "ready_for_diagnostic" else "not_ready",
        checks = list(
            spatial_record_available = record_available,
            available_analyses = available_analyses,
            n_spatial_records = n_spatial_records,
            has_spatial_coordinates = has_coordinates,
            has_counts_assay = has_counts,
            requested_analysis = requested,
            requested_analysis_available = requested_available,
            raw_values_included = FALSE
        ),
        blocked_actions = if (record_available) character() else "spatial",
        questions = questions,
        notes = paste(
            "No spatial method is recommended, ranked or preferred here: deconvolution, colocalization and",
            "niche detection answer different questions, and the choice belongs to the user. Spatial",
            "deconvolution proportions and niche assignments are model estimates rather than measured",
            "cell-type composition or ground-truth tissue architecture, and no cell is removed or filtered",
            "by this diagnostic."
        )
    )
}

#' Summarize the stored spatial analyses
#'
#' Reads the spatial state records and returns a bounded summary of what they
#' hold. Feature columns, per-spot values and input identifiers exist in the
#' underlying records but never leave this function: only record counts,
#' algorithm names, anonymized analysis codes and fixed caveats are returned.
#'
#' @param object A \code{SingleCellExperiment} object.
#' @param analysis Optional analysis key to summarize. One of
#'   \code{"spatial_deconv"}, \code{"colocalization"} or
#'   \code{"spatial_niche"}. When \code{NULL}, every analysis present is
#'   summarized.
#' @return A read-only summary list with \code{status} (\code{"available"} or
#'   \code{"not_available"}), \code{reason}, \code{n_analyses},
#'   \code{analyses} (per-analysis bounded counts keyed by anonymized codes),
#'   \code{claim_level}, \code{caveats} and \code{raw_values_included = FALSE}.
#' @export
summarize_spatial_evidence <- function(object, analysis = NULL) {
    if (!inherits(object, "SingleCellExperiment")) {
        stop("object must be a SingleCellExperiment", call. = FALSE)
    }
    if (!is.null(analysis) && (!is.character(analysis) || length(analysis) != 1L || !nzchar(analysis))) {
        stop("analysis must be NULL or a single non-empty character value", call. = FALSE)
    }

    all_records <- sclet_ai_spatial_records(object)
    if (!length(all_records)) {
        return(list(
            status = "not_available",
            reason = "no_spatial_record",
            n_analyses = 0L,
            analyses = list(),
            claim_level = "consistent_with",
            caveats = sclet_ai_spatial_caveats(),
            raw_values_included = FALSE
        ))
    }

    by_key <- sclet_ai_spatial_by_key(object, analysis = analysis)
    if (!length(by_key)) {
        return(list(
            status = "not_available",
            reason = "analysis_not_available",
            requested_analysis = if (is.null(analysis)) NULL else as.character(analysis),
            n_analyses = 0L,
            analyses = list(),
            claim_level = "consistent_with",
            caveats = sclet_ai_spatial_caveats(),
            raw_values_included = FALSE
        ))
    }

    available_keys <- names(sclet_ai_spatial_by_key(object))
    summaries <- lapply(names(by_key), function(key) {
        entries <- by_key[[key]]
        n_records <- as.integer(length(entries))
        counts <- lapply(entries, function(record) sclet_ai_spatial_counts(object, record))
        features <- vapply(counts, function(entry) entry$n_features, integer(1L))
        columns <- vapply(counts, function(entry) entry$n_columns, integer(1L))
        cells <- vapply(counts, function(entry) entry$n_cells, integer(1L))
        methods <- unique(vapply(entries, function(record) {
            as.character(sclet_ai_spatial_method(record) %||% NA_character_)
        }, character(1L)))
        methods <- methods[!is.na(methods) & nzchar(methods)]
        list(
            analysis_code = sclet_ai_spatial_code(key, available_keys),
            n_records = n_records,
            n_features = as.integer(max(features, 0L)),
            n_columns = as.integer(max(columns, 0L)),
            n_cells = as.integer(max(cells, 0L)),
            methods = as.character(methods),
            claim_level = "consistent_with",
            raw_values_included = FALSE
        )
    })
    names(summaries) <- vapply(summaries, function(entry) entry$analysis_code, character(1L))

    list(
        status = "available",
        n_analyses = as.integer(length(summaries)),
        analyses = summaries,
        claim_level = "consistent_with",
        caveats = sclet_ai_spatial_caveats(),
        raw_values_included = FALSE
    )
}

#' Compare the available spatial analyses
#'
#' Observational comparison across the spatial analyses stored on the object.
#' The function reports how many analyses are present, how many records and
#' feature columns each carries, and how those bounded counts differ. It never
#' selects a best analysis, never ranks one above another and never reports
#' which analysis is \dQuote{right}: deconvolution, colocalization and niche
#' detection are different algorithms answering different questions, so their
#' outputs are not competing answers to one question.
#'
#' @param object A \code{SingleCellExperiment} object.
#' @return A read-only comparison list with \code{status}, \code{n_analyses},
#'   \code{analysis_codes} (anonymized \code{analysis_N} codes), per-analysis
#'   bounded counts keyed by those codes, \code{claim_level}, \code{caveats},
#'   \code{recommendation = NULL} and \code{raw_values_included = FALSE}.
#' @export
compare_spatial_analyses <- function(object) {
    if (!inherits(object, "SingleCellExperiment")) {
        stop("object must be a SingleCellExperiment", call. = FALSE)
    }

    by_key <- sclet_ai_spatial_by_key(object)
    available_keys <- names(by_key)
    if (length(available_keys) < 2L) {
        return(list(
            status = "not_available",
            reason = "at_least_two_analyses_required",
            n_analyses = as.integer(length(available_keys)),
            analysis_codes = character(),
            analyses = list(),
            claim_level = "consistent_with",
            caveats = c(sclet_ai_spatial_caveats(), paste(
                "At least two spatial analyses are required for a comparison; none is selected or preferred",
                "when they are missing."
            )),
            recommendation = NULL,
            raw_values_included = FALSE
        ))
    }

    codes <- vapply(available_keys, function(key) {
        as.character(sclet_ai_spatial_code(key, available_keys))
    }, character(1L))
    analyses <- lapply(seq_along(available_keys), function(index) {
        entries <- by_key[[available_keys[[index]]]]
        counts <- lapply(entries, function(record) sclet_ai_spatial_counts(object, record))
        methods <- unique(vapply(entries, function(record) {
            as.character(sclet_ai_spatial_method(record) %||% NA_character_)
        }, character(1L)))
        methods <- methods[!is.na(methods) & nzchar(methods)]
        list(
            analysis_code = unname(codes[[index]]),
            n_records = as.integer(length(entries)),
            n_features = as.integer(max(vapply(counts, function(entry) entry$n_features, integer(1L)), 0L)),
            n_columns = as.integer(max(vapply(counts, function(entry) entry$n_columns, integer(1L)), 0L)),
            n_cells = as.integer(max(vapply(counts, function(entry) entry$n_cells, integer(1L)), 0L)),
            methods = as.character(methods),
            raw_values_included = FALSE
        )
    })
    names(analyses) <- unname(codes)

    list(
        status = "available",
        n_analyses = as.integer(length(analyses)),
        analysis_codes = unname(codes),
        analyses = analyses,
        claim_level = "consistent_with",
        caveats = c(sclet_ai_spatial_caveats(), paste(
            "This compares model OUTPUTS of different algorithms that answer different questions:",
            "deconvolution, colocalization and niche detection are not competing answers to one question."
        ), paste(
            "No analysis is preferred, ranked or selected as best by this function, and n_analyses alone says",
            "nothing about the quality of any of them."
        )),
        recommendation = NULL,
        raw_values_included = FALSE
    )
}
