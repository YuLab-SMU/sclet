# Read-only AI contract for cell-cell interaction analysis (batch B).
#
# CCI answers "which cell populations appear to talk to each other". It is the
# most leakage-prone layer in the package: the payload of an interaction table
# IS gene identity -- the ligand and receptor columns, and CellPhoneDB's
# `pathway` column, which is a "LINGAND_RECEPTOR" pair string. So this file
# reads only the columns that describe flow between groups and the strength
# columns, and never touches the identity columns at all.
#
# The interaction table lives in the ANALYSIS layer, not the state record:
# RunCellChat()/RunCellPhoneDB()/RunNicheNet() write the table through
# sclet_set_analysis(..., interactions = ...) and keep only group/interaction
# counts in the state record, so everything here reaches the table through
# get_cci().
#
# The second reason this layer needs its own contract: corroboration between two
# CCI engines is not independent evidence. CellChat and CellPhoneDB run on the
# same expression matrix and the same group labels, so agreement between them
# re-describes one measurement twice rather than measuring it twice. The
# recorded nodes therefore share an explicit dependency_group derived from the
# grouping structure, which keeps check_claim_ceiling() counting them as a
# single line of support.

#' Check whether an object is ready for cell-cell interaction analysis
#'
#' Reports whether a resolvable group column and an expression source are
#' present, and which of the three CCI engines have their packages available.
#' It deliberately recommends no method: choosing between CellChat, CellPhoneDB
#' and NicheNet is a modelling decision about the database, the statistical
#' model and the sender/receiver framing, and it belongs to the user. When more
#' than one engine is possible and none was named, the questions say the user
#' must choose and report only the count of available engines.
#'
#' @param object A \code{SingleCellExperiment} object.
#' @param method Optional engine name, one of \code{"CellChat"},
#'   \code{"CellPhoneDB"}, \code{"NicheNet"}. When supplied, that engine's
#'   package must be available.
#' @return A read-only readiness summary list.
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' readiness <- check_communication_readiness(sce)
#' readiness$status
#' readiness$questions
#' @export
check_communication_readiness <- function(object, method = NULL) {
    if (!inherits(object, "SingleCellExperiment")) {
        stop("object must be a SingleCellExperiment", call. = FALSE)
    }
    engines <- sclet_ai_clarification_engines()
    required_packages <- c(CellChat = "CellChat", CellPhoneDB = "basilisk",
        NicheNet = "nichenetr")

    columns <- colnames(SummarizedExperiment::colData(object))
    group_key <- tryCatch(ActiveIdent(object), error = function(e) NULL)
    if (!is.null(group_key) && !as.character(group_key) %in% columns &&
            !identical(as.character(group_key), "colLabels")) {
        group_key <- NULL
    }
    if (is.null(group_key)) {
        candidates <- grep("cluster|leiden|louvain|ident", columns,
            ignore.case = TRUE, value = TRUE)
        if (length(candidates)) group_key <- candidates[[1L]]
    }
    group_resolvable <- !is.null(group_key)

    assay_names <- SummarizedExperiment::assayNames(object)
    expression_available <- length(assay_names) > 0L

    packages_available <- vapply(engines, function(engine) {
        base::requireNamespace(required_packages[[engine]], quietly = TRUE)
    }, logical(1L))
    available_engines <- engines[packages_available]

    records <- tryCatch(sclet_get_state_records(object, "communication"),
        error = function(e) list())
    if (is.null(records)) records <- list()
    record_ids <- names(records)

    requested <- NULL
    if (!is.null(method)) {
        requested <- as.character(method)[[1L]]
    }
    unknown_method <- !is.null(requested) && !requested %in% engines

    questions <- character()
    if (unknown_method) {
        questions <- c(questions, paste0(
            "Unknown CCI engine '", requested, "'. Supported engines are: ",
            paste(engines, collapse = ", "), ". No engine is chosen on the user's behalf."
        ))
    }
    if (!group_resolvable) {
        questions <- c(questions, paste0(
            "Cell-cell interaction analysis groups cells into sender and receiver populations, and no group ",
            "column is resolvable; set an active identity or add a cluster column. Available colData columns: ",
            paste(columns, collapse = ", "), "."
        ))
    }
    if (!expression_available) {
        questions <- c(questions, paste0(
            "No assay is available on the object; CCI needs an expression source. ",
            "Available assays: ", paste(assay_names, collapse = ", "), "."
        ))
    }
    if (!length(available_engines)) {
        questions <- c(questions,
            "No CCI engine has its package installed (CellChat, basilisk for CellPhoneDB, nichenetr). None is recommended here.")
    } else if (is.null(requested) && length(available_engines) > 1L) {
        questions <- c(questions, paste0(
            "Interaction inference is engine-dependent: each engine uses a different ligand-receptor database, ",
            "expression test and sender/receiver framing, so results are not interchangeable. ",
            length(available_engines), " engines are available and none was named; the user must choose one. ",
            "Existing communication records: ",
            if (length(record_ids)) paste(record_ids, collapse = ", ") else "none", "."
        ))
    } else if (!is.null(requested) && requested %in% engines && !isTRUE(packages_available[[requested]])) {
        questions <- c(questions, paste0(
            "The requested engine '", requested, "' needs the '",
            required_packages[[requested]], "' package, which is not installed."
        ))
    }

    blocked <- !group_resolvable || !expression_available || !length(available_engines) ||
        unknown_method || (!is.null(requested) && !isTRUE(packages_available[[requested]]))
    list(
        status = if (blocked) "not_ready" else "ready_for_diagnostic",
        checks = list(
            group_column_resolvable = group_resolvable,
            group_column_resolved = group_key,
            expression_source_available = expression_available,
            available_assays = assay_names,
            available_methods = available_engines,
            method_packages_available = packages_available,
            communication_record_available = length(record_ids) > 0L,
            available_communication_ids = record_ids,
            raw_values_included = FALSE
        ),
        blocked_actions = if (blocked) "communication" else character(),
        questions = questions,
        notes = paste(
            "No CCI engine is recommended or ranked here: CellChat, CellPhoneDB and NicheNet differ in their",
            "ligand-receptor database, their statistical test and their sender/receiver framing, so the choice",
            "is a modelling decision reserved for the user. Inferred interactions are model estimates derived",
            "from expression plus a curated database, not observed signalling events. Two engines agreeing on",
            "the same expression matrix and the same group labels are not independent measurements, and the",
            "evidence layer enforces that by sharing a dependency group across communication runs."
        )
    )
}

# Resolve the group column the same way the CCI runners do, so anonymized codes
# line up with what the engine actually grouped on.
sclet_ai_communication_group_key <- function(object, record) {
    recorded <- as.character(record$inputs$group %||% "")
    columns <- colnames(SummarizedExperiment::colData(object))
    if (length(recorded) == 1L && nzchar(recorded) && recorded %in% columns) {
        return(recorded)
    }
    if (length(recorded) == 1L && identical(recorded, "sclet_active_ident")) {
        return(recorded)
    }
    readiness <- check_communication_readiness(object)
    readiness$checks$group_column_resolved
}

# Stable across engines on the same object: the multiset of group sizes. It is
# what makes two CCI runs share a dependency_group without hashing user-supplied
# column names, labels or any identity-bearing string.
sclet_ai_communication_dependency_group <- function(object, group_key) {
    counts <- NULL
    if (!is.null(group_key)) {
        values <- tryCatch(
            if (identical(group_key, "colLabels") || identical(group_key, "sclet_active_ident"))
                as.character(SingleCellExperiment::colLabels(object))
            else as.character(SummarizedExperiment::colData(object)[[group_key]]),
            error = function(e) NULL
        )
        if (!is.null(values)) counts <- as.integer(sort(table(values)))
    }
    if (is.null(counts)) counts <- as.integer(ncol(object))
    payload <- paste0(counts, collapse = ",")
    if (base::requireNamespace("digest", quietly = TRUE)) {
        return(paste0("cci_", substring(digest::digest(payload, algo = "sha256"), 1L, 32L)))
    }
    paste0("cci_", substring(payload, 1L, 32L))
}

# Map population labels onto anonymized codes. Falls back to positional indexes
# when a label is not one of the object's cluster levels, so a code is always
# well formed and never leaks the original string.
sclet_ai_communication_codes <- function(labels, levels) {
    labels <- as.character(labels)
    labels <- labels[!is.na(labels)]
    if (!length(labels)) return(labels[0L])
    index <- match(labels, levels)
    if (any(is.na(index))) {
        positional <- seq_along(unique(labels))
        fallback <- match(labels, unique(labels))
        index[is.na(index)] <- positional[fallback[is.na(index)]]
    }
    paste0("cluster_", index)
}

# The only columns this layer is allowed to read from an interaction table.
# ligand, receptor and pathway are never selected: they are the identity payload.
sclet_ai_communication_table_facts <- function(table) {
    if (is.null(table) || !length(dim(table)) || !nrow(table)) {
        return(list(
            n_interactions = 0L, source_group_codes = character(),
            target_group_codes = character(), pair_flow_counts = numeric(),
            probability_quantile_25 = NA_real_, probability_median = NA_real_,
            probability_quantile_75 = NA_real_, probability_max = NA_real_,
            p_value_quantile_25 = NA_real_, p_value_median = NA_real_,
            p_value_quantile_75 = NA_real_, n_significant_pairs = 0L
        ))
    }
    columns <- tolower(colnames(table))
    source_column <- colnames(table)[columns %in% "source_group"][1L]
    target_column <- colnames(table)[columns %in% "target_group"][1L]
    probability_column <- colnames(table)[columns %in% "probability"][1L]
    p_column <- colnames(table)[columns %in% "p_value"][1L]
    if (is.na(source_column) || is.na(target_column)) {
        return(NULL)
    }
    source_levels <- unique(as.character(table[[source_column]]))
    target_levels <- unique(as.character(table[[target_column]]))
    source_codes <- sclet_ai_communication_codes(table[[source_column]], source_levels)
    target_codes <- sclet_ai_communication_codes(table[[target_column]], target_levels)
    unique_source <- source_codes[!duplicated(source_codes)]
    unique_target <- target_codes[!duplicated(target_codes)]
    # Full (source, target) grid so the count vector has a position for every
    # possible edge, including empty ones; a reader realigns it from the two
    # code vectors. The join keys are built from anonymized codes, never from
    # the population labels themselves.
    pairs <- paste(source_codes, target_codes, sep = "|")
    grid_keys <- as.vector(outer(unique_source, unique_target,
        function(a, b) paste(a, b, sep = "|")))
    flow <- as.numeric(tabulate(match(pairs, grid_keys), length(grid_keys)))

    probability <- if (!is.na(probability_column))
        suppressWarnings(as.numeric(table[[probability_column]])) else numeric()
    probability <- probability[is.finite(probability)]
    p_value <- if (!is.na(p_column))
        suppressWarnings(as.numeric(table[[p_column]])) else numeric()
    p_value <- p_value[is.finite(p_value)]
    probability_quantiles <- if (length(probability))
        round(unname(stats::quantile(probability, c(0.25, 0.5, 0.75))), 6L) else rep(NA_real_, 3L)
    p_quantiles <- if (length(p_value))
        round(unname(stats::quantile(p_value, c(0.25, 0.5, 0.75))), 8L) else rep(NA_real_, 3L)
    list(
        n_interactions = as.integer(nrow(table)),
        source_group_codes = unique_source,
        target_group_codes = unique_target,
        # Row-major over (source, target) so a reader can rebuild the alignment
        # from the two code vectors; unnamed, because names would be labels.
        pair_flow_counts_rowmajor = flow,
        probability_quantile_25 = as.numeric(probability_quantiles[[1L]]),
        probability_median = as.numeric(probability_quantiles[[2L]]),
        probability_quantile_75 = as.numeric(probability_quantiles[[3L]]),
        probability_max = if (length(probability)) round(max(probability), 6L) else NA_real_,
        p_value_quantile_25 = as.numeric(p_quantiles[[1L]]),
        p_value_median = as.numeric(p_quantiles[[2L]]),
        p_value_quantile_75 = as.numeric(p_quantiles[[3L]]),
        # `n_significant_pairs` counts p < 0.05 rows IN THIS TABLE. When the
        # engine already filtered the table on its own test (CellChat does:
        # measured 733 of 733 rows with max p = 0.04 on PBMC 10k, 2026-10-20),
        # that count restates the row count. `p_value_max` is what lets a
        # reader tell the two situations apart without trusting a note.
        p_value_max = if (length(p_value)) round(max(p_value), 8L) else NA_real_,
        n_significant_pairs = as.integer(sum(p_value < 0.05))
    )
}

# The engine is reported as a code, not a string: a character value in evidence
# must match the anonymized-code pattern, and "CellChat" does not.
sclet_ai_communication_method_code <- function(method) {
    engines <- sclet_ai_clarification_engines()
    index <- match(as.character(method %||% "")[[1L]], engines)
    if (is.na(index)) 0L else as.integer(index)
}

#' Record a bounded aggregate evidence node for a CCI analysis
#'
#' Reads the interaction table through \code{get_cci()} and registers one
#' evidence node describing group-to-group flow and strength aggregates. The
#' ligand, receptor and pathway columns are never read, because in this layer
#' those columns are gene identity and an evidence node is eligible to leave the
#' process as provider context. Every communication node on an object is given
#' the same dependency group, so agreement between two engines over the same
#' expression matrix and the same labels is not counted as two independent lines
#' of support by \code{check_claim_ceiling()}.
#'
#' @param object A \code{SingleCellExperiment} object.
#' @param name Communication record id.
#' @param method Optional engine name, used only to derive an integer code.
#' @return The object, with one \code{ev:communication_<name>} node appended
#'   when a communication record exists; otherwise unchanged.
#' @keywords internal
#' @noRd
sclet_ai_record_communication_evidence <- function(object, name, method = NULL) {
    record <- tryCatch(sclet_get_state_record(object, "communication", name),
        error = function(e) NULL)
    if (is.null(record)) return(object)

    table <- tryCatch(get_cci(object, id = name), error = function(e) NULL)
    facts <- sclet_ai_communication_table_facts(table)
    if (is.null(facts)) return(object)

    group_key <- sclet_ai_communication_group_key(object, record)
    dependency_group <- sclet_ai_communication_dependency_group(object, group_key)
    resolved_method <- method %||% record$method

    values <- list(
        n_interactions = facts$n_interactions,
        n_source_groups = as.integer(length(facts$source_group_codes)),
        n_target_groups = as.integer(length(facts$target_group_codes)),
        source_group_codes = facts$source_group_codes,
        target_group_codes = facts$target_group_codes,
        pair_flow_counts_rowmajor = facts$pair_flow_counts_rowmajor,
        probability_quantile_25 = facts$probability_quantile_25,
        probability_median = facts$probability_median,
        probability_quantile_75 = facts$probability_quantile_75,
        probability_max = facts$probability_max,
        p_value_quantile_25 = facts$p_value_quantile_25,
        p_value_median = facts$p_value_median,
        p_value_quantile_75 = facts$p_value_quantile_75,
        p_value_max = facts$p_value_max,
        n_significant_pairs = facts$n_significant_pairs,
        method_code = sclet_ai_communication_method_code(resolved_method),
        communication_is_model_estimate = TRUE,
        # Records the honest fact of who chose the engine: the AI path passes a
        # method explicitly, while a record read back without one was resolved
        # by the runner's own default. It is not a claim that the user was asked.
        engine_was_named_by_caller = !is.null(method),
        engines_on_shared_inputs_are_not_independent = TRUE,
        raw_values_included = FALSE
    )
    evidence <- list(
        id = paste0("ev:communication_", name),
        kind = "deterministic_summary",
        values = values,
        claim_level = "consistent_with",
        dependency_group = dependency_group
    )
    sclet_ai_record_evidence_for_run(object, evidence, source = name,
        parents = character(), scope = NULL)
}

#' Summarize CCI results in bounded form without claiming a conclusion
#'
#' Read-only companion to \code{RunCCI()}: it reports the same group-flow
#' aggregates the evidence node carries and writes nothing. It never ranks
#' populations or pathways, never names a ligand or receptor, and never
#' recommends an engine.
#'
#' @param object A \code{SingleCellExperiment} object.
#' @param id Optional communication record id. When \code{NULL}, the active one.
#' @return A read-only summary list.
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' # No communication run is on record, so nothing is reported.
#' summary <- summarize_communication_evidence(sce)
#' summary$status
#' summary$reason
#' @export
summarize_communication_evidence <- function(object, id = NULL) {
    if (!inherits(object, "SingleCellExperiment")) {
        stop("object must be a SingleCellExperiment", call. = FALSE)
    }
    unavailable <- function(reason) {
        list(
            status = "not_available", reason = reason, n_records = 0L,
            recommendation = NULL, raw_values_included = FALSE,
            caveats = "No interaction is reported, no population is ranked and no ligand or receptor is named."
        )
    }
    resolved_id <- if (!is.null(id) && nzchar(id)) as.character(id)[[1L]] else {
        tryCatch(sclet_get_active_state(object, "communication"), error = function(e) NULL)
    }
    if (is.null(resolved_id)) return(unavailable("no_communication_record"))
    record <- tryCatch(sclet_get_state_record(object, "communication", resolved_id),
        error = function(e) NULL)
    if (is.null(record)) return(unavailable("no_communication_record"))

    table <- tryCatch(get_cci(object, id = resolved_id), error = function(e) NULL)
    facts <- sclet_ai_communication_table_facts(table)
    if (is.null(facts)) return(unavailable("no_interaction_table"))

    list(
        status = "available",
        communication_id = resolved_id,
        method_code = sclet_ai_communication_method_code(record$method),
        group_column_present = !is.null(sclet_ai_communication_group_key(object, record)),
        n_interactions = facts$n_interactions,
        n_source_groups = as.integer(length(facts$source_group_codes)),
        n_target_groups = as.integer(length(facts$target_group_codes)),
        source_group_codes = facts$source_group_codes,
        target_group_codes = facts$target_group_codes,
        pair_flow_counts_rowmajor = facts$pair_flow_counts_rowmajor,
        probability_median = facts$probability_median,
        probability_max = facts$probability_max,
        p_value_median = facts$p_value_median,
        n_significant_pairs = facts$n_significant_pairs,
        evidence_node_recorded = paste0("ev:communication_", resolved_id) %in%
            names(tryCatch(sclet_ai_evidence_get_all(object), error = function(e) list())),
        claim_level = "consistent_with",
        recommendation = NULL,
        raw_values_included = FALSE,
        caveats = paste(
            "Interactions are model-inferred ligand-receptor events computed from expression levels and a",
            "curated database, not observed signalling. A predicted pair is not evidence that the cells",
            "signalled, and the absence of a pair is uninformative because each engine filters against its",
            "own database and complex model. Only group-flow counts and strength quantiles are reported: no",
            "ligand, receptor or pathway identity is included, and no per-cell value is read. Engines run on",
            "the same expression matrix and labels share a dependency group, so their agreement is not two",
            "independent lines of support."
        )
    )
}
