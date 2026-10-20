# Read-only AI contract for enrichment / gene-program interpretation (batch A).
#
# Enrichment is the interpretation layer over differential expression: without
# it an AI can report that a cluster differs, but never what that cluster IS.
# This file exposes that layer to the AI under two non-negotiable properties.
#
# 1. NO GENE OR TERM IDENTITY ENTERS EVIDENCE. `RunEnrichment()` stores the
#    full compareCluster result in artifacts$result, whose columns include
#    `geneID` (slash-joined gene symbols), `Description` (GO/KEGG term names)
#    and `ID` (term accessions). Evidence carries only counts, anonymized
#    cluster_N codes and adjusted-p quantiles. The complete table stays
#    reachable through get_enrichment() for R-side work; it is never copied
#    into an evidence node, because evidence nodes feed provider payloads.
#
# 2. ENRICHMENT DECLARES ITS DE LINEAGE. Enrichment is computed FROM the DE
#    marker table, so it is not a second independent line of support. The
#    recorded node sets `parents` to the DE evidence node(s) for that detest
#    source, which sclet_ai_claim_dependency_reasons() resolves through the
#    ancestor chain, so check_claim_ceiling() keeps counting DE + enrichment
#    as ONE line. Without this, merely re-describing one measurement would
#    look like corroboration and would inflate every downstream claim.
#
# Neither function here ranks, suggests or recommends anything: db, ont,
# keyType and the significance cutoffs are biological modelling choices that
# belong to the user.

#' Check whether an object is ready for enrichment analysis
#'
#' Reports whether a completed differential-expression test and a usable
#' cluster assignment are present for \code{RunEnrichment()}. This diagnostic
#' deliberately proposes no database, ontology, key type or significance
#' cutoff: those are modelling choices only the user can make. It also guards a
#' specific silent failure -- \code{RunEnrichment()} falls back to labelling
#' every marker \code{"Cluster1"} when the marker table carries no cluster
#' column, so it warns instead of erroring and produces a single-group
#' comparison that means nothing. That hazard is reported rather than absorbed.
#'
#' @param object A \code{SingleCellExperiment} object.
#' @param detest_id Optional detest record id. When supplied, that exact
#'   record must exist.
#' @return A read-only readiness summary.
#' @export
check_enrichment_readiness <- function(object, detest_id = NULL) {
    if (!inherits(object, "SingleCellExperiment")) {
        stop("object must be a SingleCellExperiment", call. = FALSE)
    }
    detest <- tryCatch(get_detest(object, detest_id), error = function(e) NULL)
    detest_available <- !is.null(detest)
    resolved_detest_id <- if (detest_available) as.character(detest$id %||% "") else NULL

    # The marker table is what RunEnrichment() actually groups on, so the
    # degeneracy check reads the table rather than the object's colData.
    markers <- if (detest_available) detest$artifacts$result else NULL
    marker_columns <- if (is.data.frame(markers) || is.matrix(markers)) {
        colnames(as.data.frame(markers))
    } else {
        character()
    }
    marker_has_cluster <- "cluster" %in% tolower(marker_columns)

    columns <- colnames(SummarizedExperiment::colData(object))
    cluster_key <- tryCatch(ActiveIdent(object), error = function(e) NULL)
    if (!is.null(cluster_key) && !as.character(cluster_key) %in% columns &&
            !identical(as.character(cluster_key), "colLabels")) {
        cluster_key <- NULL
    }
    if (is.null(cluster_key)) {
        candidates <- grep("cluster|leiden|louvain|ident", columns,
            ignore.case = TRUE, value = TRUE)
        if (length(candidates)) cluster_key <- candidates[[1L]]
    }
    cluster_resolvable <- !is.null(cluster_key)

    has_cluster_profiler <- base::requireNamespace("clusterProfiler", quietly = TRUE)
    default_org_db <- "org.Hs.eg.db"
    has_org_db <- base::requireNamespace(default_org_db, quietly = TRUE)
    packages_available <- has_cluster_profiler && has_org_db

    questions <- character()
    if (!detest_available) {
        questions <- c(questions, paste0(
            "Enrichment is computed from a differential-expression marker table, and no completed ",
            "detest record is available. Run RunDEtest() first",
            if (!is.null(detest_id)) {
                paste0(" (requested detest id '", as.character(detest_id), "' is not a completed record)")
            } else {
                ""
            }, "."
        ))
    } else if (!marker_has_cluster) {
        questions <- c(questions, paste0(
            "The detest marker table carries no cluster column, so RunEnrichment() would label every ",
            "marker as a single group and compare one pooled set against itself: the result would be a ",
            "warning, not an error, and would carry no biological meaning. Re-run RunDEtest() so the ",
            "marker table retains the cluster each marker was found in.",
            " Available marker columns: ", paste(marker_columns, collapse = ", "), "."
        ))
    }
    if (!cluster_resolvable) {
        questions <- c(questions, paste0(
            "No cluster assignment is available on the object; enrichment groups markers by cluster, so ",
            "an explicit cluster column is required. Available colData columns: ",
            paste(columns, collapse = ", "), "."
        ))
    }
    if (!packages_available) {
        missing_packages <- c(
            if (!has_cluster_profiler) "clusterProfiler" else NULL,
            if (!has_org_db) default_org_db else NULL
        )
        questions <- c(questions, paste0(
            "Enrichment needs ", paste(missing_packages, collapse = " and "),
            " installed; install it with BiocManager::install().",
            " No database, ontology or cutoff is suggested here."
        ))
    }

    # The degeneracy gate blocks, because a run that silently returns a
    # meaningless comparison is worse for a planning AI than a blocked one.
    ready <- detest_available && marker_has_cluster && cluster_resolvable &&
        packages_available
    list(
        status = if (ready) "ready_for_diagnostic" else "not_ready",
        checks = list(
            detest_state_available = detest_available,
            detest_id_resolved = resolved_detest_id,
            detest_marker_column_has_cluster = marker_has_cluster,
            cluster_identity_resolvable = cluster_resolvable,
            cluster_key_resolved = cluster_key,
            enrichment_packages_available = packages_available,
            available_dbs = c("GO", "KEGG"),
            raw_values_included = FALSE
        ),
        blocked_actions = if (ready) character() else "enrichment",
        questions = questions,
        notes = paste(
            "No database, ontology, key type or significance cutoff is recommended here: each is a",
            "modelling choice reserved for the user. Enrichment categories are derived from the",
            "differential-expression marker table and are not independent measurements of anything; they",
            "inherit that table's thresholds and confounders, and enrichment by itself does not",
            "constitute a second line of independent evidence for a claim ceiling."
        )
    )
}

# Resolve the enrichment record the AI layer should read, or NULL.
sclet_ai_enrichment_record <- function(object, name = NULL) {
    if (!is.null(name) && nzchar(name)) {
        return(tryCatch(sclet_get_state_record(object, "enrichment", name),
            error = function(e) NULL))
    }
    tryCatch(get_enrichment(object), error = function(e) NULL)
}

# The detest record id an enrichment run was derived from. RunEnrichment()
# stores it in inputs$detest_id, so lineage is read from the record itself
# rather than guessed.
sclet_ai_enrichment_parent_source <- function(record) {
    as.character(record$inputs$detest_id %||% "")
}

# Evidence ids already recorded against a given analysis source.
sclet_ai_enrichment_parent_nodes <- function(object, source_id) {
    if (!length(source_id) || !nzchar(source_id)) return(character())
    nodes <- tryCatch(sclet_ai_evidence_get_all(object), error = function(e) list())
    if (!length(nodes)) return(character())
    find_source <- function(id) {
        names(nodes)[vapply(nodes, function(node) {
            identical(as.character(node$source %||% ""), id)
        }, logical(1L))]
    }
    hits <- find_source(source_id)
    if (length(hits)) return(as.character(hits))
    # A marker table registered as a DERIVED record -- the id-converted copy
    # KEGG needs, for instance -- carries no evidence node of its own, so the
    # lookup above returns nothing and the enrichment node ends up with no
    # parents. That is not cosmetic: the claim ceiling counts nodes without a
    # shared parent as INDEPENDENT lines of support, so two annotations of one
    # marker table read as corroboration. Measured on the M1 report
    # (2026-10-20): markers + GO + KEGG came out as two independent lines when
    # they are one. Follow the recorded `parent_detest` one hop.
    derived <- tryCatch(
        sclet_get_state_record(object, "detest", source_id),
        error = function(e) NULL
    )
    parent_detest <- as.character(derived$inputs$parent_detest %||% "")
    if (!nzchar(parent_detest) || identical(parent_detest, source_id)) {
        return(character())
    }
    as.character(find_source(parent_detest))
}

# Map cluster labels onto anonymized cluster_N codes. Labels come from the
# marker table and are user- or pipeline-chosen, so a label may already be a
# code, may be an arbitrary string, or may not appear among the object's
# cluster levels. Order is preserved so a companion numeric vector stays
# aligned; unresolvable labels fall back to a positional index.
sclet_ai_enrichment_cluster_codes <- function(labels, cluster_levels) {
    labels <- as.character(labels)
    labels <- labels[!is.na(labels) & nzchar(labels)]
    if (!length(labels)) return(character())
    index <- match(labels, cluster_levels)
    if (any(is.na(index))) {
        pattern_ok <- grepl("^[0-9]+$", sub("^(group|cluster|condition|sample|route)_", "", labels))
        position <- seq_along(labels)
        index[is.na(index) & pattern_ok] <- as.integer(
            sub("^(group|cluster|condition|sample|route)_", "", labels[is.na(index) & pattern_ok])
        )
        index[is.na(index)] <- position[is.na(index)]
    }
    paste0("cluster_", index)
}

# Cluster label universe used for anonymization: the object's own cluster
# levels when present, otherwise the marker table's labels.
sclet_ai_enrichment_cluster_levels <- function(object, markers) {
    columns <- colnames(SummarizedExperiment::colData(object))
    cluster_key <- tryCatch(ActiveIdent(object), error = function(e) NULL)
    if (!is.null(cluster_key) && as.character(cluster_key) %in% columns) {
        levels <- levels(as.factor(as.character(SummarizedExperiment::colData(object)[[cluster_key]])))
        if (length(levels)) return(levels)
    }
    candidates <- grep("cluster|leiden|louvain|ident", columns, ignore.case = TRUE, value = TRUE)
    for (candidate in candidates) {
        levels <- levels(as.factor(as.character(SummarizedExperiment::colData(object)[[candidate]])))
        if (length(levels)) return(levels)
    }
    if (is.data.frame(markers) || is.matrix(markers)) {
        marker_frame <- as.data.frame(markers)
        column <- tolower(colnames(marker_frame))
        hit <- colnames(marker_frame)[column %in% "cluster"][1L]
        if (!is.na(hit)) return(as.character(unique(marker_frame[[hit]])))
    }
    character()
}

# Read the enrichment result table's bounded facts. Returns NULL when the
# record carries nothing usable. geneID / Description / ID are never read out.
sclet_ai_enrichment_bounded_facts <- function(record, cluster_levels) {
    result <- record$artifacts$result
    if (is.null(result)) return(NULL)
    table <- tryCatch(as.data.frame(result), error = function(e) NULL)
    if (is.null(table) || !nrow(table)) return(NULL)
    columns <- tolower(colnames(table))
    cluster_column <- colnames(table)[columns %in% "cluster"][1L]
    padj_column <- colnames(table)[columns %in% "p.adjust"][1L]
    if (is.na(cluster_column) || is.na(padj_column)) return(NULL)

    labels <- as.character(table[[cluster_column]])
    codes <- sclet_ai_enrichment_cluster_codes(labels, cluster_levels)
    unique_codes <- codes[!duplicated(codes)]
    per_cluster <- as.numeric(vapply(unique_codes, function(code) sum(codes == code), numeric(1L)))
    padj <- suppressWarnings(as.numeric(table[[padj_column]]))
    padj <- padj[is.finite(padj)]
    quantiles <- if (length(padj)) {
        round(unname(stats::quantile(padj, c(0.25, 0.5, 0.75))), 8L)
    } else {
        NULL
    }
    list(
        n_enriched = as.integer(nrow(table)),
        n_clusters_enriched = as.integer(length(unique_codes)),
        cluster_codes = unique_codes,
        enriched_per_cluster = per_cluster,
        padj_quantile_25 = if (length(quantiles)) as.numeric(quantiles[[1L]]) else NA_real_,
        padj_median = if (length(quantiles)) as.numeric(quantiles[[2L]]) else NA_real_,
        padj_quantile_75 = if (length(quantiles)) as.numeric(quantiles[[3L]]) else NA_real_
    )
}

#' Record a bounded aggregate evidence node for an enrichment analysis
#'
#' Reads back the enrichment state record written by \code{RunEnrichment()} and
#' registers one evidence node describing it. Only aggregate counts, anonymized
#' cluster codes and adjusted-p-value quantiles are recorded: the marker genes
#' and the enriched category names stay in the state record, because an
#' evidence node is a claim-support artefact and is eligible to leave the
#' process as provider context. The node declares the detest evidence it was
#' derived from as its parent, so a claim ceiling still counts the DE test and
#' its enrichment as a single line of support.
#'
#' @param object A \code{SingleCellExperiment} object.
#' @param name Enrichment record id.
#' @param detest_id Optional detest record id used to resolve the lineage.
#'   When NULL it is read from the enrichment record itself.
#' @return The object, with one \code{ev:enrichment_<name>} node appended when
#'   an enrichment record exists; otherwise unchanged.
#' @keywords internal
#' @noRd
sclet_ai_record_enrichment_evidence <- function(object, name, detest_id = NULL) {
    record <- sclet_ai_enrichment_record(object, name)
    if (is.null(record)) return(object)

    marker_table <- if (!is.null(detest_id)) {
        tryCatch(get_detest(object, detest_id), error = function(e) NULL)
    } else {
        tryCatch(get_detest(object, NULL), error = function(e) NULL)
    }
    markers <- if (!is.null(marker_table)) marker_table$artifacts$result else NULL
    cluster_levels <- sclet_ai_enrichment_cluster_levels(object, markers)
    facts <- sclet_ai_enrichment_bounded_facts(record, cluster_levels)
    if (is.null(facts)) return(object)

    # Lineage: prefer the id RunEnrichment() recorded, fall back to the id the
    # caller supplied, so the DE node can be located under either.
    source_detest_id <- sclet_ai_enrichment_parent_source(record)
    if (!nzchar(source_detest_id)) {
        source_detest_id <- as.character(detest_id %||% "")
    }
    parents <- sclet_ai_enrichment_parent_nodes(object, source_detest_id)

    values <- list(
        n_enriched = facts$n_enriched,
        n_clusters_enriched = facts$n_clusters_enriched,
        cluster_codes = facts$cluster_codes,
        # Unnamed and aligned with cluster_codes by position: attaching label
        # names here would reintroduce the identity the codes exist to hide.
        enriched_per_cluster = facts$enriched_per_cluster,
        padj_quantile_25 = facts$padj_quantile_25,
        padj_median = facts$padj_median,
        padj_quantile_75 = facts$padj_quantile_75,
        enrichment_is_derived_from_markers = TRUE,
        db_is_user_selected = TRUE,
        raw_values_included = FALSE
    )
    evidence <- list(
        id = paste0("ev:enrichment_", name),
        kind = "deterministic_summary",
        values = values,
        claim_level = "consistent_with"
    )
    sclet_ai_record_evidence_for_run(object, evidence,
        source = as.character(record$id %||% name),
        parents = parents, scope = NULL)
}

#' Summarize recorded enrichment evidence without claiming a result
#'
#' Read-only companion to \code{RunEnrichment()}: it reports what an
#' enrichment run recorded, in the same bounded shape the evidence node uses,
#' and writes nothing. It never ranks categories, never recommends a database
#' or an ontology, and never surfaces gene or term identity.
#'
#' @param object A \code{SingleCellExperiment} object.
#' @param name Optional enrichment record id. When NULL, the active record.
#' @return A read-only summary list.
#' @export
summarize_enrichment_evidence <- function(object, name = NULL) {
    unavailable <- function(reason) {
        list(
            status = "not_available", reason = reason, n_records = 0L,
            claim_level = "consistent_with", recommendation = NULL,
            raw_values_included = FALSE,
            caveats = "No enrichment evidence is reported and no category, gene or term is named."
        )
    }
    if (!inherits(object, "SingleCellExperiment")) {
        stop("object must be a SingleCellExperiment", call. = FALSE)
    }
    record <- sclet_ai_enrichment_record(object, name)
    if (is.null(record)) return(unavailable("no_enrichment_record"))

    markers <- tryCatch({
        detest <- get_detest(object, sclet_ai_enrichment_parent_source(record))
        detest$artifacts$result
    }, error = function(e) NULL)
    if (is.null(markers)) {
        markers <- tryCatch(get_detest(object, NULL)$artifacts$result, error = function(e) NULL)
    }
    cluster_levels <- sclet_ai_enrichment_cluster_levels(object, markers)
    facts <- sclet_ai_enrichment_bounded_facts(record, cluster_levels)
    if (is.null(facts)) return(unavailable("no_enrichment_result"))

    source_detest_id <- sclet_ai_enrichment_parent_source(record)
    parents <- sclet_ai_enrichment_parent_nodes(object, source_detest_id)

    list(
        status = "available",
        enrichment_id = as.character(record$id %||% name %||% ""),
        detest_id = if (nzchar(source_detest_id)) source_detest_id else NA_character_,
        method = as.character(record$method %||% ""),
        n_enriched = facts$n_enriched,
        n_clusters_enriched = facts$n_clusters_enriched,
        cluster_codes = facts$cluster_codes,
        enriched_per_cluster = facts$enriched_per_cluster,
        padj_quantile_25 = facts$padj_quantile_25,
        padj_median = facts$padj_median,
        padj_quantile_75 = facts$padj_quantile_75,
        derived_from_evidence_node = paste0("ev:enrichment_", record$id %||% name %||% "") %in%
            names(tryCatch(sclet_ai_evidence_get_all(object), error = function(e) list())),
        lineage_parents = parents,
        claim_level = "consistent_with",
        recommendation = NULL,
        raw_values_included = FALSE,
        caveats = paste(
            "Enrichment categories are derived from the differential-expression marker table and are",
            "model output, not independent measurements; they inherit that table's thresholds and",
            "confounders. Because enrichment re-describes a single DE measurement it does not add a",
            "second independent line of support to a claim ceiling. No gene symbol, category name or",
            "category accession is included here; the full result stays in the state record for R-side",
            "access through get_enrichment()."
        )
    )
}

# The marker-table shape gate, shared by the action's prerequisite and its
# handler.
#
# A plan that chains run_de_test -> run_enrichment cannot be checked at
# validation time: validation runs before anything executes, so the marker
# table does not exist yet and the promise of it is all that can be known. This
# is what makes the execution-time re-check load-bearing rather than
# redundant -- a table without a cluster column makes RunEnrichment() pool
# every marker into one group, warn, and return a meaningless comparison that
# still looks like a successful run.
sclet_ai_enrichment_marker_problem <- function(object, detest_id = NULL) {
    detest <- tryCatch(get_detest(object, detest_id), error = function(e) NULL)
    if (is.null(detest)) {
        return(NULL)
    }
    markers <- detest$artifacts$result
    columns <- if (is.data.frame(markers) || is.matrix(markers)) {
        colnames(as.data.frame(markers))
    } else {
        character()
    }
    if (!"cluster" %in% tolower(columns)) {
        return(paste0(
            "cluster_identity_missing: run_enrichment groups markers by cluster, but the ",
            "detest marker table carries no cluster column, so RunEnrichment() would label ",
            "every marker as one pooled group and compare it against itself, warning rather ",
            "than erroring. Re-run RunDEtest() so the table retains each marker's cluster. ",
            "Available marker columns: ", paste(columns, collapse = ", "), "."
        ))
    }
    NULL
}
