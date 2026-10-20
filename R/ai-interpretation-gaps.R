# Batch C: closing interpretation gaps for three single-cell stages that could
# previously run without leaving anything the AI layer can reason over.
#
#  - run_doublet_detection called scDblFinder and wrote per-cell colData, but
#    registered no evidence node at all, so a planning agent could execute it
#    and then find nothing in the ledger to cite.
#  - RunDEtest had no readiness check, so a single-group object (or one whose
#    smallest group has one cell) reached the DE engine and produced an empty
#    marker table, which then surfaced far away as an enrichment no-result.
#  - RunMilo takes sample and condition columns, i.e. design semantics, and the
#    existing confirmation gate is keyed on the batch role only.
#
# Shared discipline: aggregates only. No per-cell vector, no cell barcode, no
# class label, no sample or condition label enters evidence, because an evidence
# node is eligible to leave the process as provider context.

# ---------------------------------------------------------------------------
# C1: doublet evidence
# ---------------------------------------------------------------------------

# Bounded aggregates from scDblFinder's two colData columns. Returns NULL when
# the columns are absent, so a recorder never invents a verdict.
sclet_ai_doublet_facts <- function(object) {
    columns <- colnames(SummarizedExperiment::colData(object))
    score_column <- "scDblFinder.score"
    class_column <- "scDblFinder.class"
    if (!all(c(score_column, class_column) %in% columns)) return(NULL)
    coldata <- SummarizedExperiment::colData(object)
    calls <- as.character(coldata[[class_column]])
    calls <- calls[!is.na(calls)]
    if (!length(calls)) return(NULL)
    scores <- suppressWarnings(as.numeric(coldata[[score_column]]))
    scores <- scores[is.finite(scores)]
    quantiles <- if (length(scores)) {
        round(unname(stats::quantile(scores, c(0.25, 0.5, 0.75))), 6L)
    } else {
        rep(NA_real_, 3L)
    }
    n_cells <- length(calls)
    n_doublets <- sum(calls == "doublet")
    list(
        n_cells = as.integer(n_cells),
        n_called_doublets = as.integer(n_doublets),
        doublet_fraction = round(n_doublets / n_cells, 6L),
        score_quantile_25 = as.numeric(quantiles[[1L]]),
        score_median = as.numeric(quantiles[[2L]]),
        score_quantile_75 = as.numeric(quantiles[[3L]]),
        score_max = if (length(scores)) round(max(scores), 6L) else NA_real_
    )
}

#' Record a bounded evidence node for a doublet call
#'
#' Registers what scDblFinder concluded as aggregate counts and score quantiles
#' only. The per-cell class labels and per-cell scores stay in colData: a score
#' per cell is cell-level identity plus a model output, and an evidence node can
#' be sent to a provider as context. The claim level is deliberately
#' \code{"associated"} rather than \code{"observed"} or
#' \code{"consistent_with"}: one classifier's call, with no independent signal
#' corroborating it, is not a measurement of doubleness.
#'
#' @param object A \code{SingleCellExperiment} object.
#' @param name Preprocess record id carrying the doublet call.
#' @return The object with one \code{ev:doublet_<name>} node appended, or
#'   unchanged when no doublet columns exist or the node already exists.
#' @keywords internal
#' @noRd
sclet_ai_record_doublet_evidence <- function(object, name = "doublet_finder",
                                             replace = FALSE) {
    facts <- sclet_ai_doublet_facts(object)
    if (is.null(facts)) return(object)
    evidence <- list(
        id = paste0("ev:doublet_", name),
        kind = "deterministic_summary",
        values = c(facts, list(
            doublet_call_is_model_estimate = TRUE,
            no_cell_removed_or_filtered = TRUE,
            raw_values_included = FALSE
        )),
        claim_level = "associated"
    )
    # Routed through the re-run-aware recorder: a second doublet call on the
    # same object reproduces the node id, and swallowing the duplicate-id error
    # left the FIRST run's counts in the ledger while colData showed the second
    # (measured 2026-10-21: colData 0 doublets, node still 150).
    sclet_ai_record_evidence_for_run(object, evidence, source = name,
        parents = character(), scope = NULL, replace = replace)
}

# Record the bounded aggregate evidence node for a completed CNV run.
#
# Extracted from the action handler so the node can be RE-RECORDED from the
# record without repeating the inferCNV run -- which matters because that run
# costs about an hour (measured 2026-10-20: 61.8 minutes on 1,400 cells with
# HMM on). The operating point belongs in the node: a user report on GSE150290
# ran the same line with hmm = FALSE after the HMM stage OOMed at 32 GB, and
# the old node could not tell that run apart from an HMM run even though the
# burden comes from a different profile.
sclet_ai_record_cnv_evidence <- function(object, name = "cnv", replace = FALSE) {
    records <- sclet_get_analysis(object, "malignant_cnv")
    record <- records[[name]]
    if (is.null(record) || !identical(record$status, "completed")) {
        # A recorded `unavailable` state is an honest answer, not a failure; no
        # evidence node is written for it.
        return(object)
    }
    burden <- record$burden$burden
    values <- list(
        n_cells = as.integer(record$n_cells_submitted),
        n_cells_retained = as.integer(record$n_cells_retained),
        n_cells_dropped = as.integer(length(record$dropped_cells)),
        n_reference_cells = as.integer(record$n_reference_cells),
        n_genes = as.integer(record$n_genes),
        burden_q10 = as.numeric(stats::quantile(burden, 0.1, na.rm = TRUE, names = FALSE)),
        burden_q50 = as.numeric(stats::quantile(burden, 0.5, na.rm = TRUE, names = FALSE)),
        burden_q90 = as.numeric(stats::quantile(burden, 0.9, na.rm = TRUE, names = FALSE)),
        # `parameter_` keeps these inside the privacy allowlist.
        parameter_hmm = isTRUE(record$params$hmm),
        parameter_cutoff = as.numeric(record$params$cutoff),
        parameter_n_genes = as.integer(record$n_genes)
    )
    sclet_ai_record_evidence_for_run(object, list(
            id = paste0("ev:cnv_", sclet_ai_evidence_id_token(name)),
            kind = "deterministic_summary",
            values = values,
            claim_level = "associated"
        ),
        source = name,
        parents = character(),
        scope = NULL,
        replace = replace)
}

# ---------------------------------------------------------------------------
# C2: differential-expression readiness
# ---------------------------------------------------------------------------

#' Check whether an object is ready for differential-expression testing
#'
#' Reports whether an active group identity resolves to at least two
#' populations with enough cells to test, and whether a detest record already
#' exists. It deliberately does not choose a group column, a clustering
#' resolution or a marker threshold: those are user decisions.
#'
#' @param object A \code{SingleCellExperiment} object.
#' @param group Optional group column name. When NULL it resolves from the
#'   active identity.
#' @return A read-only readiness summary list.
#' @export
check_de_readiness <- function(object, group = NULL) {
    if (!inherits(object, "SingleCellExperiment")) {
        stop("object must be a SingleCellExperiment", call. = FALSE)
    }
    assay_names <- SummarizedExperiment::assayNames(object)
    counts_available <- "counts" %in% assay_names
    columns <- colnames(SummarizedExperiment::colData(object))

    group_key <- group %||% tryCatch(ActiveIdent(object), error = function(e) NULL)
    if (!is.null(group_key) && !as.character(group_key) %in% columns &&
            !identical(as.character(group_key), "colLabels")) {
        group_key <- NULL
    }
    if (is.null(group_key)) {
        candidates <- grep("cluster|leiden|louvain|ident", columns,
            ignore.case = TRUE, value = TRUE)
        if (length(candidates)) group_key <- candidates[[1L]]
    }

    n_groups <- 0L
    smallest <- 0L
    if (!is.null(group_key)) {
        values <- tryCatch(
            if (identical(as.character(group_key), "colLabels"))
                as.character(SingleCellExperiment::colLabels(object))
            else as.character(SummarizedExperiment::colData(object)[[group_key]]),
            error = function(e) NULL
        )
        if (!is.null(values)) {
            counts <- as.integer(table(values[!is.na(values)]))
            n_groups <- as.integer(length(counts))
            smallest <- if (length(counts)) min(counts) else 0L
        }
    }
    group_resolvable <- !is.null(group_key) && n_groups > 0L
    enough_groups <- n_groups >= 2L
    enough_cells <- smallest >= 2L

    records <- tryCatch(sclet_get_state_records(object, "detest"), error = function(e) list())
    if (is.null(records)) records <- list()
    detest_ids <- names(records)
    de_package_available <- base::requireNamespace("limma", quietly = TRUE)

    questions <- character()
    if (!group_resolvable) {
        questions <- c(questions, paste0(
            "No group identity is resolvable, and a marker test compares groups against each other; ",
            "cluster the data or set an active identity first. Available colData columns: ",
            paste(columns, collapse = ", "), "."
        ))
    } else if (!enough_groups) {
        questions <- c(questions, paste0(
            "The resolved grouping '", group_key, "' has a single population, so every marker test would ",
            "compare one group against nothing and return an empty table. Re-cluster so at least two ",
            "populations exist; no resolution or alternative column is suggested here."
        ))
    }
    if (group_resolvable && enough_groups && !enough_cells) {
        questions <- c(questions, paste0(
            "The smallest population holds ", smallest, " cell(s); a two-group comparison with fewer than ",
            "two observations per group has no within-group variance to test, so a marker table can come ",
            "back empty or degenerate rather than erroring."
        ))
    }
    if (!counts_available) {
        questions <- c(questions, paste0(
            "No 'counts' assay is present; available assays: ",
            if (length(assay_names)) paste(assay_names, collapse = ", ") else "none", "."
        ))
    }
    if (!de_package_available) {
        questions <- c(questions,
            "The 'limma' package is not installed, which the marker test relies on.")
    }

    ready <- group_resolvable && enough_groups && enough_cells && de_package_available
    list(
        status = if (ready) "ready_for_diagnostic" else "not_ready",
        checks = list(
            counts_assay_available = counts_available,
            available_assays = assay_names,
            group_column_resolvable = group_resolvable,
            group_column_resolved = group_key,
            n_groups = n_groups,
            smallest_group_n = smallest,
            detest_record_available = length(detest_ids) > 0L,
            available_detest_ids = detest_ids,
            de_package_available = de_package_available,
            raw_values_included = FALSE
        ),
        blocked_actions = if (ready) character() else "de",
        questions = questions,
        notes = paste(
            "No group column, clustering resolution or significance threshold is chosen here. Marker",
            "results are statistical associations under the selected test and thresholds, not measured",
            "quantities, and they are the input to enrichment rather than a second independent line of",
            "evidence for it. No cell is removed or filtered by this diagnostic."
        )
    )
}

# ---------------------------------------------------------------------------
# C3: Milo (differential abundance) readiness and evidence
# ---------------------------------------------------------------------------

# A design confirmation the user already made, restricted to the roles an
# action needs. Returns the role->column map or NULL when nothing on this
# object both covers those roles and still matches its recorded value
# signature. This is the user's own declaration being read back, never a guess
# about which column "is the sample".
sclet_ai_design_declaration <- function(object, roles) {
    records <- tryCatch(sclet_get_state_records(object, "ai_design_confirmation"),
        error = function(e) list())
    if (!length(records)) return(NULL)
    design_fp <- as.character(sclet_ai_design_fingerprint(object))
    for (id in names(records)) {
        summary <- records[[id]]$summary %||% list()
        if (!identical(as.character(summary$object_fingerprint %||% ""), design_fp)) next
        design_map <- summary$design %||% list()
        if (!all(roles %in% names(design_map))) next
        if (!is.null(summary$design_value_key)) {
            current <- sclet_ai_design_value_key(object, design_map)
            if (is.null(current) || !identical(current, summary$design_value_key)) next
        }
        return(design_map[roles])
    }
    NULL
}

# Milo's design roles come from the caller or from a confirmation the user
# already recorded -- never from us pattern-matching column names.
sclet_ai_milo_design_roles <- function(object, sample_col = NULL, condition_col = NULL) {
    columns <- colnames(SummarizedExperiment::colData(object))
    clean <- function(value) {
        if (is.character(value) && length(value) == 1L && nzchar(value) && value %in% columns) {
            value
        }
    }
    caller_sample <- clean(sample_col)
    caller_condition <- clean(condition_col)
    declared <- sclet_ai_design_declaration(object, c("sample", "condition"))
    sample <- caller_sample %||% declared$sample
    condition <- caller_condition %||% declared$condition
    confirmed <- if (!is.null(sample) && !is.null(condition)) {
        !is.null(tryCatch(
            sclet_ai_find_design_confirmation(object, design = list(
                sample = sample, condition = condition)),
            error = function(e) NULL))
    } else {
        FALSE
    }
    list(
        sample = sample,
        condition = condition,
        caller_supplied = !is.null(caller_sample) && !is.null(caller_condition),
        confirmed = confirmed
    )
}

#' Check whether an object is ready for differential-abundance testing
#'
#' Reports whether miloR is installed, whether a PCA embedding exists, and
#' whether the sample and condition design roles are known -- supplied here or
#' declared through an earlier \code{ConfirmAIDesignSemantics()} call. It never
#' guesses which colData column plays which design role and never recommends
#' one: that mapping is the user's own statement about the experiment, which is
#' what the confirmation gate exists to pin down.
#'
#' @param object A \code{SingleCellExperiment} object.
#' @param sample_col Optional column name the caller intends as the sample.
#' @param condition_col Optional column name the caller intends as the
#'   condition.
#' @return A read-only readiness summary list.
#' @export
check_milo_readiness <- function(object, sample_col = NULL, condition_col = NULL) {
    if (!inherits(object, "SingleCellExperiment")) {
        stop("object must be a SingleCellExperiment", call. = FALSE)
    }
    columns <- colnames(SummarizedExperiment::colData(object))
    reductions <- tryCatch(SingleCellExperiment::reducedDimNames(object),
        error = function(e) character())
    pca_available <- "PCA" %in% reductions
    milo_available <- base::requireNamespace("miloR", quietly = TRUE)

    roles <- sclet_ai_milo_design_roles(object, sample_col, condition_col)
    design_known <- !is.null(roles$sample) && !is.null(roles$condition)
    n_samples <- 0L
    n_conditions <- 0L
    if (design_known) {
        coldata <- SummarizedExperiment::colData(object)
        n_samples <- length(unique(as.character(coldata[[roles$sample]])))
        n_conditions <- length(unique(as.character(coldata[[roles$condition]])))
    }
    # Milo attributes neighbourhood counts to samples, so fewer than two of
    # either factor leaves nothing to contrast. Two is a floor, not advice.
    design_testable <- design_known && n_samples >= 2L && n_conditions >= 2L

    records <- tryCatch(sclet_get_state_records(object, "milo"), error = function(e) list())
    if (is.null(records)) records <- list()
    milo_ids <- names(records)

    questions <- character()
    if (!milo_available) {
        questions <- c(questions,
            "The 'miloR' package is not installed; install it via BiocManager::install('miloR').")
    }
    if (!pca_available) {
        questions <- c(questions, paste0(
            "No PCA embedding is available and Milo builds its neighbourhood graph on one; run RunPCA() ",
            "first. Available reductions: ",
            if (length(reductions)) paste(reductions, collapse = ", ") else "none", "."
        ))
    }
    if (!design_known) {
        missing_roles <- c(if (is.null(roles$sample)) "sample" else NULL,
            if (is.null(roles$condition)) "condition" else NULL)
        questions <- c(questions, paste0(
            "The design role(s) ", paste(missing_roles, collapse = " and "), " are not declared for this ",
            "object. Declare them yourself with ConfirmAIDesignSemantics(object, design = list(...)), or ",
            "pass the column names explicitly. No column is guessed or recommended here. Available ",
            "colData columns: ", paste(columns, collapse = ", "), "."
        ))
    } else if (!design_testable) {
        questions <- c(questions, paste0(
            "The declared design has ", n_samples, " sample(s) and ", n_conditions,
            " condition(s); differential abundance needs at least two of each, or there is no ",
            "replication to attribute neighbourhood counts to."
        ))
    } else if (roles$caller_supplied && !roles$confirmed) {
        # The columns are known and testable, but run_milo still refuses
        # without a confirmation, so readiness must refuse to say "ready": a
        # diagnostic that green-lights an action the gate will block is worse
        # than one that blocks, because the planner acts on the green light.
        questions <- c(questions, paste0(
            "The sample and condition roles resolve to real columns but were never confirmed as design ",
            "roles; call ConfirmAIDesignSemantics(object, design = list(sample = '", roles$sample,
            "', condition = '", roles$condition, "')) so that mapping is your own declaration."
        ))
    }

    ready <- milo_available && pca_available && design_testable && isTRUE(roles$confirmed)
    list(
        status = if (ready) "ready_for_diagnostic" else "not_ready",
        checks = list(
            milo_package_available = milo_available,
            pca_reduction_available = pca_available,
            available_reductions = reductions,
            sample_column_resolved = roles$sample,
            condition_column_resolved = roles$condition,
            design_roles_confirmed = isTRUE(roles$confirmed),
            design_roles_supplied_by_caller = isTRUE(roles$caller_supplied),
            n_samples = as.integer(n_samples),
            n_conditions = as.integer(n_conditions),
            milo_record_available = length(milo_ids) > 0L,
            available_milo_ids = milo_ids,
            raw_values_included = FALSE
        ),
        blocked_actions = if (ready) character() else "milo",
        questions = questions,
        notes = paste(
            "Neither the sample column nor the condition column is guessed: design roles come from the",
            "caller or from a prior ConfirmAIDesignSemantics() declaration, because that mapping is the",
            "user's own statement about the experiment. Differential abundance results are associations",
            "between neighbourhoods and conditions under a chosen graph and test, not causal or measured",
            "effects, and no per-neighbourhood or per-cell value is returned here."
        )
    )
}

# A count that cannot be determined is reported as NA, not 0: a zero reads as
# the measured verdict "no samples observed", while NA says "this was not
# recorded", which is the honest statement when a record omits it.
sclet_ai_milo_count_or_na <- function(value) {
    if (is.null(value) || !length(value)) return(NA_integer_)
    number <- suppressWarnings(as.integer(value[[1L]]))
    if (is.na(number)) NA_integer_ else number
}

# Distinct levels of one design column, or NA when the column is unknown. The
# labels themselves never leave here -- only the count does.
sclet_ai_milo_level_count <- function(object, column) {
    if (is.null(column) || !length(column)) return(NA_integer_)
    coldata <- tryCatch(SummarizedExperiment::colData(object), error = function(e) NULL)
    if (is.null(coldata) || !column %in% colnames(coldata)) return(NA_integer_)
    values <- tryCatch(as.character(coldata[[column]]), error = function(e) NULL)
    if (is.null(values)) return(NA_integer_)
    as.integer(length(unique(values[!is.na(values)])))
}

# Every Milo node on an object shares a dependency group: the runs use one
# neighbourhood graph over one design, so a second contrast is not a second
# independent line of support.
sclet_ai_milo_dependency_group <- function(object, roles) {
    coldata <- SummarizedExperiment::colData(object)
    sizes <- integer(0)
    for (column in c(roles$sample, roles$condition)) {
        if (!is.null(column) && column %in% colnames(coldata)) {
            sizes <- c(sizes, as.integer(sort(table(as.character(coldata[[column]])))))
        }
    }
    if (!length(sizes)) sizes <- as.integer(ncol(object))
    payload <- paste0(sizes, collapse = ",")
    if (base::requireNamespace("digest", quietly = TRUE)) {
        return(paste0("milo_", substring(digest::digest(payload, algo = "sha256"), 1L, 32L)))
    }
    paste0("milo_", substring(payload, 1L, 32L))
}

# Bounded facts from a Milo record. Only numeric test columns of the DA table
# are read; design_df, the formula and every sample/condition label stay put.
sclet_ai_milo_table_facts <- function(record) {
    da <- record$artifacts$da_results
    if (is.null(da)) return(NULL)
    table <- tryCatch(as.data.frame(da), error = function(e) NULL)
    if (is.null(table)) return(NULL)
    columns <- tolower(colnames(table))
    fdr_column <- colnames(table)[columns %in% c("fdr", "adj.p.value", "padj")][1L]
    p_column <- colnames(table)[columns %in% "pvalue"][1L]
    fdr <- if (!is.na(fdr_column)) suppressWarnings(as.numeric(table[[fdr_column]])) else numeric()
    fdr <- fdr[is.finite(fdr)]
    p_value <- if (!is.na(p_column)) suppressWarnings(as.numeric(table[[p_column]])) else numeric()
    p_value <- p_value[is.finite(p_value)]
    fdr_quantiles <- if (length(fdr)) {
        round(unname(stats::quantile(fdr, c(0.25, 0.5, 0.75))), 8L)
    } else {
        rep(NA_real_, 3L)
    }
    list(
        n_da_results = as.integer(nrow(table)),
        n_significant = as.integer(sum(fdr < 0.05)),
        fdr_quantile_25 = as.numeric(fdr_quantiles[[1L]]),
        fdr_median = as.numeric(fdr_quantiles[[2L]]),
        fdr_quantile_75 = as.numeric(fdr_quantiles[[3L]]),
        p_value_median = if (length(p_value)) round(stats::median(p_value), 8L) else NA_real_
    )
}

#' Record a bounded evidence node for a differential-abundance run
#'
#' Registers neighbourhood counts and FDR quantiles only. The design frame, the
#' model formula and any sample or condition label stay in the state record:
#' those are the user's experimental metadata, and an evidence node is eligible
#' to leave the process as provider context. Every Milo node on an object shares
#' a dependency group, so a second contrast over the same graph is not treated
#' as a second independent line of support.
#'
#' @param object A \code{SingleCellExperiment} object.
#' @param name Milo record id.
#' @return The object with one \code{ev:milo_<name>} node appended, or
#'   unchanged when no Milo record exists.
#' @keywords internal
#' @noRd
sclet_ai_record_milo_evidence <- function(object, name, replace = FALSE) {
    record <- tryCatch(sclet_get_state_record(object, "milo", name), error = function(e) NULL)
    if (is.null(record)) return(object)
    facts <- sclet_ai_milo_table_facts(record)
    if (is.null(facts)) return(object)
    roles <- sclet_ai_milo_design_roles(object,
        sample_col = record$inputs$sample_col,
        condition_col = record$inputs$condition_col)
    record_summary <- record$summary %||% list()
    covariates <- record$inputs$covariates
    evidence <- list(
        id = paste0("ev:milo_", name),
        kind = "deterministic_summary",
        values = c(facts, list(
            # NA when the role is unknown: 0 would read as a measured verdict
            # of "no conditions" instead of "this could not be determined".
            n_samples = sclet_ai_milo_count_or_na(record_summary$n_samples),
            n_conditions = sclet_ai_milo_level_count(object, roles$condition),
            n_contrasts = sclet_ai_milo_count_or_na(length(record$artifacts$contrasts %||% character())),
            has_covariates = length(as.character(covariates %||% character())) > 0L,
            neighbourhood_estimates_are_model_based = TRUE,
            raw_values_included = FALSE
        )),
        claim_level = "associated",
        dependency_group = sclet_ai_milo_dependency_group(object, roles)
    )
    # Same stale-node trap as the doublet recorder: route the re-run through
    # the shared recorder, which refreshes on request and otherwise refuses.
    sclet_ai_record_evidence_for_run(object, evidence, source = name,
        parents = character(), scope = NULL, replace = replace)
}

#' Summarize Milo results in bounded form without claiming a conclusion
#'
#' Read-only companion to \code{RunMilo()}: it reports the same neighbourhood
#' counts and FDR quantiles the evidence node carries and writes nothing, so an
#' agent can reason about a differential-abundance run, or about whether one is
#' even possible, without executing anything. It never ranks neighbourhoods,
#' never names a sample or condition label, and never recommends a design column,
#' a contrast or a graph parameter.
#'
#' @param object A \code{SingleCellExperiment} object.
#' @param id Optional Milo record id. When \code{NULL}, the active one.
#' @return A read-only summary list.
#' @export
summarize_milo_evidence <- function(object, id = NULL) {
    if (!inherits(object, "SingleCellExperiment")) {
        stop("object must be a SingleCellExperiment", call. = FALSE)
    }
    unavailable <- function(reason) {
        list(
            status = "not_available", reason = reason, n_records = 0L,
            claim_level = "associated", recommendation = NULL,
            raw_values_included = FALSE,
            caveats = "No differential abundance is reported and no design label is named."
        )
    }
    records <- tryCatch(sclet_get_state_records(object, "milo"), error = function(e) list())
    if (is.null(records)) records <- list()
    if (!length(records)) return(unavailable("no_milo_record"))

    resolved_id <- if (!is.null(id) && nzchar(id)) as.character(id)[[1L]] else {
        tryCatch(sclet_get_active_state(object, "milo"), error = function(e) NULL)
    }
    if (is.null(resolved_id) || !resolved_id %in% names(records)) {
        return(unavailable("no_milo_record"))
    }
    record <- records[[resolved_id]]
    facts <- sclet_ai_milo_table_facts(record)
    if (is.null(facts)) return(unavailable("no_milo_results"))

    roles <- sclet_ai_milo_design_roles(object,
        sample_col = record$inputs$sample_col,
        condition_col = record$inputs$condition_col)
    record_summary <- record$summary %||% list()
    list(
        status = "available",
        milo_id = resolved_id,
        available_milo_ids = names(records),
        n_da_results = facts$n_da_results,
        n_significant = facts$n_significant,
        fdr_quantile_25 = facts$fdr_quantile_25,
        fdr_median = facts$fdr_median,
        fdr_quantile_75 = facts$fdr_quantile_75,
        p_value_median = facts$p_value_median,
        n_samples = sclet_ai_milo_count_or_na(record_summary$n_samples),
        n_conditions = sclet_ai_milo_level_count(object, roles$condition),
        n_contrasts = sclet_ai_milo_count_or_na(length(record$artifacts$contrasts %||% character())),
        evidence_node_recorded = paste0("ev:milo_", resolved_id) %in%
            names(tryCatch(sclet_ai_evidence_get_all(object), error = function(e) list())),
        claim_level = "associated",
        recommendation = NULL,
        raw_values_included = FALSE,
        caveats = paste(
            "Differential abundance results are associations between neighbourhoods and conditions",
            "under a chosen graph, library-size normalisation and test; they are not causal effects and",
            "not measured quantities. A non-significant neighbourhood is not evidence of no difference,",
            "because the graph and the FDR procedure both depend on the design the user declared. Only",
            "counts and quantiles are reported: no per-neighbourhood value, no design frame contents",
            "and no sample or condition label is included, and the full DA table stays in the state",
            "record for R-side access."
        )
    )
}
