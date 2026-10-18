# Combining evidence into a malignant-cell call
#
# Neither line of evidence in sclet is sufficient alone:
#
#   * expression alone finds the largest transcriptional split, which on tumour
#     tissue can be a metaplastic or intestinal identity programme rather than
#     malignancy (the failure documented in YuLab-SMU/ProjectYulab#26);
#   * copy number alone misses copy-number-quiet but hypermutated tumours and
#     is noisy in low-depth cells.
#
# The two are close to orthogonal -- one reads a transcriptional programme, the
# other reads genomic structure -- so requiring both, and naming the cases where
# they disagree, is a stronger statement than either alone. `CallMalignantCells()`
# is the single entry point that does this and it is the package's answer to the
# issue: not "one side of a split", but /which cells carry two independent kinds
# of evidence, and which are undecided/.

#' Call malignant cells from independent lines of evidence
#'
#' Combines a copy-number measurement (see [RunCNVInference()]) with a
#' transcriptional signature score and reports, for every cell, a four-state
#' malignancy call with its supporting evidence.
#'
#' @param object A `SingleCellExperiment` object.
#' @param cnv_id Analysis record id of a completed [RunCNVInference()] run
#'   (default `"cnv"`). When the run reported `status = "unavailable"`, the call
#'   degrades to expression-only and every cell that cannot be decided becomes
#'   `"undetermined"` rather than being called malignant on one line of evidence.
#' @param signature A `c(malignant = ...)` list or the output of
#'   [FindMalignantSeedGenes()] whose `tumor_genes` is used.
#' @param signature_score Optional `colData` column already holding a malignant
#'   signature score. When supplied, `signature` is not scored again.
#' @param transcriptomic_cut Optional numeric. Cells scoring at or above it are
#'   "signature positive". `NULL` uses the score's median.
#' @param cnv_cut Optional numeric cutoff on CNV burden. `NULL` puts the cutoff
#'   at the reference cells' upper quantile plus three robust standard
#'   deviations, computed from `reference_quantile`.
#' @param reference_quantile Numeric in (0, 1). Reference cells' quantile used as
#'   the CNV baseline. Defaults to 0.99.
#' @param identity_genes Optional character vector naming a competing identity
#'   programme (for a gastric tumour, the intestinal-metaplasia profile). When
#'   supplied, the malignant side's mean identity score is compared to the
#'   non-malignant side's and the call is demoted if the programme tracks the
#'   malignant side.
#' @param identity_demote Numeric. Mean difference on `identity_genes` above
#'   which a cell found malignant is demoted to `"undetermined"`.
#' @param assay Assay for scoring. Defaults to `"logcounts"`.
#' @param name Analysis record id. Defaults to `"malignant_call"`.
#'
#' @return `object` with `colData` columns `<name>_state`, `<name>_evidence`,
#'   `<name>_signature_score` and `<name>_identity_score`, and a state record of
#'   type `annotation` under `name`. `<name>_state` is a factor with levels
#'   `malignant`, `likely_malignant`, `undetermined` and `non_malignant`.
#'
#' @section The four states:
#' \describe{
#'   \item{`malignant`}{CNV-positive and signature-positive. Two independent
#'     lines agree.}
#'   \item{`likely_malignant`}{Exactly one line is positive, usually because the
#'     other could not be measured. Reported separately so that a one-line call is
#'     never read as a two-line one.}
#'   \item{`undetermined`}{Evidence conflicts (CNV-positive but
#'     signature-negative, or the reverse), the CNV layer did not run, or a
#'     competing identity programme explains the split.}
#'   \item{`non_malignant`}{CNV-negative and signature-negative.}
#' }
#'
#' @section What this does not establish:
#' No per-cell ground truth exists for malignant identity, so this is a
#' classification from two observations, not a verification. `cnv_cut` and
#' `transcriptomic_cut` are analysis parameters, not discovered thresholds, and
#' moving them moves the answer. Report the state distribution and the
#' disagreement cells alongside any conclusion.
#'
#' @seealso [RunCNVInference()], [check_cnv_readiness()], [FindMalignantSeedGenes()]
#' @export
CallMalignantCells <- function(object,
                               cnv_id = "cnv",
                               signature = NULL,
                               signature_score = NULL,
                               transcriptomic_cut = NULL,
                               cnv_cut = NULL,
                               reference_quantile = 0.99,
                               identity_genes = NULL,
                               identity_demote = 0.1,
                               assay = NULL,
                               name = "malignant_call") {
    if (!nzchar(name)) {
        stop("`name` must be a non-empty prefix.")
    }
    cdata <- SummarizedExperiment::colData(object)
    col <- function(suffix) paste0(name, suffix)

    # ---- CNV line ------------------------------------------------------------
    cnv_rec <- sclet_get_analysis(object, "malignant_cnv")
    if (is.null(cnv_rec) || is.null(cnv_rec[[cnv_id]])) {
        stop(
            "No CNV record '", cnv_id, "' found. Run RunCNVInference(object, name = \"",
            cnv_id, "\") first, or pass `cnv_id` for the record you ran."
        )
    }
    cnv_rec <- cnv_rec[[cnv_id]]
    cnv_available <- identical(cnv_rec$status, "completed") &&
        paste0(cnv_id, "_burden") %in% colnames(cdata)
    burden <- if (cnv_available) cdata[[paste0(cnv_id, "_burden")]] else
        rep(NA_real_, nrow(cdata))
    is_ref <- if (cnv_available) cdata[[paste0(cnv_id, "_is_reference")]] else
        rep(FALSE, nrow(cdata))

    # ---- signature score -----------------------------------------------------
    if (!is.null(signature_score)) {
        if (!signature_score %in% colnames(cdata)) {
            stop("`signature_score` column '", signature_score, "' not found in colData.")
        }
        sig_score <- as.numeric(cdata[[signature_score]])
        signature_used <- list(status = "supplied", column = signature_score)
    } else {
        genes <- sclet_call_resolve_signature(signature)
        if (is.null(assay)) {
            assay <- if ("logcounts" %in% SummarizedExperiment::assayNames(object)) "logcounts" else NULL
        }
        if (is.null(assay)) {
            stop("No `logcounts` assay; pass `assay` or provide `signature_score`.")
        }
        expr <- sclet_malignancy_dense(SummarizedExperiment::assay(object, assay))
        bins <- sclet_gene_expression_bins(expr, n_bins = 25)
        sig_score <- sclet_add_module_score(expr, genes, bins = bins, n_ctrl = 100,
                                            exclude = genes)
        signature_used <- list(status = "scored", assay = assay, genes = length(genes))
    }

    # ---- identity line -------------------------------------------------------
    identity_score <- NULL
    if (!is.null(identity_genes) && length(identity_genes) > 0) {
        if (is.null(assay)) {
            assay <- if ("logcounts" %in% SummarizedExperiment::assayNames(object)) "logcounts" else NULL
        }
        if (!is.null(assay)) {
            expr <- sclet_malignancy_dense(SummarizedExperiment::assay(object, assay))
            present <- intersect(identity_genes, rownames(object))
            bins <- sclet_gene_expression_bins(expr, n_bins = 25)
            identity_score <- sclet_add_module_score(expr, present, bins = bins,
                                                     n_ctrl = 100, exclude = present)
        }
    }

    # ---- thresholds ----------------------------------------------------------
    sig_pos <- if (is.null(transcriptomic_cut)) {
        stats::median(sig_score, na.rm = TRUE)
    } else {
        transcriptomic_cut
    }
    if (cnv_available) {
        ref_burden <- burden[!is.na(is_ref) & is_ref]
        ref_burden <- ref_burden[!is.na(ref_burden)]
        if (length(ref_burden) >= 5) {
            cnv_pos <- if (!is.null(cnv_cut)) burden > cnv_cut else burden > sclet_cnv_reference_cut(ref_burden, reference_quantile)
        } else {
            cnv_pos <- rep(NA, length(burden))
        }
    } else {
        cnv_pos <- rep(NA, length(burden))
    }
    sig_ok <- !is.na(sig_score) & sig_score >= sig_pos

    # ---- identity conflict ---------------------------------------------------
    identity_flag <- NULL
    if (!is.null(identity_score) && cnv_available) {
        on <- !is.na(cnv_pos) & cnv_pos & sig_ok
        off <- !is.na(cnv_pos) & !cnv_pos
        if (sum(on) > 2 && sum(off) > 2) {
            diff_id <- mean(identity_score[on]) - mean(identity_score[off])
            identity_flag <- list(
                difference = diff_id,
                demote = abs(diff_id) > identity_demote,
                n_on = sum(on),
                n_off = sum(off)
            )
        }
    }

    # ---- four-state call ----------------------------------------------------
    state <- rep("undetermined", length(burden))
    if (!cnv_available) {
        # One line only: a positive signature cannot reach "malignant".
        state[sig_ok] <- "likely_malignant"
        state[!sig_ok] <- "undetermined"
    } else {
        both <- !is.na(cnv_pos) & cnv_pos & sig_ok
        neither <- !is.na(cnv_pos) & !cnv_pos & !sig_ok
        one <- !is.na(cnv_pos) & xor(cnv_pos, sig_ok)
        state[both] <- "malignant"
        state[neither] <- "non_malignant"
        state[one] <- "likely_malignant"
        state[is.na(cnv_pos)] <- "undetermined"
        if (!is.null(identity_flag) && isTRUE(identity_flag$demote)) {
            state[state == "malignant"] <- "undetermined"
        }
    }
    state_factor <- factor(state, levels = c("malignant", "likely_malignant",
                                             "undetermined", "non_malignant"))
    evidence <- ifelse(
        is.na(cnv_pos),
        ifelse(sig_ok, "signature_only", "none"),
        ifelse(cnv_pos & sig_ok, "cnv_signature",
               ifelse(cnv_pos, "cnv_only", ifelse(sig_ok, "signature_only", "neither")))
    )

    cdata[[col("_state")]] <- state_factor
    cdata[[col("_evidence")]] <- evidence
    cdata[[col("_signature_score")]] <- sig_score
    if (!is.null(identity_score)) {
        cdata[[col("_identity_score")]] <- identity_score
    }
    SummarizedExperiment::colData(object) <- cdata

    state_counts <- table(state_factor)
    result <- list(
        state = stats::setNames(as.list(as.integer(state_counts)), names(state_counts)),
        n_cells = length(state),
        cnv_available = cnv_available,
        cnv_status = if (cnv_available) "completed" else "unavailable",
        thresholds = list(
            transcriptomic = sig_pos,
            cnv_cut = if (!cnv_available) {
                NA_real_
            } else if (!is.null(cnv_cut)) {
                cnv_cut
            } else {
                ref_burden <- burden[!is.na(is_ref) & is_ref]
                if (length(ref_burden) >= 5) {
                    sclet_cnv_reference_cut(ref_burden, reference_quantile)
                } else {
                    NA_real_
                }
            },
            reference_quantile = reference_quantile
        ),
        signature = signature_used,
        identity = identity_flag,
        params = list(cnv_id = cnv_id, transcriptomic_cut = transcriptomic_cut,
                      cnv_cut = cnv_cut, reference_quantile = reference_quantile,
                      identity_demote = identity_demote, name = name)
    )

    existing <- sclet_get_analysis(object, "malignant_call")
    if (is.null(existing)) existing <- list()
    existing[[name]] <- result
    object <- sclet_set_analysis(object, "malignant_call", existing)
    md <- S4Vectors::metadata(object)
    md$malignant_call <- if (is.null(md$malignant_call)) list() else md$malignant_call
    md$malignant_call[[name]] <- result
    S4Vectors::metadata(object) <- md

    object <- sclet_set_state_record(
        object = object,
        type = "annotation",
        id = name,
        value = list(
            method = "four-state malignant-cell call from CNV + signature evidence",
            inputs = list(cnv_id = cnv_id, cnv_status = result$cnv_status,
                          signature = signature_used$status),
            artifacts = list(
                columns = paste0(name, c("_state", "_evidence", "_signature_score")),
                metadata = paste0("malignant_call$", name)
            ),
            params = result$params,
            summary = list(
                n_malignant = sclet_state_count(state_counts, "malignant"),
                n_likely_malignant = sclet_state_count(state_counts, "likely_malignant"),
                n_undetermined = sclet_state_count(state_counts, "undetermined"),
                n_non_malignant = sclet_state_count(state_counts, "non_malignant"),
                cnv_available = cnv_available,
                identity_conflict = !is.null(identity_flag) && isTRUE(identity_flag$demote)
            ),
            created_at = Sys.time()
        )
    )
    object <- sclet_log_command(
        object,
        "CallMalignantCells",
        params = result$params,
        outputs = list(columns = paste0(name, c("_state", "_evidence")),
                       metadata = paste0("malignant_call$", name))
    )

    if (!cnv_available) {
        warning(
            "CNV evidence is unavailable, so no cell is called `malignant`: ",
            "signature-positive cells are reported as `likely_malignant` and the ",
            "rest as `undetermined`. Fix the CNV layer to obtain a two-line call.",
            call. = FALSE
        )
    }
    if (!is.null(identity_flag) && isTRUE(identity_flag$demote)) {
        warning(sprintf(
            "The supplied `identity_genes` programme scores %.3f higher on the CNV+signature side than on the rest; cells that were both evidence-positive are downgraded to `undetermined` because the split may track that identity rather than malignancy.",
            identity_flag$difference
        ), call. = FALSE)
    }
    object
}

# ---------------------------------------------------------------------------
# internals
# ---------------------------------------------------------------------------

# CNV cut for a reference population.
#
# The baseline is the reference cells' `reference_quantile`, widened by a robust
# spread estimate AND by a floor strictly above the reference maximum. The floor
# matters: an exactly measured reference (or a reference whose values are equal)
# has mad = 0, and without the floor the reference cells themselves sit ON the
# cut and read as CNV-positive -- a threshold that cannot separate anything.
sclet_cnv_reference_cut <- function(ref_burden, reference_quantile = 0.99) {
    ref_burden <- ref_burden[!is.na(ref_burden)]
    if (length(ref_burden) == 0) {
        return(NA_real_)
    }
    base <- stats::quantile(ref_burden, reference_quantile, names = FALSE)
    spread <- stats::mad(ref_burden)
    candidate <- base + 3 * spread
    ceiling <- max(ref_burden)
    if (candidate > ceiling) {
        return(candidate)
    }
    # Not above every observed reference cell: step out by a small relative
    # margin above the reference maximum instead.
    ceiling + max(abs(ceiling) * 0.05, .Machine$double.eps)
}

sclet_state_count <- function(counts, level) {
    v <- unname(counts[[level]])
    if (length(v) == 0 || is.na(v)) 0L else as.integer(v)
}

sclet_call_resolve_signature <- function(signature) {
    if (is.null(signature)) {
        stop("Supply `signature` (gene vector, named list) or `signature_score`.")
    }
    if (is.character(signature)) {
        return(signature)
    }
    if (is.list(signature) && !is.null(signature$tumor_genes)) {
        return(signature$tumor_genes)
    }
    if (is.list(signature) && !is.null(names(signature))) {
        idx <- match(c("malignant", "tumor"), tolower(names(signature)))
        idx <- idx[!is.na(idx)]
        if (length(idx) == 0) {
            stop("A named `signature` list needs a 'malignant' or 'tumor' element.")
        }
        return(signature[[idx[[1]]]])
    }
    stop("`signature` must be a gene vector, a named list, or FindMalignantSeedGenes() output.")
}
