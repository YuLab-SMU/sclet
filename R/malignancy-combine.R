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
#' The malignant-cell special case of [CallPhenotypeCells()]: a two-line call
#' where the copy-number line (see [RunCNVInference()]) is treated as independent
#' of the transcriptional signature, because genomic instability is not a
#' transcriptional identity programme.
#'
#' @param object A `SingleCellExperiment` object.
#' @param cnv_id Analysis record id of a completed [RunCNVInference()] run
#'   (default `"cnv"`).
#' @param signature A `c(malignant = ...)` list or the output of
#'   [FindMalignantSeedGenes()] whose `tumor_genes` is used.
#' @param signature_score Optional `colData` column holding a malignant
#'   signature score. When supplied, `signature` is not scored again.
#' @param transcriptomic_cut Optional numeric. Cells scoring at or above it are
#'   "signature positive". `NULL` uses the score's median.
#' @param cnv_cut Optional numeric cutoff on CNV burden. `NULL` derives it from
#'   the reference cells (see [sclet_cnv_reference_cut()]).
#' @param reference_quantile Numeric in (0, 1). Reference cells' quantile used as
#'   the CNV baseline. Defaults to 0.99.
#' @param identity_genes Optional character vector naming a competing identity
#'   programme (for a gastric tumour, the intestinal-metaplasia profile). It is
#'   passed to [CallPhenotypeCells()] as a competing hypothesis.
#' @param identity_demote Numeric. Separation above which a call is demoted.
#' @param assay Assay for scoring. Defaults to `"logcounts"`.
#' @param name Analysis record id. Defaults to `"malignant_call"`.
#'
#' @return `object` with `colData` columns `<name>_state`, `<name>_evidence`,
#'   `<name>_signature_score`, `<name>_identity_score` (when `identity_genes` is
#'   supplied), `<name>_cnv_state`, plus the generic columns written by
#'   [CallPhenotypeCells()]. The state vocabulary is the four malignant states,
#'   mapped from the generic call:
#'
#' \describe{
#'   \item{`malignant`}{`supported`: CNV and signature both positive.}
#'   \item{`likely_malignant`}{`probable`: exactly one line positive.}
#'   \item{`undetermined`}{`undetermined`: conflict, CNV unavailable, or a
#'     competing identity programme explains the split.}
#'   \item{`non_malignant`}{`refuted` or both lines negative.}
#' }
#'
#' @section What this does not establish:
#' No per-cell ground truth exists for malignant identity, so this is a
#' classification from two observations, not a verification. `cnv_cut` and
#' `transcriptomic_cut` are analysis parameters, not discovered thresholds.
#'
#' @seealso [CallPhenotypeCells()], [RunCNVInference()], [check_cnv_readiness()]
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
    if (!is.null(signature_score) && !signature_score %in% colnames(cdata)) {
        stop("`signature_score` column '", signature_score, "' not found in colData.")
    }
    # The malignant call is defined as a TWO-line call. Running it with no CNV
    # record would silently degrade it to a signature-only call, so say so here
    # rather than producing a weaker answer under the stronger name.
    cnv_rec <- sclet_get_analysis(object, "malignant_cnv")
    if (is.null(cnv_rec) || is.null(cnv_rec[[cnv_id]])) {
        stop(
            "No CNV record '", cnv_id, "' found. Run RunCNVInference(object, name = \"",
            cnv_id, "\") first, or pass `cnv_id` for the record you ran."
        )
    }

    competing <- NULL
    if (!is.null(identity_genes) && length(identity_genes) > 0) {
        competing <- stats::setNames(list(identity_genes), "identity")
    }

    object <- CallPhenotypeCells(
        object = object,
        signature = signature,
        signature_score = signature_score,
        # Only the CNV line: `module` shares the expression matrix with the
        # signature and would not be a second observation (see
        # PhenotypeEvidenceProviders()).
        evidence = "cnv",
        transcriptomic_cut = transcriptomic_cut,
        competing = competing,
        competing_demote = identity_demote,
        assay = assay,
        name = name,
        cnv_id = cnv_id,
        cnv_cut = cnv_cut,
        reference_quantile = reference_quantile
    )

    # ---- map the generic call onto the malignant vocabulary ------------------
    cdata <- SummarizedExperiment::colData(object)
    # One mapping, stated once. The generic layer already guarantees that
    # `refuted` means "an explicit negative and no positive", so it cannot be
    # reached with zero lines -- there is no second rule to overlay.
    generic <- as.character(cdata[[paste0(name, "_call")]])
    state <- ifelse(generic == "supported", "malignant",
                    ifelse(generic == "probable", "likely_malignant",
                           ifelse(generic == "refuted", "non_malignant", "undetermined")))

    evidence <- paste0(name, "_evidence")
    cnv_state <- cdata[[paste0(name, "_cnv_state")]]
    # Use the cutoff the generic layer actually applied, not a fresh median:
    # recomputing here would disagree with the `call` column whenever the caller
    # supplied a cutoff.
    sig_cut <- sclet_get_analysis(object, "phenotype_call")[[name]]$signature_cut$value
    sig_score <- cdata[[paste0(name, "_signature_score")]]
    sig_pos <- !is.na(sig_score) & sig_score >= sig_cut
    cdata[[paste0(name, "_state")]] <- factor(state, levels = c("malignant", "likely_malignant",
                                                                "undetermined", "non_malignant"))
    cdata[[evidence]] <- ifelse(is.na(cnv_state) | cnv_state == "unavailable",
                                ifelse(sig_pos, "signature_only", "none"),
                                ifelse(cnv_state == "positive" & sig_pos, "cnv_signature",
                                       ifelse(cnv_state == "positive", "cnv_only",
                                              ifelse(sig_pos, "signature_only", "neither"))))
    SummarizedExperiment::colData(object) <- cdata

    counts <- table(cdata[[paste0(name, "_state")]])
    pheno <- sclet_get_analysis(object, "phenotype_call")
    pheno_rec <- if (!is.null(pheno)) pheno[[name]] else NULL
    demoted <- if (!is.null(pheno_rec$competing)) {
        any(vapply(pheno_rec$competing, function(x) isTRUE(x$demote), logical(1)))
    } else {
        FALSE
    }
    object <- sclet_set_state_record(
        object = object,
        type = "annotation",
        id = name,
        value = list(
            method = "four-state malignant-cell call from CNV + signature evidence",
            inputs = list(cnv_id = cnv_id),
            artifacts = list(
                columns = paste0(name, c("_state", "_evidence", "_signature_score"))
            ),
            summary = list(
                n_malignant = sclet_state_count(counts, "malignant"),
                n_likely_malignant = sclet_state_count(counts, "likely_malignant"),
                n_undetermined = sclet_state_count(counts, "undetermined"),
                n_non_malignant = sclet_state_count(counts, "non_malignant"),
                cnv_available = !all(cnv_state == "unavailable"),
                identity_conflict = demoted
            ),
            created_at = Sys.time()
        )
    )
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
