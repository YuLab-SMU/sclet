# Malignant-epithelium classification
#
# Reference method: Zhang P. et al., "Dissecting transcriptional heterogeneity in
# primary gastric adenocarcinoma by single cell RNA sequencing", Gut 2020,
# doi:10.1136/gutjnl-2019-320368 (issue YuLab-SMU/ProjectYulab#26).
#
# Step 1 (FindMalignantSeedGenes) derives tumour-up / normal-up seed signatures
# from paired bulk data.  Step 2 (ClassifyMalignantEpithelium) scores single
# cells with those seeds, splits them by k-means, re-derives signatures from the
# split (pseudobulk DE) and repeats until the labels stop moving.
#
# The method finds the largest transcriptional split in the cells supplied to
# it.  That split is not automatically malignancy: metaplastic or intestinal
# identity can dominate it, and after the first round the non-malignant
# reference is no longer "the normal samples" but "the other side of the split".
# The functions therefore report diagnostics (agreement history, malignant
# fraction, seed overlap, side-of-origin composition, stage markers) and never
# present the result as ground truth; outcomes that need interpretation raise a
# warning instead of being silently returned.

#' Bulk tumour versus normal signatures for malignant-cell scoring
#'
#' Derive the seed gene sets used by [ClassifyMalignantEpithelium()] from bulk
#' transcriptomes: the top `n_top` genes most highly expressed in tumour relative
#' to normal, and the top `n_top` most highly expressed in normal relative to
#' tumour.
#'
#' @param expr A gene x sample expression matrix on a log scale, with gene ids
#'   in `rownames` and sample ids in `colnames`.
#' @param sample_info A data frame with one row per sample. Must contain a
#'   column holding the sample ids that appear in `colnames(expr)` (`sample_col`,
#'   default `"sample"`), a `condition` column and a `patient` column. The
#'   `condition` column accepts `"Tumor"`/`"Normal"` or `"Tumour"`;
#'   `case_values` and `control_values` name the levels explicitly when the data
#'   use other words.
#' @param assay_note How the caller states that `expr` was built, e.g.
#'   `"log2(TPM+1)"` or `"logCPM"`. Recorded in the result. If `assay_note`
#'   says TPM, the log offset is detected and undone before filtering so that a
#'   `min_expr` cutoff is compared on the TPM scale.
#' @param n_top Integer. Signature size per direction.
#' @param padj Numeric. Adjusted p-value cutoff (Benjamini-Hochberg). The
#'   reference method uses 0.01 and applies no fold-change cutoff.
#' @param lfc Numeric or `NULL`. Optional |logFC| cutoff. `NULL` (default) keeps
#'   the reference method's behaviour of padj only.
#' @param min_expr,min_frac Numeric. Low-expression filter: keep genes above
#'   `min_expr` in at least `min_frac` of the paired samples. Both `NULL`
#'   disables the filter.
#' @param paired Logical. Use a paired design (`~ patient + condition`). Defaults
#'   to `TRUE` when at least two patients contribute both conditions, which is
#'   what the reference analysis was paired on.
#' @param sample_col Name of the `sample_info` column holding the sample ids that
#'   appear in `colnames(expr)`. Defaults to `"sample"`.
#' @param case_values Character vector naming the tumour levels of `condition`.
#'   Defaults to `"Tumor"`/`"Tumour"`/`"Cancer"`.
#' @param control_values Character vector naming the normal levels of
#'   `condition`. Defaults to `"Normal"`/`"Control"`/`"Adjacent"`.
#' @param min_samples Integer. Fewest samples per side required to run.
#' @param n_ctrl `AddModuleScore` control pool size, stored with the signature
#'   so that scoring reproduces the reference's control-binned offsets.
#' @param n_bins `AddModuleScore` expression bins, stored with the signature.
#'
#' @return A list with `tumor_genes`, `normal_genes` (character vectors of the
#'   same length as `n_top` unless the data are thin), `deg` (the full limma
#'   table, ordered by adjusted p-value) and `info` (counts, cutoffs, design,
#'   transform and any warnings raised). It is accepted as `seed` by
#'   [ClassifyMalignantEpithelium()].
#'
#' @seealso [ClassifyMalignantEpithelium()]
#' @export
FindMalignantSeedGenes <- function(
    expr,
    sample_info,
    assay_note = c("log2(TPM+1)", "logCPM", "log2(FPKM+1)", "other"),
    n_top = 50,
    padj = 0.01,
    lfc = NULL,
    min_expr = 1,
    min_frac = 0.2,
    paired = NULL,
    min_samples = 3,
    sample_col = "sample",
    case_values = NULL,
    control_values = NULL,
    n_ctrl = 100,
    n_bins = 25
) {
    if (missing(assay_note) || is.null(assay_note)) {
        assay_note <- "other"
    }
    if (length(assay_note) > 1) {
        assay_note <- assay_note[[1]]
    }
    if (!is.matrix(expr)) {
        expr <- as.matrix(expr)
    }
    if (!is.numeric(expr)) {
        stop("`expr` must be a numeric gene x sample matrix.")
    }
    if (is.null(rownames(expr)) || is.null(colnames(expr))) {
        stop("`expr` needs gene ids in rownames and sample ids in colnames.")
    }
    if (is.null(sample_info) || !is.data.frame(sample_info) && !is.matrix(sample_info)) {
        stop("`sample_info` must be a data frame with one row per bulk sample.")
    }
    sample_info <- as.data.frame(sample_info, stringsAsFactors = FALSE)
    if (!sample_col %in% colnames(sample_info)) {
        stop("`sample_info` has no column '", sample_col, "' holding the sample ids.")
    }
    for (col in c("condition", "patient")) {
        if (!col %in% colnames(sample_info)) {
            stop("`sample_info` has no '", col, "' column.")
        }
    }

    # ---- align samples -------------------------------------------------------
    keep <- intersect(sample_info[[sample_col]], colnames(expr))
    dropped <- setdiff(sample_info[[sample_col]], colnames(expr))
    if (length(keep) < min_samples) {
        stop(
            "Only ", length(keep), " bulk samples match between `expr` and ",
            "`sample_info` (need at least ", min_samples, ")."
        )
    }
    st <- sample_info[match(keep, sample_info[[sample_col]]), , drop = FALSE]
    expr <- expr[, keep, drop = FALSE]

    # ---- case / control levels ----------------------------------------------
    cond_raw <- as.character(st$condition)
    cond_key <- tolower(cond_raw)
    if (is.null(case_values)) {
        case_values <- unique(cond_raw[cond_key %in% c("tumor", "tumour", "cancer")])
    }
    if (is.null(control_values)) {
        control_values <- unique(cond_raw[cond_key %in% c("normal", "control", "adjacent")])
    }
    if (length(case_values) == 0 || length(control_values) == 0) {
        stop(
            "Could not tell tumour from normal in `condition` (values: ",
            paste(unique(cond_raw), collapse = ", "),
            "). Name them with `case_values` and `control_values`."
        )
    }
    is_case <- cond_raw %in% case_values
    is_ctrl <- cond_raw %in% control_values
    if (any(!is_case & !is_ctrl)) {
        stop(
            "`condition` has values that are neither case nor control: ",
            paste(unique(cond_raw[!is_case & !is_ctrl]), collapse = ", "), "."
        )
    }
    if (sum(is_case) < min_samples || sum(is_ctrl) < min_samples) {
        stop(
            "Need at least ", min_samples, " samples per side; have ",
            sum(is_case), " case and ", sum(is_ctrl), " control."
        )
    }

    patient <- as.character(st$patient)
    n_paired <- length(intersect(unique(patient[is_case]), unique(patient[is_ctrl])))
    if (is.null(paired)) {
        paired <- n_paired >= 2
    }

    # ---- transform / filter --------------------------------------------------
    notes <- character(0)
    tpm <- NULL
    offset <- sclet_bulk_log_offset(assay_note)
    if (!is.null(offset)) {
        tpm <- 2^expr - offset
        # A matrix already log2(TPM+1) re-entered through this branch would
        # carry many negative values; recover by trying the other offset.
        if (mean(tpm < 0, na.rm = TRUE) > 0.01) {
            alt <- if (offset == 0.001) 1 else 0.001
            alt_tpm <- 2^expr - alt
            if (mean(alt_tpm < 0, na.rm = TRUE) < mean(tpm < 0, na.rm = TRUE)) {
                tpm <- alt_tpm
                notes <- c(notes, sprintf(
                    "`assay_note` suggests a log offset of %s but %s fits the data better.", offset, alt
                ))
            }
        }
    }

    if (!is.null(tpm) && !is.null(min_expr) && !is.null(min_frac)) {
        keep_gene <- rowSums(tpm > min_expr, na.rm = TRUE) >= ceiling(min_frac * ncol(tpm))
        n_genes_in <- nrow(expr)
        expr <- expr[keep_gene, , drop = FALSE]
        filter_note <- sprintf("TPM > %g in >= %d%% of matched samples", min_expr, round(min_frac * 100))
    } else {
        n_genes_in <- nrow(expr)
        keep_gene <- rep(TRUE, nrow(expr))
        filter_note <- "none (no expression-scale cutoff applied)"
        if (!is.null(min_expr) || !is.null(min_frac)) {
            notes <- c(notes, "Low-expression filter skipped: it needs `assay_note` on the TPM scale.")
        }
    }
    if (nrow(expr) < 100) {
        stop("Only ", nrow(expr), " genes left after filtering; too few to build a signature.")
    }

    # ---- limma ---------------------------------------------------------------
    cond <- factor(ifelse(is_case, "Tumor", "Normal"), levels = c("Normal", "Tumor"))
    design <- if (isTRUE(paired)) {
        pat <- factor(patient, levels = unique(patient))
        stats::model.matrix(~ pat + cond)
    } else {
        stats::model.matrix(~ cond)
    }
    if (qr(design)$rank < ncol(design)) {
        notes <- c(notes, "Design is rank deficient; switching to a paired-free design.")
        design <- stats::model.matrix(~ cond)
        paired <- FALSE
    }
    fit <- limma::lmFit(expr, design)
    fit <- limma::eBayes(fit)
    deg <- limma::topTable(fit, coef = ncol(design), number = Inf, sort.by = "P")
    deg$gene <- rownames(deg)
    deg <- deg[, c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]

    # ---- select signatures ---------------------------------------------------
    ok <- !is.na(deg$adj.P.Val) & deg$adj.P.Val < padj
    if (!is.null(lfc)) {
        ok <- ok & abs(deg$logFC) >= lfc
    }
    sig <- deg[ok, , drop = FALSE]
    up <- sig$gene[sig$logFC > 0]
    dn <- sig$gene[sig$logFC < 0]
    if (length(up) < n_top || length(dn) < n_top) {
        notes <- c(notes, sprintf(
            "Only %d tumour-up and %d normal-up genes pass padj < %g%s; signatures are smaller than n_top = %d.",
            length(up), length(dn), padj,
            if (is.null(lfc)) "" else sprintf(" and |logFC| >= %g", lfc), n_top
        ))
    }
    tumor_genes <- utils::head(up, n_top)
    normal_genes <- utils::head(dn, n_top) # already ordered by -logFC
    if (length(tumor_genes) == 0 || length(normal_genes) == 0) {
        stop(
            "No gene passed padj < ", padj,
            if (is.null(lfc)) "" else paste0(" and |logFC| >= ", lfc),
            " in one direction; lower the cutoffs or check the condition assignment."
        )
    }

    info <- list(
        signature = list(
            tumor = tumor_genes,
            normal = normal_genes,
            n_ctrl = n_ctrl,
            n_bins = n_bins,
            assay_note = assay_note,
            transform = if (is.null(offset)) NULL else paste0("undo log offset ", offset, " for filtering only")
        ),
        design = list(
            paired = isTRUE(paired),
            n_case = sum(is_case), n_control = sum(is_ctrl),
            n_patients = length(unique(patient)),
            n_patients_paired = n_paired,
            formula = paste("~", paste(attr(design, "term.labels"), collapse = " + ")),
            n_samples_dropped = length(dropped)
        ),
        cutoffs = list(padj = padj, lfc = lfc, n_top = n_top,
                       min_expr = if (is.null(tpm)) NULL else min_expr,
                       min_frac = if (is.null(tpm)) NULL else min_frac),
        filter = filter_note,
        counts = list(
            n_genes_in = n_genes_in,
            n_genes_used = nrow(expr),
            n_pass_padj = nrow(sig),
            n_up = length(up), n_down = length(dn)
        ),
        warnings = notes
    )
    if (length(notes) > 0) {
        for (note in notes) {
            warning(note, call. = FALSE)
        }
    }

    list(
        tumor_genes = tumor_genes,
        normal_genes = normal_genes,
        deg = deg,
        info = info
    )
}

#' Classify malignant versus non-malignant epithelial cells by iterative scoring
#'
#' Implements the iterative malignant-epithelium classification of Zhang et al.
#' (Gut 2020) on a single-cell object: score every cell in `cells` with the
#' malignant and the non-malignant signature (an `AddModuleScore`-compatible
#' control-binned score), split
#' the cells into two groups with k-means on the two scores, re-derive each
#' signature from the split by pseudobulk differential expression, and repeat
#' until the labels stop changing. The first round is seeded from bulk data
#' (see [FindMalignantSeedGenes()]); the reference method used the initial bulk
#' scores because bulk data mix in non-epithelial cells.
#'
#' @param object A `SingleCellExperiment` object.
#' @param seed Either the output of [FindMalignantSeedGenes()] or a
#'   `c(tumor = ..., normal = ...)`-style list of two gene vectors. The
#'   `tumor` side is the malignant signature.
#' @param method `"iterative"` (default) runs the reference method ported here.
#'   `"cnv"` delegates to [CallMalignantCells()], sclet's own two-line call, so
#'   the two can be compared from one entry point.
#' @param cells Character vector of cell names to classify. Defaults to every
#'   cell. Supply the epithelial compartment; signatures are re-derived from
#'   these cells, so non-epithelial cells would enter the contrast.
#' @param sample_col Column of `colData(object)` naming the biological sample.
#'   The reference method uses the patient, because pseudobulk pairing and the
#'   side-of-origin diagnostic need one donor per column.
#' @param assay Assay scored. Defaults to `"logcounts"` when present, otherwise
#'   the active layer.
#' @param n_top Integer. Signature size per direction after the first round.
#' @param padj Numeric. Adjusted p-value cutoff for the re-derived signatures.
#' @param lfc Numeric or `NULL`. Optional |logFC| cutoff for re-derivation.
#' @param convergence Numeric in (0, 1]. Label agreement with the previous round
#'   at or above this value stops the iteration. The reference method calls the
#'   result stable; issue #26 suggested 0.9.
#' @param max_iter Integer. Hard cap on rounds.
#' @param min_cells_per_group Integer. Fewest cells per side for a split to be
#'   used.
#' @param min_patients Integer. Fewest paired patients with cells on both sides
#'   for the pseudobulk test.
#' @param n_ctrl,n_bins `AddModuleScore` control pool size and bin count.
#'   Taken from `seed$info$signature` when present.
#' @param name Character. Column prefix in `colData` and the analysis record id.
#' @param nstart Integer. k-means restarts.
#' @param seed_genes Optional character vector. Signature genes to track across
#'   rounds (e.g. an intestinal metaplasia profile) so that the result reports
#'   how much of the malignant side expresses them.
#' @param verbose Logical. Print one line per round.
#'
#' @return `object` with `colData` columns `<name>_malignant_score`,
#'   `<name>_nonmalignant_score`, `<name>_cluster` and `<name>_classification`
#'   (a factor with levels `malignant` and `non_malignant`), the analysis
#'   registered under `analyses[["malignancy"]]` and a state record of type
#'   `annotation` with id `name`. The full result (iteration log, final
#'   signatures, diagnostics) is at
#'   `metadata(object)$malignancy[[name]]`.
#'
#' @section What the result is and is not:
#' Each round splits the cells at the direction of largest transcriptional
#' difference between the two score distributions. On tumour tissue that
#' direction can be a metaplastic or intestinal programme rather than a
#' malignant one, and from the second round on the non-malignant side is "the
#' other side of the split", not the normal samples. A converged run therefore
#' means "the split stopped moving", not "the cells were proven malignant".
#' Validation needs an independent line of evidence (copy-number inference,
#' histology, or scoring against a signature derived elsewhere).
#'
#' @seealso [FindMalignantSeedGenes()]
#' @export
ClassifyMalignantEpithelium <- function(
    object,
    seed,
    method = c("iterative", "cnv"),
    cells = NULL,
    sample_col = "sample",
    assay = NULL,
    n_top = 50,
    padj = 0.01,
    lfc = NULL,
    convergence = 0.9,
    max_iter = 10,
    min_cells_per_group = 10,
    min_patients = 3,
    n_ctrl = NULL,
    n_bins = NULL,
    name = "malignant",
    nstart = 25,
    seed_genes = NULL,
    verbose = TRUE
) {
    method <- match.arg(method)
    if (identical(method, "cnv")) {
        # Delegates to the two-line call. Kept as a dispatch so a caller can
        # compare the reference method against sclet's own without changing
        # their surrounding script.
        warning(
            "method = \"cnv\" runs the CNV + signature call; the iterative ",
            "signature method is method = \"iterative\" (the default).",
            call. = FALSE
        )
        return(CallMalignantCells(
            object = object,
            signature = seed,
            name = name
        ))
    }
    if (!requireNamespace("limma", quietly = TRUE)) {
        stop("Package 'limma' is needed. Please install it.")
    }
    if (!is.numeric(convergence) || length(convergence) != 1 || convergence <= 0 || convergence > 1) {
        stop("`convergence` must be a number in (0, 1].")
    }
    if (!nzchar(name)) {
        stop("`name` must be a non-empty prefix.")
    }

    # ---- seeds ---------------------------------------------------------------
    if (is.list(seed) && all(c("tumor_genes", "normal_genes") %in% names(seed))) {
        tumor_genes <- seed$tumor_genes
        normal_genes <- seed$normal_genes
        if (is.null(n_ctrl) && !is.null(seed$info$signature$n_ctrl)) {
            n_ctrl <- seed$info$signature$n_ctrl
        }
        if (is.null(n_bins) && !is.null(seed$info$signature$n_bins)) {
            n_bins <- seed$info$signature$n_bins
        }
        seed_source <- "FindMalignantSeedGenes"
    } else if (is.list(seed) && length(seed) >= 2 && !is.null(names(seed))) {
        idx <- match(c("tumor", "malignant"), tolower(names(seed)))
        idx <- idx[!is.na(idx)]
        if (length(idx) == 0) {
            stop("A named `seed` list needs a 'tumor' or 'malignant' element.")
        }
        idx2 <- match(c("normal", "nonmalignant", "non_malignant"), tolower(names(seed)))
        idx2 <- idx2[!is.na(idx2)]
        if (length(idx2) == 0) {
            stop("A named `seed` list needs a 'normal' element.")
        }
        tumor_genes <- seed[[idx[[1]]]]
        normal_genes <- seed[[idx2[[1]]]]
        seed_source <- "user-supplied"
    } else {
        stop("`seed` must be the output of FindMalignantSeedGenes() or a named list of two gene vectors.")
    }
    if (is.null(n_ctrl)) n_ctrl <- 100
    if (is.null(n_bins)) n_bins <- 25

    # ---- cells ---------------------------------------------------------------
    if (is.null(cells)) {
        cells <- colnames(object)
    }
    cells <- intersect(as.character(cells), colnames(object))
    if (length(cells) < 2 * min_cells_per_group) {
        stop(
            "Need at least ", 2 * min_cells_per_group, " cells to split; got ",
            length(cells), "."
        )
    }
    sub <- object[, cells, drop = FALSE]

    if (!sample_col %in% colnames(SummarizedExperiment::colData(sub))) {
        stop(
            "colData has no '", sample_col, "'. The pseudobulk test needs a ",
            "sample (patient) column; set one with colData(object)$", sample_col,
            " <- ... ."
        )
    }
    sample <- as.character(SummarizedExperiment::colData(sub)[[sample_col]])

    if (is.null(assay)) {
        assay <- resolve_marker_assay(sub, assay = NULL, layer = NULL)
    }
    if (!assay %in% SummarizedExperiment::assayNames(sub)) {
        stop("Assay '", assay, "' not found.")
    }
    mat <- SummarizedExperiment::assay(sub, assay)

    tumor_genes <- intersect(tumor_genes, rownames(sub))
    normal_genes <- intersect(normal_genes, rownames(sub))
    if (length(tumor_genes) < 2 || length(normal_genes) < 2) {
        stop(
            "Too few seed genes are present in the object (tumour: ", length(tumor_genes),
            ", normal: ", length(normal_genes), "). Check the gene id style."
        )
    }
    seed_genes <- intersect(seed_genes, rownames(sub))

    # ---- gene binning --------------------------------------------------------
    # Scoring must see the WHOLE feature set: AddModuleScore subtracts, for each
    # signature gene, the mean of same-expression-bin control genes. If the
    # expression matrix handed to the scorer held only the signature genes they
    # would be each other's controls and the score would collapse to zero.
    expr_full <- sclet_malignancy_dense(mat)
    expr_use <- expr_full[unique(c(tumor_genes, normal_genes, seed_genes)), , drop = FALSE]
    bins <- sclet_gene_expression_bins(expr_full, n_bins = n_bins)

    # ---- iteration -----------------------------------------------------------
    prev_labels <- NULL
    iteration_log <- list()
    labels <- NULL
    agreement <- NA_real_
    stop_reason <- "max_iter"
    converged <- FALSE
    seen <- character(0)
    seed_history <- list()

    for (iter in seq_len(max_iter)) {
        # Both signatures are excluded from each other's control pool: a
        # signature's own high-expression genes would otherwise land in the top
        # expression bin, use each other as controls, and cancel to zero.
        all_sig <- unique(c(tumor_genes, normal_genes))
        mal_score <- sclet_add_module_score(expr_full, tumor_genes, bins = bins,
                                            n_ctrl = n_ctrl, cells = cells,
                                            exclude = all_sig)
        non_score <- sclet_add_module_score(expr_full, normal_genes, bins = bins,
                                            n_ctrl = n_ctrl, cells = cells,
                                            exclude = all_sig)

        split <- sclet_kmeans_two_scores(mal_score, non_score, nstart = nstart)
        if (is.null(split)) {
            stop_reason <- "scores_constant"
            break
        }

        n_mal <- sum(split$labels == "malignant")
        n_non <- sum(split$labels == "non_malignant")
        if (n_mal < min_cells_per_group || n_non < min_cells_per_group) {
            labels <- split$labels
            stop_reason <- "split_degenerate"
            iteration_log[[iter]] <- sclet_malignancy_log_entry(
                iter, tumor_genes, normal_genes, n_mal, n_non, NA_real_, FALSE
            )
            break
        }

        # Label agreement, allowing either k-means id to be "the match": a
        # flipped assignment is the same partition, so the gain is the number of
        # cells on the correct side rather than the diagonal of the raw table.
        agreement <- NA_real_
        if (!is.null(prev_labels)) {
            same <- sum(split$labels == prev_labels)
            swapped <- sum(split$labels != prev_labels)
            # max(), not the raw diagonal: a flipped k-means assignment is the
            # same partition, and the raw diagonal would read it as 0 agreement.
            agreement <- max(same, swapped) / length(split$labels)
        }
        converged <- !is.na(agreement) && agreement >= convergence

        # Canonical key for the partition, so the SAME split reached again with
        # the opposite k-means ids still matches.
        keys <- sclet_malignancy_split_key(split$cluster)

        labels <- split$labels
        iteration_log[[iter]] <- sclet_malignancy_log_entry(
            iter, tumor_genes, normal_genes, n_mal, n_non, agreement, converged
        )
        seed_history[[iter]] <- list(tumor = tumor_genes, normal = normal_genes)

        if (isTRUE(verbose)) {
            message(sprintf(
                "[%s] round %d: %s signature %d/%d genes, %d malignant / %d non-malignant, agreement %s",
                name, iter, seed_source, length(tumor_genes), length(normal_genes),
                n_mal, n_non,
                if (is.na(agreement)) "n/a" else sprintf("%.3f", agreement)
            ))
        }

        if (converged) {
            stop_reason <- "converged"
            break
        }
        prev_labels <- labels

        # Re-derive signatures from the split taken THIS round. The reference
        # method's re-scores the cells with the previous round's labels, which
        # makes each round lag one split behind.
        deg <- sclet_malignancy_pseudobulk_de(
            mat, labels, sample, n_top = n_top, padj = padj, lfc = lfc,
            min_patients = min_patients
        )
        if (is.null(deg) || length(deg$malignant) == 0 || length(deg$non_malignant) == 0) {
            stop_reason <- "signature_empty"
            if (verbose) message("  no gene passes the cutoffs on the split; stopping")
            break
        }
        tumor_genes <- deg$malignant
        normal_genes <- deg$non_malignant

        # A split configuration revisited twice is a genuine cycle, not slow
        # convergence: the next round would reproduce the same split.
        if (keys %in% seen) {
            stop_reason <- "cycle"
            if (verbose) message("  the same split was produced before; stopping")
            break
        }
        seen <- c(seen, keys)
    }

    if (is.null(labels)) {
        stop(
            "The classifier produced no labels (reason: ", stop_reason,
            "). Check the expression values and the seed genes."
        )
    }

    # ---- diagnostics ---------------------------------------------------------
    n_mal <- sum(labels == "malignant")
    n_non <- sum(labels == "non_malignant")
    malignant_fraction <- n_mal / length(labels)

    # Side-of-origin: how much of each side comes from each sample. The
    # non-malignant side holding many tumour-sample cells is the "control pool
    # changed" situation, which is expected by design and not an error.
    by_sample <- table(sample, factor(labels, levels = c("malignant", "non_malignant")))
    side_of_origin <- list(
        table = by_sample,
        non_malignant_from_samples = sort(table(sample[labels == "non_malignant"]), decreasing = TRUE),
        malignant_from_samples = sort(table(sample[labels == "malignant"]), decreasing = TRUE)
    )

    # Seed tracking: a known metaplasia / intestinal profile expressed on the
    # malignant side means the split may be an identity axis, not malignancy.
    seed_tracking <- NULL
    if (length(seed_genes) > 0) {
        sg_score <- sclet_add_module_score(expr_full, seed_genes, bins = bins, n_ctrl = n_ctrl,
                                          cells = cells, exclude = seed_genes)
        seed_tracking <- list(
            genes = seed_genes,
            mean_malignant = mean(sg_score[labels == "malignant"]),
            mean_non_malignant = mean(sg_score[labels == "non_malignant"]),
            difference = mean(sg_score[labels == "malignant"]) - mean(sg_score[labels == "non_malignant"])
        )
    }

    # Gene stability: genes that recur in the final signatures.
    gene_frequency <- sort(
        table(unlist(lapply(seed_history, unlist))),
        decreasing = TRUE
    )

    warnings_out <- character(0)
    if (!isTRUE(converged)) {
        warnings_out <- c(warnings_out, sprintf(
            "The labels had not reached %.2f agreement after %d round(s) (stopped: %s); the split is not stable and the labels should not be read as a final classification.",
            convergence, length(iteration_log), stop_reason
        ))
    }
    # k-means splits ANY data into two parts, including a population with no
    # real two-way structure. A lopsided split, or one where a side is too thin
    # to re-derive a signature from, is the signature of that failure mode.
    degenerate <- identical(stop_reason, "split_degenerate")
    lopsided <- malignant_fraction < 0.05 || malignant_fraction > 0.95 ||
        identical(stop_reason, "scores_constant")
    if (lopsided || degenerate) {
        warnings_out <- c(warnings_out, sprintf(
            "k-means split the cells %.1f%% / %.1f%%: a two-way split of a population with no clear two-way structure produces this kind of division, and it is not a meaningful grouping at this size.",
            100 * malignant_fraction, 100 * (1 - malignant_fraction)
        ))
    }
    if (!is.null(seed_tracking) && abs(seed_tracking$difference) > 0.1) {
        warnings_out <- c(warnings_out, sprintf(
            "The supplied `seed_genes` score differ on the malignant side (difference %.3f). A signature that separates the two sides this sharply may mean the split tracks that programme rather than malignancy; compare with an independent malignant-cell signature.",
            seed_tracking$difference
        ))
    }
    warnings_out <- c(warnings_out,
        "The label means 'the side of the largest transcriptional split'. Converged labels are consistent, not validated; malignant identity needs an independent line of evidence (CNV inference, histology, or a signature from another dataset)."
    )

    result <- list(
        labels = stats::setNames(factor(labels, levels = c("malignant", "non_malignant")), cells),
        scores = data.frame(
            cell = cells,
            malignant_score = as.numeric(mal_score),
            non_malignant_score = as.numeric(non_score),
            stringsAsFactors = FALSE
        ),
        iteration_log = iteration_log,
        iterations = length(iteration_log),
        converged = isTRUE(converged),
        agreement = agreement,
        stop_reason = stop_reason,
        final_signature = list(malignant = tumor_genes, non_malignant = normal_genes),
        seed_history = seed_history,
        gene_frequency = gene_frequency,
        diagnostics = list(
            n_cells = length(labels),
            n_malignant = n_mal,
            n_non_malignant = n_non,
            malignant_fraction = malignant_fraction,
            side_of_origin = side_of_origin,
            seed_tracking = seed_tracking
        ),
        params = list(
            cells = length(cells), sample_col = sample_col, assay = assay,
            n_top = n_top, padj = padj, lfc = lfc, convergence = convergence,
            max_iter = max_iter, min_cells_per_group = min_cells_per_group,
            min_patients = min_patients, n_ctrl = n_ctrl, n_bins = n_bins,
            nstart = nstart, name = name, seed_source = seed_source
        ),
        warnings = warnings_out
    )

    # ---- write back ----------------------------------------------------------
    cdata <- SummarizedExperiment::colData(object)
    cols <- list(
        malignant_score = rep(NA_real_, nrow(cdata)),
        nonmalignant_score = rep(NA_real_, nrow(cdata)),
        cluster = rep(NA_character_, nrow(cdata)),
        classification = factor(rep(NA_character_, nrow(cdata)),
                                levels = c("malignant", "non_malignant"))
    )
    pos <- match(cells, rownames(cdata))
    cols$malignant_score[pos] <- result$scores$malignant_score
    cols$nonmalignant_score[pos] <- result$scores$non_malignant_score
    cols$cluster[pos] <- as.character(labels)
    cols$classification[pos] <- as.character(labels)
    for (nm in names(cols)) {
        cdata[[paste0(name, "_", nm)]] <- cols[[nm]]
    }
    SummarizedExperiment::colData(object) <- cdata

    object <- sclet_set_analysis(object, "malignancy", result)
    md <- S4Vectors::metadata(object)
    md$malignancy <- if (is.null(md$malignancy)) list() else md$malignancy
    md$malignancy[[name]] <- result
    S4Vectors::metadata(object) <- md

    object <- sclet_set_state_record(
        object = object,
        type = "annotation",
        id = name,
        value = list(
            method = "iterative AddModuleScore + k-means + pseudobulk limma (Zhang et al., Gut 2020)",
            inputs = list(
                assay = assay,
                cells = length(cells),
                sample_col = sample_col,
                seed = seed_source
            ),
            artifacts = list(
                columns = paste0(name, c("_malignant_score", "_nonmalignant_score",
                                         "_cluster", "_classification")),
                metadata = paste0("malignancy$", name)
            ),
            params = result$params,
            summary = list(
                n_malignant = n_mal,
                n_non_malignant = n_non,
                malignant_fraction = malignant_fraction,
                iterations = length(iteration_log),
                converged = isTRUE(converged),
                agreement = agreement,
                stop_reason = stop_reason
            ),
            created_at = Sys.time()
        )
    )
    object <- sclet_log_command(
        object,
        "ClassifyMalignantEpithelium",
        params = result$params,
        outputs = list(
            columns = paste0(name, c("_malignant_score", "_nonmalignant_score",
                                     "_cluster", "_classification")),
            analysis = "malignancy"
        )
    )

    for (w in unique(warnings_out)) {
        warning(w, call. = FALSE)
    }
    object
}

# ---------------------------------------------------------------------------
# internal helpers
# ---------------------------------------------------------------------------

sclet_bulk_log_offset <- function(assay_note) {
    note <- tolower(assay_note)
    if (!grepl("tpm|fpkm", note)) {
        return(NULL)
    }
    m <- regmatches(note, regexpr("\\+\\s*[0-9.]+", note))
    if (length(m) == 0) {
        return(0)
    }
    as.numeric(sub("\\+", "", m))
}

sclet_malignancy_dense <- function(x) {
    if (inherits(x, c("DelayedArray", "Matrix"))) {
        x <- as.matrix(x)
    }
    if (is.null(dim(x))) {
        x <- matrix(x, nrow = 1)
    }
    x
}

# Expression bins over the whole feature set, as AddModuleScore bins control
# genes by mean expression.
sclet_gene_expression_bins <- function(mat, n_bins = 25) {
    means <- if (inherits(mat, "DelayedMatrix")) {
        DelayedMatrixStats::rowMeans2(mat)
    } else {
        MatrixGenerics::rowMeans2(sclet_malignancy_dense(mat))
    }
    names(means) <- rownames(mat)
    if (length(unique(means)) < 2) {
        return(stats::setNames(rep(1L, length(means)), names(means)))
    }
    b <- as.integer(cut(means, breaks = max(1L, min(n_bins, length(unique(means)))),
                        labels = FALSE))
    b[is.na(b)] <- 1L
    stats::setNames(b, names(means))
}

# AddModuleScore-compatible scoring: for each signature gene, subtract the mean
# expression of `n_ctrl` genes drawn from the same expression bin, then average
# over signature genes. Signature genes absent from `expr` are dropped; empty
# signatures score 0.
sclet_add_module_score <- function(expr, genes, bins = NULL, n_ctrl = 100, cells = NULL,
                                    exclude = NULL) {
    if (is.null(dim(expr))) {
        expr <- matrix(expr, nrow = 1)
    }
    ids <- colnames(expr)
    if (is.null(ids) && !is.null(cells) && length(cells) == ncol(expr)) {
        ids <- cells
    }
    genes <- intersect(unique(genes), rownames(expr))
    if (length(genes) == 0) {
        return(stats::setNames(rep(0, ncol(expr)), ids))
    }
    if (is.null(bins)) {
        means <- rowMeans(expr)
        bins <- sclet_gene_expression_bins(
            matrix(means, ncol = 1, dimnames = list(names(means), "x"))
        )
    }
    # Row indices are resolved by NAME against `bins`, which covers the full
    # feature set; indexing it with the submatrix's row numbers would silently
    # map to the wrong genes (or go out of bounds).
    bin_of <- bins[rownames(expr)]
    bin_of[is.na(bin_of)] <- NA_integer_
    scores <- rep(0, ncol(expr))
    row_idx <- seq_len(nrow(expr))
    for (g in genes) {
        b <- bin_of[[g]]
        same <- !is.na(bin_of) & bin_of == b
        if (is.na(b)) same <- rep(FALSE, length(bin_of))
        same[rownames(expr) == g] <- FALSE
        if (length(exclude) > 0) {
            same[rownames(expr) %in% exclude] <- FALSE
        }
        pool <- row_idx[same]
        if (length(pool) == 0) {
            pool <- row_idx[!rownames(expr) %in% c(g, exclude)]
        }
        ctrl <- if (length(pool) > n_ctrl) sample(pool, n_ctrl) else pool
        if (length(ctrl) == 0) {
            next
        }
        if (is.matrix(expr)) {
            ctrl_mean <- colMeans(expr[ctrl, , drop = FALSE])
            scores <- scores + expr[g, ] - ctrl_mean
        } else {
            ctrl_mean <- Matrix::colMeans(expr[ctrl, , drop = FALSE])
            scores <- scores + as.numeric(expr[g, ]) - as.numeric(ctrl_mean)
        }
    }
    scores <- scores / length(genes)
    stats::setNames(as.numeric(scores), ids)
}

# k-means on the two scaled scores, k = 2. The cluster whose cells have the
# higher malignant score is "malignant". Returns NULL when the scores are
# constant (k-means would split identical points arbitrarily).
sclet_kmeans_two_scores <- function(mal_score, non_score, nstart = 25) {
    x <- cbind(malignant = as.numeric(mal_score), non = as.numeric(non_score))
    if (any(!is.finite(x))) {
        stop("Non-finite scores encountered; check for NA expression values in the scored genes.")
    }
    if (all(x == x[1, ])) {
        return(NULL)
    }
    sd_x <- apply(x, 2, stats::sd)
    sd_x[sd_x == 0] <- 1
    z <- sweep(x, 2, sd_x, "/")
    km <- stats::kmeans(z, centers = 2, nstart = nstart)
    mal_cluster <- which.max(tapply(x[, "malignant"], km$cluster, mean))
    list(
        cluster = km$cluster,
        labels = ifelse(km$cluster == mal_cluster, "malignant", "non_malignant"),
        sizes = tabulate(km$cluster)
    )
}

# Canonical one-string key for a two-way split, independent of which k-means
# cluster id was called "malignant".
sclet_malignancy_split_key <- function(cluster) {
    ord <- order(cluster)
    n_a <- sum(cluster == cluster[[1]])
    group_a <- if (n_a <= length(cluster) / 2) ord[seq_len(n_a)] else ord[-seq_len(n_a)]
    if (length(group_a) > length(cluster) / 2) {
        group_a <- setdiff(seq_along(cluster), group_a)
    }
    paste(group_a, collapse = ",")
}

sclet_malignancy_log_entry <- function(iter, tumor_genes, normal_genes, n_mal, n_non,
                                       agreement, converged) {
    list(
        iteration = iter,
        n_malignant_genes = length(tumor_genes),
        n_non_malignant_genes = length(normal_genes),
        n_malignant = n_mal,
        n_non_malignant = n_non,
        agreement = agreement,
        converged = converged
    )
}

# Pseudobulk differential expression between the two sides: mean expression per
# (sample, side), then limma with a sample-blocked design (paired by donor) or a
# one-way design when too few donors have cells on both sides.
# Returns NULL when the test cannot run or nothing passes the cutoffs.
sclet_malignancy_pseudobulk_de <- function(mat, labels, sample, n_top = 50, padj = 0.01,
                                           lfc = NULL, min_patients = 3) {
    side <- as.character(labels)
    samples <- unique(sample)
    pb_side <- character(0)
    pb_patient <- character(0)
    pb_ids <- character(0)
    pb_cells <- list()
    for (s in samples) {
        for (k in c("malignant", "non_malignant")) {
            idx <- which(sample == s & side == k)
            if (length(idx) == 0) {
                next
            }
            pb_ids <- c(pb_ids, sprintf("%s__%s", s, k))
            pb_side <- c(pb_side, k)
            pb_patient <- c(pb_patient, s)
            pb_cells[[length(pb_cells) + 1L]] <- idx
        }
    }
    # Aggregate by column-subset means, never transposing the full matrix.
    pb <- matrix(0, nrow = nrow(mat), ncol = length(pb_ids),
                 dimnames = list(rownames(mat), pb_ids))
    for (j in seq_along(pb_ids)) {
        pb[, j] <- MatrixGenerics::rowMeans(mat[, pb_cells[[j]], drop = FALSE])
    }

    cond <- factor(pb_side, levels = c("non_malignant", "malignant"))
    pat <- factor(pb_patient)
    paired <- sum(table(pb_patient, pb_side)[, "malignant"] > 0 &
                      table(pb_patient, pb_side)[, "non_malignant"] > 0) >= min_patients
    design <- if (paired) {
        stats::model.matrix(~ pat + cond)
    } else {
        stats::model.matrix(~ cond)
    }
    if (qr(design)$rank < ncol(design)) {
        design <- stats::model.matrix(~ cond)
        paired <- FALSE
    }
    if (ncol(pb) < 4 || qr(design)$rank < 2) {
        return(NULL)
    }
    fit <- limma::lmFit(pb, design)
    fit <- limma::eBayes(fit)
    tab <- limma::topTable(fit, coef = ncol(design), number = Inf, sort.by = "P")
    tab$gene <- rownames(tab)

    ok <- !is.na(tab$adj.P.Val) & tab$adj.P.Val < padj
    if (!is.null(lfc)) {
        ok <- ok & abs(tab$logFC) >= lfc
    }
    hygiene <- grepl("^MT-|^MTRNR|^RPL[0-9]|^RPS[0-9]|^MRPL|^MRPS|^MALAT1$|^HBA[12]$|^HBB$",
                     tab$gene)
    sig <- tab[ok & !hygiene, , drop = FALSE]
    list(
        malignant = utils::head(sig$gene[sig$logFC > 0], n_top),
        non_malignant = utils::head(sig$gene[sig$logFC < 0], n_top),
        deg = sig,
        paired = paired,
        n_pseudobulk_columns = ncol(pb)
    )
}
