# Copy-number inference as an independent line of evidence for malignant cells
#
# The iterative classification ported in R/malignancy.R decides malignancy from
# expression alone: it finds the largest transcriptional split and calls one side
# "malignant". Issue YuLab-SMU/ProjectYulab#26 documents, correctly, that on
# tumour tissue that split is often metaplastic or intestinal identity rather
# than malignancy, and that the non-malignant reference stops being the normal
# samples after the first round.
#
# The defining feature of a malignant cell is genomic instability, which is
# orthogonal to any transcriptional identity programme and does not depend on a
# self-referential reference. This layer infers it with inferCNV and returns a
# per-cell CNV burden that a caller can combine with expression evidence
# without either one being able to fake the other.
#
# Design points that matter:
#
#   * The reference is EXOGENOUS. Reference cells are normal cells or normal
#     samples named by the caller. They never come from the iterative split, so
#     the evidence does not drift when the split moves.
#   * inferCNV is a Suggests dependency and needs a JAGS system library. When it
#     is unavailable this layer returns status = "unavailable" with the reason
#     rather than quietly degrading to "no CNV found" -- absence of a
#     measurement must not read as a negative result.
#   * Gene order needs chromosome and position. It is taken from rowData when
#     present, from a caller-supplied file otherwise, and only then from an
#     optional TxDb. Which source was used is recorded.

# A reference from another lineage can look like copy-number change.
#
# Measured on GSE144735 (2026-10-19): an immune/stromal reference gave CNV
# burden AUROC 0.709 for tumour vs normal epithelium and 0.6% tumour cells above
# the q99 reference cut; using normal epithelium as the reference raised AUROC to
# 0.773 and q99 sensitivity to 13.3% (4.7% held-out normal false positives). The
# difference is lineage expression, not copy number: reference-centring against
# immune cells leaves epithelial identity itself as a genomic pattern. The
# package cannot choose the reference for the caller, but it can say when the
# labels suggest this situation.
sclet_cnv_lineage_note <- function(object, ref) {
    cd <- SummarizedExperiment::colData(object)
    cols <- grep("^cell[._]?type|celltype|^annotation$|^cell_type", colnames(cd),
        ignore.case = TRUE, value = TRUE)
    if (!length(cols)) return(NULL)
    is_ref <- colnames(object) %in% ref
    for (col in cols) {
        labels <- as.character(cd[[col]])
        ref_labels <- unique(labels[is_ref])
        obs_labels <- unique(labels[!is_ref])
        ref_labels <- ref_labels[!is.na(ref_labels) & nzchar(ref_labels)]
        obs_labels <- obs_labels[!is.na(obs_labels) & nzchar(obs_labels)]
        if (length(ref_labels) && length(obs_labels) &&
            !length(intersect(ref_labels, obs_labels))) {
            return(sprintf(paste0(
                "Reference cells and observation cells share no '%s' label ",
                "(reference: %s; observation: %s). A reference from another ",
                "lineage can be measured as copy-number burden, because ",
                "lineage expression survives reference-centring; prefer normal ",
                "cells of the same lineage when such cells exist."),
                col,
                paste(utils::head(ref_labels, 3), collapse = "/"),
                paste(utils::head(obs_labels, 3), collapse = "/")))
        }
    }
    NULL
}

#' Readiness check for CNV-based malignant-cell inference
#'
#' Reports whether a copy-number measurement can run on this object: inferCNV
#' and its JAGS dependency present, gene positions resolvable, and reference
#' cells available. Read-only; runs nothing.
#'
#' @param object A `SingleCellExperiment` object.
#' @param reference_cells Character vector of cell names to use as the normal
#'   reference. Required unless `reference_col` names a `colData` column whose
#'   values mark reference cells.
#' @param reference_col Optional `colData` column; its non-`FALSE` values select
#'   reference cells.
#' @param gene_order Optional data frame or path with `gene`, `chr`, `start`,
#'   `end`.
#' @param species `"human"` or `"mouse"`, used only to name the TxDb fallback.
#'
#' @return A list with `ready` (logical), `status`, `issues` (character),
#'   `n_reference_cells`, `n_genes_with_position`, and `gene_order_source`.
#'
#' @seealso [RunCNVInference()]
#' @export
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' readiness <- check_cnv_readiness(sce, reference_cells = colnames(sce)[1:5])
#' readiness$ready
#' readiness$issues
check_cnv_readiness <- function(object,
                                reference_cells = NULL,
                                reference_col = NULL,
                                gene_order = NULL,
                                species = c("human", "mouse")) {
    species <- match.arg(species)
    issues <- character(0)

    has_infercnv <- requireNamespace("infercnv", quietly = TRUE)
    has_rjags <- requireNamespace("rjags", quietly = TRUE)
    if (!has_infercnv) {
        issues <- c(issues, paste(
            "Package 'infercnv' is not installed. Install it with",
            "BiocManager::install(\"infercnv\"); it also needs the JAGS system",
            "library (apt: jags / brew: jags)."
        ))
    } else if (!has_rjags) {
        issues <- c(issues, paste(
            "Package 'infercnv' is installed but 'rjags' is not loadable.",
            "infercnv imports rjags at load time and needs JAGS 4.x on the",
            "system (apt: jags)."
        ))
    }

    ref <- sclet_cnv_resolve_reference(object, reference_cells = reference_cells,
                                       reference_col = reference_col)
    lineage_note <- if (length(ref)) sclet_cnv_lineage_note(object, ref) else NULL
    if (length(ref) == 0) {
        issues <- c(issues, paste(
            "No reference cells. Name normal cells or normal samples via",
            "`reference_cells` or `reference_col`; CNV inference is only",
            "meaningful relative to a normal reference."
        ))
    }

    # read.table rather than the helper: a gene order with NO overlap must be
    # REPORTED here, and the helper stops on that condition by design.
    order <- tryCatch(
        sclet_cnv_gene_order(object, gene_order = gene_order, species = species),
        error = function(e) e
    )
    if (inherits(order, "error")) {
        issues <- c(issues, conditionMessage(order))
        order <- sclet_cnv_empty_order()
    }
    n_chr <- if (nrow(order) > 0) length(unique(order$chr)) else 0L
    if (length(issues) == 0 && nrow(order) < 100) {
        issues <- c(issues, sprintf(
            "Only %d genes got a chromosome position. inferCNV needs a genome-wide ordered gene list; supply `gene_order`.",
            nrow(order)
        ))
    }
    # inferCNV builds a sliding window ALONG each chromosome, so a gene order
    # that piles every gene onto one chromosome produces a method that cannot
    # see a copy-number profile at all. A gene-count gate alone passed such an
    # input, which is why it is checked separately here.
    if (length(issues) == 0 && n_chr < 2) {
        issues <- c(issues, sprintf(
            "The gene order covers %d chromosome(s) (%d genes). inferCNV slides a window along each chromosome, so a single chromosome cannot carry a copy-number profile; supply a genome-wide ordered gene list.",
            n_chr, nrow(order)
        ))
    }
    # Merging duplicated rows is routine (GTF/Ensembl export one row per
    # transcript), so it is reported as a note and never blocks the run.
    notes <- character(0)
    if (!is.null(lineage_note)) notes <- c(notes, lineage_note)
    if (!is.null(attr(order, "n_duplicate_genes")) &&
        attr(order, "n_duplicate_genes") > 0) {
        notes <- c(notes, sprintf(
            "%d duplicated gene rows in `gene_order` were merged into one interval per gene.",
            attr(order, "n_duplicate_genes")
        ))
    }

    list(
        ready = length(issues) == 0,
        status = if (length(issues) == 0) "ready" else "not_ready",
        issues = issues,
        n_reference_cells = length(ref),
        n_genes_with_position = nrow(order),
        n_chromosomes = n_chr,
        notes = notes,
        gene_order_source = attr(order, "source"),
        infercnv_version = if (has_infercnv) as.character(utils::packageVersion("infercnv")) else NA_character_
    )
}

#' Infer copy-number burden per cell with inferCNV
#'
#' Runs inferCNV against caller-supplied reference cells and returns a per-cell
#' CNV burden that [CallMalignantCells()] can combine with expression evidence.
#'
#' @param object A `SingleCellExperiment` object.
#' @param reference_cells Character vector of cells to use as the normal
#'   reference. Mutually exclusive with `reference_col`.
#' @param reference_col A `colData` column: rows with a non-`FALSE`, non-`NA`
#'   value are reference cells. Mutually exclusive with `reference_cells`.
#' @param gene_order Optional data frame or path with columns `gene`, `chr`,
#'   `start`, `end`. Required when `rowData` carries no positions.
#' @param species `"human"` or `"mouse"`; names the TxDb fallback.
#' @param assay Assay to use. Defaults to `"counts"` because inferCNV models
#'    raw counts. A non-count assay (non-integer, or with negative values) is
#'    refused rather than truncated, and the refusal happens BEFORE the
#'    readiness gate so that a malformed matrix is never masked by an
#'    unrelated prerequisite message.
#' @param hmm Logical. Run inferCNV's HMM post-processing so that per-cell CNV
#'   state probabilities are reported. Defaults to `TRUE`.
#' @param cutoff,window_length,cluster_by_groups,analysis_mode,k_obs_groups
#'   Passed to [infercnv::run()].
#' @param num_threads Threads for inferCNV.
#' @param out_dir Directory for inferCNV's intermediate files. A temporary
#'   directory is used and removed unless `keep_files = TRUE`.
#' @param keep_files Logical. Keep the inferCNV output directory.
#' @param name Analysis record id. Defaults to `"cnv"`.
#'
#'
#' @return `object` with `colData` columns `<name>_burden`, `<name>_n_cnv_regions`
#'   and `<name>_is_reference`, plus a state record of type `annotation` under
#'   `name`. When the measurement cannot run, the columns are filled with `NA`
#'   and the record carries `status = "unavailable"` with the reason; the
#'   function does not error and does not write a zero, so a downstream call can
#'   tell "not measured" from "measured negative".
#'
#' @section What this measures:
#' CNV burden is the mean absolute deviation of a cell's smoothed, reference-
#' centred expression across the genome. Large deletions or amplifications shared
#' by a clone raise it. It is evidence of genomic instability, which is why it is
#' treated as an independent line from any transcriptional signature.
#'
#' The limit is resolution: sub-chromosomal events below inferCNV's window are
#' not seen, low-depth cells carry noise that masks real events, and a
#' hypermutated but copy-number-quiet tumour can score low while still being
#' malignant. CNV-positivity is therefore strong positive evidence, and
#' CNV-negativity is weak negative evidence.
#'
#' @seealso [check_cnv_readiness()], [CallMalignantCells()]
#' @export
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' \dontrun{
#' # needs the infercnv package and the JAGS system library
#' sce <- RunCNVInference(sce, reference_cells = colnames(sce)[1:5])
#' }
RunCNVInference <- function(object,
                            reference_cells = NULL,
                            reference_col = NULL,
                            gene_order = NULL,
                            species = c("human", "mouse"),
                            assay = "counts",
                            hmm = TRUE,
                            cutoff = 1,
                            window_length = 101,
                            cluster_by_groups = TRUE,
                            analysis_mode = "subclusters",
                            k_obs_groups = 1,
                            num_threads = 1,
                            out_dir = NULL,
                            keep_files = FALSE,
                            name = "cnv") {
    species <- match.arg(species)
    if (!nzchar(name)) {
        stop("`name` must be a non-empty prefix.")
    }
    # Validate the matrix BEFORE the readiness gate. Readiness now refuses a
    # single-chromosome gene order, which is correct, but that must not mask a
    # malformed assay -- a more basic error the caller needs to hear first.
    if (!assay %in% SummarizedExperiment::assayNames(object)) {
        stop("Assay '", assay, "' not found; inferCNV models raw counts.")
    }
    assay_matrix <- SummarizedExperiment::assay(object, assay)
    # A sparse dgCMatrix is numeric counts, but `is.numeric()` is FALSE for an
    # S4 Matrix -- so every ordinary SingleCellExperiment was refused with
    # "not numeric" on real data (2026-10-18).
    if (!is.numeric(assay_matrix) && !inherits(assay_matrix, c("Matrix", "DelayedArray"))) {
        stop("Assay '", assay, "' is not numeric.")
    }
    # A non-count matrix is an input error, not a missing prerequisite: report it
    # rather than letting the readiness gate return `unavailable`, which would
    # bury the real problem under a gene-order or reference message.
    probe <- sclet_malignancy_dense(assay_matrix)
    if (any(probe < 0, na.rm = TRUE)) {
        stop(
            "Assay '", assay, "' has negative values; inferCNV needs raw counts.",
            if (assay %in% c("logcounts", "scaled")) {
                " This looks like a transformed assay -- use assay = \"counts\"."
            } else ""
        )
    }
    if (any(probe != round(probe), na.rm = TRUE)) {
        stop(
            "Assay '", assay, "' holds non-integer values; inferCNV needs raw counts. ",
            "Coercing would truncate every value silently. Use assay = \"counts\", ",
            "or round explicitly if you are certain the values are counts."
        )
    }

    # The CNV cutoff is estimated from the reference population, so a handful of
    # reference cells cannot support it and the line stays unavailable. The call
    # then degrades to signature-only, which is correct but silent -- so say so
    # before the readiness gate, which may return early for another reason.
    n_ref <- length(sclet_cnv_resolve_reference(object,
                                                 reference_cells = reference_cells,
                                                 reference_col = reference_col))
    if (n_ref > 0 && n_ref < 5) {
        warning(sprintf(
            "Only %d reference cell(s) supplied. The CNV cutoff is estimated from the reference population, so with fewer than 5 the CNV line will be unavailable and the call will degrade to signature-only. Use more normal cells.",
            n_ref
        ), call. = FALSE)
    }

    readiness <- check_cnv_readiness(
        object,
        reference_cells = reference_cells,
        reference_col = reference_col,
        gene_order = gene_order,
        species = species
    )
    if (!isTRUE(readiness$ready)) {
        return(sclet_cnv_record_unavailable(object, name, readiness))
    }

    ref <- sclet_cnv_resolve_reference(object, reference_cells = reference_cells,
                                       reference_col = reference_col)
    counts <- SummarizedExperiment::assay(object, assay)
    order <- sclet_cnv_gene_order(object, gene_order = gene_order, species = species)

    common <- intersect(rownames(object), order$gene)
    if (length(common) < 100) {
        stop("Fewer than 100 genes have both expression and a chromosome position.")
    }

    own_dir <- is.null(out_dir)
    if (own_dir) {
        out_dir <- tempfile("sclet-infercnv-")
    }
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    if (!isTRUE(keep_files)) {
        on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)
    }

    cnv_counts <- sclet_malignancy_dense(counts[common, , drop = FALSE])
    mode(cnv_counts) <- "integer"

    order_file <- file.path(out_dir, "gene_order.tsv")
    utils::write.table(order[order$gene %in% common, c("gene", "chr", "start", "end")],
                       file = order_file, sep = "\t", quote = FALSE,
                       row.names = FALSE, col.names = FALSE)

    annot <- data.frame(
        cell = colnames(cnv_counts),
        group = ifelse(colnames(cnv_counts) %in% ref, "reference", "observation"),
        stringsAsFactors = FALSE
    )
    annot_file <- file.path(out_dir, "annotations.tsv")
    utils::write.table(annot, file = annot_file, sep = "\t", quote = FALSE,
                       row.names = FALSE, col.names = FALSE)

    cnv_obj <- infercnv::CreateInfercnvObject(
        raw_counts_matrix = cnv_counts,
        annotations_file = annot_file,
        gene_order_file = order_file,
        ref_group_names = "reference"
    )

    cnv_obj <- infercnv::run(
        cnv_obj,
        cutoff = cutoff,
        out_dir = out_dir,
        window_length = window_length,
        cluster_by_groups = cluster_by_groups,
        analysis_mode = analysis_mode,
        k_obs_groups = k_obs_groups,
        HMM = isTRUE(hmm),
        HMM_report_by = "cell",
        num_threads = num_threads,
        no_plot = TRUE,
        save_rds = FALSE,
        save_final_rds = FALSE,
        write_expr_matrix = FALSE,
        resume_mode = FALSE
    )

    burden <- sclet_cnv_burden(cnv_obj, cells = colnames(object), ref = ref)

    # inferCNV removes cells below its minimum-count threshold, so the returned
    # matrix routinely covers fewer cells than were supplied. Those cells would
    # land in colData as unexplained NA; say how many and why instead.
    dropped <- setdiff(colnames(object), burden$cell)
    if (length(dropped) > 0) {
        message(sprintf(
            "inferCNV kept %d of %d cells; %d were dropped by its own filters (e.g. min_max_counts_per_cell) and have NA burden. Supply fewer cells or lower inferCNV's cutoff if that is not expected.",
            length(burden$cell), ncol(object), length(dropped)
        ))
    }

    cdata <- SummarizedExperiment::colData(object)
    cdata[[paste0(name, "_burden")]] <- burden$burden[match(rownames(cdata), burden$cell)]
    cdata[[paste0(name, "_is_reference")]] <- rownames(cdata) %in% ref
    SummarizedExperiment::colData(object) <- cdata

    result <- list(
        status = "completed",
        burden = burden,
        n_cells_submitted = ncol(object),
        n_cells_retained = length(burden$cell),
        dropped_cells = dropped,
        reference_cells = ref,
        n_reference_cells = length(ref),
        n_genes = length(common),
        gene_order_source = readiness$gene_order_source,
        reference_notes = readiness$notes %||% character(),
        params = list(
            assay = assay, hmm = isTRUE(hmm), cutoff = cutoff,
            window_length = window_length, cluster_by_groups = cluster_by_groups,
            analysis_mode = analysis_mode, k_obs_groups = k_obs_groups,
            num_threads = num_threads
        ),
        infercnv_version = readiness$infercnv_version
    )

    # `analyses[["malignant_cnv"]]` is a MAP keyed by run id, mirroring
    # `metadata$malignant_cnv[[name]]`, so two runs do not overwrite each other
    # and the reader can address one by name.
    existing <- sclet_get_analysis(object, "malignant_cnv")
    if (is.null(existing)) existing <- list()
    existing[[name]] <- result
    object <- sclet_set_analysis(object, "malignant_cnv", existing)
    md <- S4Vectors::metadata(object)
    md$malignant_cnv <- if (is.null(md$malignant_cnv)) list() else md$malignant_cnv
    md$malignant_cnv[[name]] <- result
    S4Vectors::metadata(object) <- md

    # Its own `cnv` record, not an `annotation` one: the annotation bucket has
    # a single "active" slot, so this made `get_annotation()` return a CNV
    # burden run instead of the object's annotation (measured 2026-10-22). The
    # payload is unchanged -- only the bucket.
    object <- sclet_set_state_record(
        object = object,
        type = "cnv",
        id = name,
        value = list(
            method = "inferCNV per-cell copy-number burden",
            inputs = list(
                assay = assay,
                n_reference_cells = length(ref),
                gene_order_source = readiness$gene_order_source
            ),
            artifacts = list(
                columns = paste0(name, c("_burden", "_is_reference")),
                metadata = paste0("malignant_cnv$", name)
            ),
            params = result$params,
            summary = list(
                status = "completed",
                n_cells = ncol(object),
                n_cells_retained = length(burden$cell),
                n_cells_dropped = length(dropped),
                n_reference_cells = length(ref),
                burden_q = unname(stats::quantile(burden$burden, c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE))
            ),
            created_at = Sys.time()
        )
    )
    object <- sclet_log_command(
        object,
        "RunCNVInference",
        params = result$params,
        outputs = list(metadata = paste0("malignant_cnv$", name))
    )
    object
}

# ---------------------------------------------------------------------------
# inference internals
# ---------------------------------------------------------------------------

sclet_cnv_resolve_reference <- function(object, reference_cells = NULL, reference_col = NULL) {
    if (!is.null(reference_cells) && !is.null(reference_col)) {
        stop("Supply only one of `reference_cells` and `reference_col`.")
    }
    if (!is.null(reference_cells)) {
        ref <- intersect(as.character(reference_cells), colnames(object))
        return(ref)
    }
    if (!is.null(reference_col)) {
        if (!reference_col %in% colnames(SummarizedExperiment::colData(object))) {
            stop("colData has no '", reference_col, "' column.")
        }
        v <- SummarizedExperiment::colData(object)[[reference_col]]
        keep <- !is.na(v) & !(is.logical(v) & !v)
        return(colnames(object)[keep])
    }
    # No explicit reference: fall back to a colData column literally called
    # "reference" when one exists, but never guess from clustering.
    cdata <- SummarizedExperiment::colData(object)
    if ("reference" %in% colnames(cdata)) {
        v <- cdata[["reference"]]
        return(colnames(object)[!is.na(v) & !(is.logical(v) & !v)])
    }
    character(0)
}

sclet_cnv_empty_order <- function() {
    empty <- data.frame(gene = character(0), chr = character(0),
                        start = numeric(0), end = numeric(0), stringsAsFactors = FALSE)
    attr(empty, "source") <- "none"
    empty
}

sclet_cnv_gene_order <- function(object, gene_order = NULL, species = "human") {
    # 1. Caller-supplied data frame or file.
    if (!is.null(gene_order)) {
        ord <- if (is.character(gene_order) && length(gene_order) == 1) {
            utils::read.table(gene_order, sep = "\t", header = grepl("gene", readLines(gene_order, n = 1)),
                              stringsAsFactors = FALSE)
        } else if (is.data.frame(gene_order)) {
            gene_order
        } else {
            stop("`gene_order` must be a file path or a data frame.")
        }
        if (!all(c("gene", "chr", "start", "end") %in% colnames(ord))) {
            stop("`gene_order` needs columns `gene`, `chr`, `start`, `end`.")
        }
        # Normalise BEFORE subsetting: assigning a column to a zero-row data frame
        # raises `replacement has 1 row, data has 0`, which used to surface here
        # whenever the gene ids did not overlap the object at all (ENSG vs symbol
        # being the usual cause). An empty result is a reportable condition, not
        # an internal error.
        ord$chr <- sub("^chr", "", as.character(ord$chr))
        ord$chr_num <- suppressWarnings(as.integer(ord$chr))
        ord$start <- suppressWarnings(as.numeric(ord$start))
        ord <- ord[!is.na(ord$chr_num) & !is.na(ord$start), , drop = FALSE]

        matched <- ord$gene %in% rownames(object)
        if (!any(matched)) {
            stop(sprintf(
                "No gene in `gene_order` is present in the object (%d ids checked, %d object features). Check the gene id style: `gene_order` has %s, the object has %s.",
                length(ord$gene), nrow(object),
                paste(utils::head(unique(ord$gene), 3L), collapse = ", "),
                paste(utils::head(rownames(object), 3L), collapse = ", ")
            ))
        }
        # Duplicated gene rows are normal in GTF/Ensembl exports. inferCNV pairs
        # expression rows with order rows positionally after its own match, so
        # duplicated ids must not be forwarded: keep the widest interval, which is
        # the union of the annotated transcripts for that gene.
        dup <- duplicated(ord$gene) | duplicated(ord$gene, fromLast = TRUE)
        if (any(dup)) {
            ord <- ord[order(ord$gene, ord$start), , drop = FALSE]
            agg <- stats::aggregate(
                cbind(chr_num, start, end) ~ gene,
                data = ord,
                FUN = function(v) range(v, na.rm = TRUE)
            )
            # `chr_num` must survive: the sort below orders on it.
            ord <- data.frame(
                gene = agg$gene,
                chr = paste0("chr", agg$chr_num[, 1L]),
                chr_num = agg$chr_num[, 1L],
                start = agg$start[, 1L],
                end = agg$end[, 2L],
                stringsAsFactors = FALSE
            )
        }
        ord <- ord[ord$gene %in% rownames(object), , drop = FALSE]
        ord <- ord[order(ord$chr_num, ord$start), , drop = FALSE]
        out <- ord[, c("gene", "chr", "start", "end"), drop = FALSE]
        attr(out, "source") <- "caller-supplied"
        attr(out, "n_duplicate_genes") <- sum(dup)
        return(out)
    }

    # 2. Positions already in rowData (e.g. from a genome annotation).
    rd <- SummarizedExperiment::rowData(object)
    need <- c("chr", "start", "end")
    if (all(need %in% colnames(rd))) {
        ord <- data.frame(
            gene = rownames(object)[!is.na(rd$chr)],
            chr = as.character(rd$chr[!is.na(rd$chr)]),
            start = rd$start[!is.na(rd$chr)],
            end = rd$end[!is.na(rd$chr)],
            stringsAsFactors = FALSE
        )
        ord$chr <- paste0("chr", sub("^chr", "", ord$chr))
        ord <- ord[order(suppressWarnings(as.integer(sub("^chr", "", ord$chr))), ord$start), ]
        out <- ord[, c("gene", "chr", "start", "end"), drop = FALSE]
        attr(out, "source") <- "rowData"
        return(out)
    }

    # 3. TxDb fallback, only when the package happens to be installed.
    txdb_pkg <- switch(species,
                       human = "TxDb.Hsapiens.UCSC.hg38.knownGene",
                       mouse = "TxDb.Mmusculus.UCSC.mm10.knownGene")
    if (requireNamespace(txdb_pkg, quietly = TRUE) &&
        requireNamespace("GenomicFeatures", quietly = TRUE)) {
        txdb <- getExportedValue(txdb_pkg, txdb_pkg)
        g <- tryCatch(
            suppressMessages(GenomicFeatures::genes(txdb, columns = "gene_id", single.strand.genes.only = TRUE)),
            error = function(e) NULL
        )
        if (!is.null(g)) {
            symbols <- tryCatch(
                suppressMessages(AnnotationDbi::mapIds(
                    getExportedValue("org.Hs.eg.db", "org.Hs.eg.db"),
                    keys = as.character(names(g)),
                    column = "SYMBOL", keytype = "ENTREZID", multiVals = "first"
                )),
                error = function(e) NULL
            )
            if (!is.null(symbols)) {
                ord <- data.frame(
                    gene = unname(symbols),
                    chr = paste0("chr", as.character(GenomeInfoDb::seqnames(g))),
                    start = as.integer(GenomicRanges::start(g)),
                    end = as.integer(GenomicRanges::end(g)),
                    stringsAsFactors = FALSE
                )
                ord <- ord[!is.na(ord$gene) & ord$gene %in% rownames(object), , drop = FALSE]
                out <- ord[, c("gene", "chr", "start", "end"), drop = FALSE]
                attr(out, "source") <- paste0("TxDb:", txdb_pkg)
                return(out)
            }
        }
    }

    empty <- data.frame(gene = character(0), chr = character(0),
                        start = integer(0), end = integer(0), stringsAsFactors = FALSE)
    attr(empty, "source") <- "none"
    empty
}

# Per-cell burden: mean absolute deviation of the smoothed expression from the
# CELL'S OWN MEDIAN across genes.
#
# `expr.data` holds reference-centred, smoothed values whose level hovers around
# 1 for every cell; taking abs() of it measures the baseline, not the deviation,
# and on inferCNV's own bundled example the malignant and normal groups differ by
# -0.005 (measured). Subtracting the per-cell median first removes each cell's
# offset and leaves the shape of the profile, which separates the groups by
# +0.038 on that same data. Use a median rather than the reference mean so a
# large amplified or deleted tract cannot pull the baseline toward itself.
sclet_cnv_burden <- function(cnv_obj, cells, ref) {
    # Both `infercnv` and `MCMC_inferCNV` carry the processed matrix in
    # `expr.data`; the BayesNet subclass inherits the slot.
    expr <- methods::slot(cnv_obj, "expr.data")
    baseline <- apply(expr, 2, stats::median, na.rm = TRUE)
    centred <- sweep(expr, 2, baseline, "-")
    burden <- colMeans(abs(centred), na.rm = TRUE)
    obs_cells <- intersect(colnames(expr), cells)
    data.frame(
        cell = obs_cells,
        burden = as.numeric(burden[obs_cells]),
        is_reference = obs_cells %in% ref,
        stringsAsFactors = FALSE
    )
}


sclet_cnv_record_unavailable <- function(object, name, readiness) {
    cdata <- SummarizedExperiment::colData(object)
    cdata[[paste0(name, "_burden")]] <- rep(NA_real_, nrow(cdata))
    cdata[[paste0(name, "_is_reference")]] <- rep(NA, nrow(cdata))
    SummarizedExperiment::colData(object) <- cdata

    result <- list(
        status = "unavailable",
        reason = readiness$issues,
        readiness = readiness
    )
    # `analyses[["malignant_cnv"]]` is a MAP keyed by run id, mirroring
    # `metadata$malignant_cnv[[name]]`, so two runs do not overwrite each other
    # and the reader can address one by name.
    existing <- sclet_get_analysis(object, "malignant_cnv")
    if (is.null(existing)) existing <- list()
    existing[[name]] <- result
    object <- sclet_set_analysis(object, "malignant_cnv", existing)
    md <- S4Vectors::metadata(object)
    md$malignant_cnv <- if (is.null(md$malignant_cnv)) list() else md$malignant_cnv
    md$malignant_cnv[[name]] <- result
    S4Vectors::metadata(object) <- md

    # Its own `cnv` record for the failed run too, for the same reason: the
    # annotation bucket must keep answering "what is this object annotated as?".
    object <- sclet_set_state_record(
        object = object,
        type = "cnv",
        id = name,
        value = list(
            method = "inferCNV per-cell copy-number burden",
            status = "unavailable",
            reason = readiness$issues,
            artifacts = list(columns = paste0(name, c("_burden", "_is_reference"))),
            summary = list(status = "unavailable", n_issues = length(readiness$issues)),
            # A stamp is not optional: the AI trace diff detects a same-id
            # rewrite by `created_at`, so a record without one makes a re-run
            # invisible to the graph (measured 2026-10-23).
            created_at = Sys.time()
        )
    )
    warning(
        paste0(
            "CNV inference did not run; ", name, "_burden is NA (not zero).\n  - ",
            paste(readiness$issues, collapse = "\n  - ")
        ),
        call. = FALSE
    )
    object
}
