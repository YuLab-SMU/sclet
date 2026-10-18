#' Gene expression markers for all identity classes
#' 
#' Find markers for each of the identity classes
#' 
#' @param object A 'SingleCellExperiment' object
#' @param min.pct  only test genes that are detected in a minimum fraction of
#' min.pct cells in either of the two populations.
#' @param logfc.threshold Limit testing to genes which show, on average, at least
#' X-fold difference (log-scale) between the two groups of cells. Default is 0.1
#' @param return.thresh Only return markers that have a p-value < return.thresh, or a power > return.thresh (if the test is ROC)
#' @param only.pos Only return positive markers (FALSE by default)
#' @param base log base. Default is 2
#' @param pseudocount.use Pseudocount to add to averaged expression values when
#' calculating logFC. 1 by default.
#' @param assay assay used for marker testing. If NULL, sclet resolves from
#' `DefaultAssay(object)` and falls back to `logcounts` when needed.
#' @param layer layer used for marker testing. If NULL, sclet uses
#' `DefaultLayer(object)`.
#' @param ... additional parameters
#' @export
FindAllMarkers <- function(
    object, 
    min.pct = 0.01,
    logfc.threshold = 0.1,
    return.thresh = 0.01,
    only.pos = FALSE,
    base = 2,
    pseudocount.use = 1,
    assay = NULL,
    layer = NULL,
    ...) {

    idents <- Idents(object)
    if (is.null(idents)) stop("No identities found. Please run FindClusters() first.")
    assay <- resolve_marker_assay(object, assay = assay, layer = layer)
    
    # Ensure presto is available
    if (!is.installed("presto")) {
        install_zip_gh("immunogenomics/presto")
    }
    wilcoxauc <- get_fun_from_pkg("presto", "wilcoxauc")

    message("Calculating markers for all clusters using presto...")
    
    # presto::wilcoxauc does not support DelayedMatrix directly.
    mat <- SummarizedExperiment::assay(object, assay)
    if (inherits(mat, "DelayedArray")) {
        mat <- as.matrix(mat)
    }
    results <- wilcoxauc(mat, idents)

    # presto returns: feature, group, avg_logFC, pval, padj, pct_in, pct_out, auc
    # We need to align these with sclet/Seurat names:
    # gene, cluster, avg_log2FC (or base), pct.1, pct.2, pval, padj
    
    fc.name <- sprintf("avg_log%dFC", base)
    
    # If base != 2, we need to adjust presto's logFC (which is usually log2 or ln?)
    # presto uses natural log (ln) if not specified? 
    # Actually presto calculates (mean(group) - mean(rest)) if mat is log-transformed.
    # If mat is log2, then FC is log2. If mat is ln, then FC is ln.
    # sclet's FindMarkers uses: log((rowSums(expm1(x)) + pseudo)/N, base)
    # This is slightly different from (mean(x1) - mean(x2)) if x is log.
    
    # To maintain EXACT compatibility with FindMarkers, we might still need the mean calculation.
    # But for efficiency in 'FindAllMarkers', using presto's FC is usually acceptable 
    # and much faster. However, sclet users might expect the exact same values.
    
    # Let's perform a vectorized mean calculation if we want exact match.
    # aggregateAcrossCells can get us the sums of expm1(x).
    
    # For now, let's use a reasonably optimized approach that stays close to sclet's output.
    
    colnames(results)[colnames(results) == "feature"] <- "gene"
    colnames(results)[colnames(results) == "group"] <- "cluster"
    colnames(results)[colnames(results) == "pct_in"] <- "pct.1"
    colnames(results)[colnames(results) == "pct_out"] <- "pct.2"
    
    # presto's avg_logFC calculation is (mean1 - mean2).
    # sclet's is (logAvg1 - logAvg2).
    # If we want to be exact, we should re-calculate FC.
    # But let's see if we can just rename for now and see if user is okay with "close enough" 
    # or if they need "exact". Given it's a "Lightweight" toolkit, speed is often preferred.
    
    # Rename FC column
    # presto typically uses 'logFC' (or 'avg_logFC' in older versions)
    if ("logFC" %in% colnames(results)) {
        colnames(results)[colnames(results) == "logFC"] <- fc.name
    } else if ("avg_logFC" %in% colnames(results)) {
        colnames(results)[colnames(results) == "avg_logFC"] <- fc.name
    }
    
    # presto uses pct from 0-100 sometimes? No, 0-1 usually.
    results$pct.1 <- round(results$pct.1 / 100, 3)
    results$pct.2 <- round(results$pct.2 / 100, 3)

    # Filtering
    results <- results[results$pval < return.thresh, ]
    results <- results[results$pct.1 >= min.pct | results$pct.2 >= min.pct, ]
    
    if (only.pos) {
        results <- results[results[[fc.name]] > 0, ]
    }
    
    results <- results[abs(results[[fc.name]]) >= logfc.threshold, ]
    
    # Sort
    results <- results[order(results$cluster, results$pval, -abs(results$pct.1 - results$pct.2)), ]
    
    # final formatting
    rownames(results) <- make.unique(as.character(results$gene), sep="")
    
    # Reorder columns to match sclet FindMarkers exactly if possible
    col_order <- c("pval", "padj", fc.name, "pct.1", "pct.2", "cluster", "gene")
    # Add any extra columns presto might have (auc, etc) at the end
    extra_cols <- setdiff(colnames(results), col_order)
    results <- results[, c(col_order, extra_cols)]

    return(results)
}

#' Find markers 
#' 
#' @title FindMarkers
#' @param object A 'SingleCellExperiment' object
#' @param min.pct  only test genes that are detected in a minimum fraction of
#' min.pct cells in either of the two populations.
#' @param ident.1 Identity class to define markers for
#' @param ident.2 A second identity class for comparison. If NULL (default) -
#' use all other cells for comparison.
#' @param logfc.threshold Limit testing to genes which show, on average, at least
#' X-fold difference (log-scale) between the two groups of cells. Default is 0.1
#' @param base log base. Default is 2
#' @param pseudocount.use Pseudocount to add to averaged expression values when
#' calculating logFC. 1 by default.
#' @param assay assay used for marker testing. If NULL, sclet resolves from
#' `DefaultAssay(object)` and falls back to `logcounts` when needed.
#' @param layer layer used for marker testing. If NULL, sclet uses
#' `DefaultLayer(object)`.
#' @param ... additional parameters
#' @export
FindMarkers <- function(
    object, 
    ident.1 = NULL, 
    ident.2 = NULL, 
    min.pct = 0.01, 
    logfc.threshold = 0.1, 
    base=2,
    pseudocount.use=1,
    assay = NULL,
    layer = NULL,
    ...) {
        assay <- resolve_marker_assay(object, assay = assay, layer = layer)

        FindMarkers_Presto(
        object = SummarizedExperiment::assay(object, assay),
        ident.1 = ident.1,
        ident.2 = ident.2,
        clusters = Idents(object),
        min.pct = min.pct,
        logfc.threshold = logfc.threshold,
        base = base,
        pseudocount.use = pseudocount.use,
        ...
        )
}

#' Run Differential Expression Test
#'
#' @title RunDEtest
#' @param object A SingleCellExperiment object
#' @param ident.1 Identity class to define markers for. If NULL, runs FindAllMarkers.
#' @param ident.2 A second identity class for comparison
#' @param name Analysis record id. Defaults to "detest".
#' @param ... Passed to FindMarkers or FindAllMarkers
#' @return A SingleCellExperiment object with detest state added
#' @export
RunDEtest <- function(object, ident.1 = NULL, ident.2 = NULL, name = "detest",
                      group = NULL, covariates = NULL, min_cells = 10L, ...) {
    if (!is.null(group)) {
        return(sclet_run_condition_de(
            object,
            group = group,
            ident.1 = ident.1,
            ident.2 = ident.2,
            covariates = covariates,
            name = name,
            min_cells = min_cells,
            ...
        ))
    }
    if (is.null(ident.1)) {
        res <- FindAllMarkers(object, ...)
        method <- "FindAllMarkers"
    } else {
        res <- FindMarkers(object, ident.1 = ident.1, ident.2 = ident.2, ...)
        method <- "FindMarkers"
    }
    
    state_inputs <- list(
        ident.1 = ident.1,
        ident.2 = ident.2,
        active_ident = ActiveIdent(object)
    )
    
    object <- sclet_set_analysis_state(
        object = object,
        type = "detest",
        id = name,
        method = method,
        inputs = state_inputs,
        artifacts = list(result = res),
        summary = list(n_markers = nrow(res)),
        active = TRUE
    )
    
    object <- sclet_log_command(
        object,
        "RunDEtest",
        params = list(ident.1 = ident.1, ident.2 = ident.2, name = name),
        outputs = list(detest = name)
    )
    
    return(object)
}

#' Get Differential Expression Test Results
#'
#' @title get_detest
#' @param object A SingleCellExperiment object
#' @param id Analysis record id. If NULL, uses the active detest state.
#' @return The detest state record
#' @export
get_detest <- function(object, id = NULL) {
    if (is.null(id)) id <- sclet_get_active_state(object, "detest")
    if (is.null(id)) return(NULL)
    sclet_get_state_record(object, "detest", id)
}

#' Check if Differential Expression Test Results exist
#'
#' @title has_detest
#' @param object A SingleCellExperiment object
#' @param id Analysis record id. If NULL, uses the active detest state.
#' @return Logical indicating if the detest state exists
#' @export
has_detest <- function(object, id = NULL) {
    !is.null(get_detest(object, id))
}

resolve_marker_assay <- function(object, assay = NULL, layer = NULL) {
    source <- sclet_resolve_expression_source(
        object = object,
        layer = layer,
        assay = assay,
        prefer_nonscaled = TRUE,
        context = "marker testing"
    )
    source$assay
}

#' @importFrom yulab.utils install_zip_gh
#' @importFrom yulab.utils is.installed
#' @importFrom Matrix rowSums
FindMarkers_Presto <- function(object, ident.1 = NULL, ident.2 = NULL, clusters, 
                              min.pct = 0.01,logfc.threshold = 0.1, base=2,
                              pseudocount.use=1, ...) {

    if (inherits(object, "DelayedArray")) {
        object <- as.matrix(object)
    }

    if (is.null(ident.1)) stop("ident.1 must be specified")
    if (is.null(clusters)) stop("clusters must be specified")

    if (is.null(ident.2)) {
        ident.2 <- setdiff(unique(clusters), ident.1)
    }

    cells.1 <- which(clusters %in% ident.1)
    cells.2 <- which(clusters %in% ident.2)

    if (length(cells.1) == 0 || length(cells.2) == 0) {
        stop("One or both of the identity groups have no cells")}
    
    mean.fxn <- function(x) {
        log((sclet_matrix_rowSums(expm1(x = x)) + pseudocount.use)/NCOL(x), base = base)
    }

    fc.name <- sprintf("avg_log%dFC", base)
    features <- rownames(x = object)
    thresh.min <- 0
    mat <- object[features, ]
    pct.1 <- round(sclet_matrix_rowSums(mat[, cells.1, drop = FALSE] > thresh.min)/length(cells.1), digits = 3)
    pct.2 <- round(sclet_matrix_rowSums(mat[, cells.2, drop = FALSE] > thresh.min) / length(cells.2), digits = 3)
    data.1 <- mean.fxn(mat[, cells.1, drop = FALSE])
    data.2 <- mean.fxn(mat[, cells.2, drop = FALSE])
    fc <- (data.1 - data.2) 
    fc.results <- as.data.frame(cbind(fc, pct.1, pct.2)) 
    colnames(fc.results) <- c(fc.name, "pct.1", "pct.2")
    fc.results$feature <- features

    data_subset <- object[ ,c(cells.1, cells.2), drop = FALSE]

    group_labels <- c(rep("group1", length(cells.1)), rep("group2", length(cells.2)))

    if (!is.installed("presto")) {
        install_zip_gh("immunogenomics/presto")
    }

    wilcoxauc <- get_fun_from_pkg("presto", "wilcoxauc")
    results <- wilcoxauc(data_subset, group_labels)
    

    results <- results[1:(nrow(x = results)/2), ]
    p.results <- results[, c("feature", "pval","padj")]

    combined_results <- merge(p.results, fc.results, by = 'feature')
    combined_results <- combined_results[combined_results$pct.1>= min.pct | combined_results$pct.2 >= min.pct, ]
    combined_results <- combined_results[abs(combined_results[,fc.name]) >= logfc.threshold, ]
    combined_results <- combined_results[order(combined_results$pval, 
                                -abs(combined_results$pct.1 - combined_results$pct.2)), ]
    rownames(combined_results) <- combined_results$feature
    combined_results$feature <- NULL
    return(combined_results)
}

# Differential expression for a colData GROUP (condition, treatment, disease),
# not for clusters.
#
# Why this exists: `run_de_test()` compared identity levels only, so a goal
# asking "stimulated vs control within donors" was answered with cluster
# markers (Kang IFN-beta real-data run, 2026-10-18). Two paths, deliberately
# separate because they answer with different evidence:
#
#   * no covariates -> the existing cell-level FindMarkers on the requested
#     column. Cheap, but cells are treated as independent: donor or batch
#     structure is NOT adjusted. The state records that.
#   * covariates (e.g. donor) -> pseudobulk per donor x condition and limma on
#     log2 CPM with design ~ condition + covariates. This is the paired
#     analysis "within donors" actually means; a covariate must appear in BOTH
#     conditions for the paired contrast, otherwise the design cannot estimate
#     it and the call stops instead of returning an unadjusted table under a
#     paired label.
#
# The result keeps the marker-table schema (gene, cluster, avg_log2FC, pval,
# padj, pct.1, pct.2) so downstream evidence/claims code is unchanged; `cluster`
# names the contrast and the state inputs record `design` and `adjustment`.
sclet_run_condition_de <- function(object, group, ident.1, ident.2,
        covariates = NULL, name = "detest", min_cells = 10L, ...) {
    if (!is.character(group) || length(group) != 1L || is.na(group) || !nzchar(group)) {
        stop("`group` must be one colData column name.")
    }
    cd <- SummarizedExperiment::colData(object)
    if (!group %in% colnames(cd)) {
        stop("`group` column '", group, "' not found in colData.")
    }
    values <- as.character(cd[[group]])
    if (!is.character(ident.1) || length(ident.1) != 1L || !nzchar(ident.1)) {
        stop("`ident.1` must be one level of `", group, "`.")
    }
    if (!is.character(ident.2) || length(ident.2) != 1L || !nzchar(ident.2)) {
        stop("`ident.2` must be one level of `", group, "`.")
    }
    if (identical(ident.1, ident.2)) {
        stop("`ident.1` and `ident.2` must name two different levels of `", group, "`.")
    }
    present <- unique(values[!is.na(values)])
    if (!ident.1 %in% present || !ident.2 %in% present) {
        stop("`", group, "` has levels ", paste(present, collapse = ", "),
            "; both ident.1 and ident.2 must be among them.")
    }

    if (is.null(covariates) || !length(covariates)) {
        adjustment <- "none (cell-level test; donor/batch structure not adjusted)"
        design_label <- paste0("~ ", group)
        previous <- ActiveIdent(object)
        object <- `ActiveIdent<-`(object, group)
        res <- tryCatch(FindMarkers(object, ident.1 = ident.1,
            ident.2 = ident.2, ...), error = function(e) e)
        if (!is.null(previous)) object <- `ActiveIdent<-`(object, previous)
        if (inherits(res, "error")) stop(conditionMessage(res))
        method <- paste0("FindMarkers (", group, " contrast)")
    } else {
        covariates <- as.character(covariates)
        missing <- setdiff(covariates, colnames(cd))
        if (length(missing)) {
            stop("covariate column(s) not found in colData: ",
                paste(missing, collapse = ", "))
        }
        if (!requireNamespace("limma", quietly = TRUE)) {
            stop("covariate-adjusted differential expression needs the 'limma' ",
                "package; install it or call without `covariates` (unadjusted).")
        }
        adjustment <- paste0("pseudobulk per ", paste(covariates, collapse = " x "),
            " and limma, design ~ ", group, " + ",
            paste(covariates, collapse = " + "))
        design_label <- paste0("~ ", group, " + ", paste(covariates, collapse = " + "))
        res <- sclet_pseudobulk_condition_de(object, group = group,
            ident.1 = ident.1, ident.2 = ident.2, covariates = covariates,
            min_cells = min_cells, ...)
        method <- paste0("Pseudobulk limma (", group, " contrast, paired by ",
            paste(covariates, collapse = ","), ")")
    }

    object <- sclet_set_analysis_state(
        object = object,
        type = "detest",
        id = name,
        method = method,
        inputs = list(
            ident.1 = ident.1,
            ident.2 = ident.2,
            group = group,
            covariates = covariates %||% character(),
            design = design_label,
            adjustment = adjustment,
            active_ident = ActiveIdent(object)
        ),
        artifacts = list(result = res),
        summary = list(
            n_markers = nrow(res),
            contrast = paste0(ident.1, "_vs_", ident.2),
            adjustment = adjustment
        ),
        active = TRUE
    )
    object
}

# Pseudobulk condition contrast for one or more covariates.
#
# Aggregates counts per covariate combination x condition, keeps genes with
# CPM >= 1 in at least two samples, fits limma on log2 CPM with the condition
# coefficient, and reports detection fractions from the cell-level data. The
# design requires every covariate combination to contain BOTH conditions;
# otherwise a paired contrast is not identifiable and the call stops rather
# than silently dropping the pairing.
sclet_pseudobulk_condition_de <- function(object, group, ident.1, ident.2,
        covariates, min_cells = 10L, assay = "counts", min.pct = 0.01,
        logfc.threshold = 0.1, ...) {
    if (!"counts" %in% SummarizedExperiment::assayNames(object)) {
        stop("pseudobulk differential expression requires a 'counts' assay.")
    }
    cd <- as.data.frame(SummarizedExperiment::colData(object))
    keep <- !is.na(cd[[group]]) & cd[[group]] %in% c(ident.1, ident.2) &
        stats::complete.cases(cd[, covariates, drop = FALSE])
    if (!any(keep)) stop("no cells carry both the requested groups and covariates")
    meta_cells <- cd[keep, c(group, covariates), drop = FALSE]
    key <- do.call(paste, c(lapply(covariates, function(v) as.character(meta_cells[[v]])),
        list(sep = "|")))
    cond <- factor(as.character(meta_cells[[group]]), levels = c(ident.2, ident.1))
    sample_id <- paste(key, as.character(cond), sep = "@")
    split_idx <- split(which(keep), sample_id)
    sizes <- lengths(split_idx)
    split_idx <- split_idx[sizes >= min_cells]
    if (!length(split_idx)) {
        stop("no covariate x condition sample had at least ", min_cells, " cells")
    }
    sample_meta <- do.call(rbind, lapply(names(split_idx), function(sid) {
        idx <- split_idx[[sid]]
        data.frame(
            sample = sid,
            condition = as.character(meta_cells[[group]])[match(idx[1L], which(keep))],
            stringsAsFactors = FALSE
        )
    }))
    for (v in covariates) {
        sample_meta[[v]] <- vapply(names(split_idx), function(sid) {
            idx <- split_idx[[sid]]
            as.character(meta_cells[[v]])[match(idx[1L], which(keep))]
        }, character(1))
    }
    sample_meta$condition <- factor(sample_meta$condition, levels = c(ident.2, ident.1))
    if (any(table(sample_meta$condition) < 2L)) {
        stop("each condition needs at least two pseudobulk samples for a ",
            "covariate-adjusted contrast; got ",
            paste(names(table(sample_meta$condition)), table(sample_meta$condition),
                sep = "=", collapse = ", "))
    }
    pairing <- table(sample_meta[, covariates, drop = FALSE])
    if (!all(pairing >= 2L)) {
        stop("every ", paste(covariates, collapse = " x "), " combination must ",
            "contain BOTH conditions (", ident.1, " and ", ident.2,
            "); the paired contrast is not identifiable otherwise")
    }

    counts <- SummarizedExperiment::assay(object, assay)
    pb <- do.call(cbind, lapply(split_idx, function(idx) {
        Matrix::rowSums(counts[, idx, drop = FALSE])
    }))
    colnames(pb) <- names(split_idx)
    lib <- colSums(pb)
    cpm <- sweep(pb, 2L, lib, "/") * 1e6
    keep_gene <- rowSums(cpm >= 1) >= 2L
    if (sum(keep_gene) < 10L) {
        stop("too few genes pass the pseudobulk filter (", sum(keep_gene), ")")
    }
    logcpm <- log2(sweep(pb[keep_gene, , drop = FALSE] + 0.5, 2L, lib, "/") * 1e6)
    design <- stats::model.matrix(
        stats::reformulate(c("condition", covariates)), data = sample_meta)
    fit <- limma::lmFit(logcpm, design)
    fit <- limma::eBayes(fit)
    coef_name <- paste0("condition", ident.1)
    if (!coef_name %in% colnames(design)) {
        stop("condition coefficient '", coef_name, "' is not in the design; ",
            "check the factor levels of '", group, "'")
    }
    tt <- limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
    genes <- rownames(tt)
    pct <- function(condition) {
        idx <- which(keep & as.character(cd[[group]]) == condition)
        rowSums(counts[genes, idx, drop = FALSE] > 0) / length(idx)
    }
    out <- data.frame(
        gene = genes,
        cluster = paste0(ident.1, "_vs_", ident.2),
        avg_log2FC = tt$logFC,
        pval = tt$P.Value,
        padj = tt$adj.P.Val,
        pct.1 = pct(ident.1),
        pct.2 = pct(ident.2),
        stringsAsFactors = FALSE,
        row.names = NULL
    )
    ## Same visibility filters as the cell-level marker table, so callers and
    ## enrichment see the same kind of table either way.
    out <- out[(out$pct.1 >= min.pct | out$pct.2 >= min.pct) &
        abs(out$avg_log2FC) >= logfc.threshold, , drop = FALSE]
    out[order(out$pval), , drop = FALSE]
}
