#' add QC metrics
#' 
#' @title QCMetrics
#' @param object a 'SingleCellExperiment' object
#' @return update object with two QC metrics ('nFeature_RNA' and 'nCount_RNA')
#' @export
QCMetrics <- function(object) {
    qc_metrics <- scrapper::computeRnaQcMetrics(
        SummarizedExperiment::assay(object, "counts"),
        subsets = list()
    )
    SummarizedExperiment::colData(object)$nFeature_RNA <- qc_metrics$detected
    SummarizedExperiment::colData(object)$nCount_RNA <- qc_metrics$sum
    return(object)
}

#' calculate percentage of genes that matched a pattern
#' 
#' @title PercentageFeatureSet
#' @param object a SingleCellExperiment object
#' @param pattern the pattern to search for
#' @param feature search pattern from feature of the object. If NULL, use the `rownames(object)`, otherwise, use `rowData(object)[[feature]]`
#' @return the percentage of each cell
#' @export
PercentageFeatureSet <- function(object, pattern = NULL, feature=NULL) {
    if (is.null(feature)) {
        features <- rownames(object)
    } else {
        features <- rowData(object)[[feature]]
    }

    has_pattern <- grep(pattern, features)
    counts_mat <- SummarizedExperiment::assay(object, "counts")
    total_counts <- MatrixGenerics::colSums(counts_mat)
    if (length(has_pattern) == 0) {
        return(rep(0, ncol(object)))
    }
    subset_counts <- MatrixGenerics::colSums(counts_mat[has_pattern, , drop = FALSE])
    percent <- (subset_counts / total_counts) * 100
    percent[is.na(percent)] <- 0
    return(percent)
}

#' normalize data
#' 
#' @title NormalizeData
#' @param object a SingleCellExperiment object
#' @param scale.factor scale factor
#' @return a SingleCellExperiment object with 'logcounts' assay created/updated.
#' @importFrom SummarizedExperiment 'assay<-'
#' @importFrom yulab.utils get_fun_from_pkg
#' @importFrom SummarizedExperiment assay
#' @importFrom methods as
#' @export
NormalizeData <- function(object, scale.factor = 10000) {
    libsize <- MatrixGenerics::colSums(SummarizedExperiment::assay(object, "counts"))
    size_factors <- libsize / scale.factor
    prev_state <- sclet_get_state(object)
    counts_mat <- SummarizedExperiment::assay(object, "counts")
    if (is.matrix(counts_mat) && !isS4(counts_mat)) {
        logcounts <- t(t(counts_mat) / size_factors)
    } else if (methods::is(counts_mat, "Matrix")) {
        logcounts <- counts_mat %*% Matrix::Diagonal(x = 1 / size_factors)
    } else {
        logcounts <- DelayedArray::sweep(
            DelayedArray::DelayedArray(counts_mat),
            MARGIN = 2L,
            STATS = size_factors,
            FUN = "/"
        )
    }
    logcounts <- log1p(logcounts)
    dimnames(logcounts) <- dimnames(object)
    SummarizedExperiment::assay(object, "logcounts") <- logcounts
    object <- sclet_restore_state(object, prev_state)
    object <- sclet_set_layer(object, name = "counts", assay = "counts", role = "counts", active = FALSE)
    object <- sclet_set_layer(
        object,
        name = "logcounts",
        assay = "logcounts",
        role = "normalized",
        params = list(scale.factor = scale.factor)
    )
    object <- sclet_set_active_assay(object, "logcounts")
    object <- sclet_set_analysis_state(
        object = object,
        type = "preprocess",
        id = "normalize_logcounts",
        method = "scuttle::calculateCPM + log1p",
        inputs = list(
            assay = "counts"
        ),
        artifacts = list(
            assay = "logcounts"
        ),
        params = list(
            scale.factor = scale.factor
        ),
        summary = list(
            n_cells = ncol(object),
            n_features = nrow(object)
        )
    )
    object <- sclet_log_command(
        object,
        "NormalizeData",
        params = list(scale.factor = scale.factor)
    )
    
    return(object)
}

#' identify variable features
#' 
#' @title FindVariableFeatures
#' @param object a SingleCellExperiment object
#' @param nfeatures number of features to be selected as highly variable features
#' @param method one of 'scran' or 'scrapper'. The legacy value `'seurat'`
#'   is accepted for transition but is deprecated and redirected to `'scran'`.
#' @param ... additional parameters for variance modelling. For `method = "scran"`,
#'   they are passed to `scran::modelGeneVar()`. For `method = "scrapper"`,
#'   function-level arguments such as `block`, `num.threads`, `assay.type`,
#'   `include.per.block`, `more.var.args` and `more.choose.args` are supported;
#'   any remaining named arguments are forwarded into `more.var.args`.
#' @return an updated SingleCellExperiment object with identified highly variable features
#' @importFrom SingleCellExperiment counts
#' @importFrom SingleCellExperiment logcounts
#' @importFrom SummarizedExperiment 'rowData<-'
#' @importFrom stats loess
#' @importFrom scran modelGeneVar
#' @export
FindVariableFeatures <- function(object, nfeatures = 2000, method = "scran", ...) {
    method <- match.arg(method, c("seurat", "scran", "scrapper"))
    if (identical(method, "seurat")) {
        warning(
            "'method = \"seurat\"' is deprecated in sclet and now redirects to 'scran'.",
            call. = FALSE
        )
        method <- "scran"
    }
    prev_state <- sclet_get_state(object)

    rn <- rownames(object)
    if (is.null(rn)) {
        rn <- as.character(seq_len(nrow(object)))
        rownames(object) <- rn
    }

    if (method == "scrapper") {
        if (!requireNamespace("scrapper", quietly = TRUE)) {
            warning("scrapper is not installed. Falling back to scran method.")
            hvf.info <- sclet_model_gene_var(object, ...)
            method <- "scran"
        } else {
            hvf.info <- sclet_choose_hvgs_scrapper(object, nfeatures = nfeatures, ...)
        }
    } else {
        hvf.info <- sclet_model_gene_var(object, ...)
    }

    hvf.info <- hvf.info[rn, , drop = FALSE]
    rd <- rowData(object)
    SummarizedExperiment::rowData(object) <- cbind(rd, hvf.info)
    object <- sclet_restore_state(object, prev_state)
    object <- sclet_set_hvg_state(
        object,
        nfeatures = nfeatures,
        method = method,
        hvgcols = names(hvf.info)
    )
    object <- sclet_log_command(
        object,
        "FindVariableFeatures",
        params = list(nfeatures = nfeatures, method = method)
    )

    return(object)
}


#' get variable features
#' 
#' @title VariableFeatures
#' @param object a SingleCellExperiment object
#' @param method one of 'scran' or 'scrapper'. The legacy value `'seurat'`
#'   remains available when requested explicitly for older objects, but default
#'   resolution prefers Bioconductor-native HVG metadata.
#' @param ... additional parameters for ranking the `method = "scran"` HVG table
#' @return highly variable features
#' @export
VariableFeatures <- function(object, method = "scran", ...) {
    if (missing(method)) {
        method <- sclet_resolve_hvg_method(object)
    } else {
        method <- sclet_resolve_hvg_method(object, method = method, allow_legacy = TRUE)
    }
    if (is.null(method)) {
        stop(
            "No Bioconductor-native HVG statistics found. Re-run `FindVariableFeatures()` ",
            "with `method = \"scran\"` or `method = \"scrapper\"`, or request ",
            "`method = \"seurat\"` explicitly for legacy objects."
        )
    }

    nfeatures <- sclet_get_hvg_nfeatures(object)
    if (is.null(nfeatures)) {
        stop("You should run 'FindVariableFeatures' first.")
    }

    if (method == "scran") {
        # Pass rowData instead of object to avoid assay access issues 
        # when assays are missing (e.g. after BatchRemover)
        res <- choose_hvgs_by_variance(SummarizedExperiment::rowData(object), n = nfeatures, ...)
    } else if (method == "scrapper") {
        res <- getTopHVGs_scrapper(object, n = nfeatures)
    } else {
        res <- getTopHVGs_seurat(object, n=nfeatures)
    }

    return(res)
}

getTopHVGs_seurat <- function(object, n) {
    if (is(object, "SingleCellExperiment")) {
        d <- rowData(object)
        d <- d[rownames(d) %in% rownames(object), ]
    } else {
        d <- object
    }
    i <- order(d$variance.standardized, decreasing = TRUE)
    d <- d[i, ]
    res <- rownames(d)[1:n]
    return(res)
}

sclet_choose_hvgs_scrapper <- function(object, nfeatures, ...) {
    dots <- list(...)
    wrapper_args <- c("block", "num.threads", "more.var.args", "more.choose.args",
        "assay.type", "output.prefix", "include.per.block")
    direct <- dots[intersect(names(dots), wrapper_args)]
    passthrough <- dots[setdiff(names(dots), c(wrapper_args, "top"))]

    if (is.null(direct$more.var.args)) {
        direct$more.var.args <- list()
    }
    direct$more.var.args <- c(direct$more.var.args, passthrough)
    direct$top <- nfeatures

    choose_fun <- get_fun_from_pkg("scrapper", "chooseRnaHvgs.se")
    updated <- do.call(choose_fun, c(list(x = object), direct))
    info <- SummarizedExperiment::rowData(updated)

    output.prefix <- direct$output.prefix
    if (is.null(output.prefix)) {
        output.prefix <- ""
    }
    wanted <- paste0(output.prefix, c("means", "variances", "fitted", "residuals", "hvg", "per.block"))
    wanted <- intersect(wanted, colnames(info))
    if (!any(grepl("hvg$", wanted))) {
        stop("scrapper::chooseRnaHvgs.se() did not produce an HVG indicator column.")
    }

    info[, wanted, drop = FALSE]
}

choose_hvgs_by_variance <- function(stats, var.field = "bio", n = NULL, prop = NULL,
                                    var.threshold = 0, fdr.field = "FDR",
                                    fdr.threshold = NULL,
                                    row.names = !is.null(rownames(stats))) {
    survivors <- seq_len(nrow(stats))
    if (!is.null(fdr.threshold)) {
        fdr <- stats[[fdr.field]]
        keep <- !is.na(fdr) & fdr <= fdr.threshold
        survivors <- survivors[keep]
        stats <- stats[keep, , drop = FALSE]
    }
    if (!is.null(var.threshold)) {
        var <- stats[[var.field]]
        keep <- !is.na(var) & var > var.threshold
        survivors <- survivors[keep]
        stats <- stats[keep, , drop = FALSE]
    }
    o <- order(stats[[var.field]], decreasing = TRUE)
    if (!is.null(n) || !is.null(prop)) {
        n <- max(c(n, round(prop * nrow(stats))), na.rm = TRUE)
        o <- utils::head(o, n)
    }
    if (row.names) {
        return(rownames(stats)[o])
    }
    survivors[o]
}

sclet_find_hvg_column <- function(stats, suffix, hvg_cols = NULL) {
    if (is.null(hvg_cols) || !length(hvg_cols)) {
        hvg_cols <- colnames(stats)
    }
    matches <- hvg_cols[grepl(paste0(suffix, "$"), hvg_cols)]
    if (!length(matches)) {
        return("")
    }
    matches[[1]]
}

getTopHVGs_scrapper <- function(object, n) {
    if (is(object, "SingleCellExperiment")) {
        stats <- SummarizedExperiment::rowData(object)
        hvg_cols <- sclet_get_hvg_cols(object)
    } else {
        stats <- object
        hvg_cols <- colnames(stats)
    }
    hvg_col <- sclet_find_hvg_column(stats, "hvg", hvg_cols = hvg_cols)
    if (!nzchar(hvg_col)) {
        stop("scrapper HVG indicator column not found. Please run `FindVariableFeatures(method = \"scrapper\")` first.")
    }
    residual_col <- sclet_find_hvg_column(stats, "residuals", hvg_cols = hvg_cols)
    keep <- which(stats[[hvg_col]])
    if (!length(keep)) {
        return(character())
    }
    if (nzchar(residual_col)) {
        keep <- keep[order(stats[[residual_col]][keep], decreasing = TRUE)]
    }
    rownames(stats)[utils::head(keep, n)]
}

#' @importFrom utils getFromNamespace
sclet_fit_trend_var <- function(means, vars, min.mean = 0.1, parametric = TRUE,
                                lowess = TRUE, density.weights = TRUE,
                                nls.args = list(), ...) {
    inverse_density_weights <- getFromNamespace(".inverse_density_weights", "scran")
    setup_nls_args <- getFromNamespace(".setup_nls_args", "scran")
    correct_logged_expectation <- getFromNamespace(".correct_logged_expectation", "scran")
    weighted_lowess <- getFromNamespace("weightedLowess", "limma")

    is.okay <- !is.na(vars) & vars > 1e-08 & means >= min.mean
    v <- vars[is.okay]
    m <- means[is.okay]

    if (density.weights) {
        w <- inverse_density_weights(m, adjust = 1)
    } else {
        w <- rep(1, length(m))
    }

    if (length(v) < 2L) {
        stop("need at least 2 points for non-parametric curve fitting")
    }

    to.fit <- log(v)
    left.edge <- min(m)
    param_fun <- function(x) {
        pmin(1, x / left.edge)
    }

    if (parametric) {
        if (length(v) <= 3L) {
            stop("need at least 4 points for non-linear curve fitting")
        }

        nls.args <- setup_nls_args(nls.args, start.args = list(vars = v, means = m))
        nls.args$formula <- v ~ (exp(A) * m)/(m^(1 + exp(N)) + exp(B))
        nls.args$weights <- w
        nls.args$control$warnOnly <- FALSE

        init.fit <- try(do.call(stats::nls, nls.args), silent = TRUE)
        if (is(init.fit, "try-error")) {
            aest <- exp(nls.args$start$A)
            best <- exp(nls.args$start$B)
            nest <- exp(nls.args$start$N) + 1
            param_fun <- function(x) {
                aest * x / (x^nest + best)
            }
            to.fit <- to.fit - log(param_fun(m))
        } else {
            to.fit <- to.fit - log(stats::fitted(init.fit))
            param_fun <- function(x) {
                stats::predict(init.fit, data.frame(m = x))
            }
        }
    } else if (!lowess) {
        stop("at least one of 'lowess' or 'parametric' must be 'TRUE'")
    }

    if (lowess) {
        lfit <- weighted_lowess(m, to.fit, weights = w, ...)
        lowess_fun <- stats::approxfun(m, lfit$fitted, rule = 2)
        unscaled_fun <- function(x) {
            exp(lowess_fun(x)) * param_fun(x)
        }
    } else {
        unscaled_fun <- param_fun
    }

    correct_logged_expectation(m, v, w, unscaled_fun)
}

sclet_combine_blocks <- function(blocks, ave.fields, pval.field, method,
                                 geometric, equiweight, weights, valid) {
    combine_parallel_pvalues <- getFromNamespace("combineParallelPValues", "metapod")

    if (length(blocks) == 1L) {
        return(blocks[[1]])
    }

    rn <- unique(lapply(blocks, rownames))
    if (length(rn) != 1L) {
        stop("gene identities should be the same")
    }

    if (equiweight) {
        weights <- rep(1, length(blocks))
    } else if (is.null(weights)) {
        stop("'weights' must be specified if 'equiweight=FALSE'")
    }

    original <- blocks
    if (length(unique(vapply(original, nrow, 0L))) != 1L) {
        stop("not all 'blocks' have the same number of rows")
    }
    if (!any(valid)) {
        stop("no entry of 'blocks' has positive weights")
    }

    blocks <- blocks[valid]
    weights <- weights[valid]
    combined <- list()

    for (field in ave.fields) {
        extracted <- lapply(blocks, "[[", i = field)
        if (geometric) {
            extracted <- lapply(extracted, log)
        }
        extracted <- mapply(
            "*",
            extracted,
            weights,
            SIMPLIFY = FALSE,
            USE.NAMES = FALSE
        )
        averaged <- Reduce("+", extracted) / sum(weights)
        if (geometric) {
            averaged <- exp(averaged)
        }
        combined[[field]] <- averaged
    }

    extracted <- lapply(blocks, "[[", i = pval.field)
    if (identical(method, "z")) {
        method <- "stouffer"
    } else if (identical(method, "holm-middle")) {
        method <- "holm-min"
    }

    combined$p.value <- combine_parallel_pvalues(
        extracted,
        method = method,
        weights = weights
    )$p.value
    combined$FDR <- stats::p.adjust(combined$p.value, method = "BH")

    output <- S4Vectors::DataFrame(combined, row.names = rn[[1]])
    S4Vectors::metadata(output)$per.block <- original
    output
}

sclet_combine_blocked_statistics <- function(collected, method, equiweight, ncells,
                                             geometric = FALSE,
                                             fields = c("mean", "total", "tech", "bio"),
                                             pval = "p.value") {
    sclet_combine_blocks(
        collected,
        method = method,
        equiweight = equiweight,
        weights = ncells,
        valid = ncells >= 2L,
        geometric = geometric,
        ave.fields = fields,
        pval.field = pval
    )
}

sclet_combine_var <- function(..., method = "fisher", pval.field = "p.value",
                              other.fields = NULL, equiweight = TRUE,
                              ncells = NULL) {
    unpack_lists <- getFromNamespace(".unpackLists", "scuttle")
    find_other_fields <- getFromNamespace(".find_other_fields", "scran")

    collected <- unpack_lists(...)
    if (is.null(ncells)) {
        ncells <- rep(10L, length(collected))
    }
    if (is.null(other.fields)) {
        other.fields <- find_other_fields(collected, c(pval.field, "FDR"))
    }

    sclet_combine_blocked_statistics(
        collected,
        method = method,
        equiweight = equiweight,
        ncells = ncells,
        geometric = FALSE,
        fields = other.fields,
        pval = pval.field
    )
}

sclet_subset2index <- function(subset.row, x) {
    if (is.null(subset.row)) {
        return(seq_len(nrow(x)))
    }
    if (is.logical(subset.row)) {
        if (length(subset.row) != nrow(x)) {
            stop("logical 'subset.row' should have length equal to nrow(x)")
        }
        return(which(subset.row))
    }
    if (is.numeric(subset.row)) {
        return(as.integer(subset.row))
    }
    if (is.character(subset.row)) {
        rn <- rownames(x)
        if (is.null(rn)) {
            stop("cannot resolve character 'subset.row' without rownames")
        }
        idx <- match(subset.row, rn)
        if (anyNA(idx)) {
            stop("some entries of 'subset.row' are not present in rownames(x)")
        }
        return(idx)
    }
    stop("unsupported 'subset.row' type")
}

sclet_muffle_known_warnings <- function(expr, patterns) {
    withCallingHandlers(
        expr,
        warning = function(w) {
            msg <- conditionMessage(w)
            if (any(vapply(patterns, grepl, logical(1), x = msg, fixed = TRUE))) {
                invokeRestart("muffleWarning")
            }
        }
    )
}

sclet_decompose_log_exprs <- function(x.means, x.vars, fit.means, fit.vars,
                                      ncells, ...) {
    dummy_trend_fit <- getFromNamespace("dummy.trend.fit", "scran")

    collected <- vector("list", ncol(x.means))
    for (i in seq_along(collected)) {
        fm <- fit.means[, i]
        fv <- fit.vars[, i]
        if (ncells[i] >= 2L) {
            fit <- sclet_fit_trend_var(fm, fv, ...)
        } else {
            fit <- dummy_trend_fit
        }

        xm <- unname(x.means[, i])
        xv <- unname(x.vars[, i])
        output <- S4Vectors::DataFrame(mean = xm, total = xv, tech = fit$trend(xm))
        output$bio <- output$total - output$tech
        output$p.value <- stats::pnorm(
            output$bio / output$tech,
            sd = fit$std.dev,
            lower.tail = FALSE
        )
        output$FDR <- stats::p.adjust(output$p.value, method = "BH")
        rownames(output) <- rownames(x.means)
        S4Vectors::metadata(output) <- c(list(mean = fm, var = fv), fit)
        collected[[i]] <- output
    }

    names(collected) <- colnames(x.means)
    collected
}

sclet_model_gene_var <- function(object, ..., assay.type = "logcounts",
                                 block = NULL, design = NULL,
                                 subset.row = NULL, subset.fit = NULL,
                                 equiweight = TRUE, method = "fisher",
                                 BPPARAM = NULL) {
    compute_mean_var <- getFromNamespace(".compute_mean_var", "scran")
    compute_blocked_stats_none <- getFromNamespace("compute_blocked_stats_none", "scran")
    compute_residual_stats_none <- getFromNamespace("compute_residual_stats_none", "scran")

    if (is.null(BPPARAM)) {
        BPPARAM <- getFromNamespace("SerialParam", "BiocParallel")()
    }

    mat <- SummarizedExperiment::assay(object, i = assay.type)
    dots <- list(...)

    compute_stats <- function(selected_rows) {
        do.call(
            compute_mean_var,
            c(
                list(
                    x = mat,
                    block = block,
                    design = design,
                    subset.row = selected_rows,
                    block.FUN = compute_blocked_stats_none,
                    residual.FUN = compute_residual_stats_none,
                    BPPARAM = BPPARAM
                ),
                dots
            )
        )
    }

    x.stats <- tryCatch(
        compute_stats(subset.row),
        error = function(e) {
            if (!grepl("no residual d\\.f\\.", conditionMessage(e))) {
                stop(e)
            }

            warning(
                "scran variance modelling failed for this dataset; falling back to simple logcounts variance ranking.",
                call. = FALSE
            )

            idx <- sclet_subset2index(subset.row, mat)
            submat <- mat[idx, , drop = FALSE]
            means <- MatrixGenerics::rowMeans(submat)
            vars <- MatrixGenerics::rowVars(submat)
            out <- S4Vectors::DataFrame(
                mean = means,
                total = vars,
                tech = NA_real_,
                bio = vars,
                p.value = NA_real_,
                FDR = NA_real_
            )
            rownames(out) <- rownames(mat)[idx]
            return(out)
        }
    )
    if (inherits(x.stats, "DFrame")) {
        return(x.stats)
    }
    if (is.null(subset.fit)) {
        fit.stats <- x.stats
    } else {
        fit.stats <- compute_stats(subset.fit)
    }

    collected <- do.call(
        sclet_decompose_log_exprs,
        c(
            list(
                x.means = x.stats$means,
                x.vars = x.stats$vars,
                fit.means = fit.stats$means,
                fit.vars = fit.stats$vars,
                ncells = x.stats$ncells
            ),
            dots
        )
    )

    output <- sclet_combine_blocked_statistics(
        collected,
        method = method,
        equiweight = equiweight,
        ncells = x.stats$ncells
    )

    rownames(output) <- rownames(mat)[sclet_subset2index(subset.row, mat)]
    output
}

#' scale data
#' 
#' @title ScaleData
#' @param object a SingleCellExperiment object
#' @param features selected features to be scaled, default is all features
#' @param assay selected assay to be scaled
#' @return an updated SingleCellExperiment object with 'scaled' assay
#' @export
#' @importFrom stats sd
ScaleData <- function(object, features = NULL, assay = "logcounts") {
    prev_state <- sclet_get_state(object)
    if (is.null(features)) {
        # all genes
        features <- rownames(object)
    } 

    # For DelayedArray, we must avoid densifying the matrix.
    # scater::logNormCounts already provides logcounts.
    # Standard scaling (subtract mean, divide sd) on sparse/HDF5 matrix results in a DENSE matrix.
    # This is bad for memory.
    # However, 'scaled' slot in Seurat is typically dense. 
    # If we want to support big data, we should use DelayedMatrix to represent the scaled matrix lazily.
    
    mat <- assay(object, assay)[features, , drop=FALSE]
    
    # Calculate row means and sds efficiently
    # DelayedMatrixStats or MatrixGenerics should be used
    if (is(mat, "DelayedMatrix")) {
        # DelayedArray logic
        rm <- MatrixGenerics::rowMeans(mat)
        rs <- MatrixGenerics::rowSds(mat)
    } else {
        # In-memory logic (could be sparse or dense)
        if (is(mat, "dgCMatrix")) {
             rm <- Matrix::rowMeans(mat)
             # Sparse matrix rowSds calculation can be tricky if not careful, 
             # but MatrixGenerics handles it or we use sparseMatrixStats if available.
             # Fallback to apply for now if not available, but let's assume MatrixGenerics works.
             rs <- MatrixGenerics::rowSds(mat)
        } else {
             rm <- rowMeans(mat)
             rs <- apply(mat, 1, stats::sd)
        }
    }
    
    rs[rs == 0] <- 0.01
    
    # Create a DelayedMatrix that represents the scaled data
    # (mat - rm) / rs
    # In DelayedArray, we can do arithmetic directly.
    
    if (is(mat, "DelayedMatrix")) {
        # This creates a DelayedMatrix (lazy evaluation)
        # We need to broadcast the subtraction and division.
        # DelayedArray supports standard broadcasting for vector-matrix ops if dimensions match?
        # Usually it matches on columns? No, R matrices match on columns (recycling).
        # mat is genes x cells. rm is genes.
        # (mat - rm) needs rm to be repeated for each cell.
        
        # DelayedArray arithmetic:
        # We can use sweep-like operations or just rely on broadcasting if implemented correctly.
        # Alternatively, create a DelayedMatrix backend.
        
        scaled <- (mat - rm) / rs
        
    } else {
        # For memory matrices, we might still want to return a DelayedMatrix if it's too big?
        # But if it's in memory, maybe just compute it.
        # NOTE: Seurat ScaleData forces dense. 
        # Here we try to keep it compatible but lazy if possible.
        
        # If the user explicitly wants HDF5/BigData, they should have used Read10X(backend="HDF5").
        # So here, if input is DelayedMatrix, output is DelayedMatrix.
        # If input is memory, output is memory (dense, as scaling destroys sparsity).
        
        scaled <- (mat - rm) / rs
    }
    
    if (length(features) == nrow(object)) {
        scaled_assay <- scaled
    } else {
        scaled_assay <- SummarizedExperiment::assay(object, assay)
        scaled_assay[features, ] <- scaled
    }
    SummarizedExperiment::assay(object, "scaled") <- scaled_assay
    object <- sclet_restore_state(object, prev_state)
    object <- sclet_set_layer(object, name = assay, assay = assay, role = assay, active = FALSE)
    object <- sclet_set_layer(
        object,
        name = "scaled",
        assay = "scaled",
        role = "scaled",
        params = list(
            source = assay,
            features = features
        )
    )
    object <- sclet_set_active_assay(object, "scaled")
    object <- sclet_log_command(
        object,
        "ScaleData",
        params = list(
            features = features,
            assay = assay
        )
    )
    return(object)
}
