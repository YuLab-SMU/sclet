#' Run Milo Differential Abundance Analysis
#'
#' @param sce A SingleCellExperiment object containing single-cell data. Must have PCA in \code{reducedDims}.
#' @param sample_col The column name in \code{colData(sce)} that indicates the sample label.
#' @param condition_col The column name in \code{colData(sce)} that indicates the experimental condition to test.
#' @param covariates Optional character vector of additional columns as covariates in \code{colData(sce)}.
#' @param contrasts Character vector of contrast expressions.
#' @param k1 Number of neighbors used in \code{buildGraph()}, default = 10.
#' @param k2 Number of neighbors used in \code{makeNhoods()}, default = 10.
#' @param d1 Number of PCA dimensions used in \code{buildGraph()}, default = 30.
#' @param d2 Number of PCA dimensions used in \code{makeNhoods()}, default = 30.
#' @param d3 Number of PCA dimensions used in \code{calcNhoodDistance()}, default = 30.
#' @param prop Proportion of cells to sample as neighborhood centers, default = 0.1.
#' @param overlap Integer, minimum percentage of overlap between neighbourhoods to draw an edge. Default = 1.
#' @param per_contrast Logical; test contrasts separately (\code{TRUE}) or jointly (\code{FALSE}).
#' @param is_refined Logical; whether to apply refinement algorithm when sampling neighborhoods.
#' @param refinement_scheme Refinement scheme passed to \code{makeNhoods}.
#' @param adjustment Method for multiple testing correction within neighborhoods.
#' @param glmm Logical; use NB-GLMM with a random intercept.
#' @param random_effect Character; column in colData(sce) for the random intercept (GLMM only).
#' @param glmm_solver Solver for GLMM variance components: "Fisher", "HE", or "HE-NNLS".
#' @param REML Logical; restricted ML for GLMM variance components (GLMM only).
#' @param max.iters Integer; maximum iterations for GLMM fitting.
#' @param fail.on.error Logical; GLMM continues on errors if FALSE.
#' @param BPPARAM BiocParallel param for fitting.
#' @param seed Random seed for reproducibility.
#' @param ... Additional arguments passed to underlying miloR functions.
#'
#' @return A list with:
#' \item{sce}{Input \code{SingleCellExperiment} with DA results stored in metadata.}
#' \item{milo}{Milo object with graph, neighborhoods, and DA results.}
#' \item{da_results}{Data.frame of DA testing results.}
#'
#' @export
runMilo <- function(
    sce,
    sample_col,
    condition_col,
    covariates = NULL,
    contrasts = NULL,
    k1 = 10,
    k2 = 10,
    d1 = 30,
    d2 = 30,
    d3 = 30,
    prop = 0.1,
    overlap = 1,
    per_contrast = TRUE,
    is_refined = TRUE,
    refinement_scheme = "graph",
    adjustment = "BH",
    glmm = FALSE,
    random_effect = NULL,
    glmm_solver = "Fisher",
    REML = TRUE,
    max.iters = 50,
    fail.on.error = FALSE,
    BPPARAM = NULL,
    seed = 2025,
    ...
) {
    if (!requireNamespace("miloR", quietly = TRUE)) {
        stop("Package 'miloR' is needed for this function to work. Please install it.")
    }

    if (!is(sce, "SingleCellExperiment")) {
        stop("The input must be a SingleCellExperiment object.")
    }

    set.seed(seed)

    if (!("PCA" %in% reducedDimNames(sce))) stop("PCA not found in reducedDims(sce).")
    if (!(sample_col %in% colnames(colData(sce)))) stop("`sample_col` not found in colData(sce).")
    if (!(condition_col %in% colnames(colData(sce)))) stop("`condition_col` not found in colData(sce).")

    if (!is.null(covariates) && length(covariates) > 0) {
        miss <- setdiff(covariates, colnames(colData(sce)))
        if (length(miss) > 0) stop("`covariates` not found in colData(sce): ", paste(miss, collapse = ", "))
    }

    if (is.null(contrasts) || length(contrasts) == 0) {
        stop("A valid contrasts vector should be supplied, e.g. 'ConditionB - ConditionA'.")
    }

    if (glmm) {
        if (is.null(random_effect)) stop("For GLMM, `random_effect` must be provided.")
        if (!(random_effect %in% colnames(colData(sce)))) stop("`random_effect` not found in colData(sce).")
    }

    milo <- miloR::Milo(sce)
    if ("UMAP" %in% reducedDimNames(sce)) {
        reducedDim(milo, "UMAP") <- reducedDim(sce, "UMAP")
    }

    milo <- miloR::buildGraph(milo, k = k1, d = d1, ...)
    milo <- miloR::makeNhoods(
        milo,
        prop = prop,
        k = k2,
        d = d2,
        refined = is_refined,
        refinement_scheme = refinement_scheme,
        ...
    )

    milo <- miloR::countCells(milo, meta.data = data.frame(colData(milo)), samples = sample_col, ...)
    milo <- miloR::calcNhoodDistance(milo, d = d3, ...)
    milo <- miloR::buildNhoodGraph(milo, overlap = overlap, ...)

    need_cols <- unique(c(sample_col, condition_col, covariates, if (glmm) random_effect))
    design_df <- data.frame(colData(milo))[, need_cols, drop = FALSE]

    by_sample <- split(design_df, design_df[[sample_col]])
    var_ok <- vapply(by_sample, function(x) {
        cols <- setdiff(colnames(x), sample_col)
        all(vapply(cols, function(colnm) {
            length(unique(x[[colnm]])) == 1
        }, logical(1)))
    }, logical(1))
    if (!all(var_ok)) {
        bad <- names(var_ok)[!var_ok]
        stop("Non-unique condition/covariates within sample: ", paste(bad, collapse = ", "))
    }

    design_df <- design_df[!duplicated(design_df[[sample_col]]), , drop = FALSE]
    rownames(design_df) <- design_df[[sample_col]]
    design_df <- design_df[colnames(miloR::nhoodCounts(milo)), , drop = FALSE]

    design_df[[condition_col]] <- droplevels(factor(design_df[[condition_col]]))
    if (glmm) design_df[[random_effect]] <- droplevels(factor(design_df[[random_effect]]))

    if (!is.null(covariates) && length(covariates) > 0) {
        for (cv in covariates) {
            if (is.character(design_df[[cv]]) || is.logical(design_df[[cv]])) {
                design_df[[cv]] <- factor(design_df[[cv]])
            }
        }
    }

    if (!glmm) {
        rhs <- c(paste0("0 + ", condition_col), covariates)
        fml <- stats::as.formula(paste("~", paste(rhs, collapse = " + ")))
    } else {
        rhs <- c(condition_col, covariates, paste0("(1|", random_effect, ")"))
        fml <- stats::as.formula(paste("~", paste(rhs, collapse = " + ")))
    }

    test_args <- list(
        milo,
        design = fml,
        design.df = design_df,
        fdr.weighting = "graph-overlap",
        norm.method = "TMM",
        ...
    )
    if (!is.null(BPPARAM)) test_args$BPPARAM <- BPPARAM
    if (glmm) {
        test_args$glmm.solver <- glmm_solver
        test_args$REML <- REML
        test_args$max.iters <- max.iters
        test_args$fail.on.error <- fail.on.error
        test_args$force <- TRUE # Override the small N warning
    }

    if (!per_contrast) {
        test_args$model.contrasts <- contrasts
        da_results <- do.call(miloR::testNhoods, test_args)
    } else {
        res_list <- lapply(contrasts, function(ct) {
            test_args$model.contrasts <- ct
            r <- do.call(miloR::testNhoods, test_args)
            r$contrast <- ct
            r
        })
        da_results <- do.call(rbind, res_list)

        if ("PValue" %in% colnames(da_results)) {
            da_results$PValue_adj_withinNhood <- ave(
                da_results$PValue,
                da_results$Nhood,
                FUN = function(x) stats::p.adjust(x, method = adjustment)
            )
        }
    }

    sce@metadata$milo_results <- list(
        da_results = da_results,
        design_df = design_df,
        formula = deparse(fml),
        contrasts = contrasts,
        params = list(
            k1 = k1,
            k2 = k2,
            d1 = d1,
            d2 = d2,
            d3 = d3,
            prop = prop,
            overlap = overlap,
            is_refined = is_refined,
            seed = seed,
            glmm = glmm,
            random_effect = random_effect,
            glmm_solver = glmm_solver,
            REML = REML,
            max.iters = max.iters,
            fail.on.error = fail.on.error
        )
    )

    return(list(
        sce = sce,
        milo = milo,
        da_results = da_results
    ))
}
