#' Refit Milo DA after dropping separated neighbourhoods
#'
#' @param milo A Milo object containing neighbourhood counts and metadata.
#' @param design_df Data.frame of experimental design information (samples x covariates).
#' @param bad_hood Integer vector of indices of nhoods to drop.
#' @param formula Model formula used in \code{testNhoods}.
#' @param contrasts Character string or vector of model contrasts.
#' @param glmm Logical; if TRUE fit a NB-GLMM, otherwise fit a NB-GLM.
#' @param glmm.solver Solver for GLMM variance components.
#' @param REML Logical; use restricted maximum likelihood for GLMM.
#' @param max.iters Maximum number of iterations for GLMM fitting.
#' @param fail.on.error Logical; continue on errors if FALSE.
#' @param fdr.weighting Method for FDR weighting, passed to \code{testNhoods}.
#' @param norm.method Normalisation method.
#' @param BPPARAM BiocParallel parameter.
#' @param old_results Optional; previous DA results for comparison of significant nhood counts.
#' @param seed Random seed for reproducibility.
#' @param ... Additional arguments passed to \code{testNhoods}.
#'
#' @return A list with:
#' \item{milo_clean}{Milo object with separated nhoods removed.}
#' \item{da_clean}{Data.frame of refitted DA results.}
#' \item{dropped}{Indices of dropped nhoods.}
#'
#' @export
refit_milo <- function(
    milo,
    design_df,
    bad_hood,
    formula,
    contrasts,
    glmm = TRUE,
    glmm.solver = "Fisher",
    REML = TRUE,
    max.iters = 50,
    fail.on.error = FALSE,
    fdr.weighting = "graph-overlap",
    norm.method = "TMM",
    BPPARAM = NULL,
    old_results = NULL,
    seed = 2025,
    ...
) {
    if (!requireNamespace("miloR", quietly = TRUE)) {
        stop("Package 'miloR' is needed for this function to work. Please install it.")
    }

    set.seed(seed)

    milo2 <- if (length(bad_hood) > 0) milo[-bad_hood, ] else milo

    test_args <- list(
        x = milo2,
        design = stats::as.formula(formula),
        design.df = design_df,
        model.contrasts = contrasts,
        fdr.weighting = fdr.weighting,
        norm.method = norm.method,
        ...
    )
    if (!is.null(BPPARAM)) test_args$BPPARAM <- BPPARAM

    if (glmm) {
        test_args$glmm.solver <- glmm.solver
        test_args$REML <- REML
        test_args$max.iters <- max.iters
        test_args$fail.on.error <- fail.on.error
        test_args$force <- TRUE
    }
    da_new <- do.call(miloR::testNhoods, test_args)

    if (!is.null(old_results) && "SpatialFDR" %in% names(old_results)) {
        before <- sum(old_results$SpatialFDR < 0.1, na.rm = TRUE)
        after <- sum(da_new$SpatialFDR < 0.1, na.rm = TRUE)
        message(sprintf("Significant nhoods (FDR<0.1): before=%s, after=%s", before, after))
    }

    return(list(milo_clean = milo2, da_clean = da_new, dropped = bad_hood))
}
