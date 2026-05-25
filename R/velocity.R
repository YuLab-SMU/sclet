#' Run RNA Velocity Analysis
#'
#' This function runs RNA velocity analysis using `velociraptor::scvelo`.
#' It requires the input SingleCellExperiment to have "spliced" and "unspliced" assays.
#' The velocity results will be stored in the analysis-state contract under `velocity`.
#'
#' @param sce A SingleCellExperiment object.
#' @param mode String specifying the mode for velocity estimation. Can be "deterministic", "stochastic", or "dynamical".
#' @param use.dimred String specifying the dimensionality reduction to use. If NULL, defaults to `DefaultReduction(sce)`.
#' @param subset.row Optional row subset passed to `velociraptor::scvelo()`.
#' @param subset_row Legacy alias of `subset.row`, kept for compatibility.
#' @param ... Additional arguments passed to `velociraptor::scvelo`.
#' 
#' @return A SingleCellExperiment object with velocity state updated.
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom SummarizedExperiment colData colData<-
#' @export
RunVelocity <- function(sce, mode = c("deterministic", "stochastic", "dynamical"), use.dimred = NULL,
                        subset.row = NULL, subset_row = NULL, ...) {
    if (!requireNamespace("velociraptor", quietly = TRUE)) {
        stop("Please install 'velociraptor' to run RNA Velocity.")
    }
    
    mode <- match.arg(mode)
    
    if (is.null(use.dimred)) {
        use.dimred <- DefaultReduction(sce)
    }
    
    if (is.null(use.dimred) || !use.dimred %in% reducedDimNames(sce)) {
        stop("Valid dimensionality reduction is required for velocity projection.")
    }
    
    if (!all(c("spliced", "unspliced") %in% assayNames(sce))) {
        stop("Assays 'spliced' and 'unspliced' are required for RNA velocity.")
    }

    if (!is.null(subset.row) && !is.null(subset_row)) {
        stop("Please supply only one of 'subset.row' or 'subset_row'.")
    }
    if (is.null(subset.row) && !is.null(subset_row)) {
        subset.row <- subset_row
    }
    
    # Run scvelo via velociraptor
    # Note: velociraptor uses its own basilisk environment
    message("Running scVelo using velociraptor...")
    vel_res <- velociraptor::scvelo(
        sce,
        mode = mode,
        use.dimred = use.dimred,
        subset.row = subset.row,
        ...
    )
    
    # vel_res is a SingleCellExperiment containing the velocity outputs
    # We want to extract the velocity results and attach them to the original sce
    # velociraptor typically adds 'velocity' and 'velocity_pseudotime' etc.
    
    # Store velocity reduction
    if ("velocity_umap" %in% reducedDimNames(vel_res)) {
        reducedDim(sce, "velocity_umap") <- reducedDim(vel_res, "velocity_umap")
    }
    
    # Register state
    metadata(sce)$sclet$analyses$velocity <- list(
        mode = mode,
        use.dimred = use.dimred,
        subset.row = subset.row,
        timestamp = Sys.time()
    )
    
    # Mark as active trajectory/velocity
    metadata(sce)$sclet$active$velocity <- "velocity"
    
    # Move velocity pseudotime if available
    if ("velocity_pseudotime" %in% colnames(colData(vel_res))) {
        sce$velocity_pseudotime <- vel_res$velocity_pseudotime
    }
    
    # We also keep the full vel_res in analyses for advanced access
    metadata(sce)$sclet$analyses$velocity$results <- vel_res
    
    return(sce)
}

#' Plot RNA Velocity
#'
#' Plots RNA velocity stream or arrows on a reduced dimension plot.
#'
#' @param sce A SingleCellExperiment object with velocity calculated.
#' @param reduction String specifying the reduction to plot on. Defaults to "UMAP".
#' @param group.by String specifying the colData column for cell grouping.
#' @param ... Additional arguments passed to plotting function.
#' 
#' @return A ggplot object.
#' @importFrom rlang .data
#' @importFrom S4Vectors metadata
#' @importFrom SummarizedExperiment colData
#' @export
VelocityPlot <- function(sce, reduction = "UMAP", group.by = ActiveIdent(sce), ...) {
    if (!requireNamespace("velociraptor", quietly = TRUE)) {
        stop("Please install 'velociraptor' to plot RNA Velocity.")
    }
    
    vel_state <- S4Vectors::metadata(sce)$sclet$analyses$velocity
    if (is.null(vel_state) || is.null(vel_state$results)) {
        stop("No velocity results found. Please run RunVelocity() first.")
    }
    
    vel_res <- vel_state$results
    
    # Currently velociraptor does not have a built-in ggplot grid, it uses base R
    # or scater. But there are custom wrappers.
    # We can use velociraptor::plotVelocity if available, but it doesn't exist?
    # Actually velociraptor::scvelo returns SCE with velocity_umap. We can use
    # scater::plotReducedDim with grid/stream. 
    # velociraptor has gridVectors and embedVelocity.
    
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 is required.")
    }
    
    # Basic implementation using velociraptor::gridVectors
    if (!reduction %in% reducedDimNames(sce)) {
        stop(sprintf("Reduction %s not found in sce.", reduction))
    }
    
    # Wait, embedVelocity needs the original sce and vel_res
    emb <- reducedDim(sce, reduction)
    
    # Actually velociraptor provides `embedVelocity`
    # Let's extract velocity vectors for the specific reduction
    vel_emb <- velociraptor::embedVelocity(emb, vel_res)
    grid_df <- velociraptor::gridVectors(emb, vel_emb)
    
    # Prepare base plot
    df <- data.frame(
        x = emb[, 1],
        y = emb[, 2],
        group = colData(sce)[[group.by]]
    )
    
    p <- ggplot2::ggplot() +
        ggplot2::geom_point(data = df, ggplot2::aes(x = .data$x, y = .data$y, color = .data$group), size = 0.5, alpha = 0.8) +
        ggplot2::geom_segment(data = data.frame(grid_df), 
                              ggplot2::aes(x = .data$start.1, y = .data$start.2, xend = .data$end.1, yend = .data$end.2),
                              arrow = ggplot2::arrow(length = ggplot2::unit(0.05, "inches")),
                              linewidth = 0.3, color = "black") +
        ggplot2::theme_classic() +
        ggplot2::labs(x = paste0(reduction, " 1"), y = paste0(reduction, " 2"), color = group.by)
    
    return(p)
}
