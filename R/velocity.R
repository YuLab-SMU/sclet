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
#' @param name velocity record id. Defaults to `"velocity"`.
#' @param ... Additional arguments passed to `velociraptor::scvelo`.
#' 
#' @return A SingleCellExperiment object with velocity state updated.
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom SummarizedExperiment colData colData<-
#' @export
RunVelocity <- function(sce, mode = c("deterministic", "stochastic", "dynamical"), use.dimred = NULL,
                        subset.row = NULL, subset_row = NULL, name = "velocity", ...) {
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
    
    prev_state <- sclet_get_state(sce)

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
    sce <- sclet_restore_state(sce, prev_state)
    
    # vel_res is a SingleCellExperiment containing the velocity outputs
    # We want to extract the velocity results and attach them to the original sce
    # velociraptor typically adds 'velocity' and 'velocity_pseudotime' etc.
    
    velocity_reductions <- grep("^velocity_", SingleCellExperiment::reducedDimNames(vel_res), value = TRUE)
    for (rd_name in velocity_reductions) {
        SingleCellExperiment::reducedDim(sce, rd_name) <- SingleCellExperiment::reducedDim(vel_res, rd_name)
    }

    velocity_coldata <- intersect(
        c("velocity_pseudotime", "velocity_confidence", "root_cells", "end_points"),
        colnames(SummarizedExperiment::colData(vel_res))
    )
    for (field in velocity_coldata) {
        SummarizedExperiment::colData(sce)[[field]] <- SummarizedExperiment::colData(vel_res)[[field]]
    }

    velocity_analysis <- list(
        id = name,
        mode = mode,
        use.dimred = use.dimred,
        subset.row = subset.row,
        results = vel_res,
        reductions = velocity_reductions,
        coldata_fields = velocity_coldata,
        timestamp = Sys.time()
    )
    sce <- sclet_set_analysis(sce, "velocity", velocity_analysis)
    sce <- sclet_set_state_record(
        object = sce,
        type = "velocity",
        id = name,
        active = TRUE,
        value = list(
            method = "velociraptor::scvelo",
            inputs = list(
                assays = c("spliced", "unspliced"),
                reduction = use.dimred
            ),
            artifacts = list(
                analysis_key = name,
                colData = velocity_coldata,
                reductions = velocity_reductions
            ),
            params = list(
                mode = mode,
                subset.row = subset.row
            ),
            summary = list(
                n_features = if (is.null(subset.row)) nrow(sce) else length(subset.row),
                n_cells = ncol(sce)
            ),
            results = vel_res,
            reductions = velocity_reductions,
            coldata_fields = velocity_coldata,
            created_at = Sys.time()
        )
    )
    sce <- sclet_log_command(
        sce,
        "RunVelocity",
        params = list(
            mode = mode,
            use.dimred = use.dimred,
            subset.row = subset.row,
            name = name
        ),
        outputs = list(
            analysis = "velocity",
            velocity = name
        )
    )

    return(sce)
}

#' Plot RNA Velocity
#'
#' Plots RNA velocity stream or arrows on a reduced dimension plot.
#'
#' @param sce A SingleCellExperiment object with velocity calculated.
#' @param reduction String specifying the reduction to plot on. If `NULL`,
#'   prefers the reduction recorded in the velocity analysis and otherwise
#'   falls back to the active reduction.
#' @param group.by String specifying the colData column for cell grouping. If NULL,
#'   uses `ActiveIdent(sce)`, including the special `"colLabels"` identity source.
#' @param id Optional velocity record id.
#' @param ... Additional arguments passed to plotting function.
#' 
#' @return A ggplot object.
#' @importFrom rlang .data
#' @importFrom S4Vectors metadata
#' @importFrom SummarizedExperiment colData
#' @export
VelocityPlot <- function(sce, reduction = NULL, group.by = NULL, id = NULL, ...) {
    if (!requireNamespace("velociraptor", quietly = TRUE)) {
        stop("Please install 'velociraptor' to plot RNA Velocity.")
    }
    
    vel_state <- get_velocity(sce, id = id)
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
    if (is.null(reduction)) {
        reduction <- vel_state$inputs$reduction
    }
    if (is.null(reduction)) {
        reduction <- DefaultReduction(sce)
    }
    if (is.null(reduction) || !reduction %in% reducedDimNames(sce)) {
        stop(sprintf("Reduction %s not found in sce.", reduction))
    }
    
    # Wait, embedVelocity needs the original sce and vel_res
    emb <- reducedDim(sce, reduction)
    
    # Actually velociraptor provides `embedVelocity`
    # Let's extract velocity vectors for the specific reduction
    vel_emb <- velociraptor::embedVelocity(emb, vel_res)
    grid_df <- velociraptor::gridVectors(emb, vel_emb, as.data.frame = FALSE)
    if (is.null(grid_df$start) || is.null(grid_df$end)) {
        stop("`velociraptor::gridVectors()` did not return the expected `start`/`end` components.")
    }
    grid_df <- data.frame(
        x = grid_df$start[, 1],
        y = grid_df$start[, 2],
        xend = grid_df$end[, 1],
        yend = grid_df$end[, 2]
    )

    if (is.null(group.by)) {
        group.by <- ActiveIdent(sce)
    }
    if (is.null(group.by)) {
        stop("No active identity found. Please provide 'group.by'.")
    }

    if (identical(group.by, "colLabels")) {
        groups <- SingleCellExperiment::colLabels(sce)
        if (is.null(groups)) {
            stop("Active identity is 'colLabels' but `colLabels(sce)` is empty.")
        }
    } else {
        meta <- as.data.frame(SummarizedExperiment::colData(sce))
        if (!group.by %in% colnames(meta)) {
            stop("Column '", group.by, "' not found in colData(sce).")
        }
        groups <- meta[[group.by]]
    }
    if (length(groups) != nrow(emb)) {
        stop(
            "Grouping vector derived from '", group.by, "' has length ", length(groups),
            ", but reduction '", reduction, "' has ", nrow(emb), " cells."
        )
    }
    
    # Prepare base plot
    df <- data.frame(
        x = emb[, 1],
        y = emb[, 2],
        group = groups
    )
    
    p <- ggplot2::ggplot() +
        ggplot2::geom_point(data = df, ggplot2::aes(x = .data$x, y = .data$y, color = .data$group), size = 0.5, alpha = 0.8) +
        ggplot2::geom_segment(data = grid_df, 
                              ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
                              arrow = ggplot2::arrow(length = ggplot2::unit(0.05, "inches")),
                              linewidth = 0.3, color = "black") +
        ggplot2::theme_classic() +
        ggplot2::labs(x = paste0(reduction, " 1"), y = paste0(reduction, " 2"), color = group.by)
    
    return(p)
}

#' @name VelocityPlot
#' @aliases plot_velocity
#' @export
plot_velocity <- VelocityPlot
