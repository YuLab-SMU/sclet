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

#' Plot velocity latent time on an embedding
#'
#' Colors cells on a dimensional reduction by velocity pseudotime (latent time),
#' showing the predicted temporal ordering inferred by scVelo.
#'
#' @param object A SingleCellExperiment object with velocity results.
#' @param reduction Character. Dimensional reduction to use. If `NULL`, uses
#'   `DefaultReduction(object)`.
#' @param id Optional velocity record id.
#' @param point_size Numeric. Point size. Defaults to 0.6.
#' @param low,high Colors for the pseudotime gradient.
#'
#' @return A ggplot object.
#' @export
plot_velocity_latent_time <- function(object, reduction = NULL, id = NULL,
                                      point_size = 0.6, low = "grey90", high = "firebrick") {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is needed for this function to work. Please install it.")
    }

    vel <- get_velocity(object, id = id)
    if (is.null(vel)) {
        stop("No velocity results found. Please run RunVelocity() first.")
    }

    if (is.null(reduction)) {
        reduction <- vel$artifacts$reduction
    }
    if (is.null(reduction)) {
        reduction <- vel$reduction
    }
    if (is.null(reduction)) {
        reduction <- DefaultReduction(object)
    }
    if (is.null(reduction) || !reduction %in% SingleCellExperiment::reducedDimNames(object)) {
        stop("A valid dimensional reduction is required.")
    }

    object_coldata <- colnames(SummarizedExperiment::colData(object))
    velocity_coldata <- vel$artifacts$colData
    if (is.null(velocity_coldata)) {
        velocity_coldata <- vel$artifacts$velocity_coldata
    }
    if (is.null(velocity_coldata)) {
        velocity_coldata <- vel$velocity_coldata
    }
    candidate_cols <- unique(c(
        velocity_coldata[grepl("latent_time|pseudotime", velocity_coldata)],
        "velocity_pseudotime",
        "latent_time_regvelo",
        "regvelo_latent_time"
    ))
    candidate_cols <- candidate_cols[candidate_cols %in% object_coldata]
    if (!length(candidate_cols)) {
        stop("No velocity latent-time or pseudotime column found in colData.")
    }
    pseudotime_col <- candidate_cols[[1]]

    emb <- SingleCellExperiment::reducedDim(object, reduction)
    df <- data.frame(
        dim1 = emb[, 1],
        dim2 = emb[, 2],
        pseudotime = as.numeric(SummarizedExperiment::colData(object)[[pseudotime_col]])
    )

    ggplot2::ggplot(df, ggplot2::aes(x = .data$dim1, y = .data$dim2, color = .data$pseudotime)) +
        ggplot2::geom_point(size = point_size) +
        ggplot2::scale_color_gradient(low = low, high = high) +
        ggplot2::labs(
            x = paste0(reduction, " 1"),
            y = paste0(reduction, " 2"),
            color = "Latent\ntime",
            title = "Velocity Latent Time"
        ) +
        ggplot2::theme_classic() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold", size = 13)
        )
}

sclet_velocity_assay_name <- function(object, id = NULL, assay = NULL) {
    if (!is.null(assay)) {
        if (!assay %in% SummarizedExperiment::assayNames(object)) {
            stop("Velocity assay '", assay, "' not found in object.")
        }
        return(assay)
    }

    vel <- get_velocity(object, id = id)
    if (is.null(vel)) {
        stop("No velocity results found. Please run RunVelocity() or RunRegVelo() first.")
    }

    assay <- vel$artifacts$velocity_assay
    if (is.null(assay) && !is.null(vel$id)) {
        candidate <- paste0(vel$id, "_velocity")
        if (candidate %in% SummarizedExperiment::assayNames(object)) {
            assay <- candidate
        }
    }
    if (is.null(assay) || !assay %in% SummarizedExperiment::assayNames(object)) {
        label <- vel$id
        if (is.null(label)) {
            label <- id
        }
        if (is.null(label)) {
            label <- "<active>"
        }
        stop("No velocity assay recorded for velocity result '", label, "'.")
    }
    assay
}

#' Compute Velocity Magnitude
#'
#' Computes the Euclidean norm of a stored velocity matrix per cell or per gene.
#' For RegVelo results this uses the assay recorded in
#' `get_velocity(object, id)$artifacts$velocity_assay`.
#'
#' @param object A SingleCellExperiment object with velocity results.
#' @param id Optional velocity record id. If `NULL`, uses the active velocity
#'   result.
#' @param assay Optional velocity assay name. If supplied, overrides `id`.
#' @param margin One of `"cell"` or `"gene"`.
#' @param name Optional column name used when `store = TRUE`.
#' @param store Logical. If `TRUE`, stores cell magnitudes in `colData` or gene
#'   magnitudes in `rowData` and returns the updated object.
#'
#' @return A numeric vector, or an updated SingleCellExperiment when
#'   `store = TRUE`.
#' @export
VelocityMagnitude <- function(object, id = NULL, assay = NULL,
                              margin = c("cell", "gene"), name = NULL,
                              store = FALSE) {
    margin <- match.arg(margin)
    assay <- sclet_velocity_assay_name(object, id = id, assay = assay)
    velocity <- SummarizedExperiment::assay(object, assay)
    velocity <- Matrix::Matrix(velocity, sparse = TRUE)

    values <- switch(
        margin,
        cell = sqrt(Matrix::colSums(velocity ^ 2)),
        gene = sqrt(Matrix::rowSums(velocity ^ 2))
    )
    values <- as.numeric(values)
    names(values) <- if (identical(margin, "cell")) colnames(object) else rownames(object)

    if (!isTRUE(store)) {
        return(values)
    }

    if (is.null(name)) {
        name <- paste0(assay, "_magnitude")
    }
    if (identical(margin, "cell")) {
        SummarizedExperiment::colData(object)[[name]] <- values
    } else {
        SummarizedExperiment::rowData(object)[[name]] <- values
    }
    object
}

#' Plot Velocity Magnitude on an Embedding
#'
#' Colors cells on a dimensional reduction by the per-cell magnitude of a stored
#' velocity matrix.
#'
#' @param object A SingleCellExperiment object with velocity results.
#' @param reduction Character. Dimensional reduction to use. If `NULL`, uses
#'   the reduction recorded by the velocity result, then `DefaultReduction()`.
#' @param id Optional velocity record id.
#' @param assay Optional velocity assay name. If supplied, overrides `id`.
#' @param point_size Numeric point size.
#' @param low,high Colors for the magnitude gradient.
#'
#' @return A ggplot object.
#' @export
plot_velocity_magnitude <- function(object, reduction = NULL, id = NULL,
                                    assay = NULL, point_size = 0.6,
                                    low = "grey90", high = "navy") {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is needed for this function to work. Please install it.")
    }

    vel <- get_velocity(object, id = id)
    if (is.null(reduction) && !is.null(vel)) {
        reduction <- vel$inputs$reduction
    }
    if (is.null(reduction)) {
        reduction <- DefaultReduction(object)
    }
    if (is.null(reduction) || !reduction %in% SingleCellExperiment::reducedDimNames(object)) {
        stop("A valid dimensional reduction is required.")
    }

    assay <- sclet_velocity_assay_name(object, id = id, assay = assay)
    magnitude <- VelocityMagnitude(object, assay = assay, margin = "cell")
    emb <- SingleCellExperiment::reducedDim(object, reduction)
    df <- data.frame(
        dim1 = emb[, 1],
        dim2 = emb[, 2],
        magnitude = magnitude
    )

    ggplot2::ggplot(df, ggplot2::aes(x = .data$dim1, y = .data$dim2, color = .data$magnitude)) +
        ggplot2::geom_point(size = point_size) +
        ggplot2::scale_color_gradient(low = low, high = high) +
        ggplot2::labs(
            x = paste0(reduction, " 1"),
            y = paste0(reduction, " 2"),
            color = "Velocity\nmagnitude",
            title = "Velocity Magnitude"
        ) +
        ggplot2::theme_classic() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold", size = 13)
        )
}

#' Top Velocity Magnitude Genes
#'
#' Returns genes ranked by the magnitude of their velocity vector across cells.
#'
#' @param object A SingleCellExperiment object with velocity results.
#' @param id Optional velocity record id.
#' @param assay Optional velocity assay name. If supplied, overrides `id`.
#' @param n Number of genes to return.
#' @param decreasing Logical. If `TRUE`, returns the largest magnitudes.
#'
#' @return A data.frame with `gene`, `velocity_magnitude`, and `rank`.
#' @export
TopVelocityGenes <- function(object, id = NULL, assay = NULL, n = 20,
                             decreasing = TRUE) {
    magnitude <- VelocityMagnitude(object, id = id, assay = assay, margin = "gene")
    ord <- order(magnitude, decreasing = decreasing, na.last = NA)
    ord <- ord[seq_len(min(as.integer(n), length(ord)))]
    data.frame(
        gene = names(magnitude)[ord],
        velocity_magnitude = as.numeric(magnitude[ord]),
        rank = seq_along(ord),
        stringsAsFactors = FALSE
    )
}

#' Plot Top Velocity Magnitude Genes
#'
#' Shows genes with the largest velocity magnitude across cells.
#'
#' @param object A SingleCellExperiment object with velocity results.
#' @param id Optional velocity record id.
#' @param assay Optional velocity assay name. If supplied, overrides `id`.
#' @param n Number of genes to show.
#' @param fill Bar fill color.
#'
#' @return A ggplot object.
#' @export
plot_top_velocity_genes <- function(object, id = NULL, assay = NULL, n = 20,
                                    fill = "steelblue") {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is needed for this function to work. Please install it.")
    }

    df <- TopVelocityGenes(object, id = id, assay = assay, n = n)
    df$gene <- factor(df$gene, levels = rev(df$gene))
    ggplot2::ggplot(df, ggplot2::aes(x = .data$gene, y = .data$velocity_magnitude)) +
        ggplot2::geom_col(fill = fill, width = 0.75) +
        ggplot2::coord_flip() +
        ggplot2::labs(
            x = NULL,
            y = "Velocity magnitude",
            title = "Top Velocity Genes"
        ) +
        ggplot2::theme_classic() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold", size = 13)
        )
}
