#' Plot Cell-Cell Communication Bubble Chart
#'
#' @title sc_cci_bubble
#' @param object A SingleCellExperiment object.
#' @param sources A character vector of source groups to include. If `NULL`, all are included.
#' @param targets A character vector of target groups to include. If `NULL`, all are included.
#' @param id The communication record id to use.
#' @param pval_thresh The p-value threshold for filtering interactions.
#' @return A ggplot object
#' @importFrom ggplot2 ggplot aes geom_point scale_color_gradientn scale_size_continuous theme_minimal theme element_text labs
#' @importFrom rlang .data
#' @export
sc_cci_bubble <- function(object, sources = NULL, targets = NULL, id = NULL, pval_thresh = 0.05) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is needed for this function to work. Please install it.")
    }
    
    df <- get_cci(object, id = id)
    if (is.null(df)) {
        stop("No CCI data found in the object.")
    }
    
    # Filter by p-value
    if (!is.null(pval_thresh)) {
        df <- df[df$p_value <= pval_thresh, , drop = FALSE]
    }
    
    if (!is.null(sources)) {
        df <- df[df$source_group %in% sources, , drop = FALSE]
    }
    
    if (!is.null(targets)) {
        df <- df[df$target_group %in% targets, , drop = FALSE]
    }
    
    if (nrow(df) == 0) {
        stop("No interactions left after filtering.")
    }
    
    df$interaction <- paste(df$ligand, df$receptor, sep = " - ")
    df$cell_pair <- paste(df$source_group, df$target_group, sep = " -> ")
    
    # -log10(p_value) for size or color
    # Add a small epsilon to avoid -log10(0)
    df$mlog10p <- -log10(df$p_value + 1e-10)
    
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$cell_pair, y = .data$interaction)) +
        ggplot2::geom_point(ggplot2::aes(size = .data$mlog10p, color = .data$probability)) +
        ggplot2::scale_color_gradientn(colors = c("blue", "white", "red")) +
        ggplot2::scale_size_continuous(range = c(1, 6)) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
        ggplot2::labs(x = "Cell Pair (Source -> Target)", y = "Interaction (Ligand - Receptor)", 
                      color = "Probability", size = "-log10(P-value)")
    
    return(p)
}

#' Plot Cell-Cell Communication Chord Diagram
#'
#' @title sc_cci_chord
#' @param object A SingleCellExperiment object.
#' @param signaling A character vector of pathways or interactions to include.
#' @param id The communication record id to use.
#' @param pval_thresh The p-value threshold for filtering interactions.
#' @param ... Additional arguments passed to `circlize::chordDiagram`.
#' @return NULL (Plots the diagram)
#' @export
sc_cci_chord <- function(object, signaling = NULL, id = NULL, pval_thresh = 0.05, ...) {
    if (!requireNamespace("circlize", quietly = TRUE)) {
        stop("Package 'circlize' is needed for this function to work. Please install it.")
    }
    
    df <- get_cci(object, id = id)
    if (is.null(df)) {
        stop("No CCI data found in the object.")
    }
    
    if (!is.null(pval_thresh)) {
        df <- df[df$p_value <= pval_thresh, , drop = FALSE]
    }
    
    if (!is.null(signaling)) {
        df <- df[df$pathway %in% signaling | paste(df$ligand, df$receptor, sep=" - ") %in% signaling, , drop = FALSE]
    }
    
    if (nrow(df) == 0) {
        stop("No interactions left after filtering.")
    }
    
    # Aggregate probabilities between cell groups
    agg_df <- stats::aggregate(probability ~ source_group + target_group, data = df, sum)
    
    circlize::chordDiagram(agg_df, ...)
}
