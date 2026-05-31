#' Plot Cell-Cell Communication Network
#'
#' @title sc_cci_network
#' @param object A SingleCellExperiment object.
#' @param signaling A character vector of pathways or interactions to include.
#' @param id The communication record id to use.
#' @param pval_thresh The p-value threshold for filtering interactions.
#' @param layout The layout algorithm to use for the network (default: "circle").
#' @param edge_width_scale A scaling factor for edge width.
#' @return A ggplot object
#' @export
sc_cci_network <- function(object, signaling = NULL, id = NULL, pval_thresh = 0.05, layout = "circle", edge_width_scale = 1) {
    if (!requireNamespace("igraph", quietly = TRUE)) {
        stop("Package 'igraph' is needed for this function to work. Please install it.")
    }
    if (!requireNamespace("ggraph", quietly = TRUE)) {
        stop("Package 'ggraph' is needed for this function to work. Please install it.")
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
    
    # Create igraph object
    g <- igraph::graph_from_data_frame(agg_df, directed = TRUE)

    # Plot with ggraph
    p <- ggraph::ggraph(g, layout = layout) +
        ggraph::geom_edge_link(
            ggplot2::aes(width = .data$probability, alpha = .data$probability),
            arrow = ggplot2::arrow(length = ggplot2::unit(4, 'mm')),
            end_cap = ggraph::circle(5, 'mm'),
            color = "grey50"
        ) +
        ggraph::geom_node_point(size = 8, color = "lightblue") +
        ggraph::geom_node_text(ggplot2::aes(label = .data$name), repel = TRUE, size = 5) +
        ggraph::scale_edge_width(range = c(0.5, 3) * edge_width_scale) +
        ggplot2::theme_void() +
        ggplot2::theme(legend.position = "none")

    return(p)
}
#' Plot Cell-Cell Communication Bubble Chart
#'
#' @title sc_cci_bubble
#' @param object A SingleCellExperiment object.
#' @param sources A character vector of source groups to include. If `NULL`, all are included.
#' @param targets A character vector of target groups to include. If `NULL`, all are included.
#' @param id The communication record id to use.
#' @param pval_thresh The p-value threshold for filtering interactions.
#' @param type The type of bubble chart to plot: `"interaction"` (default), `"summary"`, or `"facet"`.
#' @param top_n Integer. If provided, limits the plot to the top N interactions (by probability) per cell pair.
#' @return A ggplot object
#' @importFrom ggplot2 ggplot aes geom_point scale_color_viridis_c scale_size_continuous theme_minimal theme_bw theme element_text element_line element_rect labs facet_grid
#' @importFrom rlang .data
#' @export
sc_cci_bubble <- function(object, sources = NULL, targets = NULL, id = NULL, pval_thresh = 0.05, type = c("interaction", "summary", "facet"), top_n = NULL) {
    type <- match.arg(type)
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
    
    # Handle "summary" type
    if (type == "summary") {
        # Count interactions per source-target pair
        agg_df <- as.data.frame(table(source_group = df$source_group, target_group = df$target_group))
        agg_df <- agg_df[agg_df$Freq > 0, ]
        
        p <- ggplot2::ggplot(agg_df, ggplot2::aes(x = .data$target_group, y = .data$source_group)) +
            ggplot2::geom_point(ggplot2::aes(size = .data$Freq, color = .data$Freq)) +
            ggplot2::scale_color_viridis_c(option = "C", direction = -1) +
            ggplot2::scale_size_continuous(range = c(2, 8)) +
            ggplot2::theme_minimal() +
            ggplot2::theme(
                axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                panel.grid.major = ggplot2::element_line(color = "grey90")
            ) +
            ggplot2::labs(x = "Receiver cell type", y = "Sender cell type", 
                          color = "No. of\ninteractions", size = "No. of\ninteractions")
        return(p)
    }
    
    # For interaction and facet, filter top N if requested
    if (!is.null(top_n)) {
        df <- df[order(df$probability, decreasing = TRUE), ]
        int_names <- paste(df$ligand, df$receptor, sep = "|")
        top_ints <- unique(int_names)[1:min(top_n, length(unique(int_names)))]
        df <- df[int_names %in% top_ints, ]
    }

    df$interaction <- paste(df$ligand, df$receptor, sep = " | ")
    df$cell_pair <- paste(df$source_group, df$target_group, sep = " -> ")
    
    # -log10(p_value) for size or color
    # Add a small epsilon to avoid -log10(0)
    df$mlog10p <- -log10(df$p_value + 1e-10)
    
    if (type == "interaction") {
        p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$cell_pair, y = .data$interaction)) +
            ggplot2::geom_point(ggplot2::aes(size = .data$mlog10p, color = .data$probability)) +
            ggplot2::scale_color_viridis_c(option = "magma", direction = -1) +
            ggplot2::scale_size_continuous(range = c(1, 6)) +
            ggplot2::theme_minimal() +
            ggplot2::theme(
                axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
                panel.grid.major = ggplot2::element_line(color = "grey90")
            ) +
            ggplot2::labs(x = "Sender -> Receiver", y = "Ligand | Receptor", 
                          color = "Probability\n/Score", size = "-log10(P)")
        return(p)
    } else if (type == "facet") {
        p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$target_group, y = .data$interaction)) +
            ggplot2::geom_point(ggplot2::aes(size = .data$mlog10p, color = .data$probability)) +
            ggplot2::scale_color_viridis_c(option = "viridis", direction = -1) +
            ggplot2::scale_size_continuous(range = c(1, 6)) +
            ggplot2::facet_grid(~ source_group, scales = "free_x", space = "free_x") +
            ggplot2::theme_bw() +
            ggplot2::theme(
                axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                strip.background = ggplot2::element_rect(fill = "white"),
                strip.text = ggplot2::element_text(face = "bold")
            ) +
            ggplot2::labs(x = "Target", y = "Interactions (Ligand -> Receptor)", 
                          color = "Score", size = "-log10(P)")
        return(p)
    }
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
    
    # Set default arguments for circlize::chordDiagram if not provided
    args <- list(...)
    if (!"annotationTrack" %in% names(args)) {
        args$annotationTrack <- c("grid", "name") # Remove default 'axis' to hide tiny scientific notation numbers
    }
    if (!"preAllocateTracks" %in% names(args)) {
        args$preAllocateTracks <- list(track.height = max(graphics::strwidth(unlist(dimnames(agg_df)))))
    }
    
    do.call(circlize::chordDiagram, c(list(x = agg_df), args))
}

#' Plot Cell-Cell Communication Sigmoid Diagram
#'
#' @title sc_cci_sigmoid
#' @param object A SingleCellExperiment object.
#' @param id The communication record id to use.
#' @param pval_thresh The p-value threshold for filtering interactions.
#' @param top_n Integer. If provided, limits the plot to the top N interactions (by probability).
#' @return A ggplot object
#' @export
sc_cci_sigmoid <- function(object, id = NULL, pval_thresh = 0.05, top_n = NULL) {
    if (!requireNamespace("ggbump", quietly = TRUE)) {
        stop("Package 'ggbump' is needed for this function to work. Please install it.")
    }
    
    df <- get_cci(object, id = id)
    if (is.null(df)) {
        stop("No CCI data found in the object.")
    }
    
    if (!is.null(pval_thresh)) {
        df <- df[df$p_value <= pval_thresh, , drop = FALSE]
    }
    
    if (!is.null(top_n)) {
        df <- df[order(df$probability, decreasing = TRUE), ]
        df <- utils::head(df, top_n)
    }
    
    if (nrow(df) == 0) {
        stop("No interactions left after filtering.")
    }
    
    # Extract senders (ligands) and receivers (receptors)
    senders <- unique(df[, c("ligand", "source_group")])
    senders <- senders[order(senders$source_group, senders$ligand), ]
    receivers <- unique(df[, c("receptor", "target_group")])
    receivers <- receivers[order(receivers$target_group, receivers$receptor), ]
    
    # Calculate Y coordinates
    max_y <- max(nrow(senders), nrow(receivers))
    if (nrow(senders) > 1) {
        senders$y <- (seq_len(nrow(senders)) - 1) / (nrow(senders) - 1) * max_y
    } else {
        senders$y <- max_y / 2
    }
    senders$x <- 0
    
    if (nrow(receivers) > 1) {
        receivers$y <- (seq_len(nrow(receivers)) - 1) / (nrow(receivers) - 1) * max_y
    } else {
        receivers$y <- max_y / 2
    }
    receivers$x <- 1
    
    # Merge coords back to edges
    df$y_start <- senders$y[match(paste(df$ligand, df$source_group), paste(senders$ligand, senders$source_group))]
    df$y_end <- receivers$y[match(paste(df$receptor, df$target_group), paste(receivers$receptor, receivers$target_group))]
    
    # Nodes data frame
    nodes <- rbind(
        data.frame(name = senders$ligand, group = senders$source_group, x = senders$x, y = senders$y, type = "Sender"),
        data.frame(name = receivers$receptor, group = receivers$target_group, x = receivers$x, y = receivers$y, type = "Receiver")
    )
    
    p <- ggplot2::ggplot() +
        ggbump::geom_sigmoid(data = df, ggplot2::aes(x = 0, xend = 1, y = .data$y_start, yend = .data$y_end, linewidth = .data$probability), 
                             alpha = 0.6, color = "grey30") +
        ggplot2::geom_point(data = nodes, ggplot2::aes(x = .data$x, y = .data$y, color = .data$group), shape = 15, size = 6) +
        ggplot2::geom_text(data = nodes[nodes$x == 0, ], ggplot2::aes(x = .data$x - 0.05, y = .data$y, label = .data$name), hjust = 1) +
        ggplot2::geom_text(data = nodes[nodes$x == 1, ], ggplot2::aes(x = .data$x + 0.05, y = .data$y, label = .data$name), hjust = 0) +
        ggplot2::annotate("text", x = 0, y = max_y + max_y * 0.05, label = "Sender (Ligand)", size = 5, fontface = "bold") +
        ggplot2::annotate("text", x = 1, y = max_y + max_y * 0.05, label = "Receiver (Receptor)", size = 5, fontface = "bold") +
        ggplot2::scale_linewidth_continuous(range = c(0.5, 2)) +
        ggplot2::theme_void() +
        ggplot2::expand_limits(x = c(-0.6, 1.6), y = c(-max_y * 0.05, max_y + max_y * 0.1)) +
        ggplot2::labs(color = "Cell type", linewidth = "Probability") +
        ggplot2::theme(legend.position = "bottom")
        
    return(p)
}

#' Plot Cell-Cell Communication Arrow Diagram
#'
#' @title sc_cci_arrow
#' @param object A SingleCellExperiment object.
#' @param cell_a The first cell type (left side).
#' @param cell_b The second cell type (right side).
#' @param group A character string specifying the column in `colData(object)` containing cell type labels.
#' @param id The communication record id to use.
#' @param pval_thresh The p-value threshold for filtering interactions.
#' @param top_n Integer. Limit to top N interactions.
#' @param layer The assay layer to use for expression calculation (default: "logcounts").
#' @return A ggplot object
#' @export
sc_cci_arrow <- function(object, cell_a, cell_b, group, id = NULL, pval_thresh = 0.05, top_n = NULL, layer = "logcounts") {
    if (!layer %in% SummarizedExperiment::assayNames(object)) {
        stop(paste("Layer", layer, "not found in the object assays."))
    }
    
    if (!group %in% colnames(SummarizedExperiment::colData(object))) {
        stop(paste("Group", group, "not found in colData."))
    }
    
    df <- get_cci(object, id = id)
    if (is.null(df)) {
        stop("No CCI data found in the object.")
    }
    
    # Filter for cell_a <-> cell_b
    df <- df[(df$source_group == cell_a & df$target_group == cell_b) | 
             (df$source_group == cell_b & df$target_group == cell_a), , drop = FALSE]
             
    if (!is.null(pval_thresh)) {
        df <- df[df$p_value <= pval_thresh, , drop = FALSE]
    }
    
    if (!is.null(top_n)) {
        df <- df[order(df$probability, decreasing = TRUE), ]
        df <- utils::head(df, top_n)
    }
    
    if (nrow(df) == 0) {
        stop("No interactions found between these two cell types after filtering.")
    }
    
    edges <- data.frame(
        gene_A = ifelse(df$source_group == cell_a, df$ligand, df$receptor),
        gene_B = ifelse(df$source_group == cell_a, df$receptor, df$ligand),
        direction = ifelse(df$source_group == cell_a, "A_to_B", "B_to_A"),
        probability = df$probability
    )
    
    genes_A <- unique(edges$gene_A)
    genes_B <- unique(edges$gene_B)
    
    # Function to get mean expression safely
    get_mean_exp <- function(genes, cell_type) {
        idx <- which(object[[group]] == cell_type)
        res <- numeric(length(genes))
        names(res) <- genes
        if (length(idx) == 0) return(res)
        
        valid <- genes %in% rownames(object)
        if (any(valid)) {
            mat <- SummarizedExperiment::assay(object, layer)[genes[valid], idx, drop = FALSE]
            res[valid] <- Matrix::rowMeans(mat)
        }
        return(res)
    }
    
    exp_A <- get_mean_exp(genes_A, cell_a)
    exp_B <- get_mean_exp(genes_B, cell_b)
    
    nodes_A <- data.frame(gene = genes_A, exp = exp_A, x = 0, y = seq_along(genes_A))
    nodes_B <- data.frame(gene = genes_B, exp = exp_B, x = 1, y = seq_along(genes_B))
    
    max_y <- max(nrow(nodes_A), nrow(nodes_B))
    if (nrow(nodes_A) > 1) nodes_A$y <- (seq_len(nrow(nodes_A)) - 1) / (nrow(nodes_A) - 1) * max_y else nodes_A$y <- max_y / 2
    if (nrow(nodes_B) > 1) nodes_B$y <- (seq_len(nrow(nodes_B)) - 1) / (nrow(nodes_B) - 1) * max_y else nodes_B$y <- max_y / 2
    
    edges$y_A <- nodes_A$y[match(edges$gene_A, nodes_A$gene)]
    edges$y_B <- nodes_B$y[match(edges$gene_B, nodes_B$gene)]
    
    edges$x_start <- ifelse(edges$direction == "A_to_B", 0.05, 0.95)
    edges$x_end   <- ifelse(edges$direction == "A_to_B", 0.95, 0.05)
    edges$y_start <- ifelse(edges$direction == "A_to_B", edges$y_A, edges$y_B)
    edges$y_end   <- ifelse(edges$direction == "A_to_B", edges$y_B, edges$y_A)
    
    nodes <- rbind(nodes_A, nodes_B)
    
    p <- ggplot2::ggplot() +
        ggplot2::geom_segment(data = edges, ggplot2::aes(x = .data$x_start, xend = .data$x_end, y = .data$y_start, yend = .data$y_end, linewidth = .data$probability),
                              arrow = ggplot2::arrow(length = ggplot2::unit(3, "mm"), type = "closed"), 
                              color = "black") +
        ggplot2::geom_point(data = nodes, ggplot2::aes(x = .data$x, y = .data$y, fill = .data$exp), shape = 22, size = 6, color = "black") +
        ggplot2::geom_text(data = nodes_A, ggplot2::aes(x = .data$x - 0.05, y = .data$y, label = .data$gene), hjust = 1) +
        ggplot2::geom_text(data = nodes_B, ggplot2::aes(x = .data$x + 0.05, y = .data$y, label = .data$gene), hjust = 0) +
        ggplot2::annotate("text", x = -0.3, y = max_y / 2, label = cell_a, angle = 90, fontface = "bold", size = 6) +
        ggplot2::annotate("text", x = 1.3,  y = max_y / 2, label = cell_b, angle = 270, fontface = "bold", size = 6) +
        ggplot2::scale_fill_gradientn(colors = c("#FFF5EB", "#FDB668", "#7F2704")) +
        ggplot2::scale_linewidth_continuous(range = c(0.5, 2)) +
        ggplot2::theme_void() +
        ggplot2::expand_limits(x = c(-0.4, 1.4), y = c(-max_y * 0.05, max_y + max_y * 0.05)) +
        ggplot2::labs(fill = "Mean expression", linewidth = "Probability") +
        ggplot2::theme(legend.position = "bottom")
        
    return(p)
}
