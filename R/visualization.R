#' violine plot of selected features
#' 
#' @title VlnPlot
#' @param object a SingleCellExperiment object
#' @param features cell features (e.g., nCounts)
#' @param ncol number of columns if multiple features selected
#' @return violin plot
#' @export
#' @importFrom aplot plot_list
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 geom_violin
#' @importFrom ggplot2 geom_jitter
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 theme_classic
#' @importFrom ggplot2 element_blank
VlnPlot <- function(object, features, ncol="auto") {
    pp <- lapply(features, function(y) {
        ## scater::plotColData(object, y=y) +
        ggplot(colData(object), aes(x=factor(1), y=.data[[y]])) +
        ggtitle(y) +
        geom_violin(aes(fill=y), color="black") +
        geom_jitter(color='black', size=.1, width=.4) +
        theme_classic() + 
        theme(legend.position="none",
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()) +
        xlab(NULL) +
        ylab(NULL)
    })

    if (ncol == "auto") ncol <- length(features)
    aplot::plot_list(gglist=pp, ncol=ncol)
}

#' scatter plot to compare 2 features
#' 
#' @title FeatureScatter 
#' @param object a SingleCellExperiment object
#' @param feature1 selected feature 
#' @param feature2 selected feature
#' @return a scatter plot
#' @importFrom ggplot2 geom_point
#' @export
FeatureScatter <- function(object, feature1, feature2) {
    # scater::plotColData(object, x=feature1, y=feature2) +
        #scale_x_log10() + 
     ggplot(colData(object), aes(x=.data[[feature1]], y=.data[[feature2]])) +
        geom_point(color = "red",size=0.5) +
    theme_classic()
}

#' variable feature plot
#' 
#' @title VariableFeaturePlot
#' @param object a SingleCellExperiment object
#' @param label display selected features
#' @param method one of `'scran'` or legacy `'seurat'`. If `NULL`, sclet
#'   prefers Bioconductor-native HVG metadata and does not silently fall back
#'   to legacy Seurat-style statistics.
#' @return scatter plot
#' @export
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggplot2 scale_x_log10
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 scale_color_manual
VariableFeaturePlot <- function(object, label = NULL, method = NULL) {
    d <- as.data.frame(rowData(object))
    gene <- if ("Symbol" %in% names(d)) d$Symbol else rownames(object)
    d$gene <- gene

    if (is.null(method)) method <- sclet_resolve_hvg_method(object)
    if (is.null(method)) {
        stop(
            "No Bioconductor-native HVG statistics found. Run `FindVariableFeatures()` ",
            "with `method = \"scran\"` or `method = \"scrapper\"`, or set ",
            "`method = \"seurat\"` explicitly for legacy objects."
        )
    }
    method <- match.arg(method, c("seurat", "scran"))

    if (method == "seurat") {
        ycol <- "variance.standardized"
        ylab_txt <- "Standardized Variance"
    } else {
        if ("bio" %in% names(d)) {
            ycol <- "bio"
            ylab_txt <- "Biological Variance"
        } else if ("total" %in% names(d)) {
            ycol <- "total"
            ylab_txt <- "Total Variance"
        } else {
            stop("scran HVG stats not found in rowData(): expect column 'bio' or 'total'.")
        }
    }

    hvgs <- VariableFeatures(object, method = method)
    d$type <- "Non-variable"
    d$type[d$gene %in% hvgs] <- "Variable"
    tab <- table(d$type)
    legend <- sprintf("%s count: %d", names(tab), as.integer(tab))

    p <- ggplot(d, aes(.data$mean, .data[[ycol]])) + 
        geom_point(aes(color = .data$type)) + 
        scale_x_log10() +
        xlab("Average Expression") + 
        ylab(ylab_txt) +
        scale_color_manual(values = c("black", "red"), name="", labels=legend) +
        theme_classic()

    if (is.null(label)) return(p)

    d2 <- d[d$gene %in% label, ]
    p + ggrepel::geom_text_repel(
            aes(label = .data$gene), 
            data = d2
        )
}

#' Dimensional reduction plot
#' 
#' @title DimPlot
#' @param object a SingleCellExperiment object
#' @param reduction reduction method. If NULL, use the active reduction tracked by sclet.
#' @param ... additional parameters passed to ggsc::sc_dim
#' @return ggplot object
#' @importFrom ggsc sc_dim
#' @export
DimPlot <- function(object, reduction = NULL, ...) {
    if (is.null(reduction)) {
        reduction <- DefaultReduction(object)
    }
    if (is.null(reduction)) {
        stop("No dimensional reduction found in object.")
    }
    ggsc::sc_dim(object, reduction = reduction, ...)
}

#' Cell Dimensional Reduction Plot
#' 
#' Dimensional reduction plot that defaults to the active reduction and active identity.
#' 
#' @title CellDimPlot
#' @param object a SingleCellExperiment object
#' @param group.by cell identity group. If NULL, uses `ActiveIdent(object)`.
#' @param reduction reduction method. If NULL, uses `DefaultReduction(object)`.
#' @param ... additional parameters passed to `ggsc::sc_dim`
#' @return ggplot object
#' @export
CellDimPlot <- function(object, group.by = NULL, reduction = NULL, ...) {
    if (is.null(reduction)) {
        reduction <- DefaultReduction(object)
    }
    if (is.null(reduction)) {
        stop("No dimensional reduction found in object.")
    }
    if (is.null(group.by)) {
        group.by <- ActiveIdent(object)
    }
    
    if (!is.null(group.by)) {
        valid_cols <- colnames(SummarizedExperiment::colData(object))
        if (group.by == "colLabels" || group.by %in% valid_cols) {
            if (group.by == "colLabels") {
                 group_vec <- SingleCellExperiment::colLabels(object)
            } else {
                 group_vec <- SummarizedExperiment::colData(object)[[group.by]]
            }
            
            p <- ggsc::sc_dim(object, reduction = reduction, ...)
            p$data[["_CellDimPlot_Group"]] <- group_vec
            p$mapping <- utils::modifyList(p$mapping, aes(color = .data[["_CellDimPlot_Group"]]))
            
            if (length(p$layers) > 0) {
                p$layers[[1]]$mapping <- utils::modifyList(
                    if(is.null(p$layers[[1]]$mapping)) aes() else p$layers[[1]]$mapping,
                    aes(color = .data[["_CellDimPlot_Group"]])
                )
            }
            return(p)
        } else {
             warning(paste("group.by column", group.by, "not found. Plotting without grouping."))
        }
    }
    ggsc::sc_dim(object, reduction = reduction, ...)
}

#' Feature Dimensional Reduction Plot
#' 
#' Plot feature expression on a dimensional reduction, defaulting to the active reduction and active layer.
#' 
#' @title FeatureDimPlot
#' @param object a SingleCellExperiment object
#' @param features features to plot
#' @param layer layer to extract expression from. If NULL, uses `DefaultLayer(object)`.
#' @param reduction reduction method. If NULL, uses `DefaultReduction(object)`.
#' @param ... additional parameters passed to `ggsc::sc_feature`
#' @return ggplot object
#' @importFrom ggsc sc_feature
#' @export
FeatureDimPlot <- function(object, features, layer = NULL, reduction = NULL, ...) {
    if (is.null(reduction)) {
        reduction <- DefaultReduction(object)
    }
    if (is.null(reduction)) {
        stop("No dimensional reduction found in object.")
    }
    
    source <- sclet_resolve_expression_source(
        object = object,
        layer = layer,
        assay = NULL,
        context = "FeatureDimPlot"
    )
    
    ggsc::sc_feature(object, features = features, slot = source$assay, reduction = reduction, ...)
}

#' Group Heatmap Plot
#' 
#' Plot a heatmap of features across cell groups.
#' If `ggsc` handles it, we use it, otherwise fall back to `scater`.
#' 
#' @title GroupHeatmap
#' @param object a SingleCellExperiment object
#' @param features features to plot
#' @param group.by cell identity group. If NULL, uses `ActiveIdent(object)`.
#' @param layer layer to extract expression from. If NULL, uses `DefaultLayer(object)`.
#' @param ... additional parameters passed to underlying plotting function.
#' @return ggplot or ComplexHeatmap object
#' @export
GroupHeatmap <- function(object, features, group.by = NULL, layer = NULL, ...) {
    if (is.null(group.by)) {
        group.by <- ActiveIdent(object)
    }
    
    if (is.null(group.by)) {
        stop("No active identity found. Please provide 'group.by'.")
    }
    
    source <- sclet_resolve_expression_source(
        object = object,
        layer = layer,
        assay = NULL,
        context = "GroupHeatmap"
    )
    
    if (!requireNamespace("scater", quietly = TRUE)) {
        stop("Package 'scater' is required for GroupHeatmap. Please install it.")
    }
    
    group_col <- group.by
    if (group.by == "colLabels") {
        group_col <- "label"
    }
    
    scater::plotGroupedHeatmap(object, features = features, group = group_col, exprs_values = source$assay, ...)
}

#' Cell Statistics Plot
#' 
#' Plot the proportion or count of cells across groups (e.g. clusters across batches).
#' 
#' @title CellStatPlot
#' @param object a SingleCellExperiment object
#' @param group.by the primary grouping variable (e.g. cluster)
#' @param split.by the secondary splitting variable (e.g. batch or condition). If NULL, plots overall counts of `group.by`.
#' @param position position adjustment for ggplot2 geom_bar (e.g. "fill" for proportion, "stack" or "dodge" for counts)
#' @param ... additional parameters passed to ggplot2
#' @return ggplot object
#' @importFrom ggplot2 ggplot aes geom_bar theme_classic ylab xlab theme
#' @export
CellStatPlot <- function(object, group.by = NULL, split.by = NULL, position = "fill", ...) {
    if (is.null(group.by)) {
        group.by <- ActiveIdent(object)
    }
    if (is.null(group.by)) {
        stop("No active identity found. Please provide 'group.by'.")
    }
    
    d <- as.data.frame(SummarizedExperiment::colData(object))
    
    if (group.by == "colLabels") {
        d[[group.by]] <- SingleCellExperiment::colLabels(object)
    }
    
    if (!group.by %in% colnames(d)) {
        stop(paste("Column", group.by, "not found in colData."))
    }
    
    if (!is.null(split.by)) {
        if (split.by == "colLabels") {
            d[[split.by]] <- SingleCellExperiment::colLabels(object)
        }
        if (!split.by %in% colnames(d)) {
            stop(paste("Column", split.by, "not found in colData."))
        }
        
        p <- ggplot(d, aes(x = .data[[split.by]], fill = .data[[group.by]])) +
            geom_bar(position = position, ...) +
            theme_classic() +
            xlab(split.by)
            
        if (position == "fill") {
            p <- p + ylab("Proportion")
        } else {
            p <- p + ylab("Count")
        }
    } else {
        p <- ggplot(d, aes(x = .data[[group.by]], fill = .data[[group.by]])) +
            geom_bar(stat = "count", ...) +
            theme_classic() +
            theme(legend.position = "none") +
            xlab(group.by) + ylab("Count")
    }
    return(p)
}

#' Differential Expression Test Volcano Plot
#' 
#' Plot a volcano plot from detest state results.
#' 
#' @title DEtestPlot
#' @param object A SingleCellExperiment object
#' @param id Analysis record id. If NULL, uses the active detest state.
#' @param p_val_cutoff p-value threshold for significance (default 0.05)
#' @param logfc_cutoff logFC threshold for significance (default 0.25)
#' @param ... additional parameters passed to ggplot2
#' @return ggplot object
#' @importFrom ggplot2 ggplot aes geom_point theme_classic scale_color_manual geom_vline geom_hline xlab ylab
#' @export
DEtestPlot <- function(object, id = NULL, p_val_cutoff = 0.05, logfc_cutoff = 0.25, ...) {
    detest_state <- get_detest(object, id)
    if (is.null(detest_state)) {
        stop("No detest state found. Please run RunDEtest() first.")
    }
    
    res <- detest_state$artifacts$result
    if (!"pval" %in% colnames(res)) {
        if ("p_val" %in% colnames(res)) colnames(res)[colnames(res) == "p_val"] <- "pval"
    }
    if (!"padj" %in% colnames(res)) {
        if ("p_val_adj" %in% colnames(res)) colnames(res)[colnames(res) == "p_val_adj"] <- "padj"
    }
    
    fc_col <- grep("log.*FC", colnames(res), value = TRUE, ignore.case = TRUE)
    if (length(fc_col) == 0) {
        stop("Could not find logFC column in detest results.")
    }
    fc_col <- fc_col[1]
    
    res$log10pval <- -log10(res$pval)
    
    # Define significance
    res$Significance <- "Not Sig"
    res$Significance[res$pval < p_val_cutoff & res[[fc_col]] > logfc_cutoff] <- "Up"
    res$Significance[res$pval < p_val_cutoff & res[[fc_col]] < -logfc_cutoff] <- "Down"
    
    p <- ggplot(res, aes(x = .data[[fc_col]], y = .data$log10pval, color = .data$Significance)) +
        geom_point(alpha = 0.6, size = 1.2) +
        scale_color_manual(values = c("Down" = "#2b8cbe", "Not Sig" = "grey80", "Up" = "#e34a33")) +
        geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), linetype = "dashed", color = "grey50") +
        geom_hline(yintercept = -log10(p_val_cutoff), linetype = "dashed", color = "grey50") +
        theme_classic() +
        xlab("logFC") +
        ylab("-log10(p-value)")
        
    return(p)
}

#' Enrichment Plot
#' 
#' Plot enrichment results using dotplot from enrichplot.
#' 
#' @title EnrichmentPlot
#' @param object A SingleCellExperiment object
#' @param id Analysis record id. If NULL, uses the active enrichment state.
#' @param showCategory Number of categories to show
#' @param ... additional parameters passed to enrichplot::dotplot
#' @return ggplot object
#' @export
EnrichmentPlot <- function(object, id = NULL, showCategory = 10, ...) {
    enrich_state <- get_enrichment(object, id)
    if (is.null(enrich_state)) {
        stop("No enrichment state found. Please run RunEnrichment() first.")
    }
    
    res <- enrich_state$artifacts$result
    
    if (!requireNamespace("enrichplot", quietly = TRUE)) {
        stop("Package 'enrichplot' is needed for this function to work.")
    }
    
    p <- enrichplot::dotplot(res, showCategory = showCategory, ...)
    return(p)
}

#' Projection Plot for Reference Mapping
#' 
#' Plot reference and query cells in the same dimensional reduction space.
#' 
#' @title ProjectionPlot
#' @param query A SingleCellExperiment object representing the query dataset
#' @param ref A SingleCellExperiment object representing the reference dataset
#' @param reduction The name of the dimensional reduction to use (must be present in both objects)
#' @param group.by The metadata column name containing the labels. If NULL, tries to use the active identity or mapping state from query.
#' @param ... additional parameters passed to ggplot2
#' @return ggplot object
#' @importFrom ggplot2 ggplot aes geom_point theme_classic scale_color_manual guides guide_legend xlab ylab
#' @export
ProjectionPlot <- function(query, ref, reduction = "UMAP", group.by = NULL, ...) {
    if (!reduction %in% SingleCellExperiment::reducedDimNames(query)) {
        stop(paste("Reduction", reduction, "not found in query object."))
    }
    if (!reduction %in% SingleCellExperiment::reducedDimNames(ref)) {
        stop(paste("Reduction", reduction, "not found in reference object."))
    }
    
    # Try to resolve group.by
    query_group_col <- group.by
    if (is.null(query_group_col)) {
        # Check if there is an active mapping state
        map_state <- get_mapping(query)
        if (!is.null(map_state) && !is.null(map_state$artifacts$labels_col)) {
            query_group_col <- map_state$artifacts$labels_col
        } else {
            query_group_col <- ActiveIdent(query)
        }
    }
    
    if (is.null(query_group_col) || !query_group_col %in% colnames(SummarizedExperiment::colData(query))) {
        stop("Could not resolve grouping column for query. Please specify 'group.by'.")
    }
    
    ref_group_col <- group.by
    if (is.null(ref_group_col) || !ref_group_col %in% colnames(SummarizedExperiment::colData(ref))) {
        # Try some common names
        for (cname in c("label.main", "label", "CellType", "celltype")) {
            if (cname %in% colnames(SummarizedExperiment::colData(ref))) {
                ref_group_col <- cname
                break
            }
        }
        if (is.null(ref_group_col)) {
             stop("Could not resolve grouping column for reference. Please specify 'group.by'.")
        }
    }
    
    # Extract coordinates
    query_coords <- as.data.frame(SingleCellExperiment::reducedDim(query, reduction)[, 1:2])
    colnames(query_coords) <- c("Dim1", "Dim2")
    query_coords$Label <- SummarizedExperiment::colData(query)[[query_group_col]]
    query_coords$Dataset <- "Query"
    
    ref_coords <- as.data.frame(SingleCellExperiment::reducedDim(ref, reduction)[, 1:2])
    colnames(ref_coords) <- c("Dim1", "Dim2")
    ref_coords$Label <- SummarizedExperiment::colData(ref)[[ref_group_col]]
    ref_coords$Dataset <- "Reference"
    
    # Combine
    plot_df <- rbind(ref_coords, query_coords)
    
    # We plot Reference first (larger, lower alpha), then Query on top (smaller, opaque)
    # Using ggplot2 directly as we are merging two datasets
    p <- ggplot() +
        geom_point(data = plot_df[plot_df$Dataset == "Reference", ], 
                   aes(x = .data$Dim1, y = .data$Dim2, color = .data$Label), 
                   alpha = 0.3, size = 1.5) +
        geom_point(data = plot_df[plot_df$Dataset == "Query", ], 
                   aes(x = .data$Dim1, y = .data$Dim2, color = .data$Label), 
                   alpha = 1.0, size = 0.8, shape = 17) + # shape 17 is triangle
        theme_classic() +
        xlab(paste0(reduction, "_1")) +
        ylab(paste0(reduction, "_2")) +
        guides(color = guide_legend(override.aes = list(alpha = 1, size = 2, shape = 16)))
        
    return(p)
}

#' @name DimPlot
#' @aliases plot_dim_scatter
#' @export
plot_dim_scatter <- DimPlot

#' @name CellDimPlot
#' @aliases plot_cell_dim
#' @export
plot_cell_dim <- CellDimPlot

#' @name FeatureDimPlot
#' @aliases plot_feature_dim
#' @export
plot_feature_dim <- FeatureDimPlot

#' @name GroupHeatmap
#' @aliases plot_group_heatmap
#' @export
plot_group_heatmap <- GroupHeatmap

#' @name CellStatPlot
#' @aliases plot_cell_stat
#' @export
plot_cell_stat <- CellStatPlot

#' @name DEtestPlot
#' @aliases plot_detest
#' @export
plot_detest <- DEtestPlot

#' @name EnrichmentPlot
#' @aliases plot_enrichment
#' @export
plot_enrichment <- EnrichmentPlot

#' @name ProjectionPlot
#' @aliases plot_projection
#' @export
plot_projection <- ProjectionPlot
