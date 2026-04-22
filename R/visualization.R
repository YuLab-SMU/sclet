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
#' @param method one of 'seurat' or 'scran'. If NULL, use stored method or infer from rowData.
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

    if (is.null(method)) method <- sclet_get_hvg_method(object)
    if (is.null(method)) {
        method <- if ("variance.standardized" %in% names(d)) "seurat" else "scran"
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
