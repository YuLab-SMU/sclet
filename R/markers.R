#' Gene expression markers for all identity classes
#' 
#' Find markers for each of the identity classes
#' 
#' @param object A 'SingleCellExperiment' object
#' @param min.pct  only test genes that are detected in a minimum fraction of
#' min.pct cells in either of the two populations.
#' @param logfc.threshold Limit testing to genes which show, on average, at least
#' X-fold difference (log-scale) between the two groups of cells. Default is 0.1
#' @param return.thresh Only return markers that have a p-value < return.thresh, or a power > return.thresh (if the test is ROC)
#' @param only.pos Only return positive markers (FALSE by default)
#' @param ... additional parameters
#' @export
FindAllMarkers <- function(
    object, 
    min.pct = 0.01,
    logfc.threshold = 0.1,
    return.thresh=0.01,
    only.pos = FALSE,
    ...) {

    idents.all <- levels(Idents(object))
                        
    markers_list <- list()

    for (i in seq_along(idents.all)) {
        message("Calculating cluster ", idents.all[i])
        
        markers <- FindMarkers(object = object, 
        ident.1 = idents.all[i],
        ident.2=NULL,
        min.pct = min.pct, 
        logfc.threshold = logfc.threshold, 
        ...)
        
        markers$cluster <- idents.all[i]
        markers$gene <- rownames(markers)
        
        markers <- markers[markers$pval < return.thresh, ]
        
        if (only.pos) {
            j <- grep("avg_log", names(markers))
            markers <- markers[markers[ ,j] > 0, ] 
        }  
        
        markers_list[[i]] <- markers
    }

    combined_markers <- yulab.utils::rbindlist(markers_list)
    
    rownames(combined_markers) <- make.unique(names=combined_markers$gene, sep = "")
    
    return(combined_markers)
}

#' Find markers 
#' 
#' @title FindMarkers
#' @param object A 'SingleCellExperiment' object
#' @param min.pct  only test genes that are detected in a minimum fraction of
#' min.pct cells in either of the two populations.
#' @param ident.1 Identity class to define markers for
#' @param ident.2 A second identity class for comparison. If NULL (default) -
#' use all other cells for comparison.
#' @param logfc.threshold Limit testing to genes which show, on average, at least
#' X-fold difference (log-scale) between the two groups of cells. Default is 0.1
#' @param base log base. Default is 2
#' @param pseudocount.use Pseudocount to add to averaged expression values when
#' calculating logFC. 1 by default.
#' @param ... additional parameters
#' @export
FindMarkers <- function(
    object, 
    ident.1 = NULL, 
    ident.2 = NULL, 
    min.pct = 0.01, 
    logfc.threshold = 0.1, 
    base=2,
    pseudocount.use=1,
    ...) {

        FindMarkers_Presto(
        object = object@assays@data$logcounts,
        ident.1 = ident.1,
        ident.2 = ident.2,
        clusters = Idents(object),
        min.pct = min.pct,
        logfc.threshold = logfc.threshold,
        base = base,
        pseudocount.use = pseudocount.use,
        ...
        )
}

#' @importFrom yulab.utils install_zip_gh
#' @importFrom yulab.utils is.installed
#' @importFrom Matrix rowSums
FindMarkers_Presto <- function(object, ident.1 = NULL, ident.2 = NULL, clusters, 
                              min.pct = 0.01,logfc.threshold = 0.1, base=2,
                              pseudocount.use=1, ...) {

    if (is.null(ident.1)) stop("ident.1 must be specified")
    if (is.null(clusters)) stop("clusters must be specified")

    if (is.null(ident.2)) {
        ident.2 <- setdiff(unique(clusters), ident.1)
    }

    cells.1 <- which(clusters %in% ident.1)
    cells.2 <- which(clusters %in% ident.2)

    if (length(cells.1) == 0 || length(cells.2) == 0) {
        stop("One or both of the identity groups have no cells")}
    
    mean.fxn <- function(x) {
        log((Matrix::rowSums(expm1(x = x)) + pseudocount.use)/NCOL(x), base = base)
    }

    fc.name <- sprintf("avg_log%dFC", base)
    features <- rownames(x = object)
    thresh.min <- 0
    mat <- object[features, ]
    pct.1 <- round(Matrix::rowSums(mat[, cells.1, drop = FALSE] > thresh.min)/length(cells.1), digits = 3)
    pct.2 <- round(Matrix::rowSums(mat[, cells.2, drop = FALSE] > thresh.min)/length(cells.2), digits = 3)
    data.1 <- mean.fxn(mat[, cells.1, drop = FALSE])
    data.2 <- mean.fxn(mat[, cells.2, drop = FALSE])
    fc <- (data.1 - data.2) 
    fc.results <- as.data.frame(cbind(fc, pct.1, pct.2)) 
    colnames(fc.results) <- c(fc.name, "pct.1", "pct.2")
    fc.results$feature <- features

    data_subset <- object[ ,c(cells.1, cells.2), drop = FALSE]

    group_labels <- c(rep("group1", length(cells.1)), rep("group2", length(cells.2)))

    if (!is.installed("presto")) {
        install_zip_gh("immunogenomics/presto")
    }

    wilcoxauc <- get_fun_from_pkg("presto", "wilcoxauc")
    results <- wilcoxauc(data_subset, group_labels)
    

    results <- results[1:(nrow(x = results)/2), ]
    p.results <- results[, c("feature", "pval","padj")]

    combined_results <- merge(p.results, fc.results, by = 'feature')
    combined_results <- combined_results[combined_results$pct.1>= min.pct | combined_results$pct.2 >= min.pct, ]
    combined_results <- combined_results[abs(combined_results[,fc.name]) >= logfc.threshold, ]
    combined_results <- combined_results[order(combined_results$pval, 
                                -abs(combined_results$pct.1 - combined_results$pct.2)), ]
    rownames(combined_results) <- combined_results$feature
    combined_results$feature <- NULL
    return(combined_results)
}
