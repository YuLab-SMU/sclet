#' read 10x single cell expression matrix folder
#' 
#' @title Read10X
#' @param data.dir directory of the data
#' @param ... additional parameters passed to 'DropletUtils::read10xCounts()'
#' @return a 'SingleCellExperiment' object
#' @importFrom DropletUtils read10xCounts
#' @importFrom scuttle uniquifyFeatureNames
#' @importFrom SummarizedExperiment rowData
#' @importFrom SummarizedExperiment colData
#' @export
Read10X <- function(data.dir, ...) {
    sce <- DropletUtils::read10xCounts(samples = data.dir, ...)
    rownames(sce) <- scuttle::uniquifyFeatureNames(
        rowData(sce)$ID, rowData(sce)$Symbol)

    sce <- QCMetrics(sce)

    return(sce)
}

#' add QC metrics
#' 
#' @title QCMetrics
#' @param object a 'SingleCellExperiment' object
#' @return update object with two QC metrics ('nFeature_RNA' and 'nCount_RNA')
#' @importFrom scuttle perCellQCMetrics
#' @export
QCMetrics <- function(object) {
    qc_metrics <- scuttle::perCellQCMetrics(object)
    SummarizedExperiment::colData(object)$nFeature_RNA <- qc_metrics$detected
    SummarizedExperiment::colData(object)$nCount_RNA <- qc_metrics$sum
    return(object)
}

#' calculate percentage of genes that matched a pattern
#' 
#' @title PercentageFeatureSet
#' @param object a SingleCellExperiment object
#' @param pattern the pattern to search for
#' @param feature search pattern from feature of the object. If NULL, use the `rownames(object)`, otherwise, use `rowData(object)[[feature]]`
#' @return the percentage of each cell
#' @export
PercentageFeatureSet <- function(object, pattern = NULL, feature=NULL) {
    if (is.null(feature)) {
        features <- rownames(object)
    } else {
        features <- rowData(object)[[feature]]
    }

    has_pattern <- grep(pattern, features)
    qc_metrics <- scuttle::perCellQCMetrics(object, subsets = list(pattern=has_pattern))
    return(qc_metrics$subsets_pattern_percent)
}


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



#' normalize data
#' 
#' @title NormalizeData
#' @param object a SingleCellExperiment object
#' @param scale.factor scale factor
#' @return a SingleCellExperiment object with 'logcounts' assay created/updated.
#' @importFrom SummarizedExperiment 'assay<-'
#' @importFrom yulab.utils get_fun_from_pkg
#' @export
NormalizeData <- function(object, scale.factor = 10000) {
    counts_matrix <- counts(object)
    LogNorm <- get_fun_from_pkg('Seurat', 'LogNorm')
    logcounts <- LogNorm(counts_matrix, scale_factor = scale.factor)
    colnames(logcounts) <- colnames(counts_matrix)
    rownames(logcounts) <- rownames(counts_matrix)

    SummarizedExperiment::assay(object, 'logcounts') <- logcounts

    return(object)
}

#' identify variable features
#' 
#' @title FindVariableFeatures
#' @param object a SingleCellExperiment object
#' @param nfeatures number of features to be selected as highly variable features
#' @return an updated SingleCellExperiment object with identified highly variable features
#' @importFrom SingleCellExperiment counts
#' @importFrom SingleCellExperiment logcounts
#' @importFrom SummarizedExperiment 'rowData<-'
#' @importFrom stats loess
#' @export
FindVariableFeatures <- function(object, nfeatures = 2000) {
    sce <- object
    if (!is.null(sce@metadata$nVariableFeatures)) {
        sce@metadata$nVariableFeatures <- nfeatures
        return(sce)        
    }

    object <- counts(sce)
    SparseRowVar2 <- get_fun_from_pkg('Seurat', "SparseRowVar2")
    SparseRowVarStd <- get_fun_from_pkg('Seurat', "SparseRowVarStd")

    ## taken from Seurat
    clip.max <- sqrt(x = ncol(x = object))
    hvf.info <- data.frame(mean = Matrix::rowMeans(object))
    hvf.info$variance <- SparseRowVar2(
      mat = object,
      mu = hvf.info$mean,
      display_progress = TRUE
    )
    hvf.info$variance.expected <- 0
    hvf.info$variance.standardized <- 0
    not.const <- hvf.info$variance > 0
    fit <- loess(
      formula = log10(x = variance) ~ log10(x = mean),
      data = hvf.info[not.const, ],
      span = .3
    )
    hvf.info$variance.expected[not.const] <- 10 ^ fit$fitted
    hvf.info$variance.standardized <- SparseRowVarStd(
      mat = object,
      mu = hvf.info$mean,
      sd = sqrt(hvf.info$variance.expected),
      vmax = clip.max,
      display_progress = TRUE
    )
    ## end of seurat code

    rd <- rowData(sce)
    SummarizedExperiment::rowData(sce) <- cbind(rd, hvf.info[rownames(rd),])
    sce@metadata$nVariableFeatures <-nfeatures
    return(sce)
}

#' get variable features
#' 
#' @title VariableFeatures
#' @param object a SingleCellExperiment object
#' @return highly variable features
#' @export
VariableFeatures <- function(object) {
    nfeatures <- object@metadata$nVariableFeatures
    if (is.null(nfeatures)) {
        stop("You should run 'FindVariableFeatures' first.")
    }

    d <- rowData(object)
    d <- d[rownames(d) %in% rownames(object), ]
    i <- order(d$variance.standardized, decreasing = TRUE)
    d <- d[i, ]
    res <- rownames(d)[1:nfeatures]
    return(res)
}

#' variable feature plot
#' 
#' @title VariableFeaturePlot
#' @param object a SingleCellExperiment object
#' @param label display selected features
#' @return scatter plot
#' @export
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggplot2 scale_x_log10
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 scale_color_manual
VariableFeaturePlot <- function(object, label = NULL) {
    d <- as.data.frame(rowData(object))
    d$type <- "Non-variable"
    d$type[d$Symbol %in% VariableFeatures(object)] <- "Variable"
    legend <- sprintf("%s count: %d", unique(d$type), table(d$type))

    p <- ggplot(d, aes(.data$mean, .data$variance.standardized)) + 
        geom_point(aes(color = .data$type)) + 
        scale_x_log10() +
        xlab("Average Expression") + 
        ylab("Standardized Variance") +
        scale_color_manual(values = c("black", "red"), name="", labels=legend) +
        theme_classic()

    if (is.null(label)) return(p)

    d2 <- d[d$Symbol %in% label, ]
    p + ggrepel::geom_text_repel(
            aes(label = .data$Symbol), 
            data = d2
        )
}

#' scale data
#' 
#' @title ScaleData
#' @param object a SingleCellExperiment object
#' @param features selected features to be scaled, default is all features
#' @param assay selected assay to be scaled
#' @return an updated SingleCellExperiment object with 'scaled' assay
#' @export
#' @importFrom stats sd
ScaleData <- function(object, features = NULL, assay = "logcounts") {
    if (is.null(features)) {
        # all genes
        features <- rownames(object)
    } 

    # mat <- t(logcounts(object))
    mat <- Matrix::t(assay(object, assay)[features, ])
    m <- Matrix::rowMeans(mat)
    s <- apply(mat, 1, stats::sd)
    s[s == 0] <- 0.01
    SummarizedExperiment::assay(object, "scaled") <- Matrix::t((mat - m) / s) 
    return(object)
}

#' elbow plot 
#' 
#' @title ElbowPlot
#' @param object a SingleCellExperiment object
#' @return elbow plot
#' @export
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 theme_classic
#' @importFrom rlang .data
ElbowPlot <- function(object) {
    pca_results <- reducedDim(object, "PCA")

    v <- attr(pca_results, "percentVar")
    d <- data.frame(PC = seq_along(v), var = v)

    ggplot(d, aes(.data$PC, .data$var)) + 
        geom_point() +
        ylab("Standard Deviation") +
        theme_classic()      
}

#' @importFrom SingleCellExperiment 'reducedDim<-'
set_dimred <- function(object, dims) {
    reducedDim(object, ".dimred") <- reducedDim(object, "PCA")[, dims]
    return(object)
}

#' build KNN graph
#' 
#' @title FindNeighbors
#' @param object a SingleCellExperiment object
#' @param dims number of dimensions to be used to build the KNN graph
#' @param k number of neighbors
#' @return updated SingleCellExperiment object with KNN graph
#' @importFrom scran buildSNNGraph
#' @export
FindNeighbors <- function(object, dims, k = 10) {
    object <- set_dimred(object, dims)    
    object@metadata$knn_graph <- scran::buildSNNGraph(
        object, use.dimred = ".dimred", k = k, type = "rank")
    return(object)
}

#' identify clusters
#' 
#' @title FindClusters
#' @param object a SingleCellExperiment object
#' @param resolution resolution in clustering
#' @return updated SingleCellExperiment object with clustering memberships
#' @importFrom igraph cluster_louvain
# @importFrom SingleCellExperiment 'colLabels<-'
#' @export
FindClusters <- function(object, resolution = 1) {
    g <- object@metadata$knn_graph
    clusters <- igraph::cluster_louvain(g, resolution=resolution)
    
    SingleCellExperiment::colLabels(object) <- factor(clusters$membership)
    return(object)
}

#' access the clustering memberships
#' 
#' @title Idents
#' @param object a SingleCellExperiment object
#' @return clustering memberships
#' @export
#' @importFrom SummarizedExperiment colData
#' @importFrom stats setNames
Idents <- function(object) {
    d <- colData(object)
    setNames(d$label, d$Barcode)
}

#' run umap
#' 
#' @title RunUMAP
#' @param object a SingleCellExperiment object
#' @param dims dimensions used in UMAP
#' @return updated SingleCellExperiment object with UMAP dimension reduction 
#' @export
#' @importFrom scater runUMAP
RunUMAP <- function(object, dims) {
    object <- set_dimred(object, dims)
    scater::runUMAP(object, dimred = '.dimred')    
}

#' rename cluster ids
#' 
#' @title RenameIdents
#' @param object a SingleCellExperiment object
#' @param new_ids new ID labels
#' @return updated SingleCellExperiment object
#' @importFrom SingleCellExperiment colLabels
#' @importFrom SingleCellExperiment "colLabels<-"
#' @export 
RenameIdents <- function(object, new_ids) {
    old_ids <- colLabels(object)
    lv <- levels(old_ids)

    if (!is.null(names(new_ids))) {
        new_ids <- new_ids[lv]
    }

    colLabels(object) <- factor(
        old_ids, 
        levels = lv, 
        labels = new_ids
    )
    
    return(object)
}

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
        object = assay(object),
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

    cells.1 <- which(clusters == ident.1)
    cells.2 <- which(clusters == ident.2)

    if (length(cells.1) == 0 || length(cells.2) == 0) {
        stop("One or both of the identity groups have no cells")}
    
    mean.fxn <- function(x) {
        log((Matrix::rowSums(expm1(x = x)) + pseudocount.use)/NCOL(x), base = base)
    }

    fc.name <- sprintf("avg_log%dFC", base)
    features <- rownames(x = object)
    thresh.min <- 0
    mat <- object[features, ]
    pct.1 <- round(Matrix::rowSums(mat[, cells.1, drop = FALSE]) > thresh.min/length(cells.1), digits = 3)
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

    combined_results <- cbind(p.results,fc.results)
    combined_results <- merge(p.results, fc.results, by = 'feature')
    combined_results <- combined_results[combined_results$pct.1>= min.pct | combined_results$pct.2 >= min.pct, ]
    combined_results <- combined_results[abs(combined_results[,fc.name]) >= logfc.threshold, ]
    combined_results <- combined_results[order(combined_results$pval, 
                                -abs(combined_results$pct.1 - combined_results$pct.2)), ]
    rownames(combined_results) <- combined_results$feature
    combined_results$feature <- NULL
    return(combined_results)
}


