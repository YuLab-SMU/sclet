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

#' run pca
#' 
#' @title runPCA
#' @param object a SingleCellExperiment object
#' @param subset_row subset of rows used for PCA
#' @param exprs_values assay used for PCA
#' @param ncomponents number of components
#' @param ... additional parameters passed to 'scater::runPCA'
#' @return an updated SingleCellExperiment object with PCA dimension reduction 
#' @export
#' @importFrom scater runPCA
runPCA <- function(object, subset_row = NULL, exprs_values = "logcounts", ncomponents = 50, ...) {
    scater::runPCA(object, subset_row = subset_row, exprs_values = exprs_values, ncomponents = ncomponents, ...)
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
