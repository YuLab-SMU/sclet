#' Run standard single-cell analysis workflow
#' 
#' @description This function wraps the standard preprocessing and dimensional reduction steps into a single pipeline.
#' It sequentially runs NormalizeData, FindVariableFeatures, ScaleData, RunPCA, FindNeighbors, FindClusters, and RunUMAP.
#' 
#' @title RunStandardPipeline
#' @param object a SingleCellExperiment object
#' @param nfeatures number of highly variable features to select, default is 2000
#' @param npcs number of principal components to compute, default is 50
#' @param dims dimensions of PCA to use for neighbor finding and UMAP, default is 1:30
#' @param resolution clustering resolution, default is 0.8
#' @param ... additional parameters passed to individual steps (currently not fully forwarded)
#' @return an updated SingleCellExperiment object
#' @export
RunStandardPipeline <- function(object, nfeatures = 2000, npcs = 50, dims = 1:30, resolution = 0.8, ...) {
    object <- NormalizeData(object)
    object <- FindVariableFeatures(object, nfeatures = nfeatures)
    object <- ScaleData(object)
    object <- RunPCA(object, ncomponents = npcs)
    object <- FindNeighbors(object, dims = dims)
    object <- FindClusters(object, resolution = resolution)
    object <- RunUMAP(object, dims = dims)
    return(object)
}
