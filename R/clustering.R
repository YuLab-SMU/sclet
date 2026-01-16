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
    if (is.null(colLabels(object))) {
        if (!is.null(colData(object)$label)) {
            colLabels(object) <- colData(object)$label
        } else {
            return(NULL)
        }
    }
    setNames(colLabels(object), colnames(object))
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
