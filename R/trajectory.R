#' Run Slingshot trajectory inference
#' 
#' Wrapper for slingshot trajectory inference.
#' 
#' @title RunSlingshot
#' @param object A SingleCellExperiment object
#' @param reduction Dimensional reduction to use (default: "UMAP")
#' @param cluster.labels Cluster labels to use (default: Idents(object))
#' @param start.clus Cluster to start trajectory from
#' @param ... Additional arguments passed to slingshot::slingshot
#' @return A SingleCellExperiment object with slingshot results in metadata and colData
#' @export
RunSlingshot <- function(object, reduction = "UMAP", cluster.labels = NULL, start.clus = NULL, ...) {
    if (!requireNamespace("slingshot", quietly = TRUE)) {
        stop("Package 'slingshot' is needed for this function to work. Please install it.")
    }
    
    if (is.null(cluster.labels)) {
        cluster.labels <- Idents(object)
        if (is.null(cluster.labels)) {
            stop("No cluster labels found. Please run FindClusters or provide cluster.labels.")
        }
    }
    
    # Get reduced dims
    rd <- reducedDim(object, reduction)
    if (is.null(rd)) {
        stop(paste0("Dimensional reduction '", reduction, "' not found."))
    }
    
    # Run slingshot
    sds <- slingshot::slingshot(data = rd, clusterLabels = cluster.labels, start.clus = start.clus, ...)
    
    # Store results
    # We can store the SlingshotDataSet object in metadata
    object@metadata$slingshot <- sds
    
    # And maybe put pseudotime in colData for easy plotting
    pt <- slingshot::slingPseudotime(sds)
    colnames(pt) <- paste0("pseudotime_", colnames(pt))
    colData(object) <- cbind(colData(object), pt)
    
    return(object)
}
