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
    # Store in 'slingshot_info' to match sling_plot.R expectations
    object@metadata$slingshot_info <- sds
    
    # Store pseudotime as a single DataFrame column named 'slingPseudotime'
    # This matches the expectation of pseudo_plot and genecurve_plot
    pt <- as.data.frame(slingshot::slingPseudotime(sds))
    colData(object)$slingPseudotime <- pt
    
    return(object)
}
