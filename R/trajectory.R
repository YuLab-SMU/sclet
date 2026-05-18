RunSlingshot_trajectory <- function(object, reduction = "UMAP", cluster.labels = NULL, start.clus = NULL, name = "slingshot_trajectory", ...) {
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
    
    prev_state <- sclet_get_state(object)
    pt <- as.data.frame(slingshot::slingPseudotime(sds))
    colData(object)$slingPseudotime <- pt
    object <- sclet_restore_state(object, prev_state)
    object <- sclet_set_analysis(
        object,
        "trajectory",
        list(
            method = "slingshot",
            reduction = reduction,
            start_cluster = start.clus,
            dataset = sds
        )
    )
    object <- sclet_set_analysis_state(
        object = object,
        type = "trajectory",
        id = name,
        method = "slingshot",
        inputs = list(
            reduction = reduction,
            cluster_labels = if (is.character(cluster.labels) && length(cluster.labels) == 1) cluster.labels else NULL
        ),
        artifacts = list(
            analysis_key = "trajectory",
            pseudotime_col = "slingPseudotime"
        ),
        params = list(
            start.clus = start.clus
        ),
        summary = list(
            n_lineages = ncol(pt)
        ),
        active = TRUE
    )
    object <- sclet_log_command(
        object,
        "RunSlingshot_trajectory",
        params = list(
            reduction = reduction,
            start.clus = start.clus,
            name = name
        ),
        outputs = list(
            analysis = "trajectory",
            trajectory = name
        )
    )
    
    return(object)
}
