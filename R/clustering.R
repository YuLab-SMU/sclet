#' build KNN graph
#' 
#' @title FindNeighbors
#' @param object a SingleCellExperiment object
#' @param dims number of dimensions to be used to build the KNN graph
#' @param reduction reduction to use. If NULL, use `DefaultReduction(object)` or fallback to "PCA".
#' @param k number of neighbors
#' @return updated SingleCellExperiment object with KNN graph
#' @importFrom bluster makeSNNGraph
#' @export
FindNeighbors <- function(object, dims, reduction = NULL, k = 10) {
    if (is.null(reduction)) {
        reduction <- DefaultReduction(object)
    }
    if (is.null(reduction)) {
        reduction <- "PCA"
    }
    source_reduction <- reduction
    reduction_id <- tolower(reduction)
    reduction_state <- sclet_get_state_record(object, "reduction", reduction_id)
    integration <- NULL
    layer <- NULL
    if (!is.null(reduction_state)) {
        layer <- reduction_state$inputs$layer
        integration <- reduction_state$inputs$integration
    }
    object <- set_dimred(object, dims, reduction = source_reduction)    
    graph <- bluster::makeSNNGraph(
        SingleCellExperiment::reducedDim(object, ".dimred"),
        k = k,
        type = "rank"
    )
    object <- sclet_set_graph(
        object,
        graph = graph,
        name = "knn_graph",
        params = list(dims = dims, reduction = source_reduction, k = k, type = "rank")
    )
    object <- sclet_set_analysis_state(
        object = object,
        type = "graph",
        id = "knn_graph",
        method = "bluster::makeSNNGraph",
        inputs = list(
            dims = dims,
            reduction = source_reduction,
            layer = layer,
            integration = integration
        ),
        artifacts = list(
            graph = "knn_graph"
        ),
        params = list(
            k = k,
            type = "rank"
        ),
        summary = list(
            n_vertices = igraph::vcount(graph),
            n_edges = igraph::ecount(graph)
        )
    )
    object <- sclet_log_command(
        object,
        "FindNeighbors",
        params = list(dims = dims, reduction = source_reduction, k = k)
    )
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
FindClusters <- function(object, resolution = 0.5) {
    graph_name <- DefaultGraph(object)
    if (is.null(graph_name)) {
        graph_name <- "knn_graph"
    }
    g <- get_graph(object, graph_name, "object")
    if (is.null(g)) {
        stop("No KNN graph found. Please run FindNeighbors() first.")
    }
    graph_state <- sclet_get_state_record(object, "graph", graph_name)
    clusters <- igraph::cluster_louvain(g, resolution=resolution)
    
    SingleCellExperiment::colLabels(object) <- factor(clusters$membership)
    object <- sclet_set_active_ident(object, "colLabels")
    object <- sclet_set_analysis_state(
        object = object,
        type = "clustering",
        id = "louvain_clusters",
        method = "igraph::cluster_louvain",
        inputs = list(
            graph = graph_name,
            reduction = if (!is.null(graph_state)) graph_state$inputs$reduction else NULL,
            layer = if (!is.null(graph_state)) graph_state$inputs$layer else NULL,
            integration = if (!is.null(graph_state)) graph_state$inputs$integration else NULL
        ),
        artifacts = list(
            ident = "colLabels"
        ),
        params = list(
            resolution = resolution
        ),
        summary = list(
            n_clusters = length(unique(clusters$membership))
        )
    )
    object <- sclet_log_command(
        object,
        "FindClusters",
        params = list(resolution = resolution)
    )
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
    active_ident <- ActiveIdent(object)
    if (is.null(active_ident) || identical(active_ident, "colLabels")) {
        if (is.null(colLabels(object))) {
            if (!is.null(colData(object)$label)) {
                colLabels(object) <- colData(object)$label
            } else {
                return(NULL)
            }
        }
        return(setNames(colLabels(object), colnames(object)))
    }
    values <- colData(object)[[active_ident]]
    if (is.null(values)) {
        return(NULL)
    }
    setNames(values, colnames(object))
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

    if (is.null(old_ids)) {
        stop("No identities found in 'object'. Please run FindClusters() first.")
    }

    if (is.null(names(new_ids))) {
        if (length(new_ids) == 1) {
            labels <- rep(new_ids, length(lv))
        } else {
            if (length(new_ids) != length(lv)) {
                stop("`new_ids` must have length equal to number of identity levels, or be a named vector mapping old->new.")
            }
            labels <- new_ids
        }
    } else {
        labels <- lv
        matched <- intersect(names(new_ids), lv)
        labels[match(matched, lv)] <- unname(new_ids[matched])
    }

    colLabels(object) <- factor(
        old_ids, 
        levels = lv, 
        labels = labels
    )
    object <- sclet_set_active_ident(object, "colLabels")
    
    return(object)
}
