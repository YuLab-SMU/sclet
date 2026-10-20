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
    # Hand the dimension subset straight to makeSNNGraph. Writing it into a
    # scratch reduction first (`set_dimred()` did: reducedDim(object, ".dimred"))
    # left that scratch on the object permanently -- after a plain PBMC 10k
    # pipeline on 2026-10-19, reducedDimNames() was c("PCA", ".dimred") and the
    # stray entry was an 11,769 x 30 duplicate of PCA[, 1:30], visible to every
    # later step that enumerates reductions. The subset is a local matrix.
    dimred <- SingleCellExperiment::reducedDim(object, source_reduction)[, dims, drop = FALSE]
    graph <- bluster::makeSNNGraph(
        dimred,
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
#' @param seed integer seed for the Louvain clustering. `igraph::cluster_louvain()`
#'   consumes random numbers, so the same graph and resolution otherwise return
#'   different labels on each call; the seed is recorded with the clustering
#'   state so a plan's clusters can actually be reproduced.
#' @return updated SingleCellExperiment object with clustering memberships
#' @importFrom igraph cluster_louvain
# @importFrom SingleCellExperiment 'colLabels<-'
#' @export
FindClusters <- function(object, resolution = 0.5, seed = 1L) {
    graph_name <- DefaultGraph(object)
    if (is.null(graph_name)) {
        graph_name <- "knn_graph"
    }
    g <- get_graph(object, graph_name, "object")
    if (is.null(g)) {
        stop("No KNN graph found. Please run FindNeighbors() first.")
    }
    graph_state <- sclet_get_state_record(object, "graph", graph_name)
    if (is.null(seed) || !length(seed) || is.na(seed[[1L]])) {
        seed <- 1L
    }
    seed <- as.integer(seed[[1L]])
    # Seed locally and restore afterwards: a clustering call must not move the
    # caller's random stream. Without any seed (the behaviour until
    # 2026-10-19) two runs of the same plan returned 13 and 12 clusters on
    # PBMC 10k, so "re-run the plan and get the same clusters" was not true.
    had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    previous_seed <- if (had_seed) get(".Random.seed", envir = globalenv(), inherits = FALSE) else NULL
    on.exit({
        if (had_seed) {
            assign(".Random.seed", previous_seed, envir = globalenv())
        } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
            rm(".Random.seed", envir = globalenv())
        }
    }, add = TRUE)
    set.seed(seed)
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
            resolution = resolution,
            seed = seed
        ),
        summary = list(
            n_clusters = length(unique(clusters$membership))
        )
    )
    object <- sclet_log_command(
        object,
        "FindClusters",
        params = list(resolution = resolution, seed = seed)
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
    if (is.null(active_ident)) {
        return(NULL)
    }
    if (identical(active_ident, "colLabels")) {
        if (is.null(colLabels(object))) {
            return(NULL)
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
