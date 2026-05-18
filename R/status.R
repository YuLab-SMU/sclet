#' Inspect the current sclet object status
#'
#' @title Status
#' @param object a SingleCellExperiment object
#' @param details logical, whether to include detailed records and command log
#' @return a structured summary of the current active view, available resources,
#'   and recorded analysis context
#' @export
Status <- function(object, details = FALSE) {
    context <- get_analysis_context(object)
    commands <- CommandLog(object, details = isTRUE(details))
    state <- sclet_get_state(object)

    graph_names <- names(state$graphs)
    if (has_graph(object, "knn_graph")) {
        graph_names <- unique(c(graph_names, "knn_graph"))
    }

    analyses <- c(
        if (has_hvg(object)) "hvg",
        if (has_integration(object)) "integration",
        if (has_annotation(object)) "annotation",
        if (has_mapping(object)) "mapping",
        if (has_trajectory(object)) "trajectory",
        if (has_cellchat(object)) "communication",
        if (has_milo(object)) "milo",
        if (has_supercell(object)) "aggregation"
    )

    result <- list(
        active = context$active,
        available = list(
            assays = SummarizedExperiment::assayNames(object),
            layers = Layers(object),
            reductions = SingleCellExperiment::reducedDimNames(object),
            graphs = graph_names,
            analyses = analyses
        ),
        n_commands = nrow(commands),
        last_command = if (nrow(commands)) utils::tail(commands$command, 1) else NULL
    )

    if (isTRUE(details)) {
        result$records <- context$records
        result$hvg <- get_hvg(object)
        result$commands <- commands
    }

    result
}
