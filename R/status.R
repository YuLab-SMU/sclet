#' Inspect the current sclet object status
#'
#' @title Status
#' @param object a SingleCellExperiment object
#' @param details logical, whether to include detailed records and command log
#' @return a structured summary of the current active view, available resources,
#'   and recorded analysis context. In addition to `available$analyses`, a light
#'   `health` block reports mainline readiness flags (velocity inputs, fate
#'   probability, perturbation priority / rare cells), the currently active
#'   velocity / trajectory / perturbation records, and `mainline_missing`, the
#'   ordered names of the velocity -> fate -> perturbation stages that are not
#'   yet satisfied.
#' @export
Status <- function(object, details = FALSE) {
    context <- get_analysis_context(object)
    commands <- CommandLog(object, details = isTRUE(details))
    state <- sclet_get_state(object)

    graph_names <- names(state$graphs)
    if (has_graph(object, "knn_graph")) {
        graph_names <- unique(c(graph_names, "knn_graph"))
    }

    assay_names <- SummarizedExperiment::assayNames(object)
    fate_cols <- grep(
        "^cellrank_fate_",
        colnames(SummarizedExperiment::colData(object)),
        value = TRUE
    )

    has_velocity_inputs <- all(c("spliced", "unspliced") %in% assay_names)
    has_vel <- has_velocity(object)
    has_traj <- has_trajectory(object)
    has_cr <- has_cellrank(object)
    has_fate <- length(fate_cols) > 0L
    has_pert <- has_perturbation(object)
    has_prio <- has_perturbation_priority(object)
    has_rare <- has_rare_cells(object)

    mainline <- c(
        velocity_inputs = has_velocity_inputs,
        velocity = has_vel,
        fate = has_fate,
        perturbation = has_pert
    )

    analyses <- c(
        if (has_hvg(object)) "hvg",
        if (has_integration(object)) "integration",
        if (has_annotation(object)) "annotation",
        if (has_mapping(object)) "mapping",
        if (has_traj) "trajectory",
        if (has_cellchat(object)) "communication",
        if (has_milo(object)) "milo",
        if (has_supercell(object)) "aggregation",
        if (has_vel) "velocity",
        if (has_scenic(object)) "scenic",
        if (has_geneset_scoring(object)) "geneset_scoring",
        if (has_cr) "cellrank",
        if (has_spatial(object)) "spatial",
        if (has_pert) "perturbation",
        if (has_prio) "priority",
        if (has_rare) "rare_cells"
    )

    result <- list(
        active = context$active,
        available = list(
            assays = assay_names,
            layers = Layers(object),
            reductions = SingleCellExperiment::reducedDimNames(object),
            graphs = graph_names,
            analyses = analyses
        ),
        health = list(
            has_spliced_assay = "spliced" %in% assay_names,
            has_unspliced_assay = "unspliced" %in% assay_names,
            has_velocity = has_vel,
            has_trajectory = has_traj,
            has_cellrank = has_cr,
            has_fate = has_fate,
            has_perturbation = has_pert,
            has_perturbation_priority = has_prio,
            has_rare_cells = has_rare,
            active_velocity = sclet_get_active_state(object, "velocity"),
            active_trajectory = sclet_get_active_state(object, "trajectory"),
            active_perturbation = sclet_get_active_state(object, "perturbation"),
            mainline_missing = names(mainline)[!mainline]
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
