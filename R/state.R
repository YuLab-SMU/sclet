sclet_state_template <- function() {
    list(
        schema_version = "0.1.0",
        active = list(),
        features = list(),
        graphs = list(),
        analyses = list(),
        commands = list()
    )
}

sclet_merge_state <- function(base, override = NULL) {
    if (is.null(override)) {
        return(base)
    }
    merged <- utils::modifyList(base, override)
    if (!is.null(override$commands)) {
        merged$commands <- override$commands
    }
    merged
}

sclet_default_command_outputs <- function(command) {
    switch(
        command,
        NormalizeData = list(assay = "logcounts"),
        FindVariableFeatures = list(feature_set = "hvg"),
        ScaleData = list(assay = "scaled"),
        RunPCA = list(reduction = "PCA"),
        runPCA = list(reduction = "PCA"),
        RunUMAP = list(reduction = "UMAP"),
        FindNeighbors = list(graph = "knn_graph"),
        FindClusters = list(ident = "colLabels"),
        RunSlingshot = list(analysis = "trajectory"),
        RunSlingshot_trajectory = list(analysis = "trajectory"),
        RunMilo = list(analysis = "milo"),
        runMilo = list(analysis = "milo"),
        RunCellChat = list(analysis = "cellchat"),
        BatchRemover = list(analysis = "batch"),
        list()
    )
}

sclet_init_state <- function(object) {
    md <- S4Vectors::metadata(object)
    if (is.null(md$sclet)) {
        md$sclet <- sclet_state_template()
        S4Vectors::metadata(object) <- md
    }
    object
}

sclet_get_state <- function(object) {
    md <- S4Vectors::metadata(object)
    state <- md$sclet
    if (is.null(state)) {
        return(sclet_state_template())
    }
    sclet_merge_state(sclet_state_template(), state)
}

sclet_set_state <- function(object, state) {
    md <- S4Vectors::metadata(object)
    md$sclet <- sclet_merge_state(sclet_state_template(), state)
    S4Vectors::metadata(object) <- md
    object
}

sclet_restore_state <- function(object, previous_state = NULL) {
    md <- S4Vectors::metadata(object)
    current_state <- md$sclet
    if (is.null(previous_state)) {
        previous_state <- sclet_state_template()
    }
    if (is.null(current_state)) {
        current_state <- list()
    }
    md$sclet <- sclet_merge_state(sclet_state_template(), sclet_merge_state(previous_state, current_state))
    S4Vectors::metadata(object) <- md
    object
}

sclet_set_active_reduction <- function(object, reduction) {
    object <- sclet_init_state(object)
    md <- S4Vectors::metadata(object)
    state <- sclet_get_state(object)
    state$active$reduction <- reduction
    md$sclet <- state
    S4Vectors::metadata(object) <- md
    object
}

sclet_set_active_assay <- function(object, assay) {
    object <- sclet_init_state(object)
    md <- S4Vectors::metadata(object)
    state <- sclet_get_state(object)
    state$active$assay <- assay
    md$sclet <- state
    S4Vectors::metadata(object) <- md
    object
}

sclet_get_active_assay <- function(object) {
    state <- sclet_get_state(object)
    assay <- state$active$assay
    if (!is.null(assay)) {
        return(assay)
    }
    assay_names <- SummarizedExperiment::assayNames(object)
    for (candidate in c("scaled", "logcounts", "reconstructed", "counts")) {
        if (candidate %in% assay_names) {
            return(candidate)
        }
    }
    if (length(assay_names) > 0) {
        return(assay_names[[1]])
    }
    NULL
}

sclet_get_active_reduction <- function(object) {
    state <- sclet_get_state(object)
    reduction <- state$active$reduction
    if (!is.null(reduction)) {
        return(reduction)
    }
    rd_names <- SingleCellExperiment::reducedDimNames(object)
    for (candidate in c("UMAP", "PCA", "tSNE")) {
        if (candidate %in% rd_names) {
            return(candidate)
        }
    }
    if (length(rd_names) > 0) {
        return(rd_names[[1]])
    }
    NULL
}

sclet_set_active_ident <- function(object, ident) {
    object <- sclet_init_state(object)
    md <- S4Vectors::metadata(object)
    state <- sclet_get_state(object)
    state$active$ident <- ident
    md$sclet <- state
    S4Vectors::metadata(object) <- md
    object
}

sclet_get_active_ident <- function(object) {
    state <- sclet_get_state(object)
    ident <- state$active$ident
    if (!is.null(ident)) {
        return(ident)
    }
    if (!is.null(SingleCellExperiment::colLabels(object))) {
        return("colLabels")
    }
    cd_names <- colnames(SummarizedExperiment::colData(object))
    for (candidate in c("label", "cluster", "ident")) {
        if (candidate %in% cd_names) {
            return(candidate)
        }
    }
    NULL
}

sclet_log_command <- function(object, command, params = list(), outputs = NULL) {
    object <- sclet_init_state(object)
    md <- S4Vectors::metadata(object)
    state <- sclet_get_state(object)
    if (is.null(outputs)) {
        outputs <- sclet_default_command_outputs(command)
    }
    entry <- list(
        command = command,
        params = params,
        outputs = outputs,
        timestamp = Sys.time()
    )
    state$commands <- c(state$commands, list(entry))
    md$sclet <- state
    S4Vectors::metadata(object) <- md
    object
}

sclet_get_commands <- function(object) {
    state <- sclet_get_state(object)
    cmds <- state$commands
    if (is.null(cmds)) {
        return(list())
    }
    cmds
}

sclet_get_hvg_cols <- function(object) {
    state <- sclet_get_state(object)
    cols <- state$features$hvg$rowData_cols
    if (!is.null(cols)) {
        return(cols)
    }
    md <- S4Vectors::metadata(object)
    if (!is.null(md$hvgcols)) {
        return(md$hvgcols)
    }
    rd_names <- colnames(SummarizedExperiment::rowData(object))
    if (all(c("mean", "variance", "variance.expected", "variance.standardized") %in% rd_names)) {
        return(c("mean", "variance", "variance.expected", "variance.standardized"))
    }
    cols <- intersect(c("mean", "total", "tech", "bio", "p.value", "FDR"), rd_names)
    if (length(cols) == 0) {
        return(NULL)
    }
    cols
}

sclet_get_hvg_method <- function(object) {
    state <- sclet_get_state(object)
    method <- state$features$hvg$method
    if (!is.null(method)) {
        return(method)
    }
    md <- S4Vectors::metadata(object)
    if (!is.null(md$hvgmethod)) {
        return(md$hvgmethod)
    }
    rd_names <- colnames(SummarizedExperiment::rowData(object))
    if ("variance.standardized" %in% rd_names) {
        return("seurat")
    }
    if (any(c("bio", "total") %in% rd_names)) {
        return("scran")
    }
    NULL
}

sclet_get_hvg_nfeatures <- function(object) {
    state <- sclet_get_state(object)
    nfeatures <- state$features$hvg$n
    if (!is.null(nfeatures)) {
        return(nfeatures)
    }
    md <- S4Vectors::metadata(object)
    md$nVariableFeatures
}

sclet_set_hvg_state <- function(object, nfeatures = NULL, method = NULL, hvgcols = NULL, selected = NULL) {
    object <- sclet_init_state(object)
    md <- S4Vectors::metadata(object)
    state <- sclet_get_state(object)
    current <- state$features$hvg
    if (is.null(current)) {
        current <- list()
    }
    if (!is.null(nfeatures)) {
        current$n <- nfeatures
    }
    if (!is.null(method)) {
        current$method <- method
    }
    if (!is.null(hvgcols)) {
        current$rowData_cols <- hvgcols
    }
    if (!is.null(selected)) {
        current$selected <- selected
    }
    state$features$hvg <- current
    md$sclet <- state
    S4Vectors::metadata(object) <- md
    object
}

sclet_set_graph <- function(object, graph, name = "knn_graph", params = list(), active = TRUE) {
    object <- sclet_init_state(object)
    md <- S4Vectors::metadata(object)
    state <- sclet_get_state(object)
    state$graphs[[name]] <- list(
        object = graph,
        params = params
    )
    if (isTRUE(active)) {
        state$active$graph <- name
    }
    md$sclet <- state
    S4Vectors::metadata(object) <- md
    object
}

sclet_get_graph <- function(object, name = NULL) {
    state <- sclet_get_state(object)
    if (is.null(name)) {
        name <- state$active$graph
    }
    if (!is.null(name)) {
        graph_entry <- state$graphs[[name]]
        if (!is.null(graph_entry$object)) {
            return(graph_entry$object)
        }
    }
    md <- S4Vectors::metadata(object)
    if (!is.null(md$knn_graph)) {
        return(md$knn_graph)
    }
    NULL
}

sclet_set_active_graph <- function(object, name) {
    object <- sclet_init_state(object)
    md <- S4Vectors::metadata(object)
    state <- sclet_get_state(object)
    state$active$graph <- name
    md$sclet <- state
    S4Vectors::metadata(object) <- md
    object
}

sclet_get_active_graph <- function(object) {
    state <- sclet_get_state(object)
    graph <- state$active$graph
    if (!is.null(graph)) {
        return(graph)
    }
    if (!is.null(S4Vectors::metadata(object)$knn_graph)) {
        return("knn_graph")
    }
    if (length(state$graphs) > 0) {
        return(names(state$graphs)[[1]])
    }
    NULL
}

sclet_set_analysis <- function(object, name, value) {
    object <- sclet_init_state(object)
    md <- S4Vectors::metadata(object)
    state <- sclet_get_state(object)
    state$analyses[[name]] <- value
    md$sclet <- state
    S4Vectors::metadata(object) <- md
    object
}

sclet_get_analysis <- function(object, name, legacy_key = NULL) {
    state <- sclet_get_state(object)
    value <- state$analyses[[name]]
    if (!is.null(value)) {
        return(value)
    }
    if (is.null(legacy_key)) {
        return(NULL)
    }
    md <- S4Vectors::metadata(object)
    legacy_value <- md[[legacy_key]]
    if (is.null(legacy_value)) {
        return(NULL)
    }
    if (is.list(legacy_value)) {
        return(legacy_value)
    }
    list(dataset = legacy_value)
}
