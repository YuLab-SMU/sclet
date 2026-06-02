sclet_state_template <- function() {
    list(
        schema_version = "0.1.2",
        active = list(),
        layers = list(),
        states = list(
            active = list(),
            records = list()
        ),
        features = list(),
        graphs = list(),
        analyses = list(),
        commands = list()
    )
}

sclet_modify_list <- function(x, val) {
    if (!is.list(x) || !is.list(val)) {
        return(val)
    }
    vnames <- names(val)
    vnames <- vnames[nzchar(vnames)]
    for (v in vnames) {
        xv <- x[[v]]
        vv <- val[[v]]
        if (is.list(xv) && is.list(vv) && 
            (is.null(class(xv)) || identical(class(xv), "list")) &&
            (is.null(class(vv)) || identical(class(vv), "list"))) {
            x[[v]] <- sclet_modify_list(xv, vv)
        } else {
            x[[v]] <- vv
        }
    }
    x
}

sclet_merge_state <- function(base, override = NULL) {
    if (is.null(override)) {
        return(base)
    }
    merged <- sclet_modify_list(base, override)
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
        RunSingleR = list(analysis = "annotation", state = "annotation"),
        RunSlingshot = list(analysis = "trajectory"),
        RunSlingshot_trajectory = list(analysis = "trajectory"),
        RunVelocity = list(analysis = "velocity", state = "velocity"),
        RunCellRank = list(analysis = "cellrank", state = "trajectory"),
        RunCellFate = list(analysis = "cellrank", state = "trajectory"),
        RunSuperCell = list(analysis = "supercell"),
        RunMilo = list(analysis = "milo"),
        runMilo = list(analysis = "milo"),
        RunCellChat = list(analysis = "cellchat"),
        BatchRemover = list(analysis = "batch", layer = "corrected", state = "integration"),
        RunPerturbationPriority = list(analysis = "augur", state = "priority"),
        RunRareCellDetection = list(analysis = "rare_cells", state = "rare_cells"),
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

sclet_internal_metadata_keys <- function() {
    c(
        "sclet",
        "hvgcols",
        "hvgmethod",
        "nVariableFeatures",
        "knn_graph",
        "slingshot_info",
        "milo_results",
        "SuperCell",
        "batch_correction"
    )
}

sclet_strip_internal_metadata <- function(metadata) {
    if (is.null(metadata) || length(metadata) == 0) {
        return(list())
    }
    metadata[setdiff(names(metadata), sclet_internal_metadata_keys())]
}

sclet_merge_external_metadata <- function(...) {
    inputs <- list(...)
    merged <- list()
    for (x in inputs) {
        if (inherits(x, "SingleCellExperiment")) {
            x <- S4Vectors::metadata(x)
        }
        merged <- sclet_modify_list(merged, sclet_strip_internal_metadata(x))
    }
    merged
}

sclet_rebuild_internal_state <- function(object, hvg = NULL, commands = NULL, active_assay = NULL) {
    external_metadata <- sclet_strip_internal_metadata(S4Vectors::metadata(object))
    state <- sclet_state_template()
    if (!is.null(hvg)) {
        state$features$hvg <- hvg
    }
    if (!is.null(commands)) {
        state$commands <- commands
    }
    if (!is.null(active_assay) && active_assay %in% SummarizedExperiment::assayNames(object)) {
        state$active$assay <- active_assay
    }
    external_metadata$sclet <- state
    S4Vectors::metadata(object) <- external_metadata
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

sclet_get_layer_registry <- function(object) {
    state <- sclet_get_state(object)
    layers <- state$layers
    if (is.null(layers)) {
        layers <- list()
    }
    assay_names <- SummarizedExperiment::assayNames(object)
    for (assay_name in assay_names) {
        if (is.null(layers[[assay_name]])) {
            layers[[assay_name]] <- list(
                assay = assay_name,
                role = assay_name,
                params = list()
            )
        }
    }
    layers
}

sclet_set_layer <- function(object, name, assay = NULL, role = NULL, params = list(), active = TRUE) {
    object <- sclet_init_state(object)
    md <- S4Vectors::metadata(object)
    state <- sclet_get_state(object)

    if (is.null(assay)) {
        assay <- name
    }
    if (!assay %in% SummarizedExperiment::assayNames(object)) {
        stop("Assay '", assay, "' not found in assays(object).")
    }
    if (is.null(role)) {
        role <- name
    }

    state$layers[[name]] <- list(
        assay = assay,
        role = role,
        params = params
    )
    if (isTRUE(active)) {
        state$active$layer <- name
    }
    md$sclet <- state
    S4Vectors::metadata(object) <- md
    object
}

sclet_get_active_layer <- function(object) {
    state <- sclet_get_state(object)
    layer <- state$active$layer
    if (!is.null(layer)) {
        return(layer)
    }
    layers <- sclet_get_layer_registry(object)
    for (candidate in c("scaled", "logcounts", "corrected", "reconstructed", "counts")) {
        if (candidate %in% names(layers)) {
            return(candidate)
        }
    }
    active_assay <- sclet_get_active_assay(object)
    if (!is.null(active_assay)) {
        return(active_assay)
    }
    if (length(layers) > 0) {
        return(names(layers)[[1]])
    }
    NULL
}

sclet_set_active_layer <- function(object, layer) {
    object <- sclet_init_state(object)
    md <- S4Vectors::metadata(object)
    state <- sclet_get_state(object)
    layers <- sclet_get_layer_registry(object)
    if (!layer %in% names(layers)) {
        stop("Layer '", layer, "' not found in sclet layer registry.")
    }
    state$active$layer <- layer
    md$sclet <- state
    S4Vectors::metadata(object) <- md
    object
}

sclet_get_layer <- function(object, name = NULL) {
    layers <- sclet_get_layer_registry(object)
    if (is.null(name)) {
        name <- sclet_get_active_layer(object)
    }
    layers[[name]]
}

sclet_resolve_layer <- function(object, layer = NULL) {
    layers <- sclet_get_layer_registry(object)
    if (is.null(layer)) {
        layer <- sclet_get_active_layer(object)
    }
    if (is.null(layer)) {
        return(NULL)
    }
    entry <- layers[[layer]]
    if (!is.null(entry)) {
        return(list(
            layer = layer,
            assay = entry$assay,
            role = entry$role,
            params = entry$params
        ))
    }
    if (layer %in% SummarizedExperiment::assayNames(object)) {
        return(list(
            layer = layer,
            assay = layer,
            role = layer,
            params = list()
        ))
    }
    stop("Layer '", layer, "' not found.")
}

sclet_resolve_expression_source <- function(
    object,
    layer = NULL,
    assay = NULL,
    prefer_nonscaled = FALSE,
    fallback_layers = c("logcounts"),
    context = "analysis"
) {
    assay_names <- SummarizedExperiment::assayNames(object)
    if (!is.null(assay)) {
        if (!assay %in% assay_names) {
            stop("Assay '", assay, "' not found in object.")
        }
        return(list(
            layer = assay,
            assay = assay,
            role = assay,
            params = list()
        ))
    }

    resolved <- sclet_resolve_layer(object, layer = layer)
    if (is.null(resolved)) {
        stop("No valid layer found in object.")
    }

    if (isTRUE(prefer_nonscaled) && identical(resolved$role, "scaled")) {
        layer_names <- Layers(object)
        for (candidate in fallback_layers) {
            if (candidate %in% layer_names) {
                fallback <- sclet_resolve_layer(object, candidate)
                message(
                    "Using '", fallback$assay, "' for ", context,
                    " instead of active layer '", resolved$layer, "'."
                )
                return(fallback)
            }
        }
    }

    resolved
}

sclet_resolve_integration_dependency <- function(object, layer = NULL) {
    integration_id <- sclet_get_active_state(object, "integration")
    if (is.null(integration_id)) {
        return(NULL)
    }
    integration <- sclet_get_state_record(object, "integration", integration_id)
    if (is.null(integration)) {
        return(NULL)
    }
    if (is.null(layer)) {
        layer <- sclet_get_active_layer(object)
    }
    artifact_layer <- integration$artifacts$layer
    if (!is.null(layer) && !is.null(artifact_layer) && !identical(layer, artifact_layer)) {
        return(NULL)
    }
    list(
        type = "integration",
        id = integration_id,
        method = integration$method,
        layer = artifact_layer,
        assay = integration$artifacts$assay
    )
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
    if (any(c("bio", "total") %in% rd_names)) {
        return("scran")
    }
    if ("variance.standardized" %in% rd_names) {
        return("seurat")
    }
    NULL
}

sclet_resolve_hvg_method <- function(object, method = NULL, allow_legacy = FALSE) {
    if (!is.null(method)) {
        return(match.arg(method, c("seurat", "scran", "scrapper")))
    }

    inferred <- sclet_get_hvg_method(object)
    if (!is.null(inferred) && !identical(inferred, "seurat")) {
        return(inferred)
    }

    rd_names <- colnames(SummarizedExperiment::rowData(object))
    if (any(c("bio", "total") %in% rd_names)) {
        return("scran")
    }

    if (isTRUE(allow_legacy) && ("variance.standardized" %in% rd_names || identical(inferred, "seurat"))) {
        return("seurat")
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

sclet_get_legacy_graph_entry <- function(object, name) {
    if (!identical(name, "knn_graph")) {
        return(NULL)
    }
    legacy <- S4Vectors::metadata(object)$knn_graph
    if (is.null(legacy)) {
        return(NULL)
    }
    list(object = legacy)
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
    legacy_entry <- sclet_get_legacy_graph_entry(object, name)
    if (!is.null(legacy_entry$object)) {
        return(legacy_entry$object)
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
    if (!is.null(sclet_get_legacy_graph_entry(object, "knn_graph"))) {
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

sclet_get_legacy_analysis_record <- function(object, type) {
    md <- S4Vectors::metadata(object)
    switch(
        type,
        trajectory = {
            legacy <- md$slingshot_info
            if (is.null(legacy)) {
                return(NULL)
            }
            list(
                method = "slingshot",
                dataset = legacy
            )
        },
        milo = md$milo_results,
        aggregation = {
            legacy <- md$SuperCell
            if (is.null(legacy)) {
                return(NULL)
            }
            list(object = legacy)
        },
        batch = md$batch_correction,
        NULL
    )
}

sclet_state_types <- function() {
    c(
        "preprocess",
        "reduction",
        "graph",
        "clustering",
        "annotation",
        "integration",
        "trajectory",
        "aggregation",
        "milo",
        "da",
        "communication",
        "mapping",
        "detest",
        "enrichment",
        "velocity",
        "scenic",
        "geneset_scoring",
        "spatial",
        "perturbation",
        "priority",
        "rare_cells"
    )
}

sclet_validate_state_type <- function(type) {
    if (!is.character(type) || length(type) != 1 || is.na(type) || !nzchar(type)) {
        stop("`type` must be a single non-empty character string.")
    }
    if (!type %in% sclet_state_types()) {
        stop(
            "Unsupported state type '", type, "'. Supported types are: ",
            paste(sclet_state_types(), collapse = ", "), "."
        )
    }
    invisible(type)
}

sclet_get_state_records <- function(object, type = NULL) {
    state <- sclet_get_state(object)
    records <- state$states$records
    if (is.null(records)) {
        records <- list()
    }
    if (is.null(type)) {
        return(records)
    }
    sclet_validate_state_type(type)
    type_records <- records[[type]]
    if (is.null(type_records)) {
        return(list())
    }
    type_records
}

sclet_get_active_state <- function(object, type = NULL) {
    state <- sclet_get_state(object)
    active <- state$states$active
    if (is.null(active)) {
        active <- list()
    }
    if (is.null(type)) {
        return(active)
    }
    sclet_validate_state_type(type)
    active[[type]]
}

sclet_set_active_state <- function(object, type, id) {
    sclet_validate_state_type(type)
    object <- sclet_init_state(object)
    state <- sclet_get_state(object)
    if (!is.null(id) && is.null(state$states$records[[type]][[id]])) {
        stop("State record '", id, "' does not exist for type '", type, "'.")
    }
    state$states$active[[type]] <- id
    sclet_set_state(object, state)
}

sclet_get_state_record <- function(object, type, id = NULL) {
    sclet_validate_state_type(type)
    records <- sclet_get_state_records(object, type)
    if (is.null(id)) {
        id <- sclet_get_active_state(object, type)
    }
    if (!is.null(id)) {
        return(records[[id]])
    }
    if (length(records) == 1) {
        return(records[[1]])
    }
    NULL
}

sclet_set_state_record <- function(object, type, value, id = NULL, active = TRUE) {
    sclet_validate_state_type(type)
    if (!is.list(value)) {
        stop("`value` must be a list.")
    }

    object <- sclet_init_state(object)
    state <- sclet_get_state(object)
    records <- state$states$records
    if (is.null(records)) {
        records <- list()
    }
    type_records <- records[[type]]
    if (is.null(type_records)) {
        type_records <- list()
    }

    if (is.null(id)) {
        id <- value$id
    }
    if (is.null(id) || !nzchar(id)) {
        id <- paste0(type, "_", length(type_records) + 1L)
    }

    record <- utils::modifyList(
        list(
            id = id,
            type = type,
            status = "completed"
        ),
        value
    )
    type_records[[id]] <- record
    records[[type]] <- type_records
    state$states$records <- records
    if (isTRUE(active)) {
        state$states$active[[type]] <- id
    }
    sclet_set_state(object, state)
}

sclet_set_analysis_state <- function(
    object,
    type,
    id,
    method,
    inputs = list(),
    artifacts = list(),
    params = list(),
    summary = list(),
    active = TRUE
) {
    sclet_set_state_record(
        object = object,
        type = type,
        id = id,
        active = active,
        value = list(
            method = method,
            inputs = inputs,
            artifacts = artifacts,
            params = params,
            summary = summary,
            created_at = Sys.time()
        )
    )
}
