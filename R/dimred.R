#' elbow plot 
#' 
#' @title ElbowPlot
#' @param object a SingleCellExperiment object
#' @return elbow plot
#' @export
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 theme_classic
#' @importFrom rlang .data
ElbowPlot <- function(object) {
    pca_results <- reducedDim(object, "PCA")

    v <- attr(pca_results, "percentVar")
    d <- data.frame(PC = seq_along(v), var = v)

    ggplot(d, aes(.data$PC, .data$var)) + 
        geom_point() +
        ylab("Standard Deviation") +
        theme_classic()      
}

#' get or set the default dimensional reduction
#'
#' @title DefaultReduction
#' @param object a SingleCellExperiment object
#' @param value reduction name
#' @return current default reduction or updated object
#' @export
DefaultReduction <- function(object) {
    state <- S4Vectors::metadata(object)$sclet
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

#' @rdname DefaultReduction
#' @export
"DefaultReduction<-" <- function(object, value) {
    if (!value %in% SingleCellExperiment::reducedDimNames(object)) {
        stop("Reduction '", value, "' not found in reducedDims(object).")
    }
    sclet_set_active_reduction(object, value)
}

#' get or set the default assay
#'
#' @title DefaultAssay
#' @param object a SingleCellExperiment object
#' @param value assay name
#' @return current default assay or updated object
#' @export
DefaultAssay <- function(object) {
    state <- S4Vectors::metadata(object)$sclet
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

#' @rdname DefaultAssay
#' @export
"DefaultAssay<-" <- function(object, value) {
    if (!value %in% SummarizedExperiment::assayNames(object)) {
        stop("Assay '", value, "' not found in assays(object).")
    }
    sclet_set_active_assay(object, value)
}

#' get or set the default graph
#'
#' @title DefaultGraph
#' @param object a SingleCellExperiment object
#' @param value graph name
#' @return current default graph or updated object
#' @export
DefaultGraph <- function(object) {
    state <- S4Vectors::metadata(object)$sclet
    graph <- state$active$graph
    if (!is.null(graph)) {
        return(graph)
    }
    if (!is.null(sclet_get_legacy_graph_entry(object, "knn_graph")$object)) {
        return("knn_graph")
    }
    if (length(state$graphs) > 0) {
        return(names(state$graphs)[[1]])
    }
    NULL
}

#' @rdname DefaultGraph
#' @export
"DefaultGraph<-" <- function(object, value) {
    graph_names <- names(sclet_get_state(object)$graphs)
    legacy_graphs <- if (!is.null(sclet_get_legacy_graph_entry(object, "knn_graph"))) "knn_graph" else character()
    valid_graphs <- unique(c(graph_names, legacy_graphs))
    if (!value %in% valid_graphs) {
        stop("Graph '", value, "' not found in sclet graph registry.")
    }
    sclet_set_active_graph(object, value)
}

#' get or set the active identity field
#'
#' @title ActiveIdent
#' @param object a SingleCellExperiment object
#' @param value identity source, e.g. `colLabels` or a `colData` column name
#' @return current active identity source or updated object
#' @export
ActiveIdent <- function(object) {
    state <- S4Vectors::metadata(object)$sclet
    ident <- state$active$ident
    if (!is.null(ident)) {
        return(ident)
    }
    if (!is.null(SingleCellExperiment::colLabels(object))) {
        return("colLabels")
    }
    NULL
}

#' @rdname ActiveIdent
#' @export
"ActiveIdent<-" <- function(object, value) {
    valid <- c("colLabels", colnames(SummarizedExperiment::colData(object)))
    if (!value %in% valid) {
        stop("Identity source '", value, "' not found.")
    }
    sclet_set_active_ident(object, value)
}

#' inspect sclet command log
#'
#' @title CommandLog
#' @param object a SingleCellExperiment object
#' @param details logical, whether to include raw `params` and `outputs` list columns
#' @return a data.frame summarizing recorded commands
#' @export
CommandLog <- function(object, details = FALSE) {
    commands <- sclet_get_commands(object)
    if (length(commands) == 0) {
        out <- data.frame(
            command = character(),
            timestamp = as.POSIXct(character()),
            params_summary = character(),
            outputs_summary = character(),
            stringsAsFactors = FALSE
        )
        if (isTRUE(details)) {
            out$params <- I(list())
            out$outputs <- I(list())
        }
        return(out)
    }
    summarize_field <- function(x) {
        if (is.null(x) || length(x) == 0) {
            return("")
        }
        parts <- vapply(seq_along(x), function(i) {
            key <- names(x)[i]
            value <- x[[i]]
            value_txt <- if (length(value) == 1 && !is.list(value)) {
                as.character(value)
            } else if (is.null(dim(value))) {
                paste0("[", length(value), "]")
            } else {
                paste0(dim(value)[1], "x", dim(value)[2])
            }
            paste0(key, "=", value_txt)
        }, character(1))
        paste(parts, collapse = "; ")
    }
    out <- data.frame(
        command = vapply(commands, `[[`, character(1), "command"),
        timestamp = as.POSIXct(vapply(commands, function(x) as.character(x$timestamp), character(1))),
        params_summary = vapply(commands, function(x) summarize_field(x$params), character(1)),
        outputs_summary = vapply(commands, function(x) summarize_field(x$outputs), character(1)),
        stringsAsFactors = FALSE
    )
    if (isTRUE(details)) {
        out$params <- I(lapply(commands, `[[`, "params"))
        out$outputs <- I(lapply(commands, `[[`, "outputs"))
    }
    out
}

#' @importFrom SingleCellExperiment 'reducedDim<-'
set_dimred <- function(object, dims, reduction = "PCA") {
    if (!reduction %in% SingleCellExperiment::reducedDimNames(object)) {
        stop(paste("Reduction", reduction, "not found in object."))
    }
    reducedDim(object, ".dimred") <- reducedDim(object, reduction)[, dims, drop = FALSE]
    return(object)
}

#' Run PCA
#' 
#' @title RunPCA
#' @param object a SingleCellExperiment object
#' @param subset_row subset of rows used for PCA
#' @param exprs_values assay used for PCA. If NULL, sclet resolves it from the
#' selected `layer`.
#' @param layer layer used for PCA. If NULL, use `DefaultLayer(object)`.
#' @param ncomponents number of components
#' @param ... additional parameters passed to 'scater::runPCA'
#' @return an updated SingleCellExperiment object with PCA dimension reduction 
#' @export
#' @importFrom scater runPCA
RunPCA <- function(object, subset_row = NULL, exprs_values = NULL, layer = NULL, ncomponents = 50, ...) {
    prev_state <- sclet_get_state(object)
    source <- sclet_resolve_expression_source(
        object = object,
        layer = layer,
        assay = exprs_values,
        context = "PCA"
    )
    integration <- sclet_resolve_integration_dependency(object, layer = source$layer)
    exprs_values <- source$assay
    object <- sclet_muffle_known_warnings(
        scater::runPCA(
            object,
            subset_row = subset_row,
            exprs_values = exprs_values,
            ncomponents = ncomponents,
            ...
        ),
        patterns = c(
            "You're computing too large a percentage of total singular values, use a standard svd instead.",
            "more singular values/vectors requested than available"
        )
    )
    object <- sclet_restore_state(object, prev_state)
    object <- sclet_set_active_reduction(object, "PCA")
    object <- sclet_set_analysis_state(
        object = object,
        type = "reduction",
        id = "pca",
        method = "scater::runPCA",
        inputs = list(
            layer = source$layer,
            assay = exprs_values,
            subset_row = subset_row,
            integration = integration
        ),
        artifacts = list(
            reduction = "PCA"
        ),
        params = list(
            ncomponents = ncomponents
        ),
        summary = list(
            n_components = ncol(SingleCellExperiment::reducedDim(object, "PCA"))
        )
    )
    object <- sclet_log_command(
        object,
        "RunPCA",
        params = list(
            layer = source$layer,
            subset_row = subset_row,
            exprs_values = exprs_values,
            ncomponents = ncomponents
        )
    )
    object
}

#' run umap
#' 
#' @title RunUMAP
#' @param object a SingleCellExperiment object
#' @param dims dimensions used in UMAP. If NULL, use all dimensions from the
#' selected source reduction.
#' @param reduction source reduction used to build UMAP. If NULL, sclet first
#' tries `DefaultReduction(object)`, then falls back to a non-UMAP reduction.
#' @param layer source layer annotation used for provenance. If NULL, sclet
#' tries to infer it from the source reduction state, then falls back to
#' `DefaultLayer(object)`.
#' @return updated SingleCellExperiment object with UMAP dimension reduction 
#' @export
#' @importFrom scater runUMAP
RunUMAP <- function(object, dims = NULL, reduction = NULL, layer = NULL) {
    prev_state <- sclet_get_state(object)
    available_reductions <- SingleCellExperiment::reducedDimNames(object)
    if (is.null(reduction)) {
        reduction <- DefaultReduction(object)
    }
    if (is.null(reduction) || identical(reduction, "UMAP")) {
        fallback <- setdiff(c("PCA", available_reductions), "UMAP")
        reduction <- if (length(fallback) > 0) fallback[[1]] else NULL
    }
    if (is.null(reduction) || !reduction %in% available_reductions) {
        stop("No valid source reduction available for UMAP. Please run PCA or specify `reduction`.")
    }
    if (is.null(dims)) {
        dims <- seq_len(ncol(SingleCellExperiment::reducedDim(object, reduction)))
    }
    reducedDim(object, ".dimred") <- reducedDim(object, reduction)[, dims, drop = FALSE]
    object <- scater::runUMAP(object, dimred = '.dimred')
    object <- sclet_restore_state(object, prev_state)
    object <- sclet_set_active_reduction(object, "UMAP")
    reduction_id <- tolower(reduction)
    reduction_state <- sclet_get_state_record(object, "reduction", reduction_id)
    if (is.null(layer) && !is.null(reduction_state) && !is.null(reduction_state$inputs$layer)) {
        layer <- reduction_state$inputs$layer
    }
    if (is.null(layer)) {
        layer <- DefaultLayer(object)
    }
    integration <- NULL
    if (!is.null(reduction_state) && !is.null(reduction_state$inputs$integration)) {
        integration <- reduction_state$inputs$integration
    }
    if (is.null(integration)) {
        integration <- sclet_resolve_integration_dependency(object, layer = layer)
    }
    object <- sclet_set_analysis_state(
        object = object,
        type = "reduction",
        id = "umap",
        method = "scater::runUMAP",
        inputs = list(
            reduction = reduction,
            layer = layer,
            dims = dims,
            integration = integration
        ),
        artifacts = list(
            reduction = "UMAP"
        ),
        summary = list(
            n_components = ncol(SingleCellExperiment::reducedDim(object, "UMAP"))
        )
    )
    object <- sclet_log_command(
        object,
        "RunUMAP",
        params = list(dims = dims, reduction = reduction, layer = layer)
    )
    object
}

#' Run Diffusion Map for complex trajectory inference
#'
#' @title RunDiffusionMap
#' @param sce a SingleCellExperiment object
#' @param n_components number of diffusion components to compute. Defaults to 50.
#' @param reduction the name of the dimensional reduction to use as input (e.g. 'PCA'). Defaults to NULL (uses the active reduction).
#' @param ... additional arguments passed to \code{destiny::DiffusionMap}
#' @return an updated SingleCellExperiment object with a 'DM' reduction
#' @export
#' @importFrom SingleCellExperiment reducedDim<- reducedDim
#' @importFrom SingleCellExperiment reducedDimNames
RunDiffusionMap <- function(sce, n_components = 50, reduction = NULL, ...) {
    if (!requireNamespace("destiny", quietly = TRUE)) {
        stop("Please install 'destiny' via BiocManager::install('destiny')")
    }
    
    if (is.null(reduction)) {
        reduction <- sclet_get_active_reduction(sce)
    }
    if (is.null(reduction) || !reduction %in% SingleCellExperiment::reducedDimNames(sce)) {
        stop("No valid reduction found for DiffusionMap. Please run RunPCA first or specify a valid 'reduction'.")
    }
    
    rd <- SingleCellExperiment::reducedDim(sce, reduction)
    
    # Run destiny::DiffusionMap
    # destiny uses the data matrix to compute distances. Since we provide PCA, 
    # we can use the PCA matrix as the 'data' for DiffusionMap
    dm <- destiny::DiffusionMap(rd, n_pcs = NA, ...)
    
    # Extract the eigenvectors
    dm_coords <- dm@eigenvectors[, seq_len(min(n_components, ncol(dm@eigenvectors))), drop = FALSE]
    rownames(dm_coords) <- rownames(rd)
    colnames(dm_coords) <- paste0("DC_", seq_len(ncol(dm_coords)))
    
    # Store in SCE
    SingleCellExperiment::reducedDim(sce, "DM") <- dm_coords
    
    # Log and update state
    sce <- sclet_log_command(
        sce,
        "RunDiffusionMap",
        params = list(n_components = n_components, reduction = reduction, ...)
    )
    
    sce <- sclet_set_active_reduction(sce, "DM")
    
    return(sce)
}
