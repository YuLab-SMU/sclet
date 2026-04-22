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
    sclet_get_active_reduction(object)
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
    sclet_get_active_assay(object)
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
    sclet_get_active_graph(object)
}

#' @rdname DefaultGraph
#' @export
"DefaultGraph<-" <- function(object, value) {
    graph_names <- names(sclet_get_state(object)$graphs)
    legacy_graphs <- if (!is.null(S4Vectors::metadata(object)$knn_graph)) "knn_graph" else character()
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
    sclet_get_active_ident(object)
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
set_dimred <- function(object, dims) {
    reducedDim(object, ".dimred") <- reducedDim(object, "PCA")[, dims]
    return(object)
}

#' Run PCA
#' 
#' @title RunPCA
#' @param object a SingleCellExperiment object
#' @param subset_row subset of rows used for PCA
#' @param exprs_values assay used for PCA
#' @param ncomponents number of components
#' @param ... additional parameters passed to 'scater::runPCA'
#' @return an updated SingleCellExperiment object with PCA dimension reduction 
#' @export
#' @importFrom scater runPCA
RunPCA <- function(object, subset_row = NULL, exprs_values = "logcounts", ncomponents = 50, ...) {
    prev_state <- sclet_get_state(object)
    object <- scater::runPCA(
        object,
        subset_row = subset_row,
        exprs_values = exprs_values,
        ncomponents = ncomponents,
        ...
    )
    object <- sclet_restore_state(object, prev_state)
    object <- sclet_set_active_reduction(object, "PCA")
    object <- sclet_log_command(
        object,
        "RunPCA",
        params = list(
            subset_row = subset_row,
            exprs_values = exprs_values,
            ncomponents = ncomponents
        )
    )
    object
}

#' @rdname RunPCA
#' @export
runPCA <- function(object, subset_row = NULL, exprs_values = "logcounts", ncomponents = 50, ...) {
    RunPCA(
        object = object,
        subset_row = subset_row,
        exprs_values = exprs_values,
        ncomponents = ncomponents,
        ...
    )
}

#' run umap
#' 
#' @title RunUMAP
#' @param object a SingleCellExperiment object
#' @param dims dimensions used in UMAP. If NULL, use all dimensions from the
#' selected source reduction.
#' @param reduction source reduction used to build UMAP. If NULL, sclet first
#' tries `DefaultReduction(object)`, then falls back to a non-UMAP reduction.
#' @return updated SingleCellExperiment object with UMAP dimension reduction 
#' @export
#' @importFrom scater runUMAP
RunUMAP <- function(object, dims = NULL, reduction = NULL) {
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
    object <- sclet_log_command(
        object,
        "RunUMAP",
        params = list(dims = dims, reduction = reduction)
    )
    object
}
