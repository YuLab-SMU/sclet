#' Run the standard basic single-cell workflow
#'
#' A beginner-oriented convenience wrapper that chains the standard sclet
#' preprocessing, dimensionality reduction, graph and clustering steps in the
#' conventional order. It is purely deterministic: no AI call is made and no
#' action registry is involved. This is the "do not think about the details"
#' entry point; expert users should call the individual functions directly.
#'
#' Every step delegates to the corresponding existing sclet function, which is
#' responsible for its own state registration. This function only sequences
#' the calls and never re-implements any analysis logic.
#'
#' @title RunBasicWorkflow
#' @param object a SingleCellExperiment object.
#' @param n_features number of highly variable features requested by
#'   `FindVariableFeatures()`.
#' @param n_pcs number of principal components requested by `RunPCA()` and the
#'   number of leading dimensions used by `FindNeighbors()` and `RunUMAP()`.
#' @param cluster_resolution clustering resolution passed to `FindClusters()`.
#' @param steps Character vector naming the steps to run, in execution order.
#'   Supported names are `"normalize"`, `"variable_features"`, `"scale"`,
#'   `"pca"`, `"neighbors"`, `"clusters"` and `"umap"`.
#' @param verbose Logical. Report progress with `cli::cli_alert_info()`, matching
#'   the logging style used elsewhere in the package.
#' @return the updated SingleCellExperiment object.
#' @export
RunBasicWorkflow <- function(
    object,
    n_features = 2000,
    n_pcs = 30,
    cluster_resolution = 0.5,
    steps = c("normalize", "variable_features", "scale", "pca", "neighbors", "clusters", "umap"),
    verbose = TRUE) {
    if (!inherits(object, "SingleCellExperiment")) {
        stop("object must be a SingleCellExperiment", call. = FALSE)
    }
    supported <- c("normalize", "variable_features", "scale", "pca",
        "neighbors", "clusters", "umap")
    steps <- unique(as.character(steps))
    steps <- steps[nzchar(steps)]
    unknown <- setdiff(steps, supported)
    if (length(unknown)) {
        stop("unknown step(s): ", paste(unknown, collapse = ", "),
            ". Supported steps: ", paste(supported, collapse = ", "), call. = FALSE)
    }
    steps <- supported[supported %in% steps]
    if (!length(steps)) {
        stop("`steps` must name at least one step", call. = FALSE)
    }
    n_features <- as.numeric(n_features)
    n_pcs <- as.numeric(n_pcs)
    if (length(n_features) != 1L || !is.finite(n_features) || n_features < 1) {
        stop("n_features must be one positive number", call. = FALSE)
    }
    if (length(n_pcs) != 1L || !is.finite(n_pcs) || n_pcs < 1) {
        stop("n_pcs must be one positive number", call. = FALSE)
    }
    # a principal component count can never exceed what the data supports
    n_pcs <- max(1L, min(as.integer(n_pcs), nrow(object), ncol(object)))

    announce <- function(step, label) {
        if (isTRUE(verbose)) {
            cli::cli_alert_info("Running basic workflow step {.val {step}}: {label}")
        }
    }

    if ("normalize" %in% steps) {
        announce("normalize", "NormalizeData()")
        object <- NormalizeData(object)
    }
    if ("variable_features" %in% steps) {
        announce("variable_features", "FindVariableFeatures()")
        object <- FindVariableFeatures(object, nfeatures = n_features)
    }
    if ("scale" %in% steps) {
        announce("scale", "ScaleData()")
        object <- ScaleData(object)
    }
    if ("pca" %in% steps) {
        announce("pca", "RunPCA()")
        object <- RunPCA(object, ncomponents = n_pcs)
    }
    # neighbors, clusters and UMAP all consume leading dimensions of a reduction
    use_dims <- if ("PCA" %in% SingleCellExperiment::reducedDimNames(object)) {
        seq_len(min(n_pcs, ncol(SingleCellExperiment::reducedDim(object, "PCA"))))
    } else {
        NULL
    }
    if ("neighbors" %in% steps) {
        if (is.null(use_dims)) {
            stop("step 'neighbors' requires a PCA reduction; include the 'pca' step", call. = FALSE)
        }
        announce("neighbors", "FindNeighbors()")
        object <- FindNeighbors(object, dims = use_dims, reduction = "PCA")
    }
    if ("clusters" %in% steps) {
        announce("clusters", "FindClusters()")
        object <- FindClusters(object, resolution = cluster_resolution)
    }
    if ("umap" %in% steps) {
        if (is.null(use_dims)) {
            stop("step 'umap' requires a PCA reduction; include the 'pca' step", call. = FALSE)
        }
        announce("umap", "RunUMAP()")
        object <- RunUMAP(object, dims = use_dims, reduction = "PCA")
    }
    object
}
