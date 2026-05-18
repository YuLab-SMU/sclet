#' Run SingleR for cell type annotation
#' 
#' Wrapper for SingleR to perform automatic cell type annotation.
#' 
#' @title RunSingleR
#' @param object A SingleCellExperiment object
#' @param ref Reference dataset (e.g. from celldex). If NULL, will try to use HumanPrimaryCellAtlasData.
#' @param labels Labels in the reference dataset. If NULL, tries to use 'label.main'.
#' @param assay.type Assay to use for annotation. If NULL, resolved from `layer`.
#' @param layer Layer to use for annotation. If NULL, use `DefaultLayer(object)`.
#' @param ... Additional arguments passed to SingleR::SingleR
#' @return A SingleCellExperiment object with annotation added to colData
#' @export
RunSingleR <- function(object, ref = NULL, labels = NULL, assay.type = NULL, layer = NULL, ...) {
    if (!requireNamespace("SingleR", quietly = TRUE)) {
        stop("Package 'SingleR' is needed for this function to work. Please install it.")
    }
    
    if (is.null(ref)) {
        if (!requireNamespace("celldex", quietly = TRUE)) {
            stop("Package 'celldex' is needed to load default reference. Please install it.")
        }
        message("Loading HumanPrimaryCellAtlasData from celldex...")
        ref <- celldex::HumanPrimaryCellAtlasData()
    }
    
    if (is.null(labels)) {
        if ("label.main" %in% colnames(SummarizedExperiment::colData(ref))) {
            labels <- ref$label.main
        } else {
            stop("Please specify 'labels' for the reference dataset.")
        }
    }
    
    source <- sclet_resolve_expression_source(
        object = object,
        layer = layer,
        assay = assay.type,
        prefer_nonscaled = TRUE,
        context = "SingleR"
    )
    
    pred <- SingleR::SingleR(test = object, ref = ref, labels = labels, assay.type.test = source$assay, ...)
    
    # Add predictions to colData
    SummarizedExperiment::colData(object)$SingleR_labels <- pred$labels
    SummarizedExperiment::colData(object)$SingleR_pruned.labels <- pred$pruned.labels
    
    # Check for integration dependency
    integration_id <- NULL
    integration_dep <- sclet_resolve_integration_dependency(object, layer = source$layer)
    if (!is.null(integration_dep)) {
        integration_id <- integration_dep$id
    }
    
    state_inputs <- list(
        assay = source$assay,
        layer = source$layer,
        reference_class = class(ref)[[1]]
    )
    
    if (!is.null(integration_id)) {
        state_inputs$integration <- list(id = integration_id)
    }
    
    state_artifacts <- list(
        labels_col = "SingleR_labels",
        pruned_labels_col = "SingleR_pruned.labels"
    )
    
    object <- sclet_set_analysis_state(
        object = object,
        type = "annotation",
        id = "singler",
        method = "SingleR",
        inputs = state_inputs,
        artifacts = state_artifacts,
        active = TRUE
    )
    
    # Also register as mapping for cross-compatibility
    mapping_artifacts <- list(
        labels_col = "SingleR_labels",
        mapping_type = "reference_mapping"
    )
    object <- sclet_set_analysis_state(
        object = object,
        type = "mapping",
        id = "singler",
        method = "SingleR",
        inputs = utils::modifyList(state_inputs, list(mode = "label_transfer")),
        artifacts = mapping_artifacts,
        active = TRUE
    )
    
    object <- sclet_log_command(
        object,
        "RunSingleR",
        params = list(
            layer = source$layer,
            assay.type = source$assay
        ),
        outputs = list(
            annotation = "singler",
            mapping = "singler"
        )
    )
    
    message("Annotation added to colData columns: 'SingleR_labels' and 'SingleR_pruned.labels'")
    
    return(object)
}

#' Run KNN-based reference mapping for label transfer
#' 
#' A lightweight alternative to SingleR that maps cells to a reference dataset
#' using KNN in a shared feature space (e.g., highly variable genes).
#' 
#' @title RunKNNPredict
#' @param object A SingleCellExperiment object
#' @param ref Reference dataset (SingleCellExperiment or similar with colData)
#' @param labels Labels in the reference dataset to transfer. If NULL, tries to use 'label.main'.
#' @param features Features to use for mapping. If NULL, uses the intersection of 
#' highly variable features from `object` and available features in `ref`.
#' @param layer Layer to use for mapping. If NULL, uses `DefaultLayer(object)`.
#' @param k Number of nearest neighbors to use for prediction.
#' @param name Mapping record id. Defaults to `"knn_map"`.
#' @param ... Additional arguments
#' @return A SingleCellExperiment object with predicted labels added to colData
#' @export
RunKNNPredict <- function(object, ref, labels = NULL, features = NULL, layer = NULL, k = 5, name = "knn_map", ...) {
    if (!requireNamespace("BiocNeighbors", quietly = TRUE)) {
        stop("Package 'BiocNeighbors' is needed for this function to work. Please install it.")
    }
    
    if (is.null(labels)) {
        if ("label.main" %in% colnames(SummarizedExperiment::colData(ref))) {
            labels <- ref$label.main
        } else if ("label" %in% colnames(SummarizedExperiment::colData(ref))) {
            labels <- ref$label
        } else {
            stop("Please specify 'labels' for the reference dataset.")
        }
    } else {
        if (length(labels) == 1 && is.character(labels) && labels %in% colnames(SummarizedExperiment::colData(ref))) {
            labels <- SummarizedExperiment::colData(ref)[[labels]]
        }
    }
    
    source <- sclet_resolve_expression_source(
        object = object,
        layer = layer,
        assay = NULL,
        prefer_nonscaled = FALSE,
        context = "KNN Mapping"
    )
    
    if (is.null(features)) {
        features <- VariableFeatures(object)
        if (length(features) == 0) {
            features <- rownames(object)
        }
        features <- intersect(features, rownames(ref))
    }
    if (length(features) == 0) {
        stop("No overlapping features found between object and reference.")
    }
    
    # Extract data matrices
    # For a simple KNN, we need a dense or sparse matrix. We'll transpose for BiocNeighbors
    query_mat <- t(as.matrix(SummarizedExperiment::assay(object, source$assay)[features, , drop = FALSE]))
    
    # Find matching assay in ref, prefer logcounts or the same name
    ref_assay <- "logcounts"
    if (!ref_assay %in% SummarizedExperiment::assayNames(ref)) {
        ref_assay <- SummarizedExperiment::assayNames(ref)[1]
    }
    ref_mat <- t(as.matrix(SummarizedExperiment::assay(ref, ref_assay)[features, , drop = FALSE]))
    
    # Run KNN
    # By default queryKNN returns indices.
    knn_res <- BiocNeighbors::queryKNN(X = ref_mat, query = query_mat, k = k)
    
    # Predict label by majority vote among K neighbors
    pred_labels <- apply(knn_res$index, 1, function(idx) {
        neighbor_labels <- labels[idx]
        # Get mode
        tab <- table(neighbor_labels)
        names(tab)[which.max(tab)]
    })
    
    # Calculate a simple confidence score: proportion of the majority label among K neighbors
    conf_scores <- apply(knn_res$index, 1, function(idx) {
        neighbor_labels <- labels[idx]
        tab <- table(neighbor_labels)
        max(tab) / k
    })
    
    col_name <- paste0(name, "_labels")
    conf_col_name <- paste0(name, "_score")
    SummarizedExperiment::colData(object)[[col_name]] <- pred_labels
    SummarizedExperiment::colData(object)[[conf_col_name]] <- conf_scores
    
    # Register state
    state_inputs <- list(
        assay = source$assay,
        layer = source$layer,
        features = length(features),
        k = k,
        reference_class = class(ref)[[1]]
    )
    
    state_artifacts <- list(
        labels_col = col_name,
        score_col = conf_col_name,
        mapping_type = "knn_reference_mapping"
    )
    
    state_summary <- list(
        n_labels = length(unique(pred_labels)),
        n_reference_labels = length(unique(labels))
    )
    
    object <- sclet_set_analysis_state(
        object = object,
        type = "mapping",
        id = name,
        method = "KNNPredict",
        inputs = utils::modifyList(state_inputs, list(mode = "label_transfer")),
        artifacts = state_artifacts,
        summary = state_summary,
        active = TRUE
    )
    
    object <- sclet_log_command(
        object,
        "RunKNNPredict",
        params = list(
            layer = source$layer,
            k = k,
            name = name
        ),
        outputs = list(
            mapping = name
        )
    )
    
    message(sprintf("KNN mapping added to colData column: '%s'", col_name))
    return(object)
}
