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
#' @param name Annotation/mapping record id. Defaults to `"singler"`.
#' @param ... Additional arguments passed to SingleR::SingleR
#' @return A SingleCellExperiment object with annotation added to colData
#' @export
RunSingleR <- function(object, ref = NULL, labels = NULL, assay.type = NULL,
                       layer = NULL, name = "singler", ...) {
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
    
    labels_col <- if (identical(name, "singler")) "SingleR_labels" else paste0(name, "_labels")
    pruned_labels_col <- if (identical(name, "singler")) "SingleR_pruned.labels" else paste0(name, "_pruned.labels")
    score_col <- paste0(name, "_score")

    # Add predictions to colData
    SummarizedExperiment::colData(object)[[labels_col]] <- pred$labels
    SummarizedExperiment::colData(object)[[pruned_labels_col]] <- pred$pruned.labels
    if ("delta.next" %in% colnames(pred)) {
        SummarizedExperiment::colData(object)[[score_col]] <- pred$delta.next
    } else {
        score_col <- NULL
    }
    
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
        labels_col = labels_col,
        pruned_labels_col = pruned_labels_col,
        score_col = score_col
    )
    state_summary <- list(
        n_labels = length(unique(stats::na.omit(pred$labels))),
        n_pruned_labels = length(unique(stats::na.omit(pred$pruned.labels))),
        mean_score = if (!is.null(score_col)) mean(SummarizedExperiment::colData(object)[[score_col]], na.rm = TRUE) else NULL
    )
    
    object <- sclet_set_analysis_state(
        object = object,
        type = "annotation",
        id = name,
        method = "SingleR",
        inputs = state_inputs,
        artifacts = state_artifacts,
        summary = state_summary,
        active = TRUE
    )
    
    # Also register as mapping for cross-compatibility
    mapping_artifacts <- list(
        labels_col = labels_col,
        pruned_labels_col = pruned_labels_col,
        score_col = score_col,
        mapping_type = "reference_mapping"
    )
    object <- sclet_set_analysis_state(
        object = object,
        type = "mapping",
        id = name,
        method = "SingleR",
        inputs = utils::modifyList(state_inputs, list(mode = "label_transfer")),
        artifacts = mapping_artifacts,
        summary = state_summary,
        active = TRUE
    )
    
    object <- sclet_log_command(
        object,
        "RunSingleR",
        params = list(
            layer = source$layer,
            assay.type = source$assay,
            name = name
        ),
        outputs = list(
            annotation = name,
            mapping = name
        )
    )
    
    message(sprintf(
        "Annotation added to colData columns: '%s' and '%s'",
        labels_col,
        pruned_labels_col
    ))
    
    return(object)
}

#' Run reference mapping through a unified semantic entry point
#'
#' This wrapper keeps "reference mapping" as the stable user-facing concept,
#' while routing to a concrete backend such as `SingleR` or the lightweight
#' KNN mapper.
#'
#' @title RunReferenceMapping
#' @param object A SingleCellExperiment object
#' @param ref Reference dataset
#' @param labels Labels in the reference dataset
#' @param method Backend to use. One of `"SingleR"`, `"KNN"`, or `"Symphony"`.
#' @param assay.type Assay to use for `SingleR`. Ignored for other methods.
#' @param layer Layer to use for mapping.
#' @param features Features to use for `method = "KNN"`. Ignored for other methods.
#' @param k Number of neighbors for `method = "KNN"` or kNN prediction in Symphony.
#' @param vars Character vector of batch variables for Symphony reference building.
#' @param name Mapping record id.
#' @param ... Additional arguments passed to the selected backend.
#' @return A SingleCellExperiment object with mapping results added and recorded.
#' @export
RunReferenceMapping <- function(object, ref, labels = NULL,
                                method = c("SingleR", "KNN", "Symphony"),
                                assay.type = NULL, layer = NULL,
                                features = NULL, k = 5, vars = NULL,
                                name = "mapping", ...) {
    method <- match.arg(method)

    if (identical(method, "SingleR")) {
        return(RunSingleR(
            object = object,
            ref = ref,
            labels = labels,
            assay.type = assay.type,
            layer = layer,
            name = name,
            ...
        ))
    }

    if (identical(method, "Symphony")) {
        return(RunSymphonyMapping(
            object = object,
            ref = ref,
            labels = labels,
            vars = vars,
            k = k,
            layer = layer,
            name = name,
            ...
        ))
    }

    RunKNNPredict(
        object = object,
        ref = ref,
        labels = labels,
        features = features,
        layer = layer,
        k = k,
        name = name,
        ...
    )
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
    
    # Keep matrices sparse to avoid materializing large dense cell-by-gene blocks.
    query_mat <- sclet_extract_cell_feature_matrix(
        object = object,
        assay_name = source$assay,
        features = features
    )
    
    # Find matching assay in ref, prefer logcounts or the same name
    ref_assay <- "logcounts"
    if (!ref_assay %in% SummarizedExperiment::assayNames(ref)) {
        ref_assay <- SummarizedExperiment::assayNames(ref)[1]
    }
    ref_mat <- sclet_extract_cell_feature_matrix(
        object = ref,
        assay_name = ref_assay,
        features = features
    )
    
    # Run KNN
    # By default queryKNN returns indices.
    # BiocNeighbors 2.0.x-2.4.x cannot dispatch on dgCMatrix inputs (see #25);
    # sclet_as_knn_matrix() falls back to ordinary matrices on those versions.
    knn_res <- BiocNeighbors::queryKNN(
        X = sclet_as_knn_matrix(ref_mat),
        query = sclet_as_knn_matrix(query_mat),
        k = k
    )
    
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

#' Run Symphony reference mapping
#'
#' Symphony builds a harmonized reference atlas and then maps query cells
#' into the same space, predicting cell type labels via k-NN.
#'
#' @param object A SingleCellExperiment query object.
#' @param ref A SingleCellExperiment reference object with cell type labels.
#' @param labels Character. Column in `colData(ref)` with reference labels, or a
#'   character vector of labels.
#' @param vars Character vector of batch variables for Symphony reference building.
#' @param k Integer. Number of neighbors for kNN prediction. Defaults to 5.
#' @param layer Character. Expression layer to use. Defaults to "logcounts".
#' @param name Character. Mapping record id. Defaults to `"symphony_map"`.
#' @param ... Additional arguments passed to `symphony::buildReference()`.
#'
#' @return Updated SingleCellExperiment with Symphony mapping in colData.
#' @export
RunSymphonyMapping <- function(object, ref, labels = NULL, vars = NULL,
                               k = 5, layer = NULL, name = "symphony_map", ...) {
    if (!requireNamespace("symphony", quietly = TRUE)) {
        stop(
            "Package 'symphony' is needed for this function to work. ",
            "It is no longer available on CRAN; please install it from ",
            "GitHub with: remotes::install_github('immunogenomics/symphony')",
            call. = FALSE
        )
    }
    if (!is(object, "SingleCellExperiment")) {
        stop("`object` must be a SingleCellExperiment object.")
    }
    if (!is(ref, "SingleCellExperiment")) {
        stop("`ref` must be a SingleCellExperiment object.")
    }

    if (is.null(labels)) {
        if ("label.main" %in% colnames(SummarizedExperiment::colData(ref))) {
            label_col <- "label.main"
        } else if ("label" %in% colnames(SummarizedExperiment::colData(ref))) {
            label_col <- "label"
        } else {
            stop("Please specify `labels` for the reference dataset.")
        }
    } else if (length(labels) == 1 && is.character(labels) &&
               labels %in% colnames(SummarizedExperiment::colData(ref))) {
        label_col <- labels
    } else {
        label_col <- "symphony_label_temp"
        ref[[label_col]] <- as.character(labels)
    }
    ref_labels <- as.character(SummarizedExperiment::colData(ref)[[label_col]])

    common_genes <- intersect(rownames(object), rownames(ref))
    if (length(common_genes) == 0) {
        stop("No common genes between query and reference objects.")
    }
    object <- object[common_genes, , drop = FALSE]
    ref <- ref[common_genes, , drop = FALSE]

    src <- if (is.null(layer)) {
        sclet_resolve_expression_source(
            object, prefer_nonscaled = TRUE,
            fallback_layers = c("logcounts"), context = "Symphony mapping"
        )
    } else {
        sclet_resolve_expression_source(
            object, layer = layer, context = "Symphony mapping"
        )
    }
    ref_src <- if (!is.null(layer) && layer %in% Layers(ref)) {
        sclet_resolve_expression_source(
            ref, layer = layer, context = "Symphony reference"
        )
    } else {
        sclet_resolve_expression_source(
            ref, prefer_nonscaled = TRUE,
            fallback_layers = c("logcounts"), context = "Symphony reference"
        )
    }

    ref_exp <- as.matrix(SummarizedExperiment::assay(ref, ref_src$assay))
    ref_meta <- as.data.frame(SummarizedExperiment::colData(ref))

    cli::cli_alert_info("Building Symphony reference with {.val {ncol(ref_exp)}} cells")
    sym_ref <- symphony::buildReference(
        exp_ref = ref_exp,
        metadata_ref = ref_meta,
        vars = vars,
        verbose = FALSE,
        ...
    )

    query_exp <- as.matrix(SummarizedExperiment::assay(object, src$assay))
    query_meta <- as.data.frame(SummarizedExperiment::colData(object))

    cli::cli_alert_info("Mapping {.val {ncol(query_exp)}} query cells onto reference")
    sym_query <- symphony::mapQuery(
        exp_query = query_exp,
        metadata_query = query_meta,
        ref_obj = sym_ref,
        vars = vars,
        verbose = FALSE,
        do_normalize = FALSE
    )

    cli::cli_alert_info("Predicting labels via k-NN (k = {.val {k}})")
    sym_query <- symphony::knnPredict(
        sym_query,
        sym_ref,
        train_labels = ref_labels,
        k = k,
        save_as = "symphony_predicted",
        confidence = TRUE
    )

    predicted_labels <- sym_query$meta_data$symphony_predicted
    pred_confidence <- sym_query$meta_data$symphony_predicted_confidence
    if (is.null(pred_confidence)) {
        pred_confidence <- rep(NA_real_, length(predicted_labels))
    }

    object$symphony_predicted <- as.character(predicted_labels)
    object$symphony_confidence <- as.numeric(pred_confidence)

    object <- sclet_set_analysis(
        object,
        "mapping",
        list(
            id = name,
            method = "Symphony",
            label_col = "symphony_predicted",
            confidence_col = "symphony_confidence",
            k = k,
            vars = vars,
            n_common_genes = length(common_genes),
            timestamp = Sys.time()
        )
    )
    object <- sclet_set_analysis_state(
        object = object,
        type = "mapping",
        id = name,
        method = "Symphony",
        inputs = list(
            layer = src$layer,
            assay = src$assay,
            vars = vars,
            k = k
        ),
        artifacts = list(
            analysis_key = "mapping",
            label_col = "symphony_predicted",
            confidence_col = "symphony_confidence"
        ),
        summary = list(
            n_query_cells = ncol(object),
            n_ref_cells = ncol(ref),
            n_common_genes = length(common_genes)
        )
    )
    object <- sclet_log_command(
        object,
        "RunReferenceMapping",
        params = list(
            method = "Symphony",
            vars = vars,
            k = k,
            layer = src$layer,
            name = name
        ),
        outputs = list(
            analysis = "mapping",
            state = "mapping"
        )
    )

    object
}

#' Plot query groups versus transferred reference labels
#'
#' This visualization summarizes label transfer results as a heatmap whose rows
#' are query groups (for example clusters) and whose columns are predicted
#' reference labels.
#'
#' @title plot_reference_label_transfer_heatmap
#' @param object A SingleCellExperiment object.
#' @param group.by Query grouping variable. If `NULL`, uses `ActiveIdent(object)`.
#'   The special value `"colLabels"` uses `SingleCellExperiment::colLabels(object)`.
#' @param id Optional mapping record id. If `NULL`, uses the active mapping record.
#' @param labels_col Optional predicted-label column in `colData(object)`. If
#'   `NULL`, resolves it from the selected mapping record.
#' @param normalize One of `"row"`, `"column"`, or `"none"`.
#' @param show_text Logical. If `TRUE`, draws values on the tiles.
#' @param digits Integer. Number of digits to show when normalized values are plotted.
#' @param low,high Colors for the fill gradient.
#' @return A ggplot object.
#' @export
plot_reference_label_transfer_heatmap <- function(
    object,
    group.by = NULL,
    id = NULL,
    labels_col = NULL,
    normalize = c("row", "column", "none"),
    show_text = TRUE,
    digits = 2,
    low = "white",
    high = "#B2182B"
) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is needed for this function to work. Please install it.")
    }

    normalize <- match.arg(normalize)
    mapping_record <- get_mapping(object, id = id)
    if (is.null(mapping_record)) {
        stop("No mapping results found. Please run RunReferenceMapping() first.")
    }

    if (is.null(labels_col)) {
        labels_col <- mapping_record$artifacts$labels_col
    }
    if (is.null(labels_col) || !labels_col %in% colnames(SummarizedExperiment::colData(object))) {
        stop("Could not resolve a valid transferred-label column from the mapping record.")
    }

    if (is.null(group.by)) {
        group.by <- ActiveIdent(object)
    }
    if (is.null(group.by)) {
        stop("No active identity found. Please provide `group.by`.")
    }

    if (identical(group.by, "colLabels")) {
        query_group <- SingleCellExperiment::colLabels(object)
        if (is.null(query_group)) {
            stop("group.by is 'colLabels' but `colLabels(object)` is empty.")
        }
    } else if (group.by %in% colnames(SummarizedExperiment::colData(object))) {
        query_group <- SummarizedExperiment::colData(object)[[group.by]]
    } else {
        stop(sprintf("Grouping column '%s' not found in colData(object).", group.by))
    }

    reference_label <- SummarizedExperiment::colData(object)[[labels_col]]
    keep <- !is.na(query_group) & !is.na(reference_label)
    if (!any(keep)) {
        stop("No non-missing query groups and transferred labels available to plot.")
    }

    df <- data.frame(
        query_group = as.character(query_group[keep]),
        reference_label = as.character(reference_label[keep]),
        stringsAsFactors = FALSE
    )

    heatmap_df <- as.data.frame(
        table(
            query_group = df$query_group,
            reference_label = df$reference_label
        ),
        stringsAsFactors = FALSE
    )

    heatmap_df$query_group <- factor(heatmap_df$query_group, levels = unique(df$query_group))
    heatmap_df$reference_label <- factor(heatmap_df$reference_label, levels = unique(df$reference_label))

    if (identical(normalize, "row")) {
        totals <- ave(heatmap_df$Freq, heatmap_df$query_group, FUN = sum)
        heatmap_df$value <- ifelse(totals > 0, heatmap_df$Freq / totals, NA_real_)
        heatmap_df$label <- formatC(heatmap_df$value, digits = digits, format = "f")
        fill_name <- "Row proportion"
    } else if (identical(normalize, "column")) {
        totals <- ave(heatmap_df$Freq, heatmap_df$reference_label, FUN = sum)
        heatmap_df$value <- ifelse(totals > 0, heatmap_df$Freq / totals, NA_real_)
        heatmap_df$label <- formatC(heatmap_df$value, digits = digits, format = "f")
        fill_name <- "Column proportion"
    } else {
        heatmap_df$value <- heatmap_df$Freq
        heatmap_df$label <- as.character(heatmap_df$Freq)
        fill_name <- "Count"
    }

    p <- ggplot2::ggplot(
        heatmap_df,
        ggplot2::aes(x = .data$reference_label, y = .data$query_group, fill = .data$value)
    ) +
        ggplot2::geom_tile(color = "white", linewidth = 0.3) +
        ggplot2::scale_fill_gradient(low = low, high = high, name = fill_name) +
        ggplot2::labs(
            x = "Predicted reference label",
            y = if (identical(group.by, "colLabels")) "Query group (colLabels)" else group.by
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
            panel.grid = ggplot2::element_blank()
        )

    if (isTRUE(show_text)) {
        p <- p + ggplot2::geom_text(ggplot2::aes(label = .data$label), size = 3)
    }

    p
}

#' Plot mapping confidence by transferred label
#'
#' This plot summarizes the confidence or score associated with transferred
#' reference labels. For `KNN` backends this usually corresponds to the
#' majority-neighbor proportion; for `SingleR` it uses the stored delta-based
#' score when available.
#'
#' @title plot_reference_label_confidence
#' @param object A SingleCellExperiment object.
#' @param id Optional mapping record id. If `NULL`, uses the active mapping record.
#' @param labels_col Optional transferred-label column in `colData(object)`. If
#'   `NULL`, resolves it from the selected mapping record.
#' @param score_col Optional confidence-score column in `colData(object)`. If
#'   `NULL`, resolves it from the selected mapping record.
#' @param type One of `"violin"` or `"boxplot"`.
#' @param sort_by_score Logical. If `TRUE`, orders labels by median score.
#' @param point_size Numeric. Point size for jittered observations.
#' @param jitter_width Numeric. Jitter width for overlaid points.
#' @return A ggplot object.
#' @export
plot_reference_label_confidence <- function(
    object,
    id = NULL,
    labels_col = NULL,
    score_col = NULL,
    type = c("violin", "boxplot"),
    sort_by_score = TRUE,
    point_size = 0.5,
    jitter_width = 0.15
) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is needed for this function to work. Please install it.")
    }

    type <- match.arg(type)
    mapping_record <- get_mapping(object, id = id)
    if (is.null(mapping_record)) {
        stop("No mapping results found. Please run RunReferenceMapping() first.")
    }

    if (is.null(labels_col)) {
        labels_col <- mapping_record$artifacts$labels_col
    }
    if (is.null(score_col)) {
        score_col <- mapping_record$artifacts$score_col
    }
    if (is.null(labels_col) || !labels_col %in% colnames(SummarizedExperiment::colData(object))) {
        stop("Could not resolve a valid transferred-label column from the mapping record.")
    }
    if (is.null(score_col) || !score_col %in% colnames(SummarizedExperiment::colData(object))) {
        stop("Could not resolve a valid confidence-score column from the mapping record.")
    }

    df <- data.frame(
        reference_label = SummarizedExperiment::colData(object)[[labels_col]],
        score = SummarizedExperiment::colData(object)[[score_col]],
        stringsAsFactors = FALSE
    )
    df <- df[!is.na(df$reference_label) & !is.na(df$score), , drop = FALSE]
    if (nrow(df) == 0) {
        stop("No non-missing transferred labels and confidence scores available to plot.")
    }

    if (isTRUE(sort_by_score)) {
        med <- stats::aggregate(score ~ reference_label, data = df, FUN = stats::median)
        med <- med[order(med$score, decreasing = TRUE), , drop = FALSE]
        df$reference_label <- factor(df$reference_label, levels = med$reference_label)
    } else {
        df$reference_label <- factor(df$reference_label, levels = unique(df$reference_label))
    }

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$reference_label, y = .data$score))
    if (identical(type, "violin")) {
        p <- p + ggplot2::geom_violin(
            ggplot2::aes(fill = .data$reference_label),
            color = NA,
            alpha = 0.7,
            scale = "width"
        )
    } else {
        p <- p + ggplot2::geom_boxplot(
            ggplot2::aes(fill = .data$reference_label),
            outlier.shape = NA,
            alpha = 0.7
        )
    }

    p +
        ggplot2::geom_jitter(width = jitter_width, size = point_size, alpha = 0.35) +
        ggplot2::labs(
            x = "Predicted reference label",
            y = "Confidence / score"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
            legend.position = "none",
            panel.grid.minor = ggplot2::element_blank()
        )
}
