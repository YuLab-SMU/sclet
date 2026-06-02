#' Run perturbation priority analysis with Augur
#'
#' Prioritize cell types by their responsiveness to experimental perturbation
#' using the Augur framework. Augur trains a machine-learning model for each
#' cell type to predict experimental condition from molecular measurements, then
#' ranks cell types by cross-validated AUC.
#'
#' @param sce A SingleCellExperiment object.
#' @param condition_col Character. Column in `colData(sce)` that contains
#'   condition labels (e.g. treatment vs control).
#' @param group_col Character. Column in `colData(sce)` that contains cell
#'   type or cluster labels. If `NULL`, defaults to `ActiveIdent(sce)`.
#'   The special value `"colLabels"` uses `SingleCellExperiment::colLabels(sce)`.
#' @param layer Character. Expression layer to use for the Augur analysis.
#'   If `NULL`, uses the active layer and avoids scaled layers.
#' @param assay Character. Assay to use if `layer` is not specified.
#' @param name Character. Analysis record id. Defaults to `"augur"`.
#' @param ... Additional arguments passed to `Augur::calculate_auc()`.
#'
#' @return A `SingleCellExperiment` object with perturbation priority results
#'   registered in the analysis-state layer. The raw Augur result list is also
#'   stored under `metadata(sce)$sclet$analyses$augur` for backwards
#'   compatibility.
#' @export
RunPerturbationPriority <- function(sce, condition_col, group_col = NULL,
    layer = NULL, assay = NULL, name = "augur", ...) {
    if (!requireNamespace("Augur", quietly = TRUE)) {
        stop(
            "Package 'Augur' is needed for this function to work. ",
            "Please install it with: ",
            'remotes::install_github("neurorestore/Augur")',
            call. = FALSE
        )
    }
    if (!is(sce, "SingleCellExperiment")) {
        stop("`sce` must be a SingleCellExperiment object.")
    }
    if (is.null(condition_col) || !nzchar(condition_col) ||
        !condition_col %in% colnames(SummarizedExperiment::colData(sce))) {
        stop(
            "A valid `condition_col` must be provided that exists in colData(sce)."
        )
    }
    cd <- SummarizedExperiment::colData(sce)

    if (is.null(group_col)) {
        group_col <- ActiveIdent(sce)
    }
    if (!is.null(group_col) && identical(group_col, "colLabels")) {
        labels <- SingleCellExperiment::colLabels(sce)
        sce$sclet_priority_group_temp <- as.character(labels)
        group_col <- "sclet_priority_group_temp"
    }
    if (is.null(group_col) || !nzchar(group_col) ||
        !group_col %in% colnames(cd)) {
        stop(
            "A valid `group_col` must be provided that exists in colData(sce)."
        )
    }
    condition_values <- unique(as.character(cd[[condition_col]]))
    if (length(condition_values) < 2L) {
        stop("`condition_col` must contain at least two distinct condition levels.")
    }

    src <- sclet_resolve_expression_source(
        sce, layer = layer, assay = assay,
        prefer_nonscaled = TRUE, fallback_layers = c("logcounts"),
        context = "perturbation priority"
    )

    expr_matrix <- SummarizedExperiment::assay(sce, src$assay)
    meta <- as.data.frame(cd)

    cli::cli_alert_info("Running Augur perturbation priority on {.val {group_col}} by {.val {condition_col}}")
    augur_result <- Augur::calculate_auc(
        input = expr_matrix,
        meta = meta,
        cell_type_col = group_col,
        label_col = condition_col,
        ...
    )

    sce <- sclet_set_analysis(
        sce,
        "augur",
        list(
            id = name,
            condition_col = condition_col,
            group_col = group_col,
            assay = src$assay,
            AUC = augur_result$AUC,
            result = augur_result,
            timestamp = Sys.time()
        )
    )
    sce <- sclet_set_analysis_state(
        object = sce,
        type = "priority",
        id = name,
        method = "Augur",
        inputs = list(
            condition_col = condition_col,
            group_col = group_col,
            assay = src$assay,
            layer = src$layer
        ),
        artifacts = list(
            analysis_key = "augur",
            AUC = augur_result$AUC
        ),
        params = list(...)
    )
    sce <- sclet_log_command(
        sce,
        "RunPerturbationPriority",
        params = list(
            condition_col = condition_col,
            group_col = group_col,
            assay = src$assay,
            name = name
        ),
        outputs = list(
            analysis = "augur",
            state = "priority"
        )
    )

    sce
}

#' Plot perturbation priority ranking from Augur results
#'
#' @param object A SingleCellExperiment object with Augur priority results.
#' @param n Integer. Number of top cell types to display. Defaults to 20.
#' @param fill_color Character. Bar fill color. Defaults to `"steelblue"`.
#'
#' @return A ggplot object ranking cell types by perturbation AUC.
#' @export
plot_perturbation_ranking <- function(object, n = 20, fill_color = "steelblue") {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is needed for this function to work. Please install it.")
    }

    priority <- get_perturbation_priority(object)
    if (is.null(priority)) {
        stop(
            "No perturbation priority results found. ",
            "Run RunPerturbationPriority() or RunStatePriorityWorkflow() first."
        )
    }
    auc_table <- priority$AUC
    if (is.null(auc_table)) {
        stop("No AUC table found in perturbation priority results.")
    }

    cell_type_col <- if ("cell_type" %in% colnames(auc_table)) "cell_type" else colnames(auc_table)[1]
    auc_col <- if ("auc" %in% colnames(auc_table)) "auc" else colnames(auc_table)[2]

    df <- as.data.frame(auc_table)
    df <- df[order(df[[auc_col]], decreasing = TRUE), , drop = FALSE]
    df <- utils::head(df, n)
    df[[cell_type_col]] <- factor(df[[cell_type_col]], levels = rev(df[[cell_type_col]]))

    ggplot2::ggplot(df, ggplot2::aes(x = .data[[cell_type_col]], y = .data[[auc_col]])) +
        ggplot2::geom_col(fill = fill_color, width = 0.7) +
        ggplot2::coord_flip() +
        ggplot2::labs(
            x = "Cell Type",
            y = "Perturbation AUC",
            title = "Perturbation Priority Ranking"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold", size = 13)
        )
}

#' Run rare cell detection
#'
#' Identifies rare and major cell populations using the RareQ algorithm.
#' RareQ operates on the kNN graph, quantifying local topological structure
#' via a neighborhood connectivity metric (Q). High-Q cells represent locally
#' dense micro-clusters suggestive of rare states.
#'
#' @param sce A SingleCellExperiment object with PCA reduction.
#' @param method Character. Rare cell detection backend. Currently only
#'   `"RareQ"` is supported.
#' @param reduction Character. Reduction containing PCA. Defaults to "PCA".
#' @param dims Integer vector. PCA dimensions to use. Defaults to 1:50.
#' @param k_param Integer. k parameter for kNN graph construction. Defaults to 20.
#' @param k Integer. Neighborhood size for Q value computation in RareQ.
#'   Smaller values enhance sensitivity to extremely rare populations.
#'   Defaults to 6.
#' @param Q_cut Numeric. Minimum average Q for a cluster to be retained as a
#'   rare group. Range 0.5-1.0. Defaults to 0.6.
#' @param ratio Numeric. Threshold for merging clusters. Defaults to 0.2.
#' @param max_iter Integer. Maximum iterations for label propagation. Defaults to 100.
#' @param name Analysis record id. Defaults to `"rareq"`.
#' @param ... Additional arguments passed to `RareQ::FindRare()`.
#'
#' @return Updated SingleCellExperiment with `rare_cluster` in colData.
#' @export
RunRareCellDetection <- function(sce, method = c("RareQ"),
    reduction = "PCA", dims = 1:50, k_param = 20, k = 6,
    Q_cut = 0.6, ratio = 0.2, max_iter = 100, name = "rareq", ...) {
    method <- match.arg(method)
    stopifnot(is(sce, "SingleCellExperiment"))
    stopifnot(is.character(name) && nzchar(name))

    if (identical(method, "RareQ")) {
        if (!requireNamespace("RareQ", quietly = TRUE)) {
            stop(
                "Package 'RareQ' is needed for rare cell detection. ",
                "Please install it from GitHub. ",
                'See: https://xiaolab-xjtu.github.io/RareQ/',
                call. = FALSE
            )
        }
        if (!requireNamespace("Seurat", quietly = TRUE)) {
            stop(
                "Package 'Seurat' is needed for SCE-to-Seurat conversion ",
                "during RareQ analysis. Please install it.",
                call. = FALSE
            )
        }

        if (!reduction %in% SingleCellExperiment::reducedDimNames(sce)) {
            stop(sprintf(
                "Reduction '%s' not found in reducedDimNames(sce). ",
                reduction
            ), "Please run RunPCA() first with at least 50 PCs.")
        }

        pca_mat <- SingleCellExperiment::reducedDim(sce, reduction)
        max_dim <- min(max(dims), ncol(pca_mat))
        dims <- seq_len(max_dim)

        seurat_obj <- SeuratObject::CreateSeuratObject(
            counts = SummarizedExperiment::assay(sce, "counts"),
            meta.data = as.data.frame(SummarizedExperiment::colData(sce))
        )
        colnames(pca_mat) <- paste0("PC_", seq_len(ncol(pca_mat)))
        seurat_obj[["pca"]] <- SeuratObject::CreateDimReducObject(
            embeddings = pca_mat,
            key = "PC_",
            assay = "RNA"
        )

        seurat_obj <- Seurat::FindNeighbors(
            seurat_obj,
            k.param = k_param,
            compute.SNN = FALSE,
            prune.SNN = 0,
            reduction = "pca",
            dims = dims,
            force.recalc = TRUE,
            return.neighbor = TRUE
        )

        cli::cli_alert_info("Running RareQ::FindRare() with k={.val {k}}, Q_cut={.val {Q_cut}}")
        cluster <- RareQ::FindRare(
            seurat_obj,
            k = k,
            Q_cut = Q_cut,
            ratio = ratio,
            max_iter = max_iter,
            ...
        )

        sce$rare_cluster <- as.character(cluster)

        cluster_counts <- sort(table(cluster))
        n_clusters <- length(cluster_counts)
        rare_threshold <- 5
        rare_clusters <- names(cluster_counts[cluster_counts <= rare_threshold])

        sce <- sclet_set_analysis(
            sce,
            "rare_cells",
            list(
                id = name,
                method = "RareQ",
                label_col = "rare_cluster",
                n_clusters = n_clusters,
                n_rare_clusters = length(rare_clusters),
                rare_clusters = rare_clusters,
                params = list(k = k, Q_cut = Q_cut, ratio = ratio),
                timestamp = Sys.time()
            )
        )
        sce <- sclet_set_analysis_state(
            object = sce,
            type = "rare_cells",
            id = name,
            method = "RareQ",
            inputs = list(
                reduction = reduction,
                dims = dims,
                k_param = k_param,
                k = k,
                Q_cut = Q_cut,
                ratio = ratio
            ),
            artifacts = list(
                analysis_key = "rare_cells",
                label_col = "rare_cluster"
            ),
            summary = list(
                n_clusters = n_clusters,
                n_rare_clusters = length(rare_clusters),
                n_cells = ncol(sce)
            )
        )
        sce <- sclet_log_command(
            sce,
            "RunRareCellDetection",
            params = list(
                method = "RareQ",
                reduction = reduction,
                dims = paste0(range(dims), collapse = ":"),
                k_param = k_param,
                k = k,
                Q_cut = Q_cut,
                ratio = ratio,
                name = name
            ),
            outputs = list(
                analysis = "rare_cells",
                state = "rare_cells"
            )
        )

        return(sce)
    }

    stopifnot(is(sce, "SingleCellExperiment"))
    stopifnot(is.character(name) && nzchar(name))

    message(
        "RunRareCellDetection() is currently a reserved interface. ",
        "Rare cell detection backends will be added in a future release."
    )

    placeholder <- list(
        id = name,
        method = if (is.null(method)) "reserved" else method,
        status = "reserved",
        timestamp = Sys.time()
    )

    sce <- sclet_set_analysis(sce, "rare_cells", placeholder)
    sce <- sclet_set_analysis_state(
        object = sce,
        type = "rare_cells",
        id = name,
        method = if (is.null(method)) "reserved" else method,
        inputs = list(),
        artifacts = list(analysis_key = "rare_cells"),
        summary = list(status = "reserved")
    )
    sce <- sclet_log_command(
        sce,
        "RunRareCellDetection",
        params = list(name = name),
        outputs = list(analysis = "rare_cells", state = "rare_cells")
    )

    sce
}

#' Run a state priority and perturbation sensitivity workflow
#'
#' This is the semantic shell over the state-priority mainline. It currently
#' dispatches perturbation priority analysis through Augur, and other
#' prioritization modules will be added as backends mature.
#'
#' @param sce A SingleCellExperiment object.
#' @param condition_col Character. Column in colData containing condition labels.
#' @param group_col Character. Column in colData containing cell type labels.
#' @param run_perturbation_priority Logical. Whether to run Augur priority.
#' @param run_rare_cell_detection Logical. Whether to run rare cell detection.
#' @param name Workflow analysis id. Defaults to `"state_priority"`.
#' @param ... Additional arguments passed to the selected backends.
#'
#' @return Updated `SingleCellExperiment` with state priority workflow record.
#' @export
RunStatePriorityWorkflow <- function(sce, condition_col = NULL,
    group_col = NULL, run_perturbation_priority = TRUE,
    run_rare_cell_detection = FALSE, name = "state_priority", ...) {

    component_ids <- list()

    if (isTRUE(run_perturbation_priority)) {
        priority_id <- paste0(name, "_priority")
        sce <- RunPerturbationPriority(
            sce,
            condition_col = condition_col,
            group_col = group_col,
            name = priority_id,
            ...
        )
        component_ids$perturbation_priority <- priority_id
    }

    if (isTRUE(run_rare_cell_detection)) {
        rare_id <- paste0(name, "_rare")
        sce <- RunRareCellDetection(sce, name = rare_id)
        component_ids$rare_cells <- rare_id
    }

    if (!length(component_ids)) {
        stop(
            "RunStatePriorityWorkflow() needs at least one ",
            "prioritization task enabled."
        )
    }

    workflow_analysis <- list(
        id = name,
        method = "state_priority_mainline",
        inputs = list(
            condition_col = condition_col,
            group_col = group_col
        ),
        artifacts = component_ids,
        summary = list(
            has_perturbation_priority = !is.null(component_ids$perturbation_priority),
            has_rare_cells = !is.null(component_ids$rare_cells)
        ),
        created_at = Sys.time()
    )
    sce <- sclet_set_analysis(sce, "state_priority", workflow_analysis)
    sce <- sclet_log_command(
        sce,
        "RunStatePriorityWorkflow",
        params = list(
            condition_col = condition_col,
            group_col = group_col,
            run_perturbation_priority = run_perturbation_priority,
            run_rare_cell_detection = run_rare_cell_detection,
            name = name
        ),
        outputs = list(
            workflow = name,
            perturbation_priority = component_ids$perturbation_priority,
            rare_cells = component_ids$rare_cells
        )
    )

    sce
}
