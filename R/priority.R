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
        # `cd` was read from `sce` above, so it does not yet see the column
        # this branch just wrote. Re-read it, or the check below refuses the
        # very column the special value exists to create: `group_col =
        # "colLabels"` died with "A valid `group_col` must be provided" on
        # every object, and no test covered the branch, which is why it
        # survived. Measured 2026-10-21.
        cd <- SummarizedExperiment::colData(sce)
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

    # Augur::calculate_auc() accepts a Seurat object, a monocle object, a sparse
    # matrix, a numeric matrix or a numeric data frame -- and nothing else. A
    # `SingleCellExperiment` assay is frequently a `DelayedMatrix` (that is what
    # TENxPBMCData and HDF5Array-backed objects hand back), which is none of
    # those, so the call dies with "input must be Seurat, monocle, sparse
    # matrix, numeric matrix, or numeric data frame" before it ever looks at the
    # data. Materialise it as a dense matrix here, once, at the boundary.
    #
    # The cost is bounded by what the caller already has in memory: the assay is
    # fully realised either way for a non-delayed backend, and for a delayed one
    # this is the first time it is touched. `sclet` does not keep the dense copy.
    expr_matrix <- sclet_malignancy_dense(expr_matrix)

    cli::cli_alert_info("Running Augur perturbation priority on {.val {group_col} by {.val {condition_col}}")
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
#' Identifies rare and major cell populations using local density analysis
#' on the kNN graph. The algorithm computes a mutual connectivity metric
#' (Q-value) from each cell's k-nearest neighbors, identifies locally dense
#' micro-clusters as seeds, propagates labels outward, and flags small
#' clusters as rare populations.
#'
#' @param sce A SingleCellExperiment object with a reduction (e.g. PCA).
#' @param method Character. Currently only `"density"` is supported.
#' @param reduction Character. Reduction to use. Defaults to "PCA".
#' @param dims Integer vector. Dimensions to use. Defaults to 1:20.
#' @param k Integer. Number of nearest neighbors for the kNN graph. Defaults to 20.
#' @param q_threshold Numeric. Mutual connectivity threshold above which a cell
#'   is considered a high-density seed. Range 0-1. Defaults to 0.25.
#' @param rare_threshold Integer. Clusters with ≤ this many cells are flagged
#'   as rare. Defaults to 10.
#' @param min_cluster_size Integer. Seed groups smaller than this are dissolved
#'   into their densest neighbouring seed group rather than reported as
#'   populations. A one-cell "population" is a density artefact, not a
#'   population, and reporting thousands of them -- one per high-connectivity
#'   cell, which is what the previous seeding rule did -- made the record
#'   uninterpretable. Defaults to 3.
#' @param name Analysis record id. Defaults to `"rareq"`.
#'
#' @return Updated SingleCellExperiment with `rare_cluster` in colData.
#' @export
RunRareCellDetection <- function(sce, method = c("density"),
    reduction = "PCA", dims = 1:20, k = 20,
    q_threshold = 0.25, rare_threshold = 10, min_cluster_size = 3L,
    name = "rareq") {
    method <- match.arg(method)
    stopifnot(is(sce, "SingleCellExperiment"))
    stopifnot(is.character(name) && nzchar(name))

    if (identical(method, "density")) {
        if (!requireNamespace("BiocNeighbors", quietly = TRUE)) {
            stop(
                "Package 'BiocNeighbors' is needed for density-based ",
                "rare cell detection. Please install it.",
                call. = FALSE
            )
        }
        if (!reduction %in% SingleCellExperiment::reducedDimNames(sce)) {
            stop(sprintf(
                "Reduction '%s' not found in reducedDimNames(sce). ",
                reduction
            ), "Please run RunPCA() first.")
        }

        emb <- SingleCellExperiment::reducedDim(sce, reduction)
        max_dim <- min(max(dims), ncol(emb))
        dims_use <- seq_len(max_dim)
        emb_sub <- emb[, dims_use, drop = FALSE]

        cli::cli_alert_info(
            "Building kNN graph (k={.val {k}}) on {.val {max_dim}} dims"
        )
        knn <- BiocNeighbors::findKNN(
            sclet_as_knn_matrix(emb_sub),
            k = k,
            BNPARAM = BiocNeighbors::KmknnParam()
        )

        n_cells <- nrow(emb_sub)
        mutual_counts <- integer(n_cells)
        for (i in seq_len(n_cells)) {
            for (j in knn$index[i, ]) {
                if (j > 0 && i %in% knn$index[j, ]) {
                    mutual_counts[i] <- mutual_counts[i] + 1L
                }
            }
        }
        q_values <- mutual_counts / k

        # Seeding rule, rewritten 2026-10-17 after the external-agent run on
        # pbmc3k produced 2419 "populations" for 2700 cells -- 75% of them
        # one cell. Two compounding defects:
        #
        # 1. `q_values >= q_threshold` is a THRESHOLD, but with k=10-20 on a
        #    compact PCA embedding the mutual-connectivity distribution sits
        #    near 1 (median 0.6, 99th pct 1.0), so thousands of cells clear
        #    any lenient threshold. Every one of them became a seed.
        # 2. Each seed was labelled with its OWN row index
        #    (`cluster[seeds] <- as.character(seeds)`), so n_seeds seeds meant
        #    n_seeds populations before any merging -- adjacent seeds were
        #    never consolidated, though by construction they are dense cells
        #    sitting next to other dense cells.
        #
        # The corrected rule: a seed is a LOCAL density maximum -- denser than
        # the median of its own neighbours -- which is the property "rare
        # population" actually needs, and adjacent seeds are merged into one
        # group via their mutual kNN links before any label propagates. The
        # old threshold survives as a FLOOR (a local maximum below q_threshold
        # is not a seed), so the parameter keeps its meaning and its default.
        n_cells_use <- length(q_values)
        # q_values[matrix] returns a FLAT vector on R >= 4.4 (matrix indexing
        # drops dimensions for atomic vectors), so rebuild the neighbours-by-
        # cell layout explicitly: row i is cell i's neighbours' Q values, and
        # the row median is the neighbourhood baseline cell i is compared
        # against.
        neighbour_q <- matrix(
            q_values[as.vector(knn$index[seq_len(n_cells_use), , drop = FALSE])],
            nrow = n_cells_use, ncol = ncol(knn$index), byrow = TRUE)
        local_max <- q_values > apply(neighbour_q, 1L, stats::median)
        seeds <- which(local_max & q_values >= q_threshold)
        n_seeds <- length(seeds)

        cli::cli_alert_info(
            "Found {.val {n_seeds}} density local maxima (Q >= {.val {q_threshold}})"
        )

        if (n_seeds == 0) {
            warning(
                "No density local maxima found with Q >= ", q_threshold,
                ". All cells assigned to cluster 1."
            )
            cluster <- rep("1", n_cells)
        } else {
            # Merge adjacent seeds: two seeds that are kNN neighbours of each
            # other belong to the same density peak, not two populations.
            if (n_seeds > 1L) {
                seed_idx <- knn$index[seeds, , drop = FALSE]
                seed_pos <- matrix(
                    match(as.vector(seed_idx), as.character(seeds)),
                    nrow = n_seeds, ncol = ncol(seed_idx))
                seed_adj <- matrix(0L, n_seeds, n_seeds)
                for (i in seq_len(n_seeds)) {
                    hit <- seed_pos[i, ]
                    hit <- hit[!is.na(hit) & hit != i]
                    seed_adj[i, hit] <- 1L
                }
                seed_groups <- igraph::components(
                    igraph::graph_from_adjacency_matrix(seed_adj,
                        mode = "max"))$membership
            } else {
                seed_groups <- setNames(1L, seeds)
            }
            group_ids <- as.integer(factor(seed_groups))
            cluster <- rep(NA_character_, n_cells)
            cluster[seeds] <- as.character(group_ids)

            for (iter in seq_len(5)) {
                for (i in seq_len(n_cells)) {
                    if (!is.na(cluster[i])) next
                    nb_labels <- cluster[knn$index[i, ]]
                    nb_labels <- nb_labels[!is.na(nb_labels)]
                    if (length(nb_labels) > 0) {
                        tbl <- sort(table(nb_labels), decreasing = TRUE)
                        cluster[i] <- names(tbl)[1]
                    }
                }
            }
            cluster[is.na(cluster)] <- "1"

            cluster <- as.integer(as.factor(cluster))
            cluster <- as.character(cluster)
        }

        cluster_counts <- sort(table(cluster))
        n_clusters <- length(cluster_counts)

        # Dissolve sub-min_cluster_size populations into their densest
        # neighbouring population. A label with one cell is where the cell
        # happened to sit, not a population; propagating thousands of them
        # into the evidence layer (each graded independently!) is how the
        # 2026-10-17 run recorded 2419 rare nodes for 2700 cells.
        if (length(min_cluster_size) == 1L && min_cluster_size > 1L) {
            small <- names(cluster_counts)[cluster_counts < min_cluster_size]
            if (length(small)) {
                for (s in small) {
                    members <- which(cluster == s)
                    if (!length(members)) next
                    nb <- unique(as.vector(knn$index[members, , drop = FALSE]))
                    nb <- nb[cluster[nb] %in% names(cluster_counts)[
                        cluster_counts >= min_cluster_size]]
                    if (!length(nb)) next
                    sink <- names(sort(table(cluster[nb]), decreasing = TRUE))[[1L]]
                    cluster[members] <- sink
                }
                cluster <- as.character(as.factor(cluster))
                cluster_counts <- sort(table(cluster))
                n_clusters <- length(cluster_counts)
                cli::cli_alert_info(
                    "Dissolved {.val {length(small)}} group{?s} below
                    min_cluster_size={.val {min_cluster_size}} into their
                    densest neighbours"
                )
            }
        }

        # Written AFTER dissolution: the local `cluster` is the final labelling,
        # and assigning it before the dissolve loop is how the first version of
        # this fix printed "Dissolved 4 groups" while the column kept the
        # dissolved labels -- the message and the artifact disagreed.
        sce$rare_cluster <- cluster

        rare_clusters <- names(cluster_counts[cluster_counts <= rare_threshold])

        sce <- sclet_set_analysis(
            sce,
            "rare_cells",
            list(
                id = name,
                method = "density",
                label_col = "rare_cluster",
                n_clusters = n_clusters,
                n_rare_clusters = length(rare_clusters),
                rare_clusters = rare_clusters,
                params = list(k = k, q_threshold = q_threshold),
                timestamp = Sys.time()
            )
        )
        sce <- sclet_set_analysis_state(
            object = sce,
            type = "rare_cells",
            id = name,
            method = "density",
            inputs = list(
                reduction = reduction,
                dims = dims,
                k = k,
                q_threshold = q_threshold
            ),
            artifacts = list(
                analysis_key = "rare_cells",
                label_col = "rare_cluster"
            ),
            summary = list(
                n_clusters = n_clusters,
                n_rare_clusters = length(rare_clusters),
                n_cells = ncol(sce),
                n_seeds = n_seeds
            )
        )
        sce <- sclet_log_command(
            sce,
            "RunRareCellDetection",
            params = list(
                method = "density",
                reduction = reduction,
                dims = paste0(range(dims_use), collapse = ":"),
                k = k,
                q_threshold = q_threshold,
                name = name
            ),
            outputs = list(
                analysis = "rare_cells",
                state = "rare_cells"
            )
        )

        return(sce)
    }
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
