#' Run Spatial Deconvolution
#'
#' This function runs spatial transcriptomics deconvolution using `cell2location` via `basilisk`.
#' It requires a spatial SingleCellExperiment (or SpatialExperiment) and a reference scRNA-seq SCE.
#'
#' @param sce_spatial A SingleCellExperiment object containing spatial data.
#' @param sce_ref A SingleCellExperiment object containing reference scRNA-seq data.
#' @param ref_group_key String specifying the column in `colData(sce_ref)` containing cell type labels.
#' @param ref_batch_key Optional string specifying the batch column in `sce_ref`.
#' @param max_epochs Integer. Number of training epochs. Defaults to 250.
#' @param name Character. Spatial analysis record id. Defaults to `"spatial_deconv"`.
#' @param ... Additional arguments.
#'
#' @return A SingleCellExperiment object with deconvolution proportions added to `colData`.
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom SummarizedExperiment assay colData colData<-
#' @export
RunSpatialDeconvolution <- function(sce_spatial, sce_ref, ref_group_key,
    ref_batch_key = NULL, max_epochs = 250L, name = "spatial_deconv", ...) {
    if (!requireNamespace("basilisk", quietly = TRUE)) {
        stop("Please install 'basilisk' to run cell2location.")
    }

    if (!"counts" %in% SummarizedExperiment::assayNames(sce_spatial)) {
        stop("Assay 'counts' is required in spatial data.")
    }
    if (!"counts" %in% SummarizedExperiment::assayNames(sce_ref)) {
        stop("Assay 'counts' is required in reference data.")
    }

    common_genes <- intersect(rownames(sce_spatial), rownames(sce_ref))
    if (length(common_genes) == 0) {
        stop("No common genes found between spatial and reference data.")
    }
    sp_mat <- sclet_extract_cell_feature_matrix(sce_spatial, "counts", common_genes)
    ref_mat <- sclet_extract_cell_feature_matrix(sce_ref, "counts", common_genes)

    ref_groups <- as.character(SummarizedExperiment::colData(sce_ref)[[ref_group_key]])

    ref_batches <- NULL
    if (!is.null(ref_batch_key)) {
        ref_batches <- as.character(SummarizedExperiment::colData(sce_ref)[[ref_batch_key]])
    }

    message("Running Cell2location via basilisk (this is computationally intensive)...")
    message("(First run will take several minutes to set up the Python environment)")

    c2l_res <- basilisk::basiliskRun(env = sclet_cell2location_env, fun = function(sp, ref, grp, btc, max_ep, ...) {
        c2l <- reticulate::import("cell2location")
        ad <- reticulate::import("anndata")
        pd <- reticulate::import("pandas")
        scvi <- reticulate::import("scvi")
        sp_sparse <- reticulate::import("scipy.sparse")

        ref <- reticulate::r_to_py(ref)
        sp <- reticulate::r_to_py(sp)
        if (sp_sparse$issparse(ref)) {
            ref <- ref$tocsr()
        }
        if (sp_sparse$issparse(sp)) {
            sp <- sp$tocsr()
        }

        obs_ref <- pd$DataFrame(list(labels = grp), index = rownames(ref))
        if (!is.null(btc)) {
            obs_ref$batch <- btc
        }
        adata_ref <- ad$AnnData(X = ref, obs = obs_ref)

        adata_sp <- ad$AnnData(X = sp)

        if (!is.null(btc)) {
            c2l$models$RegressionModel$setup_anndata(adata_ref, batch_key = "batch", labels_key = "labels")
        } else {
            c2l$models$RegressionModel$setup_anndata(adata_ref, labels_key = "labels")
        }

        mod_ref <- c2l$models$RegressionModel(adata_ref)
        mod_ref$train(max_epochs = as.integer(max_ep))

        inf_dict <- mod_ref$export_posterior(adata_ref, sample_kwargs=reticulate::dict(num_samples=1000L, batch_size=2500L))
        inf_df <- adata_ref$varm['means_per_cluster_mu_fg']

        c2l$models$Cell2location$setup_anndata(adata_sp)
        mod_sp <- c2l$models$Cell2location(adata_sp, cell_state_df = inf_df,
                                           N_cells_per_location = 30L, detection_alpha = 20L)
        mod_sp$train(max_epochs = as.integer(max_ep * 10))

        sp_post <- mod_sp$export_posterior(adata_sp, sample_kwargs=reticulate::dict(num_samples=1000L, batch_size=2500L))

        w_df <- adata_sp$obsm['q05_cell_abundance_w_sf']

        return(as.matrix(w_df))
    }, sp = sp_mat, ref = ref_mat, grp = ref_groups, btc = ref_batches, max_ep = max_epochs, ...)

    cnames <- gsub("^q05cell_abundance_w_sf_", "", colnames(c2l_res))
    colnames(c2l_res) <- paste0("c2l_", cnames)

    for (cn in colnames(c2l_res)) {
        SummarizedExperiment::colData(sce_spatial)[[cn]] <- c2l_res[, cn]
    }

    sce_spatial <- sclet_set_analysis(
        sce_spatial,
        "spatial_deconv",
        list(
            id = name,
            method = "cell2location",
            ref_group_key = ref_group_key,
            columns = colnames(c2l_res),
            timestamp = Sys.time()
        )
    )
    sce_spatial <- sclet_set_analysis_state(
        object = sce_spatial,
        type = "spatial",
        id = name,
        method = "cell2location",
        inputs = list(
            ref_group_key = ref_group_key,
            ref_batch_key = ref_batch_key,
            n_common_genes = length(common_genes)
        ),
        artifacts = list(
            analysis_key = "spatial_deconv",
            columns = colnames(c2l_res)
        ),
        params = list(
            max_epochs = max_epochs
        )
    )
    sce_spatial <- sclet_log_command(
        sce_spatial,
        "RunSpatialDeconvolution",
        params = list(
            ref_group_key = ref_group_key,
            ref_batch_key = ref_batch_key,
            max_epochs = max_epochs,
            name = name
        ),
        outputs = list(
            analysis = "spatial_deconv",
            state = "spatial"
        )
    )

    sce_spatial
}

#' Run spatial colocalization analysis
#'
#' Detects spatial co-distribution patterns between features using SVP's
#' global bivariate spatial autocorrelation analysis. Given a set of
#' spatially-resolved feature columns (e.g., cell type proportions from
#' deconvolution), `runGLOBALBV()` computes the spatial correlation between
#' every pair, identifying which features tend to co-occur in space.
#'
#' @param sce A SingleCellExperiment or SpatialExperiment object.
#' @param features Character vector. Column names in `colData(sce)` to test
#'   for spatial colocalization. If `NULL`, uses the columns from the most
#'   recent `RunSpatialDeconvolution()` result.
#' @param method Character. Backend. Currently only `"SVP"` is supported.
#' @param name Analysis record id. Defaults to `"colocalization"`.
#' @param ... Additional arguments passed to `SVP::runGLOBALBV()`.
#'
#' @return Updated SingleCellExperiment with colocalization results in
#'   the analysis state and a pairwise correlation matrix in
#'   `metadata(sce)$sclet$analyses$colocalization$matrix`.
#' @export
RunSpatialColocalization <- function(sce, features = NULL,
    method = c("SVP"), name = "colocalization", ...) {
    method <- match.arg(method)
    stopifnot(is(sce, "SingleCellExperiment"))
    stopifnot(is.character(name) && nzchar(name))

    if (identical(method, "SVP")) {
        if (!requireNamespace("SVP", quietly = TRUE)) {
            stop(
                "Package 'SVP' is needed for spatial colocalization. ",
                "Please install it with: BiocManager::install('SVP')",
                call. = FALSE
            )
        }

        if (is.null(features)) {
            deconv <- get_spatial(sce)
            if (!is.null(deconv)) {
                cols <- deconv$artifacts$columns
                if (is.null(cols)) cols <- deconv$columns
                features <- intersect(cols, colnames(SummarizedExperiment::colData(sce)))
            }
            if (is.null(features) || length(features) < 2) {
                stop(
                    "`features` must be provided, or spatial deconvolution ",
                    "results must be available with at least 2 feature columns."
                )
            }
        }
        if (!all(features %in% colnames(SummarizedExperiment::colData(sce)))) {
            missing <- setdiff(features, colnames(SummarizedExperiment::colData(sce)))
            stop(sprintf(
                "Feature column(s) not found in colData: %s",
                paste(missing, collapse = ", ")
            ))
        }

        cli::cli_alert_info(
            "Running SVP global bivariate analysis on {.val {length(features)}} features"
        )
        sce <- SVP::runGLOBALBV(sce, features1 = features, ...)

        gbv_result <- S4Vectors::metadata(sce)$SVP$GLOBALBV
        corr_matrix <- if (!is.null(gbv_result)) gbv_result$corr else NULL

        sce <- sclet_set_analysis(sce, "colocalization", list(
            id = name,
            method = "SVP",
            features = features,
            n_features = length(features),
            matrix = corr_matrix,
            timestamp = Sys.time()
        ))
        sce <- sclet_set_analysis_state(
            object = sce,
            type = "spatial",
            id = name,
            method = "SVP_GLOBALBV",
            inputs = list(
                features = features,
                n_features = length(features)
            ),
            artifacts = list(
                analysis_key = "colocalization",
                n_features = length(features)
            ),
            summary = list(
                n_features = length(features)
            )
        )
        sce <- sclet_log_command(
            sce,
            "RunSpatialColocalization",
            params = list(
                method = "SVP",
                features = features,
                name = name
            ),
            outputs = list(
                analysis = "colocalization",
                state = "spatial"
            )
        )

        return(sce)
    }

    stopifnot(is(sce, "SingleCellExperiment"))
    stopifnot(is.character(name) && nzchar(name))

    message(
        "RunSpatialColocalization() is currently a reserved interface. ",
        "Spatial colocalization backends will be added in a future release."
    )

    sce <- sclet_set_analysis(sce, "colocalization", list(
        id = name,
        method = "reserved",
        status = "reserved",
        timestamp = Sys.time()
    ))
    sce <- sclet_set_analysis_state(
        object = sce,
        type = "spatial",
        id = name,
        method = "reserved_colocalization",
        inputs = list(),
        artifacts = list(analysis_key = "colocalization"),
        summary = list(status = "reserved")
    )
    sce <- sclet_log_command(
        sce,
        "RunSpatialColocalization",
        params = list(name = name),
        outputs = list(analysis = "colocalization", state = "spatial")
    )

    sce
}

#' Plot spatial deconvolution cell type proportions
#'
#' Summarizes the mean estimated cell type proportion across all spatial spots
#' as a ranked bar chart.
#'
#' @param object A SingleCellExperiment object with spatial deconvolution results.
#' @param n Integer. Number of top cell types to display. Defaults to 20.
#' @param fill_color Character. Bar fill color. Defaults to `"darkgreen"`.
#'
#' @return A ggplot object.
#' @export
plot_spatial_deconvolution <- function(object, n = 20, fill_color = "darkgreen") {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is needed for this function to work. Please install it.")
    }

    deconv <- get_spatial(object)
    if (is.null(deconv)) {
        stop(
            "No spatial deconvolution results found. ",
            "Run RunSpatialDeconvolution() or RunSpatialWorkflow() first."
        )
    }
    columns <- deconv$artifacts$columns
    if (is.null(columns)) {
        columns <- deconv$columns
    }
    if (is.null(columns) || length(columns) == 0) {
        stop("No deconvolution columns recorded in spatial state.")
    }
    cd <- as.data.frame(SummarizedExperiment::colData(object))
    present <- intersect(columns, colnames(cd))
    if (length(present) == 0) {
        stop("None of the recorded deconvolution columns are present in colData.")
    }

    means <- colMeans(cd[, present, drop = FALSE], na.rm = TRUE)
    labels <- gsub("^c2l_", "", present)
    df <- data.frame(cell_type = labels, proportion = means, stringsAsFactors = FALSE)
    df <- df[order(df$proportion, decreasing = TRUE), , drop = FALSE]
    df <- utils::head(df, n)
    df$cell_type <- factor(df$cell_type, levels = rev(df$cell_type))

    ggplot2::ggplot(df, ggplot2::aes(x = .data$cell_type, y = .data$proportion)) +
        ggplot2::geom_col(fill = fill_color, width = 0.7) +
        ggplot2::coord_flip() +
        ggplot2::labs(
            x = "Cell Type",
            y = "Mean Proportion",
            title = "Spatial Deconvolution - Cell Type Composition"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold", size = 13)
        )
}

#' Plot spatial cell type composition heatmap
#'
#' Displays cell type proportions across spatial spots as a heatmap.
#'
#' @param object A SingleCellExperiment object with spatial deconvolution results.
#' @param n Integer. Number of top cell types to display. Defaults to 15.
#' @param scale_rows Logical. If TRUE, z-score scales each cell type across spots.
#' @param low,high Colors for the fill gradient.
#'
#' @return A ggplot object.
#' @export
plot_spatial_composition <- function(object, n = 15, scale_rows = TRUE,
                                     low = "#2166AC", high = "#B2182B") {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is needed for this function to work. Please install it.")
    }

    deconv <- get_spatial(object)
    if (is.null(deconv)) {
        stop(
            "No spatial deconvolution results found. ",
            "Run RunSpatialDeconvolution() or RunSpatialWorkflow() first."
        )
    }
    columns <- deconv$artifacts$columns
    if (is.null(columns)) {
        columns <- deconv$columns
    }
    if (is.null(columns) || length(columns) == 0) {
        stop("No deconvolution columns recorded in spatial state.")
    }
    cd <- as.data.frame(SummarizedExperiment::colData(object))
    present <- intersect(columns, colnames(cd))
    if (length(present) == 0) {
        stop("None of the recorded deconvolution columns are present in colData.")
    }

    means <- colMeans(cd[, present, drop = FALSE], na.rm = TRUE)
    top_idx <- order(means, decreasing = TRUE)[seq_len(min(n, length(means)))]
    keep_cols <- present[top_idx]
    labels <- gsub("^c2l_", "", keep_cols)

    mat <- as.matrix(cd[, keep_cols, drop = FALSE])
    if (scale_rows) {
        mat <- t(scale(t(mat)))
    }

    df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
    colnames(df) <- c("spot", "cell_type", "proportion")
    df$cell_type <- factor(rep(labels, each = nrow(mat)), levels = labels)

    ggplot2::ggplot(df, ggplot2::aes(x = .data$spot, y = .data$cell_type, fill = .data$proportion)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradient2(low = low, mid = "white", high = high) +
        ggplot2::labs(
            x = "Spot",
            y = "Cell Type",
            fill = if (scale_rows) "Scaled\nproportion" else "Proportion",
            title = "Spatial Composition Heatmap"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            plot.title = ggplot2::element_text(face = "bold", size = 13)
        )
}

#' Run spatial niche analysis
#'
#' Run spatial niche detection
#'
#' Detects local spatial hotspots (niches) for specified features using
#' SVP's Local Indicators of Spatial Association (LISA). For each feature,
#' `runLISA()` categorizes every spot as High-High (hotspot), Low-Low,
#' High-Low, Low-High, or not significant, identifying regions where a
#' cell state is spatially enriched above random expectation.
#'
#' @param sce A SingleCellExperiment or SpatialExperiment object.
#' @param features Character vector. Column names in `colData(sce)` to run
#'   LISA analysis on. If `NULL`, uses columns from the most recent
#'   `RunSpatialDeconvolution()` result.
#' @param method Character. Backend. Currently only `"SVP"` is supported.
#' @param name Analysis record id. Defaults to `"spatial_niche"`.
#' @param ... Additional arguments passed to `SVP::runLISA()`.
#'
#' @return Updated SingleCellExperiment with LISA hotspot categories in
#'   `colData` (columns prefixed with `lisa_`) and a niche state record.
#' @export
RunSpatialNiche <- function(sce, features = NULL,
    method = c("SVP"), name = "spatial_niche", ...) {
    method <- match.arg(method)
    stopifnot(is(sce, "SingleCellExperiment"))
    stopifnot(is.character(name) && nzchar(name))

    if (identical(method, "SVP")) {
        if (!requireNamespace("SVP", quietly = TRUE)) {
            stop(
                "Package 'SVP' is needed for spatial niche detection. ",
                "Please install it with: BiocManager::install('SVP')",
                call. = FALSE
            )
        }

        if (is.null(features)) {
            deconv <- get_spatial(sce)
            if (!is.null(deconv)) {
                cols <- deconv$artifacts$columns
                if (is.null(cols)) cols <- deconv$columns
                features <- intersect(cols, colnames(SummarizedExperiment::colData(sce)))
            }
            if (is.null(features) || length(features) == 0) {
                stop(
                    "`features` must be provided, or spatial deconvolution ",
                    "results must be available with feature columns."
                )
            }
        }
        missing <- setdiff(features, colnames(SummarizedExperiment::colData(sce)))
        if (length(missing) > 0) {
            stop(sprintf(
                "Feature column(s) not found in colData: %s",
                paste(missing, collapse = ", ")
            ))
        }

        lisa_features <- list()
        for (feat in features) {
            cli::cli_alert_info("Running SVP LISA niche detection for {.val {feat}}")
            sce <- SVP::runLISA(sce, feature = feat, ...)
            lisa_features[[feat]] <- feat
        }

        sce <- sclet_set_analysis(sce, "spatial_niche", list(
            id = name,
            method = "SVP_LISA",
            features = features,
            n_features = length(features),
            timestamp = Sys.time()
        ))
        sce <- sclet_set_analysis_state(
            object = sce,
            type = "spatial",
            id = name,
            method = "SVP_LISA",
            inputs = list(
                features = features,
                n_features = length(features)
            ),
            artifacts = list(
                analysis_key = "spatial_niche",
                features = features,
                n_features = length(features)
            ),
            summary = list(
                n_features = length(features)
            )
        )
        sce <- sclet_log_command(
            sce,
            "RunSpatialNiche",
            params = list(
                method = "SVP",
                features = features,
                name = name
            ),
            outputs = list(
                analysis = "spatial_niche",
                state = "spatial"
            )
        )

        return(sce)
    }

    stopifnot(is(sce, "SingleCellExperiment"))
    stopifnot(is.character(name) && nzchar(name))

    message(
        "RunSpatialNiche() is currently a reserved interface. ",
        "Spatial niche detection backends will be added in a future release."
    )

    sce <- sclet_set_analysis(sce, "spatial_niche", list(
        id = name,
        method = "reserved",
        status = "reserved",
        timestamp = Sys.time()
    ))
    sce <- sclet_set_analysis_state(
        object = sce,
        type = "spatial",
        id = name,
        method = "reserved_niche",
        inputs = list(),
        artifacts = list(analysis_key = "spatial_niche"),
        summary = list(status = "reserved")
    )
    sce <- sclet_log_command(
        sce,
        "RunSpatialNiche",
        params = list(name = name),
        outputs = list(analysis = "spatial_niche", state = "spatial")
    )

    sce
}

#' Run a spatial context and niche analysis workflow
#'
#' This is the semantic shell over the spatial-context mainline. It currently
#' dispatches spatial deconvolution and reserves slots for colocalization and
#' niche detection as backends mature.
#'
#' @param sce A SingleCellExperiment object with spatial data.
#' @param sce_ref Optional reference SingleCellExperiment for deconvolution.
#' @param ref_group_key Column in ref colData with cell type labels.
#' @param run_deconvolution Logical. Whether to run cell2location deconvolution.
#' @param run_colocalization Logical. Whether to run reserved colocalization.
#' @param run_niche Logical. Whether to run reserved niche detection.
#' @param name Workflow analysis id. Defaults to `"spatial_workflow"`.
#' @param ... Additional arguments passed to the selected backends.
#'
#' @return Updated `SingleCellExperiment` with spatial workflow record.
#' @export
RunSpatialWorkflow <- function(sce, sce_ref = NULL, ref_group_key = NULL,
    run_deconvolution = TRUE, run_colocalization = FALSE,
    run_niche = FALSE, name = "spatial_workflow", ...) {

    component_ids <- list()

    if (isTRUE(run_deconvolution)) {
        if (is.null(sce_ref) || is.null(ref_group_key)) {
            stop(
                "`sce_ref` and `ref_group_key` are required ",
                "when `run_deconvolution = TRUE`."
            )
        }
        deconv_id <- paste0(name, "_deconv")
        sce <- RunSpatialDeconvolution(
            sce,
            sce_ref = sce_ref,
            ref_group_key = ref_group_key,
            name = deconv_id,
            ...
        )
        component_ids$deconvolution <- deconv_id
    }

    if (isTRUE(run_colocalization)) {
        coloc_id <- paste0(name, "_coloc")
        sce <- RunSpatialColocalization(sce, name = coloc_id)
        component_ids$colocalization <- coloc_id
    }

    if (isTRUE(run_niche)) {
        niche_id <- paste0(name, "_niche")
        sce <- RunSpatialNiche(sce, name = niche_id)
        component_ids$niche <- niche_id
    }

    if (!length(component_ids)) {
        stop(
            "RunSpatialWorkflow() needs at least one ",
            "spatial task enabled."
        )
    }

    workflow_analysis <- list(
        id = name,
        method = "spatial_mainline",
        inputs = list(
            ref_group_key = ref_group_key
        ),
        artifacts = component_ids,
        summary = list(
            has_deconvolution = !is.null(component_ids$deconvolution),
            has_colocalization = !is.null(component_ids$colocalization),
            has_niche = !is.null(component_ids$niche)
        ),
        created_at = Sys.time()
    )
    sce <- sclet_set_analysis(sce, "spatial_workflow", workflow_analysis)
    sce <- sclet_log_command(
        sce,
        "RunSpatialWorkflow",
        params = list(
            ref_group_key = ref_group_key,
            run_deconvolution = run_deconvolution,
            run_colocalization = run_colocalization,
            run_niche = run_niche,
            name = name
        ),
        outputs = list(
            workflow = name,
            deconvolution = component_ids$deconvolution,
            colocalization = component_ids$colocalization,
            niche = component_ids$niche
        )
    )

    sce
}
