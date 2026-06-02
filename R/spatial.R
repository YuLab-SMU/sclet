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

    c2l_res <- basilisk::basiliskRun(env = sclet_cell2location_env, fun = function(sp, ref, grp, btc, max_ep, ...) {
        c2l <- reticulate::import("cell2location")
        ad <- reticulate::import("anndata")
        pd <- reticulate::import("pandas")
        scvi <- reticulate::import("scvi")
        sp_sparse <- reticulate::import("scipy.sparse")

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
#' Reserve interface for spatial colocalization workflows.
#' Records the analysis intent as a typed state so that downstream
#' accessors and visualizations can be designed around this layer before
#' a specific backend is finalized.
#'
#' @param sce A SingleCellExperiment object with spatial data.
#' @param name Analysis record id. Defaults to `"colocalization"`.
#' @param ... Reserved for future parameters.
#'
#' @return Updated `SingleCellExperiment` with a reserved colocalization record.
#' @export
RunSpatialColocalization <- function(sce, name = "colocalization", ...) {
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

#' Run spatial niche analysis
#'
#' Reserve interface for spatial niche detection workflows.
#' Records the analysis intent so that niche-level state and visualization
#' conventions can be designed before backend integration.
#'
#' @param sce A SingleCellExperiment object with spatial data.
#' @param name Analysis record id. Defaults to `"spatial_niche"`.
#' @param ... Reserved for future parameters.
#'
#' @return Updated `SingleCellExperiment` with a reserved niche record.
#' @export
RunSpatialNiche <- function(sce, name = "spatial_niche", ...) {
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
