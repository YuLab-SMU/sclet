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
#' @param ... Additional arguments.
#' 
#' @return A SingleCellExperiment object with deconvolution proportions added to `colData`.
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom SummarizedExperiment assay colData colData<-
#' @export
RunSpatialDeconvolution <- function(sce_spatial, sce_ref, ref_group_key, ref_batch_key = NULL, max_epochs = 250L, ...) {
    if (!requireNamespace("basilisk", quietly = TRUE)) {
        stop("Please install 'basilisk' to run cell2location.")
    }
    
    if (!"counts" %in% SummarizedExperiment::assayNames(sce_spatial)) {
        stop("Assay 'counts' is required in spatial data.")
    }
    if (!"counts" %in% SummarizedExperiment::assayNames(sce_ref)) {
        stop("Assay 'counts' is required in reference data.")
    }
    
    sp_mat <- t(as.matrix(SummarizedExperiment::assay(sce_spatial, "counts")))
    ref_mat <- t(as.matrix(SummarizedExperiment::assay(sce_ref, "counts")))
    
    # Common genes
    common_genes <- intersect(colnames(sp_mat), colnames(ref_mat))
    if (length(common_genes) == 0) {
        stop("No common genes found between spatial and reference data.")
    }
    
    sp_mat <- sp_mat[, common_genes]
    ref_mat <- ref_mat[, common_genes]
    
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
        
        # Reference AnnData
        obs_ref <- pd$DataFrame(list(labels = grp), index = rownames(ref))
        if (!is.null(btc)) {
            obs_ref$batch <- btc
        }
        adata_ref <- ad$AnnData(X = ref, obs = obs_ref)
        
        # Spatial AnnData
        adata_sp <- ad$AnnData(X = sp)
        
        # Train reference model
        if (!is.null(btc)) {
            c2l$models$RegressionModel$setup_anndata(adata_ref, batch_key = "batch", labels_key = "labels")
        } else {
            c2l$models$RegressionModel$setup_anndata(adata_ref, labels_key = "labels")
        }
        
        mod_ref <- c2l$models$RegressionModel(adata_ref)
        mod_ref$train(max_epochs = as.integer(max_ep))
        
        # Export reference signatures
        inf_dict <- mod_ref$export_posterior(adata_ref, sample_kwargs=reticulate::dict(num_samples=1000L, batch_size=2500L))
        inf_df <- adata_ref$varm['means_per_cluster_mu_fg']
        
        # Train spatial model
        c2l$models$Cell2location$setup_anndata(adata_sp)
        mod_sp <- c2l$models$Cell2location(adata_sp, cell_state_df = inf_df, 
                                           N_cells_per_location = 30L, detection_alpha = 20L)
        mod_sp$train(max_epochs = as.integer(max_ep * 10))
        
        sp_post <- mod_sp$export_posterior(adata_sp, sample_kwargs=reticulate::dict(num_samples=1000L, batch_size=2500L))
        
        # Extract proportions
        w_df <- adata_sp$obsm['q05_cell_abundance_w_sf']
        
        return(as.matrix(w_df))
    }, sp = sp_mat, ref = ref_mat, grp = ref_groups, btc = ref_batches, max_ep = max_epochs, ...)
    
    # Store results in colData
    # Remove the prefix added by cell2location to column names
    cnames <- gsub("^q05cell_abundance_w_sf_", "", colnames(c2l_res))
    colnames(c2l_res) <- paste0("c2l_", cnames)
    
    for (cn in colnames(c2l_res)) {
        SummarizedExperiment::colData(sce_spatial)[[cn]] <- c2l_res[, cn]
    }
    
    # Register state
    S4Vectors::metadata(sce_spatial)$sclet$analyses$spatial_deconv <- list(
        method = "cell2location",
        ref_group_key = ref_group_key,
        columns = colnames(c2l_res),
        timestamp = Sys.time()
    )
    
    S4Vectors::metadata(sce_spatial)$sclet$active$spatial <- "spatial_deconv"
    
    return(sce_spatial)
}
