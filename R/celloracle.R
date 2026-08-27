#' Run In Silico Perturbation
#'
#' This function performs in silico gene perturbation (e.g., virtual knockout) using 
#' `celloracle` via `basilisk`. It predicts the shift in cell states upon perturbation.
#'
#' @param sce A SingleCellExperiment object.
#' @param target_gene String. The name of the gene to perturb.
#' @param perturbation_value Numeric. The expression value to set for the target gene (e.g., 0 for knockout).
#' @param base_grn_path Character. Path to the base GRN file (e.g., from motif scanning).
#' @param cluster_key String. The column in `colData(sce)` containing cluster labels.
#' @param reduction String. Reduction to use for kNN graph. Defaults to "PCA".
#' @param ... Additional arguments.
#' 
#' @return A SingleCellExperiment object with predicted cell state shifts recorded.
#' @importFrom cli cli_abort
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom SummarizedExperiment colData colData<- assay
#' @export
RunInSilicoPerturbation <- function(sce, target_gene, perturbation_value = 0.0, base_grn_path, cluster_key = ActiveIdent(sce), reduction = "PCA", ...) {
    if (!requireNamespace("basilisk", quietly = TRUE)) {
        stop("Please install 'basilisk' to run CellOracle.")
    }
    
    if (!target_gene %in% rownames(sce)) {
        stop(sprintf("Target gene '%s' not found in sce.", target_gene))
    }
    
    if (!file.exists(base_grn_path)) {
        stop(sprintf("Base GRN file '%s' not found.", base_grn_path))
    }
    
    if (!"counts" %in% SummarizedExperiment::assayNames(sce)) {
        stop("Assay 'counts' is required.")
    }
    
    counts_mat <- sclet_extract_cell_feature_matrix(sce, "counts")
    clusters <- as.character(SummarizedExperiment::colData(sce)[[cluster_key]])
    
    if (!reduction %in% SingleCellExperiment::reducedDimNames(sce)) {
        cli::cli_abort(
            "Reduction {.val {reduction}} not found in {.cls SingleCellExperiment}.",
            class = "sclet_missing_reduction"
        )
    }
    emb <- SingleCellExperiment::reducedDim(sce, reduction)
    
    message("Running CellOracle via basilisk...")
    message("(First run will take several minutes to set up the Python environment)")
    
    co_res <- sclet_basilisk_run(
        sclet_celloracle_env,
        function(counts, cluster_labels, emb, tgt, pval, grn_f, reduction_name, ...) {
        ad <- reticulate::import("anndata")
        pd <- reticulate::import("pandas")
        co <- reticulate::import("celloracle")
        sc <- reticulate::import("scanpy")
        sp <- reticulate::import("scipy.sparse")

        counts <- reticulate::r_to_py(counts)
        if (sp$issparse(counts)) {
            counts <- counts$tocsr()
        }
        
        obs <- pd$DataFrame(list(clusters = cluster_labels), index = rownames(counts))
        adata <- ad$AnnData(X = counts, obs = obs)
        adata$obsm[paste0("X_", tolower(reduction_name))] <- emb
        
        # Preprocessing inside Python
        sc$pp$normalize_total(adata, target_sum = 10000)
        sc$pp$log1p(adata)
        sc$pp$neighbors(adata, use_rep = paste0("X_", tolower(reduction_name)))
        
        # Load base GRN
        base_grn <- pd$read_csv(grn_f, index_col = 0L)
        
        # Oracle setup
        oracle <- co$Oracle()
        oracle$import_anndata_as_raw_count(adata, cluster_column_name = "clusters", ndims = as.integer(ncol(emb)))
        oracle$import_TF_data(TF_info_matrix = base_grn)
        oracle$perform_PCA()
        oracle$knn_imputation(n_pca_dims = as.integer(ncol(emb)), k = 50L, balanced = TRUE, b_sight = 3000L, b_maxl = 1500L, n_jobs = 1L)
        
        # GRN Inference
        links <- oracle$get_links(cluster_name_for_GRN_unit = "clusters", alpha = 10L, verbose_level = 0L)
        links$filter_links(p = 0.001, weight = "coef_abs", threshold_number = 2000L)
        oracle$get_cluster_specific_TFdict_from_Links(links_object = links)
        oracle$fit_GRN_for_simulation(alpha = 10L, use_cluster_specific_TFdict = TRUE)
        
        # Perturbation
        oracle$simulate_shift(perturb_condition = reticulate::dict(tgt = pval), n_propagation = 3L)
        oracle$estimate_transition_prob(n_neighbors = 200L, knn_random = TRUE, sampled_fraction = 1L)
        oracle$calculate_embedding_shift(sigma_corr = 0.05)
        
        shift <- oracle$adata$obsm["delta_X"]
        
        return(as.matrix(shift))
    },
        counts = counts_mat,
        cluster_labels = clusters,
        emb = emb,
        tgt = target_gene,
        pval = perturbation_value,
        grn_f = base_grn_path,
        reduction_name = reduction,
        ...
    )
    
    # Store the shift vector in reducedDim, e.g., "co_shift_<gene>"
    rownames(co_res) <- colnames(sce)
    colnames(co_res) <- paste0("shift_", seq_len(ncol(co_res)))

    shift_reduction <- paste0("co_shift_", target_gene)
    SingleCellExperiment::reducedDim(sce, shift_reduction) <- co_res

    # Register state via the standard analysis-state contract
    perturbation_id <- paste0("celloracle_", target_gene)
    sce <- sclet_set_analysis(
        sce,
        "celloracle",
        list(
            id = perturbation_id,
            method = "CellOracle",
            target_gene = target_gene,
            perturbation_value = perturbation_value,
            reduction = reduction,
            cluster_key = cluster_key,
            shift_reduction = shift_reduction,
            timestamp = Sys.time()
        )
    )
    sce <- sclet_set_state_record(
        object = sce,
        type = "perturbation",
        id = perturbation_id,
        active = TRUE,
        value = list(
            method = "celloracle::Oracle",
            inputs = list(
                assay = "counts",
                cluster_key = cluster_key,
                reduction = reduction,
                base_grn_path = base_grn_path
            ),
            artifacts = list(shift_reduction = shift_reduction),
            params = list(
                target_gene = target_gene,
                perturbation_value = perturbation_value
            ),
            summary = list(n_cells = ncol(sce)),
            shift_reduction = shift_reduction
        )
    )
    sce <- sclet_log_command(
        sce,
        "RunInSilicoPerturbation",
        params = list(
            target_gene = target_gene,
            perturbation_value = perturbation_value,
            id = perturbation_id
        ),
        outputs = list(
            analysis = "celloracle",
            perturbation = perturbation_id
        )
    )

    return(sce)
}
