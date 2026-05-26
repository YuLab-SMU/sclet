#' Run CellRank for Advanced Trajectory Inference
#'
#' This function runs CellRank via `basilisk` to infer complex cellular trajectories
#' and terminal states based on RNA velocity results.
#'
#' @param sce A SingleCellExperiment object with velocity computed (via `RunVelocity`).
#' @param reduction String specifying the reduction to use for kNN graph. Defaults to "PCA".
#' @param cluster_key String specifying the column in `colData` containing cluster labels.
#' @param ... Additional arguments.
#' 
#' @return A SingleCellExperiment object with CellRank states and transition probabilities recorded.
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom SummarizedExperiment colData colData<-
#' @export
RunCellRank <- function(sce, reduction = "PCA", cluster_key = ActiveIdent(sce), ...) {
    if (!requireNamespace("basilisk", quietly = TRUE)) {
        stop("Please install 'basilisk' to run CellRank.")
    }
    
    vel_state <- S4Vectors::metadata(sce)$sclet$analyses$velocity
    if (is.null(vel_state) || is.null(vel_state$results)) {
        stop("No velocity results found. Please run RunVelocity() first.")
    }
    
    vel_res <- vel_state$results
    
    # We need to construct an AnnData for CellRank.
    # We pass the velocity data to Python.
    # For simplicity, we use zellkonverter to convert vel_res to AnnData, 
    # but zellkonverter itself uses reticulate/basilisk.
    # Alternatively, we can construct AnnData inside basiliskRun.
    
    if (!reduction %in% SingleCellExperiment::reducedDimNames(sce)) {
        stop(sprintf("Reduction %s not found.", reduction))
    }
    
    emb <- SingleCellExperiment::reducedDim(sce, reduction)
    
    # Check if velocity matrix exists
    if (!"velocity" %in% SummarizedExperiment::assayNames(vel_res)) {
        stop("Assay 'velocity' not found in velocity results.")
    }
    
    v_mat <- t(as.matrix(SummarizedExperiment::assay(vel_res, "velocity")))
    s_mat <- t(as.matrix(SummarizedExperiment::assay(vel_res, "spliced")))
    u_mat <- t(as.matrix(SummarizedExperiment::assay(vel_res, "unspliced")))
    
    clusters <- as.character(SummarizedExperiment::colData(sce)[[cluster_key]])
    
    message("Running CellRank via basilisk...")
    
    cr_res <- basilisk::basiliskRun(
        env = sclet_cellrank_env,
        fun = function(v, s, u, emb, cluster_labels, reduction_name, ...) {
        ad <- reticulate::import("anndata")
        cr <- reticulate::import("cellrank")
        pd <- reticulate::import("pandas")
        
        obs <- pd$DataFrame(list(clusters = cluster_labels), index = rownames(v))
        adata <- ad$AnnData(X = s, obs = obs)
        adata$layers["spliced"] <- s
        adata$layers["unspliced"] <- u
        adata$layers["velocity"] <- v
        adata$obsm[paste0("X_", tolower(reduction_name))] <- emb
        
        # Compute the kNN graph with scanpy to avoid the deprecated/removed
        # write_knn_indices path in scvelo's neighbors wrapper.
        sc <- reticulate::import("scanpy")
        scv <- reticulate::import("scvelo")
        sc$pp$neighbors(adata, use_rep = paste0("X_", tolower(reduction_name)))
        scv$pp$moments(adata, n_pcs = NULL, n_neighbors = NULL)
        scv$tl$velocity_graph(adata)
        
        # CellRank pipeline
        vk <- cr$kernels$VelocityKernel(adata)
        vk$compute_transition_matrix()
        
        # Estimator
        g <- cr$estimators$GPCCA(vk)
        g$compute_eigendecomposition()
        g$compute_schur(method = "krylov")
        g$compute_macrostates(cluster_key = "clusters")
        g$compute_terminal_states()
        g$compute_absorption_probabilities()
        
        term_states <- g$terminal_states
        abs_probs <- g$absorption_probabilities
        
        return(list(
            terminal_states = as.character(term_states),
            absorption_probs = as.matrix(abs_probs$X)
        ))
    },
    v = v_mat,
    s = s_mat,
    u = u_mat,
    emb = emb,
    cluster_labels = clusters,
    reduction_name = reduction,
    ...)
    
    # Store results
    sce$cellrank_terminal_states <- cr_res$terminal_states
    
    # Store absorption probs as a matrix in metadata
    colnames(cr_res$absorption_probs) <- paste0("prob_", seq_len(ncol(cr_res$absorption_probs)))
    
    S4Vectors::metadata(sce)$sclet$analyses$cellrank <- list(
        method = "CellRank",
        cluster_key = cluster_key,
        terminal_states = cr_res$terminal_states,
        absorption_probs = cr_res$absorption_probs,
        timestamp = Sys.time()
    )
    
    S4Vectors::metadata(sce)$sclet$active$trajectory <- "cellrank"
    
    return(sce)
}
