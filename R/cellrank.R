#' Run CellRank for Advanced Trajectory Inference
#'
#' This function runs CellRank via `basilisk` to infer complex cellular trajectories
#' and terminal states based on RNA velocity results.
#'
#' @param sce A SingleCellExperiment object with velocity computed (via `RunVelocity`).
#' @param reduction String specifying the reduction to use for kNN graph. Defaults to "PCA".
#' @param cluster_key String specifying the column in `colData` containing cluster labels.
#' @param ... Additional arguments (currently unused).
#' 
#' @return A SingleCellExperiment object with CellRank states and transition probabilities recorded.
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom SummarizedExperiment colData colData<-
#' @export
RunCellRank <- function(sce, reduction = "PCA", cluster_key = ActiveIdent(sce), ...) {
    if (!requireNamespace("basilisk", quietly = TRUE)) {
        cli::cli_abort(
            "Please install {.pkg basilisk} to run CellRank.",
            class = "sclet_package_missing"
        )
    }

    vel_state <- S4Vectors::metadata(sce)$sclet$analyses$velocity
    if (is.null(vel_state) || is.null(vel_state$results)) {
        cli::cli_abort(
            c("No velocity results found.",
              i = "Please run {.fun RunVelocity} first."),
            class = "sclet_missing_velocity"
        )
    }

    vel_res <- vel_state$results

    if (!reduction %in% SingleCellExperiment::reducedDimNames(sce)) {
        cli::cli_abort(
            "Reduction {.val {reduction}} not found in {.cls SingleCellExperiment}.",
            class = "sclet_missing_reduction"
        )
    }

    emb <- SingleCellExperiment::reducedDim(sce, reduction)

    if (!"velocity" %in% SummarizedExperiment::assayNames(vel_res)) {
        cli::cli_abort(
            "Assay {.val velocity} not found in velocity results.",
            class = "sclet_missing_velocity_assay"
        )
    }

    v_mat <- t(as.matrix(SummarizedExperiment::assay(vel_res, "velocity")))
    s_mat <- t(as.matrix(SummarizedExperiment::assay(vel_res, "spliced")))
    u_mat <- t(as.matrix(SummarizedExperiment::assay(vel_res, "unspliced")))

    clusters <- as.character(SummarizedExperiment::colData(sce)[[cluster_key]])

    cli::cli_alert_info("Running CellRank via basilisk...")

    cr_res <- basilisk::basiliskRun(
        env = sclet_cellrank_env,
        fun = function(v, s, u, emb, cluster_labels, reduction_name) {
        ad <- reticulate::import("anndata")
        cr <- reticulate::import("cellrank")
        pd <- reticulate::import("pandas")

        obs <- pd$DataFrame(
            list(clusters = pd$Categorical(cluster_labels)),
            index = rownames(v)
        )
        adata <- ad$AnnData(X = s, obs = obs)
        adata$layers["spliced"] <- s
        adata$layers["unspliced"] <- u
        adata$layers["velocity"] <- v
        adata$obsm[paste0("X_", tolower(reduction_name))] <- emb

        sc <- reticulate::import("scanpy")
        scv <- reticulate::import("scvelo")
        sc$pp$neighbors(adata, use_rep = paste0("X_", tolower(reduction_name)))
        scv$pp$moments(adata, n_pcs = NULL, n_neighbors = NULL)
        scv$tl$velocity_graph(adata)

        vk <- cr$kernels$VelocityKernel(adata)
        vk$compute_transition_matrix()

        g <- cr$estimators$GPCCA(vk)
        g$compute_eigendecomposition()
        g$compute_schur(method = "krylov")
        g$compute_macrostates(cluster_key = "clusters")
        g$compute_terminal_states()
        g$compute_absorption_probabilities()

        term_states <- g$terminal_states
        abs_probs <- g$absorption_probabilities

        abs_mat <- tryCatch(
            as.matrix(abs_probs$X),
            error = function(e) {
                tryCatch(
                    as.matrix(abs_probs),
                    error = function(e2) {
                        cli::cli_abort(
                            c("Failed to extract absorption probabilities.",
                              i = "CellRank returned: {.val {class(abs_probs)}}"),
                            parent = e2
                        )
                    }
                )
            }
        )

        return(list(
            terminal_states = as.character(term_states),
            absorption_probs = abs_mat
        ))
    },
    v = v_mat,
    s = s_mat,
    u = u_mat,
    emb = emb,
    cluster_labels = clusters,
    reduction_name = reduction)
    
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
