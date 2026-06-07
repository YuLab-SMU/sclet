#' Run pySCENIC Pipeline
#'
#' This function runs the pySCENIC pipeline (GRN inference, motif enrichment, and AUCell scoring)
#' via `basilisk` to isolate the Python environment. The resulting AUC matrix is stored as an 
#' alternative experiment (`altExp`) in the SingleCellExperiment object, and the analysis state 
#' is updated to track the SCENIC results.
#'
#' @param sce A SingleCellExperiment object.
#' @param tfs_path Character. Path to the transcription factors list file (e.g., a `.txt` file with one TF per line).
#' @param motif_annotations_path Character. Path to the motif annotations file (e.g., a `.tbl` file).
#' @param database_paths Character vector. Paths to the cisTarget ranking databases (`.feather` files).
#' @param assay_use Character. The name of the assay to use for inference. Defaults to "counts".
#' @param num_workers Integer. Number of workers for parallel processing. Defaults to 1.
#' @param seed Integer. Random seed for reproducibility. Defaults to 123.
#' @param name Character. Record id for this SCENIC run. Defaults to "scenic".
#' 
#' @return A SingleCellExperiment object with the SCENIC AUC matrix stored in `altExp(sce, "SCENIC_AUC")`
#'   and the state registered in the analysis-state contract.
#' 
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom SummarizedExperiment assay
#' @importFrom SingleCellExperiment altExp altExp<-
#' @export
RunSCENIC <- function(sce, tfs_path, motif_annotations_path, database_paths, 
                      assay_use = "counts", num_workers = 1L, seed = 123L,
                      name = "scenic") {
    if (!requireNamespace("basilisk", quietly = TRUE)) {
        cli::cli_abort(
            "Please install {.pkg basilisk} to run pySCENIC.",
            class = "sclet_package_missing"
        )
    }
    
    if (!assay_use %in% assayNames(sce)) {
        cli::cli_abort(
            "Assay {.val {assay_use}} not found in the {.cls SingleCellExperiment} object.",
            class = "sclet_missing_assay"
        )
    }
    
    if (!file.exists(tfs_path)) {
        cli::cli_abort(
            "TFs file not found: {.val {tfs_path}}",
            class = "sclet_file_not_found"
        )
    }
    if (!file.exists(motif_annotations_path)) {
        cli::cli_abort(
            "Motif annotations file not found: {.val {motif_annotations_path}}",
            class = "sclet_file_not_found"
        )
    }
    for (db in database_paths) {
        if (!file.exists(db)) {
            cli::cli_abort(
                "Database file not found: {.val {db}}",
                class = "sclet_file_not_found"
            )
        }
    }
    
    message("Extracting expression matrix...")
    # pySCENIC supports sparse inputs, so keep the matrix sparse across the R/Python boundary.
    expr_mat <- sclet_extract_cell_feature_matrix(sce, assay_use)
    gene_names <- colnames(expr_mat)
    cell_names <- rownames(expr_mat)
    
    message("Running pySCENIC via basilisk (this may take a while)...")
    message("(First run will take several minutes to set up the Python environment)")
    
    auc_matrix <- basilisk::basiliskRun(env = sclet_scenic_env, fun = function(expr, genes, cells, tfs_f, motif_f, dbs, n_workers, s) {
        # Import Python modules
        pd <- reticulate::import("pandas")
        sp <- reticulate::import("scipy.sparse")
        grnboost2 <- reticulate::import("arboreto.algo")$grnboost2
        pyscenic_prune <- reticulate::import("pyscenic.prune")
        pyscenic_aucell <- reticulate::import("pyscenic.aucell")
        
        # 1. Prepare DataFrame
        expr <- reticulate::r_to_py(expr)
        if (sp$issparse(expr)) {
            expr <- expr$tocsr()
            ex_df <- pd$DataFrame$sparse$from_spmatrix(expr, index = cells, columns = genes)
        } else {
            ex_df <- pd$DataFrame(expr, index = cells, columns = genes)
        }
        
        # Load TFs
        tfs <- readLines(tfs_f)
        tfs <- intersect(tfs, genes)
        
        # 2. GRN inference
        message("Step 1: Inferring Gene Regulatory Networks (GRNBoost2)...")
        adj <- grnboost2(expression_data = ex_df, tf_names = tfs, seed = as.integer(s))
        
        # 3. Motif enrichment
        message("Step 2: Motif enrichment and regulon pruning...")
        regulons <- pyscenic_prune$prune2df(
            dbs, 
            motif_f, 
            adj, 
            num_workers = as.integer(n_workers)
        )
        
        # Convert df to regulons
        regs <- pyscenic_prune$df2regulons(regulons)
        
        # 4. AUCell scoring
        message("Step 3: AUCell scoring...")
        auc_mtx <- pyscenic_aucell$aucell(
            ex_df, 
            regs, 
            num_workers = as.integer(n_workers),
            seed = as.integer(s)
        )
        
        return(auc_mtx)
        
    }, expr = expr_mat, genes = gene_names, cells = cell_names, 
       tfs_f = tfs_path, motif_f = motif_annotations_path, dbs = database_paths, 
       n_workers = num_workers, s = seed)
    
    # auc_matrix is a pandas DataFrame returned as R data.frame (cells x regulons)
    # Convert to matrix (regulons x cells)
    auc_mat <- t(as.matrix(auc_matrix))
    
    # Store as altExp
    auc_sce <- SingleCellExperiment(assays = list(AUC = auc_mat))
    altExp(sce, "SCENIC_AUC") <- auc_sce

    scenic_analysis <- list(
        id = name,
        method = "pySCENIC",
        assay_use = assay_use,
        num_workers = num_workers,
        seed = seed,
        artifacts = list(
            altExp = "SCENIC_AUC"
        ),
        regulon_names = rownames(auc_mat),
        timestamp = Sys.time()
    )
    sce <- sclet_set_analysis(sce, "scenic", scenic_analysis)
    sce <- sclet_set_state_record(
        object = sce,
        type = "scenic",
        id = name,
        active = TRUE,
        value = list(
            method = "pySCENIC",
            inputs = list(
                assay = assay_use,
                tfs_path = tfs_path,
                motif_annotations_path = motif_annotations_path,
                database_paths = database_paths
            ),
            artifacts = list(
                analysis_key = name,
                altExp = "SCENIC_AUC",
                assay = "AUC"
            ),
            params = list(
                num_workers = num_workers,
                seed = seed
            ),
            summary = list(
                n_regulons = nrow(auc_mat),
                n_cells = ncol(auc_mat)
            ),
            regulon_names = rownames(auc_mat),
            created_at = Sys.time()
        )
    )
    
    cli::cli_alert_success("pySCENIC completed. AUC matrix stored in altExp(sce, 'SCENIC_AUC').")
    return(sce)
}

#' Plot SCENIC regulon activity on an embedding
#'
#' @title plot_scenic_activity
#' @param object A SingleCellExperiment object.
#' @param regulon Character. Regulon name to display.
#' @param reduction Character. Dimensional reduction to use. If `NULL`, uses
#'   `DefaultReduction(object)`.
#' @param id Optional scenic record id.
#' @param assay Character. Assay name inside the SCENIC altExp. Defaults to "AUC".
#' @param point_size Numeric. Point size. Defaults to 0.6.
#' @param low,high Colors for the activity gradient.
#' @return A ggplot object.
#' @export
plot_scenic_activity <- function(object, regulon, reduction = NULL, id = NULL,
                                 assay = "AUC", point_size = 0.6,
                                 low = "grey90", high = "firebrick") {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is needed for this function to work. Please install it.")
    }

    scenic_record <- get_scenic(object, id = id)
    if (is.null(scenic_record)) {
        stop("No SCENIC results found. Please run RunSCENIC() first.")
    }

    alt_exp_name <- scenic_record$artifacts$altExp
    if (is.null(alt_exp_name)) {
        alt_exp_name <- "SCENIC_AUC"
    }
    if (!alt_exp_name %in% SingleCellExperiment::altExpNames(object)) {
        stop(sprintf("SCENIC altExp '%s' not found in object.", alt_exp_name))
    }

    if (is.null(reduction)) {
        reduction <- DefaultReduction(object)
    }
    if (is.null(reduction) || !reduction %in% SingleCellExperiment::reducedDimNames(object)) {
        stop("A valid dimensional reduction is required to plot SCENIC activity.")
    }

    auc_sce <- SingleCellExperiment::altExp(object, alt_exp_name)
    auc_mat <- SummarizedExperiment::assay(auc_sce, assay)
    if (!regulon %in% rownames(auc_mat)) {
        stop(sprintf("Regulon '%s' not found in altExp '%s' assay '%s'.", regulon, alt_exp_name, assay))
    }

    emb <- SingleCellExperiment::reducedDim(object, reduction)
    df <- data.frame(
        dim1 = emb[, 1],
        dim2 = emb[, 2],
        activity = as.numeric(auc_mat[regulon, ]),
        stringsAsFactors = FALSE
    )
    df <- df[order(df$activity), , drop = FALSE]

    ggplot2::ggplot(df, ggplot2::aes(x = .data$dim1, y = .data$dim2, color = .data$activity)) +
        ggplot2::geom_point(size = point_size) +
        ggplot2::scale_color_gradient(low = low, high = high) +
        ggplot2::labs(
            x = paste0(reduction, "_1"),
            y = paste0(reduction, "_2"),
            color = "AUC",
            title = regulon
        ) +
        ggplot2::theme_classic()
}

#' Run Gene Regulatory Network analysis
#'
#' Unified semantic entry point for gene regulatory network inference. Currently
#' delegates to `RunSCENIC()` as the primary backend, keeping "GRN analysis" as
#' the user-facing concept rather than "pySCENIC" as the implementation detail.
#' Future backends (e.g. SCENIC+, decoupleR, FigR) can be routed through the
#' same interface.
#'
#' @param sce A SingleCellExperiment object.
#' @param tfs_path Character. Path to the transcription factors list file.
#' @param motif_annotations_path Character. Path to the motif annotations file.
#' @param database_paths Character vector. Paths to cisTarget ranking databases.
#' @param assay_use Character. Assay for inference. Defaults to "counts".
#' @param method Character. GRN backend. Currently only `"SCENIC"` is supported.
#' @param name Character. Record id. Defaults to `"grn"`.
#' @param ... Additional arguments passed to the backend.
#'
#' @return Updated SingleCellExperiment with GRN results.
#' @export
RunGRN <- function(sce, tfs_path, motif_annotations_path, database_paths,
                   assay_use = "counts", method = c("SCENIC"), name = "grn", ...) {
    method <- match.arg(method)
    if (identical(method, "SCENIC")) {
        sce <- RunSCENIC(
            sce,
            tfs_path = tfs_path,
            motif_annotations_path = motif_annotations_path,
            database_paths = database_paths,
            assay_use = assay_use,
            name = name,
            ...
        )
    }
    sce
}
