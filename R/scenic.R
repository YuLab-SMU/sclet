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
#' 
#' @return A SingleCellExperiment object with the SCENIC AUC matrix stored in `altExp(sce, "SCENIC_AUC")`
#'   and the state registered under `metadata(sce)$sclet$analyses$scenic`.
#' 
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom SummarizedExperiment assay
#' @importFrom SingleCellExperiment altExp altExp<-
#' @export
RunSCENIC <- function(sce, tfs_path, motif_annotations_path, database_paths, 
                      assay_use = "counts", num_workers = 1L, seed = 123L) {
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
    # pySCENIC expects cells as rows and genes as columns
    expr_mat <- t(as.matrix(assay(sce, assay_use)))
    gene_names <- colnames(expr_mat)
    cell_names <- rownames(expr_mat)
    
    message("Running pySCENIC via basilisk (this may take a while)...")
    
    auc_matrix <- basilisk::basiliskRun(env = sclet_scenic_env, fun = function(expr, genes, cells, tfs_f, motif_f, dbs, n_workers, s) {
        # Import Python modules
        pd <- reticulate::import("pandas")
        grnboost2 <- reticulate::import("arboreto.algo")$grnboost2
        pyscenic_prune <- reticulate::import("pyscenic.prune")
        pyscenic_aucell <- reticulate::import("pyscenic.aucell")
        
        # 1. Prepare DataFrame
        ex_df <- pd$DataFrame(expr, index = cells, columns = genes)
        
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
    
    # Register state
    metadata(sce)$sclet$analyses$scenic <- list(
        method = "pySCENIC",
        assay_use = assay_use,
        artifacts = list(
            altExp = "SCENIC_AUC"
        ),
        timestamp = Sys.time()
    )
    
    metadata(sce)$sclet$active$scenic <- "scenic"
    
    cli::cli_alert_success("pySCENIC completed. AUC matrix stored in altExp(sce, 'SCENIC_AUC').")
    return(sce)
}
