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
#' @param output_dir Character. Directory for intermediate checkpoints. Required
#'   when `save_intermediate` or `resume` is `TRUE`. Defaults to `NULL`.
#' @param save_intermediate Logical. Whether to persist intermediate steps to
#'   `output_dir` as parquet checkpoints (GRNBoost2 adjacency and the pruned
#'   motif data frame) so a failed downstream step does not force a full rerun.
#'   Defaults to `FALSE`.
#' @param resume Logical. Whether to reuse checkpoints already present in
#'   `output_dir`, skipping any completed step (GRNBoost2 and/or motif pruning).
#'   AUCell is always recomputed from the SCE. Defaults to `FALSE`.
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
                      name = "scenic", output_dir = NULL,
                      save_intermediate = FALSE, resume = FALSE) {
    if (!requireNamespace("basilisk", quietly = TRUE)) {
        cli::cli_abort(
            "Please install {.pkg basilisk} to run pySCENIC.",
            class = "sclet_package_missing"
        )
    }
    
    if (!assay_use %in% SummarizedExperiment::assayNames(sce)) {
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

    if ((isTRUE(save_intermediate) || isTRUE(resume)) &&
        (is.null(output_dir) || !nzchar(output_dir))) {
        cli::cli_abort(
            c(
                "{.arg output_dir} must be supplied when {.arg save_intermediate} or {.arg resume} is {.val TRUE}.",
                i = "Set {.arg output_dir} to a writable directory for SCENIC checkpoints."
            ),
            class = "sclet_missing_output_dir"
        )
    }
    if (!is.null(output_dir) && nzchar(output_dir) && !dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    message("Extracting expression matrix...")
    # pySCENIC supports sparse inputs, so keep the matrix sparse across the R/Python boundary.
    expr_mat <- sclet_extract_cell_feature_matrix(sce, assay_use)
    gene_names <- colnames(expr_mat)
    cell_names <- rownames(expr_mat)
    
    message("Running pySCENIC via basilisk (this may take a while)...")
    message("(First run will take several minutes to set up the Python environment)")
    
    auc_matrix <- sclet_basilisk_run(
        sclet_scenic_env,
        function(expr, genes, cells, tfs_f, motif_f, db_paths, n_workers, s,
                 out_dir, save_mid, resume) {
        # Import Python modules. `convert = FALSE` keeps pandas/DataFrame objects
        # as Python objects so we can pass lists (modules, ranking DBs) through
        # the pySCENIC API without R auto-conversion mangling them.
        builtins <- reticulate::import_builtins(convert = FALSE)
        pd <- reticulate::import("pandas", convert = FALSE)
        sp <- reticulate::import("scipy.sparse", convert = FALSE)
        arboreto <- reticulate::import("arboreto.algo", convert = FALSE)
        scenic_utils <- reticulate::import("pyscenic.utils", convert = FALSE)
        scenic_prune <- reticulate::import("pyscenic.prune", convert = FALSE)
        scenic_aucell <- reticulate::import("pyscenic.aucell", convert = FALSE)
        rnkdb <- reticulate::import("ctxcore.rnkdb", convert = FALSE)
        os <- reticulate::import("os", convert = FALSE)

        # Build a thread-based dask client when num_workers > 1 so GRNBoost2 gets
        # a real parallel backend without launching the Nanny subprocess that
        # newer dask/distributed turns into `Nanny failed to start`. Fall back to
        # arboreto's default `'local'` scheduler when dask.distributed is absent.
        make_client <- function(nw) {
            tryCatch({
                ddist <- reticulate::import("dask.distributed", convert = FALSE)
                cluster <- ddist$LocalCluster(
                    n_workers = 1L,
                    threads_per_worker = as.integer(nw),
                    processes = FALSE,
                    dashboard_address = NULL
                )
                ddist$Client(cluster)
            }, error = function(e) NULL)
        }

        # 1. Prepare a cells x genes DataFrame (sparse-aware)
        expr <- reticulate::r_to_py(expr, convert = FALSE)
        genes_py <- reticulate::r_to_py(genes, convert = FALSE)
        cells_py <- reticulate::r_to_py(cells, convert = FALSE)
        if (reticulate::py_to_r(sp$issparse(expr))) {
            expr <- expr$tocsr()
            ex_df <- pd$DataFrame$sparse$from_spmatrix(expr, index = cells_py, columns = genes_py)
        } else {
            ex_df <- pd$DataFrame(expr, index = cells_py, columns = genes_py)
        }

        # Load TFs and keep only those present in the expression matrix
        tfs <- readLines(tfs_f)
        tfs <- intersect(tfs, genes)
        if (!length(tfs)) {
            stop("No TFs overlap with the genes in the expression matrix.")
        }

        # Intermediate checkpoint paths (only valid when output_dir is supplied)
        adj_path <- NULL
        motif_path <- NULL
        if (!is.null(out_dir) && nzchar(out_dir)) {
            os$makedirs(out_dir, exist_ok = TRUE)
            adj_path <- file.path(out_dir, "grn_adjacency.parquet")
            motif_path <- file.path(out_dir, "motifs.parquet")
        }

        # 2. GRNBoost2 inference (resumable)
        adj <- NULL
        if (isTRUE(resume) && !is.null(adj_path) && file.exists(adj_path)) {
            adj <- pd$read_parquet(adj_path)
            message("Resuming SCENIC: reused cached GRNBoost2 adjacency from ", adj_path)
        }
        if (is.null(adj)) {
            message("Step 1: Inferring Gene Regulatory Networks (GRNBoost2)...")
            client <- make_client(n_workers)
            if (is.null(client)) {
                adj <- arboreto$grnboost2(
                    expression_data = ex_df,
                    tf_names = reticulate::r_to_py(tfs, convert = FALSE),
                    seed = as.integer(s)
                )
            } else {
                adj <- arboreto$grnboost2(
                    expression_data = ex_df,
                    tf_names = reticulate::r_to_py(tfs, convert = FALSE),
                    seed = as.integer(s),
                    client_or_address = client
                )
            }
            if (isTRUE(save_mid) && !is.null(adj_path)) {
                adj$to_parquet(adj_path)
            }
        }

        # 3. Build co-expression modules from the adjacency, then (from scratch or
        #    from a checkpoint) prune to regulons via the ranking databases.
        message("Step 2: Motif enrichment and regulon pruning...")
        modules <- builtins$list(
            scenic_utils$modules_from_adjacencies(adj, ex_df)
        )

        motif_df <- NULL
        if (isTRUE(resume) && !is.null(motif_path) && file.exists(motif_path)) {
            motif_df <- pd$read_parquet(motif_path)
            message("Resuming SCENIC: reused cached motif enrichment from ", motif_path)
        }
        if (is.null(motif_df)) {
            # Wrap each feather cisTarget ranking database so prune2df() gets the
            # ranking-database objects it expects (not raw file paths).
            py_dbs <- builtins$list()
            for (f in db_paths) {
                db_name <- tools::file_path_sans_ext(basename(f))
                py_dbs$append(rnkdb$FeatherRankingDatabase(fname = f, name = db_name))
            }
            motif_df <- scenic_prune$prune2df(
                py_dbs,
                modules,
                motif_f,
                num_workers = as.integer(n_workers),
                client_or_address = "custom_multiprocessing"
            )
            if (isTRUE(save_mid) && !is.null(motif_path)) {
                motif_df$to_parquet(motif_path)
            }
        }

        # 4. Convert the pruned motif data frame into regulons.
        regs <- scenic_prune$df2regulons(motif_df)

        # 5. AUCell scoring
        message("Step 3: AUCell scoring...")
        auc_mtx <- scenic_aucell$aucell(
            ex_df,
            regs,
            num_workers = as.integer(n_workers),
            seed = as.integer(s)
        )

        reticulate::py_to_r(auc_mtx)
    },
        expr = expr_mat,
        genes = gene_names,
        cells = cell_names,
        tfs_f = tfs_path,
        motif_f = motif_annotations_path,
        db_paths = database_paths,
        n_workers = num_workers,
        s = seed,
        out_dir = output_dir,
        save_mid = save_intermediate,
        resume = resume
    )
    
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
            altExp = "SCENIC_AUC",
            output_dir = if (is.null(output_dir)) NULL else output_dir
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
                assay = "AUC",
                output_dir = if (is.null(output_dir)) NULL else output_dir
            ),
            params = list(
                num_workers = num_workers,
                seed = seed,
                save_intermediate = isTRUE(save_intermediate),
                resume = isTRUE(resume)
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
