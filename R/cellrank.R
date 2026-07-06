#' Run CellRank for Advanced Trajectory Inference
#'
#' This function runs CellRank via `basilisk` to infer complex cellular trajectories
#' and terminal states based on RNA velocity results.
#'
#' @param sce A SingleCellExperiment object with velocity computed via
#'   `RunVelocity()` or `RunRegVelo()`.
#' @param reduction String specifying the reduction to use for kNN graph. Defaults to "PCA".
#' @param cluster_key String specifying the column in `colData` containing cluster labels.
#' @param velocity_id Optional velocity record id. If `NULL`, uses the active
#'   velocity result.
#' @param backend Python execution backend. `"basilisk"` creates the packaged
#'   CellRank environment; `"reticulate"` uses the caller-configured Python.
#' @param python Optional Python executable or conda environment name used when
#'   `backend = "reticulate"`.
#' @param name trajectory record id. Defaults to `"cellrank"`.
#' @param ... Additional arguments (currently unused).
#' 
#' @return A SingleCellExperiment object with CellRank states and transition probabilities recorded.
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom utils write.csv
#' @export
RunCellRank <- function(sce, reduction = "PCA", cluster_key = ActiveIdent(sce),
                        velocity_id = NULL, backend = c("basilisk", "reticulate"),
                        python = NULL, name = "cellrank", ...) {
    backend <- match.arg(backend)
    if (identical(backend, "basilisk") && !requireNamespace("basilisk", quietly = TRUE)) {
        cli::cli_abort(
            "Please install {.pkg basilisk} to run CellRank with {.val backend = 'basilisk'}.",
            class = "sclet_package_missing"
        )
    }
    if (identical(backend, "reticulate") && !requireNamespace("reticulate", quietly = TRUE)) {
        cli::cli_abort(
            "Please install {.pkg reticulate} to run CellRank with {.val backend = 'reticulate'}.",
            class = "sclet_package_missing"
        )
    }

    vel_state <- get_velocity(sce, id = velocity_id)
    if (is.null(vel_state)) {
        cli::cli_abort(
            c("No velocity results found.",
              i = "Please run {.fun RunVelocity} or {.fun RunRegVelo} first."),
            class = "sclet_missing_velocity"
        )
    }

    if (!reduction %in% SingleCellExperiment::reducedDimNames(sce)) {
        cli::cli_abort(
            "Reduction {.val {reduction}} not found in {.cls SingleCellExperiment}.",
            class = "sclet_missing_reduction"
        )
    }
    if (is.null(cluster_key) || !nzchar(cluster_key) ||
        !cluster_key %in% colnames(SummarizedExperiment::colData(sce))) {
        cli::cli_abort(
            c("A valid {.arg cluster_key} is required for {.fun RunCellRank}.",
              i = "It must refer to an existing column in {.fun colData}."),
            class = "sclet_missing_cluster_key"
        )
    }

    emb <- as.matrix(SingleCellExperiment::reducedDim(sce, reduction))

    cellrank_input <- sclet_prepare_cellrank_velocity_input(
        sce = sce,
        vel_state = vel_state,
        velocity_id = velocity_id
    )
    sce <- cellrank_input$sce
    emb <- as.matrix(SingleCellExperiment::reducedDim(sce, reduction))
    v_mat <- cellrank_input$velocity
    s_mat <- cellrank_input$spliced
    u_mat <- cellrank_input$unspliced
    velocity_source <- cellrank_input$velocity_source
    velocity_assay <- cellrank_input$velocity_assay
    cell_names <- cellrank_input$cell_names
    gene_names <- cellrank_input$gene_names

    # Write sparse matrices to .mtx files to bypass all reticulate conversion issues
    tmp_dir <- tempfile("cellrank_")
    dir.create(tmp_dir)
    on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
    Matrix::writeMM(v_mat, file.path(tmp_dir, "velocity.mtx"))
    Matrix::writeMM(s_mat, file.path(tmp_dir, "spliced.mtx"))
    Matrix::writeMM(u_mat, file.path(tmp_dir, "unspliced.mtx"))
    # Write all auxiliary data to disk so basilisk fun receives only one string
    clusters <- as.character(SummarizedExperiment::colData(sce)[[cluster_key]])
    writeLines(clusters, file.path(tmp_dir, "clusters.txt"))
    writeLines(cell_names, file.path(tmp_dir, "cell_names.txt"))
    writeLines(gene_names, file.path(tmp_dir, "gene_names.txt"))
    writeLines(reduction, file.path(tmp_dir, "reduction.txt"))
    write.csv(emb, file.path(tmp_dir, "emb.csv"), row.names = FALSE)

    cli::cli_alert_info("Running CellRank via {.val {backend}}...")
    if (identical(backend, "basilisk")) {
        cli::cli_alert_info("(First run will take several minutes to set up the Python environment)")
    }
    prev_state <- sclet_get_state(sce)

    cr_res <- sclet_run_cellrank_backend(
        backend = backend,
        python = python,
        mtx_dir = tmp_dir
    )
    sce <- sclet_restore_state(sce, prev_state)
    
    # Store results
    terminal_state_col <- "cellrank_terminal_states"
    sce[[terminal_state_col]] <- cr_res$terminal_states

    fate_prob_names <- cr_res$absorption_prob_names
    if (is.null(fate_prob_names) || length(fate_prob_names) != ncol(cr_res$absorption_probs)) {
        fate_prob_names <- as.character(seq_len(ncol(cr_res$absorption_probs)))
    }
    fate_prob_names[is.na(fate_prob_names) | !nzchar(fate_prob_names)] <- as.character(which(is.na(fate_prob_names) | !nzchar(fate_prob_names)))
    fate_prob_cols <- paste0("cellrank_fate_", make.names(fate_prob_names, unique = TRUE))
    colnames(cr_res$absorption_probs) <- fate_prob_cols
    fate_prob_df <- as.data.frame(cr_res$absorption_probs)
    SummarizedExperiment::colData(sce)[, fate_prob_cols] <- fate_prob_df

    lineage_drivers <- sclet_normalize_cellrank_driver_table(cr_res$lineage_drivers)
    recorded_velocity_id <- velocity_id
    if (is.null(recorded_velocity_id)) {
        recorded_velocity_id <- sclet_get_active_state(sce, "velocity")
    }

    cellrank_analysis <- list(
        id = name,
        method = "CellRank",
        reduction = reduction,
        cluster_key = cluster_key,
        terminal_states = cr_res$terminal_states,
        absorption_probs = cr_res$absorption_probs,
        terminal_state_col = terminal_state_col,
        fate_probability_cols = fate_prob_cols,
        fate_probability_names = fate_prob_names,
        lineage_drivers = lineage_drivers,
        timestamp = Sys.time()
    )
    sce <- sclet_set_analysis(sce, "cellrank", cellrank_analysis)
    sce <- sclet_set_state_record(
        object = sce,
        type = "trajectory",
        id = name,
        active = TRUE,
        value = list(
            method = "CellRank",
            inputs = list(
                reduction = reduction,
                cluster_key = cluster_key,
                velocity_id = recorded_velocity_id,
                velocity_source = velocity_source,
                velocity_assay = velocity_assay,
                backend = backend,
                python = python
            ),
            artifacts = list(
                analysis_key = name,
                terminal_state_col = terminal_state_col,
                fate_probability_cols = fate_prob_cols,
                fate_probability_names = fate_prob_names,
                has_lineage_drivers = !is.null(lineage_drivers)
            ),
            params = list(),
            summary = list(
                n_terminal_states = length(unique(stats::na.omit(cr_res$terminal_states))),
                n_fate_dimensions = ncol(cr_res$absorption_probs),
                n_lineage_drivers = if (is.null(lineage_drivers)) 0L else nrow(lineage_drivers)
            ),
            reduction = reduction,
            cluster_key = cluster_key,
            terminal_states = cr_res$terminal_states,
            absorption_probs = cr_res$absorption_probs,
            terminal_state_col = terminal_state_col,
            fate_probability_cols = fate_prob_cols,
            fate_probability_names = fate_prob_names,
            lineage_drivers = lineage_drivers,
            created_at = Sys.time()
        )
    )
    sce <- sclet_log_command(
        sce,
        "RunCellRank",
        params = list(
            reduction = reduction,
            cluster_key = cluster_key,
            velocity_id = velocity_id,
            backend = backend,
            python = python,
            name = name
        ),
        outputs = list(
            analysis = "cellrank",
            trajectory = name
        )
    )

    return(sce)
}

sclet_prepare_cellrank_velocity_input <- function(sce, vel_state, velocity_id = NULL) {
    velocity_assay <- tryCatch(
        sclet_velocity_assay_name(sce, id = velocity_id, assay = NULL),
        error = function(e) NULL
    )

    if (!is.null(velocity_assay)) {
        required <- c("spliced", "unspliced", velocity_assay)
        missing <- setdiff(required, SummarizedExperiment::assayNames(sce))
        if (length(missing)) {
            cli::cli_abort(
                "Required assay(s) not found for CellRank input: {.val {missing}}.",
                class = "sclet_missing_assay"
            )
        }

        return(sclet_build_cellrank_velocity_input(
            sce = sce,
            source_object = sce,
            velocity_assay = velocity_assay,
            velocity_source = "object"
        ))
    }

    if (!is.null(vel_state$results) && methods::is(vel_state$results, "SingleCellExperiment")) {
        vel_res <- vel_state$results
        common_cells <- intersect(colnames(sce), colnames(vel_res))
        if (length(common_cells) == 0) {
            cli::cli_abort(
                "No common cells between the query object and the velocity results.",
                class = "sclet_cell_mismatch"
            )
        }

        sce <- sce[, common_cells, drop = FALSE]
        vel_res <- vel_res[, common_cells, drop = FALSE]
        if (!"velocity" %in% SummarizedExperiment::assayNames(vel_res)) {
            cli::cli_abort(
                "Assay {.val velocity} not found in velocity results.",
                class = "sclet_missing_velocity_assay"
            )
        }

        return(sclet_build_cellrank_velocity_input(
            sce = sce,
            source_object = vel_res,
            velocity_assay = "velocity",
            velocity_source = "results"
        ))
    }

    cli::cli_abort(
        "No velocity assay was recorded and no compatible velocity result object is available.",
        class = "sclet_missing_velocity_assay"
    )
}

sclet_build_cellrank_velocity_input <- function(sce, source_object, velocity_assay, velocity_source) {
    required <- c("spliced", "unspliced", velocity_assay)
    missing <- setdiff(required, SummarizedExperiment::assayNames(source_object))
    if (length(missing)) {
        cli::cli_abort(
            "Required assay(s) not found for CellRank input: {.val {missing}}.",
            class = "sclet_missing_assay"
        )
    }

    v_mat <- sclet_extract_cell_feature_matrix(source_object, velocity_assay)
    s_mat <- sclet_extract_cell_feature_matrix(source_object, "spliced")
    u_mat <- sclet_extract_cell_feature_matrix(source_object, "unspliced")
    dims <- list(
        velocity = dim(v_mat),
        spliced = dim(s_mat),
        unspliced = dim(u_mat)
    )
    if (!identical(dims$velocity, dims$spliced) || !identical(dims$velocity, dims$unspliced)) {
        cli::cli_abort(
            c("Velocity, spliced, and unspliced matrices must have matching dimensions for CellRank.",
              i = "Observed dimensions are velocity={paste(dims$velocity, collapse = ' x ')}, spliced={paste(dims$spliced, collapse = ' x ')}, unspliced={paste(dims$unspliced, collapse = ' x ')}."),
            class = "sclet_cellrank_matrix_mismatch"
        )
    }

    list(
        sce = sce,
        velocity = v_mat,
        spliced = s_mat,
        unspliced = u_mat,
        velocity_source = velocity_source,
        velocity_assay = velocity_assay,
        cell_names = colnames(source_object),
        gene_names = rownames(source_object)
    )
}

sclet_run_cellrank_backend <- function(backend, python, mtx_dir) {
    if (identical(backend, "basilisk")) {
        proc <- basilisk::basiliskStart(sclet_cellrank_env)
        on.exit(basilisk::basiliskStop(proc), add = TRUE)
        return(basilisk::basiliskRun(
            proc,
            sclet_run_cellrank_python,
            mtx_dir = mtx_dir
        ))
    }

    if (!is.null(python) && !reticulate::py_available(initialize = FALSE)) {
        if (file.exists(python)) {
            reticulate::use_python(python, required = TRUE)
        } else {
            reticulate::use_condaenv(python, required = TRUE)
        }
    }
    sclet_run_cellrank_python(mtx_dir = mtx_dir)
}

sclet_run_cellrank_python <- function(mtx_dir) {
    ad <- reticulate::import("anndata", convert = FALSE)
    cr <- reticulate::import("cellrank", convert = FALSE)
    pd <- reticulate::import("pandas", convert = FALSE)
    scio <- reticulate::import("scipy.io", convert = FALSE)

    py_attr <- function(x, name, default = NULL) {
        if (reticulate::py_has_attr(x, name)) {
            reticulate::py_get_attr(x, name)
        } else {
            default
        }
    }
    py_call <- function(x, name, ...) {
        reticulate::py_call(py_attr(x, name), ...)
    }
    py_to_r <- function(x) {
        reticulate::py_to_r(x)
    }

    as_cell_feature <- function(path, n_cells, n_genes) {
        mat <- scio$mmread(path)
        shape <- as.integer(py_to_r(py_attr(mat, "shape")))
        if (identical(shape, c(n_genes, n_cells))) {
            mat <- py_attr(mat, "T")
        } else if (!identical(shape, c(n_cells, n_genes))) {
            stop(
                "Unexpected CellRank matrix dimensions for ", basename(path),
                ": got ", paste(shape, collapse = " x "),
                ", expected ", n_cells, " x ", n_genes,
                " or ", n_genes, " x ", n_cells, "."
            )
        }
        py_call(mat, "tocsr")
    }

    cell_names <- readLines(file.path(mtx_dir, "cell_names.txt"))
    gene_names <- readLines(file.path(mtx_dir, "gene_names.txt"))
    n_cells <- length(cell_names)
    n_genes <- length(gene_names)

    v <- as_cell_feature(file.path(mtx_dir, "velocity.mtx"), n_cells, n_genes)
    s <- as_cell_feature(file.path(mtx_dir, "spliced.mtx"), n_cells, n_genes)
    u <- as_cell_feature(file.path(mtx_dir, "unspliced.mtx"), n_cells, n_genes)
    clusters <- readLines(file.path(mtx_dir, "clusters.txt"))
    emb <- as.matrix(utils::read.csv(file.path(mtx_dir, "emb.csv")))
    reduction_name <- readLines(file.path(mtx_dir, "reduction.txt"))

    obs <- pd$DataFrame(index = reticulate::r_to_py(cell_names))
    reticulate::py_set_item(
        obs,
        "clusters",
        pd$Series(
            reticulate::r_to_py(clusters),
            dtype = "category",
            index = reticulate::r_to_py(cell_names)
        )
    )
    var <- pd$DataFrame(index = reticulate::r_to_py(gene_names))

    adata <- ad$AnnData(X = s, obs = obs, var = var)
    reticulate::py_set_item(py_attr(adata, "layers"), "spliced", s)
    reticulate::py_set_item(py_attr(adata, "layers"), "unspliced", u)
    reticulate::py_set_item(py_attr(adata, "layers"), "velocity", v)

    reduction_key <- paste0("X_", tolower(reduction_name))
    reticulate::py_set_item(
        py_attr(adata, "obsm"),
        reduction_key,
        reticulate::r_to_py(emb)
    )

    sc <- reticulate::import("scanpy", convert = FALSE)
    scv <- reticulate::import("scvelo", convert = FALSE)
    sc$pp$neighbors(adata, use_rep = reduction_key)
    scv$pp$moments(adata, n_pcs = NULL, n_neighbors = NULL)
    scv$tl$velocity_graph(adata)

    vk <- cr$kernels$VelocityKernel(adata)
    py_call(vk, "compute_transition_matrix")

    g <- cr$estimators$GPCCA(vk)
    py_call(g, "compute_eigendecomposition")
    py_call(g, "compute_schur", method = "krylov")
    py_call(g, "compute_macrostates")
    py_call(g, "predict_terminal_states")
    py_call(g, "compute_fate_probabilities")

    term_states <- py_attr(g, "terminal_states")
    abs_probs <- py_attr(g, "fate_probabilities")

    extract_array <- function(x) {
        x_mat <- py_attr(x, "X", default = x)
        if (reticulate::py_module_available("scipy.sparse")) {
            sp_sparse <- reticulate::import("scipy.sparse", convert = FALSE)
            if (py_to_r(sp_sparse$issparse(x_mat))) {
                x_mat <- py_call(x_mat, "toarray")
            }
        }
        as.matrix(py_to_r(x_mat))
    }

    abs_mat <- tryCatch(
        extract_array(abs_probs),
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

    abs_prob_names <- colnames(abs_mat)
    if (is.null(abs_prob_names) && !is.null(abs_probs)) {
        abs_prob_names <- tryCatch({
            names_attr <- py_attr(abs_probs, "names")
            if (!is.null(names_attr)) {
                as.character(py_to_r(names_attr))
            } else {
                NULL
            }
        }, error = function(e) NULL)
    }
    lineage_drivers <- NULL
    tryCatch({
        py_call(g, "compute_lineage_drivers")
        lineage_drivers <- as.data.frame(py_to_r(py_attr(g, "lineage_drivers")))
        if (!"gene" %in% colnames(lineage_drivers)) {
            lineage_drivers$gene <- rownames(lineage_drivers)
        }
        lineage_drivers <- lineage_drivers[, c("gene", setdiff(colnames(lineage_drivers), "gene")), drop = FALSE]
        rownames(lineage_drivers) <- NULL
    }, error = function(e) {
        NULL
    })

    list(
        terminal_states = as.character(py_to_r(term_states)),
        absorption_probs = abs_mat,
        absorption_prob_names = abs_prob_names,
        lineage_drivers = lineage_drivers
    )
}

#' Run Cell Fate Inference
#'
#' Thin frontend for cell fate inference. Currently, `CellRank` is the only
#' supported backend and is used to estimate terminal states and fate
#' probabilities from RNA velocity results.
#'
#' @param sce A SingleCellExperiment object with velocity computed.
#' @param method Backend used for fate inference. Currently only `"CellRank"`.
#' @param reduction String specifying the reduction to use for kNN graph. Defaults to `"PCA"`.
#' @param cluster_key String specifying the column in `colData` containing cluster labels.
#' @param velocity_id Optional velocity record id. If `NULL`, uses the active
#'   velocity result.
#' @param backend Python execution backend passed to `RunCellRank()`.
#' @param python Optional Python executable or conda environment name used when
#'   `backend = "reticulate"`.
#' @param name trajectory record id. Defaults to `"cellfate"`.
#' @param ... Additional arguments passed to the selected backend.
#'
#' @return A SingleCellExperiment object with cell fate results recorded.
#' @export
RunCellFate <- function(
    sce,
    method = c("CellRank"),
    reduction = "PCA",
    cluster_key = ActiveIdent(sce),
    velocity_id = NULL,
    backend = c("basilisk", "reticulate"),
    python = NULL,
    name = "cellfate",
    ...
) {
    method <- match.arg(method)
    backend <- match.arg(backend)

    switch(
        method,
        CellRank = RunCellRank(
            sce = sce,
            reduction = reduction,
            cluster_key = cluster_key,
            velocity_id = velocity_id,
            backend = backend,
            python = python,
            name = name,
            ...
        )
    )
}

sclet_normalize_cellrank_driver_table <- function(drivers) {
    if (is.null(drivers) || !nrow(drivers)) {
        return(NULL)
    }

    if (!"gene" %in% colnames(drivers)) {
        drivers$gene <- rownames(drivers)
    }
    drivers <- drivers[, c("gene", setdiff(colnames(drivers), "gene")), drop = FALSE]
    if ("lineage" %in% colnames(drivers)) {
        return(drivers)
    }

    candidate_cols <- setdiff(colnames(drivers), "gene")
    split_candidates <- list(
        "__" = function(x) strsplit(x, "__", fixed = TRUE)[[1]],
        "." = function(x) strsplit(x, ".", fixed = TRUE)[[1]],
        "_" = function(x) strsplit(x, "_", fixed = TRUE)[[1]]
    )
    parsed_rows <- list()

    for (col in candidate_cols) {
        parts <- NULL
        for (parser in split_candidates) {
            candidate <- parser(col)
            if (length(candidate) >= 2) {
                parts <- candidate
                break
            }
        }
        if (is.null(parts)) {
            next
        }
        parsed_rows[[length(parsed_rows) + 1L]] <- data.frame(
            gene = drivers$gene,
            lineage = parts[[1]],
            metric = paste(parts[-1], collapse = "_"),
            value = drivers[[col]],
            stringsAsFactors = FALSE
        )
    }

    if (!length(parsed_rows)) {
        return(drivers)
    }

    parsed_long <- do.call(rbind, parsed_rows)
    tidyr::pivot_wider(parsed_long, names_from = "metric", values_from = "value")
}

sclet_resolve_cellrank_lineage <- function(cellrank, lineage = NULL) {
    drivers <- cellrank$lineage_drivers
    available_lineages <- NULL
    if (!is.null(drivers) && "lineage" %in% colnames(drivers)) {
        available_lineages <- unique(as.character(drivers$lineage))
    } else {
        available_lineages <- cellrank$fate_probability_names
    }
    available_lineages <- available_lineages[!is.na(available_lineages) & nzchar(available_lineages)]
    if (is.null(lineage)) {
        if (length(available_lineages) == 1) {
            return(available_lineages[[1]])
        }
        return(NULL)
    }
    lineage <- as.character(lineage)
    if (lineage %in% available_lineages) {
        lineage
    } else {
        NULL
    }
}

sclet_pick_cellrank_driver_metric <- function(drivers) {
    metric_cols <- setdiff(colnames(drivers), c("gene", "lineage"))
    metric_cols <- metric_cols[vapply(drivers[metric_cols], is.numeric, logical(1))]
    preferred <- c("corr", "correlation", "score", "stat", "pearson", "spearman")
    hit <- preferred[preferred %in% metric_cols]
    if (length(hit)) {
        return(hit[[1]])
    }
    if (length(metric_cols)) {
        return(metric_cols[[1]])
    }
    NULL
}

sclet_cellrank_fate_columns <- function(cellrank) {
    fate_cols <- cellrank$artifacts$fate_probability_cols
    if (is.null(fate_cols)) {
        fate_cols <- cellrank$fate_probability_cols
    }
    fate_cols
}

sclet_cellrank_fate_names <- function(cellrank) {
    fate_prob_names <- cellrank$artifacts$fate_probability_names
    if (is.null(fate_prob_names)) {
        fate_prob_names <- cellrank$fate_probability_names
    }
    fate_prob_names
}

sclet_resolve_cellrank_fate_col <- function(cellrank, fate = NULL) {
    fate_cols <- sclet_cellrank_fate_columns(cellrank)
    if (is.null(fate_cols) || length(fate_cols) == 0) {
        stop("No fate probability columns found. Please run RunCellRank() or RunCellFate() first.")
    }

    fate_prob_names <- sclet_cellrank_fate_names(cellrank)
    fate_name_map <- stats::setNames(fate_cols, fate_prob_names)
    fate_col <- fate
    if (is.null(fate_col)) {
        fate_col <- fate_cols[[1]]
    } else if (is.numeric(fate_col)) {
        if (length(fate_col) != 1 || is.na(fate_col) ||
            fate_col < 1 || fate_col > length(fate_cols)) {
            stop("Numeric `fate` must be a valid 1-based index into stored fate probability columns.")
        }
        fate_col <- fate_cols[[fate_col]]
    } else if (fate_col %in% names(fate_name_map)) {
        fate_col <- fate_name_map[[fate_col]]
    }
    if (!fate_col %in% fate_cols) {
        stop(
            "Unknown fate column '", fate_col, "'. Available columns are: ",
            paste(c(fate_cols, names(fate_name_map)), collapse = ", "), "."
        )
    }
    fate_col
}

#' Summarize CellRank Results
#'
#' Reports key diagnostics for stored CellRank results, including terminal-state
#' counts and fate-probability ranges. This is useful for detecting smoke-test
#' cases where CellRank returns a single terminal state or near-constant fate
#' probabilities.
#'
#' @param object A SingleCellExperiment object with CellRank results.
#' @param id Optional trajectory record id.
#'
#' @return A list with `summary`, `terminal_states`, and `fate_probabilities`.
#' @export
CellRankSummary <- function(object, id = NULL) {
    cellrank <- get_cellrank(object, id = id)
    if (is.null(cellrank)) {
        stop("No CellRank results found in the object.")
    }

    terminal_state_col <- cellrank$artifacts$terminal_state_col
    if (is.null(terminal_state_col)) {
        terminal_state_col <- cellrank$terminal_state_col
    }
    fate_cols <- sclet_cellrank_fate_columns(cellrank)

    terminal_df <- data.frame(
        terminal_state = character(),
        n_cells = integer(),
        stringsAsFactors = FALSE
    )
    if (!is.null(terminal_state_col) &&
        terminal_state_col %in% colnames(SummarizedExperiment::colData(object))) {
        terminal_state <- as.character(SummarizedExperiment::colData(object)[[terminal_state_col]])
        terminal_state[is.na(terminal_state) | !nzchar(terminal_state)] <- "Unassigned"
        terminal_tab <- table(terminal_state, useNA = "ifany")
        terminal_df <- data.frame(
            terminal_state = names(terminal_tab),
            n_cells = as.integer(terminal_tab),
            stringsAsFactors = FALSE
        )
    }

    fate_df <- data.frame(
        fate_probability_col = character(),
        min = numeric(),
        median = numeric(),
        max = numeric(),
        sd = numeric(),
        stringsAsFactors = FALSE
    )
    if (!is.null(fate_cols) && length(fate_cols)) {
        fate_df <- do.call(rbind, lapply(fate_cols, function(col) {
            values <- as.numeric(SummarizedExperiment::colData(object)[[col]])
            data.frame(
                fate_probability_col = col,
                min = min(values, na.rm = TRUE),
                median = stats::median(values, na.rm = TRUE),
                max = max(values, na.rm = TRUE),
                sd = stats::sd(values, na.rm = TRUE),
                stringsAsFactors = FALSE
            )
        }))
    }

    terminal_state_col_value <- if (is.null(terminal_state_col)) "" else terminal_state_col
    fate_cols_value <- if (is.null(fate_cols)) character() else fate_cols

    list(
        summary = data.frame(
            metric = c(
                "n_cells",
                "n_terminal_states",
                "n_fate_dimensions",
                "has_lineage_drivers",
                "terminal_state_col",
                "fate_probability_cols"
            ),
            value = c(
                ncol(object),
                length(setdiff(terminal_df$terminal_state, "Unassigned")),
                if (is.null(fate_cols)) 0L else length(fate_cols),
                isTRUE(cellrank$artifacts$has_lineage_drivers),
                terminal_state_col_value,
                paste(fate_cols_value, collapse = ";")
            ),
            stringsAsFactors = FALSE
        ),
        terminal_states = terminal_df,
        fate_probabilities = fate_df
    )
}

#' Correlate Velocity Magnitude with CellRank Fate Probabilities
#'
#' Computes correlations between per-cell velocity magnitude and CellRank fate
#' probabilities. This is a lightweight diagnostic for connecting RegVelo or
#' scVelo velocity strength to downstream fate assignments.
#'
#' @param object A SingleCellExperiment object with velocity and CellRank
#'   results.
#' @param velocity_id Optional velocity record id.
#' @param trajectory_id Optional CellRank trajectory record id.
#' @param assay Optional velocity assay name. If supplied, overrides
#'   `velocity_id`.
#' @param method Correlation method passed to `stats::cor.test()`.
#'
#' @return A data.frame with one row per fate probability column.
#' @export
VelocityFateCorrelation <- function(
    object,
    velocity_id = NULL,
    trajectory_id = NULL,
    assay = NULL,
    method = c("spearman", "pearson", "kendall")
) {
    method <- match.arg(method)
    cellrank <- get_cellrank(object, id = trajectory_id)
    if (is.null(cellrank)) {
        stop("No CellRank results found in the object.")
    }
    fate_cols <- sclet_cellrank_fate_columns(cellrank)
    if (is.null(fate_cols) || !length(fate_cols)) {
        stop("No fate probability columns found. Please run RunCellRank() or RunCellFate() first.")
    }

    velocity_magnitude <- VelocityMagnitude(
        object,
        id = velocity_id,
        assay = assay,
        margin = "cell"
    )
    velocity_assay <- sclet_velocity_assay_name(object, id = velocity_id, assay = assay)
    do.call(rbind, lapply(fate_cols, function(fate_col) {
        fate <- as.numeric(SummarizedExperiment::colData(object)[[fate_col]])
        ok <- stats::complete.cases(velocity_magnitude, fate)
        estimate <- NA_real_
        p_value <- NA_real_
        if (sum(ok) >= 3 &&
            stats::sd(velocity_magnitude[ok]) > 0 &&
            stats::sd(fate[ok]) > 0) {
            test <- stats::cor.test(velocity_magnitude[ok], fate[ok], method = method)
            estimate <- unname(test$estimate)
            p_value <- test$p.value
        }
        data.frame(
            fate_probability_col = fate_col,
            velocity_assay = velocity_assay,
            method = method,
            n_cells = sum(ok),
            estimate = estimate,
            p_value = p_value,
            stringsAsFactors = FALSE
        )
    }))
}

#' Plot Velocity Magnitude Against CellRank Fate Probability
#'
#' Draws a per-cell scatter plot linking velocity magnitude to one CellRank fate
#' probability. Cells can be colored by terminal state when available.
#'
#' @param object A SingleCellExperiment object with velocity and CellRank
#'   results.
#' @param fate Fate branch to plot. Can be a stored fate probability column, a
#'   lineage name, or a 1-based column index. If `NULL`, uses the first one.
#' @param velocity_id Optional velocity record id.
#' @param trajectory_id Optional CellRank trajectory record id.
#' @param assay Optional velocity assay name. If supplied, overrides
#'   `velocity_id`.
#' @param point_size Point size for cells.
#' @param alpha Point alpha.
#'
#' @return A ggplot object.
#' @importFrom rlang .data
#' @export
plot_velocity_fate_correlation <- function(
    object,
    fate = NULL,
    velocity_id = NULL,
    trajectory_id = NULL,
    assay = NULL,
    point_size = 0.8,
    alpha = 0.8
) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is needed for this function to work. Please install it.")
    }

    cellrank <- get_cellrank(object, id = trajectory_id)
    if (is.null(cellrank)) {
        stop("No CellRank results found in the object.")
    }
    fate_col <- sclet_resolve_cellrank_fate_col(cellrank, fate = fate)
    velocity_magnitude <- VelocityMagnitude(
        object,
        id = velocity_id,
        assay = assay,
        margin = "cell"
    )
    terminal_state_col <- cellrank$artifacts$terminal_state_col
    if (is.null(terminal_state_col)) {
        terminal_state_col <- cellrank$terminal_state_col
    }
    terminal_state <- rep("Unassigned", ncol(object))
    if (!is.null(terminal_state_col) &&
        terminal_state_col %in% colnames(SummarizedExperiment::colData(object))) {
        terminal_state <- as.character(SummarizedExperiment::colData(object)[[terminal_state_col]])
        terminal_state[is.na(terminal_state) | !nzchar(terminal_state)] <- "Unassigned"
    }

    plot_df <- data.frame(
        velocity_magnitude = as.numeric(velocity_magnitude),
        fate_probability = as.numeric(SummarizedExperiment::colData(object)[[fate_col]]),
        terminal_state = terminal_state,
        stringsAsFactors = FALSE
    )

    ggplot2::ggplot(
        plot_df,
        ggplot2::aes(
            x = .data$velocity_magnitude,
            y = .data$fate_probability,
            color = .data$terminal_state
        )
    ) +
        ggplot2::geom_point(size = point_size, alpha = alpha) +
        ggplot2::theme_classic() +
        ggplot2::labs(
            x = "Velocity magnitude",
            y = fate_col,
            color = "Terminal state",
            title = "Velocity Magnitude vs Fate Probability"
        )
}

#' Plot Cell Fate Probability
#'
#' Plot fate probabilities inferred by CellRank on a reduced-dimensional
#' embedding using the `plot_module_graphic` naming convention.
#'
#' @param object A SingleCellExperiment object with CellRank results.
#' @param fate Fate branch to plot. Can be a stored fate probability column
#'   name, a lineage name, or a 1-based column index. If `NULL`, uses the
#'   first one.
#' @param reduction Reduction to plot on. If `NULL`, prefers the reduction
#'   recorded in CellRank results and otherwise falls back to the active
#'   reduction.
#' @param id Optional trajectory record id.
#' @param point_size Point size for cells.
#' @param alpha Point alpha.
#'
#' @return A ggplot object.
#' @importFrom rlang .data
#' @export
plot_fate_probability <- function(
    object,
    fate = NULL,
    reduction = NULL,
    id = NULL,
    point_size = 0.5,
    alpha = 0.9
) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is needed for this function to work. Please install it.")
    }

    cellrank <- get_cellrank(object, id = id)
    if (is.null(cellrank)) {
        stop("No CellRank results found in the object.")
    }

    fate_cols <- sclet_cellrank_fate_columns(cellrank)
    if (is.null(fate_cols) || length(fate_cols) == 0) {
        stop("No fate probability columns found. Please run RunCellRank() or RunCellFate() first.")
    }

    if (is.null(reduction)) {
        reduction <- cellrank$inputs$reduction
    }
    if (is.null(reduction)) {
        reduction <- DefaultReduction(object)
    }
    if (is.null(reduction) || !reduction %in% SingleCellExperiment::reducedDimNames(object)) {
        stop("A valid reduction is required for plotting fate probabilities.")
    }

    fate_col <- sclet_resolve_cellrank_fate_col(cellrank, fate = fate)

    emb <- SingleCellExperiment::reducedDim(object, reduction)
    plot_df <- data.frame(
        x = emb[, 1],
        y = emb[, 2],
        fate_probability = SummarizedExperiment::colData(object)[[fate_col]]
    )

    ggplot2::ggplot(
        plot_df,
        ggplot2::aes(x = .data$x, y = .data$y, color = .data$fate_probability)
    ) +
        ggplot2::geom_point(size = point_size, alpha = alpha) +
        ggplot2::scale_color_viridis_c(option = "C") +
        ggplot2::theme_classic() +
        ggplot2::labs(
            x = paste0(reduction, " 1"),
            y = paste0(reduction, " 2"),
            color = fate_col,
            title = paste("Fate Probability:", fate_col)
        )
}

#' Plot Cell Fate Terminal States
#'
#' Plot terminal state assignments inferred by CellRank on a reduced-dimensional
#' embedding using the `plot_module_graphic` naming convention.
#'
#' @param object A SingleCellExperiment object with CellRank results.
#' @param reduction Reduction to plot on. If `NULL`, prefers the reduction
#'   recorded in CellRank results and otherwise falls back to the active
#'   reduction.
#' @param id Optional trajectory record id.
#' @param point_size Point size for cells.
#' @param alpha Point alpha.
#' @param na_color Color used for cells without terminal state labels.
#'
#' @return A ggplot object.
#' @importFrom rlang .data
#' @export
plot_fate_terminal_states <- function(
    object,
    reduction = NULL,
    id = NULL,
    point_size = 0.5,
    alpha = 0.9,
    na_color = "grey80"
) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is needed for this function to work. Please install it.")
    }

    cellrank <- get_cellrank(object, id = id)
    if (is.null(cellrank)) {
        stop("No CellRank results found in the object.")
    }

    terminal_state_col <- cellrank$artifacts$terminal_state_col
    if (is.null(terminal_state_col)) {
        terminal_state_col <- cellrank$terminal_state_col
    }
    if (is.null(terminal_state_col) || !terminal_state_col %in% colnames(SummarizedExperiment::colData(object))) {
        stop("No terminal state column found. Please run RunCellRank() or RunCellFate() first.")
    }

    if (is.null(reduction)) {
        reduction <- cellrank$inputs$reduction
    }
    if (is.null(reduction)) {
        reduction <- DefaultReduction(object)
    }
    if (is.null(reduction) || !reduction %in% SingleCellExperiment::reducedDimNames(object)) {
        stop("A valid reduction is required for plotting terminal states.")
    }

    emb <- SingleCellExperiment::reducedDim(object, reduction)
    terminal_states <- SummarizedExperiment::colData(object)[[terminal_state_col]]
    plot_df <- data.frame(
        x = emb[, 1],
        y = emb[, 2],
        terminal_state = terminal_states
    )
    plot_df$terminal_state <- as.character(plot_df$terminal_state)
    plot_df$terminal_state[is.na(plot_df$terminal_state) | !nzchar(plot_df$terminal_state)] <- "Unassigned"
    plot_df$terminal_state <- factor(plot_df$terminal_state)

    ggplot2::ggplot(
        plot_df,
        ggplot2::aes(x = .data$x, y = .data$y, color = .data$terminal_state)
    ) +
        ggplot2::geom_point(size = point_size, alpha = alpha) +
        ggplot2::scale_color_discrete(na.value = na_color) +
        ggplot2::theme_classic() +
        ggplot2::labs(
            x = paste0(reduction, " 1"),
            y = paste0(reduction, " 2"),
            color = "Terminal state",
            title = "Cell Fate Terminal States"
        )
}

#' Plot Cell Fate Driver Trends
#'
#' Plot expression trends of CellRank driver genes against pseudotime or fate
#' probability using the `plot_module_graphic` naming convention.
#'
#' @param object A SingleCellExperiment object with CellRank results.
#' @param lineage Optional lineage name. If `NULL`, uses the only lineage
#'   available; otherwise must match one of the stored CellRank lineages.
#' @param features Optional driver genes to plot. If `NULL`, picks the top
#'   genes for the selected lineage from stored CellRank driver statistics.
#' @param top_n Number of top driver genes to plot when `features = NULL`.
#' @param layer Expression layer to use. If `NULL`, prefers an active
#'   non-scaled layer.
#' @param x Axis used for trend plotting: `"auto"`, `"pseudotime"` or
#'   `"fate_probability"`.
#' @param id Optional trajectory record id.
#' @param alpha Point alpha.
#' @param point_size Point size.
#' @param line_size Smoothing line size.
#' @param se Whether to show confidence ribbons for smooth curves.
#' @param ncol Number of facet columns.
#'
#' @return A ggplot object.
#' @importFrom rlang .data
#' @export
plot_fate_driver_trends <- function(
    object,
    lineage = NULL,
    features = NULL,
    top_n = 6,
    layer = NULL,
    x = c("auto", "pseudotime", "fate_probability"),
    id = NULL,
    alpha = 0.15,
    point_size = 0.4,
    line_size = 0.6,
    se = FALSE,
    ncol = NULL
) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is needed for this function to work. Please install it.")
    }

    x <- match.arg(x)
    cellrank <- get_cellrank(object, id = id)
    if (is.null(cellrank)) {
        stop("No CellRank results found in the object.")
    }

    drivers <- cellrank$lineage_drivers
    if (is.null(drivers) || !nrow(drivers)) {
        stop("No lineage driver table found. Please rerun RunCellRank() with driver computation available.")
    }
    if (!"lineage" %in% colnames(drivers)) {
        stop("Stored lineage driver table does not contain a normalized 'lineage' column.")
    }

    lineage <- sclet_resolve_cellrank_lineage(cellrank, lineage = lineage)
    if (is.null(lineage)) {
        stop("Please provide a valid `lineage`. Multiple lineages are available in the stored CellRank results.")
    }

    lineage_drivers <- drivers[as.character(drivers$lineage) == lineage, , drop = FALSE]
    if (!nrow(lineage_drivers)) {
        stop("No driver genes found for lineage '", lineage, "'.")
    }

    if (is.null(features)) {
        rank_col <- sclet_pick_cellrank_driver_metric(lineage_drivers)
        if (is.null(rank_col)) {
            features <- utils::head(unique(lineage_drivers$gene), top_n)
        } else {
            ord <- order(abs(lineage_drivers[[rank_col]]), decreasing = TRUE, na.last = NA)
            features <- utils::head(unique(lineage_drivers$gene[ord]), top_n)
        }
    }
    features <- unique(as.character(features))
    features <- intersect(features, rownames(object))
    if (!length(features)) {
        stop("None of the requested driver genes are present in `rownames(object)`.")
    }

    source <- sclet_resolve_expression_source(
        object = object,
        layer = layer,
        assay = NULL,
        prefer_nonscaled = TRUE,
        fallback_layers = c("logcounts", "counts"),
        context = "plot_fate_driver_trends"
    )
    expr_mat <- SummarizedExperiment::assay(object, source$assay)[features, , drop = FALSE]
    expr_df <- as.data.frame(t(as.matrix(expr_mat)))
    colnames(expr_df) <- features

    fate_prob_cols <- cellrank$artifacts$fate_probability_cols
    if (is.null(fate_prob_cols)) {
        fate_prob_cols <- cellrank$fate_probability_cols
    }
    fate_prob_names <- cellrank$artifacts$fate_probability_names
    if (is.null(fate_prob_names)) {
        fate_prob_names <- cellrank$fate_probability_names
    }
    fate_name_map <- stats::setNames(fate_prob_cols, fate_prob_names)
    fate_col <- fate_name_map[[lineage]]
    if (is.null(fate_col) || !fate_col %in% colnames(SummarizedExperiment::colData(object))) {
        stop("Could not map lineage '", lineage, "' to a stored fate probability column.")
    }

    x_var <- NULL
    x_label <- NULL
    if (x == "auto") {
        x <- if ("velocity_pseudotime" %in% colnames(SummarizedExperiment::colData(object))) "pseudotime" else "fate_probability"
    }
    if (x == "pseudotime") {
        if (!"velocity_pseudotime" %in% colnames(SummarizedExperiment::colData(object))) {
            stop("`velocity_pseudotime` not found in colData(object); use `x = \"fate_probability\"` instead.")
        }
        x_var <- SummarizedExperiment::colData(object)[["velocity_pseudotime"]]
        x_label <- "Velocity pseudotime"
    } else {
        x_var <- SummarizedExperiment::colData(object)[[fate_col]]
        x_label <- paste("Fate probability:", lineage)
    }

    plot_df <- data.frame(
        .x = x_var,
        .fate = SummarizedExperiment::colData(object)[[fate_col]],
        expr_df,
        check.names = FALSE
    )
    plot_df$cell_id <- colnames(object)
    plot_long <- tidyr::pivot_longer(
        plot_df,
        cols = dplyr::all_of(features),
        names_to = "gene",
        values_to = "expression"
    )

    ggplot2::ggplot(
        plot_long,
        ggplot2::aes(x = .data$.x, y = .data$expression)
    ) +
        ggplot2::geom_point(
            ggplot2::aes(color = .data$.fate),
            alpha = alpha,
            size = point_size
        ) +
        ggplot2::geom_smooth(
            color = "black",
            se = se,
            linewidth = line_size,
            method = "loess",
            formula = y ~ x
        ) +
        ggplot2::scale_color_viridis_c(option = "C") +
        ggplot2::theme_classic() +
        ggplot2::facet_wrap(~ gene, scales = "free_y", ncol = ncol) +
        ggplot2::labs(
            x = x_label,
            y = paste0("Expression (", source$assay, ")"),
            color = paste("Fate:", lineage),
            title = paste("Driver Trends:", lineage)
        )
}
