#' Run RegVelo RNA velocity analysis
#'
#' Run RegVelo through either an isolated `basilisk` Python environment or the
#' active `reticulate` Python environment and register its outputs as a sclet
#' velocity state. RegVelo combines spliced/unspliced counts with a prior gene
#' regulatory network.
#'
#' @param sce A SingleCellExperiment object.
#' @param grn Prior gene regulatory network. Either a numeric matrix with
#'   targets in rows and regulators in columns, or a data.frame with
#'   `target`, `regulator`, and optional `weight` columns.
#' @param regulators Optional character vector of regulator genes. If `NULL`,
#'   regulators are inferred from `grn`.
#' @param spliced_assay,unspliced_assay Assays containing spliced and unspliced
#'   counts.
#' @param reduction Optional dimensionality reduction copied to AnnData obsm.
#'   Defaults to `DefaultReduction(sce)`.
#' @param name Velocity record id. Defaults to `"regvelo"`.
#' @param max_epochs Maximum training epochs passed to RegVelo.
#' @param lam,lam2 RegVelo regularization parameters.
#' @param batch_size Optional training batch size.
#' @param seed Optional random seed.
#' @param backend Python execution backend. `"basilisk"` creates the packaged
#'   CPU environment; `"reticulate"` uses the caller-configured Python, which is
#'   preferred on GPU clusters.
#' @param python Optional Python executable or conda environment name used when
#'   `backend = "reticulate"`.
#' @param ... Additional parameters passed to the RegVelo model constructor.
#'
#' @return A SingleCellExperiment object with RegVelo velocity state updated.
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom SummarizedExperiment colData colData<- rowData rowData<-
#' @export
RunRegVelo <- function(
    sce,
    grn,
    regulators = NULL,
    spliced_assay = "spliced",
    unspliced_assay = "unspliced",
    reduction = NULL,
    name = "regvelo",
    max_epochs = 1500L,
    lam = 1,
    lam2 = 0,
    batch_size = NULL,
    seed = 0L,
    backend = c("basilisk", "reticulate"),
    python = NULL,
    ...
) {
    backend <- match.arg(backend)
    if (identical(backend, "basilisk") && !requireNamespace("basilisk", quietly = TRUE)) {
        cli::cli_abort(
            "Please install {.pkg basilisk} to run RegVelo with {.val backend = 'basilisk'}.",
            class = "sclet_package_missing"
        )
    }
    if (identical(backend, "reticulate") && !requireNamespace("reticulate", quietly = TRUE)) {
        cli::cli_abort(
            "Please install {.pkg reticulate} to run RegVelo with {.val backend = 'reticulate'}.",
            class = "sclet_package_missing"
        )
    }

    assay_names <- SummarizedExperiment::assayNames(sce)
    missing_assays <- setdiff(c(spliced_assay, unspliced_assay), assay_names)
    if (length(missing_assays)) {
        cli::cli_abort(
            "Required assay(s) not found: {.val {missing_assays}}.",
            class = "sclet_missing_assay"
        )
    }

    if (is.null(reduction)) {
        reduction <- DefaultReduction(sce)
    }
    if (!is.null(reduction) && !reduction %in% SingleCellExperiment::reducedDimNames(sce)) {
        cli::cli_abort(
            "Reduction {.val {reduction}} not found in {.cls SingleCellExperiment}.",
            class = "sclet_missing_reduction"
        )
    }

    prior <- sclet_prepare_regvelo_prior(grn, rownames(sce), regulators = regulators)
    features <- prior$genes
    if (length(features) == 0) {
        cli::cli_abort(
            "No overlapping genes between {.arg grn} and {.arg sce}.",
            class = "sclet_regvelo_empty_prior"
        )
    }
    if (length(prior$regulators) == 0) {
        cli::cli_abort(
            "No valid regulator genes found in {.arg grn}.",
            class = "sclet_regvelo_empty_regulators"
        )
    }

    spliced <- sclet_extract_cell_feature_matrix(sce, spliced_assay, features = features)
    unspliced <- sclet_extract_cell_feature_matrix(sce, unspliced_assay, features = features)

    tmp_dir <- tempfile("regvelo_")
    dir.create(tmp_dir)
    keep_workdir <- isTRUE(as.logical(Sys.getenv("SCLET_REGVELO_KEEP_WORKDIR", "FALSE"))) ||
        identical(Sys.getenv("SCLET_REGVELO_KEEP_WORKDIR", ""), "1")
    if (!keep_workdir) {
        on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
    } else {
        cli::cli_alert_info("Keeping RegVelo workdir for debugging: {.path {tmp_dir}}")
    }

    Matrix::writeMM(spliced, file.path(tmp_dir, "spliced.mtx"))
    Matrix::writeMM(unspliced, file.path(tmp_dir, "unspliced.mtx"))
    utils::write.csv(prior$W, file.path(tmp_dir, "prior_grn.csv"), row.names = TRUE)
    writeLines(colnames(sce), file.path(tmp_dir, "cell_names.txt"))
    writeLines(features, file.path(tmp_dir, "gene_names.txt"))
    writeLines(prior$regulators, file.path(tmp_dir, "regulators.txt"))
    writeLines(if (is.null(reduction)) "" else reduction, file.path(tmp_dir, "reduction.txt"))
    if (!is.null(reduction)) {
        utils::write.csv(
            as.matrix(SingleCellExperiment::reducedDim(sce, reduction)),
            file.path(tmp_dir, "embedding.csv"),
            row.names = FALSE
        )
    }

    model_params <- list(...)

    cli::cli_alert_info("Running RegVelo via {.val {backend}}...")
    if (identical(backend, "basilisk")) {
        cli::cli_alert_info("(First run will take several minutes to set up the Python environment)")
    }
    prev_state <- sclet_get_state(sce)

    rv_res <- sclet_run_regvelo_backend(
        backend = backend,
        python = python,
        workdir = tmp_dir,
        max_epochs = as.integer(max_epochs),
        lam = lam,
        lam2 = lam2,
        batch_size = batch_size,
        seed = if (is.null(seed)) NULL else as.integer(seed),
        model_params = model_params
    )
    sce <- sclet_restore_state(sce, prev_state)

    velocity_assay <- paste0(name, "_velocity")
    SummarizedExperiment::assay(sce, velocity_assay) <- sclet_gene_cell_matrix(
        rv_res$velocity,
        rownames_target = rownames(sce),
        colnames_target = colnames(sce),
        features = features
    )

    coldata_fields <- character()
    if (!is.null(rv_res$latent_time)) {
        latent_col <- paste0(name, "_latent_time")
        SummarizedExperiment::colData(sce)[[latent_col]] <- rv_res$latent_time
        coldata_fields <- c(coldata_fields, latent_col)
    }
    if (!is.null(rv_res$velocity_confidence)) {
        confidence_col <- paste0(name, "_velocity_confidence")
        SummarizedExperiment::colData(sce)[[confidence_col]] <- rv_res$velocity_confidence
        coldata_fields <- c(coldata_fields, confidence_col)
    }
    if (!is.null(rv_res$directional_uncertainty)) {
        uncertainty_col <- paste0(name, "_directional_uncertainty")
        SummarizedExperiment::colData(sce)[[uncertainty_col]] <- rv_res$directional_uncertainty
        coldata_fields <- c(coldata_fields, uncertainty_col)
    }

    rowdata_fields <- character()
    for (field in names(rv_res$rowdata)) {
        row_col <- paste0(name, "_", field)
        values <- rep(NA_real_, nrow(sce))
        names(values) <- rownames(sce)
        values[features] <- rv_res$rowdata[[field]]
        SummarizedExperiment::rowData(sce)[[row_col]] <- values
        rowdata_fields <- c(rowdata_fields, row_col)
    }

    regvelo_analysis <- list(
        id = name,
        method = "regvelo::REGVELOVI",
        inputs = list(
            assays = c(spliced = spliced_assay, unspliced = unspliced_assay),
            reduction = reduction,
            grn = list(
                n_targets = length(features),
                n_regulators = length(prior$regulators),
                n_edges = sum(prior$W != 0)
            )
        ),
        artifacts = list(
            velocity_assay = velocity_assay,
            colData = coldata_fields,
            rowData = rowdata_fields,
            analysis_key = name
        ),
        params = list(
            max_epochs = max_epochs,
            lam = lam,
            lam2 = lam2,
            batch_size = batch_size,
            seed = seed,
            backend = backend,
            python = python
        ),
        summary = list(
            n_features = length(features),
            n_cells = ncol(sce),
            n_regulators = length(prior$regulators)
        ),
        results = rv_res,
        created_at = Sys.time()
    )
    sce <- sclet_set_analysis(sce, "velocity", regvelo_analysis)
    sce <- sclet_set_state_record(
        object = sce,
        type = "velocity",
        id = name,
        active = TRUE,
        value = regvelo_analysis
    )
    sce <- sclet_log_command(
        sce,
        "RunRegVelo",
        params = list(
            spliced_assay = spliced_assay,
            unspliced_assay = unspliced_assay,
            reduction = reduction,
            name = name,
            max_epochs = max_epochs,
            lam = lam,
            lam2 = lam2,
            batch_size = batch_size,
            seed = seed,
            backend = backend,
            python = python
        ),
        outputs = list(
            analysis = "velocity",
            velocity = name,
            assay = velocity_assay
        )
    )

    sce
}

sclet_prepare_regvelo_prior <- function(grn, genes, regulators = NULL) {
    genes <- as.character(genes)

    if (is.matrix(grn) || methods::is(grn, "Matrix")) {
        W <- as.matrix(grn)
        if (is.null(rownames(W)) || is.null(colnames(W))) {
            stop("Matrix `grn` must have rownames as targets and colnames as regulators.")
        }
        W[is.na(W)] <- 0
        keep_targets <- intersect(genes, rownames(W))
        if (is.null(regulators)) {
            regulators <- colnames(W)
        }
        regulators <- intersect(as.character(regulators), genes)
        regulators <- intersect(regulators, colnames(W))
        W <- W[keep_targets, regulators, drop = FALSE]
        return(list(W = W, genes = rownames(W), regulators = colnames(W)))
    }

    if (!is.data.frame(grn)) {
        stop("`grn` must be a matrix-like object or a data.frame edge table.")
    }

    required <- c("target", "regulator")
    if (!all(required %in% colnames(grn))) {
        stop("Data-frame `grn` must contain columns: target, regulator.")
    }
    if (!"weight" %in% colnames(grn)) {
        grn$weight <- 1
    }

    grn$target <- as.character(grn$target)
    grn$regulator <- as.character(grn$regulator)
    grn$weight <- as.numeric(grn$weight)
    grn$weight[is.na(grn$weight)] <- 0
    grn <- grn[grn$target %in% genes & grn$regulator %in% genes, , drop = FALSE]

    targets <- intersect(genes, unique(grn$target))
    if (is.null(regulators)) {
        regulators <- unique(grn$regulator)
    }
    regulators <- intersect(as.character(regulators), genes)
    regulators <- intersect(regulators, unique(grn$regulator))

    W <- matrix(0, nrow = length(targets), ncol = length(regulators),
                dimnames = list(targets, regulators))
    if (nrow(grn) && length(targets) && length(regulators)) {
        grn <- grn[grn$target %in% targets & grn$regulator %in% regulators, , drop = FALSE]
        idx <- cbind(match(grn$target, targets), match(grn$regulator, regulators))
        W[idx] <- grn$weight
    }

    list(W = W, genes = rownames(W), regulators = colnames(W))
}

sclet_gene_cell_matrix <- function(mat, rownames_target, colnames_target, features) {
    out <- Matrix::Matrix(0, nrow = length(rownames_target), ncol = length(colnames_target), sparse = TRUE)
    rownames(out) <- rownames_target
    colnames(out) <- colnames_target
    if (is.null(mat)) {
        return(out)
    }
    mat <- sclet_as_dgCMatrix(mat)
    expected_cell_gene <- c(length(colnames_target), length(features))
    expected_gene_cell <- rev(expected_cell_gene)
    if (identical(dim(mat), expected_cell_gene)) {
        mat <- Matrix::t(mat)
    } else if (!identical(dim(mat), expected_gene_cell)) {
        stop(
            "Unexpected RegVelo velocity matrix dimensions: got ",
            paste(dim(mat), collapse = " x "),
            ", expected cells x genes (", paste(expected_cell_gene, collapse = " x "),
            ") or genes x cells (", paste(expected_gene_cell, collapse = " x "), ")."
        )
    }
    rownames(mat) <- features
    colnames(mat) <- colnames_target
    out[features, colnames_target] <- mat
    out
}

sclet_run_regvelo_backend <- function(
    backend,
    python,
    workdir,
    max_epochs,
    lam,
    lam2,
    batch_size,
    seed,
    model_params
) {
    if (identical(backend, "basilisk")) {
        proc <- basilisk::basiliskStart(sclet_regvelo_env)
        on.exit(basilisk::basiliskStop(proc), add = TRUE)
        return(basilisk::basiliskRun(
            proc,
            sclet_run_regvelo_python,
            workdir = workdir,
            max_epochs = max_epochs,
            lam = lam,
            lam2 = lam2,
            batch_size = batch_size,
            seed = seed,
            model_params = model_params
        ))
    }

    if (!is.null(python) && !reticulate::py_available(initialize = FALSE)) {
        if (file.exists(python)) {
            reticulate::use_python(python, required = TRUE)
        } else {
            reticulate::use_condaenv(python, required = TRUE)
        }
    }
    sclet_run_regvelo_python(
        workdir = workdir,
        max_epochs = max_epochs,
        lam = lam,
        lam2 = lam2,
        batch_size = batch_size,
        seed = seed,
        model_params = model_params
    )
}

sclet_run_regvelo_python <- function(workdir, max_epochs, lam, lam2, batch_size, seed, model_params) {
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
    py_get_item <- function(x, key, default = NULL) {
        tryCatch(reticulate::py_get_item(x, key), error = function(e) default)
    }
    py_set_item <- function(x, key, value) {
        reticulate::py_set_item(x, key, value)
    }
    reticulate::py_run_string(
        paste(
            "def sclet_py_keys(x):",
            "    try:",
            "        return list(x.keys())",
            "    except Exception:",
            "        return []",
            "",
            "def sclet_py_velocity_attrs(x):",
            "    return [name for name in dir(x) if 'vel' in name.lower() or 'latent' in name.lower() or 'time' in name.lower()]",
            "",
            "def sclet_py_columns(x):",
            "    try:",
            "        return x.columns.astype(str).tolist()",
            "    except Exception:",
            "        return []",
            sep = "\n"
        )
    )
    py_keys <- function(x) {
        as.character(py_to_r(reticulate::py$sclet_py_keys(x)))
    }
    py_velocity_attrs <- function(x) {
        as.character(py_to_r(reticulate::py$sclet_py_velocity_attrs(x)))
    }
    py_columns <- function(x) {
        as.character(py_to_r(reticulate::py$sclet_py_columns(x)))
    }
    py_call_any <- function(fun, attempts) {
        last_error <- NULL
        for (args in attempts) {
            result <- tryCatch(
                do.call(fun, args),
                error = function(e) {
                    last_error <<- e
                    NULL
                }
            )
            if (!is.null(result)) {
                return(result)
            }
        }
        stop(last_error)
    }
    py_import_regvelo_model <- function() {
        rv <- reticulate::import("regvelo", convert = FALSE)
        model <- py_attr(rv, "REGVELOVI")
        if (!is.null(model)) {
            return(model)
        }
        rv_model <- reticulate::import("regvelo._model", convert = FALSE)
        py_attr(rv_model, "REGVELOVI")
    }

    if (!is.null(seed)) {
        reticulate::import("numpy", convert = FALSE)$random$seed(seed)
        torch <- reticulate::import("torch", convert = FALSE)
        torch$manual_seed(seed)
    }

    ad <- reticulate::import("anndata", convert = FALSE)
    pd <- reticulate::import("pandas", convert = FALSE)
    np <- reticulate::import("numpy", convert = FALSE)
    sc <- reticulate::import("scanpy", convert = FALSE)
    scv <- reticulate::import("scvelo", convert = FALSE)
    scio <- reticulate::import("scipy.io", convert = FALSE)
    REGVELOVI <- py_import_regvelo_model()

    cell_names <- readLines(file.path(workdir, "cell_names.txt"))
    gene_names <- readLines(file.path(workdir, "gene_names.txt"))
    regulators <- readLines(file.path(workdir, "regulators.txt"))
    reduction <- readLines(file.path(workdir, "reduction.txt"))
    if (!length(reduction) || !nzchar(reduction[[1]])) {
        reduction <- NULL
    } else {
        reduction <- reduction[[1]]
    }

    spliced <- py_call(scio$mmread(file.path(workdir, "spliced.mtx")), "tocsr")
    unspliced <- py_call(scio$mmread(file.path(workdir, "unspliced.mtx")), "tocsr")
    obs <- pd$DataFrame(index = reticulate::r_to_py(cell_names))
    var <- pd$DataFrame(index = reticulate::r_to_py(gene_names))
    adata <- ad$AnnData(X = spliced, obs = obs, var = var)
    py_set_item(py_attr(adata, "layers"), "spliced", spliced)
    py_set_item(py_attr(adata, "layers"), "unspliced", unspliced)
    if (!is.null(reduction)) {
        emb <- as.matrix(utils::read.csv(file.path(workdir, "embedding.csv")))
        py_set_item(py_attr(adata, "obsm"), paste0("X_", tolower(reduction)), reticulate::r_to_py(emb))
    }

    prior <- utils::read.csv(file.path(workdir, "prior_grn.csv"), row.names = 1, check.names = FALSE)
    W <- np$array(as.matrix(prior), dtype = "float32")

    sc$pp$normalize_total(adata)
    sc$pp$log1p(adata)
    if (!is.null(reduction)) {
        sc$pp$neighbors(adata, use_rep = paste0("X_", tolower(reduction)))
    } else {
        sc$pp$pca(adata)
        sc$pp$neighbors(adata)
    }
    scv$pp$moments(adata, n_pcs = NULL, n_neighbors = NULL)

    setup_fun <- function(...) REGVELOVI$setup_anndata(adata, ...)
    py_call_any(
        setup_fun,
        list(
            list(spliced_layer = "Ms", unspliced_layer = "Mu"),
            list(spliced_layer = "spliced", unspliced_layer = "unspliced"),
            list()
        )
    )

    base_params <- utils::modifyList(list(
        W = W,
        regulators = reticulate::r_to_py(regulators),
        lam = lam,
        lam2 = lam2
    ), model_params)
    model <- py_call_any(
        function(...) do.call(REGVELOVI, list(...)),
        list(
            c(list(adata), base_params),
            c(list(adata), base_params[setdiff(names(base_params), "lam2")]),
            c(list(adata), base_params[setdiff(names(base_params), c("lam", "lam2"))]),
            c(list(adata), list(W = W)),
            list(adata)
        )
    )
    train_args <- list(max_epochs = max_epochs)
    if (!is.null(batch_size)) {
        train_args$batch_size <- as.integer(batch_size)
    }
    do.call(model$train, train_args)

    add_outputs <- py_attr(model, "add_regvelo_outputs_to_adata")
    if (!is.null(add_outputs)) {
        updated_adata <- py_call_any(
            function(...) reticulate::py_call(add_outputs, ...),
            list(
                list(),
                list(adata),
                list(adata = adata)
            )
        )
        if (!is.null(updated_adata) && reticulate::py_has_attr(updated_adata, "layers")) {
            adata <- updated_adata
        }
    } else {
        get_velocity <- py_attr(model, "get_velocity")
        if (is.null(get_velocity)) {
            stop("RegVelo model does not expose add_regvelo_outputs_to_adata() or get_velocity().")
        }
        py_set_item(
            py_attr(adata, "layers"),
            "velocity",
            py_call_any(
                function(...) reticulate::py_call(get_velocity, ...),
                list(list(), list(adata), list(adata = adata))
            )
        )
    }

    layers <- py_attr(adata, "layers")
    obs <- py_attr(adata, "obs")
    var <- py_attr(adata, "var")
    obsm <- py_attr(adata, "obsm")

    layer_keys <- py_keys(layers)
    velocity_candidates <- c(
        "velocity",
        "velocity_regvelo",
        "regvelo_velocity",
        "velocities",
        "velocity_u",
        "velocity_unspliced"
    )
    velocity <- NULL
    for (candidate in velocity_candidates) {
        velocity <- py_get_item(layers, candidate)
        if (!is.null(velocity)) {
            break
        }
    }
    if (is.null(velocity)) {
        get_velocity <- py_attr(model, "get_velocity")
        if (!is.null(get_velocity)) {
            velocity <- tryCatch(
                py_call_any(
                    function(...) reticulate::py_call(get_velocity, ...),
                    list(list(), list(adata), list(adata = adata))
                ),
                error = function(e) NULL
            )
            if (!is.null(velocity)) {
                py_set_item(layers, "velocity", velocity)
            }
        }
    }
    if (is.null(velocity)) {
        debug_lines <- c(
            "SCLET_REGVELO_DEBUG_VERSION=2026-07-01-velocity-keys",
            paste0("workdir=", workdir),
            paste0("layers=", paste(layer_keys, collapse = ",")),
            paste0("obs_columns=", paste(py_columns(obs), collapse = ",")),
            paste0("var_columns=", paste(py_columns(var), collapse = ",")),
            paste0("obsm_keys=", paste(py_keys(obsm), collapse = ",")),
            paste0("model_velocity_attrs=", paste(py_velocity_attrs(model), collapse = ","))
        )
        writeLines(debug_lines, file.path(workdir, "regvelo_debug.txt"))
        stop(
            "RegVelo did not return a velocity layer. ",
            "SCLET_REGVELO_DEBUG_VERSION=2026-07-01-velocity-keys. ",
            "Debug file: ", file.path(workdir, "regvelo_debug.txt"), "\n",
            "Available AnnData layers: ", paste(layer_keys, collapse = ", "), "\n",
            "Available obs columns: ", paste(py_columns(obs), collapse = ", "), "\n",
            "Available var columns: ", paste(py_columns(var), collapse = ", "), "\n",
            "Available obsm keys: ", paste(py_keys(obsm), collapse = ", "), "\n",
            "Velocity-like model attributes: ", paste(py_velocity_attrs(model), collapse = ", ")
        )
    }
    if (reticulate::py_has_attr(velocity, "toarray")) {
        velocity <- py_call(velocity, "toarray")
    }
    velocity <- as.matrix(py_to_r(velocity))

    obs_names <- as.character(py_to_r(py_attr(obs, "columns")))
    obs_vec <- function(candidates) {
        hit <- candidates[candidates %in% obs_names]
        if (!length(hit)) {
            return(NULL)
        }
        as.numeric(py_to_r(py_get_item(obs, hit[[1]])))
    }
    var_names <- as.character(py_to_r(py_attr(var, "columns")))
    var_vec <- function(candidates) {
        hit <- candidates[candidates %in% var_names]
        if (!length(hit)) {
            return(NULL)
        }
        as.numeric(py_to_r(py_get_item(var, hit[[1]])))
    }

    rowdata <- list()
    for (nm in c("alpha", "alpha_1", "beta", "gamma", "fit_t")) {
        val <- var_vec(c(nm, paste0(nm, "_regvelo"), paste0("regvelo_", nm)))
        if (!is.null(val)) {
            rowdata[[nm]] <- val
        }
    }

    list(
        velocity = velocity,
        latent_time = obs_vec(c("latent_time_regvelo", "regvelo_latent_time", "fit_t")),
        velocity_confidence = obs_vec(c("velocity_confidence", "velocity_confidence_regvelo")),
        directional_uncertainty = obs_vec(c("directional_uncertainty", "directional_uncertainty_regvelo")),
        rowdata = rowdata,
        regulators = regulators,
        genes = gene_names
    )
}
