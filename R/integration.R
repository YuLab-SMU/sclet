#' Run Integration
#' 
#' Unified entry point for batch correction and integration methods.
#' 
#' @title RunIntegration
#' @param object A SingleCellExperiment object
#' @param method Integration method to use. Currently supports "fastMNN", "Harmony", and "scVI".
#' @param batch The column name in colData representing the batch variable.
#' @param features Features to use for integration. If NULL, automatically uses VariableFeatures.
#' @param layer Layer to use for expression-based integration (e.g. fastMNN). If NULL, uses DefaultLayer.
#' @param reduction Reduction to use for reduction-based integration (e.g. Harmony). If NULL, uses DefaultReduction.
#' @param name Integration record id. Defaults to the method name.
#' @param seed integer seed for the stochastic backends. `harmony::RunHarmony()`
#' @param device Which device the torch-backed methods use: `"cpu"` (the
#'   default), `"cuda"`, or `"auto"`. Only `method = "scVI"` has a torch backend.
#'   A session-wide default can be set with `options(sclet.python.device = )`
#'   or `SCLET_PYTHON_DEVICE`.
#' and scVI's training both initialise from random numbers, so the same input
#' otherwise returns a different embedding on each call; the seed is recorded
#' with the integration state so the route can actually be reproduced.
#' fastMNN is deterministic and ignores it.
#' @param ... Additional arguments passed to the underlying integration functions.
#' @return A SingleCellExperiment object with integration results recorded.
#' @export
#' @examples
#' sce <- scuttle::mockSCE(ncells = 20, ngenes = 50)
#' sce <- NormalizeData(sce)
#' sce <- RunPCA(sce, ncomponents = 5)
#' colData(sce)$batch <- rep(c("b1", "b2"), length.out = ncol(sce))
#' \dontrun{
#' # needs the 'harmony' Suggests package (scVI also needs a Python environment)
#' sce <- RunIntegration(sce, method = "Harmony", batch = "batch")
#' SingleCellExperiment::reducedDimNames(sce)
#' }
RunIntegration <- function(object, method = c("fastMNN", "Harmony", "scVI"), batch = NULL, features = NULL, layer = NULL, reduction = NULL, name = NULL, seed = 1L, device = NULL, ...) {
    method <- match.arg(method)
    
    if (is.null(batch)) {
        # Batch is a design role, and the AI layer refuses design-gated actions
        # until the role is declared (ConfirmAIDesignSemantics). Point at that
        # rather than only saying a name is missing, so the caller learns the
        # channel it has to go through.
        stop(
            "`batch` must name the colData column that is the technical batch. ",
            "If you are driving this through the AI layer, declare the role first ",
            "with ConfirmAIDesignSemantics(object, design = list(batch = \"<column>\")), ",
            "or see check_integration_readiness() for what is missing.",
            call. = FALSE
        )
    }
    
    if (!batch %in% colnames(SummarizedExperiment::colData(object))) {
        stop(paste("Batch column", batch, "not found in colData."))
    }
    
    if (is.null(name)) {
        name <- tolower(method)
    }
    
    if (method == "fastMNN") {
        # fastMNN uses BatchRemover internally
        # Prepare inputs
        source <- sclet_resolve_expression_source(object, layer = layer, context = "RunIntegration(fastMNN)")
        
        # BatchRemover uses 'sce$batch' if batch is a vector, or looks for it
        # Actually BatchRemover expects `batch` as a vector or NULL. If NULL, it uses sce$batch.
        batch_vec <- SummarizedExperiment::colData(object)[[batch]]
        
        # Run BatchRemover
        object <- BatchRemover(sce = object, batch = batch_vec, HVG = features, assay.type = source$assay, name = name, ...)
        
    } else if (method == "Harmony") {
        if (!requireNamespace("harmony", quietly = TRUE)) {
            stop("Package 'harmony' is required. Please install it.")
        }
        
        if (is.null(reduction)) {
            reduction <- DefaultReduction(object)
        }
        if (is.null(reduction) || !reduction %in% SingleCellExperiment::reducedDimNames(object)) {
            reduction <- "PCA"
        }
        if (!reduction %in% SingleCellExperiment::reducedDimNames(object)) {
            object <- RunPCA(
                object,
                subset_row = features,
                layer = layer
            )
            reduction <- "PCA"
        }
        
        reduction_state <- sclet_get_state_record(object, "reduction", tolower(reduction))

        # RunHarmony.SingleCellExperiment operates on the "PCA" reducedDim and
        # does not accept Seurat-style `reduction=` arguments. If the user asks
        # for another source reduction, temporarily mirror it into "PCA".
        had_pca <- "PCA" %in% SingleCellExperiment::reducedDimNames(object)
        original_pca <- if (had_pca) SingleCellExperiment::reducedDim(object, "PCA") else NULL
        if (!identical(reduction, "PCA")) {
            SingleCellExperiment::reducedDim(object, "PCA") <- SingleCellExperiment::reducedDim(object, reduction)
        }

        object <- sclet_with_seed(seed, {
            harmony::RunHarmony(
                object,
                group.by.vars = batch,
                reduction.save = "HARMONY",
                ...
            )
        })

        if (!identical(reduction, "PCA")) {
            if (had_pca) {
                SingleCellExperiment::reducedDim(object, "PCA") <- original_pca
            } else {
                SingleCellExperiment::reducedDim(object, "PCA") <- NULL
            }
        }
        
        # Register integration state
        # Harmony outputs a reduction, not a layer
        object <- sclet_set_analysis_state(
            object = object,
            type = "integration",
            id = name,
            method = "harmony::RunHarmony",
            inputs = list(
                reduction = reduction,
                layer = if (!is.null(reduction_state)) reduction_state$inputs$layer else NULL,
                batch = batch
            ),
            artifacts = list(
                reduction = "HARMONY"
            ),
            params = list(
                seed = seed
            ),
            summary = list(
                batch_var = batch,
                n_components = ncol(SingleCellExperiment::reducedDim(object, "HARMONY"))
            )
        )
        
        # Register it as a reduction state as well, so downstream tools like RunUMAP can trace it
        object <- sclet_set_analysis_state(
            object = object,
            type = "reduction",
            id = "harmony",
            method = "harmony::RunHarmony",
            inputs = list(
                reduction = reduction,
                layer = if (!is.null(reduction_state)) reduction_state$inputs$layer else NULL,
                integration = name
            ),
            artifacts = list(
                reduction = "HARMONY"
            ),
            summary = list(
                n_components = ncol(SingleCellExperiment::reducedDim(object, "HARMONY"))
            )
        )
        
        object <- sclet_set_active_reduction(object, "HARMONY")
        
        object <- sclet_log_command(
            object,
            "RunIntegration",
            params = list(method = method, batch = batch, reduction = reduction, seed = seed),
            outputs = list(integration = name, reduction = "HARMONY")
        )
    } else if (method == "scVI") {
        if (!requireNamespace("basilisk", quietly = TRUE)) {
            stop("Package 'basilisk' is required for scVI.")
        }
        if (!requireNamespace("zellkonverter", quietly = TRUE)) {
            stop("Package 'zellkonverter' is required for scVI.")
        }
        
        # Determine features
        if (is.null(features)) {
            features <- VariableFeatures(object)
        }
        if (length(features) == 0) {
            stop("No VariableFeatures found. Please provide 'features' or run FindVariableFeatures.")
        }
        
        # scVI typically uses raw counts.
        if (!"counts" %in% SummarizedExperiment::assayNames(object)) {
            stop("Assay 'counts' is required for scVI.")
        }

        batch_vec <- as.character(SummarizedExperiment::colData(object)[[batch]])

        # Subset to HVG genes and keep only counts assay
        sce_sub <- object[features, ]
        SummarizedExperiment::assays(sce_sub) <- SummarizedExperiment::assays(sce_sub)["counts"]
        # Strip colData / rowData columns to avoid h5py type errors on write
        SummarizedExperiment::colData(sce_sub) <- SummarizedExperiment::colData(sce_sub)[, 0, drop = FALSE]
        SummarizedExperiment::rowData(sce_sub) <- SummarizedExperiment::rowData(sce_sub)[, 0, drop = FALSE]
        # Also strip metadata (sclet state records may contain datetime objects)
        S4Vectors::metadata(sce_sub) <- list()

        # Write to H5AD to bypass dgCMatrix -> scipy conversion issues
        tmp_h5ad <- tempfile(fileext = ".h5ad")
        on.exit(unlink(tmp_h5ad), add = TRUE)
        message("Writing H5AD for scVI...")
        zellkonverter::writeH5AD(sce_sub, tmp_h5ad)

        message("Running scVI via basilisk...")
        message("(First run will take several minutes to set up the Python environment)")

        # The CPU environment is the default; the GPU one is only used when the
        # caller asks for it and it can actually see a card. Both run the same
        # callback -- scvi-tools takes the accelerator as an argument, so the
        # only difference is which interpreter runs it.
        use_gpu <- identical(sclet_ai_device(device), "cuda")
        if (use_gpu && !sclet_gpu_available()) {
            sclet_ai_error(
                "sclet_gpu_unavailable",
                paste(
                    "device = \"cuda\" was requested but the GPU environment cannot see a device.",
                    "Run sclet_gpu_status() to see what it found.",
                    "Build the environment once with sclet_env_path(\"sclet_gpu_env\", create = TRUE),",
                    "or use device = \"cpu\"."
                )
            )
        }
        env <- if (use_gpu) sclet_gpu_env else sclet_scvi_env
        accelerator <- if (use_gpu) "gpu" else "cpu"

        latent <- sclet_basilisk_run(
            env,
            function(h5ad_path, batches, seed = NULL, accelerator = "cpu", ...) {
            scvi <- reticulate::import("scvi")
            ad <- reticulate::import("anndata")
            pl <- reticulate::import("pytorch_lightning")

            adata <- ad$read_h5ad(h5ad_path)
            adata$obs[["batch"]] <- batches

            scvi$model$SCVI$setup_anndata(adata, batch_key = "batch")
            # Without a seed the training is random, so the same input returns a
            # different latent space on each call and the route cannot be
            # reproduced. `seed` is NOT a constructor argument (measured on
            # scvi-tools 1.4.3 and 1.5.1: neither declares it, and 1.5 rejects it
            # even through **kwargs), so it is set where the library reads it --
            # lightning's global seed -- before the model is constructed.
            if (!is.null(seed)) {
                pl$seed_everything(as.integer(seed))
            }
            model <- scvi$model$SCVI(adata)
            # The accelerator is a training argument, not a global: a GPU model
            # trained on one device is reproducible the same way a CPU one is.
            model$train(accelerator = accelerator)

            return(model$get_latent_representation())
        },
            h5ad_path = tmp_h5ad,
            batches = batch_vec,
            seed = seed,
            accelerator = accelerator
        )
        
        rownames(latent) <- colnames(object)
        colnames(latent) <- paste0("scVI_", seq_len(ncol(latent)))
        
        SingleCellExperiment::reducedDim(object, "scVI") <- latent
        
        # Register integration state
        object <- sclet_set_analysis_state(
            object = object,
            type = "integration",
            id = name,
            method = "scvi-tools",
            inputs = list(
                layer = "counts",
                batch = batch
            ),
            artifacts = list(
                reduction = "scVI"
            ),
            params = list(
                seed = seed
            ),
            summary = list(
                batch_var = batch,
                n_components = ncol(latent)
            )
        )
        
        object <- sclet_set_analysis_state(
            object = object,
            type = "reduction",
            id = "scvi",
            method = "scvi-tools",
            inputs = list(
                layer = "counts",
                integration = name
            ),
            artifacts = list(
                reduction = "scVI"
            ),
            summary = list(
                n_components = ncol(latent)
            )
        )
        
        object <- sclet_set_active_reduction(object, "scVI")
        
        object <- sclet_log_command(
            object,
            "RunIntegration",
            params = list(method = method, batch = batch, seed = seed),
            outputs = list(integration = name, reduction = "scVI")
        )
    }
    
    return(object)
}
