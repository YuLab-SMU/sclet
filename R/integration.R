#' Run Integration
#' 
#' A universal wrapper for batch correction and integration methods.
#' 
#' @title RunIntegration
#' @param object A SingleCellExperiment object
#' @param method Integration method to use. Currently supports "fastMNN", "Harmony", and "scVI".
#' @param batch The column name in colData representing the batch variable.
#' @param features Features to use for integration. If NULL, automatically uses VariableFeatures.
#' @param layer Layer to use for expression-based integration (e.g. fastMNN). If NULL, uses DefaultLayer.
#' @param reduction Reduction to use for reduction-based integration (e.g. Harmony). If NULL, uses DefaultReduction.
#' @param name Integration record id. Defaults to the method name.
#' @param ... Additional arguments passed to the underlying integration functions.
#' @return A SingleCellExperiment object with integration results recorded.
#' @export
RunIntegration <- function(object, method = c("fastMNN", "Harmony", "scVI"), batch = NULL, features = NULL, layer = NULL, reduction = NULL, name = NULL, ...) {
    method <- match.arg(method)
    
    if (is.null(batch)) {
        stop("Please specify 'batch' column.")
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
            stop(paste("Reduction", reduction, "not found. Please run PCA first."))
        }
        
        reduction_state <- sclet_get_state_record(object, "reduction", tolower(reduction))
        
        # harmony::RunHarmony modifies the object and adds a new reduction
        # RunHarmony supports SCE out of the box if Seurat is not loaded, but to be safe we can use the matrix
        # Or just pass the SCE directly
        object <- harmony::RunHarmony(
            object, 
            group.by.vars = batch, 
            reduction = reduction, 
            reduction.save = "HARMONY", 
            ...
        )
        
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
            params = list(method = method, batch = batch, reduction = reduction),
            outputs = list(integration = name, reduction = "HARMONY")
        )
    } else if (method == "scVI") {
        if (!requireNamespace("basilisk", quietly = TRUE)) {
            stop("Package 'basilisk' is required for scVI.")
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
        
        counts_mat <- t(as.matrix(SummarizedExperiment::assay(object, "counts")[features, ]))
        batch_vec <- as.character(SummarizedExperiment::colData(object)[[batch]])
        
        message("Running scVI via basilisk...")
        
        latent <- basilisk::basiliskRun(env = sclet_env, fun = function(counts, batches, ...) {
            scvi <- reticulate::import("scvi")
            ad <- reticulate::import("anndata")
            pd <- reticulate::import("pandas")
            
            obs <- pd$DataFrame(list(batch = batches), index = rownames(counts))
            adata <- ad$AnnData(X = counts, obs = obs)
            
            scvi$model$SCVI$setup_anndata(adata, batch_key = "batch")
            model <- scvi$model$SCVI(adata)
            model$train(...)
            
            return(model$get_latent_representation())
        }, counts = counts_mat, batches = batch_vec, ...)
        
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
            params = list(method = method, batch = batch),
            outputs = list(integration = name, reduction = "scVI")
        )
    }
    
    return(object)
}
